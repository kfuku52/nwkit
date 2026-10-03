"""Experimental conditional allopolyploid MSC CLI with a separate likelihood schema."""

import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor

import pandas as pd
from ete4.parser.newick import PARSERS

from nwkit import __version__
from nwkit.mul_msc_model import dated_candidates, score_gene, validate_sampling
from nwkit.output_transaction import output_transaction, validate_output_targets
from nwkit.species_parser import get_species_parser
from nwkit.util import (
    _serialize_newick_node_name,
    read_tree,
    read_trees,
    validate_outputs_do_not_replace_inputs,
)


def dated_text(tree):
    # Match RADTE's dated writer: ETE's default six digits lose shared ages.
    parser = {
        key: [dict(item) for item in fields] for key, fields in PARSERS[1].items()
    }
    for fields in parser.values():
        fields[1]["write"] = lambda value: format(float(value), ".17g")
    parser["leaf"][0]["write"] = lambda value: _serialize_newick_node_name(
        value, is_internal=False
    )
    parser["internal"][0]["write"] = lambda value: _serialize_newick_node_name(
        value, is_internal=True
    )
    return "[&R]" + tree.write(parser=parser, format_root_node=True, props=[])


def _time_scale(args, *, require_age=True):
    if args.species_time_unit is None or (
        require_age and args.hybridization_age is None
    ):
        raise ValueError("MSC requires --species-time-unit and --hybridization-age.")
    if args.species_time_unit == "coalescent":
        if args.effective_population_size is not None:
            raise ValueError(
                "Coalescent units already include Ne; omit --effective-population-size."
            )
        return 1.0
    size = args.effective_population_size
    if (
        size is None
        or not math.isfinite(size)
        or size <= 0
        or not math.isfinite(2 * size)
    ):
        raise ValueError(
            "Generations require a finite positive --effective-population-size."
        )
    return 2 * size


def _score_candidate(task):
    candidate, genes, parser, scale, max_states, max_assignments, collect = task
    rows, values = [], []
    for number, gene in enumerate(genes, 1):
        likelihood, count, states = score_gene(
            gene,
            candidate,
            parser,
            time_scale=scale,
            max_states=max_states,
            max_assignments=max_assignments,
        )
        values.append(likelihood)
        if not collect:
            continue
        rows.append(
            {
                "mul.tree": candidate.id,
                "gene.tree": number,
                "log_likelihood": likelihood,
                "num_assignments": count,
                "coalescent_states": states,
            }
        )
    total = math.fsum(values)
    if not math.isfinite(total):
        raise ArithmeticError("MSC total log likelihood must be finite.")
    return candidate.id, total, rows


def _search(candidates, genes, parser, args, scale):
    valid = [candidate for candidate in candidates if candidate.status == "evaluated"]
    collect = any(
        getattr(args, name) is not None for name in ("report", "check_out", "model_out")
    )
    tasks = (
        (
            candidate,
            genes,
            parser,
            scale,
            args.max_coalescent_states,
            args.max_coalescent_assignments,
            collect,
        )
        for candidate in valid
    )
    if args.cpus == 1:
        return tuple(map(_score_candidate, tasks))
    with ProcessPoolExecutor(max_workers=min(args.cpus, len(valid))) as pool:
        return tuple(pool.map(_score_candidate, tasks, chunksize=1))


def mul_msc_main(args):
    if getattr(args, "msc_fit", None) not in (None, "fixed"):
        from nwkit.mul_msc_fit_cli import mul_msc_fit_main

        return mul_msc_fit_main(args)
    if any(
        getattr(args, name, None) is not None
        for name in (
            "hybridization_age_bounds",
            "population_size_bounds",
            "msc_grid_points",
            "msc_fit_starts",
            "msc_maxiter",
            "msc_max_evaluations",
            "msc_profile_out",
        )
    ):
        raise ValueError("MSC fit options require --msc-fit age, ne or joint.")
    scale = _time_scale(args)
    species, genes, parser, paths = read_msc_inputs(args)
    candidates, polyploid_species = dated_candidates(
        species,
        args.h1,
        args.h2,
        args.hybridization_age,
        max_candidates=args.max_candidates,
    )
    validate_sampling(genes, parser, polyploid_species, species.leaf_names())
    results = _search(candidates, genes, parser, args, scale)
    write_fixed_results(
        args, candidates, polyploid_species, genes, results, scale, paths
    )


def read_msc_inputs(args, *, extra_outputs=()):
    for name, default in (
        ("max_coalescent_states", 100000),
        ("max_coalescent_assignments", 10000),
    ):
        if getattr(args, name) is None:
            setattr(args, name, default)
    for name in (
        "cpus",
        "max_candidates",
        "max_coalescent_states",
        "max_coalescent_assignments",
    ):
        if getattr(args, name) < 1:
            raise ValueError(f"--{name.replace('_', '-')} must be positive.")
    if args.multree != "no":
        raise ValueError(
            "MSC prototype requires H1/H2 search, not arbitrary supplied MUL-trees."
        )
    if args.max_state_pairs != 10000000 or args.max_maps != 100000:
        raise ValueError(
            "--max-state-pairs/--max-maps are DL-only; use MSC work/assignment limits."
        )
    outputs = {
        name: getattr(args, name)
        for name in (
            "outfile",
            "report",
            "check_out",
            "tree_out",
            "model_out",
            *extra_outputs,
        )
        if getattr(args, name) is not None
    }
    if any(
        not path or (path == "-" and name != "outfile")
        for name, path in outputs.items()
    ):
        raise ValueError("Only the primary MSC likelihood table may use stdout.")
    paths = {name: path for name, path in outputs.items() if path != "-"}
    validate_output_targets(paths.values())
    validate_outputs_do_not_replace_inputs(
        [
            ("--infile", args.infile),
            ("--species-tree", args.species_tree),
            ("--species-map-tsv", args.species_map_tsv),
        ],
        list(paths.items()),
    )
    species = read_tree(
        args.species_tree,
        "auto",
        args.quoted_node_names,
        rooted=args.species_tree_rooted,
    )
    if species.dist not in (None, 0):
        raise ValueError(
            "MSC requires an omitted/zero species-root stem; the ancestral population is infinite."
        )
    genes = read_trees(
        args.infile,
        args.format,
        args.quoted_node_names,
        rooted=args.input_rooted,
        quiet=True,
    )
    if not genes:
        raise ValueError("At least one gene tree is required.")
    parser = get_species_parser(args)
    return species, genes, parser, paths


def write_fixed_results(
    args, candidates, polyploid_species, genes, results, scale, paths
):
    likelihoods = {number: value for number, value, _ in results}
    valid = sorted(
        (candidate for candidate in candidates if candidate.status == "evaluated"),
        key=lambda candidate: (-likelihoods[candidate.id], candidate.id),
    )
    best = valid[0]
    rows = []
    for candidate in (
        *valid,
        *(candidate for candidate in candidates if candidate.status == "excluded"),
    ):
        value = likelihoods.get(candidate.id)
        rows.append(
            {
                "mul.tree": candidate.id,
                "h1.node": candidate.h1,
                "h2.node": candidate.h2,
                "log_likelihood": value,
                "delta_log_likelihood": None
                if value is None
                else likelihoods[best.id] - value,
                "dated.tree": dated_text(candidate.tree) if value is not None else None,
                "status": candidate.status,
                "reason": candidate.reason,
            }
        )
    details = [row for _, _, genes_rows in results for row in genes_rows]
    metadata = {
        "schema_version": 1,
        "method": "conditional-direct-parent-MUL-MSC-topology-v1",
        "nwkit_version": __version__,
        "score_model": "msc",
        "score_schema": "nwkit-mul-msc-likelihood-v1",
        "objective": "sum of gene-topology log probabilities, conditional on observed copies and fixed population parameters",
        "num_gene_trees": len(genes),
        "polyploid_species": list(polyploid_species),
        "species_time_unit": args.species_time_unit,
        "effective_population_size": args.effective_population_size,
        "time_scale": scale,
        "hybridization_age": args.hybridization_age,
        "estimated_parameters": [],
        "best_hypotheses": [
            candidate.id
            for candidate in valid
            if likelihoods[candidate.id] == likelihoods[best.id]
        ],
        "reported_hypothesis": best.id,
        "assignment_prior": "uniform over injective species-preserving assignments, normalized separately within each family",
        "gene_filtering": "none; every evaluated hypothesis uses every input gene tree",
        "root_population": "infinite ancestral population with the same fixed Ne scale",
        "numerics": "ancestral configurations and compatible-history counts; log-space probabilities; uniformization/long-time relative remainder bounds 1e-15",
        "limitations": [
            "Conditional parent-candidate comparison assumes one disomic allotetraploid event; not a test of WGD or its mode.",
            "H1 first-parent stem is supplied; the second parent directly donates at the supplied hybridization age, with no ghost-parent divergence parameter.",
            "One shared fixed Ne; duplicated descendant populations share ages and the population scale.",
            "At most one sampled copy per diploid/subgenome; no SSD, allelic replicates, gene loss or detection/ascertainment model.",
            "Missing copies are conditioned on, not scored as losses; observed homoeolog identity is latent.",
            "Fixed rooted binary gene topologies; no sequence likelihood, gene-tree/root uncertainty, subgenome exchange, or linked-family model.",
            "No parameter optimization, calibrated support, posterior candidate probabilities, or inferred hybridization time.",
        ],
        "scores": rows,
        "gene_likelihoods": details,
    }
    tables = {
        "outfile": pd.DataFrame(rows),
        "report": pd.DataFrame(details),
        "check_out": pd.DataFrame(details),
    }
    with output_transaction(paths.values()) as staged:
        for name, path in paths.items():
            if name in tables:
                staged.write_text(
                    path,
                    lambda handle, name=name: tables[name].to_csv(
                        handle, sep="\t", index=False
                    ),
                )
            elif name == "tree_out":
                staged.write_text(
                    path, lambda handle: handle.write(dated_text(best.tree) + "\n")
                )
            else:
                staged.write_text(
                    path,
                    lambda handle: json.dump(
                        metadata, handle, indent=2, allow_nan=False
                    ),
                )
    if args.outfile == "-":
        tables["outfile"].to_csv(sys.stdout, sep="\t", index=False)
    sys.stderr.write(
        f"Conditional MUL MSC: {len(genes)} gene trees, {len(valid)} evaluated parent candidates; best={best.id}, log_likelihood={likelihoods[best.id]}. Not a WGD test.\n"
    )
