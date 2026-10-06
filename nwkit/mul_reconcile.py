"""MUL-tree search with separate D+L, conditional MSC and locus MC contracts."""

import csv
import json
import re
import sys
from concurrent.futures import ProcessPoolExecutor
from io import StringIO
from types import SimpleNamespace
from urllib.parse import quote

import pandas as pd

from nwkit import __version__
from nwkit.mul_reconcile_model import Reconciliation, hypotheses, validate_binary
from nwkit.output_transaction import output_transaction, validate_output_targets
from nwkit.species_parser import get_species_parser
from nwkit.util import (
    _serialize_newick_node_name,
    copy_tree_iteratively,
    read_tree,
    read_trees,
    validate_outputs_do_not_replace_inputs,
    write_tree,
)


def topology_text(tree, *, internal_names=True):
    output = StringIO()
    write_tree(
        tree,
        SimpleNamespace(outfile=output),
        8 if internal_names else 9,
        quiet=True,
        props=[],
    )
    return output.getvalue()


def annotated_mapping_text(gene, mapping):
    """Retain GRAMPA's node[species-map-duplication] map serialization."""
    labeled = copy_tree_iteratively(gene)
    replacements = {}
    internal_number = 0
    for index, node in enumerate(labeled.traverse("postorder")):
        if node.is_leaf:
            label = node.name
        else:
            internal_number += 1
            label = f"<{internal_number}>"
        mapped = mapping[index]["mul_label"]
        if not mapped.endswith(("+", "*")):
            mapped += "+"
        mapped = quote(mapped, safe="._-+*<>")
        name = _serialize_newick_node_name(
            label, is_internal=not node.is_leaf, quote_style="auto"
        )
        replacements[str(index)] = f"{name}[{mapped}-{mapping[index]['duplication']}]"
        node.name = f"NWKITMULMAP{index}"
    # Replace writer-generated tokens, never user labels or Newick delimiters.
    return re.sub(
        r"NWKITMULMAP(\d+)",
        lambda match: replacements[match[1]],
        topology_text(labeled),
    ).rstrip(";")


def _score_task(task):
    candidate, genes, parser, limit, collect_checks = task
    checks, score = [], 0
    for number, gene in enumerate(genes, 1):
        result = Reconciliation(gene, candidate.tree, parser, limit)
        score += result.score
        if not collect_checks:
            continue
        checks.append(
            {
                "mul.tree": candidate.id,
                "gene.tree": number,
                "groups": result.ambiguous_tips,
                "fixed": 0,
                "combinations": str(result.combinations),
                "over.cap.filtered": "N",
                "optimal.mappings": str(result.num_maps),
                "state.pairs": result.state_pairs,
            }
        )
    return candidate.id, score, checks


def run_search(
    candidates, genes, parser, *, cpus=1, max_state_pairs=10000000, collect_checks=True
):
    tasks = (
        (candidate, genes, parser, max_state_pairs, collect_checks)
        for candidate in candidates
    )
    if cpus == 1:
        return tuple(map(_score_task, tasks))
    with ProcessPoolExecutor(max_workers=min(cpus, len(candidates))) as pool:
        return tuple(pool.map(_score_task, tasks, chunksize=1))


def mul_reconcile_main(args):
    if (
        getattr(args, "node_out", None) is not None
        and getattr(args, "score_model", "dl") != "dl"
    ):
        raise ValueError("--node-out requires --score-model dl.")
    if getattr(args, "score_model", "dl") == "locus-mc":
        from nwkit.mul_locus_cli import mul_locus_main

        return mul_locus_main(args)
    if any(
        getattr(args, name, None) is not None
        for name in (
            "locus_model",
            "locus_bootstrap",
            "locus_calibration_out",
            "locus_null_calibration",
        )
    ):
        raise ValueError("Locus options require --score-model locus-mc.")
    if getattr(args, "score_model", "dl") == "msc":
        from nwkit.mul_msc import mul_msc_main

        return mul_msc_main(args)
    if any(
        getattr(args, name, None) is not None
        for name in (
            "species_time_unit",
            "effective_population_size",
            "hybridization_age",
            "max_coalescent_states",
            "max_coalescent_assignments",
            "msc_fit",
            "hybridization_age_bounds",
            "population_size_bounds",
            "msc_grid_points",
            "msc_fit_starts",
            "msc_maxiter",
            "msc_max_evaluations",
            "msc_profile_out",
        )
    ):
        raise ValueError("MSC time/population options require --score-model msc.")
    for name in ("cpus", "max_candidates", "max_state_pairs", "max_maps"):
        if getattr(args, name) < 1:
            raise ValueError(f"--{name.replace('_', '-')} must be positive.")
    outputs = {
        name: getattr(args, name)
        for name in (
            "outfile",
            "report",
            "check_out",
            "tree_out",
            "model_out",
            "node_out",
        )
        if getattr(args, name, None) is not None
    }
    if any(
        not value or (value == "-" and name != "outfile")
        for name, value in outputs.items()
    ):
        raise ValueError(
            "Only the primary score table may use stdout; paths must be nonempty."
        )
    paths = {name: value for name, value in outputs.items() if value != "-"}
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
    genes = read_trees(
        args.infile,
        args.format,
        args.quoted_node_names,
        rooted=args.input_rooted,
        quiet=True,
    )
    if not genes:
        raise ValueError("At least one gene tree is required.")
    for gene in genes:
        validate_binary(gene, "Gene tree")
    parser = get_species_parser(args)
    candidates = hypotheses(
        species,
        args.h1,
        args.h2,
        multree=args.multree == "yes",
        max_candidates=args.max_candidates,
    )
    scored = run_search(
        candidates,
        genes,
        parser,
        cpus=args.cpus,
        max_state_pairs=args.max_state_pairs,
        collect_checks=args.check_out is not None,
    )
    scores = {number: value for number, value, _ in scored}
    ordered = sorted(
        candidates, key=lambda candidate: (scores[candidate.id], candidate.id)
    )
    best = ordered[0]
    rows = [
        {
            "mul.tree": candidate.id,
            "h1.node": candidate.h1,
            "h2.node": candidate.h2,
            "score": scores[candidate.id],
            "labeled.tree": topology_text(candidate.tree),
            "hypothesis.kind": candidate.kind,
        }
        for candidate in ordered
    ]

    def write_details(handle):
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "mul.tree",
                "gene.tree",
                "dups",
                "losses",
                "total.score",
                "maps",
                "node.maps",
                "maps.format",
                "node.maps.format",
            ],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for number, gene in enumerate(genes, 1):
            result = Reconciliation(gene, best.tree, parser, args.max_state_pairs)
            for duplications, losses, mapping in result.mappings(args.max_maps):
                writer.writerow(
                    {
                        "mul.tree": best.id,
                        "gene.tree": number,
                        "dups": duplications,
                        "losses": losses,
                        "total.score": duplications + losses,
                        "maps": annotated_mapping_text(gene, mapping),
                        "node.maps": json.dumps(mapping, separators=(",", ":")),
                        "maps.format": "grampa-annotated-newick",
                        "node.maps.format": "nwkit-node-map-json-v1",
                    }
                )

    metadata = {
        "schema_version": 1,
        "method": "exact-MUL-LCA-DL-parsimony-v1",
        "nwkit_version": __version__,
        "num_gene_trees": len(genes),
        "best_hypotheses": [
            candidate.id
            for candidate in ordered
            if scores[candidate.id] == scores[best.id]
        ],
        "reported_hypothesis": best.id,
        "objective": "duplications + losses, including species-root-to-gene-root losses",
        "gene_filtering": "none; every hypothesis uses every input gene tree",
        "internal_numbering": "one-based closing-parenthesis/postorder, independent of input names",
        "mapping_format": "grampa-annotated-newick; all optimal leaf-occurrence assignments for reported hypothesis",
        "node_mapping_format": "nwkit-node-map-json-v1 in detailed node.maps column",
        "limitations": [
            "Parsimony ranking is not posterior probability or a calibrated polyploidy test.",
            "Fixed rooted binary trees; no ILS, HGT, rooting uncertainty, or sequence likelihood.",
            "Search represents one duplicated clade and one second-parent placement; supplied MUL-trees may be more general.",
        ],
        "scores": rows,
    }
    if getattr(args, "node_out", None) is not None:
        metadata["node_diagnostics"] = {
            "schema": "nwkit-mul-node-assignments-v1",
            "scope": "all optimal mappings for every tied best global D+L hypothesis",
            "identity": "rooted topology and descendant-tip clades; lengths/support ignored",
            "meaning": "enumerated assignments, not probabilities or WGD/SSD origin labels",
        }
    tables = {
        "outfile": pd.DataFrame(rows),
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
            elif name == "report":
                staged.write_text(path, write_details)
            elif name == "check_out":

                def write_checks(handle):
                    writer = csv.DictWriter(
                        handle,
                        fieldnames=[
                            "mul.tree",
                            "gene.tree",
                            "groups",
                            "fixed",
                            "combinations",
                            "over.cap.filtered",
                            "optimal.mappings",
                            "state.pairs",
                        ],
                        delimiter="\t",
                        lineterminator="\n",
                    )
                    writer.writeheader()
                    for _, _, checks in scored:
                        writer.writerows(checks)

                staged.write_text(path, write_checks)
            elif name == "tree_out":
                staged.write_text(
                    path, lambda handle: handle.write(topology_text(best.tree) + "\n")
                )
            elif name == "node_out":
                from nwkit.mul_reconcile_nodes import write_node_diagnostics

                staged.write_text(
                    path,
                    lambda handle: write_node_diagnostics(
                        handle,
                        [c for c in ordered if scores[c.id] == scores[best.id]],
                        genes,
                        parser,
                        species,
                        max_state_pairs=args.max_state_pairs,
                        max_maps=args.max_maps,
                    ),
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
        f"MUL reconciliation: {len(genes)} gene trees, {len(candidates)} hypotheses; best={best.id}, score={scores[best.id]}.\n"
    )
