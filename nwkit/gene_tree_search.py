"""Auditable targeted search; candidate detection never labels a gene erroneous."""

import hashlib
import json
import os
import sys
from pathlib import Path

import pandas as pd

from nwkit import __version__
from nwkit.gene_tree_search_generax import (
    evaluate_candidates,
    read_alignment,
    tree_text,
)
from nwkit.gene_tree_search_model import (
    DLContext,
    discover_proposals,
    generate_candidates,
)
from nwkit.output_transaction import output_transaction, validate_output_targets
from nwkit.reconcile import _validate_rooted_binary_tree
from nwkit.species_parser import get_species_parser
from nwkit.util import (
    read_input_text,
    read_tree,
    validate_outputs_do_not_replace_inputs,
)


def _validate_options(args):
    for name in (
        "max_moved_tips",
        "max_proposals",
        "max_set_states",
        "beam_width",
        "max_candidates",
        "max_evaluations",
        "timeout",
        "eval_rounds",
    ):
        if getattr(args, name) < 1:
            raise ValueError(f"--{name.replace('_', '-')} must be positive.")
    if args.evaluation == "generax" and not all(
        (args.alignment, args.subst_model, args.workdir)
    ):
        raise ValueError(
            "--evaluation generax requires --alignment, --subst-model and a new --workdir."
        )
    if args.evaluation == "none" and any(
        (args.alignment, args.subst_model, args.workdir)
    ):
        raise ValueError(
            "Alignment/model/workdir options require --evaluation generax."
        )
    outputs = {
        n: getattr(args, n)
        for n in ("outfile", "tree_out", "candidates_out", "sets_out", "report_out")
        if getattr(args, n) is not None
    }
    if any(
        not value or (value == "-" and name != "outfile")
        for name, value in outputs.items()
    ):
        raise ValueError(
            "Only the primary score table may use stdout; output paths must be nonempty."
        )
    paths = {name: path for name, path in outputs.items() if path != "-"}
    validate_output_targets(paths.values())
    for path in paths.values():
        if not Path(path).parent.is_dir():
            raise ValueError(
                f"Output directory must already exist: {Path(path).parent}"
            )
    inputs = [
        ("--infile", args.infile),
        ("--species-tree", args.species_tree),
        ("--species-map-tsv", args.species_map_tsv),
        ("--alignment", args.alignment),
    ]
    validate_outputs_do_not_replace_inputs(inputs, list(paths.items()))
    if args.workdir:
        workdir = Path(args.workdir).resolve()
        if workdir.exists():
            raise ValueError("--workdir must be a new directory.")
        for label, path in [*inputs, *paths.items()]:
            if path not in (None, "-") and Path(path).resolve().is_relative_to(workdir):
                raise ValueError(f"{label} must be outside --workdir.")
    return paths


def _input_provenance(args):
    records = {}
    for name in ("species_map_tsv", "alignment"):
        path = getattr(args, name)
        if path == "-":
            records[name] = {
                "source": "stdin",
                "hash_status": "see --audit stdin record",
            }
        elif path is not None:
            source = Path(path).resolve()
            records[name] = {
                "path": str(source),
                "sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
            }
    return records


def _tree_source(source):
    # Read once: hashes and parsed trees describe the same snapshot, including
    # inline Newick and stdin. File hashes retain the original bytes/BOM.
    if source != "-" and os.path.isfile(source):
        path = Path(source).resolve()
        data = path.read_bytes()
        text = data.decode("utf-8-sig")
        record = {"path": str(path)}
    else:
        text = read_input_text(source)
        data = text.encode("utf-8")
        record = {"source": "stdin" if source == "-" else "inline"}
    record["sha256"] = hashlib.sha256(data).hexdigest()
    return text, record


def _best_evaluation(evaluations):
    baseline = evaluations["baseline"]
    improved = {
        name: fit
        for name, fit in evaluations.items()
        if fit.joint - fit.rounding_bound > baseline.joint + baseline.rounding_bound
    }
    return (
        max(improved, key=lambda name: improved[name].joint) if improved else "baseline"
    )


def _evaluation_subset(candidates, budget):
    # Balanced across detected sets, including larger coupled sets from the outset.
    selected = [candidates[0]]
    seen = set()
    for candidate in candidates[1:]:
        if candidate.moved_tips not in seen:
            selected.append(candidate)
            seen.add(candidate.moved_tips)
    selected_ids = {c.id for c in selected}
    selected.extend(c for c in candidates if c.id not in selected_ids)
    return selected[:budget]


def _tables(candidates, proposals, evaluations, best_id, context):
    baseline = evaluations.get("baseline")
    rows = []
    for candidate in candidates:
        fitted = evaluations.get(candidate.id)
        fitted_cost = (
            context.cost(read_tree(fitted.optimized_tree, 1, True, quiet=True))
            if fitted
            else None
        )
        rows.append(
            {
                "candidate_id": candidate.id,
                "moved_tips": json.dumps(candidate.moved_tips, ensure_ascii=False),
                "num_moved_tips": len(candidate.moved_tips),
                "num_pruned_components": candidate.components,
                "proposal_sources": ";".join(candidate.reasons),
                "lca_duplications": candidate.cost.duplications,
                "lca_losses": candidate.cost.losses,
                "species_overlap_duplications": candidate.cost.overlaps,
                "fitted_lca_duplications": fitted_cost.duplications
                if fitted_cost
                else None,
                "fitted_lca_losses": fitted_cost.losses if fitted_cost else None,
                "fitted_species_overlap_duplications": fitted_cost.overlaps
                if fitted_cost
                else None,
                "evaluation_status": "evaluated" if fitted else "not_evaluated",
                "sequence_log_likelihood": fitted.sequence_log_likelihood
                if fitted
                else None,
                "reconciliation_log_likelihood": fitted.reconciliation_log_likelihood
                if fitted
                else None,
                "joint_log_likelihood": fitted.joint if fitted else None,
                "likelihood_rounding_bound": fitted.rounding_bound if fitted else None,
                "delta_joint_log_likelihood": fitted.joint - baseline.joint
                if fitted and baseline
                else None,
                "selected": candidate.id == best_id,
            }
        )
    sets = pd.DataFrame(
        [
            {
                "set_id": f"set_{n:06d}",
                "tips": json.dumps(p.tips, ensure_ascii=False),
                "num_tips": len(p.tips),
                "sources": ";".join(sorted(p.reasons)),
                "constraint_clade_ids": ";".join(sorted(p.constraints)),
                "diagnostic_pruning_gain": p.diagnostic_gain,
            }
            for n, p in enumerate(proposals, 1)
        ],
        columns=[
            "set_id",
            "tips",
            "num_tips",
            "sources",
            "constraint_clade_ids",
            "diagnostic_pruning_gain",
        ],
    )
    trees = pd.DataFrame(
        [{"candidate_id": c.id, "newick": tree_text(c.tree)} for c in candidates]
    )
    return {"outfile": pd.DataFrame(rows), "sets_out": sets, "candidates_out": trees}


def gene_tree_search_main(args):
    paths = _validate_options(args)
    input_records = _input_provenance(args)
    gene_input, input_records["infile"] = _tree_source(args.infile)
    species_input, input_records["species_tree"] = _tree_source(args.species_tree)
    gene = read_tree(
        gene_input, args.format, args.quoted_node_names, rooted=args.input_rooted
    )
    species = read_tree(
        species_input,
        "auto",
        args.quoted_node_names,
        rooted=args.species_tree_rooted,
    )
    _validate_rooted_binary_tree(gene, "--infile")
    parser = get_species_parser(args)
    mapping = {str(n.name): parser.parse(n.name).species_label for n in gene.leaves()}
    if any(not value for value in mapping.values()):
        raise ValueError("Every gene tip needs an explicit or parsed species label.")
    context = DLContext(species, mapping)
    if args.evaluation == "generax":
        read_alignment(args.alignment, mapping, subst_model=args.subst_model)
    proposals, detection_coverage = discover_proposals(
        gene,
        context,
        max_moved_tips=args.max_moved_tips,
        max_proposals=args.max_proposals,
        max_set_states=args.max_set_states,
    )
    candidates, search_coverage = generate_candidates(
        gene,
        context,
        proposals,
        beam_width=args.beam_width,
        max_candidates=args.max_candidates,
        rooted=args.evaluation == "none" or args.root_policy == "keep",
    )
    evaluations = {}
    evaluation_metadata = {}
    best_id = "baseline"
    if args.evaluation == "generax":
        selected = _evaluation_subset(candidates, args.max_evaluations)
        sys.stderr.write(
            f"GeneRax EVAL: fitting {len(selected)} complete-tip topologies including baseline.\n"
        )
        evaluations, evaluation_metadata = evaluate_candidates(
            selected,
            species,
            mapping,
            args.alignment,
            args.workdir,
            command=args.generax_command,
            subst_model=args.subst_model,
            rec_model=args.rec_model,
            root_policy=args.root_policy,
            seed=args.seed,
            timeout=args.timeout,
            rounds=args.eval_rounds,
        )
        best_id = _best_evaluation(evaluations)
    best = next(c for c in candidates if c.id == best_id)
    tree = evaluations[best_id].optimized_tree if evaluations else tree_text(gene)
    tables = _tables(candidates, proposals, evaluations, best_id, context)
    metadata = {
        "schema_version": 1,
        "nwkit_version": __version__,
        "command": "gene-tree-search",
        "inputs": input_records,
        "gene_tips": len(mapping),
        "species_tips": len(list(species.leaves())),
        "detection_coverage": detection_coverage,
        "search_coverage": search_coverage,
        "budgets": {
            name: getattr(args, name)
            for name in (
                "max_moved_tips",
                "max_proposals",
                "max_set_states",
                "beam_width",
                "max_candidates",
                "max_evaluations",
            )
        },
        "evaluated_topologies": len(evaluations),
        "unevaluated_topologies": len(candidates) - len(evaluations),
        "evaluation": evaluation_metadata,
        "selected_candidate": best.id,
        "interpretation": "Structural anomaly candidates, not automatically diagnosed misplaced genes. D+L ranks proposals only. The output remains unchanged without likelihood evaluation. Selection is conditional on the model and searched/evaluated subset, not a global optimum or statistical significance test. Partial fits are beam-ranked without requiring single-move improvement; inferred sets above budget are excluded.",
    }
    with output_transaction(paths.values()) as staged:
        for name, path in paths.items():
            if name in tables:
                staged.write_text(
                    path,
                    lambda handle, n=name: tables[n].to_csv(
                        handle, sep="\t", index=False
                    ),
                )
            elif name == "tree_out":
                staged.write_text(path, lambda handle: handle.write(tree + "\n"))
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
        f"Targeted search: {len(proposals)} detected sets, {len(candidates)} retained topologies, best={best_id}.\n"
    )
