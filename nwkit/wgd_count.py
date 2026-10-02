"""Validated count-model inputs and transactional WGM candidate outputs."""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from nwkit import __version__
from nwkit.clade_index import CladeIndex
from nwkit.gene_family_input import read_gene_family_table
from nwkit.output_transaction import output_transaction, validate_output_targets
from nwkit.rooting_state import require_rooted
from nwkit.util import (
    assign_branch_ids,
    is_missing_table_value,
    read_tree,
    validate_outputs_do_not_replace_inputs,
    validate_unique_named_leaves,
)
from nwkit.wgd_count_fit import _burst_model, calibrate_scan, scan_counts
from nwkit.wgd_count_model import CountLikelihood, CountTree, rate_categories


def read_counts(path, names, missing_values):
    header, rows = read_gene_family_table(path, ["family_id", *names])
    if set(header) != {"family_id", *names}:
        raise ValueError(
            "Count columns must exactly match the species-tree tips plus family_id."
        )
    families = [row["family_id"] for row in rows]
    if any(not family for family in families) or len(set(families)) != len(families):
        raise ValueError("family_id values must be nonempty and unique.")
    counts = np.empty((len(rows), len(names)))
    for i, row in enumerate(rows):
        for j, name in enumerate(names):
            value = row[name]
            if is_missing_table_value(value, missing_values):
                counts[i, j] = np.nan
                continue
            try:
                number = float(value)
            except ValueError as exc:
                raise ValueError(
                    f"{families[i]}/{name}: count must be a nonnegative integer or missing."
                ) from exc
            if not np.isfinite(number) or number < 0 or number != np.floor(number):
                raise ValueError(
                    f"{families[i]}/{name}: count must be a nonnegative integer or missing."
                )
            counts[i, j] = number
    return families, counts


def _count_tree(tree):
    require_rooted(tree, "WGM count inference requires a rooted species tree.")
    validate_unique_named_leaves(tree, "--infile")
    nodes = tuple(tree.traverse("preorder"))
    if any(len(node.children) == 1 for node in nodes):
        raise ValueError(
            "WGM count inference does not accept unary species-tree nodes."
        )
    indices = {node: index for index, node in enumerate(nodes)}
    identifiers = assign_branch_ids(tree)
    clades = CladeIndex(tree)
    names = tuple(sorted(leaf.name for leaf in tree.leaves()))
    named = {node.name: node for node in tree.leaves()}
    if any(node.dist is None for node in nodes[1:]):
        raise ValueError("WGM count inference requires every non-root branch length.")
    design = CountTree(
        tuple(-1 if node.is_root else indices[node.up] for node in nodes),
        tuple(0.0 if node.is_root else float(node.dist) for node in nodes),
        tuple(indices[named[name]] for name in names),
        names,
        tuple(identifiers[node] for node in nodes),
        tuple(clades.clade_id_for_node(node) for node in nodes),
    )
    return design, nodes, clades


def _detection(path, names):
    if not path:
        return np.ones(len(names))
    _, rows = read_gene_family_table(path, ["leaf_name", "detection_probability"])
    keys = [row["leaf_name"] for row in rows]
    if len(set(keys)) != len(keys) or set(keys) != set(names):
        raise ValueError(
            "Detection TSV must exactly and uniquely cover species-tree tips."
        )
    by_name = {row["leaf_name"]: float(row["detection_probability"]) for row in rows}
    return np.array([by_name[name] for name in names])


def _groups(args, design):
    if args.rate_groups_tsv:
        _, rows = read_gene_family_table(args.rate_groups_tsv, ["branch_id", "regime"])
        mapping = {}
        for row in rows:
            branch = int(row["branch_id"])
            if branch in mapping or not row["regime"]:
                raise ValueError(
                    "Rate-group rows need unique branch IDs and nonempty regimes."
                )
            mapping[branch] = row["regime"]
        if set(mapping) != set(design.branch_ids[1:]):
            raise ValueError(
                "Rate-group TSV must exactly cover non-root input branch IDs."
            )
        labels = ["root"] + [mapping[branch] for branch in design.branch_ids[1:]]
    elif args.rate_model == "homogeneous":
        labels = ["root"] + ["background"] * (len(design.parents) - 1)
    else:
        tips = set(design.tip_nodes)
        labels = ["root"] + [
            "terminal" if node in tips else "internal"
            for node in range(1, len(design.parents))
        ]
    unique = sorted(set(labels[1:]))
    indices = {label: index for index, label in enumerate(unique)}
    return tuple(
        0 if node == 0 else indices[label] for node, label in enumerate(labels)
    ), unique


def _fit_json(fit):
    return {
        "rates": fit.rates.tolist(),
        "rate_column_order": ["duplication", "loss"],
        "root_mean": fit.root_mean,
        "log_likelihood": fit.log_likelihood,
        "num_parameters": fit.num_parameters,
        "aic": fit.aic,
        "max_count": fit.max_count,
        "max_family_log_likelihood_change_on_doubling": fit.state_error,
        "state_converged": fit.converged,
        "nuisance_bound_reached": fit.boundary,
        "optimizer_message": fit.message,
        "event": None
        if fit.event is None
        else {
            "node": fit.event.node,
            "retention": fit.event.retention,
            "fraction": fit.event.fraction,
            "multiplicity": fit.event.multiplicity,
        },
    }


def _paths(args):
    inputs = [
        ("--infile", args.infile),
        ("--counts", args.counts),
        ("--detection-tsv", args.detection_tsv),
        ("--rate-groups-tsv", args.rate_groups_tsv),
    ]
    paths = {}
    for role in ("outfile", "model_out"):
        path = getattr(args, role)
        if path == "" or (path == "-" and role != "outfile"):
            raise ValueError(
                f"--{role.replace('_', '-')} requires a nonempty file path."
            )
        if path not in (None, "-"):
            paths[role] = path
    validate_outputs_do_not_replace_inputs(
        inputs, [(role, path) for role, path in paths.items()]
    )
    validate_output_targets(paths.values())
    return paths


def wgd_count_main(args):
    paths = _paths(args)
    if (
        args.bootstrap < 0
        or args.seed < 0
        or args.max_iterations < 1
        or args.max_states < 16
    ):
        raise ValueError(
            "Bootstrap/seed must be nonnegative; iterations positive; max-states >= 16."
        )
    if not np.isfinite(args.alpha) or not 0 < args.alpha < 1:
        raise ValueError("--alpha must be in (0, 1).")
    if args.multiplicity < 2:
        raise ValueError("--multiplicity must be >= 2.")
    fractions = tuple(float(value) for value in args.event_fractions.split(","))
    shape = (
        None if args.family_gamma_shape == "none" else float(args.family_gamma_shape)
    )
    scales = rate_categories(shape, args.family_rate_categories)
    tree = read_tree(
        args.infile, args.format, args.quoted_node_names, rooted=args.input_rooted
    )
    design, nodes, clades = _count_tree(tree)
    families, counts = read_counts(args.counts, design.tip_names, args.missing_values)
    groups, group_names = _groups(args, design)
    model = CountLikelihood(
        design,
        counts,
        detection=_detection(args.detection_tsv, design.tip_names),
        rate_scales=scales,
        branch_groups=groups,
        ascertainment=args.ascertainment,
    )
    selected = None
    if args.candidate_branches:
        requested = [int(value) for value in args.candidate_branches.split(",")]
        if len(set(requested)) != len(requested) or any(
            branch not in design.branch_ids[1:] for branch in requested
        ):
            raise ValueError(
                "--candidate-branches needs unique non-root input branch IDs."
            )
        selected = tuple(design.branch_ids.index(branch) for branch in requested)
    options = dict(
        fractions=fractions,
        multiplicity=args.multiplicity,
        max_states=args.max_states,
        state_tolerance=args.state_tolerance,
        max_iterations=args.max_iterations,
    )
    result = scan_counts(model, nodes=selected, **options)
    if args.bootstrap:

        def progress(done, total):
            print(f"WGM count bootstrap: {done}/{total}", file=sys.stderr)

        result = calibrate_scan(
            model, result, args.bootstrap, args.seed, progress=progress, **options
        )
    rows = []
    for rank, candidate in enumerate(result.candidates, 1):
        fit = candidate.event_fit
        event = fit.event
        if event is None:
            raise ValueError("WGM scan returned a candidate without an event.")
        if candidate.p_value is None:
            status = "not_calibrated"
        elif candidate.p_value > args.alpha:
            status = "background_compatible"
        elif candidate.burst_aic_difference <= 0:
            status = "branch_burst_preferred"
        else:
            status = "count_supported_conditional"
        rows.append(
            {
                "rank": rank,
                "event_id": f"wgm:{design.clade_ids[candidate.node]}:m{args.multiplicity}",
                "branch_id": design.branch_ids[candidate.node],
                "species_event_id": design.clade_ids[candidate.node],
                "parent_species_event_id": design.clade_ids[
                    design.parents[candidate.node]
                ],
                "descendant_taxa": clades.csv_for_node(nodes[candidate.node]),
                "branch_length": design.lengths[candidate.node],
                "multiplicity": args.multiplicity,
                "retention": event.retention,
                "event_fraction": event.fraction,
                "log_likelihood": fit.log_likelihood,
                "background_log_likelihood": result.background.log_likelihood,
                "branch_burst_log_likelihood": candidate.burst_fit.log_likelihood,
                "likelihood_ratio_statistic": candidate.improvement,
                "burst_minus_event_aic": candidate.burst_aic_difference,
                "p_value": candidate.p_value,
                "p_value_mc_se": candidate.p_value_mc_se,
                "p_value_method": result.calibration,
                "num_bootstrap": args.bootstrap,
                "state_error": fit.state_error,
                "nuisance_bound_reached": fit.boundary,
                "background_nuisance_bound_reached": result.background.boundary,
                "branch_burst_nuisance_bound_reached": candidate.burst_fit.boundary,
                "num_families": len(families),
                "count_support": status,
            }
        )
    table = pd.DataFrame(rows).to_csv(
        sep="\t", index=False, float_format="%.12g", na_rep="NA"
    )
    metadata = {
        "schema_version": 1,
        "method": "native-linear-DL-retained-WGM-count-v1",
        "nwkit_version": __version__,
        "tip_names": list(design.tip_names),
        "num_families": len(families),
        "num_missing_cells": int(np.isnan(counts).sum()),
        "detection_probabilities": model.detection.tolist(),
        "family_rate_scales": model.rate_scales.tolist(),
        "family_gamma_shape": shape,
        "branch_rate_groups": list(groups),
        "branch_rate_group_names": group_names,
        "branch_ids": list(design.branch_ids),
        "species_event_ids": list(design.clade_ids),
        "species_tree": {
            "parents": list(design.parents),
            "lengths": list(design.lengths),
        },
        "event_fractions": list(fractions),
        "multiplicity": args.multiplicity,
        "background": _fit_json(result.background),
        "candidates": [
            {
                "branch_id": design.branch_ids[c.node],
                "event": _fit_json(c.event_fit),
                "branch_burst": _fit_json(c.burst_fit),
                "branch_burst_rate_groups": list(
                    _burst_model(model, c.node).branch_groups
                ),
            }
            for c in result.candidates
        ],
        "calibration": result.calibration,
        "bootstrap_max_statistics": list(result.bootstrap_statistics),
        "seed": args.seed,
        "alpha": args.alpha,
        "ascertainment": f"{args.ascertainment}; conditional on observation, separately for each missingness mask, after mixing family rate categories",
        "root_prior": "positive geometric with estimated mean; omitted tail is not renormalized",
        "uncertainty": "Search-wide plug-in bootstrap conditional on this species tree, observation process and finite candidate grid; not posterior probabilities or guaranteed composite-null error control.",
        "limitations": [
            "Families are present at the species root; later family origins and horizontal transfer are not modeled.",
            "One genome event is searched at a time; simultaneous events, allopolyploid donor trees, and gene-tree-node attribution are not inferred.",
            "Gamma quantile categories approximate family rate heterogeneity; gamma shape and detection probabilities are supplied, not estimated.",
            "Branch-burst AIC is a conditional diagnostic, not a calibrated test against arbitrary branch-specific SSD rates.",
            "Event fraction is selected on a grid, not a dated event estimate or confidence interval.",
            "Rates have numerical bounds exp(-16)..20 per median positive branch length; root mean upper bound is max(10, twice maximum observed count). Inspect nuisance_bound_reached.",
        ],
    }
    model_text = (
        json.dumps(metadata, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    )
    with output_transaction(paths.values()) as staged:
        for role, path in paths.items():
            Path(staged[path]).write_text(
                model_text if role == "model_out" else table, encoding="utf-8"
            )
    if args.outfile == "-":
        sys.stdout.write(table)
    print(
        f"WGM count: {len(families)} families, {len(result.candidates)} candidate branches; {result.calibration}.",
        file=sys.stderr,
    )
