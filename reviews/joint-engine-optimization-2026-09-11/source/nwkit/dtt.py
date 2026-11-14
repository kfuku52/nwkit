"""Validated DTT inputs, result tables and transactional publication."""

import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from nwkit import __version__
from nwkit.disparity import (
    clade_disparities,
    curve_area,
    fit_brownian,
    make_design,
    simulate_curves,
    transform_traits,
    validate_work,
)
from nwkit.evolution import build_evolutionary_process, tree_depths
from nwkit.output_transaction import output_transaction, validate_output_targets
from nwkit.rooting_state import require_rooted, set_rooting_info
from nwkit.trait_input import numeric_trait_value, parse_trait_columns
from nwkit.util import (
    assign_branch_ids,
    get_node_class,
    read_tip_table,
    read_tree,
    validate_outputs_do_not_replace_inputs,
    validate_unique_named_leaves,
)

_OUTPUTS = (
    "outfile",
    "summary_out",
    "clades_out",
    "simulations_out",
    "model_out",
    "figure_out",
)
_ID = "_nwkit_dtt_input_branch_id"


def _options(args):
    if args.figure_columns is not None:
        displayed = parse_trait_columns(args.figure_columns)
        calculated = parse_trait_columns(args.columns)
        if any(column not in calculated for column in displayed):
            raise ValueError("--figure-columns must be a subset of --columns.")
    if args.n_sim != 0 and not 2 <= args.n_sim <= 10000:
        raise ValueError("--n-sim must be 0 or between 2 and 10,000.")
    if not 1 <= args.threads <= 32:
        raise ValueError("--threads must be between 1 and 32.")
    if args.seed < 0:
        raise ValueError("--seed must be nonnegative.")
    if not np.isfinite(args.ci_level) or not 0 < args.ci_level < 1:
        raise ValueError("--ci-level must be finite and strictly between 0 and 1.")
    try:
        interval = tuple(float(part) for part in args.mdi_range.split(","))
    except ValueError as exc:
        raise ValueError(
            "--mdi-range must be two numbers satisfying 0 <= start < end <= 1."
        ) from exc
    if len(interval) != 2 or not 0 <= interval[0] < interval[1] <= 1:
        raise ValueError("--mdi-range must satisfy 0 <= start < end <= 1.")
    if args.n_sim == 0 and args.simulations_out:
        raise ValueError("--simulations-out requires BM simulations (--n-sim >= 2).")
    paths = {}
    for role in _OUTPUTS:
        path = getattr(args, role)
        if path == "" or (path == "-" and role != "outfile"):
            raise ValueError(
                f"--{role.replace('_', '-')} requires a nonempty file path."
            )
        if role == "outfile" and path is None:
            raise ValueError("--outfile must be a path or '-'.")
        if path not in (None, "-"):
            paths[role] = path
    validate_outputs_do_not_replace_inputs(
        [("--infile", args.infile), ("--trait", args.trait)],
        [("--" + role.replace("_", "-"), path) for role, path in paths.items()],
    )
    validate_output_targets(paths.values())
    if args.figure_out and Path(args.figure_out).suffix.lower() not in {
        ".png",
        ".pdf",
        ".svg",
    }:
        raise ValueError("--figure-out requires PNG, PDF or SVG.")
    return paths, interval


def _ultrametric_depths(tree):
    depths = tree_depths(tree, allow_zero=True)
    if not all(np.isfinite(value) for value in depths.values()):
        raise ValueError("DTT tree depths overflow; rescale branch lengths.")
    tips = np.asarray([depths[node] for node in tree.leaves()])
    height = float(tips.max())
    if height <= 0:
        raise ValueError("DTT requires a positive crown age.")
    if np.ptp(tips / height) > 1e-8:
        raise ValueError(
            "DTT requires an ultrametric tree; relative tip-depth spread must be <= 1e-8."
        )
    return depths, height


def _prune_to_crown(tree, names):
    identifiers = assign_branch_ids(tree)
    crown = tree.common_ancestor(names)
    retained = set(names)
    keep: dict[Any, bool] = {}
    for node in crown.traverse("postorder"):
        keep[node] = (
            node.name in retained
            if node.is_leaf
            else any(keep[c] for c in node.children)
        )
    # ETE's deepcopy recursively copies child links and fails on valid deep trees.
    # Copy the retained topology iteratively; node property dictionaries are distinct.
    copies = {}
    for node in crown.traverse("preorder"):
        if not keep[node]:
            continue
        copy = type(tree)()
        copy.props.update(node.props)
        copy.add_prop(_ID, identifiers[node])
        copies[node] = copy
        if node is not crown:
            copies[node.up].add_child(copy)
    root = copies[crown]
    set_rooting_info(root, True, source="operation")
    return root


def _read_data(args):
    columns = parse_trait_columns(args.columns)
    if len(columns) > 64:
        raise ValueError("DTT supports at most 64 selected traits.")
    tree = read_tree(
        args.infile, args.format, args.quoted_node_names, rooted=args.input_rooted
    )
    require_rooted(tree, "DTT requires a rooted tree.")
    validate_unique_named_leaves(tree, "--infile")
    original_depths, _ = _ultrametric_depths(tree)
    if len(original_depths) > 4000:
        raise ValueError("DTT supports at most 4,000 input nodes.")
    names = sorted(node.name for node in tree.leaves())
    if len(names) > 2000:
        raise ValueError("DTT supports at most 2,000 input tips.")
    table, _, _ = read_tip_table(
        args.trait,
        tree_leaf_names=names,
        required_columns=columns,
        missing_values=args.missing_values,
        unmatched=args.unmatched,
    )
    table = table.set_index("leaf_name").reindex(names)
    values = np.asarray(
        [
            [
                numeric_trait_value(
                    table.loc[name, column], column, args.missing_values
                )
                for column in columns
            ]
            for name in names
        ]
    )
    complete = np.isfinite(values).all(axis=1)
    excluded = [name for name, valid in zip(names, complete, strict=True) if not valid]
    if excluded and args.missing == "error":
        raise ValueError(
            "DTT requires complete tips; use --missing drop to omit: "
            + ", ".join(excluded)
        )
    used = [name for name, valid in zip(names, complete, strict=True) if valid]
    if len(used) < 3:
        raise ValueError("DTT requires at least three complete tips.")
    retained = _prune_to_crown(tree, used)
    original_depths, height = _ultrametric_depths(retained)
    for node in retained.traverse():
        if not node.is_root:
            length = float(node.dist) / height
            if node.dist > 0 and length == 0:
                raise ValueError(
                    "DTT normalized branch lengths underflow; reduce tree scale disparity."
                )
            node.dist = length
    # Divide already accumulated depths once. Re-summing normalized edges can
    # round a present-day internal split above one, making event times unsorted.
    depths = {node: depth / height for node, depth in original_depths.items()}
    # Tips can differ within the declared ultrametric roundoff tolerance.
    for node in retained.leaves():
        depths[node] = 1.0
    return retained, columns, used, excluded, values[complete], depths, height


def _result_tables(tree, design, values, simulations, height, interval, args):
    disparities, counts = clade_disparities(design, values)
    relative = disparities / disparities[0]
    observed = np.asarray(design.weights @ relative)
    table = pd.DataFrame(
        {
            "time_index": np.arange(1, len(design.times) + 1),
            "phase": [
                "root_before_split",
                *(["after_split"] * (len(design.times) - 2)),
                "present",
            ],
            "relative_time": design.times,
            "time_from_crown": design.times * height,
            "relative_disparity": observed,
            "num_clades": design.clade_counts,
        }
    )
    summary = {
        "num_tips": int(counts[0]),
        "num_traits": values.shape[1],
        "num_simulations": len(simulations),
        "scale": args.scale,
        "metric": "avg.sq",
        "crown_age": height,
        "mdi_start": interval[0],
        "mdi_end": interval[1],
        "envelope_level": args.ci_level if len(simulations) else None,
        "observed_area": float(curve_area(design.times, observed, interval)),
        "bm_median_area": None,
        "mdi": None,
        "bm_mdi_q_lower": None,
        "bm_mdi_q_upper": None,
        "status": "ok"
        if np.any((design.times > 0) & (design.times < 1))
        else "uninformative_topology",
    }
    null_mdi = np.empty(0)
    if len(simulations):
        tail = (1 - args.ci_level) / 2
        median = np.median(simulations, axis=0)
        table["bm_mean"] = simulations.mean(axis=0)
        table["bm_median"] = median
        table["bm_lower"], table["bm_upper"] = np.quantile(
            simulations, [tail, 1 - tail], axis=0
        )
        null_mdi = curve_area(design.times, simulations - median, interval)
        summary.update(
            bm_median_area=float(curve_area(design.times, median, interval)),
            mdi=float(curve_area(design.times, observed - median, interval)),
            bm_mdi_q_lower=float(np.quantile(null_mdi, tail)),
            bm_mdi_q_upper=float(np.quantile(null_mdi, 1 - tail)),
        )
    else:
        for column in ("bm_mean", "bm_median", "bm_lower", "bm_upper"):
            table[column] = ""
    table["num_simulations"] = len(simulations)
    clades = pd.DataFrame(
        [
            {
                "branch_id": node.props[_ID],
                "parent_branch_id": -1 if node.is_root else node.up.props[_ID],
                "name": node.name or "",
                "node_class": get_node_class(node),
                "num_tips": int(counts[i]),
                "relative_disparity": relative[i],
            }
            for i, node in enumerate(tree.traverse("preorder"))
        ]
    )
    return table, pd.DataFrame([summary]), clades, null_mdi


def _simulation_table(times, simulations):
    return pd.DataFrame(
        {
            "simulation": np.repeat(np.arange(1, len(simulations) + 1), len(times)),
            "time_index": np.tile(np.arange(1, len(times) + 1), len(simulations)),
            "relative_time": np.tile(times, len(simulations)),
            "relative_disparity": simulations.ravel(),
        }
    )


def dtt_main(args):
    paths, interval = _options(args)
    tree, columns, names, excluded, original, depths, height = _read_data(args)
    design = make_design(tree, names, depths)
    validate_work(design, len(columns), args.n_sim, bool(args.simulations_out))
    values, center, scales = transform_traits(original, args.scale)
    simulations = np.empty((0, len(design.times)))
    rate = bm_center = None
    if args.n_sim:
        covariance = build_evolutionary_process(tree, allow_zero=True).tip_covariance(
            names
        )
        bm_center, rate, factor = fit_brownian(covariance, values)
        simulations = simulate_curves(
            design, factor, args.n_sim, args.seed, args.threads
        )
    curve, summary, clades, null_mdi = _result_tables(
        tree, design, values, simulations, height, interval, args
    )
    tables = {"outfile": curve, "summary_out": summary, "clades_out": clades}
    if args.simulations_out:
        tables["simulations_out"] = _simulation_table(design.times, simulations)
    metadata = {
        "schema_version": 1,
        "nwkit_version": __version__,
        "metric": "avg.sq",
        "trait_names": columns,
        "scale": args.scale,
        "used_taxa": names,
        "excluded_taxa": excluded,
        "crown_root_branch_id": tree.props[_ID],
        "crown_age": height,
        "data_center": center.tolist(),
        "data_scale": scales.tolist(),
        "data_transform": "(X - data_center) / data_scale",
        "bm_center_transformed": None if bm_center is None else bm_center.tolist(),
        "bm_rate_covariance_transformed": None if rate is None else rate.tolist(),
        "bm_rate_time_units": "relative crown time",
        "bm_rate_divisor": len(names) - 1,
        "num_simulations": args.n_sim,
        "seed": args.seed,
        "mdi_range": list(interval),
        "mdi": summary.mdi.iloc[0],
        "mdi_reference": "pointwise BM median; linear interpolation and trapezoidal integration",
        "clade_average": "unweighted mean of active clades with at least two retained tips; singleton lineages excluded",
        "uncertainty": "Pointwise BM null simulation envelope conditional on fitted covariance, scaling and retained tree; not a confidence band. Parameter/tree uncertainty excluded.",
    }
    serialized = {
        role: table.to_csv(sep="\t", index=False, float_format="%.12g")
        for role, table in tables.items()
    }
    model_json = (
        json.dumps(metadata, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    )
    with output_transaction(paths.values()) as staged:
        for role, path in paths.items():
            if role == "figure_out":
                from nwkit.dtt_figure import draw_dtt

                draw_dtt(
                    curve,
                    summary.iloc[0],
                    null_mdi,
                    staged[path],
                    Path(path).suffix[1:],
                    columns,
                    tree,
                    names,
                    original,
                    figure_layout=args.figure_layout,
                    figure_columns=args.figure_columns,
                    figure_scale=args.figure_scale,
                )
            else:
                Path(staged[path]).write_text(
                    model_json if role == "model_out" else serialized[role],
                    encoding="utf-8",
                )
    if args.outfile == "-":
        sys.stdout.write(serialized["outfile"])
    sys.stderr.write(
        f"DTT: {len(names)} tips, {len(columns)} traits, {args.n_sim} BM simulations; singleton lineages excluded.\n"
    )
