"""Tip-table PCA, stable input-tree node IDs and shared-ASR projection."""

import json
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import norm

from nwkit import __version__
from nwkit.continuous_asr import compute_bm_marginals
from nwkit.evolution import build_evolutionary_process
from nwkit.file_paths import validate_outputs_do_not_replace_inputs
from nwkit.output_transaction import output_transaction, validate_output_targets
from nwkit.phylogenetic_pca import fit_pca
from nwkit.rooting_state import require_rooted, set_rooting_info
from nwkit.trait_input import numeric_trait_value, parse_trait_columns
from nwkit.util import (
    assign_branch_ids,
    get_node_class,
    read_tip_table,
    read_tree,
    validate_unique_named_leaves,
)

_OUTPUTS = (
    "outfile",
    "loadings_out",
    "eigenvalues_out",
    "model_out",
    "ancestral_out",
    "figure_out",
)
_ID = "_nwkit_pca_input_branch_id"


def _validate_outputs(args):
    paths = {}
    for role in _OUTPUTS:
        path = getattr(args, role)
        if path == "-" and role != "outfile":
            raise ValueError(f"--{role.replace('_', '-')} requires a file path.")
        if path not in (None, "-"):
            paths[role] = path
    validate_outputs_do_not_replace_inputs(
        [("--infile", args.infile), ("--trait", args.trait)],
        [("--" + role.replace("_", "-"), path) for role, path in paths.items()],
    )
    validate_output_targets(paths.values())
    if not np.isfinite(args.ci_level) or not 0 < args.ci_level < 1:
        raise ValueError("--ci-level must be finite and between zero and one.")
    if args.figure_out and Path(args.figure_out).suffix.lower() not in {
        ".pdf",
        ".svg",
        ".png",
    }:
        raise ValueError("--figure-out requires PDF, SVG or PNG.")
    return paths


def _retained_tree(tree, names):
    for node, identifier in assign_branch_ids(tree).items():
        node.add_prop(_ID, identifier)
    retained = set(names)
    root = tree.common_ancestor(names).copy(method="deepcopy")
    keep: dict[Any, bool] = {}
    for node in root.traverse("postorder"):
        keep[node] = (
            node.name in retained
            if node.is_leaf
            else any(keep[child] for child in node.children)
        )
    for node in list(root.traverse("postorder")):
        if not keep[node] and not node.is_root:
            node.detach()
    set_rooting_info(root, True, source="operation")
    return root


def _read_data(args):
    columns = parse_trait_columns(args.columns)
    if len(columns) < 2:
        raise ValueError("PCA requires at least two trait columns.")
    tree = read_tree(
        args.infile, args.format, args.quoted_node_names, rooted=args.input_rooted
    )
    require_rooted(tree, "Phylogenetic PCA requires a rooted tree.")
    validate_unique_named_leaves(tree, "--infile")
    # Validate every original branch, including branches of subsequently dropped tips.
    build_evolutionary_process(tree, allow_zero=True)
    names = sorted(node.name for node in tree.leaves())
    table, _, _ = read_tip_table(
        args.trait,
        tree_leaf_names=names,
        required_columns=columns,
        missing_values=args.missing_values,
        unmatched=args.unmatched,
    )
    table = table.set_index("leaf_name").reindex(names)
    values = np.array(
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
    complete = np.all(np.isfinite(values), axis=1)
    excluded = [name for name, valid in zip(names, complete, strict=True) if not valid]
    if excluded and args.missing == "error":
        raise ValueError(
            "PCA requires complete trait vectors; use --missing drop to omit: "
            + ", ".join(excluded)
        )
    if int(complete.sum()) < 3:
        raise ValueError("PCA requires at least three complete tips.")
    used = [name for name, valid in zip(names, complete, strict=True) if valid]
    return _retained_tree(tree, used), columns, used, excluded, values[complete]


def _tables(fit, columns, names):
    components = [f"PC{i + 1}" for i in range(len(fit.eigenvalues))]
    scores = pd.DataFrame(fit.scores, columns=components)
    scores.insert(0, "leaf_name", names)
    # Normalize before summing to avoid overflow of individually finite eigenvalues.
    fractions = fit.eigenvalues / fit.eigenvalues[0]
    fractions /= fractions.sum()
    eigenvalues = pd.DataFrame(
        {
            "component": components,
            "eigenvalue": fit.eigenvalues,
            "explained_variance_ratio": fractions,
            "cumulative_variance_ratio": np.cumsum(fractions),
        }
    )
    loadings = pd.DataFrame(
        [
            {
                "trait": name,
                "component": component,
                "rotation": fit.rotation[i, j],
                "loading": fit.loadings[i, j],
                "center": fit.center[i],
                "scale": fit.scale[i],
            }
            for i, name in enumerate(columns)
            for j, component in enumerate(components)
        ]
    )
    return scores, eigenvalues, loadings


def ancestral_scores(tree, names, fit, level):
    process = build_evolutionary_process(
        tree,
        model="lambda",
        parameter=fit.lambda_value,
        root_mode="flat",
        allow_zero=True,
    )
    rows = []
    z = norm.ppf((1 + level) / 2)
    for j, eigenvalue in enumerate(fit.eigenvalues):
        posterior, _ = compute_bm_marginals(
            tree,
            dict(zip(names, fit.scores[:, j], strict=True)),
            sigma2=float(eigenvalue),
            _process=process,
        )
        for node, marginal in posterior.items():
            interval = z * np.sqrt(marginal.variance)
            rows.append(
                {
                    "branch_id": node.props[_ID],
                    "parent_branch_id": -1 if node.is_root else node.up.props[_ID],
                    "name": node.name or "",
                    "node_class": get_node_class(node),
                    "component": f"PC{j + 1}",
                    "mean": marginal.mean,
                    "variance": marginal.variance,
                    "ci_level": level,
                    "ci_lower": marginal.mean - interval,
                    "ci_upper": marginal.mean + interval,
                }
            )
    return pd.DataFrame(rows).sort_values(["branch_id", "component"])


def _metadata(tree, columns, names, excluded, fit, args):
    return {
        "schema_version": 1,
        "nwkit_version": __version__,
        "model": args.model,
        "mode": args.mode,
        "lambda": fit.lambda_value,
        "lambda_estimated": fit.lambda_estimated,
        "lambda_bounds": [0, 1],
        "log_likelihood": fit.log_likelihood,
        "likelihood_method": "ML",
        "pca_covariance_divisor": len(names) - 1,
        "status": fit.status,
        "rank": len(fit.eigenvalues),
        "repeated_eigenvalues": fit.repeated_eigenvalues,
        "trait_names": columns,
        "used_taxa": names,
        "excluded_taxa": excluded,
        "root_branch_id": tree.props[_ID],
        "center": fit.center.tolist(),
        "scale": fit.scale.tolist(),
        "rotation": fit.rotation.tolist(),
        "loadings": fit.loadings.tolist(),
        "eigenvalues": fit.eigenvalues.tolist(),
        "score_transform": "((X - center) / scale) @ rotation",
        "ancestral_uncertainty": "Conditional on the fitted tree, lambda, PCA axes and evolutionary variances; includes root-mean uncertainty. Connectors are not sampled evolutionary paths.",
    }


def pca_main(args):
    paths = _validate_outputs(args)
    tree, columns, names, excluded, values = _read_data(args)
    covariance = build_evolutionary_process(tree, allow_zero=True).tip_covariance(names)
    fit = fit_pca(
        covariance,
        values,
        model=args.model,
        mode=args.mode,
        lambda_value=args.lambda_value,
    )
    scores, eigenvalues, loadings = _tables(fit, columns, names)
    ancestors = (
        ancestral_scores(tree, names, fit, args.ci_level)
        if args.ancestral_out or args.figure_out
        else None
    )
    tables = {
        "outfile": scores,
        "loadings_out": loadings,
        "eigenvalues_out": eigenvalues,
        "ancestral_out": ancestors,
    }
    serialized = {
        role: table.to_csv(sep="\t", index=False, float_format="%.12g")
        for role, table in tables.items()
        if table is not None
    }
    metadata = (
        json.dumps(
            _metadata(tree, columns, names, excluded, fit, args),
            indent=2,
            ensure_ascii=False,
            allow_nan=False,
        )
        + "\n"
    )
    with output_transaction(paths.values()) as staged:
        for role, path in paths.items():
            if role == "figure_out":
                from nwkit.pca_figure import draw_pca

                draw_pca(ancestors, loadings, eigenvalues, staged[path], args)
            else:
                Path(staged[path]).write_text(
                    metadata if role == "model_out" else serialized[role],
                    encoding="utf-8",
                )
        if args.outfile == "-":
            sys.stdout.write(serialized["outfile"])
