"""Strict inputs for native multivariate shift fitting and discovery."""

import numpy as np
import pandas as pd

from nwkit.asr_regimes import read_regime_map
from nwkit.rooting_state import require_rooted
from nwkit.shift_native_fit import NativeFitOptions
from nwkit.shift_native_model import ShiftData, ShiftLayout, ShiftTree
from nwkit.util import read_tip_table, read_tree, validate_unique_named_leaves


def _columns(value, option):
    columns = tuple(item.strip() for item in str(value).split(","))
    if (
        not columns
        or any(not item or item == "leaf_name" for item in columns)
        or len(set(columns)) != len(columns)
    ):
        raise ValueError(
            f"{option} requires distinct nonempty columns different from leaf_name."
        )
    return columns


def numeric_parameters(value, count, option):
    if value is None:
        return None
    try:
        parameters = np.asarray([float(item) for item in str(value).split(",")])
    except ValueError as exc:
        raise ValueError(f"{option} requires numeric values.") from exc
    if (
        len(parameters) not in {1, count}
        or np.isnan(parameters).any()
        or np.any(parameters < 0)
    ):
        raise ValueError(f"{option} requires one nonnegative value or one per trait.")
    if option != "--alpha" and not np.isfinite(parameters).all():
        raise ValueError(f"{option} must be finite.")
    return np.broadcast_to(parameters, (count,)).copy()


def read_native_data(args):
    tree = read_tree(
        args.infile, args.format, args.quoted_node_names, rooted=args.input_rooted
    )
    require_rooted(tree, "Native shift inference requires a rooted tree.")
    validate_unique_named_leaves(tree, "--infile")
    prepared = ShiftTree.build(tree)
    columns = _columns(args.state_column, "--state-column")
    error_columns = (
        _columns(args.standard_error_column, "--standard-error-column")
        if args.standard_error_column is not None
        else ()
    )
    if error_columns and (
        len(error_columns) != len(columns) or set(columns) & set(error_columns)
    ):
        raise ValueError(
            "Native standard-error columns must match the traits and differ from their columns."
        )
    table, _, _ = read_tip_table(
        args.trait,
        tree_leaf_names=prepared.leaf_names,
        required_columns=(*columns, *error_columns),
        unmatched="error",
    )
    table = table.set_index("leaf_name").loc[list(prepared.leaf_names)]
    values = np.column_stack(
        [
            pd.to_numeric(table[column], errors="raise").to_numpy(
                dtype=float, na_value=np.nan
            )
            for column in columns
        ]
    )
    errors = np.zeros_like(values)
    if error_columns:
        errors = np.column_stack(
            [
                pd.to_numeric(table[column], errors="raise").to_numpy(
                    dtype=float, na_value=np.nan
                )
                for column in error_columns
            ]
        )
        if (
            np.isinf(errors).any()
            or np.any(errors < 0)
            or np.any(np.isnan(errors) & np.isfinite(values))
        ):
            raise ValueError(
                "Each observed trait value requires a finite, nonnegative standard error."
            )
        errors = np.where(np.isnan(errors), 0.0, errors)
    with np.errstate(over="ignore", under="ignore"):
        variances = errors**2
    if not np.isfinite(variances).all() or np.any((errors > 0) & (variances == 0)):
        raise ValueError("Squared standard errors are not representable.")
    return ShiftData.build(prepared, values, columns, variances)


def read_native_layout(path, tree):
    assignment = read_regime_map(path, tree.compiled.tree)
    groups: dict[str, list[int]] = {}
    shifts = []
    for index, node in enumerate(tree.compiled.nodes):
        label = assignment.by_node[node]
        if index == 0 or label != assignment.by_node[node.up]:
            groups.setdefault(label, []).append(tree.branch_ids[index])
            if index:
                shifts.append(tree.branch_ids[index])
    layout = ShiftLayout.build(tree, shifts, groups.values())
    names = {min(group): label for label, group in groups.items()}
    return layout, [names[min(group)] for group in layout.groups]


def native_fit_arguments(args, data):
    if args.fit_out or args.criterion is not None:
        raise ValueError(
            "Native shift inference exports JSON; --fit-out and IC criteria belong to --selection ic."
        )
    alpha = numeric_parameters(args.alpha, len(data.trait_names), "--alpha")
    if alpha is not None:
        original = alpha.copy()
        with np.errstate(over="ignore", under="ignore"):
            alpha *= data.tree.height
        if np.any(np.isfinite(original) & ~np.isfinite(alpha)) or np.any(
            (original > 0) & (alpha == 0)
        ):
            raise ValueError(
                "Alpha times tree height is not representable; rescale time units."
            )
    options = NativeFitOptions(
        root_model=args.root_model,
        estimate_measurement_error=args.estimate_measurement_error,
        optimizer_starts=args.optimizer_starts,
        maxiter=args.max_iterations,
    )
    options.validate()
    return {
        "options": options,
        "alpha_height": alpha,
        "process_variance": numeric_parameters(
            args.process_tip_variance, len(data.trait_names), "--process-tip-variance"
        ),
        "measurement_variance": numeric_parameters(
            args.measurement_variance, len(data.trait_names), "--measurement-variance"
        ),
    }
