"""Input validation and batch reporting for phylogenetic signal."""

import hashlib
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from nwkit.evolution import build_evolutionary_process
from nwkit.file_paths import validate_outputs_do_not_replace_inputs
from nwkit.output_transaction import output_transaction, validate_output_targets
from nwkit.rooting_state import require_rooted
from nwkit.signal_stats import bh_adjust, k_fit, lambda_fit, permutation_test
from nwkit.trait_input import numeric_trait_value as _number
from nwkit.trait_input import parse_trait_columns as _columns
from nwkit.util import (
    read_tip_table,
    read_tree,
    validate_unique_named_leaves,
)

COLUMNS = [
    "trait",
    "method",
    "num_taxa",
    "num_missing_taxa",
    "estimate",
    "status",
    "message",
    "measurement_error",
    "sigma2",
    "root_mean",
    "log_likelihood",
    "null_log_likelihood",
    "likelihood_ratio",
    "test_method",
    "p_value",
    "p_adjusted",
    "p_adjust",
    "num_tests",
    "num_simulations",
    "seed",
    "ci_level",
    "ci_lower",
    "ci_upper",
    "lambda_lower",
    "lambda_upper",
]


def _inputs(args):
    columns = _columns(args.columns)
    se_columns = (
        _columns(args.standard_error_column) if args.standard_error_column else []
    )
    if se_columns and (
        len(se_columns) != len(columns) or set(se_columns) & set(columns)
    ):
        raise ValueError("Provide one distinct SE column per trait, in matching order.")
    if args.n_sim < 1 or args.seed < 0:
        raise ValueError("--n-sim must be positive and --seed non-negative.")
    if not np.isfinite(args.ci_level) or not 0 < args.ci_level < 1:
        raise ValueError("--ci-level must be finite and between zero and one.")
    tree = read_tree(
        args.infile, args.format, args.quoted_node_names, rooted=args.input_rooted
    )
    require_rooted(tree, "Phylogenetic signal requires a rooted tree.")
    validate_unique_named_leaves(tree, "--infile")
    names = sorted(node.name for node in tree.leaves())
    process = build_evolutionary_process(
        tree, model="brownian", root_mode="fixed", allow_zero=True
    )
    covariance = process.tip_covariance(names)
    if not np.all(np.isfinite(covariance)):
        raise ValueError("Tree covariance must be finite.")
    table, _, _ = read_tip_table(
        args.trait,
        tree_leaf_names=names,
        required_columns=columns + se_columns,
        missing_values=args.missing_values,
        unmatched=args.unmatched,
    )
    table = table.set_index("leaf_name").reindex(names)
    data = []
    for index, column in enumerate(columns):
        values = np.array(
            [_number(x, column, args.missing_values) for x in table[column]]
        )
        errors = np.zeros(len(names))
        if se_columns:
            errors = np.array(
                [
                    _number(x, se_columns[index], args.missing_values)
                    for x in table[se_columns[index]]
                ]
            )
            observed = np.isfinite(values)
            if np.any(~np.isfinite(errors[observed])) or np.any(errors[observed] < 0):
                raise ValueError(
                    f"Observed trait '{column}' requires finite, non-negative SEs."
                )
        data.append((column, values, errors))
    return covariance, data


def _trait_rows(covariance, column, values, errors, args):
    observed = np.isfinite(values)
    n = int(observed.sum())
    methods = ["K", "lambda"] if args.method == "both" else [args.method]
    base = {
        "trait": column,
        "num_taxa": n,
        "num_missing_taxa": len(values) - n,
        "measurement_error": "yes" if args.standard_error_column else "no",
    }
    status = "insufficient_taxa" if n < 3 else None
    values, errors = values[observed], errors[observed]
    if not status and np.all(values == values[0]):
        status = "constant_trait"
    if status:
        return [dict(base, method=method, status=status) for method in methods]
    covariance = covariance[np.ix_(observed, observed)].copy()
    # Pruning to the observed tips drops the shared stem above their MRCA.
    covariance -= np.min(covariance)
    try:
        np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError:
        return [
            dict(base, method=method, status="singular_covariance")
            for method in methods
        ]
    offset = float(values[0])
    magnitude = max(float(np.max(np.abs(values))), float(np.max(errors)))
    scale = math.ldexp(1.0, math.frexp(magnitude)[1] - 1)
    values = values / scale - offset / scale
    time_scale = float(np.max(np.diag(covariance)))
    errors = errors / scale
    covariance /= time_scale
    digest = int.from_bytes(hashlib.sha256(column.encode()).digest()[:8], "little")
    rng = np.random.default_rng(np.random.SeedSequence([args.seed, digest]))
    rows = []
    for method in methods:
        try:
            result = _method_result(method, covariance, values, errors, args, rng)
            result = _restore_units(result, scale, time_scale, offset, n)
        except (
            ValueError,
            np.linalg.LinAlgError,
            OverflowError,
            FloatingPointError,
        ) as exc:
            result = {"status": "fit_failed", "message": str(exc)}
        rows.append(dict(base, method=method, **result))
    return rows


def _restore_units(result, scale, time_scale, offset, n):
    if result.get("sigma2") is not None and result["sigma2"] != 0:
        rate, rate_exponent = math.frexp(result["sigma2"])
        unit, unit_exponent = math.frexp(scale)
        time, time_exponent = math.frexp(time_scale)
        result["sigma2"] = math.ldexp(
            rate * unit**2 / time, rate_exponent + 2 * unit_exponent - time_exponent
        )
        if result["sigma2"] == 0:
            raise ValueError(
                "Estimated diffusion rate underflows in original units; rescale units."
            )
    if result.get("root_mean") is not None:
        result["root_mean"] = (offset / scale + result["root_mean"]) * scale
    for key in ("log_likelihood", "null_log_likelihood"):
        if key in result:
            result[key] -= n * np.log(scale)
    if any(
        isinstance(value, (float, np.floating)) and not np.isfinite(value)
        for value in result.values()
    ):
        raise ValueError(
            "Non-finite signal result; rescale units or inspect numerical conditioning."
        )
    return result


def _method_result(method, covariance, values, errors, args, rng):
    if method == "lambda":
        result = lambda_fit(covariance, values, errors, args.ci_level)
        result.update(lambda_lower=0.0, lambda_upper=1.0)
        if "ci_lower" in result:
            result["ci_level"] = args.ci_level
        if args.test == "no":
            result.pop("p_value", None)
        elif "p_value" in result:
            result["test_method"] = "likelihood_ratio_chi2_1"
        return result
    estimate, rate = k_fit(covariance, values, errors)
    result = {"estimate": estimate, "sigma2": rate, "status": "ok"}
    if args.test == "yes":
        result.update(
            p_value=permutation_test(
                covariance, values, errors, estimate, args.n_sim, rng
            ),
            test_method="upper_tail_permutation",
            num_simulations=args.n_sim,
            seed=args.seed,
        )
    return result


def signal_main(args):
    outputs = [] if args.outfile == "-" else [args.outfile]
    validate_outputs_do_not_replace_inputs(
        [("--infile", args.infile), ("--trait", args.trait)],
        [("--outfile", path) for path in outputs],
    )
    validate_output_targets(outputs)
    covariance, data = _inputs(args)
    rows = []
    for column, values, errors in data:
        rows.extend(_trait_rows(covariance, column, values, errors, args))
    for method in ("K", "lambda"):
        tested = [row for row in rows if row["method"] == method and "p_value" in row]
        p_values = [row["p_value"] for row in tested]
        adjusted = bh_adjust(p_values) if args.p_adjust == "bh" else p_values
        for row, value in zip(tested, adjusted, strict=True):
            row.update(p_adjusted=value, p_adjust=args.p_adjust, num_tests=len(tested))
    text = (
        pd.DataFrame(rows)
        .reindex(columns=COLUMNS)
        .to_csv(sep="\t", index=False, na_rep="", float_format="%.12g")
    )
    with output_transaction(outputs) as staged:
        if outputs:
            Path(staged[args.outfile]).write_text(text, encoding="utf-8")
        else:
            sys.stdout.write(text)
