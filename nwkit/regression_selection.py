"""Nested group CV and transactional artifacts for exploratory gene selection."""

import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd

from nwkit.file_paths import validate_outputs_do_not_replace_inputs
from nwkit.ordinary_regression import build_phylogenetic_covariance
from nwkit.output_transaction import output_transaction
from nwkit.penalized_phylogenetic import (
    SelectionFit,
    fit_elastic_net,
    prediction_loss,
    validate_inputs,
)
from nwkit.regression_selection_cli import selection_output_paths
from nwkit.rooting_state import require_rooted
from nwkit.util import read_tree


class SelectionFailure(RuntimeError):
    """A tuning grid had no candidate successful in every required fold."""


def _grid(strengths, ratios):
    if not strengths or not ratios:
        raise ValueError("The regularization grid must not be empty.")
    if any(not np.isfinite(v) or v <= 0 for v in strengths):
        raise ValueError("All strengths must be positive finite numbers.")
    if any(not np.isfinite(v) or not 0 < v <= 1 for v in ratios):
        raise ValueError("All l1 ratios must be in (0, 1].")
    if len(set(strengths)) != len(strengths) or len(set(ratios)) != len(ratios):
        raise ValueError("Duplicate grid values are not allowed.")
    return [
        (s, r)
        for r in sorted(ratios, reverse=True)
        for s in sorted(strengths, reverse=True)
    ]


def _fit_rows(
    x, y, covariance, rows, family, strength, ratio, unpenalized, initial=None
):
    return fit_elastic_net(
        x[rows],
        y[rows],
        covariance[np.ix_(rows, rows)],
        family=family,
        strength=strength,
        l1_ratio=ratio,
        unpenalized=unpenalized,
        initial=initial,
    )


def _predict_rows(fit, x, covariance, train, test, prediction):
    cross = covariance[np.ix_(test, train)] if prediction == "conditional" else None
    return fit.predict(x[test], cross_covariance=cross)


def _tune(x, y, covariance, groups, rows, grid, family, unpenalized, prediction, stage):
    records = []
    labels = sorted(set(groups[rows]))
    for fold in labels:
        train = rows[groups[rows] != fold]
        test = rows[groups[rows] == fold]
        previous: dict[float, SelectionFit] = {}
        for strength, ratio in grid:
            record = dict(
                stage=stage,
                fold=fold,
                strength=strength,
                l1_ratio=ratio,
                n_train=len(train),
                n_test=len(test),
                status="ok",
                error="",
                boundary_warning="",
                loss=np.nan,
            )
            try:
                fitted = _fit_rows(
                    x,
                    y,
                    covariance,
                    train,
                    family,
                    strength,
                    ratio,
                    unpenalized,
                    previous.get(ratio),
                )
                previous[ratio] = fitted
                predicted = _predict_rows(
                    fitted, x, covariance, train, test, prediction
                )
                record["loss"] = float(
                    np.mean(prediction_loss(y[test], predicted, family))
                )
                record["boundary_warning"] = fitted.boundary_warning
            except (
                ValueError,
                RuntimeError,
                FloatingPointError,
                np.linalg.LinAlgError,
            ) as exc:
                record["status"], record["error"] = "failed", str(exc)
            records.append(record)
    table = pd.DataFrame(records)
    candidates = []
    for strength, ratio in grid:
        subset = table[(table.strength == strength) & (table.l1_ratio == ratio)]
        if (subset.status == "ok").all() and np.isfinite(subset.loss).all():
            # Each held-out tip has equal weight, even when clades differ in size.
            score = float(np.average(subset.loss, weights=subset.n_test))
            candidates.append((score, -strength, -ratio))
    if not candidates:
        errors = "; ".join(table.loc[table.status != "ok", "error"].unique()[:3])
        raise SelectionFailure(
            f"No tuning candidate succeeded across all folds in {stage}: {errors}"
        )
    _, negative_strength, negative_ratio = min(candidates)
    return -negative_strength, -negative_ratio, table


def _baseline_prediction(
    x, y, covariance, train, test, family, unpenalized, prediction
):
    # Retain the same prespecified covariates and tree, but no selectable genes.
    reduced = x[:, list(unpenalized)]
    try:
        fit = _fit_rows(
            reduced,
            y,
            covariance,
            train,
            family,
            0.1,
            0.5,
            tuple(range(reduced.shape[1])),
        )
        values = _predict_rows(fit, reduced, covariance, train, test, prediction)
        return values, prediction_loss(y[test], values, family), "ok", ""
    except (ValueError, RuntimeError, FloatingPointError, np.linalg.LinAlgError) as exc:
        return (
            np.full(len(test), np.nan),
            np.full(len(test), np.nan),
            "failed",
            str(exc),
        )


def nested_selection(
    x,
    y,
    covariance,
    groups,
    *,
    family="gaussian",
    strengths=(1.0, 0.1, 0.01),
    l1_ratios=(1.0, 0.5),
    unpenalized=(),
    prediction="conditional",
):
    """Leave-one-group-out outer evaluation with inner group tuning.

    A user supplies phylogenetic groups. Grouping is not claimed to make species
    independent. Fold-specific scales, constants and all fitted parameters use
    training data only. Outer predictions never tune the final model.
    """
    x, y, covariance = validate_inputs(x, y, covariance, family)
    # Fix the covariance scale from the supplied tree, independent of traits.
    covariance = covariance / float(np.mean(np.diag(covariance)))
    groups = np.asarray(groups, dtype=str)
    if (
        groups.ndim != 1
        or len(groups) != len(y)
        or any(not group.strip() for group in groups)
    ):
        raise ValueError("Every tip needs one non-empty group.")
    labels = sorted(set(groups))
    if len(labels) < 3:
        raise ValueError("Nested group CV requires at least three groups.")
    if prediction not in ("conditional", "fixed"):
        raise ValueError("prediction must be conditional or fixed.")
    grid = _grid(strengths, l1_ratios)
    rows = np.arange(len(y))
    # Fail before expensive work if an inner split is intrinsically unusable.
    for outer in labels:
        for inner in labels:
            if outer == inner:
                continue
            train = rows[(groups != outer) & (groups != inner)]
            validate_inputs(
                x[train], y[train], covariance[np.ix_(train, train)], family
            )
    predictions, fits, tuning = [], [], []
    for fold in labels:
        train, test = rows[groups != fold], rows[groups == fold]
        strength, ratio, table = _tune(
            x,
            y,
            covariance,
            groups,
            train,
            grid,
            family,
            unpenalized,
            prediction,
            f"outer:{fold}",
        )
        tuning.append(table)
        fitted = _fit_rows(
            x, y, covariance, train, family, strength, ratio, unpenalized
        )
        fits.append((fold, fitted))
        predicted = _predict_rows(fitted, x, covariance, train, test, prediction)
        losses = prediction_loss(y[test], predicted, family)
        baseline, baseline_loss, baseline_status, baseline_error = _baseline_prediction(
            x, y, covariance, train, test, family, unpenalized, prediction
        )
        for row, value, loss, base, base_loss in zip(
            test, predicted, losses, baseline, baseline_loss, strict=True
        ):
            predictions.append(
                dict(
                    row=int(row),
                    fold=fold,
                    observed=float(y[row]),
                    predicted=float(value),
                    loss=float(loss),
                    baseline_predicted=float(base),
                    baseline_loss=float(base_loss),
                    baseline_status=baseline_status,
                    baseline_error=baseline_error,
                    strength=strength,
                    l1_ratio=ratio,
                    boundary_warning=fitted.boundary_warning,
                )
            )
    strength, ratio, table = _tune(
        x,
        y,
        covariance,
        groups,
        rows,
        grid,
        family,
        unpenalized,
        prediction,
        "final_tuning",
    )
    tuning.append(table)
    # Full-data paths are exploratory; they are never used for outer evaluation.
    path = []
    previous: dict[float, SelectionFit] = {}
    for point_strength, point_ratio in grid:
        fitted = _fit_rows(
            x,
            y,
            covariance,
            rows,
            family,
            point_strength,
            point_ratio,
            unpenalized,
            previous.get(point_ratio),
        )
        previous[point_ratio] = fitted
        path.append(fitted)
    final = next(f for f in path if f.strength == strength and f.l1_ratio == ratio)
    return (
        final,
        path,
        fits,
        pd.concat(tuning, ignore_index=True),
        pd.DataFrame(predictions).sort_values("row"),
    )


def _names(text):
    names = [name.strip() for name in text.split(",") if name.strip()]
    if len(names) != len(set(names)):
        raise ValueError("Column lists must not contain duplicates.")
    return names


def _read_table(path):
    with open(path, encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, [])
        if (
            not header
            or any(not name.strip() for name in header)
            or len(header) != len(set(header))
        ):
            raise ValueError(
                f"Table column names must be unique and non-empty: {path}."
            )
        for line, row in enumerate(reader, 2):
            if len(row) != len(header):
                raise ValueError(
                    f"Table row {line} has the wrong number of fields: {path}."
                )
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if (
        "leaf_name" not in table
        or table.leaf_name.duplicated().any()
        or (table.leaf_name == "").any()
    ):
        raise ValueError(f"Table requires unique non-empty leaf_name keys: {path}.")
    return table.set_index("leaf_name")


def _coefficient_rows(fit, names, unpenalized):
    rows = []
    for i, name in enumerate(["(intercept)"] + names):
        penalty = i != 0 and i - 1 not in unpenalized
        rows.append(
            dict(
                term=name,
                coefficient=float(fit.coefficients[i]),
                standardized_coefficient=float(fit.standardized_coefficients[i]),
                penalized=penalty,
                selected=bool(penalty and abs(fit.standardized_coefficients[i]) > 1e-7),
                constant_in_training=False if i == 0 else not bool(fit.active[i - 1]),
                training_center=0.0 if i == 0 else float(fit.center[i - 1]),
                training_scale=1.0 if i == 0 else float(fit.scale[i - 1]),
                strength=fit.strength,
                l1_ratio=fit.l1_ratio,
                objective=fit.objective,
                projected_gradient=fit.projected_gradient,
                random_variance=fit.random_variance,
                dispersion=fit.dispersion,
                boundary_warning=fit.boundary_warning,
                inference_status="exploratory_no_post_selection_inference",
            )
        )
    return rows


def selection_main(args):
    predictor_text = (
        args.predictors
        if args.predictors is not None
        else ",".join(
            Path(args.predictor_file).read_text(encoding="utf-8").splitlines()
        )
    )
    names, free_names = _names(predictor_text), _names(args.unpenalized)
    if (
        not names
        or "(intercept)" in names
        or args.response in names
        or not set(free_names) <= set(names)
    ):
        raise ValueError(
            "Require predictors excluding the response and reserved (intercept) term; "
            "unpenalized columns must be predictors."
        )
    prefix = Path(args.out_prefix)
    bundle_paths = selection_output_paths(prefix)
    suffixes = tuple(bundle_paths)
    outputs = list(map(Path, bundle_paths.values()))
    inputs = [("tree", args.tree), ("data", args.data), ("folds", args.folds)]
    if args.predictor_file:
        inputs.append(("predictor_file", args.predictor_file))
    validate_outputs_do_not_replace_inputs(
        inputs, list(zip(suffixes, map(str, outputs), strict=True))
    )
    tree = read_tree(args.tree, "auto", True, rooted=args.input_rooted)
    if args.evolution_model == "brownian":
        require_rooted(
            tree, "Brownian predictor selection requires a rooted species tree."
        )
    leaves = [str(leaf.name) for leaf in tree.leaves()]
    if len(leaves) != len(set(leaves)) or any(not name for name in leaves):
        raise ValueError("Tree tips must have unique non-empty names.")
    if len(leaves) > 500:
        raise ValueError(
            "Exploratory elastic net currently supports at most 500 tips (dense covariance)."
        )
    data, folds = _read_table(args.data), _read_table(args.folds)
    if set(data.index) != set(leaves) or set(folds.index) != set(leaves):
        raise ValueError(
            "Tree, data and folds must contain exactly the same tips; subset explicitly upstream."
        )
    if "fold" not in folds:
        raise ValueError("Fold table requires a fold column.")
    x = data.loc[leaves, names].to_numpy(dtype=float)
    y = data.loc[leaves, args.response].to_numpy(dtype=float)
    groups = folds.loc[leaves, "fold"].to_numpy()
    covariance = build_phylogenetic_covariance(
        tree, leaves, evolution_model=args.evolution_model
    )
    unpenalized = tuple(names.index(name) for name in free_names)
    final, path, outer, cv, predictions = nested_selection(
        x,
        y,
        covariance,
        groups,
        family=args.family,
        strengths=tuple(map(float, args.strengths.split(","))),
        l1_ratios=tuple(map(float, args.l1_ratios.split(","))),
        unpenalized=unpenalized,
        prediction=args.prediction,
    )
    coefficients = pd.DataFrame(_coefficient_rows(final, names, unpenalized))
    path_table = pd.DataFrame(
        [row for fit in path for row in _coefficient_rows(fit, names, unpenalized)]
    )
    outer_rows = pd.DataFrame(
        [
            dict(row, fold=fold)
            for fold, fit in outer
            for row in _coefficient_rows(fit, names, unpenalized)
        ]
    )
    stability = (
        outer_rows.groupby("term", sort=False)
        .agg(
            outer_fit_count=("selected", "size"),
            selection_frequency=("selected", "mean"),
            penalized=("penalized", "first"),
        )
        .reset_index()
    )
    stability["interpretation"] = "descriptive_outer_fold_frequency_not_FDR_control"
    predictions.insert(0, "leaf_name", [leaves[i] for i in predictions.pop("row")])
    metadata = dict(
        response=args.response,
        response_family=args.family,
        evolution_model=args.evolution_model,
        prediction=args.prediction,
        prediction_estimand="plug_in_conditional_mode"
        if args.prediction == "conditional"
        else "inverse_link_fixed_effect_only",
        objective="GLS_squared_error_per_tip"
        if args.family == "gaussian"
        else "Laplace_negative_log_likelihood_per_tip",
        loss="binary_log_loss" if args.family == "binomial" else "squared_error",
        nested_cv_mean_loss=float(predictions.loss.mean()),
        nested_cv_baseline_mean_loss=float(predictions.baseline_loss.mean())
        if (predictions.baseline_status == "ok").all()
        else None,
        baseline="same_tree_and_unpenalized_covariates_without_selected_predictors",
        n_tips=len(leaves),
        n_predictors=len(names),
        n_folds=len(set(groups)),
        strengths=list(map(float, args.strengths.split(","))),
        l1_ratios=list(map(float, args.l1_ratios.split(","))),
        selected_strength=final.strength,
        selected_l1_ratio=final.l1_ratio,
        inference_status="exploratory_no_post_selection_inference",
        fold_policy="user_supplied_groups_do_not_imply_independence",
    )
    with output_transaction(outputs, create_parents=True) as staged:
        for target, table in zip(
            outputs[:-1],
            (coefficients, path_table, cv, predictions, stability),
            strict=True,
        ):
            table.to_csv(staged[target], sep="\t", index=False, na_rep="NA")
        Path(staged[outputs[-1]]).write_text(json.dumps(metadata, indent=2) + "\n")
    return 0
