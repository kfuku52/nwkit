"""Conventional tip-level phylogenetic generalized least squares."""

import math
import os
import warnings
from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import t as student_t

from nwkit.contrast import (
    _read_numeric_traits,
    _read_replicate_traits,
    _validate_replicate_options,
)
from nwkit.conventions import DEFAULT_TABLE_MISSING_VALUES_CSV
from nwkit.evolution import (
    EVOLUTION_MODELS,
    build_evolutionary_covariance,
    evolution_model_spec,
    optimization_parameterization,
    parameter_near_boundary,
    read_custom_covariance,
    validate_evolution_parameter,
)
from nwkit.measurement_error import (
    fit_conditional_eiv_gaussian,
    fit_latent_predictor,
)
from nwkit.pgls import _profile_covariance_fit, _solve_positive_definite
from nwkit.replicates import TIP_SUMMARY_COLUMNS
from nwkit.util import (
    is_rooted,
    normalized_missing_path_key,
    read_tip_table,
    read_tree,
    validate_distinct_output_paths,
    validate_unique_named_leaves,
)

ORDINARY_RESULT_COLUMNS = [
    "model_id",
    "tree_id",
    "response",
    "term",
    "coefficient",
    "standard_error",
    "statistic",
    "degrees_of_freedom",
    "p_value",
    "confidence_level",
    "confidence_interval_lower",
    "confidence_interval_upper",
    "n_species",
    "n_predictors",
    "num_parameters",
    "matrix_rank",
    "condition_number",
    "generalized_residual_sum_squares",
    "evolutionary_rate",
    "r_squared",
    "intercept",
    "evolution_model",
    "evolution_parameter_name",
    "evolution_parameter",
    "evolution_parameter_status",
    "branch_length_mode",
    "covariance_estimator",
    "inference_method",
    "reml",
    "mean_sampling_variance",
    "sampling_variance_fraction",
    "mean_predictor_sampling_variance",
    "mean_latent_predictor_variance",
    "predictor_uncertainty_fraction",
    "predictor_evolutionary_rate",
    "predictor_evolution_model",
    "predictor_evolution_parameter",
    "predictor_evolution_parameter_status",
    "predictor_evolution_optimizer_converged",
    "predictor_evolution_optimizer_message",
    "predictor_evolution_boundary_warning",
    "predictor_evolution_log_likelihood",
    "measurement_error_model",
    "log_likelihood",
    "optimizer_converged",
    "optimizer_message",
    "boundary_warning",
    "small_sample_warning",
    "inference_status",
    "model",
]

ORDINARY_MODEL_COMPARISON_COLUMNS = [
    "response",
    "evolution_model",
    "evolution_parameter_name",
    "evolution_parameter",
    "evolution_parameter_status",
    "branch_length_mode",
    "n_species",
    "n_coefficients",
    "n_likelihood_parameters",
    "log_likelihood",
    "aic",
    "delta_aic",
    "akaike_weight",
    "aicc",
    "delta_aicc",
    "aicc_weight",
    "bic",
    "optimizer_converged",
    "optimizer_message",
    "boundary_warning",
]

ORDINARY_SAMPLING_COVARIANCE_COLUMNS = [
    "tree_id",
    "trait",
    "leaf_name_1",
    "leaf_name_2",
    "sampling_covariance",
]


@dataclass
class OrdinaryPglsArtifacts:
    results: pd.DataFrame
    response_sampling_covariance: pd.DataFrame
    response_tip_summary: pd.DataFrame
    model_comparison: pd.DataFrame
    predictor_sampling_covariance: pd.DataFrame | None = None
    predictor_tip_summary: pd.DataFrame | None = None


def _validate_ordinary_tree(tree, branch_length: str) -> None:
    validate_unique_named_leaves(
        tree,
        option_name="--tree",
        context=" for ordinary PGLS",
    )
    if len(list(tree.leaves())) < 3:
        raise ValueError("'--tree' must contain at least three tips for PGLS.")
    if not is_rooted(tree) or len(tree.children) != 2:
        raise ValueError("'--tree' must be rooted with two root descendants.")
    if branch_length not in {"original", "unit"}:
        raise ValueError("Unsupported ordinary PGLS branch-length mode.")
    for node in tree.traverse():
        if node.is_root:
            continue
        if branch_length == "unit":
            continue
        try:
            distance = float(node.dist)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("PGLS tree branch lengths must be numeric.") from exc
        if not math.isfinite(distance) or distance <= 0.0:
            raise ValueError(
                "PGLS requires positive finite non-root branch lengths; "
                "use '--branch-length unit' to ignore input lengths."
            )


def build_phylogenetic_covariance(
    tree,
    leaf_names,
    *,
    evolution_model="brownian",
    parameter=None,
    branch_length="original",
    custom_covariance=None,
):
    """Build a tip covariance matrix in the requested tree-tip order."""
    return build_evolutionary_covariance(
        tree,
        leaf_names,
        model=evolution_model,
        parameter=parameter,
        branch_length=branch_length,
        custom_covariance=custom_covariance,
    )


def _bounded_scalar_minimize(function, bounds):
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="invalid value encountered in scalar subtract",
            category=RuntimeWarning,
            module=r"scipy\.optimize\._optimize",
        )
        return minimize_scalar(
            function,
            bounds=bounds,
            method="bounded",
            options={"xatol": 1e-7, "maxiter": 500},
        )


def _global_bounded_scalar_minimize(function, bounds, *, grid_size=17):
    """Search a bounded 1-D objective without assuming it is globally unimodal."""
    lower, upper = (float(bounds[0]), float(bounds[1]))
    grid = np.linspace(lower, upper, grid_size)
    grid_values = np.asarray([float(function(value)) for value in grid], dtype=float)
    finite_indices = np.flatnonzero(np.isfinite(grid_values))
    if not len(finite_indices):
        return SimpleNamespace(
            x=(lower + upper) / 2.0,
            fun=float("inf"),
            success=False,
            message="global grid search found no finite objective",
        )

    local_indices = []
    for index in finite_indices:
        if index == 0 or index == grid_size - 1:
            continue
        left = grid_values[index - 1]
        right = grid_values[index + 1]
        if grid_values[index] <= left and grid_values[index] <= right:
            local_indices.append(int(index))
    ranked_indices = sorted(
        (int(index) for index in finite_indices if 0 < index < grid_size - 1),
        key=lambda index: grid_values[index],
    )
    refinement_indices = []
    for index in local_indices + ranked_indices:
        if index not in refinement_indices:
            refinement_indices.append(index)
        if len(refinement_indices) == 4:
            break

    candidates = [
        SimpleNamespace(
            x=float(grid[index]),
            fun=float(grid_values[index]),
            success=index in {0, grid_size - 1},
            message="global grid point",
        )
        for index in finite_indices
    ]
    successful_refinements = 0
    for index in refinement_indices:
        refined = _bounded_scalar_minimize(
            function,
            (float(grid[index - 1]), float(grid[index + 1])),
        )
        if math.isfinite(float(refined.fun)):
            candidates.append(refined)
            successful_refinements += int(bool(refined.success))
    best = min(candidates, key=lambda candidate: float(candidate.fun))
    best_is_boundary = math.isclose(float(best.x), lower) or math.isclose(
        float(best.x), upper
    )
    return SimpleNamespace(
        x=float(best.x),
        fun=float(best.fun),
        success=bool(best_is_boundary or successful_refinements),
        message=(
            "global grid search ({} points; {} successful local refinement(s))"
        ).format(grid_size, successful_refinements),
    )


def _fit_ordinary_gaussian(
    y,
    design,
    fixed_covariance,
    tree,
    leaf_names,
    *,
    evolution_model,
    evolution_parameter,
    branch_length,
    custom_covariance,
    reml,
    predictor_uncertainties=(),
    predictor_columns=(),
):
    parameter_status = "not-applicable"
    outer_converged = True
    outer_message = "fixed evolutionary model"
    spec = evolution_model_spec(evolution_model)

    def fit_at(parameter):
        phylogenetic_covariance = build_phylogenetic_covariance(
            tree,
            leaf_names,
            evolution_model=evolution_model,
            parameter=parameter,
            branch_length=branch_length,
            custom_covariance=custom_covariance,
        )
        components = [("evolutionary_rate", phylogenetic_covariance)]
        if predictor_uncertainties:
            fit = fit_conditional_eiv_gaussian(
                y,
                design,
                predictor_uncertainties,
                predictor_columns,
                fixed_covariance,
                components,
                reml=reml,
            )
        else:
            fit = _profile_covariance_fit(
                y,
                design,
                fixed_covariance,
                components,
                reml=reml,
            )
        fit["phylogenetic_covariance"] = phylogenetic_covariance
        fit["evolution_parameter"] = parameter
        return fit

    if spec.parameter_name is None:
        fit = fit_at(None)
    elif evolution_parameter is not None:
        parameter_status = "fixed"
        fit = fit_at(validate_evolution_parameter(evolution_model, evolution_parameter))
    else:
        parameter_status = "estimated"
        parameter_bounds, decode = optimization_parameterization(
            tree,
            evolution_model,
            branch_length=branch_length,
        )

        cache = {}

        def cached_fit(value):
            decoded = decode(value)
            key = float(decoded)
            if key not in cache:
                try:
                    cache[key] = fit_at(decoded)
                except ValueError:
                    cache[key] = None
            return cache[key]

        def objective(value):
            candidate = cached_fit(value)
            return float("inf") if candidate is None else float(candidate["objective"])

        optimized = _global_bounded_scalar_minimize(objective, parameter_bounds)
        candidate_values = [parameter_bounds[0], parameter_bounds[1]]
        if math.isfinite(float(optimized.fun)):
            candidate_values.append(float(optimized.x))
        candidates = [cached_fit(value) for value in candidate_values]
        candidates = [candidate for candidate in candidates if candidate is not None]
        if not candidates:
            raise ValueError(
                "Evolution-parameter optimization failed to find a finite PGLS fit."
            )
        fit = min(candidates, key=lambda candidate: float(candidate["objective"]))
        outer_converged = bool(optimized.success)
        outer_message = str(optimized.message)
    parameter = fit["evolution_parameter"]
    parameter_boundary = False
    if parameter_status == "estimated" and parameter is not None:
        parameter_boundary = parameter_near_boundary(
            tree,
            evolution_model,
            float(parameter),
            branch_length=branch_length,
        )
    fit["evolution_parameter_status"] = parameter_status
    fit["optimizer_converged"] = bool(fit["optimizer_converged"] and outer_converged)
    fit["optimizer_message"] = "{}; {}".format(
        outer_message,
        fit["optimizer_message"],
    )
    fit["boundary_warning"] = bool(fit["boundary_warning"] or parameter_boundary)
    return fit


def _ordinary_bootstrap_coefficients(
    fit,
    design,
    fixed_covariance,
    tree,
    leaf_names,
    *,
    evolution_model,
    evolution_parameter,
    branch_length,
    custom_covariance,
    reml,
    replicates,
    seed,
    predictor_uncertainties=(),
    predictor_columns=(),
):
    rng = np.random.default_rng(seed)
    coefficients: list[np.ndarray] = []
    mean = design @ fit["beta"]
    maximum_attempts = max(replicates * 3, replicates + 10)
    attempts = 0
    while len(coefficients) < replicates and attempts < maximum_attempts:
        attempts += 1
        simulated = mean + fit["cholesky"] @ rng.standard_normal(len(mean))
        try:
            bootstrap_fit = _fit_ordinary_gaussian(
                simulated,
                design,
                fixed_covariance,
                tree,
                leaf_names,
                evolution_model=evolution_model,
                evolution_parameter=evolution_parameter,
                branch_length=branch_length,
                custom_covariance=custom_covariance,
                reml=reml,
                predictor_uncertainties=predictor_uncertainties,
                predictor_columns=predictor_columns,
            )
        except ValueError:
            continue
        coefficients.append(bootstrap_fit["beta"])
    if len(coefficients) < replicates:
        raise ValueError(
            "Parametric bootstrap produced only {} successful ordinary PGLS fits in {} attempts.".format(
                len(coefficients), attempts
            )
        )
    return np.asarray(coefficients, dtype=float)


def _ordered_values(values_by_leaf, leaf_names, trait):
    missing = sorted(set(leaf_names) - set(values_by_leaf))
    if missing:
        raise ValueError(
            "Trait '{}' is missing tree tips: {}.".format(trait, ", ".join(missing))
        )
    try:
        values = np.asarray([values_by_leaf[name] for name in leaf_names], dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "Trait '{}' must contain numeric values.".format(trait)
        ) from exc
    if not np.isfinite(values).all():
        raise ValueError("Trait '{}' must contain finite values.".format(trait))
    return values


def estimate_marginal_evolution_parameter(
    tree,
    values_by_leaf,
    trait,
    *,
    evolution_model,
    evolution_parameter=None,
    branch_length="original",
    sampling_covariance=None,
):
    """Estimate one trait's evolutionary shape parameter by marginal ML."""
    spec = evolution_model_spec(evolution_model)
    if not spec.contrast_supported:
        raise ValueError(
            "Evolutionary model '{}' cannot produce reconciled contrasts.".format(
                evolution_model
            )
        )
    _validate_ordinary_tree(
        tree,
        branch_length if spec.branch_lengths_used else "unit",
    )
    leaf_names = [str(leaf.name) for leaf in tree.leaves()]
    response = _ordered_values(values_by_leaf, leaf_names, trait)
    design = np.ones((len(leaf_names), 1), dtype=float)
    fixed_covariance = (
        np.zeros((len(leaf_names), len(leaf_names)), dtype=float)
        if sampling_covariance is None
        else _coerce_response_sampling_covariance(
            {trait: sampling_covariance},
            trait,
            leaf_names,
            label="Predictor sampling covariance",
        )
    )
    fit = _fit_ordinary_gaussian(
        response,
        design,
        fixed_covariance,
        tree,
        leaf_names,
        evolution_model=evolution_model,
        evolution_parameter=evolution_parameter,
        branch_length=branch_length,
        custom_covariance=None,
        reml=False,
    )
    return {
        "parameter": fit["evolution_parameter"],
        "parameter_status": fit["evolution_parameter_status"],
        "log_likelihood": -float(fit["objective"]),
        "optimizer_converged": bool(fit["optimizer_converged"]),
        "optimizer_message": str(fit["optimizer_message"]),
        "boundary_warning": bool(fit["boundary_warning"]),
    }


def _validate_ordinary_fit_settings(
    tree,
    *,
    evolution_model,
    evolution_parameter,
    branch_length,
    custom_covariance,
    intercept,
    confidence_level,
    inference,
    bootstrap_replicates,
    seed,
    reml,
):
    if branch_length not in {"original", "unit"}:
        raise ValueError("Unsupported ordinary PGLS branch-length mode.")
    spec = evolution_model_spec(evolution_model)
    _validate_ordinary_tree(
        tree,
        branch_length if spec.branch_lengths_used else "unit",
    )
    if evolution_parameter is not None:
        validate_evolution_parameter(evolution_model, evolution_parameter)
    if evolution_model == "custom" and custom_covariance is None:
        raise ValueError("Custom evolution model requires a covariance matrix.")
    if evolution_model != "custom" and custom_covariance is not None:
        raise ValueError("A custom covariance is only valid with model 'custom'.")
    if spec.parameter_name is None and evolution_parameter is not None:
        raise ValueError(
            "Evolutionary model '{}' does not take a parameter.".format(evolution_model)
        )
    if not 0.0 < confidence_level < 1.0:
        raise ValueError("confidence_level must be between zero and one.")
    if inference not in {"wald", "parametric-bootstrap"}:
        raise ValueError("Unsupported inference method: {}.".format(inference))
    if not isinstance(intercept, bool) or not isinstance(reml, bool):
        raise ValueError("intercept and reml must be booleans.")
    if (
        not isinstance(bootstrap_replicates, int)
        or isinstance(bootstrap_replicates, bool)
        or bootstrap_replicates < 2
    ):
        raise ValueError("bootstrap_replicates must be an integer of at least two.")
    if not isinstance(seed, int) or isinstance(seed, bool) or seed < 0:
        raise ValueError("Bootstrap seed must be a non-negative integer.")


def _normalize_ordinary_traits(
    predictor_values_by_trait,
    responses,
    predictors,
):
    responses = list(responses)
    predictors = list(predictors)
    if not responses or len(responses) != len(set(responses)):
        raise ValueError("responses must be a non-empty sequence of unique names.")
    if not predictors or len(predictors) != len(set(predictors)):
        raise ValueError("predictors must be a non-empty sequence of unique names.")
    missing_predictors = [
        predictor
        for predictor in predictors
        if predictor not in predictor_values_by_trait
    ]
    if missing_predictors:
        raise ValueError(
            "Predictor trait(s) are absent: {}.".format(", ".join(missing_predictors))
        )
    return responses, predictors


def _ordinary_design_matrix(
    predictor_values_by_trait,
    predictors,
    leaf_names,
    *,
    intercept,
):
    predictor_columns = [
        _ordered_values(predictor_values_by_trait[predictor], leaf_names, predictor)
        for predictor in predictors
    ]
    design = np.column_stack(predictor_columns)
    term_names = list(predictors)
    if intercept:
        design = np.column_stack([np.ones(len(leaf_names)), design])
        term_names = ["(intercept)"] + term_names
    num_parameters = design.shape[1]
    if len(leaf_names) <= num_parameters:
        raise ValueError(
            "Ordinary PGLS needs more species than fitted coefficients "
            "(species={}; coefficients={}).".format(len(leaf_names), num_parameters)
        )
    matrix_rank = int(np.linalg.matrix_rank(design))
    if matrix_rank != num_parameters:
        raise ValueError("Ordinary PGLS design matrix is rank deficient.")
    return design, term_names, num_parameters, matrix_rank


def _prepare_ordinary_design(
    tree,
    predictor_values_by_trait,
    responses,
    predictors,
    *,
    evolution_model,
    evolution_parameter,
    branch_length,
    custom_covariance,
    intercept,
    confidence_level,
    inference,
    bootstrap_replicates,
    seed,
    reml,
):
    _validate_ordinary_fit_settings(
        tree,
        evolution_model=evolution_model,
        evolution_parameter=evolution_parameter,
        branch_length=branch_length,
        custom_covariance=custom_covariance,
        intercept=intercept,
        confidence_level=confidence_level,
        inference=inference,
        bootstrap_replicates=bootstrap_replicates,
        seed=seed,
        reml=reml,
    )
    responses, predictors = _normalize_ordinary_traits(
        predictor_values_by_trait,
        responses,
        predictors,
    )
    leaf_names = [str(leaf.name) for leaf in tree.leaves()]
    design, term_names, num_parameters, matrix_rank = _ordinary_design_matrix(
        predictor_values_by_trait,
        predictors,
        leaf_names,
        intercept=intercept,
    )
    return (
        responses,
        predictors,
        leaf_names,
        design,
        term_names,
        num_parameters,
        matrix_rank,
    )


def _coerce_response_sampling_covariance(
    covariance_by_trait,
    response,
    leaf_names,
    *,
    label="Response sampling covariance",
):
    covariance = covariance_by_trait.get(response)
    if covariance is None:
        covariance = np.zeros((len(leaf_names), len(leaf_names)))
    elif isinstance(covariance, pd.DataFrame):
        covariance = covariance.copy()
        covariance.index = covariance.index.map(str)
        covariance.columns = covariance.columns.map(str)
        if covariance.index.duplicated().any() or covariance.columns.duplicated().any():
            raise ValueError("{} contains duplicated tree-tip names.".format(label))
        missing_rows = sorted(set(leaf_names) - set(covariance.index))
        missing_columns = sorted(set(leaf_names) - set(covariance.columns))
        extra_rows = sorted(set(covariance.index) - set(leaf_names))
        extra_columns = sorted(set(covariance.columns) - set(leaf_names))
        if missing_rows or missing_columns or extra_rows or extra_columns:
            raise ValueError("{} species must exactly match tree tips.".format(label))
        covariance = covariance.loc[leaf_names, leaf_names].to_numpy(dtype=float)
    else:
        covariance = np.asarray(covariance, dtype=float)
    if covariance.shape != (len(leaf_names), len(leaf_names)):
        raise ValueError("{} has the wrong dimensions.".format(label))
    return covariance


def _ordinary_inference_samples(
    fit,
    design,
    fixed_covariance,
    tree,
    leaf_names,
    *,
    evolution_model,
    evolution_parameter,
    branch_length,
    custom_covariance,
    inference,
    bootstrap_replicates,
    seed,
    reml,
    predictor_uncertainties=(),
    predictor_columns=(),
):
    if inference != "parametric-bootstrap":
        standard_errors = np.sqrt(np.maximum(np.diag(fit["beta_covariance"]), 0.0))
        return None, standard_errors
    coefficients = _ordinary_bootstrap_coefficients(
        fit,
        design,
        fixed_covariance,
        tree,
        leaf_names,
        evolution_model=evolution_model,
        evolution_parameter=evolution_parameter,
        branch_length=branch_length,
        custom_covariance=custom_covariance,
        reml=reml,
        replicates=bootstrap_replicates,
        seed=seed,
        predictor_uncertainties=predictor_uncertainties,
        predictor_columns=predictor_columns,
    )
    return coefficients, np.std(coefficients, axis=0, ddof=1)


def _ordinary_response_statistics(y, fit, fixed_covariance, *, intercept):
    if intercept:
        null_design = np.ones((len(y), 1))
        inverse_null = _solve_positive_definite(fit["cholesky"], null_design)
        null_beta = np.linalg.solve(
            null_design.T @ inverse_null,
            null_design.T @ _solve_positive_definite(fit["cholesky"], y),
        )
        null_residual = y - null_design @ null_beta
    else:
        null_residual = y
    total_quadratic = float(
        null_residual @ _solve_positive_definite(fit["cholesky"], null_residual)
    )
    r_squared = (
        float("nan")
        if total_quadratic == 0.0
        else 1.0 - float(fit["quadratic"]) / total_quadratic
    )
    evolutionary_rate = float(fit["component_variances"]["evolutionary_rate"])
    mean_sampling_variance = float(np.mean(np.diag(fixed_covariance)))
    mean_fitted_variance = float(np.mean(np.diag(fit["covariance"])))
    sampling_fraction = (
        0.0
        if mean_fitted_variance == 0.0
        else mean_sampling_variance / mean_fitted_variance
    )
    return {
        "condition_number": float(np.linalg.cond(fit["beta_covariance"])),
        "evolutionary_rate": evolutionary_rate,
        "mean_sampling_variance": mean_sampling_variance,
        "r_squared": r_squared,
        "sampling_fraction": sampling_fraction,
    }


def _ordinary_coefficient_rows(
    response,
    fit,
    standard_errors,
    bootstrap_coefficients,
    term_names,
    statistics,
    *,
    confidence_level,
    degrees_of_freedom,
    n_species,
    n_predictors,
    num_parameters,
    matrix_rank,
    intercept,
    evolution_model,
    branch_length,
    inference,
    reml,
    predictor_diagnostics,
    mean_predictor_sampling_variance,
):
    critical = float(student_t.ppf(0.5 + confidence_level / 2.0, degrees_of_freedom))
    rows = []
    for term_index, term in enumerate(term_names):
        coefficient = float(fit["beta"][term_index])
        standard_error = float(standard_errors[term_index])
        if bootstrap_coefficients is None:
            lower = coefficient - critical * standard_error
            upper = coefficient + critical * standard_error
        else:
            alpha = (1.0 - confidence_level) / 2.0
            lower, upper = np.quantile(
                bootstrap_coefficients[:, term_index],
                [alpha, 1.0 - alpha],
            )
            lower = float(lower)
            upper = float(upper)
        if standard_error == 0.0:
            statistic: float | str = ""
            p_value: float | str = ""
            inference_status = "zero-model-variance"
        else:
            statistic_value = coefficient / standard_error
            statistic = statistic_value
            if bootstrap_coefficients is None:
                p_value = float(
                    2.0 * student_t.sf(abs(statistic_value), degrees_of_freedom)
                )
            else:
                centered = bootstrap_coefficients[:, term_index] - coefficient
                p_value = float(
                    (1 + np.sum(np.abs(centered) >= abs(coefficient)))
                    / (len(centered) + 1)
                )
            inference_status = "ok"
        parameter = fit["evolution_parameter"]
        spec = evolution_model_spec(evolution_model)
        predictor_diagnostic = predictor_diagnostics.get(term)
        rows.append(
            {
                "model_id": "ordinary:{}".format(response),
                "tree_id": "species",
                "response": response,
                "term": term,
                "coefficient": coefficient,
                "standard_error": standard_error,
                "statistic": statistic,
                "degrees_of_freedom": degrees_of_freedom,
                "p_value": p_value,
                "confidence_level": confidence_level,
                "confidence_interval_lower": lower,
                "confidence_interval_upper": upper,
                "n_species": n_species,
                "n_predictors": n_predictors,
                "num_parameters": num_parameters,
                "matrix_rank": matrix_rank,
                "condition_number": statistics["condition_number"],
                "generalized_residual_sum_squares": float(fit["quadratic"]),
                "evolutionary_rate": statistics["evolutionary_rate"],
                "r_squared": statistics["r_squared"],
                "intercept": "yes" if intercept else "no",
                "evolution_model": evolution_model,
                "evolution_parameter_name": spec.parameter_name or "",
                "evolution_parameter": parameter if parameter is not None else "",
                "evolution_parameter_status": fit["evolution_parameter_status"],
                "branch_length_mode": (
                    branch_length if spec.branch_lengths_used else "not-applicable"
                ),
                "covariance_estimator": "gaussian{}-{}".format(
                    "-eiv" if predictor_diagnostics else "",
                    "REML" if reml else "ML",
                ),
                "inference_method": inference,
                "reml": "yes" if reml else "no",
                "mean_sampling_variance": statistics["mean_sampling_variance"],
                "sampling_variance_fraction": statistics["sampling_fraction"],
                "mean_predictor_sampling_variance": (
                    mean_predictor_sampling_variance.get(term, "")
                ),
                "mean_latent_predictor_variance": (
                    ""
                    if predictor_diagnostic is None
                    else predictor_diagnostic["mean_posterior_variance"]
                ),
                "predictor_uncertainty_fraction": (
                    ""
                    if predictor_diagnostic is None
                    else predictor_diagnostic["uncertainty_fraction"]
                ),
                "predictor_evolutionary_rate": (
                    ""
                    if predictor_diagnostic is None
                    else predictor_diagnostic["evolutionary_rate"]
                ),
                "predictor_evolution_model": (
                    ""
                    if predictor_diagnostic is None
                    else predictor_diagnostic["model"]
                ),
                "predictor_evolution_parameter": (
                    ""
                    if predictor_diagnostic is None
                    or predictor_diagnostic["parameter"] is None
                    else predictor_diagnostic["parameter"]
                ),
                "predictor_evolution_parameter_status": (
                    ""
                    if predictor_diagnostic is None
                    else predictor_diagnostic["parameter_status"]
                ),
                "predictor_evolution_optimizer_converged": (
                    ""
                    if predictor_diagnostic is None
                    else "yes"
                    if predictor_diagnostic["optimizer_converged"]
                    else "no"
                ),
                "predictor_evolution_optimizer_message": (
                    ""
                    if predictor_diagnostic is None
                    else predictor_diagnostic["optimizer_message"]
                ),
                "predictor_evolution_boundary_warning": (
                    ""
                    if predictor_diagnostic is None
                    else "yes"
                    if predictor_diagnostic["boundary_warning"]
                    else "no"
                ),
                "predictor_evolution_log_likelihood": (
                    ""
                    if predictor_diagnostic is None
                    else predictor_diagnostic["log_likelihood"]
                ),
                "measurement_error_model": (
                    "latent-predictor"
                    if predictor_diagnostics
                    else "response-only"
                    if statistics["mean_sampling_variance"] > 0.0
                    else "none"
                ),
                "log_likelihood": -float(fit["objective"]),
                "optimizer_converged": ("yes" if fit["optimizer_converged"] else "no"),
                "optimizer_message": fit["optimizer_message"],
                "boundary_warning": "yes" if fit["boundary_warning"] else "no",
                "small_sample_warning": "yes" if n_species < 20 else "no",
                "inference_status": inference_status,
                "model": "ordinary",
            }
        )
    return rows


def _prepare_latent_ordinary_predictors(
    tree,
    predictor_values_by_trait,
    predictors,
    predictor_sampling_covariance,
    *,
    evolution_model,
    evolution_parameter,
    branch_length,
    custom_covariance,
):
    leaf_names = [str(leaf.name) for leaf in tree.leaves()]
    predictor_covariance_by_trait = predictor_sampling_covariance or {}
    if not predictor_covariance_by_trait:
        return predictor_values_by_trait, [], {}
    missing_covariance = [
        predictor
        for predictor in predictors
        if predictor not in predictor_covariance_by_trait
    ]
    if missing_covariance:
        raise ValueError(
            "Predictor sampling covariance is missing trait(s): {}.".format(
                ", ".join(missing_covariance)
            )
        )
    posterior_predictor_values = dict(predictor_values_by_trait)
    predictor_uncertainties = []
    diagnostics = {}
    for predictor in predictors:
        observed = _ordered_values(
            predictor_values_by_trait[predictor], leaf_names, predictor
        )
        sampling = _coerce_response_sampling_covariance(
            predictor_covariance_by_trait,
            predictor,
            leaf_names,
            label="Predictor sampling covariance",
        )
        marginal_fit = _fit_ordinary_gaussian(
            observed,
            np.ones((len(leaf_names), 1), dtype=float),
            sampling,
            tree,
            leaf_names,
            evolution_model=evolution_model,
            evolution_parameter=evolution_parameter,
            branch_length=branch_length,
            custom_covariance=(
                custom_covariance if evolution_model == "custom" else None
            ),
            reml=False,
        )
        posterior = fit_latent_predictor(
            observed,
            marginal_fit["phylogenetic_covariance"],
            sampling,
            include_intercept=True,
        )
        posterior_predictor_values[predictor] = dict(
            zip(leaf_names, posterior.mean, strict=True)
        )
        predictor_uncertainties.append(posterior.covariance)
        diagnostics[predictor] = {
            "evolutionary_rate": posterior.evolutionary_rate,
            "model": evolution_model,
            "parameter": marginal_fit["evolution_parameter"],
            "parameter_status": marginal_fit["evolution_parameter_status"],
            "optimizer_converged": bool(
                posterior.optimizer_converged and marginal_fit["optimizer_converged"]
            ),
            "optimizer_message": "{}; {}".format(
                marginal_fit["optimizer_message"],
                posterior.optimizer_message,
            ),
            "boundary_warning": bool(
                posterior.boundary_warning or marginal_fit["boundary_warning"]
            ),
            "log_likelihood": posterior.log_likelihood,
            "mean_sampling_variance": float(np.mean(np.diag(sampling))),
            "mean_posterior_variance": float(np.mean(np.diag(posterior.covariance))),
            "uncertainty_fraction": float(
                np.mean(np.diag(posterior.covariance))
                / np.mean(
                    np.diag(
                        posterior.evolutionary_rate
                        * marginal_fit["phylogenetic_covariance"]
                    )
                )
            ),
        }
    return posterior_predictor_values, predictor_uncertainties, diagnostics


def fit_ordinary_pgls(
    tree,
    response_values_by_trait,
    predictor_values_by_trait,
    responses,
    predictors,
    *,
    response_sampling_covariance=None,
    predictor_sampling_covariance=None,
    evolution_model="brownian",
    evolution_parameter=None,
    branch_length="original",
    custom_covariance=None,
    intercept=True,
    confidence_level=0.95,
    inference="wald",
    bootstrap_replicates=1000,
    seed=1,
    reml=True,
    predictor_evolution_model=None,
    predictor_evolution_parameter=None,
    predictor_branch_length=None,
):
    """Fit conventional tip-level PGLS models, one per response trait."""
    leaf_names = [str(leaf.name) for leaf in tree.leaves()]
    predictor_model = predictor_evolution_model or evolution_model
    predictor_length = predictor_branch_length or branch_length
    (
        posterior_predictor_values,
        predictor_uncertainties,
        predictor_diagnostics,
    ) = _prepare_latent_ordinary_predictors(
        tree,
        predictor_values_by_trait,
        predictors,
        predictor_sampling_covariance,
        evolution_model=predictor_model,
        evolution_parameter=predictor_evolution_parameter,
        branch_length=predictor_length,
        custom_covariance=custom_covariance,
    )
    (
        responses,
        predictors,
        leaf_names,
        design,
        term_names,
        num_parameters,
        matrix_rank,
    ) = _prepare_ordinary_design(
        tree,
        posterior_predictor_values,
        responses,
        predictors,
        evolution_model=evolution_model,
        evolution_parameter=evolution_parameter,
        branch_length=branch_length,
        custom_covariance=custom_covariance,
        intercept=intercept,
        confidence_level=confidence_level,
        inference=inference,
        bootstrap_replicates=bootstrap_replicates,
        seed=seed,
        reml=reml,
    )
    covariance_by_trait = response_sampling_covariance or {}
    predictor_columns = list(
        range(1 if intercept else 0, (1 if intercept else 0) + len(predictors))
    )
    mean_predictor_sampling_variance = {
        predictor: diagnostics["mean_sampling_variance"]
        for predictor, diagnostics in predictor_diagnostics.items()
    }
    rows = []
    for response_index, response in enumerate(responses):
        if response not in response_values_by_trait:
            raise ValueError("Response trait '{}' is absent.".format(response))
        y = _ordered_values(response_values_by_trait[response], leaf_names, response)
        fixed_covariance = _coerce_response_sampling_covariance(
            covariance_by_trait,
            response,
            leaf_names,
        )
        fit = _fit_ordinary_gaussian(
            y,
            design,
            fixed_covariance,
            tree,
            leaf_names,
            evolution_model=evolution_model,
            evolution_parameter=evolution_parameter,
            branch_length=branch_length,
            custom_covariance=custom_covariance,
            reml=reml,
            predictor_uncertainties=predictor_uncertainties,
            predictor_columns=predictor_columns,
        )
        bootstrap_coefficients, standard_errors = _ordinary_inference_samples(
            fit,
            design,
            fixed_covariance,
            tree,
            leaf_names,
            evolution_model=evolution_model,
            evolution_parameter=evolution_parameter,
            branch_length=branch_length,
            custom_covariance=custom_covariance,
            inference=inference,
            bootstrap_replicates=bootstrap_replicates,
            seed=seed + response_index,
            reml=reml,
            predictor_uncertainties=predictor_uncertainties,
            predictor_columns=predictor_columns,
        )
        degrees_of_freedom = len(leaf_names) - num_parameters
        statistics = _ordinary_response_statistics(
            y,
            fit,
            fixed_covariance,
            intercept=intercept,
        )
        rows.extend(
            _ordinary_coefficient_rows(
                response,
                fit,
                standard_errors,
                bootstrap_coefficients,
                term_names,
                statistics,
                confidence_level=confidence_level,
                degrees_of_freedom=degrees_of_freedom,
                n_species=len(leaf_names),
                n_predictors=len(predictors),
                num_parameters=num_parameters,
                matrix_rank=matrix_rank,
                intercept=intercept,
                evolution_model=evolution_model,
                branch_length=branch_length,
                inference=inference,
                reml=reml,
                predictor_diagnostics=predictor_diagnostics,
                mean_predictor_sampling_variance=mean_predictor_sampling_variance,
            )
        )
    return pd.DataFrame(rows, columns=ORDINARY_RESULT_COLUMNS)


def _parse_comparison_models(value) -> list[str]:
    if value in (None, ""):
        return []
    models = [item.strip() for item in str(value).split(",")]
    if any(model == "" for model in models):
        raise ValueError("'--compare-evolution-models' contains an empty model.")
    if len(models) != len(set(models)):
        raise ValueError("'--compare-evolution-models' contains duplicated models.")
    unsupported = [model for model in models if model not in EVOLUTION_MODELS]
    if unsupported:
        raise ValueError(
            "'--compare-evolution-models' contains unsupported model(s): {}.".format(
                ", ".join(unsupported)
            )
        )
    return models


def fit_ordinary_model_comparison(
    tree,
    response_values_by_trait,
    predictor_values_by_trait,
    responses,
    predictors,
    evolution_models,
    *,
    response_sampling_covariance=None,
    predictor_sampling_covariance=None,
    branch_length="original",
    custom_covariance=None,
    intercept=True,
    predictor_evolution_model="brownian",
    predictor_evolution_parameter=None,
    predictor_branch_length=None,
):
    """Compare evolutionary covariance models using maximum likelihood."""
    if isinstance(evolution_models, (str, bytes)):
        raise ValueError("evolution_models must be a non-empty unique sequence.")
    models = list(evolution_models)
    if not models or len(models) != len(set(models)):
        raise ValueError("evolution_models must be a non-empty unique sequence.")
    specs = [evolution_model_spec(model) for model in models]
    validation_model = (
        "brownian" if any(spec.branch_lengths_used for spec in specs) else "independent"
    )
    (
        predictor_values_by_trait,
        predictor_uncertainties,
        _,
    ) = _prepare_latent_ordinary_predictors(
        tree,
        predictor_values_by_trait,
        predictors,
        predictor_sampling_covariance,
        evolution_model=predictor_evolution_model,
        evolution_parameter=predictor_evolution_parameter,
        branch_length=predictor_branch_length or branch_length,
        custom_covariance=custom_covariance,
    )
    (
        responses,
        predictors,
        leaf_names,
        design,
        _,
        num_coefficients,
        _,
    ) = _prepare_ordinary_design(
        tree,
        predictor_values_by_trait,
        responses,
        predictors,
        evolution_model=validation_model,
        evolution_parameter=None,
        branch_length=branch_length,
        custom_covariance=None,
        intercept=intercept,
        confidence_level=0.95,
        inference="wald",
        bootstrap_replicates=2,
        seed=0,
        reml=False,
    )
    if "custom" in models and custom_covariance is None:
        raise ValueError(
            "Model comparison containing 'custom' requires '--evolution-covariance'."
        )
    if (
        "custom" not in models
        and predictor_evolution_model != "custom"
        and custom_covariance is not None
    ):
        raise ValueError(
            "A custom covariance is only valid when model comparison includes 'custom'."
        )
    covariance_by_trait = response_sampling_covariance or {}
    predictor_columns = list(
        range(1 if intercept else 0, (1 if intercept else 0) + len(predictors))
    )
    rows = []
    for response in responses:
        if response not in response_values_by_trait:
            raise ValueError("Response trait '{}' is absent.".format(response))
        y = _ordered_values(response_values_by_trait[response], leaf_names, response)
        fixed_covariance = _coerce_response_sampling_covariance(
            covariance_by_trait,
            response,
            leaf_names,
        )
        for model in models:
            fit = _fit_ordinary_gaussian(
                y,
                design,
                fixed_covariance,
                tree,
                leaf_names,
                evolution_model=model,
                evolution_parameter=None,
                branch_length=branch_length,
                custom_covariance=(custom_covariance if model == "custom" else None),
                reml=False,
                predictor_uncertainties=predictor_uncertainties,
                predictor_columns=predictor_columns,
            )
            parameter = fit["evolution_parameter"]
            parameter_count = 1 if evolution_model_spec(model).parameter_name else 0
            likelihood_parameters = num_coefficients + 1 + parameter_count
            log_likelihood = -float(fit["objective"])
            aic = 2.0 * likelihood_parameters - 2.0 * log_likelihood
            denominator = len(leaf_names) - likelihood_parameters - 1
            aicc = (
                aic
                + 2.0
                * likelihood_parameters
                * (likelihood_parameters + 1)
                / denominator
                if denominator > 0
                else float("nan")
            )
            bic = (
                math.log(len(leaf_names)) * likelihood_parameters - 2.0 * log_likelihood
            )
            rows.append(
                {
                    "response": response,
                    "evolution_model": model,
                    "evolution_parameter_name": evolution_model_spec(
                        model
                    ).parameter_name
                    or "",
                    "evolution_parameter": parameter if parameter is not None else "",
                    "evolution_parameter_status": fit["evolution_parameter_status"],
                    "branch_length_mode": (
                        branch_length
                        if evolution_model_spec(model).branch_lengths_used
                        else "not-applicable"
                    ),
                    "n_species": len(leaf_names),
                    "n_coefficients": num_coefficients,
                    "n_likelihood_parameters": likelihood_parameters,
                    "log_likelihood": log_likelihood,
                    "aic": aic,
                    "delta_aic": "",
                    "akaike_weight": "",
                    "aicc": aicc,
                    "delta_aicc": "",
                    "aicc_weight": "",
                    "bic": bic,
                    "optimizer_converged": (
                        "yes" if fit["optimizer_converged"] else "no"
                    ),
                    "optimizer_message": fit["optimizer_message"],
                    "boundary_warning": "yes" if fit["boundary_warning"] else "no",
                }
            )
    result = pd.DataFrame(rows, columns=ORDINARY_MODEL_COMPARISON_COLUMNS)
    for column in ["delta_aic", "akaike_weight", "delta_aicc", "aicc_weight"]:
        result[column] = result[column].astype(object)
    for _response, response_rows in result.groupby("response", sort=False):
        aic_values = response_rows["aic"].to_numpy(dtype=float)
        delta_aic = aic_values - float(np.min(aic_values))
        aic_relative = np.exp(-0.5 * delta_aic)
        result.loc[response_rows.index, "delta_aic"] = delta_aic
        result.loc[response_rows.index, "akaike_weight"] = aic_relative / float(
            np.sum(aic_relative)
        )
        aicc_values = response_rows["aicc"].to_numpy(dtype=float)
        if not np.isfinite(aicc_values).all():
            continue
        delta_aicc = aicc_values - float(np.min(aicc_values))
        aicc_relative = np.exp(-0.5 * delta_aicc)
        result.loc[response_rows.index, "delta_aicc"] = delta_aicc
        result.loc[response_rows.index, "aicc_weight"] = aicc_relative / float(
            np.sum(aicc_relative)
        )
    return result


def _effective_ordinary_args(args) -> SimpleNamespace:
    values = vars(args).copy()
    defaults = {
        "batch": None,
        "biological_id": None,
        "branch_length": "original",
        "evolution_model": "brownian",
        "intercept": True,
        "missing_values": DEFAULT_TABLE_MISSING_VALUES_CSV,
        "predictor_batch": None,
        "predictor_biological_id": None,
        "predictor_sample_size_columns": None,
        "predictor_standard_error_columns": None,
        "predictor_technical_aggregation": "error",
        "predictor_technical_id": None,
        "predictor_within_variance": "pooled",
        "quoted_node_names": True,
        "sample_size_columns": None,
        "standard_error_columns": None,
        "technical_aggregation": "error",
        "technical_id": None,
        "tree_format": "auto",
        "unmatched": "error",
        "within_variance": "pooled",
    }
    for name, default in defaults.items():
        if values.get(name) is None:
            values[name] = default
    values["trait"] = values["data"]
    return SimpleNamespace(**values)


def _predictor_replicate_args(args):
    values = vars(args).copy()
    values.update(
        {
            "batch": args.predictor_batch,
            "biological_id": args.predictor_biological_id,
            "sample_size_columns": args.predictor_sample_size_columns,
            "sampling_covariance_out": getattr(
                args, "predictor_sampling_covariance_out", None
            ),
            "standard_error_columns": args.predictor_standard_error_columns,
            "technical_aggregation": args.predictor_technical_aggregation,
            "technical_id": args.predictor_technical_id,
            "tip_summary_out": getattr(args, "predictor_tip_summary_out", None),
            "within_variance": args.predictor_within_variance,
        }
    )
    return SimpleNamespace(**values)


def _ordinary_predictor_values(args, tree, predictors, duplicate_leaf_names):
    dataframe, _, _ = read_tip_table(
        args.data,
        option_name="--data",
        tree_leaf_names=list(tree.leaf_names()),
        required_columns=predictors,
        unmatched=args.unmatched,
        missing_values=args.missing_values,
        duplicate_leaf_names=duplicate_leaf_names,
    )
    dataframe = dataframe[dataframe["leaf_name"].isin(set(tree.leaf_names()))]
    values_by_trait = {}
    for predictor in predictors:
        values_by_leaf = {}
        for leaf_name in tree.leaf_names():
            selected = dataframe.loc[
                dataframe["leaf_name"] == str(leaf_name), predictor
            ]
            if selected.empty or selected.isna().any():
                raise ValueError(
                    "Predictor '{}' has missing values for tree tip '{}'.".format(
                        predictor, leaf_name
                    )
                )
            try:
                values = pd.to_numeric(selected, errors="raise").to_numpy(float)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    "Predictor '{}' must be numeric for tree tip '{}'.".format(
                        predictor, leaf_name
                    )
                ) from exc
            if not np.isfinite(values).all():
                raise ValueError(
                    "Predictor '{}' must be finite for tree tip '{}'.".format(
                        predictor, leaf_name
                    )
                )
            if np.any(values != values[0]):
                raise ValueError(
                    "Species-level predictor '{}' differs among rows for tree tip '{}'.".format(
                        predictor, leaf_name
                    )
                )
            values_by_leaf[str(leaf_name)] = float(values[0])
        values_by_trait[predictor] = values_by_leaf
    return values_by_trait


def _sampling_covariance_table(covariance_by_trait, leaf_names):
    rows = []
    for trait, covariance in covariance_by_trait.items():
        matrix = covariance.loc[leaf_names, leaf_names].to_numpy(dtype=float)
        for first, first_name in enumerate(leaf_names):
            for second in range(first, len(leaf_names)):
                rows.append(
                    {
                        "tree_id": "species",
                        "trait": trait,
                        "leaf_name_1": first_name,
                        "leaf_name_2": leaf_names[second],
                        "sampling_covariance": matrix[first, second],
                    }
                )
    return pd.DataFrame(rows, columns=ORDINARY_SAMPLING_COVARIANCE_COLUMNS)


def build_ordinary_pgls(args, responses, predictors) -> OrdinaryPglsArtifacts:
    effective = _effective_ordinary_args(args)
    tree = read_tree(
        effective.tree,
        effective.tree_format,
        effective.quoted_node_names,
    )
    evolution_spec = evolution_model_spec(effective.evolution_model)
    _validate_ordinary_tree(
        tree,
        effective.branch_length if evolution_spec.branch_lengths_used else "unit",
    )
    response_replicate_requested = _validate_replicate_options(effective)
    predictor_args = _predictor_replicate_args(effective)
    predictor_replicate_requested = _validate_replicate_options(predictor_args)
    if not predictor_replicate_requested:
        irrelevant_predictor_model_options = [
            option
            for name, option in [
                ("predictor_evolution_model", "--predictor-evolution-model"),
                ("predictor_evolution_parameter", "--predictor-evolution-parameter"),
                ("predictor_branch_length", "--predictor-branch-length"),
            ]
            if getattr(effective, name, None) is not None
        ]
        if irrelevant_predictor_model_options:
            raise ValueError(
                "Predictor evolution option(s) require predictor replicates or "
                "known standard errors: {}.".format(
                    ", ".join(irrelevant_predictor_model_options)
                )
            )
    if response_replicate_requested:
        replicate_estimates = _read_replicate_traits(
            effective,
            tree,
            responses,
            "species",
            option_name="--data",
        )
        response_values = replicate_estimates.values_by_trait
        covariance_by_trait = replicate_estimates.sampling_covariance_by_trait
        tip_summary = replicate_estimates.tip_summary
    else:
        replicate_estimates = None
        response_values = (
            _ordinary_predictor_values(
                effective,
                tree,
                responses,
                duplicate_leaf_names="allow",
            )
            if predictor_replicate_requested
            else _read_numeric_traits(
                effective,
                tree,
                responses,
                option_name="--data",
            )
        )
        covariance_by_trait = {}
        tip_summary = pd.DataFrame(columns=TIP_SUMMARY_COLUMNS)
    if predictor_replicate_requested:
        predictor_estimates = _read_replicate_traits(
            predictor_args,
            tree,
            predictors,
            "species",
            option_name="--data",
        )
        predictor_values = predictor_estimates.values_by_trait
        predictor_covariance_by_trait = predictor_estimates.sampling_covariance_by_trait
        predictor_tip_summary = predictor_estimates.tip_summary
    else:
        predictor_values = _ordinary_predictor_values(
            effective,
            tree,
            predictors,
            duplicate_leaf_names=("allow" if response_replicate_requested else "error"),
        )
        predictor_covariance_by_trait = {}
        predictor_tip_summary = pd.DataFrame(columns=TIP_SUMMARY_COLUMNS)
    leaf_names = [str(leaf.name) for leaf in tree.leaves()]
    comparison_models = _parse_comparison_models(
        getattr(effective, "compare_evolution_models", None)
    )
    covariance_path = getattr(effective, "evolution_covariance", None)
    custom_covariance = (
        None
        if covariance_path is None
        else read_custom_covariance(covariance_path, leaf_names)
    )
    results = fit_ordinary_pgls(
        tree,
        response_values,
        predictor_values,
        responses,
        predictors,
        response_sampling_covariance=covariance_by_trait,
        predictor_sampling_covariance=predictor_covariance_by_trait,
        evolution_model=effective.evolution_model,
        evolution_parameter=getattr(effective, "evolution_parameter", None),
        branch_length=effective.branch_length,
        custom_covariance=(
            custom_covariance if effective.evolution_model == "custom" else None
        ),
        intercept=effective.intercept,
        confidence_level=effective.confidence_level,
        inference=effective.inference,
        bootstrap_replicates=effective.bootstrap_replicates,
        seed=effective.seed,
        reml=effective.reml,
        predictor_evolution_model=getattr(effective, "predictor_evolution_model", None),
        predictor_evolution_parameter=getattr(
            effective, "predictor_evolution_parameter", None
        ),
        predictor_branch_length=getattr(effective, "predictor_branch_length", None),
    )
    comparison = (
        fit_ordinary_model_comparison(
            tree,
            response_values,
            predictor_values,
            responses,
            predictors,
            comparison_models,
            response_sampling_covariance=covariance_by_trait,
            predictor_sampling_covariance=predictor_covariance_by_trait,
            branch_length=effective.branch_length,
            custom_covariance=custom_covariance,
            intercept=effective.intercept,
            predictor_evolution_model=(
                getattr(effective, "predictor_evolution_model", None)
                or effective.evolution_model
            ),
            predictor_evolution_parameter=getattr(
                effective, "predictor_evolution_parameter", None
            ),
            predictor_branch_length=getattr(effective, "predictor_branch_length", None),
        )
        if comparison_models
        else pd.DataFrame(columns=ORDINARY_MODEL_COMPARISON_COLUMNS)
    )
    return OrdinaryPglsArtifacts(
        results=results,
        response_sampling_covariance=_sampling_covariance_table(
            covariance_by_trait,
            leaf_names,
        ),
        response_tip_summary=tip_summary,
        predictor_sampling_covariance=_sampling_covariance_table(
            predictor_covariance_by_trait,
            leaf_names,
        ),
        predictor_tip_summary=predictor_tip_summary,
        model_comparison=comparison,
    )


def _paths_identify_same_file(input_path, output_path):
    if normalized_missing_path_key(input_path) == normalized_missing_path_key(
        output_path
    ):
        return True
    if not os.path.exists(input_path) or not os.path.exists(output_path):
        return False
    try:
        return os.path.samefile(input_path, output_path)
    except OSError:
        return False


def _validate_outputs_do_not_replace_inputs(input_paths, output_paths):
    for input_path in input_paths:
        if input_path in (None, "", "-"):
            continue
        for output_path in output_paths:
            if output_path in (None, "", "-"):
                continue
            if _paths_identify_same_file(input_path, output_path):
                raise ValueError(
                    "Conventional PGLS output must not overwrite an input file: '{}'.".format(
                        os.path.realpath(output_path)
                    )
                )


def validate_ordinary_pgls_output_paths(args) -> None:
    """Validate conventional PGLS outputs before fitting or writing."""
    sampling_path = getattr(args, "sampling_covariance_out", None)
    tip_summary_path = getattr(args, "tip_summary_out", None)
    predictor_sampling_path = getattr(args, "predictor_sampling_covariance_out", None)
    predictor_tip_summary_path = getattr(args, "predictor_tip_summary_out", None)
    comparison_path = getattr(args, "model_comparison_out", None)
    for option, path in [
        ("--sampling-covariance-out", sampling_path),
        ("--tip-summary-out", tip_summary_path),
        ("--predictor-sampling-covariance-out", predictor_sampling_path),
        ("--predictor-tip-summary-out", predictor_tip_summary_path),
        ("--model-comparison-out", comparison_path),
    ]:
        if path == "-":
            raise ValueError("'{}' cannot be STDOUT.".format(option))
    validate_distinct_output_paths(
        [
            ("--outfile", args.outfile),
            ("--sampling-covariance-out", sampling_path),
            ("--tip-summary-out", tip_summary_path),
            ("--predictor-sampling-covariance-out", predictor_sampling_path),
            ("--predictor-tip-summary-out", predictor_tip_summary_path),
            ("--model-comparison-out", comparison_path),
        ]
    )
    output_paths = [
        args.outfile,
        sampling_path,
        tip_summary_path,
        predictor_sampling_path,
        predictor_tip_summary_path,
        comparison_path,
    ]
    _validate_outputs_do_not_replace_inputs(
        [
            getattr(args, "tree", None),
            getattr(args, "data", None),
            getattr(args, "evolution_covariance", None),
        ],
        output_paths,
    )


def write_ordinary_pgls_outputs(args, artifacts: OrdinaryPglsArtifacts) -> None:
    validate_ordinary_pgls_output_paths(args)
    sampling_path = getattr(args, "sampling_covariance_out", None)
    tip_summary_path = getattr(args, "tip_summary_out", None)
    predictor_sampling_path = getattr(args, "predictor_sampling_covariance_out", None)
    predictor_tip_summary_path = getattr(args, "predictor_tip_summary_out", None)
    comparison_path = getattr(args, "model_comparison_out", None)
    file_outputs = []
    if args.outfile != "-":
        file_outputs.append((args.outfile, artifacts.results))
    if sampling_path is not None:
        file_outputs.append((sampling_path, artifacts.response_sampling_covariance))
    if tip_summary_path is not None:
        file_outputs.append((tip_summary_path, artifacts.response_tip_summary))
    if predictor_sampling_path is not None:
        file_outputs.append(
            (predictor_sampling_path, artifacts.predictor_sampling_covariance)
        )
    if predictor_tip_summary_path is not None:
        file_outputs.append(
            (predictor_tip_summary_path, artifacts.predictor_tip_summary)
        )
    if comparison_path is not None:
        file_outputs.append((comparison_path, artifacts.model_comparison))
    if file_outputs:
        from nwkit.pgls_pipeline import _write_dataframes_transactionally

        _write_dataframes_transactionally(file_outputs)
    if args.outfile == "-":
        print(artifacts.results.to_csv(sep="\t", index=False), end="")
