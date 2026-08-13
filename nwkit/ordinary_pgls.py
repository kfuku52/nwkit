"""Conventional tip-level phylogenetic generalized least squares."""

import math
import os
import warnings
from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import minimize_scalar
from scipy.stats import chi2, norm
from scipy.stats import t as student_t

from nwkit.contrast import (
    _read_mixed_replicate_traits,
    _validate_replicate_options,
)
from nwkit.conventions import DEFAULT_TABLE_MISSING_VALUES_CSV
from nwkit.evolution import (
    EVOLUTION_MODELS,
    build_evolutionary_covariance,
    encoded_evolution_parameter,
    evolution_model_spec,
    evolutionary_covariance_factory,
    optimization_parameterization,
    parameter_near_boundary,
    read_custom_covariance,
    validate_evolution_parameter,
)
from nwkit.gaussian import DiagonalLowRankCovariance, draw_from_factor
from nwkit.measurement_error import (
    fit_conditional_eiv_gaussian,
    fit_latent_predictor,
)
from nwkit.model_matrix import (
    PredictorTerm,
    ReplicatedObservation,
    encode_predictors,
    parse_key_values,
    parse_name_list,
    parse_numeric_key_values,
    parse_ordered_levels,
    resolve_response_specs,
    response_family_link,
    validate_response_auxiliaries,
)
from nwkit.multivariate_pgls import fit_multivariate_pgls
from nwkit.pgls import _profile_covariance_fit, _solve_positive_definite
from nwkit.phylogenetic_glmm import (
    SCALAR_RESPONSE_FAMILIES,
    fit_phylogenetic_glmm,
    summarize_glmm_coefficient,
    summarize_glmm_omnibus,
    summarize_glmm_threshold,
)
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
    "response_type",
    "response_family",
    "response_level",
    "response_reference",
    "link_function",
    "response_dispersion",
    "zero_probability",
    "coefficient_penalty",
    "coefficient_prior_sd",
    "separation_warning",
    "term",
    "source_term",
    "predictor_type",
    "predictor_level",
    "predictor_reference",
    "factor_coding",
    "term_test",
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
    allow_large_dense=False,
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
                reml=False,
                allow_large_dense=allow_large_dense,
            )
        else:
            fit = _profile_covariance_fit(
                y,
                design,
                fixed_covariance,
                components,
                reml=reml,
                allow_large_dense=allow_large_dense,
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
    allow_large_dense=False,
):
    rng = np.random.default_rng(seed)
    coefficients: list[np.ndarray] = []
    mean = design @ fit["beta"]
    maximum_attempts = max(replicates * 3, replicates + 10)
    attempts = 0
    while len(coefficients) < replicates and attempts < maximum_attempts:
        attempts += 1
        simulated = mean + draw_from_factor(
            fit["cholesky"], rng.standard_normal(len(mean)), rng=rng
        )
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
                allow_large_dense=allow_large_dense,
            )
        except ValueError:
            continue
        if not bootstrap_fit["optimizer_converged"]:
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
    if inference not in {
        "wald",
        "parametric-bootstrap",
        "likelihood-ratio",
        "profile-likelihood",
    }:
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
    encoded_predictors,
    leaf_names,
    *,
    intercept,
):
    predictor_columns = [
        _ordered_values(
            encoded_predictors.values_by_trait[term.name], leaf_names, term.name
        )
        for term in encoded_predictors.terms
    ]
    design = np.column_stack(predictor_columns)
    term_names = encoded_predictors.term_names
    term_metadata = list(encoded_predictors.terms)
    if intercept:
        design = np.column_stack([np.ones(len(leaf_names)), design])
        term_names = ["(intercept)"] + term_names
        term_metadata = [
            PredictorTerm("(intercept)", "(intercept)", "intercept")
        ] + term_metadata
    num_parameters = design.shape[1]
    if len(leaf_names) <= num_parameters:
        raise ValueError(
            "Ordinary PGLS needs more species than fitted coefficients "
            "(species={}; coefficients={}).".format(len(leaf_names), num_parameters)
        )
    matrix_rank = int(np.linalg.matrix_rank(design))
    if matrix_rank != num_parameters:
        raise ValueError("Ordinary PGLS design matrix is rank deficient.")
    return design, term_names, term_metadata, num_parameters, matrix_rank


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
    categorical_predictors=(),
    ordered_predictors=None,
    factor_references=None,
    factor_coding="treatment",
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
    encoded_predictors = encode_predictors(
        predictor_values_by_trait,
        predictors,
        leaf_names,
        categorical=categorical_predictors,
        ordered_levels=ordered_predictors,
        factor_references=factor_references,
        factor_coding=factor_coding,
    )
    design, term_names, term_metadata, num_parameters, matrix_rank = (
        _ordinary_design_matrix(
            encoded_predictors,
            leaf_names,
            intercept=intercept,
        )
    )
    return (
        responses,
        predictors,
        leaf_names,
        design,
        term_names,
        term_metadata,
        encoded_predictors,
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
    if covariance.shape not in {
        (len(leaf_names),),
        (len(leaf_names), len(leaf_names)),
    }:
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
    allow_large_dense=False,
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
        allow_large_dense=allow_large_dense,
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
    fixed_array = np.asarray(fixed_covariance, dtype=float)
    mean_sampling_variance = float(
        np.mean(fixed_array if fixed_array.ndim == 1 else np.diag(fixed_array))
    )
    fitted_covariance = fit["covariance"]
    if isinstance(fitted_covariance, DiagonalLowRankCovariance):
        update = fitted_covariance.low_rank
        update_diagonal = (
            np.asarray(update.multiply(update).sum(axis=1)).reshape(-1)
            if sparse.issparse(update)
            else np.sum(np.square(update), axis=1)
        )
        mean_fitted_variance = float(
            np.mean(fitted_covariance.diagonal + update_diagonal)
        )
    else:
        fitted_array = np.asarray(fitted_covariance, dtype=float)
        mean_fitted_variance = float(
            np.mean(fitted_array if fitted_array.ndim == 1 else np.diag(fitted_array))
        )
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


def _ordinary_coefficient_interval(
    coefficient,
    standard_error,
    critical,
    bootstrap_coefficients,
    term_index,
    confidence_level,
):
    if bootstrap_coefficients is None:
        return (
            coefficient - critical * standard_error,
            coefficient + critical * standard_error,
        )
    alpha = (1.0 - confidence_level) / 2.0
    lower, upper = np.quantile(
        bootstrap_coefficients[:, term_index], [alpha, 1.0 - alpha]
    )
    return float(lower), float(upper)


def _ordinary_coefficient_test(
    coefficient,
    standard_error,
    bootstrap_coefficients,
    term_index,
    degrees_of_freedom,
):
    if standard_error == 0.0:
        return "", "", "zero-model-variance"
    statistic = coefficient / standard_error
    if bootstrap_coefficients is None:
        p_value = float(2.0 * student_t.sf(abs(statistic), degrees_of_freedom))
    else:
        centered = bootstrap_coefficients[:, term_index] - coefficient
        p_value = float(
            (1 + np.sum(np.abs(centered) >= abs(coefficient))) / (len(centered) + 1)
        )
    return statistic, p_value, "ok"


def _ordinary_omnibus_rows(rows, fit, term_metadata):
    metadata_by_source: dict[str, list[int]] = {}
    for index, metadata in enumerate(term_metadata):
        if metadata.kind in {"categorical", "ordered"}:
            metadata_by_source.setdefault(metadata.source, []).append(index)
    omnibus_rows = []
    for source, indices in metadata_by_source.items():
        coefficients = np.asarray(fit["beta"], dtype=float)[indices]
        covariance = np.asarray(fit["beta_covariance"], dtype=float)[
            np.ix_(indices, indices)
        ]
        try:
            statistic = float(coefficients @ np.linalg.solve(covariance, coefficients))
        except np.linalg.LinAlgError:
            statistic = float(coefficients @ np.linalg.pinv(covariance) @ coefficients)
        omnibus = rows[indices[0]].copy()
        omnibus.update(
            {
                "term": source,
                "source_term": source,
                "predictor_level": "",
                "term_test": "omnibus",
                "coefficient": "",
                "standard_error": "",
                "statistic": statistic,
                "degrees_of_freedom": len(indices),
                "p_value": float(chi2.sf(statistic, len(indices))),
                "confidence_interval_lower": "",
                "confidence_interval_upper": "",
                "inference_status": "ok",
                "inference_method": "wald",
            }
        )
        omnibus_rows.append(omnibus)
    return omnibus_rows


def _ordinary_coefficient_rows(
    response,
    fit,
    standard_errors,
    bootstrap_coefficients,
    term_names,
    term_metadata,
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
    has_predictor_uncertainty,
):
    critical = float(student_t.ppf(0.5 + confidence_level / 2.0, degrees_of_freedom))
    rows = []
    for term_index, (term, metadata) in enumerate(
        zip(term_names, term_metadata, strict=True)
    ):
        coefficient = float(fit["beta"][term_index])
        standard_error = float(standard_errors[term_index])
        lower, upper = _ordinary_coefficient_interval(
            coefficient,
            standard_error,
            critical,
            bootstrap_coefficients,
            term_index,
            confidence_level,
        )
        statistic, p_value, inference_status = _ordinary_coefficient_test(
            coefficient,
            standard_error,
            bootstrap_coefficients,
            term_index,
            degrees_of_freedom,
        )
        parameter = fit["evolution_parameter"]
        spec = evolution_model_spec(evolution_model)
        predictor_diagnostic = predictor_diagnostics.get(metadata.source)
        rows.append(
            {
                "model_id": "ordinary:{}".format(response),
                "tree_id": "species",
                "response": response,
                "response_type": "continuous",
                "response_family": "gaussian",
                "response_level": "",
                "response_reference": "",
                "link_function": "identity",
                "term": term,
                "source_term": metadata.source,
                "predictor_type": metadata.kind,
                "predictor_level": metadata.level,
                "predictor_reference": metadata.reference,
                "factor_coding": metadata.coding,
                "term_test": "coefficient",
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
                    mean_predictor_sampling_variance.get(metadata.source, "")
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
                    if has_predictor_uncertainty
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
    return rows + _ordinary_omnibus_rows(rows, fit, term_metadata)


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
        return predictor_values_by_trait, {}, {}
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
    predictor_uncertainties = {}
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
        predictor_uncertainties[predictor] = posterior.covariance
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
            "mean_sampling_variance": float(
                np.mean(sampling if sampling.ndim == 1 else np.diag(sampling))
            ),
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


def _categorical_shape_settings(tree, model, parameter, branch_length):
    spec = evolution_model_spec(model)
    if spec.parameter_name is None or parameter is not None:
        return None, float, None
    bounds, decoder = optimization_parameterization(
        tree, model, branch_length=branch_length
    )
    encoded_defaults = {
        "lambda": 0.5,
        "ou": (bounds[0] + bounds[1]) / 2.0,
        "kappa": 1.0,
        "delta": encoded_evolution_parameter("delta", 1.0),
        "eb": bounds[1] - 0.1 * (bounds[1] - bounds[0]),
        "acdc": 0.0,
    }
    return bounds, decoder, encoded_defaults[model]


def _categorical_common_row(
    response,
    response_spec,
    fit,
    *,
    n_species,
    n_predictors,
    num_parameters,
    matrix_rank,
    intercept,
    evolution_model,
    branch_length,
    inference_status,
    measurement_error_model,
):
    spec = evolution_model_spec(evolution_model)
    parameter_status = (
        "not-applicable"
        if spec.parameter_name is None
        else fit.evolution_parameter_status
    )
    return {
        "model_id": "ordinary:{}".format(response),
        "tree_id": "species",
        "response": response,
        "response_type": response_spec.kind,
        "response_family": response_spec.family,
        "response_level": "",
        "response_reference": fit.reference,
        "link_function": response_family_link(response_spec.family),
        "response_dispersion": "" if fit.dispersion is None else fit.dispersion,
        "zero_probability": (
            "" if fit.zero_probability is None else fit.zero_probability
        ),
        "coefficient_penalty": fit.coefficient_penalty,
        "coefficient_prior_sd": (
            "" if fit.coefficient_prior_sd is None else fit.coefficient_prior_sd
        ),
        "separation_warning": "yes" if fit.separation_warning else "no",
        "degrees_of_freedom": 1,
        "n_species": n_species,
        "n_predictors": n_predictors,
        "num_parameters": num_parameters,
        "matrix_rank": matrix_rank,
        "condition_number": (
            float(np.linalg.cond(fit.coefficient_covariance))
            if np.isfinite(fit.coefficient_covariance).all()
            else ""
        ),
        "generalized_residual_sum_squares": "",
        "evolutionary_rate": fit.component_variances["phylogenetic"],
        "r_squared": "",
        "intercept": "yes" if intercept else "no",
        "evolution_model": evolution_model,
        "evolution_parameter_name": spec.parameter_name or "",
        "evolution_parameter": (
            "" if fit.evolution_parameter is None else fit.evolution_parameter
        ),
        "evolution_parameter_status": parameter_status,
        "branch_length_mode": (
            branch_length if spec.branch_lengths_used else "not-applicable"
        ),
        "covariance_estimator": "laplace-ML",
        "inference_method": fit.coefficient_inference,
        "reml": "no",
        "mean_sampling_variance": "",
        "sampling_variance_fraction": "",
        "measurement_error_model": measurement_error_model,
        "log_likelihood": fit.log_likelihood,
        "optimizer_converged": "yes" if fit.optimizer_converged else "no",
        "optimizer_message": fit.optimizer_message,
        "boundary_warning": "yes" if fit.boundary_warning else "no",
        "small_sample_warning": "yes" if n_species < 20 else "no",
        "inference_status": inference_status,
        "model": "ordinary-pglmm",
    }


def _categorical_coefficient_rows(
    response,
    response_spec,
    fit,
    term_names,
    term_metadata,
    *,
    confidence_level,
    n_species,
    n_predictors,
    matrix_rank,
    intercept,
    evolution_model,
    branch_length,
    predictor_uncertainties,
    mean_predictor_sampling_variance,
):
    coefficient_count = fit.coefficients.size
    num_parameters = (
        coefficient_count
        + len(fit.thresholds)
        + len(fit.component_variances)
        + (1 if fit.evolution_parameter_status == "estimated" else 0)
        + (1 if fit.dispersion is not None else 0)
        + (1 if fit.zero_probability is not None else 0)
    )
    common = _categorical_common_row(
        response,
        response_spec,
        fit,
        n_species=n_species,
        n_predictors=n_predictors,
        num_parameters=num_parameters,
        matrix_rank=matrix_rank,
        intercept=intercept,
        evolution_model=evolution_model,
        branch_length=branch_length,
        inference_status="ok",
        measurement_error_model=(
            "latent-predictor" if predictor_uncertainties else "none"
        ),
    )
    critical = float(norm.ppf(0.5 + confidence_level / 2.0))
    response_levels = (
        [""]
        if response_spec.family == "ordinal"
        or response_spec.family in SCALAR_RESPONSE_FAMILIES
        else list(fit.levels[:-1])
    )
    rows = []
    coefficient_indices_by_source: dict[str, list[int]] = {}
    dimension = fit.coefficients.shape[1]
    for term_index, (term, metadata) in enumerate(
        zip(term_names, term_metadata, strict=True)
    ):
        for level_index, level in enumerate(response_levels):
            flat_index = term_index * dimension + level_index
            coefficient = float(fit.coefficients[term_index, level_index])
            (
                standard_error,
                statistic,
                p_value,
                lower,
                upper,
                inference_status,
            ) = summarize_glmm_coefficient(fit, flat_index, coefficient, critical)
            row = common.copy()
            row.update(
                {
                    "response_level": level,
                    "term": term,
                    "source_term": metadata.source,
                    "predictor_type": metadata.kind,
                    "predictor_level": metadata.level,
                    "predictor_reference": metadata.reference,
                    "factor_coding": metadata.coding,
                    "term_test": "coefficient",
                    "coefficient": coefficient,
                    "standard_error": standard_error,
                    "statistic": statistic,
                    "p_value": p_value,
                    "confidence_level": confidence_level,
                    "confidence_interval_lower": lower,
                    "confidence_interval_upper": upper,
                    "inference_status": inference_status,
                }
            )
            if metadata.kind != "intercept":
                row["mean_predictor_sampling_variance"] = (
                    mean_predictor_sampling_variance.get(metadata.source, "")
                )
            rows.append(row)
            if metadata.kind != "intercept" and (
                dimension > 1 or metadata.kind in {"categorical", "ordered"}
            ):
                coefficient_indices_by_source.setdefault(metadata.source, []).append(
                    flat_index
                )
    for source, indices in coefficient_indices_by_source.items():
        omnibus_statistic, omnibus_p_value, omnibus_status = summarize_glmm_omnibus(
            fit, indices
        )
        template_index = (
            next(
                index
                for index, metadata in enumerate(term_metadata)
                if metadata.source == source
            )
            * dimension
        )
        template_row = rows[template_index].copy()
        template_row.update(
            {
                "response_level": "",
                "term": source,
                "source_term": source,
                "predictor_level": "",
                "term_test": "omnibus",
                "coefficient": "",
                "standard_error": "",
                "statistic": omnibus_statistic,
                "degrees_of_freedom": len(indices),
                "p_value": omnibus_p_value,
                "confidence_interval_lower": "",
                "confidence_interval_upper": "",
                "inference_status": omnibus_status,
                "inference_method": "wald",
            }
        )
        rows.append(template_row)
    for threshold_index, threshold in enumerate(fit.thresholds):
        (
            threshold_standard_error,
            lower,
            upper,
            threshold_status,
        ) = summarize_glmm_threshold(fit, threshold_index, critical)
        row = common.copy()
        row.update(
            {
                "response_level": fit.levels[threshold_index],
                "term": "threshold[{}]".format(fit.levels[threshold_index]),
                "source_term": "(threshold)",
                "predictor_type": "threshold",
                "term_test": "threshold",
                "coefficient": float(threshold),
                "standard_error": threshold_standard_error,
                "statistic": "",
                "p_value": "",
                "confidence_level": confidence_level,
                "confidence_interval_lower": lower,
                "confidence_interval_upper": upper,
                "inference_status": threshold_status,
            }
        )
        rows.append(row)
    return rows


def _ordinary_multivariate_rows(
    responses,
    fit,
    term_names,
    term_metadata,
    *,
    confidence_level,
    n_species,
    n_predictors,
    intercept,
    evolution_model,
    branch_length,
):
    critical = float(norm.ppf(0.5 + confidence_level / 2.0))
    spec = evolution_model_spec(evolution_model)
    evolutionary_covariance = fit.component_trait_covariances["phylogenetic"]
    num_parameters = (
        fit.coefficients.size
        + sum(
            covariance.shape[0] * (covariance.shape[0] + 1) // 2
            for covariance in fit.component_trait_covariances.values()
        )
        + (1 if fit.evolution_parameter_status == "estimated" else 0)
    )
    common = {
        "tree_id": "species",
        "response_type": "multivariate-continuous",
        "response_family": "multivariate-gaussian",
        "response_level": "",
        "response_reference": "",
        "link_function": "identity",
        "response_dispersion": "",
        "zero_probability": "",
        "coefficient_penalty": "none",
        "coefficient_prior_sd": "",
        "separation_warning": "no",
        "degrees_of_freedom": fit.n_observations - fit.coefficients.size,
        "n_species": n_species,
        "n_predictors": n_predictors,
        "num_parameters": num_parameters,
        "matrix_rank": int(np.linalg.matrix_rank(fit.coefficient_covariance)),
        "condition_number": float(np.linalg.cond(fit.coefficient_covariance)),
        "generalized_residual_sum_squares": "",
        "r_squared": "",
        "intercept": "yes" if intercept else "no",
        "evolution_model": evolution_model,
        "evolution_parameter_name": spec.parameter_name or "",
        "evolution_parameter": (
            "" if fit.evolution_parameter is None else fit.evolution_parameter
        ),
        "evolution_parameter_status": (
            "not-applicable"
            if spec.parameter_name is None
            else fit.evolution_parameter_status
        ),
        "branch_length_mode": (
            branch_length if spec.branch_lengths_used else "not-applicable"
        ),
        "covariance_estimator": (
            "multivariate-gaussian-REML" if fit.reml else "multivariate-gaussian-ML"
        ),
        "inference_method": "wald",
        "reml": "yes" if fit.reml else "no",
        "measurement_error_model": "response-covariance",
        "log_likelihood": fit.log_likelihood,
        "optimizer_converged": "yes" if fit.optimizer_converged else "no",
        "optimizer_message": fit.optimizer_message,
        "boundary_warning": "yes" if fit.boundary_warning else "no",
        "small_sample_warning": "yes" if n_species < 20 else "no",
        "inference_status": "ok",
        "model": "ordinary-multivariate-pgls",
    }
    rows = []
    coefficient_count = len(term_names)
    for response_index, response in enumerate(responses):
        for term_index, (term, metadata) in enumerate(
            zip(term_names, term_metadata, strict=True)
        ):
            flat_index = response_index * coefficient_count + term_index
            coefficient = float(fit.coefficients[response_index, term_index])
            standard_error = math.sqrt(
                max(float(fit.coefficient_covariance[flat_index, flat_index]), 0.0)
            )
            statistic = "" if standard_error == 0.0 else coefficient / standard_error
            row = common.copy()
            row.update(
                {
                    "model_id": "ordinary:multivariate",
                    "response": response,
                    "term": term,
                    "source_term": metadata.source,
                    "predictor_type": metadata.kind,
                    "predictor_level": metadata.level,
                    "predictor_reference": metadata.reference,
                    "factor_coding": metadata.coding,
                    "term_test": "coefficient",
                    "coefficient": coefficient,
                    "standard_error": standard_error,
                    "statistic": statistic,
                    "p_value": (
                        ""
                        if standard_error == 0.0
                        else float(2.0 * norm.sf(abs(float(statistic))))
                    ),
                    "confidence_level": confidence_level,
                    "confidence_interval_lower": coefficient
                    - critical * standard_error,
                    "confidence_interval_upper": coefficient
                    + critical * standard_error,
                    "evolutionary_rate": evolutionary_covariance[
                        response_index, response_index
                    ],
                }
            )
            rows.append(row)
    for first, first_response in enumerate(responses):
        for second in range(first, len(responses)):
            row = common.copy()
            second_response = responses[second]
            row.update(
                {
                    "model_id": "ordinary:multivariate",
                    "response": "{}|{}".format(first_response, second_response),
                    "term": "covariance[{},{}]".format(first_response, second_response),
                    "source_term": "(response-covariance)",
                    "predictor_type": "response-covariance",
                    "term_test": "response-covariance",
                    "coefficient": evolutionary_covariance[first, second],
                    "confidence_level": confidence_level,
                    "evolutionary_rate": evolutionary_covariance[first, second],
                }
            )
            rows.append(row)
    return rows


def _ordinary_auxiliary_sequence(mapping, trait, leaf_names):
    if trait not in mapping:
        return None
    return [mapping[trait][name] for name in leaf_names]


def _fit_ordinary_multivariate_response(
    tree,
    response_values_by_trait,
    responses,
    predictors,
    leaf_names,
    design,
    term_names,
    term_metadata,
    response_specs,
    covariance_by_trait,
    predictor_uncertainty_values,
    *,
    allow_missing_responses,
    inference,
    evolution_model,
    evolution_parameter,
    branch_length,
    custom_covariance,
    reml,
    confidence_level,
    intercept,
    allow_large_dense,
):
    if len(responses) < 2:
        raise ValueError("Multivariate PGLS requires at least two responses.")
    if any(spec.family != "gaussian" for spec in response_specs.values()):
        raise ValueError("Multivariate PGLS currently requires Gaussian responses.")
    if inference != "wald":
        raise ValueError("Multivariate PGLS currently supports Wald inference.")
    if predictor_uncertainty_values:
        raise ValueError(
            "Multivariate PGLS does not yet combine predictor measurement error; "
            "fit separate response models for that combination."
        )
    response_matrix = np.asarray(
        [
            [response_values_by_trait[response][name] for response in responses]
            for name in leaf_names
        ],
        dtype=float,
    )
    if not allow_missing_responses and not np.isfinite(response_matrix).all():
        raise ValueError(
            "Missing multivariate responses require allow_missing_responses=True."
        )
    fixed_diagonal = np.zeros(len(leaf_names) * len(responses), dtype=float)
    dense_sampling_blocks: list[tuple[int, np.ndarray]] = []
    for response_index, response in enumerate(responses):
        covariance = _coerce_response_sampling_covariance(
            covariance_by_trait, response, leaf_names
        )
        start = response_index * len(leaf_names)
        if covariance.ndim == 1 or np.array_equal(
            covariance, np.diag(np.diag(covariance))
        ):
            fixed_diagonal[start : start + len(leaf_names)] = (
                covariance if covariance.ndim == 1 else np.diag(covariance)
            )
        else:
            dense_sampling_blocks.append((start, covariance))
    if dense_sampling_blocks:
        fixed_covariance = np.diag(fixed_diagonal)
        for start, covariance in dense_sampling_blocks:
            fixed_covariance[
                start : start + len(leaf_names), start : start + len(leaf_names)
            ] = covariance
    else:
        fixed_covariance = fixed_diagonal
    shape_bounds, shape_decoder, shape_initial = _categorical_shape_settings(
        tree, evolution_model, evolution_parameter, branch_length
    )

    covariance_factory = evolutionary_covariance_factory(
        tree,
        leaf_names,
        model=evolution_model,
        branch_length=branch_length,
        custom_covariance=custom_covariance,
    )

    multivariate_fit = fit_multivariate_pgls(
        response_matrix,
        design,
        {"phylogenetic": covariance_factory},
        fixed_covariance=fixed_covariance,
        evolution_parameter=evolution_parameter,
        evolution_parameter_bounds=shape_bounds,
        evolution_parameter_decoder=shape_decoder,
        evolution_parameter_initial=shape_initial,
        reml=reml,
        allow_large_dense=allow_large_dense,
    )
    return pd.DataFrame(
        _ordinary_multivariate_rows(
            responses,
            multivariate_fit,
            term_names,
            term_metadata,
            confidence_level=confidence_level,
            n_species=len(leaf_names),
            n_predictors=len(predictors),
            intercept=intercept,
            evolution_model=evolution_model,
            branch_length=branch_length,
        ),
        columns=ORDINARY_RESULT_COLUMNS,
    )


def _ordinal_glmm_design(
    response_spec,
    design,
    term_names,
    term_metadata,
    predictor_columns,
    intercept,
):
    if response_spec.family != "ordinal" or term_metadata[0].kind != "intercept":
        return design, term_names, term_metadata, predictor_columns, intercept
    adjusted_columns = [
        (
            column - 1
            if isinstance(column, int)
            else tuple(index - 1 for index in column)
        )
        for column in predictor_columns
    ]
    return design[:, 1:], term_names[1:], term_metadata[1:], adjusted_columns, False


def _fit_ordinary_non_gaussian_response(
    tree,
    response,
    response_index,
    response_spec,
    response_values_by_trait,
    predictors,
    leaf_names,
    design,
    term_names,
    term_metadata,
    predictor_uncertainty_values,
    predictor_columns,
    covariance_by_trait,
    response_offsets,
    response_trials,
    response_censor_lower,
    response_censor_upper,
    response_dispersions,
    response_zero_probabilities,
    mean_predictor_sampling_variance,
    *,
    evolution_model,
    evolution_parameter,
    branch_length,
    custom_covariance,
    coefficient_penalty,
    coefficient_prior_sd,
    inference,
    confidence_level,
    bootstrap_replicates,
    seed,
    intercept,
    allow_large_dense,
):
    if response in covariance_by_trait:
        raise ValueError(
            "Gaussian known-SE response covariance does not apply to "
            "non-Gaussian likelihood responses."
        )
    shape_bounds, shape_decoder, shape_initial = _categorical_shape_settings(
        tree, evolution_model, evolution_parameter, branch_length
    )
    (
        glmm_design,
        glmm_term_names,
        glmm_term_metadata,
        glmm_predictor_columns,
        glmm_intercept,
    ) = _ordinal_glmm_design(
        response_spec,
        design,
        term_names,
        term_metadata,
        predictor_columns,
        intercept,
    )

    covariance_factory = evolutionary_covariance_factory(
        tree,
        leaf_names,
        model=evolution_model,
        branch_length=branch_length,
        custom_covariance=custom_covariance,
    )

    fitted = fit_phylogenetic_glmm(
        [response_values_by_trait[response][name] for name in leaf_names],
        glmm_design,
        covariance_factory,
        family=response_spec.family,
        levels=response_spec.levels,
        reference=response_spec.reference,
        evolution_parameter=evolution_parameter,
        evolution_parameter_bounds=shape_bounds,
        evolution_parameter_decoder=shape_decoder,
        evolution_parameter_initial=shape_initial,
        predictor_uncertainties=predictor_uncertainty_values,
        predictor_columns=glmm_predictor_columns,
        offset=_ordinary_auxiliary_sequence(response_offsets, response, leaf_names),
        trials=_ordinary_auxiliary_sequence(response_trials, response, leaf_names),
        censor_lower=_ordinary_auxiliary_sequence(
            response_censor_lower, response, leaf_names
        ),
        censor_upper=_ordinary_auxiliary_sequence(
            response_censor_upper, response, leaf_names
        ),
        dispersion=response_dispersions.get(response),
        zero_probability=response_zero_probabilities.get(response),
        coefficient_penalty=coefficient_penalty,
        coefficient_prior_sd=coefficient_prior_sd,
        inference=inference,
        confidence_level=confidence_level,
        bootstrap_replicates=bootstrap_replicates,
        seed=seed + response_index,
        allow_large_dense=allow_large_dense,
    )
    return _categorical_coefficient_rows(
        response,
        response_spec,
        fitted,
        glmm_term_names,
        glmm_term_metadata,
        confidence_level=confidence_level,
        n_species=len(leaf_names),
        n_predictors=len(predictors),
        matrix_rank=int(np.linalg.matrix_rank(glmm_design)),
        intercept=glmm_intercept,
        evolution_model=evolution_model,
        branch_length=branch_length,
        predictor_uncertainties=predictor_uncertainty_values,
        mean_predictor_sampling_variance=mean_predictor_sampling_variance,
    )


def _fit_ordinary_gaussian_response(
    tree,
    response,
    response_index,
    response_values_by_trait,
    predictors,
    leaf_names,
    design,
    term_names,
    term_metadata,
    covariance_by_trait,
    predictor_uncertainty_values,
    predictor_columns,
    predictor_diagnostics,
    mean_predictor_sampling_variance,
    *,
    evolution_model,
    evolution_parameter,
    branch_length,
    custom_covariance,
    inference,
    bootstrap_replicates,
    seed,
    reml,
    intercept,
    confidence_level,
    num_parameters,
    matrix_rank,
    allow_large_dense,
):
    if inference in {"likelihood-ratio", "profile-likelihood"}:
        raise ValueError(
            "Gaussian PGLS supports Wald or parametric-bootstrap inference."
        )
    y = _ordered_values(response_values_by_trait[response], leaf_names, response)
    fixed_covariance = _coerce_response_sampling_covariance(
        covariance_by_trait, response, leaf_names
    )
    fitted = _fit_ordinary_gaussian(
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
        predictor_uncertainties=predictor_uncertainty_values,
        predictor_columns=predictor_columns,
        allow_large_dense=allow_large_dense,
    )
    effective_reml = bool(fitted.get("reml", reml))
    bootstrap_coefficients, standard_errors = _ordinary_inference_samples(
        fitted,
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
        reml=effective_reml,
        predictor_uncertainties=predictor_uncertainty_values,
        predictor_columns=predictor_columns,
        allow_large_dense=allow_large_dense,
    )
    statistics = _ordinary_response_statistics(
        y, fitted, fixed_covariance, intercept=intercept
    )
    return _ordinary_coefficient_rows(
        response,
        fitted,
        standard_errors,
        bootstrap_coefficients,
        term_names,
        term_metadata,
        statistics,
        confidence_level=confidence_level,
        degrees_of_freedom=len(leaf_names) - num_parameters,
        n_species=len(leaf_names),
        n_predictors=len(predictors),
        num_parameters=num_parameters,
        matrix_rank=matrix_rank,
        intercept=intercept,
        evolution_model=evolution_model,
        branch_length=branch_length,
        inference=inference,
        reml=effective_reml,
        predictor_diagnostics=predictor_diagnostics,
        mean_predictor_sampling_variance=mean_predictor_sampling_variance,
        has_predictor_uncertainty=bool(predictor_uncertainty_values),
    )


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
    categorical_predictors=(),
    ordered_predictors=None,
    factor_references=None,
    factor_coding="treatment",
    categorical_responses=(),
    ordered_responses=None,
    response_references=None,
    response_families=None,
    response_offsets=None,
    response_trials=None,
    response_censor_lower=None,
    response_censor_upper=None,
    response_dispersions=None,
    response_zero_probabilities=None,
    coefficient_penalty="student-t",
    coefficient_prior_sd=2.5,
    multivariate_responses=False,
    allow_missing_responses=False,
    allow_large_dense=False,
):
    """Fit conventional tip-level PGLS models, one per response trait."""
    if allow_missing_responses and not multivariate_responses:
        raise ValueError(
            "allow_missing_responses requires multivariate_responses=True."
        )
    leaf_names = [str(leaf.name) for leaf in tree.leaves()]
    initial_encoding = encode_predictors(
        predictor_values_by_trait,
        predictors,
        leaf_names,
        categorical=categorical_predictors,
        ordered_levels=ordered_predictors,
        factor_references=factor_references,
        factor_coding=factor_coding,
    )
    continuous_predictors = [
        term.source for term in initial_encoding.terms if term.kind == "continuous"
    ]
    predictor_model = predictor_evolution_model or evolution_model
    predictor_length = predictor_branch_length or branch_length
    (
        posterior_predictor_values,
        predictor_uncertainties,
        predictor_diagnostics,
    ) = _prepare_latent_ordinary_predictors(
        tree,
        predictor_values_by_trait,
        continuous_predictors,
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
        term_metadata,
        encoded_predictors,
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
        categorical_predictors=categorical_predictors,
        ordered_predictors=ordered_predictors,
        factor_references=factor_references,
        factor_coding=factor_coding,
    )
    covariance_by_trait = response_sampling_covariance or {}
    offset = 1 if intercept else 0
    term_index = {
        term.name: index + offset for index, term in enumerate(encoded_predictors.terms)
    }
    predictor_uncertainty_values = []
    predictor_columns: list[int | tuple[int, ...]] = []
    for predictor, covariance in predictor_uncertainties.items():
        predictor_uncertainty_values.append(covariance)
        predictor_columns.append(term_index[predictor])
    for uncertainty in encoded_predictors.uncertainties:
        predictor_uncertainty_values.append(uncertainty.covariance_by_observation)
        predictor_columns.append(
            tuple(term_index[term] for term in uncertainty.term_names)
        )
    mean_predictor_sampling_variance = {
        predictor: diagnostics["mean_sampling_variance"]
        for predictor, diagnostics in predictor_diagnostics.items()
    }
    mean_predictor_sampling_variance.update(
        {
            uncertainty.source: float(
                np.mean(
                    np.diagonal(uncertainty.covariance_by_observation, axis1=1, axis2=2)
                )
            )
            for uncertainty in encoded_predictors.uncertainties
        }
    )
    response_specs = resolve_response_specs(
        response_values_by_trait,
        responses,
        leaf_names,
        categorical=categorical_responses,
        ordered_levels=ordered_responses,
        references=response_references,
        families=response_families,
        allow_missing=(
            allow_missing_responses
            or any(
                family == "censored-gaussian"
                for family in (response_families or {}).values()
            )
        ),
    )
    response_offsets = {} if response_offsets is None else response_offsets
    response_trials = {} if response_trials is None else response_trials
    response_censor_lower = (
        {} if response_censor_lower is None else response_censor_lower
    )
    response_censor_upper = (
        {} if response_censor_upper is None else response_censor_upper
    )
    response_dispersions = {} if response_dispersions is None else response_dispersions
    response_zero_probabilities = (
        {} if response_zero_probabilities is None else response_zero_probabilities
    )
    validate_response_auxiliaries(
        response_specs,
        offsets=response_offsets,
        trials=response_trials,
        censor_lower=response_censor_lower,
        censor_upper=response_censor_upper,
        dispersions=response_dispersions,
        zero_probabilities=response_zero_probabilities,
    )

    if multivariate_responses:
        return _fit_ordinary_multivariate_response(
            tree,
            response_values_by_trait,
            responses,
            predictors,
            leaf_names,
            design,
            term_names,
            term_metadata,
            response_specs,
            covariance_by_trait,
            predictor_uncertainty_values,
            allow_missing_responses=allow_missing_responses,
            inference=inference,
            evolution_model=evolution_model,
            evolution_parameter=evolution_parameter,
            branch_length=branch_length,
            custom_covariance=custom_covariance,
            reml=reml,
            confidence_level=confidence_level,
            intercept=intercept,
            allow_large_dense=allow_large_dense,
        )

    rows = []
    for response_index, response in enumerate(responses):
        if response not in response_values_by_trait:
            raise ValueError("Response trait '{}' is absent.".format(response))
        response_spec = response_specs[response]
        if response_spec.family != "gaussian":
            rows.extend(
                _fit_ordinary_non_gaussian_response(
                    tree,
                    response,
                    response_index,
                    response_spec,
                    response_values_by_trait,
                    predictors,
                    leaf_names,
                    design,
                    term_names,
                    term_metadata,
                    predictor_uncertainty_values,
                    predictor_columns,
                    covariance_by_trait,
                    response_offsets,
                    response_trials,
                    response_censor_lower,
                    response_censor_upper,
                    response_dispersions,
                    response_zero_probabilities,
                    mean_predictor_sampling_variance,
                    evolution_model=evolution_model,
                    evolution_parameter=evolution_parameter,
                    branch_length=branch_length,
                    custom_covariance=custom_covariance,
                    coefficient_penalty=coefficient_penalty,
                    coefficient_prior_sd=coefficient_prior_sd,
                    inference=inference,
                    confidence_level=confidence_level,
                    bootstrap_replicates=bootstrap_replicates,
                    seed=seed,
                    intercept=intercept,
                    allow_large_dense=allow_large_dense,
                )
            )
            continue
        rows.extend(
            _fit_ordinary_gaussian_response(
                tree,
                response,
                response_index,
                response_values_by_trait,
                predictors,
                leaf_names,
                design,
                term_names,
                term_metadata,
                covariance_by_trait,
                predictor_uncertainty_values,
                predictor_columns,
                predictor_diagnostics,
                mean_predictor_sampling_variance,
                evolution_model=evolution_model,
                evolution_parameter=evolution_parameter,
                branch_length=branch_length,
                custom_covariance=custom_covariance,
                inference=inference,
                bootstrap_replicates=bootstrap_replicates,
                seed=seed,
                reml=reml,
                intercept=intercept,
                confidence_level=confidence_level,
                num_parameters=num_parameters,
                matrix_rank=matrix_rank,
                allow_large_dense=allow_large_dense,
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
    categorical_predictors=(),
    ordered_predictors=None,
    factor_references=None,
    factor_coding="treatment",
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
    leaf_names = [str(leaf.name) for leaf in tree.leaves()]
    initial_encoding = encode_predictors(
        predictor_values_by_trait,
        predictors,
        leaf_names,
        categorical=categorical_predictors,
        ordered_levels=ordered_predictors,
        factor_references=factor_references,
        factor_coding=factor_coding,
    )
    continuous_predictors = [
        term.source for term in initial_encoding.terms if term.kind == "continuous"
    ]
    (
        predictor_values_by_trait,
        predictor_uncertainties,
        _,
    ) = _prepare_latent_ordinary_predictors(
        tree,
        predictor_values_by_trait,
        continuous_predictors,
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
        _,
        encoded_predictors,
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
        categorical_predictors=categorical_predictors,
        ordered_predictors=ordered_predictors,
        factor_references=factor_references,
        factor_coding=factor_coding,
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
    offset = 1 if intercept else 0
    term_index = {
        term.name: index + offset for index, term in enumerate(encoded_predictors.terms)
    }
    predictor_uncertainty_values = []
    predictor_columns: list[int | tuple[int, ...]] = []
    for predictor, covariance in predictor_uncertainties.items():
        predictor_uncertainty_values.append(covariance)
        predictor_columns.append(term_index[predictor])
    for uncertainty in encoded_predictors.uncertainties:
        predictor_uncertainty_values.append(uncertainty.covariance_by_observation)
        predictor_columns.append(
            tuple(term_index[term] for term in uncertainty.term_names)
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
                predictor_uncertainties=predictor_uncertainty_values,
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
        "categorical_predictors": None,
        "categorical_responses": None,
        "categorical_replicate_policy": "error",
        "factor_coding": "treatment",
        "factor_reference": None,
        "ordered_predictors": None,
        "ordered_responses": None,
        "response_reference": None,
        "response_family": None,
        "response_offset": None,
        "response_trials": None,
        "response_censor_lower": None,
        "response_censor_upper": None,
        "response_dispersion": None,
        "response_zero_probability": None,
        "coefficient_penalty": "student-t",
        "coefficient_prior_sd": 2.5,
        "multivariate_responses": False,
        "allow_missing_responses": False,
        "allow_large_dense": False,
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


def _ordinary_predictor_values(
    args,
    tree,
    predictors,
    duplicate_leaf_names,
    *,
    categorical=(),
    ordered=(),
    allow_missing=False,
):
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
        observed = dataframe[predictor].dropna()
        numeric_observed = pd.to_numeric(observed, errors="coerce")
        categorical_predictor = (
            predictor in set(categorical)
            or predictor in set(ordered)
            or numeric_observed.isna().any()
        )
        values_by_leaf = {}
        for leaf_name in tree.leaf_names():
            selected = dataframe.loc[
                dataframe["leaf_name"] == str(leaf_name), predictor
            ]
            if selected.empty or (selected.isna().any() and not allow_missing):
                raise ValueError(
                    "Predictor '{}' has missing values for tree tip '{}'.".format(
                        predictor, leaf_name
                    )
                )
            if categorical_predictor:
                if selected.isna().any():
                    raise ValueError(
                        "Categorical trait '{}' cannot contain missing values.".format(
                            predictor
                        )
                    )
                values = selected.astype(str).to_numpy()
            else:
                try:
                    values = pd.to_numeric(selected, errors="raise").to_numpy(float)
                except (TypeError, ValueError) as exc:
                    raise ValueError(
                        "Predictor '{}' must be numeric for tree tip '{}'.".format(
                            predictor, leaf_name
                        )
                    ) from exc
                if not allow_missing and not np.isfinite(values).all():
                    raise ValueError(
                        "Predictor '{}' must be finite for tree tip '{}'.".format(
                            predictor, leaf_name
                        )
                    )
            finite_values = values[
                np.asarray([not pd.isna(value) for value in values], dtype=bool)
            ]
            if len(finite_values) > 1 and np.any(finite_values != finite_values[0]):
                raise ValueError(
                    "Species-level predictor '{}' differs among rows for tree tip '{}'.".format(
                        predictor, leaf_name
                    )
                )
            values_by_leaf[str(leaf_name)] = (
                str(finite_values[0])
                if categorical_predictor
                else float("nan")
                if not len(finite_values)
                else float(finite_values[0])
            )
        values_by_trait[predictor] = values_by_leaf
    return values_by_trait


def _ordinary_auxiliary_values(
    args,
    tree,
    columns,
    duplicate_leaf_names,
    *,
    allow_missing=False,
):
    if not columns:
        return {}
    dataframe, _, _ = read_tip_table(
        args.data,
        option_name="--data",
        tree_leaf_names=list(tree.leaf_names()),
        required_columns=columns,
        unmatched=args.unmatched,
        missing_values=args.missing_values,
        duplicate_leaf_names=duplicate_leaf_names,
    )
    dataframe = dataframe[dataframe["leaf_name"].isin(set(tree.leaf_names()))]
    values_by_column = {}
    for column in columns:
        values_by_leaf = {}
        for leaf_name in tree.leaf_names():
            selected_rows = dataframe.loc[dataframe["leaf_name"] == str(leaf_name)]
            biological_id = getattr(args, "biological_id", None)
            if biological_id is None:
                numeric = pd.to_numeric(
                    selected_rows[column], errors="coerce"
                ).to_numpy(float)
            else:
                observations = []
                for _identifier, group in selected_rows.groupby(
                    biological_id, sort=False, dropna=False
                ):
                    group_values = pd.to_numeric(
                        group[column], errors="coerce"
                    ).to_numpy(float)
                    finite_group = group_values[np.isfinite(group_values)]
                    if len(finite_group) and np.any(finite_group != finite_group[0]):
                        raise ValueError(
                            "Auxiliary response column '{}' differs among technical "
                            "rows for tree tip '{}'.".format(column, leaf_name)
                        )
                    observations.append(
                        float("nan")
                        if not len(finite_group)
                        else float(finite_group[0])
                    )
                numeric = np.asarray(observations, dtype=float)
            if not allow_missing and (
                not len(numeric) or not np.isfinite(numeric).all()
            ):
                raise ValueError(
                    "Auxiliary response column '{}' must be finite for tree tip '{}'.".format(
                        column, leaf_name
                    )
                )
            finite = numeric[np.isfinite(numeric)]
            if biological_id is None and len(finite) and np.any(finite != finite[0]):
                raise ValueError(
                    "Auxiliary response column '{}' differs among rows for tree tip '{}'.".format(
                        column, leaf_name
                    )
                )
            values_by_leaf[str(leaf_name)] = (
                ReplicatedObservation(tuple(numeric))
                if biological_id is not None
                else float("nan")
                if not len(finite)
                else float(finite[0])
            )
        values_by_column[column] = values_by_leaf
    return values_by_column


def _response_auxiliary_mapping(column_by_response, values_by_column):
    return {
        response: values_by_column[column]
        for response, column in column_by_response.items()
    }


def _sampling_covariance_table(covariance_by_trait, leaf_names):
    rows = []
    for trait, covariance in covariance_by_trait.items():
        matrix = (
            covariance.loc[leaf_names, leaf_names].to_numpy(dtype=float)
            if isinstance(covariance, pd.DataFrame)
            else np.asarray(covariance, dtype=float)
        )
        for first, first_name in enumerate(leaf_names):
            second_indices = (
                range(first, first + 1)
                if matrix.ndim == 1
                else range(first, len(leaf_names))
            )
            for second in second_indices:
                rows.append(
                    {
                        "tree_id": "species",
                        "trait": trait,
                        "leaf_name_1": first_name,
                        "leaf_name_2": leaf_names[second],
                        "sampling_covariance": (
                            matrix[first] if matrix.ndim == 1 else matrix[first, second]
                        ),
                    }
                )
    return pd.DataFrame(rows, columns=ORDINARY_SAMPLING_COVARIANCE_COLUMNS)


def build_ordinary_pgls(args, responses, predictors) -> OrdinaryPglsArtifacts:
    effective = _effective_ordinary_args(args)
    if effective.allow_missing_responses and not effective.multivariate_responses:
        raise ValueError(
            "--allow-missing-responses requires '--multivariate-responses yes'."
        )
    categorical_responses = parse_name_list(effective.categorical_responses)
    ordered_responses = parse_ordered_levels(
        effective.ordered_responses, "--ordered-responses"
    )
    response_references = parse_key_values(
        effective.response_reference, "--response-reference"
    )
    response_families = parse_key_values(effective.response_family, "--response-family")
    response_offset_columns = parse_key_values(
        effective.response_offset, "--response-offset"
    )
    response_trial_columns = parse_key_values(
        effective.response_trials, "--response-trials"
    )
    response_lower_columns = parse_key_values(
        effective.response_censor_lower, "--response-censor-lower"
    )
    response_upper_columns = parse_key_values(
        effective.response_censor_upper, "--response-censor-upper"
    )
    response_dispersions = parse_numeric_key_values(
        effective.response_dispersion, "--response-dispersion", lower=0.0
    )
    response_zero_probabilities = parse_numeric_key_values(
        effective.response_zero_probability,
        "--response-zero-probability",
        lower=0.0,
        upper=1.0,
    )
    categorical_predictors = parse_name_list(effective.categorical_predictors)
    ordered_predictors = parse_ordered_levels(effective.ordered_predictors)
    factor_references = parse_key_values(
        effective.factor_reference, "--factor-reference"
    )
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
        replicate_estimates = _read_mixed_replicate_traits(
            effective,
            tree,
            responses,
            set(categorical_responses) | set(ordered_responses),
            "species",
            categorical_policy="counts",
            non_gaussian_columns={
                response
                for response, family in response_families.items()
                if family in SCALAR_RESPONSE_FAMILIES
            },
            allow_missing_columns=(
                (
                    {response for response in responses}
                    if effective.multivariate_responses
                    and effective.allow_missing_responses
                    else set()
                )
                | {
                    response
                    for response, family in response_families.items()
                    if family == "censored-gaussian"
                }
            ),
            option_name="--data",
        )
        response_values = replicate_estimates.values_by_trait
        covariance_by_trait = replicate_estimates.sampling_covariance_by_trait
        tip_summary = replicate_estimates.tip_summary
    else:
        replicate_estimates = None
        response_values = _ordinary_predictor_values(
            effective,
            tree,
            responses,
            duplicate_leaf_names=(
                "allow" if predictor_replicate_requested else "error"
            ),
            categorical=categorical_responses,
            ordered=ordered_responses,
            allow_missing=(
                effective.allow_missing_responses
                or any(
                    family == "censored-gaussian"
                    for family in response_families.values()
                )
            ),
        )
        covariance_by_trait = {}
        tip_summary = pd.DataFrame(columns=TIP_SUMMARY_COLUMNS)
    if predictor_replicate_requested:
        predictor_estimates = _read_mixed_replicate_traits(
            predictor_args,
            tree,
            predictors,
            set(categorical_predictors) | set(ordered_predictors),
            "species",
            categorical_policy=effective.categorical_replicate_policy,
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
            categorical=categorical_predictors,
            ordered=ordered_predictors,
        )
        predictor_covariance_by_trait = {}
        predictor_tip_summary = pd.DataFrame(columns=TIP_SUMMARY_COLUMNS)
    leaf_names = [str(leaf.name) for leaf in tree.leaves()]
    comparison_models = _parse_comparison_models(
        getattr(effective, "compare_evolution_models", None)
    )
    if effective.multivariate_responses and comparison_models:
        raise ValueError(
            "--compare-evolution-models is not available with --multivariate-responses."
        )
    covariance_path = getattr(effective, "evolution_covariance", None)
    custom_covariance = (
        None
        if covariance_path is None
        else read_custom_covariance(covariance_path, leaf_names)
    )
    response_specs = resolve_response_specs(
        response_values,
        responses,
        leaf_names,
        categorical=categorical_responses,
        ordered_levels=ordered_responses,
        references=response_references,
        families=response_families,
        allow_missing=(
            effective.allow_missing_responses
            or any(
                family == "censored-gaussian" for family in response_families.values()
            )
        ),
    )
    validate_response_auxiliaries(
        response_specs,
        offsets=response_offset_columns,
        trials=response_trial_columns,
        censor_lower=response_lower_columns,
        censor_upper=response_upper_columns,
        dispersions=response_dispersions,
        zero_probabilities=response_zero_probabilities,
    )
    if effective.coefficient_prior_sd <= 0.0:
        raise ValueError("--coefficient-prior-sd must be positive.")
    duplicate_policy = (
        "allow"
        if response_replicate_requested or predictor_replicate_requested
        else "error"
    )
    regular_auxiliary = _ordinary_auxiliary_values(
        effective,
        tree,
        sorted(
            set(response_offset_columns.values()) | set(response_trial_columns.values())
        ),
        duplicate_policy,
    )
    censor_auxiliary = _ordinary_auxiliary_values(
        effective,
        tree,
        sorted(
            set(response_lower_columns.values()) | set(response_upper_columns.values())
        ),
        duplicate_policy,
        allow_missing=True,
    )
    response_offsets = _response_auxiliary_mapping(
        response_offset_columns, regular_auxiliary
    )
    response_trials = _response_auxiliary_mapping(
        response_trial_columns, regular_auxiliary
    )
    response_censor_lower = _response_auxiliary_mapping(
        response_lower_columns, censor_auxiliary
    )
    response_censor_upper = _response_auxiliary_mapping(
        response_upper_columns, censor_auxiliary
    )
    if comparison_models and any(
        spec.family != "gaussian" for spec in response_specs.values()
    ):
        raise ValueError(
            "--compare-evolution-models is currently available for Gaussian responses only."
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
        categorical_predictors=categorical_predictors,
        ordered_predictors=ordered_predictors,
        factor_references=factor_references,
        factor_coding=effective.factor_coding,
        categorical_responses=categorical_responses,
        ordered_responses=ordered_responses,
        response_references=response_references,
        response_families=response_families,
        response_offsets=response_offsets,
        response_trials=response_trials,
        response_censor_lower=response_censor_lower,
        response_censor_upper=response_censor_upper,
        response_dispersions=response_dispersions,
        response_zero_probabilities=response_zero_probabilities,
        coefficient_penalty=effective.coefficient_penalty,
        coefficient_prior_sd=effective.coefficient_prior_sd,
        multivariate_responses=effective.multivariate_responses,
        allow_missing_responses=effective.allow_missing_responses,
        allow_large_dense=effective.allow_large_dense,
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
            categorical_predictors=categorical_predictors,
            ordered_predictors=ordered_predictors,
            factor_references=factor_references,
            factor_coding=effective.factor_coding,
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
