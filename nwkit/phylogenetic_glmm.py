"""Laplace-approximated phylogenetic generalized linear mixed models."""

import math
import warnings
from dataclasses import dataclass, replace
from typing import Callable, Mapping, Sequence

import numpy as np
from scipy import sparse
from scipy.optimize import brentq, linprog, minimize
from scipy.special import (
    betaln,
    digamma,
    expit,
    gammaln,
    log_ndtr,
    logsumexp,
    polygamma,
)
from scipy.stats import chi2
from scipy.stats import nbinom as nbinom_distribution
from scipy.stats import poisson as poisson_distribution

from nwkit.measurement_error import _finite_difference_hessian
from nwkit.model_matrix import CategoricalObservation, ReplicatedObservation
from nwkit.sparse_laplace import (
    ContinuousPredictorUncertainty,
    GmrfPredictorUncertainty,
    GroupedPredictorUncertainty,
    JointPredictorUncertainty,
    SparseCovarianceModel,
    SparseLatentModel,
    append_identity_latent_components,
    append_latent_component,
    combine_sparse_covariance_models,
    continuous_predictor_loading,
    factor_sparse_positive_definite,
    gmrf_predictor_loading,
    grouped_predictor_loading,
    joint_predictor_loading,
    prepare_sparse_latent_sampler,
)


@dataclass(frozen=True)
class PhylogeneticGlmmFit:
    """Fitted categorical phylogenetic GLMM and its inferential state."""

    family: str
    levels: tuple[str, ...]
    reference: str
    coefficients: np.ndarray
    coefficient_covariance: np.ndarray
    thresholds: np.ndarray
    threshold_covariance: np.ndarray
    random_modes: np.ndarray
    component_random_modes: Mapping[str, np.ndarray]
    component_variances: Mapping[str, float]
    evolution_parameter: float | None
    evolution_parameter_status: str
    log_likelihood: float
    optimizer_converged: bool
    optimizer_message: str
    boundary_warning: bool
    dispersion: float | None = None
    zero_probability: float | None = None
    coefficient_penalty: str = "none"
    coefficient_prior_sd: float | None = None
    separation_warning: bool = False
    coefficient_statistics: np.ndarray | None = None
    coefficient_p_values: np.ndarray | None = None
    coefficient_confidence_lower: np.ndarray | None = None
    coefficient_confidence_upper: np.ndarray | None = None
    coefficient_inference: str = "wald"
    coefficient_covariance_status: str = "ok"


SCALAR_RESPONSE_FAMILIES = {
    "poisson",
    "negative-binomial",
    "zero-inflated-poisson",
    "zero-inflated-negative-binomial",
    "hurdle-poisson",
    "hurdle-negative-binomial",
    "gamma",
    "lognormal",
    "beta",
    "beta-binomial",
    "censored-gaussian",
}

DISPERSION_FAMILIES = {
    "negative-binomial",
    "zero-inflated-negative-binomial",
    "hurdle-negative-binomial",
    "gamma",
    "lognormal",
    "beta",
    "beta-binomial",
    "censored-gaussian",
}

ZERO_COMPONENT_FAMILIES = {
    "zero-inflated-poisson",
    "zero-inflated-negative-binomial",
    "hurdle-poisson",
    "hurdle-negative-binomial",
}

MAX_DENSE_GLMM_TIPS = 500
MAX_SPARSE_GLMM_TIPS = 5000
MAX_SPARSE_GLMM_LINEAR_PREDICTORS = 20000


def _estimated_dense_glmm_gib(n_tips: int, random_dimension: int) -> float:
    """Conservative working-memory estimate for dense GLMM factorizations."""
    dimension = int(n_tips) * int(random_dimension)
    return 8.0 * 8.0 * dimension * dimension / 1024.0**3


def _select_sparse_glmm_backend(
    n_tips: int,
    random_dimension: int,
    sparse_capable: bool,
    allow_large_dense: bool,
) -> bool:
    """Select only the numerical backend; both paths use the same likelihood."""
    if n_tips <= MAX_DENSE_GLMM_TIPS:
        return False
    if sparse_capable:
        if n_tips > MAX_SPARSE_GLMM_TIPS:
            warnings.warn(
                "Sparse phylogenetic GLMM fitting above {} tips is not yet "
                "included in the routine validation suite (received {}); "
                "attempting the fit with the same sparse likelihood backend.".format(
                    MAX_SPARSE_GLMM_TIPS, n_tips
                ),
                RuntimeWarning,
                stacklevel=3,
            )
        return True
    estimated_gib = _estimated_dense_glmm_gib(n_tips, random_dimension)
    message = (
        "This GLMM cannot use the sparse backend and its dense latent dimension "
        "is {}. Dense fitting is estimated to need about {:.2f} GiB of working "
        "memory plus cubic-time factorizations."
    ).format(n_tips * random_dimension, estimated_gib)
    if not allow_large_dense:
        raise ValueError(
            message
            + " Pass allow_large_dense=True (CLI: '--allow-large-dense yes') to "
            "attempt the allocation."
        )
    warnings.warn(
        message + " Large dense allocation was explicitly enabled; attempting it.",
        RuntimeWarning,
        stacklevel=3,
    )
    return False


def _dense_glmm_memory_error(n_tips: int, random_dimension: int) -> MemoryError:
    estimated_gib = _estimated_dense_glmm_gib(n_tips, random_dimension)
    return MemoryError(
        "Dense phylogenetic GLMM allocation failed for latent dimension {}. "
        "The estimated working-memory requirement was about {:.2f} GiB; actual "
        "requirements depend on the numerical library and model components.".format(
            n_tips * random_dimension, estimated_gib
        )
    )


def _wald_coefficient_summary(fit, flat_index, coefficient, critical):
    variance = float(fit.coefficient_covariance[flat_index, flat_index])
    available = np.isfinite(variance) and variance >= 0.0
    if available and variance > 0.0:
        standard_error: float | str = math.sqrt(variance)
        statistic_value = coefficient / standard_error
        statistic: float | str = statistic_value
        p_value: float | str = float(chi2.sf(statistic_value**2, 1))
        lower: float | str = coefficient - critical * standard_error
        upper: float | str = coefficient + critical * standard_error
        status = "ok"
    else:
        standard_error = 0.0 if available else ""
        statistic = p_value = lower = upper = ""
        status = (
            "zero-model-variance" if available else fit.coefficient_covariance_status
        )
    return standard_error, statistic, p_value, lower, upper, status


def summarize_glmm_coefficient(fit, flat_index, coefficient, critical):
    """Return one coefficient's inference fields with singularity propagation."""
    standard_error, statistic, p_value, lower, upper, status = (
        _wald_coefficient_summary(fit, flat_index, coefficient, critical)
    )
    if fit.coefficient_statistics is not None:
        statistic = float(fit.coefficient_statistics[flat_index])
        assert fit.coefficient_p_values is not None
        p_value = float(fit.coefficient_p_values[flat_index])
    if fit.coefficient_confidence_lower is not None and np.isfinite(
        fit.coefficient_confidence_lower[flat_index]
    ):
        lower = float(fit.coefficient_confidence_lower[flat_index])
    if fit.coefficient_confidence_upper is not None and np.isfinite(
        fit.coefficient_confidence_upper[flat_index]
    ):
        upper = float(fit.coefficient_confidence_upper[flat_index])
    return standard_error, statistic, p_value, lower, upper, status


def summarize_glmm_omnibus(fit, indices):
    """Return a Wald omnibus result or an explicit unavailable status."""
    selected = fit.coefficients.reshape(-1)[indices]
    covariance = fit.coefficient_covariance[np.ix_(indices, indices)]
    if not np.isfinite(covariance).all():
        return "", "", fit.coefficient_covariance_status
    statistic = float(selected @ np.linalg.pinv(covariance) @ selected)
    return statistic, float(chi2.sf(statistic, len(indices))), "ok"


def summarize_glmm_threshold(fit, threshold_index, critical):
    """Return one ordinal threshold's uncertainty fields."""
    threshold = float(fit.thresholds[threshold_index])
    variance = float(fit.threshold_covariance[threshold_index, threshold_index])
    if not np.isfinite(variance) or variance < 0.0:
        return "", "", "", fit.coefficient_covariance_status
    standard_error = math.sqrt(variance)
    return (
        standard_error,
        threshold - critical * standard_error,
        threshold + critical * standard_error,
        "ok",
    )


def observed_response_levels(values: Sequence[object]) -> list[str]:
    """Return globally sorted response levels, including replicate counts."""
    levels: set[str] = set()
    for value in values:
        if isinstance(value, CategoricalObservation):
            levels.update(str(level) for level in value.probabilities)
        else:
            levels.add(str(value))
    return sorted(levels)


def response_count_matrix(
    values: Sequence[object], levels: Sequence[str]
) -> np.ndarray:
    """Encode exact states or aggregated biological replicates as counts."""
    level_index = {level: index for index, level in enumerate(levels)}
    counts = np.zeros((len(values), len(levels)), dtype=float)
    for row, value in enumerate(values):
        if isinstance(value, CategoricalObservation):
            for level, probability in value.probabilities.items():
                label = str(level)
                if label not in level_index:
                    raise ValueError(
                        "Categorical response contains an undeclared level '{}'.".format(
                            label
                        )
                    )
                count = float(probability) * int(value.n_observations)
                if not np.isfinite(count) or count < 0.0:
                    raise ValueError(
                        "Categorical response counts must be non-negative."
                    )
                counts[row, level_index[label]] += count
        else:
            label = str(value)
            if label not in level_index:
                raise ValueError(
                    "Categorical response contains an undeclared level '{}'.".format(
                        label
                    )
                )
            counts[row, level_index[label]] = 1.0
    if np.any(counts.sum(axis=1) <= 0.0):
        raise ValueError("Every tree tip needs a categorical response observation.")
    return counts


def _positive_definite_cholesky(matrix: np.ndarray) -> np.ndarray:
    symmetric = (matrix + matrix.T) / 2.0
    scale = max(1.0, float(np.max(np.diag(symmetric))))
    for multiplier in (0.0, 1e-12, 1e-10, 1e-8, 1e-6):
        try:
            return np.linalg.cholesky(
                symmetric + np.eye(len(symmetric)) * scale * multiplier
            )
        except np.linalg.LinAlgError:
            continue
    raise np.linalg.LinAlgError("Matrix is not positive definite.")


def _valid_optimizer_result(result) -> bool:
    return bool(result.success and np.isfinite(result.fun) and float(result.fun) < 1e99)


def _bounded_perturbed_starts(initial, first, bounds):
    starts = [np.asarray(initial, dtype=float), np.asarray(first.x, dtype=float)]
    for direction in (-1.0, 1.0):
        perturbed = np.asarray(initial, dtype=float).copy()
        perturbed += direction * 0.25 * np.maximum(1.0, np.abs(perturbed))
        for index, (lower, upper) in enumerate(bounds):
            if lower is not None:
                perturbed[index] = max(perturbed[index], float(lower))
            if upper is not None:
                perturbed[index] = min(perturbed[index], float(upper))
        starts.append(perturbed)
    return starts


def _powell_fallback_candidates(objective, starts, bounds, maxiter):
    results = [
        minimize(
            objective,
            start,
            method="Powell",
            bounds=bounds,
            options={"maxiter": maxiter * 2, "xtol": 1e-7, "ftol": 1e-9},
        )
        for start in starts
    ]
    return results, [result for result in results if _valid_optimizer_result(result)]


def _lbfgs_multistart_candidates(objective, starts, bounds, maxiter):
    results = [
        minimize(
            objective,
            start,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": maxiter, "ftol": 1e-10, "gtol": 1e-6},
        )
        for start in starts
    ]
    return results, [result for result in results if _valid_optimizer_result(result)]


def _robust_minimize(
    objective: Callable[[np.ndarray], float],
    initial: np.ndarray,
    bounds: Sequence[tuple[float | None, float | None]],
    *,
    label: str,
    maxiter: int = 700,
    force_multistart: bool = False,
):
    """Minimize a Laplace objective and never return a failed optimization."""
    candidates = []
    fallback_results = []
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="invalid value encountered in subtract",
            category=RuntimeWarning,
            module=r"scipy\.optimize\._numdiff",
        )
        first = minimize(
            objective,
            np.asarray(initial, dtype=float),
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": maxiter, "ftol": 1e-10, "gtol": 1e-6},
        )
        initial_value = float(objective(np.asarray(initial, dtype=float)))
        if _valid_optimizer_result(first) and first.fun <= initial_value + 1e-8 * max(
            1.0, abs(initial_value)
        ):
            candidates.append(first)
        if force_multistart:
            multistart_results, multistart_candidates = _lbfgs_multistart_candidates(
                objective,
                _bounded_perturbed_starts(initial, first, bounds)[2:],
                bounds,
                maxiter,
            )
            fallback_results.extend(multistart_results)
            candidates.extend(multistart_candidates)
        if not candidates:
            fallback_results, fallback_candidates = _powell_fallback_candidates(
                objective,
                _bounded_perturbed_starts(initial, first, bounds),
                bounds,
                maxiter,
            )
            candidates.extend(fallback_candidates)
            if fallback_candidates:
                fallback = min(candidates, key=lambda candidate: float(candidate.fun))
                refined = minimize(
                    objective,
                    fallback.x,
                    method="L-BFGS-B",
                    bounds=bounds,
                    options={"maxiter": maxiter, "ftol": 1e-10, "gtol": 1e-6},
                )
                if _valid_optimizer_result(refined):
                    candidates.append(refined)
    if not candidates:
        messages = "; ".join(
            str(result.message) for result in [first, *fallback_results]
        )
        raise RuntimeError("{} optimization failed: {}".format(label, messages))
    return min(candidates, key=lambda candidate: float(candidate.fun))


def _normalize_covariance(
    matrix: np.ndarray, size: int, label: str, *, allow_semidefinite: bool = False
) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=float)
    if matrix.shape != (size, size) or not np.isfinite(matrix).all():
        raise ValueError(
            "{} covariance has invalid dimensions or values.".format(label)
        )
    matrix = (matrix + matrix.T) / 2.0
    diagonal_mean = float(np.mean(np.diag(matrix)))
    if diagonal_mean <= 0.0:
        raise ValueError("{} covariance needs a positive diagonal.".format(label))
    normalized = matrix / diagonal_mean
    if allow_semidefinite:
        eigenvalues = np.linalg.eigvalsh(normalized)
        if float(np.min(eigenvalues)) < -1e-9 * max(1.0, float(np.max(eigenvalues))):
            raise ValueError(
                "{} covariance must be positive semidefinite.".format(label)
            )
    else:
        _positive_definite_cholesky(normalized)
    return normalized


def _multinomial_terms(
    counts: np.ndarray, linear: np.ndarray
) -> tuple[float, np.ndarray, np.ndarray]:
    zeros = np.zeros((len(linear), 1), dtype=float)
    all_linear = np.column_stack([linear, zeros])
    normalizer = logsumexp(all_linear, axis=1)
    probabilities = np.exp(all_linear - normalizer[:, None])
    totals = counts.sum(axis=1)
    log_likelihood = float(
        np.sum(counts[:, :-1] * linear) - np.sum(totals * normalizer)
    )
    gradient = totals[:, None] * probabilities[:, :-1] - counts[:, :-1]
    dimension = linear.shape[1]
    weights = np.empty((len(linear), dimension, dimension), dtype=float)
    for row, (total, probability) in enumerate(
        zip(totals, probabilities[:, :-1], strict=True)
    ):
        block = total * (np.diag(probability) - np.outer(probability, probability))
        weights[row] = block
    return log_likelihood, gradient.reshape(-1), weights


def _ordinal_thresholds(parameters: np.ndarray) -> np.ndarray:
    thresholds = np.empty(len(parameters), dtype=float)
    thresholds[0] = parameters[0]
    if len(parameters) > 1:
        thresholds[1:] = parameters[0] + np.cumsum(np.exp(parameters[1:]))
    return thresholds


def _ordinal_terms(
    counts: np.ndarray, linear: np.ndarray, thresholds: np.ndarray
) -> tuple[float, np.ndarray, np.ndarray]:
    negative_infinity = np.full(len(linear), -np.inf)
    positive_infinity = np.full(len(linear), np.inf)
    lower = np.column_stack([negative_infinity, thresholds[None, :] - linear[:, None]])
    upper = np.column_stack([thresholds[None, :] - linear[:, None], positive_infinity])
    lower_cdf = expit(lower)
    upper_cdf = expit(upper)
    probabilities = np.maximum(upper_cdf - lower_cdf, np.finfo(float).tiny)
    lower_density = lower_cdf * (1.0 - lower_cdf)
    upper_density = upper_cdf * (1.0 - upper_cdf)
    first_probability = lower_density - upper_density
    lower_derivative = lower_density * (1.0 - 2.0 * lower_cdf)
    upper_derivative = upper_density * (1.0 - 2.0 * upper_cdf)
    second_probability = upper_derivative - lower_derivative
    ratios = first_probability / probabilities
    gradient = -np.sum(counts * ratios, axis=1)
    curvature = np.sum(
        counts * (ratios**2 - second_probability / probabilities), axis=1
    )
    curvature = np.maximum(curvature, 1e-9)
    log_likelihood = float(np.sum(counts * np.log(probabilities)))
    return log_likelihood, gradient, curvature


def _dense_weight_matrix(weights: np.ndarray) -> np.ndarray:
    if weights.ndim == 1:
        return np.diag(weights)
    dimension = weights.shape[1]
    matrix = np.zeros((len(weights) * dimension,) * 2, dtype=float)
    for row, block in enumerate(weights):
        start = row * dimension
        matrix[start : start + dimension, start : start + dimension] = block
    return matrix


def _sparse_weight_matrix(weights: np.ndarray) -> sparse.csr_matrix:
    if weights.ndim == 1:
        return sparse.diags(weights, format="csr")
    n_blocks, dimension, _ = weights.shape
    block_offsets = np.arange(n_blocks, dtype=int) * dimension
    local_rows = np.repeat(np.arange(dimension, dtype=int), dimension)
    local_columns = np.tile(np.arange(dimension, dtype=int), dimension)
    rows = (block_offsets[:, None] + local_rows[None, :]).reshape(-1)
    columns = (block_offsets[:, None] + local_columns[None, :]).reshape(-1)
    matrix = sparse.csr_matrix(
        (weights.reshape(-1), (rows, columns)),
        shape=(n_blocks * dimension, n_blocks * dimension),
    )
    matrix.eliminate_zeros()
    return matrix


def _sparse_covariance_from_factory(factory, parameter):
    builder = getattr(factory, "sparse_model", None)
    return None if not callable(builder) else builder(parameter)


def _sparse_glmm_latent_model(
    phylogenetic_covariance,
    parameter,
    component_variances,
    group_covariance,
    random_dimension,
    predictor_uncertainties=(),
    predictor_columns=(),
    coefficients=None,
):
    phylogenetic = _sparse_covariance_from_factory(phylogenetic_covariance, parameter)
    if phylogenetic is None:
        return None
    components = {"phylogenetic": (component_variances["phylogenetic"], phylogenetic)}
    if group_covariance is not None:
        if not isinstance(group_covariance, SparseCovarianceModel):
            return None
        components["group"] = (component_variances["group"], group_covariance)
    latent = combine_sparse_covariance_models(
        components, random_dimension=random_dimension
    )
    if predictor_uncertainties:
        if coefficients is None:
            raise ValueError("Sparse predictor uncertainty requires coefficients.")
        for index, (uncertainty, selected_columns) in enumerate(
            zip(predictor_uncertainties, predictor_columns, strict=True)
        ):
            name = "predictor_uncertainty_{}".format(index)
            loading = _structured_predictor_loading(
                uncertainty, selected_columns, coefficients
            )
            if isinstance(uncertainty, GmrfPredictorUncertainty):
                latent = append_latent_component(
                    latent,
                    name,
                    uncertainty.model.precision,
                    loading,
                    uncertainty.model.logdet_covariance,
                    uncertainty.model.precision_factor(),
                )
            else:
                latent = append_identity_latent_components(latent, {name: loading})
    return latent


def _sparse_latent_mode(
    fixed_linear,
    loading,
    precision,
    likelihood_terms,
    *,
    label,
):
    mode = np.zeros(precision.shape[0], dtype=float)

    def terms(candidate):
        linear = fixed_linear + np.asarray(loading @ candidate).reshape(-1)
        return likelihood_terms(linear)

    def state(candidate):
        log_likelihood, likelihood_gradient, weights = terms(candidate)
        gradient = np.asarray(loading.T @ likelihood_gradient).reshape(-1)
        gradient += precision @ candidate
        weight_matrix = _sparse_weight_matrix(weights)
        hessian = precision + loading.T @ weight_matrix @ loading
        factor = factor_sparse_positive_definite(hessian)
        value = -log_likelihood + 0.5 * float(candidate @ (precision @ candidate))
        return value, gradient, weights, factor

    converged = False
    for _iteration in range(80):
        value, gradient, _weights, factor = state(mode)
        if not np.isfinite(gradient).all():
            break
        if float(np.max(np.abs(gradient))) < 1e-7:
            converged = True
            break
        step = factor.solve(gradient)
        accepted = False
        step_scale = 1.0
        for _backtrack in range(24):
            candidate = mode - step_scale * step
            candidate_log_likelihood, _, _ = terms(candidate)
            candidate_value = -candidate_log_likelihood + 0.5 * float(
                candidate @ (precision @ candidate)
            )
            if np.isfinite(candidate_value) and candidate_value <= value:
                mode = candidate
                accepted = True
                break
            step_scale *= 0.5
        if not accepted:
            break
        if float(np.max(np.abs(step_scale * step))) < 1e-8:
            converged = True
            break
    if not converged:

        def value_and_gradient(candidate):
            log_likelihood, likelihood_gradient, _weights = terms(candidate)
            gradient = np.asarray(loading.T @ likelihood_gradient).reshape(-1)
            gradient += precision @ candidate
            value = -log_likelihood + 0.5 * float(candidate @ (precision @ candidate))
            return value, gradient

        fallback = minimize(
            value_and_gradient,
            mode,
            method="L-BFGS-B",
            jac=True,
            options={"maxiter": 1000, "ftol": 1e-12, "gtol": 1e-7},
        )
        if not _valid_optimizer_result(fallback):
            raise RuntimeError(
                "{} random-effect mode optimization failed: {}".format(
                    label, fallback.message
                )
            )
        mode = np.asarray(fallback.x, dtype=float)
    value, gradient, weights, factor = state(mode)
    if not np.isfinite(value) or float(np.max(np.abs(gradient))) >= 1e-5:
        raise RuntimeError("{} random-effect mode did not converge.".format(label))
    return mode, value, weights, factor


def _sparse_marginal_coefficient_covariance(design, weights, loading, posterior_factor):
    weight_matrix = _sparse_weight_matrix(weights)
    weighted_design = weight_matrix @ design
    joint_fixed = design.T @ weighted_design
    fixed_random = np.asarray(design.T @ weight_matrix @ loading)
    adjustment = fixed_random @ posterior_factor.solve(fixed_random.T)
    information = np.asarray(joint_fixed - adjustment, dtype=float)
    information = (information + information.T) / 2.0
    return np.linalg.pinv(information, hermitian=True)


def _sparse_component_modes(
    latent: SparseLatentModel,
    state: np.ndarray,
    n_tips: int,
    random_dimension: int,
    component_names,
):
    modes = {}
    for name in component_names:
        selected = latent.component_slices[name]
        values = latent.component_loadings[name] @ state[selected]
        modes[name] = np.asarray(values).reshape(n_tips, random_dimension)
    return modes


def _random_mode(
    counts: np.ndarray,
    fixed_linear: np.ndarray,
    precision: np.ndarray,
    *,
    family: str,
    thresholds: np.ndarray,
    initial: np.ndarray | None = None,
) -> tuple[np.ndarray, float, np.ndarray]:
    mode = (
        np.zeros(precision.shape[0], dtype=float)
        if initial is None
        else np.asarray(initial, dtype=float).copy()
    )

    def terms(candidate: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
        if family == "ordinal":
            return _ordinal_terms(counts, fixed_linear + candidate, thresholds)
        dimension = fixed_linear.shape[1]
        return _multinomial_terms(
            counts, fixed_linear + candidate.reshape(len(counts), dimension)
        )

    converged = False
    for _iteration in range(80):
        log_likelihood, likelihood_gradient, weights = terms(mode)
        gradient = likelihood_gradient + precision @ mode
        if not np.isfinite(gradient).all():
            break
        if float(np.max(np.abs(gradient))) < 1e-7:
            converged = True
            break
        hessian = _dense_weight_matrix(weights) + precision
        try:
            cholesky = _positive_definite_cholesky(hessian)
        except np.linalg.LinAlgError:
            break
        step = np.linalg.solve(cholesky.T, np.linalg.solve(cholesky, gradient))
        current = -log_likelihood + 0.5 * float(mode @ precision @ mode)
        step_scale = 1.0
        accepted = False
        for _backtrack in range(24):
            candidate = mode - step_scale * step
            candidate_log_likelihood, _, _ = terms(candidate)
            candidate_objective = -candidate_log_likelihood + 0.5 * float(
                candidate @ precision @ candidate
            )
            if np.isfinite(candidate_objective) and candidate_objective <= current:
                mode = candidate
                accepted = True
                break
            step_scale *= 0.5
        if not accepted:
            break
        if float(np.max(np.abs(step_scale * step))) < 1e-8:
            converged = True
            break
    if converged:
        _log_likelihood, likelihood_gradient, _weights = terms(mode)
        converged = bool(
            np.isfinite(likelihood_gradient).all()
            and float(np.max(np.abs(likelihood_gradient + precision @ mode))) < 1e-5
        )
    if not converged:

        def value_and_gradient(candidate: np.ndarray):
            log_likelihood, likelihood_gradient, _weights = terms(candidate)
            value = -log_likelihood + 0.5 * float(candidate @ precision @ candidate)
            gradient = likelihood_gradient + precision @ candidate
            return value, gradient

        fallback = minimize(
            value_and_gradient,
            mode,
            method="L-BFGS-B",
            jac=True,
            options={"maxiter": 1000, "ftol": 1e-12, "gtol": 1e-7},
        )
        if (
            not fallback.success
            or not np.isfinite(fallback.fun)
            or float(fallback.fun) >= 1e99
        ):
            raise RuntimeError(
                "Categorical random-effect mode optimization failed: {}".format(
                    fallback.message
                )
            )
        mode = np.asarray(fallback.x, dtype=float)
    log_likelihood, _, weights = terms(mode)
    objective = -log_likelihood + 0.5 * float(mode @ precision @ mode)
    _positive_definite_cholesky(_dense_weight_matrix(weights) + precision)
    if not np.isfinite(objective):
        raise RuntimeError("Categorical random-effect mode is non-finite.")
    return mode, objective, weights


def _initial_coefficients(
    counts: np.ndarray, design: np.ndarray, family: str
) -> np.ndarray:
    columns = 1 if family == "ordinal" else counts.shape[1] - 1
    coefficients = np.zeros((design.shape[1], columns), dtype=float)
    if not np.allclose(design[:, 0], 1.0):
        return coefficients
    totals = counts.sum(axis=0) + 0.5
    if family == "ordinal":
        return coefficients
    coefficients[0] = np.log(totals[:-1] / totals[-1])
    return coefficients


def _initial_threshold_parameters(counts: np.ndarray) -> np.ndarray:
    cumulative = np.cumsum(counts.sum(axis=0))[:-1]
    total = float(np.sum(counts))
    probabilities = np.clip(cumulative / total, 1e-4, 1.0 - 1e-4)
    thresholds = np.log(probabilities / (1.0 - probabilities))
    parameters = np.empty_like(thresholds)
    parameters[0] = thresholds[0]
    if len(thresholds) > 1:
        parameters[1:] = np.log(np.maximum(np.diff(thresholds), 1e-3))
    return parameters


def _project_predictor_uncertainty(
    coefficients: np.ndarray,
    uncertainties: Sequence[np.ndarray],
    columns: Sequence[int | tuple[int, ...]],
    n_observations: int,
) -> np.ndarray:
    """Project predictor uncertainty into the multivariate linear predictor."""
    dimension = coefficients.shape[1]
    projected = np.zeros(
        (n_observations * dimension, n_observations * dimension), dtype=float
    )
    if len(uncertainties) != len(columns):
        raise ValueError("Predictor uncertainty and coefficient-column counts differ.")
    for uncertainty, selected_columns in zip(uncertainties, columns, strict=True):
        if isinstance(
            uncertainty,
            (
                ContinuousPredictorUncertainty,
                GmrfPredictorUncertainty,
                GroupedPredictorUncertainty,
                JointPredictorUncertainty,
            ),
        ):
            loading = _structured_predictor_loading(
                uncertainty, selected_columns, coefficients
            )
            if isinstance(uncertainty, GmrfPredictorUncertainty):
                factor = factor_sparse_positive_definite(uncertainty.model.precision)
                projected += np.asarray(loading @ factor.solve(loading.T.toarray()))
            else:
                projected += np.asarray((loading @ loading.T).toarray())
            continue
        covariance = np.asarray(uncertainty, dtype=float)
        if isinstance(selected_columns, int):
            if covariance.shape != (n_observations, n_observations):
                raise ValueError(
                    "Continuous predictor uncertainty has wrong dimensions."
                )
            slopes = coefficients[selected_columns]
            projected += np.kron(covariance, np.outer(slopes, slopes))
            continue
        selected = tuple(selected_columns)
        slopes = coefficients[np.asarray(selected, dtype=int), :]
        if covariance.shape == (n_observations, len(selected), len(selected)):
            for observation in range(n_observations):
                block = slopes.T @ covariance[observation] @ slopes
                start = observation * dimension
                projected[start : start + dimension, start : start + dimension] += block
            continue
        expected = (
            len(selected),
            len(selected),
            n_observations,
            n_observations,
        )
        if covariance.shape != expected:
            raise ValueError("Grouped predictor uncertainty has wrong dimensions.")
        for first in range(n_observations):
            first_start = first * dimension
            for second in range(n_observations):
                second_start = second * dimension
                block = slopes.T @ covariance[:, :, first, second] @ slopes
                projected[
                    first_start : first_start + dimension,
                    second_start : second_start + dimension,
                ] += block
    return (projected + projected.T) / 2.0


def _structured_predictor_loading(uncertainty, selected_columns, coefficients):
    if isinstance(uncertainty, ContinuousPredictorUncertainty):
        if not isinstance(selected_columns, int):
            raise ValueError("Continuous predictor uncertainty requires one column.")
        return continuous_predictor_loading(uncertainty, coefficients[selected_columns])
    if isinstance(uncertainty, GmrfPredictorUncertainty):
        if not isinstance(selected_columns, int):
            raise ValueError("GMRF predictor uncertainty requires one column.")
        return gmrf_predictor_loading(uncertainty, coefficients[selected_columns])
    if isinstance(uncertainty, GroupedPredictorUncertainty):
        if isinstance(selected_columns, int):
            raise ValueError("Grouped predictor uncertainty requires column indices.")
        selected = np.asarray(tuple(selected_columns), dtype=int)
        return grouped_predictor_loading(uncertainty, coefficients[selected])
    if isinstance(uncertainty, JointPredictorUncertainty):
        if isinstance(selected_columns, int):
            raise ValueError("Joint predictor uncertainty requires column indices.")
        selected = np.asarray(tuple(selected_columns), dtype=int)
        return joint_predictor_loading(uncertainty, coefficients[selected])
    raise TypeError("Unsupported structured predictor uncertainty.")


def _sparse_predictor_uncertainties(uncertainties):
    return all(
        isinstance(
            uncertainty,
            (
                ContinuousPredictorUncertainty,
                GmrfPredictorUncertainty,
                GroupedPredictorUncertainty,
                JointPredictorUncertainty,
            ),
        )
        for uncertainty in uncertainties
    )


def _coefficient_penalty_value(
    coefficients: np.ndarray, penalty: str, prior_sd: float | None
) -> float:
    if penalty == "none":
        return 0.0
    if prior_sd is None or not np.isfinite(prior_sd) or prior_sd <= 0.0:
        raise ValueError(
            "A positive coefficient prior SD is required for penalization."
        )
    scaled = np.asarray(coefficients, dtype=float).reshape(-1) / prior_sd
    if penalty == "gaussian":
        return 0.5 * float(scaled @ scaled)
    if penalty == "student-t":
        degrees_of_freedom = 3.0
        return (
            0.5
            * (degrees_of_freedom + 1.0)
            * float(np.sum(np.log1p(scaled**2 / degrees_of_freedom)))
        )
    raise ValueError("Unsupported coefficient penalty: {}.".format(penalty))


def _reported_log_likelihood(
    fitted_objective: float,
    coefficients: np.ndarray,
    penalty: str,
    prior_sd: float | None,
) -> float:
    """Remove coefficient regularization from the fitted Laplace objective."""
    penalty_value = _coefficient_penalty_value(coefficients, penalty, prior_sd)
    return -float(fitted_objective - penalty_value)


def _count_log_pmf(
    values: np.ndarray, means: np.ndarray, family: str, dispersion: float | None
) -> tuple[np.ndarray, np.ndarray]:
    if "negative-binomial" not in family:
        log_pmf = values * np.log(means) - means - gammaln(values + 1.0)
        return log_pmf, -means
    if dispersion is None or dispersion <= 0.0:
        raise ValueError("Negative-binomial dispersion must be positive.")
    size = 1.0 / dispersion
    log_pmf = (
        gammaln(values + size)
        - gammaln(size)
        - gammaln(values + 1.0)
        + size * (np.log(size) - np.log(size + means))
        + values * (np.log(means) - np.log(size + means))
    )
    log_zero = size * (np.log(size) - np.log(size + means))
    return log_pmf, log_zero


def _count_log_likelihood(values, linear, family, dispersion, zero_probability, offset):
    log_means = linear + offset
    with np.errstate(over="ignore", invalid="ignore"):
        means = np.exp(log_means)
    if "negative-binomial" not in family:
        log_pmf = values * log_means - means - gammaln(values + 1.0)
        log_zero = -means
    else:
        if dispersion is None or dispersion <= 0.0:
            raise ValueError("Negative-binomial dispersion must be positive.")
        size = 1.0 / dispersion
        log_size = math.log(size)
        log_denominator = np.logaddexp(log_size, log_means)
        log_pmf = (
            gammaln(values + size)
            - gammaln(size)
            - gammaln(values + 1.0)
            + size * (log_size - log_denominator)
            + values * (log_means - log_denominator)
        )
        log_zero = size * (log_size - log_denominator)
    if family in {"poisson", "negative-binomial"}:
        return log_pmf
    if zero_probability is None or not 0.0 < zero_probability < 1.0:
        raise ValueError("Zero-component probability must lie strictly in (0, 1).")
    log_pi = math.log(zero_probability)
    log_one_minus_pi = math.log1p(-zero_probability)
    is_zero = values == 0.0
    if family.startswith("zero-inflated"):
        return np.where(
            is_zero,
            np.logaddexp(log_pi, log_one_minus_pi + log_zero),
            log_one_minus_pi + log_pmf,
        )
    positive_normalizer = np.log(np.maximum(-np.expm1(log_zero), 1e-300))
    return np.where(
        is_zero,
        log_pi,
        log_one_minus_pi + log_pmf - positive_normalizer,
    )


def _gamma_log_likelihood(values, linear, dispersion):
    if dispersion is None or dispersion <= 0.0:
        raise ValueError("Gamma shape must be positive.")
    with np.errstate(over="ignore", invalid="ignore"):
        inverse_means = np.exp(-linear)
        return (
            dispersion * np.log(dispersion)
            - gammaln(dispersion)
            + (dispersion - 1.0) * np.log(values)
            - dispersion * linear
            - dispersion * values * inverse_means
        )


def _lognormal_log_likelihood(values, linear, dispersion):
    if dispersion is None or dispersion <= 0.0:
        raise ValueError("Lognormal SD must be positive.")
    standardized = (np.log(values) - linear) / dispersion
    return (
        -0.5 * standardized**2
        - np.log(values)
        - math.log(dispersion)
        - 0.5 * math.log(2.0 * math.pi)
    )


def _beta_parameters(linear, dispersion, label):
    if dispersion is None or dispersion <= 0.0:
        raise ValueError("{} precision must be positive.".format(label))
    means = np.clip(expit(linear), 1e-10, 1.0 - 1e-10)
    return means * dispersion, (1.0 - means) * dispersion


def _beta_log_likelihood(values, linear, dispersion):
    alpha, beta = _beta_parameters(linear, dispersion, "Beta")
    return (
        -betaln(alpha, beta)
        + (alpha - 1.0) * np.log(values)
        + (beta - 1.0) * np.log1p(-values)
    )


def _beta_binomial_log_likelihood(values, linear, dispersion, trials):
    if trials is None:
        raise ValueError("Beta-binomial responses require trial counts.")
    alpha, beta = _beta_parameters(linear, dispersion, "Beta-binomial")
    return (
        gammaln(trials + 1.0)
        - gammaln(values + 1.0)
        - gammaln(trials - values + 1.0)
        + betaln(values + alpha, trials - values + beta)
        - betaln(alpha, beta)
    )


def _censored_gaussian_log_likelihood(
    values, linear, dispersion, censor_lower, censor_upper
):
    if dispersion is None or dispersion <= 0.0:
        raise ValueError("Censored-Gaussian SD must be positive.")
    lower = np.full(len(values), np.nan) if censor_lower is None else censor_lower
    upper = np.full(len(values), np.nan) if censor_upper is None else censor_upper
    exact = np.isnan(lower) & np.isnan(upper)
    left = np.isnan(lower) & ~np.isnan(upper)
    right = ~np.isnan(lower) & np.isnan(upper)
    interval = ~np.isnan(lower) & ~np.isnan(upper)
    result = np.empty(len(values), dtype=float)
    standardized = (values - linear) / dispersion
    result[exact] = (
        -0.5 * standardized[exact] ** 2
        - math.log(dispersion)
        - 0.5 * math.log(2.0 * math.pi)
    )
    result[left] = log_ndtr((upper[left] - linear[left]) / dispersion)
    result[right] = log_ndtr((linear[right] - lower[right]) / dispersion)
    if np.any(interval):
        standardized_lower = (lower[interval] - linear[interval]) / dispersion
        standardized_upper = (upper[interval] - linear[interval]) / dispersion
        interval_result = np.empty(len(standardized_lower), dtype=float)
        negative = standardized_upper <= 0.0
        positive = standardized_lower >= 0.0
        central = ~(negative | positive)

        def log_difference(larger: np.ndarray, smaller: np.ndarray) -> np.ndarray:
            ratio = np.exp(np.minimum(smaller - larger, 0.0))
            return larger + np.log1p(-ratio)

        if np.any(negative):
            interval_result[negative] = log_difference(
                log_ndtr(standardized_upper[negative]),
                log_ndtr(standardized_lower[negative]),
            )
        if np.any(positive):
            interval_result[positive] = log_difference(
                log_ndtr(-standardized_lower[positive]),
                log_ndtr(-standardized_upper[positive]),
            )
        if np.any(central):
            interval_result[central] = log_difference(
                log_ndtr(standardized_upper[central]),
                log_ndtr(standardized_lower[central]),
            )
        result[interval] = interval_result
    return result


def _scalar_log_likelihood_contributions(
    values: np.ndarray,
    linear: np.ndarray,
    *,
    family: str,
    dispersion: float | None,
    zero_probability: float | None,
    offset: np.ndarray,
    trials: np.ndarray | None,
    censor_lower: np.ndarray | None,
    censor_upper: np.ndarray | None,
) -> np.ndarray:
    if family in {
        "poisson",
        "negative-binomial",
        "zero-inflated-poisson",
        "zero-inflated-negative-binomial",
        "hurdle-poisson",
        "hurdle-negative-binomial",
    }:
        return _count_log_likelihood(
            values, linear, family, dispersion, zero_probability, offset
        )
    if family == "gamma":
        return _gamma_log_likelihood(values, linear, dispersion)
    if family == "lognormal":
        return _lognormal_log_likelihood(values, linear, dispersion)
    if family == "beta":
        return _beta_log_likelihood(values, linear, dispersion)
    if family == "beta-binomial":
        return _beta_binomial_log_likelihood(values, linear, dispersion, trials)
    if family == "censored-gaussian":
        return _censored_gaussian_log_likelihood(
            values, linear, dispersion, censor_lower, censor_upper
        )
    raise ValueError("Unsupported scalar response family: {}.".format(family))


def _poisson_scalar_terms(values, linear, offset):
    with np.errstate(over="ignore"):
        means = np.exp(linear + offset)
    return means - values, means


def _negative_binomial_scalar_terms(values, linear, offset, dispersion):
    assert dispersion is not None
    size = 1.0 / dispersion
    log_means = linear + offset
    log_denominator = np.logaddexp(math.log(size), log_means)
    mean_fraction = expit(log_means - math.log(size))
    inverse_denominator = np.exp(-log_denominator)
    gradient = size * (mean_fraction - values * inverse_denominator)
    curvature = size * mean_fraction * (size + values) * inverse_denominator
    return gradient, curvature


def _softplus_extended(values):
    """Extended-precision softplus, including sub-epsilon exponential tails.

    Some NumPy/platform long-double log1p implementations lose the entire
    leading term for tiny positive arguments (logaddexp uses it too). Evaluate
    log(1+t) by its alternating series there, rather than rounding 1+t to 1
    or reducing the rest of the calculation to double precision.
    """
    with np.errstate(under="ignore"):
        tail = np.exp(-np.abs(values))
    small = tail < np.longdouble("1e-4")
    log_tail = np.empty_like(tail)
    if np.any(small):
        small_tail = tail[small]
        series = np.zeros_like(small_tail)
        # Twelve terms bound relative truncation error below 1e-48 here,
        # including platforms where longdouble is IEEE quadruple precision.
        for order in range(12, 0, -1):
            coefficient = np.longdouble(1 if order % 2 else -1) / order
            series = coefficient + small_tail * series
        log_tail[small] = small_tail * series
    log_tail[~small] = np.log1p(tail[~small])
    return np.maximum(values, 0.0) + log_tail


def _zero_truncation_terms(eta, dispersion):
    """Stable zero-truncation score correction and its derivative.

    Evaluate q / expm1(x) through exp(-x) for large x. This avoids overflow
    while retaining expm1 accuracy at the small-mean limit, where t tends to 1.
    """
    eta = np.asarray(eta, dtype=np.longdouble)
    if dispersion is None:
        with np.errstate(over="ignore", under="ignore"):
            q = np.exp(eta)
        exponent = q
        leading_derivative = np.ones_like(q)
    else:
        size = np.longdouble(1.0) / np.longdouble(dispersion)
        log_relative_mean = eta - np.log(size)
        fraction = np.exp(-_softplus_extended(-log_relative_mean))
        q = size * fraction
        exponent = size * _softplus_extended(log_relative_mean)
        leading_derivative = 1.0 - fraction
    correction = np.zeros_like(q)
    small = exponent <= 1.0
    denominator = np.expm1(exponent[small])
    correction[small] = np.divide(
        q[small], denominator, out=np.ones_like(denominator), where=denominator != 0.0
    )
    large = ~small & np.isfinite(q)
    correction[large] = (
        q[large] * np.exp(-exponent[large]) / -np.expm1(-exponent[large])
    )
    derivative = np.zeros_like(q)
    positive = correction > 0.0
    derivative[positive] = correction[positive] * (
        leading_derivative[positive] - q[positive] - correction[positive]
    )
    return correction, derivative


def _zero_component_count_scalar_terms(
    values,
    linear,
    offset,
    dispersion,
    zero_probability,
    *,
    hurdle,
):
    """Exact score and curvature for zero-modified Poisson/NB responses."""
    assert zero_probability is not None
    eta = np.asarray(linear + offset, dtype=float)
    negative_binomial = dispersion is not None
    if negative_binomial:
        size = float(1.0 / dispersion)
        fraction = expit(eta - math.log(size))
        zero_log_probability = -size * np.logaddexp(0.0, eta - math.log(size))
        zero_score = -size * fraction
        zero_second = -size * fraction * (1.0 - fraction)
        base_gradient, base_curvature = _negative_binomial_scalar_terms(
            values, linear, offset, dispersion
        )
    else:
        with np.errstate(over="ignore"):
            means = np.exp(eta)
        zero_log_probability = -means
        zero_score = -means
        zero_second = -means
        base_gradient, base_curvature = _poisson_scalar_terms(values, linear, offset)

    zeros = values == 0.0
    gradient = np.asarray(base_gradient, dtype=float).copy()
    curvature = np.asarray(base_curvature, dtype=float).copy()
    if hurdle:
        gradient[zeros] = 0.0
        curvature[zeros] = 0.0
        positive = ~zeros
        if np.any(positive):
            # If t = -p0'/(1-p0), the truncated likelihood has gradient
            # base_gradient + t and curvature base_curvature + t'.  This form
            # and extended precision avoid subtracting two order-one terms
            # when the mean approaches zero.
            t, t_derivative = _zero_truncation_terms(eta[positive], dispersion)
            selected_gradient = (
                np.asarray(base_gradient[positive], dtype=np.longdouble) + t
            )
            selected_curvature = (
                np.asarray(base_curvature[positive], dtype=np.longdouble) + t_derivative
            )
            gradient[positive] = np.asarray(selected_gradient, dtype=float)
            curvature[positive] = np.asarray(selected_curvature, dtype=float)
            # A zero-truncated Poisson likelihood is convex in its natural
            # parameter.  Roundoff in the limiting eta -> -inf cancellation
            # must not manufacture a negative Hessian.
            if not negative_binomial:
                curvature[positive] = np.maximum(curvature[positive], 0.0)
        return gradient, curvature

    log_odds = math.log1p(-zero_probability) - math.log(zero_probability)
    selected_log_zero = zero_log_probability[zeros]
    posterior_count_zero = expit(log_odds + selected_log_zero)
    score_zero = zero_score[zeros]
    second_zero = zero_second[zeros]
    with np.errstate(over="ignore", invalid="ignore"):
        selected_gradient = -posterior_count_zero * score_zero
        selected_curvature = -(
            posterior_count_zero * second_zero
            + posterior_count_zero * (1.0 - posterior_count_zero) * score_zero**2
        )
    infinite_tail = np.isneginf(selected_log_zero)
    gradient[zeros] = np.where(infinite_tail, 0.0, selected_gradient)
    curvature[zeros] = np.where(infinite_tail, 0.0, selected_curvature)
    return gradient, curvature


def _gamma_scalar_terms(values, linear, dispersion):
    assert dispersion is not None
    with np.errstate(over="ignore"):
        scaled = values * np.exp(-linear)
    return dispersion * (1.0 - scaled), dispersion * scaled


def _lognormal_scalar_terms(values, linear, dispersion):
    assert dispersion is not None
    variance = dispersion**2
    return (linear - np.log(values)) / variance, np.full(len(values), 1.0 / variance)


def _beta_scalar_terms(values, linear, dispersion):
    assert dispersion is not None
    means = np.clip(expit(linear), 1e-10, 1.0 - 1e-10)
    alpha = means * dispersion
    beta = (1.0 - means) * dispersion
    likelihood_derivative = dispersion * (
        -digamma(alpha) + digamma(beta) + np.log(values) - np.log1p(-values)
    )
    likelihood_second = -(dispersion**2) * (polygamma(1, alpha) + polygamma(1, beta))
    first = means * (1.0 - means)
    second = first * (1.0 - 2.0 * means)
    gradient = -likelihood_derivative * first
    curvature = -(likelihood_second * first**2 + likelihood_derivative * second)
    return gradient, curvature


def _beta_binomial_scalar_terms(values, linear, dispersion, trials):
    assert dispersion is not None
    assert trials is not None
    means = np.clip(expit(linear), 1e-10, 1.0 - 1e-10)
    alpha = means * dispersion
    beta = (1.0 - means) * dispersion
    likelihood_mean_derivative = dispersion * (
        digamma(values + alpha)
        - digamma(alpha)
        - digamma(trials - values + beta)
        + digamma(beta)
    )
    likelihood_mean_second = dispersion**2 * (
        polygamma(1, values + alpha)
        - polygamma(1, alpha)
        + polygamma(1, trials - values + beta)
        - polygamma(1, beta)
    )
    first = means * (1.0 - means)
    second = first * (1.0 - 2.0 * means)
    gradient = -likelihood_mean_derivative * first
    curvature = -(
        likelihood_mean_second * first**2 + likelihood_mean_derivative * second
    )
    return gradient, curvature


def _numerical_scalar_terms(contributions, current, linear):
    step = 1e-5 * np.maximum(1.0, np.abs(linear))
    plus = contributions(linear + step)
    minus = contributions(linear - step)
    with np.errstate(invalid="ignore", over="ignore", divide="ignore"):
        gradient = -(plus - minus) / (2.0 * step)
        curvature = -(plus - 2.0 * current + minus) / (step**2)
    return gradient, curvature


def _scalar_terms(
    values: np.ndarray,
    linear: np.ndarray,
    *,
    family: str,
    dispersion: float | None,
    zero_probability: float | None,
    offset: np.ndarray,
    trials: np.ndarray | None,
    censor_lower: np.ndarray | None,
    censor_upper: np.ndarray | None,
) -> tuple[float, np.ndarray, np.ndarray]:
    def contributions(candidate: np.ndarray) -> np.ndarray:
        return _scalar_log_likelihood_contributions(
            values,
            candidate,
            family=family,
            dispersion=dispersion,
            zero_probability=zero_probability,
            offset=offset,
            trials=trials,
            censor_lower=censor_lower,
            censor_upper=censor_upper,
        )

    current = contributions(linear)
    if family == "poisson":
        derivative = _poisson_scalar_terms(values, linear, offset)
    elif family == "negative-binomial":
        derivative = _negative_binomial_scalar_terms(values, linear, offset, dispersion)
    elif family in ZERO_COMPONENT_FAMILIES:
        derivative = _zero_component_count_scalar_terms(
            values,
            linear,
            offset,
            dispersion if "negative-binomial" in family else None,
            zero_probability,
            hurdle=family.startswith("hurdle"),
        )
    elif family == "gamma":
        derivative = _gamma_scalar_terms(values, linear, dispersion)
    elif family == "lognormal":
        derivative = _lognormal_scalar_terms(values, linear, dispersion)
    elif family == "beta":
        derivative = _beta_scalar_terms(values, linear, dispersion)
    elif family == "beta-binomial":
        derivative = _beta_binomial_scalar_terms(values, linear, dispersion, trials)
    else:
        derivative = _numerical_scalar_terms(contributions, current, linear)
    gradient, curvature = derivative
    return float(np.sum(current)), gradient, curvature


def _coerce_scalar_auxiliary(values, size, label, default=None):
    if values is None:
        return None if default is None else np.full(size, default, dtype=float)
    result = np.asarray(values, dtype=float)
    if result.shape != (size,):
        raise ValueError(
            "Response {} length differs from response length.".format(label)
        )
    return result


def _validate_count_response(response, family):
    count_family = family in {
        "poisson",
        "negative-binomial",
        "zero-inflated-poisson",
        "zero-inflated-negative-binomial",
        "hurdle-poisson",
        "hurdle-negative-binomial",
        "beta-binomial",
    }
    if count_family and (
        np.any(response < 0.0) or np.any(response != np.floor(response))
    ):
        raise ValueError("Count responses must be non-negative integers.")


def _validate_scalar_range(response, family):
    if family in {"gamma", "lognormal"} and np.any(response <= 0.0):
        raise ValueError("Positive response families require values greater than zero.")
    if family == "beta" and np.any((response <= 0.0) | (response >= 1.0)):
        raise ValueError("Beta responses must lie strictly between zero and one.")


def _validate_beta_binomial_response(response, family, trial_values):
    if family != "beta-binomial":
        return
    if trial_values is None:
        raise ValueError("Beta-binomial responses require --response-trials.")
    invalid = (
        not np.isfinite(trial_values).all()
        or np.any(trial_values <= 0.0)
        or np.any(trial_values != np.floor(trial_values))
        or np.any(response > trial_values)
    )
    if invalid:
        raise ValueError("Beta-binomial trials must be positive integers >= response.")


def _validate_censored_response(response, family, lower, upper):
    if family != "censored-gaussian":
        if not np.isfinite(response).all():
            raise ValueError("Non-Gaussian response values must be finite.")
        return
    size = len(response)
    if np.any(np.isinf(response)):
        raise ValueError("Censored-Gaussian response values must be finite or missing.")
    lower_values = np.full(size, np.nan) if lower is None else lower
    upper_values = np.full(size, np.nan) if upper is None else upper
    if np.any(np.isinf(lower_values)) or np.any(np.isinf(upper_values)):
        raise ValueError("Censor bounds must be finite or missing.")
    exact = np.isnan(lower_values) & np.isnan(upper_values)
    if np.any(exact & ~np.isfinite(response)):
        raise ValueError("Uncensored Gaussian observations must be finite.")
    if np.any(~exact & np.isfinite(response)):
        raise ValueError(
            "Censored Gaussian observations must have a missing response value."
        )
    invalid_interval = (
        np.isfinite(lower_values)
        & np.isfinite(upper_values)
        & (lower_values >= upper_values)
    )
    if np.any(invalid_interval):
        raise ValueError("Censor lower bounds must be smaller than upper bounds.")


def _validate_scalar_response(
    values: Sequence[object],
    family: str,
    *,
    offset: Sequence[float] | None,
    trials: Sequence[float] | None,
    censor_lower: Sequence[float] | None,
    censor_upper: Sequence[float] | None,
) -> tuple[
    np.ndarray, np.ndarray, np.ndarray | None, np.ndarray | None, np.ndarray | None
]:
    response = np.asarray(values, dtype=float)
    size = len(response)
    offsets = _coerce_scalar_auxiliary(offset, size, "offset", 0.0)
    assert offsets is not None
    trial_values = _coerce_scalar_auxiliary(trials, size, "trials")
    lower = _coerce_scalar_auxiliary(censor_lower, size, "censor lower bound")
    upper = _coerce_scalar_auxiliary(censor_upper, size, "censor upper bound")
    if not np.isfinite(offsets).all():
        raise ValueError("Response offsets must be finite.")
    _validate_censored_response(response, family, lower, upper)
    _validate_count_response(response, family)
    _validate_scalar_range(response, family)
    _validate_beta_binomial_response(response, family, trial_values)
    return response, offsets, trial_values, lower, upper


def _expand_replicated_scalar_inputs(
    response_values: Sequence[object],
    design: np.ndarray,
    offset: Sequence[float] | None,
    trials: Sequence[float] | None,
    censor_lower: Sequence[float] | None,
    censor_upper: Sequence[float] | None,
):
    counts = [
        len(value.values) if isinstance(value, ReplicatedObservation) else 1
        for value in response_values
    ]
    if any(count <= 0 for count in counts):
        raise ValueError("Replicated responses must contain at least one observation.")
    mapping = np.repeat(np.arange(len(response_values), dtype=int), counts)
    expanded_values = [
        observation
        for value in response_values
        for observation in (
            value.values if isinstance(value, ReplicatedObservation) else (value,)
        )
    ]
    expanded_design = design[mapping]

    def expand_auxiliary(values, default):
        if values is None:
            return (
                None if default is None else np.full(len(mapping), default, dtype=float)
            )
        if len(values) != len(response_values):
            raise ValueError("Response auxiliary length differs from tree-tip count.")
        expanded: list[object] = []
        for value, count in zip(values, counts, strict=True):
            if isinstance(value, ReplicatedObservation):
                if len(value.values) != count:
                    raise ValueError("Replicated response auxiliaries are misaligned.")
                expanded.extend(value.values)
            else:
                expanded.extend([value] * count)
        return np.asarray(expanded, dtype=float)

    return (
        expanded_values,
        expanded_design,
        mapping,
        expand_auxiliary(offset, 0.0),
        expand_auxiliary(trials, None),
        expand_auxiliary(censor_lower, None),
        expand_auxiliary(censor_upper, None),
    )


def _validate_glmm_inputs(
    response_values: Sequence[object], design: np.ndarray, family: str
) -> tuple[list[object], np.ndarray]:
    if family not in {"binomial", "multinomial", "ordinal"}:
        raise ValueError("Unsupported categorical response family: {}.".format(family))
    design = np.asarray(design, dtype=float)
    if design.ndim != 2 or not np.isfinite(design).all():
        raise ValueError("Categorical PGLMM design must be a finite matrix.")
    if family == "ordinal" and np.any(np.ptp(design, axis=0) <= 1e-14):
        raise ValueError(
            "Ordinal PGLMM design must not contain an intercept or constant column; "
            "the thresholds supply the location parameters."
        )
    values = list(response_values)
    if len(values) != len(design):
        raise ValueError("Categorical response and design lengths differ.")
    return values, design


def _resolve_glmm_levels(
    values: Sequence[object],
    family: str,
    levels: Sequence[str] | None,
    reference: str | None,
) -> tuple[list[str], str]:
    resolved = (
        observed_response_levels(values) if levels is None else [str(x) for x in levels]
    )
    observed = set(observed_response_levels(values))
    if (
        len(resolved) < 2
        or len(resolved) != len(set(resolved))
        or observed != set(resolved)
    ):
        raise ValueError(
            "Categorical response levels must be unique and cover every observed state."
        )
    if family == "binomial" and len(resolved) != 2:
        raise ValueError("Binomial responses must have exactly two levels.")
    if family == "multinomial" and len(resolved) < 3:
        raise ValueError("Multinomial responses must have at least three levels.")
    if family == "ordinal":
        return resolved, ""
    resolved_reference = str(reference or resolved[0])
    if resolved_reference not in resolved:
        raise ValueError("Response reference is not an observed level.")
    ordered = [level for level in resolved if level != resolved_reference]
    return ordered + [resolved_reference], resolved_reference


def _glmm_initial_parameters(
    counts: np.ndarray,
    design: np.ndarray,
    family: str,
    component_count: int,
    evolution_parameter_bounds: tuple[float, float] | None,
    evolution_parameter_initial: float | None,
) -> tuple[np.ndarray, list[tuple[float | None, float | None]]]:
    random_dimension = 1 if family == "ordinal" else counts.shape[1] - 1
    coefficient_count = design.shape[1] * random_dimension
    threshold_count = counts.shape[1] - 1 if family == "ordinal" else 0
    initial_thresholds = (
        _initial_threshold_parameters(counts)
        if family == "ordinal"
        else np.empty(0, dtype=float)
    )
    initial = np.concatenate(
        [
            _initial_coefficients(counts, design, family).reshape(-1),
            initial_thresholds,
            np.full(component_count, -1.0),
        ]
    )
    bounds: list[tuple[float | None, float | None]] = [(None, None)] * (
        coefficient_count + threshold_count
    ) + [(-12.0, 6.0)] * component_count
    if evolution_parameter_bounds is not None:
        lower, upper = evolution_parameter_bounds
        shape_initial = (
            (lower + upper) / 2.0
            if evolution_parameter_initial is None
            else float(evolution_parameter_initial)
        )
        initial = np.append(initial, shape_initial)
        bounds.append((float(lower), float(upper)))
    return initial, bounds


def _marginal_coefficient_covariance(
    design: np.ndarray,
    random_dimension: int,
    weights: np.ndarray,
    precision: np.ndarray,
) -> np.ndarray:
    expanded_design = np.kron(design, np.eye(random_dimension))
    weight_matrix = _dense_weight_matrix(weights)
    joint_fixed = expanded_design.T @ weight_matrix @ expanded_design
    fixed_random = expanded_design.T @ weight_matrix
    posterior_cholesky = _positive_definite_cholesky(weight_matrix + precision)
    adjustment = fixed_random @ np.linalg.solve(
        posterior_cholesky.T,
        np.linalg.solve(posterior_cholesky, fixed_random.T),
    )
    information = joint_fixed - adjustment
    information = (information + information.T) / 2.0
    return np.linalg.pinv(information, hermitian=True)


def _fixed_parameter_covariances(
    objective: Callable[[np.ndarray], float],
    optimized_parameters: np.ndarray,
    coefficient_count: int,
    threshold_count: int,
    coefficient_covariance: np.ndarray,
    *,
    numerical: bool,
) -> tuple[np.ndarray, np.ndarray, str]:
    if not numerical:
        return coefficient_covariance, np.empty((0, 0), dtype=float), "conditional"
    full_hessian = _finite_difference_hessian(objective, optimized_parameters)
    full_hessian = (full_hessian + full_hessian.T) / 2.0
    if not np.isfinite(full_hessian).all():
        return (
            np.full((coefficient_count, coefficient_count), np.nan),
            np.full((threshold_count, threshold_count), np.nan),
            "joint-information-non-finite",
        )
    fixed_parameter_count = coefficient_count + threshold_count
    fixed_information = full_hessian[:fixed_parameter_count, :fixed_parameter_count]
    if fixed_parameter_count < len(full_hessian):
        cross_information = full_hessian[:fixed_parameter_count, fixed_parameter_count:]
        nuisance_information = full_hessian[
            fixed_parameter_count:, fixed_parameter_count:
        ]
        fixed_information = (
            fixed_information
            - cross_information
            @ np.linalg.pinv(nuisance_information, hermitian=True)
            @ cross_information.T
        )
    fixed_information = (fixed_information + fixed_information.T) / 2.0
    eigenvalues = np.linalg.eigvalsh(fixed_information)
    tolerance = (
        np.finfo(float).eps
        * max(1.0, float(np.max(np.abs(eigenvalues))))
        * max(1, fixed_parameter_count)
        * 1000.0
    )
    if float(np.min(eigenvalues)) <= tolerance:
        return (
            np.full((coefficient_count, coefficient_count), np.nan),
            np.full((threshold_count, threshold_count), np.nan),
            "nuisance-information-singular",
        )
    fixed_covariance = np.linalg.inv(fixed_information)
    if not np.isfinite(fixed_covariance).all():
        return (
            np.full((coefficient_count, coefficient_count), np.nan),
            np.full((threshold_count, threshold_count), np.nan),
            "joint-covariance-non-finite",
        )
    coefficient_covariance = fixed_covariance[:coefficient_count, :coefficient_count]
    if not threshold_count:
        return coefficient_covariance, np.empty((0, 0), dtype=float), "ok"
    threshold_parameter_covariance = fixed_covariance[
        coefficient_count:, coefficient_count:
    ]
    threshold_parameters = optimized_parameters[coefficient_count:fixed_parameter_count]
    jacobian = np.zeros((threshold_count, threshold_count), dtype=float)
    jacobian[:, 0] = 1.0
    for row in range(1, threshold_count):
        jacobian[row, 1 : row + 1] = np.exp(threshold_parameters[1 : row + 1])
    return (
        coefficient_covariance,
        jacobian @ threshold_parameter_covariance @ jacobian.T,
        "ok",
    )


def _glmm_component_modes(
    components: Mapping[str, float],
    phylogenetic_base: np.ndarray,
    group_base: np.ndarray | None,
    precision: np.ndarray,
    mode: np.ndarray,
    random_dimension: int,
) -> dict[str, np.ndarray]:
    latent_score = precision @ mode
    n_observations = len(phylogenetic_base)

    def component_mode(name: str, covariance: np.ndarray) -> np.ndarray:
        expanded = np.kron(components[name] * covariance, np.eye(random_dimension))
        return (expanded @ latent_score).reshape(n_observations, random_dimension)

    modes = {"phylogenetic": component_mode("phylogenetic", phylogenetic_base)}
    if group_base is not None:
        modes["group"] = component_mode("group", group_base)
    return modes


def _glmm_boundary_warning(
    log_variances: np.ndarray,
    optimized_parameters: np.ndarray,
    evolution_parameter_bounds: tuple[float, float] | None,
) -> bool:
    variance_boundary = bool(
        np.any(log_variances < -11.99) or np.any(log_variances > 5.99)
    )
    if evolution_parameter_bounds is None:
        return variance_boundary
    lower, upper = evolution_parameter_bounds
    width = upper - lower
    shape_value = float(optimized_parameters[-1])
    shape_boundary = (
        shape_value <= lower + width * 1e-5 or shape_value >= upper - width * 1e-5
    )
    return variance_boundary or shape_boundary


def _binary_separation_detected(counts: np.ndarray, design: np.ndarray) -> bool:
    if counts.shape[1] != 2:
        return False
    first_only = (counts[:, 0] > 0.0) & (counts[:, 1] == 0.0)
    second_only = (counts[:, 1] > 0.0) & (counts[:, 0] == 0.0)
    informative = first_only | second_only
    if not np.any(first_only) or not np.any(second_only) or not np.all(informative):
        return False
    signs = np.where(first_only, 1.0, -1.0)
    constraints = -(signs[:, None] * design)
    result = linprog(
        np.zeros(design.shape[1], dtype=float),
        A_ub=constraints,
        b_ub=-np.ones(len(design), dtype=float),
        bounds=[(None, None)] * design.shape[1],
        method="highs",
    )
    return bool(result.success)


def _optimize_with_fixed_coefficient(
    objective: Callable[[np.ndarray], float],
    optimized: np.ndarray,
    bounds: Sequence[tuple[float | None, float | None]],
    coefficient_index: int,
    coefficient_value: float,
) -> float:
    free_indices = [
        index for index in range(len(optimized)) if index != coefficient_index
    ]

    def reduced_objective(free: np.ndarray) -> float:
        parameters = optimized.copy()
        parameters[coefficient_index] = coefficient_value
        parameters[free_indices] = free
        return objective(parameters)

    result = _robust_minimize(
        reduced_objective,
        optimized[free_indices],
        bounds=[bounds[index] for index in free_indices],
        label="Fixed-coefficient profile",
        maxiter=500,
    )
    return float(result.fun)


def _profile_endpoint(
    objective: Callable[[np.ndarray], float],
    optimized: np.ndarray,
    bounds: Sequence[tuple[float | None, float | None]],
    coefficient_index: int,
    standard_error: float,
    target: float,
    direction: float,
) -> float:
    estimate = float(optimized[coefficient_index])

    def difference(value: float) -> float:
        return (
            _optimize_with_fixed_coefficient(
                objective, optimized, bounds, coefficient_index, value
            )
            - target
        )

    inner = estimate
    outer = estimate + direction * max(standard_error, 0.1)
    for _attempt in range(12):
        value = difference(outer)
        if np.isfinite(value) and value >= 0.0:
            lower, upper = sorted((inner, outer))
            return float(brentq(difference, lower, upper, maxiter=60))
        outer = estimate + 2.0 * (outer - estimate)
    return float("nan")


def _coefficient_likelihood_inference(
    objective: Callable[[np.ndarray], float],
    optimized: np.ndarray,
    bounds: Sequence[tuple[float | None, float | None]],
    coefficient_count: int,
    coefficient_covariance: np.ndarray,
    *,
    inference: str,
    confidence_level: float,
):
    if inference == "wald":
        return None, None, None, None
    fitted_objective = objective(optimized)
    statistics = np.empty(coefficient_count, dtype=float)
    p_values = np.empty(coefficient_count, dtype=float)
    lower = np.full(coefficient_count, np.nan, dtype=float)
    upper = np.full(coefficient_count, np.nan, dtype=float)
    profile_target = fitted_objective + 0.5 * float(chi2.ppf(confidence_level, 1))
    for index in range(coefficient_count):
        null_objective = _optimize_with_fixed_coefficient(
            objective, optimized, bounds, index, 0.0
        )
        statistics[index] = max(0.0, 2.0 * (null_objective - fitted_objective))
        p_values[index] = float(chi2.sf(statistics[index], 1))
        if inference == "profile-likelihood":
            standard_error = math.sqrt(
                max(float(coefficient_covariance[index, index]), 0.0)
            )
            lower[index] = _profile_endpoint(
                objective,
                optimized,
                bounds,
                index,
                standard_error,
                profile_target,
                -1.0,
            )
            upper[index] = _profile_endpoint(
                objective,
                optimized,
                bounds,
                index,
                standard_error,
                profile_target,
                1.0,
            )
    return statistics, p_values, lower, upper


def _fitted_random_covariance(
    fit: PhylogeneticGlmmFit,
    design: np.ndarray,
    phylogenetic_covariance: Callable[[float | None], np.ndarray],
    group_covariance: np.ndarray | None,
    predictor_uncertainties: Sequence[np.ndarray],
    predictor_columns: Sequence[int | tuple[int, ...]],
) -> np.ndarray:
    """Reconstruct the fitted latent covariance used for simulation."""
    n_tips = len(design)
    phylogenetic_base = _normalize_covariance(
        phylogenetic_covariance(fit.evolution_parameter), n_tips, "Phylogenetic"
    )
    random_dimension = fit.coefficients.shape[1]
    tip_covariance = fit.component_variances["phylogenetic"] * phylogenetic_base
    if group_covariance is not None:
        group_value = (
            group_covariance.materialize()
            if isinstance(group_covariance, SparseCovarianceModel)
            else group_covariance
        )
        group_base = _normalize_covariance(
            group_value,
            n_tips,
            "Grouping random-effect",
            allow_semidefinite=True,
        )
        tip_covariance += fit.component_variances["group"] * group_base
    return np.kron(tip_covariance, np.eye(random_dimension)) + (
        _project_predictor_uncertainty(
            fit.coefficients,
            predictor_uncertainties,
            predictor_columns,
            n_tips,
        )
    )


def _draw_categorical_responses(
    rng: np.random.Generator,
    response_values: Sequence[object],
    design: np.ndarray,
    fit: PhylogeneticGlmmFit,
    latent: np.ndarray,
) -> list[CategoricalObservation]:
    counts = response_count_matrix(response_values, fit.levels)
    totals = np.rint(counts.sum(axis=1)).astype(int)
    fixed = design @ fit.coefficients
    if fit.family == "ordinal":
        linear = fixed[:, 0] + latent.reshape(len(design))
        lower = np.column_stack(
            [np.full(len(design), -np.inf), fit.thresholds[None, :] - linear[:, None]]
        )
        upper = np.column_stack(
            [fit.thresholds[None, :] - linear[:, None], np.full(len(design), np.inf)]
        )
        probabilities = expit(upper) - expit(lower)
    else:
        linear = fixed + latent.reshape(fixed.shape)
        logits = np.column_stack([linear, np.zeros(len(design), dtype=float)])
        probabilities = np.exp(logits - logsumexp(logits, axis=1)[:, None])
    simulated = []
    for total, probability in zip(totals, probabilities, strict=True):
        probability = np.maximum(probability, 0.0)
        probability /= probability.sum()
        sampled = rng.multinomial(int(total), probability)
        simulated.append(
            CategoricalObservation(
                {
                    level: float(count) / int(total)
                    for level, count in zip(fit.levels, sampled, strict=True)
                },
                int(total),
            )
        )
    return simulated


def _draw_count(
    rng: np.random.Generator,
    means: np.ndarray,
    family: str,
    dispersion: float | None,
) -> np.ndarray:
    if "negative-binomial" in family:
        assert dispersion is not None
        size = 1.0 / dispersion
        probability = size / (size + means)
        return rng.negative_binomial(size, probability).astype(float)
    return rng.poisson(means).astype(float)


def _draw_zero_truncated_count(
    rng: np.random.Generator,
    means: np.ndarray,
    family: str,
    dispersion: float | None,
) -> np.ndarray:
    if "negative-binomial" in family:
        assert dispersion is not None
        size = 1.0 / dispersion
        probability = size / (size + means)
        zero_cdf = nbinom_distribution.cdf(0, size, probability)
        quantiles = zero_cdf + rng.random(len(means)) * (1.0 - zero_cdf)
        return nbinom_distribution.ppf(
            np.minimum(quantiles, np.nextafter(1.0, 0.0)), size, probability
        ).astype(float)
    zero_cdf = poisson_distribution.cdf(0, means)
    quantiles = zero_cdf + rng.random(len(means)) * (1.0 - zero_cdf)
    return poisson_distribution.ppf(
        np.minimum(quantiles, np.nextafter(1.0, 0.0)), means
    ).astype(float)


def _draw_mixture_count(rng, linear, family, dispersion, zero_probability, offset):
    means = np.exp(np.clip(linear + offset, -30.0, 30.0))
    if family.startswith("hurdle"):
        assert zero_probability is not None
        values = _draw_zero_truncated_count(rng, means, family, dispersion)
        values[rng.random(len(values)) < zero_probability] = 0.0
        return values
    values = _draw_count(rng, means, family, dispersion)
    if family.startswith("zero-inflated"):
        assert zero_probability is not None
        values[rng.random(len(values)) < zero_probability] = 0.0
    return values


def _draw_gamma(rng, linear, dispersion):
    assert dispersion is not None
    means = np.exp(np.clip(linear, -30.0, 30.0))
    return rng.gamma(dispersion, means / dispersion)


def _draw_lognormal(rng, linear, dispersion):
    assert dispersion is not None
    return rng.lognormal(linear, dispersion)


def _draw_beta(rng, linear, dispersion):
    assert dispersion is not None
    means = np.clip(expit(linear), 1e-10, 1.0 - 1e-10)
    return rng.beta(means * dispersion, (1.0 - means) * dispersion)


def _draw_beta_binomial(rng, linear, dispersion, trials):
    assert dispersion is not None and trials is not None
    means = np.clip(expit(linear), 1e-10, 1.0 - 1e-10)
    probabilities = rng.beta(means * dispersion, (1.0 - means) * dispersion)
    return rng.binomial(trials.astype(int), probabilities).astype(float)


def _draw_censored_gaussian(rng, linear, dispersion, censor_lower, censor_upper):
    assert dispersion is not None
    values = rng.normal(linear, dispersion)
    lower = np.full(len(values), np.nan) if censor_lower is None else censor_lower
    upper = np.full(len(values), np.nan) if censor_upper is None else censor_upper
    values[~(np.isnan(lower) & np.isnan(upper))] = np.nan
    return values


def _draw_scalar_observations(
    rng: np.random.Generator,
    linear: np.ndarray,
    *,
    family: str,
    dispersion: float | None,
    zero_probability: float | None,
    offset: np.ndarray,
    trials: np.ndarray | None,
    censor_lower: np.ndarray | None,
    censor_upper: np.ndarray | None,
) -> np.ndarray:
    if family in {
        "poisson",
        "negative-binomial",
        "zero-inflated-poisson",
        "zero-inflated-negative-binomial",
        "hurdle-poisson",
        "hurdle-negative-binomial",
    }:
        return _draw_mixture_count(
            rng, linear, family, dispersion, zero_probability, offset
        )
    if family == "gamma":
        return _draw_gamma(rng, linear, dispersion)
    if family == "lognormal":
        return _draw_lognormal(rng, linear, dispersion)
    if family == "beta":
        return _draw_beta(rng, linear, dispersion)
    if family == "beta-binomial":
        return _draw_beta_binomial(rng, linear, dispersion, trials)
    if family == "censored-gaussian":
        return _draw_censored_gaussian(
            rng, linear, dispersion, censor_lower, censor_upper
        )
    raise ValueError("Unsupported scalar response family: {}.".format(family))


def _draw_scalar_responses(
    rng: np.random.Generator,
    response_values: Sequence[object],
    design: np.ndarray,
    fit: PhylogeneticGlmmFit,
    latent: np.ndarray,
    *,
    offset: Sequence[float] | None,
    trials: Sequence[float] | None,
    censor_lower: Sequence[float] | None,
    censor_upper: Sequence[float] | None,
) -> list[object]:
    (
        expanded_values,
        expanded_design,
        observation_to_tip,
        expanded_offset,
        expanded_trials,
        expanded_lower,
        expanded_upper,
    ) = _expand_replicated_scalar_inputs(
        response_values, design, offset, trials, censor_lower, censor_upper
    )
    _values, offsets, trial_values, lower, upper = _validate_scalar_response(
        expanded_values,
        fit.family,
        offset=expanded_offset,
        trials=expanded_trials,
        censor_lower=expanded_lower,
        censor_upper=expanded_upper,
    )
    linear = (expanded_design @ fit.coefficients)[:, 0] + latent.reshape(len(design))[
        observation_to_tip
    ]
    simulated = _draw_scalar_observations(
        rng,
        linear,
        family=fit.family,
        dispersion=fit.dispersion,
        zero_probability=fit.zero_probability,
        offset=offsets,
        trials=trial_values,
        censor_lower=lower,
        censor_upper=upper,
    )
    packed: list[object] = []
    position = 0
    for original in response_values:
        count = (
            len(original.values) if isinstance(original, ReplicatedObservation) else 1
        )
        selected = tuple(
            float(value) for value in simulated[position : position + count]
        )
        packed.append(
            ReplicatedObservation(selected)
            if isinstance(original, ReplicatedObservation)
            else selected[0]
        )
        position += count
    return packed


def _bootstrap_coefficient_inference(
    fit: PhylogeneticGlmmFit,
    coefficient_samples: np.ndarray,
    confidence_level: float,
) -> PhylogeneticGlmmFit:
    estimates = fit.coefficients.reshape(-1)
    centered = coefficient_samples - estimates[None, :]
    standard_errors = np.std(coefficient_samples, axis=0, ddof=1)
    statistics = np.divide(
        estimates,
        standard_errors,
        out=np.full_like(estimates, np.nan),
        where=standard_errors > 0.0,
    )
    exceedances = np.sum(np.abs(centered) >= np.abs(estimates)[None, :], axis=0)
    p_values = (exceedances + 1.0) / (len(coefficient_samples) + 1.0)
    alpha = 1.0 - confidence_level
    lower, upper = np.quantile(
        coefficient_samples, [alpha / 2.0, 1.0 - alpha / 2.0], axis=0
    )
    covariance = np.atleast_2d(np.cov(coefficient_samples, rowvar=False, ddof=1))
    return replace(
        fit,
        coefficient_covariance=covariance,
        coefficient_statistics=statistics,
        coefficient_p_values=p_values,
        coefficient_confidence_lower=lower,
        coefficient_confidence_upper=upper,
        coefficient_inference="parametric-bootstrap",
        coefficient_covariance_status="ok",
    )


def _scalar_random_mode(
    values: np.ndarray,
    fixed_linear: np.ndarray,
    precision: np.ndarray,
    observation_to_tip: np.ndarray,
    *,
    family: str,
    dispersion: float | None,
    zero_probability: float | None,
    offset: np.ndarray,
    trials: np.ndarray | None,
    censor_lower: np.ndarray | None,
    censor_upper: np.ndarray | None,
) -> tuple[np.ndarray, float, np.ndarray]:
    mode = np.zeros(len(precision), dtype=float)

    def terms(candidate: np.ndarray):
        return _scalar_terms(
            values,
            fixed_linear + candidate[observation_to_tip],
            family=family,
            dispersion=dispersion,
            zero_probability=zero_probability,
            offset=offset,
            trials=trials,
            censor_lower=censor_lower,
            censor_upper=censor_upper,
        )

    converged = False
    for _iteration in range(80):
        log_likelihood, likelihood_gradient, weights = terms(mode)
        mapped_gradient = np.zeros(len(mode), dtype=float)
        np.add.at(mapped_gradient, observation_to_tip, likelihood_gradient)
        mapped_weights = np.bincount(
            observation_to_tip,
            weights=weights,
            minlength=len(mode),
        )
        gradient = mapped_gradient + precision @ mode
        if not np.isfinite(gradient).all():
            break
        if float(np.max(np.abs(gradient))) < 1e-7:
            converged = True
            break
        hessian = np.diag(mapped_weights) + precision
        try:
            cholesky = _positive_definite_cholesky(hessian)
        except np.linalg.LinAlgError:
            break
        step = np.linalg.solve(cholesky.T, np.linalg.solve(cholesky, gradient))
        current = -log_likelihood + 0.5 * float(mode @ precision @ mode)
        accepted = False
        step_scale = 1.0
        for _backtrack in range(24):
            candidate = mode - step_scale * step
            candidate_log_likelihood, _, _ = terms(candidate)
            candidate_objective = -candidate_log_likelihood + 0.5 * float(
                candidate @ precision @ candidate
            )
            if np.isfinite(candidate_objective) and candidate_objective <= current:
                mode = candidate
                accepted = True
                break
            step_scale *= 0.5
        if not accepted:
            break
        if float(np.max(np.abs(step_scale * step))) < 1e-8:
            converged = True
            break
    if converged:
        _log_likelihood, likelihood_gradient, _weights = terms(mode)
        mapped_gradient = np.zeros(len(mode), dtype=float)
        np.add.at(mapped_gradient, observation_to_tip, likelihood_gradient)
        final_gradient = mapped_gradient + precision @ mode
        converged = bool(
            np.isfinite(final_gradient).all()
            and float(np.max(np.abs(final_gradient))) < 1e-5
        )
    if not converged:

        def value_and_gradient(candidate: np.ndarray):
            log_likelihood, likelihood_gradient, _weights = terms(candidate)
            mapped_gradient = np.zeros(len(candidate), dtype=float)
            np.add.at(mapped_gradient, observation_to_tip, likelihood_gradient)
            value = -log_likelihood + 0.5 * float(candidate @ precision @ candidate)
            return value, mapped_gradient + precision @ candidate

        fallback = minimize(
            value_and_gradient,
            mode,
            method="L-BFGS-B",
            jac=True,
            options={"maxiter": 1000, "ftol": 1e-12, "gtol": 1e-7},
        )
        if (
            not fallback.success
            or not np.isfinite(fallback.fun)
            or float(fallback.fun) >= 1e99
        ):
            raise RuntimeError(
                "Scalar random-effect mode optimization failed: {}".format(
                    fallback.message
                )
            )
        mode = np.asarray(fallback.x, dtype=float)
    log_likelihood, _, weights = terms(mode)
    objective = -log_likelihood + 0.5 * float(mode @ precision @ mode)
    mapped_weights = np.bincount(
        observation_to_tip,
        weights=weights,
        minlength=len(mode),
    )
    _positive_definite_cholesky(np.diag(mapped_weights) + precision)
    if not np.isfinite(objective):
        raise RuntimeError("Scalar random-effect mode is non-finite.")
    return mode, objective, weights


def _marginal_scalar_coefficient_covariance(
    design: np.ndarray,
    observation_to_tip: np.ndarray,
    weights: np.ndarray,
    precision: np.ndarray,
) -> np.ndarray:
    diagonal_weights = weights
    weighted_design = diagonal_weights[:, None] * design
    joint_fixed = design.T @ weighted_design
    fixed_random = np.zeros((design.shape[1], len(precision)), dtype=float)
    for observation, tip_index in enumerate(observation_to_tip):
        fixed_random[:, tip_index] += weighted_design[observation]
    mapped_weights = np.bincount(
        observation_to_tip,
        weights=diagonal_weights,
        minlength=len(precision),
    )
    posterior_cholesky = _positive_definite_cholesky(
        np.diag(mapped_weights) + precision
    )
    adjustment = fixed_random @ np.linalg.solve(
        posterior_cholesky.T,
        np.linalg.solve(posterior_cholesky, fixed_random.T),
    )
    information = joint_fixed - adjustment
    information = (information + information.T) / 2.0
    return np.linalg.pinv(information, hermitian=True)


def _scalar_initial_coefficient(
    values: np.ndarray,
    design: np.ndarray,
    family: str,
    offset: np.ndarray,
    trials: np.ndarray | None,
) -> np.ndarray:
    coefficients = np.zeros((design.shape[1], 1), dtype=float)
    if design.shape[1] == 0 or not np.allclose(design[:, 0], 1.0):
        return coefficients
    if family in {
        "poisson",
        "negative-binomial",
        "zero-inflated-poisson",
        "zero-inflated-negative-binomial",
        "hurdle-poisson",
        "hurdle-negative-binomial",
    }:
        mean = float(np.mean(values / np.exp(offset)))
        coefficients[0, 0] = math.log(max(mean, 1e-3))
    elif family in {"gamma", "lognormal"}:
        coefficients[0, 0] = float(np.mean(np.log(values)))
    elif family == "beta":
        mean = float(np.clip(np.mean(values), 1e-4, 1.0 - 1e-4))
        coefficients[0, 0] = math.log(mean / (1.0 - mean))
    elif family == "beta-binomial":
        assert trials is not None
        mean = float(np.clip(np.sum(values) / np.sum(trials), 1e-4, 1.0 - 1e-4))
        coefficients[0, 0] = math.log(mean / (1.0 - mean))
    elif family == "censored-gaussian":
        finite = values[np.isfinite(values)]
        coefficients[0, 0] = 0.0 if not len(finite) else float(np.mean(finite))
    return coefficients


def _scalar_auxiliary_initial(values: np.ndarray, family: str) -> float:
    if family == "lognormal":
        return math.log(max(float(np.std(np.log(values))), 0.5))
    if family == "censored-gaussian":
        finite = values[np.isfinite(values)]
        return math.log(max(float(np.std(finite)) if len(finite) else 1.0, 0.5))
    if family in {"gamma", "beta", "beta-binomial"}:
        return math.log(10.0)
    return math.log(0.2)


def _fit_scalar_phylogenetic_glmm(
    response_values: Sequence[object],
    design: np.ndarray,
    phylogenetic_covariance: Callable[[float | None], np.ndarray],
    *,
    family: str,
    group_covariance: np.ndarray | None,
    evolution_parameter: float | None,
    evolution_parameter_bounds: tuple[float, float] | None,
    evolution_parameter_decoder: Callable[[float], float],
    evolution_parameter_initial: float | None,
    predictor_uncertainties: Sequence[np.ndarray],
    predictor_columns: Sequence[int | tuple[int, ...]],
    offset: Sequence[float] | None,
    trials: Sequence[float] | None,
    censor_lower: Sequence[float] | None,
    censor_upper: Sequence[float] | None,
    dispersion: float | None,
    zero_probability: float | None,
    coefficient_penalty: str,
    coefficient_prior_sd: float | None,
    inference: str,
    confidence_level: float,
    allow_large_dense: bool,
) -> PhylogeneticGlmmFit:
    tip_design = np.asarray(design, dtype=float)
    if tip_design.ndim != 2 or not np.isfinite(tip_design).all():
        raise ValueError("Phylogenetic GLMM design must be a finite matrix.")
    n_tips = len(response_values)
    sparse_builder = getattr(phylogenetic_covariance, "sparse_model", None)
    sparse_capable = (
        callable(sparse_builder)
        and (
            group_covariance is None
            or isinstance(group_covariance, SparseCovarianceModel)
        )
        and _sparse_predictor_uncertainties(predictor_uncertainties)
    )
    use_sparse = _select_sparse_glmm_backend(
        n_tips, 1, sparse_capable, allow_large_dense
    )
    if n_tips != len(tip_design):
        raise ValueError("Response and design lengths differ.")
    (
        expanded_values,
        design,
        observation_to_tip,
        expanded_offset,
        expanded_trials,
        expanded_lower,
        expanded_upper,
    ) = _expand_replicated_scalar_inputs(
        response_values,
        tip_design,
        offset,
        trials,
        censor_lower,
        censor_upper,
    )
    values, offsets, trial_values, lower, upper = _validate_scalar_response(
        expanded_values,
        family,
        offset=expanded_offset,
        trials=expanded_trials,
        censor_lower=expanded_lower,
        censor_upper=expanded_upper,
    )
    group_base = None
    component_names = ["phylogenetic"]
    if group_covariance is not None:
        if isinstance(group_covariance, SparseCovarianceModel):
            if not use_sparse:
                group_base = _normalize_covariance(
                    group_covariance.materialize(),
                    n_tips,
                    "Grouping random-effect",
                    allow_semidefinite=True,
                )
        else:
            group_base = _normalize_covariance(
                group_covariance,
                n_tips,
                "Grouping random-effect",
                allow_semidefinite=True,
            )
        component_names.append("group")
    coefficient_count = design.shape[1]
    estimate_dispersion = family in DISPERSION_FAMILIES and dispersion is None
    estimate_zero = family in ZERO_COMPONENT_FAMILIES and zero_probability is None
    estimate_shape = evolution_parameter_bounds is not None
    initial = _scalar_initial_coefficient(
        values, design, family, offsets, trial_values
    ).reshape(-1)
    bounds: list[tuple[float | None, float | None]] = [(None, None)] * coefficient_count
    if estimate_dispersion:
        initial = np.append(initial, _scalar_auxiliary_initial(values, family))
        bounds.append((-10.0, 10.0))
    if estimate_zero:
        observed_zero = float(np.mean(values == 0.0))
        observed_zero = float(np.clip(observed_zero, 0.01, 0.99))
        initial = np.append(initial, math.log(observed_zero / (1.0 - observed_zero)))
        bounds.append((-10.0, 10.0))
    initial = np.append(initial, np.full(len(component_names), -1.0))
    bounds.extend([(-12.0, 6.0)] * len(component_names))
    if estimate_shape:
        assert evolution_parameter_bounds is not None
        shape_initial = (
            sum(evolution_parameter_bounds) / 2.0
            if evolution_parameter_initial is None
            else float(evolution_parameter_initial)
        )
        initial = np.append(initial, shape_initial)
        bounds.append(
            (
                float(evolution_parameter_bounds[0]),
                float(evolution_parameter_bounds[1]),
            )
        )

    def unpack(parameters: np.ndarray):
        position = coefficient_count
        coefficients = parameters[:position].reshape(design.shape[1], 1)
        fitted_dispersion = dispersion
        if estimate_dispersion:
            fitted_dispersion = math.exp(float(parameters[position]))
            position += 1
        fitted_zero = zero_probability
        if estimate_zero:
            fitted_zero = float(expit(parameters[position]))
            position += 1
        log_variances = parameters[position : position + len(component_names)]
        position += len(component_names)
        decoded_parameter = (
            evolution_parameter_decoder(float(parameters[position]))
            if estimate_shape
            else evolution_parameter
        )
        return (
            coefficients,
            fitted_dispersion,
            fitted_zero,
            log_variances,
            decoded_parameter,
        )

    def state(parameters: np.ndarray):
        coefficients, fitted_dispersion, fitted_zero, log_variances, decoded = unpack(
            parameters
        )
        components = dict(zip(component_names, np.exp(log_variances), strict=True))
        if use_sparse:
            latent = _sparse_glmm_latent_model(
                phylogenetic_covariance,
                decoded,
                components,
                group_covariance,
                1,
                predictor_uncertainties,
                predictor_columns,
                coefficients,
            )
            if latent is None:
                raise ValueError("Sparse phylogenetic covariance is unavailable.")
            observation_loading = latent.loading[observation_to_tip]

            def likelihood_terms(linear):
                return _scalar_terms(
                    values,
                    linear,
                    family=family,
                    dispersion=fitted_dispersion,
                    zero_probability=fitted_zero,
                    offset=offsets,
                    trials=trial_values,
                    censor_lower=lower,
                    censor_upper=upper,
                )

            mode, posterior_objective, weights, posterior_factor = _sparse_latent_mode(
                (design @ coefficients)[:, 0],
                observation_loading,
                latent.precision,
                likelihood_terms,
                label="Scalar",
            )
            objective_value = (
                posterior_objective
                + 0.5 * latent.prior_logdet
                + 0.5 * posterior_factor.logdet
                + _coefficient_penalty_value(
                    coefficients, coefficient_penalty, coefficient_prior_sd
                )
            )
            return (
                objective_value,
                mode,
                weights,
                latent,
                observation_loading,
                posterior_factor,
            )
        phylo_base = _normalize_covariance(
            phylogenetic_covariance(decoded), n_tips, "Phylogenetic"
        )
        random_covariance = math.exp(float(log_variances[0])) * phylo_base
        if group_base is not None:
            random_covariance += math.exp(float(log_variances[1])) * group_base
        expanded_covariance = random_covariance + _project_predictor_uncertainty(
            coefficients,
            predictor_uncertainties,
            predictor_columns,
            n_tips,
        )
        covariance_cholesky = _positive_definite_cholesky(expanded_covariance)
        precision = np.linalg.solve(
            covariance_cholesky.T,
            np.linalg.solve(covariance_cholesky, np.eye(n_tips)),
        )
        mode, posterior_objective, weights = _scalar_random_mode(
            values,
            (design @ coefficients)[:, 0],
            precision,
            observation_to_tip,
            family=family,
            dispersion=fitted_dispersion,
            zero_probability=fitted_zero,
            offset=offsets,
            trials=trial_values,
            censor_lower=lower,
            censor_upper=upper,
        )
        mapped_weights = np.bincount(
            observation_to_tip,
            weights=weights,
            minlength=n_tips,
        )
        hessian_cholesky = _positive_definite_cholesky(
            np.diag(mapped_weights) + precision
        )
        logdet_random = 2.0 * float(np.sum(np.log(np.diag(covariance_cholesky))))
        logdet_hessian = 2.0 * float(np.sum(np.log(np.diag(hessian_cholesky))))
        objective = (
            posterior_objective
            + 0.5 * logdet_random
            + 0.5 * logdet_hessian
            + _coefficient_penalty_value(
                coefficients, coefficient_penalty, coefficient_prior_sd
            )
        )
        return objective, mode, weights, precision, phylo_base, None

    def objective(parameters: np.ndarray) -> float:
        try:
            value = float(state(parameters)[0])
        except MemoryError as exc:
            raise _dense_glmm_memory_error(n_tips, 1) from exc
        except (
            ValueError,
            RuntimeError,
            np.linalg.LinAlgError,
            FloatingPointError,
            OverflowError,
        ):
            return 1e100
        return value if np.isfinite(value) else 1e100

    result = _robust_minimize(
        objective,
        initial,
        bounds=bounds,
        label="Scalar phylogenetic GLMM",
        force_multistart=family
        in ZERO_COMPONENT_FAMILIES | {"beta-binomial", "censored-gaussian"},
    )
    coefficients, fitted_dispersion, fitted_zero, log_variances, decoded = unpack(
        result.x
    )
    fitted_state = state(result.x)
    fitted_objective, mode, weights = fitted_state[:3]
    if use_sparse:
        latent, observation_loading, posterior_factor = fitted_state[3:]
        approximate_covariance = _sparse_marginal_coefficient_covariance(
            design, weights, observation_loading, posterior_factor
        )
    else:
        precision, _phylo_base, _unused = fitted_state[3:]
        approximate_covariance = _marginal_scalar_coefficient_covariance(
            design, observation_to_tip, weights, precision
        )
    coefficient_covariance, _unused, covariance_status = _fixed_parameter_covariances(
        objective,
        result.x,
        coefficient_count,
        0,
        approximate_covariance,
        numerical=True,
    )
    components = dict(zip(component_names, np.exp(log_variances), strict=True))
    if use_sparse:
        tip_mode = np.asarray(latent.loading @ mode).reshape(n_tips, 1)
        component_modes = _sparse_component_modes(
            latent, mode, n_tips, 1, component_names
        )
    else:
        tip_mode = mode.reshape(n_tips, 1)
        component_modes = _glmm_component_modes(
            components, _phylo_base, group_base, precision, mode, 1
        )
    condition = (
        float(np.linalg.cond(coefficient_covariance))
        if np.isfinite(coefficient_covariance).all()
        else float("inf")
    )
    separation = bool(
        np.max(np.abs(coefficients)) > 10.0
        or not np.isfinite(condition)
        or condition > 1e10
    )
    statistics, p_values, lower_limits, upper_limits = (
        _coefficient_likelihood_inference(
            objective,
            result.x,
            bounds,
            coefficient_count,
            coefficient_covariance,
            inference=inference,
            confidence_level=confidence_level,
        )
    )
    return PhylogeneticGlmmFit(
        family=family,
        levels=(),
        reference="",
        coefficients=coefficients,
        coefficient_covariance=coefficient_covariance,
        thresholds=np.empty(0, dtype=float),
        threshold_covariance=np.empty((0, 0), dtype=float),
        random_modes=tip_mode,
        component_random_modes=component_modes,
        component_variances=components,
        evolution_parameter=decoded,
        evolution_parameter_status="estimated" if estimate_shape else "fixed",
        log_likelihood=_reported_log_likelihood(
            fitted_objective,
            coefficients,
            coefficient_penalty,
            coefficient_prior_sd,
        ),
        optimizer_converged=True,
        optimizer_message=str(result.message),
        boundary_warning=_glmm_boundary_warning(
            log_variances, result.x, evolution_parameter_bounds
        ),
        dispersion=fitted_dispersion,
        zero_probability=fitted_zero,
        coefficient_penalty=coefficient_penalty,
        coefficient_prior_sd=coefficient_prior_sd,
        separation_warning=separation,
        coefficient_statistics=statistics,
        coefficient_p_values=p_values,
        coefficient_confidence_lower=lower_limits,
        coefficient_confidence_upper=upper_limits,
        coefficient_inference=inference,
        coefficient_covariance_status=covariance_status,
    )


@dataclass(frozen=True)
class _GlmmCallOptions:
    family: str
    levels: Sequence[str] | None
    reference: str | None
    group_covariance: np.ndarray | None
    evolution_parameter: float | None
    evolution_parameter_bounds: tuple[float, float] | None
    evolution_parameter_decoder: Callable[[float], float]
    evolution_parameter_initial: float | None
    predictor_uncertainties: Sequence[np.ndarray]
    predictor_columns: Sequence[int | tuple[int, ...]]
    offset: Sequence[float] | None
    trials: Sequence[float] | None
    censor_lower: Sequence[float] | None
    censor_upper: Sequence[float] | None
    dispersion: float | None
    zero_probability: float | None
    coefficient_penalty: str
    coefficient_prior_sd: float | None
    inference: str
    confidence_level: float
    bootstrap_replicates: int
    seed: int
    allow_large_dense: bool


def _validate_fixed_dispersion(value: float | None) -> None:
    if value is not None and (not np.isfinite(value) or value <= 0.0):
        raise ValueError("A fixed response dispersion must be positive and finite.")


def _validate_fixed_zero_probability(value: float | None) -> None:
    if value is not None and (not np.isfinite(value) or not 0.0 < value < 1.0):
        raise ValueError("A fixed zero probability must lie strictly in (0, 1).")


def _validate_coefficient_prior(penalty: str, prior_sd: float | None) -> None:
    if penalty != "none" and (
        prior_sd is None or not np.isfinite(prior_sd) or prior_sd <= 0.0
    ):
        raise ValueError("A positive finite coefficient prior SD is required.")


def _validate_glmm_call_options(options: _GlmmCallOptions) -> None:
    if not isinstance(options.allow_large_dense, (bool, np.bool_)):
        raise ValueError("allow_large_dense must be boolean.")
    count_families = {
        "poisson",
        "negative-binomial",
        "zero-inflated-poisson",
        "zero-inflated-negative-binomial",
        "hurdle-poisson",
        "hurdle-negative-binomial",
    }
    if options.offset is not None and options.family not in count_families:
        raise ValueError("Response offsets apply only to count families.")
    if options.trials is not None and options.family != "beta-binomial":
        raise ValueError("Response trials apply only to beta-binomial families.")
    if (
        options.censor_lower is not None or options.censor_upper is not None
    ) and options.family != "censored-gaussian":
        raise ValueError("Censor bounds apply only to the censored-gaussian family.")
    if options.dispersion is not None and options.family not in DISPERSION_FAMILIES:
        raise ValueError(
            "A fixed dispersion does not apply to family '{}'.".format(options.family)
        )
    if (
        options.zero_probability is not None
        and options.family not in ZERO_COMPONENT_FAMILIES
    ):
        raise ValueError(
            "A fixed zero probability does not apply to family '{}'.".format(
                options.family
            )
        )
    if options.coefficient_penalty not in {"none", "gaussian", "student-t"}:
        raise ValueError(
            "Unsupported coefficient penalty: {}.".format(options.coefficient_penalty)
        )
    _validate_fixed_dispersion(options.dispersion)
    _validate_fixed_zero_probability(options.zero_probability)
    _validate_coefficient_prior(
        options.coefficient_penalty, options.coefficient_prior_sd
    )
    if options.inference not in {
        "wald",
        "parametric-bootstrap",
        "likelihood-ratio",
        "profile-likelihood",
    }:
        raise ValueError(
            "Unsupported non-Gaussian inference: {}.".format(options.inference)
        )
    if not 0.0 < options.confidence_level < 1.0:
        raise ValueError("Confidence level must lie strictly between zero and one.")
    if (
        not isinstance(options.bootstrap_replicates, int)
        or isinstance(options.bootstrap_replicates, bool)
        or options.bootstrap_replicates < 2
    ):
        raise ValueError("bootstrap_replicates must be an integer of at least two.")
    if (
        not isinstance(options.seed, int)
        or isinstance(options.seed, bool)
        or options.seed < 0
    ):
        raise ValueError("seed must be a non-negative integer.")


def _call_phylogenetic_glmm(
    response_values,
    design,
    phylogenetic_covariance,
    options: _GlmmCallOptions,
) -> PhylogeneticGlmmFit:
    return fit_phylogenetic_glmm(
        response_values,
        design,
        phylogenetic_covariance,
        family=options.family,
        levels=options.levels,
        reference=options.reference,
        group_covariance=options.group_covariance,
        evolution_parameter=options.evolution_parameter,
        evolution_parameter_bounds=options.evolution_parameter_bounds,
        evolution_parameter_decoder=options.evolution_parameter_decoder,
        evolution_parameter_initial=options.evolution_parameter_initial,
        predictor_uncertainties=options.predictor_uncertainties,
        predictor_columns=options.predictor_columns,
        offset=options.offset,
        trials=options.trials,
        censor_lower=options.censor_lower,
        censor_upper=options.censor_upper,
        dispersion=options.dispersion,
        zero_probability=options.zero_probability,
        coefficient_penalty=options.coefficient_penalty,
        coefficient_prior_sd=options.coefficient_prior_sd,
        inference=options.inference,
        confidence_level=options.confidence_level,
        bootstrap_replicates=options.bootstrap_replicates,
        seed=options.seed,
        allow_large_dense=options.allow_large_dense,
    )


def _draw_bootstrap_responses(
    rng,
    response_values,
    design,
    fit,
    latent,
    options: _GlmmCallOptions,
) -> Sequence[object]:
    if options.family not in SCALAR_RESPONSE_FAMILIES:
        return _draw_categorical_responses(rng, response_values, design, fit, latent)
    return _draw_scalar_responses(
        rng,
        response_values,
        design,
        fit,
        latent,
        offset=options.offset,
        trials=options.trials,
        censor_lower=options.censor_lower,
        censor_upper=options.censor_upper,
    )


def _fit_parametric_bootstrap_glmm(
    response_values,
    design,
    phylogenetic_covariance,
    options: _GlmmCallOptions,
) -> PhylogeneticGlmmFit:
    wald_options = replace(options, inference="wald")
    fit = _call_phylogenetic_glmm(
        response_values, design, phylogenetic_covariance, wald_options
    )
    design_array = np.asarray(design, dtype=float)
    sparse_builder = getattr(phylogenetic_covariance, "sparse_model", None)
    sparse_capable = (
        callable(sparse_builder)
        and (
            options.group_covariance is None
            or isinstance(options.group_covariance, SparseCovarianceModel)
        )
        and _sparse_predictor_uncertainties(options.predictor_uncertainties)
    )
    sparse_sampler = None
    random_cholesky = None
    if sparse_capable:
        latent_model = _sparse_glmm_latent_model(
            phylogenetic_covariance,
            fit.evolution_parameter,
            fit.component_variances,
            options.group_covariance,
            fit.coefficients.shape[1],
            options.predictor_uncertainties,
            options.predictor_columns,
            fit.coefficients,
        )
        if latent_model is None:
            raise ValueError("Sparse phylogenetic covariance is unavailable.")
        sparse_sampler = prepare_sparse_latent_sampler(latent_model)
    else:
        try:
            random_covariance = _fitted_random_covariance(
                fit,
                design_array,
                phylogenetic_covariance,
                options.group_covariance,
                options.predictor_uncertainties,
                options.predictor_columns,
            )
            random_cholesky = _positive_definite_cholesky(random_covariance)
        except MemoryError as exc:
            raise _dense_glmm_memory_error(
                len(design_array), fit.coefficients.shape[1]
            ) from exc
    rng = np.random.default_rng(options.seed)
    samples: list[np.ndarray] = []
    maximum_attempts = max(
        options.bootstrap_replicates * 3, options.bootstrap_replicates + 10
    )
    attempts = 0
    refit_options = replace(
        wald_options,
        levels=fit.levels or options.levels,
        reference=fit.reference or options.reference,
    )
    while len(samples) < options.bootstrap_replicates and attempts < maximum_attempts:
        attempts += 1
        if sparse_sampler is not None:
            latent = sparse_sampler.sample(rng)
        else:
            assert random_cholesky is not None
            latent = random_cholesky @ rng.normal(size=len(random_cholesky))
        simulated = _draw_bootstrap_responses(
            rng, response_values, design_array, fit, latent, options
        )
        try:
            refit = _call_phylogenetic_glmm(
                simulated, design_array, phylogenetic_covariance, refit_options
            )
        except (ValueError, RuntimeError, np.linalg.LinAlgError):
            continue
        if not refit.optimizer_converged or not np.isfinite(refit.log_likelihood):
            continue
        coefficients = refit.coefficients.reshape(-1)
        if np.isfinite(coefficients).all():
            samples.append(coefficients)
    if len(samples) < options.bootstrap_replicates:
        raise RuntimeError(
            "Non-Gaussian parametric bootstrap produced only {} of {} "
            "successful refits after {} attempts.".format(
                len(samples), options.bootstrap_replicates, attempts
            )
        )
    return _bootstrap_coefficient_inference(
        fit, np.asarray(samples, dtype=float), options.confidence_level
    )


def fit_phylogenetic_glmm(
    response_values: Sequence[object],
    design: np.ndarray,
    phylogenetic_covariance: Callable[[float | None], np.ndarray],
    *,
    family: str,
    levels: Sequence[str] | None = None,
    reference: str | None = None,
    group_covariance: np.ndarray | None = None,
    evolution_parameter: float | None = None,
    evolution_parameter_bounds: tuple[float, float] | None = None,
    evolution_parameter_decoder: Callable[[float], float] = float,
    evolution_parameter_initial: float | None = None,
    predictor_uncertainties: Sequence[np.ndarray] = (),
    predictor_columns: Sequence[int | tuple[int, ...]] = (),
    offset: Sequence[float] | None = None,
    trials: Sequence[float] | None = None,
    censor_lower: Sequence[float] | None = None,
    censor_upper: Sequence[float] | None = None,
    dispersion: float | None = None,
    zero_probability: float | None = None,
    coefficient_penalty: str = "student-t",
    coefficient_prior_sd: float | None = 2.5,
    inference: str = "wald",
    confidence_level: float = 0.95,
    bootstrap_replicates: int = 1000,
    seed: int = 1,
    allow_large_dense: bool = False,
) -> PhylogeneticGlmmFit:
    """Fit a categorical, count, positive, or proportion phylogenetic GLMM."""
    options = _GlmmCallOptions(
        family=family,
        levels=levels,
        reference=reference,
        group_covariance=group_covariance,
        evolution_parameter=evolution_parameter,
        evolution_parameter_bounds=evolution_parameter_bounds,
        evolution_parameter_decoder=evolution_parameter_decoder,
        evolution_parameter_initial=evolution_parameter_initial,
        predictor_uncertainties=predictor_uncertainties,
        predictor_columns=predictor_columns,
        offset=offset,
        trials=trials,
        censor_lower=censor_lower,
        censor_upper=censor_upper,
        dispersion=dispersion,
        zero_probability=zero_probability,
        coefficient_penalty=coefficient_penalty,
        coefficient_prior_sd=coefficient_prior_sd,
        inference=inference,
        confidence_level=confidence_level,
        bootstrap_replicates=bootstrap_replicates,
        seed=seed,
        allow_large_dense=allow_large_dense,
    )
    _validate_glmm_call_options(options)
    if inference == "parametric-bootstrap":
        return _fit_parametric_bootstrap_glmm(
            response_values, design, phylogenetic_covariance, options
        )
    if family in SCALAR_RESPONSE_FAMILIES:
        return _fit_scalar_phylogenetic_glmm(
            response_values,
            design,
            phylogenetic_covariance,
            family=family,
            group_covariance=group_covariance,
            evolution_parameter=evolution_parameter,
            evolution_parameter_bounds=evolution_parameter_bounds,
            evolution_parameter_decoder=evolution_parameter_decoder,
            evolution_parameter_initial=evolution_parameter_initial,
            predictor_uncertainties=predictor_uncertainties,
            predictor_columns=predictor_columns,
            offset=offset,
            trials=trials,
            censor_lower=censor_lower,
            censor_upper=censor_upper,
            dispersion=dispersion,
            zero_probability=zero_probability,
            coefficient_penalty=coefficient_penalty,
            coefficient_prior_sd=coefficient_prior_sd,
            inference=inference,
            confidence_level=confidence_level,
            allow_large_dense=allow_large_dense,
        )
    values, design = _validate_glmm_inputs(response_values, design, family)
    sparse_builder = getattr(phylogenetic_covariance, "sparse_model", None)
    sparse_capable = (
        callable(sparse_builder)
        and (
            group_covariance is None
            or isinstance(group_covariance, SparseCovarianceModel)
        )
        and _sparse_predictor_uncertainties(predictor_uncertainties)
    )
    ordered_levels, resolved_reference = _resolve_glmm_levels(
        values, family, levels, reference
    )
    counts = response_count_matrix(values, ordered_levels)
    random_dimension = 1 if family == "ordinal" else len(ordered_levels) - 1
    use_sparse = _select_sparse_glmm_backend(
        len(design), random_dimension, sparse_capable, allow_large_dense
    )
    if (
        use_sparse
        and len(design) * random_dimension > MAX_SPARSE_GLMM_LINEAR_PREDICTORS
    ):
        warnings.warn(
            "Sparse categorical phylogenetic GLMM fitting above {} tip-level "
            "linear predictors is not yet included in the routine validation "
            "suite (received {}); attempting the fit.".format(
                MAX_SPARSE_GLMM_LINEAR_PREDICTORS,
                len(design) * random_dimension,
            ),
            RuntimeWarning,
            stacklevel=2,
        )
    coefficient_count = design.shape[1] * random_dimension
    threshold_count = len(ordered_levels) - 1 if family == "ordinal" else 0
    component_names = ["phylogenetic"]
    group_base = None
    if group_covariance is not None:
        if isinstance(group_covariance, SparseCovarianceModel):
            if not use_sparse:
                group_base = _normalize_covariance(
                    group_covariance.materialize(),
                    len(design),
                    "Grouping random-effect",
                    allow_semidefinite=True,
                )
        else:
            group_base = _normalize_covariance(
                group_covariance,
                len(design),
                "Grouping random-effect",
                allow_semidefinite=True,
            )
        component_names.append("group")
    component_count = len(component_names)
    estimate_shape = evolution_parameter_bounds is not None
    initial, bounds = _glmm_initial_parameters(
        counts,
        design,
        family,
        component_count,
        evolution_parameter_bounds,
        evolution_parameter_initial,
    )

    def unpack(parameters: np.ndarray):
        coefficients = parameters[:coefficient_count].reshape(
            design.shape[1], random_dimension
        )
        threshold_parameters = parameters[
            coefficient_count : coefficient_count + threshold_count
        ]
        variance_start = coefficient_count + threshold_count
        log_variances = parameters[variance_start : variance_start + component_count]
        decoded_parameter: float | None
        if estimate_shape:
            decoded_parameter = evolution_parameter_decoder(
                float(parameters[variance_start + component_count])
            )
        else:
            decoded_parameter = evolution_parameter
        thresholds = (
            _ordinal_thresholds(threshold_parameters)
            if family == "ordinal"
            else np.empty(0, dtype=float)
        )
        return coefficients, thresholds, log_variances, decoded_parameter

    def state(parameters: np.ndarray):
        coefficients, thresholds, log_variances, decoded_parameter = unpack(parameters)
        components = dict(zip(component_names, np.exp(log_variances), strict=True))
        if use_sparse:
            latent = _sparse_glmm_latent_model(
                phylogenetic_covariance,
                decoded_parameter,
                components,
                group_covariance,
                random_dimension,
                predictor_uncertainties,
                predictor_columns,
                coefficients,
            )
            if latent is None:
                raise ValueError("Sparse phylogenetic covariance is unavailable.")
            fixed_linear = design @ coefficients
            expanded_design = (
                design
                if family == "ordinal"
                else np.kron(design, np.eye(random_dimension))
            )

            def likelihood_terms(linear):
                if family == "ordinal":
                    return _ordinal_terms(counts, linear, thresholds)
                return _multinomial_terms(
                    counts, linear.reshape(len(design), random_dimension)
                )

            mode, posterior_objective, weights, posterior_factor = _sparse_latent_mode(
                fixed_linear.reshape(-1),
                latent.loading,
                latent.precision,
                likelihood_terms,
                label="Categorical",
            )
            objective_value = (
                posterior_objective
                + 0.5 * latent.prior_logdet
                + 0.5 * posterior_factor.logdet
                + _coefficient_penalty_value(
                    coefficients, coefficient_penalty, coefficient_prior_sd
                )
            )
            return (
                objective_value,
                mode,
                weights,
                latent,
                expanded_design,
                posterior_factor,
            )
        phylo_base = _normalize_covariance(
            phylogenetic_covariance(decoded_parameter),
            len(design),
            "Phylogenetic",
        )
        random_covariance = np.exp(log_variances[0]) * phylo_base
        if group_base is not None:
            random_covariance = (
                random_covariance + np.exp(log_variances[1]) * group_base
            )
        expanded_covariance = np.kron(
            random_covariance, np.eye(random_dimension)
        ) + _project_predictor_uncertainty(
            coefficients,
            predictor_uncertainties,
            predictor_columns,
            len(design),
        )
        covariance_cholesky = _positive_definite_cholesky(expanded_covariance)
        precision = np.linalg.solve(
            covariance_cholesky.T,
            np.linalg.solve(
                covariance_cholesky, np.eye(len(design) * random_dimension)
            ),
        )
        fixed_linear = design @ coefficients
        if family == "ordinal":
            fixed_linear = fixed_linear[:, 0]
        mode, posterior_objective, weights = _random_mode(
            counts,
            fixed_linear,
            precision,
            family=family,
            thresholds=thresholds,
        )
        posterior_hessian = _dense_weight_matrix(weights) + precision
        hessian_cholesky = _positive_definite_cholesky(posterior_hessian)
        logdet_random = 2.0 * float(np.sum(np.log(np.diag(covariance_cholesky))))
        logdet_hessian = 2.0 * float(np.sum(np.log(np.diag(hessian_cholesky))))
        objective = (
            posterior_objective
            + 0.5 * logdet_random
            + 0.5 * logdet_hessian
            + _coefficient_penalty_value(
                coefficients, coefficient_penalty, coefficient_prior_sd
            )
        )
        return objective, mode, weights, precision, random_covariance, None

    def objective(parameters: np.ndarray) -> float:
        try:
            value = float(state(parameters)[0])
        except MemoryError as exc:
            raise _dense_glmm_memory_error(len(design), random_dimension) from exc
        except (
            ValueError,
            RuntimeError,
            np.linalg.LinAlgError,
            FloatingPointError,
            OverflowError,
        ):
            return 1e100
        return value if np.isfinite(value) else 1e100

    result = _robust_minimize(
        objective,
        initial,
        bounds=bounds,
        label="Categorical phylogenetic GLMM",
    )
    coefficients, thresholds, log_variances, decoded_parameter = unpack(result.x)
    fitted_state = state(result.x)
    fitted_objective, mode, weights = fitted_state[:3]
    if use_sparse:
        latent, expanded_design, posterior_factor = fitted_state[3:]
        coefficient_covariance = _sparse_marginal_coefficient_covariance(
            expanded_design, weights, latent.loading, posterior_factor
        )
    else:
        precision, _random_covariance, _unused = fitted_state[3:]
        coefficient_covariance = _marginal_coefficient_covariance(
            design, random_dimension, weights, precision
        )
    (
        coefficient_covariance,
        threshold_covariance,
        covariance_status,
    ) = _fixed_parameter_covariances(
        objective,
        result.x,
        coefficient_count,
        threshold_count,
        coefficient_covariance,
        numerical=True,
    )
    variance_values = np.exp(log_variances)
    components = dict(zip(component_names, variance_values, strict=True))
    if use_sparse:
        tip_mode = np.asarray(latent.loading @ mode).reshape(
            len(design), random_dimension
        )
        component_random_modes = _sparse_component_modes(
            latent, mode, len(design), random_dimension, component_names
        )
    else:
        tip_mode = mode.reshape(len(design), random_dimension)
        fitted_phylogenetic = _normalize_covariance(
            phylogenetic_covariance(decoded_parameter),
            len(design),
            "Phylogenetic",
        )
        component_random_modes = _glmm_component_modes(
            components,
            fitted_phylogenetic,
            group_base,
            precision,
            mode,
            random_dimension,
        )
    condition = (
        float(np.linalg.cond(coefficient_covariance))
        if np.isfinite(coefficient_covariance).all()
        else float("inf")
    )
    separation = bool(
        np.max(np.abs(coefficients)) > 10.0
        or not np.isfinite(condition)
        or condition > 1e10
        or (family == "binomial" and _binary_separation_detected(counts, design))
    )
    statistics, p_values, lower_limits, upper_limits = (
        _coefficient_likelihood_inference(
            objective,
            result.x,
            bounds,
            coefficient_count,
            coefficient_covariance,
            inference=inference,
            confidence_level=confidence_level,
        )
    )
    return PhylogeneticGlmmFit(
        family=family,
        levels=tuple(ordered_levels),
        reference=resolved_reference,
        coefficients=coefficients,
        coefficient_covariance=coefficient_covariance,
        thresholds=thresholds,
        threshold_covariance=threshold_covariance,
        random_modes=tip_mode,
        component_random_modes=component_random_modes,
        component_variances=components,
        evolution_parameter=decoded_parameter,
        evolution_parameter_status="estimated" if estimate_shape else "fixed",
        log_likelihood=_reported_log_likelihood(
            fitted_objective,
            coefficients,
            coefficient_penalty,
            coefficient_prior_sd,
        ),
        optimizer_converged=True,
        optimizer_message=str(result.message),
        boundary_warning=_glmm_boundary_warning(
            log_variances, result.x, evolution_parameter_bounds
        ),
        coefficient_penalty=coefficient_penalty,
        coefficient_prior_sd=coefficient_prior_sd,
        separation_warning=separation,
        coefficient_statistics=statistics,
        coefficient_p_values=p_values,
        coefficient_confidence_lower=lower_limits,
        coefficient_confidence_upper=upper_limits,
        coefficient_inference=inference,
        coefficient_covariance_status=covariance_status,
    )
