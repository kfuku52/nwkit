"""Laplace-approximated phylogenetic generalized linear mixed models."""

import math
import warnings
from dataclasses import dataclass, replace
from decimal import Decimal, localcontext
from functools import lru_cache
from typing import Callable, Mapping, Sequence, cast

import numpy as np
from scipy import sparse
from scipy.linalg import cho_solve
from scipy.optimize import brentq, linprog, minimize
from scipy.special import (
    betainc,
    betaln,
    digamma,
    erfcx,
    expit,
    gammaln,
    log_ndtr,
    logsumexp,
    pdtrc,
    polygamma,
)
from scipy.stats import chi2
from scipy.stats import nbinom as nbinom_distribution
from scipy.stats import poisson as poisson_distribution

from nwkit.evolution import EvolutionaryCovarianceFactory
from nwkit.gaussian import finite_mean, symmetric_average
from nwkit.measurement_error import _finite_difference_hessian
from nwkit.model_matrix import CategoricalObservation, ReplicatedObservation
from nwkit.sparse_laplace import (
    ContinuousPredictorUncertainty,
    GmrfPredictorUncertainty,
    GroupedPredictorUncertainty,
    JointPredictorUncertainty,
    SparseCovarianceModel,
    SparseLatentModel,
    SparsePositiveDefiniteFactor,
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

_SMALL_NB_DISPERSION = 1e-6
_POISSON_LIMIT_DISPERSION = (np.finfo(float).eps / 2.0) ** 2
# Coefficients of ((1+z)*log1p(z)/z - 1) for the Stirling log-gamma ratio.
_NB_LOG_GAMMA_RATIO_SERIES = (0.0,) + tuple(
    (-1.0) ** (order + 1) / (order * (order + 1)) for order in range(1, 9)
)
_SMALL_BETA_PRECISION = 1e-20
_LARGE_BETA_PRECISION = 1e6
_LARGE_GAMMA_SHAPE = 128.0
_LARGE_POISSON_COUNT = 128.0
_LARGE_NEGATIVE_BINOMIAL_COUNT = 1e6
_EXPONENTIAL_DEVIANCE_SERIES = tuple(
    1.0 / math.factorial(order) for order in range(2, 11)
)
_LOG1P_REMAINDER_SERIES = tuple((-1.0) ** (order + 1) / order for order in range(2, 15))


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
    symmetric = symmetric_average(matrix)
    try:
        return np.linalg.cholesky(symmetric)
    except np.linalg.LinAlgError:
        return _jittered_cholesky(symmetric)


def _jittered_cholesky(symmetric):
    diagonal = np.diag_indices_from(symmetric)
    original = symmetric[diagonal]
    scale = max(1.0, float(np.max(original)))
    for multiplier in (1e-12, 1e-10, 1e-8, 1e-6):
        symmetric[diagonal] = original + scale * multiplier
        try:
            return np.linalg.cholesky(symmetric)
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


def _scale_covariance_by_mean_variance(matrix, size, label):
    matrix = np.asarray(matrix, dtype=float)
    if not size or matrix.shape != (size, size) or not np.isfinite(matrix).all():
        raise ValueError(
            "{} covariance has invalid dimensions or values.".format(label)
        )
    matrix = symmetric_average(matrix)
    diagonal_mean = finite_mean(np.diag(matrix))
    if diagonal_mean < np.finfo(float).tiny:
        # Normalize a subnormal mean before rounding it to float precision.
        # Conversely, normal means can divide the original matrix directly,
        # preserving entries that would underflow in an intermediate rescaling.
        _, exponent = math.frexp(float(np.max(np.diag(matrix))))
        with np.errstate(over="ignore"):
            matrix = np.ldexp(matrix, -exponent)
        diagonal_mean = finite_mean(np.diag(matrix))
    if diagonal_mean <= 0.0:
        raise ValueError("{} covariance needs a positive diagonal.".format(label))
    with np.errstate(over="ignore"):
        normalized = matrix / diagonal_mean
    if not np.isfinite(normalized).all():
        raise ValueError(
            "{} covariance has non-finite normalized values.".format(label)
        )
    return normalized


def _normalize_covariance(
    matrix: np.ndarray, size: int, label: str, *, allow_semidefinite: bool = False
) -> np.ndarray:
    normalized = _scale_covariance_by_mean_variance(matrix, size, label)
    if allow_semidefinite:
        eigenvalues = np.linalg.eigvalsh(normalized)
        if float(np.min(eigenvalues)) < -1e-9 * max(1.0, float(np.max(eigenvalues))):
            raise ValueError(
                "{} covariance must be positive semidefinite.".format(label)
            )
    else:
        _positive_definite_cholesky(normalized)
    return normalized


def _reuse_tree_covariance(factory, covariance):
    if type(factory) is EvolutionaryCovarianceFactory:
        return lru_cache(maxsize=1)(covariance)
    return covariance


def _prepare_dense_phylogenetic_covariance(factory, n_tips):
    def covariance(parameter):
        return _normalize_covariance(factory(parameter), n_tips, "Phylogenetic")

    return _reuse_tree_covariance(factory, covariance)


def _multinomial_terms(
    counts: np.ndarray, linear: np.ndarray
) -> tuple[float, np.ndarray, np.ndarray]:
    log_probabilities = _multinomial_log_probabilities(linear)
    probabilities = np.exp(log_probabilities)
    totals = counts.sum(axis=1)
    log_likelihood = float(
        np.sum(
            np.multiply(
                counts,
                log_probabilities,
                out=np.zeros_like(log_probabilities),
                where=counts != 0.0,
            )
        )
    )
    nonreference = probabilities[:, :-1]
    # Obtain the complement without subtracting a probability rounded to one.
    remaining = _multinomial_competing_sum(probabilities)
    competing_counts = _multinomial_competing_sum(counts)
    gradient = competing_counts * nonreference
    gradient -= counts[:, :-1] * remaining
    weights = -nonreference[:, :, None] * nonreference[:, None, :]
    diagonal = np.arange(linear.shape[1])
    weights[:, diagonal, diagonal] = nonreference * remaining
    weights *= totals[:, None, None]
    return log_likelihood, gradient.reshape(-1), weights


def _multinomial_log_probabilities(linear):
    logits = np.column_stack([linear, np.zeros(len(linear), dtype=float)])
    with np.errstate(over="ignore"):
        logits -= np.max(logits, axis=1, keepdims=True)
    logits -= np.logaddexp.reduce(logits, axis=1, keepdims=True)
    return logits


def _multinomial_competing_sum(values):
    remaining = np.cumsum(values[:, :0:-1], axis=1)[:, ::-1].copy()
    remaining[:, 1:] += np.cumsum(values[:, :-2], axis=1)
    return remaining


def _ordinal_thresholds(parameters: np.ndarray) -> np.ndarray:
    thresholds = np.empty(len(parameters), dtype=float)
    thresholds[0] = parameters[0]
    with np.errstate(over="ignore"):
        thresholds[1:] = parameters[0] + np.cumsum(np.exp(parameters[1:]))
    if not np.isfinite(thresholds).all():
        raise ValueError("Ordinal response thresholds must be finite.")
    return thresholds


def _ordinal_log_probabilities(
    centered_thresholds: np.ndarray, interval_log_widths: np.ndarray
) -> np.ndarray:
    log_probabilities = np.zeros(
        (len(centered_thresholds), centered_thresholds.shape[1] + 1), dtype=float
    )
    # The logistic interval mass is sigmoid(upper) * sigmoid(-lower)
    # times (1 - exp(lower - upper)). Its width is independent of the predictor.
    log_probabilities[:, :-1] = -np.logaddexp(0.0, -centered_thresholds)
    log_probabilities[:, 1:] -= np.logaddexp(0.0, centered_thresholds)
    log_probabilities[:, 1:-1] += interval_log_widths
    return log_probabilities


def _ordinal_interval_log_widths(thresholds):
    with np.errstate(divide="ignore", over="ignore"):
        return np.log(-np.expm1(-np.diff(thresholds)))


def _ordinal_terms(
    counts: np.ndarray,
    linear: np.ndarray,
    thresholds: np.ndarray,
) -> tuple[float, np.ndarray, np.ndarray]:
    return _ordinal_terms_with_widths(
        counts, linear, thresholds, _ordinal_interval_log_widths(thresholds)
    )


def _ordinal_terms_with_widths(
    counts: np.ndarray,
    linear: np.ndarray,
    thresholds: np.ndarray,
    interval_log_widths: np.ndarray,
) -> tuple[float, np.ndarray, np.ndarray]:
    with np.errstate(over="ignore"):
        centered = thresholds[None, :] - linear[:, None]
    log_probabilities = _ordinal_log_probabilities(centered, interval_log_widths)
    below = expit(centered)
    above = expit(-centered)
    gradient = np.sum(counts[:, :-1] * above - counts[:, 1:] * below, axis=1)
    curvature = np.sum((counts[:, :-1] + counts[:, 1:]) * below * above, axis=1)
    log_likelihood = float(
        np.sum(
            np.multiply(
                counts,
                log_probabilities,
                out=np.zeros_like(log_probabilities),
                where=counts != 0.0,
            )
        )
    )
    return log_likelihood, gradient, curvature


def _weighted_dense_hessian(precision, weights):
    if weights.ndim == 1:
        return precision + np.diag(weights)
    matrix = precision.copy()
    dimension = weights.shape[1]
    indices = np.arange(len(weights) * dimension).reshape(-1, dimension)
    matrix[indices[:, :, None], indices[:, None, :]] += weights
    return matrix


def _sparse_weight_matrix(weights: np.ndarray) -> sparse.csr_matrix:
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


def _weighted_sparse_hessian(loading, loading_transpose, weights):
    if weights.ndim == 1:
        weighted = loading.multiply(weights[:, None]).tocsr()
        weighted.eliminate_zeros()
        return loading_transpose @ weighted
    return loading_transpose @ _sparse_weight_matrix(weights) @ loading


def _sparse_glmm_latent_model(
    phylogenetic,
    component_variances,
    group_covariance,
    random_dimension,
    predictor_uncertainties=(),
    predictor_columns=(),
    coefficients=None,
    *,
    include_sampling_factor=False,
):
    if phylogenetic is None:
        return None
    components = {"phylogenetic": (component_variances["phylogenetic"], phylogenetic)}
    if isinstance(group_covariance, SparseCovarianceModel):
        components["group"] = (component_variances["group"], group_covariance)
    latent = combine_sparse_covariance_models(
        components,
        random_dimension=random_dimension,
        include_sampling_factor=include_sampling_factor,
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
                    (
                        uncertainty.model.precision_factor()
                        if include_sampling_factor
                        else None
                    ),
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
    loading_transpose = loading.T

    def evaluate(candidate):
        linear = fixed_linear + np.asarray(loading @ candidate).reshape(-1)
        log_likelihood, likelihood_gradient, weights = likelihood_terms(linear)
        prior_score = precision @ candidate
        gradient = np.asarray(loading_transpose @ likelihood_gradient).reshape(-1)
        gradient += prior_score
        value = -log_likelihood + 0.5 * float(candidate @ prior_score)
        return value, gradient, weights

    def factor_hessian(weights):
        hessian = precision + _weighted_sparse_hessian(
            loading, loading_transpose, weights
        )
        return factor_sparse_positive_definite(hessian)

    state = evaluate(mode)
    converged = False
    for _iteration in range(80):
        value, gradient, weights = state
        if not np.isfinite(gradient).all():
            break
        if float(np.max(np.abs(gradient))) < 1e-7:
            converged = True
            break
        step = factor_hessian(weights).solve(gradient)
        accepted = False
        step_scale = 1.0
        for _backtrack in range(24):
            candidate = mode - step_scale * step
            candidate_state = evaluate(candidate)
            candidate_value = candidate_state[0]
            if np.isfinite(candidate_value) and candidate_value <= value:
                mode = candidate
                state = candidate_state
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
            value, gradient, _weights = evaluate(candidate)
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
        state = evaluate(mode)
    value, gradient, weights = state
    if not np.isfinite(value) or float(np.max(np.abs(gradient))) >= 1e-5:
        raise RuntimeError("{} random-effect mode did not converge.".format(label))
    return mode, value, factor_hessian(weights)


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
    interval_log_widths = _ordinal_interval_log_widths(thresholds)

    def terms(candidate: np.ndarray) -> tuple[float, np.ndarray, np.ndarray]:
        if family == "ordinal":
            return _ordinal_terms_with_widths(
                counts, fixed_linear + candidate, thresholds, interval_log_widths
            )
        dimension = fixed_linear.shape[1]
        return _multinomial_terms(
            counts, fixed_linear + candidate.reshape(len(counts), dimension)
        )

    return _dense_latent_mode(precision, terms, "Categorical", mode)


def _dense_latent_mode(precision, likelihood_terms, label, mode):
    def evaluate(candidate):
        log_likelihood, likelihood_gradient, weights = likelihood_terms(candidate)
        value = -log_likelihood + 0.5 * float(candidate @ precision @ candidate)
        return value, likelihood_gradient, weights

    def factor_hessian(weights):
        return _positive_definite_cholesky(_weighted_dense_hessian(precision, weights))

    state = evaluate(mode)
    converged = False
    for _iteration in range(80):
        value, likelihood_gradient, weights = state
        gradient = likelihood_gradient + precision @ mode
        if not np.isfinite(gradient).all():
            break
        if float(np.max(np.abs(gradient))) < 1e-7:
            converged = True
            break
        try:
            cholesky = factor_hessian(weights)
        except np.linalg.LinAlgError:
            break
        step = cho_solve((cholesky, True), gradient, check_finite=False)
        accepted = False
        step_scale = 1.0
        for _backtrack in range(24):
            candidate = mode - step_scale * step
            candidate_state = evaluate(candidate)
            candidate_value = candidate_state[0]
            if np.isfinite(candidate_value) and candidate_value <= value:
                mode = candidate
                state = candidate_state
                accepted = True
                break
            step_scale *= 0.5
        if not accepted:
            break
        if float(np.max(np.abs(step_scale * step))) < 1e-8:
            converged = _dense_mode_gradient_is_small(state[1] + precision @ mode)
            break
    if not converged:

        def value_and_gradient(candidate):
            value, likelihood_gradient, _weights = evaluate(candidate)
            return value, likelihood_gradient + precision @ candidate

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
        state = evaluate(mode)
    value, _gradient, weights = state
    if not np.isfinite(value):
        raise RuntimeError("{} random-effect mode is non-finite.".format(label))
    return mode, value, factor_hessian(weights)


def _dense_mode_gradient_is_small(gradient):
    return np.isfinite(gradient).all() and float(np.max(np.abs(gradient))) < 1e-5


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
    *,
    gmrf_factors=None,
) -> np.ndarray:
    """Project predictor uncertainty into the multivariate linear predictor."""
    dimension = coefficients.shape[1]
    projected = np.zeros(
        (n_observations * dimension, n_observations * dimension), dtype=float
    )
    if len(uncertainties) != len(columns):
        raise ValueError("Predictor uncertainty and coefficient-column counts differ.")
    for index, (uncertainty, selected_columns) in enumerate(
        zip(uncertainties, columns, strict=True)
    ):
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
                factor = _gmrf_projection_factor(uncertainty, index, gmrf_factors)
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
            blocks = (slopes.T @ covariance) @ slopes
            indices = np.arange(n_observations * dimension).reshape(-1, dimension)
            projected[indices[:, :, None], indices[:, None, :]] += blocks
            continue
        expected = (
            len(selected),
            len(selected),
            n_observations,
            n_observations,
        )
        if covariance.shape != expected:
            raise ValueError("Grouped predictor uncertainty has wrong dimensions.")
        _add_joint_dense_predictor_covariance(projected, slopes, covariance)
    return (projected + projected.T) / 2.0


def _gmrf_projection_factor(uncertainty, index, cached):
    if cached is None:
        return factor_sparse_positive_definite(uncertainty.model.precision)
    if index not in cached:
        cached[index] = factor_sparse_positive_definite(uncertainty.model.precision)
    return cached[index]


def _add_joint_dense_predictor_covariance(projected, slopes, covariance):
    n_observations = covariance.shape[2]
    dimension = slopes.shape[1]
    # Batch over source observations to bound temporary covariance blocks.
    row_size = n_observations * dimension * max(slopes.shape)
    batch = max(1, 262_144 // max(1, row_size))
    for start in range(0, n_observations, batch):
        stop = min(start + batch, n_observations)
        source = covariance[:, :, start:stop, :].transpose(2, 3, 0, 1)
        blocks = (slopes.T @ source) @ slopes
        projected[start * dimension : stop * dimension] += blocks.transpose(
            0, 2, 1, 3
        ).reshape((stop - start) * dimension, -1)


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


def _log_positive_count_probability(log_zero, eta, log_relative_mean):
    # For small mean and relative mean, Pr(Y > 0) = exp(eta) to
    # leading order, even if the count probability underflows.
    tiny_mean = (eta < -36.0) & (log_relative_mean < -36.0)
    with np.errstate(divide="ignore"):
        return np.where(tiny_mean, eta, np.log(-np.expm1(log_zero)))


def _negative_binomial_log_zero(log_means, log_relative_mean, dispersion):
    with np.errstate(over="ignore", under="ignore"):
        tiny = log_relative_mean < -36.0
        log_zero = np.empty_like(log_means)
        log_zero[tiny] = -np.exp(log_means[tiny])
        log_zero[~tiny] = -np.logaddexp(0.0, log_relative_mean[~tiny]) / dispersion
    return log_zero


def _small_dispersion_log_rising_factorial(values, dispersion):
    """Return log(Gamma(size + y) / (Gamma(size) * size**y))."""
    relative_count = values * dispersion
    small = relative_count < 1e-3
    correction = np.empty_like(values)
    z = relative_count[small]
    correction[small] = values[small] * np.polynomial.polynomial.polyval(
        z, _NB_LOG_GAMMA_RATIO_SERIES
    )
    z = relative_count[~small]
    correction[~small] = values[~small] * ((1.0 + z) * (np.log1p(z) / z) - 1.0)
    correction -= 0.5 * np.log1p(relative_count)
    # DLMF 5.11.1. The next Stirling term is bounded by dispersion**3 / 360.
    correction -= (dispersion / 12.0) * (relative_count / (1.0 + relative_count))
    return np.where(values <= 1.0, 0.0, correction)


def _negative_binomial_log_pmf(values, log_means, dispersion):
    log_dispersion = math.log(dispersion)
    log_relative_mean = log_means + log_dispersion
    log_zero = _negative_binomial_log_zero(log_means, log_relative_mean, dispersion)
    if np.max(values, initial=0.0) >= _LARGE_NEGATIVE_BINOMIAL_COUNT:
        return (
            _negative_binomial_with_large_counts(
                values,
                log_means,
                log_relative_mean,
                log_zero,
                dispersion,
                log_dispersion,
            ),
            log_zero,
            log_relative_mean,
        )
    return (
        _ordinary_negative_binomial_log_mass(
            values, log_means, log_relative_mean, log_zero, dispersion
        ),
        log_zero,
        log_relative_mean,
    )


def _ordinary_negative_binomial_log_mass(
    values, log_means, log_relative_mean, log_zero, dispersion
):
    if dispersion < _SMALL_NB_DISPERSION:
        return (
            _small_dispersion_log_rising_factorial(values, dispersion)
            - gammaln(values + 1.0)
            + values * (log_means - np.logaddexp(0.0, log_relative_mean))
            + log_zero
        )
    size = 1.0 / dispersion
    return (
        gammaln(values + size)
        - gammaln(size)
        - gammaln(values + 1.0)
        + log_zero
        - values * np.logaddexp(0.0, -log_relative_mean)
    )


def _negative_binomial_with_large_counts(
    values, log_means, log_relative_mean, log_zero, dispersion, log_dispersion
):
    large = values >= _LARGE_NEGATIVE_BINOMIAL_COUNT
    small = ~large
    result = np.empty_like(values)
    result[small] = _ordinary_negative_binomial_log_mass(
        values[small],
        log_means[small],
        log_relative_mean[small],
        log_zero[small],
        dispersion,
    )
    counts = values[large]
    log_count = np.log(counts)
    log_ratio = log_count + log_dispersion
    # Evaluate the mass at its count-specific mean using Stirling corrections.
    # The remaining log likelihood is the drop from that mean.
    peak = (
        -0.5 * (math.log(math.tau) + log_count + np.logaddexp(0.0, log_ratio))
        - _negative_binomial_size_correction(dispersion)
        - _stirling_log_gamma_correction(counts)
    )
    inverse_total = np.exp(-np.logaddexp(log_count, -log_dispersion))
    peak += _stirling_log_gamma_correction_from_inverse(inverse_total)
    result[large] = peak - _negative_binomial_deviance(
        log_means[large] - log_count, log_count, log_ratio, log_dispersion
    )
    return result


def _negative_binomial_size_correction(dispersion):
    if dispersion <= 1.0 / _LARGE_GAMMA_SHAPE:
        return _stirling_log_gamma_correction_from_inverse(dispersion)
    size = 1.0 / dispersion
    return (
        gammaln(1.0 + size)
        + (size + 0.5) * math.log(dispersion)
        + size
        - 0.5 * math.log(math.tau)
    )


def _negative_binomial_deviance(relative, log_count, log_ratio, log_dispersion):
    with np.errstate(over="ignore", invalid="ignore", under="ignore", divide="ignore"):
        log_q = -np.logaddexp(0.0, log_ratio)
        log_p = -np.logaddexp(0.0, -log_ratio)
        deviance = np.empty_like(relative)
        nearby = np.abs(relative) < 0.1
        deviance[nearby] = _nearby_negative_binomial_deviance(
            relative[nearby], log_count[nearby], log_q[nearby], log_p[nearby]
        )
        far = ~nearby
        change = relative[far]
        previous = log_ratio[far]
        new = previous + change
        deviance[far] = _far_negative_binomial_deviance(
            change,
            previous,
            new,
            log_count[far],
            np.broadcast_to(log_dispersion, relative.shape)[far],
        )
        return np.where(np.isinf(relative), np.inf, deviance)


def _nearby_negative_binomial_deviance(relative, log_count, log_q, log_p):
    # Split log1p(p * expm1(t)) - p*t into two quadratic remainders.
    # Use the smaller of p and q so neither remainder is lost to saturation.
    small_probability = np.exp(np.minimum(log_p, log_q))
    direction = np.where(log_p <= log_q, relative, -relative)
    z = small_probability * np.expm1(direction)
    log_weight = log_count + np.where(log_p <= log_q, 0.0, log_q - log_p)
    remainder = _exponential_deviance(direction, np.exp(log_weight))
    ratio = np.expm1(direction) / np.where(direction == 0.0, 1.0, direction)
    correction = np.polynomial.polynomial.polyval(z, _LOG1P_REMAINDER_SERIES)
    scale = np.exp(log_weight + np.minimum(log_p, log_q))
    return remainder + ((scale * direction) * direction) * ratio**2 * correction


def _far_negative_binomial_deviance(change, previous, new, log_count, log_dispersion):
    # Orient the softplus difference so that its arguments have a nonpositive
    # average. Keep the multiplication by the shape in logarithmic units.
    mostly_failures = (previous + 0.5 * change) > 0.0
    direction = np.where(mostly_failures, -1.0, 1.0)
    first, last = previous * direction, new * direction
    high, low = np.maximum(first, last), np.minimum(first, last)
    softplus_high = np.logaddexp(0.0, high)
    difference = softplus_high - np.logaddexp(0.0, low)
    # log(softplus(high) - softplus(low)). In the extreme tail, both
    # arguments are negative and log(softplus(high)) tends to high.
    tail = high < -36.0
    log_difference = np.log(difference)
    log_difference[tail] = high[tail] + np.log(-np.expm1(low[tail] - high[tail]))
    coefficient = np.logaddexp(log_count, -log_dispersion)
    log_term = coefficient + log_difference
    log_linear = np.log(np.abs(change)) + np.where(
        mostly_failures, -log_dispersion, log_count
    )
    larger, smaller = np.maximum(log_term, log_linear), np.minimum(log_term, log_linear)
    log_deviance = larger + np.log(-np.expm1(smaller - larger))
    return np.exp(log_deviance)


def _poisson_log_mass(values, log_means):
    with np.errstate(over="ignore", invalid="ignore"):
        log_zero = -np.exp(log_means)
        if np.max(values, initial=0.0) < _LARGE_POISSON_COUNT:
            return values * log_means + log_zero - gammaln(values + 1.0), log_zero
        result = _poisson_log_mass_with_large_counts(values, log_means, log_zero)
    return result, log_zero


def _poisson_log_mass_with_large_counts(values, log_means, log_zero):
    large = values >= _LARGE_POISSON_COUNT
    small = ~large
    result = np.empty_like(values)
    y, eta = values[small], log_means[small]
    result[small] = y * eta + log_zero[small] - gammaln(y + 1.0)
    counts = values[large]
    log_counts = np.log(counts)
    # A Stirling correction plus the Poisson deviance avoids cancellation
    # between the factorial and the power of the mean. See Loader (2002).
    result[large] = (
        -0.5 * (math.log(math.tau) + log_counts)
        - _stirling_log_gamma_correction(counts)
        - _exponential_deviance(log_means[large] - log_counts, counts)
    )
    return result


def _count_log_likelihood(values, linear, family, dispersion, zero_probability, offset):
    log_means = linear + offset
    if "negative-binomial" not in family:
        log_pmf, log_zero = _poisson_log_mass(values, log_means)
        log_relative_mean = log_means
    else:
        if dispersion is None or dispersion <= 0.0:
            raise ValueError("Negative-binomial dispersion must be positive.")
        log_pmf, log_zero, log_relative_mean = _negative_binomial_log_pmf(
            values, log_means, dispersion
        )
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
    positive_normalizer = _log_positive_count_probability(
        log_zero, log_means, log_relative_mean
    )
    return np.where(
        is_zero,
        log_pi,
        log_one_minus_pi + log_pmf - positive_normalizer,
    )


def _gamma_log_likelihood(values, linear, dispersion):
    if dispersion is None or dispersion <= 0.0:
        raise ValueError("Gamma shape must be positive.")
    log_values = np.log(values)
    log_ratio = log_values - linear
    return _gamma_log_density(log_values, log_ratio, dispersion)


def _gamma_log_density(log_values, log_ratio, shape):
    if shape >= _LARGE_GAMMA_SHAPE:
        return _large_shape_gamma_log_density(log_values, log_ratio, shape)
    log_shape = math.log(shape)
    log_gamma = gammaln(1.0 + shape) - log_shape if shape < 1.0 else gammaln(shape)
    with np.errstate(over="ignore"):
        shape_scaled_ratio = np.exp(log_shape + log_ratio)
    return (
        shape * log_shape
        - log_gamma
        - log_values
        + shape * log_ratio
        - shape_scaled_ratio
    )


def _stirling_log_gamma_correction(shape):
    return _stirling_log_gamma_correction_from_inverse(1.0 / shape)


def _stirling_log_gamma_correction_from_inverse(inverse):
    inverse_squared = inverse * inverse
    return inverse * (1.0 / 12.0 - inverse_squared * (1.0 / 360.0))


def _large_shape_gamma_log_density(log_values, log_ratio, shape):
    # Cancel the shape-sized terms from the log density using DLMF 5.11.1.
    correction = _stirling_log_gamma_correction(shape)
    normalization = 0.5 * (math.log(shape) - math.log(math.tau)) - correction
    return normalization - log_values - _exponential_deviance(log_ratio, shape)


def _exponential_deviance(log_ratio, scale):
    scale = np.broadcast_to(scale, log_ratio.shape)
    close = np.abs(log_ratio) < 0.1
    relative = log_ratio[close]
    deviance = np.empty_like(log_ratio)
    deviance[close] = (
        (scale[close] * relative) * relative
    ) * np.polynomial.polynomial.polyval(relative, _EXPONENTIAL_DEVIANCE_SERIES)
    relative = log_ratio[~close]
    with np.errstate(over="ignore", invalid="ignore"):
        divergence = np.expm1(relative) - relative
        deviance[~close] = scale[~close] * np.where(
            np.isinf(relative), np.inf, divergence
        )
    return deviance


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


def _beta_mean(linear):
    return np.clip(expit(linear), 1e-10, 1.0 - 1e-10)


def _beta_parameters(linear, dispersion, label):
    if dispersion is None or dispersion <= 0.0:
        raise ValueError("{} precision must be positive.".format(label))
    mean = _beta_mean(linear)
    complement = 1.0 - mean
    return mean, complement, mean * dispersion, complement * dispersion


def _beta_log_likelihood(values, linear, dispersion):
    _, _, alpha, beta = _beta_parameters(linear, dispersion, "Beta")
    return (
        -betaln(alpha, beta)
        + (alpha - 1.0) * np.log(values)
        + (beta - 1.0) * np.log1p(-values)
    )


def _beta_binomial_log_likelihood(values, linear, dispersion, trials):
    if trials is None:
        raise ValueError("Beta-binomial responses require trial counts.")
    _, _, alpha, beta = _beta_parameters(linear, dispersion, "Beta-binomial")
    return (
        gammaln(trials + 1.0)
        - gammaln(values + 1.0)
        - gammaln(trials - values + 1.0)
        + betaln(values + alpha, trials - values + beta)
        - betaln(alpha, beta)
    )


def _small_beta_binomial_log_mass(values, trials, dispersion, mean, complement):
    log_inner_mass = math.log(dispersion) + np.log(mean) + np.log(complement)
    # For two positive counts the limiting probability is the beta
    # interior mass times n / (y * (n-y)).
    result = np.where(values == trials, np.log(mean), np.log(complement))
    interior = (values > 0.0) & (values < trials)
    positive = values[interior]
    totals = trials[interior]
    result[interior] = (
        log_inner_mass[interior]
        + np.log(totals)
        - np.log(positive)
        - np.log(totals - positive)
    )
    result[trials == 0.0] = 0.0
    return result


def _small_beta_binomial_terms(values, linear, dispersion, trials):
    mean = _beta_mean(linear)
    complement = 1.0 - mean
    has_successes = (values > 0.0).astype(float)
    has_failures = (values < trials).astype(float)
    contributions = _small_beta_binomial_log_mass(
        values, trials, dispersion, mean, complement
    )
    interior_mean = (mean > 1e-10) & (mean < 1.0 - 1e-10)
    gradient = (mean * has_failures - complement * has_successes) * interior_mean
    curvature = mean * complement * (has_successes + has_failures) * interior_mean
    return float(np.sum(contributions)), gradient, curvature


def _small_beta_terms(values, linear, dispersion):
    mean = _beta_mean(linear)
    complement = 1.0 - mean
    # Gamma(1 + shape) is one to working precision at these shapes.
    # The beta density is then precision * mean * (1 - mean) / (y * (1 - y)).
    contributions = (
        math.log(dispersion)
        + np.log(mean)
        + np.log(complement)
        - np.log(values)
        - np.log1p(-values)
    )
    interior_mean = (mean > 1e-10) & (mean < 1.0 - 1e-10)
    gradient = (mean - complement) * interior_mean
    curvature = (2.0 * mean * complement) * interior_mean
    return float(np.sum(contributions)), gradient, curvature


def _beta_stirling_terms(shape):
    """Scaled log-gamma, digamma, and trigamma remainders."""
    large = shape >= _LARGE_GAMMA_SHAPE
    inverse = 1.0 / shape[large]
    inverse2 = inverse * inverse
    normalizer = np.empty_like(shape)
    score = np.empty_like(shape)
    curvature = np.empty_like(shape)
    normalizer[large] = _stirling_log_gamma_correction_from_inverse(inverse)
    score[large] = -0.5 - inverse * (1.0 / 12.0 - inverse2 / 120.0)
    curvature[large] = 0.5 + inverse * (1.0 / 6.0 - inverse2 / 30.0)
    small = shape[~large]
    log_shape = np.log(small)
    normalizer[~large] = (
        gammaln(small) - (small - 0.5) * log_shape + small - 0.5 * math.log(math.tau)
    )
    score[~large] = small * (digamma(small) - log_shape)
    curvature[~large] = small * (small * polygamma(1, small) - 1.0)
    return normalizer, score, curvature


@lru_cache(maxsize=1024)
def _precise_log_odds(value, total=1.0):
    with localcontext() as context:
        context.prec = 80
        observed = Decimal.from_float(value)
        logit = observed.ln() - (Decimal.from_float(total) - observed).ln()
        leading = float(logit)
        return leading, float(logit - Decimal.from_float(leading))


def _beta_weighted_divergence(observed, fitted, difference, precision):
    relative = difference / fitted
    close = np.abs(relative) < 0.1
    change = relative[close]
    result = np.empty_like(observed)
    weight = precision * fitted
    result[close] = -((weight[close] * change) * change) * (
        np.polynomial.polynomial.polyval(change, _LOG1P_REMAINDER_SERIES)
    )
    far = ~close
    log_ratio = np.log(observed[far]) - np.log(fitted[far])
    result[far] = weight[far] * (relative[far] - log_ratio)
    return result


def _beta_logit_difference(values, eta, log_values, log_complement):
    observed_logit = log_values - log_complement
    difference = eta - observed_logit
    # Highly concentrated densities resolve differences much smaller than
    # an ulp of the observed logit. Keep the logit's residual near the mode.
    near = np.abs(difference) < 1e-3 * np.abs(observed_logit)
    for index in np.flatnonzero(near):
        leading, trailing = _precise_log_odds(float(values[index]))
        difference[index] = (eta[index] - leading) - trailing
    return difference


def _beta_logit_probabilities(linear):
    logit_lower = math.log(1e-10) - math.log1p(-1e-10)
    max_mean = 1.0 - 1e-10
    logit_upper = math.log(max_mean) - math.log1p(-max_mean)
    eta = np.clip(linear, logit_lower, logit_upper)
    active = (linear > logit_lower) & (linear < logit_upper)
    return eta, expit(eta), expit(-eta), active


def _large_beta_terms(values, linear, precision):
    eta, mean, complement, active = _beta_logit_probabilities(linear)
    log_values, log_complement = np.log(values), np.log1p(-values)
    eta_difference = _beta_logit_difference(values, eta, log_values, log_complement)
    near = np.abs(eta_difference) < 0.1
    # Subtract probabilities through their log odds where the mean is close
    # to the observation, even when expit(eta) rounds to the observation.
    change = np.expm1(eta_difference[near])
    near_values = values[near]
    proportion = 1.0 - values
    difference = np.where(eta > 0.0, proportion - complement, mean - values)
    difference[near] = (near_values * (1.0 - near_values)) * (
        change / (1.0 + near_values * change)
    )
    a_correction, a_score, a_curvature = _beta_stirling_terms(precision * mean)
    b_correction, b_score, b_curvature = _beta_stirling_terms(precision * complement)
    slope = -np.tanh(0.5 * eta)
    information = (precision * mean) * complement
    score_correction = complement * a_score - mean * b_score
    with np.errstate(over="ignore"):
        divergence = _beta_weighted_divergence(values, mean, -difference, precision)
        divergence += _beta_weighted_divergence(
            proportion, complement, difference, precision
        )
        gradient = information * eta_difference + score_correction
        curvature = information * (1.0 + slope * eta_difference) + (
            complement**2 * a_curvature
            + mean**2 * b_curvature
            + slope * score_correction
        )
        density = (
            0.5
            * (
                math.log(precision)
                + np.log(mean)
                + np.log(complement)
                - math.log(math.tau)
            )
            - log_values
            - log_complement
            + _stirling_log_gamma_correction(precision)
            - a_correction
            - b_correction
            - divergence
        )
        likelihood = float(np.sum(density))
    return (
        likelihood,
        np.where(active, gradient, 0.0),
        np.where(active, curvature, 0.0),
    )


def _large_beta_binomial_components(counts, trials, mean, shape, residual):
    """Negative-binomial mass and scaled shape derivative remainders."""
    inverse_shape = 1.0 / shape
    # Both count components have the same NB success probability. Choose it
    # to put their sum at the observed trial count.
    log_trials = np.log(trials)
    log_means = log_trials + np.log(mean)
    log_shape = np.log(shape)
    log_relative_mean = log_trials - (log_shape - np.log(mean))
    log_zero = -shape * np.logaddexp(0.0, log_relative_mean)
    tiny = log_relative_mean < -36.0
    log_zero[tiny] = -(trials[tiny] * mean[tiny])
    large = counts >= _LARGE_POISSON_COUNT
    small = ~large & (shape >= _LARGE_GAMMA_SHAPE)
    mass = np.empty_like(counts)
    small_counts = counts[small]
    mass[small] = (
        _small_dispersion_log_rising_factorial(small_counts, inverse_shape[small])
        - gammaln(small_counts + 1.0)
        + small_counts
        * (log_means[small] - np.logaddexp(0.0, log_relative_mean[small]))
        + log_zero[small]
    )
    ordinary = ~large & ~small
    mass[ordinary] = _beta_binomial_small_component_mass(
        counts[ordinary],
        shape[ordinary],
        log_relative_mean[ordinary],
        log_zero[ordinary],
    )
    selected = counts[large]
    log_count = np.log(selected)
    log_ratio = log_count - log_shape[large]
    log_total = np.logaddexp(log_count, log_shape[large])
    log_difference = log_means[large] - log_count
    with np.errstate(over="ignore"):
        relative = residual[large] / selected
    close = np.abs(relative) < 0.1
    log_difference[close] = np.log1p(relative[close])
    mass[large] = (
        -0.5 * (math.log(math.tau) + log_count + np.logaddexp(0.0, log_ratio))
        - _beta_binomial_shape_correction(shape[large], inverse_shape[large])
        - _stirling_log_gamma_correction(selected)
        + _stirling_log_gamma_correction_from_inverse(np.exp(-log_total))
        - _negative_binomial_deviance(
            log_difference, log_count, log_ratio, -log_shape[large]
        )
    )
    log_factor = _log_count_shape_ratio(counts, shape)
    score = -0.5 * np.expm1(-log_factor) - (inverse_shape / 12.0) * np.expm1(
        -2.0 * log_factor
    )
    curvature = -0.5 * np.expm1(-2.0 * log_factor) - (inverse_shape / 6.0) * np.expm1(
        -3.0 * log_factor
    )
    small_shape = (shape < _LARGE_GAMMA_SHAPE) & (counts > 0.0)
    if np.any(small_shape):
        score[small_shape], curvature[small_shape] = _beta_binomial_small_shape_terms(
            counts[small_shape], shape[small_shape], log_factor[small_shape]
        )
    return mass, score, curvature


def _beta_binomial_shape_correction(shape, inverse_shape):
    correction = np.empty_like(shape)
    large = shape >= _LARGE_GAMMA_SHAPE
    correction[large] = _stirling_log_gamma_correction_from_inverse(
        inverse_shape[large]
    )
    small = shape[~large]
    correction[~large] = (
        gammaln(1.0 + small)
        - (small + 0.5) * np.log(small)
        + small
        - 0.5 * math.log(math.tau)
    )
    return correction


def _beta_binomial_small_component_mass(counts, shape, log_relative_mean, log_zero):
    # Use Gamma(1+a) and one factor of a explicitly when a is tiny.
    positive = counts > 0.0
    log_rising = np.zeros_like(counts)
    log_rising[positive] = (
        gammaln(counts[positive] + shape[positive])
        - gammaln(1.0 + shape[positive])
        + np.log(shape[positive])
    )
    return (
        log_rising
        - gammaln(counts + 1.0)
        - counts * np.logaddexp(0.0, -log_relative_mean)
        + log_zero
    )


def _beta_binomial_small_shape_terms(counts, shape, log_factor):
    total = counts + shape
    # Shift the special-function arguments away from the pole at zero.
    score = 1.0 + shape * (digamma(total) - digamma(1.0 + shape) - log_factor)
    curvature = 1.0 + shape * (
        shape * (polygamma(1, 1.0 + shape) - polygamma(1, total)) - counts / total
    )
    return score, curvature


def _log_count_shape_ratio(counts, shape):
    return np.logaddexp(0.0, np.log(counts) - np.log(shape))


def _large_beta_binomial_terms(values, linear, precision, trials):
    eta, mean, complement, active = _beta_logit_probabilities(linear)
    a, b = precision * mean, precision * complement
    failures = trials - values
    residual = _beta_binomial_residual(values, failures, trials, eta, mean, complement)
    with np.errstate(over="ignore", divide="ignore", under="ignore"):
        x_mass, x_score, x_curvature = _large_beta_binomial_components(
            values, trials, mean, a, -residual
        )
        y_mass, y_score, y_curvature = _large_beta_binomial_components(
            failures, trials, complement, b, residual
        )
        log_trial = np.log(trials)
        relative_scale = np.logaddexp(0.0, log_trial - math.log(precision))
        trial_correction = _stirling_log_gamma_correction(trials)
        small_trials = trials < _LARGE_POISSON_COUNT
        small_n = trials[small_trials]
        log_n = log_trial[small_trials]
        trial_correction[small_trials] = (
            gammaln(small_n + 1.0)
            - (small_n + 0.5) * log_n
            + small_n
            - 0.5 * math.log(math.tau)
        )
        total_mass = (
            -0.5 * (math.log(math.tau) + log_trial + relative_scale)
            - _negative_binomial_size_correction(1.0 / precision)
            - trial_correction
            + _stirling_log_gamma_correction_from_inverse(
                np.exp(-np.logaddexp(log_trial, math.log(precision)))
            )
        )
        mass = x_mass + y_mass - total_mass
        relative_success = values / a
        relative_failure = failures / b
        log_success = _log_count_shape_ratio(values, a)
        log_failure = _log_count_shape_ratio(failures, b)
        difference = log_success - log_failure
        nearby = np.abs(difference) < 0.1
        # The scaled ratio keeps the denominator finite near the upper limit.
        bounded = nearby & np.isfinite(relative_failure)
        difference[bounded] = np.log1p(
            (residual[bounded] / (1.0 + relative_failure[bounded]))
            / (a[bounded] * complement[bounded])
        )
        huge = nearby & ~bounded
        difference[huge] = np.log1p((residual[huge] / failures[huge]) / mean[huge])
        gradient = (
            -(a * complement) * difference - complement * x_score + mean * y_score
        )
        slope = complement - mean
        curvature_correction = (
            complement**2 * x_curvature
            + mean**2 * y_curvature
            + slope * (-complement * x_score + mean * y_score)
        )
        curvature = (
            complement**2 * (a * -np.expm1(-log_success))
            + mean**2 * (b * -np.expm1(-log_failure))
            - slope * (a * complement) * difference
        ) + curvature_correction
        # In the binomial limit the leading terms combine to n*m*(1-m).
        # Expand log1p(s)-log1p(f) around s-f before this combination, so
        # subtracting a saturated score cannot erase the remaining curvature.
        small = np.maximum(relative_success, relative_failure) < 0.1
        ratio = np.expm1(difference[small])
        log_remainder = (ratio * ratio) * np.polynomial.polynomial.polyval(
            ratio, _LOG1P_REMAINDER_SERIES
        )
        s, f = relative_success[small], relative_failure[small]
        scaled_success = values[small] * (s / (1.0 + s))
        scaled_failure = failures[small] * (f / (1.0 + f))
        score_adjustment = (
            residual[small] * (f / (1.0 + f))
            - (a[small] * complement[small]) * log_remainder
        )
        curvature[small] = (
            (trials[small] * mean[small]) * complement[small]
            - complement[small] ** 2 * scaled_success
            - mean[small] ** 2 * scaled_failure
            + slope[small] * score_adjustment
            + curvature_correction[small]
        )
        likelihood = float(np.sum(mass))
    return (
        likelihood,
        np.where(active, gradient, 0.0),
        np.where(active, curvature, 0.0),
    )


def _beta_binomial_residual(values, failures, trials, eta, mean, complement):
    residual = np.where(
        eta > 0.0,
        trials * complement - failures,
        values - trials * mean,
    )
    interior = (values > 0.0) & (failures > 0.0)
    near = interior & (np.abs(residual) < 1e-3 * np.minimum(values, failures))
    for index in np.flatnonzero(near):
        leading, trailing = _precise_log_odds(
            float(values[index]), float(trials[index])
        )
        change = math.expm1((eta[index] - leading) - trailing)
        p = values[index] / trials[index]
        residual[index] = -(p * failures[index]) * (change / (1.0 + p * change))
    return residual


def _standardized_censor_bounds(bounds, center, dispersion):
    difference = bounds - center
    result = difference / dispersion
    if dispersion > 1.0:
        overflowing = np.isinf(difference)
        result[overflowing] = bounds[overflowing] / dispersion - (
            center[overflowing] / dispersion
        )
    return result


def _normal_censor_interval_terms(lower, upper, center, dispersion, *, derivatives):
    with np.errstate(over="ignore", under="ignore", invalid="ignore", divide="ignore"):
        standardized_lower = _standardized_censor_bounds(lower, center, dispersion)
        standardized_upper = _standardized_censor_bounds(upper, center, dispersion)
        reflect = standardized_lower >= 0.0
        upper_cdf = log_ndtr(np.where(reflect, -standardized_lower, standardized_upper))
        lower_cdf = log_ndtr(np.where(reflect, -standardized_upper, standardized_lower))
        difference = np.minimum(lower_cdf - upper_cdf, 0.0)
        log_mass = upper_cdf + np.log(-np.expm1(difference))

        width = upper - lower
        log_width = np.log(width) - math.log(dispersion)
        midpoint = 0.5 * standardized_lower + 0.5 * standardized_upper
        # Integrate the quadratic expansion of the density for short intervals.
        # The original endpoints preserve width after centering and scaling.
        small_interval = log_width + np.log1p(np.abs(midpoint)) < math.log(1e-3)
        standardized_width = np.exp(np.where(small_interval, log_width, -np.inf))
        scaled_midpoint = np.where(small_interval, midpoint, 0.0) * standardized_width
        correction = np.log1p((scaled_midpoint**2 - standardized_width**2) / 24.0)
        log_density = -(0.5 * midpoint) * midpoint - 0.5 * math.log(math.tau)
        midpoint_mass = log_width + log_density
        both_outside_range = np.isneginf(upper_cdf)
        log_mass[both_outside_range] = -np.inf
        midpoint_rule = small_interval | (
            (lower_cdf == upper_cdf) & ~both_outside_range
        )
        log_mass = np.where(midpoint_rule, midpoint_mass + correction, log_mass)
        if not derivatives:
            return log_mass, None, None
        gradient, curvature = _normal_interval_derivatives(
            standardized_lower,
            standardized_upper,
            log_mass,
            midpoint_rule,
            midpoint,
            standardized_width,
        )
        return log_mass, gradient, curvature


def _normal_interval_derivatives(
    lower, upper, log_mass, midpoint_rule, midpoint, width
):
    upper_tail = (lower > 8.0) & np.isfinite(lower) & ~midpoint_rule
    lower_tail = (upper < -8.0) & np.isfinite(upper) & ~midpoint_rule
    regular = ~(midpoint_rule | upper_tail | lower_tail)
    a, b, log_p = lower[regular], upper[regular], log_mass[regular]
    log_a = -(0.5 * a) * a - 0.5 * math.log(math.tau)
    log_b = -(0.5 * b) * b - 0.5 * math.log(math.tau)
    ratio_a, ratio_b = np.exp(log_a - log_p), np.exp(log_b - log_p)
    finite_a, finite_b = (
        np.where(np.isfinite(a), a, 0.0),
        np.where(np.isfinite(b), b, 0.0),
    )
    mean = ratio_a - ratio_b
    gradient = np.empty_like(log_mass)
    curvature = np.empty_like(log_mass)
    gradient[regular] = -mean
    curvature[regular] = np.clip(
        mean**2 - finite_a * ratio_a + finite_b * ratio_b, 0.0, 1.0
    )
    mean, curvature[upper_tail] = _normal_interval_tail_moments(
        lower[upper_tail], upper[upper_tail]
    )
    gradient[upper_tail] = -mean
    mean, curvature[lower_tail] = _normal_interval_tail_moments(
        -upper[lower_tail], -lower[lower_tail]
    )
    gradient[lower_tail] = mean
    _set_small_interval_derivatives(gradient, curvature, midpoint_rule, midpoint, width)
    return gradient, curvature


def _normal_interval_tail_moments(lower, upper):
    """Standard normal mean and one minus variance in a right-tail interval."""
    width = upper - lower
    log_density_ratio = -width * (lower + 0.5 * width)
    density_ratio = np.exp(log_density_ratio)
    mills_lower = math.sqrt(math.pi / 2.0) * erfcx(lower / math.sqrt(2.0))
    mills_upper = math.sqrt(math.pi / 2.0) * erfcx(upper / math.sqrt(2.0))
    mass_ratio = mills_lower - density_ratio * mills_upper
    near = _inverse_mills_fraction(lower)
    excess = 1.0 / (lower + near)
    q = density_ratio * mills_upper / mass_ratio
    delta = np.where(np.isfinite(width), width, 0.0)
    upper_excess = 1.0 / (upper + _inverse_mills_fraction(upper))
    excess_difference = upper_excess - excess
    adjustment = q * (delta + excess_difference)
    shift = excess - adjustment
    mean = lower + shift
    tail_weight = 1.0 + excess * (excess - near)
    # Compute the variance around the lower endpoint. The one-sided moments
    # at the upper endpoint are shifted by width before their mass is removed.
    curvature = (
        tail_weight
        + (q * delta) * (delta + upper_excess - 2.0 * excess)
        - (q * excess_difference) * (lower + 2.0 * excess)
        + adjustment**2
    )
    return mean, np.clip(curvature, 0.0, 1.0)


def _set_small_interval_derivatives(gradient, curvature, mask, midpoint, width):
    m, w = midpoint[mask], width[mask]
    # Differentiate the same midpoint log mass used for the likelihood.
    normalizer = 1.0 + ((m * w) ** 2 - w**2) / 24.0
    adjustment = (w**2 / 12.0) / normalizer
    gradient[mask] = -(m - m * adjustment)
    curvature[mask] = 1.0 - adjustment + (m * adjustment) ** 2


def _censor_observations(size, lower, upper):
    lower = np.full(size, np.nan) if lower is None else lower
    upper = np.full(size, np.nan) if upper is None else upper
    lower_missing, upper_missing = np.isnan(lower), np.isnan(upper)
    return (
        lower,
        upper,
        lower_missing & upper_missing,
        lower_missing & ~upper_missing,
        ~lower_missing & upper_missing,
        ~lower_missing & ~upper_missing,
    )


def _censored_gaussian_log_likelihood(
    values, linear, dispersion, censor_lower, censor_upper
):
    if dispersion is None or dispersion <= 0.0:
        raise ValueError("Censored-Gaussian SD must be positive.")
    lower, upper, exact, left, right, interval = _censor_observations(
        len(values), censor_lower, censor_upper
    )
    result, _, _ = _censored_gaussian_log_contributions(
        values,
        linear,
        dispersion,
        lower,
        upper,
        exact,
        left,
        right,
        interval,
        derivatives=False,
    )
    return result


def _censored_gaussian_log_contributions(
    values,
    linear,
    dispersion,
    lower,
    upper,
    exact,
    left,
    right,
    interval,
    *,
    derivatives,
):
    result = np.empty(len(values), dtype=float)
    with np.errstate(over="ignore"):
        standardized_exact = _standardized_censor_bounds(
            values[exact], linear[exact], dispersion
        )
        result[exact] = (
            (-0.5 * standardized_exact) * standardized_exact
            - math.log(dispersion)
            - 0.5 * math.log(math.tau)
        )
        result[left] = log_ndtr(
            _standardized_censor_bounds(upper[left], linear[left], dispersion)
        )
        result[right] = log_ndtr(
            _standardized_censor_bounds(linear[right], lower[right], dispersion)
        )
    interval_score, interval_curvature = None, None
    if np.any(interval):
        interval_mass, interval_score, interval_curvature = (
            _normal_censor_interval_terms(
                lower[interval],
                upper[interval],
                linear[interval],
                dispersion,
                derivatives=derivatives,
            )
        )
        result[interval] = interval_mass
    return result, interval_score, interval_curvature


def _scalar_log_likelihood_contributions(
    values: np.ndarray,
    linear: np.ndarray,
    *,
    family: str,
    dispersion: float | None,
    zero_probability: float | None,
    offset: np.ndarray,
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
    raise ValueError("Unsupported scalar response family: {}.".format(family))


def _poisson_scalar_terms(values, linear, offset):
    with np.errstate(over="ignore"):
        means = np.exp(linear + offset)
    return means - values, means


def _negative_binomial_score_terms(log_means, log_relative_mean, log_dispersion):
    log_denominator = np.logaddexp(0.0, log_relative_mean)
    with np.errstate(over="ignore", under="ignore"):
        log_scaled_mean = np.where(
            log_relative_mean <= 0.0,
            log_means - log_denominator,
            -log_dispersion - np.logaddexp(0.0, -log_relative_mean),
        )
        scaled_mean = np.exp(log_scaled_mean)
        complement = np.exp(-log_denominator)
        zero_curvature = np.exp(log_scaled_mean - log_denominator)
    return scaled_mean, complement, zero_curvature


def _negative_binomial_scalar_terms(values, linear, offset, dispersion):
    assert dispersion is not None
    log_means = linear + offset
    log_dispersion = math.log(dispersion)
    log_relative_mean = log_means + log_dispersion
    scaled_mean, complement, zero_curvature = _negative_binomial_score_terms(
        log_means, log_relative_mean, log_dispersion
    )
    gradient = scaled_mean - values * complement
    curvature = zero_curvature + (values * complement) * expit(log_relative_mean)
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
        log_relative_mean = eta + np.log(np.longdouble(dispersion))
        log_denominator = _softplus_extended(log_relative_mean)
        with np.errstate(over="ignore", under="ignore"):
            fraction = np.exp(-_softplus_extended(-log_relative_mean))
            q = np.exp(eta - log_denominator)
            tiny = log_relative_mean < -36.0
            exponent = np.empty_like(eta)
            exponent[tiny] = np.exp(eta[tiny])
            exponent[~tiny] = log_denominator[~tiny] / dispersion
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
        log_dispersion = math.log(dispersion)
        log_relative_mean = eta + log_dispersion
        zero_log_probability = _negative_binomial_log_zero(
            eta, log_relative_mean, dispersion
        )
        scaled_mean, complement, zero_curvature = _negative_binomial_score_terms(
            eta, log_relative_mean, log_dispersion
        )
        zero_score = -scaled_mean
        zero_second = -zero_curvature
        base_gradient = scaled_mean - values * complement
        base_curvature = zero_curvature + (values * complement) * expit(
            log_relative_mean
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
    log_ratio = np.log(values) - linear
    with np.errstate(over="ignore"):
        curvature = np.exp(math.log(dispersion) + log_ratio)
        gradient = np.where(
            np.abs(log_ratio) < 0.5,
            -dispersion * np.expm1(log_ratio),
            dispersion - curvature,
        )
    return gradient, curvature


def _lognormal_scalar_terms(values, linear, dispersion):
    assert dispersion is not None
    difference = linear - np.log(values)
    with np.errstate(over="ignore"):
        gradient = (difference / dispersion) / dispersion
        curvature = np.full(len(values), np.divide(1.0, dispersion) / dispersion)
    return gradient, curvature


def _beta_scalar_terms(values, linear, dispersion):
    assert dispersion is not None
    means, complement, alpha, beta = _beta_parameters(linear, dispersion, "Beta")
    likelihood_derivative = dispersion * (
        -digamma(alpha) + digamma(beta) + np.log(values) - np.log1p(-values)
    )
    likelihood_second = -(dispersion**2) * (polygamma(1, alpha) + polygamma(1, beta))
    first = means * complement
    second = first * (1.0 - 2.0 * means)
    gradient = -likelihood_derivative * first
    curvature = -(likelihood_second * first**2 + likelihood_derivative * second)
    return gradient, curvature


def _beta_binomial_scalar_terms(values, linear, dispersion, trials):
    assert dispersion is not None
    assert trials is not None
    means, complement, alpha, beta = _beta_parameters(
        linear, dispersion, "Beta-binomial"
    )
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
    first = means * complement
    second = first * (1.0 - 2.0 * means)
    gradient = -likelihood_mean_derivative * first
    curvature = -(
        likelihood_mean_second * first**2 + likelihood_mean_derivative * second
    )
    return gradient, curvature


def _normal_one_sided_score_and_curvature(standardized):
    """Location terms in units of SD for the event Z >= standardized."""
    lower_mask = standardized < 0.0
    far_mask = standardized > 8.0
    center = ~(lower_mask | far_mask)
    lower = standardized[lower_mask]
    middle = standardized[center]
    with np.errstate(over="ignore", under="ignore"):
        lower_ratio = np.exp(
            -(0.5 * lower) * lower - 0.5 * math.log(math.tau) - log_ndtr(-lower)
        )
    middle_ratio = math.sqrt(2.0 / math.pi) / erfcx(middle / math.sqrt(2.0))
    # Evaluate the inverse Mills ratio's excess over z directly in the far
    # tail to avoid cancellation in the Hessian (DLMF 7.9.1).
    far = standardized[far_mask]
    nested = _inverse_mills_fraction(far)
    excess = 1.0 / (far + nested)
    ratio = np.empty_like(standardized)
    curvature = np.empty_like(standardized)
    ratio[lower_mask] = lower_ratio
    ratio[center] = middle_ratio
    ratio[far_mask] = far + excess
    with np.errstate(invalid="ignore"):
        curvature[lower_mask] = lower_ratio * (lower_ratio - lower)
    curvature[center] = middle_ratio * (middle_ratio - middle)
    curvature[far_mask] = 1.0 + excess * (excess - nested)
    curvature[np.isneginf(standardized)] = 0.0
    return ratio, curvature


def _inverse_mills_fraction(standardized):
    fraction = np.zeros_like(standardized)
    for numerator in range(16, 1, -1):
        fraction = numerator / (standardized + fraction)
    return fraction


def _censored_gaussian_scalar_terms(values, linear, dispersion, lower, upper):
    if dispersion is None or dispersion <= 0.0:
        raise ValueError("Censored-Gaussian SD must be positive.")
    lower, upper, exact, left, right, interval = _censor_observations(
        len(values), lower, upper
    )
    current, interval_gradient, interval_curvature = (
        _censored_gaussian_log_contributions(
            values,
            linear,
            dispersion,
            lower,
            upper,
            exact,
            left,
            right,
            interval,
            derivatives=True,
        )
    )
    gradient = np.empty_like(current)
    curvature = np.empty_like(current)
    with np.errstate(over="ignore", divide="ignore"):
        gradient[exact] = (
            _standardized_censor_bounds(linear[exact], values[exact], dispersion)
            / dispersion
        )
        curvature[exact] = np.divide(1.0, dispersion) / dispersion
        _set_one_sided_censor_terms(
            gradient, curvature, linear, dispersion, upper, left, 1.0
        )
        _set_one_sided_censor_terms(
            gradient, curvature, linear, dispersion, lower, right, -1.0
        )
        if interval_gradient is not None:
            gradient[interval] = interval_gradient / dispersion
            curvature[interval] = (interval_curvature / dispersion) / dispersion
    return float(np.sum(current)), gradient, curvature


def _set_one_sided_censor_terms(
    gradient, curvature, linear, dispersion, bounds, mask, sign
):
    standardized = sign * _standardized_censor_bounds(
        linear[mask], bounds[mask], dispersion
    )
    ratio, weight = _normal_one_sided_score_and_curvature(standardized)
    gradient[mask] = (sign * ratio) / dispersion
    curvature[mask] = (weight / dispersion) / dispersion


def _beta_family_terms(values, linear, dispersion, trials, family):
    if family == "beta":
        if dispersion is not None and 0.0 < dispersion < _SMALL_BETA_PRECISION:
            return _small_beta_terms(values, linear, dispersion)
        if dispersion is not None and dispersion >= _LARGE_BETA_PRECISION:
            return _large_beta_terms(values, linear, dispersion)
        contributions = _beta_log_likelihood(values, linear, dispersion)
        gradient, curvature = _beta_scalar_terms(values, linear, dispersion)
        return float(np.sum(contributions)), gradient, curvature
    if trials is None:
        raise ValueError("Beta-binomial responses require trial counts.")
    if dispersion is not None and 0.0 < dispersion < _SMALL_BETA_PRECISION:
        return _small_beta_binomial_terms(values, linear, dispersion, trials)
    if dispersion is not None and dispersion >= _LARGE_GAMMA_SHAPE / 1e-10:
        return _large_beta_binomial_terms(values, linear, dispersion, trials)
    large_trials = trials >= _LARGE_NEGATIVE_BINOMIAL_COUNT
    if dispersion is not None and np.any(large_trials):
        return _mixed_beta_binomial_terms(
            values, linear, dispersion, trials, large_trials
        )
    contributions = _beta_binomial_log_likelihood(values, linear, dispersion, trials)
    gradient, curvature = _beta_binomial_scalar_terms(
        values, linear, dispersion, trials
    )
    return float(np.sum(contributions)), gradient, curvature


def _mixed_beta_binomial_terms(values, linear, precision, trials, large_trials):
    small = ~large_trials
    ordinary_mass = _beta_binomial_log_likelihood(
        values[small], linear[small], precision, trials[small]
    )
    ordinary_gradient, ordinary_curvature = _beta_binomial_scalar_terms(
        values[small], linear[small], precision, trials[small]
    )
    large_mass, large_gradient, large_curvature = _large_beta_binomial_terms(
        values[large_trials], linear[large_trials], precision, trials[large_trials]
    )
    gradient, curvature = np.empty_like(values), np.empty_like(values)
    gradient[small], curvature[small] = ordinary_gradient, ordinary_curvature
    gradient[large_trials], curvature[large_trials] = large_gradient, large_curvature
    return float(np.sum(ordinary_mass)) + large_mass, gradient, curvature


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
    if family in {"beta", "beta-binomial"}:
        return _beta_family_terms(values, linear, dispersion, trials, family)
    if family == "censored-gaussian":
        return _censored_gaussian_scalar_terms(
            values, linear, dispersion, censor_lower, censor_upper
        )

    current = _scalar_log_likelihood_contributions(
        values,
        linear,
        family=family,
        dispersion=dispersion,
        zero_probability=zero_probability,
        offset=offset,
    )
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
    else:
        derivative = _lognormal_scalar_terms(values, linear, dispersion)
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


def _fixed_parameter_covariances(
    objective: Callable[[np.ndarray], float],
    optimized_parameters: np.ndarray,
    coefficient_count: int,
    threshold_count: int,
) -> tuple[np.ndarray, np.ndarray, str]:
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
    n_observations = len(phylogenetic_base)
    latent_score = (precision @ mode).reshape(n_observations, random_dimension)

    def component_mode(name: str, covariance: np.ndarray) -> np.ndarray:
        return (components[name] * covariance) @ latent_score

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


def _coefficient_profile_scale(estimate, standard_error):
    if standard_error > 0.0:
        return max(standard_error, math.ulp(estimate))
    return max(abs(estimate), 0.1)


def _profile_endpoint(
    objective: Callable[[np.ndarray], float],
    optimized: np.ndarray,
    bounds: Sequence[tuple[float | None, float | None]],
    coefficient_index: int,
    standard_error: float,
    fitted_objective: float,
    likelihood_drop: float,
    direction: float,
) -> float:
    estimate = float(optimized[coefficient_index])
    scale = _coefficient_profile_scale(estimate, standard_error)

    def difference(offset: float) -> float:
        value = estimate + direction * scale * offset
        return (
            _optimize_with_fixed_coefficient(
                objective, optimized, bounds, coefficient_index, value
            )
            - fitted_objective
            - likelihood_drop
        )

    outer = 1.0
    for _attempt in range(12):
        value = difference(outer)
        if 0.0 <= value < math.inf:
            endpoint = float(
                brentq(
                    lambda offset: difference(offset) if offset else -likelihood_drop,
                    0.0,
                    outer,
                    xtol=1e-12,
                    maxiter=100,
                )
            )
            return estimate + direction * scale * endpoint
        outer *= 2.0
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
    likelihood_drop = 0.5 * float(chi2.ppf(confidence_level, 1))
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
                fitted_objective,
                likelihood_drop,
                -1.0,
            )
            upper[index] = _profile_endpoint(
                objective,
                optimized,
                bounds,
                index,
                standard_error,
                fitted_objective,
                likelihood_drop,
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
    return _expanded_glmm_covariance(
        tip_covariance,
        random_dimension,
        fit.coefficients,
        predictor_uncertainties,
        predictor_columns,
    )


def _expanded_glmm_covariance(
    tip_covariance,
    random_dimension,
    coefficients,
    predictor_uncertainties,
    predictor_columns,
    *,
    gmrf_factors=None,
):
    covariance = (
        tip_covariance
        if random_dimension == 1
        else np.kron(tip_covariance, np.eye(random_dimension))
    )
    if predictor_uncertainties:
        covariance += _project_predictor_uncertainty(
            coefficients,
            predictor_uncertainties,
            predictor_columns,
            len(tip_covariance),
            gmrf_factors=gmrf_factors,
        )
    return covariance


def _draw_categorical_responses(
    rng: np.random.Generator,
    totals: np.ndarray,
    fixed: np.ndarray,
    fit: PhylogeneticGlmmFit,
    latent: np.ndarray,
) -> list[CategoricalObservation]:
    if fit.family == "ordinal":
        linear = fixed[:, 0] + latent.reshape(len(fixed))
        with np.errstate(over="ignore"):
            centered = fit.thresholds[None, :] - linear[:, None]
        log_widths = _ordinal_interval_log_widths(fit.thresholds)
        probabilities = np.exp(_ordinal_log_probabilities(centered, log_widths))
    else:
        linear = fixed + latent.reshape(fixed.shape)
        probabilities = np.exp(_multinomial_log_probabilities(linear))
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


def _negative_binomial_sampling_dispersion(family, dispersion):
    # At this limit the relative Gamma mixing deviation is below a rounding unit.
    if (
        "negative-binomial" in family
        and dispersion is not None
        and dispersion > _POISSON_LIMIT_DISPERSION
    ):
        return dispersion
    return None


def _draw_negative_binomial(rng, means, dispersion):
    size = 1.0 / dispersion
    if dispersion < _SMALL_NB_DISPERSION:
        # Gamma-Poisson mixing preserves the rate when the NB success
        # probability would round to one.
        mixing = rng.gamma(size, size=means.shape) * dispersion
        return rng.poisson(means * mixing).astype(float)
    probability = size / (size + means)
    return rng.negative_binomial(size, probability).astype(float)


def _draw_count(
    rng: np.random.Generator,
    means: np.ndarray,
    family: str,
    dispersion: float | None,
) -> np.ndarray:
    dispersion = _negative_binomial_sampling_dispersion(family, dispersion)
    if dispersion is not None:
        return _draw_negative_binomial(rng, means, dispersion)
    return rng.poisson(means).astype(float)


def _zero_truncated_negative_binomial(means, dispersion):
    size = 1.0 / dispersion
    with np.errstate(divide="ignore"):
        log_means = np.log(means)
    log_relative_mean = log_means + math.log(dispersion)
    parameters: tuple[float | np.ndarray, ...]
    if dispersion < _SMALL_NB_DISPERSION:
        distribution, parameters = poisson_distribution, (means,)
    else:
        distribution, parameters = nbinom_distribution, (size, size / (size + means))
    log_zero = _negative_binomial_log_zero(log_means, log_relative_mean, dispersion)
    with np.errstate(under="ignore"):
        failure_probability = np.exp(-np.logaddexp(0.0, -log_relative_mean))

    def survival(counts, indices):
        return betainc(counts + 1.0, size, failure_probability[indices])

    return distribution, parameters, -np.expm1(log_zero), survival


def _draw_zero_truncated_count(
    rng: np.random.Generator,
    means: np.ndarray,
    family: str,
    dispersion: float | None,
) -> np.ndarray:
    dispersion = _negative_binomial_sampling_dispersion(family, dispersion)
    if np.any(means < np.finfo(float).tiny):
        return _draw_zero_truncated_count_with_underflow(rng, means, dispersion)
    return _draw_representable_positive_counts(
        rng.random(len(means)), means, dispersion
    )


def _draw_zero_truncated_count_with_underflow(rng, means, dispersion):
    # P(Y > 1 | Y > 0) is below a rounding unit for underflowed means.
    active = means >= np.finfo(float).tiny
    uniforms = rng.random(len(means))
    values = np.ones(len(means), dtype=float)
    if np.any(active):
        values[active] = _draw_representable_positive_counts(
            uniforms[active], means[active], dispersion
        )
    return values


def _draw_representable_positive_counts(uniforms, means, dispersion):
    if dispersion is not None:
        distribution, parameters, positive_probability, survival = (
            _zero_truncated_negative_binomial(means, dispersion)
        )
    else:
        distribution = poisson_distribution
        parameters = (means,)
        positive_probability = -np.expm1(-means)

        def survival(counts, indices):
            return pdtrc(counts, means[indices])

    tail = (1.0 - uniforms) * positive_probability
    quantiles = np.minimum(1.0 - tail, np.nextafter(1.0, 0.0))
    values = distribution.ppf(quantiles, *parameters).astype(float)
    return _verify_positive_count_quantiles(values, tail, survival)


def _verify_positive_count_quantiles(values, tail, survival):
    # Direct CDF quantiles can round to the excluded zero mass or to the wrong
    # count in a narrow tail. Verify the crossing in survival coordinates.
    indices = np.arange(len(values))
    valid = (
        np.isfinite(values)
        & (values >= 1.0)
        & (values <= 2**53)
        & (survival(values, indices) <= tail)
        & ((values == 1.0) | (survival(values - 1.0, indices) > tail))
    )
    for index in indices[~valid]:
        target = tail[index]
        lower, upper = 0, 1
        while survival(float(upper), index) > target:
            lower, upper = upper, 2 * upper
            if upper > 2**53:
                raise ValueError(
                    "A zero-truncated count draw exceeds integer precision."
                )
        while upper - lower > 1:
            middle = (lower + upper) // 2
            if survival(float(middle), index) > target:
                lower = middle
            else:
                upper = middle
        values[index] = upper
    return values


def _draw_mixture_count(rng, linear, family, dispersion, zero_probability, offset):
    means = _count_bootstrap_means(linear, offset)
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


def _count_bootstrap_means(linear, offset):
    with np.errstate(over="ignore", under="ignore"):
        means = np.exp(linear + offset)
    if not np.isfinite(means).all():
        raise RuntimeError("A count bootstrap mean exceeds the floating-point range.")
    return means


def _draw_gamma(rng, linear, dispersion):
    with np.errstate(over="ignore", under="ignore"):
        scale = np.exp(linear) / dispersion
    if _gamma_draw_uses_finite_scale(dispersion, scale):
        values = rng.gamma(dispersion, scale)
    else:
        values = _draw_gamma_in_log_units(rng, linear, dispersion)
    return _check_gamma_draw_range(values)


def _gamma_draw_uses_finite_scale(shape, scale):
    return (
        shape >= 0.1
        and np.isfinite(scale).all()
        and np.all(scale >= np.finfo(float).tiny)
    )


def _draw_gamma_in_log_units(rng, linear, shape):
    log_scale = linear - math.log(shape)
    with np.errstate(over="ignore", under="ignore", divide="ignore"):
        if shape < 1.0:
            # Gamma(shape) = Gamma(shape + 1) * U ** (1 / shape).
            unit_log = np.log(rng.standard_gamma(shape + 1.0, size=linear.shape))
            unit_log += np.log(rng.random(linear.shape)) / shape
        else:
            unit_log = np.log(rng.standard_gamma(shape, size=linear.shape))
        return np.exp(log_scale + unit_log)


def _check_gamma_draw_range(values):
    if not np.isfinite(values).all() or np.any(values <= 0.0):
        raise RuntimeError(
            "Gamma bootstrap draw exceeds the positive floating-point range."
        )
    return values


def _draw_lognormal(rng, linear, dispersion):
    assert dispersion is not None
    return rng.lognormal(linear, dispersion)


def _draw_beta(rng, linear, dispersion):
    _, _, alpha, beta = _beta_parameters(linear, dispersion, "Beta")
    if dispersion < _SMALL_BETA_PRECISION:
        raise RuntimeError(
            "Beta bootstrap draws round to 0 or 1 at the fitted precision, "
            "outside the open response interval."
        )
    return rng.beta(alpha, beta)


def _draw_beta_binomial(rng, linear, dispersion, trials):
    assert dispersion is not None and trials is not None
    return _draw_beta_binomial_with_precision(rng, linear, dispersion, trials)


def _draw_beta_binomial_with_precision(rng, linear, dispersion, trials):
    means = _beta_mean(linear)
    if dispersion < _SMALL_BETA_PRECISION:
        return np.where(rng.random(len(linear)) < means, trials, 0.0)
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


def _prepare_scalar_bootstrap_responses(
    response_values: Sequence[object],
    design: np.ndarray,
    fit: PhylogeneticGlmmFit,
    *,
    offset: Sequence[float] | None,
    trials: Sequence[float] | None,
    censor_lower: Sequence[float] | None,
    censor_upper: Sequence[float] | None,
) -> Callable[[np.random.Generator, np.ndarray], list[object]]:
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
    fixed = (expanded_design @ fit.coefficients)[:, 0]

    def draw(rng, latent):
        linear = fixed + latent.reshape(len(design))[observation_to_tip]
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
        return _pack_scalar_responses(response_values, simulated)

    return draw


def _pack_scalar_responses(response_values, simulated):
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
        log_likelihood, gradient, weights = _scalar_terms(
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
        mapped_gradient = np.bincount(
            observation_to_tip, weights=gradient, minlength=len(mode)
        )
        mapped_weights = np.bincount(
            observation_to_tip, weights=weights, minlength=len(mode)
        )
        return log_likelihood, mapped_gradient, mapped_weights

    return _dense_latent_mode(precision, terms, "Scalar", mode)


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
        with np.errstate(divide="ignore"):
            log_mean = logsumexp(np.log(values) - offset) - math.log(len(values))
        coefficients[0, 0] = max(log_mean, math.log(1e-3))
    elif family in {"gamma", "lognormal"}:
        coefficients[0, 0] = float(np.mean(np.log(values)))
    elif family == "beta":
        mean = float(np.clip(np.mean(values), 1e-4, 1.0 - 1e-4))
        coefficients[0, 0] = math.log(mean / (1.0 - mean))
    elif family == "beta-binomial":
        assert trials is not None
        mean = _initial_beta_binomial_mean(values, trials)
        coefficients[0, 0] = math.log(mean / (1.0 - mean))
    elif family == "censored-gaussian":
        finite = values[np.isfinite(values)]
        coefficients[0, 0] = 0.0 if not len(finite) else float(np.mean(finite))
    return coefficients


def _initial_beta_binomial_mean(values, trials):
    scale = float(np.max(trials))
    if scale > np.finfo(float).max / len(trials):
        values, trials = values / scale, trials / scale
    return float(np.clip(np.sum(values) / np.sum(trials), 1e-4, 1.0 - 1e-4))


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
    response_values, design, phylogenetic_covariance, options: "_GlmmCallOptions"
) -> "PhylogeneticGlmmFit | _GlmmRefitEstimates":
    family = options.family
    group_covariance = options.group_covariance
    predictor_uncertainties, predictor_columns = (
        options.predictor_uncertainties,
        options.predictor_columns,
    )
    evolution_parameter = options.evolution_parameter
    evolution_parameter_initial = options.evolution_parameter_initial
    evolution_parameter_bounds = options.evolution_parameter_bounds
    evolution_parameter_decoder = options.evolution_parameter_decoder
    coefficient_penalty, coefficient_prior_sd = (
        options.coefficient_penalty,
        options.coefficient_prior_sd,
    )
    inference, confidence_level = options.inference, options.confidence_level
    allow_large_dense = options.allow_large_dense
    offset, trials = options.offset, options.trials
    censor_lower, censor_upper = options.censor_lower, options.censor_upper
    dispersion, zero_probability = options.dispersion, options.zero_probability

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
    dense_phylogenetic = _prepare_dense_phylogenetic_covariance(
        phylogenetic_covariance, n_tips
    )
    sparse_phylogenetic = _reuse_tree_covariance(
        phylogenetic_covariance, sparse_builder
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

    gmrf_factors: dict[int, SparsePositiveDefiniteFactor] = {}

    def state(parameters: np.ndarray):
        coefficients, fitted_dispersion, fitted_zero, log_variances, decoded = unpack(
            parameters
        )
        components = dict(zip(component_names, np.exp(log_variances), strict=True))
        if use_sparse:
            latent = _sparse_glmm_latent_model(
                sparse_phylogenetic(decoded),
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

            mode, posterior_objective, posterior_factor = _sparse_latent_mode(
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
            return objective_value, mode, latent, None
        phylo_base = dense_phylogenetic(decoded)
        random_covariance = math.exp(float(log_variances[0])) * phylo_base
        if group_base is not None:
            random_covariance += math.exp(float(log_variances[1])) * group_base
        expanded_covariance = _expanded_glmm_covariance(
            random_covariance,
            1,
            coefficients,
            predictor_uncertainties,
            predictor_columns,
            gmrf_factors=gmrf_factors,
        )
        covariance_cholesky = _positive_definite_cholesky(expanded_covariance)
        precision = cho_solve(
            (covariance_cholesky, True), np.eye(n_tips), check_finite=False
        )
        mode, posterior_objective, hessian_cholesky = _scalar_random_mode(
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
        return objective, mode, precision, phylo_base

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
    if not options.compute_uncertainty:
        return _bootstrap_refit_estimates(
            coefficients, result, coefficient_penalty, coefficient_prior_sd
        )
    fitted_objective, mode, random_effect_state, phylo_base = state(result.x)
    coefficient_covariance, _unused, covariance_status = _fixed_parameter_covariances(
        objective,
        result.x,
        coefficient_count,
        0,
    )
    components = dict(zip(component_names, np.exp(log_variances), strict=True))
    if use_sparse:
        latent = random_effect_state
        tip_mode = np.asarray(latent.loading @ mode).reshape(n_tips, 1)
        component_modes = _sparse_component_modes(
            latent, mode, n_tips, 1, component_names
        )
    else:
        tip_mode = mode.reshape(n_tips, 1)
        component_modes = _glmm_component_modes(
            components, phylo_base, group_base, random_effect_state, mode, 1
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
    compute_uncertainty: bool = True


@dataclass(frozen=True)
class _GlmmRefitEstimates:
    coefficients: np.ndarray
    log_likelihood: float
    optimizer_converged: bool


def _bootstrap_refit_estimates(coefficients, result, penalty, prior_sd):
    return _GlmmRefitEstimates(
        coefficients=coefficients,
        log_likelihood=_reported_log_likelihood(
            result.fun, coefficients, penalty, prior_sd
        ),
        optimizer_converged=bool(result.success),
    )


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


def _prepare_bootstrap_responses(
    response_values,
    design,
    fit,
    options: _GlmmCallOptions,
) -> Callable[[np.random.Generator, np.ndarray], Sequence[object]]:
    if options.family not in SCALAR_RESPONSE_FAMILIES:
        counts = response_count_matrix(response_values, fit.levels)
        totals = np.rint(counts.sum(axis=1)).astype(int)
        fixed = design @ fit.coefficients
        return lambda rng, latent: _draw_categorical_responses(
            rng, totals, fixed, fit, latent
        )
    return _prepare_scalar_bootstrap_responses(
        response_values,
        design,
        fit,
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
    fit = cast(
        PhylogeneticGlmmFit,
        _fit_phylogenetic_glmm(
            response_values, design, phylogenetic_covariance, wald_options
        ),
    )
    design_array = np.asarray(design, dtype=float)
    sparse_builder = getattr(phylogenetic_covariance, "sparse_model", None)
    sparse_group = options.group_covariance is None or isinstance(
        options.group_covariance, SparseCovarianceModel
    )
    sparse_sampler = None
    random_cholesky = None
    if (
        callable(sparse_builder)
        and sparse_group
        and _sparse_predictor_uncertainties(options.predictor_uncertainties)
    ):
        latent_model = _sparse_glmm_latent_model(
            sparse_builder(fit.evolution_parameter),
            fit.component_variances,
            options.group_covariance,
            fit.coefficients.shape[1],
            options.predictor_uncertainties,
            options.predictor_columns,
            fit.coefficients,
            include_sampling_factor=True,
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
        compute_uncertainty=False,
    )
    draw_responses = _prepare_bootstrap_responses(
        response_values, design_array, fit, options
    )
    while len(samples) < options.bootstrap_replicates and attempts < maximum_attempts:
        attempts += 1
        if sparse_sampler is not None:
            latent = sparse_sampler.sample(rng)
        else:
            assert random_cholesky is not None
            latent = random_cholesky @ rng.normal(size=len(random_cholesky))
        simulated = draw_responses(rng, latent)
        try:
            refit = _fit_phylogenetic_glmm(
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
    result = _fit_phylogenetic_glmm(
        response_values, design, phylogenetic_covariance, options
    )
    assert isinstance(result, PhylogeneticGlmmFit)
    return result


def _fit_phylogenetic_glmm(
    response_values, design, phylogenetic_covariance, options: _GlmmCallOptions
) -> PhylogeneticGlmmFit | _GlmmRefitEstimates:
    _validate_glmm_call_options(options)
    if options.inference == "parametric-bootstrap":
        return _fit_parametric_bootstrap_glmm(
            response_values, design, phylogenetic_covariance, options
        )
    if options.family in SCALAR_RESPONSE_FAMILIES:
        return _fit_scalar_phylogenetic_glmm(
            response_values, design, phylogenetic_covariance, options
        )
    return _fit_categorical_phylogenetic_glmm(
        response_values, design, phylogenetic_covariance, options
    )


def _fit_categorical_phylogenetic_glmm(
    response_values, design, phylogenetic_covariance, options: _GlmmCallOptions
) -> PhylogeneticGlmmFit | _GlmmRefitEstimates:
    family = options.family
    levels, reference = options.levels, options.reference
    group_covariance = options.group_covariance
    predictor_uncertainties, predictor_columns = (
        options.predictor_uncertainties,
        options.predictor_columns,
    )
    evolution_parameter = options.evolution_parameter
    evolution_parameter_initial = options.evolution_parameter_initial
    evolution_parameter_bounds = options.evolution_parameter_bounds
    evolution_parameter_decoder = options.evolution_parameter_decoder
    coefficient_penalty, coefficient_prior_sd = (
        options.coefficient_penalty,
        options.coefficient_prior_sd,
    )
    inference, confidence_level = options.inference, options.confidence_level
    allow_large_dense = options.allow_large_dense
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
    dense_phylogenetic = _prepare_dense_phylogenetic_covariance(
        phylogenetic_covariance, len(design)
    )
    sparse_phylogenetic = _reuse_tree_covariance(
        phylogenetic_covariance, sparse_builder
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

    gmrf_factors: dict[int, SparsePositiveDefiniteFactor] = {}

    def state(parameters: np.ndarray):
        coefficients, thresholds, log_variances, decoded_parameter = unpack(parameters)
        components = dict(zip(component_names, np.exp(log_variances), strict=True))
        if use_sparse:
            latent = _sparse_glmm_latent_model(
                sparse_phylogenetic(decoded_parameter),
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
            interval_log_widths = _ordinal_interval_log_widths(thresholds)

            def likelihood_terms(linear):
                if family == "ordinal":
                    return _ordinal_terms_with_widths(
                        counts, linear, thresholds, interval_log_widths
                    )
                return _multinomial_terms(
                    counts, linear.reshape(len(design), random_dimension)
                )

            mode, posterior_objective, posterior_factor = _sparse_latent_mode(
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
            return objective_value, mode, latent, None
        phylo_base = dense_phylogenetic(decoded_parameter)
        random_covariance = np.exp(log_variances[0]) * phylo_base
        if group_base is not None:
            random_covariance = (
                random_covariance + np.exp(log_variances[1]) * group_base
            )
        expanded_covariance = _expanded_glmm_covariance(
            random_covariance,
            random_dimension,
            coefficients,
            predictor_uncertainties,
            predictor_columns,
            gmrf_factors=gmrf_factors,
        )
        covariance_cholesky = _positive_definite_cholesky(expanded_covariance)
        precision = cho_solve(
            (covariance_cholesky, True),
            np.eye(len(design) * random_dimension),
            check_finite=False,
        )
        fixed_linear = design @ coefficients
        if family == "ordinal":
            fixed_linear = fixed_linear[:, 0]
        mode, posterior_objective, hessian_cholesky = _random_mode(
            counts,
            fixed_linear,
            precision,
            family=family,
            thresholds=thresholds,
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
        return objective, mode, precision, phylo_base

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
    if not options.compute_uncertainty:
        return _bootstrap_refit_estimates(
            coefficients, result, coefficient_penalty, coefficient_prior_sd
        )
    fitted_objective, mode, random_effect_state, phylo_base = state(result.x)
    (
        coefficient_covariance,
        threshold_covariance,
        covariance_status,
    ) = _fixed_parameter_covariances(
        objective,
        result.x,
        coefficient_count,
        threshold_count,
    )
    variance_values = np.exp(log_variances)
    components = dict(zip(component_names, variance_values, strict=True))
    if use_sparse:
        latent = random_effect_state
        tip_mode = np.asarray(latent.loading @ mode).reshape(
            len(design), random_dimension
        )
        component_random_modes = _sparse_component_modes(
            latent, mode, len(design), random_dimension, component_names
        )
    else:
        tip_mode = mode.reshape(len(design), random_dimension)
        component_random_modes = _glmm_component_modes(
            components,
            phylo_base,
            group_base,
            random_effect_state,
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
