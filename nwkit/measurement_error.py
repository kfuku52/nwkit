"""Latent-predictor Gaussian models for errors-in-variables PGLS."""

import math
import warnings
from dataclasses import dataclass

import numpy as np
from scipy import sparse
from scipy.optimize import minimize, minimize_scalar

from nwkit.gaussian import (
    DiagonalLowRankCovariance,
    DiagonalSparsePrecisionCovariance,
    center_response,
    covariance_marginal_diagonal,
    effective_likelihood_settings,
    factor_covariance,
    factor_diagonal_low_rank_updates,
    factor_diagonal_sparse_precision_updates,
    factor_logdet,
    grouped_average_marginal_logdet,
    grouped_mean_covariance_diagonal,
    is_diagonal,
    materialize_covariance,
    residual_variance_scale,
    solve_factor,
    sparse_precision_update_diagonal,
)
from nwkit.optimization import FitResourceError
from nwkit.sparse_laplace import (
    ContinuousPredictorUncertainty,
    GmrfPredictorUncertainty,
    GroupedPredictorUncertainty,
    JointPredictorUncertainty,
    SparseCovarianceModel,
    continuous_predictor_loading,
    factor_sparse_nonsingular,
    grouped_predictor_loading,
    joint_predictor_loading,
)

MAX_DENSE_EIV_OBSERVATIONS = 2000
MAX_STRUCTURED_EIV_OBSERVATIONS = 5000


@dataclass(frozen=True)
class LatentPredictorPosterior:
    """Conditional distribution of a latent predictor given its noisy estimate."""

    mean: np.ndarray
    covariance: np.ndarray
    prior_mean: float
    evolutionary_rate: float
    log_likelihood: float
    optimizer_converged: bool
    optimizer_message: str
    boundary_warning: bool


@dataclass(frozen=True)
class SparseLatentPredictorPosterior:
    """Sparse counterpart retaining posterior uncertainty as a tree GMRF."""

    mean: np.ndarray
    covariance_model: SparseCovarianceModel
    mean_posterior_variance: float
    prior_mean: float
    evolutionary_rate: float
    log_likelihood: float
    optimizer_converged: bool
    optimizer_message: str
    boundary_warning: bool


def _positive_semidefinite(matrix, label):
    matrix = np.asarray(matrix, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("{} must be a square matrix.".format(label))
    if not np.isfinite(matrix).all():
        raise ValueError("{} must contain finite values.".format(label))
    symmetric = (matrix + matrix.T) / 2.0
    scale = max(1.0, float(np.max(np.abs(symmetric))) if matrix.size else 1.0)
    tolerance = np.finfo(float).eps * scale * max(1, len(matrix)) * 100.0
    if matrix.size and float(np.max(np.abs(matrix - matrix.T))) > tolerance:
        raise ValueError("{} must be symmetric.".format(label))
    eigenvalues = (
        np.diag(symmetric) if is_diagonal(symmetric) else np.linalg.eigvalsh(symmetric)
    )
    if eigenvalues.size and float(eigenvalues.min()) < -tolerance:
        raise ValueError("{} must be positive semidefinite.".format(label))
    symmetric[np.abs(symmetric) < tolerance] = 0.0
    return symmetric


def _solve(cholesky, values):
    return solve_factor(cholesky, values)


def _latent_predictor_likelihood(
    observed,
    evolutionary_covariance,
    sampling_covariance,
    rate,
    *,
    include_intercept,
    return_details=False,
):
    covariance = rate * evolutionary_covariance + sampling_covariance
    covariance = (covariance + covariance.T) / 2.0
    try:
        cholesky = factor_covariance(covariance)
        inverse_observed = _solve(cholesky, observed)
        if include_intercept:
            ones = np.ones(len(observed), dtype=float)
            inverse_ones = _solve(cholesky, ones)
            denominator = float(ones @ inverse_ones)
            if denominator <= 0.0:
                return float("inf")
            prior_mean = float(ones @ inverse_observed / denominator)
        else:
            prior_mean = 0.0
        residual = observed - prior_mean
        quadratic = float(residual @ _solve(cholesky, residual))
        logdet = factor_logdet(cholesky)
    except np.linalg.LinAlgError:
        return float("inf")
    objective = 0.5 * (len(observed) * math.log(2.0 * math.pi) + logdet + quadratic)
    if not return_details:
        return objective
    return objective, prior_mean, cholesky


def fit_latent_predictor(
    observed,
    evolutionary_covariance,
    sampling_covariance,
    *,
    include_intercept,
):
    """Estimate predictor evolutionary rate and condition on its noisy values."""
    observed = np.asarray(observed, dtype=float)
    if observed.ndim != 1 or not len(observed) or not np.isfinite(observed).all():
        raise ValueError("Observed predictor values must be a finite non-empty vector.")
    evolutionary_covariance = _positive_semidefinite(
        evolutionary_covariance, "Predictor evolutionary covariance"
    )
    if isinstance(sampling_covariance, DiagonalLowRankCovariance):
        sampling_covariance = materialize_covariance(sampling_covariance)
    sampling_covariance = np.asarray(sampling_covariance, dtype=float)
    if sampling_covariance.ndim == 1:
        if (
            sampling_covariance.shape != observed.shape
            or not np.isfinite(sampling_covariance).all()
            or np.any(sampling_covariance < 0.0)
        ):
            raise ValueError("Predictor sampling covariance diagonal is invalid.")
        sampling_covariance = np.diag(sampling_covariance)
    else:
        sampling_covariance = _positive_semidefinite(
            sampling_covariance, "Predictor sampling covariance"
        )
    expected_shape = (len(observed), len(observed))
    if evolutionary_covariance.shape != expected_shape:
        raise ValueError("Predictor evolutionary covariance has the wrong dimensions.")
    if sampling_covariance.shape != expected_shape:
        raise ValueError("Predictor sampling covariance has the wrong dimensions.")
    evolutionary_diagonal = np.diag(evolutionary_covariance)
    if not np.any(evolutionary_diagonal > 0.0):
        raise ValueError("Predictor evolutionary covariance has zero diagonal.")
    observed_offset = float(np.mean(observed)) if include_intercept else 0.0
    observed = observed - observed_offset
    observed_scale = max(
        float(np.mean(observed**2)),
        float(np.mean(np.diag(sampling_covariance))),
        np.finfo(float).tiny,
    )
    evolutionary_scale = float(
        np.mean(evolutionary_diagonal[evolutionary_diagonal > 0])
    )
    rate_scale = observed_scale / evolutionary_scale
    lower_rate = max(rate_scale * 1e-12, np.finfo(float).tiny)
    upper_rate = max(rate_scale * 1e6, lower_rate * 1e6)
    log_bounds = (math.log(lower_rate), math.log(upper_rate))

    def objective(log_rate):
        return _latent_predictor_likelihood(
            observed,
            evolutionary_covariance,
            sampling_covariance,
            math.exp(float(log_rate)),
            include_intercept=include_intercept,
        )

    grid = np.linspace(log_bounds[0], log_bounds[1], 21)
    candidates = [(float(objective(value)), float(value), "grid") for value in grid]
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="invalid value encountered in scalar subtract",
            category=RuntimeWarning,
            module=r"scipy\.optimize\._optimize",
        )
        optimized = minimize_scalar(
            objective,
            bounds=log_bounds,
            method="bounded",
            options={"xatol": 1e-8, "maxiter": 1000},
        )
    if math.isfinite(float(optimized.fun)):
        candidates.append(
            (float(optimized.fun), float(optimized.x), str(optimized.message))
        )
    finite = [candidate for candidate in candidates if math.isfinite(candidate[0])]
    if not finite:
        raise ValueError(
            "Predictor evolutionary-rate optimization found no finite fit."
        )
    objective_value, log_rate, message = min(finite, key=lambda candidate: candidate[0])
    rate = math.exp(log_rate)
    details = _latent_predictor_likelihood(
        observed,
        evolutionary_covariance,
        sampling_covariance,
        rate,
        include_intercept=include_intercept,
        return_details=True,
    )
    if not isinstance(details, tuple):
        raise ValueError(
            "Predictor evolutionary-rate optimization produced an invalid fit."
        )
    _, prior_mean, observation_cholesky = details
    prior_covariance = rate * evolutionary_covariance
    residual = observed - prior_mean
    if observation_cholesky.ndim == 1 and is_diagonal(prior_covariance):
        prior_diagonal = np.diag(prior_covariance)
        observation_diagonal = np.square(observation_cholesky)
        gain = prior_diagonal / observation_diagonal
        posterior_mean = prior_mean + gain * residual
        posterior_covariance = np.diag(
            prior_diagonal - np.square(prior_diagonal) / observation_diagonal
        )
    else:
        gain_residual = prior_covariance @ _solve(observation_cholesky, residual)
        posterior_mean = prior_mean + gain_residual
        posterior_covariance = prior_covariance - (
            prior_covariance @ _solve(observation_cholesky, prior_covariance.T)
        )
    posterior_covariance = _positive_semidefinite(
        posterior_covariance, "Latent predictor posterior covariance"
    )
    boundary = bool(rate <= lower_rate * 10.0 or rate >= upper_rate / 10.0)
    return LatentPredictorPosterior(
        mean=np.asarray(posterior_mean + observed_offset, dtype=float),
        covariance=posterior_covariance,
        prior_mean=float(prior_mean + observed_offset),
        evolutionary_rate=float(rate),
        log_likelihood=-float(objective_value),
        optimizer_converged=bool(optimized.success),
        optimizer_message=message,
        boundary_warning=boundary,
    )


def fit_factor_latent_predictor(
    observed,
    evolutionary_variance,
    sampling_covariance: DiagonalLowRankCovariance,
    *,
    include_intercept,
) -> SparseLatentPredictorPosterior:
    """Condition independent evolutionary contrasts on sparse factor errors.

    Reconciled predictor contrasts have a diagonal evolutionary covariance and
    a sampling covariance supplied as ``U @ U.T``.  This routine retains ``U``
    and the posterior latent precision instead of constructing either dense
    covariance.  A nonzero sampling diagonal is intentionally handled by the
    general sparse tree conditioner, not by this exact-factor specialization.
    """
    observed = np.asarray(observed, dtype=float)
    evolutionary_variance = np.asarray(evolutionary_variance, dtype=float)
    sampling_diagonal = np.asarray(sampling_covariance.diagonal, dtype=float)
    loading = sparse.csr_matrix(sampling_covariance.low_rank, dtype=float)
    if (
        observed.ndim != 1
        or not len(observed)
        or observed.shape != evolutionary_variance.shape
        or sampling_diagonal.shape != observed.shape
        or loading.shape[0] != len(observed)
        or not np.isfinite(observed).all()
        or not np.isfinite(evolutionary_variance).all()
        or not np.isfinite(sampling_diagonal).all()
        or not np.isfinite(loading.data).all()
        or np.any(evolutionary_variance <= 0.0)
        or np.any(sampling_diagonal != 0.0)
    ):
        raise ValueError(
            "Factor predictor conditioning requires positive evolutionary "
            "variances and an exact zero-diagonal sampling factor."
        )
    loading.eliminate_zeros()
    if loading.shape[1]:
        nonzero_columns = np.asarray(loading.getnnz(axis=0)).reshape(-1) > 0
        loading = loading[:, nonzero_columns]
    observed_offset = float(np.mean(observed)) if include_intercept else 0.0
    observed = observed - observed_offset
    sampling_variance = np.asarray(loading.multiply(loading).sum(axis=1)).reshape(-1)
    observed_scale = max(
        float(np.mean(observed**2)),
        float(np.mean(sampling_variance)),
        np.finfo(float).tiny,
    )
    evolutionary_scale = float(np.mean(evolutionary_variance))
    rate_scale = observed_scale / evolutionary_scale
    lower_rate = max(rate_scale * 1e-12, np.finfo(float).tiny)
    upper_rate = max(rate_scale * 1e6, lower_rate * 1e6)
    log_bounds = (math.log(lower_rate), math.log(upper_rate))

    def state(log_rate, *, return_details=False):
        rate = math.exp(float(log_rate))
        diagonal = rate * evolutionary_variance
        try:
            factor = factor_diagonal_low_rank_updates(diagonal, [loading])
            inverse_observed = _solve(factor, observed)
            if include_intercept:
                ones = np.ones(len(observed), dtype=float)
                inverse_ones = _solve(factor, ones)
                denominator = float(ones @ inverse_ones)
                if denominator <= 0.0:
                    return float("inf")
                prior_mean = float(ones @ inverse_observed / denominator)
            else:
                prior_mean = 0.0
            residual = observed - prior_mean
            inverse_residual = _solve(factor, residual)
            quadratic = float(residual @ inverse_residual)
            logdet = factor_logdet(factor)
        except (ValueError, np.linalg.LinAlgError):
            return float("inf")
        objective = 0.5 * (len(observed) * math.log(2.0 * math.pi) + logdet + quadratic)
        if not math.isfinite(objective):
            return float("inf")
        if not return_details:
            return objective
        posterior_mean = prior_mean + diagonal * inverse_residual
        return objective, prior_mean, posterior_mean

    grid = np.linspace(log_bounds[0], log_bounds[1], 21)
    candidates = [(float(state(value)), float(value), "grid") for value in grid]
    optimized = minimize_scalar(
        state,
        bounds=log_bounds,
        method="bounded",
        options={"xatol": 1e-8, "maxiter": 1000},
    )
    if math.isfinite(float(optimized.fun)):
        candidates.append(
            (float(optimized.fun), float(optimized.x), str(optimized.message))
        )
    finite = [candidate for candidate in candidates if math.isfinite(candidate[0])]
    if not finite:
        raise ValueError(
            "Factor predictor evolutionary-rate optimization found no finite fit."
        )
    objective_value, log_rate, message = min(finite, key=lambda candidate: candidate[0])
    details = state(log_rate, return_details=True)
    if not isinstance(details, tuple):
        raise ValueError("Factor predictor conditioning produced an invalid fit.")
    _, prior_mean, posterior_mean = details
    rate = math.exp(log_rate)
    prior_diagonal = rate * evolutionary_variance
    if loading.shape[1] == 0:
        posterior_loading = sparse.csr_matrix((len(observed), 1), dtype=float)
        posterior_precision = sparse.eye(1, format="csc")
        posterior_precision_factor = sparse.eye(1, format="csr")
        precision_logdet = 0.0
    else:
        weighted_loading = loading.multiply((1.0 / np.sqrt(prior_diagonal))[:, None])
        posterior_precision = (
            sparse.eye(loading.shape[1], format="csc")
            + weighted_loading.T @ weighted_loading
        ).tocsc()
        precision_factor = factor_sparse_nonsingular(posterior_precision)
        posterior_loading = loading
        posterior_precision_factor = sparse.vstack(
            [sparse.eye(loading.shape[1], format="csr"), weighted_loading],
            format="csr",
        )
        precision_logdet = precision_factor.logdet
    posterior_model = SparseCovarianceModel(
        precision=posterior_precision,
        tip_loading=posterior_loading,
        logdet_covariance=-float(precision_logdet),
        sampling_parent=np.empty(0, dtype=int),
        sampling_transition=np.empty(0, dtype=float),
        sampling_variance=np.empty(0, dtype=float),
        sampling_precision_factor=posterior_precision_factor,
    )
    factor = factor_sparse_nonsingular(posterior_precision)
    posterior_marginal_variance = sparse_precision_update_diagonal(
        posterior_loading,
        posterior_precision,
        solver=factor,
    )
    mean_posterior_variance = float(np.mean(posterior_marginal_variance))
    boundary = bool(rate <= lower_rate * 10.0 or rate >= upper_rate / 10.0)
    return SparseLatentPredictorPosterior(
        mean=np.asarray(posterior_mean + observed_offset, dtype=float),
        covariance_model=posterior_model,
        mean_posterior_variance=mean_posterior_variance,
        prior_mean=float(prior_mean + observed_offset),
        evolutionary_rate=float(rate),
        log_likelihood=-float(objective_value),
        optimizer_converged=bool(optimized.success),
        optimizer_message=message,
        boundary_warning=boundary,
    )


def _sparse_sampling_parts(observed, evolutionary_model, sampling_variance):
    if isinstance(sampling_variance, DiagonalLowRankCovariance):
        sampling = np.asarray(sampling_variance.diagonal, dtype=float)
        sampling_loading = sparse.csr_matrix(sampling_variance.low_rank, dtype=float)
    else:
        sampling = np.asarray(sampling_variance, dtype=float)
        sampling_loading = sparse.csr_matrix((len(observed), 0), dtype=float)
    if (
        observed.shape != (evolutionary_model.n_tips,)
        or not np.isfinite(observed).all()
    ):
        raise ValueError("Observed predictor values have invalid dimensions or values.")
    if (
        sampling.shape != observed.shape
        or not np.isfinite(sampling).all()
        or np.any(sampling <= 0.0)
        or sampling_loading.shape[0] != len(observed)
        or not np.isfinite(sampling_loading.data).all()
    ):
        raise ValueError(
            "Sparse predictor conditioning requires positive diagonal sampling "
            "variance."
        )
    return sampling, sampling_loading


def _sparse_predictor_precision(
    evolutionary_model, variance, nuisance_count, observation_information
):
    prior_precision = sparse.block_diag(
        [
            evolutionary_model.precision / variance,
            sparse.eye(nuisance_count, format="csc"),
        ],
        format="csc",
    )
    return (prior_precision + observation_information).tocsc()


def _sparse_predictor_covariance_model(
    evolutionary_model,
    latent_variance,
    sampling,
    observation_loading,
    observation_information,
    posterior_factor,
):
    nuisance_count = observation_loading.shape[1] - evolutionary_model.n_states
    posterior_precision = _sparse_predictor_precision(
        evolutionary_model,
        latent_variance,
        nuisance_count,
        observation_information,
    )
    augmented_loading = sparse.hstack(
        [
            evolutionary_model.tip_loading,
            sparse.csr_matrix((evolutionary_model.n_tips, nuisance_count), dtype=float),
        ],
        format="csr",
    )
    prior_factor = sparse.hstack(
        [
            evolutionary_model.precision_factor() / math.sqrt(latent_variance),
            sparse.csr_matrix(
                (evolutionary_model.n_states, nuisance_count), dtype=float
            ),
        ],
        format="csr",
    )
    nuisance_factor = sparse.hstack(
        [
            sparse.csr_matrix(
                (nuisance_count, evolutionary_model.n_states), dtype=float
            ),
            sparse.eye(nuisance_count, format="csr"),
        ],
        format="csr",
    )
    observation_factor = (
        sparse.diags(1.0 / np.sqrt(sampling), format="csr") @ observation_loading
    )
    covariance_model = SparseCovarianceModel(
        precision=posterior_precision,
        tip_loading=augmented_loading,
        logdet_covariance=-posterior_factor.logdet,
        sampling_parent=np.empty(0, dtype=int),
        sampling_transition=np.empty(0, dtype=float),
        sampling_variance=np.empty(0, dtype=float),
        sampling_precision_factor=sparse.vstack(
            [prior_factor, nuisance_factor, observation_factor], format="csr"
        ),
    )
    return covariance_model, augmented_loading


def fit_sparse_latent_predictor(
    observed,
    evolutionary_model: SparseCovarianceModel,
    sampling_variance,
    *,
    include_intercept,
) -> SparseLatentPredictorPosterior:
    """Fit and condition a predictor using sparse tree precision throughout."""
    observed = np.asarray(observed, dtype=float)
    sampling, sampling_loading = _sparse_sampling_parts(
        observed, evolutionary_model, sampling_variance
    )
    loading = evolutionary_model.tip_loading
    nuisance_count = sampling_loading.shape[1]
    observation_loading = sparse.hstack([loading, sampling_loading], format="csr")
    inverse_sampling = 1.0 / sampling
    observation_information = (
        observation_loading.T @ sparse.diags(inverse_sampling) @ observation_loading
    )
    observed_offset = float(np.mean(observed)) if include_intercept else 0.0
    observed = observed - observed_offset
    observed_scale = max(
        float(np.mean(observed**2)),
        float(np.mean(sampling)),
        np.finfo(float).tiny,
    )
    lower_variance = max(observed_scale * 1e-12, np.finfo(float).tiny)
    upper_variance = max(observed_scale * 1e6, lower_variance * 1e6)
    log_bounds = (math.log(lower_variance), math.log(upper_variance))

    def state(log_variance, *, return_details=False):
        variance = math.exp(float(log_variance))
        posterior_precision = _sparse_predictor_precision(
            evolutionary_model,
            variance,
            nuisance_count,
            observation_information,
        )
        try:
            factor = factor_sparse_nonsingular(posterior_precision)
        except (RuntimeError, np.linalg.LinAlgError):
            return float("inf")

        def inverse(values):
            values = np.asarray(values, dtype=float)
            weighted = inverse_sampling * values
            latent_rhs = np.asarray(observation_loading.T @ weighted).reshape(-1)
            correction = np.asarray(
                observation_loading @ factor.solve(latent_rhs)
            ).reshape(-1)
            return weighted - inverse_sampling * correction

        inverse_observed = inverse(observed)
        if include_intercept:
            ones = np.ones(len(observed), dtype=float)
            inverse_ones = inverse(ones)
            denominator = float(ones @ inverse_ones)
            if denominator <= 0.0:
                return float("inf")
            prior_mean = float(ones @ inverse_observed / denominator)
        else:
            prior_mean = 0.0
        residual = observed - prior_mean
        quadratic = float(residual @ inverse(residual))
        logdet = (
            float(np.sum(np.log(sampling)))
            + factor.logdet
            + evolutionary_model.logdet_covariance
            + evolutionary_model.n_states * math.log(variance)
        )
        objective = 0.5 * (len(observed) * math.log(2.0 * math.pi) + logdet + quadratic)
        if not math.isfinite(objective):
            return float("inf")
        if not return_details:
            return objective
        latent_rhs = np.asarray(
            observation_loading.T @ (inverse_sampling * residual)
        ).reshape(-1)
        posterior_mean = prior_mean + np.asarray(
            loading @ factor.solve(latent_rhs)[: evolutionary_model.n_states]
        ).reshape(-1)
        return objective, prior_mean, posterior_mean, factor

    grid = np.linspace(log_bounds[0], log_bounds[1], 21)
    candidates = [(float(state(value)), float(value), "grid") for value in grid]
    optimized = minimize_scalar(
        state,
        bounds=log_bounds,
        method="bounded",
        options={"xatol": 1e-8, "maxiter": 1000},
    )
    if math.isfinite(float(optimized.fun)):
        candidates.append(
            (float(optimized.fun), float(optimized.x), str(optimized.message))
        )
    finite = [candidate for candidate in candidates if math.isfinite(candidate[0])]
    if not finite:
        raise ValueError(
            "Sparse predictor evolutionary-rate optimization found no finite fit."
        )
    objective_value, log_variance, message = min(
        finite, key=lambda candidate: candidate[0]
    )
    details = state(log_variance, return_details=True)
    if not isinstance(details, tuple):
        raise ValueError("Sparse predictor conditioning produced an invalid fit.")
    _, prior_mean, posterior_mean, posterior_factor = details
    latent_variance = math.exp(log_variance)
    covariance_model, augmented_loading = _sparse_predictor_covariance_model(
        evolutionary_model,
        latent_variance,
        sampling,
        observation_loading,
        observation_information,
        posterior_factor,
    )
    # Only the mean marginal posterior variance is needed for diagnostics.
    # A deterministic Hutchinson trace estimate avoids n sparse solves and an
    # n-by-n covariance while remaining accurate enough for this diagnostic.
    probe_count = min(64, max(16, int(math.ceil(math.log2(len(observed) + 1))) * 4))
    rng = np.random.default_rng(0)
    probes = rng.choice((-1.0, 1.0), size=(len(observed), probe_count))
    latent_rhs = augmented_loading.T @ probes
    projected = augmented_loading @ posterior_factor.solve(np.asarray(latent_rhs))
    mean_posterior_variance = float(np.mean(probes * projected))
    mean_posterior_variance = max(0.0, mean_posterior_variance)
    boundary = bool(
        latent_variance <= lower_variance * 10.0
        or latent_variance >= upper_variance / 10.0
    )
    return SparseLatentPredictorPosterior(
        mean=np.asarray(posterior_mean + observed_offset, dtype=float),
        covariance_model=covariance_model,
        mean_posterior_variance=mean_posterior_variance,
        prior_mean=float(prior_mean + observed_offset),
        evolutionary_rate=float(latent_variance / evolutionary_model.covariance_scale),
        log_likelihood=-float(objective_value),
        optimizer_converged=bool(optimized.success),
        optimizer_message=message,
        boundary_warning=boundary,
    )


def _finite_difference_hessian(function, point):
    point = np.asarray(point, dtype=float)
    size = len(point)
    steps = np.finfo(float).eps ** 0.25 * np.maximum(1.0, np.abs(point))
    hessian = np.zeros((size, size), dtype=float)
    center = float(function(point))
    for first in range(size):
        first_step = np.zeros(size, dtype=float)
        first_step[first] = steps[first]
        hessian[first, first] = (
            float(function(point + first_step))
            - 2.0 * center
            + float(function(point - first_step))
        ) / (steps[first] ** 2)
        for second in range(first + 1, size):
            second_step = np.zeros(size, dtype=float)
            second_step[second] = steps[second]
            value = (
                float(function(point + first_step + second_step))
                - float(function(point + first_step - second_step))
                - float(function(point - first_step + second_step))
                + float(function(point - first_step - second_step))
            ) / (4.0 * steps[first] * steps[second])
            hessian[first, second] = value
            hessian[second, first] = value
    return (hessian + hessian.T) / 2.0


def _predictor_error_covariance(beta, uncertainty, columns, n_observations):
    """Project scalar or multicolumn predictor uncertainty onto the response."""
    if isinstance(uncertainty, ContinuousPredictorUncertainty):
        if not isinstance(columns, (int, np.integer)):
            raise ValueError("Continuous predictor uncertainty requires one column.")
        loading = continuous_predictor_loading(
            uncertainty, np.asarray([beta[int(columns)]])
        )
        return np.asarray((loading @ loading.T).toarray())
    if isinstance(uncertainty, GroupedPredictorUncertainty):
        if isinstance(columns, (int, np.integer)):
            raise ValueError("Grouped predictor uncertainty requires column indices.")
        selected = np.asarray(tuple(columns), dtype=int)
        loading = grouped_predictor_loading(uncertainty, beta[selected, None])
        return np.asarray((loading @ loading.T).toarray())
    if isinstance(uncertainty, JointPredictorUncertainty):
        if isinstance(columns, (int, np.integer)):
            raise ValueError("Joint predictor uncertainty requires column indices.")
        selected = np.asarray(tuple(columns), dtype=int)
        loading = joint_predictor_loading(uncertainty, beta[selected])
        return np.asarray((loading @ loading.T).toarray())
    if isinstance(uncertainty, GmrfPredictorUncertainty):
        if not isinstance(columns, (int, np.integer)):
            raise ValueError("GMRF predictor uncertainty requires one column.")
        loading, precision, _ = _gmrf_uncertainty_update(
            beta, int(columns), uncertainty
        )
        factor = factor_sparse_nonsingular(precision)
        solved = factor.solve(loading.T.toarray())
        covariance = np.asarray(loading @ solved)
        return (covariance + covariance.T) / 2.0
    if isinstance(columns, (int, np.integer)):
        matrix = _positive_semidefinite(uncertainty, "Predictor posterior covariance")
        if matrix.shape != (n_observations, n_observations):
            raise ValueError("Predictor posterior covariance has the wrong dimensions.")
        return float(beta[int(columns)]) ** 2 * matrix
    columns = tuple(int(column) for column in columns)
    coefficients = beta[np.asarray(columns, dtype=int)]
    uncertainty = np.asarray(uncertainty, dtype=float)
    if uncertainty.shape == (n_observations, len(columns), len(columns)):
        diagonal = np.einsum(
            "j,ijk,k->i", coefficients, uncertainty, coefficients, optimize=True
        )
        covariance = np.diag(diagonal)
    elif uncertainty.shape == (
        len(columns),
        len(columns),
        n_observations,
        n_observations,
    ):
        covariance = np.einsum(
            "j,jkil,k->il",
            coefficients,
            uncertainty,
            coefficients,
            optimize=True,
        )
    else:
        raise ValueError("Grouped predictor uncertainty has the wrong dimensions.")
    return _positive_semidefinite(covariance, "Grouped predictor covariance")


def _predictor_error_variance_diagonal(beta, uncertainty, columns, n_observations):
    """Return only the predictor-error marginal variances."""
    if isinstance(uncertainty, ContinuousPredictorUncertainty):
        assert isinstance(columns, (int, np.integer))
        loading = continuous_predictor_loading(
            uncertainty, np.asarray([beta[int(columns)]])
        )
        return np.asarray(loading.multiply(loading).sum(axis=1)).reshape(-1)
    if isinstance(uncertainty, GroupedPredictorUncertainty):
        assert not isinstance(columns, (int, np.integer))
        selected = np.asarray(tuple(columns), dtype=int)
        loading = grouped_predictor_loading(uncertainty, beta[selected, None])
        return np.asarray(loading.multiply(loading).sum(axis=1)).reshape(-1)
    if isinstance(uncertainty, JointPredictorUncertainty):
        assert not isinstance(columns, (int, np.integer))
        selected = np.asarray(tuple(columns), dtype=int)
        loading = joint_predictor_loading(uncertainty, beta[selected])
        return np.asarray(loading.multiply(loading).sum(axis=1)).reshape(-1)
    if isinstance(uncertainty, GmrfPredictorUncertainty):
        assert isinstance(columns, (int, np.integer))
        loading, precision, _ = _gmrf_uncertainty_update(
            beta, int(columns), uncertainty
        )
        return _sparse_precision_marginal_diagonal(loading, precision)
    return np.diag(
        _predictor_error_covariance(beta, uncertainty, columns, n_observations)
    )


def _validate_conditional_eiv_size(n_observations, allow_large_dense):
    if n_observations > MAX_DENSE_EIV_OBSERVATIONS and not allow_large_dense:
        raise FitResourceError(
            "Conditional errors-in-variables Gaussian fitting is limited to {} "
            "observations (received {}) because its predictor-dependent "
            "covariance is dense; pass allow_large_dense=True to attempt the "
            "allocation.".format(MAX_DENSE_EIV_OBSERVATIONS, n_observations)
        )
    if n_observations > MAX_DENSE_EIV_OBSERVATIONS:
        warnings.warn(
            "Large dense conditional errors-in-variables fitting was explicitly "
            "enabled for {} observations; attempting the allocation.".format(
                n_observations
            ),
            RuntimeWarning,
            stacklevel=3,
        )


def _uses_structured_eiv_covariance(
    fixed_covariance, components, component_factors, predictor_uncertainties
):
    fixed_input = (
        None
        if isinstance(fixed_covariance, DiagonalLowRankCovariance)
        else np.asarray(fixed_covariance, dtype=float)
    )
    return (
        (fixed_input is None or fixed_input.ndim == 1)
        and all(
            matrix is None or np.asarray(matrix).ndim == 1 or name in component_factors
            for name, matrix in components
        )
        and all(
            isinstance(
                uncertainty,
                (
                    ContinuousPredictorUncertainty,
                    GroupedPredictorUncertainty,
                    JointPredictorUncertainty,
                    GmrfPredictorUncertainty,
                ),
            )
            for uncertainty in predictor_uncertainties
        )
    )


def _prepare_eiv_fixed_covariance(fixed_covariance, n_observations):
    if isinstance(fixed_covariance, DiagonalLowRankCovariance):
        diagonal = np.asarray(fixed_covariance.diagonal, dtype=float)
        loading = fixed_covariance.low_rank
        if (
            diagonal.shape != (n_observations,)
            or np.any(diagonal < 0.0)
            or not np.isfinite(diagonal).all()
            or loading.shape[0] != n_observations
        ):
            raise ValueError("Fixed covariance factor is malformed.")
        return fixed_covariance
    fixed_input = np.asarray(fixed_covariance, dtype=float)
    if fixed_input.ndim != 1:
        return _positive_semidefinite(fixed_input, "Fixed covariance")
    if (
        fixed_input.shape != (n_observations,)
        or (not np.isfinite(fixed_input).all())
        or np.any(fixed_input < 0.0)
    ):
        raise ValueError("Fixed covariance diagonal is invalid.")
    return fixed_input


def _prepare_eiv_components(components, component_factors, n_observations):
    normalized_components = []
    component_scales = []
    for name, matrix in components:
        factor = component_factors.get(name)
        if matrix is None:
            if factor is None:
                raise ValueError(
                    "Implicit EIV variance component '{}' requires a factor.".format(
                        name
                    )
                )
            diagonal = (
                np.asarray(factor.multiply(factor).sum(axis=1)).reshape(-1)
                if sparse.issparse(factor)
                else np.sum(np.square(np.asarray(factor, dtype=float)), axis=1)
            )
            prepared = None
        elif np.asarray(matrix).ndim == 1:
            prepared = np.asarray(matrix, dtype=float)
            if (
                prepared.shape != (n_observations,)
                or (not np.isfinite(prepared).all())
                or np.any(prepared < 0.0)
            ):
                raise ValueError(
                    "Variance component '{}' has an invalid diagonal.".format(name)
                )
            diagonal = prepared
        else:
            prepared = _positive_semidefinite(
                matrix, "Variance component '{}'".format(name)
            )
            diagonal = np.diag(prepared)
        positive = diagonal[diagonal > 0.0]
        if not len(positive):
            raise ValueError("Variance component '{}' has zero diagonal.".format(name))
        scale = float(np.mean(positive))
        normalized_components.append(
            (name, None if prepared is None else prepared / scale)
        )
        if factor is not None:
            component_factors[name] = (
                sparse.csr_matrix(factor, dtype=float) / math.sqrt(scale)
                if sparse.issparse(factor)
                else np.asarray(factor, dtype=float) / math.sqrt(scale)
            )
        component_scales.append(scale)
    return normalized_components, component_scales


def _structured_eiv_base(fixed_covariance):
    if isinstance(fixed_covariance, DiagonalLowRankCovariance):
        return np.asarray(fixed_covariance.diagonal, dtype=float).copy(), [
            fixed_covariance.low_rank
        ]
    return np.asarray(fixed_covariance, dtype=float).copy(), []


def _add_diagonal_continuous_uncertainty(diagonal, beta, columns, uncertainty):
    factor = uncertainty.factor
    indices = np.asarray(uncertainty.observation_index, dtype=int)
    factor_array = None if sparse.issparse(factor) else np.asarray(factor, dtype=float)
    if (
        factor_array is None
        or factor_array.ndim != 1
        or len(np.unique(indices)) != len(indices)
    ):
        return False
    if (
        indices.ndim != 1
        or np.any(indices < 0)
        or np.any(indices >= len(factor_array))
        or not np.isfinite(factor_array).all()
    ):
        raise ValueError("Continuous predictor factor is malformed.")
    row_scale = (
        np.ones(len(indices), dtype=float)
        if uncertainty.row_scale is None
        else np.asarray(uncertainty.row_scale, dtype=float)
    )
    if row_scale.shape != indices.shape or not np.isfinite(row_scale).all():
        raise ValueError("Continuous predictor row scale is malformed.")
    diagonal += np.square(float(beta[columns]) * factor_array[indices] * row_scale)
    return True


def _structured_uncertainty_loading(beta, columns, uncertainty):
    if isinstance(uncertainty, ContinuousPredictorUncertainty):
        assert isinstance(columns, int)
        return continuous_predictor_loading(uncertainty, np.asarray([beta[columns]]))
    assert not isinstance(columns, int)
    selected = beta[np.asarray(tuple(columns), dtype=int)]
    if isinstance(uncertainty, GroupedPredictorUncertainty):
        return grouped_predictor_loading(uncertainty, selected[:, None])
    assert isinstance(uncertainty, JointPredictorUncertainty)
    return joint_predictor_loading(uncertainty, selected)


def _gmrf_uncertainty_update(beta, column, uncertainty):
    indices = np.asarray(uncertainty.observation_index, dtype=int)
    model = uncertainty.model
    if (
        indices.ndim != 1
        or np.any(indices < 0)
        or np.any(indices >= model.n_tips)
        or not np.isfinite(float(beta[column]))
    ):
        raise ValueError("GMRF predictor uncertainty is malformed.")
    row_scale = (
        np.ones(len(indices), dtype=float)
        if uncertainty.row_scale is None
        else np.asarray(uncertainty.row_scale, dtype=float)
    )
    if row_scale.shape != indices.shape or not np.isfinite(row_scale).all():
        raise ValueError("GMRF predictor row scale is malformed.")
    loading = sparse.csr_matrix(model.tip_loading[indices], dtype=float).multiply(
        (float(beta[column]) * row_scale)[:, None]
    )
    return loading.tocsr(), model.precision, model.precision_factor()


def _sparse_precision_marginal_diagonal(loading, precision):
    return sparse_precision_update_diagonal(loading, precision)


def _gmrf_likelihood_marginal_profiles(
    predictor_uncertainties,
    predictor_columns,
    n_coefficients,
    likelihood_groups,
    *,
    use_row_marginals,
):
    """Precompute exact unit-slope GMRF marginals for an EIV objective."""
    if not any(
        isinstance(uncertainty, GmrfPredictorUncertainty)
        for uncertainty in predictor_uncertainties
    ):
        return {}
    profiles = {}
    groups = None if likelihood_groups is None else np.asarray(likelihood_groups)
    if groups is not None:
        unique_groups, inverse = np.unique(groups, return_inverse=True)
        if len(unique_groups) == 0 or not np.array_equal(
            unique_groups, np.arange(len(unique_groups))
        ):
            raise ValueError("Likelihood groups must be non-empty and contiguous.")
        counts = np.bincount(inverse, minlength=len(unique_groups))
        aggregation = sparse.csr_matrix(
            (
                1.0 / counts[inverse].astype(float),
                (inverse, np.arange(len(inverse))),
            ),
            shape=(len(unique_groups), len(inverse)),
        )
    else:
        aggregation = None
    for index, (uncertainty, columns) in enumerate(
        zip(predictor_uncertainties, predictor_columns, strict=True)
    ):
        if not isinstance(uncertainty, GmrfPredictorUncertainty):
            continue
        if not isinstance(columns, (int, np.integer)):
            raise ValueError("GMRF predictor uncertainty requires one column.")
        column = int(columns)
        unit_beta = np.zeros(n_coefficients, dtype=float)
        unit_beta[column] = 1.0
        loading, precision, _ = _gmrf_uncertainty_update(unit_beta, column, uncertainty)
        profiles[index] = {
            "column": column,
            "row": (
                sparse_precision_update_diagonal(loading, precision)
                if use_row_marginals
                else None
            ),
            "grouped": (
                sparse_precision_update_diagonal(
                    sparse.csr_matrix(aggregation @ loading), precision
                )
                if aggregation is not None and not use_row_marginals
                else None
            ),
        }
    return profiles


def _factor_sparse_precision_eiv(
    diagonal,
    updates,
    precision_updates,
    *,
    compute_marginal,
    marginal_diagonal=None,
    grouped_marginal_diagonal=None,
):
    latent_updates = list(precision_updates)
    for update in updates:
        loading = sparse.csr_matrix(update, dtype=float)
        identity = sparse.eye(loading.shape[1], format="csc")
        latent_updates.append((loading, identity, identity.tocsr()))
    cholesky = factor_diagonal_sparse_precision_updates(diagonal, latent_updates)
    loading = sparse.hstack(
        [sparse.csr_matrix(update[0], dtype=float) for update in latent_updates],
        format="csr",
    )
    precision = sparse.block_diag(
        [sparse.csc_matrix(update[1], dtype=float) for update in latent_updates],
        format="csc",
    )
    if marginal_diagonal is not None:
        marginal_diagonal = np.asarray(marginal_diagonal, dtype=float)
    elif compute_marginal:
        marginal_diagonal = np.asarray(diagonal, dtype=float).copy()
        start = 0
        for loading_value, precision_value, _ in latent_updates:
            width = loading_value.shape[1]
            marginal_diagonal += _sparse_precision_marginal_diagonal(
                loading[:, start : start + width], precision_value
            )
            start += width
    else:
        marginal_diagonal = None
    covariance = DiagonalSparsePrecisionCovariance(
        diagonal=np.asarray(diagonal, dtype=float),
        loading=loading,
        precision=precision,
        marginal_diagonal=marginal_diagonal,
        grouped_marginal_diagonal=(
            None
            if grouped_marginal_diagonal is None
            else np.asarray(grouped_marginal_diagonal, dtype=float)
        ),
    )
    return covariance, cholesky


def _add_structured_eiv_components(
    diagonal, updates, variances, normalized_components, component_factors
):
    for variance, (name, component) in zip(
        variances, normalized_components, strict=True
    ):
        if component is not None and np.asarray(component).ndim == 1:
            diagonal += variance * component
        elif name in component_factors:
            updates.append(math.sqrt(float(variance)) * component_factors[name])
        else:
            raise ValueError("Structured EIV covariance received a dense component.")


def _factor_structured_eiv(diagonal, updates, n_observations):
    cholesky = factor_diagonal_low_rank_updates(diagonal, updates)
    low_rank = (
        sparse.hstack([sparse.csr_matrix(update) for update in updates], format="csr")
        if updates
        else sparse.csr_matrix((n_observations, 0), dtype=float)
    )
    return DiagonalLowRankCovariance(diagonal, low_rank), cholesky


def _dense_eiv_covariance(
    beta,
    variances,
    fixed_covariance,
    normalized_components,
    predictor_uncertainties,
    predictor_columns,
    n_observations,
):
    fixed_array = np.asarray(fixed_covariance, dtype=float)
    covariance = np.diag(fixed_array) if fixed_array.ndim == 1 else fixed_array.copy()
    for columns, uncertainty in zip(
        predictor_columns, predictor_uncertainties, strict=True
    ):
        covariance += _predictor_error_covariance(
            beta, uncertainty, columns, n_observations
        )
    for variance, (_, component) in zip(variances, normalized_components, strict=True):
        assert component is not None
        component_array = np.asarray(component)
        covariance += variance * (
            np.diag(component_array) if component_array.ndim == 1 else component_array
        )
    return (covariance + covariance.T) / 2.0


def _conditional_eiv_covariance(
    beta,
    variances,
    fixed_covariance,
    normalized_components,
    component_factors,
    predictor_uncertainties,
    predictor_columns,
    n_observations,
    structured,
    *,
    compute_marginal=False,
    gmrf_marginal_profiles=None,
    likelihood_groups=None,
):
    if structured:
        diagonal, updates = _structured_eiv_base(fixed_covariance)
        precision_updates = []
        row_precision_marginal = None
        grouped_precision_marginal = None
        for uncertainty_index, (columns, uncertainty) in enumerate(
            zip(predictor_columns, predictor_uncertainties, strict=True)
        ):
            if isinstance(uncertainty, GmrfPredictorUncertainty):
                assert isinstance(columns, (int, np.integer))
                precision_updates.append(
                    _gmrf_uncertainty_update(beta, int(columns), uncertainty)
                )
                profile = (gmrf_marginal_profiles or {}).get(uncertainty_index)
                if profile is not None:
                    scale = float(beta[int(columns)]) ** 2
                    if profile["row"] is not None:
                        contribution = scale * np.asarray(profile["row"], dtype=float)
                        row_precision_marginal = (
                            contribution
                            if row_precision_marginal is None
                            else row_precision_marginal + contribution
                        )
                    if profile["grouped"] is not None:
                        contribution = scale * np.asarray(
                            profile["grouped"], dtype=float
                        )
                        grouped_precision_marginal = (
                            contribution
                            if grouped_precision_marginal is None
                            else grouped_precision_marginal + contribution
                        )
                continue
            if not (
                isinstance(uncertainty, ContinuousPredictorUncertainty)
                and _add_diagonal_continuous_uncertainty(
                    diagonal, beta, columns, uncertainty
                )
            ):
                updates.append(
                    _structured_uncertainty_loading(beta, columns, uncertainty)
                )
        _add_structured_eiv_components(
            diagonal,
            updates,
            variances,
            normalized_components,
            component_factors,
        )
        if precision_updates:
            nonprecision = DiagonalLowRankCovariance(
                np.asarray(diagonal, dtype=float),
                (
                    sparse.hstack(
                        [sparse.csr_matrix(update) for update in updates],
                        format="csr",
                    )
                    if updates
                    else sparse.csr_matrix((n_observations, 0), dtype=float)
                ),
            )
            marginal_diagonal = None
            if row_precision_marginal is not None:
                marginal_diagonal = (
                    np.asarray(covariance_marginal_diagonal(nonprecision), dtype=float)
                    + row_precision_marginal
                )
            grouped_marginal_diagonal = None
            if grouped_precision_marginal is not None:
                grouped_marginal_diagonal = (
                    np.asarray(
                        grouped_mean_covariance_diagonal(
                            nonprecision, likelihood_groups
                        ),
                        dtype=float,
                    )
                    + grouped_precision_marginal
                )
            return _factor_sparse_precision_eiv(
                diagonal,
                updates,
                precision_updates,
                compute_marginal=compute_marginal,
                marginal_diagonal=marginal_diagonal,
                grouped_marginal_diagonal=grouped_marginal_diagonal,
            )
        return _factor_structured_eiv(diagonal, updates, n_observations)
    covariance = _dense_eiv_covariance(
        beta,
        variances,
        fixed_covariance,
        normalized_components,
        predictor_uncertainties,
        predictor_columns,
        n_observations,
    )
    return covariance, factor_covariance(covariance)


def fit_conditional_eiv_gaussian(
    response,
    design,
    predictor_uncertainties,
    predictor_columns,
    fixed_covariance,
    components,
    *,
    reml,
    starting_parameters=None,
    component_factors=None,
    allow_large_dense=False,
    likelihood_observations=None,
    likelihood_logdet_offset=0.0,
    likelihood_groups=None,
):
    """Fit y | x-hat when latent predictor uncertainty depends on beta."""
    if reml:
        raise ValueError(
            "Predictor-dependent covariance has no standard REML objective; "
            "fit the conditional errors-in-variables model with ML."
        )
    response = np.asarray(response, dtype=float)
    design = np.asarray(design, dtype=float)
    effective_likelihood_count, logdet_weight, likelihood_logdet_offset = (
        effective_likelihood_settings(
            len(response),
            design.shape[1],
            reml,
            likelihood_observations,
            likelihood_logdet_offset,
        )
    )
    if likelihood_groups is not None:
        likelihood_groups = np.asarray(likelihood_groups)
        grouped_count = len(np.unique(likelihood_groups))
        expected_count = (
            len(response)
            if likelihood_observations is None
            else float(likelihood_observations)
        )
        if not math.isclose(
            float(grouped_count), expected_count, rel_tol=1e-12, abs_tol=1e-12
        ):
            raise ValueError("Likelihood groups do not match likelihood_observations.")
    component_factors = {} if component_factors is None else dict(component_factors)
    structured = _uses_structured_eiv_covariance(
        fixed_covariance, components, component_factors, predictor_uncertainties
    )
    if len(response) > MAX_STRUCTURED_EIV_OBSERVATIONS:
        warnings.warn(
            "Structured conditional errors-in-variables Gaussian fitting above "
            "{} observations is outside the routine validation range (received "
            "{}); attempting the fit.".format(
                MAX_STRUCTURED_EIV_OBSERVATIONS, len(response)
            ),
            RuntimeWarning,
            stacklevel=2,
        )
    if len(response) > MAX_DENSE_EIV_OBSERVATIONS and not structured:
        _validate_conditional_eiv_size(len(response), allow_large_dense)
    fixed_covariance = _prepare_eiv_fixed_covariance(fixed_covariance, len(response))
    if len(predictor_uncertainties) != len(predictor_columns):
        raise ValueError("Each uncertain predictor requires coefficient columns.")
    normalized_components, component_scales = _prepare_eiv_components(
        components, component_factors, len(response)
    )
    centered_response, beta_offset = center_response(
        response,
        design,
        excluded_columns={
            int(column)
            for columns in predictor_columns
            for column in np.atleast_1d(columns)
        },
    )
    ordinary_beta = np.linalg.lstsq(design, centered_response, rcond=None)[0]
    ordinary_residual = centered_response - design @ ordinary_beta
    response_scale = residual_variance_scale(
        centered_response,
        ordinary_residual,
        covariance_marginal_diagonal(fixed_covariance),
    )
    lower_variance = max(response_scale * 1e-12, np.finfo(float).tiny)
    upper_variance = max(response_scale * 1e6, lower_variance * 1e6)
    num_coefficients = design.shape[1]
    gmrf_marginal_profiles = _gmrf_likelihood_marginal_profiles(
        predictor_uncertainties,
        predictor_columns,
        num_coefficients,
        likelihood_groups,
        use_row_marginals=(
            likelihood_groups is not None
            and any(
                name == "species_event_variance"
                for name, _matrix in normalized_components
            )
        ),
    )
    variance_bounds = [(math.log(lower_variance), math.log(upper_variance))] * len(
        normalized_components
    )
    bounds = [(None, None)] * num_coefficients + variance_bounds

    def evaluate(parameters, return_details=False):
        parameters = np.asarray(parameters, dtype=float)
        working_beta = parameters[:num_coefficients]
        beta = working_beta + beta_offset
        log_variances = parameters[num_coefficients:]
        variances = np.exp(log_variances)
        try:
            covariance, cholesky = _conditional_eiv_covariance(
                beta,
                variances,
                fixed_covariance,
                normalized_components,
                component_factors,
                predictor_uncertainties,
                predictor_columns,
                len(response),
                structured,
                compute_marginal=return_details,
                gmrf_marginal_profiles=gmrf_marginal_profiles,
                likelihood_groups=likelihood_groups,
            )
        except (ValueError, np.linalg.LinAlgError):
            return float("inf")
        residual = centered_response - design @ working_beta
        try:
            inverse_residual = _solve(cholesky, residual)
            quadratic = float(residual @ inverse_residual)
            covariance_logdet = factor_logdet(cholesky)
            if reml:
                gram = design.T @ _solve(cholesky, design)
                gram_sign, gram_logdet = np.linalg.slogdet(gram)
                if gram_sign <= 0.0:
                    return float("inf")
            else:
                gram_logdet = 0.0
        except np.linalg.LinAlgError:
            return float("inf")
        if likelihood_groups is None:
            determinant_term = logdet_weight * (
                covariance_logdet - likelihood_logdet_offset
            )
        else:
            if any(name == "species_event_variance" for name, _ in components):
                determinant_term = grouped_average_marginal_logdet(
                    covariance,
                    likelihood_groups,
                    precision_factor=cholesky,
                )
            else:
                grouped_variances = grouped_mean_covariance_diagonal(
                    covariance,
                    likelihood_groups,
                    precision_factor=cholesky,
                )
                if (
                    np.any(grouped_variances <= 0.0)
                    or not np.isfinite(grouped_variances).all()
                ):
                    return float("inf")
                determinant_term = float(np.log(grouped_variances).sum())
            determinant_term -= likelihood_logdet_offset
        objective = 0.5 * (
            effective_likelihood_count * math.log(2.0 * math.pi)
            + determinant_term
            + quadratic
            + gram_logdet
        )
        if not math.isfinite(objective):
            return float("inf")
        if not return_details:
            return objective
        component_variances = {
            name: float(variance / scale)
            for variance, (name, _), scale in zip(
                variances,
                normalized_components,
                component_scales,
                strict=True,
            )
        }
        return {
            "objective": objective,
            "beta": beta,
            "residual": residual,
            "inverse_residual": inverse_residual,
            "covariance": covariance,
            "cholesky": cholesky,
            "quadratic": quadratic,
            "component_variances": component_variances,
            "log_variances": log_variances,
            "parameters": np.concatenate([beta, log_variances]),
            "lower_variance": lower_variance,
            "upper_variance": upper_variance,
        }

    if starting_parameters is None:
        starts = [
            np.concatenate(
                [ordinary_beta, np.log(np.repeat(response_scale, len(components)))]
            )
        ]
    else:
        start = np.asarray(starting_parameters, dtype=float).copy()
        start[:num_coefficients] -= beta_offset
        starts = [start]
    candidates = []
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="invalid value encountered in subtract",
            category=RuntimeWarning,
            module=r"scipy\.optimize\._numdiff",
        )
        for start in starts:
            result = minimize(
                evaluate,
                start,
                method="L-BFGS-B",
                bounds=bounds,
                options={"maxiter": 3000, "ftol": 1e-11},
            )
            if math.isfinite(float(result.fun)):
                candidates.append(result)
    if not candidates:
        raise ValueError("Errors-in-variables optimization found no finite fit.")
    result = min(candidates, key=lambda candidate: float(candidate.fun))
    if not result.success:
        fallback = minimize(
            evaluate,
            result.x,
            method="Powell",
            bounds=bounds,
            options={"maxiter": 10000, "xtol": 1e-8, "ftol": 1e-10},
        )
        if math.isfinite(float(fallback.fun)) and float(fallback.fun) <= float(
            result.fun
        ):
            result = fallback
    if not result.success:
        raise ValueError(
            "Errors-in-variables optimization did not converge: {}".format(
                result.message
            )
        )
    details = evaluate(result.x, return_details=True)
    if not isinstance(details, dict):
        raise ValueError("Errors-in-variables optimization produced an invalid fit.")
    hessian = _finite_difference_hessian(evaluate, result.x)
    try:
        parameter_covariance = np.linalg.inv(hessian)
    except np.linalg.LinAlgError:
        parameter_covariance = np.linalg.pinv(hessian)
    beta_covariance = parameter_covariance[:num_coefficients, :num_coefficients]
    beta_covariance = (beta_covariance + beta_covariance.T) / 2.0
    if not np.isfinite(beta_covariance).all():
        raise ValueError("Errors-in-variables coefficient covariance is non-finite.")
    covariance_scale = max(1.0, float(np.max(np.abs(beta_covariance))))
    covariance_tolerance = (
        np.finfo(float).eps * covariance_scale * max(1, num_coefficients) * 1000.0
    )
    minimum_eigenvalue = float(np.linalg.eigvalsh(beta_covariance).min())
    if minimum_eigenvalue < -covariance_tolerance:
        raise ValueError(
            "Errors-in-variables coefficient information is not positive definite."
        )
    if minimum_eigenvalue < 0.0:
        eigenvalues, eigenvectors = np.linalg.eigh(beta_covariance)
        beta_covariance = (eigenvectors * np.maximum(eigenvalues, 0.0)) @ eigenvectors.T
    details["beta_covariance"] = beta_covariance
    details["parameter_covariance"] = parameter_covariance

    fitted_variances = np.exp(details["log_variances"])

    def covariance_for_beta(beta):
        return _conditional_eiv_covariance(
            np.asarray(beta, dtype=float),
            fitted_variances,
            fixed_covariance,
            normalized_components,
            component_factors,
            predictor_uncertainties,
            predictor_columns,
            len(response),
            structured,
            compute_marginal=False,
            gmrf_marginal_profiles=gmrf_marginal_profiles,
            likelihood_groups=likelihood_groups,
        )

    details["covariance_for_beta"] = covariance_for_beta
    details["optimizer_converged"] = bool(result.success)
    details["optimizer_message"] = str(result.message)
    details["reml"] = False
    details["boundary_warning"] = bool(
        np.any(np.exp(details["log_variances"]) <= lower_variance * 10.0)
        or np.any(np.exp(details["log_variances"]) >= upper_variance / 10.0)
    )
    return details
