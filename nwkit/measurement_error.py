"""Latent-predictor Gaussian models for errors-in-variables PGLS."""

import math
import warnings
from dataclasses import dataclass, replace

import numpy as np
from scipy import sparse
from scipy.linalg import solve_triangular
from scipy.optimize import minimize, minimize_scalar

from nwkit.gaussian import (
    DiagonalLowRankCovariance,
    DiagonalSparsePrecisionCovariance,
    ScaledCovarianceFactor,
    center_response,
    covariance_marginal_diagonal,
    effective_likelihood_settings,
    factor_covariance,
    factor_diagonal_low_rank_updates,
    factor_diagonal_sparse_precision_updates,
    factor_logdet,
    finite_mean,
    gaussian_response_units,
    grouped_mean_covariance_diagonal,
    is_diagonal,
    materialize_covariance,
    materialize_low_rank_covariance,
    rescale_gaussian_covariance,
    residual_variance_scale,
    restore_gaussian_response_units,
    scale_covariance_by_slope,
    solve_factor,
    sparse_precision_update_diagonal,
    symmetric_average,
    validate_psd_covariance,
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
    factor_sparse_positive_definite,
    grouped_predictor_loading,
    joint_predictor_loading,
    materialize_sparse_covariance,
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
    return validate_psd_covariance(matrix, label, clip_negative_diagonal=True)


def _solve(cholesky, values):
    return solve_factor(cholesky, values)


def _latent_predictor_units(observed, sampling_covariance):
    # Scaling each standard deviation also scales both covariance dimensions.
    if not np.isfinite(sampling_covariance).all():
        raise ValueError("Predictor sampling covariance must contain finite values.")
    largest = max(
        float(np.max(np.abs(observed))),
        math.sqrt(float(np.max(np.abs(sampling_covariance), initial=0.0))),
    )
    _, exponent = math.frexp(largest)
    scaled = np.ldexp(observed, -exponent)
    sampling = np.ldexp(sampling_covariance, -2 * exponent)
    return scaled, sampling, exponent


def _sparse_predictor_units(observed, sampling, loading):
    largest = max(
        float(np.max(np.abs(observed), initial=0.0)),
        math.sqrt(float(np.max(sampling, initial=0.0))),
        float(np.max(np.abs(loading.data), initial=0.0)),
    )
    _, exponent = math.frexp(largest)
    loading = loading.copy()
    loading.data = np.ldexp(loading.data, -exponent)
    return (
        np.ldexp(observed, -exponent),
        np.ldexp(sampling, -2 * exponent),
        loading,
        exponent,
    )


def _restore_sparse_predictor_covariance(model, exponent):
    if exponent == 0:
        return model
    loading = model.tip_loading.copy()
    loading.data = np.ldexp(loading.data, exponent)
    return replace(model, tip_loading=loading)


def _restore_predictor_rate(rate, exponent):
    try:
        restored = math.ldexp(rate, exponent)
    except OverflowError as exc:
        raise ValueError(
            "The fitted predictor evolutionary rate exceeds floating-point range."
        ) from exc
    if restored == 0.0:
        raise ValueError(
            "The fitted predictor evolutionary rate is below floating-point range."
        )
    return restored


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
    try:
        cholesky = factor_covariance(covariance)
        prior_mean, quadratic = _latent_predictor_profile_moments(
            observed, cholesky, include_intercept
        )
        logdet = factor_logdet(cholesky)
    except np.linalg.LinAlgError:
        return float("inf")
    objective = 0.5 * (len(observed) * math.log(2.0 * math.pi) + logdet + quadratic)
    if not return_details:
        return objective
    return objective, prior_mean, cholesky


def _latent_predictor_profile_moments(observed, cholesky, include_intercept):
    rhs = (
        np.column_stack((observed, np.ones(len(observed))))
        if include_intercept
        else observed[:, None]
    )
    whitened = (
        rhs / cholesky[:, None]
        if cholesky.ndim == 1
        else solve_triangular(cholesky, rhs, lower=True, check_finite=False)
    )
    whitened_values = whitened[:, 0]
    prior_mean = 0.0
    if include_intercept:
        whitened_ones = whitened[:, 1]
        prior_mean = _whitened_predictor_intercept(whitened_ones, whitened_values)
        whitened_values = whitened_values - prior_mean * whitened_ones
    return prior_mean, float(whitened_values @ whitened_values)


def _whitened_predictor_intercept(whitened_ones, whitened_values):
    with np.errstate(over="ignore"):
        denominator = float(whitened_ones @ whitened_ones)
        numerator = float(whitened_ones @ whitened_values)
    if math.isfinite(denominator) and math.isfinite(numerator) and denominator > 0.0:
        return numerator / denominator
    scale = float(np.max(np.abs(whitened_ones)))
    if not math.isfinite(scale) or scale == 0.0:
        raise np.linalg.LinAlgError("Predictor intercept information is invalid.")
    weights = whitened_ones / scale
    return float(weights @ (whitened_values / scale) / (weights @ weights))


def _condition_dense_predictor(prior, sampling, factor, residual, prior_mean):
    weighted_residual = _solve(factor, residual)
    mean = prior_mean + prior @ weighted_residual
    # The conditional covariance is R (P + R)^-1 P and its transpose,
    # with P denoting prior covariance and R sampling covariance. Place
    # the larger matrix inside the solve to preserve small variances.
    larger, smaller = prior, sampling
    if float(np.max(np.diag(sampling))) > float(np.max(np.diag(prior))):
        larger, smaller = sampling, prior
    covariance = smaller @ _solve(factor, larger)
    return mean, symmetric_average(covariance)


def _center_dense_predictor(observed, include_intercept):
    if not include_intercept:
        return observed, 0.0, float(np.mean(observed**2))
    lower, upper = float(np.min(observed)), float(np.max(observed))
    # An unweighted mean offset can erase near-zero observations that
    # dominate the intercept because their measurement error is small.
    offset = max(lower, min(0.0, upper))
    return observed - offset, offset, float(np.var(observed))


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
    observed, sampling_covariance, trait_exponent = _latent_predictor_units(
        observed, sampling_covariance
    )
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
    _, covariance_exponent = math.frexp(float(np.max(evolutionary_diagonal)))
    evolutionary_covariance = np.ldexp(evolutionary_covariance, -covariance_exponent)
    evolutionary_diagonal = np.diag(evolutionary_covariance)
    observed, observed_offset, observed_variance = _center_dense_predictor(
        observed, include_intercept
    )
    observed_scale = max(
        observed_variance,
        float(np.mean(np.diag(sampling_covariance))),
        np.finfo(float).tiny,
    )
    evolutionary_scale = finite_mean(evolutionary_diagonal[evolutionary_diagonal > 0])
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
        posterior_covariance = np.diag(gain * np.diag(sampling_covariance))
    else:
        posterior_mean, posterior_covariance = _condition_dense_predictor(
            prior_covariance,
            sampling_covariance,
            observation_cholesky,
            residual,
            prior_mean,
        )
    posterior_covariance = _positive_semidefinite(
        posterior_covariance, "Latent predictor posterior covariance"
    )
    boundary = bool(rate <= lower_rate * 10.0 or rate >= upper_rate / 10.0)
    return LatentPredictorPosterior(
        mean=np.ldexp(posterior_mean + observed_offset, trait_exponent),
        covariance=np.ldexp(posterior_covariance, 2 * trait_exponent),
        prior_mean=math.ldexp(prior_mean + observed_offset, trait_exponent),
        evolutionary_rate=_restore_predictor_rate(
            rate, 2 * trait_exponent - covariance_exponent
        ),
        log_likelihood=-float(objective_value)
        - len(observed) * trait_exponent * math.log(2.0),
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
    observed, _, loading, trait_exponent = _sparse_predictor_units(
        observed, sampling_diagonal, loading
    )
    observed_offset = float(np.mean(observed)) if include_intercept else 0.0
    observed = observed - observed_offset
    sampling_variance = np.asarray(loading.multiply(loading).sum(axis=1)).reshape(-1)
    observed_scale = max(
        float(np.mean(observed**2)),
        float(np.mean(sampling_variance)),
        np.finfo(float).tiny,
    )
    _, covariance_exponent = math.frexp(float(np.max(evolutionary_variance)))
    evolutionary_variance = np.ldexp(evolutionary_variance, -covariance_exponent)
    evolutionary_scale = finite_mean(evolutionary_variance)
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
        precision_factor = factor_sparse_positive_definite(posterior_precision)
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
    factor = factor_sparse_positive_definite(posterior_precision)
    posterior_marginal_variance = sparse_precision_update_diagonal(
        posterior_loading,
        posterior_precision,
        solver=factor,
    )
    mean_posterior_variance = float(np.mean(posterior_marginal_variance))
    boundary = bool(rate <= lower_rate * 10.0 or rate >= upper_rate / 10.0)
    return SparseLatentPredictorPosterior(
        mean=np.ldexp(posterior_mean + observed_offset, trait_exponent),
        covariance_model=_restore_sparse_predictor_covariance(
            posterior_model, trait_exponent
        ),
        mean_posterior_variance=math.ldexp(mean_posterior_variance, 2 * trait_exponent),
        prior_mean=math.ldexp(prior_mean + observed_offset, trait_exponent),
        evolutionary_rate=_restore_predictor_rate(
            rate, 2 * trait_exponent - covariance_exponent
        ),
        log_likelihood=-float(objective_value)
        - len(observed) * trait_exponent * math.log(2.0),
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
):
    nuisance_count = observation_loading.shape[1] - evolutionary_model.n_states
    inverse_sampling = 1.0 / sampling
    observation_information = (
        observation_loading.T @ sparse.diags(inverse_sampling) @ observation_loading
    )
    posterior_precision = _sparse_predictor_precision(
        evolutionary_model,
        latent_variance,
        nuisance_count,
        observation_information,
    )
    posterior_factor = factor_sparse_positive_definite(posterior_precision)
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
    return covariance_model, posterior_factor


def _sparse_predictor_reparameterized_covariance(
    evolutionary_model, variance, sampling, sampling_loading
):
    """Use scaled observation residuals as latent coordinates at distinct tips."""
    tip_loading = evolutionary_model.tip_loading
    if not np.all(np.diff(tip_loading.indptr) == 1):
        return None
    indices = tip_loading.indices
    if len(np.unique(indices)) != len(indices) or np.any(tip_loading.data == 0.0):
        return None
    n_states = evolutionary_model.n_states
    nuisance_count = sampling_loading.shape[1]
    observation_sd = np.sqrt(sampling)
    with np.errstate(over="ignore", divide="ignore"):
        inverse_tip_loading = 1.0 / tip_loading.data
        scaled_sd = observation_sd / tip_loading.data
    if not np.all(np.isfinite(inverse_tip_loading) & np.isfinite(scaled_sd)):
        return None
    diagonal = np.ones(n_states, dtype=float)
    diagonal[indices] = scaled_sd
    tip_selector = sparse.csr_matrix(
        (inverse_tip_loading, (indices, np.arange(len(indices)))),
        shape=(n_states, len(indices)),
    )
    transform = sparse.hstack(
        [
            sparse.diags(diagonal, format="csr"),
            -tip_selector @ sampling_loading,
        ],
        format="csr",
    )
    prior_factor = (
        evolutionary_model.precision_factor() @ transform / math.sqrt(variance)
    )
    observation_factor = sparse.eye(n_states, format="csr")[indices]
    independent = sparse.block_diag(
        [observation_factor, sparse.eye(nuisance_count, format="csr")], format="csr"
    )
    sampling_factor = sparse.vstack([prior_factor, independent], format="csr")
    precision = (sampling_factor.T @ sampling_factor).tocsc()
    posterior_factor = factor_sparse_positive_definite(precision)
    observed_loading = sparse.hstack(
        [
            observation_factor.multiply(observation_sd[:, None]),
            -sampling_loading,
        ],
        format="csr",
    )
    return SparseCovarianceModel(
        precision=precision,
        tip_loading=observed_loading,
        logdet_covariance=-posterior_factor.logdet,
        sampling_parent=np.empty(0, dtype=int),
        sampling_transition=np.empty(0),
        sampling_variance=np.empty(0),
        sampling_precision_factor=sampling_factor,
    ), posterior_factor


def _fitted_sparse_predictor_covariance(
    model, variance, sampling, sampling_loading, observation_loading
):
    reparameterized = _sparse_predictor_reparameterized_covariance(
        model, variance, sampling, sampling_loading
    )
    if reparameterized is not None:
        return reparameterized
    return _sparse_predictor_covariance_model(
        model, variance, sampling, observation_loading
    )


def _sparse_predictor_observation_solver(
    evolutionary_model,
    observation_loading,
    sampling,
    nuisance_count,
    observed,
    include_intercept,
):
    n_states = evolutionary_model.n_states
    latent_count = n_states + nuisance_count
    observation_count = len(sampling)
    independent = sparse.block_diag(
        [
            sparse.csc_matrix((n_states, n_states)),
            sparse.eye(nuisance_count, format="csc"),
        ],
        format="csc",
    )
    observation_variance = sparse.diags(-sampling, format="csc")
    coupling = sparse.bmat(
        [
            [independent, observation_loading.T],
            [observation_loading, observation_variance],
        ],
        format="csc",
    )
    prior = sparse.block_diag(
        [
            evolutionary_model.precision,
            sparse.csc_matrix((nuisance_count + observation_count,) * 2),
        ],
        format="csc",
    )
    system = (prior + coupling).tocsc()
    prior_indices = _augmented_prior_entries(system, prior)

    response = (
        np.column_stack((observed, np.ones_like(observed)))
        if include_intercept
        else observed[:, None]
    )
    right = np.vstack((np.zeros((latent_count, response.shape[1])), response))

    def solve(variance):
        # Schur complement is -(R + L Q^-1 L.T). The augmented solve
        # keeps inverse-observation scores free of Woodbury cancellation.
        system.data[prior_indices] = prior.data * (1.0 / variance)
        # The matrix is owned by this synchronous optimizer. Every call
        # replaces all variable entries before factoring the current trial.
        factor = factor_sparse_nonsingular(system)
        solved = factor.solve(right)
        logdet = (
            factor.logdet
            + evolutionary_model.logdet_covariance
            + n_states * math.log(variance)
        )
        return solved[:latent_count], -solved[latent_count:], logdet

    return solve


def _augmented_prior_entries(system, prior):
    # Prior coordinates occupy the top left block and are disjoint from
    # sampling and observation terms. Locate them once in the CSC structure.
    def entries(matrix):
        columns = np.repeat(
            np.arange(matrix.shape[1], dtype=np.int64), np.diff(matrix.indptr)
        )
        return columns * matrix.shape[0] + matrix.indices

    system.sort_indices()
    prior.sum_duplicates()
    prior.eliminate_zeros()
    return np.searchsorted(entries(system), entries(prior))


def _sparse_predictor_reported_rate(variance, covariance_scale, trait_exponent):
    scale, exponent = math.frexp(covariance_scale)
    return _restore_predictor_rate(variance / scale, 2 * trait_exponent - exponent)


def _sparse_predictor_observation_fit(observed, latent, inverse, include_intercept):
    if include_intercept:
        denominator = float(np.sum(inverse[:, 1]))
        if denominator <= 0.0:
            raise np.linalg.LinAlgError(
                "Predictor intercept information is nonpositive."
            )
        prior_mean = float(np.sum(inverse[:, 0])) / denominator
        residual = observed - prior_mean
        inverse_residual = inverse[:, 0] - prior_mean * inverse[:, 1]
        latent_residual = latent[:, 0] - prior_mean * latent[:, 1]
        return prior_mean, float(residual @ inverse_residual), latent_residual
    return 0.0, float(observed @ inverse[:, 0]), latent[:, 0]


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
    observed, sampling, sampling_loading, trait_exponent = _sparse_predictor_units(
        observed, sampling, sampling_loading
    )
    loading = evolutionary_model.tip_loading
    nuisance_count = sampling_loading.shape[1]
    observation_loading = sparse.hstack([loading, sampling_loading], format="csr")
    observed_offset = float(np.mean(observed)) if include_intercept else 0.0
    observed = observed - observed_offset
    observation_solver = _sparse_predictor_observation_solver(
        evolutionary_model,
        observation_loading,
        sampling,
        nuisance_count,
        observed,
        include_intercept,
    )
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
        try:
            latent, inverse, logdet = observation_solver(variance)
            prior_mean, quadratic, conditional_latent = (
                _sparse_predictor_observation_fit(
                    observed, latent, inverse, include_intercept
                )
            )
        except (RuntimeError, np.linalg.LinAlgError):
            return float("inf")
        objective = 0.5 * (len(observed) * math.log(2.0 * math.pi) + logdet + quadratic)
        if not math.isfinite(objective):
            return float("inf")
        if not return_details:
            return objective
        posterior_mean = prior_mean + np.asarray(
            loading @ conditional_latent[: evolutionary_model.n_states]
        ).reshape(-1)
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
            "Sparse predictor evolutionary-rate optimization found no finite fit."
        )
    objective_value, log_variance, message = min(
        finite, key=lambda candidate: candidate[0]
    )
    details = state(log_variance, return_details=True)
    if not isinstance(details, tuple):
        raise ValueError("Sparse predictor conditioning produced an invalid fit.")
    _, prior_mean, posterior_mean = details
    latent_variance = math.exp(log_variance)
    covariance_model, posterior_factor = _fitted_sparse_predictor_covariance(
        evolutionary_model,
        latent_variance,
        sampling,
        sampling_loading,
        observation_loading,
    )
    augmented_loading = covariance_model.tip_loading
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
        mean=np.ldexp(posterior_mean + observed_offset, trait_exponent),
        covariance_model=_restore_sparse_predictor_covariance(
            covariance_model, trait_exponent
        ),
        mean_posterior_variance=math.ldexp(mean_posterior_variance, 2 * trait_exponent),
        prior_mean=math.ldexp(prior_mean + observed_offset, trait_exponent),
        evolutionary_rate=_sparse_predictor_reported_rate(
            latent_variance, evolutionary_model.covariance_scale, trait_exponent
        ),
        log_likelihood=-float(objective_value)
        - len(observed) * trait_exponent * math.log(2.0),
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


def _predictor_error_covariance(
    beta, uncertainty, columns, n_observations, *, continuous_loading=None
):
    """Project scalar or multicolumn predictor uncertainty onto the response."""
    if isinstance(uncertainty, ContinuousPredictorUncertainty):
        return _continuous_predictor_error_covariance(
            beta, uncertainty, columns, continuous_loading
        )
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
        loading = _gmrf_uncertainty_loading(beta, int(columns), uncertainty)
        precision = uncertainty.model.precision
        return materialize_sparse_covariance(precision, loading)
    if isinstance(columns, (int, np.integer)):
        matrix = _positive_semidefinite(uncertainty, "Predictor posterior covariance")
        if matrix.shape != (n_observations, n_observations):
            raise ValueError("Predictor posterior covariance has the wrong dimensions.")
        return scale_covariance_by_slope(matrix, float(beta[int(columns)]))
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


def _continuous_predictor_error_covariance(beta, uncertainty, column, cached):
    if not isinstance(column, (int, np.integer)):
        raise ValueError("Continuous predictor uncertainty requires one column.")
    slope = float(beta[column])
    if cached is not None and math.isfinite(slope):
        loading = slope * cached
        return loading @ loading.T
    loading = continuous_predictor_loading(uncertainty, np.asarray([slope]))
    return (loading @ loading.T).toarray()


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
        loading = _gmrf_uncertainty_loading(beta, int(columns), uncertainty)
        precision = uncertainty.model.precision
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


def _eiv_working_response(response, fixed_covariance, structured):
    fixed = _prepare_eiv_fixed_covariance(fixed_covariance, len(response))
    response, sampling_diagonal, exponent = gaussian_response_units(
        response, covariance_marginal_diagonal(fixed)
    )
    if exponent:
        fixed = rescale_gaussian_covariance(fixed, -exponent)
    return (
        response,
        sampling_diagonal,
        exponent,
        fixed if structured else materialize_covariance(fixed),
    )


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
        scale = finite_mean(positive)
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


def _eiv_working_components(components, component_factors, n_observations, structured):
    normalized, scales = _prepare_eiv_components(
        components, component_factors, n_observations
    )
    if structured:
        return normalized, scales
    for index, (name, matrix) in enumerate(normalized):
        if matrix is None:
            normalized[index] = (
                name,
                materialize_low_rank_covariance(component_factors[name]),
            )
    return normalized, scales


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


def _gmrf_uncertainty_loading(beta, column, uncertainty):
    indices = np.asarray(uncertainty.observation_index, dtype=int)
    model = uncertainty.model
    if (
        indices.ndim != 1
        or np.any(indices < 0)
        or np.any(indices >= model.n_tips)
        or not np.isfinite(float(beta[column]))
    ):
        raise ValueError("GMRF predictor uncertainty is malformed.")
    loading = sparse.csr_matrix(model.tip_loading[indices], dtype=float)
    slope = float(beta[column])
    if uncertainty.row_scale is None:
        return loading * slope
    row_scale = np.asarray(uncertainty.row_scale, dtype=float)
    if row_scale.shape != indices.shape or not np.isfinite(row_scale).all():
        raise ValueError("GMRF predictor row scale is malformed.")
    return loading.multiply((slope * row_scale)[:, None]).tocsr()


def _sparse_precision_marginal_diagonal(loading, precision):
    return sparse_precision_update_diagonal(loading, precision)


def _prepare_eiv_gmrf_components(
    predictor_uncertainties,
    predictor_columns,
    n_coefficients,
    likelihood_groups,
    *,
    use_row_marginals,
    structured,
):
    """Prepare fixed precision factors and exact unit-slope GMRF marginals."""
    if not structured or not any(
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
        loading = _gmrf_uncertainty_loading(unit_beta, column, uncertainty)
        precision = uncertainty.model.precision
        profiles[index] = {
            "loading": loading,
            "precision_factor": uncertainty.model.precision_factor(),
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
    loading = cholesky.loading
    precision = cholesky.prior_precision
    if marginal_diagonal is not None:
        marginal_diagonal = np.asarray(marginal_diagonal, dtype=float)
    elif compute_marginal:
        marginal_diagonal = diagonal + sparse_precision_update_diagonal(
            loading, precision, solver=cholesky.prior_factor
        )
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
    low_rank = (
        sparse.hstack([sparse.csr_matrix(update) for update in updates], format="csr")
        if updates
        else sparse.csr_matrix((n_observations, 0), dtype=float)
    )
    if n_observations <= 64 and low_rank.shape[1] <= 512:
        if np.any(np.asarray(diagonal) <= 0) or not np.isfinite(diagonal).all():
            raise np.linalg.LinAlgError(
                "Covariance diagonal is not positive and finite."
            )
        if not np.isfinite(low_rank.data).all():
            raise ValueError("Low-rank covariance update is non-finite.")
        # At variance boundaries Woodbury subtraction can lose the likelihood
        # curvature. Factor the rectangular covariance root directly for small
        # fits; QR avoids both subtraction and forming an ill-conditioned Gram
        # matrix. Larger fits retain their bounded-memory structured route.
        root = np.column_stack([np.diag(np.sqrt(diagonal)), low_rank.toarray()])
        _, triangular = np.linalg.qr(root.T, mode="reduced")
        signs = np.where(np.diag(triangular) < 0, -1.0, 1.0)
        cholesky = (signs[:, None] * triangular).T
    else:
        cholesky = factor_diagonal_low_rank_updates(diagonal, updates)
    return DiagonalLowRankCovariance(diagonal, low_rank), cholesky


def _dense_eiv_covariance(
    beta,
    variances,
    fixed_covariance,
    normalized_components,
    predictor_uncertainties,
    predictor_columns,
    n_observations,
    prepared_uncertainties,
):
    covariance = fixed_covariance.copy()
    diagonal = np.diag_indices_from(covariance)
    loadings, fixed_covariances = prepared_uncertainties
    for columns, uncertainty, cached, fixed in zip(
        predictor_columns,
        predictor_uncertainties,
        loadings,
        fixed_covariances,
        strict=True,
    ):
        covariance += (
            _predictor_error_covariance(
                beta, uncertainty, columns, n_observations, continuous_loading=cached
            )
            if fixed is None
            else scale_covariance_by_slope(fixed, float(beta[columns]))
        )
    for variance, (_, component) in zip(variances, normalized_components, strict=True):
        assert component is not None
        component_array = np.asarray(component)
        if component_array.ndim == 1:
            covariance[diagonal] += variance * component_array
        else:
            covariance += variance * component_array
    return symmetric_average(covariance)


def _prepare_dense_predictor_uncertainties(uncertainties, columns, size, structured):
    """Validate fixed predictor covariances and prepare reusable loadings."""
    if structured or size > MAX_DENSE_EIV_OBSERVATIONS:
        return (None,) * len(uncertainties), (None,) * len(uncertainties)
    prepared = []
    raw_covariances = []
    for uncertainty, column in zip(uncertainties, columns, strict=True):
        if isinstance(uncertainty, ContinuousPredictorUncertainty):
            if not isinstance(column, (int, np.integer)):
                raise ValueError(
                    "Continuous predictor uncertainty requires one column."
                )
            loading = continuous_predictor_loading(uncertainty, np.ones(1))
            if loading.shape[0] != size:
                raise ValueError(
                    "Predictor posterior covariance has the wrong dimensions."
                )
            prepared.append(
                loading.toarray()
                if loading.shape[1] <= MAX_DENSE_EIV_OBSERVATIONS
                else None
            )
        else:
            prepared.append(None)
        raw_covariances.append(
            _dense_fixed_predictor_covariance(uncertainty, column, size)
        )
    return tuple(prepared), tuple(raw_covariances)


def _dense_fixed_predictor_covariance(uncertainty, column, size):
    if not isinstance(column, (int, np.integer)) or isinstance(
        uncertainty,
        (
            ContinuousPredictorUncertainty,
            GroupedPredictorUncertainty,
            JointPredictorUncertainty,
        ),
    ):
        return None
    if isinstance(uncertainty, GmrfPredictorUncertainty):
        covariance = _predictor_error_covariance(np.ones(1), uncertainty, 0, size)
    else:
        covariance = _positive_semidefinite(
            uncertainty, "Predictor posterior covariance"
        )
    if covariance.shape != (size, size):
        raise ValueError("Predictor posterior covariance has the wrong dimensions.")
    return covariance


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
    gmrf_marginal_profiles,
    likelihood_groups=None,
    dense_predictor_uncertainties=((), ()),
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
                model = uncertainty.model
                profile = gmrf_marginal_profiles[uncertainty_index]
                slope = float(beta[int(columns)])
                precision_updates.append(
                    (
                        profile["loading"] * slope,
                        model.precision,
                        profile["precision_factor"],
                    )
                )
                row_precision_marginal = _add_gmrf_marginal(
                    row_precision_marginal, profile["row"], slope
                )
                grouped_precision_marginal = _add_gmrf_marginal(
                    grouped_precision_marginal, profile["grouped"], slope
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
            marginal_diagonal, grouped_marginal_diagonal = _eiv_precision_marginals(
                diagonal,
                updates,
                row_precision_marginal,
                grouped_precision_marginal,
                likelihood_groups,
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
        dense_predictor_uncertainties,
    )
    return covariance, factor_covariance(covariance)


def _eiv_precision_marginals(diagonal, updates, row, grouped, groups):
    if row is None and grouped is None:
        return None, None
    loading = (
        sparse.hstack(
            [sparse.csr_matrix(update, dtype=float) for update in updates], format="csr"
        )
        if updates
        else sparse.csr_matrix((len(diagonal), 0), dtype=float)
    )
    nonprecision = DiagonalLowRankCovariance(diagonal, loading)
    if row is not None:
        row = covariance_marginal_diagonal(nonprecision) + row
    if grouped is not None:
        grouped = grouped_mean_covariance_diagonal(nonprecision, groups) + grouped
    return row, grouped


def _add_gmrf_marginal(current, unit_marginal, slope):
    if unit_marginal is None:
        return current
    contribution = scale_covariance_by_slope(unit_marginal, slope)
    return contribution if current is None else current + contribution


def _refine_eiv_variance_boundaries(function, result, bounds, num_coefficients):
    """Profile feasible zero-variance boundaries without vanishing log scores.

    A small derivative with respect to log(variance) need not mean that the
    variance is at its optimum. Test its lower bound directly and optimize
    the remaining free parameters before accepting a better boundary fit.
    """
    for index in range(num_coefficients, len(bounds)):
        point = np.asarray(result.x, dtype=float).copy()
        # Only profile a regular zero-variance model. A singular covariance
        # at zero is not made reliable by a tiny positive diagonal floor.
        point[index] = -np.inf
        if not math.isfinite(float(function(point))):
            continue
        point[index] = bounds[index][0]
        boundary_value = float(function(point))
        if not math.isfinite(boundary_value) or boundary_value >= float(result.fun):
            continue
        free = np.arange(len(point)) != index

        def expand(values, point=point, free=free):
            candidate = point.copy()
            candidate[free] = values
            return candidate

        def objective(values, expand=expand):
            return function(expand(values))

        free_bounds = [bound for keep, bound in zip(free, bounds, strict=True) if keep]
        candidate = minimize(
            objective,
            point[free],
            method="L-BFGS-B",
            bounds=free_bounds,
            options={"maxiter": 3000, "ftol": 1e-11},
        )
        if not candidate.success:
            candidate = minimize(
                objective,
                point[free],
                method="Powell",
                bounds=free_bounds,
                options={"maxiter": 10000, "xtol": 1e-8, "ftol": 1e-10},
            )
        if not candidate.success or not math.isfinite(float(candidate.fun)):
            raise ValueError(
                "Errors-in-variables variance-boundary fit did not converge."
            )
        if float(candidate.fun) <= float(result.fun):
            candidate.x = expand(candidate.x)
            result = candidate
    return result


def _eiv_parameter_covariance(function, point, bounds):
    """Condition coefficient information on variance bounds that are active."""
    point = np.asarray(point, dtype=float)
    free = np.asarray(
        [
            value != lower and value != upper
            for value, (lower, upper) in zip(point, bounds, strict=True)
        ]
    )

    def objective(values):
        candidate = point.copy()
        candidate[free] = values
        return function(candidate)

    hessian = _finite_difference_hessian(objective, point[free])
    try:
        free_covariance = np.linalg.inv(hessian)
    except np.linalg.LinAlgError:
        free_covariance = np.linalg.pinv(hessian)
    covariance = np.zeros((len(point), len(point)), dtype=float)
    covariance[np.ix_(free, free)] = free_covariance
    return covariance


def _restore_eiv_response_units(details, exponent, effective_count):
    if not exponent:
        return details
    coefficient_count = len(details["beta"])
    covariance = details["parameter_covariance"]
    exponents = np.zeros(len(covariance), dtype=int)
    exponents[:coefficient_count] = exponent
    try:
        with np.errstate(over="raise", under="ignore"):
            details["parameter_covariance"] = np.ldexp(
                covariance, exponents[:, None] + exponents
            )
    except FloatingPointError as exc:
        raise ValueError(
            "Errors-in-variables coefficient covariance exceeds floating-point range."
        ) from exc
    details = restore_gaussian_response_units(details, exponent, effective_count)
    details["parameters"] = np.concatenate((details["beta"], details["log_variances"]))
    return details


def _eiv_covariance_in_response_units(covariance, factor, exponent, size):
    if not exponent:
        return covariance, factor
    return (
        rescale_gaussian_covariance(covariance, exponent),
        ScaledCovarianceFactor(factor, exponent, size),
    )


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
    effective_likelihood_count, _logdet_weight, likelihood_logdet_offset = (
        effective_likelihood_settings(
            len(response),
            design.shape[1],
            reml,
            likelihood_observations,
            likelihood_logdet_offset,
        )
    )
    if likelihood_groups is not None:
        raise ValueError(
            "Grouped pseudo-likelihood is unsupported; use event-average estimation."
        )
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
    response, sampling_diagonal, response_exponent, fixed_covariance = (
        _eiv_working_response(response, fixed_covariance, structured)
    )
    if len(predictor_uncertainties) != len(predictor_columns):
        raise ValueError("Each uncertain predictor requires coefficient columns.")
    normalized_components, component_scales = _eiv_working_components(
        components, component_factors, len(response), structured
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
        sampling_diagonal,
    )
    lower_variance = max(response_scale * 1e-12, np.finfo(float).tiny)
    upper_variance = max(response_scale * 1e6, lower_variance * 1e6)
    num_coefficients = design.shape[1]
    gmrf_marginal_profiles = _prepare_eiv_gmrf_components(
        predictor_uncertainties,
        predictor_columns,
        num_coefficients,
        likelihood_groups,
        structured=structured,
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
    dense_predictor_uncertainties = _prepare_dense_predictor_uncertainties(
        predictor_uncertainties, predictor_columns, len(response), structured
    )

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
                dense_predictor_uncertainties=dense_predictor_uncertainties,
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
        determinant_term = covariance_logdet - likelihood_logdet_offset
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

    starts = [
        np.concatenate(
            [ordinary_beta, np.log(np.repeat(response_scale, len(components)))]
        )
    ]
    if starting_parameters is not None:
        start = np.asarray(starting_parameters, dtype=float).copy()
        start[:num_coefficients] = (
            np.ldexp(start[:num_coefficients], -response_exponent) - beta_offset
        )
        start[num_coefficients:] -= 2 * response_exponent * math.log(2.0)
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
    result = _refine_eiv_variance_boundaries(evaluate, result, bounds, num_coefficients)
    details = evaluate(result.x, return_details=True)
    if not isinstance(details, dict):
        raise ValueError("Errors-in-variables optimization produced an invalid fit.")
    parameter_covariance = _eiv_parameter_covariance(evaluate, result.x, bounds)
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
            "Errors-in-variables coefficient information is not positive definite; "
            "check whether the fitted latent-predictor variance is near zero."
        )
    if minimum_eigenvalue < 0.0:
        eigenvalues, eigenvectors = np.linalg.eigh(beta_covariance)
        beta_covariance = (eigenvectors * np.maximum(eigenvalues, 0.0)) @ eigenvectors.T
    details["beta_covariance"] = beta_covariance
    details["parameter_covariance"] = parameter_covariance

    fitted_variances = np.exp(details["log_variances"])

    def covariance_for_beta(beta):
        covariance, factor = _conditional_eiv_covariance(
            np.ldexp(np.asarray(beta, dtype=float), -response_exponent),
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
            dense_predictor_uncertainties=dense_predictor_uncertainties,
        )
        return _eiv_covariance_in_response_units(
            covariance, factor, response_exponent, len(response)
        )

    details["covariance_for_beta"] = covariance_for_beta
    details["optimizer_converged"] = bool(result.success)
    details["optimizer_message"] = str(result.message)
    details["reml"] = False
    details["boundary_warning"] = bool(
        np.any(np.exp(details["log_variances"]) <= lower_variance * 10.0)
        or np.any(np.exp(details["log_variances"]) >= upper_variance / 10.0)
    )
    return _restore_eiv_response_units(
        details, response_exponent, effective_likelihood_count
    )


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
