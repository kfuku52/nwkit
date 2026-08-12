"""Latent-predictor Gaussian models for errors-in-variables PGLS."""

import math
import warnings
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize, minimize_scalar


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
    if matrix.size and float(np.linalg.eigvalsh(symmetric).min()) < -tolerance:
        raise ValueError("{} must be positive semidefinite.".format(label))
    symmetric[np.abs(symmetric) < tolerance] = 0.0
    return symmetric


def _solve(cholesky, values):
    return np.linalg.solve(cholesky.T, np.linalg.solve(cholesky, values))


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
        cholesky = np.linalg.cholesky(covariance)
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
        logdet = 2.0 * float(np.log(np.diag(cholesky)).sum())
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
    centered = observed - (float(np.mean(observed)) if include_intercept else 0.0)
    observed_scale = max(
        float(np.mean(centered**2)),
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
        mean=np.asarray(posterior_mean, dtype=float),
        covariance=posterior_covariance,
        prior_mean=float(prior_mean),
        evolutionary_rate=float(rate),
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
):
    """Fit y | x-hat when latent predictor uncertainty depends on beta."""
    response = np.asarray(response, dtype=float)
    design = np.asarray(design, dtype=float)
    fixed_covariance = _positive_semidefinite(fixed_covariance, "Fixed covariance")
    predictor_uncertainties = [
        _positive_semidefinite(matrix, "Predictor posterior covariance")
        for matrix in predictor_uncertainties
    ]
    if len(predictor_uncertainties) != len(predictor_columns):
        raise ValueError("Each uncertain predictor requires one coefficient column.")
    if any(
        matrix.shape != fixed_covariance.shape for matrix in predictor_uncertainties
    ):
        raise ValueError("Predictor posterior covariance has the wrong dimensions.")
    normalized_components = []
    component_scales = []
    for name, matrix in components:
        matrix = _positive_semidefinite(matrix, "Variance component '{}'".format(name))
        positive = np.diag(matrix)
        positive = positive[positive > 0.0]
        if not len(positive):
            raise ValueError("Variance component '{}' has zero diagonal.".format(name))
        scale = float(np.mean(positive))
        normalized_components.append((name, matrix / scale))
        component_scales.append(scale)
    ordinary_beta = np.linalg.lstsq(design, response, rcond=None)[0]
    ordinary_residual = response - design @ ordinary_beta
    response_scale = max(
        float(np.mean(response**2)),
        float(np.mean(ordinary_residual**2)),
        float(np.mean(np.diag(fixed_covariance))),
        np.finfo(float).tiny,
    )
    lower_variance = max(response_scale * 1e-12, np.finfo(float).tiny)
    upper_variance = max(response_scale * 1e6, lower_variance * 1e6)
    num_coefficients = design.shape[1]
    variance_bounds = [(math.log(lower_variance), math.log(upper_variance))] * len(
        normalized_components
    )
    bounds = [(None, None)] * num_coefficients + variance_bounds

    def evaluate(parameters, return_details=False):
        parameters = np.asarray(parameters, dtype=float)
        beta = parameters[:num_coefficients]
        log_variances = parameters[num_coefficients:]
        covariance = fixed_covariance.copy()
        for column, uncertainty in zip(
            predictor_columns, predictor_uncertainties, strict=True
        ):
            covariance += float(beta[column]) ** 2 * uncertainty
        variances = np.exp(log_variances)
        for variance, (_, component) in zip(
            variances, normalized_components, strict=True
        ):
            covariance += variance * component
        covariance = (covariance + covariance.T) / 2.0
        residual = response - design @ beta
        try:
            cholesky = np.linalg.cholesky(covariance)
            inverse_residual = _solve(cholesky, residual)
            quadratic = float(residual @ inverse_residual)
            covariance_logdet = 2.0 * float(np.log(np.diag(cholesky)).sum())
            if reml:
                gram = design.T @ _solve(cholesky, design)
                gram_sign, gram_logdet = np.linalg.slogdet(gram)
                if gram_sign <= 0.0:
                    return float("inf")
            else:
                gram_logdet = 0.0
        except np.linalg.LinAlgError:
            return float("inf")
        effective_n = len(response) - num_coefficients if reml else len(response)
        objective = 0.5 * (
            effective_n * math.log(2.0 * math.pi)
            + covariance_logdet
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
            "parameters": parameters,
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
        starts = [np.asarray(starting_parameters, dtype=float)]
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
    details["optimizer_converged"] = bool(result.success)
    details["optimizer_message"] = str(result.message)
    details["boundary_warning"] = bool(
        np.any(np.exp(details["log_variances"]) <= lower_variance * 10.0)
        or np.any(np.exp(details["log_variances"]) >= upper_variance / 10.0)
    )
    return details
