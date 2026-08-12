"""Joint multivariate Gaussian phylogenetic regression with missing responses."""

from dataclasses import dataclass
from typing import Callable, Mapping

import numpy as np
from scipy.optimize import minimize


@dataclass(frozen=True)
class MultivariatePglsFit:
    """Profiled multivariate Gaussian fit."""

    coefficients: np.ndarray
    coefficient_covariance: np.ndarray
    component_trait_covariances: Mapping[str, np.ndarray]
    fitted_covariance: np.ndarray
    evolution_parameter: float | None
    evolution_parameter_status: str
    log_likelihood: float
    optimizer_converged: bool
    optimizer_message: str
    boundary_warning: bool
    n_observations: int
    reml: bool


def _positive_cholesky(matrix: np.ndarray) -> np.ndarray:
    symmetric = (matrix + matrix.T) / 2.0
    scale = max(1.0, float(np.max(np.diag(symmetric))))
    for multiplier in (0.0, 1e-12, 1e-10, 1e-8, 1e-6):
        try:
            return np.linalg.cholesky(
                symmetric + np.eye(len(matrix)) * multiplier * scale
            )
        except np.linalg.LinAlgError:
            continue
    raise np.linalg.LinAlgError("Multivariate covariance is not positive definite.")


def _trait_cholesky(parameters: np.ndarray, dimension: int) -> np.ndarray:
    lower = np.zeros((dimension, dimension), dtype=float)
    position = 0
    for row in range(dimension):
        for column in range(row + 1):
            value = float(parameters[position])
            lower[row, column] = np.exp(value) if row == column else value
            position += 1
    return lower


def _initial_trait_parameters(values: np.ndarray) -> np.ndarray:
    dimension = values.shape[1]
    variances = []
    for column in range(dimension):
        observed = values[np.isfinite(values[:, column]), column]
        variance = float(np.var(observed)) if len(observed) > 1 else 1.0
        variances.append(max(variance, 1e-3))
    parameters = []
    for row in range(dimension):
        for column in range(row + 1):
            parameters.append(0.5 * np.log(variances[row]) if row == column else 0.0)
    return np.asarray(parameters, dtype=float)


def _validate_multivariate_inputs(
    responses: np.ndarray,
    design: np.ndarray,
    components: Mapping[str, np.ndarray | Callable[[float | None], np.ndarray]],
    fixed_covariance: np.ndarray | None,
):
    responses = np.asarray(responses, dtype=float)
    design = np.asarray(design, dtype=float)
    if responses.ndim != 2 or responses.shape[1] < 2:
        raise ValueError("Multivariate PGLS requires at least two response traits.")
    if design.ndim != 2 or len(design) != len(responses):
        raise ValueError("Multivariate response and design dimensions differ.")
    if not np.isfinite(design).all():
        raise ValueError("Multivariate PGLS design must be finite.")
    observed = np.isfinite(responses)
    if np.any(observed.sum(axis=0) <= design.shape[1]):
        raise ValueError(
            "Each multivariate response needs more observed tips than coefficients."
        )
    if np.any(observed.sum(axis=1) == 0):
        raise ValueError("Every retained tip needs at least one observed response.")
    if not components:
        raise ValueError("Multivariate PGLS requires a covariance component.")
    for name, component in components.items():
        if callable(component):
            continue
        matrix = np.asarray(component, dtype=float)
        if matrix.shape != (len(responses), len(responses)):
            raise ValueError(
                "Covariance component '{}' has wrong dimensions.".format(name)
            )
    total_size = responses.size
    if fixed_covariance is not None:
        fixed_covariance = np.asarray(fixed_covariance, dtype=float)
        if fixed_covariance.shape != (total_size, total_size):
            raise ValueError("Multivariate fixed covariance has wrong dimensions.")
    return responses, design, observed, fixed_covariance


def fit_multivariate_pgls(
    responses: np.ndarray,
    design: np.ndarray,
    covariance_components: Mapping[
        str, np.ndarray | Callable[[float | None], np.ndarray]
    ],
    *,
    fixed_covariance: np.ndarray | None = None,
    evolution_parameter: float | None = None,
    evolution_parameter_bounds: tuple[float, float] | None = None,
    evolution_parameter_decoder: Callable[[float], float] = float,
    evolution_parameter_initial: float | None = None,
    reml: bool = True,
) -> MultivariatePglsFit:
    """Fit separate fixed effects and full trait covariance for each component."""
    if not isinstance(reml, bool):
        raise ValueError("reml must be a boolean.")
    responses, design, observed, fixed_covariance = _validate_multivariate_inputs(
        responses, design, covariance_components, fixed_covariance
    )
    n_tips, n_traits = responses.shape
    trait_parameter_count = n_traits * (n_traits + 1) // 2
    initial_trait = _initial_trait_parameters(responses)
    initial = np.tile(initial_trait, len(covariance_components))
    bounds: list[tuple[float | None, float | None]] = []
    for _component in covariance_components:
        for row in range(n_traits):
            for column in range(row + 1):
                bounds.append((-12.0, 12.0) if row == column else (None, None))
    estimate_shape = evolution_parameter_bounds is not None
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
    response_vector = responses.reshape(-1, order="F")
    observed_indices = np.flatnonzero(observed.reshape(-1, order="F"))
    observed_response = response_vector[observed_indices]
    full_design = np.kron(np.eye(n_traits), design)
    observed_design = full_design[observed_indices]

    def unpack(parameters: np.ndarray):
        component_covariances = {}
        for index, name in enumerate(covariance_components):
            start = index * trait_parameter_count
            lower = _trait_cholesky(
                parameters[start : start + trait_parameter_count], n_traits
            )
            component_covariances[name] = lower @ lower.T
        decoded = (
            evolution_parameter_decoder(float(parameters[-1]))
            if estimate_shape
            else evolution_parameter
        )
        return component_covariances, decoded

    def state(parameters: np.ndarray):
        trait_covariances, decoded = unpack(parameters)
        full_covariance = np.zeros((n_tips * n_traits,) * 2, dtype=float)
        for name, component in covariance_components.items():
            tip_covariance = component(decoded) if callable(component) else component
            tip_covariance = np.asarray(tip_covariance, dtype=float)
            if tip_covariance.shape != (n_tips, n_tips):
                raise ValueError(
                    "Covariance component '{}' has wrong dimensions.".format(name)
                )
            full_covariance += np.kron(trait_covariances[name], tip_covariance)
        if fixed_covariance is not None:
            full_covariance += fixed_covariance
        observed_covariance = full_covariance[
            np.ix_(observed_indices, observed_indices)
        ]
        cholesky = _positive_cholesky(observed_covariance)
        inverse_design = np.linalg.solve(
            cholesky.T, np.linalg.solve(cholesky, observed_design)
        )
        information = observed_design.T @ inverse_design
        beta_covariance = np.linalg.pinv(information, hermitian=True)
        coefficients = beta_covariance @ (
            observed_design.T
            @ np.linalg.solve(cholesky.T, np.linalg.solve(cholesky, observed_response))
        )
        residual = observed_response - observed_design @ coefficients
        whitened = np.linalg.solve(cholesky, residual)
        degrees_of_freedom = len(observed_response) - observed_design.shape[1]
        if reml and degrees_of_freedom <= 0:
            raise ValueError(
                "Multivariate REML requires more observations than coefficients."
            )
        normalizing_count = degrees_of_freedom if reml else len(observed_response)
        objective = 0.5 * (
            float(whitened @ whitened)
            + 2.0 * float(np.sum(np.log(np.diag(cholesky))))
            + normalizing_count * np.log(2.0 * np.pi)
        )
        if reml:
            information_sign, information_logdet = np.linalg.slogdet(information)
            if information_sign <= 0.0:
                raise np.linalg.LinAlgError(
                    "Multivariate fixed-effect information is singular."
                )
            objective += 0.5 * information_logdet
        return objective, coefficients, beta_covariance, full_covariance

    def objective(parameters: np.ndarray) -> float:
        try:
            value = float(state(parameters)[0])
        except (ValueError, np.linalg.LinAlgError, FloatingPointError, OverflowError):
            return 1e100
        return value if np.isfinite(value) else 1e100

    result = minimize(
        objective,
        initial,
        method="L-BFGS-B",
        bounds=bounds,
        options={"maxiter": 1000, "ftol": 1e-11, "gtol": 1e-6},
    )
    fitted_objective, coefficients, coefficient_covariance, fitted_covariance = state(
        result.x
    )
    component_trait_covariances, decoded = unpack(result.x)
    shape_boundary = False
    if estimate_shape:
        assert evolution_parameter_bounds is not None
        width = evolution_parameter_bounds[1] - evolution_parameter_bounds[0]
        shape_boundary = bool(
            result.x[-1] <= evolution_parameter_bounds[0] + width * 1e-5
            or result.x[-1] >= evolution_parameter_bounds[1] - width * 1e-5
        )
    variance_boundary = any(
        np.min(np.linalg.eigvalsh(covariance)) < 1e-10
        for covariance in component_trait_covariances.values()
    )
    return MultivariatePglsFit(
        coefficients=coefficients.reshape(n_traits, design.shape[1]),
        coefficient_covariance=coefficient_covariance,
        component_trait_covariances=component_trait_covariances,
        fitted_covariance=fitted_covariance,
        evolution_parameter=decoded,
        evolution_parameter_status="estimated" if estimate_shape else "fixed",
        log_likelihood=-float(fitted_objective),
        optimizer_converged=bool(result.success),
        optimizer_message=str(result.message),
        boundary_warning=variance_boundary or shape_boundary,
        n_observations=len(observed_response),
        reml=reml,
    )
