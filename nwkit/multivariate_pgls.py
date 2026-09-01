"""Joint multivariate Gaussian phylogenetic regression with missing responses."""

import warnings
from dataclasses import dataclass
from typing import Callable, Mapping

import numpy as np
from scipy import sparse
from scipy.optimize import minimize

from nwkit.gaussian import DiagonalLowRankCovariance, materialize_covariance
from nwkit.sparse_laplace import (
    SparseCovarianceModel,
    factor_sparse_nonsingular,
    factor_sparse_positive_definite,
)

MAX_DENSE_MULTIVARIATE_DIMENSION = 2000
MAX_SPARSE_MULTIVARIATE_TIPS = 5000
MAX_SPARSE_MULTIVARIATE_DIMENSION = 20000


@dataclass(frozen=True)
class SparseFittedCovariance:
    """Unmaterialized covariance ``loading @ precision^-1 @ loading.T``."""

    precision: sparse.csc_matrix
    loading: sparse.csr_matrix

    @property
    def shape(self) -> tuple[int, int]:
        size = int(self.loading.shape[0])
        return size, size

    def materialize(self) -> np.ndarray:
        """Materialize the covariance when a caller explicitly requests it."""
        factor = factor_sparse_positive_definite(self.precision)
        solved = factor.solve(self.loading.T.toarray())
        covariance = np.asarray(self.loading @ solved, dtype=float)
        return (covariance + covariance.T) / 2.0


@dataclass(frozen=True)
class MultivariatePglsFit:
    """Profiled multivariate Gaussian fit."""

    coefficients: np.ndarray
    coefficient_covariance: np.ndarray
    component_trait_covariances: Mapping[str, np.ndarray]
    fitted_covariance: np.ndarray | SparseFittedCovariance
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


def _validate_psd_covariance(matrix, size, label):
    matrix = np.asarray(matrix, dtype=float)
    if matrix.shape != (size, size) or not np.isfinite(matrix).all():
        raise ValueError("{} has invalid dimensions or values.".format(label))
    scale = max(1.0, float(np.max(np.abs(matrix))))
    tolerance = np.finfo(float).eps * scale * max(1, size) * 100.0
    if float(np.max(np.abs(matrix - matrix.T))) > tolerance:
        raise ValueError("{} must be symmetric.".format(label))
    matrix = (matrix + matrix.T) / 2.0
    if float(np.min(np.linalg.eigvalsh(matrix))) < -tolerance:
        raise ValueError("{} must be positive semidefinite.".format(label))
    return matrix


def _fixed_covariance_is_sparse_capable(fixed_covariance):
    if fixed_covariance is None or isinstance(
        fixed_covariance, DiagonalLowRankCovariance
    ):
        return True
    fixed_array = np.asarray(fixed_covariance)
    return fixed_array.ndim == 1 or np.array_equal(
        fixed_array, np.diag(np.diag(fixed_array))
    )


def _select_multivariate_sparse_path(
    responses, components, fixed_covariance, allow_large_dense
):
    total_size = responses.size
    sparse_capable = all(
        isinstance(component, SparseCovarianceModel)
        or callable(getattr(component, "sparse_model", None))
        for component in components.values()
    )
    use_sparse = total_size > MAX_DENSE_MULTIVARIATE_DIMENSION
    if use_sparse and (
        not sparse_capable or not _fixed_covariance_is_sparse_capable(fixed_covariance)
    ):
        message = (
            "Multivariate PGLS has a dense joint covariance with {} tip-trait "
            "observations; the validated automatic dense range ends at {}."
        ).format(total_size, MAX_DENSE_MULTIVARIATE_DIMENSION)
        if not allow_large_dense:
            raise ValueError(
                message
                + " Pass allow_large_dense=True to attempt the allocation, or use "
                "tree-structured components and diagonal sampling covariance."
            )
        warnings.warn(
            message + " Large dense allocation was explicitly enabled; attempting it.",
            RuntimeWarning,
            stacklevel=4,
        )
        return False
    return use_sparse


def _warn_sparse_multivariate_range(responses):
    if len(responses) > MAX_SPARSE_MULTIVARIATE_TIPS:
        warnings.warn(
            "Sparse multivariate PGLS above {} tips is outside the routine "
            "validation range (received {}); attempting the fit.".format(
                MAX_SPARSE_MULTIVARIATE_TIPS, len(responses)
            ),
            RuntimeWarning,
            stacklevel=4,
        )
    if responses.size > MAX_SPARSE_MULTIVARIATE_DIMENSION:
        warnings.warn(
            "Sparse multivariate PGLS above {} tip-trait cells is outside the "
            "routine validation range (received {}); attempting the fit.".format(
                MAX_SPARSE_MULTIVARIATE_DIMENSION, responses.size
            ),
            RuntimeWarning,
            stacklevel=4,
        )


def _validate_multivariate_fixed_covariance(fixed_covariance, size):
    if fixed_covariance is None:
        return None
    if isinstance(fixed_covariance, DiagonalLowRankCovariance):
        diagonal = np.asarray(fixed_covariance.diagonal, dtype=float)
        loading = fixed_covariance.low_rank
        finite_loading = (
            np.isfinite(loading.data).all()
            if sparse.issparse(loading)
            else np.isfinite(np.asarray(loading, dtype=float)).all()
        )
        if (
            diagonal.shape != (size,)
            or np.any(diagonal < 0.0)
            or not np.isfinite(diagonal).all()
            or loading.shape[0] != size
            or not finite_loading
        ):
            raise ValueError("Multivariate fixed covariance factor is malformed.")
        return fixed_covariance
    fixed_array = np.asarray(fixed_covariance)
    if fixed_array.ndim != 1:
        return _validate_psd_covariance(
            fixed_array, size, "Multivariate fixed covariance"
        )
    if (
        fixed_array.shape != (size,)
        or not np.isfinite(fixed_array).all()
        or np.any(fixed_array < 0.0)
    ):
        raise ValueError("Multivariate fixed covariance diagonal is invalid.")
    return fixed_array.astype(float)


def _validate_multivariate_inputs(
    responses: np.ndarray,
    design: np.ndarray,
    components: Mapping[str, np.ndarray | Callable[[float | None], np.ndarray]],
    fixed_covariance: np.ndarray | DiagonalLowRankCovariance | None,
    allow_large_dense: bool,
):
    responses = np.asarray(responses, dtype=float)
    design = np.asarray(design, dtype=float)
    if responses.ndim != 2 or responses.shape[1] < 1:
        raise ValueError("PGLS requires at least one response trait.")
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
        if callable(component) or isinstance(component, SparseCovarianceModel):
            continue
        _validate_psd_covariance(
            component,
            len(responses),
            "Covariance component '{}'".format(name),
        )
    use_sparse = _select_multivariate_sparse_path(
        responses, components, fixed_covariance, allow_large_dense
    )
    if use_sparse:
        _warn_sparse_multivariate_range(responses)
    fixed_covariance = _validate_multivariate_fixed_covariance(
        fixed_covariance, responses.size
    )
    return responses, design, observed, fixed_covariance, use_sparse


def _sparse_component_model(component, parameter):
    if isinstance(component, SparseCovarianceModel):
        return component
    builder = getattr(component, "sparse_model", None)
    if not callable(builder):
        raise ValueError("Multivariate covariance component is not sparse-capable.")
    model = builder(parameter)
    if not isinstance(model, SparseCovarianceModel):
        raise ValueError("Sparse covariance factory returned an invalid model.")
    return model


def _sparse_multivariate_latent_model(
    covariance_components,
    trait_covariances,
    decoded,
    n_traits,
    fixed_covariance,
):
    precision_blocks = []
    loading_blocks = []
    prior_logdet = 0.0
    trait_identity = sparse.eye(n_traits, format="csr")
    for name, component in covariance_components.items():
        model = _sparse_component_model(component, decoded)
        trait_covariance = model.covariance_scale * np.asarray(
            trait_covariances[name], dtype=float
        )
        trait_sign, trait_logdet = np.linalg.slogdet(trait_covariance)
        if trait_sign <= 0.0:
            raise np.linalg.LinAlgError(
                "Multivariate trait covariance is not positive definite."
            )
        trait_precision = np.linalg.inv(trait_covariance)
        precision_blocks.append(
            sparse.kron(trait_precision, model.precision, format="csc")
        )
        loading_blocks.append(
            sparse.kron(trait_identity, model.tip_loading, format="csr")
        )
        prior_logdet += (
            model.n_states * trait_logdet + n_traits * model.logdet_covariance
        )
    if fixed_covariance is not None:
        if isinstance(fixed_covariance, DiagonalLowRankCovariance):
            diagonal = np.asarray(fixed_covariance.diagonal, dtype=float)
            fixed_loading = sparse.csr_matrix(fixed_covariance.low_rank, dtype=float)
        else:
            fixed = np.asarray(fixed_covariance, dtype=float)
            diagonal = fixed if fixed.ndim == 1 else np.diag(fixed)
            fixed_loading = sparse.csr_matrix((len(diagonal), 0), dtype=float)
        positive = np.flatnonzero(diagonal > 0.0)
        if len(positive):
            precision_blocks.append(sparse.eye(len(positive), format="csc"))
            loading_blocks.append(
                sparse.csr_matrix(
                    (
                        np.sqrt(diagonal[positive]),
                        (positive, np.arange(len(positive), dtype=int)),
                    ),
                    shape=(len(diagonal), len(positive)),
                )
            )
        if fixed_loading.shape[1]:
            precision_blocks.append(sparse.eye(fixed_loading.shape[1], format="csc"))
            loading_blocks.append(fixed_loading)
    precision = sparse.block_diag(precision_blocks, format="csc")
    loading = sparse.hstack(loading_blocks, format="csr")
    return precision, loading, float(prior_logdet)


def _sparse_multivariate_state(
    covariance_components,
    trait_covariances,
    decoded,
    n_traits,
    fixed_covariance,
    observed_indices,
    observed_design,
    observed_response,
    reml,
):
    precision, loading, prior_logdet = _sparse_multivariate_latent_model(
        covariance_components,
        trait_covariances,
        decoded,
        n_traits,
        fixed_covariance,
    )
    observed_loading = loading[observed_indices]
    zero = sparse.csc_matrix((len(observed_indices), len(observed_indices)))
    kkt = sparse.bmat(
        [[precision, observed_loading.T], [observed_loading, zero]],
        format="csc",
    )
    # The covariance KKT matrix is deliberately indefinite; only
    # nonsingularity, not positive definiteness, is required here.
    factor = factor_sparse_nonsingular(kkt)
    n_states = precision.shape[0]

    def inverse(values):
        values = np.asarray(values, dtype=float)
        vector = values.ndim == 1
        if vector:
            values = values[:, None]
        right = np.vstack([np.zeros((n_states, values.shape[1])), values])
        solved = factor.solve(right)
        result = -solved[n_states:]
        return result[:, 0] if vector else result

    inverse_design = inverse(observed_design)
    information = observed_design.T @ inverse_design
    beta_covariance = np.linalg.pinv(information, hermitian=True)
    coefficients = beta_covariance @ (observed_design.T @ inverse(observed_response))
    residual = observed_response - observed_design @ coefficients
    quadratic = float(residual @ inverse(residual))
    logdet_covariance = factor.logdet + prior_logdet
    degrees_of_freedom = len(observed_response) - observed_design.shape[1]
    if reml and degrees_of_freedom <= 0:
        raise ValueError(
            "Multivariate REML requires more observations than coefficients."
        )
    normalizing_count = degrees_of_freedom if reml else len(observed_response)
    objective = 0.5 * (
        quadratic + logdet_covariance + normalizing_count * np.log(2.0 * np.pi)
    )
    if reml:
        information_sign, information_logdet = np.linalg.slogdet(information)
        if information_sign <= 0.0:
            raise np.linalg.LinAlgError(
                "Multivariate fixed-effect information is singular."
            )
        objective += 0.5 * information_logdet
    return (
        objective,
        coefficients,
        beta_covariance,
        SparseFittedCovariance(precision=precision, loading=loading),
    )


def fit_multivariate_pgls(
    responses: np.ndarray,
    design: np.ndarray,
    covariance_components: Mapping[
        str,
        np.ndarray | SparseCovarianceModel | Callable[[float | None], np.ndarray],
    ],
    *,
    fixed_covariance: np.ndarray | DiagonalLowRankCovariance | None = None,
    evolution_parameter: float | None = None,
    evolution_parameter_bounds: tuple[float, float] | None = None,
    evolution_parameter_decoder: Callable[[float], float] = float,
    evolution_parameter_initial: float | None = None,
    reml: bool = True,
    allow_large_dense: bool = False,
) -> MultivariatePglsFit:
    """Fit separate fixed effects and full trait covariance for each component."""
    if not isinstance(reml, bool):
        raise ValueError("reml must be a boolean.")
    (
        responses,
        design,
        observed,
        fixed_covariance,
        use_sparse,
    ) = _validate_multivariate_inputs(
        responses,
        design,
        covariance_components,
        fixed_covariance,
        allow_large_dense,
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
        if use_sparse:
            return _sparse_multivariate_state(
                covariance_components,
                trait_covariances,
                decoded,
                n_traits,
                fixed_covariance,
                observed_indices,
                observed_design,
                observed_response,
                reml,
            )
        full_covariance = np.zeros((n_tips * n_traits,) * 2, dtype=float)
        for name, component in covariance_components.items():
            if isinstance(component, SparseCovarianceModel):
                tip_covariance = component.materialize()
            else:
                tip_covariance = (
                    component(decoded) if callable(component) else component
                )
            tip_covariance = np.asarray(tip_covariance, dtype=float)
            if (
                tip_covariance.shape != (n_tips, n_tips)
                or not np.isfinite(tip_covariance).all()
            ):
                raise ValueError(
                    "Covariance component '{}' has invalid dimensions or values.".format(
                        name
                    )
                )
            asymmetry = float(np.max(np.abs(tip_covariance - tip_covariance.T)))
            scale = max(1.0, float(np.max(np.abs(tip_covariance))))
            tolerance = np.finfo(float).eps * scale * max(1, n_tips) * 100.0
            if asymmetry > tolerance:
                raise ValueError(
                    "Covariance component '{}' must be symmetric.".format(name)
                )
            tip_covariance = (tip_covariance + tip_covariance.T) / 2.0
            if float(np.min(np.linalg.eigvalsh(tip_covariance))) < -tolerance:
                raise ValueError(
                    "Covariance component '{}' must be positive semidefinite.".format(
                        name
                    )
                )
            full_covariance += np.kron(trait_covariances[name], tip_covariance)
        if fixed_covariance is not None:
            full_covariance += materialize_covariance(fixed_covariance)
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
    if not result.success or not np.isfinite(result.fun):
        fallback_start = result.x if np.isfinite(result.x).all() else initial
        fallback = minimize(
            objective,
            fallback_start,
            method="Powell",
            bounds=bounds,
            options={"maxiter": 10000, "xtol": 1e-7, "ftol": 1e-9},
        )
        if fallback.success and np.isfinite(fallback.fun):
            result = fallback
    if not result.success or not np.isfinite(result.fun):
        raise RuntimeError(
            "Multivariate PGLS optimization failed: {}".format(result.message)
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
