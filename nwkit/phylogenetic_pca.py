"""GLS-centered phylogenetic PCA and shared-lambda matrix-normal fitting."""

from dataclasses import dataclass

import numpy as np
from scipy.linalg import cho_solve, eigh, solve_triangular

from nwkit.optimization import global_bounded_scalar_minimize


@dataclass(frozen=True)
class PhylogeneticPCA:
    center: np.ndarray
    scale: np.ndarray
    rotation: np.ndarray
    loadings: np.ndarray
    eigenvalues: np.ndarray
    scores: np.ndarray
    covariance: np.ndarray
    lambda_value: float
    lambda_estimated: bool
    log_likelihood: float | None
    status: str
    repeated_eigenvalues: bool


def _finite(array, description):
    if not np.all(np.isfinite(array)):
        raise ValueError(
            f"{description} exceeds floating-point range; rescale trait or tree units."
        )
    return array


def _centered(covariance, values):
    factor = np.linalg.cholesky(covariance)
    weight = cho_solve((factor, True), np.ones(len(values)))
    # Subtract a nearby reference before computing the GLS mean.
    differences = _finite(values - values[0], "Trait range")
    delta = weight @ differences / weight.sum()
    residual = differences - delta
    return factor, values[0] + delta, residual


def _profile_likelihood(covariance, values):
    n, p = values.shape
    factor, _, residual = _centered(covariance, values)
    whitened = solve_triangular(factor, residual, lower=True)
    return _likelihood_from_white(factor, whitened, n, p)


def _likelihood_from_white(factor, whitened, n, p):
    """Profile the ML trait covariance without repeating tree whitening."""
    logdet_tree = float(2 * np.log(np.diag(factor)).sum())
    return _matrix_normal_likelihood(logdet_tree, whitened, n, p)


def _matrix_normal_likelihood(logdet_tree, whitened, n, p):
    units = np.max(np.abs(whitened), axis=0)
    if np.any(units == 0):
        raise ValueError("Constant traits have no full-rank Gaussian likelihood.")
    singular = np.linalg.svd(whitened / units, compute_uv=False)
    tolerance = np.finfo(float).eps * max(n, p) * singular[0]
    if len(singular) < p or singular[-1] <= tolerance:
        raise ValueError(
            "Estimating lambda requires a full-rank trait matrix and more tips than traits."
        )
    logdet_trait = float(
        2 * np.log(singular).sum() + 2 * np.log(units).sum() - p * np.log(n)
    )
    return -0.5 * (n * p * (np.log(2 * np.pi) + 1) + p * logdet_tree + n * logdet_trait)


def _spectral_lambda_profile(covariance, diagonal, values):
    """Reuse the covariance spectrum throughout a well-conditioned search."""
    if np.any(diagonal <= 0):
        return None
    sd = np.sqrt(diagonal)
    correlation = np.empty(covariance.shape, dtype=float, order="F")
    np.divide(covariance, sd[:, None], out=correlation)
    correlation /= sd[None, :]
    eigenvalues, vectors = eigh(correlation, driver="evd", overwrite_a=True)
    n, p = values.shape
    # Near singularity retain Cholesky, whose rounding is more accurate.
    accuracy_floor = np.sqrt(np.finfo(float).eps) * n * eigenvalues[-1]
    if eigenvalues[0] <= accuracy_floor:
        return None
    projected_one = vectors.T @ (1 / sd)
    projected_y = vectors.T @ (values / sd[:, None])
    diagonal_logdet = float(2 * np.log(sd).sum())

    def profile(lam):
        variance = (1 - lam) + lam * eigenvalues
        root = np.sqrt(variance)
        design = projected_one / root
        whitened = projected_y / root[:, None]
        whitened -= np.outer(design, (design @ whitened) / (design @ design))
        logdet = diagonal_logdet + float(np.log(variance).sum())
        return -_matrix_normal_likelihood(logdet, whitened, n, p)

    return profile


def _check_singular_lambda_limit(covariance, diagonal, values):
    """Reject a likelihood that diverges as lambda approaches one."""
    sd = np.sqrt(diagonal)
    correlation = covariance / sd[:, None] / sd[None, :]
    eigenvalues, eigenvectors = np.linalg.eigh(correlation)
    n, p = values.shape
    epsilon = np.finfo(float).eps * n
    tolerance = epsilon * eigenvalues[-1]
    if eigenvalues[0] < -tolerance:
        raise ValueError("PCA Brownian covariance must be positive semidefinite.")
    null = eigenvectors[:, eigenvalues <= tolerance]
    if not null.shape[1]:
        return
    # In null directions, C(lambda) has eigenvalues proportional to 1-lambda.
    # Only independent trait residuals in those directions inflate the profiled
    # trait determinant. The mean can absorb the projected intercept direction.
    whitened = values / sd[:, None]
    contrasts = null.T @ whitened
    inverse_sd = 1 / sd
    intercept = null.T @ inverse_sd
    length = np.linalg.norm(intercept)
    if length > epsilon * np.linalg.norm(inverse_sd):
        direction = intercept / length
        contrasts -= np.outer(direction, direction @ contrasts)
    singular = np.linalg.svd(contrasts, compute_uv=False)
    rank = np.count_nonzero(singular > epsilon * np.linalg.norm(whitened))
    if n * rank < p * null.shape[1]:
        raise ValueError("PCA lambda likelihood is unbounded near singular lambda=1.")


def _lambda(covariance, values, model, fixed):
    if model not in {"BM", "LAMBDA"}:
        raise ValueError("PCA model must be BM or LAMBDA.")
    if fixed is not None and (
        model != "LAMBDA" or not np.isfinite(fixed) or not 0 <= fixed <= 1
    ):
        raise ValueError("--lambda-value requires LAMBDA and a finite value in [0,1].")
    if model == "BM" or fixed is not None:
        return 1.0 if model == "BM" else float(fixed), False
    diagonal = np.diag(covariance)
    if np.count_nonzero(covariance) == np.count_nonzero(diagonal):
        raise ValueError(
            "Lambda is unidentifiable on this tree; use BM or fix --lambda-value."
        )
    # Column normalization preserves the ML lambda optimum and avoids treating
    # different measurement units as numerical rank loss in its likelihood.
    differences = _finite(values - values[0], "Trait range")
    scaled = differences / np.max(np.abs(differences), axis=0)
    singular_endpoint = False
    trial = None

    def cholesky_profile(value):
        nonlocal singular_endpoint, trial
        if trial is None:
            trial = np.empty_like(covariance)
        np.multiply(covariance, value, out=trial)
        np.fill_diagonal(trial, diagonal)
        try:
            return -_profile_likelihood(trial, scaled)
        except np.linalg.LinAlgError:
            # A zero-length terminal split can make the BM endpoint singular
            # while the diagonal contribution keeps interior fits well defined.
            singular_endpoint |= value == 1
            return float("inf")

    spectral_profile = None
    if len(values) >= 64:
        spectral_profile = _spectral_lambda_profile(covariance, diagonal, scaled)
    objective = spectral_profile or cholesky_profile
    fit = global_bounded_scalar_minimize(objective, (0, 1))
    if not fit.success or not np.isfinite(fit.fun):
        raise ValueError("PCA lambda optimization failed.")
    if singular_endpoint:
        _check_singular_lambda_limit(covariance, diagonal, scaled)

        # Search the endpoint on a log distance scale so a finite optimum
        # closer to one than the linear search tolerance can still be found.
        def endpoint_objective(log_distance):
            return objective(-np.expm1(log_distance))

        near_one = np.nextafter(1.0, 0.0)
        endpoint = global_bounded_scalar_minimize(
            endpoint_objective, (np.log(1 - near_one), np.log(1 / 16))
        )
        if endpoint.success and endpoint.fun < fit.fun:
            return float(-np.expm1(endpoint.x)), True
    return fit.x, True


def _scaled_covariance(covariance, lam):
    scaled = np.multiply(covariance, lam)
    diagonal = scaled.diagonal() + (1 - lam) * covariance.diagonal()
    # Preserve the sign of zeros from adding the diagonal matrix.
    np.add(scaled, 0.0, out=scaled)
    np.fill_diagonal(scaled, diagonal)
    return scaled


def fit_pca(covariance, values, *, model="BM", mode="cov", lambda_value=None):
    values = np.asarray(values, dtype=float)
    covariance = np.asarray(covariance, dtype=float)
    if values.ndim != 2 or values.shape[0] < 3 or values.shape[1] < 2:
        raise ValueError(
            "Phylogenetic PCA requires at least three complete tips and two traits."
        )
    if mode not in {"cov", "corr"}:
        raise ValueError("PCA mode must be cov or corr.")
    _finite(values, "Traits")
    _finite(covariance, "Tree covariance")
    if covariance.shape != (len(values), len(values)) or not np.allclose(
        covariance, covariance.T
    ):
        raise ValueError("Tree covariance must be symmetric and match trait rows.")
    if np.any(np.all(values == values[0], axis=0)):
        raise ValueError(
            "PCA requires nonconstant trait columns; remove constant columns."
        )
    # Correlation PCA is invariant to positive column rescaling. Keep the
    # intermediate differences and GLS products in range for extreme units.
    exponents = None
    if mode == "corr":
        _, powers = np.frexp(np.max(np.abs(values), axis=0))
        if np.any(np.abs(powers) > 450):
            exponents = np.where(np.abs(powers) > 450, powers, 0)
            values = np.ldexp(values, -exponents)
    lam, estimated = _lambda(covariance, values, model, lambda_value)
    covariance = _scaled_covariance(covariance, lam)
    factor, center, residual = _centered(covariance, values)
    white = _finite(solve_triangular(factor, residual, lower=True), "Whitened traits")
    n, p = values.shape
    sizes = np.max(np.abs(white), axis=0)
    sd = sizes * np.sqrt(np.sum((white / sizes) ** 2, axis=0) / (n - 1))
    scale = sd if mode == "corr" else np.ones(p)
    if np.any(sd <= 0) or not np.all(np.isfinite(sd)):
        raise ValueError("Trait variances cannot be represented; rescale units.")
    unscaled_white = white
    white = white / scale
    _, singular, right = np.linalg.svd(white, full_matrices=False)
    rank = int(np.sum(singular > np.finfo(float).eps * max(n, p) * singular[0]))
    singular, rotation = singular[:rank], right[:rank].T
    eigenvalues = _finite((singular / np.sqrt(n - 1)) ** 2, "PCA eigenvalues")
    if np.any(eigenvalues <= 0):
        raise ValueError("PCA eigenvalues underflow; rescale units.")
    # Resolve the arbitrary sign by making the largest absolute coefficient positive.
    for column in range(rank):
        absolute = np.abs(rotation[:, column])
        pivot = np.flatnonzero(absolute >= absolute.max() * (1 - 1e-12))[0]
        if rotation[pivot, column] < 0:
            rotation[:, column] *= -1
    scores = _finite((residual / scale) @ rotation, "PC scores")
    loadings = _finite(
        rotation * np.sqrt(eigenvalues) / (sd / scale)[:, None], "PC loadings"
    )
    try:
        likelihood = float(_likelihood_from_white(factor, unscaled_white, n, p))
    except ValueError:
        likelihood = None
    if exponents is not None:
        with np.errstate(over="ignore", under="ignore"):
            center = np.ldexp(center, exponents)
            scale = np.ldexp(scale, exponents)
        if not np.isfinite(scale).all() or np.any(scale <= 0):
            raise ValueError("Trait variances cannot be represented; rescale units.")
        if likelihood is not None:
            likelihood -= n * int(exponents.sum()) * np.log(2)
    return PhylogeneticPCA(
        center=_finite(center, "GLS center"),
        scale=scale,
        rotation=rotation,
        loadings=loadings,
        eigenvalues=eigenvalues,
        scores=scores,
        covariance=covariance,
        lambda_value=lam,
        lambda_estimated=estimated,
        log_likelihood=likelihood,
        status="rank_deficient"
        if rank < p
        else "boundary"
        if estimated and lam in (0, 1)
        else "ok",
        repeated_eigenvalues=bool(
            np.any(np.abs(np.diff(eigenvalues)) <= eigenvalues[0] * 1e-8)
        ),
    )
