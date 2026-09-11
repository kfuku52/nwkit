"""GLS-centered phylogenetic PCA and shared-lambda matrix-normal fitting."""

from dataclasses import dataclass

import numpy as np
from scipy.linalg import cho_solve, solve_triangular

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
    logdet_tree = float(2 * np.log(np.diag(factor)).sum())
    return -0.5 * (n * p * (np.log(2 * np.pi) + 1) + p * logdet_tree + n * logdet_trait)


def _lambda(covariance, values, model, fixed):
    if model not in {"BM", "LAMBDA"}:
        raise ValueError("PCA model must be BM or LAMBDA.")
    if fixed is not None and (
        model != "LAMBDA" or not np.isfinite(fixed) or not 0 <= fixed <= 1
    ):
        raise ValueError("--lambda-value requires LAMBDA and a finite value in [0,1].")
    if model == "BM" or fixed is not None:
        return 1.0 if model == "BM" else float(fixed), False
    diagonal = np.diag(np.diag(covariance))
    shared = covariance - diagonal
    if not np.any(shared):
        raise ValueError(
            "Lambda is unidentifiable on this tree; use BM or fix --lambda-value."
        )
    # Column normalization preserves the ML lambda optimum and avoids treating
    # different measurement units as numerical rank loss in its likelihood.
    differences = _finite(values - values[0], "Trait range")
    scaled = differences / np.max(np.abs(differences), axis=0)
    fit = global_bounded_scalar_minimize(
        lambda value: -_profile_likelihood(diagonal + value * shared, scaled), (0, 1)
    )
    if not fit.success or not np.isfinite(fit.fun):
        raise ValueError("PCA lambda optimization failed.")
    return fit.x, True


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
    lam, estimated = _lambda(covariance, values, model, lambda_value)
    covariance = lam * covariance + (1 - lam) * np.diag(np.diag(covariance))
    factor, center, residual = _centered(covariance, values)
    white = _finite(solve_triangular(factor, residual, lower=True), "Whitened traits")
    n, p = values.shape
    sizes = np.max(np.abs(white), axis=0)
    sd = sizes * np.sqrt(np.sum((white / sizes) ** 2, axis=0) / (n - 1))
    scale = sd if mode == "corr" else np.ones(p)
    if np.any(sd <= 0) or not np.all(np.isfinite(sd)):
        raise ValueError("Trait variances cannot be represented; rescale units.")
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
        likelihood = float(_profile_likelihood(covariance, values))
    except ValueError:
        likelihood = None
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
