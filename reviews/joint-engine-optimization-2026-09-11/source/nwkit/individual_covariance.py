"""Joint BM diffusion and individual covariance, using observed coordinates only.

For observations (s,i,t), V_ab = C[s_a,s_b] Sigma[t_a,t_b]
+ 1[(s_a,i_a)==(s_b,i_b)] W[t_a,t_b]. Species intercepts are
latent; one unknown root mean per trait is profiled (ML) or integrated (REML).
"""

from dataclasses import dataclass

import numpy as np
from scipy.linalg import cho_factor, cho_solve, solve_triangular
from scipy.optimize import minimize


@dataclass(frozen=True)
class IndividualFit:
    sigma: np.ndarray
    within: np.ndarray
    root_mean: np.ndarray
    root_covariance: np.ndarray
    log_likelihood: float
    method: str
    fit_status: str
    optimizer_success: bool
    optimizer_starts: int
    optimizer_converged_starts: int
    identifiability_ratio: float
    sigma_eigenvalue_ratio: float
    within_eigenvalue_ratio: float
    scale: np.ndarray
    offset: np.ndarray
    time_scale: float
    normalized_sigma: np.ndarray
    normalized_within: np.ndarray
    profile: "GaussianProfile"


@dataclass(frozen=True)
class GaussianProfile:
    value: float
    beta: np.ndarray
    beta_covariance: np.ndarray
    alpha: np.ndarray
    inverse: np.ndarray
    precision: np.ndarray
    cholesky: np.ndarray
    inverse_design: np.ndarray


def profile_gaussian(values, design, covariance, method):
    """Exact observed-data likelihood; no dropped normalization constants."""
    factor = cho_factor(covariance, lower=True, check_finite=False)
    inverse = cho_solve(factor, np.eye(len(values)), check_finite=False)
    vi_x = cho_solve(factor, design, check_finite=False)
    information = design.T @ vi_x
    root_factor = cho_factor(information, lower=True, check_finite=False)
    beta_covariance = cho_solve(
        root_factor, np.eye(design.shape[1]), check_finite=False
    )
    beta = beta_covariance @ (vi_x.T @ values)
    residual = values - design @ beta
    alpha = cho_solve(factor, residual, check_finite=False)
    logdet = 2 * np.log(np.diag(factor[0])).sum()
    dimension = len(values)
    precision = inverse
    if method == "REML":
        logdet += 2 * np.log(np.diag(root_factor[0])).sum()
        dimension -= design.shape[1]
        precision = inverse - vi_x @ beta_covariance @ vi_x.T
    value = 0.5 * (dimension * np.log(2 * np.pi) + logdet + residual @ alpha)
    return GaussianProfile(
        float(value),
        beta,
        beta_covariance,
        alpha,
        inverse,
        precision,
        np.tril(factor[0]),
        vi_x,
    )


def covariance_components(covariance, species, individuals, traits, sigma, within):
    evolutionary = covariance[np.ix_(species, species)]
    same = individuals[:, None] == individuals[None, :]
    return (
        evolutionary * sigma[np.ix_(traits, traits)]
        + same * within[np.ix_(traits, traits)]
    )


def _resource_check(n, p):
    parameters = p * (p + 1)
    if 8 * (3 * parameters + 20) * n * n > 512 * 1024**2:
        raise ValueError(
            "Joint individual covariance fitting exceeds the 512 MiB dense workspace "
            "budget; reduce the number of observations or traits."
        )


def _identifiability(evolutionary, same, traits, design, full):
    p = design.shape[1]
    projection = np.eye(len(traits)) - design @ np.linalg.solve(
        design.T @ design, design.T
    )
    derivatives = []
    for component, matrix in enumerate((evolutionary, same)):
        for first in range(p):
            for second in range(first, p):
                if component and not full and first != second:
                    continue
                cells = (traits[:, None] == first) & (traits[None, :] == second)
                if first != second:
                    cells |= cells.T
                derivative = projection @ (matrix * cells) @ projection
                norm = np.linalg.norm(derivative)
                if norm <= 1e-12 * np.linalg.norm(matrix * cells):
                    raise ValueError(
                        "Sigma and W are not identifiable from the observed coordinates."
                    )
                derivatives.append((derivative / norm).ravel())
    # Normalizing derivative columns makes the rank check independent of time units.
    singular = np.linalg.svd(np.asarray(derivatives), compute_uv=False)
    ratio = float(singular[-1] / singular[0])
    if ratio < 1e-9:
        raise ValueError(
            "Sigma and W are not separately identifiable after removing root means."
        )
    return ratio


def _replication_check(values, species, traits, individuals, p, full):
    for trait in range(p):
        chunks = [
            values[(species == s) & (traits == trait)] for s in np.unique(species)
        ]
        chunks = [x for x in chunks if len(x) > 1]
        if sum(len(x) - 1 for x in chunks) < 2:
            raise ValueError(
                "Each trait needs at least two within-species replicate degrees of freedom."
            )
        centered = np.concatenate([x - x.mean() for x in chunks])
        if all(np.ptp(chunk) == 0 for chunk in chunks):
            raise ValueError(
                "Zero within-species variation gives an unbounded likelihood; W cannot be estimated."
            )
        if np.max(np.abs(centered)) < 1e-12:
            raise ValueError(
                "Within-species variation is numerically unresolved relative to between-species variation."
            )
    if not full:
        return
    vectors = np.full((int(individuals.max()) + 1, p), np.nan)
    vectors[individuals, traits] = values
    labels = np.full(len(vectors), -1, dtype=int)
    labels[individuals] = species
    for first in range(p):
        for second in range(first + 1, p):
            valid = np.isfinite(vectors[:, first]) & np.isfinite(vectors[:, second])
            counts = [np.sum(valid & (labels == s)) for s in np.unique(species)]
            if sum(max(0, int(n) - 1) for n in counts) < 2:
                raise ValueError(
                    "Full W needs at least two paired within-species replicate degrees of freedom for every trait pair; use diagonal W or add paired data."
                )
    complete = np.all(np.isfinite(vectors), axis=1)
    differences = []
    for s in np.unique(species):
        chunk = vectors[complete & (labels == s)]
        if len(chunk) > 1:
            differences.extend(chunk[1:] - chunk[0])
    if (
        np.all(complete)
        and len(differences)
        and np.linalg.matrix_rank(differences, tol=1e-10) < p
    ):
        raise ValueError(
            "Singular paired within-species contrasts give an unbounded full-W likelihood."
        )


def _unpack(params, p, full):
    matrices, factors = [], []
    cursor = 0
    for diagonal in (False, not full):
        factor = np.zeros((p, p))
        for row in range(p):
            for col in range(row + 1):
                if diagonal and row != col:
                    continue
                factor[row, col] = (
                    np.exp(params[cursor]) if row == col else params[cursor]
                )
                cursor += 1
        factors.append(factor)
        matrices.append(factor @ factor.T)
    return matrices, factors


def _pack(factors, full, *, gradient=False):
    result = []
    for index, factor in enumerate(factors):
        for row in range(len(factor)):
            for col in range(row + 1):
                if index and not full and row != col:
                    continue
                result.append(
                    factor[row, col]
                    if gradient or row != col
                    else np.log(factor[row, col])
                )
    return np.asarray(result)


class _Objective:
    def __init__(self, values, design, evolutionary, same, traits, method, full):
        self.values, self.design = values, design
        self.evolutionary, self.same, self.traits = evolutionary, same, traits
        self.method, self.full = method, full

    def __call__(self, params):
        p = self.design.shape[1]
        (sigma, within), factors = _unpack(params, p, self.full)
        traits = self.traits
        covariance = (
            self.evolutionary * sigma[np.ix_(traits, traits)]
            + self.same * within[np.ix_(traits, traits)]
        )
        try:
            profile = profile_gaussian(
                self.values, self.design, covariance, self.method
            )
        except np.linalg.LinAlgError:
            return 1e100, np.zeros_like(params)
        score = 0.5 * (profile.precision - np.outer(profile.alpha, profile.alpha))
        gradients = []
        for component, factor in zip(
            (self.evolutionary, self.same), factors, strict=True
        ):
            matrix_gradient = self.design.T @ (score * component) @ self.design
            factor_gradient = 2 * matrix_gradient @ factor
            factor_gradient[np.diag_indices(p)] *= np.diag(factor)
            gradients.append(factor_gradient)
        return profile.value, _pack(gradients, self.full, gradient=True)


def _optimize(objective, p, full):
    bounds = []
    for diagonal in (False, not full):
        for row in range(p):
            for col in range(row + 1):
                if not diagonal or row == col:
                    bounds.append((-18, 8) if row == col else (-1000, 1000))

    def scaled_objective(params):
        value, gradient = objective(params)
        return value / len(objective.values), gradient / len(objective.values)

    candidates = []
    for fraction in (0.2, 0.5, 0.8):
        initial = _pack(
            [np.eye(p) * np.sqrt(fraction), np.eye(p) * np.sqrt(1 - fraction)], full
        )
        candidate = minimize(
            scaled_objective,
            initial,
            jac=True,
            bounds=bounds,
            method="SLSQP",
            options={"maxiter": 1200, "ftol": 1e-12},
        )
        if np.isfinite(candidate.fun) and candidate.fun < 1e99 / len(objective.values):
            candidates.append(candidate)
    if not candidates:
        raise ValueError("Joint covariance optimization found no finite likelihood.")
    best = min(candidates, key=lambda item: item.fun)
    if not _optimizer_converged(best, bounds):
        raise ValueError(
            "Joint covariance optimization did not converge; no fitted reconstruction was produced."
        )
    boundary = any(
        abs(x - low) < 1e-4 or abs(x - high) < 1e-4
        for x, (low, high) in zip(best.x, bounds, strict=True)
    )
    return best, sum(_optimizer_converged(x, bounds) for x in candidates), boundary


def _optimizer_converged(candidate, bounds):
    projected = np.array(candidate.jac, copy=True)
    for index, (low, high) in enumerate(bounds):
        if (candidate.x[index] <= low + 1e-6 and projected[index] > 0) or (
            candidate.x[index] >= high - 1e-6 and projected[index] < 0
        ):
            projected[index] = 0
    return bool(
        candidate.success
        and np.isfinite(projected).all()
        and np.max(np.abs(projected)) <= 1e-5
    )


def fit_individual_covariance(
    covariance,
    values,
    species,
    individuals,
    traits,
    *,
    dimension,
    within="full",
    method="REML",
):
    """Fit independent biological individuals nested in species.

    ``individuals`` are globally unique integer indices (the CLI maps compound
    species/individual IDs); arrays contain only observed scalar coordinates.
    All returned covariances and likelihoods are in the supplied units.
    """
    if method not in {"ML", "REML"} or within not in {"full", "diagonal"}:
        raise ValueError("Use ML/REML and full/diagonal within-species covariance.")
    values = np.asarray(values, dtype=float)
    indices = [np.asarray(a) for a in (species, individuals, traits)]
    if any(a.dtype.kind not in "iu" for a in indices):
        raise ValueError("Observation indices must be integers.")
    species, individuals, traits = (np.asarray(a, dtype=int) for a in indices)
    covariance = np.asarray(covariance, dtype=float)
    p = dimension
    _validate_arrays(covariance, values, species, individuals, traits, p)
    _resource_check(len(values), p)
    _, individuals = np.unique(individuals, return_inverse=True)
    offset = np.array([np.min(values[traits == t]) for t in range(p)])
    differences = values - offset[traits]
    scale = np.array([np.max(np.abs(differences[traits == t])) for t in range(p)])
    if (
        not np.isfinite(differences).all()
        or not np.isfinite(scale).all()
        or np.any(scale == 0)
    ):
        raise ValueError(
            "Constant or out-of-range traits cannot identify joint covariance; rescale trait units."
        )
    normalized = differences / scale[traits]
    time_scale = float(np.max(np.diag(covariance)[species]))
    if time_scale <= 0:
        raise ValueError(
            "Evolutionary Sigma is not identifiable at zero-depth observed tips."
        )
    design = np.eye(p)[traits]
    evolutionary = covariance[np.ix_(species, species)] / time_scale
    same = individuals[:, None] == individuals[None, :]
    full = within == "full"
    _replication_check(normalized, species, traits, individuals, p, full)
    ratio = _identifiability(evolutionary, same, traits, design, full)
    objective = _Objective(normalized, design, evolutionary, same, traits, method, full)
    best, converged, boundary = _optimize(objective, p, full)
    (sigma, w), _ = _unpack(best.x, p, full)
    v = evolutionary * sigma[np.ix_(traits, traits)] + same * w[np.ix_(traits, traits)]
    profile = profile_gaussian(normalized, design, v, method)
    sigma_ratio, w_ratio = (
        float(np.linalg.eigvalsh(a)[0] / np.linalg.eigvalsh(a)[-1]) for a in (sigma, w)
    )
    status = (
        "boundary_covariance"
        if boundary
        or min(sigma_ratio, w_ratio) < 1e-6
        or min(np.linalg.eigvalsh(sigma)[0], np.linalg.eigvalsh(w)[0])
        < 1e-8 * np.linalg.eigvalsh(sigma + w)[-1]
        else "ok"
    )
    correction = np.log(scale[traits]).sum() - (
        np.log(scale).sum() if method == "REML" else 0
    )
    restored_sigma = _restore_covariance(sigma, scale, time_scale)
    restored_w = _restore_covariance(w, scale)
    beta_covariance = _restore_covariance(profile.beta_covariance, scale)
    return IndividualFit(
        restored_sigma,
        restored_w,
        offset + scale * profile.beta,
        beta_covariance,
        -profile.value - float(correction),
        method,
        status,
        True,
        3,
        converged,
        ratio,
        sigma_ratio,
        w_ratio,
        scale,
        offset,
        time_scale,
        sigma,
        w,
        profile,
    )


def _validate_arrays(c, values, species, individuals, traits, p):
    if (
        isinstance(p, (bool, np.bool_))
        or not isinstance(p, (int, np.integer))
        or p < 2
        or any(
            a.ndim != 1 or len(a) != len(values)
            for a in (values, species, individuals, traits)
        )
    ):
        raise ValueError(
            "Joint MV-BM requires aligned vectors and at least two traits."
        )
    if c.ndim != 2 or c.shape[0] != c.shape[1] or not np.isfinite(c).all():
        raise ValueError(
            "The evolutionary covariance must be finite, square and symmetric."
        )
    if not np.isfinite(values).all() or len(values) <= p:
        raise ValueError(
            "Observed values must be finite with positive residual degrees of freedom."
        )
    if (
        min(species) < 0
        or max(species) >= len(c)
        or min(individuals) < 0
        or min(traits) < 0
        or max(traits) >= p
    ):
        raise ValueError("Observation indices are outside their declared dimensions.")
    if set(traits) != set(range(p)) or np.max(np.diag(c)) <= 0:
        raise ValueError(
            "Every trait needs observations and the tree needs positive evolutionary variance."
        )
    _validate_covariance(c)
    keys = list(zip(individuals, traits, strict=True))
    if len(set(keys)) != len(keys):
        raise ValueError("An individual may have only one observation per trait.")
    for individual in np.unique(individuals):
        if len(np.unique(species[individuals == individual])) != 1:
            raise ValueError(
                "Each individual index must belong to exactly one species."
            )


def _validate_covariance(c):
    """Check in correlation units so a long unobserved tip cannot mask errors."""
    diagonal = np.diag(c)
    if np.any(diagonal < 0):
        raise ValueError("The evolutionary covariance must be positive semidefinite.")
    zero = diagonal == 0
    if np.any(c[zero] != 0) or np.any(c[:, zero] != 0):
        raise ValueError("Zero-variance tips must have zero covariance with every tip.")
    positive = ~zero
    scale = np.sqrt(diagonal[positive])
    correlation = c[np.ix_(positive, positive)] / scale[:, None] / scale[None, :]
    if not np.isfinite(correlation).all() or not np.allclose(
        correlation, correlation.T, rtol=0, atol=1e-12
    ):
        raise ValueError(
            "The evolutionary covariance must be symmetric in correlation units."
        )
    if np.linalg.eigvalsh(correlation)[0] < -1e-10:
        raise ValueError("The evolutionary covariance must be positive semidefinite.")


def _restore_covariance(matrix, scale, time_scale=1.0):
    """Restore units without intermediate under/overflow or silent zero variance."""
    mantissa, exponent = np.frexp(matrix)
    scale_mantissa, scale_exponent = np.frexp(scale)
    time_mantissa, time_exponent = np.frexp(time_scale)
    product = (
        mantissa * scale_mantissa[:, None] * scale_mantissa[None, :] / time_mantissa
    )
    power = exponent + scale_exponent[:, None] + scale_exponent[None, :] - time_exponent
    with np.errstate(over="ignore", under="ignore"):
        restored = np.ldexp(product, power)
    if not np.isfinite(restored).all():
        raise ValueError(
            "Covariance exceeds floating-point range; rescale trait units."
        )
    if np.any((np.diag(matrix) > 0) & (np.diag(restored) == 0)):
        raise ValueError("Positive variance underflows to zero; rescale trait units.")
    return restored


def conditional_vector(
    fit, traits, cross_evolution, prior_evolution, *, same_individual=None
):
    """Universal-kriging marginal; root-mean uncertainty is included."""
    p = len(fit.scale)
    sigma, within = fit.normalized_sigma, fit.normalized_within
    profile = fit.profile
    cross = (cross_evolution / fit.time_scale)[None, :] * sigma[:, traits]
    prior = (prior_evolution / fit.time_scale) * sigma
    if same_individual is not None:
        cross = cross + same_individual[None, :] * within[:, traits]
        prior = prior + within
    remainder = np.eye(p) - cross @ profile.inverse_design
    mean = profile.beta + cross @ profile.alpha
    whitened_cross = solve_triangular(
        profile.cholesky, cross.T, lower=True, check_finite=False
    )
    variance = (
        prior
        - whitened_cross.T @ whitened_cross
        + remainder @ profile.beta_covariance @ remainder.T
    )
    if not np.isfinite(variance).all() or not np.isfinite(mean).all():
        raise ValueError("Conditional prediction exceeds floating-point range.")
    variance = (variance + variance.T) / 2
    if same_individual is not None:
        exact = np.unique(traits[same_individual])
        variance[exact, :] = 0
        variance[:, exact] = 0
    eigenvalues, eigenvectors = np.linalg.eigh(variance)
    if eigenvalues[0] < -1e-10 * np.max(np.abs(prior)):
        raise ValueError("Conditional covariance lost positive semidefiniteness.")
    if eigenvalues[0] < 0:
        variance = (eigenvectors * np.maximum(eigenvalues, 0)) @ eigenvectors.T
    restored_mean = fit.offset + fit.scale * mean
    restored_variance = _restore_covariance(variance, fit.scale)
    if not np.isfinite(restored_mean).all() or not np.isfinite(restored_variance).all():
        raise ValueError(
            "Conditional prediction exceeds floating-point range; rescale trait units."
        )
    return restored_mean, restored_variance
