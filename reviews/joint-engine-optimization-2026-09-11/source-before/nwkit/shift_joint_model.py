"""Joint native shift likelihood with diagonal attraction and correlated noise.

The covariance coordinate S is a variance-scaled diffusion matrix: its diagonal
is the marginal process variance at normalized tree height 1. For shared alpha
it equals the tip covariance; for unequal alpha its off-diagonals do not.
Parameterizing S as PSD guarantees a valid diffusion and every branch innovation.
"""

import math
from dataclasses import dataclass

import numpy as np
from scipy.linalg import solve_triangular

from nwkit.gaussian_whitening import TreeWhitening
from nwkit.shift_native_model import covariance_geometry
from nwkit.vector_whitening import VectorTreeWhitening, vector_gls


def integral_decay(rate, duration):
    """Integral of exp(-rate*s) from 0 to duration, including the zero limit."""
    rate, duration = np.broadcast_arrays(
        np.asarray(rate, float), np.asarray(duration, float)
    )
    result = np.empty_like(rate)
    zero = rate == 0
    result[zero] = duration[zero]
    result[~zero] = -np.expm1(-rate[~zero] * duration[~zero]) / rate[~zero]
    return result


def joint_covariance_geometry(tree, alpha, covariance, root_model):
    """Return diagonal slopes, innovations, root, tip, and diffusion covariance."""
    alpha = np.asarray(alpha, float)
    covariance = np.asarray(covariance, float)
    p = len(alpha)
    if alpha.shape != (p,) or np.isnan(alpha).any() or (alpha < 0).any():
        raise ValueError("Joint alpha must be a nonnegative vector.")
    if covariance.shape != (p, p) or not np.isfinite(covariance).all():
        raise ValueError(
            "Joint process covariance must be a finite trait-square matrix."
        )
    if not np.allclose(covariance, covariance.T, atol=1e-13, rtol=1e-10):
        raise ValueError("Joint process covariance must be symmetric.")
    eigenvalues = np.linalg.eigvalsh(covariance)
    if eigenvalues.min() < -1e-12 * max(1.0, np.max(np.abs(eigenvalues))):
        raise ValueError("Joint process covariance must be positive semidefinite.")
    if root_model not in {"OUfixedRoot", "OUrandomRoot"}:
        raise ValueError("Unknown native root model.")
    if np.isinf(alpha).any():
        if not np.isinf(alpha).all():
            if not np.any(covariance - np.diag(np.diag(covariance))):
                entries = [
                    covariance_geometry(
                        tree, float(a), float(covariance[j, j]), root_model
                    )
                    for j, a in enumerate(alpha)
                ]
                slopes = np.column_stack([entry[0] for entry in entries])
                innovations = np.zeros((len(tree.times), p, p))
                for j, entry in enumerate(entries):
                    innovations[:, j, j] = entry[1]
                return (
                    slopes,
                    innovations,
                    np.diag([entry[2] for entry in entries]),
                    covariance.copy(),
                    None,
                )
            raise ValueError(
                "Joint covariance does not define mixed finite/infinite alpha limits."
            )
        innovations = np.broadcast_to(covariance, (len(tree.times), p, p)).copy()
        root = covariance.copy() if root_model == "OUrandomRoot" else np.zeros((p, p))
        return (
            np.zeros((len(tree.times), p)),
            innovations,
            root,
            covariance.copy(),
            None,
        )
    if root_model == "OUrandomRoot" and (alpha == 0).any():
        raise ValueError("Stationary-root OU is undefined at alpha zero.")
    marginal_integrals = (
        1 / (2 * alpha)
        if root_model == "OUrandomRoot"
        else integral_decay(2 * alpha, 1.0)
    )
    scale = np.sqrt(marginal_integrals)
    diffusion = covariance / scale[:, None] / scale[None, :]
    sums = alpha[:, None] + alpha[None, :]
    innovations = diffusion[None] * integral_decay(
        sums[None], tree.times[:, None, None]
    )
    root = diffusion / sums if root_model == "OUrandomRoot" else np.zeros((p, p))
    tip = (
        root.copy()
        if root_model == "OUrandomRoot"
        else diffusion * integral_decay(sums, 1.0)
    )
    slopes = np.exp(-tree.times[:, None] * alpha[None])
    return slopes, innovations, root, tip, diffusion


def original_diffusion_coordinate(tree, alpha, diffusion, root_model):
    """Convert a physical-time diffusion matrix to the internal S coordinate."""
    diffusion = np.asarray(diffusion, float) * tree.height
    integral = (
        1 / (2 * np.asarray(alpha))
        if root_model == "OUrandomRoot"
        else integral_decay(2 * np.asarray(alpha), 1.0)
    )
    return diffusion * np.sqrt(integral[:, None] * integral[None])


@dataclass(frozen=True)
class NativeJointFit:
    alpha_height: np.ndarray
    covariance_coordinate: np.ndarray
    process_tip_covariance: np.ndarray
    diffusion_covariance: np.ndarray | None
    measurement_variance: np.ndarray
    coefficients: np.ndarray
    coefficient_covariance: np.ndarray
    predicted: np.ndarray
    log_likelihood: float
    quadratic: float
    num_observations: int
    root_model: str
    engine: str


def joint_design(data, layout, alpha):
    """Trait-major coefficients and tip-major observed coordinates."""
    p, q = len(data.trait_names), len(layout.groups)
    matrix = np.zeros((len(data.values), p, p * q))
    node_groups = layout.node_groups(data.tree)
    for trait in range(p):
        matrix[:, trait, trait * q : (trait + 1) * q] = layout.design(
            data.tree, float(alpha[trait]), node_groups=node_groups
        )
    return matrix


def joint_factor(data, alpha, covariance, noise, *, root_model):
    slopes, innovations, root, tip, diffusion = joint_covariance_geometry(
        data.tree, alpha, covariance, root_model
    )
    variances = data.variances + noise[None]
    errors = np.zeros((*variances.shape, variances.shape[1]))
    indices = np.arange(variances.shape[1])
    errors[:, indices, indices] = variances
    factor = VectorTreeWhitening.build(
        data.tree.compiled,
        np.isfinite(data.values),
        slopes,
        innovations,
        errors,
        root_covariance=root,
    )
    return factor, tip, diffusion


def evaluate_joint(data, layout, alpha, covariance, noise, *, root_model="OUfixedRoot"):
    """General exact GLS at specified normalized covariance parameters."""
    alpha = np.broadcast_to(np.asarray(alpha, float), (len(data.trait_names),)).copy()
    noise = np.broadcast_to(np.asarray(noise, float), (len(alpha),)).copy()
    if not np.isfinite(noise).all() or (noise < 0).any():
        raise ValueError(
            "Additional measurement variance must be finite and nonnegative."
        )
    factor, tip, diffusion = joint_factor(
        data, alpha, covariance, noise, root_model=root_model
    )
    design = joint_design(data, layout, alpha)
    observed = np.isfinite(data.values)
    coefficients, loglik, quadratic, coefficient_covariance = vector_gls(
        factor, data.values[observed], design[observed]
    )
    p, q = len(alpha), len(layout.groups)
    return NativeJointFit(
        alpha,
        np.asarray(covariance).copy(),
        tip,
        diffusion,
        noise,
        coefficients.reshape(p, q).T,
        coefficient_covariance,
        np.einsum("ijk,k->ij", design, coefficients),
        loglik,
        quadratic,
        int(observed.sum()),
        root_model,
        "vector_tree_pruning",
    )


def evaluate_separable(
    data, layout, alpha, *, root_model="OUfixedRoot", diagonal=False, factor=None
):
    """Profile free covariance exactly for shared alpha, complete error-free data."""
    if not np.isfinite(data.values).all() or np.any(data.variances):
        raise ValueError(
            "Separable covariance profiling requires complete error-free observations."
        )
    n, p = data.values.shape
    design = layout.design(data.tree, float(alpha))
    q = design.shape[1]
    if n - q < (1 if diagonal else p):
        raise ValueError(
            "Insufficient residual rank for a finite full covariance ML estimate."
        )
    if factor is None:
        slopes, innovations, root = covariance_geometry(
            data.tree, alpha, 1.0, root_model
        )
        factor = TreeWhitening.build(
            data.tree.compiled,
            data.tree.compiled.leaf_indices,
            slopes,
            innovations,
            np.zeros(n),
            root_variance=root,
        )
    white = factor.apply(np.column_stack((data.values, design)))
    y, x = white[:, :p], white[:, p:]
    norms = np.linalg.norm(x, axis=0)
    if np.any(norms == 0) or np.linalg.matrix_rank(x / norms) < q:
        raise ValueError("The observed mean design is rank deficient.")
    basis, r = np.linalg.qr(x / norms, mode="reduced")
    coefficients = solve_triangular(r, basis.T @ y, check_finite=False) / norms[:, None]
    residual = y - basis @ (basis.T @ y)
    covariance = residual.T @ residual / n
    if diagonal:
        covariance = np.diag(np.diag(covariance))
    try:
        chol = np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as exc:
        raise ValueError(
            "Residual covariance is singular; no finite covariance ML estimate."
        ) from exc
    inverse = solve_triangular(r, np.eye(q), check_finite=False) / norms[:, None]
    coefficient_covariance = np.kron(covariance, inverse @ inverse.T)
    loglik = -0.5 * (
        n * (p * (math.log(2 * math.pi) + 1) + 2 * np.log(np.diag(chol)).sum())
        + p * factor.log_determinant
    )
    alphas = np.full(p, alpha)
    _, _, _, tip, diffusion = joint_covariance_geometry(
        data.tree, alphas, covariance, root_model
    )
    return NativeJointFit(
        alphas,
        covariance,
        tip,
        diffusion,
        np.zeros(p),
        coefficients,
        coefficient_covariance,
        design @ coefficients,
        float(loglik),
        float(n * p),
        n * p,
        root_model,
        "separable_profile",
    )
