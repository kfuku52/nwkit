"""Bounded observed-coordinate OU GLS and analytic profile-likelihood scores.

Only observed coordinates enter the dense covariance. Tree ages and observation
indices are prepared once per dataset. Large inputs retain tree pruning.
"""

from dataclasses import dataclass

import numpy as np
from scipy.linalg import cho_solve, cholesky, solve_triangular

from nwkit.shift_joint_model import (
    NativeJointFit,
    integral_decay,
    joint_covariance_geometry,
    joint_design,
)
from nwkit.vector_whitening import vector_gls

# Bounds cover work arrays as well as the Cholesky matrix. The coordinate cap
# limits cubic work; the byte bound also accounts for the mean design.
DENSE_MAX_COORDINATES = 384
DENSE_MAX_BYTES = 64 * 1024**2


def dense_eligible(data):
    m = int(np.isfinite(data.values).sum())
    return (
        0 < m <= DENSE_MAX_COORDINATES
        and 8 * (16 * m * m + 2 * data.values.size**2) <= DENSE_MAX_BYTES
    )


def log_integral_derivative(rate, duration):
    """d log(integral exp(-rate*t) dt) / d rate, continuous at rate=0.

    A zero duration has zero covariance and is assigned the limiting score 0.
    The series avoids cancellation; the exponential form avoids overflow.
    """
    rate, duration = np.broadcast_arrays(np.asarray(rate, float), duration)
    z = rate * duration
    small = np.abs(z) < 1e-3
    result = np.empty_like(z)
    result[small] = duration[small] * (
        -0.5 + z[small] / 12 - z[small] ** 3 / 720 + z[small] ** 5 / 30240
    )
    result[~small] = (
        duration[~small] * np.exp(-z[~small]) / -np.expm1(-z[~small]) - 1 / rate[~small]
    )
    return result


def design_derivative(tree, layout, alpha):
    """Derivative of the scaled regime design, including rounded tip ages."""
    labels = layout.node_groups(tree)
    q = len(layout.groups)
    if tree.exact_ultrametric:
        derivative = np.zeros((len(tree.leaf_names), q))
        denominator = integral_decay(alpha, 1.0)
        for branch in layout.shifts:
            node = tree.indices_by_branch[branch]
            parent = tree.compiled.parents[node]
            age = tree.remaining_times[parent]
            weight = integral_decay(alpha, age) / denominator
            score = weight * (
                log_integral_derivative(alpha, age)
                - log_integral_derivative(alpha, 1.0)
            )
            first, last = tree.tip_intervals[node]
            if labels[node]:
                derivative[first:last, labels[node]] += score
            if labels[parent]:
                derivative[first:last, labels[parent]] -= score
        return derivative
    weights = integral_decay(alpha, tree.times) / integral_decay(alpha, 1.0)
    dweights = weights * (
        log_integral_derivative(alpha, tree.times) - log_integral_derivative(alpha, 1.0)
    )
    slopes = np.exp(-alpha * tree.times)
    values = np.zeros((len(tree.times), q))
    derivative = values.copy()
    parents = np.asarray(tree.compiled.parents)
    for indices in tree.levels:
        inherited = values[parents[indices], 1:]
        derivative[indices, 1:] = slopes[indices, None] * (
            derivative[parents[indices], 1:] - tree.times[indices, None] * inherited
        )
        values[indices, 1:] = slopes[indices, None] * inherited
        shifted = indices[labels[indices] != 0]
        values[shifted, labels[shifted]] += weights[shifted]
        derivative[shifted, labels[shifted]] += dweights[shifted]
    return derivative[list(tree.compiled.leaf_indices)]


@dataclass(frozen=True)
class _DenseFactor:
    lower: np.ndarray
    log_determinant: float

    def apply(self, values):
        if not np.isfinite(values).all():
            raise ValueError("Dense GLS inputs must be finite.")
        return solve_triangular(self.lower, values, lower=True, check_finite=False)


class DenseJointContext:
    """Reusable geometry; never cache parameter-dependent matrix factorizations."""

    def __init__(self, data, root_model):
        if not dense_eligible(data):
            raise ValueError("Dense joint GLS exceeds its coordinate or memory bound.")
        self.data, self.root_model = data, root_model
        self.mask = np.isfinite(data.values)
        self.tips, self.traits = np.nonzero(self.mask)
        m = len(self.tips)
        self.shared = np.zeros((m, m))
        for node, (first, last) in enumerate(data.tree.tip_intervals):
            positions = np.flatnonzero((self.tips >= first) & (self.tips < last))
            self.shared[np.ix_(positions, positions)] += data.tree.times[node]
        self.distance = self.shared.diagonal()[:, None] - self.shared
        self.same_tip = self.tips[:, None] == self.tips[None, :]
        self.rows, self.columns = self.traits[:, None], self.traits[None, :]

    def covariance(self, alpha, covariance, noise):
        # Retain the reference geometry's parameter validation and reported
        # diffusion/tip coordinates, including stationary-root restrictions.
        _, _, _, tip, diffusion = joint_covariance_geometry(
            self.data.tree, alpha, covariance, self.root_model
        )
        if not np.isfinite(noise).all() or np.any(noise < 0):
            raise ValueError(
                "Additional measurement variance must be finite and nonnegative."
            )
        r, c = self.rows, self.columns
        if np.isinf(alpha).all():
            kernel = self.same_tip.astype(float)
            alpha_score = np.zeros_like(kernel)
        else:
            if np.isinf(alpha).any():
                raise ValueError(
                    "Dense joint GLS requires common infinite alpha limits."
                )
            sums = alpha[r] + alpha[c]
            decay = np.exp(-alpha[r] * self.distance - alpha[c] * self.distance.T)
            if self.root_model == "OUrandomRoot":
                scale = np.sqrt(2 * alpha)
                kernel = decay * scale[r] * scale[c] / sums
                alpha_score = -self.distance + 0.5 / alpha[r] - 1 / sums
            else:
                scale = np.sqrt(integral_decay(2 * alpha, 1.0))
                kernel = decay * integral_decay(sums, self.shared) / scale[r] / scale[c]
                alpha_score = (
                    -self.distance
                    + log_integral_derivative(sums, self.shared)
                    - log_integral_derivative(2 * alpha[r], 1.0)
                )
        process = covariance[r, c] * kernel
        matrix = process.copy()
        matrix.flat[:: len(matrix) + 1] += (self.data.variances + noise)[self.mask]
        if not np.isfinite(matrix).all():
            raise ValueError("Observed joint covariance must be finite.")
        lower = cholesky(matrix, lower=True, check_finite=False)
        return (
            _DenseFactor(lower, float(2 * np.log(lower.diagonal()).sum())),
            tip,
            diffusion,
            kernel,
            process * alpha_score,
        )

    def evaluate(self, layout, alpha, covariance, noise, *, gradient=False):
        p = len(self.data.trait_names)
        alpha = np.broadcast_to(np.asarray(alpha, float), (p,)).copy()
        noise = np.broadcast_to(np.asarray(noise, float), (p,)).copy()
        covariance = np.asarray(covariance, float)
        factor, tip, diffusion, kernel, alpha_rows = self.covariance(
            alpha, covariance, noise
        )
        design = joint_design(self.data, layout, alpha)
        x, y = design[self.mask], self.data.values[self.mask]
        beta, ll, quadratic, beta_covariance = vector_gls(factor, y, x)
        fit = NativeJointFit(
            alpha,
            covariance.copy(),
            tip,
            diffusion,
            noise,
            beta.reshape(p, len(layout.groups)).T,
            beta_covariance,
            np.einsum("ijk,k->ij", design, beta),
            ll,
            quadratic,
            len(y),
            self.root_model,
            "dense_observed_gls",
        )
        if not gradient:
            return fit
        # Envelope theorem: at the GLS mean optimum no derivative of beta is
        # required. d(-log L) = .5 tr((K^-1-u u')dK) - u' (dX) beta.
        inverse = cho_solve((factor.lower, True), np.eye(len(y)), check_finite=False)
        u = cho_solve((factor.lower, True), y - x @ beta, check_finite=False)
        weight = inverse - np.outer(u, u)
        covariance_score = np.zeros((p, p))
        np.add.at(covariance_score, (self.rows, self.columns), 0.5 * weight * kernel)
        noise_score = np.bincount(
            self.traits, weights=0.5 * weight.diagonal(), minlength=p
        )
        alpha_gradient = np.bincount(
            self.traits, weights=np.sum(weight * alpha_rows, axis=1), minlength=p
        )
        if not np.isinf(alpha).all():
            for j in range(p):
                dmean = (
                    design_derivative(self.data.tree, layout, float(alpha[j]))
                    @ fit.coefficients[:, j]
                )
                positions = self.traits == j
                alpha_gradient[j] -= u[positions] @ dmean[self.tips[positions]]
        return fit, (alpha_gradient, covariance_score, noise_score)
