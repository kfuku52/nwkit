"""Square-root Gaussian pruning for joint vector GLS on a tree.

Each subtree carries at most p rows involving its parent state. Orthogonal
elimination emits independent residual rows; branch integration whitens the
remaining rows. No (tips*p)-square covariance or inverse is formed. Factor
storage is O(nodes*p**2). Zero process covariance is allowed when the observed
covariance is positive definite; no jitter is added to singular models.
"""

import math
from dataclasses import dataclass

import numpy as np
from scipy.linalg import solve_triangular


def _cholesky(matrix, label):
    if not np.isfinite(matrix).all() or not np.allclose(
        matrix, matrix.T, rtol=1e-10, atol=1e-13
    ):
        raise ValueError(f"{label} must be finite and symmetric.")
    try:
        return np.linalg.cholesky((matrix + matrix.T) / 2)
    except np.linalg.LinAlgError as exc:
        raise ValueError(
            f"{label} is not positive definite; no jitter is applied."
        ) from exc


@dataclass(frozen=True)
class _Step:
    node: int
    children: tuple
    positions: np.ndarray | None
    rotation: np.ndarray | None
    retained: int
    factor: np.ndarray


@dataclass(frozen=True)
class VectorObservationPlan:
    """Parameter-independent observed subtrees and coordinate gather indices."""

    compiled: object
    mask: np.ndarray
    entries: tuple

    @classmethod
    def build(cls, compiled, mask):
        mask = np.asarray(mask, dtype=bool)
        if (
            mask.ndim != 2
            or mask.shape[0] != len(compiled.leaf_indices)
            or not mask.any()
        ):
            raise ValueError(
                "Vector observation mask must have observed tip coordinates."
            )
        mask = mask.copy()
        mask.setflags(write=False)
        positions = np.full(mask.shape, -1, dtype=int)
        positions[mask] = np.arange(mask.sum())
        leaves = {node: row for row, node in enumerate(compiled.leaf_indices)}
        active: set[int] = set()
        entries: list[tuple] = []
        for node in compiled.postorder:
            row = leaves.get(node)
            if row is not None:
                observed = np.flatnonzero(mask[row])
                if not len(observed):
                    continue
                entries.append((node, (), row, observed, positions[row, observed]))
            else:
                children = tuple(i for i in compiled.children[node] if i in active)
                if not children:
                    continue
                entries.append((node, children, None, None, None))
            active.add(node)
        return cls(compiled, mask, tuple(entries))


@dataclass(frozen=True)
class VectorTreeWhitening:
    steps: tuple
    root: int
    root_factor: np.ndarray
    observation_count: int
    log_determinant: float

    @classmethod
    def build(
        cls, compiled, mask, slopes, innovations, errors, *, root_covariance, plan=None
    ):
        """Use tip-major observed-coordinate ordering, with diagonal slopes.

        mask has shape (tips,p); slopes (nodes,p); innovations (nodes,p,p);
        errors (tips,p,p). Root state has zero mean and the supplied covariance.
        """
        mask = np.asarray(mask, dtype=bool)
        if mask.ndim != 2 or mask.shape[0] != len(compiled.leaf_indices):
            raise ValueError("Vector observation mask must have shape (tips, traits).")
        p = mask.shape[1]
        n = len(compiled.nodes)
        if not p or not mask.any():
            raise ValueError("Vector whitening needs observed coordinates.")
        slopes = np.asarray(slopes, dtype=float)
        innovations = np.asarray(innovations, dtype=float)
        errors = np.asarray(errors, dtype=float)
        root_covariance = np.asarray(root_covariance, dtype=float)
        for value, shape, label in (
            (slopes, (n, p), "Slopes"),
            (innovations, (n, p, p), "Innovations"),
            (errors, (len(mask), p, p), "Observation covariance"),
            (root_covariance, (p, p), "Root covariance"),
        ):
            if value.shape != shape or not np.isfinite(value).all():
                raise ValueError(f"{label} must be finite with shape {shape}.")
        # Validate covariance inputs, including unobserved coordinates. Tiny
        # roundoff-negative eigenvalues are tolerated for PSD validation only;
        # the matrices themselves are never clipped or modified.
        for value, label in (
            (innovations[1:], "Innovations"),
            (errors, "Observation covariance"),
            (root_covariance[None], "Root covariance"),
        ):
            if not np.allclose(
                value, np.swapaxes(value, -1, -2), rtol=1e-10, atol=1e-13
            ):
                raise ValueError(f"{label} must be symmetric.")
            eigenvalues = np.linalg.eigvalsh(value)
            tolerance = 1e-12 * np.maximum(1, np.max(np.abs(eigenvalues), axis=-1))
            if np.any(np.min(eigenvalues, axis=-1) < -tolerance):
                raise ValueError(f"{label} must be positive semidefinite.")

        if plan is None:
            plan = VectorObservationPlan.build(compiled, mask)
        elif plan.compiled is not compiled or not np.array_equal(plan.mask, mask):
            raise ValueError(
                "Vector observation plan belongs to a different tree or mask."
            )
        messages: dict[int, np.ndarray] = {}
        steps = []
        logdet = 0.0
        for node, children, row, observed, positions in plan.entries:
            if row is not None:
                local = np.eye(p)[observed]
                covariance = innovations[node][np.ix_(observed, observed)]
                covariance = covariance + errors[row][np.ix_(observed, observed)]
                factor = _cholesky(covariance, "Observed branch covariance")
                propagated = solve_triangular(
                    factor, local * slopes[node], lower=True, check_finite=False
                )
                step = _Step(node, (), positions, None, len(observed), factor)
            else:
                local = np.concatenate([messages.pop(i) for i in children])
                retained = min(p, len(local))
                rotation = None
                if len(local) > p:
                    rotation, triangular = np.linalg.qr(local, mode="complete")
                    local = triangular[:retained]
                if node == 0:
                    covariance = np.eye(retained) + local @ root_covariance @ local.T
                    factor = _cholesky(covariance, "Integrated root covariance")
                    steps.append(
                        _Step(
                            node, children, None, rotation, retained, np.eye(retained)
                        )
                    )
                    logdet += 2 * np.log(np.diag(factor)).sum()
                    return cls(
                        tuple(steps), node, factor, int(mask.sum()), float(logdet)
                    )
                covariance = np.eye(retained) + local @ innovations[node] @ local.T
                factor = _cholesky(covariance, "Integrated branch covariance")
                propagated = solve_triangular(
                    factor, local * slopes[node], lower=True, check_finite=False
                )
                step = _Step(node, children, None, rotation, retained, factor)
            messages[node] = propagated
            steps.append(step)
            logdet += 2 * np.log(np.diag(step.factor)).sum()
        raise ValueError("No observed subtree reaches the root.")

    def apply(self, values):
        """Whiten one response or multiple response/design columns together."""
        values = np.asarray(values, dtype=float)
        vector = values.ndim == 1
        if vector:
            values = values[:, None]
        if (
            values.ndim != 2
            or values.shape[0] != self.observation_count
            or not np.isfinite(values).all()
        ):
            raise ValueError(
                "Vector whitening input must be finite and observation aligned."
            )
        output = np.empty_like(values)
        messages: dict[int, np.ndarray] = {}
        offset = 0
        for step in self.steps:
            local = (
                values[step.positions]
                if step.positions is not None
                else np.concatenate([messages.pop(i) for i in step.children])
            )
            if step.rotation is not None:
                local = step.rotation.T @ local
                residual = local[step.retained :]
                output[offset : offset + len(residual)] = residual
                offset += len(residual)
                local = local[: step.retained]
            if step.node == self.root:
                output[offset:] = solve_triangular(
                    self.root_factor, local, lower=True, check_finite=False
                )
            else:
                messages[step.node] = solve_triangular(
                    step.factor, local, lower=True, check_finite=False
                )
        return output[:, 0] if vector else output


def vector_gls(factor, response, design):
    """Profile linear means using QR, retaining their full joint covariance."""
    whitened = factor.apply(np.column_stack((response, design)))
    y, x = whitened[:, 0], whitened[:, 1:]
    norms = np.linalg.norm(x, axis=0)
    if np.any(norms == 0) or np.linalg.matrix_rank(x / norms) < x.shape[1]:
        raise ValueError("The observed mean design is rank deficient.")
    if len(y) <= x.shape[1]:
        raise ValueError(
            "Joint GLS needs residual observations beyond mean parameters."
        )
    q, r = np.linalg.qr(x / norms, mode="reduced")
    coefficients = solve_triangular(r, q.T @ y, check_finite=False) / norms
    residual = y - q @ (q.T @ y)
    quadratic = float(residual @ residual)
    inverse = (
        solve_triangular(r, np.eye(len(norms)), check_finite=False) / norms[:, None]
    )
    covariance = inverse @ inverse.T
    loglik = -0.5 * (
        len(y) * math.log(2 * math.pi) + factor.log_determinant + quadratic
    )
    return coefficients, loglik, quadratic, covariance
