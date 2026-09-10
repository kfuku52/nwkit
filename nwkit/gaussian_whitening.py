"""Square-root Gaussian elimination for many right-hand sides on one tree.

The rotations whiten observed tips without forming their covariance matrix.
They can be reused for observations and mean-design columns at fixed covariance.
The root has known mean zero and an explicitly supplied variance (zero is fixed).
"""

import math
from dataclasses import dataclass

import numpy as np

from nwkit.compiled_tree import CompiledTree


def _vector(values, size, label, *, nonnegative=False):
    result = np.asarray(values, dtype=float)
    if result.shape != (size,) or not np.isfinite(result).all():
        raise ValueError(f"{label} must contain {size} finite values.")
    if nonnegative and np.any(result < 0):
        raise ValueError(f"{label} must be nonnegative.")
    return result


@dataclass(frozen=True)
class _WhiteningBatch:
    destinations: np.ndarray
    first: np.ndarray
    second: np.ndarray
    cosine: np.ndarray
    sine: np.ndarray
    residual_rows: np.ndarray
    scales: np.ndarray


def _batches(compiled, rotations, node_scales):
    heights = np.zeros(len(compiled.nodes), dtype=int)
    for i in compiled.postorder:
        if i:
            parent = compiled.parents[i]
            heights[parent] = max(heights[parent], heights[i] + 1)
    records = {}
    for destination, source, cosine, sine, row in rotations:
        if row == -1:
            records[destination] = [source, source, 1.0, 0.0, -1]
        else:
            records[destination][1:] = [source, cosine, sine, row]
    levels: dict[int, list[int]] = {}
    for destination in records:
        levels.setdefault(heights[destination], []).append(destination)
    scales = dict(node_scales)
    result = []
    for _, destinations in sorted(levels.items()):
        rows = [records[d] for d in destinations]
        result.append(
            _WhiteningBatch(
                np.array(destinations),
                np.array([r[0] for r in rows], dtype=int),
                np.array([r[1] for r in rows], dtype=int),
                np.array([r[2] for r in rows]),
                np.array([r[3] for r in rows]),
                np.array([r[4] for r in rows], dtype=int),
                np.array([scales[d] for d in destinations]),
            )
        )
    return tuple(result)


@dataclass(frozen=True)
class TreeWhitening:
    """A covariance factor with O(nodes) storage and O(nodes * columns) apply."""

    num_nodes: int
    observed_indices: tuple[int, ...]
    leaf_scales: np.ndarray
    batches: tuple[_WhiteningBatch, ...]
    log_determinant: float

    @classmethod
    def build(
        cls,
        compiled: CompiledTree,
        observed_indices,
        slopes,
        innovations,
        errors=None,
        *,
        root_variance=0.0,
    ):
        n = len(compiled.nodes)
        if any(len(children) > 2 for children in compiled.children):
            raise ValueError(
                "Square-root tree whitening currently requires at most two children per node."
            )
        indices = tuple(observed_indices)
        if not indices or len(set(indices)) != len(indices):
            raise ValueError("Whitening needs distinct observed tips.")
        leaf_indices = set(compiled.leaf_indices)
        if any(index not in leaf_indices or index == 0 for index in indices):
            raise ValueError("Whitening observations must be non-root tip indices.")
        slopes = _vector(slopes, n, "Transition slopes")
        innovations = _vector(innovations, n, "Innovation variances", nonnegative=True)
        errors = _vector(
            np.zeros(len(indices)) if errors is None else errors,
            len(indices),
            "Observation variances",
            nonnegative=True,
        )
        root_variance = float(root_variance)
        if not math.isfinite(root_variance) or root_variance < 0:
            raise ValueError("Root variance must be finite and nonnegative.")
        variances = innovations[list(indices)] + errors
        if not np.isfinite(variances).all() or np.any(variances <= 0):
            raise ValueError("Whitening requires positive observed terminal variances.")
        leaf_scales = 1 / np.sqrt(variances)
        precision_roots = np.zeros(n)
        precision_roots[list(indices)] = slopes[list(indices)] * leaf_scales
        active = set(indices)
        logdet = float(np.log(variances).sum())
        rotations = []
        node_scales = []
        row = 0
        for index in compiled.postorder:
            children = [child for child in compiled.children[index] if child in active]
            if not children:
                continue
            active.add(index)
            first = children[0]
            precision_roots[index] = precision_roots[first]
            rotations.append((index, first, 0.0, 1.0, -1))
            for child in children[1:]:
                left, right = precision_roots[index], precision_roots[child]
                norm = math.hypot(left, right)
                cosine, sine = (left / norm, right / norm) if norm else (1.0, 0.0)
                rotations.append((index, child, cosine, sine, row))
                precision_roots[index] = norm
                row += 1
            variance = root_variance if index == 0 else innovations[index]
            divisor = math.hypot(1.0, precision_roots[index] * math.sqrt(variance))
            logdet += 2 * math.log(divisor)
            node_scales.append((index, 1 / divisor))
            if index:
                precision_roots[index] *= slopes[index] / divisor
        if row != len(indices) - 1 or not math.isfinite(logdet):
            raise ValueError("Invalid or unrepresentable Gaussian elimination.")
        return cls(
            n, indices, leaf_scales, _batches(compiled, rotations, node_scales), logdet
        )

    def apply(self, values):
        matrix = np.asarray(values, dtype=float)
        vector = matrix.ndim == 1
        if vector:
            matrix = matrix[:, None]
        if matrix.ndim != 2 or matrix.shape[0] != len(self.observed_indices):
            raise ValueError("Whitening rows must match the observed-tip order.")
        if not np.isfinite(matrix).all():
            raise ValueError("Whitening requires finite values.")
        state = np.zeros((self.num_nodes, matrix.shape[1]))
        state[list(self.observed_indices)] = matrix * self.leaf_scales[:, None]
        result = np.empty_like(matrix)
        for batch in self.batches:
            if len(batch.destinations) == 1:
                left, right = state[batch.first[0]], state[batch.second[0]]
                cosine, sine = batch.cosine[0], batch.sine[0]
                state[batch.destinations[0]] = (
                    cosine * left + sine * right
                ) * batch.scales[0]
                row = batch.residual_rows[0]
                if row >= 0:
                    result[row] = sine * left - cosine * right
                continue
            left, right = state[batch.first], state[batch.second]
            state[batch.destinations] = (
                batch.cosine[:, None] * left + batch.sine[:, None] * right
            ) * batch.scales[:, None]
            paired = batch.residual_rows >= 0
            result[batch.residual_rows[paired]] = (
                batch.sine[paired, None] * left[paired]
                - batch.cosine[paired, None] * right[paired]
            )
        result[-1] = state[0]
        if not np.isfinite(result).all():
            raise ValueError("Gaussian whitening overflow; rescale the inputs.")
        return result[:, 0] if vector else result


def tree_gls(factor: TreeWhitening, values, design):
    """Profile ordinary ML mean coefficients with rank-checked QR."""
    y, x = np.asarray(values, dtype=float), np.asarray(design, dtype=float)
    if x.ndim != 2 or y.shape != (len(x),) or not 0 < x.shape[1] < len(y):
        raise ValueError(
            "GLS needs a mean design with positive residual degrees of freedom."
        )
    whitened = factor.apply(np.column_stack((y, x)))
    response, predictors = whitened[:, 0], whitened[:, 1:]
    norms = np.linalg.norm(predictors, axis=0)
    if np.any(norms == 0):
        raise ValueError("The observed mean design is rank deficient.")
    q, r = np.linalg.qr(predictors / norms)
    if np.linalg.matrix_rank(r) < x.shape[1]:
        raise ValueError("The observed mean design is rank deficient.")
    coefficients = np.linalg.solve(r, q.T @ response) / norms
    residual = response - predictors @ coefficients
    quadratic = float(residual @ residual)
    log_likelihood = -0.5 * (
        len(y) * math.log(2 * math.pi) + factor.log_determinant + quadratic
    )
    inverse = np.linalg.solve(r, np.eye(len(r))) / norms[:, None]
    return coefficients, log_likelihood, quadratic, inverse @ inverse.T
