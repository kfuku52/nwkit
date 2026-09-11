"""Square-root Gaussian elimination for many right-hand sides on one tree.

The rotations whiten observed tips without forming their covariance matrix.
They can be reused for observations and mean-design columns at fixed covariance.
The root has known mean zero and an explicitly supplied variance (zero is fixed).
"""

import math
from dataclasses import dataclass
from functools import lru_cache

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


@lru_cache(maxsize=16)
def _whitening_structure(children, postorder, parents, observed):
    """Cache topology only; covariance-dependent rotations are always rebuilt.

    Keys contain immutable topology and the ordered observation mask, never tree
    objects or numeric covariance values. The bounded cache retains no trait data.
    """
    heights = [0] * len(children)
    active = set(observed)
    steps = []
    levels: dict[int, list[tuple[int, int, int, int]]] = {}
    row = 0
    for index in postorder:
        if index:
            parent = parents[index]
            heights[parent] = max(heights[parent], heights[index] + 1)
        sources = [child for child in children[index] if child in active]
        if not sources:
            continue
        active.add(index)
        paired = len(sources) == 2
        step = (index, sources[0], sources[-1], row if paired else -1)
        steps.append(step)
        levels.setdefault(heights[index], []).append(step)
        row += paired
    batches = []
    for _, records in sorted(levels.items()):
        columns = tuple(
            np.array(column, dtype=int) for column in zip(*records, strict=True)
        )
        for column in columns:
            column.flags.writeable = False
        batches.append(columns)
    return tuple(steps), tuple(batches), row


def _level_rotations(structure, precision_roots, slopes, innovations, root_variance):
    """Eliminate independent nodes together when levels contain many nodes."""
    n = len(precision_roots)
    cosines, sines, scales, divisors = np.ones(n), np.zeros(n), np.ones(n), np.ones(n)
    for destinations, first, second, rows in structure:
        left, right = precision_roots[first], precision_roots[second]
        paired = rows >= 0
        # NumPy's hypot can round differently from the scalar implementation.
        # Keep identical rotations so finite-difference optimizers see the
        # same objective values.
        norms = np.fromiter(
            map(math.hypot, left[paired], right[paired]),
            dtype=float,
            count=int(np.count_nonzero(paired)),
        )
        cosines[destinations[paired]] = np.divide(
            left[paired], norms, out=np.ones_like(norms), where=norms != 0
        )
        sines[destinations[paired]] = np.divide(
            right[paired], norms, out=np.zeros_like(norms), where=norms != 0
        )
        roots = left.copy()
        roots[paired] = norms
        variance = np.where(destinations == 0, root_variance, innovations[destinations])
        divisor = np.fromiter(
            (math.hypot(1.0, value) for value in roots * np.sqrt(variance)),
            dtype=float,
            count=len(roots),
        )
        divisors[destinations] = divisor
        scales[destinations] = 1 / divisor
        precision_roots[destinations] = roots
        nonroot = destinations != 0
        precision_roots[destinations[nonroot]] *= (
            slopes[destinations[nonroot]] / divisor[nonroot]
        )
    return cosines, sines, scales, divisors


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
        steps, structure, rows = _whitening_structure(
            compiled.children, compiled.postorder, compiled.parents, indices
        )
        logdet = float(np.log(variances).sum())
        if len(steps) >= 32 * len(structure):
            cosines, sines, scales, divisors = _level_rotations(
                structure, precision_roots, slopes, innovations, root_variance
            )
            # Retain the original determinant accumulation order.
            for index, _, _, _ in steps:
                logdet += 2 * math.log(divisors[index])
        else:
            # Narrow/deep trees do not amortize vectorized per-level setup.
            cosines, sines, scales = np.ones(n), np.zeros(n), np.ones(n)
            for index, first, second, row in steps:
                precision_roots[index] = precision_roots[first]
                if row >= 0:
                    left, right = precision_roots[first], precision_roots[second]
                    norm = math.hypot(left, right)
                    cosines[index], sines[index] = (
                        (left / norm, right / norm) if norm else (1.0, 0.0)
                    )
                    precision_roots[index] = norm
                variance = root_variance if index == 0 else innovations[index]
                divisor = math.hypot(1.0, precision_roots[index] * math.sqrt(variance))
                logdet += 2 * math.log(divisor)
                scales[index] = 1 / divisor
                if index:
                    precision_roots[index] *= slopes[index] / divisor
        if rows != len(indices) - 1 or not math.isfinite(logdet):
            raise ValueError("Invalid or unrepresentable Gaussian elimination.")
        batches = tuple(
            _WhiteningBatch(
                destinations,
                first,
                second,
                cosines[destinations],
                sines[destinations],
                residual_rows,
                scales[destinations],
            )
            for destinations, first, second, residual_rows in structure
        )
        return cls(n, indices, leaf_scales, batches, logdet)

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
