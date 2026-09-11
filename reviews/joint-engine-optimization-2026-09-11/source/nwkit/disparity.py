"""Squared-distance disparity, crown-time DTT and a fitted Brownian null."""

import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass

import numpy as np
from scipy.integrate import trapezoid
from scipy.linalg import cho_solve, solve_triangular
from scipy.sparse import csr_matrix


@dataclass(frozen=True)
class DTTDesign:
    parents: np.ndarray
    children: tuple
    tip_rows: np.ndarray
    lengths: np.ndarray
    times: np.ndarray
    weights: csr_matrix
    clade_counts: np.ndarray


def make_design(tree, names, depths):
    nodes = list(tree.traverse("preorder"))
    index = {node: i for i, node in enumerate(nodes)}
    tip_index = {name: i for i, name in enumerate(names)}
    children = tuple(tuple(index[c] for c in node.children) for node in nodes)
    counts = np.zeros(len(nodes), dtype=int)
    for i in range(len(nodes) - 1, -1, -1):
        counts[i] = sum(counts[c] for c in children[i]) if children[i] else 1
    events = sorted({depths[node] for node in nodes if len(node.children) >= 2})
    # First root value is before its split. Simultaneous splits form one cohort.
    times = np.asarray([0.0, *events, 1.0])
    parents = np.asarray([-1 if node.is_root else index[node.up] for node in nodes])
    starts = np.asarray([0.0 if node.is_root else depths[node.up] for node in nodes])
    ends = np.asarray([depths[node] for node in nodes])
    row_indices, column_indices, weights = [0], [0], [1.0]
    clade_counts = [1]
    for row, time in enumerate(times[1:], 1):
        active = np.flatnonzero(
            (parents >= 0) & (counts > 1) & (starts <= time) & (ends > time)
        )
        clade_counts.append(len(active))
        row_indices.extend([row] * len(active))
        column_indices.extend(active)
        if len(active):
            weights.extend([1 / len(active)] * len(active))
    matrix = csr_matrix(
        (weights, (row_indices, column_indices)), shape=(len(times), len(nodes))
    )
    return DTTDesign(
        parents,
        children,
        np.asarray([tip_index[node.name] if node.is_leaf else -1 for node in nodes]),
        np.asarray([0.0 if node.is_root else node.dist for node in nodes]),
        times,
        matrix,
        np.asarray(clade_counts),
    )


def clade_disparities(design, values):
    """Mean pairwise squared Euclidean distance via stable subtree moments."""
    p = values.shape[1]
    counts = np.zeros(len(design.parents), dtype=int)
    means = np.zeros((len(counts), p))
    sums = np.zeros(len(counts))
    for node in range(len(counts) - 1, -1, -1):
        if design.tip_rows[node] >= 0:
            counts[node] = 1
            means[node] = values[design.tip_rows[node]]
            continue
        for child in design.children[node]:
            total = counts[node] + counts[child]
            delta = means[child] - means[node]
            sums[node] += (
                sums[child]
                + np.dot(delta, delta) * counts[node] * counts[child] / total
            )
            means[node] += delta * (counts[child] / total)
            counts[node] = total
    disparities = np.divide(
        2 * sums, counts - 1, out=np.zeros_like(sums), where=counts > 1
    )
    if not np.isfinite(disparities).all() or disparities[0] <= 0:
        raise ValueError(
            "DTT requires finite, positive total disparity; rescale nonconstant traits."
        )
    return disparities, counts


def dtt_curve(design, values):
    disparities, _ = clade_disparities(design, values)
    return np.asarray(design.weights @ (disparities / disparities[0]))


def transform_traits(values, scale):
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or not np.isfinite(values).all():
        raise ValueError("DTT requires a finite numeric trait matrix.")
    with np.errstate(over="ignore"):
        differences = values - values[0]
    if not np.isfinite(differences).all():
        raise ValueError("DTT trait ranges overflow; rescale input units.")
    units = np.max(np.abs(differences), axis=0)
    if np.any(units == 0):
        raise ValueError(
            "DTT requires nonconstant selected columns; remove constant traits."
        )
    if scale == "standardize":
        normalized = differences / units
        sd = np.std(normalized, axis=0, ddof=1)
        scales = units * sd
        transformed = normalized / sd
    elif scale == "raw":
        scales = np.full(values.shape[1], units.max())
        transformed = differences / units.max()
    else:
        raise ValueError("DTT scale must be raw or standardize.")
    if not np.isfinite(scales).all() or np.any(scales <= 0):
        raise ValueError("DTT trait scales cannot be represented; rescale input units.")
    return transformed, values[0].copy(), scales


def fit_brownian(covariance, values):
    n, p = values.shape
    if n <= p:
        raise ValueError(
            "BM simulations require more tips than traits; use fewer columns or --n-sim 0."
        )
    try:
        factor = np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as exc:
        raise ValueError(
            "BM simulations require nonsingular tip covariance; check zero-length terminal branches."
        ) from exc
    weight = cho_solve((factor, True), np.ones(n))
    center = weight @ values / weight.sum()
    white = solve_triangular(factor, values - center, lower=True)
    units = np.max(np.abs(white), axis=0)
    if not np.isfinite(white).all() or np.any(units <= 0):
        raise ValueError(
            "BM rate fitting is numerically unresolved; rescale the tree or traits."
        )
    _, singular, right = np.linalg.svd(white / units, full_matrices=False)
    if singular[-1] <= np.finfo(float).eps * max(n, p) * singular[0]:
        raise ValueError(
            "BM simulations require linearly independent trait columns; use --n-sim 0 for descriptive DTT."
        )
    rate_factor = units[:, None] * right.T * (singular / np.sqrt(n - 1))
    rate = rate_factor @ rate_factor.T
    if not np.isfinite(rate).all() or np.any(np.diag(rate) <= 0):
        raise ValueError(
            "BM rate covariance cannot be represented; try --scale standardize."
        )
    return center, rate, rate_factor


def _simulate_one(design, factor, seed):
    rng = np.random.default_rng(seed)
    values = np.zeros((len(design.parents), factor.shape[0]))
    innovations = rng.standard_normal((len(values) - 1, factor.shape[0])) @ factor.T
    for node in range(1, len(values)):
        values[node] = (
            values[design.parents[node]]
            + np.sqrt(design.lengths[node]) * innovations[node - 1]
        )
    tips = np.empty((int(np.sum(design.tip_rows >= 0)), factor.shape[0]))
    for node, row in enumerate(design.tip_rows):
        if row >= 0:
            tips[row] = values[node]
    return dtt_curve(design, tips)


def _simulate_chunk(payload):
    design, factor, seeds = payload
    return [(index, _simulate_one(design, factor, seed)) for index, seed in seeds]


def validate_work(design, traits, count, export_simulations=False):
    if count * len(design.times) > 2_000_000:
        raise ValueError(
            "DTT exceeds 2,000,000 retained simulation/time values; reduce --n-sim."
        )
    if export_simulations and count * len(design.times) > 500_000:
        raise ValueError("--simulations-out exceeds 500,000 rows; reduce --n-sim.")
    work = count * (len(design.parents) * traits * (traits + 1) + design.weights.nnz)
    if work > 100_000_000:
        raise ValueError(
            "DTT exceeds 100,000,000 bounded simulation operations; reduce --n-sim, tips or traits."
        )


def simulate_curves(design, factor, count, seed=1, threads=1):
    seeds = list(enumerate(np.random.SeedSequence(seed).spawn(count)))
    if threads == 1 or count <= 1:
        indexed = _simulate_chunk((design, factor, seeds))
    else:
        workers = min(threads, count)
        chunks = [seeds[i::workers] for i in range(workers)]
        with ProcessPoolExecutor(
            max_workers=workers, mp_context=multiprocessing.get_context("spawn")
        ) as pool:
            indexed = [
                item
                for chunk in pool.map(
                    _simulate_chunk, ((design, factor, chunk) for chunk in chunks)
                )
                for item in chunk
            ]
    return np.asarray([curve for _, curve in sorted(indexed, key=lambda item: item[0])])


def curve_area(times, curves, interval=(0.0, 1.0)):
    """Trapezoidal area with exact linear interpolation at interval boundaries."""
    curves = np.asarray(curves)
    lower, upper = interval

    def boundary(time):
        if time == times[-1]:
            return curves[..., -1]
        left = np.searchsorted(times, time, side="right") - 1
        fraction = (time - times[left]) / (times[left + 1] - times[left])
        return curves[..., left] * (1 - fraction) + curves[..., left + 1] * fraction

    keep = (times > lower) & (times < upper)
    x = np.concatenate(([lower], times[keep], [upper]))
    y = np.concatenate(
        (boundary(lower)[..., None], curves[..., keep], boundary(upper)[..., None]),
        axis=-1,
    )
    return trapezoid(y, x=x, axis=-1)
