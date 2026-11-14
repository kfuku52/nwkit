"""Stable full-matrix OU parameterization and local moment identifiability.

For SPD C,D and skew K, A=(D/2+K) C^-1 satisfies AC+CA'=D.
The Lyapunov theorem gives positive-real-part eigenvalues without restricting
the signs of A's individual entries or requiring real eigenvalues.
"""

import math

import numpy as np
from scipy.linalg import expm

from nwkit.clade_index import LcaIndex
from nwkit.multivariate_gaussian_asr import _decode_cholesky, _initial_cholesky


def parse_full_ou_matrix(value, dimension, label):
    """Parse semicolon-delimited rows and comma-delimited columns."""
    if value is None:
        return None
    try:
        matrix = np.asarray(
            [[float(cell) for cell in row.split(",")] for row in value.split(";")]
        )
    except (TypeError, ValueError, AttributeError) as exc:
        raise ValueError(
            f"{label} requires comma-separated entries and semicolon-separated rows."
        ) from exc
    if matrix.shape != (dimension, dimension) or not np.all(np.isfinite(matrix)):
        raise ValueError(
            f"{label} must be a finite {dimension}-by-{dimension} matrix in trait order."
        )
    return matrix


def full_ou_initial(dimension):
    initial, bounds = _initial_cholesky(dimension)
    skew_count = dimension * (dimension - 1) // 2
    return initial * 2 + [0.0] * skew_count, bounds * 2 + [(None, None)] * skew_count


def decode_full_ou(parameters, dimension, time_scale=1.0):
    """Return stable attraction, diffusion and stationary covariance."""
    covariance, offset = _decode_cholesky(parameters, dimension)
    diffusion, offset = _decode_cholesky(parameters, dimension, offset)
    skew = np.zeros((dimension, dimension))
    for row in range(dimension):
        for column in range(row):
            skew[row, column] = parameters[offset]
            skew[column, row] = -parameters[offset]
            offset += 1
    attraction = np.linalg.solve(covariance, (diffusion / 2 + skew).T).T
    return attraction / time_scale, diffusion / time_scale, covariance


def _distance_to_ancestor(node, ancestor):
    terms = []
    while node is not ancestor:
        terms.append(float(node.dist))
        node = node.up
    return math.fsum(terms)


def observation_moment_design(data, *, limit=4096):
    """Distinct observed covariance entries, capped without claiming full rank.

    A full-column-rank subset establishes local identifiability. A deficient
    subset does not establish nonidentifiability unless the design is complete.
    Known measurement covariance is constant and has zero parameter derivative.
    """
    lca = LcaIndex(data.compiled.tree)
    designs: set[tuple[float, float, int, int]] = set()
    coordinates = tuple(zip(data.node_indices, data.trait_indices, strict=True))
    examined = 0
    for first, (left_index, left_trait) in enumerate(coordinates):
        left = data.compiled.nodes[left_index]
        for right_index, right_trait in coordinates[first:]:
            if examined >= 4 * limit:
                return tuple(sorted(designs)), False
            examined += 1
            right = data.compiled.nodes[right_index]
            ancestor = lca.common_ancestor(left, right)
            designs.add(
                (
                    _distance_to_ancestor(left, ancestor),
                    _distance_to_ancestor(right, ancestor),
                    int(left_trait),
                    int(right_trait),
                )
            )
            if len(designs) > limit:
                return tuple(sorted(designs)[:limit]), False
    return tuple(sorted(designs)), True


def _moment_vector(parameters, dimension, time_scale, design):
    attraction, _, covariance = decode_full_ou(parameters, dimension, time_scale)
    distances = {distance for row in design for distance in row[:2]}
    slopes = {distance: expm(-attraction * distance) for distance in distances}
    return np.array(
        [
            (slopes[left] @ covariance @ slopes[right].T)[first, second]
            for left, right, first, second in design
        ]
    )


def full_ou_identifiability(parameters, dimension, time_scale, design, complete):
    """Finite-difference, column-scaled local covariance Jacobian diagnostic.

    This is not a proof of global identifiability or a confidence interval.
    Tolerances are deliberately conservative; rank-deficient and inconclusive
    fits must not receive ordinary regular-model information criteria.
    """
    parameters = np.asarray(parameters, dtype=float)
    columns = []
    for index, value in enumerate(parameters):
        step = 1e-4 * max(1.0, abs(float(value)))
        delta = np.zeros(len(parameters))
        delta[index] = step
        columns.append(
            (
                _moment_vector(parameters + delta, dimension, time_scale, design)
                - _moment_vector(parameters - delta, dimension, time_scale, design)
            )
            / (2 * step)
        )
    jacobian = np.column_stack(columns)
    norms = np.linalg.norm(jacobian, axis=0)
    jacobian /= np.maximum(norms, 1e-12)
    singular = np.linalg.svd(jacobian, compute_uv=False)
    rank = int(np.sum(singular > max(float(singular[0]), 1.0) * 1e-6))
    if rank == len(parameters):
        status = "local_full_rank"
    else:
        status = "local_rank_deficient" if complete else "inconclusive_design_subset"
    ratio = float(singular[-1] / singular[0]) if singular[0] > 0 else 0.0
    return status, rank, ratio
