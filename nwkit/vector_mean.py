"""GLS stationary means by tree whitening and exact-coordinate conditioning."""

import numpy as np
from scipy.linalg import solve_triangular

from nwkit.compiled_tree import CompiledTree
from nwkit.vector_gaussian import condition_vector_tree
from nwkit.vector_whitening import VectorTreeWhitening, _mean_design_qr


def stationary_mean_gls(process, observed, errors):
    """Whiten the response and all stationary-mean design columns together."""
    compiled = CompiledTree.from_tree(process.tree)
    dimension = process.dimension
    count = len(compiled.leaf_indices)
    values = np.zeros((count, dimension))
    mask = np.zeros_like(values, dtype=bool)
    noise = np.zeros((count, dimension, dimension))
    root_variance = process.root_covariance.diagonal()
    for row, index in enumerate(compiled.leaf_indices):
        node = compiled.nodes[index]
        name = str(node.name)
        observations = observed.get(name)
        if observations is None:
            continue
        mask[row] = [value is not None for value in observations]
        if not mask[row].any():
            continue
        values[row, mask[row]] = [value for value in observations if value is not None]
        if errors is not None:
            noise[row] = errors[name]
        if index == 0:
            return _stationary_mean_by_scores(process, observed, errors)
        observed_traits = np.flatnonzero(mask[row])
        selected = np.ix_(observed_traits, observed_traits)
        tip_covariance = (
            process.transitions[node].covariance[selected] + noise[row][selected]
        )
        scale = np.sqrt(root_variance[observed_traits])
        # Whitening nearly coincident observations can magnify roundoff in
        # the shared-mean design. Score pruning eliminates exact observations
        # and avoids explicitly subtracting their large whitened columns.
        normalized = (tip_covariance / scale[:, None]) / scale[None, :]
        scaled_minimum = np.linalg.eigvalsh(normalized)[0]
        if scaled_minimum <= np.finfo(float).eps ** 0.5:
            return _stationary_mean_by_scores(process, observed, errors)
    slopes = np.zeros((len(compiled.nodes), dimension, dimension))
    innovations = np.zeros_like(slopes)
    for index, node in enumerate(compiled.nodes[1:], 1):
        transition = process.transitions[node]
        slopes[index] = transition.slope
        innovations[index] = transition.covariance
    factor = VectorTreeWhitening.build(
        compiled,
        mask,
        slopes,
        innovations,
        noise,
        root_covariance=process.root_covariance,
    )
    response = values[mask]
    design = np.eye(dimension)[np.nonzero(mask)[1]]
    whitened = factor.apply(np.column_stack((response, design)))
    q, r, norms = _mean_design_qr(whitened[:, 1:])
    return solve_triangular(r, q.T @ whitened[:, 0], check_finite=False) / norms


def _stationary_mean_by_scores(process, observed, errors):
    """Solve H' V^-1 H theta = H' V^-1 y using zero-mean tree smoothing.

    Fisher's identity gives the marginal mean score from conditional latent
    means. Applying that linear score to each observed design column supplies
    the GLS information matrix. Missing coordinates and correlated errors keep
    their original observation pattern; zero-length identity edges contribute
    no mean score. Storage remains O(nodes * traits**2).
    """
    dimension = process.dimension
    identity = np.eye(dimension)

    def score(values):
        result = condition_vector_tree(process, values, error_covariances=errors)
        total = np.linalg.solve(process.root_covariance, result.means[0])
        for index, node in enumerate(result.nodes[1:], 1):
            transition = process.transitions[node]
            if not np.any(transition.covariance):
                continue
            delta = (
                result.means[index]
                - transition.slope @ result.means[result.parents[index]]
            )
            total += (identity - transition.slope).T @ np.linalg.solve(
                transition.covariance, delta
            )
        return total

    linear = score(observed)
    columns = []
    for coordinate in range(dimension):
        values = {
            name: None
            if row is None
            else [
                None if value is None else float(index == coordinate)
                for index, value in enumerate(row)
            ]
            for name, row in observed.items()
        }
        columns.append(score(values))
    information = np.column_stack(columns)
    information = (information + information.T) / 2
    # Cholesky rejects an unidentified mean instead of silently regularizing it.
    factor = np.linalg.cholesky(information)
    return np.linalg.solve(factor.T, np.linalg.solve(factor, linear))
