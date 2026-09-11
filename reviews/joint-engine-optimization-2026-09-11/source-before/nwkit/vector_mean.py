"""GLS stationary means from tree posterior scores, without a dense tip matrix."""

import numpy as np

from nwkit.vector_gaussian import condition_vector_tree


def stationary_mean_gls(process, observed, errors):
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
