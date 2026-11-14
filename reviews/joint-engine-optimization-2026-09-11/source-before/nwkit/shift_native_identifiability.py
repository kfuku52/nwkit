"""Conservative covariance-rank diagnostics, not post-selection intervals."""

import math

import numpy as np


def _observed_pair_heights(data, trait):
    tree = data.tree
    counts = np.zeros(len(tree.branch_ids), dtype=int)
    counts[list(tree.compiled.leaf_indices)] = np.isfinite(data.values[:, trait])
    depths = np.zeros(len(counts))
    for i in range(1, len(counts)):
        depths[i] = depths[tree.compiled.parents[i]] + tree.times[i]
    heights = []
    for i in tree.compiled.postorder:
        children = tree.compiled.children[i]
        if children:
            counts[i] = sum(counts[c] for c in children)
            if all(counts[c] for c in children):
                heights.append(depths[i])
    return np.unique(heights)


def _correlation_derivative(alpha, heights, stationary):
    if math.isinf(alpha):
        return np.zeros(len(heights)), np.zeros(len(heights))
    if alpha == 0:
        return heights, np.zeros(len(heights))
    correlation = np.exp(-2 * alpha * (1 - heights))
    if stationary:
        return correlation, -2 * alpha * (1 - heights) * correlation
    correlation *= -np.expm1(-2 * alpha * heights) / -np.expm1(-2 * alpha)
    if alpha < 1e-3:
        log_derivative = (
            -(1 - heights)
            + alpha * (heights**2 - 1) / 3
            - alpha**3 * (heights**4 - 1) / 45
        )
    else:
        # x/expm1(x) is well-conditioned at zero and vanishes at large x.
        def ratio(x):
            x = np.asarray(x, dtype=float)
            result = np.ones_like(x)
            middle = (x > 0) & (x < 700)
            result[middle] = x[middle] / np.expm1(x[middle])
            result[x >= 700] = 0
            return result

        log_derivative = (
            -2 * (1 - heights) + (ratio(2 * alpha * heights) - ratio(2 * alpha)) / alpha
        )
    return correlation, alpha * correlation * log_derivative


def covariance_identifiability(
    data, trait, fit, *, alpha_free, process_free, noise_free
):
    """Check local covariance information; mean information is not assumed.

    Ultrametric diagonal entries are v+eta; each distinct observed MRCA height
    supplies one off-diagonal covariance type. This avoids dense tip matrices.
    Columns are rescaled before the explicit numerical rank tolerance is applied.
    """
    heights = _observed_pair_heights(data, trait)
    correlation, derivative = _correlation_derivative(
        fit.alpha_height, heights, fit.root_model == "OUrandomRoot"
    )
    columns, names = [], []
    if process_free:
        columns.append(np.r_[1.0, correlation])
        names.append("process_variance")
    if noise_free:
        columns.append(np.r_[1.0, np.zeros(len(heights))])
        names.append("measurement_variance")
    finite_alpha = 0 < fit.alpha_height < math.inf
    if alpha_free and finite_alpha:
        columns.append(np.r_[0.0, fit.process_variance * derivative])
        names.append("alpha")
    matrix = np.column_stack(columns) if columns else np.empty((len(heights) + 1, 0))
    norms = np.linalg.norm(matrix, axis=0)
    np.divide(matrix, norms, out=matrix, where=norms > 0)

    def rank(indices):
        return (
            int(np.linalg.matrix_rank(matrix[:, indices], tol=1e-8)) if indices else 0
        )

    all_indices = list(range(len(names)))
    full_rank = rank(all_indices)
    supported = {
        name: full_rank > rank([j for j in all_indices if j != i])
        for i, name in enumerate(names)
    }
    variance_supported = all(
        supported.get(name, True)
        for name in ("process_variance", "measurement_variance")
    )
    alpha_supported = not alpha_free or (finite_alpha and supported.get("alpha", False))
    return {
        "method": "local_covariance_jacobian_rank; mean_information_not_assumed; not_a_global_identifiability_proof",
        "rank_tolerance": 1e-8,
        "rank": full_rank,
        "free_covariance_parameters": names,
        "parameter_supported": supported,
        "variance_decomposition_supported": variance_supported,
        "finite_alpha_supported": alpha_supported,
    }
