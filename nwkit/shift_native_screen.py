"""Deterministic multivariate group-lasso candidate generation.

Whitening uses the fitted null covariance. Screening and penalized coefficients
only propose locations; all retained models are fitted without this penalty.
Every bootstrap search must regenerate this data-dependent candidate set.
"""

import numpy as np

from nwkit.gaussian_whitening import TreeWhitening
from nwkit.shift_native_model import covariance_geometry


def descendant_design(tree):
    """Contiguous preorder tip intervals avoid repeated descendant traversals."""
    n = len(tree.leaf_names)
    matrix = np.zeros((n, len(tree.branch_ids) - 1))
    ranges = {node: (i, i + 1) for i, node in enumerate(tree.compiled.leaf_indices)}
    columns = {
        branch: column
        for column, branch in enumerate(sorted(set(tree.branch_ids) - {0}))
    }
    for index in tree.compiled.postorder:
        if tree.compiled.children[index]:
            intervals = [ranges[child] for child in tree.compiled.children[index]]
            ranges[index] = (intervals[0][0], intervals[-1][1])
        if index:
            first, last = ranges[index]
            matrix[first:last, columns[tree.branch_ids[index]]] = 1
    return matrix, tuple(columns)


def _whitened_matrices(data, null_fit, memory_limit):
    n, p = data.values.shape
    # Retained standardized matrices plus one working tree matrix and raw design.
    estimated = 8 * n * (len(data.tree.branch_ids) - 1) * (p + 4)
    if estimated > memory_limit:
        raise ValueError(
            f"Group-lasso screening needs approximately {estimated} bytes; increase --search-memory-mb or reduce the input."
        )
    design, branches = descendant_design(data.tree)
    matrices, responses = [], []
    for trait, fit in enumerate(null_fit["fits"]):
        mask = np.isfinite(data.values[:, trait])
        observed = tuple(
            i
            for i, keep in zip(data.tree.compiled.leaf_indices, mask, strict=True)
            if keep
        )
        slopes, innovations, root_variance = covariance_geometry(
            data.tree, fit.alpha_height, fit.process_variance, fit.root_model
        )
        factor = TreeWhitening.build(
            data.tree.compiled,
            observed,
            slopes,
            innovations,
            data.variances[mask, trait] + fit.measurement_variance,
            root_variance=root_variance,
        )
        white = factor.apply(
            np.column_stack(
                (np.ones(np.sum(mask)), data.values[mask, trait], design[mask])
            )
        )
        intercept = white[:, 0] / np.linalg.norm(white[:, 0])
        y = white[:, 1] - intercept * (intercept @ white[:, 1])
        x = white[:, 2:] - intercept[:, None] * (intercept @ white[:, 2:])[None, :]
        norms = np.linalg.norm(x, axis=0)
        informative = (
            norms > np.finfo(float).eps * max(1.0, np.linalg.norm(white[:, 2:])) * 100
        )
        x[:, informative] /= norms[informative]
        x[:, ~informative] = 0
        matrices.append(x)
        responses.append(y)
    return matrices, responses, branches


def _loss_gradient(matrices, responses, coefficients):
    residuals = [
        x @ coefficients[:, j] - y
        for j, (x, y) in enumerate(zip(matrices, responses, strict=True))
    ]
    loss = 0.5 * sum(float(r @ r) for r in residuals)
    gradient = np.column_stack(
        [x.T @ r for x, r in zip(matrices, residuals, strict=True)]
    )
    return loss, gradient


def _proximal_step(matrices, responses, coefficients, strength, step):
    loss, gradient = _loss_gradient(matrices, responses, coefficients)
    for _ in range(40):
        trial = coefficients - step * gradient
        norms = np.linalg.norm(trial, axis=1)
        multiplier = np.maximum(
            0, 1 - step * strength / np.maximum(norms, np.finfo(float).tiny)
        )
        trial *= multiplier[:, None]
        delta = trial - coefficients
        trial_loss, trial_gradient = _loss_gradient(matrices, responses, trial)
        bound = (
            loss
            + float(np.sum(gradient * delta))
            + float(np.sum(delta**2)) / (2 * step)
        )
        if trial_loss <= bound + 1e-12 * max(1, loss):
            return trial, trial_gradient, step
        step *= 0.5
    raise ValueError("Group-lasso backtracking could not verify its quadratic bound.")


def _kkt_residual(coefficients, gradient, strength):
    norms = np.linalg.norm(coefficients, axis=1)
    active = norms > 0
    residual = np.maximum(0, np.linalg.norm(gradient, axis=1) - strength)
    if np.any(active):
        residual[active] = np.linalg.norm(
            gradient[active] + strength * coefficients[active] / norms[active, None],
            axis=1,
        )
    return float(np.max(residual))


def group_lasso_screen(
    data, null_fit, *, pool_size=24, iterations=150, paths=6, memory_limit=512 * 1024**2
):
    if pool_size < 1 or iterations < 1 or paths < 1 or memory_limit < 1:
        raise ValueError("Native screening budgets must be positive.")
    matrices, responses, branches = _whitened_matrices(data, null_fit, memory_limit)
    gradients = np.column_stack(
        [x.T @ y for x, y in zip(matrices, responses, strict=True)]
    )
    correlation = np.linalg.norm(gradients, axis=1)
    preliminary_count = min(len(branches), max(64, 4 * pool_size))
    preliminary = np.argsort(-correlation, kind="stable")[:preliminary_count]
    matrices = [x[:, preliminary] for x in matrices]
    coefficients = np.zeros((preliminary_count, len(responses)))
    maximum = float(np.max(correlation))
    magnitudes = np.zeros(preliminary_count)
    records = []
    step = 1 / max(1, preliminary_count)
    for fraction in np.geomspace(0.95, 0.02, paths):
        strength = maximum * fraction
        residual = 0.0
        for _iteration in range(iterations):
            coefficients, gradient, step = _proximal_step(
                matrices, responses, coefficients, strength, step
            )
            residual = _kkt_residual(coefficients, gradient, strength)
            if residual <= 1e-5 * max(1, maximum):
                break
            step *= 1.05
        magnitudes = np.maximum(magnitudes, np.linalg.norm(coefficients, axis=1))
        records.append(
            {
                "relative_strength": float(fraction),
                "strength": strength,
                "iterations": _iteration + 1,
                "kkt_residual": residual,
                "converged": residual <= 1e-5 * max(1, maximum),
            }
        )
    ranked = sorted(
        range(preliminary_count),
        key=lambda j: (
            -magnitudes[j],
            -correlation[preliminary[j]],
            branches[preliminary[j]],
        ),
    )
    pool = [branches[preliminary[j]] for j in ranked[:pool_size]]
    return pool, {
        "method": "null_whitened_group_lasso",
        "total_branches": len(branches),
        "preliminary_branches": preliminary_count,
        "pool_size": len(pool),
        "pool": pool,
        "path": records,
        "all_paths_converged": all(record["converged"] for record in records),
        "candidate_generation_only": True,
    }
