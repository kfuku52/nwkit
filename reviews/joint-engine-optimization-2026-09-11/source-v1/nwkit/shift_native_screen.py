"""Deterministic multivariate group-lasso candidate generation.

Whitening uses the fitted null covariance. Screening and penalized coefficients
only propose locations; all retained models are fitted without this penalty.
Every bootstrap search must regenerate this data-dependent candidate set.
"""

import numpy as np

from nwkit.gaussian_whitening import TreeWhitening
from nwkit.shift_native_model import covariance_geometry

# Bound temporary whitening storage independently of the number of branches.
_SCREEN_BLOCK_COLUMNS = 512


def _screen_block_width(factor):
    # Deep trees have many small rotation batches: amortize their traversal over
    # more columns. Round up to a power of two, with bounded temporary storage.
    traversal_width = 1 << (2 * len(factor.batches) - 1).bit_length()
    return min(4096, max(_SCREEN_BLOCK_COLUMNS, traversal_width))


def _descendant_intervals(tree):
    ranges = {node: (i, i + 1) for i, node in enumerate(tree.compiled.leaf_indices)}
    for index in tree.compiled.postorder:
        children = tree.compiled.children[index]
        if children:
            ranges[index] = (ranges[children[0]][0], ranges[children[-1]][1])
    return {tree.branch_ids[index]: interval for index, interval in ranges.items()}


def descendant_design(tree, branches=None):
    """Construct only requested columns using contiguous descendant intervals."""
    branches = tuple(
        sorted(set(tree.branch_ids) - {0}) if branches is None else branches
    )
    return _interval_design(
        len(tree.leaf_names), _descendant_intervals(tree), branches
    ), branches


def _interval_design(tips, intervals, branches):
    matrix = np.zeros((tips, len(branches)))
    for column, branch in enumerate(branches):
        first, last = intervals[branch]
        matrix[first:last, column] = 1
    return matrix


def _effect_design(tree, design, alpha, branches=None):
    """Unstandardized optimum-increment columns, including near-ultrametric tips."""
    depths = tree.times.copy()
    for i in range(1, len(depths)):
        depths[i] += depths[tree.compiled.parents[i]]
    branches = sorted(set(tree.branch_ids) - {0}) if branches is None else branches
    indices = {branch: i for i, branch in enumerate(tree.branch_ids)}
    parents = [tree.compiled.parents[indices[b]] for b in branches]
    ages = np.maximum(
        0, depths[list(tree.compiled.leaf_indices), None] - depths[parents]
    )
    weights = (
        ages
        if alpha == 0
        else np.ones_like(ages)
        if np.isinf(alpha)
        else -np.expm1(-alpha * ages) / -np.expm1(-alpha)
    )
    return design * weights


def _whitened_matrices(data, null_fit, memory_limit, *, optimum_increments=False):
    n, p = data.values.shape
    # Keep the existing conservative admission budget. Blocking reduces actual
    # temporary storage, but this guard is not a whole-process peak RAM bound.
    estimated = (
        8 * n * (len(data.tree.branch_ids) - 1) * (p + (8 if optimum_increments else 4))
    )
    if estimated > memory_limit:
        raise ValueError(
            f"Group-lasso screening needs approximately {estimated} bytes; increase --search-memory-mb or reduce the input."
        )
    branches = tuple(sorted(set(data.tree.branch_ids) - {0}))
    intervals = _descendant_intervals(data.tree)
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
        x = np.empty((sum(mask), len(branches)))
        white_squared_norm = 0.0
        block_width = _screen_block_width(factor)
        for first in range(0, len(branches), block_width):
            selected = branches[first : first + block_width]
            design = _interval_design(n, intervals, selected)
            if optimum_increments:
                design = _effect_design(data.tree, design, fit.alpha_height, selected)
            if first == 0:
                # Share the first traversal with the response and intercept,
                # especially when a deep tree fits into one column block.
                white = factor.apply(
                    np.column_stack(
                        (np.ones(sum(mask)), data.values[mask, trait], design[mask])
                    )
                )
                intercept = white[:, 0] / np.linalg.norm(white[:, 0])
                y = white[:, 1] - intercept * (intercept @ white[:, 1])
                white = white[:, 2:]
            else:
                white = factor.apply(design[mask])
            white_squared_norm += float(np.sum(white * white))
            x[:, first : first + len(selected)] = (
                white - intercept[:, None] * (intercept @ white)[None, :]
            )
        norms = np.linalg.norm(x, axis=0)
        # Preserve the full-design threshold, not a block-dependent threshold.
        informative = norms > (
            np.finfo(float).eps * max(1.0, np.sqrt(white_squared_norm)) * 100
        )
        for first in range(0, len(branches), block_width):
            block = x[:, first : first + block_width]
            keep = informative[first : first + block_width]
            if not optimum_increments:
                block[:, keep] /= norms[first : first + block_width][keep]
            block[:, ~keep] = 0
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
