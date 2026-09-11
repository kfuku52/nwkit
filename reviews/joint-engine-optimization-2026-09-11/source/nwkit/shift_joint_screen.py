"""Covariance-aware candidate proposals and cached joint layout profiles.

The full joint covariance whitens both the response and every candidate effect.
No diagonal-only prefilter is used. Penalized paths only propose layouts; final
likelihoods always refit covariance and means without a screening penalty.
"""

import math

import numpy as np

from nwkit.shift_joint_model import joint_factor
from nwkit.shift_native_screen import _effect_design, descendant_design


def joint_screen_matrix(
    data, fitted, branches=None, *, normalize=True, memory_limit=512 * 1024**2
):
    joint = fitted["joint_fit"]
    n, p = data.values.shape
    branches = tuple(
        sorted(set(data.tree.branch_ids) - {0}) if branches is None else branches
    )
    mask = np.isfinite(data.values)
    observed = int(mask.sum())
    estimated = 8 * (
        3 * observed * len(branches) * p
        + 4 * observed * (p + 1)
        + len(data.tree.branch_ids) * 8 * p * p
    )
    if estimated > memory_limit:
        raise ValueError(
            f"Joint screening requires approximately {estimated} bytes, beyond the declared memory budget {memory_limit}."
        )
    factor, _, _ = joint_factor(
        data,
        joint.alpha_height,
        joint.covariance_coordinate,
        joint.measurement_variance,
        root_model=joint.root_model,
    )
    intercept = np.broadcast_to(np.eye(p), (n, p, p))[mask]
    white = factor.apply(np.column_stack((data.values[mask], intercept)))
    baseline, _ = np.linalg.qr(white[:, 1:], mode="reduced")
    response = white[:, 0] - baseline @ (baseline.T @ white[:, 0])
    matrix = np.empty((observed, len(branches) * p))
    for first in range(0, len(branches), 16):
        selected = branches[first : first + 16]
        descendant, _ = descendant_design(data.tree, selected)
        raw = np.zeros((n, p, len(selected) * p))
        for trait in range(p):
            values = _effect_design(
                data.tree, descendant, float(joint.alpha_height[trait]), selected
            )
            raw[:, trait, trait::p] = values
        block = factor.apply(raw[mask])
        block -= baseline @ (baseline.T @ block)
        matrix[:, first * p : (first + len(selected)) * p] = block
    if normalize:
        # Orthonormalize within each branch group. Group selection then depends
        # on its joint mean-effect subspace, not the choice of trait coordinates.
        for j in range(len(branches)):
            block = matrix[:, j * p : (j + 1) * p]
            u, s, _ = np.linalg.svd(block, full_matrices=False)
            keep = s > np.finfo(float).eps * max(block.shape) * max(1.0, s[0])
            block[:] = u * keep[None]
    return matrix, response, branches, factor.log_determinant


def joint_path(matrix, response, p, *, iterations, paths):
    groups = matrix.shape[1] // p
    initial_gradient = (matrix.T @ response).reshape(groups, p)
    maximum = float(np.max(np.linalg.norm(initial_gradient, axis=1)))
    coefficients = np.zeros((groups, p))
    step = 1 / max(1.0, float(np.sum(matrix * matrix)))
    for fraction in paths:
        strength = maximum * fraction
        for _iteration in range(iterations):
            residual = matrix @ coefficients.ravel() - response
            gradient = (matrix.T @ residual).reshape(groups, p)
            loss = 0.5 * float(residual @ residual)
            for _ in range(40):
                trial = coefficients - step * gradient
                norms = np.linalg.norm(trial, axis=1)
                trial *= np.maximum(
                    0.0, 1 - step * strength / np.maximum(norms, np.finfo(float).tiny)
                )[:, None]
                change = trial - coefficients
                r = matrix @ trial.ravel() - response
                if 0.5 * float(r @ r) <= loss + float(
                    np.sum(gradient * change)
                ) + float(np.sum(change * change)) / (2 * step) + 1e-12 * max(
                    1.0, loss
                ):
                    break
                step *= 0.5
            else:
                raise ValueError(
                    "Joint group-lasso backtracking failed its quadratic bound."
                )
            coefficients = trial
            gradient = (matrix.T @ r).reshape(groups, p)
            norms = np.linalg.norm(coefficients, axis=1)
            active = norms > 0
            kkt = np.maximum(0.0, np.linalg.norm(gradient, axis=1) - strength)
            kkt[active] = np.linalg.norm(
                gradient[active]
                + strength * coefficients[active] / norms[active, None],
                axis=1,
            )
            error = float(np.max(kkt))
            if error <= 1e-5 * max(1.0, maximum):
                break
            step *= 1.05
        yield (
            coefficients.copy(),
            {
                "relative_strength": float(fraction),
                "strength": strength,
                "iterations": _iteration + 1,
                "kkt_residual": error,
                "converged": error <= 1e-5 * max(1.0, maximum),
            },
        )


def joint_group_lasso_screen(
    data, null_fit, *, pool_size=24, iterations=150, paths=6, memory_limit=512 * 1024**2
):
    if min(pool_size, iterations, paths, memory_limit) < 1:
        raise ValueError("Joint screening budgets must be positive.")
    matrix, response, branches, _ = joint_screen_matrix(
        data, null_fit, memory_limit=memory_limit
    )
    p = len(data.trait_names)
    correlation = np.linalg.norm(
        (matrix.T @ response).reshape(len(branches), p), axis=1
    )
    count = min(len(branches), max(64, 4 * pool_size))
    preliminary = np.argsort(-correlation, kind="stable")[:count]
    columns = (preliminary[:, None] * p + np.arange(p)[None]).ravel()
    magnitudes = np.zeros(count)
    records = []
    for coefficients, record in joint_path(
        matrix[:, columns],
        response,
        p,
        iterations=iterations,
        paths=np.geomspace(0.95, 0.02, paths),
    ):
        magnitudes = np.maximum(magnitudes, np.linalg.norm(coefficients, axis=1))
        records.append(record)
    ranked = sorted(
        range(count),
        key=lambda j: (
            -magnitudes[j],
            -correlation[preliminary[j]],
            branches[preliminary[j]],
        ),
    )
    pool = [branches[preliminary[j]] for j in ranked[:pool_size]]
    return pool, {
        "method": "joint_covariance_whitened_group_lasso",
        "total_branches": len(branches),
        "preliminary_branches": count,
        "pool_size": len(pool),
        "pool": pool,
        "path": records,
        "all_paths_converged": all(r["converged"] for r in records),
        "candidate_generation_only": True,
    }


class JointQuickProfile:
    def __init__(
        self, data, null_fit, branches, fixed_alpha=None, *, memory_limit=512 * 1024**2
    ):
        self.data = data
        self.null_fit = null_fit
        self.matrix, self.response, branches, logdet = joint_screen_matrix(
            data, null_fit, branches, normalize=False, memory_limit=memory_limit
        )
        self.columns = {b: i for i, b in enumerate(branches)}
        self.evaluations = 0
        observed = np.isfinite(data.values)
        self.constant = -0.5 * (
            int(observed.sum()) * math.log(2 * math.pi) + logdet
        ) - sum(
            int(observed[:, j].sum()) * math.log(data.scales[j])
            for j in range(len(data.trait_names))
        )

    def score(self, layout):
        self.evaluations += 1
        p, q = len(self.data.trait_names), len(layout.groups) - 1
        if not q:
            return self.constant - 0.5 * float(self.response @ self.response)
        x = np.zeros((len(self.response), p * q))
        labels = layout.node_groups(self.data.tree)
        for branch in layout.shifts:
            i = self.data.tree.indices_by_branch[branch]
            column = self.columns[branch]
            block = self.matrix[:, column * p : (column + 1) * p]
            for sign, group in (
                (1, labels[i]),
                (-1, labels[self.data.tree.compiled.parents[i]]),
            ):
                if group:
                    x[:, (group - 1) :: q] += sign * block
        norms = np.linalg.norm(x, axis=0)
        if np.any(norms == 0) or np.linalg.matrix_rank(x / norms) < x.shape[1]:
            return -math.inf
        basis, _ = np.linalg.qr(x / norms, mode="reduced")
        residual = self.response - basis @ (basis.T @ self.response)
        return self.constant - 0.5 * float(residual @ residual)
