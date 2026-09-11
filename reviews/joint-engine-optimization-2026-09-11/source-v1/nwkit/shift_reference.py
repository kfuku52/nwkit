"""Fixed-parameter OU shift evaluation using NWKIT's Gaussian tree engine.

This is a validation evaluator, not a shift search or parameter optimizer.
The baseline uses the explicit convention root mean = baseline optimum.
At alpha=0, finite mean effects are represented as BM mean offsets, not optima.
"""

import math

import numpy as np
from scipy.linalg import cho_factor, cho_solve

from nwkit.gaussian_tree import (
    GaussianRootPrior,
    GaussianTransition,
    GaussianTreeProcess,
    ou_transition,
)
from nwkit.shift_math import observation_variances, remaining_heights


def build_shift_process(tree, *, alpha, sigma2, intercept, mean_effects, root_model):
    if root_model not in {"OUfixedRoot", "OUrandomRoot"}:
        raise ValueError("Unsupported kfl1ou root model.")
    if (
        not all(math.isfinite(x) for x in (alpha, sigma2, intercept))
        or alpha < 0
        or sigma2 <= 0
    ):
        raise ValueError("Shift parameters must be finite, alpha >= 0 and sigma2 > 0.")
    if root_model == "OUrandomRoot" and alpha == 0:
        raise ValueError(
            "Stationary/random-root OU is undefined at alpha=0; use a fixed-root BM model explicitly."
        )
    nodes = set(tree.traverse())
    if (
        tree in mean_effects
        or not set(mean_effects) <= nodes
        or not all(math.isfinite(x) for x in mean_effects.values())
    ):
        raise ValueError("Shift effects must be finite and identify non-root nodes.")
    depths = {tree: 0.0}
    for node in tree.traverse("preorder"):
        if not node.is_root:
            if not math.isfinite(node.dist) or node.dist <= 0:
                raise ValueError("Shift reference requires positive branch lengths.")
            depths[node] = depths[node.up] + node.dist
    heights = [depths[n] for n in tree.leaves()]
    height = max(heights)
    if max(heights) - min(heights) > height * 1e-8:
        raise ValueError("Shift reference requires an ultrametric tree.")
    remaining = remaining_heights(tree)
    root_variance = (
        sigma2 / (2 * alpha) if root_model == "OUrandomRoot" and alpha > 0 else 0.0
    )
    root = GaussianRootPrior(
        "stationary" if root_variance > 0 else "fixed", intercept, root_variance
    )
    optima = {tree: intercept}
    transitions = {}
    for node in tree.traverse("preorder"):
        if node.is_root:
            continue
        if alpha == 0:
            transitions[node] = GaussianTransition(
                1.0, mean_effects.get(node, 0.0), sigma2 * node.dist
            )
        else:
            denominator = -math.expm1(-alpha * remaining[node.up])
            optima[node] = optima[node.up] + mean_effects.get(node, 0.0) / denominator
            transitions[node], _ = ou_transition(node.dist, alpha, sigma2, optima[node])
    return GaussianTreeProcess(tree, transitions, root, "OU-SHIFT", alpha)


def evaluate_shift_model(tree, *, observations, standard_errors=None, **parameters):
    """Return tip mean, covariance and ordinary ML density without refitting.

    Dense covariance is intentional for small validation fixtures. Production
    shift search does not invoke this O(tips**2) reference evaluator.
    """
    process = build_shift_process(tree, **parameters)
    leaves = list(tree.leaves())
    means, _ = process.marginal_moments()
    mean = np.asarray([means[n] for n in leaves])
    covariance = process.covariance(leaves)
    observed = np.asarray(observations, dtype=float)
    if observed.shape != mean.shape or not np.isfinite(observed).all():
        raise ValueError("One finite observation per tip is required.")
    if standard_errors is not None:
        variance = np.asarray(observation_variances(standard_errors))
        if variance.shape != mean.shape:
            raise ValueError("One standard error per tip is required.")
        covariance = covariance + np.diag(variance)
    factor = cho_factor(covariance, lower=True)
    residual = observed - mean
    log_likelihood = -0.5 * (
        len(mean) * math.log(2 * math.pi)
        + 2 * np.log(np.diag(factor[0])).sum()
        + residual @ cho_solve(factor, residual)
    )
    return mean, covariance, float(log_likelihood)
