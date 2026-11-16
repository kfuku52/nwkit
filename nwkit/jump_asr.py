"""Brownian motion plus a Gaussian compound-Poisson evolutionary jump process."""

import math
from functools import lru_cache

from nwkit.gaussian_tree import (
    GaussianRootPrior,
    GaussianTransition,
    GaussianTreeProcess,
)
from nwkit.latent_gaussian import integrate_histories


def finite_parameter(value, name, *, positive=False):
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be explicitly specified and finite.") from exc
    if not math.isfinite(number) or number < 0 or (positive and number == 0):
        qualifier = "positive" if positive else "nonnegative"
        raise ValueError(f"{name} must be finite and {qualifier}.")
    return number


def fit_jump_bm(
    tree,
    observed,
    *,
    sigma2,
    jump_rate,
    jump_sd,
    standard_errors=None,
    history_samples=1000,
    seed=None,
):
    """Integrate Poisson edge counts with fixed BM/jump parameters.

    Conditional edge variance is sigma2*t + count*jump_sd**2. Jump sizes
    themselves are integrated analytically. Positive Brownian diffusion avoids
    singular atom/density mixtures when observations are exact.
    """
    sigma2 = finite_parameter(sigma2, "sigma2", positive=True)
    jump_rate = finite_parameter(jump_rate, "jump_rate")
    jump_sd = finite_parameter(jump_sd, "jump_sd")
    branches = tuple(node for node in tree.traverse() if not node.is_root)
    lengths = {node: finite_parameter(node.dist, "branch length") for node in branches}
    expected_counts = {node: jump_rate * length for node, length in lengths.items()}
    if any(count > 1e8 for count in expected_counts.values()):
        raise ValueError("Expected branch jump counts exceed the supported range.")

    @lru_cache(maxsize=1024)
    def transition(length, count):
        return GaussianTransition(1.0, 0.0, sigma2 * length + count * jump_sd**2)

    root = GaussianRootPrior("flat", variance=None)

    def propose(rng):
        counts = {
            node: int(rng.poisson(expected_count))
            for node, expected_count in expected_counts.items()
        }
        process = GaussianTreeProcess(
            tree,
            {node: transition(lengths[node], count) for node, count in counts.items()},
            root,
            "JUMP-BM",
        )
        return process, counts

    return integrate_histories(
        tree,
        observed,
        standard_errors,
        propose,
        model="JUMP-BM",
        samples=history_samples,
        seed=seed,
        parameters={"sigma2": sigma2, "jump_rate": jump_rate, "jump_sd": jump_sd},
    )
