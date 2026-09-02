"""Simulation, posterior-predictive, bootstrap, and profile diagnostics."""

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import chi2

from nwkit.gaussian_inference import (
    condition_gaussian_tree,
    simulate_gaussian_process,
)


@dataclass(frozen=True, slots=True)
class ProfileLikelihoodInterval:
    lower: float
    upper: float
    level: float
    lower_boundary_limited: bool
    upper_boundary_limited: bool


def profile_likelihood_interval(
    log_likelihood,
    mle,
    bounds,
    *,
    level=0.95,
    grid_size=65,
):
    """Invert a one-parameter profile-likelihood ratio on finite bounds."""

    mle = float(mle)
    lower, upper = float(bounds[0]), float(bounds[1])
    level = float(level)
    if not (
        math.isfinite(lower)
        and math.isfinite(mle)
        and math.isfinite(upper)
        and lower < upper
        and lower <= mle <= upper
    ):
        raise ValueError(
            "Profile likelihood requires finite bounds containing the MLE."
        )
    if isinstance(grid_size, bool) or int(grid_size) != grid_size or grid_size < 3:
        raise ValueError("Profile likelihood requires at least three grid points.")
    grid_size = int(grid_size)
    if not 0.0 < level < 1.0:
        raise ValueError("Profile-likelihood level must be between zero and one.")
    maximum = float(log_likelihood(mle))
    if not math.isfinite(maximum):
        raise ValueError("The profile likelihood is not finite at the MLE.")
    cutoff = maximum - 0.5 * float(chi2.ppf(level, 1))

    def difference(value):
        try:
            likelihood = float(log_likelihood(float(value)))
        except (ValueError, ArithmeticError, OverflowError):
            return -math.inf
        return likelihood - cutoff if math.isfinite(likelihood) else -math.inf

    left_grid = np.linspace(lower, mle, grid_size)
    right_grid = np.linspace(mle, upper, grid_size)

    def crossing(grid, reverse):
        values = [difference(value) for value in grid]
        pairs = tuple(range(len(grid) - 1))
        if reverse:
            pairs = pairs[::-1]
        for index in pairs:
            first, second = values[index], values[index + 1]
            if not (math.isfinite(first) and math.isfinite(second)):
                finite_index = None
                if math.isfinite(first) and first > 0.0:
                    finite_index = index
                elif math.isfinite(second) and second > 0.0:
                    finite_index = index + 1
                if finite_index is not None:
                    return float(grid[finite_index]), True
                continue
            if first == 0.0:
                return float(grid[index]), False
            if second == 0.0:
                return float(grid[index + 1]), False
            if (first < 0.0 < second) or (second < 0.0 < first):
                return float(
                    brentq(
                        difference,
                        float(grid[index]),
                        float(grid[index + 1]),
                        xtol=1e-10,
                    )
                ), False
        finite_indices = [
            index for index, value in enumerate(values) if math.isfinite(value)
        ]
        if not finite_indices:
            raise ValueError("The profile likelihood is finite only at the MLE.")
        boundary_index = finite_indices[0] if reverse else finite_indices[-1]
        boundary = float(grid[boundary_index])
        return boundary, True

    lower_value, lower_limited = (
        (lower, True) if mle == lower else crossing(left_grid, True)
    )
    upper_value, upper_limited = (
        (upper, True) if mle == upper else crossing(right_grid, False)
    )
    return ProfileLikelihoodInterval(
        lower_value,
        upper_value,
        level,
        lower_limited,
        upper_limited,
    )


def _summary_statistics(values):
    values = np.asarray(values, dtype=float)
    centered = values - float(np.mean(values))
    return {
        "mean": float(np.mean(values)),
        "variance": float(np.var(values, ddof=1)) if len(values) > 1 else 0.0,
        "range": float(np.ptp(values)),
        "max_abs_centered": float(np.max(np.abs(centered))),
    }


def gaussian_posterior_predictive(
    process,
    values_by_leaf,
    *,
    standard_errors=None,
    num_simulations=1000,
    seed=None,
):
    """Simulate replicated tip datasets from fitted parameters and root posterior."""

    if isinstance(num_simulations, bool) or int(num_simulations) != num_simulations:
        raise ValueError("num_simulations must be a positive integer.")
    num_simulations = int(num_simulations)
    if num_simulations <= 0:
        raise ValueError("num_simulations must be a positive integer.")
    conditioned = condition_gaussian_tree(
        process, values_by_leaf, standard_errors=standard_errors
    )
    compiled = conditioned.compiled_tree
    records = []
    for name, node_index in compiled.leaf_index_by_name.items():
        value = values_by_leaf.get(name)
        if value is None:
            continue
        error = 0.0 if standard_errors is None else float(standard_errors[name])
        records.append((node_index, float(value), error))
    if len(records) < 2:
        raise ValueError("Posterior predictive checks require two observed tips.")
    seed_sequence = np.random.SeedSequence(seed)
    root_seed, process_seed, error_seed = seed_sequence.spawn(3)
    root_rng = np.random.default_rng(root_seed)
    root = conditioned.marginals[process.tree]
    root_values = (
        np.full(num_simulations, root.mean)
        if root.variance == 0.0
        else root_rng.normal(root.mean, math.sqrt(root.variance), num_simulations)
    )
    simulations = simulate_gaussian_process(
        process,
        num_samples=num_simulations,
        seed=int(process_seed.generate_state(1)[0]),
        root_values=root_values,
        compiled_tree=compiled,
    ).values
    indices = np.asarray([record[0] for record in records], dtype=int)
    replicated = simulations[:, indices]
    errors = np.asarray([record[2] for record in records], dtype=float)
    if np.any(errors > 0.0):
        error_rng = np.random.default_rng(error_seed)
        replicated += error_rng.normal(0.0, errors, replicated.shape)
    observed = np.asarray([record[1] for record in records], dtype=float)
    observed_statistics = _summary_statistics(observed)
    replicate_statistics = {
        key: np.asarray(
            [_summary_statistics(row)[key] for row in replicated], dtype=float
        )
        for key in observed_statistics
    }
    rows = []
    for statistic, observed_value in observed_statistics.items():
        values = replicate_statistics[statistic]
        upper = (1.0 + float(np.sum(values >= observed_value))) / (
            num_simulations + 1.0
        )
        lower = (1.0 + float(np.sum(values <= observed_value))) / (
            num_simulations + 1.0
        )
        rows.append(
            {
                "statistic": statistic,
                "observed": observed_value,
                "replicate_mean": float(np.mean(values)),
                "replicate_sd": (
                    float(np.std(values, ddof=1)) if num_simulations > 1 else 0.0
                ),
                "replicate_q025": float(np.quantile(values, 0.025)),
                "replicate_q975": float(np.quantile(values, 0.975)),
                "p_lower": lower,
                "p_upper": upper,
                "p_two_sided": min(1.0, 2.0 * min(lower, upper)),
                "num_simulations": num_simulations,
            }
        )
    return pd.DataFrame(rows)


def parametric_bootstrap(
    simulator,
    fitter,
    extractor,
    *,
    num_simulations,
    seed=None,
):
    """Run a reproducible callback-based parametric bootstrap."""

    if isinstance(num_simulations, bool) or int(num_simulations) != num_simulations:
        raise ValueError("num_simulations must be a positive integer.")
    num_simulations = int(num_simulations)
    if num_simulations <= 0:
        raise ValueError("num_simulations must be a positive integer.")
    seeds = np.random.SeedSequence(seed).spawn(num_simulations)
    rows = []
    failures = []
    for replicate, child_seed in enumerate(seeds):
        replicate_seed = int(child_seed.generate_state(1)[0])
        try:
            simulated = simulator(replicate_seed)
            fitted = fitter(simulated)
            row = dict(extractor(fitted))
            row.update(replicate=replicate, seed=replicate_seed, fit_status="ok")
            rows.append(row)
        except (ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
            failures.append(
                {
                    "replicate": replicate,
                    "seed": replicate_seed,
                    "fit_status": "failed",
                    "error": str(exc),
                }
            )
    if not rows:
        raise ValueError("Every parametric-bootstrap refit failed.")
    return pd.DataFrame(rows), pd.DataFrame(failures)
