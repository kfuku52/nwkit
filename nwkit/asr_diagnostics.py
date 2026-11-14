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
from nwkit.optimization import _bounded_scalar_grid


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
    drop = 0.5 * float(chi2.ppf(level, 1))

    def difference(value):
        try:
            likelihood = float(log_likelihood(float(value)))
        except (ValueError, ArithmeticError, OverflowError):
            return -math.inf
        return (likelihood - maximum) + drop if math.isfinite(likelihood) else -math.inf

    left_grid = _bounded_scalar_grid(lower, mle, grid_size)
    right_grid = _bounded_scalar_grid(mle, upper, grid_size)

    lower_value, lower_limited = (
        (lower, True)
        if mle == lower
        else _profile_crossing(difference, left_grid[-2::-1], mle)
    )
    upper_value, upper_limited = (
        (upper, True)
        if mle == upper
        else _profile_crossing(difference, right_grid[1:], mle)
    )
    return ProfileLikelihoodInterval(
        lower_value,
        upper_value,
        level,
        lower_limited,
        upper_limited,
    )


def _profile_crossing(difference, points, mle):
    inside = mle
    for point in points:
        outside = float(point)
        value = difference(outside)
        if value > 0.0:
            inside = outside
            continue
        # Refine at a nonfinite profile value to distinguish a likelihood
        # crossing from exhaustion of the finite parameter domain.
        while not math.isfinite(value):
            middle = 0.5 * inside + 0.5 * outside
            if middle == inside or middle == outside:
                return inside, True
            midpoint_value = difference(middle)
            if midpoint_value > 0.0:
                inside = middle
            else:
                outside = middle
                value = midpoint_value
        if value == 0.0:
            return outside, False
        return float(
            brentq(
                difference,
                min(inside, outside),
                max(inside, outside),
                xtol=math.ulp(0.0),
                # A root close to zero may be hundreds of binary scales
                # closer than the nearest point of the initial profile grid.
                maxiter=2 * (np.finfo(float).maxexp - np.finfo(float).minexp),
            )
        ), False
    return inside, True


def _summary_statistics(values):
    """Summarize observation vectors along a contiguous final axis."""

    means = np.mean(values, axis=-1)
    lower = np.min(values, axis=-1)
    upper = np.max(values, axis=-1)
    return {
        "mean": means,
        "variance": np.var(values, axis=-1, ddof=1),
        "range": upper - lower,
        "max_abs_centered": np.maximum(upper - means, means - lower),
    }


def _predictive_work_units(observed, replicated):
    """Scale extreme trait values before moments and tail counts are computed."""

    size = max(
        float(np.max(observed)),
        -float(np.min(observed)),
        float(np.max(replicated)),
        -float(np.min(replicated)),
    )
    if not math.isfinite(size):
        raise ValueError(
            "Posterior predictive observations exceed floating-point range."
        )
    exponent = math.frexp(size)[1]
    if -250 <= exponent <= 250:
        return observed, replicated, 0
    return np.ldexp(observed, -exponent), np.ldexp(replicated, -exponent), exponent


def _restore_predictive_statistic(value, exponent):
    try:
        return math.ldexp(float(value), exponent)
    except OverflowError as exc:
        raise ValueError(
            "A posterior predictive statistic exceeds floating-point range."
        ) from exc


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
    replicated = simulate_gaussian_process(
        process,
        num_samples=num_simulations,
        seed=int(process_seed.generate_state(1)[0]),
        root_values=root_values,
        compiled_tree=compiled,
        _output_indices=tuple(record[0] for record in records),
    ).values
    errors = np.asarray([record[2] for record in records], dtype=float)
    if np.any(errors > 0.0):
        error_rng = np.random.default_rng(error_seed)
        replicated += error_rng.normal(0.0, errors, replicated.shape)
    observed = np.asarray([record[1] for record in records], dtype=float)
    observed, replicated, exponent = _predictive_work_units(observed, replicated)
    observed_statistics = _summary_statistics(observed)
    replicate_statistics = _summary_statistics(replicated)
    rows = []
    for statistic, observed_value in observed_statistics.items():
        values = replicate_statistics[statistic]
        upper = (1.0 + float(np.sum(values >= observed_value))) / (
            num_simulations + 1.0
        )
        lower = (1.0 + float(np.sum(values <= observed_value))) / (
            num_simulations + 1.0
        )
        statistic_exponent = 2 * exponent if statistic == "variance" else exponent
        moments = {
            "observed": observed_value,
            "replicate_mean": np.mean(values),
            "replicate_sd": np.std(values, ddof=1) if num_simulations > 1 else 0.0,
            "replicate_q025": np.quantile(values, 0.025),
            "replicate_q975": np.quantile(values, 0.975),
        }
        rows.append(
            {
                "statistic": statistic,
                **{
                    key: _restore_predictive_statistic(value, statistic_exponent)
                    for key, value in moments.items()
                },
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
