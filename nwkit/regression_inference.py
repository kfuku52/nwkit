"""Null-generated calibration and numerical inversion of coefficient tests."""

from dataclasses import dataclass
from typing import Callable

import numpy as np


@dataclass(frozen=True)
class BootstrapTest:
    statistic: float
    p_value: float
    monte_carlo_se: float
    replicates: int


def bootstrap_test(observed, simulate_statistic: Callable, *, replicates, seed):
    """Calibrate a specified statistic under its null, without failure filtering.

    A failure aborts inference: replacing a failed dataset with another one
    conditions the reference distribution on successful optimization.
    """
    if replicates < 2 or not np.isfinite(observed):
        raise ValueError("Bootstrap needs a finite statistic and at least two draws.")
    rng = np.random.default_rng(seed)
    exceedances = 0
    for index in range(replicates):
        try:
            statistic = float(simulate_statistic(rng))
            if not np.isfinite(statistic):
                raise ValueError("Non-finite bootstrap statistic.")
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
            raise RuntimeError(
                f"Null bootstrap failed at dataset {index + 1}/{replicates}; "
                "no p-value was computed and no failed dataset was replaced."
            ) from exc
        exceedances += statistic >= observed
    p_value = (1 + exceedances) / (replicates + 1)
    return BootstrapTest(
        observed,
        p_value,
        float(np.sqrt(p_value * (1.0 - p_value) / (replicates + 1))),
        replicates,
    )


def objective_difference(full_objective, null_objective):
    """A materially worse unconstrained fit is an optimization failure."""
    difference = float(null_objective - full_objective)
    tolerance = 1e-7 * max(1.0, abs(full_objective), abs(null_objective))
    if not np.isfinite(difference) or difference < -tolerance:
        raise RuntimeError("Unconstrained fit is worse than its constrained null.")
    return max(0.0, 2.0 * difference)


def invert_bootstrap_grid(test, grid, confidence_level):
    """Evaluate the *pointwise* inverted test on an explicit candidate grid.

    Returns accepted candidates, not interpolated confidence intervals. A
    finite grid cannot establish exclusion between or outside its points, and
    non-connected sets must not be collapsed into a misleading finite interval.
    Each grid candidate is calibrated under that candidate's own null model.
    """
    grid = np.asarray(grid, dtype=float)
    if (
        grid.ndim != 1
        or not len(grid)
        or not np.isfinite(grid).all()
        or np.any(np.diff(grid) <= 0)
    ):
        raise ValueError(
            "Bootstrap profile grid must be finite and strictly increasing."
        )
    if not 0.0 < confidence_level < 1.0:
        raise ValueError("Confidence level must be between zero and one.")
    evaluations = []
    for candidate in grid:
        result = test(float(candidate))
        evaluations.append(
            {
                "null_value": float(candidate),
                "p_value": result.p_value,
                "monte_carlo_se": result.monte_carlo_se,
                "accepted": bool(result.p_value > 1.0 - confidence_level),
            }
        )
    return evaluations
