"""Exact Gaussian contrast intervals when log duration is a linear parameter.

The identifiable special case has one free age, and every affected duration
equals a positive constant times (age - offset). K and the chronology are
fixed. General internal duplications or sequence data do not satisfy this.
"""

import numpy as np
from scipy.stats import norm, t

from nwkit.radte_studentized import chronology_domain


def exact_log_duration_intervals(fit, problem, level=0.95):
    """Return False outside the exact model; otherwise set fit's interval."""
    if problem.likelihood is not None or len(problem.free) != 1:
        return False
    indices = problem.observation_indices
    slope = problem.free_design.toarray()[indices, 0]
    affected = slope != 0
    if not affected.any() or affected.all() or np.any(slope[affected] <= 0):
        return False
    group = problem.free[0]
    duration = problem.chronology.durations(fit.ages)[indices]
    offsets = fit.ages[group] - duration[affected] / slope[affected]
    if not np.allclose(offsets, offsets[0], rtol=0, atol=1e-12):
        return False
    if not 0 < level < 1:
        raise ValueError("Interval level must be between zero and one.")
    offset = float(offsets[0])
    # y = mu + indicator*log(age-offset) + correlated Gaussian error.
    normalizer = duration.copy()
    normalizer[affected] = slope[affected]
    y = problem.log_observed[indices] - np.log(normalizer)
    design = np.column_stack([np.ones(len(y)), affected])
    q = problem.precision.toarray()
    information = design.T @ q @ design
    beta = np.linalg.solve(information, design.T @ q @ y)
    residual = y - design @ beta
    df = len(y) - 2
    fit.interval_lower = fit.interval_upper = None
    if problem.rate_variance_estimated and fit.log_rate_sd == 0:
        fit.interval_status = "unavailable-strict-clock-limit"
        return True
    if problem.rate_variance_estimated and df <= 0:
        fit.interval_status = "unavailable-no-residual-rate-degrees-of-freedom"
        return True
    variance = (
        float(residual @ q @ residual / df)
        if problem.rate_variance_estimated
        else problem.rate_sd**2
    )
    if variance <= 0:
        fit.interval_status = "unavailable-strict-clock-limit"
        return True
    critical = (
        t.ppf((1 + level) / 2, df)
        if problem.rate_variance_estimated
        else norm.ppf((1 + level) / 2)
    )
    width = critical * np.sqrt(variance * np.linalg.inv(information)[1, 1])
    lo, hi = chronology_domain(problem.chronology)
    # Intersect the exact contrast confidence set with the known age domain.
    # This is set intersection, not clipping a delta-method Gaussian interval.
    log_lower = max(float(beta[1] - width), float(np.log(lo[group] - offset)))
    log_upper = min(float(beta[1] + width), float(np.log(hi[group] - offset)))
    if log_lower > log_upper:
        fit.interval_status = "unavailable-empty-exact-log-duration-set"
        return True
    fit.interval_lower, fit.interval_upper = fit.ages.copy(), fit.ages.copy()
    fit.interval_lower[group] = max(lo[group], offset + np.exp(log_lower))
    fit.interval_upper[group] = min(hi[group], offset + np.exp(log_upper))
    fit.interval_status = "conditional-exact-log-duration-" + (
        "t" if problem.rate_variance_estimated else "normal"
    )
    fit.diagnostics.extend(
        [
            f"exact_log_duration_residual_df={df}",
            "exact_log_duration_conditions=fixed_K_and_branch_lengths_and_chronology",
        ]
    )
    return True
