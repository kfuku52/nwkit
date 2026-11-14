"""Small-sample adjusted curvature intervals, not a posterior calculation.

The t and residual-variance corrections are exact in a local homoscedastic
Gaussian regression. Applying them to a nonlinear, mixed sequence/rate model
is an approximation; it does not establish nominal coverage for every family.
"""

from graphlib import TopologicalSorter

import numpy as np
from scipy.special import expit, logit
from scipy.stats import norm, t

from nwkit.radte_model import curvature_covariance


def rate_variance_information(fit, problem):
    """Count rate observations and locally fitted age/mean directions.

    The exact conditional sequence fit estimates SD from non-root branches.
    Marginal sequence inference instead observes one combined root length.
    Latent fitted rates are random effects, not residual fixed-effect degrees
    of freedom. Alignment columns are not independent draws of branch rates.
    """
    c = problem.chronology
    duration = c.durations(fit.ages)
    design = problem.free_design.toarray()
    if getattr(problem, "marginal", False):
        mapping = problem.likelihood.mapping
        count = len(problem.likelihood.center)
        jacobian = np.stack(
            [
                design[mapping == row].sum(axis=0) / duration[mapping == row].sum()
                for row in range(count)
            ]
        )
    else:
        indices = (
            problem.observation_indices
            if problem.likelihood is None
            else np.array(
                [i for i, node in enumerate(c.edges) if node.up is not c.gene],
                dtype=int,
            )
        )
        count = len(indices)
        jacobian = design[indices] / duration[indices, None]
    if count == 0:
        return 0, 0
    jacobian = np.column_stack([np.ones(count), jacobian])
    # Normalize nonzero columns so parameter units do not set the rank cutoff.
    norms = np.linalg.norm(jacobian, axis=0)
    jacobian[:, norms > 0] /= norms[norms > 0]
    rank = int(np.linalg.matrix_rank(jacobian))
    return count, count - rank


def chronology_domain(chronology):
    """Project all hard age/order constraints onto each coordinate.

    These are marginal feasible ranges, not a rectangular joint confidence set.
    Difference constraints on the acyclic shared-age graph permit propagation
    instead of a separate linear program for every free age.
    """
    lower, upper = chronology.lower.copy(), chronology.upper.copy()
    children: dict[int, set[int]] = {i: set() for i in range(len(chronology.groups))}
    for parent, child in zip(
        chronology.constraint_parent, chronology.constraint_child, strict=True
    ):
        children[int(parent)].add(int(child))
    order = list(TopologicalSorter(children).static_order())
    for parent in order:
        for child in children[parent]:
            lower[parent] = max(lower[parent], lower[child] + chronology.min_duration)
    for parent in reversed(order):
        for child in children[parent]:
            upper[child] = min(upper[child], upper[parent] - chronology.min_duration)
    return lower, upper


def studentized_intervals(fit, problem, level=0.95):
    """Use a residual-df t correction in bounded age coordinates.

    For an estimated rate variance, rescale the ML curvature covariance by
    n/(n-rank(X)) and replace the normal critical value by a t critical value.
    With a supplied SD, use the normal critical value without rescaling.
    The logit delta transformation respects the feasible age domain without
    clipping endpoints or changing the fitted ages, rates, or constraints.
    """
    covariance = curvature_covariance(fit, problem, level)
    if covariance is None:
        return
    factor = 1.0
    critical = norm.ppf((1 + level) / 2)
    if problem.rate_variance_estimated:
        count, degrees = rate_variance_information(fit, problem)
        fit.diagnostics.extend(
            [
                f"studentized_rate_observations={count}",
                f"studentized_residual_df={degrees}",
            ]
        )
        if degrees <= 0:
            fit.interval_status = "unavailable-no-residual-rate-degrees-of-freedom"
            return
        factor = count / degrees
        critical = t.ppf((1 + level) / 2, degrees)
    c = problem.chronology
    lower, upper = chronology_domain(c)
    lo, hi = lower[problem.free], upper[problem.free]
    age = fit.ages[problem.free]
    if np.any(age <= lo) or np.any(age >= hi):
        fit.interval_status = "unavailable-active-bound-use-profile-or-bootstrap"
        return
    center = logit((age - lo) / (hi - lo))
    derivative = 1 / (age - lo) + 1 / (hi - age)
    width = critical * np.sqrt(factor * np.diag(covariance)) * derivative
    if not np.isfinite(width).all():
        fit.interval_status = "unavailable-nonfinite-studentized-curvature"
        return
    fit.interval_lower, fit.interval_upper = fit.ages.copy(), fit.ages.copy()
    fit.interval_lower[problem.free] = lo + (hi - lo) * expit(center - width)
    fit.interval_upper[problem.free] = lo + (hi - lo) * expit(center + width)
    fit.interval_status = (
        "conditional-studentized-curvature"
        if problem.rate_variance_estimated
        else "conditional-bounded-normal-curvature"
    )
    fit.diagnostics.extend(
        [
            f"studentized_variance_factor={factor:.12g}",
            "curvature_age_transform=logit-feasible-domain",
            "studentized_curvature_is_a_local_approximation_not_a_posterior",
        ]
    )
