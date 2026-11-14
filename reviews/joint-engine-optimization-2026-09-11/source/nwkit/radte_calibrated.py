"""Experimental constrained parametric-bootstrap profile inversion.

The nuisance plug-in and finite age grid are approximations, not an exact
confidence procedure. Fits use a marginal likelihood, never a joint MAP
penalty as a likelihood ratio. Sequence draws regenerate both rates and sites.
"""

import copy
from dataclasses import replace

import numpy as np
from scipy.optimize import linprog
from scipy.stats import binomtest

from nwkit.radte_model import solve_problem
from nwkit.radte_simulation import simulate_rate_lengths, simulate_sequence_likelihood
from nwkit.radte_studentized import chronology_domain
from nwkit.radte_uncertainty import _profile_delta, validate_profile_approximation


def _problem_like(problem, chronology, likelihood):
    settings = dict(rho=problem.rho, likelihood=likelihood, rate_sd=problem.rate_sd)
    if getattr(problem, "marginal", False):
        settings.update(
            rate_sd=problem.fixed_sd, quadrature_points=problem.quadrature_points
        )
    return type(problem)(chronology, **settings)


def _solve(problem, *, starts, maxiter, seed):
    parameters, attempts = solve_problem(
        problem, starts=starts, maxiter=maxiter, seed=seed
    )
    if getattr(problem, "marginal", False):
        from nwkit.radte_marginal import refine_quadrature

        problem, parameters, extra = refine_quadrature(
            problem, parameters, starts=starts, maxiter=maxiter, seed=seed
        )
        attempts.extend(extra)
        validate_profile_approximation(problem, parameters)
    if any(a.get("nonunique_age_solution") for a in attempts):
        raise ValueError("Calibrated profile found a nonunique age optimum.")
    bounds, _ = problem.bounds_and_constraint()
    offset = len(problem.free)
    # Zero rate variance is a scientific boundary; mean/upper variance bounds
    # are numerical protections and cannot define a bootstrap generating fit.
    nuisance = parameters[offset:]
    lower = nuisance - bounds.lb[offset:]
    if getattr(problem, "marginal", False) and problem.fixed_sd is None:
        lower = lower[:-1]
    if np.any(lower < 1e-7) or np.any(bounds.ub[offset:] - nuisance < 1e-7):
        raise ValueError("Calibrated profile nuisance reached a numerical bound.")
    return problem, parameters


def _constrained(problem, group, age):
    c = problem.chronology
    if not np.isfinite(age) or not c.lower[group] <= age <= c.upper[group]:
        raise ValueError("Candidate age is outside the hard calibration bounds.")
    lower, upper = c.lower.copy(), c.upper.copy()
    lower[group] = upper[group] = age
    feasible = linprog(
        np.zeros(len(lower)),
        A_ub=-problem.constraints,
        b_ub=np.full(problem.constraints.shape[0], -c.min_duration),
        bounds=list(zip(lower, upper, strict=True)),
        method="highs",
    )
    if not feasible.success:
        raise ValueError("Candidate age is outside the feasible chronology.")
    return _problem_like(
        problem,
        replace(c, lower=lower, upper=upper, initial=feasible.x),
        problem.likelihood,
    )


def _generating_parameters(problem, parameters):
    ages = problem.unpack_ages(parameters)
    if problem.likelihood is None:
        _, mean, _, _, sse = problem.branch_statistics(ages)
        sd = (
            problem.rate_sd
            if problem.rate_sd is not None
            else np.sqrt(sse / len(problem.observation_indices))
        )
    else:
        mean = parameters[len(problem.free)]
        sd = (
            problem.fixed_sd
            if problem.fixed_sd is not None
            else np.sqrt(parameters[-1])
        )
    return ages, float(mean), float(sd)


def _replicate(problem, null_problem, null_parameters, rng, maxiter):
    c = copy.deepcopy(problem.chronology)
    ages, mean, sd = _generating_parameters(null_problem, null_parameters)
    lengths = simulate_rate_lengths(c, ages, mean, sd, problem.rho, rng)
    for edge, length in zip(c.edges, lengths, strict=True):
        edge.dist = float(length)
    c = replace(c, initial=ages.copy())
    likelihood = None
    if problem.likelihood is not None:
        from nwkit.radte_sequence import build_quadratic
        from nwkit.radte_sequence_fit import fit_sequence_model

        exact = simulate_sequence_likelihood(
            problem.likelihood.exact, lengths, rng, chronology=c
        )
        if exact.fit_settings:
            fit_sequence_model(exact, lengths, maxiter=maxiter, **exact.fit_settings)
        likelihood = build_quadratic(
            exact, getattr(exact, "initial_lengths", lengths), maxiter
        )
    return _problem_like(problem, c, likelihood)


def calibrated_profile_test(
    fit, problem, group, age, *, replicates=199, starts=3, maxiter=2000, seed=1
):
    """Test one fixed age; failed refits count as exceeding the observed LR.

    p_upper is a one-sided 99% Monte Carlo upper bound, conditional on the
    plug-in generating parameters. It does not bound nuisance plug-in error.
    """
    if replicates < 19:
        raise ValueError("Calibrated profile requires at least 19 replicates.")
    if problem.likelihood is not None and (
        not getattr(problem, "marginal", False)
        or getattr(problem.likelihood, "exact", None) is None
    ):
        raise ValueError(
            "Calibrated sequence profile requires marginal inference and the native alignment."
        )
    if group not in problem.free:
        raise ValueError("Calibrated profile target must be a free age.")
    settings = dict(starts=starts, maxiter=maxiter, seed=seed)
    null_problem, null_parameters = _solve(
        _constrained(problem, group, age), **settings
    )
    baseline = problem.value_gradient(fit.parameters)[0]
    null_value = null_problem.value_gradient(null_parameters)[0]
    observed = _profile_delta(problem, null_value, baseline, str(group), age)
    rng = np.random.default_rng(seed)
    exceedances, failures = 0, []
    for replicate in range(replicates):
        try:
            sample = _replicate(problem, null_problem, null_parameters, rng, maxiter)
            sample, alternative = _solve(
                sample, **dict(settings, seed=seed + replicate + 1)
            )
            constrained, parameters = _solve(
                _constrained(sample, group, age),
                **dict(settings, seed=seed + replicate + 1),
            )
            statistic = _profile_delta(
                sample,
                constrained.value_gradient(parameters)[0],
                sample.value_gradient(alternative)[0],
                str(group),
                age,
            )
            exceedances += statistic >= observed - 1e-10
        except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
            # Dropping failed fits would select the bootstrap distribution.
            exceedances += 1
            failures.append(str(exc))
    upper = (
        binomtest(exceedances, replicates, alternative="less")
        .proportion_ci(confidence_level=0.99, method="exact")
        .high
    )
    _, mean, sd = _generating_parameters(null_problem, null_parameters)
    return dict(
        age=float(age * problem.chronology.scale),
        statistic=float(observed),
        p_value=(1 + exceedances) / (1 + replicates),
        p_upper=float(upper),
        replicates=replicates,
        exceedances=int(exceedances),
        failures=len(failures),
        first_failure=failures[0] if failures else None,
        null_rate_sd=sd,
        null_log_rate_mean=mean - np.log(problem.chronology.scale),
    )


def calibrated_profile_intervals(
    fit,
    problem,
    *,
    level=0.95,
    replicates=199,
    grid_points=17,
    starts=3,
    maxiter=2000,
    seed=1,
):
    """Return an explicitly experimental, padded hull of an age acceptance grid."""
    fit.interval_lower = fit.interval_upper = None
    fit.calibration_profile = []
    if not 0 < level < 1 or grid_points < 5:
        raise ValueError(
            "Interior interval level and at least five grid points required."
        )
    from nwkit.radte_exact_interval import exact_log_duration_intervals

    if exact_log_duration_intervals(fit, problem, level):
        return
    if replicates < int(np.ceil(np.log(0.01) / np.log(level))):
        raise ValueError(
            "Too few replicates for the Monte Carlo upper bound to reject at this level."
        )
    if not len(problem.free):
        fit.interval_lower, fit.interval_upper = fit.ages.copy(), fit.ages.copy()
        fit.interval_status = "fixed-ages"
        return
    if problem.likelihood is None and fit.log_rate_sd == 0:
        fit.interval_status = "unavailable-strict-clock-limit"
        return
    lower, upper = chronology_domain(problem.chronology)
    result_lower, result_upper = fit.ages.copy(), fit.ages.copy()
    for group in problem.free:
        # Include the fitted age; evaluate the entire domain, not just the
        # first threshold crossing. Pad by adjacent grid cells when reporting
        # a hull. Finite grids can still miss components: this is not exact CI.
        grid = np.unique(
            np.append(
                np.linspace(lower[group], upper[group], grid_points), fit.ages[group]
            )
        )
        accepted = []
        for index, age in enumerate(grid):
            try:
                row = calibrated_profile_test(
                    fit,
                    problem,
                    group,
                    age,
                    replicates=replicates,
                    starts=starts,
                    maxiter=maxiter,
                    seed=seed,
                )
                if row["failures"] > replicates // 10:
                    raise ValueError("More than 10% of calibration refits failed.")
            except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
                fit.calibration_profile.append(
                    dict(
                        group=problem.chronology.groups[group],
                        age=float(age * problem.chronology.scale),
                        error=str(exc),
                    )
                )
                fit.interval_status = "unavailable-calibrated-profile-fit-failure"
                return
            row["group"] = problem.chronology.groups[group]
            fit.calibration_profile.append(row)
            if row["p_upper"] > 1 - level:
                accepted.append(index)
        if not accepted:
            fit.interval_status = "unavailable-empty-calibrated-profile-grid"
            return
        result_lower[group] = grid[max(0, min(accepted) - 1)]
        result_upper[group] = grid[min(len(grid) - 1, max(accepted) + 1)]
        fit.diagnostics.append(
            f"calibrated_profile_max_grid_spacing:{problem.chronology.groups[group]}={np.max(np.diff(grid)) * problem.chronology.scale:.12g}"
        )
    fit.interval_lower, fit.interval_upper = result_lower, result_upper
    fit.interval_status = "experimental-calibrated-profile-grid-hull"
    fit.diagnostics.extend(
        [
            "calibrated_profile_nuisance=constrained-plugin",
            "calibrated_profile_mc_upper_confidence=0.99",
            "calibrated_profile_grid_hull_is_not_an_exact_confidence_interval",
        ]
    )
