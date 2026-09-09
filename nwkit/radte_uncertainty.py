"""Profile and bootstrap uncertainty for the explicitly conditional RADTE model."""

from dataclasses import replace

import numpy as np
from scipy.optimize import brentq, linprog
from scipy.stats import chi2

from nwkit.radte_model import fit_dates, solve_problem


def validate_profile_approximation(problem, parameters):
    likelihood = problem.likelihood
    if not hasattr(likelihood, "mapping"):
        return
    if likelihood.exact is None:
        lengths = problem.posterior_rates(parameters) * problem.chronology.durations(
            problem.unpack_ages(parameters)
        )
        combined = np.bincount(likelihood.mapping, weights=lengths)
        valid = np.max(abs(np.log(combined) - likelihood.center)) <= 0.5
    elif getattr(problem, "marginal", False):
        valid = all(check[0] for check in problem.approximation_checks(parameters))
    else:
        lengths = problem.posterior_rates(parameters) * problem.chronology.durations(
            problem.unpack_ages(parameters)
        )
        valid = likelihood.check(lengths)[0]
    if not valid:
        raise ValueError(
            "Profile interval left the validated quadratic likelihood region; use exact likelihood or bootstrap."
        )


def profile_intervals(fit, problem, *, level=0.95, maxiter=2000, seed=1):
    if not 0 < level < 1:
        raise ValueError("Interval level must be between zero and one.")
    c = problem.chronology
    if fit.log_rate_sd == 0:
        fit.interval_status = "unavailable-strict-clock-limit"
        return
    lower, upper = fit.ages.copy(), fit.ages.copy()
    threshold = float(chi2.ppf(level, 1))
    baseline = problem.value_gradient(fit.parameters)[0]
    # One-sided bounds are reported as the feasible endpoint; they are not
    # falsely extrapolated beyond the hard calibration domain.
    boundary_limited = False
    for group in problem.free:
        coordinate = np.zeros(len(c.groups))
        coordinate[group] = 1
        feasible_ends = []
        for sign in [1, -1]:
            extreme = linprog(
                sign * coordinate,
                A_ub=-problem.constraints,
                b_ub=np.full(problem.constraints.shape[0], -c.min_duration),
                bounds=list(zip(c.lower, c.upper, strict=True)),
                method="highs",
            )
            if not extreme.success:
                raise ValueError("Could not resolve the profile's feasible age domain.")
            feasible_ends.append(float(extreme.x[group]))

        def difference(age, group=group):
            lo, hi = c.lower.copy(), c.upper.copy()
            lo[group] = hi[group] = age
            constrained = replace(c, lower=lo, upper=hi, initial=fit.ages.copy())
            constrained.initial[group] = age
            extra = {}
            if getattr(problem, "marginal", False):
                extra["quadrature_points"] = problem.quadrature_points
            profile = type(problem)(
                constrained,
                rho=problem.rho,
                likelihood=problem.likelihood,
                rate_sd=problem.rate_sd,
                **extra,
            )
            x, _ = solve_problem(profile, starts=2, maxiter=maxiter, seed=seed)
            if getattr(profile, "marginal", False):
                from nwkit.radte_marginal import refine_quadrature

                profile, x, _ = refine_quadrature(
                    profile, x, starts=2, maxiter=maxiter, seed=seed
                )
            validate_profile_approximation(profile, x)
            value = profile.value_gradient(x)[0]
            if problem.likelihood is None and problem.rate_sd is not None:
                delta = 2 * (value - baseline) / problem.rate_sd**2
            elif problem.likelihood is None:
                if baseline <= 0:
                    raise ValueError(
                        "Profile likelihood is undefined at zero rate variance."
                    )
                delta = len(c.edges) * np.log(max(value, baseline) / baseline)
            else:
                delta = 2 * (value - baseline)
            return float(delta - threshold)

        for destination, endpoint in [
            (lower, feasible_ends[0]),
            (upper, feasible_ends[1]),
        ]:
            best_age = fit.ages[group]
            if abs(endpoint - best_age) < 1e-9:
                destination[group] = endpoint
                boundary_limited = True
                continue
            # Avoid numerical evaluation exactly on a positive-duration limit.
            inside = endpoint + (best_age - endpoint) * 1e-7
            if difference(inside) <= 0:
                destination[group] = endpoint
                boundary_limited = True
            else:
                destination[group] = brentq(
                    difference, *sorted([inside, best_age]), xtol=1e-6
                )
    fit.interval_lower, fit.interval_upper = lower, upper
    fit.interval_status = "conditional-profile"
    if boundary_limited:
        fit.interval_status += "-calibration-limited"


def bootstrap_intervals(
    fit,
    problem,
    *,
    replicates=100,
    level=0.95,
    rho=0,
    starts=1,
    maxiter=2000,
    seed=1,
    rate_sd=None,
):
    if replicates < 20:
        raise ValueError("Bootstrap intervals require at least 20 replicates.")
    if not 0 < level < 1:
        raise ValueError("Interval level must be between zero and one.")
    c = problem.chronology
    exact = getattr(problem.likelihood, "exact", problem.likelihood)
    if problem.likelihood is None and fit.log_rate_sd == 0:
        fit.interval_status = "unavailable-strict-clock-limit"
        return
    if problem.likelihood is not None and exact is None:
        raise ValueError(
            "Sequence bootstrap requires the original alignment, not only its likelihood summary."
        )
    rng = np.random.default_rng(seed)
    original = np.array([n.dist for n in c.edges])
    samples, failures = [], []
    try:
        for replicate in range(replicates):
            if exact is None:
                # Stationary Gaussian AR(1) on the actual gene genealogy.
                values = {c.gene: rng.normal(0, fit.log_rate_sd)}
                lengths = []
                duration = c.durations(fit.ages)
                for i, node in enumerate(c.edges):
                    values[node] = rho * values[node.up] + rng.normal(
                        0, fit.log_rate_sd * np.sqrt(1 - rho**2)
                    )
                    lengths.append(
                        duration[i] * c.scale * np.exp(fit.log_rate_mean + values[node])
                    )
                for node, length in zip(c.edges, lengths, strict=True):
                    node.dist = length
                likelihood = None
            else:
                likelihood = exact.bootstrap(rng)
            try:
                initial_lengths = original
                if likelihood is not None and getattr(likelihood, "fit_settings", {}):
                    from nwkit.radte_sequence_fit import fit_sequence_model

                    fit_sequence_model(
                        likelihood, original, maxiter=maxiter, **likelihood.fit_settings
                    )
                    initial_lengths = getattr(likelihood, "initial_lengths", original)
                if getattr(problem, "marginal", False):
                    from nwkit.radte_sequence import build_quadratic

                    likelihood = build_quadratic(likelihood, initial_lengths, maxiter)
                sample, sample_problem = fit_dates(
                    c,
                    rho=rho,
                    likelihood=likelihood,
                    rate_sd=rate_sd,
                    starts=starts,
                    maxiter=maxiter,
                    seed=seed + replicate + 1,
                    inference="marginal"
                    if getattr(problem, "marginal", False)
                    else "joint-map",
                )
                if getattr(problem, "marginal", False) and not all(
                    check[0]
                    for check in sample_problem.approximation_checks(sample.parameters)
                ):
                    raise ValueError(
                        "Bootstrap quadratic approximation failed exact validation."
                    )
                samples.append(sample.ages)
            except (ValueError, FloatingPointError) as exc:
                failures.append(str(exc))
    finally:
        for node, length in zip(c.edges, original, strict=True):
            node.dist = length
    fit.diagnostics.append(f"bootstrap_successes={len(samples)}/{replicates}")
    if failures:
        fit.diagnostics.append("bootstrap_first_failure=" + failures[0])
    if len(samples) < max(20, int(np.ceil(0.9 * replicates))):
        fit.interval_status = "unavailable-too-many-bootstrap-failures"
        return
    fit.samples = np.asarray(samples)
    tail = (1 - level) / 2
    fit.interval_lower, fit.interval_upper = np.quantile(
        fit.samples, [tail, 1 - tail], axis=0
    )
    fit.interval_status = (
        "site-bootstrap" if exact else "conditional-parametric-rate-bootstrap"
    )
