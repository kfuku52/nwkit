"""Profile and bootstrap uncertainty for the explicitly conditional RADTE model."""

from dataclasses import replace

import numpy as np
from scipy.optimize import brentq, linprog
from scipy.stats import chi2

from nwkit.radte_model import DatingOptimizationError, fit_dates, solve_problem


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


def _record_profile_attempts(fit, attempts, group, age, *, phase="profile"):
    successful = [
        a["objective"]
        for a in attempts
        if a["success"] and not a.get("nonunique_age_solution")
    ]
    spread = max(successful) - min(successful) if successful else None
    for attempt in attempts:
        fit.attempts.append(
            dict(
                attempt,
                phase=phase,
                profile_group=group,
                profile_age=age,
                successful_starts=len(successful),
                objective_spread=spread,
            )
        )
    flags = []
    if spread is not None and spread > 1e-4 * max(1, abs(min(successful))):
        flags.append("multiple_local_optima")
    if any(a.get("nonunique_age_solution") for a in attempts):
        flags.append("nonunique_age_optimum")
    for flag in flags:
        message = f"profile_{flag}:{group}"
        if message not in fit.diagnostics:
            fit.diagnostics.append(message)


def _profile_delta(problem, value, baseline, group, age):
    # A constrained fit cannot improve the unconstrained optimum. Stop rather
    # than build an interval around a demonstrably inferior reference fit.
    objective_scale = max(1.0, abs(baseline))
    if problem.likelihood is None:
        # This objective is a sum of squared log-rate residuals, not a log
        # likelihood. A unit-sized absolute tolerance can hide many likelihood
        # ratio units when the residuals and specified rate variance are small.
        objective_scale = max(abs(baseline), (problem.rate_sd or 0.0) ** 2)
    tolerance = 1e-7 * objective_scale
    if value < baseline - tolerance:
        raise ValueError(
            f"Profile found a better optimum for {group} at age {age:g} "
            f"(objective {value:.12g} < baseline {baseline:.12g}). "
            "Refit the point estimate with more --starts before requesting intervals."
        )
    value = max(value, baseline)  # Only roundoff-sized negative differences remain.
    if problem.likelihood is not None:
        return 2 * (value - baseline)
    if problem.rate_sd is not None:
        return 2 * (value - baseline) / problem.rate_sd**2
    if baseline <= 0:
        raise ValueError("Profile likelihood is undefined at zero rate variance.")
    return len(problem.chronology.edges) * np.log(value / baseline)


def _continue_profile(evaluate, target, cache, group, scale):
    """Retry numerical failures through closer feasible fits, with a hard budget."""
    candidate = target
    failures = 0
    while True:
        try:
            result = evaluate(candidate)
        except DatingOptimizationError as exc:
            failures += 1
            nearest = min(cache, key=lambda value: abs(value - candidate))
            halfway = nearest + (candidate - nearest) / 2
            if failures >= 32 or halfway in (nearest, candidate):
                raise DatingOptimizationError(
                    f"Profile continuation failed for {group} at age "
                    f"{target * scale:g} after {failures} step reductions: {exc}",
                    exc.attempts,
                ) from exc
            candidate = halfway
            continue
        if candidate == target:
            return result
        candidate = target


def profile_intervals(fit, problem, *, level=0.95, starts=3, maxiter=2000, seed=1):
    if not 0 < level < 1:
        raise ValueError("Interval level must be between zero and one.")
    if starts < 1:
        raise ValueError("Profile optimizer starts must be positive.")
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
            feasible_ends.append(extreme.x.copy())

        # Every cached age vector is feasible. Convex interpolation with an
        # extreme point preserves all linear chronology constraints.
        cache = {
            float(fit.ages[group]): (
                fit.ages.copy(),
                fit.parameters[len(problem.free) :].copy(),
                -threshold,
            )
        }

        def evaluate(age, group=group, cache=cache, feasible_ends=feasible_ends):
            if age in cache:
                return cache[age][2]
            nearest = min(cache, key=lambda value: abs(value - age))
            previous, nuisance, _ = cache[nearest]
            extreme = feasible_ends[0 if age < nearest else 1]
            weight = (age - nearest) / (extreme[group] - nearest)
            initial = (1 - weight) * previous + weight * extreme
            initial[group] = age
            lo, hi = c.lower.copy(), c.upper.copy()
            lo[group] = hi[group] = age
            constrained = replace(c, lower=lo, upper=hi, initial=initial)
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
            start = np.concatenate([initial[profile.free], nuisance])
            phase = "profile"
            try:
                x, attempts = solve_problem(
                    profile,
                    starts=starts,
                    maxiter=maxiter,
                    seed=seed,
                    initial_parameters=start,
                )
                _record_profile_attempts(fit, attempts, c.groups[group], age * c.scale)
                if getattr(profile, "marginal", False):
                    from nwkit.radte_marginal import refine_quadrature

                    phase = "profile-quadrature"
                    profile, x, attempts = refine_quadrature(
                        profile, x, starts=starts, maxiter=maxiter, seed=seed
                    )
                    _record_profile_attempts(
                        fit, attempts, c.groups[group], age * c.scale, phase=phase
                    )
            except DatingOptimizationError as exc:
                _record_profile_attempts(
                    fit, exc.attempts, c.groups[group], age * c.scale, phase=phase
                )
                raise
            validate_profile_approximation(profile, x)
            value = profile.value_gradient(x)[0]
            delta = _profile_delta(
                problem, value, baseline, c.groups[group], age * c.scale
            )
            result = float(delta - threshold)
            cache[age] = (
                profile.unpack_ages(x).copy(),
                x[len(profile.free) :].copy(),
                result,
            )
            return result

        def difference(age, evaluate=evaluate, cache=cache, group=group):
            return _continue_profile(evaluate, age, cache, c.groups[group], c.scale)

        for destination, extreme in zip((lower, upper), feasible_ends, strict=True):
            endpoint = float(extreme[group])
            best_age = float(fit.ages[group])
            if abs(endpoint - best_age) < 1e-9:
                destination[group] = endpoint
                boundary_limited = True
                continue
            previous = best_age
            # Continue outward from the optimum. Stop as soon as the LR
            # threshold is bracketed, without visiting a singular boundary.
            for fraction in (
                1 / 32,
                1 / 16,
                1 / 8,
                1 / 4,
                1 / 2,
                3 / 4,
                7 / 8,
                15 / 16,
                31 / 32,
                1 - 1e-7,
            ):
                candidate = best_age + fraction * (endpoint - best_age)
                if difference(candidate) > 0:
                    destination[group] = brentq(
                        difference, *sorted([previous, candidate]), xtol=1e-6
                    )
                    break
                previous = candidate
            else:
                destination[group] = endpoint
                boundary_limited = True
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
