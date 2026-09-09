"""Native shared-age relaxed clocks, with explicit conditional inference.

Tree-only mode integrates Gaussian log rates analytically and profiles their
mean/variance. Sequence mode fits a penalized likelihood conditional on the
tree-only estimate of rate variance (an empirical-Bayes MAP, not MCMC).
"""

from dataclasses import dataclass, field

import numpy as np
from scipy import sparse
from scipy.optimize import Bounds, LinearConstraint, linprog, minimize
from scipy.stats import norm


@dataclass
class DatingFit:
    ages: np.ndarray
    rates: np.ndarray
    log_rate_mean: float
    log_rate_sd: float
    objective: float
    parameters: np.ndarray
    converged: bool
    attempts: list[dict]
    diagnostics: list[str] = field(default_factory=list)
    interval_lower: np.ndarray | None = None
    interval_upper: np.ndarray | None = None
    interval_status: str = "not-requested"
    samples: np.ndarray | None = None
    sample_ids: list[int] | None = None
    sample_presence: np.ndarray | None = None
    sample_clade_presence: np.ndarray | None = None
    ensemble_metadata: dict | None = None
    conditional_intervals: list[dict] = field(default_factory=list)


def rate_precision(chronology, rho=0.0):
    """Stationary Gaussian AR(1) on log rates; integrate the root's rate.

    rho=0 gives independent lognormal branch rates. For rho>0, correlation
    decays per gene-tree edge (not per unit time). No rates are shared by species.
    """
    count = len(chronology.edges)
    if not 0 <= rho < 1:
        raise ValueError("Rate correlation must satisfy 0 <= rho < 1.")
    if rho == 0:
        return sparse.eye(count, format="csr"), 0.0
    positions = {n: i + 1 for i, n in enumerate(chronology.edges)}
    positions[chronology.gene] = 0
    diagonal = np.zeros(count + 1)
    diagonal[0] = 1
    rows, cols, values = [], [], []
    innovation = 1 - rho**2
    for node in chronology.edges:
        i, p = positions[node], positions[node.up]
        diagonal[i] += 1 / innovation
        diagonal[p] += rho**2 / innovation
        rows.extend([i, p])
        cols.extend([p, i])
        values.extend([-rho / innovation] * 2)
    rows.extend(range(count + 1))
    cols.extend(range(count + 1))
    values.extend(diagonal)
    full = sparse.csr_matrix((values, (rows, cols)), shape=(count + 1, count + 1))
    qroot = full[1:, :1]
    precision = full[1:, 1:] - qroot @ qroot.T / diagonal[0]
    logdet = count * np.log(innovation) + np.log(diagonal[0])
    return precision.tocsr(), float(logdet)


class DatingProblem:
    def __init__(
        self,
        chronology,
        *,
        rho=0.0,
        likelihood=None,
        rate_sd=None,
        omit_root_observations=False,
    ):
        self.chronology = chronology
        self.rho = rho
        self.free = np.flatnonzero(chronology.upper > chronology.lower)
        self.fixed = chronology.lower.copy()
        self.precision, self.logdet = rate_precision(chronology, rho)
        self.observation_indices = np.array(
            [
                i
                for i, node in enumerate(chronology.edges)
                if not omit_root_observations or node.up is not chronology.gene
            ]
        )
        if omit_root_observations:
            missing = np.array(
                [
                    i
                    for i, node in enumerate(chronology.edges)
                    if node.up is chronology.gene
                ]
            )
            cross = self.precision[self.observation_indices][:, missing].toarray()
            missing_precision = self.precision[missing][:, missing].toarray()
            marginal = self.precision[self.observation_indices][
                :, self.observation_indices
            ].toarray()
            marginal -= cross @ np.linalg.solve(missing_precision, cross.T)
            self.precision = sparse.csr_matrix(marginal)
            self.logdet = -float(np.linalg.slogdet(marginal)[1])
        self.ones = np.ones(len(self.observation_indices))
        self.qones = self.precision @ self.ones
        self.mean_precision = self.ones @ self.qones
        self.observed = np.array([n.dist for n in chronology.edges], dtype=float)
        self.log_observed = np.log(self.observed)
        self.likelihood = likelihood
        self.rate_sd = rate_sd
        self.rate_variance_estimated = rate_sd is None
        self.age_design = sparse.coo_matrix(
            (
                np.tile([1.0, -1.0], len(chronology.edges)),
                (
                    np.repeat(np.arange(len(chronology.edges)), 2),
                    np.column_stack([chronology.parent, chronology.child]).ravel(),
                ),
            ),
            shape=(len(chronology.edges), len(chronology.groups)),
        ).tocsr()
        self.constraints = sparse.coo_matrix(
            (
                np.tile([1.0, -1.0], len(chronology.constraint_parent)),
                (
                    np.repeat(np.arange(len(chronology.constraint_parent)), 2),
                    np.column_stack(
                        [chronology.constraint_parent, chronology.constraint_child]
                    ).ravel(),
                ),
            ),
            shape=(len(chronology.constraint_parent), len(chronology.groups)),
        ).tocsr()
        self.free_design = self.age_design[:, self.free]
        self.free_constraints = self.constraints[:, self.free]
        fixed_ages = self.fixed.copy()
        fixed_ages[self.free] = 0
        self.constraint_lower = chronology.min_duration - self.constraints @ fixed_ages
        active = np.asarray(self.free_constraints.getnnz(axis=1) > 0)
        self.free_constraints = self.free_constraints[active]
        self.constraint_lower = self.constraint_lower[active]

    def unpack_ages(self, x):
        ages = self.fixed.copy()
        ages[self.free] = x[: len(self.free)]
        return ages

    def branch_statistics(self, ages):
        duration = self.chronology.durations(ages)
        log_rates = (self.log_observed - np.log(duration))[self.observation_indices]
        mean = float(self.qones @ log_rates / self.mean_precision)
        residual = log_rates - mean
        qres = self.precision @ residual
        sse = float(residual @ qres)
        return duration, mean, residual, qres, max(0.0, sse)

    def value_gradient(self, x):
        ages = self.unpack_ages(x)
        duration = self.chronology.durations(ages)
        if np.any(duration <= 0):
            # SLSQP can evaluate outside its linear feasible region. A smooth
            # exterior penalty guides it back; it can never pass final QA.
            violation = np.minimum(duration - self.chronology.min_duration, 0)
            grad = np.zeros_like(x)
            grad[: len(self.free)] = self.free_design.T @ (2e12 * violation)
            return float(1e12 * (1 + violation @ violation)), grad
        if self.likelihood is None:
            _, _, _, qres, sse = self.branch_statistics(ages)
            gradient = self.free_design[self.observation_indices].T @ (
                -qres / duration[self.observation_indices]
            )
            # Profiling the normal mean and variance makes minimizing SSE
            # equivalent to maximizing the marginal log-rate likelihood.
            return 0.5 * sse, np.asarray(gradient)
        offset = len(self.free)
        mu = x[offset]
        deviations = x[offset + 1 :] if self.rate_sd > 0 else np.zeros_like(duration)
        lengths = np.exp(mu + deviations) * duration
        nll, branch_gradient = self.likelihood.value_gradient(lengths)
        log_length_gradient = branch_gradient * lengths
        age_gradient = self.free_design.T @ (log_length_gradient / duration)
        if self.rate_sd > 0:
            prior_gradient = self.precision @ deviations / self.rate_sd**2
            nll += 0.5 * deviations @ prior_gradient
            gradient = np.concatenate(
                [
                    age_gradient,
                    [log_length_gradient.sum()],
                    log_length_gradient + prior_gradient,
                ]
            )
        else:
            gradient = np.concatenate([age_gradient, [log_length_gradient.sum()]])
        return float(nll), gradient

    def bounds_and_constraint(self):
        lower = self.chronology.lower[self.free].tolist()
        upper = self.chronology.upper[self.free].tolist()
        extras = 0
        if self.likelihood is not None:
            # Bounds protect exponentials, not time calibration semantics.
            lower.append(-30.0)
            upper.append(30.0)
            extras = 1
            if self.rate_sd > 0:
                lower.extend([-30.0] * len(self.observed))
                upper.extend([30.0] * len(self.observed))
                extras += len(self.observed)
        matrix = sparse.hstack(
            [
                self.free_constraints,
                sparse.csr_matrix((len(self.constraint_lower), extras)),
            ]
        ).toarray()
        constraint = LinearConstraint(matrix, self.constraint_lower, np.inf)
        return Bounds(lower, upper), constraint

    def initial_parameters(self, initial_ages=None):
        ages = self.chronology.initial if initial_ages is None else initial_ages
        x = ages[self.free]
        if self.likelihood is not None:
            duration, mu, residual, _, _ = self.branch_statistics(ages)
            x = np.concatenate([x, [mu], residual if self.rate_sd > 0 else []])
        return x

    def feasible(self, x):
        ages = self.unpack_ages(x)
        c = self.chronology
        return (
            np.isfinite(x).all()
            and np.all(ages >= c.lower - 1e-10)
            and np.all(ages <= c.upper + 1e-10)
            and np.all(self.constraints @ ages >= c.min_duration * 0.5)
        )

    def posterior_rates(self, x):
        offset = len(self.free)
        deviations = (
            x[offset + 1 :] if self.rate_sd > 0 else np.zeros(len(self.observed))
        )
        return np.exp(x[offset] + deviations)


class DatingOptimizationError(ValueError):
    """Numerical solver failure with the complete attempted fits attached."""

    def __init__(self, message, attempts):
        super().__init__(message)
        self.attempts = attempts


def solve_problem(
    problem, *, starts=3, maxiter=2000, seed=1, initial=None, initial_parameters=None
):
    if starts < 1 or maxiter < 1:
        raise ValueError("Optimizer starts and max iterations must be positive.")
    x0 = (
        problem.initial_parameters(initial)
        if initial_parameters is None
        else initial_parameters.copy()
    )
    if len(x0) == 0:
        return x0, [
            dict(success=True, objective=problem.value_gradient(x0)[0], iterations=0)
        ]
    bounds, constraint = problem.bounds_and_constraint()
    rng = np.random.default_rng(seed)
    attempts, solutions = [], []
    for attempt in range(starts):
        start = x0.copy()
        if attempt and len(problem.free):
            vertex = linprog(
                rng.normal(size=len(problem.free)),
                A_ub=-problem.free_constraints,
                b_ub=-problem.constraint_lower,
                bounds=list(
                    zip(
                        bounds.lb[: len(problem.free)],
                        bounds.ub[: len(problem.free)],
                        strict=True,
                    )
                ),
                method="highs",
            )
            if vertex.success:
                start[: len(problem.free)] = (
                    0.7 * x0[: len(problem.free)] + 0.3 * vertex.x
                )
        if (
            attempt
            and problem.likelihood is not None
            and not getattr(problem, "marginal", False)
            and problem.rate_sd > 0
        ):
            start[len(problem.free) + 1 :] += rng.normal(0, 0.1, len(problem.observed))
        result = minimize(
            problem.value_gradient,
            start,
            jac=True,
            method="SLSQP",
            bounds=bounds,
            constraints=[constraint] if constraint.A.shape[0] else [],
            options={"maxiter": maxiter, "ftol": 1e-10},
        )
        success = bool(
            result.success and problem.feasible(result.x) and np.isfinite(result.fun)
        )
        attempts.append(
            dict(
                success=success,
                objective=float(result.fun),
                iterations=int(result.nit),
                message=str(result.message),
            )
        )
        if success:
            solutions.append(result)
    if not solutions:
        raise DatingOptimizationError(
            "Dating optimization failed; no constraints were dropped. "
            + "; ".join(str(a.get("message", "")) for a in attempts),
            attempts,
        )
    best = min(solutions, key=lambda res: res.fun)
    close = [
        res
        for res in solutions
        if abs(res.fun - best.fun) < 1e-7 * max(1, abs(best.fun))
    ]
    if any(
        np.max(abs(res.x[: len(problem.free)] - best.x[: len(problem.free)]), initial=0)
        > 1e-3
        for res in close
    ):
        attempts.append(
            dict(
                success=True,
                objective=float(best.fun),
                iterations=0,
                nonunique_age_solution=True,
            )
        )
    return best.x, attempts


def solution_diagnostics(problem, x, ages, attempts, starts):
    chronology = problem.chronology
    diagnostics = []
    values = [a["objective"] for a in attempts[-starts:] if a["success"]]
    if len(values) > 1 and max(values) - min(values) > 1e-4 * max(1, abs(min(values))):
        diagnostics.append("multiple_local_optima")
    if any(a.get("nonunique_age_solution", False) for a in attempts):
        diagnostics.append("nonunique_age_optimum")
    positive = ages > 0
    scale_lower = max(
        float(np.max(chronology.lower[positive] / ages[positive])),
        float(np.max(chronology.min_duration / (problem.constraints @ ages))),
    )
    scale_upper = float(np.min(chronology.upper[positive] / ages[positive]))
    if (
        scale_lower <= 1 + 1e-7
        and scale_upper >= 1 - 1e-7
        and scale_upper - scale_lower > 1e-6
    ):
        diagnostics.append("absolute_age_scale_unidentified_within_hard_bounds")
    numerical_bounds, _ = problem.bounds_and_constraint()
    offset = len(problem.free)
    if len(x) > offset and (
        np.any(x[offset:] - numerical_bounds.lb[offset:] < 1e-5)
        or np.any(numerical_bounds.ub[offset:] - x[offset:] < 1e-5)
    ):
        diagnostics.append("nuisance_parameter_at_numerical_bound")
    return diagnostics


def fit_dates(
    chronology,
    *,
    rho=0.0,
    likelihood=None,
    rate_sd=None,
    starts=3,
    maxiter=2000,
    seed=1,
    inference="auto",
):
    from nwkit.radte_sequence import QuadraticLikelihood

    if inference not in {"auto", "marginal", "joint-map"}:
        raise ValueError("Inference must be auto, marginal, or joint-map.")
    if (
        likelihood is not None
        and inference == "marginal"
        and not isinstance(likelihood, QuadraticLikelihood)
    ):
        raise ValueError(
            "Marginal sequence inference requires a validated quadratic likelihood."
        )
    omit_root = likelihood is not None and len(chronology.edges) > 2
    if likelihood is not None and len(chronology.edges) == 2 and rate_sd is None:
        raise ValueError(
            "Two-tip sequence dating cannot estimate rate variance; supply --rate-sd."
        )
    branch_problem = DatingProblem(
        chronology, rho=rho, rate_sd=rate_sd, omit_root_observations=omit_root
    )
    reference_likelihood = getattr(likelihood, "exact", likelihood)
    if isinstance(likelihood, QuadraticLikelihood):
        # Use the same identifiable ML lengths before and after serialization.
        # The root split is unobserved and omitted from variance estimation.
        multiplicity = np.bincount(likelihood.mapping)
        branch_problem.observed = (
            np.exp(likelihood.center[likelihood.mapping])
            / multiplicity[likelihood.mapping]
        )
        branch_problem.log_observed = np.log(branch_problem.observed)
    elif hasattr(reference_likelihood, "initial_lengths"):
        branch_problem.observed = reference_likelihood.initial_lengths.copy()
        branch_problem.log_observed = np.log(branch_problem.observed)
    x, attempts = solve_problem(
        branch_problem, starts=starts, maxiter=maxiter, seed=seed
    )
    final_attempts = attempts
    for attempt in attempts:
        attempt["phase"] = "branch-initialization" if likelihood else "inference"
    ages = branch_problem.unpack_ages(x)
    duration, mu, _, _, sse = branch_problem.branch_statistics(ages)
    inferred_sd = float(np.sqrt(sse / len(branch_problem.observation_indices)))
    sd = inferred_sd if rate_sd is None else float(rate_sd)
    if not np.isfinite(sd) or sd < 0:
        raise ValueError("--rate-sd must be finite and nonnegative.")
    if likelihood is None and rate_sd == 0 and inferred_sd >= 1e-7:
        raise ValueError(
            "--rate-sd 0 is incompatible with the input branch lengths: "
            "no strict-clock fit satisfies the chronology constraints."
        )
    problem = branch_problem
    diagnostics = []
    if likelihood is None and rate_sd is not None:
        diagnostics.append("rate_sd_changes_uncertainty_only_in_tree_mode")
    if sd < 1e-7:
        sd = 0.0
        diagnostics.append("strict_clock_limit")
    objective = 0.5 * sse
    rates = branch_problem.observed / duration
    if likelihood is not None:
        marginal = (
            isinstance(likelihood, QuadraticLikelihood) and inference != "joint-map"
        )
        if marginal:
            from nwkit.radte_marginal import MarginalDatingProblem

            problem = MarginalDatingProblem(
                chronology, rho=rho, likelihood=likelihood, rate_sd=rate_sd
            )
        else:
            problem = DatingProblem(
                chronology, rho=rho, likelihood=likelihood, rate_sd=sd
            )
        x, sequence_attempts = solve_problem(
            problem, starts=starts, maxiter=maxiter, seed=seed, initial=ages
        )
        final_attempts = sequence_attempts
        for attempt in sequence_attempts:
            attempt["phase"] = "inference"
        attempts.extend(sequence_attempts)
        if marginal:
            from nwkit.radte_marginal import refine_quadrature

            problem, x, quadrature_attempts = refine_quadrature(
                problem, x, starts=starts, maxiter=maxiter, seed=seed
            )
            for attempt in quadrature_attempts:
                attempt["phase"] = "quadrature-refinement"
            if quadrature_attempts:
                final_attempts = quadrature_attempts
            attempts.extend(quadrature_attempts)
        ages = problem.unpack_ages(x)
        mu = float(x[len(problem.free)])
        if marginal:
            sd = float(np.exp(x[-1])) if rate_sd is None else float(rate_sd)
        rates = problem.posterior_rates(x)
        objective = problem.value_gradient(x)[0]
        diagnostics.append(
            "sequence_rates_marginalized"
            if marginal
            else "sequence_map_conditional_on_rate_variance_and_substitution_model"
        )
    else:
        diagnostics.append("conditional_on_input_branch_lengths_and_root_split")
    diagnostics.extend(solution_diagnostics(problem, x, ages, final_attempts, starts))
    if sd > 1e-7 and "strict_clock_limit" in diagnostics:
        diagnostics.remove("strict_clock_limit")
    fit = DatingFit(
        ages,
        rates / chronology.scale,
        mu - np.log(chronology.scale),
        sd,
        objective,
        x,
        True,
        attempts,
        diagnostics,
    )
    # The exact sequence problem receives the fitted SD as a numerical input;
    # retain whether that input was estimated or supplied by the caller.
    problem.rate_variance_estimated = rate_sd is None
    return fit, problem


def curvature_covariance(fit, problem, level=0.95):
    """Return the free-age covariance, retaining all nuisance directions.

    Unavailable curvature and fixed ages set the interval status and return None.
    """
    if not 0 < level < 1:
        raise ValueError("Interval level must be between zero and one.")
    fit.interval_lower = fit.interval_upper = None
    c = problem.chronology
    lower, upper = fit.ages.copy(), fit.ages.copy()
    if len(problem.free) == 0:
        fit.interval_lower, fit.interval_upper = lower, upper
        fit.interval_status = "fixed-ages"
        return
    x = fit.parameters
    age = fit.ages[problem.free]
    numerical_bounds, _ = problem.bounds_and_constraint()
    if (
        np.any(age - c.lower[problem.free] < 1e-5)
        or np.any(c.upper[problem.free] - age < 1e-5)
        or np.any(problem.constraints @ fit.ages < 1e-5)
        or np.any(x - numerical_bounds.lb < 1e-5)
        or np.any(numerical_bounds.ub - x < 1e-5)
    ):
        fit.interval_status = "unavailable-active-bound-use-profile-or-bootstrap"
        return
    if fit.log_rate_sd == 0:
        fit.interval_status = "unavailable-strict-clock-limit"
        return
    if len(x) > 600:
        fit.interval_status = "unavailable-dense-hessian-limit-use-bootstrap"
        return
    hessian = np.empty((len(x), len(x)))
    for j in range(len(x)):
        step = 1e-5 * max(1.0, abs(x[j]))
        plus, minus = x.copy(), x.copy()
        plus[j] += step
        minus[j] -= step
        if not problem.feasible(plus) or not problem.feasible(minus):
            fit.interval_status = "unavailable-near-constraint"
            return
        hessian[:, j] = (
            problem.value_gradient(plus)[1] - problem.value_gradient(minus)[1]
        ) / (2 * step)
    hessian = (hessian + hessian.T) / 2
    if problem.likelihood is None:
        hessian /= fit.log_rate_sd**2
    eigen = np.linalg.eigvalsh(hessian)
    if eigen[0] <= max(1e-10, eigen[-1] * 1e-10):
        fit.interval_status = "unavailable-unidentified-or-nonpositive-curvature"
        return
    return np.linalg.inv(hessian)[: len(problem.free), : len(problem.free)]


def laplace_intervals(fit, problem, level=0.95):
    """Unadjusted Gaussian curvature intervals, conditional on the fitted model."""
    covariance = curvature_covariance(fit, problem, level)
    if covariance is None:
        return
    c = problem.chronology
    lower, upper = fit.ages.copy(), fit.ages.copy()
    widths = norm.ppf((1 + level) / 2) * np.sqrt(np.diag(covariance))
    lower[problem.free] -= widths
    upper[problem.free] += widths
    if np.any(lower < c.lower) or np.any(upper > c.upper):
        fit.interval_status = "unavailable-gaussian-interval-crosses-bound"
        return
    fit.interval_lower, fit.interval_upper = lower, upper
    fit.interval_status = (
        "conditional-laplace" if problem.likelihood else "conditional-profile-curvature"
    )
    if problem.rate_variance_estimated:
        fit.diagnostics.append("unadjusted_curvature_with_estimated_rate_variance")
