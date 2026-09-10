"""Profiled native OU fitting for a declared shared shift layout.

Finite alpha is bounded and searched from several deterministic starts. Exact
zero/infinity alpha limits and zero variance components are separate candidates.
The returned maximum is numerical, not a certified continuous global optimum.
"""

import math
from dataclasses import dataclass, replace

import numpy as np
from scipy.optimize import minimize

from nwkit.optimization import deterministic_multistart, global_bounded_scalar_minimize
from nwkit.shift_native_identifiability import covariance_identifiability
from nwkit.shift_native_model import evaluate_trait, restore_trait_fit


@dataclass(frozen=True)
class NativeFitOptions:
    root_model: str = "OUfixedRoot"
    estimate_measurement_error: bool = False
    alpha_height_bounds: tuple[float, float] = (1e-6, 1000.0)
    variance_bounds: tuple[float, float] = (1e-10, 1e4)
    optimizer_starts: int = 3
    maxiter: int = 300

    def validate(self):
        if self.root_model not in {"OUfixedRoot", "OUrandomRoot"}:
            raise ValueError("Unknown native root model.")
        for bounds in (self.alpha_height_bounds, self.variance_bounds):
            if (
                len(bounds) != 2
                or not all(math.isfinite(x) and x > 0 for x in bounds)
                or bounds[0] >= bounds[1]
            ):
                raise ValueError(
                    "Native optimizer bounds must be finite, positive and increasing."
                )
        if self.optimizer_starts < 1 or self.maxiter < 1:
            raise ValueError("Native optimizer starts and iterations must be positive.")


def _checked_minimize(function, initial, *, method, bounds, options):
    """Check a unit-box optimum using central, one-sided-at-boundary scores."""

    def gradient(x):
        result = []
        for j in range(len(x)):
            left, right = x.copy(), x.copy()
            left[j], right[j] = max(0, x[j] - 1e-5), min(1, x[j] + 1e-5)
            result.append((function(right) - function(left)) / (right[j] - left[j]))
        return np.asarray(result)

    result = minimize(function, initial, method=method, bounds=bounds, options=options)
    if not result.success or not np.isfinite(result.fun):
        return result
    checked = gradient(result.x)
    checked[(result.x <= 1e-7) & (checked > 0)] = 0
    checked[(result.x >= 1 - 1e-7) & (checked < 0)] = 0
    norm = float(np.max(np.abs(checked)))
    if norm > 2e-4:
        result = minimize(
            function,
            result.x,
            jac=gradient,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": options["maxiter"], "ftol": 1e-14, "gtol": 1e-7},
        )
        checked = gradient(result.x)
        checked[(result.x <= 1e-7) & (checked > 0)] = 0
        checked[(result.x >= 1 - 1e-7) & (checked < 0)] = 0
        norm = float(np.max(np.abs(checked)))
    result.success = bool(result.success and math.isfinite(norm) and norm <= 2e-4)
    result.message = f"{result.message}; checked projected gradient={norm:.6g}"
    return result


class _TraitProblem:
    def __init__(
        self, data, layout, trait, options, fixed_process=None, fixed_noise=None
    ):
        self.data, self.layout, self.trait, self.options = data, layout, trait, options
        self.node_groups = layout.node_groups(data.tree)
        self.fixed_process = fixed_process
        self.fixed_noise = fixed_noise
        self.known_error = bool(np.any(data.variances[:, trait] > 0))
        self.num_observations = int(np.sum(np.isfinite(data.values[:, trait])))
        self.evaluations = 0
        self.failed_evaluations = 0

    def evaluate(self, alpha, process, noise):
        self.evaluations += 1
        return evaluate_trait(
            self.data,
            self.layout,
            self.trait,
            alpha,
            process,
            noise,
            root_model=self.options.root_model,
            _node_groups=self.node_groups,
        )

    def analytic_scale(self, alpha):
        unit = self.evaluate(alpha, 1.0, 0.0)
        threshold = (
            100
            * np.finfo(float).eps ** 2
            * max(1.0, float(np.nansum(self.data.values[:, self.trait] ** 2)))
        )
        if unit.quadratic <= threshold:
            raise ValueError(
                "Zero residual variance: the ordinary Gaussian maximum is undefined."
            )
        variance = unit.quadratic / unit.num_observations
        likelihood = unit.log_likelihood + 0.5 * (
            unit.quadratic - unit.num_observations * (1 + math.log(variance))
        )
        return replace(
            unit,
            process_variance=variance,
            coefficient_covariance=unit.coefficient_covariance * variance,
            log_likelihood=likelihood,
            quadratic=float(unit.num_observations),
        )

    def scalar_variance(self, alpha):
        lower, upper = np.log(self.options.variance_bounds)
        candidates = []

        def objective(coordinate):
            try:
                fit = self.evaluate(
                    alpha, math.exp(coordinate), self.fixed_noise or 0.0
                )
                return -fit.log_likelihood / self.num_observations
            except (ValueError, ArithmeticError, np.linalg.LinAlgError):
                self.failed_evaluations += 1
                return math.inf

        optimized = global_bounded_scalar_minimize(
            objective, (lower, upper), grid_size=9
        )
        if not optimized.success or not math.isfinite(optimized.fun):
            raise ValueError("Positive process-variance optimization failed.")
        if optimized.success and math.isfinite(optimized.fun):
            candidates.append(
                self.evaluate(alpha, math.exp(optimized.x), self.fixed_noise or 0.0)
            )
        try:
            candidates.append(self.evaluate(alpha, 0.0, self.fixed_noise or 0.0))
        except (ValueError, ArithmeticError, np.linalg.LinAlgError):
            self.failed_evaluations += 1
        if not candidates:
            raise ValueError("No finite process-variance fit was found.")
        return max(
            candidates, key=lambda fit: (fit.log_likelihood, -fit.process_variance)
        )

    def at_alpha(self, alpha):
        if self.fixed_process is not None:
            return self.evaluate(alpha, self.fixed_process, self.fixed_noise or 0.0)
        if not self.known_error and not self.fixed_noise:
            return self.analytic_scale(alpha)
        return self.scalar_variance(alpha)


def _variance_decode(problem, coordinates, *, alpha):
    options = problem.options
    lower, upper = np.log(options.variance_bounds)
    total = math.exp(lower + coordinates[0] * (upper - lower))
    if problem.fixed_process is not None:
        return problem.fixed_process, total
    if math.isinf(alpha):
        # Only the sum of independent process and unknown observation variances
        # is estimable here. This canonical representation is labelled on export.
        return total, 0.0
    fraction = coordinates[1]
    return total * fraction, total * (1 - fraction)


def _joint_fit(problem, fixed_alpha):
    free_alpha = fixed_alpha is None
    a_lower, a_upper = np.log(problem.options.alpha_height_bounds)
    size = int(free_alpha) + (
        1
        if problem.fixed_process is not None
        or (fixed_alpha is not None and math.isinf(fixed_alpha))
        else 2
    )

    def evaluate(coordinates):
        alpha = (
            math.exp(a_lower + coordinates[0] * (a_upper - a_lower))
            if free_alpha
            else fixed_alpha
        )
        process, noise = _variance_decode(
            problem, coordinates[int(free_alpha) :], alpha=alpha
        )
        return problem.evaluate(alpha, process, noise)

    def objective(coordinates):
        try:
            return -evaluate(coordinates).log_likelihood / problem.num_observations
        except (ValueError, ArithmeticError, np.linalg.LinAlgError):
            problem.failed_evaluations += 1
            return 1e100

    fractions = tuple(
        np.linspace(0.2, 0.8, max(0, problem.options.optimizer_starts - 1))
    )
    optimized = deterministic_multistart(
        objective,
        np.full(size, 0.5),
        [(0.0, 1.0)] * size,
        maxiter=problem.options.maxiter,
        fractions=fractions,
        patterned_starts=False,
        fallback=False,
        minimizer=_checked_minimize,
    )
    fit = evaluate(optimized.x)
    return fit, {
        "success": optimized.success,
        "starts": optimized.starts,
        "converged_starts": optimized.converged_starts,
        "failed_starts": optimized.failed_starts,
        "message": optimized.message,
    }


def _profile_alpha(problem, fixed_alpha):
    if fixed_alpha is not None:
        return problem.at_alpha(fixed_alpha), {
            "success": True,
            "message": "alpha fixed",
        }
    lower, upper = np.log(problem.options.alpha_height_bounds)

    def objective(coordinate):
        try:
            return (
                -problem.at_alpha(math.exp(coordinate)).log_likelihood
                / problem.num_observations
            )
        except (ValueError, ArithmeticError, np.linalg.LinAlgError):
            problem.failed_evaluations += 1
            return math.inf

    optimized = global_bounded_scalar_minimize(objective, (lower, upper), grid_size=17)
    if not optimized.success or not math.isfinite(optimized.fun):
        raise ValueError("No finite interior-alpha profile fit was found.")
    return problem.at_alpha(math.exp(optimized.x)), {
        "success": True,
        "message": optimized.message,
    }


def _require_observed_mean_rank(problem, alpha):
    if alpha is None:
        return
    design = problem.layout.design(problem.data.tree, alpha)[
        np.isfinite(problem.data.values[:, problem.trait])
    ]
    norms = np.linalg.norm(design, axis=0)
    if np.any(norms == 0) or np.linalg.matrix_rank(design / norms) < design.shape[1]:
        raise ValueError("The observed mean design is rank deficient.")


def fit_native_trait(
    data,
    layout,
    trait,
    *,
    options=None,
    alpha_height=None,
    process_variance=None,
    measurement_variance=None,
):
    """Fit one trait in normalized units, retaining failures and limit fits."""
    options = NativeFitOptions() if options is None else options
    options.validate()
    if options.estimate_measurement_error and measurement_variance is not None:
        raise ValueError(
            "Additional measurement variance cannot be fixed and estimated together."
        )
    noise = measurement_variance if measurement_variance is not None else 0.0
    problem = _TraitProblem(data, layout, trait, options, process_variance, noise)
    alphas = [alpha_height]
    if alpha_height is None:
        alphas += ([0.0] if options.root_model == "OUfixedRoot" else []) + [math.inf]
    records, candidates = [], []
    for alpha in alphas:
        try:
            _require_observed_mean_rank(problem, alpha)
            fit, optimizer = (
                _joint_fit(problem, alpha)
                if options.estimate_measurement_error
                else _profile_alpha(problem, alpha)
            )
        except (ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
            records.append(
                {
                    "alpha_mode": "finite" if alpha is None else str(alpha),
                    "success": False,
                    "message": str(exc),
                }
            )
            continue
        records.append(
            {
                "alpha_mode": "finite" if alpha is None else str(alpha),
                "success": True,
                "log_likelihood": fit.log_likelihood,
                "optimizer": optimizer,
            }
        )
        candidates.append(fit)
    observed = np.isfinite(data.values[:, trait])
    zero_total_eligible = bool(np.all(data.variances[observed, trait] > 0))
    if options.estimate_measurement_error and (
        process_variance is not None or zero_total_eligible
    ):
        boundary_problem = (
            problem
            if process_variance is not None
            else _TraitProblem(data, layout, trait, options, 0.0, 0.0)
        )
        for alpha in alphas:
            try:
                zero_noise, optimizer = _profile_alpha(boundary_problem, alpha)
                candidates.append(zero_noise)
                records.append(
                    {
                        "alpha_mode": "finite" if alpha is None else str(alpha),
                        "measurement_variance": 0.0,
                        "process_variance": boundary_problem.fixed_process,
                        "success": True,
                        "log_likelihood": zero_noise.log_likelihood,
                        "optimizer": optimizer,
                    }
                )
            except (ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
                records.append(
                    {
                        "alpha_mode": "finite" if alpha is None else str(alpha),
                        "measurement_variance": 0.0,
                        "process_variance": boundary_problem.fixed_process,
                        "success": False,
                        "message": str(exc),
                    }
                )
        if boundary_problem is not problem:
            problem.evaluations += boundary_problem.evaluations
            problem.failed_evaluations += boundary_problem.failed_evaluations
    if not candidates:
        raise ValueError(
            "Native fixed-layout fitting failed: "
            + "; ".join(str(record["message"]) for record in records)
        )
    best_score = max(f.log_likelihood for f in candidates)
    tied = [f for f in candidates if f.log_likelihood >= best_score - 1e-9]
    fit = min(
        tied,
        key=lambda f: (
            f.alpha_height not in (0.0, math.inf),
            f.alpha_height,
            f.process_variance + f.measurement_variance,
        ),
    )
    bounded_variance = None
    if options.estimate_measurement_error:
        bounded_variance = (
            fit.measurement_variance
            if process_variance is not None
            else fit.process_variance + fit.measurement_variance
        )
    elif process_variance is None and (problem.known_error or noise):
        bounded_variance = fit.process_variance
    lower, upper = options.variance_bounds
    variance_boundary = (
        bounded_variance is not None
        and bounded_variance > 0
        and (
            bounded_variance <= lower * (1 + 1e-5)
            or bounded_variance >= upper * (1 - 1e-5)
        )
    )
    variance_mode_ambiguity = bool(
        options.estimate_measurement_error
        and process_variance is None
        and alpha_height is None
        and (
            np.ptp([candidate.process_variance for candidate in tied])
            > 1e-5 * max(1.0, fit.process_variance + fit.measurement_variance)
        )
    )
    # Partial covariance-search failure is retained, not described as full coverage.
    diagnostics = {
        "evaluations": problem.evaluations,
        "invalid_evaluations": problem.failed_evaluations,
        "alpha_candidates": records,
        "continuous_globally_optimal": False,
        "complete_alpha_modes": all(r["success"] for r in records),
        "alpha_height_bounds": list(options.alpha_height_bounds),
        "normalized_variance_bounds": list(options.variance_bounds),
        "alpha_estimated": alpha_height is None,
        "process_variance_estimated": process_variance is None,
        "measurement_variance_estimated": options.estimate_measurement_error,
        "nuisance_variance_at_numerical_bound": variance_boundary,
        "equivalent_modes_disagree_on_variance_decomposition": variance_mode_ambiguity,
    }
    return fit, diagnostics


def fit_native_layout(
    data,
    layout,
    *,
    options=None,
    alpha_height=None,
    process_variance=None,
    measurement_variance=None,
):
    """Fit independent trait covariances and shared shift/regime identities.

    Fixed variances are supplied in original trait units (scalar or per-trait).
    Alpha-height coordinates are dimensionless, scalar or per-trait.
    """
    options = NativeFitOptions() if options is None else options
    dimension = len(data.trait_names)

    def parameter(value, trait, variance=False):
        if value is None:
            return None
        array = np.broadcast_to(np.asarray(value, dtype=float), (dimension,))
        result = float(array[trait])
        normalized = (
            result / data.scales[trait] / data.scales[trait] if variance else result
        )
        if variance and (
            not math.isfinite(normalized) or (result > 0 and normalized == 0)
        ):
            raise ValueError(
                "Fixed variance is not representable in normalized trait units."
            )
        return normalized

    fits, records = [], []
    for trait in range(dimension):
        fit, diagnostics = fit_native_trait(
            data,
            layout,
            trait,
            options=options,
            alpha_height=parameter(alpha_height, trait),
            process_variance=parameter(process_variance, trait, True),
            measurement_variance=parameter(measurement_variance, trait, True),
        )
        restored = restore_trait_fit(data, trait, fit)
        restored["optimizer"] = diagnostics
        identification = covariance_identifiability(
            data,
            trait,
            fit,
            alpha_free=alpha_height is None,
            process_free=process_variance is None,
            noise_free=options.estimate_measurement_error,
        )
        identification["between_mode_variance_ambiguity"] = diagnostics[
            "equivalent_modes_disagree_on_variance_decomposition"
        ]
        if identification["between_mode_variance_ambiguity"]:
            identification["variance_decomposition_supported"] = False
        restored["identifiability_diagnostic"] = identification
        if alpha_height is None and (
            fit.process_variance == 0
            or (
                0 < fit.alpha_height < math.inf
                and not identification["finite_alpha_supported"]
            )
        ):
            restored["alpha_candidate"] = restored["alpha"]
            restored["alpha_height_candidate"] = restored["alpha_height"]
            restored["alpha_candidate_status"] = restored["alpha_status"]
            restored["alpha"] = None
            restored["alpha_height"] = None
            restored["alpha_status"] = (
                "unsupported_by_covariance_identifiability_diagnostic"
            )
        restored["covariance_components_identifiable"] = identification[
            "variance_decomposition_supported"
        ]
        if not restored["covariance_components_identifiable"]:
            total = restored["process_tip_variance"] + restored["measurement_variance"]
            restored["fitted_total_tip_variance"] = total
            if math.isinf(fit.alpha_height):
                restored["unstructured_tip_variance"] = total
            restored["process_tip_variance"] = None
            restored["measurement_variance"] = None
            restored["sigma2"] = None
        fits.append(fit)
        records.append(restored)
    return {
        "layout": layout,
        "fits": fits,
        "traits": records,
        "log_likelihood": sum(r["log_likelihood"] for r in records),
        "root_model": options.root_model,
        "trait_covariance": "diagonal",
    }
