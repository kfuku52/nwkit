"""Early-burst and directional-drift Brownian ancestral reconstruction."""

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize_scalar

from nwkit.continuous_asr import GaussianMarginal, _finite_number, compute_bm_marginals
from nwkit.evolution import build_evolutionary_process


@dataclass(frozen=True)
class EarlyBurstFit:
    sigma2: float
    sigma2_estimated: bool
    eb_rate: float
    eb_rate_estimated: bool
    eb_rate_bounds: tuple[float, float]
    restricted_log_likelihood: float | None
    num_observed: int
    num_effective_observations: int
    residual_df: int
    fit_status: str
    optimizer_success: bool
    optimizer_message: str
    optimizer_grid_evaluations: int


@dataclass(frozen=True)
class BrownianDriftFit:
    sigma2: float
    sigma2_estimated: bool
    drift: float
    drift_estimated: bool
    restricted_log_likelihood: float | None
    num_observed: int
    num_effective_observations: int
    residual_df: int
    fit_status: str
    optimizer_success: bool
    optimizer_message: str
    optimizer_grid_evaluations: int


def _node_depths(tree):
    depths = {tree: 0.0}
    for node in tree.traverse(strategy="preorder"):
        if not node.is_root:
            depths[node] = depths[node.up] + float(node.dist)
    return depths


def default_eb_rate_bounds(tree):
    depths = _node_depths(tree)
    maximum = max(depths.values(), default=0.0)
    scale = maximum if maximum > 0.0 else 1.0
    return -10.0 / scale, 10.0 / scale


def parse_eb_rate_bounds(value, tree):
    if value in (None, ""):
        return default_eb_rate_bounds(tree)
    items = [item.strip() for item in str(value).split(",")]
    if len(items) != 2:
        raise ValueError("'--eb-rate-bounds' must contain exactly two values.")
    try:
        lower, upper = (float(item) for item in items)
    except ValueError as exc:
        raise ValueError("'--eb-rate-bounds' values must be numeric.") from exc
    if not all(math.isfinite(item) for item in (lower, upper)):
        raise ValueError("'--eb-rate-bounds' values must be finite.")
    if lower >= upper:
        raise ValueError(
            "'--eb-rate-bounds' lower bound must be smaller than the upper bound."
        )
    return lower, upper


def _run_eb(tree, values_by_leaf, rate, sigma2, standard_errors, tree_validated):
    # The ASR command historically calls the two-sided exponential rate-change
    # model EB.  Internally it is the same process that regression calls ACDC,
    # including positive acceleration parameters.
    process = build_evolutionary_process(
        tree,
        model="acdc",
        parameter=rate,
        root_mode="flat",
        allow_zero=True,
    )
    posterior, fit = compute_bm_marginals(
        tree,
        values_by_leaf,
        sigma2=sigma2,
        standard_errors=standard_errors,
        _tree_validated=tree_validated,
        _process=process,
    )
    return posterior, fit


def compute_eb_marginals(
    tree,
    values_by_leaf,
    *,
    sigma2=None,
    eb_rate=None,
    eb_rate_bounds=None,
    standard_errors=None,
    _tree_validated=False,
):
    """Return flat-root marginals under exponential time-varying BM rates."""

    fixed_sigma2 = (
        None if sigma2 is None else _finite_number(sigma2, "--sigma2", nonnegative=True)
    )
    bounds = default_eb_rate_bounds(tree) if eb_rate_bounds is None else eb_rate_bounds
    bounds = float(bounds[0]), float(bounds[1])
    if not all(math.isfinite(value) for value in bounds) or bounds[0] >= bounds[1]:
        raise ValueError("EB rate bounds must be two increasing finite values.")
    if eb_rate is not None:
        fitted_rate = _finite_number(eb_rate, "--eb-rate")
        posterior, bm_fit = _run_eb(
            tree,
            values_by_leaf,
            fitted_rate,
            fixed_sigma2,
            standard_errors,
            _tree_validated,
        )
        evaluations = 1
        success = True
        message = "EB rate fixed"
    else:
        cache = {}

        def evaluate(rate):
            key = float(rate)
            if key not in cache:
                try:
                    cache[key] = _run_eb(
                        tree,
                        values_by_leaf,
                        key,
                        fixed_sigma2,
                        standard_errors,
                        _tree_validated,
                    )
                except (ValueError, ArithmeticError):
                    cache[key] = None
            result = cache[key]
            if result is None or result[1].restricted_log_likelihood is None:
                return 1e100
            return -result[1].restricted_log_likelihood

        grid = np.linspace(bounds[0], bounds[1], 41)
        objectives = np.asarray([evaluate(value) for value in grid])
        candidates = [
            (float(objectives[index]), float(grid[index]), "grid")
            for index in range(len(grid))
            if math.isfinite(float(objectives[index]))
            and float(objectives[index]) < 1e99
        ]
        for index in range(1, len(grid) - 1):
            if (
                objectives[index] <= objectives[index - 1]
                and objectives[index] <= objectives[index + 1]
            ):
                result = minimize_scalar(
                    evaluate,
                    bounds=(float(grid[index - 1]), float(grid[index + 1])),
                    method="bounded",
                    options={"xatol": 1e-10},
                )
                if result.success and math.isfinite(float(result.fun)):
                    candidates.append(
                        (float(result.fun), float(result.x), str(result.message))
                    )
        zero = 0.0
        if bounds[0] <= zero <= bounds[1]:
            candidates.append((evaluate(zero), zero, "BM boundary candidate"))
        if not candidates:
            raise ValueError("Failed to find a finite EB likelihood.")
        _, fitted_rate, message = min(candidates, key=lambda item: item[0])
        posterior, bm_fit = _run_eb(
            tree,
            values_by_leaf,
            fitted_rate,
            fixed_sigma2,
            standard_errors,
            _tree_validated,
        )
        evaluations = len(cache)
        success = True
        message = f"profile grid and bounded polishing; {message}"
        # A flat profile means EB rate and sigma2 are confounded on this tree.
        finite_objectives = objectives[np.isfinite(objectives) & (objectives < 1e99)]
        if len(finite_objectives) > 1 and np.ptp(finite_objectives) <= 1e-10 * max(
            1.0, abs(float(np.min(finite_objectives)))
        ):
            raise ValueError(
                "EB rate is not identifiable on this tree and observation pattern; fix --eb-rate."
            )
    statuses = []
    tolerance = 1e-6 * max(1.0, abs(bounds[0]), abs(bounds[1]))
    if eb_rate is None and abs(fitted_rate - bounds[0]) <= tolerance:
        statuses.append("eb_rate_lower_boundary")
    if eb_rate is None and abs(fitted_rate - bounds[1]) <= tolerance:
        statuses.append("eb_rate_upper_boundary")
    if bm_fit.fit_status != "ok":
        statuses.append(bm_fit.fit_status)
    fit = EarlyBurstFit(
        sigma2=bm_fit.sigma2,
        sigma2_estimated=bm_fit.sigma2_estimated,
        eb_rate=fitted_rate,
        eb_rate_estimated=eb_rate is None,
        eb_rate_bounds=bounds,
        restricted_log_likelihood=bm_fit.restricted_log_likelihood,
        num_observed=bm_fit.num_observed,
        num_effective_observations=bm_fit.num_effective_observations,
        residual_df=bm_fit.residual_df,
        fit_status="+".join(statuses) if statuses else "ok",
        optimizer_success=success,
        optimizer_message=message,
        optimizer_grid_evaluations=evaluations,
    )
    return posterior, fit


def _run_drift(
    tree,
    values_by_leaf,
    depths,
    leaf_by_name,
    drift,
    sigma2,
    standard_errors,
    tree_validated,
):
    residual_values = {
        name: None
        if value is None
        else float(value) - drift * depths[leaf_by_name[name]]
        for name, value in values_by_leaf.items()
    }
    residual_posterior, fit = compute_bm_marginals(
        tree,
        residual_values,
        sigma2=sigma2,
        standard_errors=standard_errors,
        _tree_validated=tree_validated,
    )
    posterior = {
        node: GaussianMarginal(marginal.mean + drift * depths[node], marginal.variance)
        for node, marginal in residual_posterior.items()
    }
    return posterior, fit


def _minimize_unbounded_profile(evaluate, center, scale, *, quadratic=False):
    """Bracket a finite profile minimum without imposing parameter bounds."""

    sampled = {float(center): float(evaluate(center))}
    exponent = 0
    while True:
        try:
            radius = math.ldexp(float(scale), exponent)
        except OverflowError:
            break
        if not math.isfinite(radius):
            break
        left = center - radius
        right = center + radius
        if not math.isfinite(left) or not math.isfinite(right):
            break
        sampled.setdefault(float(left), float(evaluate(left)))
        sampled.setdefault(float(right), float(evaluate(right)))
        finite = [
            (objective, point)
            for point, objective in sampled.items()
            if math.isfinite(objective) and objective < 1e99
        ]
        if finite:
            best_objective, best_point = min(finite)
            if (
                left < best_point < right
                and sampled[left] > best_objective
                and sampled[right] > best_objective
            ):
                candidates = [(best_objective, best_point, "adaptive grid")]
                if quadratic:
                    center_objective = sampled[float(center)]
                    curvature = sampled[left] - 2.0 * center_objective + sampled[right]
                    if math.isfinite(curvature) and curvature > 0.0:
                        vertex = center + radius * (sampled[left] - sampled[right]) / (
                            2.0 * curvature
                        )
                        if left < vertex < right:
                            vertex_objective = float(evaluate(vertex))
                            if (
                                math.isfinite(vertex_objective)
                                and vertex_objective < 1e99
                            ):
                                return (
                                    float(vertex),
                                    vertex_objective,
                                    "quadratic profile interpolation",
                                )
                result = minimize_scalar(
                    evaluate,
                    bounds=(left, right),
                    method="bounded",
                    options={"xatol": 1e-10},
                )
                if (
                    result.success
                    and math.isfinite(float(result.fun))
                    and float(result.fun) < 1e99
                ):
                    candidates.append(
                        (float(result.fun), float(result.x), str(result.message))
                    )
                objective, point, message = min(candidates)
                return point, objective, message
        exponent += 1
    raise ValueError(
        "Failed to bracket a finite BM-DRIFT trend; rescale trait or time units."
    )


def compute_bm_drift_marginals(
    tree,
    values_by_leaf,
    *,
    sigma2=None,
    drift=None,
    standard_errors=None,
    _tree_validated=False,
):
    """Return Brownian marginals with a fixed or fitted directional trend."""

    fixed_sigma2 = (
        None if sigma2 is None else _finite_number(sigma2, "--sigma2", nonnegative=True)
    )
    depths = _node_depths(tree)
    leaf_by_name = {str(leaf.name): leaf for leaf in tree.leaves()}
    observed = [
        (leaf_by_name[name], _finite_number(value, f"Trait value for '{name}'"))
        for name, value in values_by_leaf.items()
        if value is not None
    ]
    if drift is not None:
        fitted_drift = _finite_number(drift, "--drift")
        posterior, bm_fit = _run_drift(
            tree,
            values_by_leaf,
            depths,
            leaf_by_name,
            fitted_drift,
            fixed_sigma2,
            standard_errors,
            _tree_validated,
        )
        evaluations = 1
        success = True
        message = "drift fixed"
    else:
        tip_depths = np.asarray([depths[node] for node, _ in observed], dtype=float)
        if len(tip_depths) < 2 or np.ptp(tip_depths) == 0.0:
            raise ValueError(
                "BM-DRIFT drift is not identifiable when observed tips have one sampling depth; fix --drift."
            )
        responses = np.asarray([value for _, value in observed], dtype=float)
        ordinary_slope = float(
            np.dot(tip_depths - tip_depths.mean(), responses - responses.mean())
            / np.dot(tip_depths - tip_depths.mean(), tip_depths - tip_depths.mean())
        )
        if not math.isfinite(ordinary_slope):
            raise ValueError(
                "BM-DRIFT trend exceeds floating-point range; rescale trait or time units."
            )
        scale = max(
            abs(ordinary_slope),
            float(np.ptp(responses)) / float(np.ptp(tip_depths)),
            np.finfo(float).eps,
        )
        try:
            boundary_candidate = _run_drift(
                tree,
                values_by_leaf,
                depths,
                leaf_by_name,
                ordinary_slope,
                fixed_sigma2,
                standard_errors,
                _tree_validated,
            )
        except (ValueError, ArithmeticError):
            boundary_candidate = None
        if (
            boundary_candidate is not None
            and boundary_candidate[1].sigma2 == 0.0
            and boundary_candidate[1].restricted_log_likelihood is None
        ):
            fitted_drift = ordinary_slope
            posterior, bm_fit = boundary_candidate
            evaluations = 1
            success = True
            message = "exact zero-diffusion boundary"
        else:
            cache = {}
            if boundary_candidate is not None:
                cache[ordinary_slope] = boundary_candidate

            def evaluate(current_drift):
                key = float(current_drift)
                if key not in cache:
                    try:
                        cache[key] = _run_drift(
                            tree,
                            values_by_leaf,
                            depths,
                            leaf_by_name,
                            key,
                            fixed_sigma2,
                            standard_errors,
                            _tree_validated,
                        )
                    except (ValueError, ArithmeticError):
                        cache[key] = None
                result = cache[key]
                if result is None or result[1].restricted_log_likelihood is None:
                    return 1e100
                return -result[1].restricted_log_likelihood

            fitted_drift, objective, optimizer_message = _minimize_unbounded_profile(
                evaluate,
                ordinary_slope,
                scale,
                quadratic=fixed_sigma2 is not None,
            )
            if not math.isfinite(objective) or objective >= 1e99:
                raise ValueError("Failed to estimate a finite BM-DRIFT trend.")
            posterior, bm_fit = _run_drift(
                tree,
                values_by_leaf,
                depths,
                leaf_by_name,
                fitted_drift,
                fixed_sigma2,
                standard_errors,
                _tree_validated,
            )
            evaluations = len(cache)
            success = True
            message = f"adaptive unbounded profile; {optimizer_message}"
    fit = BrownianDriftFit(
        sigma2=bm_fit.sigma2,
        sigma2_estimated=bm_fit.sigma2_estimated,
        drift=fitted_drift,
        drift_estimated=drift is None,
        restricted_log_likelihood=bm_fit.restricted_log_likelihood,
        num_observed=bm_fit.num_observed,
        num_effective_observations=bm_fit.num_effective_observations,
        residual_df=bm_fit.residual_df,
        fit_status=bm_fit.fit_status,
        optimizer_success=success,
        optimizer_message=message,
        optimizer_grid_evaluations=evaluations,
    )
    return posterior, fit
