"""Linear-time single-optimum Ornstein-Uhlenbeck ancestral reconstruction."""

import math
from dataclasses import dataclass, replace
from typing import Any

import numpy as np
from scipy.optimize import minimize, minimize_scalar

from nwkit.continuous_asr import (
    GaussianMarginal,
    _finite_number,
    _ldexp,
    _observations_by_group,
    _tree_groups,
)

_LOG_2PI = math.log(2.0 * math.pi)
_LOG_2 = math.log(2.0)


@dataclass(frozen=True)
class OrnsteinUhlenbeckFit:
    alpha: float
    alpha_estimated: bool
    theta: float
    theta_estimated: bool
    sigma2: float
    sigma2_estimated: bool
    root_variance: float
    log_likelihood: float
    num_observed: int
    num_effective_observations: int
    num_observed_positions: int
    fit_status: str
    optimizer_success: bool
    optimizer_message: str
    optimizer_starts: int
    optimizer_converged_starts: int
    optimizer_failed_starts: int
    optimizer_grid_evaluations: int
    alpha_bounds: tuple[float, float]
    root_variance_bounds: tuple[float, float] | None


@dataclass(frozen=True, slots=True)
class _Factor:
    marginal: GaussianMarginal | None
    log_weight: float = 0.0


@dataclass(frozen=True, slots=True)
class _AffineFactor:
    """Gaussian factor whose mean is affine in ``theta - reference``."""

    mean_offset: float | None
    mean_slope: float = 0.0
    variance: float = 0.0
    nll_quadratic: float = 0.0
    nll_linear: float = 0.0
    nll_constant: float = 0.0


@dataclass(frozen=True)
class _OuData:
    node_groups: dict[Any, int]
    parents: tuple[int, ...]
    lengths: tuple[float, ...]
    children: tuple[tuple[int, ...], ...]
    local: tuple[_Factor, ...]
    num_observed: int
    num_effective_observations: int
    num_observed_positions: int
    observed_values: tuple[float, ...]
    observed_errors: tuple[float, ...]
    exact_values: dict[int, float]
    likelihood_dimensions: int
    trait_center: float
    trait_exponent: int


def _log_normal_density(value, mean, variance):
    if variance <= 0.0 or not math.isfinite(variance):
        raise ValueError("OU Gaussian variance must be positive and finite.")
    standardized = (value - mean) / math.sqrt(variance)
    if not math.isfinite(standardized):
        return -math.inf
    return -0.5 * (_LOG_2PI + math.log(variance) + standardized * standardized)


def _combine(first, second):
    """Multiply two normalized scalar Gaussian factors and retain the constant."""

    if first.marginal is None:
        return _Factor(second.marginal, first.log_weight + second.log_weight)
    if second.marginal is None:
        return _Factor(first.marginal, first.log_weight + second.log_weight)
    left, right = first.marginal, second.marginal
    variance_sum = left.variance + right.variance
    if not math.isfinite(variance_sum):
        raise ValueError(
            "An OU conditional variance exceeds floating-point range; rescale units."
        )
    if variance_sum == 0.0:
        if left.mean != right.mean:
            raise ValueError(
                "Conflicting exact observations have zero likelihood under the OU model."
            )
        return _Factor(left, first.log_weight + second.log_weight)
    log_overlap = _log_normal_density(left.mean, right.mean, variance_sum)
    if left.variance == 0.0:
        marginal = left
    elif right.variance == 0.0:
        marginal = right
    else:
        larger = max(left.variance, right.variance)
        ratio = min(left.variance, right.variance) / larger
        denominator = 1.0 + ratio
        left_weight = (right.variance / larger) / denominator
        right_weight = (left.variance / larger) / denominator
        variance = min(left.variance, right.variance) / denominator
        if variance <= 0.0:
            raise ValueError(
                "A positive OU conditional variance underflowed; rescale units."
            )
        marginal = GaussianMarginal(
            left_weight * left.mean + right_weight * right.mean, variance
        )
    return _Factor(marginal, first.log_weight + second.log_weight + log_overlap)


def _normal_nll_polynomial(mean_offset, mean_slope, variance):
    if variance <= 0.0 or not math.isfinite(variance):
        raise ValueError("OU Gaussian variance must be positive and finite.")
    scale = math.sqrt(variance)
    standardized_offset = mean_offset / scale
    standardized_slope = mean_slope / scale
    if not all(
        math.isfinite(value) for value in (standardized_offset, standardized_slope)
    ):
        raise ValueError("An OU profiled likelihood exceeds floating-point range.")
    return (
        0.5 * standardized_slope * standardized_slope,
        standardized_offset * standardized_slope,
        0.5
        * (_LOG_2PI + math.log(variance) + standardized_offset * standardized_offset),
    )


def _combine_affine(first, second):
    quadratic = first.nll_quadratic + second.nll_quadratic
    linear = first.nll_linear + second.nll_linear
    constant = first.nll_constant + second.nll_constant
    if first.mean_offset is None:
        return replace(
            second,
            nll_quadratic=quadratic,
            nll_linear=linear,
            nll_constant=constant,
        )
    if second.mean_offset is None:
        return replace(
            first,
            nll_quadratic=quadratic,
            nll_linear=linear,
            nll_constant=constant,
        )
    variance_sum = first.variance + second.variance
    if not math.isfinite(variance_sum):
        raise ValueError(
            "An OU conditional variance exceeds floating-point range; rescale units."
        )
    if variance_sum == 0.0:
        if (
            first.mean_offset != second.mean_offset
            or first.mean_slope != second.mean_slope
        ):
            raise ValueError(
                "Conflicting exact observations have zero likelihood under the OU model."
            )
        return replace(
            first,
            nll_quadratic=quadratic,
            nll_linear=linear,
            nll_constant=constant,
        )
    extra = _normal_nll_polynomial(
        first.mean_offset - second.mean_offset,
        first.mean_slope - second.mean_slope,
        variance_sum,
    )
    quadratic += extra[0]
    linear += extra[1]
    constant += extra[2]
    if first.variance == 0.0:
        mean_offset = first.mean_offset
        mean_slope = first.mean_slope
        variance = 0.0
    elif second.variance == 0.0:
        mean_offset = second.mean_offset
        mean_slope = second.mean_slope
        variance = 0.0
    else:
        larger = max(first.variance, second.variance)
        ratio = min(first.variance, second.variance) / larger
        denominator = 1.0 + ratio
        first_weight = (second.variance / larger) / denominator
        second_weight = (first.variance / larger) / denominator
        variance = min(first.variance, second.variance) / denominator
        if variance <= 0.0:
            raise ValueError(
                "A positive OU conditional variance underflowed; rescale units."
            )
        mean_offset = (
            first_weight * first.mean_offset + second_weight * second.mean_offset
        )
        mean_slope = first_weight * first.mean_slope + second_weight * second.mean_slope
    return _AffineFactor(
        mean_offset,
        mean_slope,
        variance,
        quadratic,
        linear,
        constant,
    )


def _prepare_data(
    tree,
    values_by_leaf,
    standard_errors,
    *,
    process_sd=0.0,
    tree_validated=False,
):
    node_groups, parents, lengths = _tree_groups(tree, validated=tree_validated)
    grouped, num_observed, exact_values = _observations_by_group(
        tree, node_groups, values_by_leaf, standard_errors
    )
    all_records = [record for items in grouped.values() for record in items]
    raw_values = [value for value, _ in all_records]
    center = min(all_records, key=lambda record: (record[1], abs(record[0])))[0]
    differences = [value - center for value in raw_values]
    if any(not math.isfinite(value) for value in differences):
        center = min(raw_values) / 2.0 + max(raw_values) / 2.0
        differences = [value - center for value in raw_values]
    if any(not math.isfinite(value) for value in differences):
        raise ValueError("OU trait range exceeds floating-point range; rescale units.")
    size = max(
        max(abs(value) for value in differences),
        max(error for _, error in all_records),
        process_sd,
    )
    trait_exponent = math.frexp(size)[1] - 1 if size > 0.0 else 0
    children: list[list[int]] = [[] for _ in parents]
    for child, parent in enumerate(parents):
        if parent >= 0:
            children[parent].append(child)
    local = [_Factor(None) for _ in parents]
    values, errors = [], []
    for group, records in grouped.items():
        factor = _Factor(None)
        for value, error in records:
            scaled_value = _ldexp(
                value - center, -trait_exponent, "OU centered trait value"
            )
            scaled_error = _ldexp(error, -trait_exponent, "OU standard error")
            variance = scaled_error * scaled_error
            if error > 0.0 and variance == 0.0:
                raise ValueError(
                    "A positive standard-error variance underflowed; rescale units."
                )
            factor = _combine(factor, _Factor(GaussianMarginal(scaled_value, variance)))
            values.append(scaled_value)
            errors.append(scaled_error)
        local[group] = factor
    return _OuData(
        node_groups=node_groups,
        parents=tuple(parents),
        lengths=tuple(lengths),
        children=tuple(tuple(items) for items in children),
        local=tuple(local),
        num_observed=num_observed,
        num_effective_observations=len(all_records),
        num_observed_positions=len(grouped),
        observed_values=tuple(values),
        observed_errors=tuple(errors),
        exact_values=exact_values,
        likelihood_dimensions=len(values),
        trait_center=center,
        trait_exponent=trait_exponent,
    )


def _transition(alpha, sigma2, theta, length):
    exponent = alpha * length
    if not math.isfinite(exponent):
        raise ValueError(
            "alpha multiplied by a branch length exceeds floating-point range."
        )
    decay = math.exp(-exponent)
    root_variance = sigma2 / (2.0 * alpha)
    if not math.isfinite(root_variance) or root_variance <= 0.0:
        raise ValueError(
            "The stationary OU variance is outside floating-point range; rescale units."
        )
    innovation_fraction = -math.expm1(-2.0 * exponent)
    innovation_variance = root_variance * innovation_fraction
    if length > 0.0 and innovation_variance <= 0.0:
        raise ValueError(
            "A positive OU branch variance underflowed; rescale trait or branch units."
        )
    intercept = (-math.expm1(-exponent)) * theta
    if not math.isfinite(intercept):
        raise ValueError("An OU transition mean exceeds floating-point range.")
    return decay, intercept, innovation_variance, root_variance


def _propagate_up(factor, alpha, sigma2, theta, length):
    if factor.marginal is None:
        return factor
    decay, intercept, innovation, _ = _transition(alpha, sigma2, theta, length)
    marginal = factor.marginal
    total_variance = marginal.variance + innovation
    if not math.isfinite(total_variance) or total_variance <= 0.0:
        raise ValueError("An OU pruning variance is not positive and finite.")
    exponent = alpha * length
    log_parent_variance = math.log(total_variance) + 2.0 * exponent
    if decay == 0.0 or log_parent_variance >= math.log(np.finfo(float).max):
        return _Factor(
            None,
            factor.log_weight
            + _log_normal_density(marginal.mean, intercept, total_variance),
        )
    parent_variance = total_variance / decay / decay
    parent_mean = (marginal.mean - intercept) / decay
    if not math.isfinite(parent_variance) or not math.isfinite(parent_mean):
        raise ValueError(
            "An OU pruning message exceeds floating-point range; reduce alpha bounds."
        )
    return _Factor(
        GaussianMarginal(parent_mean, parent_variance),
        factor.log_weight - math.log(decay),
    )


def _propagate_affine_up(factor, alpha, sigma2, theta_reference, length):
    if factor.mean_offset is None:
        return factor
    decay, intercept, innovation, _ = _transition(
        alpha, sigma2, theta_reference, length
    )
    total_variance = factor.variance + innovation
    if not math.isfinite(total_variance) or total_variance <= 0.0:
        raise ValueError("An OU pruning variance is not positive and finite.")
    exponent = alpha * length
    intercept_slope = 1.0 - decay
    log_parent_variance = math.log(total_variance) + 2.0 * exponent
    if decay == 0.0 or log_parent_variance >= math.log(np.finfo(float).max):
        extra = _normal_nll_polynomial(
            factor.mean_offset - intercept,
            factor.mean_slope - intercept_slope,
            total_variance,
        )
        return _AffineFactor(
            None,
            nll_quadratic=factor.nll_quadratic + extra[0],
            nll_linear=factor.nll_linear + extra[1],
            nll_constant=factor.nll_constant + extra[2],
        )
    return _AffineFactor(
        (factor.mean_offset - intercept) / decay,
        (factor.mean_slope - intercept_slope) / decay,
        total_variance / decay / decay,
        factor.nll_quadratic,
        factor.nll_linear,
        factor.nll_constant - exponent,
    )


def _profile_theta(data, alpha, sigma2, theta_reference):
    """Profile theta exactly from its quadratic tree likelihood in one pass."""

    inside = [
        _AffineFactor(
            None if factor.marginal is None else factor.marginal.mean,
            variance=0.0 if factor.marginal is None else factor.marginal.variance,
            nll_constant=-factor.log_weight,
        )
        for factor in data.local
    ]
    upward = [_AffineFactor(None) for _ in data.parents]
    for group in range(len(data.parents) - 1, -1, -1):
        factor = inside[group]
        for child in data.children[group]:
            factor = _combine_affine(factor, upward[child])
        inside[group] = factor
        if data.parents[group] >= 0:
            upward[group] = _propagate_affine_up(
                factor,
                alpha,
                sigma2,
                theta_reference,
                data.lengths[group],
            )
    root_variance = sigma2 / (2.0 * alpha)
    result = _combine_affine(
        _AffineFactor(theta_reference, 1.0, root_variance), inside[0]
    )
    quadratic = result.nll_quadratic
    linear = result.nll_linear
    if not math.isfinite(quadratic) or quadratic <= 0.0:
        raise ValueError("OU theta is not identifiable from the supplied observations.")
    delta = -linear / (2.0 * quadratic)
    fitted_theta = theta_reference + delta
    negative_log_likelihood = (
        quadratic * delta * delta + linear * delta + result.nll_constant
    )
    if not math.isfinite(fitted_theta) or not math.isfinite(negative_log_likelihood):
        raise ValueError("The profiled OU likelihood exceeds floating-point range.")
    return fitted_theta, negative_log_likelihood


def _predict_down(factor, alpha, sigma2, theta, length):
    decay, intercept, innovation, _ = _transition(alpha, sigma2, theta, length)
    if factor.marginal is None:
        return factor
    marginal = factor.marginal
    mean = decay * marginal.mean + intercept
    variance = decay * decay * marginal.variance + innovation
    if not math.isfinite(mean) or not math.isfinite(variance) or variance <= 0.0:
        raise ValueError("An OU smoothing message is not finite and positive.")
    return _Factor(GaussianMarginal(mean, variance), factor.log_weight)


def _prune(data, alpha, sigma2, theta):
    inside = list(data.local)
    upward = [_Factor(None) for _ in data.parents]
    for group in range(len(data.parents) - 1, -1, -1):
        factor = inside[group]
        for child in data.children[group]:
            factor = _combine(factor, upward[child])
        inside[group] = factor
        if data.parents[group] >= 0:
            upward[group] = _propagate_up(
                factor, alpha, sigma2, theta, data.lengths[group]
            )
    root_variance = sigma2 / (2.0 * alpha)
    posterior_root = _combine(
        _Factor(GaussianMarginal(theta, root_variance)), inside[0]
    )
    if not math.isfinite(posterior_root.log_weight):
        raise ValueError("The observed values have zero likelihood under the OU model.")
    return tuple(inside), tuple(upward), posterior_root.log_weight


def _smooth(data, alpha, sigma2, theta, inside, upward):
    root_variance = sigma2 / (2.0 * alpha)
    outside = [_Factor(None) for _ in data.parents]
    outside[0] = _Factor(GaussianMarginal(theta, root_variance))
    posterior = [None for _ in data.parents]
    for group in range(len(data.parents)):
        base = _combine(outside[group], data.local[group])
        children = data.children[group]
        prefix = [base]
        for child in children:
            prefix.append(_combine(prefix[-1], upward[child]))
        posterior[group] = prefix[-1].marginal
        suffix = _Factor(None)
        for index in range(len(children) - 1, -1, -1):
            child = children[index]
            cavity = _combine(prefix[index], suffix)
            outside[child] = _predict_down(
                cavity, alpha, sigma2, theta, data.lengths[child]
            )
            suffix = _combine(upward[child], suffix)
    if any(item is None for item in posterior):
        raise ValueError("Failed to calculate OU ancestral marginals.")
    return tuple(posterior)


def default_alpha_bounds(tree):
    depth_by_node = {tree: 0.0}
    maximum_depth = 0.0
    for node in tree.traverse(strategy="preorder"):
        if node.is_root:
            continue
        depth = depth_by_node[node.up] + float(node.dist)
        depth_by_node[node] = depth
        if node.is_leaf:
            maximum_depth = max(maximum_depth, depth)
    time_scale = maximum_depth if maximum_depth > 0.0 else 1.0
    return 1e-6 / time_scale, 50.0 / time_scale


def parse_alpha_bounds(value, tree):
    if value in (None, ""):
        return default_alpha_bounds(tree)
    items = [item.strip() for item in str(value).split(",")]
    if len(items) != 2:
        raise ValueError("'--alpha-bounds' must contain exactly two values.")
    try:
        lower, upper = (float(item) for item in items)
    except ValueError as exc:
        raise ValueError("'--alpha-bounds' values must be numeric.") from exc
    if not all(math.isfinite(item) and item > 0.0 for item in (lower, upper)):
        raise ValueError("'--alpha-bounds' values must be positive and finite.")
    if lower >= upper:
        raise ValueError(
            "'--alpha-bounds' lower bound must be smaller than the upper bound."
        )
    return lower, upper


def _initial_variance(data):
    values = np.asarray(data.observed_values, dtype=float)
    errors = np.asarray(data.observed_errors, dtype=float)
    variance = float(np.var(values)) + float(np.mean(errors * errors))
    scale = max(
        variance,
        float(np.ptp(values)) ** 2 if len(values) > 1 else 0.0,
        float(np.max(errors)) ** 2 if len(errors) else 0.0,
        max(1.0, abs(float(np.mean(values)))) ** 2 * np.finfo(float).eps,
    )
    return max(variance, scale * 1e-3), scale


def _one_dimensional_intervals(grid_objectives, lower_minimum, upper_minimum):
    local_indices = [
        index
        for index in range(1, len(grid_objectives) - 1)
        if grid_objectives[index] <= grid_objectives[index - 1]
        and grid_objectives[index] <= grid_objectives[index + 1]
        and grid_objectives[index] < 1e99
    ]
    if not local_indices:
        finite_indices = [
            index for index, value in enumerate(grid_objectives) if value < 1e99
        ]
        if finite_indices:
            best_index = min(finite_indices, key=grid_objectives.__getitem__)
            if 0 < best_index < len(grid_objectives) - 1:
                local_indices = [best_index]
    intervals = [
        (index - 1, index + 1, grid_objectives[index]) for index in local_indices
    ]
    if lower_minimum:
        intervals.append((0, 1, grid_objectives[0]))
    if upper_minimum:
        intervals.append(
            (len(grid_objectives) - 2, len(grid_objectives) - 1, grid_objectives[-1])
        )
    return sorted(set(intervals), key=lambda interval: interval[2])[:8]


def _search_ou_one_dimension(parameter_bounds, objective, add_candidate):
    lower, upper = parameter_bounds[0]
    grid = np.linspace(lower, upper, 49)
    grid_objectives = []
    for value in grid:
        add_candidate((value,), "grid", False, "deterministic grid evaluation")
        grid_objectives.append(objective((value,)))
    finite_values = [value for value in grid_objectives if value < 1e99]
    comparison_tolerance = (
        1e-12 * max(1.0, abs(min(finite_values))) if finite_values else 0.0
    )
    lower_minimum = (
        grid_objectives[0] < 1e99
        and grid_objectives[0] <= grid_objectives[1] + comparison_tolerance
    )
    upper_minimum = (
        grid_objectives[-1] < 1e99
        and grid_objectives[-1] <= grid_objectives[-2] + comparison_tolerance
    )
    if lower_minimum:
        add_candidate(
            (grid[0],),
            "boundary",
            True,
            "lower parameter boundary evaluated explicitly",
        )
    if upper_minimum:
        add_candidate(
            (grid[-1],),
            "boundary",
            True,
            "upper parameter boundary evaluated explicitly",
        )
    optimizer_starts = converged_starts = failed_starts = 0
    for lower_index, upper_index, _ in _one_dimensional_intervals(
        grid_objectives, lower_minimum, upper_minimum
    ):
        optimizer_starts += 1
        result = minimize_scalar(
            lambda value: objective((value,)),
            bounds=(float(grid[lower_index]), float(grid[upper_index])),
            method="bounded",
            options={"xatol": 1e-10, "maxiter": 500},
        )
        success = (
            bool(result.success)
            and math.isfinite(float(result.fun))
            and float(result.fun) < 1e99
        )
        if success:
            converged_starts += 1
            add_candidate((result.x,), "bounded", True, result.message)
        else:
            failed_starts += 1
    return optimizer_starts, converged_starts, failed_starts, len(grid)


def _two_dimensional_starts(axes, grid_objectives, initial_parameters):
    ranked_indices = [
        tuple(int(value) for value in index)
        for index in np.ndindex(grid_objectives.shape)
        if grid_objectives[index] < 1e99
    ]
    ranked_indices.sort(key=lambda index: grid_objectives[index])
    start_indices = set(ranked_indices[:4])
    for first_index in range(1, len(axes[0]) - 1):
        for second_index in range(1, len(axes[1]) - 1):
            neighborhood = grid_objectives[
                first_index - 1 : first_index + 2,
                second_index - 1 : second_index + 2,
            ]
            if grid_objectives[first_index, second_index] <= np.min(neighborhood):
                start_indices.add((first_index, second_index))
    for first_index in (0, len(axes[0]) - 1):
        second_index = int(np.argmin(grid_objectives[first_index, :]))
        start_indices.add((first_index, second_index))
    for second_index in (0, len(axes[1]) - 1):
        first_index = int(np.argmin(grid_objectives[:, second_index]))
        start_indices.add((first_index, second_index))
    nearest_initial = tuple(
        int(np.argmin(np.abs(axis - initial)))
        for axis, initial in zip(axes, initial_parameters, strict=True)
    )
    start_indices.add(nearest_initial)
    return sorted(start_indices, key=lambda index: grid_objectives[index])[:10]


def _optimize_ou_start(objective, start, parameter_bounds):
    result = minimize(
        objective,
        start,
        method="L-BFGS-B",
        bounds=parameter_bounds,
    )
    messages = [str(result.message)]
    success = (
        bool(result.success)
        and math.isfinite(float(result.fun))
        and np.all(np.isfinite(result.x))
        and float(result.fun) < 1e99
    )
    return result, success, messages


def _powell_ou_start(objective, start, parameter_bounds, messages):
    fallback = minimize(
        objective,
        start,
        method="Powell",
        bounds=parameter_bounds,
        options={"xtol": 1e-8, "ftol": 1e-10, "maxiter": 500},
    )
    messages.append(f"Powell fallback: {fallback.message}")
    if (
        bool(fallback.success)
        and math.isfinite(float(fallback.fun))
        and np.all(np.isfinite(fallback.x))
        and float(fallback.fun) < 1e99
    ):
        return fallback, messages
    return None, messages


def _search_ou_two_dimensions(
    parameter_bounds, initial_parameters, objective, add_candidate
):
    axes = [np.linspace(lower, upper, 7) for lower, upper in parameter_bounds]
    grid_objectives = np.empty((len(axes[0]), len(axes[1])), dtype=float)
    for first_index, first_value in enumerate(axes[0]):
        for second_index, second_value in enumerate(axes[1]):
            parameters = (first_value, second_value)
            add_candidate(parameters, "grid", False, "deterministic grid evaluation")
            grid_objectives[first_index, second_index] = objective(parameters)
    optimizer_starts = converged_starts = failed_starts = 0
    best_converged_objective = math.inf
    pending_fallbacks = []
    for index in _two_dimensional_starts(axes, grid_objectives, initial_parameters):
        optimizer_starts += 1
        start = np.asarray(
            [axes[dimension][index[dimension]] for dimension in range(2)],
            dtype=float,
        )
        result, success, messages = _optimize_ou_start(
            objective, start, parameter_bounds
        )
        if not success:
            if math.isfinite(float(result.fun)) and np.all(np.isfinite(result.x)):
                add_candidate(
                    result.x,
                    "l-bfgs-b-endpoint",
                    False,
                    "; ".join(messages),
                )
            pending_fallbacks.append((start, result, messages))
            continue
        converged_starts += 1
        best_converged_objective = min(best_converged_objective, float(result.fun))
        add_candidate(result.x, "multistart", True, "; ".join(messages))
    pending_fallbacks.sort(
        key=lambda item: (
            float(item[1].fun) if math.isfinite(float(item[1].fun)) else math.inf
        )
    )
    for start, failed_result, messages in pending_fallbacks:
        comparison_tolerance = 1e-10 * max(1.0, abs(best_converged_objective))
        if (
            math.isfinite(best_converged_objective)
            and float(failed_result.fun)
            >= best_converged_objective - comparison_tolerance
        ):
            failed_starts += 1
            continue
        result, messages = _powell_ou_start(
            objective, start, parameter_bounds, messages
        )
        if result is None:
            failed_starts += 1
            continue
        converged_starts += 1
        best_converged_objective = min(best_converged_objective, float(result.fun))
        add_candidate(result.x, "multistart", True, "; ".join(messages))
    return (
        optimizer_starts,
        converged_starts,
        failed_starts,
        int(grid_objectives.size),
    )


def _fit_parameters(data, alpha, sigma2, theta, alpha_bounds):
    alpha_estimated = alpha is None
    sigma2_estimated = sigma2 is None
    theta_estimated = theta is None
    if alpha is not None:
        alpha = _finite_number(alpha, "--alpha", nonnegative=True)
        if alpha <= 0.0:
            raise ValueError("--alpha must be strictly positive for stationary OU.")
    if sigma2 is not None:
        sigma2 = _finite_number(sigma2, "--sigma2", nonnegative=True)
        if sigma2 <= 0.0:
            raise ValueError("--sigma2 must be strictly positive for stationary OU.")
    if theta is not None:
        theta = _finite_number(theta, "--theta")
    estimated_names = [
        name
        for name, estimated in (
            ("alpha", alpha_estimated),
            ("sigma2", sigma2_estimated),
            ("theta", theta_estimated),
        )
        if estimated
    ]
    if data.num_observed_positions < len(estimated_names):
        raise ValueError(
            "OU free parameters are not separately identifiable: estimating "
            f"{', '.join(estimated_names)} requires at least {len(estimated_names)} "
            "distinct observed tree positions, but only "
            f"{data.num_observed_positions} are available; fix parameters explicitly."
        )

    initial_variance, variance_scale = _initial_variance(data)
    initial_theta = float(np.mean(data.observed_values)) if theta is None else theta
    root_variance_bounds = (
        variance_scale / 1e12,
        variance_scale * 1e12,
    )
    parameter_names: list[str] = []
    parameter_bounds: list[tuple[float, float]] = []
    initial_parameters: list[float] = []
    if alpha_estimated:
        parameter_names.append("log_alpha")
        parameter_bounds.append((math.log(alpha_bounds[0]), math.log(alpha_bounds[1])))
        initial_parameters.append(
            0.5 * (math.log(alpha_bounds[0]) + math.log(alpha_bounds[1]))
        )
    if sigma2_estimated:
        parameter_names.append("log_root_variance")
        parameter_bounds.append(
            (
                math.log(root_variance_bounds[0]),
                math.log(root_variance_bounds[1]),
            )
        )
        initial_parameters.append(math.log(initial_variance))

    evaluation_cache: dict[
        tuple[float, ...], tuple[float, tuple[float, float, float] | None]
    ] = {}
    evaluation_errors: set[str] = set()

    def evaluate(parameters):
        key = tuple(float(value) for value in parameters)
        if key in evaluation_cache:
            return evaluation_cache[key]
        values = dict(zip(parameter_names, key, strict=True))
        try:
            current_alpha = (
                math.exp(values["log_alpha"]) if alpha_estimated else float(alpha)
            )
            current_sigma2 = (
                2.0 * current_alpha * math.exp(values["log_root_variance"])
                if sigma2_estimated
                else float(sigma2)
            )
            if theta_estimated:
                current_theta, negative_log_likelihood = _profile_theta(
                    data, current_alpha, current_sigma2, initial_theta
                )
            else:
                current_theta = float(theta)
                negative_log_likelihood = -_prune(
                    data, current_alpha, current_sigma2, current_theta
                )[2]
            current = (current_alpha, current_sigma2, current_theta)
            if not math.isfinite(negative_log_likelihood):
                raise ValueError("non-finite likelihood")
            evaluation_cache[key] = (float(negative_log_likelihood), current)
        except (ValueError, ArithmeticError) as exc:
            evaluation_errors.add(str(exc))
            evaluation_cache[key] = (1e100, None)
        return evaluation_cache[key]

    def objective(parameters):
        return evaluate(parameters)[0]

    candidates = []

    def add_candidate(parameters, source, success, message):
        values = tuple(float(value) for value in parameters)
        objective_value, current = evaluate(values)
        if current is not None and objective_value < 1e99:
            candidates.append(
                {
                    "objective": objective_value,
                    "parameters": values,
                    "current": current,
                    "source": source,
                    "success": bool(success),
                    "message": str(message),
                }
            )

    optimizer_starts = 0
    converged_starts = 0
    failed_starts = 0
    grid_evaluations = 0
    dimensions = len(parameter_names)
    if dimensions == 0:
        message = (
            "theta profiled exactly" if theta_estimated else "all parameters fixed"
        )
        add_candidate((), "exact", True, message)
    elif dimensions == 1:
        (
            optimizer_starts,
            converged_starts,
            failed_starts,
            grid_evaluations,
        ) = _search_ou_one_dimension(parameter_bounds, objective, add_candidate)
    else:
        (
            optimizer_starts,
            converged_starts,
            failed_starts,
            grid_evaluations,
        ) = _search_ou_two_dimensions(
            parameter_bounds, initial_parameters, objective, add_candidate
        )

    if not candidates:
        failure_summary = "; ".join(sorted(evaluation_errors))
        raise ValueError(
            f"Failed to estimate OU parameters: {failure_summary or 'no finite fit'}"
        )
    if theta_estimated:
        ranked_candidates = sorted(
            range(len(candidates)),
            key=lambda candidate_index: candidates[candidate_index]["objective"],
        )
        profiled_best = candidates[ranked_candidates[0]]["objective"]
        verification_tolerance = 1e-6 * max(1.0, abs(profiled_best))
        verification_indices = set(ranked_candidates[:20])
        verification_indices.update(
            candidate_index
            for candidate_index, candidate in enumerate(candidates)
            if candidate["objective"] <= profiled_best + verification_tolerance
        )
        for candidate_index in verification_indices:
            candidate = candidates[candidate_index]
            candidate["objective"] = -_prune(data, *candidate["current"])[2]
        candidates = [
            candidates[candidate_index]
            for candidate_index in sorted(verification_indices)
        ]
    best_objective = min(candidate["objective"] for candidate in candidates)
    selection_tolerance = 1e-10 * max(1.0, abs(best_objective))
    successful_candidates = [
        candidate
        for candidate in candidates
        if candidate["success"]
        and candidate["objective"] <= best_objective + selection_tolerance
    ]
    selected = min(
        successful_candidates or candidates,
        key=lambda candidate: candidate["objective"],
    )
    fitted = selected["current"]
    fitted_alpha, fitted_sigma2, fitted_theta = fitted
    _, _, log_likelihood = _prune(data, *fitted)
    tolerance = 1e-5
    statuses = []
    if alpha_estimated:
        if fitted_alpha <= alpha_bounds[0] * (1.0 + tolerance):
            statuses.append("alpha_lower_boundary")
        elif fitted_alpha >= alpha_bounds[1] * (1.0 - tolerance):
            statuses.append("alpha_upper_boundary")
    fitted_root_variance = fitted_sigma2 / (2.0 * fitted_alpha)
    if sigma2_estimated:
        if fitted_root_variance <= root_variance_bounds[0] * (1.0 + tolerance):
            statuses.append("root_variance_lower_boundary")
        elif fitted_root_variance >= root_variance_bounds[1] * (1.0 - tolerance):
            statuses.append("root_variance_upper_boundary")
    if not selected["success"]:
        statuses.append("optimizer_grid_fallback")
    fit_status = "+".join(statuses) if statuses else "ok"
    search_summary = (
        f"selected {selected['source']}; {grid_evaluations} grid evaluations; "
        f"{converged_starts}/{optimizer_starts} optimizer starts converged"
    )
    optimizer_message = f"{search_summary}; {selected['message']}"
    return OrnsteinUhlenbeckFit(
        alpha=fitted_alpha,
        alpha_estimated=alpha_estimated,
        theta=fitted_theta,
        theta_estimated=theta_estimated,
        sigma2=fitted_sigma2,
        sigma2_estimated=sigma2_estimated,
        root_variance=fitted_root_variance,
        log_likelihood=log_likelihood,
        num_observed=data.num_observed,
        num_effective_observations=data.num_effective_observations,
        num_observed_positions=data.num_observed_positions,
        fit_status=fit_status,
        optimizer_success=bool(selected["success"]),
        optimizer_message=optimizer_message,
        optimizer_starts=optimizer_starts,
        optimizer_converged_starts=converged_starts,
        optimizer_failed_starts=failed_starts,
        optimizer_grid_evaluations=grid_evaluations,
        alpha_bounds=alpha_bounds,
        root_variance_bounds=root_variance_bounds if sigma2_estimated else None,
    )


def compute_ou_marginals(
    tree,
    values_by_leaf,
    *,
    alpha=None,
    sigma2=None,
    theta=None,
    alpha_bounds=None,
    standard_errors=None,
    _tree_validated=False,
):
    """Return all-node conditional marginals under stationary-root OU1."""

    bounds = default_alpha_bounds(tree) if alpha_bounds is None else alpha_bounds
    if len(bounds) != 2 or not all(
        math.isfinite(float(value)) and float(value) > 0.0 for value in bounds
    ):
        raise ValueError("OU alpha bounds must be two positive finite values.")
    bounds = float(bounds[0]), float(bounds[1])
    if bounds[0] >= bounds[1]:
        raise ValueError("OU alpha lower bound must be smaller than the upper bound.")
    fixed_alpha = (
        None if alpha is None else _finite_number(alpha, "--alpha", nonnegative=True)
    )
    fixed_sigma2 = (
        None if sigma2 is None else _finite_number(sigma2, "--sigma2", nonnegative=True)
    )
    fixed_theta = None if theta is None else _finite_number(theta, "--theta")
    process_sd = 0.0
    if fixed_sigma2 is not None and fixed_sigma2 > 0.0:
        scale_alpha = fixed_alpha if fixed_alpha is not None else bounds[0]
        if scale_alpha > 0.0:
            root_variance = fixed_sigma2 / (2.0 * scale_alpha)
            if math.isfinite(root_variance):
                process_sd = math.sqrt(root_variance)
    data = _prepare_data(
        tree,
        values_by_leaf,
        standard_errors,
        process_sd=process_sd,
        tree_validated=_tree_validated,
    )
    scaled_sigma2 = (
        None
        if fixed_sigma2 is None
        else _ldexp(
            fixed_sigma2,
            -2 * data.trait_exponent,
            "--sigma2",
        )
    )
    if fixed_theta is None:
        scaled_theta = None
    else:
        centered_theta = fixed_theta - data.trait_center
        if not math.isfinite(centered_theta):
            raise ValueError("--theta minus the observed trait center is not finite.")
        scaled_theta = _ldexp(centered_theta, -data.trait_exponent, "OU centered theta")
    fit = _fit_parameters(data, fixed_alpha, scaled_sigma2, scaled_theta, bounds)
    inside, upward, _ = _prune(data, fit.alpha, fit.sigma2, fit.theta)
    smoothed = _smooth(data, fit.alpha, fit.sigma2, fit.theta, inside, upward)
    restored = tuple(
        GaussianMarginal(
            data.trait_center
            + _ldexp(item.mean, data.trait_exponent, "OU ancestral mean"),
            _ldexp(
                item.variance,
                2 * data.trait_exponent,
                "OU ancestral variance",
            ),
        )
        for item in smoothed
    )
    posterior = {
        node: GaussianMarginal(data.exact_values[group], 0.0)
        if group in data.exact_values
        else restored[group]
        for node, group in data.node_groups.items()
    }
    if any(
        not math.isfinite(item.mean)
        or not math.isfinite(item.variance)
        or item.variance < 0.0
        for item in posterior.values()
    ):
        raise ValueError("An OU ancestral marginal is outside floating-point range.")
    restored_theta = data.trait_center + _ldexp(
        fit.theta, data.trait_exponent, "Fitted theta"
    )
    restored_sigma2 = _ldexp(fit.sigma2, 2 * data.trait_exponent, "Fitted sigma2")
    restored_root_variance = _ldexp(
        fit.root_variance,
        2 * data.trait_exponent,
        "Fitted OU root variance",
    )
    restored_variance_bounds = (
        None
        if fit.root_variance_bounds is None
        else (
            _ldexp(
                fit.root_variance_bounds[0],
                2 * data.trait_exponent,
                "OU root-variance bound",
            ),
            _ldexp(
                fit.root_variance_bounds[1],
                2 * data.trait_exponent,
                "OU root-variance bound",
            ),
        )
    )
    fit = replace(
        fit,
        theta=restored_theta,
        sigma2=restored_sigma2,
        root_variance=restored_root_variance,
        root_variance_bounds=restored_variance_bounds,
        log_likelihood=fit.log_likelihood
        - data.likelihood_dimensions * data.trait_exponent * _LOG_2,
    )
    if not math.isfinite(fit.theta) or not math.isfinite(fit.log_likelihood):
        raise ValueError("An OU fitted value is outside floating-point range.")
    return posterior, fit
