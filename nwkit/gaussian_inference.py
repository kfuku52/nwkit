"""Generic scalar affine-Gaussian pruning and smoothing on rooted trees.

Every branch follows ``child = slope * parent + intercept + error``.  The
engine is deliberately independent of biological model names: BM, OU,
time-transformed BM, drift, and regime models differ only in how they build a
``GaussianTreeProcess``.
"""

import math
from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np

from nwkit.compiled_tree import CompiledTree
from nwkit.gaussian_tree import GaussianTreeProcess, _scaled_gaussian_variance

_LOG_2PI = math.log(2.0 * math.pi)
_LOG_2 = math.log(2.0)
_LOG_MAX = math.log(np.finfo(float).max)


@dataclass(frozen=True, slots=True)
class GaussianMarginal:
    mean: float
    variance: float


@dataclass(frozen=True, slots=True)
class GaussianConditioningResult:
    marginals: Mapping[Any, GaussianMarginal]
    log_likelihood: float
    likelihood_rank: int
    num_observed: int
    num_observed_positions: int
    compiled_tree: CompiledTree


@dataclass(frozen=True, slots=True)
class GaussianLikelihoodResult:
    log_likelihood: float
    likelihood_rank: int
    num_observed: int
    num_observed_positions: int
    compiled_tree: CompiledTree


@dataclass(frozen=True, slots=True)
class _Factor:
    """A weighted normalized Gaussian, exact point mass, or constant factor."""

    mean: float | None = None
    variance: float = 0.0
    log_weight: float = 0.0
    density_rank: int = 0


@dataclass(frozen=True, slots=True)
class GaussianSamples:
    nodes: tuple[Any, ...]
    values: np.ndarray


@dataclass(frozen=True, slots=True)
class _MessagePass:
    compiled: CompiledTree
    records: tuple[tuple[int, float, float], ...]
    center: float
    exponent: int
    process: GaussianTreeProcess
    local: tuple[_Factor, ...]
    inside: tuple[_Factor, ...]
    upward: tuple[_Factor, ...]
    root_posterior: _Factor
    observed_positions: int


def _normal_log_density(value: float, mean: float, variance: float) -> float:
    if not math.isfinite(variance) or variance <= 0.0:
        raise ValueError("A Gaussian likelihood variance must be positive and finite.")
    standard_deviation = math.sqrt(variance)
    difference = value - mean
    standardized = (
        difference / standard_deviation
        if math.isfinite(difference)
        else value / standard_deviation - mean / standard_deviation
    )
    if not math.isfinite(standardized):
        return -math.inf
    quadratic = standardized * standardized
    if math.isfinite(quadratic):
        return -0.5 * (_LOG_2PI + math.log(variance) + quadratic)
    return -0.5 * (_LOG_2PI + math.log(variance)) - (0.5 * standardized) * standardized


def _normal_log_overlap(first, second, variance_sum):
    if math.isfinite(variance_sum):
        return _normal_log_density(first.mean, second.mean, variance_sum)
    largest = max(first.variance, second.variance)
    relative_sum = 1.0 + min(first.variance, second.variance) / largest
    standard_deviation = math.sqrt(largest) * math.sqrt(relative_sum)
    standardized = first.mean / standard_deviation - second.mean / standard_deviation
    log_variance = math.log(largest) + math.log(relative_sum)
    return -0.5 * (_LOG_2PI + log_variance) - (0.5 * standardized) * standardized


def _combine(first: _Factor, second: _Factor) -> _Factor:
    if first.mean is None:
        return _Factor(
            second.mean,
            second.variance,
            first.log_weight + second.log_weight,
            first.density_rank + second.density_rank,
        )
    if second.mean is None:
        return _Factor(
            first.mean,
            first.variance,
            first.log_weight + second.log_weight,
            first.density_rank + second.density_rank,
        )
    variance_sum = first.variance + second.variance
    base_weight = first.log_weight + second.log_weight
    base_rank = first.density_rank + second.density_rank
    if variance_sum == 0.0:
        if first.mean != second.mean:
            raise ValueError(
                "Conflicting exact observations have zero likelihood under the "
                "Gaussian tree process."
            )
        return _Factor(first.mean, 0.0, base_weight, base_rank)
    log_overlap = _normal_log_overlap(first, second, variance_sum)
    if not math.isfinite(log_overlap):
        raise ValueError("Observed values have zero Gaussian likelihood.")
    if first.variance == 0.0:
        mean, variance = first.mean, 0.0
    elif second.variance == 0.0:
        mean, variance = second.mean, 0.0
    else:
        larger = max(first.variance, second.variance)
        ratio = min(first.variance, second.variance) / larger
        denominator = 1.0 + ratio
        first_weight = (second.variance / larger) / denominator
        second_weight = (first.variance / larger) / denominator
        mean = first_weight * first.mean + second_weight * second.mean
        variance = min(first.variance, second.variance) / denominator
        if variance <= 0.0 or not math.isfinite(variance):
            raise ValueError("A positive Gaussian variance underflowed during pruning.")
    return _Factor(
        mean,
        variance,
        base_weight + log_overlap,
        base_rank + 1,
    )


def _push_up(
    factor: _Factor, slope: float, intercept: float, innovation: float
) -> _Factor:
    if factor.mean is None:
        return factor
    total = factor.variance + innovation
    if not math.isfinite(total) or total < 0.0:
        raise ValueError("A Gaussian pruning variance is invalid or non-finite.")
    absolute_slope = abs(slope)
    if absolute_slope == 0.0:
        if total == 0.0:
            if factor.mean != intercept:
                raise ValueError(
                    "An exact deterministic branch conflicts with an observation."
                )
            return _Factor(
                None,
                log_weight=factor.log_weight,
                density_rank=factor.density_rank,
            )
        log_density = _normal_log_density(factor.mean, intercept, total)
        if not math.isfinite(log_density):
            raise ValueError("Observed values have zero Gaussian likelihood.")
        return _Factor(
            None,
            log_weight=factor.log_weight + log_density,
            density_rank=factor.density_rank + 1,
        )
    if total == 0.0:
        return _Factor(
            (factor.mean - intercept) / slope,
            0.0,
            factor.log_weight - math.log(absolute_slope),
            factor.density_rank,
        )
    log_parent_variance = math.log(total) - 2.0 * math.log(absolute_slope)
    if log_parent_variance >= _LOG_MAX:
        log_density = _normal_log_density(factor.mean, intercept, total)
        if not math.isfinite(log_density):
            raise ValueError("Observed values have zero Gaussian likelihood.")
        return _Factor(
            None,
            log_weight=factor.log_weight + log_density,
            density_rank=factor.density_rank + 1,
        )
    parent_mean = (factor.mean - intercept) / slope
    parent_variance = total / slope / slope
    if not math.isfinite(parent_mean) or not math.isfinite(parent_variance):
        raise ValueError("A Gaussian pruning message exceeds floating-point range.")
    return _Factor(
        parent_mean,
        parent_variance,
        factor.log_weight - math.log(absolute_slope),
        factor.density_rank,
    )


def _push_down(
    factor: _Factor, slope: float, intercept: float, innovation: float
) -> _Factor:
    if factor.mean is None:
        if slope != 0.0:
            return factor
        return _Factor(intercept, innovation, factor.log_weight, factor.density_rank)
    variance = _scaled_gaussian_variance(slope, factor.variance) + innovation
    mean = slope * factor.mean + intercept
    if not math.isfinite(mean) or not math.isfinite(variance) or variance < 0.0:
        raise ValueError("A Gaussian smoothing message exceeds floating-point range.")
    return _Factor(mean, variance, factor.log_weight, factor.density_rank)


def _finite_number(value, label: str, *, nonnegative: bool = False) -> float:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{label} must be numeric, not boolean.")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{label} must be numeric and finite.") from exc
    if not math.isfinite(result) or (nonnegative and result < 0.0):
        qualifier = "non-negative and finite" if nonnegative else "finite"
        raise ValueError(f"{label} must be {qualifier}.")
    return result


def _scaled_process_and_observations(
    process: GaussianTreeProcess,
    compiled: CompiledTree,
    values_by_leaf,
    standard_errors,
):
    records = []
    for index in compiled.leaf_indices:
        name = str(compiled.nodes[index].name)
        value = values_by_leaf.get(name)
        if value is None:
            continue
        observed = _finite_number(value, f"Trait value for '{name}'")
        if standard_errors is not None and name not in standard_errors:
            raise ValueError(f"A standard error is required for observed tip '{name}'.")
        error = _finite_number(
            0.0 if standard_errors is None else standard_errors[name],
            f"Standard error for '{name}'",
            nonnegative=True,
        )
        records.append((index, observed, error))

    if not records:
        return records, 0.0, 0, process
    center = min(records, key=lambda item: (item[2], abs(item[1])))[1]
    differences = [value - center for _, value, _ in records]
    if any(not math.isfinite(value) for value in differences):
        midpoint = (
            min(value for _, value, _ in records) / 2.0
            + max(value for _, value, _ in records) / 2.0
        )
        center = midpoint
        differences = [value - center for _, value, _ in records]
    if any(not math.isfinite(value) for value in differences):
        raise ValueError("Trait range exceeds floating-point range; rescale units.")

    magnitudes = [abs(value) for value in differences]
    magnitudes.extend(error for _, _, error in records)
    for transition in process.transitions.values():
        if transition.variance:
            magnitudes.append(math.sqrt(transition.variance))
        offset = math.fsum(((transition.slope - 1.0) * center, transition.intercept))
        if math.isfinite(offset):
            magnitudes.append(abs(offset))
    if process.root.variance:
        magnitudes.append(math.sqrt(process.root.variance))
    root_offset = process.root.mean - center
    if math.isfinite(root_offset):
        magnitudes.append(abs(root_offset))
    size = max(magnitudes, default=0.0)
    exponent = math.frexp(size)[1] - 1 if size > 0.0 else 0

    def scaled(value, power=1):
        try:
            result = math.ldexp(value, -power * exponent)
        except OverflowError as exc:
            raise ValueError(
                "Gaussian process units exceed floating-point range; rescale units."
            ) from exc
        if not math.isfinite(result) or (value != 0.0 and result == 0.0):
            raise ValueError(
                "Gaussian process units exceed floating-point range; rescale units."
            )
        return result

    scaled_center = scaled(center)
    scaled_records = [
        (index, scaled(value) - scaled_center, scaled(error))
        for index, value, error in records
    ]
    transitions = {}
    for node, transition in process.transitions.items():
        transformed_intercept = math.fsum(
            (
                (transition.slope - 1.0) * scaled_center,
                scaled(transition.intercept),
            )
        )
        transitions[node] = type(transition)(
            transition.slope,
            transformed_intercept,
            scaled(transition.variance, power=2),
        )
    root = type(process.root)(
        process.root.mode,
        scaled(process.root.mean) - scaled_center,
        None
        if process.root.variance is None
        else scaled(process.root.variance, power=2),
    )
    scaled_process = type(process)(
        tree=process.tree,
        transitions=transitions,
        root=root,
        model=process.model,
        parameter=process.parameter,
    )
    return scaled_records, center, exponent, scaled_process


def _message_pass(
    process,
    values_by_leaf,
    standard_errors,
    compiled_tree,
    *,
    retain_inside=False,
    likelihood_only=False,
):
    compiled = (
        CompiledTree.from_tree(process.tree) if compiled_tree is None else compiled_tree
    )
    compiled.require_tree(process.tree)
    records, center, exponent, scaled_process = _scaled_process_and_observations(
        process, compiled, values_by_leaf, standard_errors
    )
    if not records and scaled_process.root.mode == "flat":
        raise ValueError(
            "A flat-root Gaussian process requires at least one observation."
        )
    local = [_Factor()] * len(compiled.nodes)
    latent_positions = list(range(len(compiled.nodes)))
    for index in range(1, len(compiled.nodes)):
        transition = scaled_process.transitions[compiled.nodes[index]]
        if (
            transition.slope == 1.0
            and transition.intercept == 0.0
            and transition.variance == 0.0
        ):
            latent_positions[index] = latent_positions[compiled.parents[index]]
    observed_positions = set()
    for index, value, error in records:
        variance = error * error
        if error > 0.0 and variance == 0.0:
            raise ValueError("A positive observation variance underflowed.")
        local[index] = _combine(local[index], _Factor(value, variance))
        observed_positions.add(latent_positions[index])
    inside = local if likelihood_only else list(local)
    upward = _prune_gaussian_messages(
        scaled_process, compiled, inside, likelihood_only=likelihood_only
    )
    root = scaled_process.root
    if root.mode == "flat":
        if inside[0].mean is None:
            raise ValueError(
                "The flat Gaussian root is not identifiable from observations."
            )
        root_posterior = inside[0]
    else:
        assert root.variance is not None
        root_posterior = _combine(_Factor(root.mean, root.variance), inside[0])
    if not math.isfinite(root_posterior.log_weight):
        raise ValueError(
            "Observed values have zero likelihood under the Gaussian process."
        )
    return _MessagePass(
        compiled,
        tuple(records),
        center,
        exponent,
        scaled_process,
        tuple(local) if not likelihood_only else (),
        tuple(inside) if retain_inside else (),
        tuple(upward) if not likelihood_only else (),
        root_posterior,
        len(observed_positions),
    )


def _prune_gaussian_messages(process, compiled, inside, *, likelihood_only):
    empty = _Factor()
    upward = inside if likelihood_only else [empty] * len(compiled.nodes)
    for index in compiled.postorder:
        factor = inside[index]
        for child in compiled.children[index]:
            factor = _combine(factor, upward[child])
        inside[index] = factor
        parent = compiled.parents[index]
        if parent >= 0:
            transition = process.transitions[compiled.nodes[index]]
            upward[index] = _push_up(
                factor,
                transition.slope,
                transition.intercept,
                transition.variance,
            )
    return upward


def condition_gaussian_tree(
    process: GaussianTreeProcess,
    values_by_leaf: Mapping[str, float | None],
    *,
    standard_errors: Mapping[str, float] | None = None,
    compiled_tree: CompiledTree | None = None,
) -> GaussianConditioningResult:
    """Condition every process node on observed tip values in linear time."""

    messages = _message_pass(process, values_by_leaf, standard_errors, compiled_tree)
    compiled = messages.compiled
    records = messages.records
    center = messages.center
    exponent = messages.exponent
    scaled_process = messages.process
    local = messages.local
    upward = messages.upward
    root_posterior = messages.root_posterior
    root = scaled_process.root

    outside = [_Factor() for _ in compiled.nodes]
    outside[0] = (
        _Factor() if root.mode == "flat" else _Factor(root.mean, root.variance or 0.0)
    )
    posterior: list[_Factor | None] = [None] * len(compiled.nodes)
    for index in range(len(compiled.nodes)):
        base = _combine(outside[index], local[index])
        children = compiled.children[index]
        prefix = [base]
        for child in children:
            prefix.append(_combine(prefix[-1], upward[child]))
        posterior[index] = prefix[-1]
        suffix = _Factor()
        for position in range(len(children) - 1, -1, -1):
            child = children[position]
            cavity = _combine(prefix[position], suffix)
            transition = scaled_process.transitions[compiled.nodes[child]]
            outside[child] = _push_down(
                cavity,
                transition.slope,
                transition.intercept,
                transition.variance,
            )
            suffix = _combine(upward[child], suffix)

    if any(factor is None or factor.mean is None for factor in posterior):
        raise ValueError(
            "Some Gaussian process nodes have improper posterior marginals."
        )
    final_posterior = tuple(
        factor for factor in posterior if factor is not None and factor.mean is not None
    )

    def restore_mean(value):
        try:
            result = center + math.ldexp(value, exponent)
        except OverflowError as exc:
            raise ValueError(
                "A posterior Gaussian mean exceeds floating-point range."
            ) from exc
        if not math.isfinite(result):
            raise ValueError("A posterior Gaussian mean exceeds floating-point range.")
        return result

    def restore_variance(value):
        try:
            result = math.ldexp(value, 2 * exponent)
        except OverflowError as exc:
            raise ValueError(
                "A posterior Gaussian variance exceeds floating-point range."
            ) from exc
        if not math.isfinite(result) or result < 0.0:
            raise ValueError(
                "A posterior Gaussian variance exceeds floating-point range."
            )
        return result

    marginals = {
        node: GaussianMarginal(
            restore_mean(final_posterior[index].mean),
            restore_variance(final_posterior[index].variance),
        )
        for index, node in enumerate(compiled.nodes)
    }
    rank = root_posterior.density_rank
    log_likelihood = root_posterior.log_weight - rank * exponent * _LOG_2
    return GaussianConditioningResult(
        marginals=marginals,
        log_likelihood=log_likelihood,
        likelihood_rank=rank,
        num_observed=len(records),
        num_observed_positions=messages.observed_positions,
        compiled_tree=compiled,
    )


def gaussian_tree_likelihood(
    process: GaussianTreeProcess,
    values_by_leaf: Mapping[str, float | None],
    *,
    standard_errors: Mapping[str, float] | None = None,
    compiled_tree: CompiledTree | None = None,
) -> GaussianLikelihoodResult:
    """Evaluate a Gaussian tree likelihood without the smoothing pass."""

    messages = _message_pass(
        process, values_by_leaf, standard_errors, compiled_tree, likelihood_only=True
    )
    rank = messages.root_posterior.density_rank
    log_likelihood = (
        messages.root_posterior.log_weight - rank * messages.exponent * _LOG_2
    )
    return GaussianLikelihoodResult(
        log_likelihood=log_likelihood,
        likelihood_rank=rank,
        num_observed=len(messages.records),
        num_observed_positions=messages.observed_positions,
        compiled_tree=messages.compiled,
    )


def gaussian_process_parameter_rank(
    process_for, parameters, observed_nodes, bounds, *, trait_exponent=0
):
    """Numerically rank free parameters through observed means and covariances."""

    parameters = np.asarray(parameters, dtype=float)
    observed_nodes = tuple(observed_nodes)
    bounds = tuple(bounds)
    if parameters.ndim != 1 or len(parameters) != len(bounds):
        raise ValueError("Gaussian parameter-rank inputs have inconsistent sizes.")
    if not len(parameters):
        return 0
    if not observed_nodes:
        return 0

    triangle = np.triu_indices(len(observed_nodes))

    def distribution_summary(values):
        process = process_for(values)
        means, _ = process.marginal_moments()
        mean_vector = np.asarray([means[node] for node in observed_nodes], dtype=float)
        covariance = process.covariance(observed_nodes)
        result = np.concatenate(
            (
                np.ldexp(mean_vector, -trait_exponent),
                np.ldexp(covariance[triangle], -2 * trait_exponent),
            )
        )
        if np.any(~np.isfinite(result)):
            raise ValueError("Gaussian parameter-rank summary is not finite.")
        return result

    derivatives = []
    for index, value in enumerate(parameters):
        step = 1e-5 * max(1.0, abs(float(value)))
        lower, upper = bounds[index]
        left = float(value) - step
        right = float(value) + step
        if lower is not None:
            left = max(left, float(lower))
        if upper is not None:
            right = min(right, float(upper))
        if not left < right:
            derivatives.append(np.zeros_like(distribution_summary(parameters)))
            continue
        left_values = parameters.copy()
        right_values = parameters.copy()
        left_values[index] = left
        right_values[index] = right
        derivatives.append(
            (distribution_summary(right_values) - distribution_summary(left_values))
            / (right - left)
        )
    jacobian = np.column_stack(derivatives)
    row_scales = np.max(np.abs(jacobian), axis=1)
    informative = np.concatenate(
        (
            _informative_rank_rows(row_scales[: len(observed_nodes)]),
            _informative_rank_rows(row_scales[len(observed_nodes) :]),
        )
    )
    if not np.any(informative):
        return 0
    normalized = jacobian[informative] / row_scales[informative, None]
    column_scales = np.linalg.norm(normalized, axis=0)
    nonzero = column_scales > 1e-10
    if not np.any(nonzero):
        return 0
    normalized[:, nonzero] /= column_scales[nonzero]
    return _relative_parameter_rank(normalized[:, nonzero])


def _informative_rank_rows(scales):
    # Means and covariances have different dimensions. Apply the numerical
    # threshold within each block so a change of units cannot erase one block.
    largest = float(np.max(scales, initial=0.0))
    return scales > largest * 1e-10


def _relative_parameter_rank(normalized):
    singular_values = np.linalg.svd(normalized, compute_uv=False)
    if not len(singular_values) or singular_values[0] == 0.0:
        return 0
    # Rows contain all observed means and covariance entries, so their count grows
    # quadratically with the number of tips.  It must not inflate the relative
    # finite-difference tolerance and make an otherwise unchanged model lose rank.
    tolerance = max(1, normalized.shape[1]) * 1e-7 * float(singular_values[0])
    return int(np.sum(singular_values > tolerance))


def _posterior_draw_parameters(transition, inside):
    """Return conditional affine weights, offset and standard deviation."""

    variance = transition.variance
    if inside.mean is None:
        return 1.0, 0.0, math.sqrt(variance)
    if inside.variance == 0.0:
        # Pruning already checked compatibility, including exact observations
        # propagated back through deterministic affine branches.
        return 0.0, inside.mean, 0.0
    if variance == 0.0:
        return 1.0, 0.0, 0.0
    larger = max(variance, inside.variance)
    smaller = min(variance, inside.variance)
    denominator = 1.0 + smaller / larger
    parent_weight = (inside.variance / larger) / denominator
    inside_weight = (variance / larger) / denominator
    conditional_variance = smaller / denominator
    if conditional_variance <= 0.0:
        raise ValueError("A positive Gaussian variance underflowed during sampling.")
    return parent_weight, inside_weight * inside.mean, math.sqrt(conditional_variance)


def _draw_gaussian_messages(messages, num_samples, rng):
    compiled = messages.compiled
    parameters = [
        _posterior_draw_parameters(
            messages.process.transitions[compiled.nodes[index]], messages.inside[index]
        )
        for index in range(1, len(compiled.nodes))
    ]
    root = messages.root_posterior
    assert root.mean is not None
    root_sd = math.sqrt(root.variance)
    stochastic_count = int(root_sd > 0.0) + sum(sd > 0.0 for _, _, sd in parameters)
    draws = np.empty((num_samples, len(compiled.nodes)), dtype=float)
    # Keep the original sample-major random-number order. Batch the normals
    # so temporary storage is bounded even for a large number of samples.
    batch_size = max(1, 131_072 // max(1, stochastic_count))
    for start in range(0, num_samples, batch_size):
        batch = draws[start : start + batch_size]
        noise = rng.standard_normal((len(batch), stochastic_count)).T
        next_noise = iter(noise)
        batch[:, 0] = root.mean
        if root_sd > 0.0:
            batch[:, 0] += root_sd * next(next_noise)
        for index, (weight, offset, sd) in enumerate(parameters, start=1):
            transition = messages.process.transitions[compiled.nodes[index]]
            parent = compiled.parents[index]
            if weight == 0.0:
                batch[:, index] = offset
            else:
                prior = transition.slope * batch[:, parent] + transition.intercept
                batch[:, index] = weight * prior + offset
            if sd > 0.0:
                batch[:, index] += sd * next(next_noise)
    return draws


def sample_gaussian_posterior(
    process: GaussianTreeProcess,
    values_by_leaf: Mapping[str, float | None],
    *,
    num_samples: int,
    standard_errors: Mapping[str, float] | None = None,
    seed: int | None = None,
    compiled_tree: CompiledTree | None = None,
) -> GaussianSamples:
    """Draw joint all-node samples using tree forward-filter sampling."""

    if isinstance(num_samples, bool) or int(num_samples) != num_samples:
        raise ValueError("num_samples must be a positive integer.")
    num_samples = int(num_samples)
    if num_samples <= 0:
        raise ValueError("num_samples must be a positive integer.")
    messages = _message_pass(
        process,
        values_by_leaf,
        standard_errors,
        compiled_tree,
        retain_inside=True,
    )
    draws = _draw_gaussian_messages(messages, num_samples, np.random.default_rng(seed))
    np.ldexp(draws, messages.exponent, out=draws)
    draws += messages.center
    if np.any(~np.isfinite(draws)):
        raise ValueError("A posterior Gaussian sample exceeds floating-point range.")
    return GaussianSamples(messages.compiled.nodes, draws)


def simulate_gaussian_process(
    process: GaussianTreeProcess,
    *,
    num_samples: int,
    seed: int | None = None,
    root_values=None,
    compiled_tree: CompiledTree | None = None,
    _output_indices: tuple[int, ...] | None = None,
) -> GaussianSamples:
    """Simulate joint process realizations in original trait units."""

    if isinstance(num_samples, bool) or int(num_samples) != num_samples:
        raise ValueError("num_samples must be a positive integer.")
    num_samples = int(num_samples)
    if num_samples <= 0:
        raise ValueError("num_samples must be a positive integer.")
    compiled = (
        CompiledTree.from_tree(process.tree) if compiled_tree is None else compiled_tree
    )
    compiled.require_tree(process.tree)
    rng = np.random.default_rng(seed)
    root_draws = _simulation_root_draws(process, root_values, num_samples, rng)
    if _output_indices is not None:
        return _simulate_selected_nodes(
            process, compiled, root_draws, rng, _output_indices
        )
    if num_samples == 1:
        return _simulate_single_draw(process, compiled, root_draws, rng)
    draws = np.empty((num_samples, len(compiled.nodes)), dtype=float)
    draws[:, 0] = root_draws
    for node_index in range(1, len(compiled.nodes)):
        parent = compiled.parents[node_index]
        transition = process.transitions[compiled.nodes[node_index]]
        _simulate_branch(transition, draws[:, parent], draws[:, node_index], rng)
    if np.any(~np.isfinite(draws)):
        raise ValueError("A Gaussian process simulation exceeds floating-point range.")
    return GaussianSamples(compiled.nodes, draws)


def _simulate_single_draw(process, compiled, root, rng):
    draws = np.empty(len(compiled.nodes), dtype=float)
    draws[0] = root[0]
    for index, node in enumerate(compiled.nodes[1:], start=1):
        transition = process.transitions[node]
        value = transition.slope * float(draws[compiled.parents[index]])
        value += transition.intercept
        if transition.variance > 0.0:
            value += rng.normal(0.0, math.sqrt(transition.variance))
        draws[index] = value
    if not np.isfinite(draws).all():
        raise ValueError("A Gaussian process simulation exceeds floating-point range.")
    return GaussianSamples(compiled.nodes, draws[None, :])


def _simulation_root_draws(process, root_values, num_samples, rng):
    if root_values is None:
        if process.root.mode == "flat":
            raise ValueError("A flat-root process simulation requires root_values.")
        assert process.root.variance is not None
        if process.root.variance == 0.0:
            return np.full(num_samples, process.root.mean)
        return rng.normal(
            process.root.mean, math.sqrt(process.root.variance), num_samples
        )
    values = np.asarray(root_values, dtype=float)
    if values.ndim == 0:
        return np.full(num_samples, float(values))
    if values.shape == (num_samples,):
        return values
    raise ValueError("root_values must be scalar or have num_samples values.")


def _simulate_branch(transition, parent_draws, child_draws, rng):
    child_draws[:] = transition.slope * parent_draws + transition.intercept
    if transition.variance > 0.0:
        child_draws += rng.normal(0.0, math.sqrt(transition.variance), len(child_draws))


def _simulate_selected_nodes(process, compiled, root, rng, output_indices):
    outputs = np.empty((len(root), len(output_indices)), dtype=float)
    output_rows = {index: row for row, index in enumerate(output_indices)}
    if len(output_rows) != len(output_indices) or any(
        not 0 <= index < len(compiled.nodes) for index in output_rows
    ):
        raise ValueError("Gaussian simulation node indices must be distinct and valid.")
    remaining = [len(children) for children in compiled.children]
    # Keep a parent's draws until its last child, and fill selected tip rows
    # directly. Traversal and random-number consumption match the full draws.
    active = {0: root}
    if 0 in output_rows:
        outputs[:, output_rows[0]] = root
    finite = bool(np.isfinite(root).all())
    for index, parent in enumerate(compiled.parents[1:], start=1):
        selected = output_rows.get(index)
        draws = np.empty_like(root) if selected is None else outputs[:, selected]
        _simulate_branch(
            process.transitions[compiled.nodes[index]], active[parent], draws, rng
        )
        finite &= bool(np.isfinite(draws).all())
        remaining[parent] -= 1
        if remaining[parent] == 0:
            del active[parent]
        if remaining[index]:
            active[index] = draws
    if not finite:
        raise ValueError("A Gaussian process simulation exceeds floating-point range.")
    return GaussianSamples(
        tuple(compiled.nodes[index] for index in output_indices), outputs
    )
