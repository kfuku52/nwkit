"""Bayesian threshold/liability ancestral reconstruction.

The latent liability follows a standardized Brownian process with
``root ~ N(0, 1)`` and diffusion variance one.  This fixes the otherwise
unidentifiable location/scale.  Binary traits use threshold zero.  For ordinal
traits the first threshold is zero and remaining thresholds are sampled, or
all thresholds may be supplied explicitly.
"""

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.special import ndtr, ndtri
from scipy.stats import truncnorm

from nwkit.compiled_tree import CompiledTree
from nwkit.gaussian_inference import GaussianMarginal
from nwkit.gaussian_tree import (
    GaussianRootPrior,
    GaussianTreeProcess,
    brownian_transition,
)
from nwkit.rooting_state import require_rooted
from nwkit.util import assign_branch_ids, get_node_class, validate_unique_named_leaves


@dataclass(frozen=True, slots=True)
class ThresholdFit:
    thresholds: tuple[float, ...]
    thresholds_estimated: bool
    num_samples: int
    burnin: int
    thin: int
    chains: int
    seed: int | None
    rhat_max: float
    ess_min: float
    fit_status: str
    liability_marginals: dict


def parse_thresholds(value, num_states):
    """Parse fixed thresholds, or return ``None`` for identified estimation."""

    if value in (None, ""):
        return None
    items = [item.strip() for item in str(value).split(",")]
    if len(items) != num_states - 1:
        raise ValueError(
            f"--thresholds requires {num_states - 1} increasing values for "
            f"{num_states} states."
        )
    try:
        thresholds = tuple(float(item) for item in items)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("--thresholds values must be numeric and finite.") from exc
    if not all(math.isfinite(item) for item in thresholds) or any(
        right <= left for left, right in zip(thresholds, thresholds[1:], strict=False)
    ):
        raise ValueError("--thresholds values must be finite and strictly increasing.")
    return thresholds


def build_threshold_process(tree):
    """Return the identified unit-diffusion Brownian liability process."""

    transitions = {}
    for node in tree.traverse(strategy="preorder"):
        if node.is_root:
            continue
        try:
            length = float(node.dist)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("Threshold-model branch lengths must be numeric.") from exc
        if not math.isfinite(length) or length <= 0.0:
            raise ValueError(
                "Threshold ASR requires positive finite non-root branch lengths."
            )
        transitions[node] = brownian_transition(length)
    return GaussianTreeProcess(
        tree=tree,
        transitions=transitions,
        root=GaussianRootPrior("gaussian", 0.0, 1.0),
        model="threshold",
    )


def _positive_integer(value, label):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{label} must be a positive integer.")
    try:
        number = float(value)
        integer = int(number)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{label} must be a positive integer.") from exc
    if not math.isfinite(number) or integer != number or integer <= 0:
        raise ValueError(f"{label} must be a positive integer.")
    return integer


def _nonnegative_integer(value, label):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{label} must be a non-negative integer.")
    try:
        number = float(value)
        integer = int(number)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{label} must be a non-negative integer.") from exc
    if not math.isfinite(number) or integer != number or integer < 0:
        raise ValueError(f"{label} must be a non-negative integer.")
    return integer


def _initial_thresholds(category_indices, num_states):
    if num_states == 2:
        return np.asarray([0.0], dtype=float)
    counts = np.bincount(category_indices, minlength=num_states).astype(float)
    cumulative = np.cumsum(counts)[:-1] / float(np.sum(counts))
    quantiles = ndtri(np.clip(cumulative, 1e-3, 1.0 - 1e-3))
    quantiles -= quantiles[0]
    result = np.empty(num_states - 1, dtype=float)
    result[0] = 0.0
    for index in range(1, len(result)):
        result[index] = max(float(quantiles[index]), result[index - 1] + 0.25)
    return result


def _category_interval(category, thresholds):
    lower = -math.inf if category == 0 else float(thresholds[category - 1])
    upper = math.inf if category == len(thresholds) else float(thresholds[category])
    return lower, upper


def _allowed_intervals(allowed, thresholds):
    intervals = []
    start = previous = int(allowed[0])
    for category in allowed[1:]:
        category = int(category)
        if category == previous + 1:
            previous = category
            continue
        lower = -math.inf if start == 0 else float(thresholds[start - 1])
        upper = math.inf if previous == len(thresholds) else float(thresholds[previous])
        intervals.append((lower, upper))
        start = previous = category
    lower = -math.inf if start == 0 else float(thresholds[start - 1])
    upper = math.inf if previous == len(thresholds) else float(thresholds[previous])
    intervals.append((lower, upper))
    return intervals


def _truncated_draw(mean, sd, allowed, thresholds, rng):
    intervals = _allowed_intervals(allowed, thresholds)
    masses = []
    for lower, upper in intervals:
        z_lower = -math.inf if lower == -math.inf else (lower - mean) / sd
        z_upper = math.inf if upper == math.inf else (upper - mean) / sd
        masses.append(max(0.0, float(ndtr(z_upper) - ndtr(z_lower))))
    total = math.fsum(masses)
    if total > 0.0 and math.isfinite(total):
        interval = intervals[
            int(rng.choice(len(intervals), p=np.asarray(masses) / total))
        ]
    else:
        interval = min(
            intervals,
            key=lambda bounds: (
                0.0
                if bounds[0] <= mean <= bounds[1]
                else min(abs(mean - bounds[0]), abs(mean - bounds[1]))
            ),
        )
    lower, upper = interval
    a = -math.inf if lower == -math.inf else (lower - mean) / sd
    b = math.inf if upper == math.inf else (upper - mean) / sd
    return float(truncnorm.rvs(a, b, loc=mean, scale=sd, random_state=rng))


def _conditional_parameters(index, values, compiled, process):
    node = compiled.nodes[index]
    if index == 0:
        assert process.root.variance is not None
        precision = 1.0 / process.root.variance
        numerator = process.root.mean / process.root.variance
    else:
        transition = process.transitions[node]
        parent = compiled.parents[index]
        precision = 1.0 / transition.variance
        numerator = (
            transition.slope * values[parent] + transition.intercept
        ) / transition.variance
    for child in compiled.children[index]:
        transition = process.transitions[compiled.nodes[child]]
        precision += transition.slope * transition.slope / transition.variance
        numerator += (
            transition.slope
            * (values[child] - transition.intercept)
            / transition.variance
        )
    variance = 1.0 / precision
    return numerator * variance, math.sqrt(variance)


def _initialize_values(compiled, constraints, thresholds):
    values = np.zeros(len(compiled.nodes), dtype=float)
    for index, allowed in constraints.items():
        candidates = []
        for category in allowed:
            lower, upper = _category_interval(int(category), thresholds)
            if math.isfinite(lower) and math.isfinite(upper):
                candidate = (lower + upper) / 2.0
            elif math.isfinite(lower):
                candidate = lower + 1.0
            elif math.isfinite(upper):
                candidate = upper - 1.0
            else:
                candidate = 0.0
            candidates.append(candidate)
        values[index] = min(candidates, key=abs)
    return values


def _update_thresholds(values, allowed_categories, thresholds, rng):
    if len(thresholds) <= 1:
        return
    scale = max(1.0, float(np.max(np.abs(values))))
    epsilon = np.finfo(float).eps * scale * 16.0
    for threshold_index in range(1, len(thresholds)):
        previous = float(thresholds[threshold_index - 1])
        following = (
            float(thresholds[threshold_index + 1])
            if threshold_index + 1 < len(thresholds)
            else math.inf
        )
        lower = previous + epsilon
        upper = following - epsilon
        for value, allowed in zip(values, allowed_categories, strict=True):
            value = float(value)
            if value <= previous or value > following:
                continue
            lower_category_allowed = threshold_index in allowed
            upper_category_allowed = threshold_index + 1 in allowed
            if lower_category_allowed and not upper_category_allowed:
                lower = max(lower, value)
            elif upper_category_allowed and not lower_category_allowed:
                upper = min(upper, value)
            elif not lower_category_allowed and not upper_category_allowed:
                raise ValueError(
                    "Ordinal threshold constraints are internally inconsistent."
                )
        if threshold_index + 1 < len(thresholds):
            upper = min(upper, float(thresholds[threshold_index + 1]) - epsilon)
        if not lower < upper:
            raise ValueError(
                "Ordinal threshold update collapsed; increase MCMC burn-in or fix --thresholds."
            )
        thresholds[threshold_index] = rng.uniform(lower, upper)


def _rhat_and_ess(traces):
    traces = np.asarray(traces, dtype=float)
    chains, samples, dimensions = traces.shape
    rhat = np.full(dimensions, np.nan, dtype=float)
    ess = np.full(dimensions, chains * samples, dtype=float)
    for dimension in range(dimensions):
        values = traces[:, :, dimension]
        if chains > 1 and samples > 1:
            within = float(np.mean(np.var(values, axis=1, ddof=1)))
            between = samples * float(np.var(np.mean(values, axis=1), ddof=1))
            if within > 0.0:
                pooled = ((samples - 1.0) / samples) * within + between / samples
                rhat[dimension] = math.sqrt(max(0.0, pooled / within))
            else:
                rhat[dimension] = 1.0 if between == 0.0 else math.inf
        if samples > 2:
            centered = values - np.mean(values, axis=1, keepdims=True)
            denominator = float(np.sum(centered * centered))
            if denominator > 0.0:
                rho = float(np.sum(centered[:, :-1] * centered[:, 1:]) / denominator)
                rho = min(0.99, max(-0.99, rho))
                ess[dimension] = min(
                    chains * samples,
                    chains * samples * (1.0 - rho) / (1.0 + rho),
                )
    return rhat, ess


def compute_threshold_marginals(
    tree,
    states,
    observed_state_by_leaf,
    likelihood_by_leaf,
    *,
    thresholds=None,
    num_samples=1000,
    burnin=500,
    thin=1,
    chains=4,
    seed=None,
    _tree_validated=False,
):
    """Sample posterior liabilities and return category probabilities."""

    states = tuple(str(state) for state in states)
    if len(states) < 2:
        raise ValueError("Threshold ASR requires at least two ordered states.")
    if len(states) != len(set(states)):
        raise ValueError("Threshold ASR states must be unique.")
    if not _tree_validated:
        validate_unique_named_leaves(tree, option_name="--infile", context=" for 'asr'")
        require_rooted(tree, "ASR requires a rooted tree.")
    num_samples = _positive_integer(num_samples, "--liability-samples")
    burnin = _nonnegative_integer(burnin, "--liability-burnin")
    thin = _positive_integer(thin, "--liability-thin")
    chains = _positive_integer(chains, "--liability-chains")
    fixed_thresholds = parse_thresholds(thresholds, len(states))
    compiled = CompiledTree.from_tree(tree)
    constraints = {}
    observed_categories = []
    for name, node_index in compiled.leaf_index_by_name.items():
        if observed_state_by_leaf.get(name) is None:
            continue
        likelihood = np.asarray(likelihood_by_leaf[name], dtype=float)
        allowed = np.flatnonzero(likelihood > 0.0)
        if not len(allowed):
            raise ValueError(f"Threshold observation for '{name}' allows no state.")
        constraints[node_index] = allowed
        if len(allowed) != 1:
            continue
        observed_categories.append(int(allowed[0]))
    if len(observed_categories) < 2:
        raise ValueError(
            "Threshold ASR requires at least two unambiguous observations."
        )
    observed_categories_array = np.asarray(observed_categories, dtype=int)
    if fixed_thresholds is None and len(states) > 2:
        missing_states = [
            states[index]
            for index in range(len(states))
            if index not in set(observed_categories)
        ]
        if missing_states:
            raise ValueError(
                "Estimating ordinal thresholds requires every state to be observed; "
                "missing: " + ", ".join(missing_states)
            )
    initial_thresholds = (
        np.asarray(fixed_thresholds, dtype=float)
        if fixed_thresholds is not None
        else _initial_thresholds(observed_categories_array, len(states))
    )
    estimate_thresholds = fixed_thresholds is None and len(states) > 2
    process = build_threshold_process(tree)
    constrained_indices = np.asarray(tuple(constraints), dtype=int)
    allowed_by_constraint = tuple(constraints[index] for index in constraints)
    seed_sequences = np.random.SeedSequence(seed).spawn(chains)
    probability_counts = np.zeros((len(compiled.nodes), len(states)), dtype=float)
    liability_sum = np.zeros(len(compiled.nodes), dtype=float)
    liability_square_sum = np.zeros(len(compiled.nodes), dtype=float)
    threshold_sum = np.zeros(len(states) - 1, dtype=float)
    diagnostic_traces = np.empty(
        (chains, num_samples, 1 + max(0, len(states) - 2)), dtype=float
    )
    total_sweeps = burnin + num_samples * thin
    for chain_index, child_seed in enumerate(seed_sequences):
        rng = np.random.default_rng(child_seed)
        current_thresholds = initial_thresholds.copy()
        values = _initialize_values(compiled, constraints, current_thresholds)
        retained = 0
        for sweep in range(total_sweeps):
            for node_index in rng.permutation(len(compiled.nodes)):
                mean, sd = _conditional_parameters(
                    int(node_index), values, compiled, process
                )
                allowed = constraints.get(int(node_index))
                values[node_index] = (
                    rng.normal(mean, sd)
                    if allowed is None
                    else _truncated_draw(mean, sd, allowed, current_thresholds, rng)
                )
            if estimate_thresholds:
                _update_thresholds(
                    values[constrained_indices],
                    allowed_by_constraint,
                    current_thresholds,
                    rng,
                )
            if sweep < burnin or (sweep - burnin) % thin:
                continue
            categories = np.searchsorted(current_thresholds, values, side="left")
            probability_counts[np.arange(len(values)), categories] += 1.0
            liability_sum += values
            liability_square_sum += values * values
            threshold_sum += current_thresholds
            diagnostic_traces[chain_index, retained, 0] = values[0]
            if len(current_thresholds) > 1:
                diagnostic_traces[chain_index, retained, 1:] = current_thresholds[1:]
            retained += 1
    total_draws = chains * num_samples
    posterior = {
        node: probability_counts[index] / total_draws
        for index, node in enumerate(compiled.nodes)
    }
    means = liability_sum / total_draws
    variances = np.maximum(0.0, liability_square_sum / total_draws - means * means)
    liabilities = {
        node: GaussianMarginal(float(means[index]), float(variances[index]))
        for index, node in enumerate(compiled.nodes)
    }
    fitted_thresholds = tuple(float(value) for value in threshold_sum / total_draws)
    rhat, ess = _rhat_and_ess(diagnostic_traces)
    available_rhat = rhat[~np.isnan(rhat)]
    rhat_max = float(np.max(available_rhat)) if len(available_rhat) else math.nan
    ess_min = float(np.min(ess))
    statuses = []
    if math.isnan(rhat_max):
        statuses.append("mcmc_rhat_unavailable")
    elif rhat_max > 1.05:
        statuses.append("mcmc_rhat")
    if ess_min < max(100.0, 0.1 * total_draws):
        statuses.append("mcmc_low_ess")
    fit = ThresholdFit(
        thresholds=fitted_thresholds,
        thresholds_estimated=estimate_thresholds,
        num_samples=num_samples,
        burnin=burnin,
        thin=thin,
        chains=chains,
        seed=seed,
        rhat_max=rhat_max,
        ess_min=ess_min,
        fit_status="+".join(statuses) if statuses else "ok",
        liability_marginals=liabilities,
    )
    return posterior, fit


def threshold_liability_table(tree, fit):
    """Return per-node posterior liability moments."""

    ids = assign_branch_ids(tree)
    rows = []
    for node in tree.traverse():
        marginal = fit.liability_marginals[node]
        rows.append(
            {
                "branch_id": ids[node],
                "parent": -1 if node.is_root else ids[node.up],
                "node_class": get_node_class(node),
                "name": "" if node.name in (None, "") else str(node.name),
                "liability_mean": marginal.mean,
                "liability_variance": marginal.variance,
                "liability_sd": math.sqrt(marginal.variance),
                "num_posterior_draws": fit.num_samples * fit.chains,
            }
        )
    return pd.DataFrame(rows)
