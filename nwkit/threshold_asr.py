"""Bayesian threshold/liability ancestral reconstruction.

The latent liability follows a standardized Brownian process with
``root ~ N(0, 1)`` and diffusion variance one.  This fixes the otherwise
unidentifiable location/scale.  Binary traits use threshold zero.  For ordinal
traits the first threshold is zero and remaining thresholds are sampled, or
all thresholds may be supplied explicitly.
"""

import math
from dataclasses import dataclass
from fractions import Fraction

import numpy as np
import pandas as pd
from scipy import special
from scipy.fft import irfft, next_fast_len, rfft
from scipy.special import log_ndtr, ndtri
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

_ndtri_exp = getattr(special, "ndtri_exp", None)
_MAX_FLOAT = float(np.finfo(float).max)


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


def _normal_interval_log_mass(lower, upper):
    # Reflect positive-tail intervals so the log CDF is evaluated in the tail.
    if lower > 0.0:
        lower, upper = -upper, -lower
    log_upper = float(log_ndtr(upper))
    log_lower = float(log_ndtr(lower))
    if log_lower < log_upper:
        return log_upper + math.log(-math.expm1(log_lower - log_upper))
    if lower < upper and math.isfinite(log_upper):
        # Near the center, distinct endpoints can round to the same CDF.
        midpoint = 0.5 * lower + 0.5 * upper
        return math.log(upper - lower) - 0.5 * (midpoint**2 + math.log(math.tau))
    return -math.inf


def _truncated_quantile(lower, upper, rng):
    if _ndtri_exp is None or not lower < upper:
        return float(truncnorm.rvs(lower, upper, random_state=rng))
    uniform = float(rng.random())
    if uniform == 0.0:
        return lower
    log_mass = _normal_interval_log_mass(lower, upper)
    if lower >= 0.0 or (upper > 0.0 and uniform > 0.5):
        log_survival = np.logaddexp(log_ndtr(-upper), math.log1p(-uniform) + log_mass)
        return -float(_ndtri_exp(min(0.0, log_survival)))
    log_cdf = np.logaddexp(log_ndtr(lower), math.log(uniform) + log_mass)
    return float(_ndtri_exp(min(0.0, log_cdf)))


def _narrow_interval_log_mass(lower, upper, a, b, sd):
    if not math.isfinite(a) or b - a > 1e-12 * (1.0 + abs(a)):
        return None
    width = upper - lower
    # In this regime the normal density changes by at most about one part in
    # 1e8 across the interval. Work in raw units when CDF inversion loses width.
    if width <= 0.0 or (width / sd) * (1.0 + abs(a)) > 1e-8:
        return None
    midpoint = 0.5 * a + 0.5 * b
    return (
        math.log(width)
        - math.log(sd)
        - (0.5 * midpoint) * midpoint
        - 0.5 * math.log(math.tau)
    )


def _interval_log_mass(lower, upper, a, b, sd):
    narrow_mass = _narrow_interval_log_mass(lower, upper, a, b, sd)
    return _normal_interval_log_mass(a, b) if narrow_mass is None else narrow_mass


def _relative_tail_log_masses(bounds, mean, sd):
    # The large-distance normal mass is proportional to exp(-z²/2)/z
    # times the integrated exponential over the finite interval. Subtract
    # the nearest interval's exponent before converting to float, keeping
    # even offsets that disappear in the standardized endpoints.
    center = Fraction(float(mean))
    variance = Fraction(float(sd)) ** 2
    distances = [
        abs(Fraction(lower if lower > mean else upper) - center)
        for lower, upper in bounds
    ]
    closest = min(distances)
    return np.asarray(
        [
            _relative_tail_log_mass(lower, upper, distance, closest, variance)
            for (lower, upper), distance in zip(bounds, distances, strict=True)
        ]
    )


def _relative_tail_log_mass(lower, upper, distance, closest, variance):
    try:
        penalty = float((distance - closest) * (distance + closest) / (2 * variance))
    except OverflowError:
        return -math.inf
    log_mass = -penalty + _log_fraction(closest / distance)
    if not math.isfinite(lower) or not math.isfinite(upper):
        return log_mass
    width = Fraction(upper) - Fraction(lower)
    rate_width = width * distance / variance
    if rate_width > 36:
        return log_mass
    if rate_width < 1e-15:
        return log_mass + _log_fraction(rate_width)
    return log_mass + math.log(-math.expm1(-float(rate_width)))


def _log_fraction(value):
    floating = float(value)
    if floating >= 2.0**-1022:
        return math.log(floating)
    return math.log(value.numerator) - math.log(value.denominator)


def _disjoint_interval_log_masses(bounds, intervals, mean, sd):
    nearest = min(max(a, -b) for a, b in intervals)
    if nearest >= 1e14:
        return _relative_tail_log_masses(bounds, mean, sd)
    return np.asarray(
        [
            _interval_log_mass(lower, upper, a, b, sd)
            for (lower, upper), (a, b) in zip(bounds, intervals, strict=True)
        ]
    )


def _interval_draw(lower, upper, a, b, mean, sd, rng):
    narrow_mass = _narrow_interval_log_mass(lower, upper, a, b, sd)
    if narrow_mass is not None:
        return lower + float(rng.random()) * (upper - lower)
    if a >= 1e5:
        return _tail_interval_draw(lower, upper, a, mean, sd, rng)
    if b <= -1e5:
        return -_tail_interval_draw(-upper, -lower, -b, -mean, sd, rng)
    if not a < b:
        raise ValueError("Could not resolve the allowed ordinal liability interval.")
    return mean + sd * _truncated_quantile(a, b, rng)


def _tail_interval_draw(lower, upper, standardized_lower, mean, sd, rng):
    # Shift the tail to the truncation boundary before drawing the excess.
    # An exponential proposal with rate a for standardized excess z accepts
    # with probability exp(-z**2/2), including for finite upper truncation.
    width = upper - lower
    reciprocal = 1.0 / standardized_lower
    if math.isinf(standardized_lower):
        distance = lower - mean
        reciprocal = (
            sd / distance
            if math.isfinite(distance)
            else ((0.5 * sd) / (0.5 * lower - 0.5 * mean))
        )
    scale = sd * reciprocal
    if scale == 0.0:
        return lower
    ratio = width / scale
    mass = -math.expm1(-ratio)
    while True:
        uniform = float(rng.random())
        proposal = -math.log1p(-uniform * mass)
        standardized_excess = proposal * reciprocal
        if float(rng.random()) < math.exp(-0.5 * standardized_excess**2):
            excess = width * uniform if ratio < 1e-12 else scale * proposal
            return lower + excess


def _truncated_draw(mean, sd, bounds, rng, *, thresholds=None):
    if not math.isfinite(mean):
        raise ValueError("An ordinal liability conditional mean is not finite.")
    if thresholds is not None:
        bounds = _allowed_intervals(bounds, thresholds)
    if len(bounds) == 1:
        lower, upper = bounds[0]
        a, b = (lower - mean) / sd, (upper - mean) / sd
    else:
        intervals = [
            ((lower - mean) / sd, (upper - mean) / sd) for lower, upper in bounds
        ]
        log_masses = _disjoint_interval_log_masses(bounds, intervals, mean, sd)
        largest = float(np.max(log_masses))
        if not math.isfinite(largest):
            raise ValueError(
                "Could not resolve the allowed ordinal liability intervals."
            )
        weights = np.exp(log_masses - largest)
        index = int(rng.choice(len(intervals), p=weights / weights.sum()))
        a, b = intervals[index]
        lower, upper = bounds[index]
    location = _interval_draw(lower, upper, a, b, mean, sd, rng)
    # Ordinal categories use (lower, upper]. A continuous tail draw can
    # round onto its excluded endpoint when projected back into trait units.
    first_allowed = math.nextafter(lower, math.inf)
    if not math.isfinite(location) or not math.isfinite(first_allowed):
        raise ValueError("A sampled ordinal liability exceeds floating-point range.")
    return min(upper, max(location, first_allowed))


def _liability_conditionals(compiled, process):
    """Cache unit-Brownian conditional weights and standard deviations."""

    assert process.root.variance is not None
    variances = (process.root.variance,) + tuple(
        process.transitions[node].variance for node in compiled.nodes[1:]
    )
    conditionals = []
    for index, parent in enumerate(compiled.parents):
        neighbors = (parent, *compiled.children[index])
        edge_variances = (
            variances[index],
            *(variances[child] for child in compiled.children[index]),
        )
        smallest = min(edge_variances)
        precisions = tuple(smallest / variance for variance in edge_variances)
        total = math.fsum(precisions)
        standard_deviation = math.sqrt(smallest) / math.sqrt(total)
        weights = tuple(weight / total for weight in precisions)
        # The identified Brownian root has mean zero and no parent value.
        if index == 0:
            neighbors = neighbors[1:]
            weights = weights[1:]
        conditionals.append(
            (tuple(zip(neighbors, weights, strict=True)), standard_deviation)
        )
    return conditionals


def _initialize_values(compiled, constraints, thresholds):
    values = np.zeros(len(compiled.nodes), dtype=float)
    for index, allowed in constraints.items():
        candidates = []
        for category in allowed:
            lower, upper = _category_interval(int(category), thresholds)
            if math.isfinite(lower) and math.isfinite(upper):
                candidate = 0.5 * lower + 0.5 * upper
            elif math.isfinite(lower):
                candidate = lower + 1.0
            elif math.isfinite(upper):
                candidate = upper - 1.0
            else:
                candidate = 0.0
            candidates.append(
                min(upper, max(candidate, math.nextafter(lower, math.inf)))
            )
        values[index] = min(candidates, key=abs)
    return values


def _update_thresholds(values, allowed_categories, thresholds, rng):
    if len(thresholds) <= 1:
        return
    for threshold_index in range(1, len(thresholds)):
        previous = float(thresholds[threshold_index - 1])
        following = (
            float(thresholds[threshold_index + 1])
            if threshold_index + 1 < len(thresholds)
            else math.inf
        )
        # Categories are (lower, upper]. Preserve every representable
        # admissible cut, including a singleton between adjacent floats.
        relevant = (values > previous) & (values <= following)
        allowed_lower = allowed_categories[threshold_index]
        allowed_upper = allowed_categories[threshold_index + 1]
        if np.any(relevant & ~allowed_lower & ~allowed_upper):
            raise ValueError(
                "Ordinal threshold constraints are internally inconsistent."
            )
        lower = float(
            np.max(
                values,
                where=relevant & allowed_lower & ~allowed_upper,
                initial=previous,
            )
        )
        upper = float(
            np.min(
                values,
                where=relevant & allowed_upper & ~allowed_lower,
                initial=following,
            )
        )
        if math.isfinite(upper):
            upper = math.nextafter(upper, -math.inf)
        lower = max(lower, math.nextafter(previous, math.inf))
        if lower > upper:
            raise ValueError(
                "Ordinal threshold update collapsed; increase MCMC burn-in or fix --thresholds."
            )
        thresholds[threshold_index] = rng.uniform(lower, upper)


def _effective_sample_size(centered, pooled):
    chains, samples = centered.shape
    total_draws = chains * samples
    if pooled == 0.0 or samples < 3:
        return float(total_draws)

    fft_size = next_fast_len(2 * samples - 1, real=True)
    frequency = rfft(centered, n=fft_size, axis=1)
    covariance = (
        irfft(np.mean(frequency.real**2 + frequency.imag**2, axis=0), n=fft_size)[
            :samples
        ]
        / samples
    )
    correlations = 1.0 - (covariance[0] - covariance) / pooled
    correlations[0] = 1.0
    paired = correlations[: 2 * (samples // 2)].reshape(-1, 2).sum(axis=1)
    # Geyer's initial positive, monotone sequence of paired correlations.
    nonpositive = np.flatnonzero(paired <= 0.0)
    if len(nonpositive):
        paired = paired[: nonpositive[0]]
    autocorrelation_time = max(
        1.0, 2.0 * float(np.minimum.accumulate(paired).sum()) - 1.0
    )
    return total_draws / autocorrelation_time


def _rhat_and_ess(traces):
    traces = np.asarray(traces, dtype=float)
    chains, samples, dimensions = traces.shape
    rhat = np.full(dimensions, np.nan, dtype=float)
    ess = np.full(dimensions, chains * samples, dtype=float)
    for dimension in range(dimensions):
        values = traces[:, :, dimension]
        if samples < 2:
            continue
        _, exponent = math.frexp(float(np.max(np.abs(values))))
        scaled = np.ldexp(values, -exponent)
        means = np.mean(scaled, axis=1, keepdims=True)
        centered = scaled - means
        within = float(np.sum(centered**2)) / (chains * (samples - 1))
        between = float(np.var(means, ddof=1)) if chains > 1 else 0.0
        pooled = (1.0 - 1.0 / samples) * within + between
        if chains > 1:
            rhat[dimension] = (
                math.sqrt(pooled / within)
                if within > 0.0
                else 1.0
                if between == 0.0
                else math.inf
            )
        ess[dimension] = _effective_sample_size(centered, pooled)
    return rhat, ess


def _accumulate_threshold_mean(mean, thresholds, draw_count, estimated):
    if estimated:
        mean += (thresholds - mean) / draw_count


def _allowed_categories_by_constraint(constraints, states):
    allowed = np.zeros((states, len(constraints)), dtype=bool)
    for index, categories in enumerate(constraints.values()):
        allowed[categories, index] = True
    return allowed


def _liability_bounds(constraints, thresholds, estimated):
    if estimated:
        return constraints
    return {
        index: _allowed_intervals(allowed, thresholds)
        for index, allowed in constraints.items()
    }


def _liability_moment_units(compiled, process, thresholds, retained_draws):
    # Ordinary trees can accumulate squared deviations in their input units.
    # For extreme branch lengths or thresholds, allow each node's units to
    # grow with the magnitude of its actual draws.
    limit = math.frexp(math.sqrt(_MAX_FLOAT / (256 * retained_draws)))[1] - 1
    branch_variance = max(
        (edge.variance for edge in process.transitions.values()), default=1.0
    )
    size = max(
        math.sqrt(branch_variance) * math.sqrt(len(compiled.nodes)),
        float(np.max(np.abs(thresholds))),
    )
    exponents = (
        None
        if size < math.ldexp(1.0, limit)
        else np.zeros(len(compiled.nodes), dtype=int)
    )
    return limit, exponents


def _scaled_liability_values(values, mean, m2, limit, exponents):
    if exponents is None:
        return values
    _, value_exponents = np.frexp(values)
    shift = np.maximum(value_exponents - limit - exponents, 0)
    if np.any(shift):
        np.ldexp(mean, -shift, out=mean)
        np.ldexp(m2, -2 * shift, out=m2)
        exponents += shift
    return np.ldexp(values, -exponents)


def _liability_moments(mean, m2, draws, exponents):
    variance = np.maximum(0.0, m2 / draws)
    with np.errstate(over="ignore"):
        if exponents is not None:
            mean = np.ldexp(mean, exponents)
            variance = np.ldexp(variance, 2 * exponents)
    if not np.isfinite(mean).all() or not np.isfinite(variance).all():
        raise ValueError("An ordinal liability summary exceeds floating-point range.")
    return mean, variance


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
    conditionals = _liability_conditionals(compiled, process)
    bounds = _liability_bounds(constraints, initial_thresholds, estimate_thresholds)
    constrained_indices = np.asarray(tuple(constraints), dtype=int)
    allowed_by_constraint = _allowed_categories_by_constraint(constraints, len(states))
    seed_sequences = np.random.SeedSequence(seed).spawn(chains)
    probability_counts = np.zeros((len(compiled.nodes), len(states)), dtype=float)
    liability_mean = np.zeros(len(compiled.nodes), dtype=float)
    liability_m2 = np.zeros(len(compiled.nodes), dtype=float)
    total_draws = chains * num_samples
    moment_limit, moment_exponents = _liability_moment_units(
        compiled, process, initial_thresholds, total_draws
    )
    threshold_mean = initial_thresholds.copy()
    diagnostic_traces = np.empty(
        (chains, num_samples, 1 + max(0, len(states) - 2)), dtype=float
    )
    total_sweeps = burnin + num_samples * thin
    retained_draws = 0
    for chain_index, child_seed in enumerate(seed_sequences):
        rng = np.random.default_rng(child_seed)
        current_thresholds = initial_thresholds.copy()
        moving_thresholds = current_thresholds if estimate_thresholds else None
        values = _initialize_values(compiled, constraints, current_thresholds)
        retained = 0
        for sweep in range(total_sweeps):
            for node_index in rng.permutation(len(compiled.nodes)):
                neighbors, sd = conditionals[node_index]
                mean = math.fsum(
                    float(values[neighbor]) * weight for neighbor, weight in neighbors
                )
                intervals = bounds.get(int(node_index))
                values[node_index] = (
                    rng.normal(mean, sd)
                    if intervals is None
                    else _truncated_draw(
                        mean, sd, intervals, rng, thresholds=moving_thresholds
                    )
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
            retained_draws += 1
            moment_values = _scaled_liability_values(
                values, liability_mean, liability_m2, moment_limit, moment_exponents
            )
            centered = moment_values - liability_mean
            liability_mean += centered / retained_draws
            liability_m2 += centered * (moment_values - liability_mean)
            _accumulate_threshold_mean(
                threshold_mean, current_thresholds, retained_draws, estimate_thresholds
            )
            diagnostic_traces[chain_index, retained, 0] = values[0]
            if len(current_thresholds) > 1:
                diagnostic_traces[chain_index, retained, 1:] = current_thresholds[1:]
            retained += 1
    posterior = {
        node: probability_counts[index] / total_draws
        for index, node in enumerate(compiled.nodes)
    }
    liability_mean, variances = _liability_moments(
        liability_mean, liability_m2, total_draws, moment_exponents
    )
    liabilities = {
        node: GaussianMarginal(float(liability_mean[index]), float(variances[index]))
        for index, node in enumerate(compiled.nodes)
    }
    fitted_thresholds = tuple(float(value) for value in threshold_mean)
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
