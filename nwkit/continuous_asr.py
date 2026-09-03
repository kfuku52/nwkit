"""Linear-time Brownian ancestral reconstruction with a flat root prior.

The upward pass integrates latent states; the downward pass smooths every
state conditional on all observations. Zero-length edges identify states,
not nearly deterministic transitions. Rates are in squared trait units per
input branch-length unit; no tip-variance normalization is applied.
"""

import math
from dataclasses import dataclass, replace
from typing import Any

import numpy as np
from scipy.optimize import minimize_scalar

from nwkit.continuous_asr_optimize import minimize_brownian_rate
from nwkit.gaussian_tree import (
    GaussianRootPrior,
    GaussianTreeProcess,
    brownian_transition,
)
from nwkit.rooting_state import require_rooted
from nwkit.util import validate_unique_named_leaves

_LOG_2PI = math.log(2.0 * math.pi)
_LOG_2 = math.log(2.0)


@dataclass(frozen=True, slots=True)
class GaussianMarginal:
    mean: float
    variance: float


@dataclass(frozen=True)
class BrownianFit:
    sigma2: float
    sigma2_estimated: bool
    restricted_log_likelihood: float | None
    num_observed: int
    num_effective_observations: int
    residual_df: int
    fit_status: str


@dataclass(frozen=True)
class _BrownianData:
    node_groups: dict[Any, int]
    parents: tuple[int, ...]
    lengths: tuple[float, ...]
    observations: tuple[GaussianMarginal | None, ...]
    exact_values: dict[int, float]
    local_log_variance: float
    local_quadratic: float
    local_df: int
    num_observed: int
    num_effective: int
    trait_center: float
    trait_exponent: int
    time_exponent: int


@dataclass(frozen=True)
class _PruningResult:
    inside: tuple[GaussianMarginal | None, ...]
    log_variance: float
    quadratic: float
    residual_df: int

    @property
    def log_likelihood(self):
        return -0.5 * (self.residual_df * _LOG_2PI + self.log_variance + self.quadratic)


class _ExactConstraintError(ValueError):
    """The data violate an exactly deterministic state constraint."""


def _finite_number(value, label, *, nonnegative=False):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{label} must be numeric, not boolean.")
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{label} must be numeric and finite.") from exc
    if not math.isfinite(number) or (nonnegative and number < 0.0):
        qualifier = "non-negative and finite" if nonnegative else "finite"
        raise ValueError(f"{label} must be {qualifier}.")
    return number


def _ldexp(value, exponent, label):
    try:
        result = math.ldexp(value, exponent)
    except OverflowError as exc:
        raise ValueError(
            f"{label} exceeds floating-point range; rescale the input units."
        ) from exc
    if not math.isfinite(result):
        raise ValueError(
            f"{label} exceeds floating-point range; rescale the input units."
        )
    if value != 0.0 and result == 0.0:
        raise ValueError(
            f"{label} underflows floating-point range; rescale the input units."
        )
    return result


def _nonnegative_sum(values):
    try:
        return math.fsum(values)
    except OverflowError:
        return math.inf


def _merge(first, second_mean, second_variance):
    """Multiply Gaussian messages without allocating an incoming wrapper."""
    if first is None:
        return GaussianMarginal(second_mean, second_variance), 0.0, 0.0, 0
    larger = max(first.variance, second_variance)
    if larger == 0.0:
        if first.mean != second_mean:
            raise _ExactConstraintError(
                "Conflicting exact observations have zero likelihood across "
                "zero-length branches or under sigma2=0."
            )
        return first, 0.0, 0.0, 0
    ratio = min(first.variance, second_variance) / larger
    denominator = 1.0 + ratio
    first_weight = (second_variance / larger) / denominator
    second_weight = (first.variance / larger) / denominator
    mean = first_weight * first.mean + second_weight * second_mean
    variance = min(first.variance, second_variance) / denominator
    if min(first.variance, second_variance) > 0.0 and variance == 0.0:
        raise ValueError(
            "A positive conditional variance underflowed; rescale the input units."
        )
    standardized = (first.mean - second_mean) / math.sqrt(larger)
    quadratic = (standardized * standardized) / denominator
    return (
        GaussianMarginal(mean, variance),
        math.log(larger) + math.log1p(ratio),
        quadratic,
        1,
    )


def _tree_groups(tree, *, validated=False, process=None):
    """Contract only exactly zero-length edges, without modifying the tree."""
    if not validated:
        validate_unique_named_leaves(tree, option_name="--infile", context=" for 'asr'")
        require_rooted(tree, "ASR requires a rooted tree.")
    raw_lengths = {}
    if process is None:
        for node in tree.traverse(strategy="preorder"):
            if node.is_root:
                continue
            if validated:
                length = float(node.dist)
            else:
                if node.dist is None:
                    raise ValueError(
                        "ASR requires branch lengths for all non-root nodes."
                    )
                length = _finite_number(node.dist, "Branch lengths", nonnegative=True)
            raw_lengths[node] = length
        process = GaussianTreeProcess(
            tree=tree,
            transitions={
                node: brownian_transition(length)
                for node, length in raw_lengths.items()
            },
            root=GaussianRootPrior("flat", 0.0, None),
            model="brownian",
        )
    if process.tree is not tree:
        raise ValueError("The Gaussian ASR process must use the input tree.")
    if process.root.mode != "flat":
        raise ValueError("Brownian ASR requires a flat-root Gaussian process.")
    node_groups = {tree: 0}
    parents, lengths = [-1], [0.0]
    for node in tree.traverse(strategy="preorder"):
        if node.is_root:
            continue
        transition = process.transitions[node]
        if transition.slope != 1.0 or transition.intercept != 0.0:
            raise ValueError(
                "Brownian ASR requires unit-slope zero-intercept branches."
            )
        length = transition.variance
        parent = node_groups[node.up]
        if length == 0.0:
            node_groups[node] = parent
            continue
        index = len(parents)
        node_groups[node] = index
        parents.append(parent)
        lengths.append(length)
    return node_groups, parents, lengths


def _observations_by_group(tree, node_groups, values_by_leaf, standard_errors):
    grouped: dict[int, list[tuple[float, float]]] = {}
    exact_by_group: dict[int, float] = {}
    num_observed = 0
    for leaf in tree.leaves():
        raw_value = values_by_leaf.get(str(leaf.name))
        if raw_value is None:
            continue
        value = _finite_number(raw_value, f"Trait value for '{leaf.name}'")
        if standard_errors is not None and str(leaf.name) not in standard_errors:
            raise ValueError(
                f"A standard error is required for observed tip '{leaf.name}'."
            )
        error = _finite_number(
            0.0 if standard_errors is None else standard_errors[str(leaf.name)],
            f"Standard error for '{leaf.name}'",
            nonnegative=True,
        )
        group = node_groups[leaf]
        num_observed += 1
        if error == 0.0 and group in exact_by_group:
            if value != exact_by_group[group]:
                raise _ExactConstraintError(
                    "Conflicting exact observations are connected by zero-length branches."
                )
            continue
        if error == 0.0:
            exact_by_group[group] = value
        grouped.setdefault(group, []).append((value, error))
    if not grouped:
        raise ValueError("Continuous ASR requires at least one observed tip value.")
    return grouped, num_observed, exact_by_group


def _scales(grouped, lengths, sigma2, trait_exponent_ceiling=None):
    records = [record for observations in grouped.values() for record in observations]
    values = [value for value, _ in records]
    # Center at a precise observation, not halfway towards an arbitrarily
    # uncertain outlier: subtraction must not erase differences between exact
    # values before they constrain either the rate or the root.
    center = min(records, key=lambda record: (record[1], abs(record[0])))[0]
    if any(not math.isfinite(value - center) for value in values):
        center = min(values) / 2.0 + max(values) / 2.0
    spread = max(abs(value - center) for value in values)
    error = max(error for observations in grouped.values() for _, error in observations)
    maximum_length = max(lengths)
    time_exponent = math.frexp(maximum_length)[1] - 1 if maximum_length > 0.0 else 0
    if time_exponent > 0:
        minimum_length = min(length for length in lengths if length > 0.0)
        # Scaling the longest edge towards one must not erase a representable
        # subnormal short edge before its compensating rate scale is applied.
        while time_exponent > 0 and math.ldexp(minimum_length, -time_exponent) == 0.0:
            time_exponent -= 1
    size = max(spread, error)
    trait_exponent = math.frexp(size)[1] - 1 if size > 0.0 else 0
    if sigma2 is not None and sigma2 > 0.0:
        process_exponent = (math.frexp(sigma2)[1] + time_exponent + 1) // 2
        trait_exponent = max(trait_exponent, process_exponent)
    if trait_exponent_ceiling is not None:
        trait_exponent = min(trait_exponent, trait_exponent_ceiling)
    return center, trait_exponent, time_exponent


def _prepare_data(
    tree,
    values_by_leaf,
    standard_errors,
    sigma2,
    *,
    trait_exponent_ceiling=None,
    tree_validated=False,
    process=None,
):
    node_groups, parents, lengths = _tree_groups(
        tree, validated=tree_validated, process=process
    )
    grouped, count, exact_values = _observations_by_group(
        tree, node_groups, values_by_leaf, standard_errors
    )
    center, trait_exponent, time_exponent = _scales(
        grouped, lengths, sigma2, trait_exponent_ceiling
    )
    observations: list[GaussianMarginal | None] = [None] * len(parents)
    log_variances, quadratics, local_df = [], [], 0
    for group, records in grouped.items():
        for value, error in records:
            mean = _ldexp(value - center, -trait_exponent, "Trait differences")
            scaled_error = _ldexp(error, -trait_exponent, "Standard errors")
            variance = scaled_error * scaled_error
            if scaled_error > 0.0 and variance == 0.0:
                raise ValueError(
                    "A positive measurement variance underflowed; rescale the input units."
                )
            message, log_variance, quadratic, df = _merge(
                observations[group], mean, variance
            )
            observations[group] = message
            log_variances.append(log_variance)
            quadratics.append(quadratic)
            local_df += df
    return _BrownianData(
        node_groups=node_groups,
        parents=tuple(parents),
        lengths=tuple(
            _ldexp(length, -time_exponent, "Branch lengths") for length in lengths
        ),
        observations=tuple(observations),
        exact_values=exact_values,
        local_log_variance=math.fsum(log_variances),
        local_quadratic=_nonnegative_sum(quadratics),
        local_df=local_df,
        num_observed=count,
        num_effective=sum(len(records) for records in grouped.values()),
        trait_center=center,
        trait_exponent=trait_exponent,
        time_exponent=time_exponent,
    )


def _edge_variance(rate, length, *, allow_underflow=False):
    variance = rate * length
    if not math.isfinite(variance):
        raise ValueError(
            "An evolutionary variance exceeds floating-point range; rescale the input units."
        )
    if rate > 0.0 and length > 0.0 and variance == 0.0 and not allow_underflow:
        raise ValueError(
            "An evolutionary variance underflowed; rescale the input units."
        )
    return variance


def _prune(data, rate, *, allow_edge_underflow=False, collect_inside=True):
    if rate == 0.0 and len(set(data.exact_values.values())) > 1:
        raise _ExactConstraintError(
            "Conflicting exact observations have zero likelihood under sigma2=0."
        )
    parents = data.parents
    lengths = data.lengths
    inside = list(data.observations)
    log_variances = [data.local_log_variance]
    quadratics = [data.local_quadratic]
    residual_df = data.local_df
    for index in range(len(parents) - 1, 0, -1):
        message = inside[index]
        if message is None:
            continue
        length = lengths[index]
        edge_variance = rate * length
        if not math.isfinite(edge_variance):
            raise ValueError(
                "An evolutionary variance exceeds floating-point range; "
                "rescale the input units."
            )
        if (
            edge_variance == 0.0
            and rate > 0.0
            and length > 0.0
            and not allow_edge_underflow
        ):
            raise ValueError(
                "An evolutionary variance underflowed; rescale the input units."
            )
        variance = message.variance + edge_variance
        if not math.isfinite(variance):
            raise ValueError("An ancestral variance exceeds floating-point range.")
        parent = parents[index]
        merged, log_variance, quadratic, df = _merge(
            inside[parent], message.mean, variance
        )
        inside[parent] = merged
        log_variances.append(log_variance)
        quadratics.append(quadratic)
        residual_df += df
    return _PruningResult(
        tuple(inside) if collect_inside else (),
        math.fsum(log_variances),
        _nonnegative_sum(quadratics),
        residual_df,
    )


def _smooth(data, rate, inside):
    root = inside[0]
    if root is None:
        raise ValueError("The flat root prior is improper without observed tip values.")
    marginals = [root]
    for index in range(1, len(data.parents)):
        parent = marginals[data.parents[index]]
        child = inside[index]
        innovation = _edge_variance(rate, data.lengths[index])
        if innovation == 0.0:
            marginal = parent
        elif child is None:
            marginal = GaussianMarginal(parent.mean, parent.variance + innovation)
        else:
            conditional, _, _, _ = _merge(child, parent.mean, innovation)
            assert conditional is not None
            scale = max(child.variance, innovation)
            gain = (child.variance / scale) / (
                child.variance / scale + innovation / scale
            )
            marginal = GaussianMarginal(
                conditional.mean, conditional.variance + gain * gain * parent.variance
            )
        if not math.isfinite(marginal.variance):
            raise ValueError("An ancestral variance exceeds floating-point range.")
        marginals.append(marginal)
    return marginals


def _zero_rate_result(data):
    try:
        return _prune(data, 0.0, collect_inside=False)
    except _ExactConstraintError:
        return None


def _rate_upper_bound(data):
    # Once each zero-edge group has one weighted observation, within-group
    # residuals are rate-independent. For its between-group contrasts z,
    # Q = z' C^-1 z bounds the REML optimum above, even with known sampling
    # covariance S: zz' <= Q C and V = rate*C + S imply a decreasing
    # likelihood for rate > Q. No arbitrary rate limits are needed.
    exact = replace(
        data,
        observations=tuple(
            None if item is None else GaussianMarginal(item.mean, 0.0)
            for item in data.observations
        ),
        local_log_variance=0.0,
        local_quadratic=0.0,
        local_df=0,
    )
    return _prune(exact, 1.0, collect_inside=False).quadratic


def _exact_rate_lower_bound(data, residual_df):
    exact = replace(
        data,
        observations=tuple(
            item if index in data.exact_values else None
            for index, item in enumerate(data.observations)
        ),
    )
    quadratic = _prune(exact, 1.0, collect_inside=False).quadratic
    lower = quadratic / residual_df
    if lower == 0.0 or math.ulp(lower) > 1e-12 * lower:
        raise ValueError(
            "Exact-observation contrasts exceed floating-point range; "
            "the input dynamic range cannot resolve Brownian sigma2."
        )
    # In whitened contrasts, exact observations supply K/r to Q(r).
    # f'(r) <= (d/r - K/r^2)/2, so no optimum lies below K/d.
    return lower


def _noise_rate_shift(data, upper):
    observed = [i for i, item in enumerate(data.observations) if item is not None]
    anchor = min(observed, key=lambda i: data.observations[i].variance)
    distances = [0.0] * len(data.parents)
    ancestors = {anchor}
    index = anchor
    while data.parents[index] >= 0:
        parent = data.parents[index]
        distances[parent] = distances[index] + data.lengths[index]
        ancestors.add(parent)
        index = parent
    for index in range(1, len(distances)):
        if index not in ancestors:
            distances[index] = distances[data.parents[index]] + data.lengths[index]
    # Contrasts to the most precise position have S >= min(other variances) I.
    # lambda_max(C) <= trace(C), whose diagonal entries are path distances.
    # The factor two leaves room for roundoff in distance accumulation.
    trace_bound = 2.0 * math.fsum(distances[i] for i in observed if i != anchor)
    minimum_noise = min(data.observations[i].variance for i in observed if i != anchor)
    shift = min(upper, minimum_noise / trace_bound)
    if shift == 0.0:
        raise ValueError("Brownian measurement variances exceed floating-point range.")
    return shift


def _optimize_noisy_rate(data, upper, zero):
    residual_df = data.num_effective - 1
    if zero is not None and zero.residual_df < residual_df:
        return 0.0  # Singular, feasible zero-variance limit; density is unbounded.
    lower = _exact_rate_lower_bound(data, residual_df) if zero is None else 0.0
    shift = _noise_rate_shift(data, upper) if zero is not None else 0.0

    def evaluate(rate):
        # A positive exploratory rate can be too small to contribute a
        # representable variance on the shortest edge even when the eventual
        # optimum is well resolved. At that probe only, the correctly rounded
        # contribution is zero. The final fit is pruned again without this
        # allowance, so an unrepresentable reported model still fails.
        result = _prune(data, rate, allow_edge_underflow=True, collect_inside=False)
        return result.log_variance, result.quadratic

    return minimize_brownian_rate(
        evaluate, lower, upper, shift, residual_df, minimize_scalar
    )


def _fit_rate(data):
    observed_groups = sum(item is not None for item in data.observations)
    if observed_groups < 2:
        raise ValueError(
            "Brownian sigma2 is not identifiable from observations at only one "
            "distinct tree position; specify --sigma2."
        )
    upper = _rate_upper_bound(data)
    if not math.isfinite(upper):
        raise ValueError(
            "Brownian rate estimation exceeds floating-point range; rescale the input units."
        )
    if upper == 0.0:
        return 0.0
    if all(item is None or item.variance == 0.0 for item in data.observations):
        return upper / float(observed_groups - 1)
    # Within-position measurement residuals are independent of the rate. Their
    # constants can dwarf the entire between-position objective; exclude them
    # before optimization, not by subtracting already-rounded likelihoods.
    fitting_data = replace(
        data,
        local_log_variance=0.0,
        local_quadratic=0.0,
        local_df=0,
        num_effective=observed_groups,
    )
    return _optimize_noisy_rate(fitting_data, upper, _zero_rate_result(fitting_data))


def compute_bm_marginals(
    tree,
    values_by_leaf,
    *,
    sigma2=None,
    standard_errors=None,
    _tree_validated=False,
    _process=None,
):
    """Return all-node Gaussian marginals and a fixed-rate or REML BM fit.

    Missing values are absent keys or None. Standard errors, when provided,
    must cover every observed tip. The root prior is flat with respect to the
    input trait units. Reported variances condition on the fitted sigma2;
    they do include uncertainty in the root value. A singular zero-rate
    boundary has no ordinary residual density, reported as None, not a
    regular finite likelihood on a different-dimensional support.
    """
    if sigma2 is not None:
        sigma2 = _finite_number(sigma2, "--sigma2", nonnegative=True)
    trait_exponent_ceiling = None
    # The private fast path is used only by the CLI after its shared ASR tree
    # validation; direct library calls always validate here.
    tree_validated = _tree_validated
    while True:
        data = _prepare_data(
            tree,
            values_by_leaf,
            standard_errors,
            sigma2,
            trait_exponent_ceiling=trait_exponent_ceiling,
            tree_validated=tree_validated,
            process=_process,
        )
        tree_validated = True
        rate = (
            _fit_rate(data)
            if sigma2 is None
            else _ldexp(
                sigma2,
                data.time_exponent - 2 * data.trait_exponent,
                "--sigma2",
            )
        )
        fitted_sigma2 = (
            _ldexp(
                rate,
                2 * data.trait_exponent - data.time_exponent,
                "Fitted sigma2",
            )
            if sigma2 is None
            else sigma2
        )
        scaled_underflow = rate > 0.0 and any(
            length > 0.0 and rate * length == 0.0 for length in data.lengths
        )
        if not scaled_underflow:
            break
        # Retry only when power-of-two trait normalization, rather than the
        # requested/fitted model in its original units, erased a variance.
        # The retry changes no mathematical units and strictly lowers the
        # exponent, so genuinely unrepresentable final models still reach the
        # explicit error in `_prune` below.
        original_underflow = any(
            node.dist > 0.0 and fitted_sigma2 * node.dist == 0.0
            for node in tree.traverse()
            if not node.is_root
        )
        if original_underflow:
            break
        candidate_rate = rate
        candidate_exponent = data.trait_exponent
        while any(
            length > 0.0 and candidate_rate * length == 0.0 for length in data.lengths
        ):
            try:
                candidate_rate = math.ldexp(candidate_rate, 2)
            except OverflowError:
                break
            if not math.isfinite(candidate_rate):
                break
            candidate_exponent -= 1
        if (
            not math.isfinite(candidate_rate)
            or candidate_exponent >= data.trait_exponent
        ):
            break
        trait_exponent_ceiling = candidate_exponent
    upward = _prune(data, rate)
    smoothed = _smooth(data, rate, upward.inside)
    restored = [
        GaussianMarginal(
            data.trait_center
            + _ldexp(item.mean, data.trait_exponent, "Ancestral mean"),
            _ldexp(item.variance, 2 * data.trait_exponent, "Ancestral variance"),
        )
        for item in smoothed
    ]
    if any(not math.isfinite(item.mean) for item in restored):
        raise ValueError("An ancestral mean exceeds floating-point range.")
    # Exact observations retain the original floating-point value, including
    # after centering/scaling; all members of a zero-edge group share it.
    posterior = {
        node: GaussianMarginal(data.exact_values[group], 0.0)
        if group in data.exact_values
        else restored[group]
        for node, group in data.node_groups.items()
    }
    residual_df = data.num_effective - 1
    singular = upward.residual_df < residual_df
    log_likelihood = (
        None
        if singular
        else upward.log_likelihood - residual_df * data.trait_exponent * _LOG_2
    )
    if log_likelihood is not None and not math.isfinite(log_likelihood):
        raise ValueError(
            "Brownian residual log-likelihood is not finite; check the rate and input scales."
        )
    fit = BrownianFit(
        sigma2=fitted_sigma2,
        sigma2_estimated=sigma2 is None,
        restricted_log_likelihood=log_likelihood,
        num_observed=data.num_observed,
        num_effective_observations=data.num_effective,
        residual_df=residual_df,
        fit_status="singular_zero_boundary"
        if singular
        else "zero_boundary"
        if rate == 0.0
        else "ok",
    )
    return posterior, fit
