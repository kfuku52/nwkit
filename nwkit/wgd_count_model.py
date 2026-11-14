"""Linear duplication/loss and retained genome-multiplication count likelihoods.

Finite matrices retain, rather than renormalize, probability outside the state
space. Callers must check convergence as the maximum count is increased.
"""

from collections import OrderedDict
from dataclasses import dataclass

import numpy as np
from scipy.special import logsumexp
from scipy.stats import binom, gamma

_LOG_PRUNING_WORKSPACE_BYTES = 1024**2
_PULSE_BATCH_ENTRIES = 4096
_TRANSITION_CACHE_BYTES = 8 * 1024**2
_TRANSITION_CACHE_ENTRIES = 4096
_EMISSION_CACHE_BYTES = 16 * 1024**2
_SELECTION_PROBABILITY_FLOOR = np.finfo(float).tiny / np.finfo(float).eps
_GEOMETRIC_RECURRENCE_MIN_COUNT = 256


def _log_pruning_contribution(matrix, values):
    # Estimate four float-sized broadcast/logsumexp temporaries. Retain the
    # row loop when even a two-family tile would exceed this workspace.
    batch = max(1, min(64, _LOG_PRUNING_WORKSPACE_BYTES // (32 * matrix.size)))
    if len(values) <= 1 or batch == 1:
        return np.array([logsumexp(matrix + row[None, :], axis=1) for row in values])
    contribution = np.empty((len(values), matrix.shape[0]))
    for start in range(0, len(values), batch):
        end = min(start + batch, len(values))
        contribution[start:end] = logsumexp(
            matrix[None, :, :] + values[start:end, None, :], axis=2
        )
    return contribution


def _nonnegative(value: float, name: str) -> float:
    value = float(value)
    if not np.isfinite(value) or value < 0:
        raise ValueError(f"{name} must be finite and nonnegative.")
    return value


def birth_death_parameters(duplication: float, loss: float, time: float):
    """Return extinction probability and surviving-lineage geometric parameter."""
    a, b, _, _ = _birth_death_probabilities(duplication, loss, time)
    return a, b


def _birth_death_probabilities(duplication, loss, time):
    """Compute both probabilities and their small complements directly."""
    duplication = _nonnegative(duplication, "Duplication rate")
    loss = _nonnegative(loss, "Loss rate")
    time = _nonnegative(time, "Time")
    if time == 0 or duplication == loss == 0:
        return 0.0, 0.0, 1.0, 1.0
    if duplication == loss:
        product = duplication * time
        value = 1.0 if np.isinf(product) else product / (1.0 + product)
        complement = 1.0 / (1.0 + product)
        return value, value, complement, complement
    delta = abs(duplication - loss)
    change = -np.expm1(-delta * time)
    denominator = delta + min(duplication, loss) * change
    complement = delta / denominator
    attenuated = complement * np.exp(-delta * time)
    survival, geometric_success = (
        (complement, attenuated) if duplication > loss else (attenuated, complement)
    )
    return (
        loss * change / denominator,
        duplication * change / denominator,
        survival,
        geometric_success,
    )


def birth_death_transition(
    duplication: float, loss: float, time: float, max_count: int
) -> np.ndarray:
    """Exact unbounded branching-process probabilities within 0..max_count.

    One ancestral copy is extinct with probability a, and otherwise has a
    geometric number of descendants. Independent copies are convolved. Paths
    that temporarily exceed the count bound and return are included here.
    """
    if isinstance(max_count, bool) or not isinstance(max_count, int) or max_count < 1:
        raise ValueError("max_count must be a positive integer.")
    a, b, survival, geometric_success = _birth_death_probabilities(
        duplication, loss, time
    )
    coefficient = survival * geometric_success
    if max_count >= _GEOMETRIC_RECURRENCE_MIN_COUNT:
        return _geometric_transition(a, b, coefficient, max_count)
    one = np.empty(max_count + 1)
    one[0] = a
    one[1:] = coefficient * b ** np.arange(max_count)
    result = np.zeros((max_count + 1, max_count + 1))
    result[0, 0] = 1.0
    for count in range(1, max_count + 1):
        result[count] = np.convolve(result[count - 1], one)[: max_count + 1]
    return result


def _geometric_transition(extinction, geometric_ratio, coefficient, max_count):
    from scipy.signal import lfilter

    # Convolution with the geometric tail obeys g[j] = c*p[j-1] + b*g[j-1].
    # Evaluate it in C for each ancestral count, without subtractive cancellation.
    result = np.zeros((max_count + 1, max_count + 1))
    result[0, 0] = 1.0
    numerator, denominator = (coefficient,), (1.0, -geometric_ratio)
    for count in range(1, max_count + 1):
        previous = result[count - 1]
        result[count] = extinction * previous
        result[count, 1:] += lfilter(numerator, denominator, previous[:-1])
    return result


def _birth_death_log_transition(duplication, loss, time, max_count):
    """The same unbounded kernel, including probabilities below float range."""
    duplication = _nonnegative(duplication, "Duplication rate")
    loss = _nonnegative(loss, "Loss rate")
    time = _nonnegative(time, "Time")
    result = np.full((max_count + 1, max_count + 1), -np.inf)
    result[0, 0] = 0.0
    if time == 0 or duplication == loss == 0:
        np.fill_diagonal(result, 0.0)
        return result
    log_a, log_b, log_survival, log_success = _birth_death_log_probabilities(
        duplication, loss, time
    )
    result[1:, 0] = np.arange(1, max_count + 1) * log_a
    log_one = log_survival + log_success
    positive = np.full(max_count, -np.inf)
    for descendants in range(1, max_count + 1):
        # For n ancestors and k descendants, the last ancestor produces
        # exactly one copy with probability d, or a positive geometric tail
        # with continuation b. Its positive contribution obeys
        # d*P[n-1,k-1] + b*positive[n,k-1]. Its extinction adds a*P[n-1,k].
        positive = np.logaddexp(
            log_one + result[:-1, descendants - 1], log_b + positive
        )
        column = positive.copy()
        # Double the number of ancestor terms accumulated at each step.
        # Terms are nonnegative, and impossible states stay at -inf.
        shift = 1
        while shift < max_count:
            column[shift:] = np.logaddexp(
                column[shift:], shift * log_a + column[:-shift]
            )
            shift *= 2
        result[1:, descendants] = column
    return result


def _birth_death_log_probabilities(duplication, loss, time):
    if time == 0 or duplication == loss == 0:
        return -np.inf, -np.inf, 0.0, 0.0
    if duplication == loss:
        log_product = np.log(duplication) + np.log(time)
        log_denominator = np.logaddexp(0.0, log_product)
        log_probability = log_product - log_denominator
        return log_probability, log_probability, -log_denominator, -log_denominator
    delta = abs(duplication - loss)
    log_delta = np.log(delta)
    product = delta * time
    with np.errstate(divide="ignore"):
        log_change = (
            log_delta + np.log(time)
            if product < np.finfo(float).tiny
            else np.log(-np.expm1(-product))
        )
        log_denominator = np.logaddexp(
            log_delta, np.log(min(duplication, loss)) + log_change
        )
        log_a = np.log(loss) + log_change - log_denominator
        log_b = np.log(duplication) + log_change - log_denominator
    complement = log_delta - log_denominator
    attenuated = complement - product
    if duplication > loss:
        return log_a, log_b, complement, attenuated
    return log_a, log_b, attenuated, complement


def _log_matrix_product(left, right):
    return np.array([logsumexp(row[:, None] + right, axis=0) for row in left])


def _log_pulse_transition(before, after, event):
    size = len(before)
    if event.retention in (0.0, 1.0):
        # The intermediate count is deterministic. Only include states whose
        # pulse output is within the bound, including for split zero retention.
        multiplier = 1 if event.retention == 0 else event.multiplicity
        limit = 1 + (size - 1) // multiplier
        if event.fraction == 0:
            result = np.full_like(after, -np.inf)
            result[:limit] = after[::multiplier]
            return result
        if event.fraction == 1:
            result = np.full_like(before, -np.inf)
            result[:, ::multiplier] = before[:, :limit]
            return result
        return _log_matrix_product(before[:, :limit], after[::multiplier])

    states = np.arange(size)
    trials = np.array([count * (event.multiplicity - 1) for count in range(size)])
    pulse = np.empty_like(before)
    batch = max(1, _PULSE_BATCH_ENTRIES // size)
    for first in range(0, size, batch):
        last = min(first + batch, size)
        extra = states[None, :] - states[first:last, None]
        pulse[first:last] = binom.logpmf(
            extra, trials[first:last, None], event.retention
        )
    if event.fraction == 0:
        return _log_matrix_product(pulse, after)
    if event.fraction == 1:
        return _log_matrix_product(before, pulse)
    return _log_matrix_product(_log_matrix_product(before, pulse), after)


def multiplication_transition(
    retention: float, multiplicity: int, max_count: int
) -> np.ndarray:
    """Keep original copies; retain each of (multiplicity-1)*n extra copies."""
    retention = float(retention)
    if not np.isfinite(retention) or not 0 <= retention <= 1:
        raise ValueError("Retention must be finite and in [0, 1].")
    if isinstance(multiplicity, bool) or not isinstance(multiplicity, int):
        raise ValueError("Multiplicity must be an integer >= 2.")
    if multiplicity < 2:
        raise ValueError("Multiplicity must be an integer >= 2.")
    if isinstance(max_count, bool) or not isinstance(max_count, int) or max_count < 1:
        raise ValueError("max_count must be a positive integer.")
    states = np.arange(max_count + 1)
    trials = np.array([count * (multiplicity - 1) for count in range(len(states))])
    result = np.empty((len(states), len(states)))
    batch = max(1, _PULSE_BATCH_ENTRIES // len(states))
    for first in range(0, len(states), batch):
        last = min(first + batch, len(states))
        extra = states[None, :] - states[first:last, None]
        result[first:last] = binom.pmf(extra, trials[first:last, None], retention)
    return result


def root_probabilities(mean: float, max_count: int) -> np.ndarray:
    """Positive geometric root prior with its omitted tail left unnormalized."""
    mean = float(mean)
    if not np.isfinite(mean) or mean < 1:
        raise ValueError("Root mean must be finite and >= 1.")
    result = np.zeros(max_count + 1)
    result[1:] = (1.0 / mean) * (1.0 - 1.0 / mean) ** np.arange(max_count)
    return result


def _union_survival(probabilities):
    result = np.zeros(probabilities.shape[0])
    for probability in probabilities.T:
        result += (1.0 - result) * probability
    return result


def _log_complement(log_probability):
    with np.errstate(divide="ignore"):
        return np.where(
            log_probability < -np.log(2),
            np.log1p(-np.exp(log_probability)),
            np.log(-np.expm1(log_probability)),
        )


def _log_union_survival(logs):
    seen = np.full(logs.shape[0], -np.inf)
    absent = np.zeros(logs.shape[0])
    for log_probability in logs.T:
        seen = np.logaddexp(seen, absent + log_probability)
        absent += _log_complement(log_probability)
    return np.minimum(seen, 0.0)


def _root_clade_coverage(logs, log_absent, log_continue, log_normalizer):
    # Each ancestral copy observes an independent subset of root clades.
    # Geometric continuation gives a triangular coverage recurrence.
    covered = {0: np.zeros(len(logs))}
    for missing in range(1, 1 << logs.shape[1]):
        columns = [i for i in range(logs.shape[1]) if missing & (1 << i)]
        terms = [logs[:, columns].sum(axis=1)]
        observed = (missing - 1) & missing
        while observed:
            present = [i for i in columns if observed & (1 << i)]
            absent = [i for i in columns if not observed & (1 << i)]
            terms.append(
                log_continue
                + logs[:, present].sum(axis=1)
                + log_absent[:, absent].sum(axis=1)
                + covered[missing ^ observed]
            )
            observed = (observed - 1) & missing
        covered[missing] = logsumexp(terms, axis=0) - log_normalizer(columns)
    return covered[(1 << logs.shape[1]) - 1]


def _root_clade_log_survival(probabilities, mean):
    """Exact geometric-root selection, without alternating subtraction."""
    p = 1.0 / mean
    continuation = 1.0 - p
    with np.errstate(divide="ignore"):
        logs = np.log(probabilities)
        if mean == 1:
            return logs.sum(axis=1)
        if probabilities.shape[1] == 1:
            survival = probabilities[:, 0]
            return logs[:, 0] - np.log(p + continuation * survival)
        if probabilities.shape[1] == 2:
            a, b = probabilities.T
            union = _union_survival(probabilities)
            factor = p * (p + 2 * continuation) + continuation**2 * union
            return (
                logs.sum(axis=1)
                + np.log(factor)
                - np.log(p + continuation * a)
                - np.log(p + continuation * b)
                - np.log(p + continuation * union)
            )

        # Remove the empty-subset self-loop using the chance of stopping or
        # observing at least one of the clades still needed.
        def log_normalizer(columns):
            union = _union_survival(probabilities[:, columns])
            return np.log(p + continuation * union)

        return _root_clade_coverage(
            logs, np.log1p(-probabilities), np.log(continuation), log_normalizer
        )


def _log_pulse_survival(log_survival, event):
    # For one ancestral copy, observe the original lineage or at least one
    # retained extra. Work with the log hazard of observing any extra copy.
    extra = np.log(event.retention) + log_survival
    cutoff = np.log(np.finfo(float).eps)
    with np.errstate(divide="ignore", over="ignore"):
        log_hazard = np.log(event.multiplicity - 1) + np.where(
            extra < cutoff, extra, np.log(-_log_complement(extra))
        )
        log_extra_survival = np.where(
            log_hazard < cutoff, log_hazard, _log_complement(-np.exp(log_hazard))
        )
    return np.minimum(
        np.logaddexp(log_survival, _log_complement(log_survival) + log_extra_survival),
        0.0,
    )


def _selection_log_from_survival(logs, mean, ascertainment):
    if ascertainment == "observed":
        logs = _log_union_survival(logs)[:, None]
    if mean == 1:
        return logs.sum(axis=1)
    log_stop = -np.log(mean)
    log_continue = np.log1p(-1 / mean)

    def log_normalizer(columns):
        union = _log_union_survival(logs[:, columns])
        return np.logaddexp(log_stop, log_continue + union)

    if logs.shape[1] == 1:
        return logs[:, 0] - log_normalizer([0])
    return _root_clade_coverage(
        logs, _log_complement(logs), log_continue, log_normalizer
    )


def rate_categories(shape: float | None, categories: int = 4) -> np.ndarray:
    """Equal-weight gamma quantile categories, normalized to mean one.

    This is a discrete approximation, not analytic gamma integration. None
    selects homogeneous families; shape and category count are recorded by CLI.
    """
    if shape is None:
        return np.ones(1)
    if not np.isfinite(shape) or shape <= 0:
        raise ValueError("Gamma shape must be finite and positive.")
    if (
        isinstance(categories, bool)
        or not isinstance(categories, int)
        or categories < 2
    ):
        raise ValueError("Gamma approximation needs at least two categories.")
    rates = gamma.ppf(
        (np.arange(categories) + 0.5) / categories, shape, scale=1 / shape
    )
    if not np.all(np.isfinite(rates)) or np.any(rates <= 0):
        raise ValueError("Gamma categories are numerically unrepresentable.")
    return rates / rates.mean()


@dataclass(frozen=True)
class CountTree:
    """Preorder tree arrays; root index zero; tip columns have explicit order."""

    parents: tuple[int, ...]
    lengths: tuple[float, ...]
    tip_nodes: tuple[int, ...]
    tip_names: tuple[str, ...]
    branch_ids: tuple[int, ...]
    clade_ids: tuple[str, ...]

    def __post_init__(self):
        n = len(self.parents)
        if n < 3 or self.parents[0] != -1:
            raise ValueError("Count tree needs a root and at least two tips.")
        if any(
            len(items) != n for items in (self.lengths, self.branch_ids, self.clade_ids)
        ):
            raise ValueError("Count tree node arrays must have equal lengths.")
        if any(
            not 0 <= parent < node for node, parent in enumerate(self.parents[1:], 1)
        ):
            raise ValueError("Count tree parents must precede their children.")
        if (
            any(not np.isfinite(t) or t < 0 for t in self.lengths)
            or self.lengths[0] != 0
        ):
            raise ValueError(
                "Count tree needs finite nonnegative lengths and root length zero."
            )
        children = {parent for parent in self.parents[1:]}
        tips = tuple(node for node in range(n) if node not in children)
        if (
            len(self.tip_nodes) != len(tips)
            or len(self.tip_names) != len(tips)
            or set(self.tip_nodes) != set(tips)
        ):
            raise ValueError("Tip columns must exactly and uniquely cover tree leaves.")
        if len(set(self.tip_nodes)) != len(tips) or len(set(self.tip_names)) != len(
            tips
        ):
            raise ValueError("Count tree tips must be unique.")
        if any(not name for name in self.tip_names):
            raise ValueError("Count tree tip names must be nonempty.")
        if len(set(self.branch_ids)) != n or len(set(self.clade_ids)) != n:
            raise ValueError("Count tree node identifiers must be unique.")

    @property
    def children(self):
        result: list[list[int]] = [[] for _ in self.parents]
        for node, parent in enumerate(self.parents[1:], 1):
            result[parent].append(node)
        return result


@dataclass(frozen=True)
class MultiplicationEvent:
    node: int
    retention: float
    fraction: float = 0.5
    multiplicity: int = 2

    def __post_init__(self):
        if (
            isinstance(self.node, bool)
            or not isinstance(self.node, int)
            or self.node < 1
        ):
            raise ValueError("Event node must be a positive integer branch index.")
        if not np.isfinite(self.fraction) or not 0 <= self.fraction <= 1:
            raise ValueError("Event fraction must be in [0, 1], measured from parent.")
        multiplication_transition(self.retention, self.multiplicity, 1)


class CountLikelihood:
    """Pruning likelihood conditional on >=1 observed copy per family.

    Missing tips contribute one, not zero copies. Binomial observation models
    use known detection probabilities. Families originate at the species root;
    later family origins and horizontal transfer are outside this model.
    """

    def __init__(
        self,
        tree: CountTree,
        counts: np.ndarray,
        *,
        detection: np.ndarray | None = None,
        rate_scales: np.ndarray | None = None,
        branch_groups: tuple[int, ...] | None = None,
        ascertainment: str = "observed",
    ):
        self.tree = tree
        self.counts = np.asarray(counts, dtype=float)
        if self.counts.ndim != 2 or self.counts.shape[1] != len(tree.tip_names):
            raise ValueError("Counts need one column per ordered species-tree tip.")
        missing = np.isnan(self.counts)
        observed = self.counts[~missing]
        if (
            np.any(~np.isfinite(observed))
            or np.any(observed < 0)
            or np.any(observed != np.floor(observed))
        ):
            raise ValueError("Observed counts must be finite nonnegative integers.")
        del observed
        if self.counts.shape[0] == 0 or np.any(np.nansum(self.counts, axis=1) == 0):
            raise ValueError("Every analyzed family needs at least one observed copy.")
        self.detection = (
            np.ones(len(tree.tip_names))
            if detection is None
            else np.asarray(detection, dtype=float)
        )
        if (
            self.detection.shape != (len(tree.tip_names),)
            or np.any(~np.isfinite(self.detection))
            or np.any(self.detection <= 0)
            or np.any(self.detection > 1)
        ):
            raise ValueError("Detection probabilities must be in (0, 1], one per tip.")
        self.rate_scales = (
            np.ones(1) if rate_scales is None else np.asarray(rate_scales, dtype=float)
        )
        if (
            self.rate_scales.ndim != 1
            or self.rate_scales.size == 0
            or np.any(~np.isfinite(self.rate_scales))
            or np.any(self.rate_scales <= 0)
        ):
            raise ValueError("Family rate scales must be a nonempty positive vector.")
        self.branch_groups = branch_groups or tuple(0 for _ in tree.parents)
        if len(self.branch_groups) != len(tree.parents) or any(
            not isinstance(group, int) or group < 0 for group in self.branch_groups
        ):
            raise ValueError(
                "Branch rate groups must be nonnegative integer node assignments."
            )
        active_groups = set(self.branch_groups[1:])
        if active_groups != set(range(max(active_groups) + 1)):
            raise ValueError(
                "Non-root branch rate groups must be contiguous from zero."
            )
        self.num_groups = len(active_groups)
        self.children = tree.children
        if ascertainment not in {"observed", "root-clades"}:
            raise ValueError("Ascertainment must be observed or root-clades.")
        self.ascertainment = ascertainment
        root_child = [0] * len(tree.parents)
        for node, parent in enumerate(tree.parents[1:], 1):
            root_child[node] = node if parent == 0 else root_child[parent]
        columns_by_clade: dict[int, list[int]] = {
            child: [] for child in self.children[0]
        }
        for column, node in enumerate(tree.tip_nodes):
            columns_by_clade[root_child[node]].append(column)
        self.root_clade_columns = list(columns_by_clade.values())
        if ascertainment == "root-clades" and any(
            np.any(np.nansum(self.counts[:, columns], axis=1) == 0)
            for columns in self.root_clade_columns
        ):
            raise ValueError(
                "Root-clade ascertainment requires observed copies in every root-child clade, for every family."
            )
        # Packing the first tip into the most significant bit preserves the
        # lexicographic order of missingness patterns and their inverse indices.
        packed = np.packbits(missing, axis=1)
        patterns, self.pattern_index = np.unique(packed, axis=0, return_inverse=True)
        self.patterns = np.unpackbits(patterns, axis=1, count=missing.shape[1]).view(
            bool
        )
        self.subtree_missing = {
            node: missing[:, column].copy()
            for column, node in enumerate(tree.tip_nodes)
        }
        for node in range(len(tree.parents) - 1, -1, -1):
            if node not in self.subtree_missing:
                self.subtree_missing[node] = np.all(
                    [self.subtree_missing[child] for child in self.children[node]],
                    axis=0,
                )
        self._transition_cache: OrderedDict[
            tuple[bool, float, float, float, int], np.ndarray
        ] = OrderedDict()
        self._transition_cache_bytes = 0
        self._emission_cache: (
            tuple[tuple[int, bytes, bytes], dict[int, np.ndarray]] | None
        ) = None

    def _cached_transition(self, duplication, loss, time, max_count, *, log=False):
        # Optimizer perturbations leave many branch groups unchanged. Cache
        # exact keys only, with bounds on both matrix payload and entry count.
        key = (log, float(duplication), float(loss), time, max_count)
        cached = self._transition_cache.get(key)
        if cached is not None:
            self._transition_cache.move_to_end(key)
            return cached
        kernel = _birth_death_log_transition if log else birth_death_transition
        matrix = kernel(duplication, loss, time, max_count)
        if matrix.nbytes <= _TRANSITION_CACHE_BYTES and _TRANSITION_CACHE_ENTRIES > 0:
            while (
                self._transition_cache_bytes + matrix.nbytes > _TRANSITION_CACHE_BYTES
                or len(self._transition_cache) >= _TRANSITION_CACHE_ENTRIES
            ):
                _, discarded = self._transition_cache.popitem(last=False)
                self._transition_cache_bytes -= discarded.nbytes
            self._transition_cache[key] = matrix
            self._transition_cache_bytes += matrix.nbytes
        return matrix

    def transitions(
        self,
        rates,
        scale,
        max_count,
        event=None,
        *,
        _split_zero_event=False,
        _kernel=None,
    ):
        rates = np.asarray(rates, dtype=float)
        if (
            rates.shape != (self.num_groups, 2)
            or np.any(~np.isfinite(rates))
            or np.any(rates < 0)
        ):
            raise ValueError(
                "Each branch group needs finite nonnegative duplication/loss rates."
            )
        if event is not None and (not 1 <= event.node < len(self.tree.parents)):
            raise ValueError("Multiplication event must identify a non-root branch.")
        kernel = birth_death_transition if _kernel is None else _kernel
        matrices = {}
        cache = {}
        for node, time in enumerate(self.tree.lengths[1:], 1):
            duplication, loss = rates[self.branch_groups[node]] * scale
            key = (float(duplication), float(loss), time)
            if (
                event is None
                or event.node != node
                or (event.retention == 0 and not _split_zero_event)
            ):
                if key not in cache:
                    cache[key] = kernel(duplication, loss, time, max_count)
                matrices[node] = cache[key]
            else:
                before = kernel(duplication, loss, time * event.fraction, max_count)
                after = kernel(
                    duplication, loss, time * (1 - event.fraction), max_count
                )
                matrices[node] = (
                    before
                    @ multiplication_transition(
                        event.retention, event.multiplicity, max_count
                    )
                    @ after
                )
        return matrices

    def _prune(self, emissions, matrices, prior):
        size = len(self.tree.parents)
        likelihoods = [None] * size
        scales = [None] * size
        rows = len(self.counts)
        with np.errstate(divide="ignore"):
            for node in range(size - 1, -1, -1):
                if node in emissions:
                    likelihoods[node] = emissions[node]
                    continue
                value = np.ones((rows, len(prior)))
                scale = np.zeros(rows)
                for child in self.children[node]:
                    if (tip := emissions.get(child)) is not None and tip.ndim == 1:
                        contribution = matrices[child].T[tip]
                    else:
                        contribution = likelihoods[child] @ matrices[child].T
                    # An unobserved subtree integrates to one over the unbounded
                    # process, not to the truncated transition's row sum.
                    missing = self.subtree_missing[child]
                    contribution[missing] = 1.0
                    value *= contribution
                    if (child_scale := scales[child]) is not None:
                        child_scale[missing] = 0.0
                        scale += child_scale
                    likelihoods[child] = scales[child] = None
                    normalizer = value.max(axis=1)
                    positive = normalizer > 0
                    value[positive] /= normalizer[positive, None]
                    scale += np.log(normalizer)
                likelihoods[node], scales[node] = value, scale
            return np.log(likelihoods[0] @ prior) + scales[0]

    def _log_transitions(
        self, rates, scale, max_count, event, *, _split_zero_event=False
    ):
        matrices, cache = {}, {}
        for node, time in enumerate(self.tree.lengths[1:], 1):
            duplication, loss = np.asarray(rates)[self.branch_groups[node]] * scale
            key = (float(duplication), float(loss), time)
            if (
                event is None
                or event.node != node
                or (event.retention == 0 and not _split_zero_event)
            ):
                if key not in cache:
                    cache[key] = self._cached_transition(
                        duplication, loss, time, max_count, log=True
                    )
                matrices[node] = cache[key]
                continue
            before = self._cached_transition(
                duplication, loss, time * event.fraction, max_count, log=True
            )
            after = self._cached_transition(
                duplication, loss, time * (1 - event.fraction), max_count, log=True
            )
            matrices[node] = _log_pulse_transition(before, after, event)
        return matrices

    def _prune_log(self, matrices, root_mean, max_count, rows):
        states = np.arange(max_count + 1)
        prior = np.full(max_count + 1, -np.inf)
        if root_mean == 1:
            prior[1] = 0.0
        else:
            prior[1:] = -np.log(root_mean) + np.arange(max_count) * np.log1p(
                -1 / root_mean
            )
        likelihoods = {}
        for column, node in enumerate(self.tree.tip_nodes):
            counts = self.counts[rows, column]
            value = binom.logpmf(
                np.nan_to_num(counts)[:, None], states[None, :], self.detection[column]
            )
            value[np.isnan(counts)] = 0.0
            likelihoods[node] = value
        for node in range(len(self.tree.parents) - 1, -1, -1):
            if node in likelihoods:
                continue
            value = np.zeros((len(rows), max_count + 1))
            for child in self.children[node]:
                contribution = _log_pruning_contribution(
                    matrices[child], likelihoods.pop(child)
                )
                contribution[self.subtree_missing[child][rows]] = 0.0
                value += contribution
            likelihoods[node] = value
        return logsumexp(likelihoods[0] + prior[None, :], axis=1)

    def _selection_logs(self, rates, scale, root_mean, event):
        rates = np.asarray(rates, dtype=float)
        if len(self.patterns) == 1:
            return self._single_pattern_selection(rates, scale, root_mean, event)
        survivals = {
            node: np.where(self.patterns[:, column], 0.0, self.detection[column])
            for column, node in enumerate(self.tree.tip_nodes)
        }
        rare = np.zeros(len(self.patterns), dtype=bool)

        def propagate(survival, duplication, loss, time):
            _, _, survives, success = _birth_death_probabilities(
                duplication, loss, time
            )
            # 1-b and b are computed separately to retain small complements;
            # their rounded sum can fall below one. Rewrite 1-b+b*s as
            # (1-b)*(1-s)+s, so s=1 stays exact and propagation cannot turn a
            # valid survival into a probability above one before a pulse.
            numerator = survives * survival
            rare[:] |= (survival > 0) & (numerator < _SELECTION_PROBABILITY_FLOOR)
            return np.divide(
                numerator,
                success * (1.0 - survival) + survival,
                out=np.zeros_like(survival),
                where=survival > 0,
            )

        for node in range(len(self.tree.parents) - 1, 0, -1):
            if node not in survivals:
                survivals[node] = _union_survival(
                    np.column_stack([survivals[child] for child in self.children[node]])
                )
            duplication, loss = rates[self.branch_groups[node]] * scale
            time = self.tree.lengths[node]
            survival = survivals[node]
            if event is not None and event.node == node and event.retention != 0:
                survival = propagate(
                    survival, duplication, loss, time * (1 - event.fraction)
                )
                with np.errstate(divide="ignore"):
                    extra = -np.expm1(
                        (event.multiplicity - 1) * np.log1p(-event.retention * survival)
                    )
                survival += (1 - survival) * extra
                time *= event.fraction
            survivals[node] = propagate(survival, duplication, loss, time)
        root_children = np.column_stack(
            [survivals[child] for child in self.children[0]]
        )

        if not np.any(rare):
            return self._condition_selection(root_children, root_mean)
        result = np.empty(len(rare))
        result[rare] = self._rare_selection_logs(rates, scale, root_mean, event, rare)
        if not np.all(rare):
            result[~rare] = self._condition_selection(root_children[~rare], root_mean)
        return result

    def _condition_selection(self, probabilities, root_mean):
        if self.ascertainment == "root-clades":
            return _root_clade_log_survival(probabilities, root_mean)
        survival = _union_survival(probabilities)
        with np.errstate(divide="ignore"):
            return np.log(survival) - np.log(
                1 / root_mean + (1 - 1 / root_mean) * survival
            )

    def _single_pattern_selection(self, rates, scale, root_mean, event):
        # A shared missingness mask needs only one survival probability per node.
        values = [0.0] * len(self.tree.parents)
        for column, node in enumerate(self.tree.tip_nodes):
            if not self.patterns[0, column]:
                values[node] = float(self.detection[column])
        rare = False

        def propagate(survival, duplication, loss, time):
            nonlocal rare
            if survival == 0:
                return 0.0
            _, _, survives, success = _birth_death_probabilities(
                duplication, loss, time
            )
            numerator = survives * survival
            rare |= numerator < _SELECTION_PROBABILITY_FLOOR
            return numerator / (success * (1.0 - survival) + survival)

        scaled_rates = rates * scale
        pulse = event.node if event is not None and event.retention != 0 else -1
        for node in range(len(values) - 1, 0, -1):
            survival = values[node]
            for child in self.children[node]:
                survival += (1.0 - survival) * values[child]
            duplication, loss = scaled_rates[self.branch_groups[node]]
            time = self.tree.lengths[node]
            if node == pulse:
                survival = propagate(
                    survival, duplication, loss, time * (1 - event.fraction)
                )
                with np.errstate(divide="ignore"):
                    extra = -np.expm1(
                        (event.multiplicity - 1) * np.log1p(-event.retention * survival)
                    )
                survival += (1 - survival) * extra
                time *= event.fraction
            values[node] = propagate(survival, duplication, loss, time)
        if rare:
            return self._rare_selection_logs(
                rates, scale, root_mean, event, np.array([True])
            )
        return self._condition_selection(
            np.array([[values[child] for child in self.children[0]]]), root_mean
        )

    def _rare_selection_logs(self, rates, scale, root_mean, event, rare):
        survivals = {
            node: np.where(self.patterns[rare, column], -np.inf, np.log(detection))
            for column, (node, detection) in enumerate(
                zip(self.tree.tip_nodes, self.detection, strict=True)
            )
        }

        def propagate(survival, duplication, loss, time):
            _, _, log_survives, log_success = _birth_death_log_probabilities(
                duplication, loss, time
            )
            result = np.full_like(survival, -np.inf)
            present = np.isfinite(survival)
            logs = survival[present]
            result[present] = log_survives - np.logaddexp(
                0, log_success + _log_complement(logs) - logs
            )
            return result

        for node in range(len(self.tree.parents) - 1, 0, -1):
            if node not in survivals:
                survivals[node] = _log_union_survival(
                    np.column_stack([survivals[child] for child in self.children[node]])
                )
            duplication, loss = rates[self.branch_groups[node]] * scale
            time = self.tree.lengths[node]
            survival = survivals[node]
            if event is not None and event.node == node and event.retention != 0:
                survival = propagate(
                    survival, duplication, loss, time * (1 - event.fraction)
                )
                survival = _log_pulse_survival(survival, event)
                time *= event.fraction
            survivals[node] = propagate(survival, duplication, loss, time)
        root_children = np.column_stack(
            [survivals[child] for child in self.children[0]]
        )
        return _selection_log_from_survival(
            root_children, root_mean, self.ascertainment
        )

    def _emissions(self, max_count):
        # Observations and the count bound determine all tip emissions. Keep
        # only the current bound and include mutable inputs in the cache key.
        if np.all(self.detection == 1):
            self._emission_cache = None
            return self._build_emissions(max_count)
        cache_key = None
        if (
            self.counts.nbytes * (max_count + 2) + self.detection.nbytes
            <= _EMISSION_CACHE_BYTES
        ):
            cache_key = (max_count, self.counts.tobytes(), self.detection.tobytes())
            if (
                self._emission_cache is not None
                and cache_key == self._emission_cache[0]
            ):
                return self._emission_cache[1]
        self._emission_cache = None
        emissions = self._build_emissions(max_count)
        if cache_key is not None:
            self._emission_cache = (cache_key, emissions)
        return emissions

    def _build_emissions(self, max_count):
        states = np.arange(max_count + 1)
        emissions = {}
        for column, node in enumerate(self.tree.tip_nodes):
            values = self.counts[:, column]
            missing = np.isnan(values)
            probability = self.detection[column]
            if probability == 1:
                # Exact tip counts select a transition column. The placeholder
                # index for missing counts is replaced by one during pruning.
                emissions[node] = np.where(missing, 0, values).astype(np.intp)
                continue
            # Families with the same observed count have identical emissions.
            counts, inverse = np.unique(np.nan_to_num(values), return_inverse=True)
            matrix = binom.pmf(counts[:, None], states[None, :], probability)[inverse]
            matrix[missing] = 1.0
            emissions[node] = matrix
        return emissions

    def family_log_likelihoods(
        self, rates, root_mean, max_count, event=None, *, _split_zero_event=False
    ):
        if max_count < np.nanmax(self.counts):
            raise ValueError("Maximum state count is below an observed count.")
        prior = root_probabilities(root_mean, max_count)
        emissions = self._emissions(max_count)
        observed_logs, selection_logs = [], []
        for scale in self.rate_scales:
            transitions = self.transitions(
                rates,
                scale,
                max_count,
                event,
                _split_zero_event=_split_zero_event,
                _kernel=self._cached_transition,
            )
            logs = self._prune(emissions, transitions, prior)
            # Scaling cannot recover a transition or latent-state contribution
            # that underflowed before normalization. Evaluate rare families in
            # the log semiring, with the same state bound and probability law.
            # Leave machine-precision headroom for accumulated underflow error
            # over every represented state and branch. This is an arithmetic
            # domain switch, not a likelihood or state-convergence tolerance.
            arithmetic_limit = (
                np.log(np.finfo(float).tiny)
                - np.log(np.finfo(float).eps)
                + np.log(len(self.tree.parents) * (max_count + 1))
            )
            rare = np.flatnonzero(logs < arithmetic_limit)
            if len(rare):
                log_transitions = self._log_transitions(
                    rates, scale, max_count, event, _split_zero_event=_split_zero_event
                )
                logs[rare] = self._prune_log(
                    log_transitions, root_mean, max_count, rare
                )
            observed_logs.append(logs)
            selection_logs.append(self._selection_logs(rates, scale, root_mean, event))
        observed_log = logsumexp(observed_logs, axis=0) - np.log(len(self.rate_scales))
        # Ascertainment is conditioned after mixing, not separately per category.
        survival_log = logsumexp(selection_logs, axis=0) - np.log(len(self.rate_scales))
        result = observed_log - survival_log[self.pattern_index]
        if np.any(np.isnan(result)) or np.any(result > 1e-8):
            raise ValueError(
                "Count likelihood is numerically undefined; increase the state bound."
            )
        return result

    def log_likelihood(self, rates, root_mean, max_count, event=None):
        return float(
            self.family_log_likelihoods(rates, root_mean, max_count, event).sum()
        )

    def convergence(self, rates, root_mean, max_count, event=None):
        """Per-family doubling error, including the q->0+ split-edge boundary."""
        # The positive-q limit has a truncated latent count at the pulse even
        # though q=0 uses the exact full-edge kernel. Evaluate the split kernel
        # at exactly zero, so no epsilon perturbation changes the probability law.
        checks = (
            (False, True) if event is not None and event.retention == 0 else (False,)
        )
        error = 0.0
        for split_zero in checks:
            lower = self.family_log_likelihoods(
                rates, root_mean, max_count, event, _split_zero_event=split_zero
            )
            upper = self.family_log_likelihoods(
                rates, root_mean, 2 * max_count, event, _split_zero_event=split_zero
            )
            if not np.all(np.isfinite(lower)) or not np.all(np.isfinite(upper)):
                return float("inf")
            error = max(error, float(np.max(np.abs(upper - lower))))
        return error

    def simulate(self, rates, root_mean, rng, event=None, *, max_attempts=10000):
        """Unbounded branching-process simulation, with original missingness masks."""
        rates = np.asarray(rates, dtype=float)
        self.transitions(rates, 1.0, 1, event)
        root_probabilities(root_mean, 1)
        result = np.full_like(self.counts, np.nan)
        tip_columns = {node: col for col, node in enumerate(self.tree.tip_nodes)}
        for row in range(len(result)):
            for _ in range(max_attempts):
                scale = float(rng.choice(self.rate_scales))
                latent = {0: int(rng.geometric(1 / root_mean))}
                for node, time in enumerate(self.tree.lengths[1:], 1):
                    duplication, loss = rates[self.branch_groups[node]] * scale
                    count = latent[self.tree.parents[node]]
                    if event is not None and event.node == node:
                        count = _simulate_edge(
                            count, duplication, loss, time * event.fraction, rng
                        )
                        count += int(
                            rng.binomial(
                                count * (event.multiplicity - 1), event.retention
                            )
                        )
                        count = _simulate_edge(
                            count, duplication, loss, time * (1 - event.fraction), rng
                        )
                    else:
                        count = _simulate_edge(count, duplication, loss, time, rng)
                    latent[node] = count
                    if node in tip_columns:
                        column = tip_columns[node]
                        if not np.isnan(self.counts[row, column]):
                            result[row, column] = rng.binomial(
                                count, self.detection[column]
                            )
                selected = np.nansum(result[row]) > 0
                if self.ascertainment == "root-clades":
                    selected = all(
                        np.nansum(result[row, columns]) > 0
                        for columns in self.root_clade_columns
                    )
                if selected:
                    break
            else:
                raise ValueError(
                    "Simulation cannot sample an ascertained family; check extinction/detection."
                )
        return result


def _simulate_edge(count, duplication, loss, time, rng):
    _, _, survival, success = _birth_death_probabilities(duplication, loss, time)
    surviving = int(rng.binomial(count, survival))
    if surviving == 0:
        return 0
    if success <= 0:
        raise ValueError("Simulated copy counts overflow; rescale rates/time.")
    return surviving + int(rng.negative_binomial(surviving, success))
