"""Linear duplication/loss and retained genome-multiplication count likelihoods.

Finite matrices retain, rather than renormalize, probability outside the state
space. Callers must check convergence as the maximum count is increased.
"""

from dataclasses import dataclass

import numpy as np
from scipy.special import gammaln, logsumexp
from scipy.stats import binom, gamma


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
    one = np.empty(max_count + 1)
    one[0] = a
    one[1:] = survival * geometric_success * b ** np.arange(max_count)
    result = np.zeros((max_count + 1, max_count + 1))
    result[0, 0] = 1.0
    for count in range(1, max_count + 1):
        result[count] = np.convolve(result[count - 1], one)[: max_count + 1]
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
    with np.errstate(divide="ignore"):
        if duplication == loss:
            log_product = np.log(duplication) + np.log(time)
            log_denominator = np.logaddexp(0.0, log_product)
            log_a = log_b = log_product - log_denominator
            log_survival = log_success = -log_denominator
        else:
            delta = abs(duplication - loss)
            log_change = np.log(-np.expm1(-delta * time))
            log_denominator = np.logaddexp(
                np.log(delta), np.log(min(duplication, loss)) + log_change
            )
            log_a = np.log(loss) + log_change - log_denominator
            log_b = np.log(duplication) + log_change - log_denominator
            complement = np.log(delta) - log_denominator
            attenuated = complement - delta * time
            log_survival, log_success = (
                (complement, attenuated)
                if duplication > loss
                else (attenuated, complement)
            )
    descendants = np.arange(1, max_count + 1)[None, :]
    for count in range(1, max_count + 1):
        surviving = np.arange(1, count + 1)[:, None]
        # Sum over the number of surviving ancestral lineages. Conditional
        # descendants follow a negative binomial, not a killed finite CTMC.
        with np.errstate(invalid="ignore"):
            extinct = np.where(count == surviving, 0.0, (count - surviving) * log_a)
            extra = np.where(
                descendants == surviving, 0.0, (descendants - surviving) * log_b
            )
            terms = (
                gammaln(count + 1)
                - gammaln(surviving + 1)
                - gammaln(count - surviving + 1)
                + extinct
                + surviving * (log_survival + log_success)
                + gammaln(descendants)
                - gammaln(surviving)
                - gammaln(descendants - surviving + 1)
                + extra
            )
        terms[descendants < surviving] = -np.inf
        result[count, 1:] = logsumexp(terms, axis=0)
        result[count, 0] = count * log_a
    return result


def _log_matrix_product(left, right):
    return np.array([logsumexp(row[:, None] + right, axis=0) for row in left])


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
    result = np.zeros((max_count + 1, max_count + 1))
    for count in range(max_count + 1):
        extra = np.arange(max_count - count + 1)
        result[count, count:] = binom.pmf(extra, count * (multiplicity - 1), retention)
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

        # Each ancestral copy observes an independent subset of root clades.
        # Geometric continuation gives a triangular coverage recurrence; remove
        # its empty-subset self-loop with p + continuation * P(any observation).
        log_absent = np.log1p(-probabilities)
        log_continue = np.log(continuation)
        covered = {0: np.zeros(len(probabilities))}
        for missing in range(1, 1 << probabilities.shape[1]):
            columns = [i for i in range(probabilities.shape[1]) if missing & (1 << i)]
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
            union = _union_survival(probabilities[:, columns])
            covered[missing] = logsumexp(terms, axis=0) - np.log(
                p + continuation * union
            )
        return covered[(1 << probabilities.shape[1]) - 1]


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
        observed = self.counts[~np.isnan(self.counts)]
        if (
            np.any(~np.isfinite(observed))
            or np.any(observed < 0)
            or np.any(observed != np.floor(observed))
        ):
            raise ValueError("Observed counts must be finite nonnegative integers.")
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
        self.root_clade_columns = []
        for child in self.children[0]:
            columns = []
            for column, node in enumerate(tree.tip_nodes):
                ancestor = node
                while tree.parents[ancestor] != 0:
                    ancestor = tree.parents[ancestor]
                if ancestor == child:
                    columns.append(column)
            self.root_clade_columns.append(columns)
        if ascertainment == "root-clades" and any(
            np.any(np.nansum(self.counts[:, columns], axis=1) == 0)
            for columns in self.root_clade_columns
        ):
            raise ValueError(
                "Root-clade ascertainment requires observed copies in every root-child clade, for every family."
            )
        self.patterns, self.pattern_index = np.unique(
            np.isnan(self.counts), axis=0, return_inverse=True
        )
        self.subtree_missing = {
            node: np.isnan(self.counts[:, column])
            for column, node in enumerate(tree.tip_nodes)
        }
        for node in range(len(tree.parents) - 1, -1, -1):
            if node not in self.subtree_missing:
                self.subtree_missing[node] = np.all(
                    [self.subtree_missing[child] for child in self.children[node]],
                    axis=0,
                )

    def transitions(
        self, rates, scale, max_count, event=None, *, _split_zero_event=False
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
        matrices = {}
        cache = {}
        for node, time in enumerate(self.tree.lengths[1:], 1):
            duplication, loss = rates[self.branch_groups[node]] * scale
            key = (float(duplication), float(loss), time)
            if key not in cache:
                cache[key] = birth_death_transition(duplication, loss, time, max_count)
            if (
                event is None
                or event.node != node
                or (event.retention == 0 and not _split_zero_event)
            ):
                matrices[node] = cache[key]
            else:
                before = birth_death_transition(
                    duplication, loss, time * event.fraction, max_count
                )
                after = birth_death_transition(
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
        likelihoods: dict[int, np.ndarray] = {}
        scales: dict[int, np.ndarray] = {}
        rows = next(iter(emissions.values())).shape[0]
        for node in range(len(self.tree.parents) - 1, -1, -1):
            if node in emissions:
                value = emissions[node].copy()
                scale = np.zeros(rows)
            else:
                value = np.ones((rows, len(prior)))
                scale = np.zeros(rows)
                for child in self.children[node]:
                    contribution = likelihoods.pop(child) @ matrices[child].T
                    child_scale = scales.pop(child)
                    # An unobserved subtree integrates to one over the unbounded
                    # process, not to the truncated transition's row sum.
                    missing = self.subtree_missing[child]
                    contribution[missing] = 1.0
                    child_scale[missing] = 0.0
                    value *= contribution
                    scale += child_scale
                    normalizer = value.max(axis=1)
                    positive = normalizer > 0
                    value[positive] /= normalizer[positive, None]
                    with np.errstate(divide="ignore"):
                        scale += np.log(normalizer)
            likelihoods[node], scales[node] = value, scale
        with np.errstate(divide="ignore"):
            return np.log(likelihoods[0] @ prior) + scales[0]

    def _log_transitions(
        self, rates, scale, max_count, event, *, _split_zero_event=False
    ):
        matrices, cache = {}, {}
        for node, time in enumerate(self.tree.lengths[1:], 1):
            duplication, loss = np.asarray(rates)[self.branch_groups[node]] * scale
            key = (float(duplication), float(loss), time)
            if key not in cache:
                cache[key] = _birth_death_log_transition(
                    duplication, loss, time, max_count
                )
            if (
                event is None
                or event.node != node
                or (event.retention == 0 and not _split_zero_event)
            ):
                matrices[node] = cache[key]
                continue
            before = _birth_death_log_transition(
                duplication, loss, time * event.fraction, max_count
            )
            after = _birth_death_log_transition(
                duplication, loss, time * (1 - event.fraction), max_count
            )
            pulse = np.full_like(before, -np.inf)
            for count in range(max_count + 1):
                pulse[count, count:] = binom.logpmf(
                    np.arange(max_count - count + 1),
                    count * (event.multiplicity - 1),
                    event.retention,
                )
            matrices[node] = _log_matrix_product(
                _log_matrix_product(before, pulse), after
            )
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
                contribution = np.array(
                    [
                        logsumexp(matrices[child] + row[None, :], axis=1)
                        for row in likelihoods.pop(child)
                    ]
                )
                contribution[self.subtree_missing[child][rows]] = 0.0
                value += contribution
            likelihoods[node] = value
        return logsumexp(likelihoods[0] + prior[None, :], axis=1)

    def _selection_logs(self, rates, scale, root_mean, event):
        survivals = {
            node: np.where(self.patterns[:, column], 0.0, self.detection[column])
            for column, node in enumerate(self.tree.tip_nodes)
        }

        def propagate(survival, duplication, loss, time):
            _, b, survives, success = _birth_death_probabilities(
                duplication, loss, time
            )
            return np.divide(
                survives * survival,
                success + b * survival,
                out=np.zeros_like(survival),
                where=survival > 0,
            )

        rates = np.asarray(rates, dtype=float)
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
        if self.ascertainment == "root-clades":
            return _root_clade_log_survival(root_children, root_mean)
        survival = _union_survival(root_children)
        with np.errstate(divide="ignore"):
            return np.log(survival) - np.log(
                1 / root_mean + (1 - 1 / root_mean) * survival
            )

    def family_log_likelihoods(
        self, rates, root_mean, max_count, event=None, *, _split_zero_event=False
    ):
        if max_count < np.nanmax(self.counts):
            raise ValueError("Maximum state count is below an observed count.")
        prior = root_probabilities(root_mean, max_count)
        states = np.arange(max_count + 1)
        emissions = {}
        for column, node in enumerate(self.tree.tip_nodes):
            values = self.counts[:, column]
            missing = np.isnan(values)
            matrix = binom.pmf(
                np.nan_to_num(values)[:, None], states[None, :], self.detection[column]
            )
            matrix[missing] = 1.0
            emissions[node] = matrix
        observed_logs, selection_logs = [], []
        for scale in self.rate_scales:
            transitions = self.transitions(
                rates, scale, max_count, event, _split_zero_event=_split_zero_event
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
