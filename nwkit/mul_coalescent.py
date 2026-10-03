"""Fixed-topology MSC probabilities, integrating all compatible coalescent histories."""

import math
from collections import defaultdict
from itertools import product
from typing import Any

import numpy as np
from scipy.linalg import expm
from scipy.special import logsumexp

from nwkit.mul_reconcile_model import validate_binary


def _rate(lineages):
    return lineages * (lineages - 1) / 2


def _small_time_transition(start, end, time):
    """Positive uniformization series; bound its omitted tail relative to the sum."""
    rate = _rate(start)
    log_mu = math.log(rate) + math.log(time)
    mu = rate * time
    states = np.full(start - end + 1, -math.inf)
    states[-1] = 0.0
    total = -math.inf
    log_poisson = -mu
    for number in range(1000):
        total = float(np.logaddexp(total, log_poisson + states[0]))
        next_poisson = log_poisson + log_mu - math.log(number + 1)
        tail = next_poisson - math.log1p(-mu / (number + 2))
        if math.isfinite(total) and tail < total + math.log(1e-15):
            return total
        advanced = np.full_like(states, -math.inf)
        for index, lineages in enumerate(range(end, start + 1)):
            holding = 1 - _rate(lineages) / rate
            if holding > 0:
                advanced[index] = states[index] + math.log(holding)
            if index + 1 < len(states):
                advanced[index] = np.logaddexp(
                    advanced[index],
                    states[index + 1] + math.log(_rate(lineages + 1) / rate),
                )
        states, log_poisson = advanced, next_poisson
    raise ArithmeticError("MSC uniformization did not meet its numerical tail bound.")


def _long_time_transition(start, end, time):
    """Use the slowest hypoexponential term only with a relative remainder bound."""
    rates = [_rate(k) for k in range(end, start + 1)]
    denominator = sum(math.log(rate - rates[0]) for rate in rates[1:])
    leading = sum(math.log(rate) for rate in rates[1:]) - denominator
    remainders = [
        denominator
        - sum(math.log(abs(other - rate)) for j, other in enumerate(rates) if i != j)
        - (rate - rates[0]) * time
        for i, rate in enumerate(rates)
        if i
    ]
    if float(logsumexp(remainders)) < math.log(1e-15):
        result = leading - rates[0] * time
        if not math.isfinite(result):
            raise ArithmeticError(
                "MSC log probability exceeds the finite numerical range."
            )
        return min(result, 0.0)
    return None


def log_lineage_transition(start, end, time):
    """Log P(k lineages become j) for a Kingman pure-death process in t units."""
    if not 1 <= end <= start or not math.isfinite(time) or time < 0:
        raise ValueError("MSC transition requires 1 <= j <= k and finite time >= 0.")
    if time == 0:
        return 0.0 if start == end else -math.inf
    if start == end:
        result = -_rate(start) * time
        if not math.isfinite(result):
            raise ArithmeticError(
                "MSC log probability exceeds the finite numerical range."
            )
        return result
    if _rate(start) * time <= 0.5:
        return _small_time_transition(start, end, time)
    if time > 40:
        asymptote = _long_time_transition(start, end, time)
        if asymptote is not None:
            return asymptote
    # Shift by the slowest exit rate before expm: retain e.g. log P=-10000.
    size = start - end + 1
    generator = np.zeros((size, size))
    for index, lineages in enumerate(range(end, start + 1)):
        generator[index, index] = _rate(end) - _rate(lineages)
        if index:
            generator[index, index - 1] = _rate(lineages)
    transition = float(expm(generator * time)[-1, 0])
    result = math.log(transition) - _rate(end) * time if transition > 0 else math.nan
    if not math.isfinite(result) or result > 1e-10:
        raise ArithmeticError("MSC lineage transition failed numerical validation.")
    return min(result, 0.0)


class CoalescentTopology:
    """Ancestral configurations are antichains of nodes of one binary gene tree."""

    def __init__(self, gene, *, max_states=100000):
        validate_binary(gene, "Gene tree")
        if max_states < 1:
            raise ValueError("MSC state limit must be positive.")
        self.nodes = tuple(gene.traverse("postorder"))
        index = {node: number for number, node in enumerate(self.nodes)}
        self.cherries = tuple(
            (index[node], tuple(index[child] for child in node.children))
            for node in self.nodes
            if not node.is_leaf
        )
        self.leaves = {node.name: index[node] for node in self.nodes if node.is_leaf}
        self.root = index[gene]
        self.max_states = max_states
        self.work = 0
        self._descendants = {}
        self._transitions = {}
        self._branches = {}

    def _charge(self, count=1):
        self.work += count
        if self.work > self.max_states:
            raise ValueError("MSC calculation exceeds --max-coalescent-states.")

    def successors(self, configuration):
        present = set(configuration)
        for parent, children in self.cherries:
            if set(children) <= present:
                yield tuple(sorted((present - set(children)) | {parent}))

    def descendants(self, configuration):
        if configuration in self._descendants:
            return self._descendants[configuration]
        levels = {configuration: 1}
        descendants = {}
        while levels:
            following: dict[tuple[int, ...], int] = defaultdict(int)
            for state, count in levels.items():
                self._charge()
                descendants[state] = count
                for successor in self.successors(state):
                    self._charge()
                    following[successor] += count
            levels = following
        self._descendants[configuration] = descendants
        return descendants

    def branch(self, configuration, time):
        key = configuration, time
        if key in self._branches:
            return self._branches[key]
        if not configuration or time == 0:
            result = {configuration: 0.0}
        else:
            start = len(configuration)
            result = {}
            for state, paths in self.descendants(configuration).items():
                end = len(state)
                transition_key = start, end, time
                if transition_key not in self._transitions:
                    self._charge((start - end + 1) ** 2)
                    self._transitions[transition_key] = log_lineage_transition(
                        start, end, time
                    ) - sum(math.log(_rate(k)) for k in range(end + 1, start + 1))
                result[state] = math.log(paths) + self._transitions[transition_key]
        self._branches[key] = result
        return result

    def _join(self, left, right):
        result: dict[tuple[int, ...], float] = {}
        for (first, log_first), (second, log_second) in product(
            left.items(), right.items()
        ):
            self._charge()
            state = tuple(sorted(first + second))
            result[state] = float(
                np.logaddexp(result.get(state, -math.inf), log_first + log_second)
            )
        return result

    def _advance(self, distribution, time):
        result: dict[tuple[int, ...], float] = {}
        for initial, log_initial in distribution.items():
            for state, log_transition in self.branch(initial, time).items():
                self._charge()
                result[state] = float(
                    np.logaddexp(
                        result.get(state, -math.inf), log_initial + log_transition
                    )
                )
        return result

    def log_probability(self, species, assignment, *, time_scale=1.0):
        """Condition on one injective gene-tip to population-tip assignment."""
        validate_binary(species, "Population tree", unique=False)
        if not math.isfinite(time_scale) or time_scale <= 0:
            raise ValueError("MSC time scale must be finite and positive.")
        tips = {node for node in species.leaves()}
        if set(assignment) != set(self.leaves) or any(
            node not in tips for node in assignment.values()
        ):
            raise ValueError(
                "MSC assignment must map every gene tip to a population tip."
            )
        if len(set(assignment.values())) != len(assignment):
            raise ValueError("MSC prototype samples at most one copy per subgenome.")
        active: dict[Any, tuple[int, ...]] = {
            node: (self.leaves[name],) for name, node in assignment.items()
        }
        distributions: dict[Any, dict[tuple[int, ...], float]] = {}
        for node in species.traverse("postorder"):
            if node.is_leaf:
                distribution = {active.get(node, ()): 0.0}
            else:
                left, right = (distributions.pop(child) for child in node.children)
                distribution = self._join(left, right)
            if not node.is_root:
                if node.dist is None or not math.isfinite(node.dist) or node.dist < 0:
                    raise ValueError(
                        "MSC needs finite, nonnegative population branch lengths."
                    )
                duration = node.dist / time_scale
                if not math.isfinite(duration):
                    raise ValueError("MSC scaled branch time must be finite.")
                distribution = self._advance(distribution, duration)
            distributions[node] = distribution
        terms = []
        for state, log_mass in distributions[species].items():
            paths = self.descendants(state).get((self.root,), 0)
            if paths:
                terms.append(
                    log_mass
                    + math.log(paths)
                    - sum(math.log(_rate(k)) for k in range(2, len(state) + 1))
                )
        result = float(logsumexp(terms))
        if not math.isfinite(result) or result > 1e-10:
            raise ArithmeticError(
                "MSC gene topology probability failed numerical validation."
            )
        return min(result, 0.0)
