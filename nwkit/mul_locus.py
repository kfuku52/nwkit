"""Linear locus birth/death followed by daughter-bounded multilocus coalescence.

The ancestral population starts with one locus on an explicit finite DL stem.
The original locus has an infinite ancestral coalescent population. Sampling
is one genome per extant locus; detection is applied AFTER genealogy sampling.
"""

import math
from collections import defaultdict
from dataclasses import dataclass, field
from functools import lru_cache
from numbers import Real
from typing import Any

import numpy as np
from scipy.special import logsumexp

from nwkit.mul_coalescent import log_lineage_transition
from nwkit.mul_reconcile_model import validate_binary
from nwkit.util import compute_node_ages


def finite_number(value):
    if isinstance(value, bool) or not isinstance(value, Real):
        return False
    try:
        return math.isfinite(value)
    except OverflowError:
        return False


@dataclass(eq=False)
class Locus:
    age: float
    kind: str
    species: str = ""
    children: list = field(default_factory=list)
    daughter: bool = False


@dataclass(frozen=True)
class LocusParameters:
    duplication: float
    loss: float
    ne: float
    hybridization_age: float

    def __post_init__(self):
        values = (self.duplication, self.loss, self.ne, self.hybridization_age)
        if any(not finite_number(v) for v in values) or min(values) < 0:
            raise ValueError("Locus parameters must be finite and nonnegative.")
        if not finite_number(self.duplication + self.loss):
            raise ValueError("Combined locus rate must be finite.")
        if self.ne <= 0 or not finite_number(2 * self.ne):
            raise ValueError("Locus Ne must be positive with finite 2*Ne.")


def sample_locus_tree(
    population, parameters, stem, rng, *, max_nodes=10000, stem_condition=None
):
    validate_binary(population, "Population tree", unique=False)
    if not finite_number(stem) or stem < 0 or max_nodes < 1:
        raise ValueError("Locus stem must be finite/nonnegative; node cap positive.")
    ages = compute_node_ages(population)
    origin_age = ages[population] + stem
    if not math.isfinite(origin_age):
        raise ValueError("Ancestral locus origin age must be finite.")
    origin = Locus(origin_age, "origin")
    pending = [(origin, population, origin.age, False)]
    nodes = 1
    rate = parameters.duplication + parameters.loss
    if not math.isfinite(rate):
        raise ValueError("Combined locus rate must be finite.")
    if stem_condition == "none":
        pending = [(origin, population, ages[population], False)]
    elif stem_condition == "birth":
        if stem <= 0 or parameters.duplication <= 0:
            raise ValueError("Ancestral birth stratum has zero prior probability.")
        wait = -math.log1p(rng.random() * math.expm1(-rate * stem)) / rate
        birth = Locus(origin.age - wait, "duplication")
        origin.children.append(birth)
        pending = [
            (birth, population, birth.age, True),
            (birth, population, birth.age, False),
        ]
        nodes += 1
    elif stem_condition is not None:
        raise ValueError("Unknown ancestral locus stratum.")
    while pending:
        parent, pop, start, daughter = pending.pop()
        wait = float(rng.exponential(1 / rate)) if rate else math.inf
        age = start - wait
        if age > ages[pop]:
            duplication = rng.random() < parameters.duplication / rate
            node = Locus(age, "duplication" if duplication else "loss")
            if duplication:
                pending.extend([(node, pop, age, True), (node, pop, age, False)])
        else:
            node = Locus(
                ages[pop],
                "tip" if pop.is_leaf else "speciation",
                pop.props.get("mul_species", pop.name) if pop.is_leaf else "",
            )
            pending.extend((node, child, node.age, False) for child in pop.children)
        node.daughter = daughter
        parent.children.append(node)
        nodes += 1
        if nodes > max_nodes:
            raise ValueError(
                "Locus simulation exceeded node cap; no history discarded."
            )
    return origin


@lru_cache(maxsize=100000)
def transition(k, j, time):
    return 0.0 if k == j == 0 else log_lineage_transition(k, j, time)


def combine_counts(distributions, charge=None):
    result = {0: 0.0}
    for distribution in distributions:
        joined: dict[int, float] = {}
        for k, first in result.items():
            for j, second in distribution.items():
                if charge is not None:
                    charge()
                joined[k + j] = float(
                    np.logaddexp(joined.get(k + j, -math.inf), first + second)
                )
        result = joined
    return result


class LocusCoalescent:
    """Count DP and backward sampling; no rejection of rare daughter bounds."""

    def __init__(self, root, ne, *, max_states=100000):
        if (
            not finite_number(ne)
            or ne <= 0
            or not finite_number(2 * ne)
            or max_states < 1
        ):
            raise ValueError("Coalescent Ne and work cap must be positive.")
        self.root, self.scale = root, 2 * ne
        self.base: dict[Locus, dict[int, float]] = {}
        self.output: dict[Locus, dict[int, float]] = {}
        self.times: dict[Locus, float] = {}
        self.bound_log_normalizers: dict[Locus, float] = {}
        self.work, self.max_states = 0, max_states
        self.order = []
        pending = [(root, None, False)]
        while pending:
            node, parent, visited = pending.pop()
            if not visited:
                pending.append((node, parent, True))
                pending.extend((child, node, False) for child in node.children)
                continue
            self.order.append(node)
            base = (
                {1 if node.kind == "tip" else 0: 0.0}
                if not node.children
                else combine_counts(
                    (self.output[child] for child in node.children), self._charge
                )
            )
            self.base[node] = base
            time = 0 if parent is None else (parent.age - node.age) / self.scale
            if not math.isfinite(time) or time < 0:
                raise ValueError("Invalid scaled locus branch duration.")
            self.times[node] = time
            output: dict[int, float] = defaultdict(lambda: -math.inf)
            for k, mass in base.items():
                for j in range(1, k + 1) if k else (0,):
                    self._charge((k - j + 1) ** 2)
                    output[j] = float(
                        np.logaddexp(output[j], mass + transition(k, j, time))
                    )
            if node.daughter and max(base) > 0:
                if not math.isfinite(output[1]):
                    raise ArithmeticError("Daughter bound has zero probability.")
                self.bound_log_normalizers[node] = output[1]
                output = {1: 0.0}
            self.output[node] = dict(output)

    def _charge(self, count=1):
        self.work += count
        if self.work > self.max_states:
            raise ValueError("Locus coalescent exceeded state cap.")

    @staticmethod
    def choose(options, rng):
        if len(options) == 1:
            key, log_weight = options[0]
            if not math.isfinite(log_weight):
                raise ArithmeticError("Invalid conditional coalescent weights.")
            # Preserve the uniform draw consumed by choice(..., p=[1]).
            rng.random()
            return key
        keys, logs = zip(*options, strict=True)
        values = np.asarray(logs)
        probabilities = np.exp(values - logsumexp(values))
        if np.any(~np.isfinite(probabilities)):
            raise ArithmeticError("Invalid conditional coalescent weights.")
        return keys[int(rng.choice(len(keys), p=probabilities))]

    def sample(self, rng):
        root_count = self.choose(list(self.base[self.root].items()), rng)
        requested = {self.root: root_count}
        # Choose all branch-end counts before constructing gene topologies.
        for node in reversed(self.order):
            j = requested[node]
            k = self.choose(
                [
                    (k, mass + transition(k, j, self.times[node]))
                    for k, mass in self.base[node].items()
                    if k >= j
                ],
                rng,
            )
            if len(node.children) == 1:
                requested[node.children[0]] = k
            elif len(node.children) == 2:
                left, right = node.children
                split = self.choose(
                    [
                        (a, mass + self.output[right][k - a])
                        for a, mass in self.output[left].items()
                        if k - a in self.output[right]
                    ],
                    rng,
                )
                requested[left], requested[right] = split, k - split
            elif node.children:
                raise ValueError("Locus nodes must have <=2 children.")
        forests: dict[Locus, list[Any]] = {}
        tip = 0
        for node in self.order:
            if node.children:
                forest = [
                    gene for child in node.children for gene in forests.pop(child)
                ]
            elif node.kind == "tip":
                forest = [("tip", node.species, tip)]
                tip += 1
            else:
                forest = []
            target = 1 if node is self.root and forest else requested[node]
            while len(forest) > target:
                a, b = sorted(
                    rng.choice(len(forest), size=2, replace=False), reverse=True
                )
                pair = sorted((forest.pop(a), forest.pop(b)))
                forest.append(("node", *pair))
            forests[node] = forest
        return forests[self.root][0] if forests[self.root] else None


def detected_signature(gene, detection, rng):
    if gene is None:
        return None
    if gene[0] == "tip":
        return ("tip", gene[1]) if rng.random() < detection[gene[1]] else None
    left, right = (detected_signature(child, detection, rng) for child in gene[1:])
    if left is None or right is None:
        return right if left is None else left
    return ("node", *sorted((left, right)))


def signature_size(signature):
    return (
        0
        if signature is None
        else 1
        if signature[0] == "tip"
        else sum(signature_size(c) for c in signature[1:])
    )


def sample_observation(population, parameters, model, rng, *, stem_condition=None):
    locus = sample_locus_tree(
        population,
        parameters,
        model["ancestral_stem"],
        rng,
        max_nodes=model["max_locus_nodes"],
        stem_condition=stem_condition,
    )
    gene = LocusCoalescent(
        locus, parameters.ne, max_states=model["max_coalescent_states"]
    ).sample(rng)
    return detected_signature(gene, model["detection"], rng)


def sample_selected(population, parameters, model, rng, *, count):
    observations: list[Any] = []
    attempts = 0
    while len(observations) < count:
        attempts += 1
        if attempts > model["max_attempts"]:
            raise ValueError("Locus selection exceeded attempt cap; no partial bank.")
        observation = sample_observation(population, parameters, model, rng)
        if 2 <= signature_size(observation) <= model["max_observed_tips"]:
            observations.append(observation)
    return observations, attempts
