"""Independent forward histories validate colored topology mass and WGD origins.

The simulator uses exponential event times, not NWKIT transition probabilities,
pruning, or simulation. Lost/undetected tips and unary ancestral nodes are removed.
"""

import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass

import numpy as np
import pytest

from nwkit.wgd_count_model import CountTree, MultiplicationEvent
from nwkit.wgd_tree_model import GeneTopology, TopologyLikelihood


@dataclass(frozen=True)
class History:
    species: str = ""
    children: tuple["History", ...] = ()
    origin: str = ""

    @property
    def shape(self):
        return (
            ("tip", self.species)
            if not self.children
            else ("node", *sorted(child.shape for child in self.children))
        )

    @property
    def wgd_shapes(self):
        result = Counter({self.shape: 1}) if self.origin == "wgd" else Counter()
        for child in self.children:
            result.update(child.wgd_shapes)
        return result


def observed_split(origin, first, second):
    if first is None:
        return second
    if second is None:
        return first
    return History(children=(first, second), origin=origin)


class ForwardHistories:
    def __init__(self, tree, rates, detection, scales, root_mean, event, seed):
        self.tree, self.rates, self.detection = tree, rates, detection
        self.scales, self.root_mean, self.event = scales, root_mean, event
        self.rng = np.random.default_rng(seed)
        self.tip_species = dict(zip(tree.tip_nodes, tree.tip_names, strict=True))
        self.tip_detection = dict(zip(tree.tip_nodes, detection, strict=True))

    def _segment(self, duration, duplication, loss, finish, origin="ssd"):
        rate = duplication + loss
        wait = self.rng.exponential(1 / rate) if rate else math.inf
        if wait >= duration:
            return finish()
        if self.rng.random() < loss / rate:
            return None
        return observed_split(
            origin,
            self._segment(duration - wait, duplication, loss, finish, origin),
            self._segment(duration - wait, duplication, loss, finish, origin),
        )

    def _species_node(self, node, scale):
        if node in self.tip_species:
            return (
                History(species=self.tip_species[node])
                if self.rng.random() < self.tip_detection[node]
                else None
            )
        left, right = self.tree.children[node]
        return observed_split(
            "speciation", self._branch(left, scale), self._branch(right, scale)
        )

    def _branch(self, node, scale):
        duplication, loss = self.rates[node] * scale
        length = self.tree.lengths[node]

        def endpoint():
            return self._species_node(node, scale)

        if self.event is None or self.event.node != node:
            return self._segment(length, duplication, loss, endpoint)
        event = self.event

        def after_event():
            return self._segment(
                length * (1 - event.fraction), duplication, loss, endpoint
            )

        def pulse():
            first = after_event()
            if self.rng.random() >= event.retention:
                return first
            return observed_split("wgd", first, after_event())

        return self._segment(length * event.fraction, duplication, loss, pulse)

    def sample(self):
        scale = self.scales[int(self.rng.integers(len(self.scales)))]
        return self._segment(
            1,
            math.log(self.root_mean),
            0,
            lambda: self._species_node(0, scale),
            origin="root_prior",
        )


def four_species_tree():
    return CountTree(
        (-1, 0, 1, 1, 0, 4, 4),
        (0.0, 0.7, 0.6, 0.6, 0.8, 0.5, 0.5),
        (2, 3, 5, 6),
        ("A", "B", "C", "D"),
        tuple(range(7)),
        ("root", "AB", "A", "B", "CD", "C", "D"),
    )


def shape_topology(shape):
    children, species = [], []

    def visit(item):
        descendants = () if item[0] == "tip" else (visit(item[1]), visit(item[2]))
        children.append(descendants)
        species.append(item[1] if item[0] == "tip" else "")
        return len(children) - 1

    visit(shape)
    return GeneTopology(
        tuple(children), tuple(species), tuple(str(i) for i in range(len(children)))
    )


def test_forward_simulator_perfect_internal_wgd_keeps_the_true_origin():
    tree = four_species_tree()
    simulator = ForwardHistories(
        tree,
        np.zeros((7, 2)),
        (1.0,) * 4,
        (1.0,),
        1.0,
        MultiplicationEvent(1, 1.0, 0.35),
        615,
    )
    a, b, c, d = (("tip", name) for name in "ABCD")
    ab, cd = ("node", a, b), ("node", c, d)
    wgd_shape = ("node", ab, ab)
    expected = ("node", *sorted((wgd_shape, cd)))
    for _ in range(20):
        history = simulator.sample()
        assert history.shape == expected
        assert history.wgd_shapes == {wgd_shape: 1}


@pytest.mark.slow
@pytest.mark.parametrize(
    "name,rates,root_mean,detection,scales,event,selection,seed",
    [
        ("null", (0.15, 0.25), 1.5, (0.85,) * 4, (1.0,), None, "observed", 611),
        (
            "internal_wgd_mixture",
            (0.15, 0.20),
            1.5,
            (0.85, 0.65, 0.9, 1.0),
            (0.2, 1.8),
            MultiplicationEvent(1, 0.65, 0.35),
            "root-clades",
            612,
        ),
        (
            "strong_loss_terminal_wgd",
            (0.2, 1.2),
            2.2,
            (0.8, 0.6, 0.9, 0.7),
            (1.0,),
            MultiplicationEvent(2, 0.8, 0.4),
            "observed",
            613,
        ),
        (
            "ssd_after_wgd",
            (0.55, 0.15),
            1.2,
            (0.9,) * 4,
            (1.0,),
            MultiplicationEvent(1, 0.5, 0.5),
            "root-clades",
            614,
        ),
        (
            "heterogeneous_background",
            ((0.12, 0.3), (0.3, 0.15)),
            1.8,
            (0.8, 0.65, 0.9, 0.7),
            (0.3, 1.7),
            None,
            "root-clades",
            617,
        ),
        (
            "critical_wgd",
            (0.35, 0.35),
            2.0,
            (0.8, 0.6, 0.9, 0.7),
            (1.0,),
            MultiplicationEvent(4, 0.4, 0.75),
            "observed",
            618,
        ),
        (
            "retained_pulse_only",
            (0.0, 0.0),
            1.0,
            (0.8, 0.65, 0.9, 1.0),
            (1.0,),
            MultiplicationEvent(1, 0.7, 0.5),
            "root-clades",
            619,
        ),
    ],
)
def test_topology_and_origin_probabilities_match_independent_forward_histories(
    name, rates, root_mean, detection, scales, event, selection, seed
):
    tree = four_species_tree()
    categories = np.atleast_2d(rates)
    groups = (0, 0, 1, 1, 0, 1, 1) if len(categories) == 2 else (0,) * 7
    branch_rates = categories[list(groups)]
    simulator = ForwardHistories(
        tree, branch_rates, detection, scales, root_mean, event, seed
    )
    observed = defaultdict(list)
    proposals, selected = 60000, 0
    for _ in range(proposals):
        history = simulator.sample()
        if history is None:
            continue
        topology = shape_topology(history.shape)
        labels = set(topology.species) - {""}
        if selection == "root-clades" and not (
            labels & {"A", "B"} and labels & {"C", "D"}
        ):
            continue
        selected += 1
        if len(topology.children) <= 11:
            observed[history.shape].append(history.wgd_shapes)

    checked, probability_errors, probability_z, origin_z = 0, [], [], []
    for shape, markers in sorted(observed.items()):
        if len(markers) < 120:
            continue
        model = TopologyLikelihood(
            tree,
            shape_topology(shape),
            detection=detection,
            branch_groups=groups,
            rate_scales=scales,
            ascertainment=selection,
        )
        prediction = model.evaluate(categories, root_mean, event)
        probability = math.exp(prediction.log_likelihood)
        frequency = len(markers) / selected
        se = math.sqrt(probability * (1 - probability) / selected)
        error = abs(frequency - probability)
        assert error <= 6 * se + 3 / selected, (name, shape, frequency, probability)
        expected_origins = defaultdict(float)
        node_multiplicities = Counter()
        signatures = model.gene.signatures
        for node, probability_origin in zip(
            model.targets, prediction.origin_probabilities, strict=True
        ):
            expected_origins[signatures[node]] += float(probability_origin)
            node_multiplicities[signatures[node]] += 1
        for subtree, expected in expected_origins.items():
            assignments = [counts[subtree] for counts in markers]
            empirical = float(np.mean(assignments))
            # A zero empirical variance is not evidence against rare events.
            # For X in [0,K], Var(X) <= E(X)*(K-E(X)); K=1 gives exact binomial variance.
            variance_bound = max(
                0.0, expected * (node_multiplicities[subtree] - expected)
            )
            marker_se = math.sqrt(variance_bound / len(markers))
            if expected == 0:
                assert empirical == 0, (name, shape, subtree, empirical)
            assert abs(empirical - expected) <= (6 * marker_se + 3 / len(markers)), (
                name,
                shape,
                subtree,
                empirical,
                expected,
            )
            origin_z.append(abs(empirical - expected) / marker_se if marker_se else 0)
        checked += 1
        probability_errors.append(error)
        probability_z.append(error / se if se else 0)
    assert checked >= 10
    print(
        json.dumps(
            {
                "study": name,
                "seed": seed,
                "proposals": proposals,
                "conditioned_histories": selected,
                "common_topologies_checked": checked,
                "maximum_probability_absolute_error": max(probability_errors),
                "maximum_probability_standard_errors": max(probability_z),
                "maximum_origin_bounded_standard_errors": max(origin_z),
            },
            sort_keys=True,
        )
    )
