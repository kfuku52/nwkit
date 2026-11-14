"""Independent all-pair forest CTMC, analytic limits, and high-precision checks."""

import itertools
import math
from collections import defaultdict
from decimal import Decimal, localcontext
from functools import lru_cache

import numpy as np
import pytest
from scipy.linalg import expm

from nwkit.mul_coalescent import CoalescentTopology, log_lineage_transition
from tests.test_mul_reconcile import tree


@lru_cache(None)
def all_topologies(labels):
    if len(labels) == 1:
        return labels
    result = set()
    for mask in range(1, 2 ** len(labels) - 1):
        if not mask & 1:
            continue
        left = tuple(label for i, label in enumerate(labels) if mask & (1 << i))
        right = tuple(label for label in labels if label not in left)
        for first, second in itertools.product(
            all_topologies(left), all_topologies(right)
        ):
            result.add("(" + ",".join(sorted((first, second))) + ")")
    return tuple(sorted(result))


def canonical_gene(gene):
    values = {}
    for node in gene.traverse("postorder"):
        values[node] = (
            node.name
            if node.is_leaf
            else "(" + ",".join(sorted(values[c] for c in node.children)) + ")"
        )
    return values[gene]


def forest_successors(forest):
    for first, second in itertools.combinations(range(len(forest)), 2):
        remaining = [
            lineage for i, lineage in enumerate(forest) if i not in (first, second)
        ]
        remaining.append("(" + ",".join(sorted((forest[first], forest[second]))) + ")")
        yield tuple(sorted(remaining))


@lru_cache(None)
def oracle_branch(initial, time):
    """Track every pair merger, including histories incompatible with the target gene."""
    forests = [initial]
    positions = {initial: 0}
    edges = {}
    for state in forests:
        successors = tuple(forest_successors(state))
        edges[state] = successors
        for successor in successors:
            if successor not in positions:
                positions[successor] = len(forests)
                forests.append(successor)
    matrix = np.zeros((len(forests), len(forests)))
    for state, successors in edges.items():
        row = positions[state]
        matrix[row, row] = -len(successors)
        for successor in successors:
            matrix[row, positions[successor]] += 1
    probabilities = expm(matrix * time)[0]
    return {
        state: float(probabilities[i])
        for i, state in enumerate(forests)
        if probabilities[i] > 0
    }


@lru_cache(None)
def oracle_ancestral(forest):
    if len(forest) == 1:
        return {forest[0]: 1.0}
    result = defaultdict(float)
    successors = tuple(forest_successors(forest))
    for successor in successors:
        for topology, probability in oracle_ancestral(successor).items():
            result[topology] += probability / len(successors)
    return dict(result)


def oracle_distribution(species, assignment):
    active = {node: (name,) for name, node in assignment.items()}
    distributions = {}
    for node in species.traverse("postorder"):
        if node.is_leaf:
            current = {active.get(node, ()): 1.0}
        else:
            current = defaultdict(float)
            left, right = (distributions[c] for c in node.children)
            for (first, p), (second, q) in itertools.product(
                left.items(), right.items()
            ):
                current[tuple(sorted(first + second))] += p * q
        if not node.is_root:
            following = defaultdict(float)
            for forest, probability in current.items():
                for state, conditional in oracle_branch(forest, node.dist).items():
                    following[state] += probability * conditional
            current = following
        distributions[node] = current
    result = defaultdict(float)
    for forest, probability in distributions[species].items():
        for topology, conditional in oracle_ancestral(forest).items():
            result[topology] += probability * conditional
    return dict(result)


def decimal_lineage_transition(start, end, time):
    """Independent hypoexponential closed form, not production uniformization/expm."""
    precision = max(80, int(max(0, -math.log10(time)) * (start - end) + 60))
    with localcontext() as context:
        context.prec = precision
        duration = Decimal(str(time))
        rates = {i: Decimal(i * (i - 1)) / 2 for i in range(end, start + 1)}
        coefficient = math.prod(rates[i] for i in range(end + 1, start + 1))
        probability = coefficient * sum(
            (-rate * duration).exp()
            / math.prod(other - rate for j, other in rates.items() if j != i)
            for i, rate in rates.items()
        )
        return float(probability.ln())


@pytest.mark.parametrize(
    "duration", [1e-100, 1e-12, 0.001, 0.1, 0.5, 2.0, 20.0, 1000.0]
)
def test_lineage_transitions_against_decimal_closed_form(duration):
    for start in range(2, 8):
        total = 0.0
        for end in range(1, start + 1):
            observed = log_lineage_transition(start, end, duration)
            expected = decimal_lineage_transition(start, end, duration)
            assert observed == pytest.approx(expected, abs=3e-11, rel=3e-13)
            total += math.exp(observed)
        assert total == pytest.approx(1.0, abs=3e-13)


@pytest.mark.parametrize("duration", [0.0, 1e-12, 0.05, 0.5, 3.0, 50.0, 10000.0])
def test_three_taxon_analytic_probability_and_long_branch_logs(duration):
    species = tree(f"((A:1,B:1):{duration},C:{1 + duration});")
    assignments = {node.name: node for node in species.leaves()}
    concordant = CoalescentTopology(tree("((A,B),C);"))
    expected = math.log1p(-2 / 3 * math.exp(-duration))
    assert concordant.log_probability(species, assignments) == pytest.approx(
        expected, abs=3e-13
    )
    for text in ("((A,C),B);", "((B,C),A);"):
        observed = CoalescentTopology(tree(text)).log_probability(species, assignments)
        assert observed == pytest.approx(-duration - math.log(3), abs=3e-12)


@pytest.mark.parametrize(
    "text",
    [
        "((A:1,B:1):0.7,(C:1,D:1):0.7);",
        "(((A:1,B:1):0.3,C:1.3):0.4,D:1.7);",
        "(((A:1,B:1):0,C:1):0,D:1);",
        "((A:1,B:1):0.2,(C:1,(D:0.4,E:0.4):0.6):0.2);",
    ],
)
def test_all_rooted_topologies_against_all_pair_forest_oracle(text):
    species = tree(text)
    assignment = {node.name: node for node in species.leaves()}
    expected = oracle_distribution(species, assignment)
    topologies = all_topologies(tuple(sorted(assignment)))
    assert len(topologies) == (15 if len(assignment) == 4 else 105)
    total = 0.0
    for topology in topologies:
        observed = math.exp(
            CoalescentTopology(tree(topology + ";")).log_probability(
                species, assignment
            )
        )
        assert observed == pytest.approx(expected[topology], abs=2e-13, rel=2e-12)
        total += observed
    assert total == pytest.approx(1.0, abs=3e-13)


def test_missing_population_tips_are_conditioned_on_not_scored_as_losses():
    species = tree("((A:1,B:1):0.5,(C:1,D:1):0.5);")
    assignment = {node.name: node for node in species.leaves() if node.name != "B"}
    expected = oracle_distribution(species, assignment)
    for topology in all_topologies(tuple(sorted(assignment))):
        observed = math.exp(
            CoalescentTopology(tree(topology + ";")).log_probability(
                species, assignment
            )
        )
        assert observed == pytest.approx(expected[topology], abs=2e-13)


def test_time_scale_and_child_order_invariance():
    gene = tree("((A,C),(B,D));")
    species = tree("((B:20,A:20):10,(D:20,C:20):10);")
    assignment = {node.name: node for node in species.leaves()}
    observed = CoalescentTopology(gene).log_probability(
        species, assignment, time_scale=20
    )
    other = tree("((A:1,B:1):0.5,(C:1,D:1):0.5);")
    expected = CoalescentTopology(gene).log_probability(
        other, {n.name: n for n in other.leaves()}
    )
    assert observed == pytest.approx(expected, abs=2e-13)


@pytest.mark.parametrize(
    "start,end,duration",
    [(1, 2, 1), (3, 0, 1), (3, 2, -1), (3, 2, math.inf), (3, 2, math.nan)],
)
def test_invalid_transitions_fail(start, end, duration):
    with pytest.raises(ValueError, match="transition"):
        log_lineage_transition(start, end, duration)


def test_zero_duration_transitions_and_invalid_assignments():
    assert log_lineage_transition(3, 3, 0) == 0
    assert log_lineage_transition(3, 2, 0) == -math.inf
    species = tree("(A:1,B:1);")
    tips = list(species.leaves())
    solver = CoalescentTopology(tree("(x,y);"))
    with pytest.raises(ValueError, match="every gene tip"):
        solver.log_probability(species, {"x": tips[0]})
    with pytest.raises(ValueError, match="one copy"):
        solver.log_probability(species, {"x": tips[0], "y": tips[0]})
    with pytest.raises(ValueError, match="time scale"):
        solver.log_probability(species, {"x": tips[0], "y": tips[1]}, time_scale=0)


def test_state_cap_is_failure_not_approximation():
    species = tree("((A:1,B:1):1,C:2);")
    with pytest.raises(ValueError, match="max-coalescent-states"):
        CoalescentTopology(tree("((A,C),B);"), max_states=1).log_probability(
            species, {n.name: n for n in species.leaves()}
        )
    with pytest.raises(ValueError, match="positive"):
        CoalescentTopology(tree("(A,B);"), max_states=0)


def test_extreme_times_do_not_overflow_the_matrix_exponential():
    assert log_lineage_transition(7, 1, 1e308) == 0.0
    assert log_lineage_transition(7, 2, 1e308) == -1e308
    with pytest.raises(ArithmeticError, match="finite numerical range"):
        log_lineage_transition(7, 7, 1e308)


def test_matrix_work_is_bounded_before_dense_allocation():
    solver = CoalescentTopology(tree("(((A,B),C),D);"), max_states=10)
    with pytest.raises(ValueError, match="max-coalescent-states"):
        solver.branch(tuple(solver.leaves.values()), 1.0)
