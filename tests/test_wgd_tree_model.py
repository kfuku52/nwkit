import math

import numpy as np
import pytest

from nwkit.wgd_count_model import (
    CountLikelihood,
    CountTree,
    MultiplicationEvent,
    birth_death_transition,
)
from nwkit.wgd_tree_model import GeneTopology, TopologyLikelihood, lineage_survival


def species():
    return CountTree(
        (-1, 0, 0), (0, 1, 1), (1, 2), ("A", "B"), (0, 1, 2), ("root", "A", "B")
    )


def shapes(n):
    if n == 1:
        return {("tip", "A")}
    return {
        ("node", *sorted((a, b)))
        for i in range(1, n)
        for a in shapes(i)
        for b in shapes(n - i)
    }


def gene(shape):
    children, labels = [], []

    def visit(item):
        descendants = () if item[0] == "tip" else (visit(item[1]), visit(item[2]))
        children.append(descendants)
        labels.append(item[1] if item[0] == "tip" else "")
        return len(children) - 1

    visit(shape)
    return GeneTopology(
        tuple(children), tuple(labels), tuple(f"g{i}" for i in range(len(children)))
    )


@pytest.mark.parametrize("retention", [None, 0, 0.1, 0.7, 1.0])
def test_likelihood_only_agrees_with_full_node_origin_calculation(retention):
    a, b = ("tip", "A"), ("tip", "B")
    topology = gene(("node", ("node", a, a), ("node", a, b)))
    model = TopologyLikelihood(
        species(),
        topology,
        detection=[0.8, 0.9],
        rate_scales=[0.2, 1.8],
        ascertainment="root-clades",
    )
    event = None if retention is None else MultiplicationEvent(1, retention, 0.25)
    args = [[0.12, 0.16]], 1.4, event
    all_nodes = model.evaluate(*args)
    only = model.evaluate(*args, origins=False)
    assert only.log_likelihood == all_nodes.log_likelihood
    assert only.log_likelihood_error == all_nodes.log_likelihood_error
    assert only.origin_node_ids == () and len(only.origin_probabilities) == 0
    assert only.origin_probability_error == 0


def test_fixed_background_doubling_never_counts_as_candidate_origin():
    a, b = ("tip", "A"), ("tip", "B")
    topology = gene(("node", ("node", a, a), ("node", b, b)))
    fixed, tested = MultiplicationEvent(1, 1.0), MultiplicationEvent(2, 1.0)
    model = TopologyLikelihood(
        species(), topology, fixed_events=(fixed,), ascertainment="root-clades"
    )
    combined = model.evaluate([[0.0, 0.0]], 1.0, tested)
    assert combined.log_likelihood == pytest.approx(0.0)
    np.testing.assert_array_equal(combined.origin_probabilities, [0.0, 1.0, 0.0])
    assert (
        model.evaluate([[0.0, 0.0]], 1.0, tested, origins=False).log_likelihood
        == combined.log_likelihood
    )
    with pytest.raises(ValueError, match="overlap"):
        model.evaluate([[0.0, 0.0]], 1.0, fixed)
    with pytest.raises(ValueError, match="all fixed background events to be doublings"):
        TopologyLikelihood(
            species(), topology, fixed_events=(MultiplicationEvent(1, 0.5, 0.5, 3),)
        )


def test_fixed_background_and_single_candidate_are_equal_without_origin_markers():
    a, b = ("tip", "A"), ("tip", "B")
    topology = gene(("node", ("node", a, a), b))
    event = MultiplicationEvent(1, 0.65, 0.25)
    settings = dict(
        detection=(0.7, 0.9), rate_scales=(0.4, 1.6), ascertainment="root-clades"
    )
    model = TopologyLikelihood(species(), topology, fixed_events=(event,), **settings)
    baseline = TopologyLikelihood(species(), topology, **settings)
    args = [[0.1, 0.22]], 1.12
    assert (
        model.evaluate(*args).log_likelihood
        == baseline.evaluate(*args, event).log_likelihood
    )
    assert not np.any(model.evaluate(*args).origin_probabilities)


@pytest.mark.parametrize("with_background", [False, True])
@pytest.mark.parametrize("fraction", [0.0, 0.4, 1.0])
def test_zero_retention_topology_equals_background_exactly(with_background, fraction):
    a, b = ("tip", "A"), ("tip", "B")
    topology = gene(("node", ("node", a, a), b))
    fixed = (MultiplicationEvent(1, 0.65, 0.25),) if with_background else ()
    model = TopologyLikelihood(
        species(),
        topology,
        fixed_events=fixed,
        detection=(0.75, 0.9),
        rate_scales=(0.4, 1.6),
        ascertainment="root-clades",
    )
    args = [[0.13, 0.22]], 1.16
    baseline = model.evaluate(*args)
    result = model.evaluate(*args, MultiplicationEvent(2, 0.0, fraction))
    assert result.log_likelihood == baseline.log_likelihood
    assert result.log_likelihood_error == baseline.log_likelihood_error
    assert result.origin_node_ids == baseline.origin_node_ids
    assert result.origin_probability_error == 0
    np.testing.assert_array_equal(result.origin_probabilities, [0.0, 0.0])


def test_one_ancestral_copy_rejects_crossed_root_before_numerical_evaluation(
    monkeypatch,
):
    a, b = ("tip", "A"), ("tip", "B")
    crossed = TopologyLikelihood(
        species(),
        gene(("node", ("node", a, b), ("node", a, b))),
        ascertainment="root-clades",
    )
    compatible = TopologyLikelihood(
        species(),
        gene(("node", ("node", a, a), ("node", b, b))),
        ascertainment="root-clades",
    )
    assert crossed.requires_multiple_root_copies
    assert not compatible.requires_multiple_root_copies
    assert np.isfinite(crossed.evaluate([[0.2, 0.3]], 1.2).log_likelihood)
    assert np.isfinite(compatible.evaluate([[0.2, 0.3]], 1.0).log_likelihood)

    def error(*args, **kwargs):
        pytest.fail("An impossible topology should not enter the pruning calculation")

    monkeypatch.setattr(crossed, "_evaluate", error)
    with pytest.raises(ValueError, match="impossible with one ancestral root copy"):
        crossed.evaluate([[0.2, 0.3]], 1.0)


@pytest.mark.parametrize("n", [1, 2, 3, 4, 5])
def test_small_colored_topology_sum_matches_independent_count_probability(n):
    duplication, loss = 0.35, 0.2
    count = birth_death_transition(duplication, loss, 1, 12)[1, n]
    survival = lineage_survival(1, duplication, loss, 1)
    likelihoods = [
        math.exp(
            TopologyLikelihood(species(), gene(shape), detection=[1, 0])
            .evaluate([[duplication, loss]], 1)
            .log_likelihood
        )
        for shape in shapes(n)
    ]
    assert sum(likelihoods) == pytest.approx(count / survival, rel=2e-6)


def test_four_tip_yule_topologies_have_the_known_one_third_two_thirds_shape_weights():
    probabilities = sorted(
        math.exp(
            TopologyLikelihood(species(), gene(shape), detection=[1, 0])
            .evaluate([[0.3, 0]], 1)
            .log_likelihood
        )
        for shape in shapes(4)
    )
    total = sum(probabilities)
    assert [value / total for value in probabilities] == pytest.approx(
        [1 / 3, 2 / 3], rel=1e-6
    )


def test_perfect_wgd_has_unit_conditional_assignment_without_ssd():
    topology = gene(next(iter(shapes(2))))
    result = TopologyLikelihood(species(), topology, detection=[1, 0]).evaluate(
        [[0, 0]], 1, MultiplicationEvent(1, 1)
    )
    assert result.log_likelihood == pytest.approx(0)
    assert result.origin_probabilities == pytest.approx([1])
    assert result.log_likelihood_error < 1e-6


def test_yule_root_count_prior_is_not_scaled_by_family_rate_categories():
    topology = gene(next(iter(shapes(2))))
    first = TopologyLikelihood(
        species(), topology, detection=[1, 0], rate_scales=[1]
    ).evaluate([[0, 0]], 2)
    second = TopologyLikelihood(
        species(), topology, detection=[1, 0], rate_scales=[0.1, 10]
    ).evaluate([[0, 0]], 2)
    assert math.exp(first.log_likelihood) == pytest.approx(0.25, rel=1e-7)
    assert first.log_likelihood == pytest.approx(second.log_likelihood, abs=1e-7)


def test_strong_loss_uses_direct_survival_and_scaled_topology_likelihood():
    assert lineage_survival(1, 0, 40, 1) == pytest.approx(math.exp(-40), rel=1e-12)
    result = TopologyLikelihood(
        species(), gene(("tip", "A")), detection=[1, 0]
    ).evaluate([[0, 40]], 1)
    assert result.log_likelihood == pytest.approx(0)
    result = TopologyLikelihood(
        species(), gene(("tip", "A")), detection=[1, 0]
    ).evaluate([[1e-12, 40]], 1)
    assert result.log_likelihood == pytest.approx(0, abs=1e-6)


def test_wgd_after_ssd_does_not_assign_the_older_node_to_wgd():
    tip = ("tip", "A")
    cherry = ("node", tip, tip)
    topology = gene(("node", cherry, cherry))
    result = TopologyLikelihood(species(), topology, detection=[1, 0]).evaluate(
        [[0.0001, 0]], 1, MultiplicationEvent(1, 1)
    )
    assert np.all(result.origin_probabilities[:2] > 0.99)
    assert result.origin_probabilities[-1] < 0.01


def test_category_assignments_are_likelihood_weighted_not_averaged():
    model = TopologyLikelihood(
        species(),
        gene(next(iter(shapes(2)))),
        detection=[0.9, 0],
        rate_scales=[0.01, 10],
    )
    rates, event = np.array([[0.3, 0.2]]), MultiplicationEvent(1, 0.8)
    categories = [
        model._category(rates, 1.5, event, scale, np.float64) for scale in model.scales
    ]
    numerator = sum(math.exp(item[0][1]) for item in categories)
    denominator = sum(math.exp(item[0][0]) for item in categories)
    result = model.evaluate(rates, 1.5, event)
    assert result.origin_probabilities[0] == pytest.approx(
        numerator / denominator, rel=1e-5
    )


def test_root_clade_selection_and_model_limits_are_explicit():
    topology = gene(("tip", "A"))
    with pytest.raises(ValueError, match="root-clade"):
        TopologyLikelihood(species(), topology, ascertainment="root-clades")
    with pytest.raises(ValueError, match="doubling"):
        TopologyLikelihood(species(), topology).evaluate(
            [[0.1, 0.2]], 1, MultiplicationEvent(1, 0.5, multiplicity=3)
        )
    with pytest.raises(ValueError, match="zero detection"):
        TopologyLikelihood(species(), topology, detection=[0, 1])


@pytest.mark.parametrize("n", [2, 3, 4])
@pytest.mark.parametrize("loss", [0.35, 0.25])
def test_wgd_and_geometric_root_topology_sum_matches_count_pruning(n, loss):
    rates = np.array([[0.25, loss]])
    event = MultiplicationEvent(1, 0.65)
    count = CountLikelihood(
        species(), np.array([[n, 0]]), detection=[0.9, 0.8], rate_scales=[0.2, 1.8]
    )
    expected = math.exp(count.log_likelihood(rates, 1.5, 96, event))
    likelihoods = [
        math.exp(
            TopologyLikelihood(
                species(), gene(shape), detection=[0.9, 0.8], rate_scales=[0.2, 1.8]
            )
            .evaluate(rates, 1.5, event)
            .log_likelihood
        )
        for shape in shapes(n)
    ]
    assert sum(likelihoods) == pytest.approx(expected, rel=3e-6)


@pytest.mark.parametrize("selection", ["observed", "root-clades"])
def test_multispecies_colored_topologies_sum_to_count_probability(selection):
    a, b = ("tip", "A"), ("tip", "B")
    topologies = [("node", ("node", a, a), b), ("node", ("node", a, b), a)]
    rates = np.array([[0.25, 0.35]])
    event = MultiplicationEvent(1, 0.65)
    count = CountLikelihood(
        species(), np.array([[2, 1]]), detection=[0.9, 0.8], ascertainment=selection
    )
    expected = math.exp(count.log_likelihood(rates, 1.5, 96, event))
    actual = sum(
        math.exp(
            TopologyLikelihood(
                species(), gene(shape), detection=[0.9, 0.8], ascertainment=selection
            )
            .evaluate(rates, 1.5, event)
            .log_likelihood
        )
        for shape in topologies
    )
    assert actual == pytest.approx(expected, rel=3e-6)
