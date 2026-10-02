"""Independent ODE transport and analytic WGD-marker references."""

import math

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from nwkit.wgd_count_model import CountTree, MultiplicationEvent
from nwkit.wgd_tree_model import GeneTopology, TopologyLikelihood


def species_tree():
    return CountTree(
        (-1, 0, 0),
        (0.0, 1.0, 1.0),
        (1, 2),
        ("A", "B"),
        (0, 1, 2),
        ("root", "A", "B"),
    )


@pytest.mark.parametrize("rates", [(0.0, 0.4), (0.3, 0.1), (0.1, 0.4), (0.2, 0.2)])
@pytest.mark.parametrize("length", [0.3, 2.0])
@pytest.mark.parametrize("initial_survival", [0.2, 0.7, 1.0])
def test_exact_marked_branch_flow_matches_independent_ode(
    rates, length, initial_survival
):
    topology = GeneTopology(
        ((), (), (0, 1), (), (), (3, 4), (2, 5)),
        ("A", "A", "", "A", "A", "", ""),
        tuple(str(node) for node in range(7)),
    )
    model = TopologyLikelihood(species_tree(), topology, detection=[0.7, 0.0])
    initial = np.array(
        [
            [0.07, 0.0],
            [0.10, 0.0],
            [0.015, 0.003],
            [0.08, 0.0],
            [0.05, 0.0],
            [0.01, 0.0],
            [0.002, 0.0002],
        ]
    )
    duplication, loss = rates

    def derivative(_time, flat):
        survival = flat[0]
        values = flat[1:].reshape(initial.shape)
        change = (duplication - loss - 2 * duplication * survival) * values
        # Every split has identical colored child shapes, so its coefficient is one.
        for node, left, right in ((2, 0, 1), (5, 3, 4), (6, 2, 5)):
            change[node, 0] += duplication * values[left, 0] * values[right, 0]
            change[node, 1] += duplication * (
                values[left, 1] * values[right, 0] + values[left, 0] * values[right, 1]
            )
        survival_change = (
            duplication - loss
        ) * survival - duplication * survival * survival
        return np.concatenate(([survival_change], change.ravel()))

    reference = solve_ivp(
        derivative,
        (0.0, length),
        np.concatenate(([initial_survival], initial.ravel())),
        method="DOP853",
        rtol=1e-12,
        atol=1e-16,
    )
    assert reference.success, reference.message
    logs = np.full_like(initial, -np.inf)
    positive = initial > 0
    logs[positive] = np.log(initial[positive])
    survival, flow = model._edge(
        initial_survival,
        logs,
        duplication,
        loss,
        length,
        np.longdouble,
        math.log(initial_survival),
    )
    expected = reference.y[1:, -1].reshape(initial.shape)
    np.testing.assert_allclose(np.exp(flow), expected, rtol=3e-10, atol=1e-12)
    assert abs(survival - reference.y[0, -1]) <= 1e-11


@pytest.mark.parametrize("loss", [0.0, 1.0, 40.0])
@pytest.mark.parametrize("root_mean", [1.0, 2.0, 10.0])
@pytest.mark.parametrize("fraction", [0.0, 0.4, 1.0])
@pytest.mark.parametrize("retention", [0.01, 0.8, 1.0])
def test_wgd_marker_matches_independent_pure_loss_geometric_root_reference(
    loss, root_mean, fraction, retention
):
    topology = GeneTopology(((), (), (0, 1)), ("A", "A", ""), ("a", "b", "root"))
    model = TopologyLikelihood(species_tree(), topology, detection=[0.7, 0.0])

    after_event = 0.7 * math.exp(-loss * (1 - fraction))
    before_event = math.exp(-loss * fraction)
    one_copy = (
        before_event * after_event * (1 + retention - 2 * retention * after_event)
    )
    two_copies = before_event * retention * after_event * after_event
    survival = before_event * after_event * (1 + retention - retention * after_event)
    geometric_ratio = 1 - 1 / root_mean
    denominator = 1 / root_mean + geometric_ratio * survival
    # Two observed copies arise either at WGD or from two surviving root copies.
    expected = two_copies / (
        two_copies + (geometric_ratio / denominator) * one_copy * one_copy
    )

    result = model.evaluate(
        [[0.0, loss]],
        root_mean,
        MultiplicationEvent(1, retention, fraction),
    )
    assert result.origin_node_ids == ("root",)
    assert abs(float(result.origin_probabilities[0]) - expected) <= 1e-7
