"""Independent regressions for fixed-topology identity and numerical tails."""

import math

import pytest

from nwkit.wgd_count_model import CountLikelihood, CountTree
from nwkit.wgd_tree import model_mapping
from nwkit.wgd_tree_model import GeneTopology, TopologyLikelihood


def species(lengths=(1.0, 1.0)):
    return CountTree(
        (-1, 0, 0),
        (0.0, *lengths),
        (1, 2),
        ("A", "B"),
        (0, 1, 2),
        ("root", "A", "B"),
    )


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


def test_full_species_tree_identity_rejects_relative_change_in_tiny_branch():
    original = species((1e-13, 1e-13))
    metadata = {
        "species_tree": {
            "parents": list(original.parents),
            "lengths": list(original.lengths),
        },
        "species_event_ids": list(original.clade_ids),
        "tip_names": list(original.tip_names),
        "branch_rate_groups": [0, 0, 0],
    }
    assert model_mapping(metadata, original) == ({0: 0, 1: 1, 2: 2}, (0, 0, 0))
    with pytest.raises(ValueError, match="branch lengths"):
        model_mapping(metadata, species((2e-13, 1e-13)))


def test_rare_extinction_retains_single_species_topology_probability():
    loss = 1e-18
    rates = [[0.0, loss]]
    extinction = -math.expm1(-loss)
    expected = -loss + math.log(extinction) - math.log1p(-(extinction**2))
    reference = CountLikelihood(species(), [[1, 0]]).log_likelihood(rates, 1.0, 16)
    assert reference == pytest.approx(expected, abs=1e-12)
    result = TopologyLikelihood(species(), gene(("tip", "A"))).evaluate(rates, 1.0)
    assert math.isfinite(result.log_likelihood)
    assert result.log_likelihood == pytest.approx(reference, abs=1e-7)


def test_root_clade_selection_stays_in_log_space_when_joint_survival_underflows():
    rates = [[0.0, 400.0]]
    reference = CountLikelihood(
        species(), [[1, 1]], ascertainment="root-clades"
    ).log_likelihood(rates, 1.0, 16)
    assert reference == pytest.approx(0.0, abs=1e-12)
    result = TopologyLikelihood(
        species(),
        gene(("node", ("tip", "A"), ("tip", "B"))),
        ascertainment="root-clades",
    ).evaluate(rates, 1.0)
    assert math.isfinite(result.log_likelihood)
    assert result.log_likelihood == pytest.approx(reference, abs=1e-7)


@pytest.mark.parametrize("tip_count", [32, 64, 128])
def test_comb_matches_exact_yule_colored_topology_probability(tip_count):
    tip = ("tip", "A")
    shape = tip
    for _ in range(tip_count - 1):
        shape = ("node", shape, tip)
    # A comb has one symmetric cherry; the other splits have coefficient two.
    shape_log_weight = (tip_count - 2) * math.log(2) - math.lgamma(tip_count)
    expected = shape_log_weight - 1 + (tip_count - 1) * math.log(-math.expm1(-1))
    if tip_count == 32:
        assert expected == pytest.approx(-72.516737644, abs=1e-9)
    result = TopologyLikelihood(species(), gene(shape), detection=[1.0, 0.0]).evaluate(
        [[1.0, 0.0]], 1.0
    )
    assert result.log_likelihood == pytest.approx(expected, abs=1e-6)
    assert result.log_likelihood_error <= 1e-6
