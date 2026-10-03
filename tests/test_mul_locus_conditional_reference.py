import importlib.util
import math
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from nwkit.mul_locus import Locus
from tests.test_mul_locus import color, nested_locus_with_loss, strip_ids
from tests.test_mul_locus_integral import joint_oracle


@pytest.fixture
def reference():
    path = Path(__file__).resolve().parents[1] / "examples/mul-locus/reference.py"
    spec = importlib.util.spec_from_file_location("conditional_reference", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_independent_conditional_topologies_match_joint_full_forest(reference):
    locus, _ = nested_locus_with_loss()
    expected = Counter()
    for text, probability in joint_oracle(locus, 0.7).items():
        expected[color(text)] += probability
    sampler = reference.ConditionalReference(locus, 0.7)
    rng, n = np.random.default_rng(20261110), 12000
    actual = Counter(strip_ids(sampler.sample(rng)) for _ in range(n))
    assert set(actual) <= set(expected)
    for signature, probability in expected.items():
        assert abs(actual[signature] / n - probability) <= (
            6 * math.sqrt(probability * (1 - probability) / n) + 6 / n
        )


@pytest.mark.parametrize("duration", [1e-100, 1e-8, 0.1])
def test_rare_daughter_is_conditioned_without_rejection(reference, duration):
    daughter = Locus(
        0,
        "speciation",
        children=[Locus(0, "tip", "A"), Locus(0, "tip", "B")],
        daughter=True,
    )
    root = Locus(duration, "origin", children=[daughter])
    sampler = reference.ConditionalReference(root, 3)
    for seed in range(10):
        gene = sampler.sample(np.random.default_rng(seed))
        assert strip_ids(gene) == ("node", ("tip", "A"), ("tip", "B"))


def test_reference_caps_and_impossible_bound_abort(reference):
    locus, _ = nested_locus_with_loss()
    with pytest.raises(ValueError, match="cap"):
        reference.ConditionalReference(locus, 1, max_work=1)
    daughter = Locus(
        0,
        "speciation",
        children=[Locus(0, "tip", "A"), Locus(0, "tip", "B")],
        daughter=True,
    )
    with pytest.raises(ArithmeticError, match="probability"):
        reference.ConditionalReference(Locus(0, "origin", children=[daughter]), 1)


def test_reference_ignores_only_pure_death_structural_roundoff(reference, monkeypatch):
    original = reference.expm

    def perturbed(generator):
        matrix = original(generator)
        if len(matrix) > 2:
            matrix[1, 2] = -1e-16
        return matrix

    monkeypatch.setattr(reference, "expm", perturbed)
    locus = Locus(
        1,
        "origin",
        children=[
            Locus(
                0, "speciation", children=[Locus(0, "tip", "A"), Locus(0, "tip", "B")]
            )
        ],
    )
    assert (
        reference.ConditionalReference(locus, 1).sample(np.random.default_rng(1))
        is not None
    )

    def invalid(generator):
        matrix = original(generator)
        if len(matrix) > 2:
            matrix[2, 1] = -1e-16
        return matrix

    monkeypatch.setattr(reference, "expm", invalid)
    with pytest.raises(ArithmeticError, match="transition"):
        reference.ConditionalReference(locus, 1)


def test_reference_detection_preserves_hidden_loci_and_extinction(reference):
    gene = ("node", ("tip", "A", 0), ("tip", "B", 1))
    assert reference.detect(gene, {"A": 1, "B": 0}, np.random.default_rng(1)) == (
        "tip",
        "A",
    )
    root = Locus(1, "origin", children=[Locus(0, "loss")])
    assert (
        reference.ConditionalReference(root, 1).sample(np.random.default_rng(1)) is None
    )
