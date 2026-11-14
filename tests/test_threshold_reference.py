"""End-to-end comparisons with independently integrated threshold posteriors."""

import numpy as np
import pytest

from nwkit.threshold_asr import compute_threshold_marginals
from tests.test_threshold_asr import likelihoods, tree_from
from tests.threshold_posterior_support import (
    binary_internal_reference,
    binary_root_reference,
    ordinal_star_reference,
)


@pytest.mark.parametrize("internal", [False, True])
def test_binary_posterior_matches_independent_quadrature(internal):
    tree = tree_from("((A:1,B:0.5)I:0.7,C:2)R;" if internal else "(A:1,B:2)R;")
    observed = {"A": "low", "B": "high", **({"C": "high"} if internal else {})}
    states = ("low", "high")
    reference = binary_internal_reference() if internal else binary_root_reference()
    node = tree["I"] if internal else tree
    probabilities, fit = compute_threshold_marginals(
        tree,
        states,
        observed,
        likelihoods(tree, states, observed),
        num_samples=2000,
        burnin=500,
        chains=4,
        seed=724,
    )
    # Fixed absolute tolerances are independent of the diagnostic implementation;
    # multi-seed error/MCSE calibration is performed by the validation tool.
    np.testing.assert_allclose(
        probabilities[node], reference["probabilities"], atol=0.025, rtol=0
    )
    assert fit.liability_marginals[node].mean == pytest.approx(
        reference["mean"], abs=0.04
    )
    assert fit.liability_marginals[node].variance == pytest.approx(
        reference["variance"], abs=0.05
    )


@pytest.mark.parametrize("ambiguous", [False, True])
def test_estimated_ordinal_posterior_matches_independent_integration(ambiguous):
    states = ("low", "mid", "high")
    tree = tree_from("(A:1,B:1,C:1,D:2)R;" if ambiguous else "(A:1,B:1,C:1)R;")
    observed = {"A": "low", "B": "mid", "C": "high"}
    weights = {
        "A": np.array([1, 0, 0]),
        "B": np.array([0, 1, 0]),
        "C": np.array([0, 0, 1]),
    }
    if ambiguous:
        observed["D"] = "low|high"
        weights["D"] = np.array([1, 0, 1])
    reference = ordinal_star_reference(ambiguous=ambiguous)
    tighter = ordinal_star_reference(ambiguous=ambiguous, tolerance=1e-10)
    for key in reference:
        np.testing.assert_allclose(reference[key], tighter[key], rtol=1e-7, atol=1e-9)
    probabilities, fit = compute_threshold_marginals(
        tree,
        states,
        observed,
        weights,
        num_samples=2500,
        burnin=800,
        chains=4,
        seed=831,
    )
    np.testing.assert_allclose(
        probabilities[tree], reference["probabilities"], atol=0.03, rtol=0
    )
    assert fit.liability_marginals[tree].mean == pytest.approx(
        reference["mean"], abs=0.05
    )
    assert fit.liability_marginals[tree].variance == pytest.approx(
        reference["variance"], abs=0.06
    )
    assert fit.thresholds[1] == pytest.approx(reference["threshold_mean"], abs=0.06)


def test_fixed_ordinal_reference():
    from tests.threshold_posterior_support import fixed_ordinal_root_reference

    tree = tree_from("(A:1,B:1,C:1)R;")
    states = ("low", "mid", "high")
    observed = dict(zip(("A", "B", "C"), states, strict=True))
    reference = fixed_ordinal_root_reference()
    probabilities, fit = compute_threshold_marginals(
        tree,
        states,
        observed,
        likelihoods(tree, states, observed),
        thresholds="0,1",
        num_samples=2000,
        burnin=500,
        chains=4,
        seed=237,
    )
    np.testing.assert_allclose(
        probabilities[tree], reference["probabilities"], atol=0.025, rtol=0
    )
    assert fit.liability_marginals[tree].mean == pytest.approx(
        reference["mean"], abs=0.04
    )
    assert fit.liability_marginals[tree].variance == pytest.approx(
        reference["variance"], abs=0.05
    )
