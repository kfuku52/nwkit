"""Independent density, Monte Carlo and failure checks for the research candidate."""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from shift_confidence_calibration import (  # noqa: E402
    NullConfidenceBank,
    restricted_envelope,
)
from shift_continuous_reference import DenseOUReference  # noqa: E402
from try_shift_confidence_calibration import summarize  # noqa: E402
from validate_shift_null_contract import tree_text  # noqa: E402

from nwkit.shift_calibration import CalibratedSearch  # noqa: E402
from nwkit.util import read_tree  # noqa: E402


@pytest.fixture(scope="module")
def bank():
    tree = read_tree(tree_text(4, "pectinate"), "auto", True, quiet=True)
    search = CalibratedSearch(tree, convergence=False, alpha_grid=[0, 0.4, np.inf])
    result = NullConfidenceBank(search)
    result.simulate(seed=149, replicates=199)
    return result


def test_confidence_profile_matches_independent_dense_gaussian(bank):
    values = np.array([0.1, -0.3, 1.2, 0.7])
    actual = bank.null_profile(bank.search.q @ values)[:, 0]
    ref = DenseOUReference(bank.search.tree, bank.search.models[bank.null_index])
    expected = [ref.at(values, a)["log_likelihood"] for a in bank.search.grid]
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-10)


def test_every_full_search_mc_probability_matches_native_replay(bank):
    values = np.array([0.1, -0.3, 1.2, 0.7])
    result = bank.evaluate(values, betas=(0, 0.005, 0.01))
    noise = np.random.default_rng(149).normal(size=(bank.search.d, 199))
    for i, item in enumerate(bank.search.cache):
        p = bank.search._probability_from_noise(
            np.zeros(bank.search.d),
            item[2],
            1.0,
            bank.search.families[0],
            result["shift_statistic"],
            noise,
        )
        assert p == result["shift_probabilities"][i]
    assert result["candidates"][0]["p_value"] == result["full_envelope_p"]
    native = bank.search.fit(values, seed=149, replicates=199)
    assert (native["tests"][0]["p_value"] <= 0.05) == result["full_envelope_reject"]


def test_confidence_statistic_is_independently_reconstructed(bank):
    values = np.array([0.1, -0.3, 1.2, 0.7])
    result = bank.evaluate(values)
    reference = DenseOUReference(bank.search.tree, bank.search.models[bank.null_index])
    ll = np.array([reference.at(values, a)["log_likelihood"] for a in bank.search.grid])
    np.testing.assert_allclose(
        result["confidence_statistics"], 2 * (ll.max() - ll), atol=1e-10
    )
    assert all(0 < p <= 1 for p in result["confidence_probabilities"])


def test_scale_and_constant_offset_do_not_change_confidence_decisions(bank):
    values = np.array([0.1, -0.3, 1.2, 0.7])
    before, after = bank.evaluate(values), bank.evaluate(7 * values + 31)
    assert before["shift_probabilities"] == after["shift_probabilities"]
    assert before["confidence_probabilities"] == after["confidence_probabilities"]
    assert before["candidates"] == after["candidates"]


def test_confidence_exclusion_is_penalized_and_empty_set_is_explicit():
    result = restricted_envelope([0.9, 0.02], [0.005, 0.8], 0.01)
    assert result["retained_indices"] == [1]
    assert result["p_value"] == pytest.approx(0.03)
    assert result["reject"]
    empty = restricted_envelope([0.9, 0.2], [0.005, 0.005], 0.01)
    assert empty["empty_confidence_set"]
    assert empty["p_value"] == 0.01


@pytest.mark.parametrize(
    "probabilities,beta", [([np.nan], 0.01), ([0], 0.01), ([1.1], 0.01), ([0.1], 0.05)]
)
def test_invalid_confidence_inputs_are_not_silently_accepted(probabilities, beta):
    with pytest.raises(ValueError):
        restricted_envelope(probabilities, [0.5], beta)


def test_known_error_is_outside_this_prototype(bank):
    search = CalibratedSearch(bank.search.tree, variances=np.ones(4), alpha_grid=[0.4])
    with pytest.raises(ValueError, match="no known errors"):
        NullConfidenceBank(search)


def test_failures_remain_in_futility_bound_and_prevent_advancement():
    result = summarize([dict(status="failed", error="numerical")], "development")
    assert result["failed_blocks"] == 1
    assert result["candidates_passing_screen"] == []
    assert all(c["gross_gain_simultaneous_95_upper"] == 1 for c in result["cells"])
    json.dumps(result, allow_nan=False)


def test_unrestricted_first_stage_contains_shared_candidates():
    tree = read_tree(tree_text(8, "balanced"), "auto", True, quiet=True)
    no = CalibratedSearch(tree, convergence=False)
    yes = CalibratedSearch(tree, convergence=True)
    values = np.random.default_rng(251).normal(size=(7, 6))
    a, _ = no.profile(values)
    b, _ = yes.profile(values)
    np.testing.assert_allclose(a.max(axis=0) - a[0], b.max(axis=0) - b[0], atol=1e-9)


@pytest.mark.parametrize(
    "field", ["full_log_likelihood", "shift_exceedances", "confidence_probabilities"]
)
def test_independent_final_audit_rejects_corrupted_evidence(field):
    import copy
    import gzip
    import json

    from verify_shift_confidence_calibration import check_condition

    with gzip.open("examples/shift/confidence-pilot/records.jsonl.gz", "rt") as stream:
        row = json.loads(next(stream))
    condition = copy.deepcopy(row["conditions"][0])
    value = condition["result"][field]
    if isinstance(value, list):
        value[0] += 1
    else:
        condition["result"][field] += 1
    tree = read_tree(row["tree"], "auto", True, quiet=True)
    search = CalibratedSearch(tree)
    with pytest.raises(ValueError):
        check_condition(tree, search, condition, row["data_seed"], 999)
