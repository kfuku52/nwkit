import numpy as np
import pytest

from nwkit.shift_native_bootstrap import (
    calibrate_native_search,
    native_selection_support,
    simulate_native_data,
)
from nwkit.shift_native_fit import NativeFitOptions, fit_native_layout
from nwkit.shift_native_model import ShiftData, ShiftLayout
from nwkit.shift_native_search import exhaustive_native_search
from tests.test_shift_native_model import dense_reference
from tests.test_shift_native_search import sample_data


@pytest.mark.parametrize(
    "root,alpha",
    [
        ("OUfixedRoot", 0),
        ("OUfixedRoot", 1),
        ("OUrandomRoot", 1),
        ("OUrandomRoot", float("inf")),
    ],
)
def test_simulator_matches_independent_dense_covariance_and_original_units(root, alpha):
    base = sample_data()
    data = ShiftData.build(
        base.tree, base.values * 3 + 7, base.trait_names, np.full((8, 2), 0.12)
    )
    result = fit_native_layout(
        data,
        ShiftLayout.build(data.tree, [1]),
        alpha_height=alpha,
        process_variance=2,
        measurement_variance=0.3,
        options=NativeFitOptions(root_model=root),
    )
    rng = np.random.default_rng(192)
    samples = []
    for _ in range(4000):
        draw = simulate_native_data(data, result, rng)
        samples.append(draw.values * draw.scales + draw.centers)
    samples = np.asarray(samples)
    for j, fit in enumerate(result["fits"]):
        _, expected = dense_reference(
            data.tree,
            result["layout"],
            alpha,
            fit.process_variance,
            data.variances[:, j] + fit.measurement_variance,
            root,
        )
        expected *= data.scales[j] ** 2
        np.testing.assert_allclose(
            samples[:, :, j].mean(axis=0),
            fit.predicted * data.scales[j] + data.centers[j],
            atol=0.09,
        )
        np.testing.assert_allclose(
            np.cov(samples[:, :, j], rowvar=False), expected, atol=0.17
        )


def test_simulation_preserves_missingness_and_known_error_units():
    base = sample_data()
    values = base.values.copy()
    values[2, 1] = np.nan
    data = ShiftData.build(base.tree, values, base.trait_names, np.full((8, 2), 0.13))
    fit = fit_native_layout(
        data, ShiftLayout.build(data.tree), alpha_height=0.7, process_variance=1
    )
    simulated = simulate_native_data(data, fit, np.random.default_rng(1))
    np.testing.assert_array_equal(np.isnan(simulated.values), np.isnan(data.values))
    np.testing.assert_allclose(simulated.variances * simulated.scales**2, 0.13)


def test_calibration_replays_search_deterministically_and_does_not_discard_failures():
    data = sample_data()
    calls = []

    def run(sample):
        calls.append(sample)
        return exhaustive_native_search(
            sample,
            max_shifts=1,
            fit_arguments={"alpha_height": 0.7, "process_variance": 1},
        )

    search = run(data)
    first, metadata = calibrate_native_search(data, search, run, replicates=19, seed=18)
    assert len(calls) == 20
    second, repeated = calibrate_native_search(
        data, search, run, replicates=19, seed=18
    )
    assert first["layout"] == second["layout"]
    assert metadata == repeated
    assert (
        metadata["tests"][0]["p_value"]
        == (1 + metadata["tests"][0]["exceedances"]) / 20
    )

    def fail(_):
        raise ValueError("optimizer failure")

    with pytest.raises(ValueError, match="draw 0; no draws discarded"):
        calibrate_native_search(data, search, fail, replicates=19, seed=18)
    with pytest.raises(ValueError, match="cannot resolve"):
        calibrate_native_search(data, search, run, replicates=18, seed=18)


def test_support_replays_selected_procedure_with_independent_seeds():
    data = sample_data()
    result = fit_native_layout(
        data, ShiftLayout.build(data.tree, [1]), alpha_height=0.7, process_variance=1
    )
    seeds = []

    def select(sample, seed):
        seeds.append(seed)
        return fit_native_layout(
            sample, result["layout"], alpha_height=0.7, process_variance=1
        )

    support = native_selection_support(data, result, select, replicates=3, seed=9)
    assert len(set(seeds)) == 3
    assert support["exact_layout_frequency"] == 1
    assert support["regime_pair_frequencies"] == [
        {
            "first_branch_id": 0,
            "second_branch_id": 1,
            "cooccurrences": 3,
            "shared_regime_count": 0,
            "shared_regime_frequency": 0,
            "shared_given_both_selected": 0,
        }
    ]
    assert (
        next(row for row in support["branch_frequencies"] if row["branch_id"] == 1)[
            "frequency"
        ]
        == 1
    )


def test_global_null_aic_gate_replays_search_and_preserves_ungated_winner():
    from nwkit.shift_native_bootstrap import gate_native_aic

    data = sample_data()
    calls = []

    def run(sample):
        calls.append(sample)
        return exhaustive_native_search(
            sample,
            max_shifts=1,
            fit_arguments={"alpha_height": 0.7, "process_variance": 1},
            criterion="AIC",
        )

    search = run(data)
    selected, metadata = gate_native_aic(data, search, run, replicates=19, seed=18)
    assert len(calls) == 20
    assert metadata["p_value"] == (1 + metadata["exceedances"]) / 20
    assert metadata["statistic"] == pytest.approx(
        max(
            0,
            search.best_by_complexity[0]["information_criterion"]["score"]
            - search.best_information["information_criterion"]["score"],
        )
    )
    expected = (
        search.best_information
        if metadata["rejected"]
        else search.best_by_complexity[0]
    )
    assert selected is expected
    _, repeated = gate_native_aic(data, search, run, replicates=19, seed=18)
    assert repeated == metadata
    assert not metadata["controls_false_branches_under_nonnull"]


def test_global_null_gate_ties_failures_and_resolution():
    from nwkit.shift_native_bootstrap import gate_native_aic

    data = sample_data()
    search = exhaustive_native_search(data, max_shifts=0, criterion="AIC")
    selected, metadata = gate_native_aic(
        data,
        search,
        lambda _: search,
        replicates=19,
        seed=5,
    )
    assert metadata["statistic"] == 0
    assert metadata["p_value"] == 1
    assert not selected["layout"].shifts
    for draws in [0, 18, True, 19.5]:
        with pytest.raises(ValueError):
            gate_native_aic(data, search, lambda _: search, replicates=draws, seed=5)

    def fail(_):
        raise ValueError("optimizer failure")

    with pytest.raises(ValueError, match="draw 0; no draws discarded"):
        gate_native_aic(data, search, fail, replicates=19, seed=5)


def test_global_null_gate_rejects_large_gain_and_accepts_small_gain(monkeypatch):
    from copy import deepcopy

    from nwkit.shift_native_bootstrap import gate_native_aic

    data = sample_data()
    observed = exhaustive_native_search(data, max_shifts=1, criterion="AIC")
    observed.best_information = deepcopy(observed.best_information)
    observed.best_information["information_criterion"]["score"] = (
        observed.best_by_complexity[0]["information_criterion"]["score"] - 10
    )
    no_gain = deepcopy(observed)
    no_gain.best_information = no_gain.best_by_complexity[0]
    selected, metadata = gate_native_aic(
        data, observed, lambda _: no_gain, replicates=19, seed=5
    )
    assert metadata["p_value"] == 0.05
    assert metadata["rejected"]
    assert selected is observed.best_information
    selected, metadata = gate_native_aic(
        data, no_gain, lambda _: observed, replicates=19, seed=5
    )
    assert metadata["p_value"] == 1
    assert selected is no_gain.best_by_complexity[0]
