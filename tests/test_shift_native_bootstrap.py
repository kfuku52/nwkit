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
