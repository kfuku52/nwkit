import math

import numpy as np
import pytest
from scipy.optimize import differential_evolution, minimize_scalar

from nwkit.shift_native_fit import NativeFitOptions, fit_native_layout
from nwkit.shift_native_model import ShiftData, ShiftLayout
from nwkit.util import read_tree
from tests.test_shift_native_model import TREE, dense_reference


def dense_fit(data, layout, alpha, variance, noise, trait=0):
    design, covariance = dense_reference(
        data.tree,
        layout,
        alpha,
        variance,
        data.variances[:, trait] + noise,
        "OUfixedRoot",
    )
    mask = np.isfinite(data.values[:, trait])
    x, y = design[mask], data.values[mask, trait]
    factor = np.linalg.cholesky(covariance[np.ix_(mask, mask)])
    xw, yw = np.linalg.solve(factor, x), np.linalg.solve(factor, y)
    beta = np.linalg.lstsq(xw, yw, rcond=None)[0]
    rss = np.sum((yw - xw @ beta) ** 2)
    likelihood = -0.5 * (
        len(y) * np.log(2 * np.pi) + 2 * np.log(np.diag(factor)).sum() + rss
    )
    return likelihood, rss, beta


def test_no_error_alpha_profile_beats_dense_grid_with_same_model():
    y = np.random.default_rng(7391).normal(size=(8, 1))
    data = ShiftData.build(read_tree(TREE, "auto", True, quiet=True), y, ["x"])
    layout = ShiftLayout.build(data.tree, [1])
    result = fit_native_layout(data, layout)
    likelihoods = []
    for alpha in np.r_[0, np.geomspace(1e-6, 1000, 150), np.inf]:
        _, rss, _ = dense_fit(data, layout, alpha, 1, 0)
        likelihoods.append(dense_fit(data, layout, alpha, rss / len(y), 0)[0])
    assert result["fits"][0].log_likelihood >= max(likelihoods) - 1e-8
    assert result["traits"][0]["optimizer"]["complete_alpha_modes"]


@pytest.mark.parametrize("known_error", [False, True])
def test_estimated_noise_with_fixed_process_matches_dense_scalar_fit(known_error):
    y = np.random.default_rng(832).normal(size=(8, 1))
    errors = np.arange(8)[:, None] / 100 if known_error else None
    data = ShiftData.build(read_tree(TREE, "auto", True, quiet=True), y, ["x"], errors)
    layout = ShiftLayout.build(data.tree, [1])
    variance = 0.03
    result = fit_native_layout(
        data,
        layout,
        alpha_height=0.8,
        process_variance=variance * data.scales[0] ** 2,
        options=NativeFitOptions(estimate_measurement_error=True),
    )

    def objective(noise):
        return -dense_fit(data, layout, 0.8, variance, noise)[0]

    optimized = minimize_scalar(
        objective, bounds=(0, 10), method="bounded", options={"xatol": 1e-11}
    )
    expected = min(optimized.fun, objective(0))
    fit = result["fits"][0]
    assert fit.log_likelihood == pytest.approx(-expected, abs=2e-7)
    assert fit.process_variance == pytest.approx(variance)
    assert result["traits"][0]["covariance_components_identifiable"]


def test_independent_limit_does_not_invent_process_noise_decomposition():
    y = np.random.default_rng(933).normal(size=(8, 2))
    data = ShiftData.build(read_tree(TREE, "auto", True, quiet=True), y, ["x", "y"])
    layout = ShiftLayout.build(data.tree, [1])
    result = fit_native_layout(
        data,
        layout,
        alpha_height=math.inf,
        options=NativeFitOptions(estimate_measurement_error=True),
    )
    for record in result["traits"]:
        assert not record["covariance_components_identifiable"]
        assert record["process_tip_variance"] is None
        assert record["measurement_variance"] is None
        assert record["unstructured_tip_variance"] > 0
        assert record["alpha_status"] == "independent_limit"


def test_fixed_variances_use_original_trait_units_and_masked_observations():
    y = np.random.default_rng(999).normal(size=(8, 2)) * [1e-9, 3e4] + [7, 200]
    y[3, 1] = np.nan
    data = ShiftData.build(read_tree(TREE, "auto", True, quiet=True), y, ["x", "y"])
    layout = ShiftLayout.build(data.tree, [1])
    result = fit_native_layout(
        data,
        layout,
        alpha_height=[0.8, 1.3],
        process_variance=[2e-18, 4e8],
        measurement_variance=[1e-20, 1e7],
    )
    for j, record in enumerate(result["traits"]):
        likelihood, _, beta = dense_fit(
            data,
            layout,
            [0.8, 1.3][j],
            [2e-18, 4e8][j] / data.scales[j] ** 2,
            [1e-20, 1e7][j] / data.scales[j] ** 2,
            j,
        )
        assert record["log_likelihood"] == pytest.approx(
            likelihood - record["num_observations"] * np.log(data.scales[j]), abs=1e-8
        )
        np.testing.assert_allclose(result["fits"][j].coefficients, beta, atol=1e-8)


@pytest.mark.parametrize("missing", [False, True])
def test_joint_process_and_noise_fit_matches_independent_dense_optimization(missing):
    values = np.random.default_rng(813).normal(size=(8, 1))
    if missing:
        values[2, 0] = np.nan
    data = ShiftData.build(
        read_tree(TREE, "auto", True, quiet=True),
        values,
        ["x"],
        np.arange(1, 9)[:, None] / 100,
    )
    layout = ShiftLayout.build(data.tree, [1])
    result = fit_native_layout(
        data,
        layout,
        alpha_height=0.8,
        options=NativeFitOptions(estimate_measurement_error=True),
    )

    def objective(parameters):
        return -dense_fit(data, layout, 0.8, parameters[0], parameters[1])[0]

    dense = differential_evolution(
        objective, [(0, 5), (0, 5)], seed=81, tol=1e-10, polish=True
    )
    assert dense.success
    assert result["fits"][0].log_likelihood == pytest.approx(-dense.fun, abs=1e-6)


def test_known_error_can_support_exact_zero_process_and_extra_measurement_variance():
    data = ShiftData.build(
        read_tree(TREE, "auto", True, quiet=True),
        np.linspace(-0.001, 0.001, 8)[:, None],
        ["x"],
        np.ones((8, 1)),
    )
    layout = ShiftLayout.build(data.tree)
    result = fit_native_layout(
        data, layout, options=NativeFitOptions(estimate_measurement_error=True)
    )
    fit = result["fits"][0]
    assert fit.process_variance == 0
    assert fit.measurement_variance == 0
    assert result["traits"][0]["alpha"] is None
    assert result["traits"][0]["alpha_candidate_status"] == "brownian_limit"
    assert not result["traits"][0]["optimizer"]["nuisance_variance_at_numerical_bound"]
    assert fit.log_likelihood == pytest.approx(
        dense_fit(data, layout, 0, 0, 0)[0], abs=1e-12
    )


def test_estimated_noise_retains_structural_rank_reason_at_independent_limit():
    values = np.random.default_rng(52).normal(size=(8, 1))
    values[2, 0] = np.nan
    data = ShiftData.build(
        read_tree(TREE, "auto", True, quiet=True), values, ["x"], np.ones((8, 1)) * 0.01
    )
    layout = ShiftLayout.build(data.tree, [4, 10], [[0, 10], [4]])
    with pytest.raises(ValueError, match="rank deficient"):
        fit_native_layout(
            data,
            layout,
            alpha_height=math.inf,
            options=NativeFitOptions(estimate_measurement_error=True),
        )


def test_equivalent_alpha_modes_do_not_invent_a_variance_decomposition():
    data = ShiftData.build(
        read_tree(TREE, "auto", True, quiet=True),
        np.tile([-1.0, 1.0], 4)[:, None],
        ["x"],
    )
    result = fit_native_layout(
        data,
        ShiftLayout.build(data.tree),
        options=NativeFitOptions(estimate_measurement_error=True),
    )
    record = result["traits"][0]
    assert record["optimizer"]["equivalent_modes_disagree_on_variance_decomposition"]
    assert not record["covariance_components_identifiable"]
    assert record["process_tip_variance"] is None
    assert record["measurement_variance"] is None
    assert record["fitted_total_tip_variance"] == pytest.approx(1.0, abs=1e-6)
