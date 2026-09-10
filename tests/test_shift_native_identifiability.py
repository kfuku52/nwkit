import numpy as np
import pytest

from nwkit.shift_native_identifiability import (
    _correlation_derivative,
    covariance_identifiability,
)
from nwkit.shift_native_model import ShiftData, ShiftLayout, evaluate_trait
from nwkit.util import read_tree
from tests.test_shift_native_search import sample_data


@pytest.mark.parametrize("stationary", [False, True])
@pytest.mark.parametrize("alpha", [1e-5, 0.2, 3, 1000])
def test_covariance_derivative_matches_numeric_differences(alpha, stationary):
    heights = np.array([0.0, 0.2, 0.6, 0.9])

    def independent(a):
        covariance = np.exp(-2 * a * (1 - heights))
        if not stationary:
            covariance *= -np.expm1(-2 * a * heights) / -np.expm1(-2 * a)
        return covariance

    value, derivative = _correlation_derivative(alpha, heights, stationary)
    delta = 1e-4
    numerical = (
        independent(alpha * np.exp(delta)) - independent(alpha * np.exp(-delta))
    ) / (2 * delta)
    np.testing.assert_allclose(value, independent(alpha), atol=1e-12)
    np.testing.assert_allclose(derivative, numerical, rtol=1e-5, atol=1e-10)


def test_balanced_four_tip_tree_cannot_separate_three_free_covariance_parameters():
    data = ShiftData.build(
        read_tree("((a:1,b:1):1,(c:1,d:1):1);", "auto", True, quiet=True),
        np.array([[0.0], [1.0], [2.0], [0.5]]),
        ["x"],
    )
    fit = evaluate_trait(data, ShiftLayout.build(data.tree), 0, 0.7, 1, 0.2)
    result = covariance_identifiability(
        data, 0, fit, alpha_free=True, process_free=True, noise_free=True
    )
    assert result["rank"] == 2
    assert not result["variance_decomposition_supported"]
    assert not result["finite_alpha_supported"]
    fixed_alpha = covariance_identifiability(
        data, 0, fit, alpha_free=False, process_free=True, noise_free=True
    )
    assert fixed_alpha["variance_decomposition_supported"]


def test_multiple_mrca_depths_separate_covariances_but_independent_limit_does_not():
    data = sample_data()
    for alpha, expected in [(0.7, True), (1000, False), (float("inf"), False)]:
        fit = evaluate_trait(data, ShiftLayout.build(data.tree), 0, alpha, 1, 0.2)
        result = covariance_identifiability(
            data, 0, fit, alpha_free=True, process_free=True, noise_free=True
        )
        assert result["variance_decomposition_supported"] == expected
