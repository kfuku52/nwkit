"""Regression cases from numerical and input-contract review of joint ASR."""

from dataclasses import replace

import numpy as np
import pytest
from scipy.optimize import OptimizeResult

from nwkit.individual_covariance import (
    _optimizer_converged,
    _restore_covariance,
    _validate_arrays,
    conditional_vector,
    fit_individual_covariance,
)
from tests import test_individual_asr as reference


@pytest.fixture(scope="module")
def problem():
    return reference.problem.__wrapped__()


@pytest.fixture(scope="module")
def fitted(problem):
    return reference.fit_problem(problem)


def fit_scaled(problem, multiplier):
    _, _, data, c = problem
    return fit_individual_covariance(
        c,
        data.values * multiplier,
        data.species,
        data.individuals,
        data.coordinates,
        dimension=2,
    )


def test_positive_covariance_cannot_silently_underflow(problem):
    with pytest.raises(ValueError, match="underflows to zero"):
        fit_scaled(problem, 1e-170)


def test_subnormal_covariance_is_retained_when_representable(problem, fitted):
    small = fit_scaled(problem, 1e-160)
    # Divide sequentially: squaring the conversion factor would underflow.
    assert small.sigma / 1e-160 / 1e-160 == pytest.approx(fitted.sigma, abs=0.001)
    assert small.within / 1e-160 / 1e-160 == pytest.approx(fitted.within, abs=0.001)
    assert (np.diag(small.within) > 0).all()


def test_covariance_restoration_avoids_intermediate_overflow_and_underflow():
    # Direct multiplication or division loses a finite, nonzero answer.
    first = _restore_covariance(np.array([[1e-200]]), np.array([1e200]), 1e200)
    second = _restore_covariance(np.array([[1e-200]]), np.array([1e-100]), 1e-300)
    assert first[0, 0] == pytest.approx(1e0)
    assert second[0, 0] == pytest.approx(1e-100)


@pytest.mark.parametrize("scale", [1e-150, 1, 1e150])
def test_asymmetry_validation_is_independent_of_units(problem, scale):
    _, _, data, c = problem
    c = c.copy() * scale
    c[0, 1] += 0.1 * scale
    with pytest.raises(ValueError, match="symmetric"):
        _validate_arrays(
            c, data.values, data.species, data.individuals, data.coordinates, 2
        )


def test_unobserved_large_variance_cannot_mask_indefinite_observed_block(problem):
    _, _, data, c = problem
    bad = np.zeros((9, 9))
    bad[:8, :8] = c
    bad[0, 1] = bad[1, 0] = 4
    bad[8, 8] = 1e20
    with pytest.raises(ValueError, match="positive semidefinite"):
        _validate_arrays(
            bad, data.values, data.species, data.individuals, data.coordinates, 2
        )


@pytest.mark.parametrize("dimension", [True, 2.0, -1])
def test_dimension_contract(problem, dimension):
    _, _, data, c = problem
    with pytest.raises(ValueError, match="at least two"):
        _validate_arrays(
            c, data.values, data.species, data.individuals, data.coordinates, dimension
        )


def test_optimizer_success_alone_is_not_convergence():
    result = OptimizeResult(success=True, x=np.array([0.0]), jac=np.array([0.1]))
    assert not _optimizer_converged(result, [(-1, 1)])
    result.jac[:] = np.nan
    assert not _optimizer_converged(result, [(-1, 1)])
    result.x[:] = -1
    result.jac[:] = 0.1
    assert _optimizer_converged(result, [(-1, 1)])


def test_exact_individual_coordinates_have_exact_zero_covariance(problem, fitted):
    _, _, data, c = problem
    mean, variance = conditional_vector(
        fitted,
        data.coordinates,
        c[0, data.species],
        c[0, 0],
        same_individual=data.individuals == 0,
    )
    assert variance == pytest.approx(np.zeros((2, 2)), abs=0)
    assert mean == pytest.approx(data.matrix[0], abs=1e-10)


def test_actual_prediction_variance_underflow_is_rejected(problem, fitted):
    _, _, data, c = problem
    tiny = replace(fitted, scale=fitted.scale * 1e-170)
    with pytest.raises(ValueError, match="underflows to zero"):
        conditional_vector(tiny, data.coordinates, c[0, data.species], c[0, 0])


def test_reordering_preserves_normalization_and_boundary_diagnostics(problem, fitted):
    _, _, data, c = problem
    reverse = np.arange(len(data.values))[::-1]
    other = fit_individual_covariance(
        c,
        data.values[reverse],
        data.species[reverse],
        data.individuals[reverse],
        data.coordinates[reverse],
        dimension=2,
    )
    assert np.array_equal(other.scale, fitted.scale)
    assert np.array_equal(other.offset, fitted.offset)
    assert other.fit_status == fitted.fit_status
    assert other.sigma_eigenvalue_ratio == pytest.approx(
        fitted.sigma_eigenvalue_ratio, abs=1e-6
    )


def test_empty_primary_output_is_rejected():
    from nwkit.cli import main

    with pytest.raises(ValueError, match="--outfile must be"):
        main(reference.arguments("--outfile", ""))


def test_sample_mean_se_retains_subnormal_individual_variance(problem, fitted):
    from nwkit.individual_asr import _species_summaries

    data = problem[2]
    variance = np.nextafter(0.0, 1.0)
    minimal = replace(fitted, within=np.eye(2) * variance)
    _, errors, _ = _species_summaries(data, minimal)
    # sqrt(W/n) underflows before the root although this SE is representable.
    assert errors["A"][0] > 0
    assert errors["A"][0] / np.sqrt(variance) == pytest.approx(1 / np.sqrt(3))
