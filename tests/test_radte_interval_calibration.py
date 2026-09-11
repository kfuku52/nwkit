import copy
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.optimize import Bounds, LinearConstraint, OptimizeResult
from scipy.stats import t

from nwkit.radte_calibrated import calibrated_profile_intervals, calibrated_profile_test
from nwkit.radte_marginal import MarginalDatingProblem
from nwkit.radte_model import _checked_zero_variance_boundary, fit_dates
from nwkit.radte_sequence import SequenceLikelihood
from nwkit.radte_simulation import simulate_rate_lengths, simulate_sequence_likelihood
from tests.test_radte import small_chronology
from tests.test_radte_studentized import quadratic_likelihood


@pytest.mark.parametrize("rho", [0.0, 0.5, 0.9])
def test_zero_variance_right_derivative_and_continuity(rho):
    c = small_chronology()
    q = quadratic_likelihood(c)
    p = MarginalDatingProblem(c, likelihood=q, rho=rho)
    x = p.initial_parameters()
    x[-1] = 0
    value, gradient = p.value_gradient(x)
    for step in [1e-6, 1e-7]:
        plus, twice = x.copy(), x.copy()
        plus[-1], twice[-1] = step, 2 * step
        # Second-order forward derivative; no native gradient in the reference.
        derivative = (
            -3 * value + 4 * p.value_gradient(plus)[0] - p.value_gradient(twice)[0]
        ) / (2 * step)
        assert gradient[-1] == pytest.approx(derivative, rel=1e-5, abs=1e-5)
        assert p.value_gradient(plus)[1][-1] == pytest.approx(
            gradient[-1], rel=2e-3, abs=2e-3
        )


def test_estimated_marginal_variance_can_reach_the_true_zero_boundary():
    c = small_chronology()
    q = quadratic_likelihood(c)
    q.center[:] = np.log(np.bincount(q.mapping, weights=c.durations(c.initial) * 0.1))
    fit, problem = fit_dates(c, likelihood=q, inference="marginal")
    assert fit.log_rate_sd == 0
    assert fit.parameters[-1] == 0
    assert problem.value_gradient(fit.parameters)[1][-1] >= -1e-6
    assert "estimated_rate_variance_at_zero_boundary" in fit.diagnostics
    assert "nuisance_parameter_at_numerical_bound" not in fit.diagnostics
    from nwkit.radte_studentized import studentized_intervals

    studentized_intervals(fit, problem)
    assert fit.interval_status == "unavailable-estimated-zero-rate-variance"


@pytest.mark.parametrize("true_variance", [0.0, 1e-8])
def test_zero_boundary_refit_preserves_genuine_small_positive_variance(true_variance):
    def objective(x):
        residual = x - np.array([0.5, true_variance])
        return float(residual @ residual), 2 * residual

    point = np.array([0.5, true_variance or 8e-17])
    best = OptimizeResult(x=point, fun=objective(point)[0], success=True)
    problem = SimpleNamespace(
        marginal=True,
        fixed_sd=None,
        value_gradient=objective,
        feasible=lambda x: x[-1] >= 0,
    )
    attempts = []
    result = _checked_zero_variance_boundary(
        problem,
        best,
        Bounds([-30.0, 0.0], [30.0, 50.0]),
        LinearConstraint(np.empty((0, 2)), [], []),
        100,
        attempts,
    )
    assert result.x[-1] == true_variance
    assert attempts[-1]["success"] == (true_variance == 0)


def test_dense_reference_matches_independent_two_group_t_formula():
    path = Path(__file__).resolve().parents[1] / "tools/radte_interval_reference.py"
    spec = importlib.util.spec_from_file_location("interval_reference", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    y = np.array([0.4, 0.8, -0.2, 0.3, -0.5, 0.0])
    x = np.column_stack([np.ones(6), [1, 1, 0, 0, 0, 0]])
    answer = module.gaussian_contrast_interval(y, x, np.eye(6), [0, 1])
    center = y[:2].mean() - y[2:].mean()
    variance = (
        ((y[:2] - y[:2].mean()) ** 2).sum() + ((y[2:] - y[2:].mean()) ** 2).sum()
    ) / 4
    width = t.ppf(0.975, 4) * np.sqrt(variance * (1 / 2 + 1 / 4))
    assert answer["lower"] == pytest.approx(center - width)
    assert answer["upper"] == pytest.approx(center + width)


def test_parametric_rates_have_genealogical_covariance():
    c = small_chronology()
    rng = np.random.default_rng(92410001)
    logs = np.array(
        [
            np.log(
                simulate_rate_lengths(c, c.initial, 0, 0.4, 0.5, rng)
                / c.durations(c.initial)
            )
            for _ in range(4000)
        ]
    )
    # The two root edges share an unobserved parent: correlation rho^2.
    assert np.cov(logs.T)[0, 1] == pytest.approx(0.4**2 * 0.5**2, abs=0.008)
    np.testing.assert_allclose(logs.var(axis=0), 0.4**2, atol=0.012)


def test_parametric_alignment_matches_jc_pair_probability():
    c = small_chronology()
    matrix = np.ones((4, 12000), dtype=np.uint64)
    exact = SequenceLikelihood(c, None, matrix=matrix, model="jc69", gamma_categories=1)
    lengths = np.full(6, 0.1)
    simulated = simulate_sequence_likelihood(
        exact, lengths, np.random.default_rng(92410002)
    )
    # The first two tips are sisters separated by .2 substitutions/site.
    expected = 0.25 + 0.75 * np.exp(-4 * 0.2 / 3)
    assert np.mean(simulated.raw_matrix[0] == simulated.raw_matrix[1]) == pytest.approx(
        expected, abs=0.012
    )
    np.testing.assert_array_equal(exact.raw_matrix, matrix)


def test_parametric_generator_does_not_reuse_informative_ambiguity_as_missing():
    c = small_chronology()
    matrix = np.ones((4, 10), dtype=np.uint64)
    matrix[0, 0] = 3  # A/C is not an uninformative missing observation.
    exact = SequenceLikelihood(c, None, matrix=matrix, model="jc69", gamma_categories=1)
    with pytest.raises(ValueError, match="observation model"):
        simulate_sequence_likelihood(exact, np.full(6, 0.1), np.random.default_rng(1))


def test_calibrated_lr_matches_independent_branch_profile_and_preserves_input():
    c = small_chronology(max_age=100)
    fit, problem = fit_dates(c)
    before = copy.deepcopy(fit)
    original = [node.dist for node in c.edges]
    group = problem.free[0]
    age = 3.0
    row = calibrated_profile_test(
        fit, problem, group, age, replicates=19, starts=1, seed=92410003
    )
    observed = np.array(original)
    durations = np.array([2, 2, 1, 1, 1, 1])
    y = np.log(observed / durations)
    null_sse = ((y - y.mean()) ** 2).sum()
    residual = np.log(observed)
    alt_sse = ((residual[:2] - residual[:2].mean()) ** 2).sum() + (
        (residual[2:] - residual[2:].mean()) ** 2
    ).sum()
    assert row["statistic"] == pytest.approx(6 * np.log(null_sse / alt_sse), abs=1e-7)
    assert row["failures"] == 0
    assert row["p_upper"] >= row["exceedances"] / row["replicates"]
    np.testing.assert_array_equal(fit.parameters, before.parameters)
    np.testing.assert_array_equal([node.dist for node in c.edges], original)


def test_calibrated_profile_rejects_joint_map_instead_of_using_penalty_lr():
    c = small_chronology()
    fit, problem = fit_dates(
        c, likelihood=quadratic_likelihood(c), inference="joint-map"
    )
    with pytest.raises(ValueError, match="marginal inference"):
        calibrated_profile_test(
            fit, problem, problem.free[0], fit.ages[problem.free[0]], replicates=19
        )


def test_calibrated_failure_clears_stale_intervals(monkeypatch):
    import nwkit.radte_calibrated as module

    c = small_chronology()
    fit, problem = fit_dates(
        c, likelihood=quadratic_likelihood(c), inference="marginal"
    )
    fit.interval_lower, fit.interval_upper = fit.ages.copy(), fit.ages.copy()

    def fail(*args, **kwargs):
        raise ValueError("reference failure")

    monkeypatch.setattr(module, "calibrated_profile_test", fail)
    calibrated_profile_intervals(fit, problem, replicates=99, grid_points=5)
    assert fit.interval_lower is None and fit.interval_upper is None
    assert fit.interval_status == "unavailable-calibrated-profile-fit-failure"


@pytest.mark.parametrize("rho", [0.0, 0.5, 0.9])
@pytest.mark.parametrize("sd", [None, 0.3])
def test_exact_log_duration_interval_matches_whitened_reference(rho, sd):
    from nwkit.radte_exact_interval import exact_log_duration_intervals

    path = Path(__file__).resolve().parents[1] / "tools/radte_interval_reference.py"
    spec = importlib.util.spec_from_file_location("interval_reference", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    c = small_chronology(max_age=100)
    fit, problem = fit_dates(c, rho=rho, rate_sd=sd)
    x = np.column_stack([np.ones(6), [1, 1, 0, 0, 0, 0]])
    # Explicit distances on this six-edge genealogy, independent of precision.
    distances = np.array(
        [
            [0, 2, 1, 1, 3, 3],
            [2, 0, 3, 3, 1, 1],
            [1, 3, 0, 2, 4, 4],
            [1, 3, 2, 0, 4, 4],
            [3, 1, 4, 4, 0, 2],
            [3, 1, 4, 4, 2, 0],
        ]
    )
    reference = module.gaussian_contrast_interval(
        np.log([n.dist for n in c.edges]), x, rho**distances, [0, 1], sd=sd
    )
    assert exact_log_duration_intervals(fit, problem)
    group = problem.free[0]
    assert fit.interval_lower[group] == pytest.approx(
        max(1 + c.min_duration, 1 + np.exp(reference["lower"]))
    )
    assert fit.interval_upper[group] == pytest.approx(
        min(10, 1 + np.exp(reference["upper"]))
    )


def test_exact_interval_does_not_turn_strict_clock_into_a_point_interval():
    from nwkit.radte_exact_interval import exact_log_duration_intervals

    c = small_chronology(gene_text="((A_1:.1,B_1:.1)S1:.1,(A_2:.1,B_2:.1)S2:.1)D;")
    fit, problem = fit_dates(c)
    assert fit.log_rate_sd == 0
    assert exact_log_duration_intervals(fit, problem)
    assert fit.interval_status == "unavailable-strict-clock-limit"
    assert fit.interval_lower is None
