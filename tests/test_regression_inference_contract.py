"""Independent small references for the scientific regression inference contract."""

from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np
import pytest
from scipy import sparse
from scipy.integrate import quad
from scipy.stats import chi2, multivariate_normal, norm

from nwkit.event_regression import (
    event_average_estimate,
    event_average_fit,
    resolve_event_weighting,
)
from nwkit.gaussian import (
    DiagonalLowRankCovariance,
    DiagonalSparsePrecisionCovariance,
    project_covariance,
)
from nwkit.model_matrix import ReplicatedObservation
from nwkit.phylogenetic_glmm import (
    _draw_bootstrap_dataset,
    fit_phylogenetic_glmm,
    summarize_glmm_coefficient,
)
from nwkit.regress import _build_covariance_components, _profile_covariance_fit
from nwkit.regression_inference import (
    BootstrapTest,
    bootstrap_test,
    invert_bootstrap_grid,
    objective_difference,
)
from nwkit.regression_observation import CensoringModel


def test_event_average_targets_equal_event_loss_not_covariance_weighted_loss():
    x = np.ones((11, 1))
    y = np.r_[np.ones(10), 3.0]
    groups = np.r_[np.zeros(10), 1]
    beta, operator = event_average_estimate(y, x, groups)
    assert beta[0] == pytest.approx(2.0)
    np.testing.assert_allclose(operator[0], np.r_[np.full(10, 0.05), 0.5])
    # Independent common-coefficient MLE under equal independent noise.
    assert float(np.mean(y)) == pytest.approx(13 / 11)


@pytest.mark.parametrize(
    "representation",
    ["dense", "diagonal-low-rank", "sparse-low-rank", "sparse-precision"],
)
def test_event_sampling_covariance_matches_independent_dense_reference(representation):
    x = np.array([[1.0, -1.0], [1.0, 0.0], [1.0, 1.0], [1.0, 2.0]])
    groups = np.array([0, 0, 1, 2])
    y = np.array([0.2, 0.4, 1.1, 2.4])
    diagonal = np.array([0.2, 0.4, 0.3, 0.1])
    loading = np.array([[1.0, 0.0], [1.0, 0.0], [0.0, 0.8], [0.7, 0.8]])
    dense = np.diag(diagonal) + loading @ loading.T
    covariance = {
        "dense": dense,
        "diagonal-low-rank": DiagonalLowRankCovariance(diagonal, loading),
        "sparse-low-rank": DiagonalLowRankCovariance(
            diagonal, sparse.csr_matrix(loading)
        ),
        "sparse-precision": DiagonalSparsePrecisionCovariance(
            diagonal, sparse.csr_matrix(loading), sparse.eye(2, format="csc")
        ),
    }[representation]
    nuisance = {"covariance": covariance, "cholesky": np.linalg.cholesky(dense)}
    result = event_average_fit(y, x, groups, nuisance)
    w = np.diag([0.5, 0.5, 1.0, 1.0])
    gram_inverse = np.linalg.inv(x.T @ w @ x)
    reference = gram_inverse @ x.T @ w @ dense @ w @ x @ gram_inverse
    np.testing.assert_allclose(result["beta_covariance"], reference, atol=1e-12)
    np.testing.assert_allclose(result["beta"], gram_inverse @ x.T @ w @ y, atol=1e-12)


def test_projection_diagonal_and_eiv_use_covariance_at_event_coefficient():
    x = np.ones((3, 1))
    y = np.array([0.0, 0.0, 6.0])
    diagonal = np.ones(3)
    fit = {
        "covariance": diagonal,
        "cholesky": diagonal,
        "covariance_for_beta": lambda beta: (
            diagonal * (1 + beta[0] ** 2),
            np.sqrt(diagonal * (1 + beta[0] ** 2)),
        ),
    }
    result = event_average_fit(y, x, [0, 0, 1], fit)
    assert result["beta"][0] == pytest.approx(3)
    assert result["beta_covariance"][0, 0] == pytest.approx(
        10 * (0.25**2 + 0.25**2 + 0.5**2)
    )
    np.testing.assert_allclose(project_covariance(diagonal, np.eye(3)), np.eye(3))


def test_gaussian_likelihood_uses_full_normalizing_constant_and_actual_n():
    x = np.ones((8, 1))
    y = np.array([-2.0, 1.0, 3.0, -1.0, 2.0, 0.0, 1.0, 4.0])
    result = _profile_covariance_fit(
        y, x, np.zeros(8), [("rate", np.ones(8))], reml=False
    )
    rate = float(np.mean((y - y.mean()) ** 2))
    assert result["component_variances"]["rate"] == pytest.approx(rate)
    assert result["objective"] == pytest.approx(
        -multivariate_normal.logpdf(y, mean=np.full(8, y.mean()), cov=rate * np.eye(8))
    )
    with pytest.raises(ValueError, match="actual observation count"):
        _profile_covariance_fit(
            y,
            x,
            np.zeros(8),
            [("rate", np.ones(8))],
            reml=False,
            likelihood_observations=4,
        )


def test_shared_event_population_example_recovers_biological_variances():
    # Four event vectors have exact second moment [[2, 1], [1, 2]]. Repeating
    # them integrates this Gaussian log likelihood exactly, without simulation.
    chol = np.linalg.cholesky(np.array([[2.0, 1.0], [1.0, 2.0]]))
    event_vectors = np.concatenate([np.sqrt(2.0) * chol.T, -np.sqrt(2.0) * chol.T])
    response = np.tile(event_vectors, (5, 1)).reshape(-1)
    groups = np.repeat(np.arange(20), 2)
    x = np.empty((40, 0))
    fixed, components, factors, *_ = _build_covariance_components(
        x,
        np.ones(40),
        np.zeros(40),
        groups,
        np.full(20, 2),
        np.arange(40),
        np.ones(40, dtype=int),
        np.arange(20),
        np.arange(40),
        [],
        n_events=20,
        num_parameters=0,
        event_weighting="event",
        model="hierarchical",
        event_random_effect="yes",
        lineage_random_slope="no",
    )
    np.testing.assert_array_equal(components[0][1], np.ones(40))
    fit = _profile_covariance_fit(
        response, x, fixed, components, reml=False, component_factors=factors
    )
    assert fit["component_variances"]["evolutionary_rate"] == pytest.approx(
        1.0, abs=1e-5
    )
    assert fit["component_variances"]["species_event_variance"] == pytest.approx(
        1.0, abs=1e-5
    )


@pytest.mark.parametrize("reml, expected", [(False, 0.95), (True, 1.0)])
def test_bootstrap_rate_expectation_matches_finite_sample_ml_or_reml(reml, expected):
    # Deterministic quadrature for a quadratic estimator: columns sqrt(n)*I
    # have mean outer product I, so their average RSS is the exact expectation.
    n = 20
    x = np.ones((n, 1))
    draws = np.sqrt(n) * np.eye(n)
    rates = [
        _profile_covariance_fit(y, x, np.zeros(n), [("rate", np.ones(n))], reml=reml)[
            "component_variances"
        ]["rate"]
        for y in draws
    ]
    assert np.mean(rates) == pytest.approx(expected, abs=1e-12)


def test_aliases_default_to_event_average_and_reject_conflicts():
    assert resolve_event_weighting() == "event"
    assert resolve_event_weighting(regression_estimand="common") == "contrast"
    assert resolve_event_weighting("event", "event-average") == "event"
    with pytest.raises(ValueError, match="conflicts"):
        resolve_event_weighting("event", "common")


def test_detection_limits_regenerate_exact_and_censored_labels():
    template = [np.nan, 0.5, 0.4]
    model = CensoringModel("detection-limits", [0.0, 0.0, 0.0], [1.0, 1.0, 1.0])
    model.validate_observed(template, [np.nan] * 3, [0.0, np.nan, np.nan])
    values, lower, upper = model.observe([0.2, -1.0, 1.4], template)
    np.testing.assert_allclose(values, [0.2, np.nan, np.nan], equal_nan=True)
    np.testing.assert_allclose(lower, [np.nan, np.nan, 1.0], equal_nan=True)
    np.testing.assert_allclose(upper, [np.nan, 0.0, np.nan], equal_nan=True)
    model.validate_observed(values, lower, upper)
    with pytest.raises(ValueError, match="inconsistent"):
        model.validate_observed([-1.0, 0.5, 0.4], None, None)


def test_interval_bins_cover_all_values_and_preserve_biological_replicates():
    template = [ReplicatedObservation((np.nan, np.nan, np.nan, np.nan))]
    model = CensoringModel("interval-bins", cutpoints=(-1.0, 1.0))
    values, lower, upper = model.observe([-2.0, 0.0, 1.0, 2.0], template)
    assert isinstance(values[0], ReplicatedObservation)
    np.testing.assert_allclose(
        lower[0].values, [np.nan, -1.0, -1.0, 1.0], equal_nan=True
    )
    np.testing.assert_allclose(
        upper[0].values, [-1.0, 1.0, 1.0, np.nan], equal_nan=True
    )
    model.validate_observed(values, lower, upper)
    with pytest.raises(ValueError, match="strictly increasing"):
        CensoringModel("interval-bins", cutpoints=(1.0, -1.0)).observe([0.0], [0.0])


def test_censoring_expected_score_cancels_density_and_probability_parts():
    score = (
        norm.cdf(0) * (-norm.pdf(0) / norm.cdf(0))
        + quad(lambda y: y * norm.pdf(y), 0.0, np.inf)[0]
    )
    assert score == pytest.approx(0.0, abs=1e-12)


def test_glmm_bootstrap_passes_fresh_bounds_to_refit():
    @dataclass(frozen=True)
    class Options:
        family: str = "censored-gaussian"
        censoring_model: object = None
        censor_lower: object = None
        censor_upper: object = None

    n = 100
    template = np.r_[np.full(50, np.nan), np.ones(50)].tolist()
    options = Options(
        censoring_model=CensoringModel("detection-limits", np.zeros(n)),
        censor_upper=np.r_[np.zeros(50), np.full(50, np.nan)],
    )
    fit = SimpleNamespace(coefficients=np.zeros((1, 1)), dispersion=1.0)
    values, generated = _draw_bootstrap_dataset(
        np.random.default_rng(93), template, np.ones((n, 1)), fit, np.zeros(n), options
    )
    assert np.isfinite(values[:50]).any()
    assert np.isnan(values[50:]).any()
    exact = np.isfinite(values)
    assert np.all(np.asarray(values)[exact] > 0)
    options.censoring_model.validate_observed(
        values, generated.censor_lower, generated.censor_upper
    )


def test_penalized_gaussian_statistic_is_calibrated_under_its_own_null():
    # Independent exact null: for precision-3 regularization the objective
    # difference is Y^2/4, so its tail at Y=2 is P(chi-square_1 >= 4).
    result = bootstrap_test(
        1.0, lambda rng: rng.normal() ** 2 / 4, replicates=4000, seed=51
    )
    assert abs(result.p_value - chi2.sf(4.0, 1)) < 4 * result.monte_carlo_se
    assert abs(result.p_value - chi2.sf(1.0, 1)) > 0.2


def test_failed_null_dataset_is_not_replaced_and_bad_optimization_is_rejected():
    calls = []

    def simulate(rng):
        calls.append(1)
        raise ValueError("failed refit")

    with pytest.raises(RuntimeError, match="dataset 1/10"):
        bootstrap_test(1.0, simulate, replicates=10, seed=1)
    assert len(calls) == 1
    with pytest.raises(RuntimeError, match="worse"):
        objective_difference(3.0, 1.0)


def test_profile_grid_retains_disconnected_acceptance_without_inventing_intervals():
    evaluated = []

    def test(candidate):
        evaluated.append(candidate)
        return BootstrapTest(1.0, 0.9 if abs(candidate) == 1 else 0.01, 0.005, 1000)

    result = invert_bootstrap_grid(test, [-2.0, -1.0, 0.0, 1.0, 2.0], 0.95)
    assert evaluated == [-2.0, -1.0, 0.0, 1.0, 2.0]
    assert [row["null_value"] for row in result if row["accepted"]] == [-1.0, 1.0]
    assert all("confidence_interval_lower" not in row for row in result)


def test_penalized_wald_does_not_report_frequentist_p_or_intervals():
    fit = SimpleNamespace(
        coefficient_inference="wald",
        coefficient_penalty="gaussian",
        coefficient_statistics=None,
        coefficient_confidence_lower=None,
        coefficient_confidence_upper=None,
    )
    assert summarize_glmm_coefficient(fit, 0, 0.5, 1.96) == (
        "",
        "",
        "",
        "",
        "",
        "penalized-point-estimate",
    )
    with pytest.raises(ValueError, match="Penalized coefficients require"):
        fit_phylogenetic_glmm(
            [0.0, 1.0, 2.0],
            np.ones((3, 1)),
            lambda _: np.eye(3),
            family="poisson",
            inference="likelihood-ratio",
        )


def test_failed_profile_endpoint_is_not_replaced_with_a_wald_interval():
    fit = SimpleNamespace(
        coefficient_inference="profile-likelihood",
        coefficient_penalty="none",
        coefficient_covariance=np.ones((1, 1)),
        coefficient_covariance_status="ok",
        coefficient_statistics=np.array([1.0]),
        coefficient_p_values=np.array([0.3]),
        coefficient_confidence_lower=np.array([np.nan]),
        coefficient_confidence_upper=np.array([3.0]),
    )
    summary = summarize_glmm_coefficient(fit, 0, 1.0, 1.96)
    assert summary[3:] == ("", 3.0, "profile-endpoint-not-found")


def test_shape_bootstrap_preserves_tip_noise_when_no_shared_effects_are_added():
    from ete4 import Tree

    from nwkit.regression_pipeline import _simulate_response_tip_values

    tree = Tree("((A:1,B:1):1,C:2);", parser=1)
    inner = tree.children[0]
    leaves = list(tree.leaves())
    base = np.array([1.0, -2.0, 4.0])
    simulator = {
        "leaf_names": [str(leaf.name) for leaf in leaves],
        "evolutionary_model": SimpleNamespace(
            covariance_scale=1.0, sample=lambda rng, variance: base.copy()
        ),
        "sampling_factor": None,
        "orientation": {tree: tuple(tree.children), inner: tuple(inner.children)},
        "weights": {tree: (0.5, 0.5), inner: (0.5, 0.5)},
        "selected_nodes": [inner],
    }
    state = {
        "evolutionary_rate": 1.0,
        "design": np.zeros((1, 1)),
        "beta": np.zeros(1),
        "sample_shared_effects": lambda rng: np.zeros(1),
    }
    simulated = _simulate_response_tip_values(
        tree, state, simulator, np.random.default_rng(1)
    )
    np.testing.assert_allclose([simulated[leaf.name] for leaf in leaves], base)


def test_penalized_poisson_null_bootstrap_refits_both_models():
    x = np.column_stack([np.ones(8), np.linspace(-1.0, 1.0, 8)])
    fit = fit_phylogenetic_glmm(
        [0.0, 1.0, 0.0, 2.0, 1.0, 3.0, 2.0, 4.0],
        x,
        lambda _: np.eye(8),
        family="poisson",
        inference="null-bootstrap",
        bootstrap_replicates=2,
        seed=4,
    )
    assert np.all(np.isfinite(fit.coefficient_p_values))
    assert np.all(np.isin(fit.coefficient_p_values, [1 / 3, 2 / 3, 1]))
    assert summarize_glmm_coefficient(fit, 1, float(fit.coefficients.flat[1]), 1.96)[
        3:5
    ] == ("", "")


@pytest.mark.parametrize("grid", ["nan|1", "1|0", "1|1", "nope"])
def test_invalid_null_profile_grid_fails_before_covariance_evaluation(grid):
    def covariance(_):
        pytest.fail("Invalid grid must be rejected before fitting.")

    with pytest.raises(ValueError, match="grid"):
        fit_phylogenetic_glmm(
            [0.0, 1.0, 2.0],
            np.ones((3, 1)),
            covariance,
            family="poisson",
            inference="null-bootstrap",
            coefficient_profile_grid=grid,
        )


def test_student_t_penalty_null_statistic_matches_independent_normal_tail():
    from scipy.optimize import minimize_scalar

    from nwkit.phylogenetic_glmm import _coefficient_penalty_value

    def statistic(y):
        # Independent one-observation normal likelihood with t_3(scale=1) penalty.
        def objective(beta):
            return (y - beta) ** 2 / 2 + 2 * np.log1p(beta**2 / 3)

        optimum = minimize_scalar(objective, bounds=(-12, 12), method="bounded")
        assert optimum.success
        assert _coefficient_penalty_value(
            np.array([optimum.x]), "student-t", 1.0
        ) == pytest.approx(2 * np.log1p(optimum.x**2 / 3))
        return objective_difference(optimum.fun, objective(0))

    # The improvement is symmetric and strictly increases with |y|, hence
    # its exact null tail at y=2 is 2*Phi(-2), despite the different statistic.
    values = [statistic(y) for y in [0.5, 1, 1.5, 2, 3]]
    assert np.all(np.diff(values) > 0)
    assert statistic(-2) == pytest.approx(statistic(2))
    result = bootstrap_test(
        statistic(2), lambda rng: statistic(rng.normal()), replicates=2000, seed=19
    )
    assert abs(result.p_value - 2 * norm.sf(2)) < 4 * result.monte_carlo_se
