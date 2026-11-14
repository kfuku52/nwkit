"""Invariance and resource-limit contracts at numerical model boundaries."""

import math
import weakref
from decimal import Decimal, localcontext

import numpy as np
import pytest
from ete4 import Tree
from scipy import sparse

from nwkit import ordinary_regression, regress
from nwkit.evolution import (
    build_evolutionary_covariance,
    build_sparse_evolutionary_model,
)
from nwkit.gaussian import (
    DiagonalLowRankCovariance,
    draw_from_factor,
    factor_logdet,
    solve_factor,
)
from nwkit.measurement_error import (
    _positive_semidefinite,
    _predictor_error_covariance,
    fit_conditional_eiv_gaussian,
    fit_factor_latent_predictor,
    fit_latent_predictor,
    fit_sparse_latent_predictor,
)
from nwkit.optimization import FitResourceError, ScalarFitCache
from nwkit.phylogenetic_glmm import (
    _count_log_likelihood,
    _draw_beta_binomial,
    _gamma_log_likelihood,
    _gamma_scalar_terms,
    _lognormal_scalar_terms,
    _negative_binomial_scalar_terms,
    _normal_censor_interval_terms,
    _scalar_initial_coefficient,
    _scalar_terms,
    _scale_covariance_by_mean_variance,
    _zero_component_count_scalar_terms,
    _zero_truncation_terms,
)
from nwkit.regression_pipeline import _positive_semidefinite_loading
from nwkit.sparse_laplace import JointPredictorUncertainty

pytestmark = pytest.mark.filterwarnings(r"error::RuntimeWarning:nwkit\.")


def _data():
    rng = np.random.default_rng(52)
    x = rng.normal(size=40)
    return np.column_stack([np.ones(40), x]), 2.0 + 0.5 * x + rng.normal(size=40)


@pytest.mark.parametrize(
    "validate", [regress._validate_covariance_matrix, _positive_semidefinite]
)
@pytest.mark.parametrize("diagonal", [False, True])
def test_gaussian_covariance_readers_retain_small_variances(validate, diagonal):
    values = np.ldexp(np.array([[1.0, 0.25], [0.25, 0.5]]), -1050)
    if diagonal:
        values = np.diag(np.diag(values))
    with np.errstate(over="raise", invalid="raise"):
        np.testing.assert_array_equal(validate(values, "Sampling covariance"), values)
        with pytest.raises(ValueError, match="positive semidefinite"):
            validate(-values, "Sampling covariance")


@pytest.mark.parametrize("diagonal", [False, True])
def test_covariance_validation_preserves_small_components_of_a_mixed_scale(diagonal):
    values = np.array([[1.0, 2e-101], [2e-101, 1e-200]])
    if diagonal:
        values = np.diag(np.diag(values))
    np.testing.assert_array_equal(_positive_semidefinite(values, "sampling"), values)


def test_glmm_covariance_normalization_preserves_subnormal_components():
    tiny = math.ulp(0.0)
    mixed = np.diag([1.0, tiny])
    normalized = _scale_covariance_by_mean_variance(mixed, 2, "Test")
    np.testing.assert_array_equal(np.diag(normalized), [2.0, 2.0 * tiny])
    entirely_small = np.diag([tiny, 2.0 * tiny])
    normalized = _scale_covariance_by_mean_variance(entirely_small, 2, "Test")
    np.testing.assert_allclose(np.diag(normalized), [2 / 3, 4 / 3])


@pytest.mark.parametrize("exponent", [-1060, 1023])
def test_categorical_sampling_loading_preserves_covariance_units(exponent):
    covariance = np.ldexp(np.array([[1.0, 0.5], [0.5, 0.25]]), exponent)
    with np.errstate(over="raise", invalid="raise"):
        loading = _positive_semidefinite_loading(covariance, 2)
    assert loading.shape == (2, 1)
    np.testing.assert_allclose(
        np.ldexp(loading @ loading.T, -exponent),
        [[1.0, 0.5], [0.5, 0.25]],
    )


@pytest.mark.parametrize("exponent", [-530, 530])
def test_dense_predictor_covariance_scales_through_extreme_slope(exponent):
    base = np.array([[1.0, 0.25], [0.25, 0.5]])
    covariance_exponent = -1050 if exponent > 0 else 1022
    covariance = np.ldexp(base, covariance_exponent)
    with np.errstate(over="raise", invalid="raise"):
        projected = _predictor_error_covariance(
            np.array([np.ldexp(1.0, exponent)]), covariance, 0, 2
        )
    np.testing.assert_array_equal(
        projected, np.ldexp(base, covariance_exponent + 2 * exponent)
    )


def test_latent_predictor_retains_small_correlated_sampling_uncertainty():
    prior = np.array([[1.0, 0.2, 0.0], [0.2, 1.0, 0.1], [0.0, 0.1, 1.2]])
    sampling = 1e-24 * np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.25], [0.0, 0.25, 1.0]])
    result = fit_latent_predictor(
        np.array([0.0, 2.0, 1.0]), prior, sampling, include_intercept=True
    )
    scale = np.max(np.diag(sampling))
    np.testing.assert_allclose(
        result.covariance / scale, sampling / scale, rtol=1e-12, atol=1e-12
    )


@pytest.mark.parametrize("last", [-0.3, -0.1, 0.2])
@pytest.mark.parametrize("correlation", [0.0, 0.15])
def test_latent_predictor_keeps_intercept_when_precision_sum_overflows(
    last, correlation
):
    values = np.array([0.0, 1e-160, 0.3, last])
    prior = np.diag([0.0, 0.0, 1.0, 1.0])
    prior[2, 3] = prior[3, 2] = correlation
    sampling = np.array([1e-310, 1e-310, 0.1, 0.1])
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        fitted = fit_latent_predictor(values, prior, sampling, include_intercept=True)
    assert np.isfinite(fitted.log_likelihood)
    assert fitted.prior_mean == pytest.approx(5e-161, rel=1e-10, abs=0.0)
    np.testing.assert_allclose(fitted.mean[:2], 5e-161, rtol=1e-10, atol=0.0)
    np.testing.assert_array_equal(fitted.covariance[:2, :2], np.zeros((2, 2)))


@pytest.mark.parametrize("sampling", [0.0, 0.1])
@pytest.mark.parametrize("reml", [False, True])
def test_gaussian_fit_is_invariant_to_an_intercept_shift(sampling, reml):
    design, response = _data()
    fits = [
        regress._profile_covariance_fit(
            response + offset,
            design,
            np.full(len(response), sampling),
            [("rate", np.ones(len(response)))],
            reml=reml,
        )
        for offset in (0.0, 1e8)
    ]
    assert fits[1]["beta"] - [1e8, 0.0] == pytest.approx(fits[0]["beta"], abs=1e-7)
    for key in ("objective", "beta_covariance", "residual"):
        assert fits[1][key] == pytest.approx(fits[0][key], rel=2e-6, abs=1e-7)
    assert fits[1]["component_variances"] == pytest.approx(
        fits[0]["component_variances"], rel=2e-6
    )
    if sampling == 0.0:
        expected_beta = np.linalg.lstsq(design, response, rcond=None)[0]
        residual = response - design @ expected_beta
        variance = (
            residual @ residual / (len(response) - (design.shape[1] if reml else 0))
        )
        assert fits[0]["beta_covariance"] == pytest.approx(
            variance * np.linalg.inv(design.T @ design)
        )


def test_gaussian_component_variance_is_invariant_to_large_covariance_units():
    design, response = _data()
    scale = 2.0**1023
    covariance = np.ones(len(response))
    kwargs = {"reml": True}
    fixed = np.zeros(len(response))
    before = regress._profile_covariance_fit(
        response, design, fixed, [("rate", covariance)], **kwargs
    )
    with np.errstate(over="raise", invalid="raise"):
        after = regress._profile_covariance_fit(
            response, design, fixed, [("rate", scale * covariance)], **kwargs
        )
    np.testing.assert_array_equal(after["beta"], before["beta"])
    np.testing.assert_array_equal(after["beta_covariance"], before["beta_covariance"])
    assert after["objective"] == before["objective"]
    assert scale * after["component_variances"]["rate"] == pytest.approx(
        before["component_variances"]["rate"]
    )


@pytest.mark.parametrize("scale", [1e-6, 1e6])
def test_gaussian_variance_follows_response_units(scale):
    design, response = _data()
    fits = [
        regress._profile_covariance_fit(
            response * unit,
            design,
            np.zeros(len(response)),
            [("rate", np.ones(len(response)))],
            reml=True,
        )
        for unit in (1.0, scale)
    ]
    assert fits[1]["beta"] / scale == pytest.approx(fits[0]["beta"])
    assert fits[1]["beta_covariance"] / scale**2 == pytest.approx(
        fits[0]["beta_covariance"]
    )


@pytest.mark.parametrize("exponent", [-500, 500])
@pytest.mark.parametrize("structured", [False, True])
def test_gaussian_regression_preserves_fits_in_extreme_response_units(
    exponent, structured
):
    design, response = _data()
    size = len(response)
    if structured:
        fixed = DiagonalLowRankCovariance(
            np.full(size, 0.04), sparse.csr_matrix(np.full((size, 1), 0.2))
        )
        scaled_fixed = DiagonalLowRankCovariance(
            np.ldexp(fixed.diagonal, 2 * exponent),
            fixed.low_rank * math.ldexp(1.0, exponent),
        )
        component = np.ones(size)
    else:
        fixed = np.full(size, 0.04)
        scaled_fixed = np.ldexp(fixed, 2 * exponent)
        distances = np.abs(np.arange(size)[:, None] - np.arange(size))
        component = 0.4**distances
    components = [("evolutionary_rate", component)]
    reference = regress._profile_covariance_fit(
        response, design, fixed, components, reml=True
    )
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        scaled = regress._profile_covariance_fit(
            np.ldexp(response, exponent), design, scaled_fixed, components, reml=True
        )
    np.testing.assert_allclose(
        np.ldexp(scaled["beta"], -exponent), reference["beta"], atol=1e-7
    )
    np.testing.assert_allclose(
        np.ldexp(scaled["beta_covariance"], -2 * exponent),
        reference["beta_covariance"],
        rtol=1e-5,
        atol=1e-10,
    )
    assert scaled["objective"] - (size - 2) * exponent * math.log(2.0) == pytest.approx(
        reference["objective"], abs=1e-9
    )
    scaled_rate = np.ldexp(
        scaled["component_variances"]["evolutionary_rate"], -2 * exponent
    )
    assert scaled_rate == pytest.approx(
        reference["component_variances"]["evolutionary_rate"], rel=1e-5
    )
    assert factor_logdet(scaled["cholesky"]) - 2 * size * exponent * math.log(
        2.0
    ) == pytest.approx(factor_logdet(reference["cholesky"]), abs=1e-4)
    scaled_response = np.ldexp(response, exponent)
    np.testing.assert_allclose(
        np.ldexp(solve_factor(scaled["cholesky"], scaled_response), exponent),
        solve_factor(reference["cholesky"], response),
        rtol=1e-5,
        atol=1e-8,
    )
    draws = np.random.default_rng(22).standard_normal(size)
    before = draw_from_factor(
        reference["cholesky"], draws, rng=np.random.default_rng(15)
    )
    after = draw_from_factor(scaled["cholesky"], draws, rng=np.random.default_rng(15))
    np.testing.assert_allclose(np.ldexp(after, -exponent), before, rtol=1e-5, atol=1e-8)


def test_no_intercept_fit_retains_the_response_mean():
    design, response = _data()
    design = design[:, 1:]
    fit = regress._profile_covariance_fit(
        response,
        design,
        np.zeros(len(response)),
        [("rate", np.ones(len(response)))],
        reml=True,
    )
    expected = np.linalg.lstsq(design, response, rcond=None)[0]
    assert fit["beta"] == pytest.approx(expected)
    assert fit["residual"] == pytest.approx(response - design @ expected)


def test_exact_constant_response_reports_a_variance_boundary():
    fit = regress._profile_covariance_fit(
        np.ones(4), np.ones((4, 1)), np.zeros(4), [("rate", np.ones(4))], reml=False
    )
    assert fit["boundary_warning"]
    assert np.isfinite(fit["objective"])
    assert fit["residual"] == pytest.approx(np.zeros(4))


def test_dense_limit_rejects_auto_lambda_before_allocating(monkeypatch):
    monkeypatch.setattr(regress, "MAX_DENSE_GAUSSIAN_OBSERVATIONS", 3)

    def unexpected_allocation(*args, **kwargs):
        pytest.fail("A covariance was allocated before the size preflight")

    monkeypatch.setattr(
        ordinary_regression, "build_phylogenetic_covariance", unexpected_allocation
    )
    with pytest.raises(FitResourceError, match="Dense Gaussian covariance fitting"):
        ordinary_regression._fit_ordinary_gaussian(
            np.arange(4.0),
            np.ones((4, 1)),
            np.zeros(4),
            Tree("((A:1,B:1):1,(C:1,D:1):1);", parser=1),
            list("ABCD"),
            evolution_model="lambda",
            evolution_parameter=None,
            branch_length="original",
            custom_covariance=None,
            reml=False,
        )


@pytest.mark.parametrize("branch_length", ["original", "unit"])
def test_reused_lambda_covariances_match_tree_transforms(monkeypatch, branch_length):
    tree = Tree("((A:1,B:2):0.7,(C:0.2,D:3):1.1);", parser=1)
    original = ordinary_regression.build_phylogenetic_covariance
    builds = []

    def counted(*args, **kwargs):
        builds.append(kwargs["parameter"])
        return original(*args, **kwargs)

    monkeypatch.setattr(ordinary_regression, "build_phylogenetic_covariance", counted)
    covariance_at = ordinary_regression._ordinary_covariance_builder(
        tree,
        list("DCBA"),
        evolution_model="lambda",
        evolution_parameter=None,
        branch_length=branch_length,
        custom_covariance=None,
    )
    for parameter in (0.0, 1e-12, 0.2, 0.999999, 1.0):
        actual = covariance_at(parameter)
        expected = original(
            tree,
            list("DCBA"),
            evolution_model="lambda",
            parameter=parameter,
            branch_length=branch_length,
        )
        assert actual == pytest.approx(expected, abs=1e-14)
        actual[:] = -1.0  # a consumer cannot corrupt later trial covariances
    assert builds == [1.0]


@pytest.mark.parametrize("kind", ["dense", "sparse", "factor"])
def test_latent_predictor_is_invariant_to_an_intercept_shift(kind):
    design, response = _data()
    names = [f"T{i}" for i in range(len(response))]
    tree = Tree("(" + ",".join(f"{name}:1" for name in names) + ");", parser=1)

    def fit(observed):
        if kind == "dense":
            return fit_latent_predictor(
                observed,
                build_evolutionary_covariance(tree, names),
                np.full(len(response), 0.1),
                include_intercept=True,
            )
        if kind == "sparse":
            return fit_sparse_latent_predictor(
                observed,
                build_sparse_evolutionary_model(tree, names),
                np.full(len(response), 0.1),
                include_intercept=True,
            )
        return fit_factor_latent_predictor(
            observed,
            np.ones(len(response)),
            DiagonalLowRankCovariance(np.zeros(len(response)), design * 0.1),
            include_intercept=True,
        )

    first, shifted = fit(response), fit(response + 1e8)
    assert shifted.mean - 1e8 == pytest.approx(first.mean, abs=1e-7)
    assert shifted.prior_mean - 1e8 == pytest.approx(first.prior_mean, abs=1e-7)
    assert shifted.evolutionary_rate == pytest.approx(first.evolutionary_rate, rel=2e-6)
    assert shifted.log_likelihood == pytest.approx(first.log_likelihood, rel=2e-7)


def test_eiv_fit_and_information_are_invariant_to_an_intercept_shift():
    design, response = _data()
    fits = [
        fit_conditional_eiv_gaussian(
            response + offset,
            design,
            [np.eye(len(response)) * 0.04],
            [1],
            np.full(len(response), 0.1),
            [("rate", np.ones(len(response)))],
            reml=False,
        )
        for offset in (0.0, 1e8)
    ]
    assert fits[1]["beta"] - [1e8, 0.0] == pytest.approx(fits[0]["beta"], abs=2e-6)
    assert fits[1]["beta_covariance"] == pytest.approx(
        fits[0]["beta_covariance"], rel=1e-4
    )
    assert fits[1]["objective"] == pytest.approx(fits[0]["objective"], rel=1e-7)


@pytest.mark.parametrize("constructor", [sparse.csr_matrix, sparse.csr_array])
def test_sparse_predictor_uncertainty_matches_a_dense_gaussian_fit(constructor):
    design, response = _data()
    loading = np.column_stack(
        [np.linspace(-0.1, 0.1, len(response)), np.full(len(response), 0.03)]
    )
    fits = []
    for uncertainty in (
        loading @ loading.T,
        JointPredictorUncertainty(factors=(constructor(loading),)),
    ):
        columns = (1,) if isinstance(uncertainty, JointPredictorUncertainty) else 1
        fits.append(
            fit_conditional_eiv_gaussian(
                response,
                design,
                [uncertainty],
                [columns],
                np.full(len(response), 0.1),
                [("rate", np.ones(len(response)))],
                reml=False,
            )
        )
    for key in ("beta", "beta_covariance", "objective"):
        assert fits[1][key] == pytest.approx(fits[0][key], rel=2e-5, abs=1e-7)


def test_scalar_fit_cache_does_not_suppress_resource_errors():
    def fit(value):
        if value:
            raise FitResourceError("too large")
        return {"objective": 1.0}

    cache = ScalarFitCache(fit, lambda result: result["objective"])
    assert cache.objective(0.0) == 1.0
    with pytest.raises(FitResourceError, match="too large"):
        cache.objective(0.5)


def test_likelihood_search_releases_evicted_arrays_and_keeps_scalar_scores():
    references = []
    evaluated = []

    def fit(value):
        covariance = np.full((32, 32), float(value))
        references.append(weakref.ref(covariance))
        evaluated.append(value)
        return {"objective": (value - 7.0) ** 2, "covariance": covariance}

    cache = ScalarFitCache(fit, lambda result: result["objective"])
    for value in range(30):
        cache.objective(value)
        assert sum(ref() is not None for ref in references) <= 2
    assert cache.objective(7) == 0.0
    assert len(evaluated) == 30  # revisiting an evicted score needs no matrix
    best = cache.best(iter(range(30)))
    assert best["objective"] == 0.0
    assert np.all(best["covariance"] == 7.0)
    assert sum(ref() is not None for ref in references) == 1


@pytest.mark.parametrize(
    "model,parameter", [("independent", None), ("lambda", 0.0), ("custom", None)]
)
def test_dense_limit_does_not_reject_diagonal_models(monkeypatch, model, parameter):
    monkeypatch.setattr(regress, "MAX_DENSE_GAUSSIAN_OBSERVATIONS", 3)
    fit = ordinary_regression._fit_ordinary_gaussian(
        np.asarray([0.0, 1.0, 3.0, 2.0]),
        np.ones((4, 1)),
        np.zeros(4),
        Tree("((A:1,B:1):1,(C:1,D:1):1);", parser=1),
        list("ABCD"),
        evolution_model=model,
        evolution_parameter=parameter,
        branch_length="original",
        custom_covariance=np.eye(4) if model == "custom" else None,
        reml=False,
    )
    assert np.isfinite(fit["objective"])


def test_explicit_dense_opt_in_keeps_the_full_lambda_search(monkeypatch):
    tree = Tree("((A:1,B:1):1,(C:1,D:1):1);", parser=1)

    def fit(**kwargs):
        return ordinary_regression._fit_ordinary_gaussian(
            np.asarray([0.0, 0.1, 2.0, 2.1]),
            np.ones((4, 1)),
            np.zeros(4),
            tree,
            list("ABCD"),
            evolution_model="lambda",
            evolution_parameter=None,
            branch_length="original",
            custom_covariance=None,
            reml=False,
            **kwargs,
        )

    expected = fit()
    monkeypatch.setattr(regress, "MAX_DENSE_GAUSSIAN_OBSERVATIONS", 3)
    with pytest.warns(
        RuntimeWarning, match="Large dense allocation was explicitly enabled"
    ):
        observed = fit(allow_large_dense=True)
    assert observed["evolution_parameter"] == pytest.approx(
        expected["evolution_parameter"]
    )
    assert observed["objective"] == pytest.approx(expected["objective"])
    assert observed["evolution_parameter"] > 0.0


@pytest.mark.parametrize("dispersion", [None, 0.2])
def test_hurdle_likelihood_retains_small_mean_limit(dispersion):
    family = "hurdle-poisson" if dispersion is None else "hurdle-negative-binomial"
    counts = np.array([0.0, 1.0, 2.0, 1.0])
    eta = np.array([-1000.0, -1000.0, -1000.0, -100.0])
    log_likelihood = _count_log_likelihood(
        counts, eta, family, dispersion, 0.2, np.zeros(len(eta))
    )

    expected = [
        np.log(0.2),
        np.log(0.8),
        np.log(0.8) - 1000.0 + np.log((1.0 + (dispersion or 0.0)) / 2.0),
        np.log(0.8),
    ]
    np.testing.assert_allclose(log_likelihood, expected, rtol=0.0, atol=1e-12)


@pytest.mark.parametrize("dispersion", [9e-7, 1e-14, 1e-20, 1e-310])
def test_negative_binomial_likelihood_preserves_near_poisson_limit(dispersion):
    counts = np.array([0.0, 1.0, 5.0, 1000.0, 10_000.0])
    means = np.array([0.5, 1.0, 4.0, 1200.0, 10_000.0])
    eta = np.log(means)
    zero_offset = np.zeros_like(eta)
    expected = []
    for count, mean, log_mean in zip(counts, means, eta, strict=True):
        expected.append(
            math.fsum(math.log1p(dispersion * j) for j in range(int(count)))
            - math.lgamma(count + 1.0)
            + count * (log_mean - math.log1p(dispersion * mean))
            - math.log1p(dispersion * mean) / dispersion
        )
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        actual = _count_log_likelihood(
            counts, eta, "negative-binomial", dispersion, None, zero_offset
        )
        gradient, curvature = _negative_binomial_scalar_terms(
            counts, eta, zero_offset, dispersion
        )
    np.testing.assert_allclose(actual, expected, rtol=0.0, atol=7e-11)
    denominator = 1.0 + dispersion * means
    np.testing.assert_allclose(gradient, (means - counts) / denominator, atol=1e-10)
    np.testing.assert_allclose(
        curvature, means * (1.0 + dispersion * counts) / denominator**2, rtol=1e-12
    )


def test_poisson_log_mass_resolves_large_counts_near_the_fitted_mean():
    counts = np.array([1e12, 1e15, 1e24, 1e100, 1e308])
    eta = np.log(counts) + np.array([0.0, 1e-7, -1e-12, 0.0, 0.0])
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        mass = _count_log_likelihood(counts, eta, "poisson", None, None, 0.0)
    # Evaluate the reference around the same representable log-count used in
    # the deviance to distinguish cancellation from rounding the input logit.
    relative = eta - np.log(counts)
    expected_mass = -0.5 * (np.log(counts) + math.log(math.tau))
    expected_mass -= (1.0 / counts) / 12.0
    expected_mass -= ((counts * relative) * relative) * (0.5 + relative / 6.0)
    np.testing.assert_allclose(mass, expected_mass, rtol=0.0, atol=5e-13)


@pytest.mark.parametrize(
    "dispersion,count,shift,expected",
    [
        (1e-8, 1e18, 0.099, -474303.9574402318),
        (1e-6, 1e6, -0.1, -2507.13229477862),
        (0.7, 1e18, 0.0, -42.24461542763349),
        (1e308, 1e18, -0.2, -750.6427403160588),
        (1e-310, 1e308, 0.0, -355.52201801971427),
        (0.7, 1e308, 1.0, -710.5198344547231),
    ],
)
def test_large_negative_binomial_counts_match_high_precision_reference(
    dispersion, count, shift, expected
):
    counts = np.array([count])
    log_mean = np.log(counts) + shift
    # Reference from the log-gamma definition at 750 decimal digits.
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        actual = _count_log_likelihood(
            counts, log_mean, "negative-binomial", dispersion, None, 0.0
        )
    assert actual[0] == pytest.approx(expected, rel=1e-11)


def test_beta_binomial_retains_endpoint_and_interior_mass_at_tiny_precision():
    dispersion = 1e-310
    values = np.array([0.0, 1.0, 2.0, 4.0, 0.0])
    trials = np.array([4.0, 4.0, 4.0, 4.0, 0.0])
    linear = np.full(len(values), math.log(3.0))
    with np.errstate(over="raise", divide="raise", invalid="raise"):
        likelihood, gradient, curvature = _scalar_terms(
            values,
            linear,
            family="beta-binomial",
            dispersion=dispersion,
            zero_probability=None,
            offset=np.zeros_like(values),
            trials=trials,
            censor_lower=None,
            censor_upper=None,
        )
    # The scaled endpoint masses are 1/4 and 3/4, and the interior
    # masses are precision * 3/16 * 4 / (x * (4-x)).
    expected = math.log(0.25 * 0.75 * 0.25 * (3 / 16)) + 2 * math.log(dispersion)
    assert likelihood == pytest.approx(expected, abs=1e-12)
    np.testing.assert_allclose(gradient, [0.75, 0.5, 0.5, -0.25, 0.0])
    np.testing.assert_allclose(curvature, [3 / 16, 3 / 8, 3 / 8, 3 / 16, 0.0])
    rng = np.random.default_rng(17)
    simulated = _draw_beta_binomial(
        rng, np.full(1200, math.log(3.0)), dispersion, np.full(1200, 4.0)
    )
    assert np.all((simulated == 0.0) | (simulated == 4.0))
    assert simulated.mean() == pytest.approx(3.0, abs=0.16)


@pytest.mark.parametrize(
    "precision, trials, value, linear, expected",
    [
        (2e12, 10.0, 0.0, 0.0, (-6.931471805576953, 4.9999999999775, 2.4999999999775)),
        (1e308, 10.0, 5.0, 0.0, (-1.4020427180880297, 0.0, 2.5)),
        (1e100, 1e200, 5e199, 0.0, (-345.6135553017516, 0.0, 2.5e99)),
        (
            1e308,
            1e308,
            1e307,
            0.0,
            (-2.0349845015839338e307, 2.118244650968009e307, 1.0119047619047619e307),
        ),
    ],
)
def test_beta_binomial_preserves_mass_and_derivatives_at_high_precision(
    precision, trials, value, linear, expected
):
    # Values from the log-gamma and polygamma definitions at 770 digits.
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        log_mass, gradient, curvature = _scalar_terms(
            np.array([value]),
            np.array([linear]),
            family="beta-binomial",
            dispersion=precision,
            zero_probability=None,
            offset=np.zeros(1),
            trials=np.array([trials]),
            censor_lower=None,
            censor_upper=None,
        )
    np.testing.assert_allclose(
        [log_mass, gradient[0], curvature[0]], expected, rtol=1e-11
    )


def test_gamma_likelihood_and_derivatives_scale_with_response_units():
    values = np.array([0.75, 2.0, 3.0])
    linear = np.array([0.3, -0.2, 0.6])
    dispersion = 2.0
    scale = 1e-310
    with np.errstate(over="raise", invalid="raise"):
        before = _gamma_log_likelihood(values, linear, dispersion)
        after = _gamma_log_likelihood(
            values * scale, linear + np.log(scale), dispersion
        )
        before_terms = _gamma_scalar_terms(values, linear, dispersion)
        after_terms = _gamma_scalar_terms(
            values * scale, linear + np.log(scale), dispersion
        )
    np.testing.assert_allclose(after + np.log(scale), before, atol=3e-12)
    np.testing.assert_allclose(after_terms, before_terms, atol=3e-12)


@pytest.mark.parametrize("shape", [1e-310, math.ulp(0.0)])
def test_gamma_likelihood_remains_finite_for_subnormal_shape(shape):
    values = np.array([0.5, 2.0])
    log_density = _gamma_log_likelihood(values, np.zeros(2), shape)
    np.testing.assert_allclose(
        log_density, math.log(shape) - np.log(values), atol=1e-13
    )


@pytest.mark.parametrize("shape", [1e18, 1e200, 1e308])
def test_high_shape_gamma_log_density_resolves_variation_near_its_mean(shape):
    # In log-response units the high-shape limit has variance 1 / shape.
    standardized = np.array([0.0, -0.25, 1.0, 3.0])
    linear = -standardized / math.sqrt(shape)
    with np.errstate(over="raise", invalid="raise"):
        result = _gamma_log_likelihood(np.ones(len(linear)), linear, shape)
    normal_limit = 0.5 * (math.log(shape) - math.log(math.tau)) - (
        0.5 * standardized**2
    )
    np.testing.assert_allclose(result, normal_limit, rtol=0.0, atol=5e-9)


def test_lognormal_score_avoids_overflow_in_squared_dispersion():
    gradient, curvature = _lognormal_scalar_terms(
        np.ones(2), np.array([1e308, -1e308]), 1e155
    )
    np.testing.assert_allclose(gradient, [0.01, -0.01])
    np.testing.assert_allclose(curvature, [1e-310, 1e-310], atol=0.0)


def test_gaussian_censoring_scores_keep_tails_and_small_curvature():
    values = np.array([1.0, np.nan, np.nan, np.nan, np.nan])
    linear = np.array([101.0, 100.0, -22.0, 8e6, -800.0])
    lower = np.array([np.nan, np.nan, 100.0, np.nan, np.nan])
    upper = np.array([np.nan, -22.0, np.nan, 0.0, 0.0])
    sd = 4.0
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        likelihood, score, weight = _scalar_terms(
            values,
            linear,
            family="censored-gaussian",
            dispersion=sd,
            zero_probability=None,
            offset=np.zeros_like(values),
            trials=None,
            censor_lower=lower,
            censor_upper=upper,
        )
    assert np.isfinite(likelihood)
    assert score[0] == 100 / sd**2
    assert weight[0] == 1 / sd**2
    assert score[1] == -score[2]
    assert score[1] == pytest.approx(7.6331791926650395, rel=1e-12)
    assert weight[1] == pytest.approx(0.06243324326357813, rel=1e-12)
    assert weight[1] == weight[2]
    assert score[3] == pytest.approx(8e6 / sd**2, rel=1e-12)
    assert weight[3] == pytest.approx(1 / sd**2, rel=1e-12)
    assert score[4] == weight[4] == 0.0


def test_interval_censoring_derivatives_match_high_precision_normal_moments():
    lower = np.array([-1.25, 8.1, 30.0, -30.2, 0.0])
    upper = np.array([0.75, 8.3, 30.2, -30.0, math.ulp(0.0)])
    center = np.array([0.0, 0.0, 0.0, 0.0, 11.0])
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        mass, score, weight = _normal_censor_interval_terms(
            lower, upper, center, 1.0, derivatives=True
        )
    # Reference moments from the normal density and upper-tail integrals
    # evaluated to 160 decimal digits, plus the limiting narrow-interval case.
    expected_mass = [
        -0.40388205087032925,
        -36.04052541639861,
        -454.32366048808217,
        -454.32366048808217,
        math.log(math.ulp(0.0)) - 60.5 - 0.5 * math.log(math.tau),
    ]
    expected_score = [
        0.17745138198386967,
        -8.173850602678082,
        -30.032776307243488,
        30.032776307243488,
        11.0,
    ]
    expected_weight = [
        0.7116580399881117,
        0.997074221430888,
        0.9989929931365967,
        0.9989929931365967,
        1.0,
    ]
    np.testing.assert_allclose(mass, expected_mass, rtol=1e-12)
    np.testing.assert_allclose(score, expected_score, rtol=1e-12)
    np.testing.assert_allclose(weight, expected_weight, rtol=1e-12)


def test_count_intercept_initialization_uses_log_exposure():
    with np.errstate(over="raise", divide="raise"):
        coefficients = _scalar_initial_coefficient(
            np.array([1.0, 0.0, 4.0]),
            np.ones((3, 1)),
            "poisson",
            np.full(3, -1000.0),
            None,
        )
    assert coefficients[0, 0] == pytest.approx(1000.0 + np.log(5.0 / 3.0))


@pytest.mark.parametrize("dispersion", [None, 0.2])
def test_hurdle_derivatives_are_stable_at_small_and_large_means(dispersion):
    eta = np.asarray([-100.0, -20.0, -2.0, 0.0, 2.0, 10.0, 100.0])
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        correction, derivative = _zero_truncation_terms(eta, dispersion)
        gradient, curvature = _zero_component_count_scalar_terms(
            np.ones(len(eta)), eta, np.zeros(len(eta)), dispersion, 0.1, hurdle=True
        )
    assert np.isfinite(gradient).all()
    assert np.isfinite(curvature).all()
    assert correction[0] == pytest.approx(1.0)
    assert derivative[0] == pytest.approx(0.0, abs=1e-15)
    assert correction[-1] == pytest.approx(0.0, abs=1e-100)
    assert derivative[-1] == pytest.approx(0.0, abs=1e-100)
    # Independent central differences check the derivative away from the limits.
    middle = eta[2:5]
    plus, _ = _zero_truncation_terms(middle + 1e-5, dispersion)
    minus, _ = _zero_truncation_terms(middle - 1e-5, dispersion)
    assert derivative[2:5] == pytest.approx((plus - minus) / 2e-5, rel=1e-7)


@pytest.mark.parametrize("dispersion", [0.01, 0.2, 1.0, 5.0])
def test_negative_binomial_zero_truncation_matches_decimal_reference(dispersion):
    eta = np.asarray([-100.0, -60.0, -40.0, -25.0, -12.0, -10.0, -5.0, 0.0, 5.0, 100.0])
    expected_correction, expected_derivative = [], []
    with localcontext() as context:
        context.prec = 100
        size = 1 / Decimal.from_float(dispersion)
        for value in eta:
            mean = Decimal.from_float(float(value)).exp()
            fraction = mean / (size + mean)
            q = size * fraction
            exponent = size * (1 + mean / size).ln()
            correction = q / (exponent.exp() - 1)
            derivative = correction * (1 - fraction - q - correction)
            expected_correction.append(float(correction))
            expected_derivative.append(float(derivative))
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        correction, derivative = _zero_truncation_terms(eta, dispersion)
    np.testing.assert_allclose(correction, expected_correction, rtol=5e-13, atol=5e-15)
    np.testing.assert_allclose(derivative, expected_derivative, rtol=5e-13, atol=5e-15)
