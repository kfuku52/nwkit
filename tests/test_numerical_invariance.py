"""Invariance and resource-limit contracts at numerical model boundaries."""

import weakref

import numpy as np
import pytest
from ete4 import Tree
from scipy import sparse

from nwkit import ordinary_regression, regress
from nwkit.evolution import (
    build_evolutionary_covariance,
    build_sparse_evolutionary_model,
)
from nwkit.gaussian import DiagonalLowRankCovariance
from nwkit.measurement_error import (
    fit_conditional_eiv_gaussian,
    fit_factor_latent_predictor,
    fit_latent_predictor,
    fit_sparse_latent_predictor,
)
from nwkit.optimization import FitResourceError, ScalarFitCache
from nwkit.phylogenetic_glmm import (
    _zero_component_count_scalar_terms,
    _zero_truncation_terms,
)
from nwkit.sparse_laplace import JointPredictorUncertainty

pytestmark = pytest.mark.filterwarnings(r"error::RuntimeWarning:nwkit\.")


def _data():
    rng = np.random.default_rng(52)
    x = rng.normal(size=40)
    return np.column_stack([np.ones(40), x]), 2.0 + 0.5 * x + rng.normal(size=40)


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
