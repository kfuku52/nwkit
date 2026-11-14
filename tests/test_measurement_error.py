import numpy as np
import pytest
from ete4 import Tree
from scipy import sparse
from scipy.optimize import minimize_scalar

from nwkit import measurement_error as measurement_error_mod
from nwkit.evolution import (
    build_evolutionary_covariance,
    build_sparse_evolutionary_model,
)
from nwkit.gaussian import (
    DiagonalLowRankCovariance,
    DiagonalSparsePrecisionCovariance,
    draw_from_factor,
    factor_diagonal_sparse_precision_updates,
    materialize_covariance,
    sparse_precision_update_diagonal,
)
from nwkit.measurement_error import (
    fit_conditional_eiv_gaussian,
    fit_factor_latent_predictor,
    fit_latent_predictor,
    fit_sparse_latent_predictor,
)
from nwkit.sparse_laplace import (
    ContinuousPredictorUncertainty,
    GmrfPredictorUncertainty,
    GroupedPredictorUncertainty,
    factor_sparse_positive_definite,
    grouped_predictor_loading,
)


def test_grouped_predictor_loading_keeps_order_and_shared_species_factors():
    factors = (np.array([[1.0], [2.0]]), np.empty((2, 0)), np.array([[3.0], [0.0]]))
    indices = np.array([2, 0, 1, 2])
    loading = grouped_predictor_loading(
        GroupedPredictorUncertainty(factors, indices), np.eye(2)
    )
    np.testing.assert_array_equal(
        loading.toarray(),
        [[0, 3], [0, 0], [1, 0], [2, 0], [0, 0], [0, 0], [0, 3], [0, 0]],
    )


def test_latent_predictor_posterior_matches_gaussian_conditioning():
    observed = np.asarray([1.0, 2.0, -0.5])
    evolutionary = np.diag([1.0, 2.0, 3.0])
    sampling = np.diag([0.2, 0.4, 0.1])

    posterior = fit_latent_predictor(
        observed,
        evolutionary,
        sampling,
        include_intercept=False,
    )

    prior = posterior.evolutionary_rate * evolutionary
    expected_mean = prior @ np.linalg.solve(prior + sampling, observed)
    expected_covariance = prior - prior @ np.linalg.solve(prior + sampling, prior)
    assert posterior.mean == pytest.approx(expected_mean)
    assert posterior.covariance == pytest.approx(expected_covariance)
    assert posterior.evolutionary_rate > 0.0
    assert np.linalg.eigvalsh(posterior.covariance).min() >= -1e-12


@pytest.mark.parametrize("scale", [1e-150, 1e150])
def test_latent_predictor_preserves_covariance_across_trait_units(scale):
    values = np.array([2.0, -1.0, 3.0, -3.0])
    evolutionary = np.diag([1.0, 0.5, 0.8, 1.5])
    sampling = np.full(len(values), 0.3)
    expected = fit_latent_predictor(
        values, evolutionary, sampling, include_intercept=True
    )

    with np.errstate(over="raise", invalid="raise"):
        scaled = fit_latent_predictor(
            values * scale, evolutionary, sampling * scale**2, include_intercept=True
        )

    np.testing.assert_allclose(scaled.mean / scale, expected.mean, rtol=1e-6)
    np.testing.assert_allclose(
        scaled.covariance / scale**2, expected.covariance, rtol=1e-6
    )
    assert scaled.evolutionary_rate / scale**2 == pytest.approx(
        expected.evolutionary_rate, rel=2e-6
    )
    assert scaled.log_likelihood + len(values) * np.log(scale) == pytest.approx(
        expected.log_likelihood, abs=1e-9
    )


def test_latent_predictor_rejects_indefinite_sampling_covariance():
    with pytest.raises(ValueError, match="positive semidefinite"):
        fit_latent_predictor(
            np.asarray([1.0, 2.0]),
            np.eye(2),
            np.asarray([[1.0, 2.0], [2.0, 1.0]]),
            include_intercept=False,
        )


def test_dense_eiv_accepts_fixed_covariance_as_a_sampling_factor():
    rng = np.random.default_rng(873)
    n = 16
    design = np.column_stack((np.ones(n), rng.normal(size=n)))
    response = design @ [0.3, -0.8] + rng.normal(size=n)
    loading = sparse.csr_matrix(rng.normal(scale=0.1, size=(n, 3)))
    fixed = DiagonalLowRankCovariance(np.full(n, 0.15), loading)
    predictor_uncertainty = [np.diag(rng.uniform(0.1, 0.3, n))]
    components = [("evolutionary_rate", np.ones(n))]
    before = fit_conditional_eiv_gaussian(
        response,
        design,
        predictor_uncertainty,
        [1],
        materialize_covariance(fixed),
        components,
        reml=False,
    )
    after = fit_conditional_eiv_gaussian(
        response, design, predictor_uncertainty, [1], fixed, components, reml=False
    )
    for key in ("beta", "beta_covariance", "covariance"):
        np.testing.assert_allclose(after[key], before[key], rtol=1e-7, atol=1e-8)
    assert after["objective"] == pytest.approx(before["objective"], abs=1e-10)


@pytest.mark.parametrize("factor_path", [False, True])
@pytest.mark.parametrize("exponent", [-1023, 1023])
def test_latent_predictor_normalizes_extreme_evolutionary_covariances(
    factor_path, exponent
):
    observed = np.array([0.4, -0.5, 0.6, 0.0, -0.25, 0.56])
    variances = np.array([0.75, 0.75, 0.5, 0.5, 0.75, 0.5])
    if factor_path:
        fit = fit_factor_latent_predictor
        evolutionary = variances
        sampling = DiagonalLowRankCovariance(
            np.zeros(len(observed)), sparse.eye(len(observed), format="csr") * 0.07
        )
    else:
        fit = fit_latent_predictor
        distance = np.abs(np.arange(len(observed))[:, None] - np.arange(len(observed)))
        evolutionary = 0.5 * 0.25**distance + np.diag(variances - 0.5)
        sampling = np.full(len(observed), 0.005)
    reference = fit(observed, evolutionary, sampling, include_intercept=True)
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        scaled = fit(
            observed, np.ldexp(evolutionary, exponent), sampling, include_intercept=True
        )
    assert np.ldexp(scaled.evolutionary_rate, exponent) == pytest.approx(
        reference.evolutionary_rate
    )
    assert scaled.log_likelihood == reference.log_likelihood
    np.testing.assert_array_equal(scaled.mean, reference.mean)
    if factor_path:
        np.testing.assert_array_equal(
            scaled.covariance_model.materialize(),
            reference.covariance_model.materialize(),
        )
    else:
        np.testing.assert_array_equal(scaled.covariance, reference.covariance)


def test_sparse_latent_predictor_matches_dense_conditioning():
    tree = Tree("(((A:1,B:1):1,C:2):1,(D:1,E:1):2);", parser=1)
    names = ["A", "B", "C", "D", "E"]
    observed = np.asarray([1.0, 1.2, 2.0, 3.0, 3.3])
    sampling = np.asarray([0.1, 0.2, 0.15, 0.3, 0.12])
    dense = fit_latent_predictor(
        observed,
        build_evolutionary_covariance(tree, names),
        np.diag(sampling),
        include_intercept=True,
    )
    sparse_fit = fit_sparse_latent_predictor(
        observed,
        build_sparse_evolutionary_model(tree, names),
        sampling,
        include_intercept=True,
    )
    assert sparse_fit.evolutionary_rate == pytest.approx(
        dense.evolutionary_rate, rel=1e-6
    )
    assert sparse_fit.log_likelihood == pytest.approx(dense.log_likelihood, rel=1e-9)
    np.testing.assert_allclose(sparse_fit.mean, dense.mean, rtol=1e-7, atol=1e-7)
    np.testing.assert_allclose(
        sparse_fit.covariance_model.materialize(),
        dense.covariance,
        rtol=1e-7,
        atol=1e-7,
    )


@pytest.mark.parametrize("scale", [1.0, 1e-150, 1e150])
def test_sparse_latent_predictor_supports_low_rank_batch_sampling_covariance(scale):
    tree = Tree("(((A:1,B:1):1,C:2):1,(D:1,E:1):2);", parser=1)
    names = ["A", "B", "C", "D", "E"]
    observed = np.asarray([1.0, 1.2, 2.0, 3.0, 3.3])
    sampling = DiagonalLowRankCovariance(
        np.asarray([0.1, 0.2, 0.15, 0.3, 0.12]),
        np.asarray([[0.1], [-0.1], [0.05], [0.2], [-0.2]]),
    )
    dense = fit_latent_predictor(
        observed,
        build_evolutionary_covariance(tree, names),
        materialize_covariance(sampling),
        include_intercept=True,
    )
    scaled_sampling = DiagonalLowRankCovariance(
        sampling.diagonal * scale**2, sampling.low_rank * scale
    )
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        sparse_fit = fit_sparse_latent_predictor(
            observed * scale,
            build_sparse_evolutionary_model(tree, names),
            scaled_sampling,
            include_intercept=True,
        )

    assert sparse_fit.evolutionary_rate / scale**2 == pytest.approx(
        dense.evolutionary_rate, rel=2e-6
    )
    assert sparse_fit.log_likelihood + len(observed) * np.log(scale) == pytest.approx(
        dense.log_likelihood, rel=1e-8
    )
    np.testing.assert_allclose(
        sparse_fit.mean / scale, dense.mean, rtol=1e-6, atol=1e-7
    )
    np.testing.assert_allclose(
        sparse_fit.covariance_model.materialize() / scale**2,
        dense.covariance,
        rtol=1e-6,
        atol=1e-7,
    )


@pytest.mark.parametrize("scale", [1.0, 1e-150, 1e150])
def test_factor_latent_predictor_matches_dense_conditioning_without_materializing(
    scale,
):
    observed = np.asarray([1.0, 1.2, 2.0, 3.0, 3.3])
    evolutionary = np.asarray([1.0, 1.5, 0.8, 2.0, 1.2])
    loading = sparse.csr_matrix(
        np.asarray(
            [
                [0.2, 0.0, 0.1],
                [0.1, 0.2, 0.0],
                [0.0, 0.1, 0.2],
                [0.3, 0.0, 0.1],
                [0.0, 0.2, 0.1],
            ]
        )
    )
    sampling = DiagonalLowRankCovariance(np.zeros(5), loading)
    dense = fit_latent_predictor(
        observed,
        np.diag(evolutionary),
        materialize_covariance(sampling),
        include_intercept=False,
    )
    scaled_sampling = DiagonalLowRankCovariance(sampling.diagonal, loading * scale)
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        structured = fit_factor_latent_predictor(
            observed * scale,
            evolutionary,
            scaled_sampling,
            include_intercept=False,
        )

    assert structured.evolutionary_rate / scale**2 == pytest.approx(
        dense.evolutionary_rate, rel=2e-6
    )
    assert structured.log_likelihood + len(observed) * np.log(scale) == pytest.approx(
        dense.log_likelihood, rel=1e-8
    )
    np.testing.assert_allclose(
        structured.mean / scale, dense.mean, rtol=1e-6, atol=1e-7
    )
    np.testing.assert_allclose(
        structured.covariance_model.materialize() / scale**2,
        dense.covariance,
        rtol=1e-6,
        atol=1e-7,
    )


def test_sparse_positive_definite_factor_rejects_indefinite_matrix():
    with pytest.raises(np.linalg.LinAlgError, match="not positive definite"):
        factor_sparse_positive_definite(sparse.diags([-1.0, 2.0]))


def test_sparse_positive_definite_factor_rejects_roundoff_rank_deficiency():
    loading = np.asarray(
        [[-20.749325931793987], [167.02490282686168], [165.85872694090048]]
    )
    low_rank = loading @ loading.T
    represented = sparse.csc_matrix(np.diag(np.full(3, 1e-16)) + low_rank)
    # The positive diagonal is below one ULP at this scale and disappears.
    # Dense eigensolvers disagree on the sign of the resulting near-zero
    # eigenvalues, so assert the portable floating-point failure directly.
    np.testing.assert_array_equal(represented.toarray(), low_rank)
    with pytest.raises(np.linalg.LinAlgError, match="singular|positive definite"):
        factor_sparse_positive_definite(represented)


def test_conditional_eiv_rejects_pseudo_reml_objective():
    with pytest.raises(ValueError, match="no standard REML objective"):
        fit_conditional_eiv_gaussian(
            np.asarray([1.0, 2.0, 3.0]),
            np.column_stack([np.ones(3), np.arange(3.0)]),
            [np.eye(3) * 0.1],
            [1],
            np.eye(3) * 0.1,
            [("evolutionary", np.eye(3))],
            reml=True,
        )


def test_large_conditional_eiv_fit_is_rejected_before_covariance_validation(
    monkeypatch,
):
    monkeypatch.setattr(measurement_error_mod, "MAX_DENSE_EIV_OBSERVATIONS", 2)
    with pytest.raises(ValueError, match="limited to 2 observations"):
        fit_conditional_eiv_gaussian(
            np.arange(3.0),
            np.ones((3, 1)),
            [],
            [],
            np.eye(3),
            [("evolutionary", np.eye(3))],
            reml=False,
        )


def test_large_structured_conditional_eiv_fit_remains_available(monkeypatch):
    monkeypatch.setattr(measurement_error_mod, "MAX_DENSE_EIV_OBSERVATIONS", 2)
    fit = fit_conditional_eiv_gaussian(
        np.asarray([1.0, 2.0, 3.2]),
        np.column_stack([np.ones(3), np.arange(3.0)]),
        [
            ContinuousPredictorUncertainty(
                np.eye(3) * 0.1,
                np.arange(3),
            )
        ],
        [1],
        np.zeros(3),
        [("evolutionary", np.ones(3))],
        reml=False,
    )
    assert np.isfinite(fit["objective"])


def test_diagonal_predictor_factor_matches_dense_eiv_covariance():
    response = np.asarray([1.0, 2.0, 3.2, 4.1])
    design = np.column_stack([np.ones(4), np.arange(4.0)])
    diagonal_factor = np.asarray([0.1, 0.2, 0.15, 0.25])
    dense = fit_conditional_eiv_gaussian(
        response,
        design,
        [np.diag(np.square(diagonal_factor))],
        [1],
        np.full(4, 0.1),
        [("evolutionary", np.ones(4))],
        reml=False,
    )
    structured = fit_conditional_eiv_gaussian(
        response,
        design,
        [ContinuousPredictorUncertainty(diagonal_factor, np.arange(4))],
        [1],
        np.full(4, 0.1),
        [("evolutionary", np.ones(4))],
        reml=False,
    )
    np.testing.assert_allclose(structured["beta"], dense["beta"], rtol=3e-5)
    assert structured["objective"] == pytest.approx(dense["objective"], rel=1e-6)
    assert structured["covariance"].low_rank.shape == (4, 0)

    # This fixture's residual variation is below its known sampling error,
    # so the additional variance has its optimum at zero. Independently
    # profile the intercept and optimize the remaining one-dimensional
    # likelihood: agreement between two prematurely stopped fits is not enough.
    def boundary_objective(slope):
        variance = 0.1 + np.square(slope * diagonal_factor)
        intercept = np.sum((response - slope * design[:, 1]) / variance) / np.sum(
            1.0 / variance
        )
        residual = response - intercept - slope * design[:, 1]
        return 0.5 * np.sum(
            np.log(2.0 * np.pi * variance) + np.square(residual) / variance
        )

    reference = minimize_scalar(
        boundary_objective,
        bounds=(0.5, 1.5),
        method="bounded",
        options={"xatol": 1e-14},
    )
    assert reference.success
    for fit in (dense, structured):
        assert fit["objective"] == pytest.approx(reference.fun, abs=1e-8, rel=0.0)
        assert fit["boundary_warning"]
        assert fit["component_variances"]["evolutionary"] == pytest.approx(
            fit["lower_variance"], rel=1e-12, abs=0.0
        )
        # Independent coefficient information for the active zero-variance
        # model, including covariance dependence on the slope.
        intercept, slope = fit["beta"]
        residual = response - intercept - slope * design[:, 1]
        variance = 0.1 + np.square(slope * diagonal_factor)
        derivative = 2.0 * slope * np.square(diagonal_factor)
        h00 = np.sum(1.0 / variance)
        h01 = np.sum(design[:, 1] / variance + residual * derivative / variance**2)
        h11 = np.sum(
            design[:, 1] ** 2 / variance
            + 2.0 * design[:, 1] * residual * derivative / variance**2
            + np.square(diagonal_factor)
            * (1.0 / variance - np.square(residual) / variance**2)
            + np.square(derivative)
            * (-0.5 / variance**2 + np.square(residual) / variance**3)
        )
        expected_covariance = np.linalg.inv([[h00, h01], [h01, h11]])
        np.testing.assert_allclose(
            fit["beta_covariance"], expected_covariance, rtol=1e-5, atol=1e-8
        )


def test_dense_eiv_preserves_correlated_predictor_loading_with_repeated_tips():
    response = np.array([1.1, 0.8, 2.6, 3.3, 4.9, 4.2])
    predictor = np.array([0.0, 0.0, 1.0, 2.0, 3.0, 3.0])
    design = np.column_stack((np.ones(len(response)), predictor))
    loading = np.array([[0.3, 0.1], [0.2, 0.4], [0.0, 0.1], [0.3, 0.5]])
    indices = np.array([0, 0, 1, 2, 3, 3])
    scales = np.array([1.0, -1.5, 0.5, 2.0, 1.0, 0.75])
    selected = loading[indices] * scales[:, None]
    fixed = np.diag(np.linspace(0.03, 0.09, len(response)))
    components = [("evolutionary", np.ones(len(response)))]
    fit = fit_conditional_eiv_gaussian(
        response,
        design,
        [ContinuousPredictorUncertainty(loading, indices, scales)],
        [1],
        fixed,
        components,
        reml=False,
    )
    reference = fit_conditional_eiv_gaussian(
        response,
        design,
        [selected @ selected.T],
        [1],
        fixed,
        components,
        reml=False,
    )
    np.testing.assert_allclose(fit["beta"], reference["beta"], rtol=1e-5)
    assert fit["objective"] == pytest.approx(reference["objective"], abs=1e-8)
    for slope in (0.0, 1.2, -0.4):
        covariance, _ = fit["covariance_for_beta"](np.array([0.0, slope]))
        expected = fixed + slope**2 * (selected @ selected.T)
        expected += np.eye(len(response)) * fit["component_variances"]["evolutionary"]
        np.testing.assert_allclose(covariance, expected, rtol=2e-14, atol=1e-15)


def test_sparse_precision_predictor_uncertainty_matches_dense_eiv_covariance():
    response = np.asarray([1.0, 2.0, 3.2, 4.1])
    design = np.column_stack([np.ones(4), np.arange(4.0)])
    loading = sparse.csr_matrix(
        np.asarray(
            [
                [0.2, 0.0],
                [0.1, 0.2],
                [0.0, 0.3],
                [0.2, 0.1],
            ]
        )
    )
    precision = sparse.csc_matrix(np.asarray([[2.0, -0.3], [-0.3, 1.5]]))
    precision_factor = np.linalg.cholesky(precision.toarray()).T
    model = measurement_error_mod.SparseCovarianceModel(
        precision=precision,
        tip_loading=loading,
        logdet_covariance=-np.linalg.slogdet(precision.toarray())[1],
        sampling_parent=np.empty(0, dtype=int),
        sampling_transition=np.empty(0),
        sampling_variance=np.empty(0),
        sampling_precision_factor=sparse.csr_matrix(precision_factor),
    )
    dense_uncertainty = model.materialize()
    common = dict(
        response=response,
        design=design,
        predictor_columns=[1],
        fixed_covariance=np.full(4, 0.1),
        components=[("evolutionary", np.ones(4))],
        reml=False,
    )
    dense = fit_conditional_eiv_gaussian(
        predictor_uncertainties=[dense_uncertainty], **common
    )
    structured = fit_conditional_eiv_gaussian(
        predictor_uncertainties=[GmrfPredictorUncertainty(model, np.arange(4))],
        **common,
    )

    np.testing.assert_allclose(structured["beta"], dense["beta"], rtol=5e-5)
    assert structured["objective"] == pytest.approx(
        dense["objective"], rel=2e-6, abs=2e-5
    )
    assert isinstance(structured["covariance"], DiagonalSparsePrecisionCovariance)
    np.testing.assert_allclose(
        materialize_covariance(structured["covariance"]),
        dense["covariance"],
        rtol=5e-5,
        atol=1e-7,
    )


def test_grouped_sparse_precision_eiv_objective_matches_dense_covariance():
    response = np.asarray([1.0, 1.9, 3.2, 4.0])
    design = np.column_stack([np.ones(4), np.arange(4.0)])
    loading = sparse.eye(4, format="csr")
    precision = sparse.csc_matrix(
        np.asarray(
            [
                [2.0, -0.3, 0.0, 0.0],
                [-0.3, 1.8, -0.2, 0.0],
                [0.0, -0.2, 1.7, -0.25],
                [0.0, 0.0, -0.25, 1.6],
            ]
        )
    )
    precision_factor = np.linalg.cholesky(precision.toarray()).T
    model = measurement_error_mod.SparseCovarianceModel(
        precision=precision,
        tip_loading=loading,
        logdet_covariance=-np.linalg.slogdet(precision.toarray())[1],
        sampling_parent=np.empty(0, dtype=int),
        sampling_transition=np.empty(0),
        sampling_variance=np.empty(0),
        sampling_precision_factor=sparse.csr_matrix(precision_factor),
    )
    common = dict(
        response=response,
        design=design,
        predictor_columns=[1],
        fixed_covariance=np.full(4, 0.1),
        components=[("evolutionary", np.ones(4))],
        reml=False,
        likelihood_observations=2,
        likelihood_groups=np.asarray([0, 0, 1, 1]),
    )

    dense = fit_conditional_eiv_gaussian(
        predictor_uncertainties=[model.materialize()], **common
    )
    structured = fit_conditional_eiv_gaussian(
        predictor_uncertainties=[GmrfPredictorUncertainty(model, np.arange(4))],
        **common,
    )

    np.testing.assert_allclose(structured["beta"], dense["beta"], rtol=5e-5)
    assert structured["objective"] == pytest.approx(
        dense["objective"], rel=2e-6, abs=2e-5
    )


def test_sparse_precision_covariance_draw_uses_retained_precision_factor():
    diagonal = np.asarray([0.2, 0.3, 0.4])
    loading = sparse.csr_matrix(np.asarray([[0.2, 0.0], [0.1, 0.2], [0.0, 0.3]]))
    precision = np.asarray([[2.0, -0.3], [-0.3, 1.5]])
    precision_factor = np.linalg.cholesky(precision).T
    factor = factor_diagonal_sparse_precision_updates(
        diagonal,
        [(loading, sparse.csc_matrix(precision), sparse.csr_matrix(precision_factor))],
    )
    observation_noise = np.asarray([0.4, -0.2, 0.1])
    latent_noise = np.asarray([0.3, -0.5])

    class FixedGenerator:
        def standard_normal(self, size):
            assert size == len(latent_noise)
            return latent_noise

    actual = draw_from_factor(
        factor,
        observation_noise,
        rng=FixedGenerator(),
    )
    latent_perturbation = precision_factor.T @ latent_noise
    expected = np.sqrt(diagonal) * observation_noise + loading @ np.linalg.solve(
        precision, latent_perturbation
    )
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)


def test_sparse_precision_marginal_diagonal_is_exact_above_legacy_threshold():
    block = sparse.csc_matrix(np.asarray([[2.0, -0.8], [-0.8, 1.5]]))
    precision = sparse.block_diag([block] * 256 + [sparse.csc_matrix([[1.2]])])
    loading = sparse.eye(513, format="csr")
    block_inverse = np.linalg.inv(block.toarray())
    expected = np.concatenate(
        [np.tile(np.diag(block_inverse), 256), np.asarray([1.0 / 1.2])]
    )

    actual = sparse_precision_update_diagonal(loading, precision)

    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)


def test_grouped_gmrf_marginals_are_precomputed_once_for_optimization(monkeypatch):
    response = np.asarray([1.0, 2.0, 3.2, 4.1])
    design = np.column_stack([np.ones(4), np.arange(4.0)])
    loading = sparse.eye(4, format="csr")
    precision = sparse.csc_matrix(
        np.asarray(
            [
                [2.0, -0.3, 0.0, 0.0],
                [-0.3, 1.8, -0.2, 0.0],
                [0.0, -0.2, 1.7, -0.25],
                [0.0, 0.0, -0.25, 1.6],
            ]
        )
    )
    model = measurement_error_mod.SparseCovarianceModel(
        precision=precision,
        tip_loading=loading,
        logdet_covariance=-np.linalg.slogdet(precision.toarray())[1],
        sampling_parent=np.empty(0, dtype=int),
        sampling_transition=np.empty(0),
        sampling_variance=np.empty(0),
        sampling_precision_factor=sparse.csr_matrix(
            np.linalg.cholesky(precision.toarray()).T
        ),
    )
    original = measurement_error_mod.sparse_precision_update_diagonal
    calls = 0

    def count_calls(*args, **kwargs):
        nonlocal calls
        calls += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(
        measurement_error_mod, "sparse_precision_update_diagonal", count_calls
    )
    fit_conditional_eiv_gaussian(
        response,
        design,
        [GmrfPredictorUncertainty(model, np.arange(4))],
        [1],
        np.full(4, 0.1),
        [("evolutionary", np.ones(4))],
        reml=False,
        likelihood_observations=2,
        likelihood_groups=np.asarray([0, 0, 1, 1]),
    )

    # One grouped profile before optimization and one row profile retained in
    # the returned diagnostic covariance; never one exact solve per objective.
    assert calls == 2


def test_structured_eiv_warns_and_attempts_above_validated_size(monkeypatch):
    monkeypatch.setattr(measurement_error_mod, "MAX_STRUCTURED_EIV_OBSERVATIONS", 2)
    with pytest.warns(RuntimeWarning, match="outside the routine validation range"):
        fit = fit_conditional_eiv_gaussian(
            np.arange(3.0),
            np.ones((3, 1)),
            [ContinuousPredictorUncertainty(np.ones(3), np.arange(3))],
            [0],
            np.ones(3),
            [("evolutionary", np.ones(3))],
            reml=False,
        )
    assert fit["optimizer_converged"]
