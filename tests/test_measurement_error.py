import numpy as np
import pytest
from ete4 import Tree
from scipy import sparse

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
    factor_sparse_positive_definite,
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


def test_latent_predictor_rejects_indefinite_sampling_covariance():
    with pytest.raises(ValueError, match="positive semidefinite"):
        fit_latent_predictor(
            np.asarray([1.0, 2.0]),
            np.eye(2),
            np.asarray([[1.0, 2.0], [2.0, 1.0]]),
            include_intercept=False,
        )


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


def test_sparse_latent_predictor_supports_low_rank_batch_sampling_covariance():
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
    sparse_fit = fit_sparse_latent_predictor(
        observed,
        build_sparse_evolutionary_model(tree, names),
        sampling,
        include_intercept=True,
    )

    assert sparse_fit.evolutionary_rate == pytest.approx(
        dense.evolutionary_rate, rel=1e-6
    )
    assert sparse_fit.log_likelihood == pytest.approx(dense.log_likelihood, rel=1e-8)
    np.testing.assert_allclose(sparse_fit.mean, dense.mean, rtol=1e-6, atol=1e-7)
    np.testing.assert_allclose(
        sparse_fit.covariance_model.materialize(),
        dense.covariance,
        rtol=1e-6,
        atol=1e-7,
    )


def test_factor_latent_predictor_matches_dense_conditioning_without_materializing():
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
    structured = fit_factor_latent_predictor(
        observed,
        evolutionary,
        sampling,
        include_intercept=False,
    )

    assert structured.evolutionary_rate == pytest.approx(
        dense.evolutionary_rate, rel=2e-6
    )
    assert structured.log_likelihood == pytest.approx(dense.log_likelihood, rel=1e-8)
    np.testing.assert_allclose(structured.mean, dense.mean, rtol=1e-6, atol=1e-7)
    np.testing.assert_allclose(
        structured.covariance_model.materialize(),
        dense.covariance,
        rtol=1e-6,
        atol=1e-7,
    )


def test_sparse_positive_definite_factor_rejects_indefinite_matrix():
    with pytest.raises(np.linalg.LinAlgError, match="not positive definite"):
        factor_sparse_positive_definite(sparse.diags([-1.0, 2.0]))


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
