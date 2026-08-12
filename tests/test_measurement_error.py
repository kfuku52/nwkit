import numpy as np
import pytest

from nwkit.measurement_error import fit_latent_predictor


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
