"""Dense GLMM factors retain their solve and information-matrix contracts."""

import numpy as np
import pytest

from nwkit import phylogenetic_glmm as glmm


@pytest.mark.parametrize("categorical", [False, True])
def test_glmm_marginal_information_uses_existing_cholesky(categorical, monkeypatch):
    rng = np.random.default_rng(451)
    n, dimension = 8, 2 if categorical else 1
    design = np.column_stack([np.ones(n), rng.normal(size=n)])
    loading = rng.normal(size=(n * dimension, n * dimension))
    precision = loading @ loading.T + np.eye(n * dimension)
    if categorical:
        weights = np.repeat(np.array([[[0.3, -0.1], [-0.1, 0.2]]]), n, axis=0)
        expanded = np.kron(design, np.eye(dimension))
        weight_matrix = glmm._dense_weight_matrix(weights)
    else:
        weights = np.linspace(0.2, 0.8, n)
        expanded = design
        weight_matrix = np.diag(weights)
    fixed_random = expanded.T @ weight_matrix
    information = (
        expanded.T @ weight_matrix @ expanded
        - fixed_random @ np.linalg.solve(weight_matrix + precision, fixed_random.T)
    )
    expected = np.linalg.pinv((information + information.T) / 2, hermitian=True)

    def forbidden(*args, **kwargs):
        raise AssertionError(
            "A Cholesky factor must not be factored as a general matrix"
        )

    monkeypatch.setattr(np.linalg, "solve", forbidden)
    if categorical:
        actual = glmm._marginal_coefficient_covariance(
            design, dimension, weights, precision
        )
    else:
        actual = glmm._marginal_scalar_coefficient_covariance(
            design, np.arange(n), weights, precision
        )
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
