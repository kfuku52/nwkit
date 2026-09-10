"""Compare tree elimination with an independently assembled dense covariance."""

import numpy as np
import pytest
from scipy.linalg import cho_solve

from nwkit.compiled_tree import CompiledTree
from nwkit.gaussian_whitening import TreeWhitening, tree_gls
from nwkit.util import read_tree

TREE = "(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,(g:1,h:1):1):1);"


def dense_covariance(compiled, slopes, innovations, root_variance, indices, errors):
    loadings = np.zeros((len(compiled.nodes), len(compiled.nodes)))
    loadings[0, 0] = 1
    for i in range(1, len(loadings)):
        loadings[i] = slopes[i] * loadings[compiled.parents[i]]
        loadings[i, i] = 1
    variances = np.array(innovations, copy=True)
    variances[0] = root_variance
    tips = loadings[list(indices)]
    return (tips * variances) @ tips.T + np.diag(errors)


@pytest.mark.parametrize("root_variance", [0, 0.7])
@pytest.mark.parametrize("slopes_kind", ["random", "independent", "brownian"])
@pytest.mark.parametrize("indices", [None, [14, 3, 10, 6], [3]])
def test_whitening_matches_dense_quadratics_and_determinant(
    root_variance, slopes_kind, indices
):
    compiled = CompiledTree.from_tree(read_tree(TREE, "auto", True, quiet=True))
    indices = tuple(compiled.leaf_indices if indices is None else indices)
    rng = np.random.default_rng(457)
    n = len(compiled.nodes)
    slopes = rng.uniform(0.01, 0.95, n)
    if slopes_kind != "random":
        slopes[:] = 0 if slopes_kind == "independent" else 1
    innovations = rng.uniform(0.05, 2, n)
    errors = rng.uniform(0, 0.2, len(indices))
    errors[::2] = 0
    values = rng.normal(size=(len(indices), 4))
    factor = TreeWhitening.build(
        compiled, indices, slopes, innovations, errors, root_variance=root_variance
    )
    covariance = dense_covariance(
        compiled, slopes, innovations, root_variance, indices, errors
    )
    expected = values.T @ np.linalg.solve(covariance, values)
    actual = factor.apply(values)
    np.testing.assert_allclose(actual.T @ actual, expected, rtol=2e-12, atol=2e-12)
    assert factor.log_determinant == pytest.approx(
        np.linalg.slogdet(covariance)[1], abs=2e-12
    )
    np.testing.assert_allclose(factor.apply(values[:, 0]), actual[:, 0])


def test_gls_matches_dense_mean_covariance_and_likelihood():
    compiled = CompiledTree.from_tree(read_tree(TREE, "auto", True, quiet=True))
    n = len(compiled.nodes)
    slopes, innovations = np.full(n, 0.9), np.full(n, 0.4)
    indices = compiled.leaf_indices
    errors = np.arange(8) / 100
    covariance = dense_covariance(compiled, slopes, innovations, 0.6, indices, errors)
    factor = TreeWhitening.build(
        compiled, indices, slopes, innovations, errors, root_variance=0.6
    )
    x = np.column_stack([np.ones(8), [0, 0, 0, 0, 1, 1, 1, 1], np.arange(8)])
    y = np.random.default_rng(12).normal(size=8)
    precision = cho_solve((np.linalg.cholesky(covariance), True), np.eye(8))
    expected_cov = np.linalg.inv(x.T @ precision @ x)
    expected_beta = expected_cov @ x.T @ precision @ y
    residual = y - x @ expected_beta
    expected_ll = -0.5 * (
        8 * np.log(2 * np.pi)
        + np.linalg.slogdet(covariance)[1]
        + residual @ precision @ residual
    )
    beta, likelihood, quadratic, beta_cov = tree_gls(factor, y, x)
    np.testing.assert_allclose(beta, expected_beta, atol=2e-12)
    np.testing.assert_allclose(beta_cov, expected_cov, atol=2e-12)
    assert likelihood == pytest.approx(expected_ll, abs=2e-12)
    assert quadratic == pytest.approx(residual @ precision @ residual, abs=2e-12)
    with pytest.raises(ValueError, match="rank deficient"):
        tree_gls(factor, y, np.column_stack([x[:, 0], x[:, 0]]))


def test_whitening_rejects_invalid_inputs():
    compiled = CompiledTree.from_tree(read_tree(TREE, "auto", True, quiet=True))
    n = len(compiled.nodes)
    with pytest.raises(ValueError, match="distinct"):
        TreeWhitening.build(compiled, [3, 3], np.ones(n), np.ones(n))
    with pytest.raises(ValueError, match="terminal variances"):
        TreeWhitening.build(compiled, [3], np.ones(n), np.zeros(n))
    with pytest.raises(ValueError, match="non-root tip"):
        TreeWhitening.build(compiled, [0], np.ones(n), np.ones(n))
    with pytest.raises(ValueError, match="nonnegative"):
        TreeWhitening.build(compiled, [3], np.ones(n), np.ones(n), [-1])
