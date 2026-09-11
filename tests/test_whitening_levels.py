"""Independent dense checks on wide trees with partial and reordered observations."""

import numpy as np
import pytest

from nwkit.compiled_tree import CompiledTree
from nwkit.gaussian_whitening import TreeWhitening, tree_gls
from nwkit.util import read_tree


def balanced_names(count):
    nodes = [f"t{i}:1" for i in range(count)]
    while len(nodes) > 1:
        nodes = [f"({nodes[i]},{nodes[i + 1]}):1" for i in range(0, len(nodes), 2)]
    return nodes[0] + ";"


@pytest.mark.parametrize("root_variance", [0.0, 0.7])
@pytest.mark.parametrize("observations", [319, 448, 512])
def test_wide_tree_gls_matches_dense_with_signed_and_zero_slopes(
    root_variance, observations
):
    tree = read_tree(balanced_names(512), "auto", True, quiet=True)
    compiled = CompiledTree.from_tree(tree)
    rng = np.random.default_rng(1907)
    size = len(compiled.nodes)
    slopes = rng.uniform(-0.95, 0.95, size)
    slopes[::7] = 0
    variances = rng.uniform(0.1, 1.5, size)
    variances[0] = root_variance
    indices = tuple(rng.permutation(compiled.leaf_indices)[:observations])
    errors = rng.uniform(0, 0.2, observations)
    errors[::3] = 0
    loadings = np.eye(size)
    for i in range(1, size):
        loadings[i] += slopes[i] * loadings[compiled.parents[i]]
    observed_loadings = loadings[list(indices)]
    covariance = (observed_loadings * variances) @ observed_loadings.T + np.diag(errors)
    factor = TreeWhitening.build(
        compiled, indices, slopes, variances, errors, root_variance=root_variance
    )
    design = np.column_stack(
        [np.ones(observations), rng.normal(size=(observations, 4))]
    )
    response = rng.normal(size=observations)
    precision_design = np.linalg.solve(covariance, design)
    expected_covariance = np.linalg.inv(design.T @ precision_design)
    expected_beta = expected_covariance @ precision_design.T @ response
    residual = response - design @ expected_beta
    expected_quadratic = residual @ np.linalg.solve(covariance, residual)
    expected_determinant = np.linalg.slogdet(covariance)[1]
    beta, likelihood, quadratic, coefficient_covariance = tree_gls(
        factor, response, design
    )
    np.testing.assert_allclose(beta, expected_beta, atol=2e-12, rtol=2e-12)
    np.testing.assert_allclose(
        coefficient_covariance, expected_covariance, atol=2e-12, rtol=2e-12
    )
    assert quadratic == pytest.approx(expected_quadratic, abs=2e-11)
    assert factor.log_determinant == pytest.approx(expected_determinant, abs=2e-11)
    assert likelihood == pytest.approx(
        -0.5
        * (
            observations * np.log(2 * np.pi) + expected_determinant + expected_quadratic
        ),
        abs=2e-11,
    )
