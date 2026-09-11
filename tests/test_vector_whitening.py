"""Independent dense recursion verifies square-root vector pruning."""

import numpy as np
import pytest
from ete4 import Tree

from nwkit.compiled_tree import CompiledTree
from nwkit.vector_whitening import VectorTreeWhitening, vector_gls


def dense_covariance(compiled, slopes, innovations, root):
    n, p = slopes.shape
    blocks = np.zeros((n, n, p, p))
    blocks[0, 0] = root
    for i in range(1, n):
        parent = compiled.parents[i]
        for j in range(i):
            blocks[i, j] = slopes[i, :, None] * blocks[parent, j]
            blocks[j, i] = blocks[i, j].T
        blocks[i, i] = (
            slopes[i, :, None] * blocks[parent, parent] * slopes[i, None, :]
            + innovations[i]
        )
    return np.block(
        [[blocks[i, j] for j in compiled.leaf_indices] for i in compiled.leaf_indices]
    )


@pytest.mark.parametrize("missing", [False, True])
@pytest.mark.parametrize("noise", [0.0, 0.2])
@pytest.mark.parametrize("random_root", [False, True])
def test_joint_whitening_matches_dense(missing, noise, random_root):
    tree = CompiledTree.from_tree(Tree("((a:1,b:1):1,(c:1,(d:0.5,e:0.5):0.5):1);"))
    rng = np.random.default_rng(4)
    p, n = 3, len(tree.nodes)
    slopes = rng.uniform(0.3, 0.95, (n, p))
    matrices = rng.normal(size=(n, p, p))
    innovations = matrices @ np.swapaxes(matrices, 1, 2) + np.eye(p) * 0.2
    root = innovations[0] if random_root else np.zeros((p, p))
    errors = np.broadcast_to(noise * (np.eye(p) + np.ones((p, p))), (5, p, p))
    mask = np.ones((5, p), dtype=bool)
    if missing:
        mask[1] = False
        mask[3, 1:] = False
    dense = dense_covariance(tree, slopes, innovations, root)
    for i in range(5):
        dense[i * p : (i + 1) * p, i * p : (i + 1) * p] += errors[i]
    dense = dense[np.ix_(mask.ravel(), mask.ravel())]
    factor = VectorTreeWhitening.build(
        tree, mask, slopes, innovations, errors, root_covariance=root
    )
    W = factor.apply(np.eye(mask.sum()))
    np.testing.assert_allclose(W @ dense @ W.T, np.eye(mask.sum()), atol=1e-12)
    assert factor.log_determinant == pytest.approx(
        np.linalg.slogdet(dense)[1], abs=1e-12
    )
    y = rng.normal(size=mask.sum())
    x = np.column_stack((np.ones(len(y)), rng.normal(size=len(y))))
    beta, ll, _, covariance = vector_gls(factor, y, x)
    precision_x = np.linalg.solve(dense, x)
    expected_cov = np.linalg.inv(x.T @ precision_x)
    expected_beta = expected_cov @ precision_x.T @ y
    residual = y - x @ expected_beta
    expected_ll = -0.5 * (
        len(y) * np.log(2 * np.pi)
        + np.linalg.slogdet(dense)[1]
        + residual @ np.linalg.solve(dense, residual)
    )
    np.testing.assert_allclose(beta, expected_beta, atol=1e-12)
    np.testing.assert_allclose(covariance, expected_cov, atol=1e-12)
    assert ll == pytest.approx(expected_ll, abs=1e-12)


def test_zero_process_with_noise_and_singular_rejection():
    tree = CompiledTree.from_tree(Tree("((a:1,b:1):1,c:2);"))
    p, n = 2, len(tree.nodes)
    args = (tree, np.ones((3, p), bool), np.ones((n, p)), np.zeros((n, p, p)))
    factor = VectorTreeWhitening.build(
        *args, np.tile(np.eye(p), (3, 1, 1)), root_covariance=np.zeros((p, p))
    )
    W = factor.apply(np.eye(6))
    np.testing.assert_allclose(W @ W.T, np.eye(6), atol=1e-12)
    with pytest.raises(ValueError, match="no jitter"):
        VectorTreeWhitening.build(
            *args, np.zeros((3, p, p)), root_covariance=np.zeros((p, p))
        )


def test_observation_plan_reuse_and_identity_checks():
    from nwkit.vector_whitening import VectorObservationPlan

    tree = CompiledTree.from_tree(Tree("((a:1,b:1):1,c:2);"))
    mask = np.array([[True, True], [False, False], [True, False]])
    plan = VectorObservationPlan.build(tree, mask)
    mask[0, 0] = False  # plan owns a stable copy
    assert plan.mask[0, 0]
    p, n = 2, len(tree.nodes)
    for scale in (0.3, 0.8):
        args = (
            tree,
            plan.mask,
            np.full((n, p), scale),
            np.tile(np.eye(p), (n, 1, 1)),
            np.tile(np.eye(p) * 0.1, (3, 1, 1)),
        )
        cached = VectorTreeWhitening.build(*args, root_covariance=np.eye(p), plan=plan)
        fresh = VectorTreeWhitening.build(*args, root_covariance=np.eye(p))
        np.testing.assert_array_equal(cached.apply(np.eye(3)), fresh.apply(np.eye(3)))
        assert cached.log_determinant == fresh.log_determinant
    with pytest.raises(ValueError, match="different tree or mask"):
        VectorTreeWhitening.build(
            tree, mask, *args[2:], root_covariance=np.eye(p), plan=plan
        )
