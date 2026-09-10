"""Independent dense Gaussian checks of multivariate fixed-layout inference."""

import numpy as np
import pytest

from nwkit.shift_native_model import (
    ShiftData,
    ShiftLayout,
    ShiftTree,
    evaluate_trait,
    restore_trait_fit,
)
from nwkit.util import read_tree

TREE = "(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,(g:1,h:1):1):1);"


def dense_reference(tree, layout, alpha, variance, noise, root):
    """Propagate independent branch innovations and regime forcing directly."""
    compiled = tree.compiled
    n, p = len(compiled.nodes), len(layout.groups)
    loading = np.zeros((n, n))
    loading[0, 0] = 1
    design = np.zeros((n, p))
    design[:, 0] = 1
    group = {b: i for i, row in enumerate(layout.groups) for b in row}
    inherited = [0] * n
    innovation = np.zeros(n)
    if root == "OUrandomRoot":
        innovation[0] = variance
    for i, node in enumerate(compiled.nodes[1:], 1):
        parent = compiled.parents[i]
        time = float(node.dist) / tree.height
        inherited[i] = group.get(tree.branch_ids[i], inherited[parent])
        if alpha == 0:
            attenuation, forcing, innovation[i] = 1, time, variance * time
        elif np.isinf(alpha):
            attenuation, forcing, innovation[i] = 0, 1, variance
        else:
            attenuation = np.exp(-alpha * time)
            forcing = -np.expm1(-alpha * time) / -np.expm1(-alpha)
            innovation[i] = (
                variance * -np.expm1(-2 * alpha * time) / -np.expm1(-2 * alpha)
            )
            if root == "OUrandomRoot":
                innovation[i] = variance * -np.expm1(-2 * alpha * time)
        loading[i] = attenuation * loading[parent]
        loading[i, i] = 1
        design[i, 1:] = attenuation * design[parent, 1:]
        if inherited[i]:
            design[i, inherited[i]] += forcing
    rows = list(compiled.leaf_indices)
    loading = loading[rows]
    return design[rows], (loading * innovation) @ loading.T + np.diag(noise)


@pytest.mark.parametrize(
    "newick",
    [TREE.replace("a:1", "a:1.000000006"), "((((a:1,b:1):1,c:2):1,d:3):1,e:4);"],
)
@pytest.mark.parametrize("alpha", [0, 0.7, 1000, np.inf])
def test_interval_mean_design_matches_branch_propagation_on_unbalanced_and_rounded_trees(
    newick, alpha
):
    tree = ShiftTree.build(read_tree(newick, "auto", True, quiet=True))
    layout = ShiftLayout.build(tree, [1, 3])
    expected, _ = dense_reference(
        tree, layout, alpha, 1, np.zeros(len(tree.leaf_names)), "OUfixedRoot"
    )
    np.testing.assert_allclose(layout.design(tree, alpha), expected, atol=2e-14)


@pytest.mark.parametrize("alpha", [0, 1e-7, 0.6, 30, np.inf])
@pytest.mark.parametrize("root", ["OUfixedRoot", "OUrandomRoot"])
@pytest.mark.parametrize(
    "layout_spec", [([1], None), ([1, 7], [[0, 7], [1]]), ([7, 14], [[0], [7, 14]])]
)
def test_multitrait_layout_matches_independent_reference(alpha, root, layout_spec):
    tree = ShiftTree.build(read_tree(TREE, "auto", True, quiet=True))
    layout = ShiftLayout.build(tree, *layout_spec)
    values = np.random.default_rng(62).normal(size=(8, 2)) * [2, 5] + [3, -1]
    values[2, 0] = np.nan
    variances = np.tile([0.04, 0.12], (8, 1))
    variances[::2] = 0
    data = ShiftData.build(tree, values, ["first", "second"], variances)
    if alpha == 0 and root == "OUrandomRoot":
        with pytest.raises(ValueError, match="Stationary-root"):
            evaluate_trait(data, layout, 0, alpha, 0.8, root_model=root)
        return
    for trait in range(2):
        fit = evaluate_trait(data, layout, trait, alpha, 0.8, 0.07, root_model=root)
        design, covariance = dense_reference(
            tree, layout, alpha, 0.8, data.variances[:, trait] + 0.07, root
        )
        mask = np.isfinite(data.values[:, trait])
        x, y = design[mask], data.values[mask, trait]
        covariance = covariance[np.ix_(mask, mask)]
        whitening = np.linalg.cholesky(covariance)
        xw, yw = np.linalg.solve(whitening, x), np.linalg.solve(whitening, y)
        beta = np.linalg.lstsq(xw, yw, rcond=None)[0]
        residual = yw - xw @ beta
        likelihood = -0.5 * (
            len(y) * np.log(2 * np.pi)
            + 2 * np.log(np.diag(whitening)).sum()
            + residual @ residual
        )
        # Stationary-root alpha~0 has a large common covariance component.
        tolerance = 2e-8 if root == "OUrandomRoot" and alpha < 1e-6 else 2e-11
        np.testing.assert_allclose(
            fit.coefficients, beta, atol=tolerance, rtol=tolerance
        )
        assert fit.log_likelihood == pytest.approx(likelihood, abs=tolerance)
        restored = restore_trait_fit(data, trait, fit)
        np.testing.assert_allclose(
            restored["predicted"],
            data.centers[trait] + data.scales[trait] * (design @ beta),
            atol=1e-6,
        )
        assert restored["log_likelihood"] == pytest.approx(
            likelihood - len(y) * np.log(data.scales[trait]), abs=tolerance
        )


def test_time_and_trait_units_preserve_fixed_model_observables():
    tree = read_tree(TREE, "auto", True, quiet=True)
    y = np.random.default_rng(815).normal(size=(8, 1))
    data = ShiftData.build(tree, y, ["x"])
    layout = ShiftLayout.build(data.tree, [1])
    first = restore_trait_fit(data, 0, evaluate_trait(data, layout, 0, 0.6, 0.7))
    for node in tree.traverse():
        if not node.is_root:
            node.dist *= 1000
    other = ShiftData.build(tree, y * 1e-12, ["x"])
    second = restore_trait_fit(other, 0, evaluate_trait(other, layout, 0, 0.6, 0.7))
    np.testing.assert_allclose(
        np.array(second["predicted"]) / 1e-12, first["predicted"], atol=1e-12
    )
    assert second["alpha"] == pytest.approx(first["alpha"] / 1000)
    assert second["sigma2"] == pytest.approx(first["sigma2"] * 1e-24 / 1000)
    assert second["log_likelihood"] + 8 * np.log(1e-12) == pytest.approx(
        first["log_likelihood"]
    )


def test_rank_deficiency_and_missing_errors_are_explicit():
    tree = ShiftTree.build(read_tree(TREE, "auto", True, quiet=True))
    data = ShiftData.build(tree, np.arange(8)[:, None], ["x"])
    # Both root children shifted leave no extant background; alpha=infinity
    # then makes the baseline and the two regime columns linearly dependent.
    layout = ShiftLayout.build(tree, [1, 2])
    with pytest.raises(ValueError, match="rank deficient"):
        evaluate_trait(data, layout, 0, np.inf, 1)
    with pytest.raises(ValueError, match="partition"):
        ShiftLayout.build(tree, [1], [[0]])
    with pytest.raises(ValueError, match="immediate parent"):
        ShiftLayout.build(tree, [1], [[0, 1]])
    with pytest.raises(ValueError, match="constant"):
        ShiftData.build(tree, np.ones((8, 1)), ["x"])
