"""Independent joint OU likelihood and ML profile checks."""

import numpy as np
import pytest
from ete4 import Tree

from nwkit.shift_joint_model import (
    evaluate_joint,
    evaluate_separable,
    joint_covariance_geometry,
    joint_design,
)
from nwkit.shift_native_fit import NativeFitOptions, fit_native_layout
from nwkit.shift_native_ic import native_information_criterion
from nwkit.shift_native_model import ShiftData, ShiftLayout, evaluate_trait


def fixture_data(noise=0.0, missing=False):
    tree = Tree(
        "(((a:0.4,b:0.4):0.3,(c:0.4,d:0.4):0.3):0.3,((e:0.4,f:0.4):0.3,(g:0.4,h:0.4):0.3):0.3);"
    )
    values = np.random.default_rng(71).normal(size=(8, 2))
    if missing:
        values[1, 0] = np.nan
        values[3, 1] = np.nan
    return ShiftData.build(tree, values, ("x", "y"), np.full((8, 2), noise))


def dense_ou(data, alpha, S, root_model):
    # Independent branch-by-branch covariance recursion, not a pruning call.
    n, p = len(data.tree.branch_ids), len(alpha)
    if np.isinf(alpha).all():
        D = None
        root = S if root_model == "OUrandomRoot" else np.zeros_like(S)
    else:
        g = np.where(
            alpha == 0, 1.0, -np.expm1(-2 * alpha) / np.where(alpha == 0, 1, 2 * alpha)
        )
        if root_model == "OUrandomRoot":
            g = 1 / (2 * alpha)
        D = S / np.sqrt(g[:, None] * g[None])
        root = (
            D / (alpha[:, None] + alpha[None])
            if root_model == "OUrandomRoot"
            else np.zeros_like(S)
        )
    cov = np.zeros((n, n, p, p))
    cov[0, 0] = root
    for i in range(1, n):
        t = data.tree.times[i]
        slope = np.exp(-alpha * t)
        if D is None:
            innovation = S
        else:
            sums = alpha[:, None] + alpha[None]
            g = np.where(
                sums == 0, t, -np.expm1(-sums * t) / np.where(sums == 0, 1, sums)
            )
            innovation = D * g
        parent = data.tree.compiled.parents[i]
        for j in range(i):
            cov[i, j] = slope[:, None] * cov[parent, j]
            cov[j, i] = cov[i, j].T
        cov[i, i] = slope[:, None] * cov[parent, parent] * slope[None] + innovation
    leaves = data.tree.compiled.leaf_indices
    return np.block([[cov[i, j] for j in leaves] for i in leaves])


@pytest.mark.parametrize("alpha", [[0.4, 2.0], [0.0, 0.0], [np.inf, np.inf]])
@pytest.mark.parametrize("missing", [False, True])
def test_general_likelihood_dense_and_diagonal_nesting(alpha, missing):
    data = fixture_data(0.02, missing)
    layout = ShiftLayout.build(data.tree, [2])
    alpha = np.asarray(alpha)
    S = np.array([[0.6, 0.2], [0.2, 0.8]])
    noise = np.array([0.05, 0.03])
    fit = evaluate_joint(data, layout, alpha, S, noise)
    dense = dense_ou(data, alpha, S, "OUfixedRoot")
    dense += np.diag((data.variances + noise).ravel())
    mask = np.isfinite(data.values)
    dense = dense[np.ix_(mask.ravel(), mask.ravel())]
    X = joint_design(data, layout, alpha)[mask]
    y = data.values[mask]
    c = np.linalg.inv(X.T @ np.linalg.solve(dense, X))
    beta = c @ X.T @ np.linalg.solve(dense, y)
    residual = y - X @ beta
    ll = -0.5 * (
        len(y) * np.log(2 * np.pi)
        + np.linalg.slogdet(dense)[1]
        + residual @ np.linalg.solve(dense, residual)
    )
    assert fit.log_likelihood == pytest.approx(ll, abs=1e-10)
    np.testing.assert_allclose(fit.coefficients.T.ravel(), beta, atol=1e-10)
    independent = evaluate_joint(data, layout, alpha, np.diag(np.diag(S)), noise)
    scalars = [
        evaluate_trait(data, layout, j, alpha[j], S[j, j], noise[j]) for j in range(2)
    ]
    assert independent.log_likelihood == pytest.approx(
        sum(f.log_likelihood for f in scalars), abs=1e-10
    )


@pytest.mark.parametrize("root", ["OUfixedRoot", "OUrandomRoot"])
def test_separable_profile_equals_general(root):
    data = fixture_data()
    layout = ShiftLayout.build(data.tree, [2])
    profiled = evaluate_separable(data, layout, 0.7, root_model=root)
    general = evaluate_joint(
        data,
        layout,
        [0.7, 0.7],
        profiled.covariance_coordinate,
        [0.0, 0.0],
        root_model=root,
    )
    assert profiled.log_likelihood == pytest.approx(general.log_likelihood, abs=1e-10)
    np.testing.assert_allclose(profiled.coefficients, general.coefficients, atol=1e-10)
    np.testing.assert_allclose(
        profiled.coefficient_covariance, general.coefficient_covariance, atol=1e-10
    )


def test_shared_alpha_full_fit_and_count():
    data = fixture_data()
    layout = ShiftLayout.build(data.tree)
    options = NativeFitOptions(trait_covariance="full", alpha_model="shared")
    result = fit_native_layout(data, layout, options=options)
    assert result["trait_covariance"] == "full"
    assert result["joint_covariance"]["engine"] == "separable_profile"
    information = native_information_criterion(data, result, "AIC")
    assert (
        information["parameter_count"] == 6
    )  # 2 means + 3 covariance + 1 shared alpha
    assert all(r["log_likelihood"] is None for r in result["traits"])
    fixed = fit_native_layout(data, layout, options=options, alpha_height=0.7)
    assert native_information_criterion(data, fixed, "AIC")["parameter_count"] == 5
    with pytest.raises(ValueError, match="pBIC"):
        native_information_criterion(data, result, "pBIC")


def test_general_fixed_alpha_fit_with_noise():
    data = fixture_data(0.04, True)
    options = NativeFitOptions(trait_covariance="full", optimizer_starts=2)
    result = fit_native_layout(
        data, ShiftLayout.build(data.tree), options=options, alpha_height=[0.5, 2.0]
    )
    assert result["joint_covariance"]["engine"] == "dense_observed_gls"
    assert np.isfinite(result["log_likelihood"])
    assert result["joint_covariance"]["optimizer"]["complete_alpha_modes"]


def test_exact_profile_and_general_optimizer_agree():
    data = fixture_data()
    layout = ShiftLayout.build(data.tree, [2])
    exact = fit_native_layout(
        data,
        layout,
        options=NativeFitOptions(trait_covariance="full", alpha_model="shared"),
        alpha_height=0.7,
    )
    general = fit_native_layout(
        data,
        layout,
        options=NativeFitOptions(
            trait_covariance="full", alpha_model="shared", joint_engine="pruning"
        ),
        alpha_height=0.7,
    )
    assert exact["log_likelihood"] == pytest.approx(general["log_likelihood"], abs=1e-7)
    np.testing.assert_allclose(
        exact["joint_covariance"]["process_tip_covariance"],
        general["joint_covariance"]["process_tip_covariance"],
        rtol=1e-5,
        atol=1e-7,
    )


def test_joint_trait_permutation_and_units():
    data = fixture_data()
    layout = ShiftLayout.build(data.tree, [2])
    options = NativeFitOptions(trait_covariance="full", alpha_model="shared")
    original = fit_native_layout(data, layout, options=options, alpha_height=0.7)
    values = data.centers + data.values * data.scales
    changed = ShiftData.build(data.tree, values[:, ::-1] * [3.0, 0.2], ("y", "x"))
    permuted = fit_native_layout(changed, layout, options=options, alpha_height=0.7)
    assert permuted["log_likelihood"] == pytest.approx(
        original["log_likelihood"] - 8 * np.log(0.6), abs=1e-10
    )
    c = np.asarray(original["joint_covariance"]["process_tip_covariance"])[::-1, ::-1]
    np.testing.assert_allclose(
        permuted["joint_covariance"]["process_tip_covariance"],
        c * np.outer([3.0, 0.2], [3.0, 0.2]),
        atol=1e-10,
    )


def test_invalid_diffusion_and_mixed_infinite_alpha_rejected():
    data = fixture_data()
    with pytest.raises(ValueError, match="positive semidefinite"):
        joint_covariance_geometry(
            data.tree, [0.3, 2.0], [[1, 2], [2, 1]], "OUfixedRoot"
        )
    with pytest.raises(ValueError, match="mixed finite/infinite"):
        joint_covariance_geometry(
            data.tree, [0.3, np.inf], [[1.0, 0.2], [0.2, 1.0]], "OUfixedRoot"
        )


def test_unbounded_full_covariance_layout_is_excluded():
    from nwkit.shift_native_search import NativeLayoutEvaluator

    data = fixture_data()
    options = NativeFitOptions(trait_covariance="full", alpha_model="shared")
    # Eight regime means for eight observations have no residual covariance.
    leaves = [data.tree.branch_ids[i] for i in data.tree.compiled.leaf_indices]
    layout = ShiftLayout.build(data.tree, leaves[:-1])
    evaluator = NativeLayoutEvaluator(data, {"options": options, "alpha_height": 0.7})
    assert evaluator.evaluate(layout) == -np.inf
    assert (
        evaluator.records[-1]["status"]
        == "structurally_excluded_unbounded_covariance_ml"
    )
    with pytest.raises(ValueError, match="residual rank"):
        fit_native_layout(data, layout, options=options, alpha_height=0.7)
