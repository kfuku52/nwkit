"""Independent likelihood, derivative, root/boundary and dispatch checks."""

from dataclasses import replace

import numpy as np
import pytest
from ete4 import Tree
from scipy.optimize._numdiff import approx_derivative

from nwkit.shift_joint_dense import DenseJointContext, dense_eligible, design_derivative
from nwkit.shift_joint_fit import _correlation, _parameter_score
from nwkit.shift_joint_model import evaluate_joint
from nwkit.shift_native_model import ShiftData, ShiftLayout
from tests.test_shift_joint_covariance import dense_ou, fixture_data


@pytest.mark.parametrize("root", ["OUfixedRoot", "OUrandomRoot"])
@pytest.mark.parametrize(
    "alpha", [[0.4, 2.0], [1e-8, 1e-6], [700, 1000], [np.inf, np.inf]]
)
@pytest.mark.parametrize("missing", [False, True])
def test_dense_matches_branch_recursion_and_pruning(root, alpha, missing):
    data = fixture_data(0.02, missing)
    layout = ShiftLayout.build(data.tree, [2])
    S = np.array([[0.6, -0.2], [-0.2, 0.8]])
    noise = np.array([0.05, 0.03])
    alpha = np.asarray(alpha)
    context = DenseJointContext(data, root)
    fit = context.evaluate(layout, alpha, S, noise)
    reference = evaluate_joint(data, layout, alpha, S, noise, root_model=root)
    factor = context.covariance(alpha, S, noise)[0].lower
    expected = dense_ou(data, alpha, S, root) + np.diag(
        (data.variances + noise).ravel()
    )
    expected = expected[np.ix_(context.mask.ravel(), context.mask.ravel())]
    np.testing.assert_allclose(factor @ factor.T, expected, atol=1e-10, rtol=1e-10)
    assert fit.log_likelihood == pytest.approx(reference.log_likelihood, abs=1e-9)
    for field in (
        "coefficients",
        "coefficient_covariance",
        "predicted",
        "process_tip_covariance",
    ):
        np.testing.assert_allclose(
            getattr(fit, field), getattr(reference, field), atol=1e-9
        )


@pytest.mark.parametrize("root", ["OUfixedRoot", "OUrandomRoot"])
@pytest.mark.parametrize("shared", [False, True])
@pytest.mark.parametrize("fixed_process", [False, True])
def test_optimizer_coordinate_gradient(root, shared, fixed_process):
    data = fixture_data(0.02, True)
    layout = ShiftLayout.build(data.tree, [2])
    context = DenseJointContext(data, root)
    ac, pc = (1 if shared else 2), (0 if fixed_process else 2)
    fixed = np.array([0.5, 0.8]) if fixed_process else None
    x = np.array(
        ([np.log(0.7)] if shared else [np.log(0.4), np.log(2)])
        + ([] if fixed_process else [0.7, 0.9])
        + [-0.3, 0.2, 0.3]
    )

    def evaluate(v):
        alpha = np.broadcast_to(np.exp(v[:ac]), (2,))
        std = v[ac : ac + pc] if pc else np.sqrt(fixed)
        corr = _correlation(v[ac + pc : ac + pc + 1], 2, False)
        fit, score = context.evaluate(
            layout, alpha, std[:, None] * corr * std[None], v[-2:], gradient=True
        )
        return -fit.log_likelihood, _parameter_score(v, score, 2, ac, pc, 1, 2, fixed)

    _, analytic = evaluate(x)
    numerical = approx_derivative(lambda v: evaluate(v)[0], x, method="3-point")
    np.testing.assert_allclose(analytic, numerical, atol=2e-7, rtol=2e-6)


@pytest.mark.parametrize(
    "alpha", [[0, 0], [0, 2], [1e-7, 2], [600, 900], [np.inf, np.inf]]
)
def test_physical_scores_at_boundaries(alpha):
    data = fixture_data(0.04, True)
    layout = ShiftLayout.build(data.tree, [2])
    context = DenseJointContext(data, "OUfixedRoot")
    S, noise = np.array([[0.6, 0.2], [0.2, 0.8]]), np.array([0.0, 0.03])
    alpha = np.array(alpha, float)
    fit, (a, s, e) = context.evaluate(layout, alpha, S, noise, gradient=True)
    for j in range(2):
        h = 1e-5
        changed = noise.copy()
        changed[j] += h
        right = -context.evaluate(layout, alpha, S, changed).log_likelihood
        changed[j] += h
        second = -context.evaluate(layout, alpha, S, changed).log_likelihood
        assert e[j] == pytest.approx(
            (-3 * -fit.log_likelihood + 4 * right - second) / (2 * h), abs=2e-6
        )
        if np.isfinite(alpha).all():
            changed = alpha.copy()
            changed[j] += h
            right = -context.evaluate(layout, changed, S, noise).log_likelihood
            changed[j] += h
            second = -context.evaluate(layout, changed, S, noise).log_likelihood
            assert a[j] == pytest.approx(
                (-3 * -fit.log_likelihood + 4 * right - second) / (2 * h), abs=2e-6
            )
    # A symmetric off-diagonal perturbation has both coordinate scores.
    delta = np.array([[0.0, 1.0], [1.0, 0.0]]) * 1e-5
    derivative = (
        -context.evaluate(layout, alpha, S + delta, noise).log_likelihood
        + context.evaluate(layout, alpha, S - delta, noise).log_likelihood
    ) / 2e-5
    assert s[0, 1] + s[1, 0] == pytest.approx(derivative, abs=2e-6)


def test_rounded_tree_and_wholly_missing_tip():
    tree = Tree(
        "(((a:0.4,b:0.400000001):0.3,(c:0.4,d:0.4):0.3):0.3,((e:0.4,f:0.4):0.3,(g:0.4,h:0.4):0.3):0.3);"
    )
    values = np.random.default_rng(29).normal(size=(8, 3))
    values[0] = np.nan
    data = ShiftData.build(tree, values, ("a", "b", "c"), np.full((8, 3), 0.1))
    layout = ShiftLayout.build(data.tree, [2])
    context = DenseJointContext(data, "OUrandomRoot")
    alpha = np.array([0.2, 0.8, 4.0])
    S = np.eye(3) + 0.3
    result = context.evaluate(layout, alpha, S, np.zeros(3))
    ref = evaluate_joint(data, layout, alpha, S, np.zeros(3), root_model="OUrandomRoot")
    assert result.log_likelihood == pytest.approx(ref.log_likelihood, abs=1e-9)
    numeric = (
        layout.design(data.tree, 0.8 + 1e-5) - layout.design(data.tree, 0.8 - 1e-5)
    ) / 2e-5
    np.testing.assert_allclose(
        design_derivative(data.tree, layout, 0.8), numeric, atol=1e-9
    )


def test_dense_preserves_rejections_and_bounds():
    data = fixture_data(0.02, True)
    layout = ShiftLayout.build(data.tree)
    context = DenseJointContext(data, "OUfixedRoot")
    for S, noise, alpha in [
        ([[1, 2], [2, 1]], [0, 0], [1, 1]),
        (np.eye(2), [-1, 0], [1, 1]),
        (np.eye(2), [0, 0], [-1, 1]),
    ]:
        with pytest.raises(ValueError):
            context.evaluate(layout, alpha, S, noise)
    with pytest.raises(ValueError, match="Stationary"):
        DenseJointContext(data, "OUrandomRoot").evaluate(
            layout, [0, 0], np.eye(2), [0, 0]
        )
    with pytest.raises(ValueError, match="rank deficient"):
        context.evaluate(
            ShiftLayout.build(
                data.tree, [data.tree.branch_ids[data.tree.compiled.leaf_indices[1]]]
            ),
            [1, 1],
            np.eye(2),
            [0, 0],
        )
    large = replace(data, values=np.ones((400, 2)))
    assert not dense_eligible(large)
    with pytest.raises(ValueError, match="bound"):
        DenseJointContext(large, "OUfixedRoot")


def test_zero_noise_has_nonzero_variance_score():
    data = fixture_data(0.001, True)
    layout = ShiftLayout.build(data.tree)
    context = DenseJointContext(data, "OUfixedRoot")
    fit, scores = context.evaluate(
        layout, [1.0, 1.0], np.eye(2) * 0.01, [0.0, 0.0], gradient=True
    )
    x = np.array([0.1, 0.1, 0.0, 0.0])
    gradient = _parameter_score(x, scores, 2, 0, 2, 0, 2, None)
    assert np.all(gradient[-2:] < 0)  # zero error is not a stationary point
    for j in range(2):
        noise = np.zeros(2)
        noise[j] = 1e-8
        finite = (
            -context.evaluate(
                layout, [1.0, 1.0], np.eye(2) * 0.01, noise
            ).log_likelihood
            + fit.log_likelihood
        ) / 1e-8
        assert gradient[2 + j] == pytest.approx(finite, rel=1e-5)


@pytest.mark.parametrize("diagonal", [False, True])
def test_three_trait_correlation_chain(diagonal):
    data = fixture_data(0.01)
    values = np.random.default_rng(456).normal(size=(8, 3))
    values[1, 1] = np.nan
    data = ShiftData.build(data.tree, values, ("x", "y", "z"), np.full((8, 3), 0.02))
    context = DenseJointContext(data, "OUfixedRoot")
    layout = ShiftLayout.build(data.tree)
    cc = 0 if diagonal else 3
    x = np.array([0.7, 0.8, 0.9] + ([] if diagonal else [0.2, -0.4, 0.3]))

    def fun(v):
        S = np.outer(v[:3], v[:3]) * _correlation(v[3:], 3, diagonal)
        fit, scores = context.evaluate(
            layout, [0.3, 0.7, 2.0], S, [0.01, 0.03, 0.02], gradient=True
        )
        return -fit.log_likelihood, _parameter_score(v, scores, 3, 0, 3, cc, 0, None)

    np.testing.assert_allclose(
        fun(x)[1], approx_derivative(lambda v: fun(v)[0], x), rtol=2e-6, atol=2e-7
    )


def test_stationarity_check_removes_step_bias_without_relaxing_tolerance():
    from nwkit.shift_joint_fit import _projected_gradient, _richardson_gradient

    # Variance profile with a known optimum at a small, positive variance.
    def objective(x):
        return 80 * (np.log(x[0]) + 0.005 / x[0]) + x[1] + x[2] ** 2

    x = np.array([0.005, 0.0, 1.0])
    bounds = [(0.0, 1.0), (0.0, 1.0), (0.0, 1.0)]
    np.testing.assert_allclose(
        _richardson_gradient(objective, x, bounds), [0.0, 1.0, 2.0], atol=2e-6
    )
    assert _projected_gradient(objective, x, bounds, accurate=True) == pytest.approx(
        2.0, abs=2e-6
    )
    x[2] = 0.0
    assert _projected_gradient(objective, x, bounds, accurate=True) < 1e-5
    assert _projected_gradient(objective, x, bounds) > 0.001


def test_zero_process_noise_only_and_singular_observations():
    data = fixture_data()
    layout = ShiftLayout.build(data.tree)
    dense = DenseJointContext(data, "OUfixedRoot")
    fit = dense.evaluate(layout, [1.0, 2.0], np.zeros((2, 2)), [0.2, 0.3])
    ref = evaluate_joint(data, layout, [1.0, 2.0], np.zeros((2, 2)), [0.2, 0.3])
    assert fit.log_likelihood == pytest.approx(ref.log_likelihood, abs=1e-10)
    with pytest.raises((ValueError, np.linalg.LinAlgError)):
        dense.evaluate(layout, [1.0, 2.0], np.zeros((2, 2)), [0.0, 0.0])


@pytest.mark.parametrize("convergent", [False, True])
def test_design_score_with_nested_and_reused_regimes(convergent):
    data = fixture_data(0.02, True)
    nodes = [
        i for i, children in enumerate(data.tree.compiled.children) if i and children
    ]
    first = nodes[0]
    second = next(i for i in nodes if data.tree.compiled.parents[i] == first)
    branches = [data.tree.branch_ids[i] for i in (first, second)]
    groups = ((0, branches[1]), (branches[0],)) if convergent else None
    layout = ShiftLayout.build(data.tree, branches, groups=groups)
    for alpha in (0.0001, 0.7, 30.0):
        numerical = (
            layout.design(data.tree, alpha + 1e-6)
            - layout.design(data.tree, alpha - 1e-6)
        ) / 2e-6
        np.testing.assert_allclose(
            design_derivative(data.tree, layout, alpha), numerical, atol=1e-9
        )


@pytest.mark.parametrize("engine", ["auto", "pruning"])
def test_other_engines_do_not_allocate_dense_geometry(monkeypatch, engine):
    import nwkit.shift_joint_fit as module
    from nwkit.shift_native_fit import NativeFitOptions, fit_native_layout

    def reject(*args, **kwargs):
        raise AssertionError("Unused dense geometry must not be allocated")

    monkeypatch.setattr(module, "DenseJointContext", reject)
    data = fixture_data()
    result = fit_native_layout(
        data,
        ShiftLayout.build(data.tree),
        options=NativeFitOptions(
            trait_covariance="full", alpha_model="shared", joint_engine=engine
        ),
        alpha_height=0.7,
    )
    assert np.isfinite(result["log_likelihood"])
