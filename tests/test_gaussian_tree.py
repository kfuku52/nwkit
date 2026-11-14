"""Cross-consumer contracts for the shared tree-Gaussian process layer."""

import math

import numpy as np
import pytest

from nwkit.continuous_asr import compute_bm_marginals
from nwkit.continuous_asr_process import fitted_scalar_process
from nwkit.evolution import (
    build_evolutionary_covariance,
    build_evolutionary_process,
    evolutionary_covariance_factory,
)
from nwkit.gaussian_tree import (
    GaussianRootPrior,
    GaussianTransition,
    GaussianTreeProcess,
    ou_transition,
)
from nwkit.nonstationary_continuous_asr import compute_eb_marginals
from nwkit.ordinary_regression import _fit_ordinary_gaussian
from nwkit.ou_asr import compute_ou_marginals
from nwkit.util import read_tree


def tree_from(source):
    return read_tree(source, "1", True, quiet=True, rooted="yes")


def test_ou_transition_preserves_finite_stationary_variance_at_extreme_alpha():
    transition, variance = ou_transition(2.0, 1e308, 1e308, 3.0)
    assert variance == 0.5
    assert transition == GaussianTransition(0.0, 3.0, 0.5)


def test_ou_innovation_and_intercept_survive_underflowing_attraction_time():
    transition, variance = ou_transition(1e-200, 1e-200, 1.0, 1e100)
    assert variance == pytest.approx(5e199)
    assert transition.slope == 1.0
    assert transition.variance == pytest.approx(1e-200, rel=1e-14, abs=0.0)
    assert transition.intercept == pytest.approx(1e-300, rel=1e-14, abs=0.0)


def test_reconstructed_stationary_ou_uses_fitted_variance():
    tree = tree_from("(A:1,B:2,C:3)R;")
    _, fit = compute_ou_marginals(
        tree, {"A": 0.5, "B": 0.6, "C": 0.2}, alpha=0.2, sigma2=1.4, theta=0.3
    )
    process = fitted_scalar_process(tree, "OU", fit)
    assert process.root.variance == fit.root_variance
    for leaf in tree.leaves():
        expected = fit.root_variance * -math.expm1(-2.0 * fit.alpha * leaf.dist)
        assert process.transitions[leaf].variance == expected


def test_brownian_asr_and_intercept_only_regression_share_reml_fit():
    tree = tree_from("((A:1,B:2)I:0.5,(C:3,D:0.25)J:0.4)R;")
    values = {"A": -1.0, "B": 2.0, "C": 4.0, "D": 1.3}
    names = [str(leaf.name) for leaf in tree.leaves()]
    posterior, asr_fit = compute_bm_marginals(tree, values)
    regression_fit = _fit_ordinary_gaussian(
        np.asarray([values[name] for name in names]),
        np.ones((len(names), 1)),
        np.zeros(len(names)),
        tree,
        names,
        evolution_model="brownian",
        evolution_parameter=None,
        branch_length="original",
        custom_covariance=None,
        reml=True,
    )
    assert regression_fit["component_variances"]["evolutionary_rate"] == pytest.approx(
        asr_fit.sigma2, rel=2e-12
    )
    assert regression_fit["objective"] == pytest.approx(
        -asr_fit.restricted_log_likelihood, abs=2e-12
    )
    assert regression_fit["beta"][0] == pytest.approx(posterior[tree].mean, abs=2e-12)
    assert regression_fit["beta_covariance"][0, 0] == pytest.approx(
        posterior[tree].variance, rel=2e-12
    )


def test_fixed_and_stationary_ou_roots_remain_explicitly_distinct():
    tree = tree_from("((A:0.4,B:0.7)I:0.3,(C:0.2,D:0.8)J:0.5)R;")
    names = [str(leaf.name) for leaf in tree.leaves()]
    alpha = 0.9
    fixed = build_evolutionary_process(tree, model="ou", parameter=alpha)
    stationary = build_evolutionary_process(
        tree,
        model="ou",
        parameter=alpha,
        root_mode="stationary",
        variance_scale=1.7,
    )
    assert fixed.root.mode == "fixed"
    assert stationary.root.mode == "stationary"
    np.testing.assert_allclose(
        fixed.tip_covariance(names),
        build_evolutionary_covariance(tree, names, model="ou", parameter=alpha),
        rtol=2e-13,
        atol=2e-13,
    )
    fixed_cross_root = fixed.tip_covariance(names)[0, 2]
    stationary_cross_root = stationary.tip_covariance(names)[0, 2]
    assert fixed_cross_root == 0.0
    assert stationary_cross_root > 0.0


def test_stationary_process_covariance_is_the_ou_asr_dense_oracle():
    tree = tree_from("((A:0.3,B:0.7)I:0.4,C:1.1,D:0.2)R;")
    values = {"A": 1.2, "B": -0.4, "C": 2.1, "D": None}
    errors = {"A": 0.2, "B": 0.0, "C": 0.4}
    alpha, sigma2, theta = 0.8, 1.7, 0.3
    process = build_evolutionary_process(
        tree,
        model="ou",
        parameter=alpha,
        root_mode="stationary",
        root_mean=theta,
        variance_scale=sigma2,
        allow_zero=True,
    )
    nodes = list(tree.traverse())
    means, _ = process.marginal_moments()
    covariance = process.covariance(nodes)
    observed = [
        node for node in tree.leaves() if values.get(str(node.name)) is not None
    ]
    observed_indices = [nodes.index(node) for node in observed]
    observed_covariance = covariance[np.ix_(observed_indices, observed_indices)]
    observed_covariance += np.diag([errors[str(node.name)] ** 2 for node in observed])
    response = np.asarray([values[str(node.name)] for node in observed])
    prior_mean = np.asarray([means[node] for node in nodes])
    centered = response - prior_mean[observed_indices]
    cross = covariance[:, observed_indices]
    conditional_mean = prior_mean + cross @ np.linalg.solve(
        observed_covariance, centered
    )
    conditional_covariance = covariance - cross @ np.linalg.solve(
        observed_covariance, cross.T
    )

    posterior, _ = compute_ou_marginals(
        tree,
        values,
        alpha=alpha,
        sigma2=sigma2,
        theta=theta,
        standard_errors=errors,
    )
    for index, node in enumerate(nodes):
        assert posterior[node].mean == pytest.approx(conditional_mean[index], abs=2e-12)
        assert posterior[node].variance == pytest.approx(
            max(0.0, conditional_covariance[index, index]), abs=2e-12
        )


@pytest.mark.parametrize("model,parameter", [("brownian", None), ("ou", 0.6)])
def test_batched_covariance_matches_pairs_with_mixed_and_repeated_nodes(
    model, parameter
):
    tree = tree_from("(((A:1,B:2)I:0.5,C:3)J:0.1,(D:0.4,E:2.5)K:0.3)R;")
    process = build_evolutionary_process(tree, model=model, parameter=parameter)
    distinct_nodes = list(tree.traverse(strategy="postorder"))
    cov = process.covariance(distinct_nodes)
    selected = [4, 3, 0, 8, 1, 2, 8, 6, 7, 5] * 4
    result = process.covariance([distinct_nodes[index] for index in selected])
    np.testing.assert_allclose(
        result, cov[np.ix_(selected, selected)], rtol=2e-15, atol=0.0
    )
    np.testing.assert_array_equal(result, result.T)


@pytest.mark.parametrize("slope,root_variance", [(1e-200, 1e300), (1e200, 1e-300)])
def test_gaussian_moments_keep_finite_variances_with_extreme_affine_slopes(
    slope, root_variance
):
    tree = tree_from("(A:1,B:1)R;")
    process = GaussianTreeProcess(
        tree,
        {
            tree["A"]: GaussianTransition(slope, 0.0, 0.0),
            tree["B"]: GaussianTransition(1.0, 0.0, 0.0),
        },
        GaussianRootPrior("gaussian", 0.0, root_variance),
        "test",
    )
    expected = slope * (slope * root_variance)
    expected_covariance = slope * root_variance
    with np.errstate(over="raise", under="raise", invalid="raise"):
        _, variances = process.marginal_moments()
        covariance = process.tip_covariance(["A", "B"])
    assert variances[tree["A"]] == expected
    assert covariance[0, 0] == expected
    assert covariance[0, 1] == pytest.approx(expected_covariance, rel=1e-12, abs=0.0)


def test_positive_asr_eb_uses_the_shared_acdc_branch_process():
    tree = tree_from("((A:0.4,B:1.2)I:0.7,C:1.5,D:0.8)R;")
    values = {"A": -1.0, "B": 2.0, "C": 4.0, "D": None}
    rate = math.log(1.5)
    process = build_evolutionary_process(
        tree,
        model="acdc",
        parameter=rate,
        root_mode="flat",
        allow_zero=True,
    )
    expected, expected_fit = compute_bm_marginals(
        tree, values, sigma2=0.8, _process=process
    )
    observed, observed_fit = compute_eb_marginals(
        tree, values, sigma2=0.8, eb_rate=rate
    )
    assert observed_fit.restricted_log_likelihood == pytest.approx(
        expected_fit.restricted_log_likelihood, abs=2e-12
    )
    for node in tree.traverse():
        assert observed[node] == pytest.approx(expected[node], abs=2e-12)


def test_process_supports_root_polytomy_and_zero_edges_before_adapter_policy():
    tree = tree_from("(A:0,B:1,C:2,D:0.5)R;")
    names = [str(leaf.name) for leaf in tree.leaves()]
    process = build_evolutionary_process(
        tree,
        model="ou",
        parameter=0.7,
        root_mode="stationary",
        allow_zero=True,
    )
    dense = process.tip_covariance(names)
    sparse_model = process.sparse_tip_model(names)
    np.testing.assert_allclose(
        sparse_model.materialize(),
        dense / np.mean(np.diag(dense)),
        rtol=2e-12,
        atol=2e-12,
    )


@pytest.mark.parametrize(
    "variances",
    [
        (np.nextafter(0.0, 1.0), 1.0),
        (np.finfo(float).max, np.finfo(float).max),
    ],
)
def test_sparse_process_handles_extreme_finite_variance_scales(variances):
    tree = tree_from("(A:1,B:1)R;")
    process = GaussianTreeProcess(
        tree=tree,
        transitions={
            tree["A"]: GaussianTransition(1.0, 0.0, variances[0]),
            tree["B"]: GaussianTransition(1.0, 0.0, variances[1]),
        },
        root=GaussianRootPrior("fixed", 0.0, 0.0),
        model="test",
    )
    dense = process.tip_covariance(["A", "B"])
    maximum = max(variances)
    normalization = maximum * np.mean(np.asarray(variances) / maximum)
    sparse_model = process.sparse_tip_model(["A", "B"])
    materialized = sparse_model.materialize()

    assert np.isfinite(sparse_model.precision.data).all()
    assert sparse_model.covariance_scale == normalization
    np.testing.assert_allclose(
        materialized,
        dense / normalization,
        rtol=5e-14,
        atol=0.0,
    )
    assert materialized[0, 0] > 0.0


def test_flat_root_process_refuses_finite_covariance_views():
    tree = tree_from("(A:1,B:2,C:3)R;")
    process = build_evolutionary_process(
        tree, model="brownian", root_mode="flat", allow_zero=True
    )
    with pytest.raises(ValueError, match="flat-root"):
        process.tip_covariance(["A", "B", "C"])
    with pytest.raises(ValueError, match="flat-root"):
        process.sparse_tip_model(["A", "B", "C"])


def test_covariance_factory_exposes_its_shared_process():
    tree = tree_from("((A:1,B:1):1,(C:1,D:1):1)R;")
    names = [str(leaf.name) for leaf in tree.leaves()]
    factory = evolutionary_covariance_factory(tree, names, model="ou")
    process = factory.process(0.6)
    np.testing.assert_allclose(process.tip_covariance(names), factory(0.6))
