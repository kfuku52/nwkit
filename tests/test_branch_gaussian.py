import math

import numpy as np
import pytest

from nwkit.branch_gaussian import (
    BranchGaussianModel,
    BrownianBranch,
    GaussianJump,
    OUBranch,
    build_branch_gaussian_process,
)
from nwkit.evolution import build_evolutionary_process
from nwkit.gaussian_inference import (
    condition_gaussian_tree,
    gaussian_tree_likelihood,
    simulate_gaussian_process,
)
from nwkit.gaussian_tree import GaussianRootPrior
from nwkit.util import assign_branch_ids, read_tree


def tree_from(source="((A:0.3,B:0.7)I:0.4,C:1.1)R;"):
    return read_tree(source, "1", True, quiet=True, rooted="yes")


def assign(tree, model):
    return {
        identifier: model
        for node, identifier in assign_branch_ids(tree).items()
        if not node.is_root
    }


@pytest.mark.parametrize("alpha", [None, 0.0, 0.7])
def test_homogeneous_reduction(alpha):
    tree = tree_from()
    diffusion = BrownianBranch(1.4) if alpha is None else OUBranch(alpha, 1.4)
    actual = build_branch_gaussian_process(
        tree,
        assign(tree, BranchGaussianModel(diffusion)),
        root=GaussianRootPrior("gaussian", 0.3, 0.9),
    )
    expected = build_evolutionary_process(
        tree,
        model="ou" if alpha else "brownian",
        parameter=alpha if alpha else None,
        variance_scale=1.4,
        root_mode="gaussian",
        root_mean=0.3,
        root_variance=0.9 / 1.4,  # Existing builder scales the root too.
        allow_zero=True,
    )
    for node in actual.transitions:
        assert actual.transitions[node].slope == pytest.approx(
            expected.transitions[node].slope
        )
        assert actual.transitions[node].variance == pytest.approx(
            expected.transitions[node].variance
        )
    np.testing.assert_allclose(
        actual.tip_covariance(["A", "B", "C"]), expected.tip_covariance(["A", "B", "C"])
    )


def mixed_process():
    tree = tree_from()
    models = {
        "I": BranchGaussianModel(BrownianBranch(0.8)),
        "A": BranchGaussianModel(OUBranch(0.6, 1.3, 0.9), GaussianJump(0.2, 0.15)),
        "B": BranchGaussianModel(jump=GaussianJump(-0.4, 0.5)),
        "C": BranchGaussianModel(BrownianBranch(0.7)),
    }
    process = build_branch_gaussian_process(
        tree,
        {
            identifier: models[node.name]
            for node, identifier in assign_branch_ids(tree).items()
            if not node.is_root
        },
        root=GaussianRootPrior("gaussian", 0.3, 0.9),
    )
    # Independent structural-equation oracle in R,I,C,A,B order, assembled
    # directly from biological parameters rather than process transitions.
    a = math.exp(-0.6 * 0.3)
    mean = np.array([0.3, 0.3, 0.3, a * 0.3 + (1 - a) * 0.9 + 0.2, -0.1])
    loading = np.array(
        [
            [1, 0, 0, 0, 0],
            [1, 1, 0, 0, 0],
            [1, 0, 1, 0, 0],
            [a, a, 0, 1, 0],
            [1, 1, 0, 0, 1],
        ]
    )
    innovation = np.array(
        [0.9, 0.8 * 0.4, 0.7 * 1.1, 1.3 / (2 * 0.6) * (1 - a * a) + 0.15, 0.5]
    )
    covariance = (loading * innovation) @ loading.T
    nodes = {node.name: node for node in tree.traverse()}
    return (
        process,
        [nodes[name] for name in ("R", "I", "C", "A", "B")],
        mean,
        covariance,
    )


def test_mixed_likelihood_smoothing_and_covariance_against_independent_oracle():
    process, nodes, mean, covariance = mixed_process()
    np.testing.assert_allclose(process.covariance(nodes), covariance, atol=1e-14)
    np.testing.assert_allclose(
        [process.marginal_moments()[0][node] for node in nodes], mean
    )
    observed = [3, 4, 2]
    values = {"A": 1.2, "B": -0.4, "C": 2.1}
    errors = {"A": 0.2, "B": 0.1, "C": 0.4}
    obs_cov = covariance[np.ix_(observed, observed)] + np.diag([0.04, 0.01, 0.16])
    residual = np.array(list(values.values())) - mean[observed]
    solve = np.linalg.solve(obs_cov, residual)
    expected_mean = mean + covariance[:, observed] @ solve
    expected_variance = np.diag(
        covariance
        - covariance[:, observed] @ np.linalg.solve(obs_cov, covariance[observed, :])
    )
    expected_ll = -0.5 * (
        3 * math.log(2 * math.pi) + np.linalg.slogdet(obs_cov)[1] + residual @ solve
    )
    result = condition_gaussian_tree(process, values, standard_errors=errors)
    assert result.log_likelihood == pytest.approx(expected_ll, abs=1e-12)
    assert gaussian_tree_likelihood(
        process, values, standard_errors=errors
    ).log_likelihood == pytest.approx(expected_ll, abs=1e-12)
    np.testing.assert_allclose(
        [result.marginals[node].mean for node in nodes], expected_mean, atol=1e-12
    )
    np.testing.assert_allclose(
        [result.marginals[node].variance for node in nodes],
        expected_variance,
        atol=1e-12,
    )


def test_mixed_simulation_moments_and_seed():
    process, nodes, mean, covariance = mixed_process()
    samples = simulate_gaussian_process(process, num_samples=80000, seed=13)
    columns = [samples.nodes.index(node) for node in nodes]
    data = samples.values[:, columns]
    np.testing.assert_allclose(data.mean(axis=0), mean, atol=0.018)
    np.testing.assert_allclose(np.cov(data, rowvar=False), covariance, atol=0.025)
    first = simulate_gaussian_process(process, num_samples=4, seed=42)
    second = simulate_gaussian_process(process, num_samples=4, seed=42)
    np.testing.assert_array_equal(first.values, second.values)


def test_jump_is_at_end_and_not_scaled_by_time():
    model = BranchGaussianModel(OUBranch(2, 0, 3), GaussianJump(4, 5))
    transition = model.transition(0.7)
    assert transition.intercept == pytest.approx((1 - math.exp(-1.4)) * 3 + 4)
    assert transition.variance == 5
    assert model.transition(0).intercept == 4
    assert model.transition(0).variance == 5
    assert BranchGaussianModel(jump=GaussianJump(4, 5)).transition(
        100
    ) == BranchGaussianModel(jump=GaussianJump(4, 5)).transition(0)
    assert BranchGaussianModel(BrownianBranch(2), GaussianJump()).transition(
        3
    ) == BranchGaussianModel(BrownianBranch(2)).transition(3)


def test_ou_limits_and_extreme_finite_time_variance():
    transition = BranchGaussianModel(OUBranch(1e-310, 2, 4)).transition(0.5)
    assert transition.variance == pytest.approx(1)
    assert transition.slope == 1
    assert BranchGaussianModel(OUBranch(1e308, 1e308)).transition(
        10
    ).variance == pytest.approx(0.5)
    assert BranchGaussianModel(OUBranch(1e-310, 1e-310)).transition(
        1e308
    ).variance == pytest.approx(-math.expm1(-0.02) / 2)
    with pytest.raises(ValueError, match="underflows"):
        BranchGaussianModel(BrownianBranch(1e-300)).transition(1e-300)
    with pytest.raises(ValueError, match="overflows"):
        BranchGaussianModel(BrownianBranch(1e300)).transition(1e300)


@pytest.mark.parametrize(
    "factory",
    [
        lambda: BrownianBranch(-1),
        lambda: BrownianBranch(True),
        lambda: OUBranch(-1),
        lambda: OUBranch(1, optimum=float("nan")),
        lambda: GaussianJump(variance=-1),
        lambda: GaussianJump(mean=float("inf")),
        lambda: BranchGaussianModel(),
        lambda: BranchGaussianModel(diffusion="bm"),
        lambda: BranchGaussianModel(jump=1),
    ],
)
def test_invalid_models(factory):
    with pytest.raises(ValueError):
        factory()


@pytest.mark.parametrize("length", [-1, float("nan"), float("inf"), True, "bad"])
def test_invalid_lengths(length):
    with pytest.raises(ValueError):
        BranchGaussianModel(BrownianBranch()).transition(length)


@pytest.mark.parametrize("key", [0, 99, "1", 1.0, True])
def test_branch_id_contract(key):
    tree = tree_from("(A:1)R;")
    with pytest.raises(ValueError):
        build_branch_gaussian_process(
            tree,
            {key: BranchGaussianModel(BrownianBranch())},
            root=GaussianRootPrior("fixed"),
        )


def test_missing_branch_invalid_model_and_root():
    tree = tree_from("(A:1)R;")
    for mapping, root in [
        ({}, GaussianRootPrior("fixed")),
        ({1: 1}, GaussianRootPrior("fixed")),
        ({1: BranchGaussianModel(BrownianBranch())}, None),
    ]:
        with pytest.raises(ValueError):
            build_branch_gaussian_process(tree, mapping, root=root)


def test_flat_root_and_missing_observation():
    tree = tree_from()
    process = build_branch_gaussian_process(
        tree,
        assign(tree, BranchGaussianModel(BrownianBranch())),
        root=GaussianRootPrior("flat", variance=None),
    )
    result = condition_gaussian_tree(process, {"A": 1.0, "B": None, "C": 2.0})
    assert result.num_observed == 2
    assert math.isfinite(result.log_likelihood)
    with pytest.raises(ValueError):
        simulate_gaussian_process(process, num_samples=2, seed=1)
    samples = simulate_gaussian_process(process, num_samples=2, seed=1, root_values=0.3)
    assert samples.values.shape == (2, 5)


def test_ou_extreme_rate_ratio_avoids_premature_rounding():
    rate = float.fromhex("0x0.0000000000003p-1022")
    transition = BranchGaussianModel(OUBranch(1e-307, rate)).transition(1e308)
    expected = (rate / 1e-307) * (-math.expm1(-20.0) / 2.0)
    assert transition.variance == pytest.approx(expected, rel=1e-13, abs=0)
    transition = BranchGaussianModel(OUBranch(0.5, 1e308)).transition(10)
    assert transition.variance == pytest.approx(1e308 * -math.expm1(-10), rel=1e-12)


@pytest.mark.parametrize(
    "factory",
    [
        lambda: BrownianBranch(np.bool_(True)),
        lambda: OUBranch(np.bool_(False)),
        lambda: GaussianJump(np.bool_(True)),
        lambda: BranchGaussianModel(BrownianBranch()).transition(np.bool_(False)),
    ],
)
def test_numpy_booleans_are_not_parameters(factory):
    with pytest.raises(ValueError, match="boolean"):
        factory()


def test_mapping_contract_and_attached_subtree():
    tree = tree_from()
    model = BranchGaussianModel(BrownianBranch())
    with pytest.raises(ValueError, match="mapping"):
        build_branch_gaussian_process(
            tree, [(1, model), (1, model)], root=GaussianRootPrior("fixed")
        )
    subtree = next(node for node in tree.traverse() if node.name == "I")
    with pytest.raises(ValueError, match="root node"):
        build_branch_gaussian_process(subtree, {}, root=GaussianRootPrior("fixed"))


def test_extreme_ou_coefficients_against_decimal_oracle():
    from decimal import Decimal, localcontext

    rng = np.random.default_rng(918)
    cases = [
        (1e-300, 1e-300, 1e300, 1e300),
        (1e-200, 1e-120, 1e300, -1e200),
        (1e308, 1e308, 1e308, -2.0),
    ]
    cases += [
        tuple(10.0 ** float(exponent) for exponent in rng.uniform(-310, 308, 4))
        for _ in range(150)
    ]
    with localcontext() as context:
        context.prec = 750
        for alpha, length, rate, optimum in cases:
            a, t, r, theta = map(Decimal.from_float, (alpha, length, rate, optimum))
            exponent = a * t
            decay = (-exponent).exp() if exponent < 1000 else Decimal(0)
            expected_variance = float(r * (1 - decay * decay) / (2 * a))
            expected_intercept = float((1 - decay) * theta)
            model = BranchGaussianModel(OUBranch(alpha, rate, optimum))
            if expected_variance == 0 or not math.isfinite(expected_variance):
                with pytest.raises(ValueError, match="floating-point range"):
                    model.transition(length)
            else:
                result = model.transition(length)
                assert result.variance == pytest.approx(
                    expected_variance, rel=2e-14, abs=5e-324
                )
                assert result.intercept == pytest.approx(
                    expected_intercept, rel=2e-14, abs=5e-324
                )
                assert result.slope == pytest.approx(
                    float(decay), rel=2e-14, abs=5e-324
                )


def test_zero_length_end_jump_is_not_contracted_by_inference():
    tree = tree_from("(A:0,B:0)R;")
    process = build_branch_gaussian_process(
        tree,
        assign(tree, BranchGaussianModel(jump=GaussianJump(0.4, 0.5))),
        root=GaussianRootPrior("fixed", 0.3),
    )
    result = condition_gaussian_tree(process, {"A": 1.0, "B": -0.1})
    expected_ll = -math.log(2 * math.pi * 0.5) - (
        (1.0 - 0.7) ** 2 + (-0.1 - 0.7) ** 2
    ) / (2 * 0.5)
    assert result.log_likelihood == pytest.approx(expected_ll)
    for node in tree.leaves():
        assert result.marginals[node].variance == 0
    assert result.marginals[tree].mean == pytest.approx(0.3)


def test_deterministic_mixed_process_and_posterior_sampling():
    from nwkit.gaussian_inference import sample_gaussian_posterior

    tree = tree_from("(A:1,B:1)R;")
    models = {
        identifier: BranchGaussianModel(OUBranch(0.5, 0, 2), GaussianJump(0.4, 0))
        for node, identifier in assign_branch_ids(tree).items()
        if not node.is_root
    }
    process = build_branch_gaussian_process(
        tree, models, root=GaussianRootPrior("fixed", 2)
    )
    result = condition_gaussian_tree(process, {"A": 2.4, "B": 2.4})
    assert all(marginal.variance == 0 for marginal in result.marginals.values())
    with pytest.raises(ValueError, match="likelihood|Conflicting"):
        condition_gaussian_tree(process, {"A": 2.4, "B": 3.0})
    samples = sample_gaussian_posterior(
        process, {"A": 2.4, "B": 2.4}, num_samples=3, seed=42
    )
    np.testing.assert_allclose(samples.values, np.array([[2, 2.4, 2.4]] * 3))
