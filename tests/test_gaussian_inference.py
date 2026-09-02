import math

import numpy as np
import pytest

from nwkit.compiled_tree import CompiledTree
from nwkit.evolution import build_evolutionary_process
from nwkit.gaussian_inference import (
    condition_gaussian_tree,
    gaussian_tree_likelihood,
    sample_gaussian_posterior,
    simulate_gaussian_process,
)
from nwkit.gaussian_tree import (
    GaussianRootPrior,
    GaussianTransition,
    GaussianTreeProcess,
)
from nwkit.util import read_tree


def tree_from(source):
    return read_tree(source, "1", True, quiet=True, rooted="yes")


def dense_condition(process, values, errors):
    nodes = list(process.tree.traverse(strategy="preorder"))
    index = {node: position for position, node in enumerate(nodes)}
    means, covariance = process.marginal_moments()[0], process.covariance(nodes)
    observed = [
        node for node in process.tree.leaves() if values.get(str(node.name)) is not None
    ]
    observed_indices = [index[node] for node in observed]
    response = np.asarray([values[str(node.name)] for node in observed], dtype=float)
    observation_covariance = covariance[np.ix_(observed_indices, observed_indices)]
    observation_covariance += np.diag(
        [errors.get(str(node.name), 0.0) ** 2 for node in observed]
    )
    mean_vector = np.asarray([means[node] for node in nodes])
    residual = response - mean_vector[observed_indices]
    solved = np.linalg.solve(observation_covariance, residual)
    cross = covariance[:, observed_indices]
    conditional_mean = mean_vector + cross @ solved
    conditional_covariance = covariance - cross @ np.linalg.solve(
        observation_covariance, cross.T
    )
    likelihood = -0.5 * (
        len(response) * math.log(2.0 * math.pi)
        + np.linalg.slogdet(observation_covariance)[1]
        + residual @ solved
    )
    return nodes, conditional_mean, np.diag(conditional_covariance), likelihood


@pytest.mark.parametrize(
    "model,parameter", [("brownian", None), ("ou", 0.8), ("acdc", 0.3)]
)
def test_proper_process_matches_dense_gaussian_conditioning(model, parameter):
    tree = tree_from("((A:0.3,B:0.7)I:0.4,C:1.1,D:0.2)R;")
    process = build_evolutionary_process(
        tree,
        model=model,
        parameter=parameter,
        root_mode="gaussian",
        root_mean=0.25,
        root_variance=0.9,
        variance_scale=1.4,
        allow_zero=True,
    )
    values = {"A": 1.2, "B": -0.4, "C": 2.1, "D": None}
    errors = {"A": 0.2, "B": 0.1, "C": 0.4}
    result = condition_gaussian_tree(process, values, standard_errors=errors)
    nodes, means, variances, likelihood = dense_condition(process, values, errors)
    assert result.log_likelihood == pytest.approx(likelihood, abs=2e-12)
    assert result.likelihood_rank == 3
    for position, node in enumerate(nodes):
        assert result.marginals[node].mean == pytest.approx(means[position], abs=2e-12)
        assert result.marginals[node].variance == pytest.approx(
            max(0.0, variances[position]), abs=2e-12
        )


def test_stationary_ou_matches_existing_dense_oracle_contract():
    tree = tree_from("((A:0.3,B:0.7)I:0.4,C:1.1,D:0.2)R;")
    process = build_evolutionary_process(
        tree,
        model="ou",
        parameter=0.8,
        root_mode="stationary",
        root_mean=0.3,
        variance_scale=1.7,
        allow_zero=True,
    )
    values = {"A": 1.2, "B": -0.4, "C": 2.1, "D": None}
    errors = {"A": 0.2, "B": 0.0, "C": 0.4}
    result = condition_gaussian_tree(process, values, standard_errors=errors)
    _, means, variances, likelihood = dense_condition(process, values, errors)
    assert result.log_likelihood == pytest.approx(likelihood, abs=2e-12)
    for position, node in enumerate(process.tree.traverse(strategy="preorder")):
        assert result.marginals[node].mean == pytest.approx(means[position], abs=2e-12)
        assert result.marginals[node].variance == pytest.approx(
            max(0.0, variances[position]), abs=2e-12
        )


def test_flat_brownian_matches_existing_pruning_with_missing_and_errors():
    from nwkit.continuous_asr import compute_bm_marginals

    tree = tree_from("((A:0.3,B:0.7)I:0.4,C:1.1,D:0.2)R;")
    process = build_evolutionary_process(
        tree, model="brownian", root_mode="flat", allow_zero=True
    ).scaled_variance(1.7)
    values = {"A": 1.2, "B": -0.4, "C": 2.1, "D": None}
    errors = {"A": 0.2, "B": 0.0, "C": 0.4}
    result = condition_gaussian_tree(process, values, standard_errors=errors)
    expected, fit = compute_bm_marginals(
        tree, values, sigma2=1.7, standard_errors=errors
    )
    assert result.log_likelihood == pytest.approx(
        fit.restricted_log_likelihood, abs=2e-12
    )
    assert result.likelihood_rank == fit.residual_df
    for node in tree.traverse():
        assert result.marginals[node].mean == pytest.approx(
            expected[node].mean, abs=2e-12
        )
        assert result.marginals[node].variance == pytest.approx(
            expected[node].variance, abs=2e-12
        )


def test_generic_affine_process_and_compiled_tree_reuse():
    tree = tree_from("((A:1,B:1)I:1,C:1)R;")
    transitions = {
        node: GaussianTransition(0.5, 1.0, 0.3)
        for node in tree.traverse()
        if not node.is_root
    }
    process = GaussianTreeProcess(
        tree,
        transitions,
        GaussianRootPrior("gaussian", -0.2, 0.7),
        "test-affine",
    )
    compiled = CompiledTree.from_tree(tree)
    first = condition_gaussian_tree(
        process, {"A": 2.0, "B": None, "C": -1.0}, compiled_tree=compiled
    )
    second = condition_gaussian_tree(
        process, {"A": 2.0, "B": None, "C": -1.0}, compiled_tree=compiled
    )
    likelihood = gaussian_tree_likelihood(
        process, {"A": 2.0, "B": None, "C": -1.0}, compiled_tree=compiled
    )
    assert first.marginals == second.marginals
    assert first.log_likelihood == second.log_likelihood
    assert likelihood.log_likelihood == first.log_likelihood
    assert likelihood.likelihood_rank == first.likelihood_rank
    assert likelihood.num_observed == first.num_observed


def test_zero_variance_identity_branch_enforces_exact_state():
    tree = tree_from("((A:0,B:1)I:0,C:1)R;")
    process = build_evolutionary_process(
        tree, model="brownian", root_mode="flat", allow_zero=True
    )
    result = condition_gaussian_tree(process, {"A": 3.0, "B": None, "C": 1.0})
    assert result.marginals[tree].mean == result.marginals[tree["I"]].mean
    assert result.marginals[tree["I"]].mean == result.marginals[tree["A"]].mean
    assert result.marginals[tree["A"]].variance == 0.0


def test_observed_position_count_contracts_deterministic_identity_edges():
    tree = tree_from("((A:0,B:0)I:0,C:1)R;")
    process = build_evolutionary_process(
        tree,
        model="brownian",
        root_mode="gaussian",
        root_variance=1.0,
        allow_zero=True,
    )
    result = condition_gaussian_tree(process, {"A": 1.0, "B": 1.0, "C": 2.0})
    assert result.num_observed == 3
    assert result.num_observed_positions == 2
    assert result.likelihood_rank == 2


def test_scaled_extreme_offset_preserves_likelihood_and_marginals():
    tree = tree_from("((A:1,B:1)I:1,C:2)R;")
    process = build_evolutionary_process(
        tree,
        model="brownian",
        root_mode="gaussian",
        root_mean=1e150,
        root_variance=4.0,
        variance_scale=2e280,
    )
    values = {"A": 1.1e150, "B": 0.9e150, "C": 1.2e150}
    result = condition_gaussian_tree(process, values)
    assert math.isfinite(result.log_likelihood)
    assert all(
        math.isfinite(marginal.mean) and math.isfinite(marginal.variance)
        for marginal in result.marginals.values()
    )


def test_joint_posterior_samples_match_marginal_moments_and_exact_tip():
    tree = tree_from("((A:0.5,B:0.8)I:0.4,C:1.1)R;")
    process = build_evolutionary_process(
        tree, model="brownian", root_mode="flat", allow_zero=True
    ).scaled_variance(0.7)
    values = {"A": 1.0, "B": -0.5, "C": 0.3}
    errors = {"A": 0.0, "B": 0.2, "C": 0.1}
    conditioned = condition_gaussian_tree(process, values, standard_errors=errors)
    samples = sample_gaussian_posterior(
        process,
        values,
        standard_errors=errors,
        num_samples=20000,
        seed=82,
    )
    index = {node: position for position, node in enumerate(samples.nodes)}
    assert np.all(samples.values[:, index[tree["A"]]] == 1.0)
    for node in (tree, tree["I"], tree["B"]):
        values_at_node = samples.values[:, index[node]]
        assert np.mean(values_at_node) == pytest.approx(
            conditioned.marginals[node].mean, abs=0.02
        )
        assert np.var(values_at_node) == pytest.approx(
            conditioned.marginals[node].variance, abs=0.02
        )


def test_prior_process_simulation_matches_affine_moments():
    tree = tree_from("((A:0.5,B:0.8)I:0.4,C:1.1)R;")
    process = build_evolutionary_process(
        tree,
        model="ou",
        parameter=0.8,
        root_mode="stationary",
        root_mean=0.3,
        variance_scale=1.2,
    )
    samples = simulate_gaussian_process(process, num_samples=30000, seed=91)
    expected_means, expected_variances = process.marginal_moments()
    for index, node in enumerate(samples.nodes):
        assert np.mean(samples.values[:, index]) == pytest.approx(
            expected_means[node], abs=0.02
        )
        assert np.var(samples.values[:, index]) == pytest.approx(
            expected_variances[node], abs=0.03
        )
