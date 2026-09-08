import numpy as np
import pytest
from scipy.stats import multivariate_normal

from nwkit.compiled_tree import CompiledTree
from nwkit.util import read_tree
from nwkit.vector_gaussian import VectorProcess, VectorTransition, condition_vector_tree


def fixture():
    tree = read_tree("((A:1,B:0.5):0.7,C:1.5)R;", "1", True, quiet=True, rooted="yes")
    nodes = CompiledTree.from_tree(tree)
    transitions = {
        node: VectorTransition(
            np.array([[0.8, 0.1], [-0.1, 0.9]]),
            np.array([0.2, -0.3]),
            np.array([[0.7, 0.2], [0.2, 0.5]]) * node.dist,
        )
        for node in nodes.nodes[1:]
    }
    process = VectorProcess(
        tree, 2, transitions, np.array([0.4, -0.2]), np.array([[1.0, 0.3], [0.3, 0.8]])
    )
    return process, nodes


def dense_joint(process, compiled):
    count = len(compiled.nodes)
    mean = np.zeros(count * 2)
    covariance = np.zeros((count * 2, count * 2))
    mean[:2] = process.root_mean
    covariance[:2, :2] = process.root_covariance
    for i in range(1, count):
        p = compiled.parents[i]
        transition = process.transitions[compiled.nodes[i]]
        index = slice(2 * i, 2 * i + 2)
        parent = slice(2 * p, 2 * p + 2)
        mean[index] = transition.slope @ mean[parent] + transition.intercept
        covariance[index, : 2 * i] = transition.slope @ covariance[parent, : 2 * i]
        covariance[: 2 * i, index] = covariance[index, : 2 * i].T
        covariance[index, index] = (
            transition.slope @ covariance[parent, parent] @ transition.slope.T
            + transition.covariance
        )
    return mean, covariance


@pytest.mark.parametrize("noisy", [False, True])
def test_vector_pruning_equals_dense_joint_conditioning(noisy):
    process, compiled = fixture()
    observed = {"A": [1.2, None], "B": [-0.3, 0.5], "C": [None, -0.8]}
    errors = (
        {name: np.array([[0.1, 0.03], [0.03, 0.2]]) for name in observed}
        if noisy
        else None
    )
    result = condition_vector_tree(process, observed, error_covariances=errors)
    mean, covariance = dense_joint(process, compiled)
    indices = [
        2 * compiled.leaf_index_by_name[name] + k
        for name, values in observed.items()
        for k, value in enumerate(values)
        if value is not None
    ]
    values = np.array(
        [value for vector in observed.values() for value in vector if value is not None]
    )
    noise = np.zeros((4, 4))
    if noisy:
        noise[np.diag_indices(4)] = [0.1, 0.1, 0.2, 0.2]
        noise[1, 2] = noise[2, 1] = 0.03
    observed_covariance = covariance[np.ix_(indices, indices)] + noise
    inverse = np.linalg.inv(observed_covariance)
    cross = covariance[:, indices]
    expected_mean = mean + cross @ inverse @ (values - mean[indices])
    expected_covariance = covariance - cross @ inverse @ cross.T
    assert result.means.reshape(-1) == pytest.approx(expected_mean, abs=1e-12)
    for i in range(len(compiled.nodes)):
        assert result.covariances[i] == pytest.approx(
            expected_covariance[2 * i : 2 * i + 2, 2 * i : 2 * i + 2], abs=1e-12
        )
    assert result.log_likelihood == pytest.approx(
        multivariate_normal.logpdf(values, mean[indices], observed_covariance),
        abs=1e-12,
    )
    draws = result.sample(20000, seed=12).reshape(20000, -1)
    assert draws.mean(axis=0) == pytest.approx(expected_mean, abs=0.025)
    assert np.cov(draws, rowvar=False) == pytest.approx(expected_covariance, abs=0.025)


def test_vector_flat_root_matches_scalar_independent_traits():
    from nwkit.evolution import build_evolutionary_process
    from nwkit.gaussian_inference import condition_gaussian_tree

    process, compiled = fixture()
    process = VectorProcess(
        process.tree,
        2,
        {
            node: VectorTransition(np.eye(2), np.zeros(2), np.eye(2) * node.dist)
            for node in compiled.nodes[1:]
        },
    )
    observed = {"A": [1.2, None], "B": [-0.3, 0.5], "C": [None, -0.8]}
    result = condition_vector_tree(process, observed)
    likelihood = 0
    for trait in range(2):
        scalar = condition_gaussian_tree(
            build_evolutionary_process(
                process.tree, model="brownian", root_mode="flat"
            ),
            {name: values[trait] for name, values in observed.items()},
        )
        likelihood += scalar.log_likelihood
        for i, node in enumerate(compiled.nodes):
            assert result.means[i, trait] == pytest.approx(scalar.marginals[node].mean)
            assert result.covariances[i, trait, trait] == pytest.approx(
                scalar.marginals[node].variance
            )
    assert result.log_likelihood == pytest.approx(likelihood)


def test_zero_length_branches_preserve_exact_observations():
    from nwkit.vector_processes import vector_brownian_process

    tree = read_tree("((A:0,B:0):1,C:1)R;", "1", True, quiet=True, rooted="yes")
    process = vector_brownian_process(tree, np.eye(2))
    result = condition_vector_tree(
        process, {"A": [1.0, None], "B": [1.0, 2.0], "C": [2.0, 3.0]}
    )
    for name in ("A", "B"):
        index = next(i for i, node in enumerate(result.nodes) if node.name == name)
        assert result.means[index] == pytest.approx([1.0, 2.0])
        assert result.covariances[index] == pytest.approx(np.zeros((2, 2)))
    with pytest.raises(ValueError, match="Conflicting exact"):
        condition_vector_tree(
            process, {"A": [1.0, None], "B": [2.0, None], "C": [2.0, 3.0]}
        )


def test_general_ou_stationarity_and_short_edges():
    from nwkit.vector_processes import vector_ou_process

    tree = read_tree(
        "(A:0.000000000001,B:1,C:1000)R;", "1", True, quiet=True, rooted="yes"
    )
    attraction = np.array([[1.0, 0.3], [-0.2, 0.5]])
    diffusion = np.array([[0.7, 0.1], [0.1, 0.8]])
    process = vector_ou_process(tree, attraction, diffusion, [2.0, -1.0])
    for node, branch in process.transitions.items():
        assert (
            branch.slope @ process.root_covariance @ branch.slope.T + branch.covariance
            == pytest.approx(process.root_covariance, abs=1e-12)
        )
        if node.name == "A":
            assert branch.covariance / node.dist == pytest.approx(diffusion, abs=1e-11)
    with pytest.raises(ValueError, match="positive real"):
        vector_ou_process(tree, -attraction, diffusion, [2.0, -1.0])


def test_nearly_singular_proper_root_does_not_create_spurious_likelihood():
    from scipy.stats import multivariate_normal

    from nwkit.vector_processes import vector_ou_process

    tree = read_tree("(A:1)R;", "1", True, quiet=True, rooted="yes")
    sigma = np.array([[1e-12, 1e-10], [1e-10, 1.000001e-8]])
    mean = np.array([2.0, -3.0])
    y = np.array([1.0, 4.0])
    error = np.diag([0.1, 0.2])
    process = vector_ou_process(tree, np.eye(2) * 0.4, 0.8 * sigma, mean)
    result = condition_vector_tree(process, {"A": y}, error_covariances={"A": error})
    assert result.log_likelihood == pytest.approx(
        multivariate_normal.logpdf(y, mean=mean, cov=sigma + error), abs=1e-8
    )
