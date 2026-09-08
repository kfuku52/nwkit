import numpy as np
import pytest

from nwkit.phylogenetic_predictive import sister_clade_contrast
from nwkit.util import read_tree
from nwkit.vector_processes import vector_brownian_process, vector_ou_process
from nwkit.vector_simulation import (
    simulate_vector_process,
    simulated_vector_observations,
)


def test_vector_simulation_preserves_joint_covariance():
    tree = read_tree("(A:1,B:2)R;", "1", True, quiet=True, rooted="yes")
    process = vector_ou_process(
        tree, np.eye(2), np.array([[2, 0.6], [0.6, 1]]), np.zeros(2)
    )
    nodes, values = simulate_vector_process(process, 30000, seed=3)
    assert np.cov(values[:, 0].T) == pytest.approx(process.root_covariance, abs=0.025)
    assert np.cov(values[:, 0, 0], values[:, 1, 1])[0, 1] == pytest.approx(
        0.3 * np.exp(-1), abs=0.02
    )
    assert len(nodes) == 3


def test_simulation_preserves_missing_mask_and_correlated_error():
    tree = read_tree("(A:0,B:1)R;", "1", True, quiet=True, rooted="yes")
    process = vector_brownian_process(tree, np.eye(2))
    observed = {"A": [1, 2], "B": [None, 3]}
    errors = {"A": np.array([[1, 0.7], [0.7, 2]]), "B": np.zeros((2, 2))}
    values = simulated_vector_observations(
        process, observed, errors, 30000, 3, root_values=np.zeros(2)
    )
    assert all(value["B"][0] is None for value in values)
    assert np.cov(np.array([value["A"] for value in values]).T) == pytest.approx(
        errors["A"], abs=0.035
    )
    with pytest.raises(ValueError, match="flat root"):
        simulate_vector_process(process)


def test_sister_clade_contrast_depends_on_topology_and_mask():
    tree = read_tree("((A:1,B:1):1,(C:1,D:1):1)R;", "1", True, quiet=True, rooted="yes")
    assert sister_clade_contrast(
        tree, {"A": 0, "B": 0, "C": 2, "D": 2}
    ) == pytest.approx(4 / 3)
    assert sister_clade_contrast(
        tree, {"A": 0, "B": 2, "C": 0, "D": 2}
    ) == pytest.approx(8 / 3)
    assert sister_clade_contrast(tree, {"A": 0, "C": 2}) == 4
