import numpy as np
import pytest

from nwkit.general_ou_asr import compute_general_ou_marginals
from nwkit.util import read_tree


def tree_from(source):
    return read_tree(source, "1", True, quiet=True, rooted="yes")


@pytest.mark.parametrize(
    "root_prior,root_variance", [("fixed", None), ("gaussian", 0.7)]
)
def test_fixed_parameters_match_dense_conditioning(root_prior, root_variance):
    tree = tree_from("((A:0.3,B:0.7)I:0.4,C:1.1,D:0.2)R;")
    values = {"A": 1.2, "B": -0.4, "C": 2.1, "D": None}
    errors = {"A": 0.2, "B": 0.1, "C": 0.4}
    posterior, fit = compute_general_ou_marginals(
        tree,
        values,
        root_prior=root_prior,
        root_mean=-0.1,
        root_variance=root_variance,
        alpha=0.8,
        sigma2=1.7,
        theta=0.3,
        standard_errors=errors,
    )
    assert fit.log_likelihood is not None
    assert fit.likelihood_rank == 3
    assert fit.root_prior == root_prior
    assert all(
        np.isfinite(item.mean) and np.isfinite(item.variance)
        for item in posterior.values()
    )


def test_nonstationary_root_parameters_can_be_estimated():
    tree = tree_from("(((A:0.3,B:0.4)I:0.2,C:0.8)K:0.5,(D:0.4,E:0.7)J:0.3)R;")
    values = {"A": -0.2, "B": 0.1, "C": 0.7, "D": 1.1, "E": 1.3}
    _, fit = compute_general_ou_marginals(
        tree,
        values,
        root_prior="gaussian",
        root_mean=0.0,
        root_variance=1.0,
        alpha_bounds=(0.05, 5.0),
    )
    assert fit.alpha_estimated
    assert fit.sigma2_estimated
    assert fit.theta_estimated
    assert fit.log_likelihood is not None


def test_fixed_root_star_rejects_joint_alpha_sigma_confounding():
    tree = tree_from("(A:1,B:1,C:1,D:1,E:1,F:1)R;")
    with pytest.raises(ValueError, match="not jointly identifiable"):
        compute_general_ou_marginals(
            tree,
            {name: float(index) for index, name in enumerate("ABCDEF")},
            root_prior="fixed",
            root_mean=0.0,
            theta=0.0,
        )


def test_gaussian_root_requires_positive_variance():
    tree = tree_from("(A:1,B:1,C:1)R;")
    with pytest.raises(ValueError, match="root-variance"):
        compute_general_ou_marginals(
            tree,
            {"A": 0.0, "B": 1.0, "C": 2.0},
            root_prior="gaussian",
            root_mean=0.0,
            alpha=0.5,
            sigma2=1.0,
            theta=0.0,
        )
