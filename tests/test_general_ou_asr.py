import numpy as np
import pytest

from nwkit.continuous_asr_process import fitted_scalar_process
from nwkit.gaussian_inference import condition_gaussian_tree
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


def test_fixed_root_ou_has_finite_brownian_limit_without_stationary_variance():
    tree = tree_from("((A:0.5,B:1.1)I:0.4,C:1.2,D:1.7)R;")
    values = {"A": -1.0, "B": 0.4, "C": 0.7, "D": 1.6}
    posterior, fit = compute_general_ou_marginals(
        tree,
        values,
        root_prior="fixed",
        root_mean=0.0,
        alpha=1e-310,
        sigma2=0.8,
        theta=0.1,
    )
    process = fitted_scalar_process(tree, "OU", fit)
    reconstructed = condition_gaussian_tree(process, values)
    for node in tree.traverse():
        assert reconstructed.marginals[node] == posterior[node]
        if not node.is_root:
            assert process.transitions[node].variance == pytest.approx(0.8 * node.dist)
    assert fit.log_likelihood == reconstructed.log_likelihood


@pytest.mark.parametrize("scale", [1e-200, 1e200])
def test_nonstationary_ou_initialization_handles_extreme_time_units(scale):
    source = "(((A:0.3,B:0.4)I:0.2,C:0.8)K:0.5,(D:0.4,E:0.7)J:0.3)R;"
    values = {"A": -0.2, "B": 0.1, "C": 0.7, "D": 1.1, "E": 1.3}
    fits = []
    for time_scale in (1.0, scale):
        tree = tree_from(source)
        for node in tree.traverse():
            if not node.is_root:
                node.dist *= time_scale
        _, fit = compute_general_ou_marginals(
            tree,
            values,
            root_prior="gaussian",
            root_mean=0.0,
            root_variance=1.0,
            alpha_bounds=(0.05 / time_scale, 5.0 / time_scale),
        )
        fits.append(fit)
    reference, scaled = fits
    assert scaled.alpha * scale == pytest.approx(reference.alpha, rel=1e-5)
    assert scaled.sigma2 * scale == pytest.approx(reference.sigma2, rel=1e-5)
    assert scaled.log_likelihood == pytest.approx(reference.log_likelihood, abs=1e-9)


@pytest.mark.parametrize("scale", [1e-95, 1e95])
def test_nonstationary_ou_theta_and_diffusion_are_trait_unit_invariant(scale):
    tree = tree_from("((A:0.5,B:1.1)I:0.4,C:1.2,D:1.7)R;")
    values = {"A": -1.0, "B": 0.4, "C": 0.7, "D": 1.6}
    options = {"root_prior": "fixed", "root_mean": 0.0, "alpha": 0.5}
    expected, fit = compute_general_ou_marginals(tree, values, **options)
    posterior, scaled = compute_general_ou_marginals(
        tree, {name: value * scale for name, value in values.items()}, **options
    )
    assert scaled.theta / scale == pytest.approx(fit.theta, rel=1e-5)
    assert scaled.sigma2 / scale**2 == pytest.approx(fit.sigma2, rel=1e-5)
    assert scaled.log_likelihood + len(values) * np.log(scale) == pytest.approx(
        fit.log_likelihood, abs=1e-9
    )
    for node in tree.traverse():
        assert posterior[node].mean / scale == pytest.approx(
            expected[node].mean, abs=1e-5
        )
        assert posterior[node].variance / scale**2 == pytest.approx(
            expected[node].variance, abs=1e-5
        )


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
