import math
from types import SimpleNamespace

import numpy as np
import pytest

from nwkit.asr import _multivariate_fit_status_message
from nwkit.multivariate_gaussian_asr import (
    _geometry,
    _prepare_observations,
    fit_dense_mvbm,
    fit_dense_mvou,
    fit_dense_mvou_diag,
    parse_alpha_by_trait,
)
from nwkit.util import read_tree


def tree_from(source):
    return read_tree(source, "1", True, quiet=True, rooted="yes")


def test_dense_mvbm_accepts_partial_vectors_and_errors():
    tree = tree_from("((A:0.5,B:0.7)I:0.4,(C:0.6,D:0.8)J:0.3,E:1.0)R;")
    values = {
        "A": (0.0, 1.0),
        "B": (0.4, None),
        "C": (1.2, 0.2),
        "D": (None, -0.3),
        "E": (1.7, 0.8),
    }
    errors = {
        "A": (0.1, 0.2),
        "B": (0.1, None),
        "C": (0.2, 0.1),
        "D": (None, 0.2),
        "E": (0.1, 0.1),
    }
    posterior, fit = fit_dense_mvbm(tree, values, ("x", "y"), standard_errors=errors)
    assert fit.restricted_log_likelihood is not None
    assert fit.num_effective_observations == 5
    assert fit.sigma.shape == (2, 2)
    assert np.all(np.linalg.eigvalsh(fit.sigma) > 0.0)
    assert posterior[tree["B"]].covariance[1, 1] > 0.0
    assert posterior[tree["D"]].covariance[0, 0] > 0.0


def test_dense_mvbm_exact_complete_case_agrees_with_fast_path_reasonably():
    from nwkit.multivariate_asr import compute_mvbm_marginals

    tree = tree_from("((A:0.5,B:0.7)I:0.4,(C:0.6,D:0.8)J:0.3,E:1.0)R;")
    values = {
        "A": (0.0, 1.0),
        "B": (0.4, 0.8),
        "C": (1.2, 0.2),
        "D": (1.5, -0.3),
        "E": (1.7, 0.8),
    }
    _, fast = compute_mvbm_marginals(tree, values, ("x", "y"))
    _, dense = fit_dense_mvbm(tree, values, ("x", "y"))
    assert dense.sigma == pytest.approx(fast.sigma, rel=2e-4, abs=2e-5)
    assert dense.restricted_log_likelihood == pytest.approx(
        fast.restricted_log_likelihood, rel=2e-5, abs=2e-5
    )


@pytest.mark.parametrize(
    "fitter,is_ou",
    [
        (fit_dense_mvbm, False),
        (fit_dense_mvou, True),
        (fit_dense_mvou_diag, True),
    ],
)
def test_dense_covariance_rank_is_invariant_to_trait_units(fitter, is_ou):
    tree = tree_from("((A:.5,B:.7):.4,(C:.6,D:.8):.3,E:1,F:1.2)R;")
    values = {
        "A": (0.0, 1.0),
        "B": (0.4, 0.8),
        "C": (1.2, 0.2),
        "D": (1.5, -0.3),
        "E": (1.7, 0.8),
        "F": (2.3, -0.6),
    }
    errors = {name: (0.05, 0.08) for name in values}
    options = {"alpha": 0.7} if is_ou else {}
    _, original = fitter(tree, values, ("x", "y"), standard_errors=errors, **options)
    multiplier = 1e10
    transformed_values = {
        name: (vector[0] * multiplier, vector[1]) for name, vector in values.items()
    }
    transformed_errors = {
        name: (vector[0] * multiplier, vector[1]) for name, vector in errors.items()
    }
    _, transformed = fitter(
        tree,
        transformed_values,
        ("x", "y"),
        standard_errors=transformed_errors,
        **options,
    )
    assert transformed.sigma_rank == original.sigma_rank == 2
    assert transformed.fit_status == original.fit_status == "ok"
    original_likelihood = (
        original.log_likelihood if is_ou else original.restricted_log_likelihood
    )
    transformed_likelihood = (
        transformed.log_likelihood if is_ou else transformed.restricted_log_likelihood
    )
    likelihood_dimensions = len(values) if is_ou else len(values) - 1
    assert transformed_likelihood == pytest.approx(
        original_likelihood - likelihood_dimensions * math.log(multiplier),
        abs=2e-4,
    )


@pytest.mark.parametrize("fitter", [fit_dense_mvou, fit_dense_mvou_diag])
def test_correlated_mvou_supports_partial_observations(fitter):
    tree = tree_from("((A:0.5,B:0.7)I:0.4,(C:0.6,D:0.8)J:0.3,E:1.0)R;")
    values = {
        "A": (0.0, 1.0),
        "B": (0.4, None),
        "C": (1.2, 0.2),
        "D": (None, -0.3),
        "E": (1.7, 0.8),
    }
    posterior, fit = fitter(tree, values, ("x", "y"), alpha=0.7)
    assert fit.log_likelihood is not None
    assert fit.alpha == 0.7
    assert fit.theta is not None and len(fit.theta) == 2
    assert fit.num_effective_observations == 5
    assert posterior[tree].covariance.shape == (2, 2)


def test_mvou_diag_with_fixed_shared_alpha_reduces_to_mvou():
    tree = tree_from("((A:.5,B:.7):.4,(C:.6,D:.8):.3,E:1,F:1.2)R;")
    values = {
        "A": (0.0, 1.0),
        "B": (0.4, 0.8),
        "C": (1.2, 0.2),
        "D": (1.5, -0.3),
        "E": (1.7, 0.8),
        "F": (2.3, -0.6),
    }
    shared_posterior, shared = fit_dense_mvou(tree, values, ("x", "y"), alpha=0.7)
    diagonal_posterior, diagonal = fit_dense_mvou_diag(
        tree, values, ("x", "y"), alpha=0.7
    )
    assert diagonal.log_likelihood == pytest.approx(shared.log_likelihood, abs=1e-9)
    assert diagonal.sigma == pytest.approx(shared.sigma, abs=1e-6)
    assert diagonal.alpha_by_trait == pytest.approx((0.7, 0.7))
    for node in tree.traverse():
        assert diagonal_posterior[node].mean == pytest.approx(
            shared_posterior[node].mean, abs=1e-8
        )
        assert diagonal_posterior[node].covariance == pytest.approx(
            shared_posterior[node].covariance, abs=1e-6
        )


def test_mvou_diag_accepts_fixed_trait_specific_alpha_and_reports_diffusion():
    tree = tree_from("((A:.5,B:.7):.4,(C:.6,D:.8):.3,E:1,F:1.2)R;")
    values = {
        "A": (0.0, 1.0),
        "B": (0.4, 0.8),
        "C": (1.2, 0.2),
        "D": (1.5, -0.3),
        "E": (1.7, 0.8),
        "F": (2.3, -0.6),
    }
    _posterior, fit = fit_dense_mvou_diag(
        tree, values, ("x", "y"), alpha_by_trait="0.3,1.2"
    )
    assert fit.alpha_by_trait == pytest.approx((0.3, 1.2))
    assert not fit.alpha_estimated
    assert fit.diffusion_sigma is not None
    assert np.all(np.linalg.eigvalsh(fit.diffusion_sigma) > 0.0)
    alpha_sum = fit.alpha_by_trait[:, None] + fit.alpha_by_trait[None, :]
    assert fit.diffusion_sigma == pytest.approx(alpha_sum * fit.sigma, abs=1e-8)


def test_mvou_diag_posterior_matches_direct_dense_conditioning():
    tree = tree_from("((A:.5,B:.7):.4,(C:.6,D:.8):.3,E:1,F:1.2)R;")
    values = {
        "A": (0.0, 1.0),
        "B": (0.4, 0.8),
        "C": (1.2, 0.2),
        "D": (1.5, -0.3),
        "E": (1.7, 0.8),
        "F": (2.3, -0.6),
    }
    posterior, fit = fit_dense_mvou_diag(
        tree, values, ("x", "y"), alpha_by_trait=(0.3, 1.2)
    )
    data = _prepare_observations(tree, values, ("x", "y"))
    geometry = _geometry(data)
    observed_depths = geometry.depths[data.node_indices]
    observed_alpha = fit.alpha_by_trait[data.trait_indices]
    left = observed_depths[:, None] - geometry.observed_shared_depth
    right = observed_depths[None, :] - geometry.observed_shared_depth
    observed_covariance = (
        np.exp(-observed_alpha[:, None] * left - observed_alpha[None, :] * right)
        * fit.sigma[data.trait_indices[:, None], data.trait_indices[None, :]]
    )
    observed_values = data.centers[data.trait_indices] + (
        data.scales[data.trait_indices] * data.values
    )
    observed_mean = fit.theta[data.trait_indices]
    solved_residual = np.linalg.solve(
        observed_covariance, observed_values - observed_mean
    )
    for node_index, node in enumerate(geometry.compiled.nodes):
        shared = geometry.node_shared_depth[node_index]
        node_distance = geometry.depths[node_index] - shared
        observed_distance = observed_depths - shared
        cross = (
            np.exp(
                -fit.alpha_by_trait[:, None] * node_distance
                - observed_alpha[None, :] * observed_distance
            )
            * fit.sigma[:, data.trait_indices]
        )
        expected_mean = fit.theta + cross @ solved_residual
        expected_covariance = fit.sigma - cross @ np.linalg.solve(
            observed_covariance, cross.T
        )
        assert posterior[node].mean == pytest.approx(expected_mean, abs=2e-7)
        assert posterior[node].covariance == pytest.approx(
            expected_covariance, abs=2e-7
        )


def test_parse_alpha_by_trait_validates_dimension_and_positivity():
    assert parse_alpha_by_trait("0.2,1.3", ("x", "y")) == pytest.approx((0.2, 1.3))
    with pytest.raises(ValueError, match="exactly one value"):
        parse_alpha_by_trait("0.2", ("x", "y"))
    with pytest.raises(ValueError, match="strictly positive"):
        parse_alpha_by_trait("0.2,0", ("x", "y"))


@pytest.mark.parametrize(
    "fitter", [fit_dense_mvbm, fit_dense_mvou, fit_dense_mvou_diag]
)
def test_dense_multivariate_models_contract_duplicate_exact_zero_length_tips(fitter):
    tree = tree_from("((A:0,B:0)I:1,C:1,D:1,E:1)R;")
    values = {
        "A": (1.0, 2.0),
        "B": (1.0, 2.0),
        "C": (2.0, 1.0),
        "D": (3.0, 0.0),
        "E": (4.0, -1.0),
    }
    options = {"alpha": 0.7} if fitter is not fit_dense_mvbm else {}
    posterior, fit = fitter(tree, values, ("x", "y"), **options)
    assert fit.num_observed == 5
    assert fit.num_effective_observations == 4
    assert posterior[tree["A"]].mean == pytest.approx(posterior[tree["B"]].mean)
    assert posterior[tree["A"]].covariance == pytest.approx(
        posterior[tree["B"]].covariance
    )


def test_dense_multivariate_model_rejects_conflicting_exact_zero_length_tips():
    tree = tree_from("((A:0,B:0)I:1,C:1,D:1,E:1)R;")
    values = {
        "A": (1.0, 2.0),
        "B": (1.1, 2.0),
        "C": (2.0, 1.0),
        "D": (3.0, 0.0),
        "E": (4.0, -1.0),
    }
    with pytest.raises(ValueError, match="Conflicting exact multivariate"):
        fit_dense_mvou(tree, values, ("x", "y"), alpha=0.7)


def test_dense_mvbm_uses_exact_root_position_to_fix_flat_trait_means():
    tree = tree_from("(A:0,B:1,C:1,D:1,E:1)R;")
    values = {
        "A": (1.0, 2.0),
        "B": (2.0, None),
        "C": (3.0, 0.0),
        "D": (4.0, -1.0),
        "E": (5.0, -2.0),
    }
    posterior, fit = fit_dense_mvbm(tree, values, ("x", "y"))
    assert fit.restricted_log_likelihood is not None
    assert posterior[tree].mean == pytest.approx((1.0, 2.0))
    assert posterior[tree].covariance == pytest.approx(np.zeros((2, 2)))


def test_dense_mvbm_rejects_covariance_without_shared_trait_branch_overlap():
    tree = tree_from("(A:1,B:1,C:1,D:1,E:1,F:1)R;")
    values = {
        "A": (0.0, None),
        "B": (1.0, None),
        "C": (2.0, None),
        "D": (None, 0.0),
        "E": (None, 1.0),
        "F": (None, 2.0),
    }
    with pytest.raises(ValueError, match=r"covariance\(x,y\)"):
        fit_dense_mvbm(tree, values, ("x", "y"))


def test_dense_mvou_rejects_free_alpha_with_constant_trait_pair_distances():
    tree = tree_from("((A:0,B:0,C:0)I:1,(D:0,E:0,F:0)J:1)R;")
    values = {
        "A": (0.0, None),
        "B": (0.2, None),
        "C": (-0.1, None),
        "D": (None, 1.0),
        "E": (None, 1.2),
        "F": (None, 0.9),
    }
    errors = {name: (0.1, 0.1) for name in values}
    with pytest.raises(ValueError, match="alpha is confounded"):
        fit_dense_mvou(
            tree,
            values,
            ("x", "y"),
            standard_errors=errors,
        )
    with pytest.raises(ValueError, match="trait-specific alpha"):
        fit_dense_mvou_diag(
            tree,
            values,
            ("x", "y"),
            standard_errors=errors,
        )


@pytest.mark.parametrize(
    "fitter", [fit_dense_mvbm, fit_dense_mvou, fit_dense_mvou_diag]
)
def test_dense_multivariate_models_reject_too_few_coordinates(fitter):
    tree = tree_from("(A:1,B:1,C:1)R;")
    values = {"A": (0.0, 1.0), "B": (1.0, None), "C": (None, 3.0)}
    options = {"alpha": 0.7} if fitter is not fit_dense_mvbm else {}
    with pytest.raises(ValueError, match="too few effective observed coordinates"):
        fitter(tree, values, ("x", "y"), **options)


def test_dense_multivariate_coordinate_cap_precedes_geometry_allocation():
    names = [f"T{index}" for index in range(501)]
    tree = tree_from("(" + ",".join(f"{name}:1" for name in names) + ")R;")
    values = {name: (float(index), float(-index)) for index, name in enumerate(names)}
    with pytest.raises(ValueError, match="larger than 1000 observed coordinates"):
        _prepare_observations(tree, values, ("x", "y"))


def test_mvou_boundary_warning_retains_available_likelihood():
    fit = SimpleNamespace(
        fit_status="alpha_upper_boundary",
        log_likelihood=-21.5,
        restricted_log_likelihood=None,
    )
    message = _multivariate_fit_status_message("MV-OU", fit)
    assert "parameter-estimation uncertainty" in message
    assert "likelihood is not" not in message


def test_multivariate_warning_treats_nonfinite_likelihood_as_unavailable():
    fit = SimpleNamespace(
        fit_status="singular_covariance",
        log_likelihood=math.nan,
        restricted_log_likelihood=None,
    )
    message = _multivariate_fit_status_message("MV-OU", fit)
    assert "likelihood is not" in message
