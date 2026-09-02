import math
from types import SimpleNamespace

import numpy as np
import pytest

from nwkit.asr import _multivariate_fit_status_message
from nwkit.multivariate_gaussian_asr import (
    _prepare_observations,
    fit_dense_mvbm,
    fit_dense_mvou,
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
    assert fit.num_effective_observations == 8
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


def test_correlated_mvou_supports_partial_observations():
    tree = tree_from("((A:0.5,B:0.7)I:0.4,(C:0.6,D:0.8)J:0.3,E:1.0)R;")
    values = {
        "A": (0.0, 1.0),
        "B": (0.4, None),
        "C": (1.2, 0.2),
        "D": (None, -0.3),
        "E": (1.7, 0.8),
    }
    posterior, fit = fit_dense_mvou(tree, values, ("x", "y"), alpha=0.7)
    assert fit.log_likelihood is not None
    assert fit.alpha == 0.7
    assert fit.theta is not None and len(fit.theta) == 2
    assert posterior[tree].covariance.shape == (2, 2)


@pytest.mark.parametrize("fitter", [fit_dense_mvbm, fit_dense_mvou])
def test_dense_multivariate_models_contract_duplicate_exact_zero_length_tips(fitter):
    tree = tree_from("((A:0,B:0)I:1,C:1,D:1,E:1)R;")
    values = {
        "A": (1.0, 2.0),
        "B": (1.0, 2.0),
        "C": (2.0, 1.0),
        "D": (3.0, 0.0),
        "E": (4.0, -1.0),
    }
    options = {"alpha": 0.7} if fitter is fit_dense_mvou else {}
    posterior, fit = fitter(tree, values, ("x", "y"), **options)
    assert fit.num_observed == 5
    assert fit.num_effective_observations == 8
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


@pytest.mark.parametrize("fitter", [fit_dense_mvbm, fit_dense_mvou])
def test_dense_multivariate_models_reject_too_few_coordinates(fitter):
    tree = tree_from("(A:1,B:1,C:1)R;")
    values = {"A": (0.0, 1.0), "B": (1.0, None), "C": (None, 3.0)}
    options = {"alpha": 0.7} if fitter is fit_dense_mvou else {}
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
