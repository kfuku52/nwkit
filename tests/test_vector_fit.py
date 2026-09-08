import numpy as np
import pytest

from nwkit.multivariate_gaussian_asr import fit_dense_mvbm
from nwkit.util import read_tree
from nwkit.vector_fit import fit_pruning_mvbm


@pytest.mark.parametrize("scale", [np.ones(2), np.array([1e-8, 1e8])])
def test_pruning_bm_matches_dense_fit_with_missing_errors(scale):
    tree = read_tree(
        "((A:1,B:2):1,(C:1,D:1):2,E:3,F:2)R;", "1", True, quiet=True, rooted="yes"
    )
    values = {
        "A": [0, 1],
        "B": [1, None],
        "C": [3, 4],
        "D": [4, 2],
        "E": [-2, 3],
        "F": [1, -1],
    }
    observed = {
        name: [None if v is None else v * scale[k] for k, v in enumerate(vector)]
        for name, vector in values.items()
    }
    errors = {name: np.array([0.1, 0.2]) * scale for name in observed}
    dense, fit = fit_dense_mvbm(tree, observed, ("x", "y"), standard_errors=errors)
    sparse, actual = fit_pruning_mvbm(
        tree, observed, ("x", "y"), standard_errors=errors
    )
    assert actual.restricted_log_likelihood == pytest.approx(
        fit.restricted_log_likelihood, abs=1e-7
    )
    assert actual.sigma / np.outer(scale, scale) == pytest.approx(
        fit.sigma / np.outer(scale, scale), abs=2e-4
    )
    for node in dense:
        assert sparse[node].mean / scale == pytest.approx(
            dense[node].mean / scale, abs=1e-4
        )
        assert sparse[node].covariance / np.outer(scale, scale) == pytest.approx(
            dense[node].covariance / np.outer(scale, scale), abs=1e-4
        )


def test_pruning_bm_rejects_unidentifiable_trait_overlap():
    tree = read_tree(
        "((A:1,B:1):1,(C:1,D:1):1,E:1,F:1)R;", "1", True, quiet=True, rooted="yes"
    )
    observed = {
        "A": [0, None],
        "B": [1, None],
        "C": [None, 2],
        "D": [None, 3],
        "E": [2, None],
        "F": [None, 4],
    }
    with pytest.raises(ValueError, match="not identifiable"):
        fit_pruning_mvbm(tree, observed, ("x", "y"))


def test_large_noisy_bm_routes_to_pruning(monkeypatch):
    from nwkit import vector_fit
    from nwkit.multivariate_asr import compute_mvbm_marginals

    tree = read_tree(
        "(" + ",".join(f"T{i}:1" for i in range(501)) + ")R;",
        "1",
        True,
        quiet=True,
        rooted="yes",
    )
    values = {f"T{i}": [float(i), float(i % 7)] for i in range(501)}
    errors = {name: [0.1, 0.2] for name in values}
    seen = []

    def record(*args, **kwargs):
        seen.append((args, kwargs))
        return {}, "pruning"

    monkeypatch.setattr(vector_fit, "fit_pruning_mvbm", record)
    assert compute_mvbm_marginals(tree, values, ("x", "y"), standard_errors=errors) == (
        {},
        "pruning",
    )
    assert seen[0][1]["standard_errors"] is errors


def test_pruning_observation_preparation_has_no_dense_limit():
    from nwkit.multivariate_gaussian_asr import _prepare_observations
    from nwkit.vector_fit import normalized_vector_observations
    from nwkit.vector_gaussian import condition_vector_tree
    from nwkit.vector_processes import vector_brownian_process

    tree = read_tree(
        "(" + ",".join(f"T{i}:1" for i in range(501)) + ")R;",
        "1",
        True,
        quiet=True,
        rooted="yes",
    )
    values = {f"T{i}": [float(i), float(i % 7)] for i in range(501)}
    data = _prepare_observations(tree, values, ("x", "y"), dense_limit=False)
    observed, errors = normalized_vector_observations(data)
    result = condition_vector_tree(
        vector_brownian_process(tree, np.eye(2)), observed, error_covariances=errors
    )
    assert len(data.values) == 1002
    assert np.isfinite(result.log_likelihood)
    assert result.means.shape == (502, 2)


@pytest.mark.parametrize("diagonal", [False, True])
def test_pruning_ou_matches_dense_fit(diagonal):
    from nwkit.multivariate_gaussian_asr import fit_dense_mvou, fit_dense_mvou_diag
    from nwkit.vector_ou_fit import fit_pruning_mvou

    tree = read_tree(
        "((A:1,B:2):1,(C:1,D:1):2,E:3,F:2)R;", "1", True, quiet=True, rooted="yes"
    )
    observed = {
        "A": [0, 1],
        "B": [1, None],
        "C": [3, 4],
        "D": [4, 2],
        "E": [-2, 3],
        "F": [1, -1],
    }
    errors = {name: [0.1, 0.2] for name in observed}
    options = {"alpha_by_trait": [0.3, 0.7]} if diagonal else {"alpha": 0.4}
    fitter = fit_dense_mvou_diag if diagonal else fit_dense_mvou
    dense, expected = fitter(
        tree, observed, ("x", "y"), standard_errors=errors, **options
    )
    sparse, actual = fit_pruning_mvou(
        tree, observed, ("x", "y"), standard_errors=errors, diagonal=diagonal, **options
    )
    assert actual.log_likelihood == pytest.approx(expected.log_likelihood, abs=1e-7)
    assert actual.sigma == pytest.approx(expected.sigma, abs=2e-4)
    assert actual.theta == pytest.approx(expected.theta, abs=1e-4)
    for node in dense:
        assert sparse[node].mean == pytest.approx(dense[node].mean, abs=1e-4)
        assert sparse[node].covariance == pytest.approx(
            dense[node].covariance, abs=2e-4
        )


@pytest.mark.parametrize("diagonal", [False, True])
def test_large_ou_routes_all_consumers_to_pruning(monkeypatch, diagonal):
    from nwkit import vector_ou_fit
    from nwkit.multivariate_gaussian_asr import fit_dense_mvou, fit_dense_mvou_diag

    seen = []

    def record(*args, **kwargs):
        seen.append(kwargs)
        return {}, "pruning"

    monkeypatch.setattr(vector_ou_fit, "fit_pruning_mvou", record)
    fitter = fit_dense_mvou_diag if diagonal else fit_dense_mvou
    assert fitter(
        None, {str(i): [1.0, 2.0] for i in range(501)}, ("x", "y"), alpha=0.4
    ) == ({}, "pruning")
    assert seen[0].get("diagonal", False) is diagonal


def test_sparse_ou_geometry_ignores_large_stem():
    from nwkit.multivariate_gaussian_asr import _prepare_observations
    from nwkit.vector_ou_fit import _ou_geometry

    tree = read_tree("((A:1,B:2):1e20)R;", "1", True, quiet=True, rooted="yes")
    data = _prepare_observations(
        tree, {"A": [0, 1], "B": [1, 2]}, ("x", "y"), dense_limit=False
    )
    assert _ou_geometry(data) == (1.5, [True, True])
