"""Sparse candidate paths against independent mean/covariance geometry."""

import numpy as np
import pytest

from nwkit.shift_native_fit import NativeFitOptions, fit_native_layout
from nwkit.shift_native_heuristic import NativeSearchOptions
from nwkit.shift_native_model import ShiftData, ShiftLayout, ShiftTree
from nwkit.shift_native_path import sparse_native_search
from nwkit.shift_native_screen import (
    _effect_design,
    _whitened_matrices,
    descendant_design,
)
from nwkit.util import read_tree
from tests.test_shift_native_model import TREE, dense_reference


@pytest.mark.parametrize("alpha", [0.0, 0.7, 1000.0, np.inf])
@pytest.mark.parametrize("newick", [TREE, TREE.replace("a:1", "a:1.000000006")])
def test_optimum_increment_columns_match_independent_branch_recursion(alpha, newick):
    tree = ShiftTree.build(read_tree(newick, "auto", True, quiet=True))
    design, branches = descendant_design(tree)
    actual = _effect_design(tree, design, alpha)
    for column, branch in enumerate(branches):
        expected, _ = dense_reference(
            tree,
            ShiftLayout.build(tree, [branch]),
            alpha,
            1.0,
            np.zeros(8),
            "OUfixedRoot",
        )
        np.testing.assert_allclose(actual[:, column], expected[:, 1], atol=1e-12)


def data_fixture():
    y = np.random.default_rng(836).normal(size=(8, 2))
    y[0, 1] = np.nan
    return ShiftData.build(
        read_tree(TREE, "auto", True, quiet=True), y, ["x", "y"], np.full((8, 2), 0.04)
    )


def test_unstandardized_whitening_matches_dense_gls_with_missingness_and_errors():
    data = data_fixture()
    fit = fit_native_layout(
        data,
        ShiftLayout.build(data.tree),
        alpha_height=[0.7, 1.3],
        process_variance=[1, 2],
    )
    matrices, responses, branches = _whitened_matrices(
        data, fit, 1024**2, optimum_increments=True
    )
    for j, trait in enumerate(fit["fits"]):
        mask = np.isfinite(data.values[:, j])
        for column, branch in enumerate(branches):
            x, v = dense_reference(
                data.tree,
                ShiftLayout.build(data.tree, [branch]),
                trait.alpha_height,
                trait.process_variance,
                data.variances[:, j],
                "OUfixedRoot",
            )
            w = np.linalg.cholesky(v[np.ix_(mask, mask)])
            white = np.linalg.solve(w, np.column_stack((x[mask], data.values[mask, j])))
            intercept = white[:, 0] / np.linalg.norm(white[:, 0])
            projected = (
                white[:, 1:] - intercept[:, None] * (intercept @ white[:, 1:])[None, :]
            )
            actual = np.column_stack((matrices[j][:, column], responses[j]))
            np.testing.assert_allclose(
                actual.T @ actual, projected.T @ projected, atol=1e-10
            )


@pytest.mark.parametrize("root", ["OUfixedRoot", "OUrandomRoot"])
def test_path_replays_and_preserves_fixed_parameters_and_budget(root):
    data = data_fixture()
    options = NativeSearchOptions(max_shifts=2, refit_budget=8, screening_budget=40)
    arguments = {
        "alpha_height": [0.7, 1.3],
        "process_variance": [1, 2],
        "options": NativeFitOptions(root_model=root),
    }
    first = sparse_native_search(data, options=options, fit_arguments=arguments)
    second = sparse_native_search(data, options=options, fit_arguments=arguments)
    assert first.records == second.records
    assert first.metadata == second.metadata
    assert first.metadata["covariance_seed_fits"] == 0
    assert len(first.records) <= options.refit_budget
    assert first.metadata["screening_evaluations"] <= options.screening_budget
    assert first.best_information["information_criterion"]["score"] == min(
        r["information_criterion"]["score"] for r in first.records
    )
    for result in first.best_by_complexity.values():
        np.testing.assert_allclose([f.alpha_height for f in result["fits"]], [0.7, 1.3])


def test_free_alpha_covariance_seed_is_counted_in_refit_budget():
    data = data_fixture()
    options = NativeSearchOptions(max_shifts=2, refit_budget=3, screening_budget=3)
    result = sparse_native_search(data, options=options)
    assert len(result.records) + result.metadata["covariance_seed_fits"] <= 3
    assert result.metadata["screening_evaluations"] <= 3
    assert result.metadata["budget_exhausted"]


def test_path_memory_guard_precedes_design_allocation():
    data = data_fixture()
    with pytest.raises(ValueError, match="approximately"):
        sparse_native_search(data, options=NativeSearchOptions(memory_limit=1))


def test_all_branch_path_does_not_require_unused_beam_pool_to_cover_cap():
    result = sparse_native_search(
        data_fixture(),
        options=NativeSearchOptions(max_shifts=2, candidate_pool=1, refit_budget=8),
        fit_arguments={"alpha_height": 0.7},
    )
    assert result.best_information is not None


def test_r_backend_rejects_native_path_before_execution():
    from types import SimpleNamespace

    from nwkit.shift import _validate_options

    with pytest.raises(ValueError, match="requires --selection native"):
        _validate_options(
            SimpleNamespace(search_strategy="native-path", selection="ic")
        )


@pytest.mark.parametrize(
    "criterion,convergence", [(None, False), ("BIC", False), ("AIC", True)]
)
def test_path_rejects_unvalidated_selection_modes(criterion, convergence):
    with pytest.raises(ValueError, match="AIC or AICc without convergence"):
        sparse_native_search(
            data_fixture(),
            options=NativeSearchOptions(convergence=convergence),
            criterion=criterion,
        )


def test_native_path_supports_aicc_and_reports_valid_minimum():
    from tests.test_shift_native_search import sample_data

    result = sparse_native_search(
        sample_data(), options=NativeSearchOptions(max_shifts=1), criterion="AICc"
    )
    scores = [
        r["information_criterion"]["score"]
        for r in result.records
        if r["information_criterion"]["score"] is not None
    ]
    assert result.best_information["information_criterion"]["score"] == min(scores)
    assert result.best_information["information_criterion"]["criterion"] == "AICc"


def test_aicc_path_rejects_a_null_with_too_few_observations():
    tree = read_tree("((a:1,b:1):1,(c:1,d:1):1);", "auto", True, quiet=True)
    data = ShiftData.build(tree, np.array([[1.0], [2.0], [4.0], [3.0]]), ["x"])
    with pytest.raises(
        ValueError, match="No candidate has finite native information criterion"
    ):
        sparse_native_search(
            data, options=NativeSearchOptions(max_shifts=0), criterion="AICc"
        )
