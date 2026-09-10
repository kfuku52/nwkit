"""Branch-recursion and dense profiles audit production contrast calculations."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from shift_continuous_reference import DenseOUReference  # noqa: E402
from validate_shift_null_contract import summarize, tree_text  # noqa: E402

from nwkit.shift_calibration import CalibratedSearch  # noqa: E402
from nwkit.util import read_tree  # noqa: E402


@pytest.mark.parametrize(
    "n,shape", [(4, "balanced"), (8, "pectinate"), (16, "balanced")]
)
@pytest.mark.parametrize("alpha", [0.0, 0.4, 30.0, np.inf])
def test_all_candidate_likelihoods_match_branch_recursion(n, shape, alpha):
    tree = read_tree(tree_text(n, shape), "auto", True, quiet=True)
    native = CalibratedSearch(tree, alpha_grid=[alpha])
    values = np.random.default_rng(502 + n).normal(size=n)
    expected, _ = native.profile(native.q @ values)
    for index, model in enumerate(native.models):
        independent = DenseOUReference(tree, model).at(values, alpha)
        assert independent["log_likelihood"] == pytest.approx(
            expected[index, 0], abs=1e-8
        )


@pytest.mark.parametrize("variance_index", [0, 20, 30])
def test_known_error_dense_likelihood_and_boundary(variance_index):
    tree = read_tree(tree_text(4, "pectinate"), "auto", True, quiet=True)
    errors = np.array([0.03, 0.04, 0.1, 0.25])
    native = CalibratedSearch(tree, variances=errors, alpha_grid=[0.7])
    values = np.array([0.1, -1, 0.7, 1.2])
    item = native.cache[variance_index]
    scores, _, _ = native._at((native.q @ values)[:, None], item)
    for index, model in enumerate(native.models):
        independent = DenseOUReference(tree, model, errors).at(values, 0.7, item[1])
        assert independent["log_likelihood"] == pytest.approx(
            scores[index, 0], abs=1e-8
        )


def test_brownian_branch_moments_match_shared_ancestry():
    tree = read_tree(tree_text(8, "pectinate"), "auto", True, quiet=True)
    reference = DenseOUReference(tree, {"groups": [[0]], "shift_branch_ids": []})
    _, covariance = reference.branch_moments(0)
    for i, left in enumerate(reference.tips):
        for j, right in enumerate(reference.tips):
            assert covariance[i, j] == pytest.approx(
                1 - tree.get_distance(left, right) / 2
            )


def test_continuous_profile_can_improve_the_finite_grid():
    tree = read_tree(tree_text(8, "balanced"), "auto", True, quiet=True)
    native = CalibratedSearch(tree)
    values = np.random.default_rng(205).normal(size=8)
    scores, _ = native.profile(native.q @ values)
    reference = DenseOUReference(tree, native.models[100])
    fit = reference.profile(values, grid_points=81)
    assert fit["log_likelihood"] >= scores[100, 0] - 1e-8
    assert not fit["continuous_global_optimum_certified"]


def test_null_summary_keeps_failures_in_worst_case_and_counts_datasets():
    specification = {
        "phase": "validation",
        "cells": [{"cell_id": 0}],
        "convergence_lanes": [False, True],
        "acceptance_upper": 0.075,
    }
    rows = [
        {
            "cell_id": 0,
            "lanes": [
                {"convergence": False, "status": "completed", "any_shift": False},
                {"convergence": True, "status": "failed"},
            ],
        },
        {
            "cell_id": 0,
            "lanes": [
                {"convergence": False, "status": "completed", "any_shift": True},
                {"convergence": True, "status": "completed", "any_shift": False},
            ],
        },
    ]
    result = summarize(rows, specification)
    assert result["independent_datasets"] == 2
    assert result["fits"] == 4
    failed = result["cell_results"][1]
    assert failed["completed_rate"] == 0
    assert failed["worst_case_rate"] == 0.5
    assert failed["failed"] == 1
    assert not result["all_cells_pass"]
