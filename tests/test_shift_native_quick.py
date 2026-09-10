"""Orthogonal candidate compression against independent full-layout fitting."""

import numpy as np
import pytest

from nwkit.shift_native_fit import fit_native_layout
from nwkit.shift_native_model import ShiftData, ShiftLayout
from nwkit.shift_native_quick import NativeQuickProfile
from nwkit.util import read_tree


def balanced_tree(size):
    def subtree(start, count):
        if count == 1:
            return f"t{start}", 0
        left, lh = subtree(start, count // 2)
        right, rh = subtree(start + count // 2, count - count // 2)
        height = max(lh, rh) + 1
        return f"({left}:{height - lh},{right}:{height - rh})", height

    return subtree(0, size)[0] + ";"


@pytest.mark.parametrize("alpha", [0, 0.7, float("inf")])
@pytest.mark.parametrize("shared", [False, True])
def test_thousand_tip_hundred_shift_compression_matches_full_fit(alpha, shared):
    values = np.random.default_rng(940).normal(size=(1000, 2))
    values[900:950, 1] = np.nan
    errors = np.linspace(0, 0.2, 1000)[:, None] * np.ones((1, 2))
    data = ShiftData.build(
        read_tree(balanced_tree(1000), "auto", True, quiet=True),
        values,
        ["x", "y"],
        errors,
    )
    branches = tuple(data.tree.branch_ids[i] for i in data.tree.compiled.leaf_indices)
    pool = branches[:128]
    shifts = pool[:100]
    groups = [[0], list(shifts[:50]), list(shifts[50:])] if shared else None
    layout = ShiftLayout.build(data.tree, shifts, groups)
    arguments = {"alpha_height": alpha, "process_variance": 1.0}
    null = fit_native_layout(data, ShiftLayout.build(data.tree), **arguments)
    quick = NativeQuickProfile(data, null, pool, alpha)
    assert all(matrix.shape[0] <= 129 for _, _, matrix, _, _ in quick.profiles[0])
    expected = fit_native_layout(data, layout, **arguments)["log_likelihood"]
    assert quick.score(layout) == pytest.approx(expected, rel=1e-12, abs=1e-8)


def test_compressed_profiles_keep_unobservable_shift_rank_deficient():
    values = np.random.default_rng(41).normal(size=(64, 1))
    values[0] = np.nan
    data = ShiftData.build(
        read_tree(balanced_tree(64), "auto", True, quiet=True), values, ["x"]
    )
    branches = tuple(data.tree.branch_ids[i] for i in data.tree.compiled.leaf_indices)
    null = fit_native_layout(
        data, ShiftLayout.build(data.tree), alpha_height=0.7, process_variance=1.0
    )
    quick = NativeQuickProfile(data, null, branches[:8], 0.7)
    assert quick.score(ShiftLayout.build(data.tree, [branches[0]])) == -float("inf")


@pytest.mark.parametrize("alpha", [0, 0.7, float("inf")])
def test_compressed_nested_return_and_shared_regime_matches_full_fit(alpha):
    values = np.random.default_rng(93).normal(size=(128, 2))
    values[100:105, 1] = np.nan
    data = ShiftData.build(
        read_tree(balanced_tree(128), "auto", True, quiet=True), values, ["x", "y"]
    )
    tips = tuple(data.tree.branch_ids[i] for i in data.tree.compiled.leaf_indices)
    branches = (1, tips[0], tips[80])
    layout = ShiftLayout.build(data.tree, branches, [[0, tips[0]], [1, tips[80]]])
    arguments = {"alpha_height": alpha, "process_variance": 1.0}
    null = fit_native_layout(data, ShiftLayout.build(data.tree), **arguments)
    quick = NativeQuickProfile(data, null, branches, alpha)
    expected = fit_native_layout(data, layout, **arguments)["log_likelihood"]
    assert quick.score(layout) == pytest.approx(expected, abs=1e-9)


def test_rank_tolerance_uses_observed_tips_after_row_reduction():
    values = np.random.default_rng(18).normal(size=(64, 1))
    data = ShiftData.build(
        read_tree(balanced_tree(64), "auto", True, quiet=True), values, ["x"]
    )
    branches = tuple(data.tree.branch_ids[i] for i in data.tree.compiled.leaf_indices)[
        :2
    ]
    null = fit_native_layout(
        data, ShiftLayout.build(data.tree), alpha_height=0.7, process_variance=1.0
    )
    quick = NativeQuickProfile(data, null, branches, 0.7)
    # An orthogonally reduced design may have only three rows although it
    # represents 64 observations. Its second direction is below the original
    # rank cutoff, but above a cutoff incorrectly computed from three rows.
    matrix = np.array([[1.0, 1.0], [0.0, 1e-14], [0.0, 0.0]])
    assert np.linalg.matrix_rank(matrix) == 2
    padded = np.pad(matrix, ((0, 61), (0, 0)))
    assert np.linalg.matrix_rank(padded) == 1
    quick.profiles = [((float("inf"), np.array([0.0, 0.0, 1.0]), matrix, 0.0, 64),)]
    assert quick.score(ShiftLayout.build(data.tree, branches)) == -float("inf")


def test_compressed_multi_alpha_ranking_matches_independent_full_fits():
    data = ShiftData.build(
        read_tree(balanced_tree(128), "auto", True, quiet=True),
        np.random.default_rng(79).normal(size=(128, 2)),
        ["x", "y"],
    )
    tips = tuple(data.tree.branch_ids[i] for i in data.tree.compiled.leaf_indices)
    branches = (1, tips[0], tips[80])
    layout = ShiftLayout.build(data.tree, branches, [[0, tips[0]], [1, tips[80]]])
    null = fit_native_layout(
        data, ShiftLayout.build(data.tree), alpha_height=0.7, process_variance=1.0
    )
    quick = NativeQuickProfile(data, null, branches)
    expected = max(
        fit_native_layout(data, layout, alpha_height=alpha, process_variance=1.0)[
            "log_likelihood"
        ]
        for alpha in (0.7, 0.3, 3.0)
    )
    assert quick.score(layout) == pytest.approx(expected, abs=1e-9)
