"""Independent dense information determinants and native criterion selection."""

import math

import numpy as np
import pytest

from nwkit.shift_native_fit import NativeFitOptions, fit_native_layout
from nwkit.shift_native_ic import native_information_criterion
from nwkit.shift_native_model import ShiftData, ShiftLayout
from nwkit.shift_native_search import exhaustive_native_search
from nwkit.util import read_tree
from tests.test_shift_native_model import TREE, dense_reference


@pytest.mark.parametrize("root", ["OUfixedRoot", "OUrandomRoot"])
@pytest.mark.parametrize("alpha", [1e-7, 0.8, float("inf")])
@pytest.mark.parametrize("groups", [None, [[0], [7, 14]]])
def test_information_scores_match_dense_optimum_information(root, alpha, groups):
    y = np.random.default_rng(29).normal(size=(8, 2)) * [2, 3] + [5, -2]
    y[2, 1] = np.nan
    data = ShiftData.build(read_tree(TREE, "auto", True, quiet=True), y, ["x", "y"])
    layout = ShiftLayout.build(data.tree, [7, 14], groups)
    result = fit_native_layout(
        data,
        layout,
        alpha_height=alpha,
        process_variance=[2, 3],
        options=NativeFitOptions(root_model=root),
    )
    penalty = 4 * math.log(13)
    for j, fit in enumerate(result["fits"]):
        x, v = dense_reference(
            data.tree, layout, alpha, fit.process_variance, np.zeros(8), root
        )
        x[:, 1:] *= -math.expm1(-alpha)
        x[:, 0] -= x[:, 1:].sum(axis=1)
        mask = np.isfinite(data.values[:, j])
        x, v = x[mask], v[np.ix_(mask, mask)]
        xw = np.linalg.solve(np.linalg.cholesky(v), x)
        # QR avoids normal-equation cancellation near alpha=0.
        logdet = 2 * np.log(np.abs(np.diag(np.linalg.qr(xw, mode="reduced")[1]))).sum()
        penalty += logdet + len(layout.groups) * np.log(
            np.var(data.values[mask, j], ddof=1)
        )
    actual = native_information_criterion(data, result, "pBIC")
    assert actual["score"] == pytest.approx(
        -2 * result["log_likelihood"] + penalty, abs=2e-7
    )
    for criterion in ("AIC", "BIC"):
        expected = (
            4 + 4 * len(layout.groups)
            if criterion == "AIC"
            else 2 * math.log(8) + len(layout.groups) * (math.log(8) + math.log(7))
        )
        assert native_information_criterion(data, result, criterion)[
            "penalty"
        ] == pytest.approx(expected)


@pytest.mark.parametrize("criterion", ["AIC", "AICc", "BIC", "pBIC"])
def test_exhaustive_selection_minimizes_all_candidate_scores(criterion):
    data = ShiftData.build(
        read_tree(TREE, "auto", True, quiet=True),
        np.random.default_rng(3).normal(size=(8, 1)),
        ["x"],
    )
    searched = exhaustive_native_search(
        data, max_shifts=1, criterion=criterion, fit_arguments={"alpha_height": 0.8}
    )
    best = searched.best_information["information_criterion"]["score"]
    assert best == min(
        row["information_criterion"]["score"] for row in searched.records
    )
    assert searched.best_information["information_criterion"][
        "parameter_count"
    ] == 1 + len(searched.best_information["layout"].groups) + len(
        searched.best_information["layout"].shifts
    )


def test_pbic_retains_better_penalized_layout_within_same_complexity():
    data = ShiftData.build(
        read_tree(TREE, "auto", True, quiet=True),
        np.random.default_rng(14).normal(size=(8, 1)),
        ["x"],
    )
    searched = exhaustive_native_search(
        data, max_shifts=1, criterion="pBIC", fit_arguments={"alpha_height": 0.8}
    )
    assert searched.best_by_complexity[2]["layout"].shifts == (4,)
    assert searched.best_information["layout"].shifts == (10,)
    assert (
        searched.best_information["information_criterion"]["score"]
        < searched.best_by_complexity[2]["information_criterion"]["score"]
    )


def test_zero_alpha_shift_has_no_finite_optimum_information():
    data = ShiftData.build(
        read_tree(TREE, "auto", True, quiet=True),
        np.random.default_rng(14).normal(size=(8, 1)),
        ["x"],
    )
    shifted = fit_native_layout(data, ShiftLayout.build(data.tree, [7]), alpha_height=0)
    null = fit_native_layout(data, ShiftLayout.build(data.tree), alpha_height=0)
    assert native_information_criterion(data, shifted, "pBIC")["score"] is None
    assert native_information_criterion(data, null, "pBIC")["score"] is not None


@pytest.mark.parametrize("shared", [False, True])
def test_aicc_counts_shared_locations_once_and_observed_coordinates(shared):
    values = np.random.default_rng(52).normal(size=(8, 2))
    values[2, 1] = np.nan
    data = ShiftData.build(
        read_tree(TREE, "auto", True, quiet=True), values, ["x", "y"]
    )
    layout = ShiftLayout.build(data.tree, [7, 14], [[0], [7, 14]] if shared else None)
    result = fit_native_layout(data, layout, alpha_height=0.8, process_variance=[1, 2])
    score = native_information_criterion(data, result, "AICc")
    parameters = 2 + 2 * (2 if shared else 3)
    correction = 2 * parameters * (parameters + 1) / (15 - parameters - 1)
    assert score["sample_size"] == 15
    assert score["parameter_count"] == parameters
    assert score["small_sample_correction"] == pytest.approx(correction)
    assert score["score"] == pytest.approx(
        -2 * result["log_likelihood"] + 2 * parameters + correction
    )


def test_aicc_excludes_nonpositive_denominator_and_counts_estimated_parameters():
    data = ShiftData.build(
        read_tree(TREE, "auto", True, quiet=True), np.arange(8)[:, None], ["x"]
    )
    layout = ShiftLayout.build(data.tree, [7, 14])
    estimated = fit_native_layout(data, layout)
    invalid = native_information_criterion(data, estimated, "AICc")
    assert invalid["parameter_count"] == 7
    assert invalid["sample_size"] == 8
    assert invalid["score"] is None
    assert invalid["status"] == "insufficient_aicc_sample_size"
    fixed = fit_native_layout(data, layout, alpha_height=0.8)
    valid = native_information_criterion(data, fixed, "AICc")
    assert valid["parameter_count"] == 6
    assert valid["small_sample_correction"] == 84
    assert valid["score"] == pytest.approx(-2 * fixed["log_likelihood"] + 12 + 84)
