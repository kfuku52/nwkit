"""Blocked screening retains the same geometry and branch ordering."""

import numpy as np
import pytest

from nwkit import shift_native_screen as screen
from nwkit.shift_native_fit import fit_native_layout
from nwkit.shift_native_model import ShiftLayout
from tests.test_shift_native_path import data_fixture


@pytest.mark.parametrize("alpha", [0.0, 0.7, 1000.0, np.inf])
@pytest.mark.parametrize("increments", [False, True])
def test_blocked_whitening_matches_whole_matrix(alpha, increments, monkeypatch):
    data = data_fixture()
    fit = fit_native_layout(
        data, ShiftLayout.build(data.tree), alpha_height=alpha, process_variance=1.0
    )
    monkeypatch.setattr(screen, "_screen_block_width", lambda factor: 100)
    expected_x, expected_y, expected_branches = screen._whitened_matrices(
        data, fit, 1024**2, optimum_increments=increments
    )
    for width in (1, 3):
        monkeypatch.setattr(
            screen, "_screen_block_width", lambda factor, width=width: width
        )
        actual_x, actual_y, branches = screen._whitened_matrices(
            data, fit, 1024**2, optimum_increments=increments
        )
        assert branches == expected_branches
        for actual, expected in zip(
            actual_x + actual_y, expected_x + expected_y, strict=True
        ):
            np.testing.assert_allclose(actual, expected, atol=1e-12, rtol=1e-12)


def test_selected_descendant_columns_keep_requested_order():
    data = data_fixture()
    whole, branches = screen.descendant_design(data.tree)
    selected = branches[::-3]
    actual, order = screen.descendant_design(data.tree, selected)
    assert order == selected
    np.testing.assert_array_equal(
        actual, whole[:, [branches.index(b) for b in selected]]
    )
