import numpy as np
import pytest

from nwkit.shift_candidates import enumerate_candidates
from nwkit.shift_native_model import ShiftData, ShiftLayout
from nwkit.shift_native_search import (
    enumerate_native_layouts,
    exhaustive_native_search,
    layout_complexity,
)
from nwkit.util import read_tree
from tests.test_shift_native_model import TREE


def sample_data():
    return ShiftData.build(
        read_tree(TREE, "auto", True, quiet=True),
        np.random.default_rng(48).normal(size=(8, 2)),
        ["x", "y"],
    )


def test_native_reference_contains_existing_two_shift_candidate_space():
    data = sample_data()
    layouts, _ = enumerate_native_layouts(data, 2, convergence=True)
    old, _ = enumerate_candidates(data.tree.compiled.tree, 2)
    expected = {
        ShiftLayout.build(data.tree, row["shift_branch_ids"], row["groups"])
        for row in old
    }
    assert expected <= set(layouts)
    assert all(layout_complexity(row) <= 4 for row in layouts)


def test_three_shift_native_candidates_include_nested_return_and_shared_groups():
    data = sample_data()
    layouts, metadata = enumerate_native_layouts(data, 3, convergence=True, limit=10000)
    expected = ShiftLayout.build(data.tree, [1, 7, 14], [[0, 7], [1, 14]])
    assert expected in layouts
    assert any(len(layout.shifts) == 3 for layout in layouts)
    assert metadata["complete_discrete_enumeration"]
    assert metadata["eligible"] == len(layouts)
    assert len(layouts) == len(set(layouts))


def test_exhaustive_limits_fail_before_any_fits():
    data = sample_data()
    with pytest.raises(ValueError, match="candidate limit"):
        enumerate_native_layouts(data, 2, convergence=True, limit=1)
    with pytest.raises(ValueError, match="traversal budget"):
        enumerate_native_layouts(data, 6, convergence=True, limit=100)


def test_exhaustive_scores_and_nested_families_are_complete_for_fixed_covariance():
    data = sample_data()
    result = exhaustive_native_search(
        data, max_shifts=1, fit_arguments={"alpha_height": 0.8, "process_variance": 1.0}
    )
    assert len(result.records) == 15
    assert result.metadata["complete_discrete_enumeration"]
    families = result.families()
    assert families[0][0] == 0 and families[-1][0] == 2
    assert families[-1][1]["log_likelihood"] >= families[0][1]["log_likelihood"]
    assert families[-1][1]["log_likelihood"] == max(
        row["log_likelihood"] for row in result.records
    )


def test_exact_independent_limit_excludes_only_rank_deficient_layouts():
    base = sample_data()
    values = base.values.copy()
    values[2, 1] = np.nan
    data = ShiftData.build(base.tree, values, base.trait_names)
    result = exhaustive_native_search(
        data,
        max_shifts=2,
        convergence=True,
        fit_arguments={"alpha_height": float("inf"), "process_variance": 1},
    )
    excluded = [row for row in result.records if row["log_likelihood"] is None]
    assert excluded
    assert all(
        row["status"] == "structurally_excluded_at_fixed_alpha" for row in excluded
    )
    assert np.isfinite(result.families()[-1][1]["log_likelihood"])
