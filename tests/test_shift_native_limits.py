"""Automatic caps retain integer semantics and explicit search-budget limits."""

from types import SimpleNamespace

import numpy as np
import pytest

from nwkit.shift import _validate_options
from nwkit.shift_cli import maximum_shift_count
from nwkit.shift_native_limits import native_shift_limit
from nwkit.shift_native_model import ShiftData
from nwkit.shift_native_search import enumerate_native_layouts
from nwkit.shift_native_selection import NativeSearchRunner
from nwkit.util import read_tree
from tests.test_shift_native_search import sample_data


def arguments(**overrides):
    return SimpleNamespace(
        **{
            "max_shifts": "auto",
            "convergence": True,
            "criterion": "AICc",
            "search_strategy": "auto",
            "exhaustive_max_configurations": 5000,
            "candidate_pool": 24,
            "refit_budget": 48,
            "screening_budget": 2000,
            "beam_width": 2,
            "lasso_iterations": 150,
            "search_memory_mb": 512,
            **overrides,
        }
    )


@pytest.mark.parametrize(
    "tips,strategy,pool,refits,expected",
    [
        (8, "auto", 24, 48, 6),
        (1000, "auto", 24, 48, 24),
        (1000, "lasso", 128, 220, 128),
        (1000, "native-path", 1, 220, 219),
        (8, "lasso", 24, 4, 3),
        (8, "lasso", 24, 1, 0),
    ],
)
def test_auto_cap_uses_only_applicable_search_constraints(
    tips, strategy, pool, refits, expected
):
    data = SimpleNamespace(tree=SimpleNamespace(leaf_names=range(tips)))
    args = arguments(search_strategy=strategy, candidate_pool=pool, refit_budget=refits)
    cap, record = native_shift_limit(data, args, budgeted=True)
    assert cap == expected
    assert record["requested"] == "auto"
    assert record["resolved"] == expected
    assert record["budget_limited"] == (expected < tips - 2)
    assert ("candidate_pool" in record["constraints"]) == (strategy != "native-path")
    assert args.max_shifts == "auto"


def test_explicit_limit_is_not_silently_clipped_to_budget():
    data = sample_data()
    args = arguments(
        max_shifts=5, candidate_pool=2, refit_budget=3, search_strategy="lasso"
    )
    assert native_shift_limit(data, args, budgeted=True)[0] == 5
    with pytest.raises(ValueError, match="candidate pool covering"):
        NativeSearchRunner(data, args, {})


@pytest.mark.parametrize("value", [-1, 7, True, 2.5, "bad"])
def test_invalid_explicit_limits_remain_errors(value):
    with pytest.raises(ValueError, match="Maximum shifts"):
        native_shift_limit(sample_data(), arguments(max_shifts=value))


def test_large_automatic_enumeration_stops_before_expensive_bell_numbers(monkeypatch):
    data = SimpleNamespace(
        tree=SimpleNamespace(leaf_names=range(1000), branch_ids=range(1999))
    )
    original = __import__(
        "nwkit.shift_native_search", fromlist=["_bell_number"]
    )._bell_number
    visited = []

    def bounded_bell(size):
        visited.append(size)
        assert size <= 3
        return original(size)

    monkeypatch.setattr("nwkit.shift_native_search._bell_number", bounded_bell)
    with pytest.raises(ValueError, match="traversal budget"):
        enumerate_native_layouts(data, 998, convergence=True)
    assert visited == [1, 2, 3]


def test_small_auto_space_is_exhaustive_even_with_small_heuristic_budgets():
    data = ShiftData.build(
        read_tree("((a:1,b:1):1,(c:1,d:1):1);", "auto", True, quiet=True),
        np.array([[1.0, 2.0], [2.0, 1.0], [4.0, 3.0], [3.0, 4.0]]),
        ["x", "y"],
    )
    runner = NativeSearchRunner(
        data,
        arguments(candidate_pool=1, refit_budget=1),
        {"alpha_height": 0.7, "process_variance": 1.0},
    )
    result = runner(data)
    assert result.metadata["strategy"] == "exhaustive"
    assert result.metadata["shift_limit"]["resolved"] == 2
    assert result.metadata["shift_limit"]["constraints"] == {"tree": 2}
    assert result.metadata["complete_discrete_enumeration"]


def test_parser_keeps_integer_values_and_accepts_auto():
    assert maximum_shift_count("auto") == "auto"
    assert maximum_shift_count("0") == 0
    assert maximum_shift_count("100") == 100


@pytest.mark.parametrize("selection", ["ic", "calibrated"])
def test_auto_is_rejected_by_legacy_backends_before_fitting(selection):
    with pytest.raises(ValueError, match="requires --selection native"):
        _validate_options(SimpleNamespace(selection=selection, max_shifts="auto"))


def test_auto_path_ignores_candidate_pool_and_respects_resolved_refit_cap():
    data = sample_data()
    runner = NativeSearchRunner(
        data,
        arguments(
            search_strategy="native-path",
            convergence=False,
            candidate_pool=1,
            refit_budget=4,
        ),
        {"alpha_height": 0.7, "process_variance": 1.0},
    )
    result = runner(data)
    assert result.metadata["shift_limit"]["resolved"] == 3
    assert "candidate_pool" not in result.metadata["shift_limit"]["constraints"]
    assert len(result.records) <= 4
    assert all(len(row["shift_branch_ids"]) <= 3 for row in result.records)
