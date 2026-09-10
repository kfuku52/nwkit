import numpy as np
import pytest

from nwkit.shift_native_fit import fit_native_layout
from nwkit.shift_native_heuristic import NativeSearchOptions, heuristic_native_search
from nwkit.shift_native_model import ShiftData, ShiftLayout
from nwkit.shift_native_quick import NativeQuickProfile
from nwkit.shift_native_screen import _kkt_residual, _proximal_step, group_lasso_screen
from nwkit.shift_native_search import enumerate_native_layouts, exhaustive_native_search
from tests.test_shift_native_search import sample_data


@pytest.mark.parametrize("alpha", [0, 0.7, float("inf")])
def test_cached_candidate_likelihood_matches_refit_for_nested_shared_regimes(alpha):
    base = sample_data()
    values = base.values.copy()
    values[2, 1] = np.nan
    data = ShiftData.build(base.tree, values, base.trait_names)
    arguments = {"alpha_height": alpha, "process_variance": 1.0}
    null = fit_native_layout(data, ShiftLayout.build(data.tree), **arguments)
    branches = tuple(b for b in data.tree.branch_ids if b)
    quick = NativeQuickProfile(data, null, branches, alpha)
    layouts, _ = enumerate_native_layouts(data, 2, convergence=True)
    for layout in layouts:
        try:
            expected = fit_native_layout(data, layout, **arguments)["log_likelihood"]
        except ValueError as exc:
            assert "rank deficient" in str(exc)
            assert quick.score(layout) == -float("inf")
            continue
        assert quick.score(layout) == pytest.approx(expected, abs=1e-9)


def test_group_proximal_solution_matches_orthogonal_closed_form():
    matrices = [np.eye(4), np.eye(4)]
    responses = [np.array([1.0, 2.0, 0.0, -3.0]), np.array([2.0, 0.0, 0.0, 1.0])]
    target = np.column_stack(responses)
    strength = 1.2
    norms = np.linalg.norm(target, axis=1)
    expected = target * np.maximum(0, 1 - strength / np.maximum(norms, 1e-100))[:, None]
    result, gradient, _ = _proximal_step(
        matrices, responses, np.zeros((4, 2)), strength, 1
    )
    np.testing.assert_allclose(result, expected, atol=1e-12)
    assert _kkt_residual(result, gradient, strength) < 1e-12


def test_screening_memory_guard_precedes_allocation():
    data = sample_data()
    with pytest.raises(ValueError, match="approximately"):
        group_lasso_screen(data, {}, memory_limit=1)


def test_one_shift_full_pool_search_agrees_with_exhaustive_and_replays():
    data = sample_data()
    arguments = {"alpha_height": 0.7, "process_variance": 1.0}
    exact = exhaustive_native_search(data, max_shifts=1, fit_arguments=arguments)
    options = NativeSearchOptions(
        max_shifts=1, candidate_pool=14, refit_budget=15, beam_width=14
    )
    first = heuristic_native_search(data, options=options, fit_arguments=arguments)
    second = heuristic_native_search(data, options=options, fit_arguments=arguments)
    assert first.records == second.records
    assert first.metadata == second.metadata
    assert first.families()[-1][1]["log_likelihood"] == pytest.approx(
        exact.families()[-1][1]["log_likelihood"]
    )
    assert not first.metadata["complete_discrete_enumeration"]
    assert len(first.records) <= options.refit_budget
