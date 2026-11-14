from types import SimpleNamespace

import numpy as np
import pytest

from nwkit import optimization


def test_scalar_search_refines_best_basin_even_after_four_other_minima():
    centers = np.array([1.0, 3.0, 5.0, 7.0, 9.25])
    offsets = np.array([0.0, 0.0, 0.0, 0.0, -10.0])
    result = optimization.global_bounded_scalar_minimize(
        lambda x: float(np.min(100 * (x - centers) ** 2 + offsets)), (0, 16)
    )
    assert result.success
    assert result.x == pytest.approx(9.25, abs=1e-6)
    assert result.fun == pytest.approx(-10.0)


def test_scalar_search_does_not_borrow_convergence_from_another_basin(monkeypatch):
    results = iter(
        [
            SimpleNamespace(x=0.51, fun=-1.0, success=False),
            *[SimpleNamespace(x=0.75, fun=1.0, success=True)] * 3,
        ]
    )
    monkeypatch.setattr(
        optimization, "_bounded_scalar_minimize", lambda *args: next(results)
    )
    result = optimization.global_bounded_scalar_minimize(
        lambda x: (x - 0.5) ** 2, (0, 1)
    )
    assert result.x == 0.51
    assert result.fun == -1.0
    assert not result.success


def test_scalar_search_accepts_a_flat_profile_confirmed_in_the_same_interval():
    result = optimization.global_bounded_scalar_minimize(
        lambda x: max(abs(x - 0.5) - 0.25, 0.0) ** 2, (0, 1)
    )
    assert result.success
    assert result.fun == 0.0
    assert 0.25 <= result.x <= 0.75


@pytest.mark.parametrize("center", [0.0, 0.5, 1.0])
def test_scalar_search_keeps_convergence_at_exact_grid_and_boundary_minima(center):
    result = optimization.global_bounded_scalar_minimize(
        lambda x: (x - center) ** 2, (0, 1)
    )
    assert result.success
    assert result.x == pytest.approx(center, abs=1e-7)
