"""Horizontal tracks shrink uniformly only when geometric clearance requires it."""

import pytest

from nwkit.reconciliation_spacing import horizontal_pitch, horizontal_tracks


def test_two_and_three_turn_groups_have_the_same_step():
    pitch = horizontal_pitch(100, [3, 5], [(2, 32, 6), (3, 44, 12)])
    assert pitch == 6
    for count in (2, 3):
        xs = horizontal_tracks(100, count, pitch)
        assert [b - a for a, b in zip(xs[:-1], xs[1:], strict=True)] == pytest.approx(
            [6] * (count - 1)
        )


def test_crowding_reduces_the_common_pitch_without_breaking_uniformity():
    pitch = horizontal_pitch(20, [12, 3], [(3, 44, 12)])
    assert 0 < pitch < 6
    for count in (3, 12):
        xs = horizontal_tracks(20, count, pitch)
        assert xs[0] > 0 and xs[-1] == 20
        assert [b - a for a, b in zip(xs[:-1], xs[1:], strict=True)] == pytest.approx(
            [pitch] * (count - 1)
        )
