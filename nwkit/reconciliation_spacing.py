"""Common horizontal tracks for topological gene joins and species turns."""

from __future__ import annotations

import math


def horizontal_pitch(unit, columns, junctions, spacing=12):
    """Use one pitch per scene, shrinking only for branch or corridor clearance.

    Junctions provide (turn count, corridor width, largest lane offset).
    The pitch is independent of each junction's own turn count.
    """
    pitch = min(spacing / 2, unit / (max(columns, default=0) + 1))
    for count, width, offset in junctions:
        if count > 1:
            radius = width / 2 - 1
            clearance = math.sqrt(max(0, radius * radius - offset * offset))
            pitch = min(pitch, clearance / (count - 1))
    if pitch <= 0:
        raise ValueError(
            "Species corridors have no room for distinct horizontal tracks."
        )
    return pitch


def horizontal_tracks(right, count, pitch):
    """Right-aligned increasing coordinates with a constant pitch."""
    return [right - pitch * (count - 1 - index) for index in range(count)]
