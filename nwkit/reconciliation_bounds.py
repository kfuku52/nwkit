"""Conservative display extents shared by vector and offline HTML exports."""

from __future__ import annotations

import copy
import textwrap


def primitive_bounds(primitive, font_size):
    points = list(primitive.points)
    if primitive.kind == "transfer":
        (x, y), (xx, yy) = points[0], points[-1]
        control = ((x + xx) / 2 - (yy - y) / 4, (y + yy) / 2 + (xx - x) / 4)
        points.append(control)  # The quadratic curve stays inside its control hull.
    xs, ys = zip(*points, strict=True)
    if primitive.kind == "text":
        lines = textwrap.wrap(primitive.text, 24) or [""]
        half_height = len(lines) * font_size * 0.65 + 2
        return (
            min(xs),
            min(ys) - half_height,
            max(xs) + max(map(len, lines)) * font_size + 4,
            max(ys) + half_height,
        )
    radius = (
        max(4, primitive.width / 2)
        if primitive.kind == "event"
        else primitive.width / 2
    )
    radius += 1 if primitive.kind == "upper" else 2
    return min(xs) - radius, min(ys) - radius, max(xs) + radius, max(ys) + radius


def fitted_scene(scene, font_size):
    """Fit source geometry without modifying the supplied scene or science."""
    if not scene.primitives:
        return scene
    bounds = [primitive_bounds(p, font_size) for p in scene.primitives]
    dx = max(0, 12 - min(b[0] for b in bounds))
    dy = max(0, 12 - min(b[1] for b in bounds))
    width = max(scene.width, max(b[2] for b in bounds) + dx + 12)
    height = max(scene.height, max(b[3] for b in bounds) + dy + 12)
    if dx == dy == 0 and width == scene.width and height == scene.height:
        return scene
    result = copy.deepcopy(scene)
    for p in result.primitives:
        p.points = [(x + dx, y + dy) for x, y in p.points]
    result.event_points = {
        key: (x + dx, y + dy) for key, (x, y) in result.event_points.items()
    }
    result.width, result.height = width, height
    return result
