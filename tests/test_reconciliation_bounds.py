"""Regressions for clipped labels/corridors and distorted nested corners."""

import copy
import math
from pathlib import Path
from xml.etree import ElementTree as ET

import pytest

from nwkit.reconciliation_bounds import fitted_scene, primitive_bounds
from nwkit.reconciliation_document import read_document, validate_document
from nwkit.reconciliation_html import _svg
from nwkit.reconciliation_nested import _offset
from nwkit.reconciliation_scene import Primitive, Scene, build_scene


@pytest.mark.parametrize("amount", [-12, 0, 5])
def test_nested_offsets_preserve_right_angles_and_duplicate_corner_points(amount):
    path = [(0, 0), (0, 0), (100, 0), (100, 0), (100, 100), (200, 100)]
    points = _offset(path, amount)
    assert all(
        a[0] == b[0] or a[1] == b[1]
        for a, b in zip(points[:-1], points[1:], strict=True)
    )
    assert points[0] == (0, -amount) and points[-1] == (200, 100 - amount)


def test_nested_reverse_and_zero_length_offsets_stay_finite():
    for path in [[(0, 0), (100, 0), (0, 0)], [(0, 0), (0, 0)]]:
        assert all(math.isfinite(v) for p in _offset(path, 5) for v in p)


def test_fitting_wide_corridors_and_multiline_labels_keeps_semantic_coordinates():
    scene = Scene(
        100,
        40,
        [
            Primitive("upper", "wide", [(5, 20), (80, 20)], width=100),
            Primitive("event", "event", [(30, 20)]),
            Primitive("text", "label", [(85, 20)], text="W" * 80),
        ],
    )
    scene.event_points["rec", "event"] = (30, 20)
    before = copy.deepcopy(scene)
    fitted = fitted_scene(scene, 16)
    assert fitted.width > scene.width and fitted.height > scene.height
    assert fitted.event_points["rec", "event"] == fitted.primitives[1].points[0]
    for primitive in fitted.primitives:
        x0, y0, x1, y1 = primitive_bounds(primitive, 16)
        assert 0 <= x0 <= x1 <= fitted.width and 0 <= y0 <= y1 <= fitted.height
    assert scene == before


def test_real_wide_labels_fit_svg_and_html_without_source_changes(
    tmp_path, monkeypatch
):
    import nwkit.reconciliation_render as render

    record = read_document(
        Path(__file__).parents[1] / "examples/reconciliation/toy.json"
    ).record
    for tree in record["trees"][1:]:
        for node in tree["nodes"]:
            if not node["children"]:
                node["label"] += "W" * 70
    doc = validate_document(record)
    before = copy.deepcopy(doc.record)
    scene = build_scene(doc, font_size=16)
    original = render._figure
    checked = []

    def inspect(*args, **kwargs):
        fig = original(*args, **kwargs)
        fig.canvas.draw()
        ax = fig.axes[0]
        renderer = fig.canvas.get_renderer()
        for text in ax.texts:
            box = text.get_window_extent(renderer)
            assert box.x0 >= ax.bbox.x0 - 1 and box.x1 <= ax.bbox.x1 + 1
            assert box.y0 >= ax.bbox.y0 - 1 and box.y1 <= ax.bbox.y1 + 1
        checked.append(True)
        return fig

    monkeypatch.setattr(render, "_figure", inspect)
    render.write_scene(scene, tmp_path / "wide.svg", font_size=16)
    assert checked
    svg = ET.fromstring(_svg(scene, 16))
    full = list(map(float, svg.get("data-full").split()))
    assert full[2] > scene.width
    assert doc.record == before
