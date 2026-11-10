import json

import pytest

from nwkit.dated_presentation import render_dated_tree


@pytest.mark.parametrize("legacy_canvas", [False, True])
def test_dated_presentation_preserves_tip_order_and_annotation_layout(
    tmp_path, monkeypatch, legacy_canvas
):
    if legacy_canvas:
        from matplotlib.figure import Figure

        monkeypatch.delattr(Figure, "draw_without_rendering", raising=False)

    tree = tmp_path / "tree.nwk"
    tree.write_text("((A:8,B:8):8,(C:8,D:8):8);\n")
    tip_order = tmp_path / "tips.tsv"
    tip_order.write_text("species_id\nD\nC\nB\nA\n")
    age_clades = tmp_path / "ages.tsv"
    age_clades.write_text("descendant_species\nA,B\n")
    events = tmp_path / "events.tsv"
    events.write_text("descendant_species\tlabel\tsymbol\nC\tEvent\tD\n")
    plot = tmp_path / "dated.svg"
    layout = tmp_path / "layout.json"

    report = render_dated_tree(
        tree,
        plot,
        tip_order=tip_order,
        age_clades=age_clades,
        branch_annotations=events,
        layout_report=layout,
        figure_width=6,
        figure_height=4,
        font_family="DejaVu Sans",
    )

    assert report == json.loads(layout.read_text())
    assert report["species_order"] == ["D", "C", "B", "A"]
    assert report["mean_age_label_count"] == 1
    assert report["branch_annotations"][0]["display_position_Ma"] == 4
    assert len(report["branch_symbol_bbox_points"]) == 1
    assert "geological-period-" in plot.read_text()


def test_dated_presentation_rejects_interleaved_tip_order(tmp_path):
    tree = tmp_path / "tree.nwk"
    tree.write_text("((A:8,B:8):8,(C:8,D:8):8);\n")
    tip_order = tmp_path / "tips.tsv"
    tip_order.write_text("species_id\nA\nC\nB\nD\n")
    plot = tmp_path / "dated.svg"

    with pytest.raises(ValueError, match="tip order is not compatible"):
        render_dated_tree(tree, plot, tip_order=tip_order)

    assert not plot.exists()


def test_dated_presentation_keeps_repeated_unary_clade_ages(tmp_path):
    tree = tmp_path / "tree.nwk"
    tree.write_text("((A:1,B:1):1);\n")
    age_clades = tmp_path / "ages.tsv"
    age_clades.write_text("descendant_species\nA,B\n")
    events = tmp_path / "events.tsv"
    events.write_text("descendant_species\tlabel\tsymbol\nB,A\tEvent\tD\n")

    report = render_dated_tree(
        tree,
        tmp_path / "dated.svg",
        age_clades=age_clades,
        branch_annotations=events,
    )

    assert report["mean_age_label_count"] == 2
    assert report["branch_annotations"][0]["display_position_Ma"] == 1.5
