"""Heatmaps show observed tips in tree order, with explicit missingness and scales."""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from ete4 import Tree
from matplotlib import colormaps
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.colors import to_rgba

from nwkit.asr_figure import build_continuous_asr_figure
from nwkit.asr_heatmap import _trait_scale, tip_trait_values
from nwkit.cli import main
from nwkit.continuous_asr import compute_bm_marginals
from nwkit.continuous_asr_io import continuous_output_table
from nwkit.rooting_state import set_rooting_info


def summary():
    # Internal name duplicates a tip; non-ultrametric input depths stay intact.
    tree = Tree("((A:1,B:2)A:1,C:1);", parser=1)
    set_rooting_info(tree, True)
    observed = {"A": 1.0, "B": 5.0}
    posterior, fit = compute_bm_marginals(tree, observed, sigma2=0.5)
    table = continuous_output_table(
        tree,
        list(tree.traverse()),
        observed,
        None,
        posterior,
        trait="value",
        ci_level=0.95,
    )
    return tree, table, fit


def test_heatmap_observations_order_colors_and_layout():
    tree, table, fit = summary()
    traits, values = tip_trait_values(table, list(tree.leaves()))
    assert traits == ["value"]
    np.testing.assert_allclose(values, [[1, 5, np.nan]], equal_nan=True)
    figure = build_continuous_asr_figure(
        tree, table, model="BM", fit=fit, tip_heatmap=True
    )
    canvas = FigureCanvasAgg(figure)
    canvas.draw()
    heatmap, key = figure.axes[0].child_axes
    assert heatmap.get_label() == "observed-tip-heatmap"
    cells = heatmap.patches
    assert [cell.get_x() + 0.5 for cell in cells] == [0, 1, 2]
    np.testing.assert_allclose(cells[0].get_facecolor(), colormaps["viridis"](0.0))
    np.testing.assert_allclose(cells[1].get_facecolor(), colormaps["viridis"](1.0))
    np.testing.assert_allclose(cells[2].get_facecolor(), to_rgba("#D0D0D0"))
    assert key.get_xticks().tolist() == [1, 5]
    assert figure.axes[0].get_ylim() == figure.axes[1].get_ylim()
    renderer = canvas.get_renderer()
    labels = [text.get_window_extent(renderer) for text in figure.axes[0].texts]
    assert max(box.y1 for box in labels) < heatmap.get_window_extent(renderer).y0
    assert min(box.y0 for box in labels) > key.get_tightbbox(renderer).y1
    figure.clear()


def test_multiple_traits_get_separate_rows_and_scales():
    tree, table, fit = summary()
    second = table.copy()
    second["trait"] = "second"
    second.loc[second.observed_value != "", "observed_value"] = [10, 20]
    combined = pd.concat([table, second], ignore_index=True)
    figure = build_continuous_asr_figure(
        tree, combined, model="BM", fit=fit, tip_heatmap=True
    )
    FigureCanvasAgg(figure).draw()
    heatmap, first_key, second_key = figure.axes[0].child_axes
    assert len(heatmap.patches) == 6
    assert first_key.get_xticks().tolist() == [1, 5]
    assert second_key.get_xticks().tolist() == [10, 20]
    assert first_key.get_title(loc="left").startswith("1.")
    assert second_key.get_title(loc="left").startswith("2.")
    figure.clear()


@pytest.mark.parametrize("values, ticks", [([2, 2], [2]), ([np.nan, np.nan], [])])
def test_constant_and_missing_scales_are_defined(values, ticks):
    norm, actual = _trait_scale(np.array(values))
    assert actual == ticks
    assert norm.vmin < norm.vmax


def test_heatmap_option_requires_figure_and_panel_layout():
    from nwkit.asr_compare_panels import validate_panel_options
    from nwkit.asr_figure import validate_figure_options

    args = SimpleNamespace(figure_tip_heatmap="yes")
    with pytest.raises(ValueError, match="requires --figure-out"):
        validate_figure_options(args)
    with pytest.raises(ValueError, match="require --figure-layout panels"):
        validate_panel_options(args, "continuous")


def test_comparison_heatmaps_have_identical_observations_and_scales(
    tmp_path, monkeypatch
):
    import nwkit.asr_compare_panels as plotting

    source = tmp_path / "traits.tsv"
    source.write_text("leaf_name\tx\nA\t1\nB\t3\nC\tNA\nD\t4\n")
    original = plotting.build_comparison_panels
    captured = []

    def build(context, table):
        figure = original(context, table)
        trees = [ax for ax in figure.axes if ax.child_axes]
        assert len(trees) == 2
        first, second = [ax.child_axes for ax in trees]
        np.testing.assert_allclose(
            [p.get_facecolor() for p in first[0].patches],
            [p.get_facecolor() for p in second[0].patches],
        )
        np.testing.assert_allclose(first[1].get_xticks(), second[1].get_xticks())
        captured.append(True)
        return figure

    monkeypatch.setattr(plotting, "build_comparison_panels", build)
    main(
        [
            "asrcompare",
            "-i",
            "[&R]((A:1,B:1):1,(C:1,D:1):1);",
            "--trait",
            str(source),
            "--state-column",
            "x",
            "--models",
            "BM,OU",
            "--alpha",
            "0.5",
            "--sigma2",
            "1",
            "--figure-layout",
            "panels",
            "--figure-tip-heatmap",
            "yes",
            "--figure-out",
            str(tmp_path / "heatmap.pdf"),
            "-o",
            str(tmp_path / "comparison.tsv"),
        ]
    )
    assert captured
