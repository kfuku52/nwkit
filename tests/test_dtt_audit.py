"""Adversarial DTT numerical, topology and plotting checks."""

import numpy as np
import pandas as pd
import pytest
from ete4 import Tree
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure

from nwkit.cli import parser
from nwkit.disparity import make_design
from nwkit.dtt import _prune_to_crown, _read_data
from nwkit.dtt_figure import _trait_heatmap, _trait_tree
from nwkit.util import assign_branch_ids
from tests.test_dtt import DATA, design_for


def test_deep_crown_copy_is_iterative_and_preserves_input():
    tree = Tree()
    node = tree
    for i in range(300):
        node.add_child(name=f"t{i}", dist=302 - i)
        node = node.add_child(dist=1)
    node.add_child(name="A", dist=2)
    node.add_child(name="B", dist=2)
    names = list(tree.leaf_names())
    identifiers = {
        tip.name: identifier
        for tip, identifier in assign_branch_ids(tree).items()
        if tip.is_leaf
    }
    result = _prune_to_crown(tree, names)
    assert len(list(result.leaves())) == 302
    for tip in result.leaves():
        assert tip.props["_nwkit_dtt_input_branch_id"] == identifiers[tip.name]
    assert all(
        "_nwkit_dtt_input_branch_id" not in node.props for node in tree.traverse()
    )


def test_normalized_event_times_cannot_exceed_present(tmp_path):
    lengths = [
        2.623505221150671,
        2.9919265227070917,
        8.14411514853686,
        0.9282402619296182,
        6.005004254396884,
    ]
    height = sum(lengths)
    inner = "(A:0,B:0)"
    for length in reversed(lengths[1:]):
        inner = f"({inner}:{length!r})"
    text = f"[&R]({inner}:{lengths[0]!r},C:{height!r},D:{height!r});"
    path = tmp_path / "traits.tsv"
    pd.DataFrame({"leaf_name": list("ABCD"), "x": [1, 2, 4, 5]}).to_csv(
        path, sep="\t", index=False
    )
    args = parser.parse_args(
        ["dtt", "-i", text, "--trait", str(path), "--columns", "x", "--n-sim", "0"]
    )
    tree, _, names, _, _, depths, _ = _read_data(args)
    times = make_design(tree, names, depths).times
    assert np.all((times >= 0) & (times <= 1))
    assert np.all(np.diff(times) >= 0)


@pytest.mark.parametrize("layout", ["trees", "heatmap"])
def test_raw_colors_are_invariant_to_extreme_unit_changes(layout):
    tree, _ = design_for()
    colors = []
    for factor in (1.0, 1e-300, 1e300):
        figure = Figure()
        FigureCanvasAgg(figure)
        left, right, color_axis = figure.subplots(1, 3)
        values = DATA[:, :1] * factor
        if layout == "trees":
            _trait_tree(left, color_axis, tree, list("ABCDEFGH"), values[:, 0], "size")
            artist = left.collections[0]
        else:
            _trait_heatmap(
                left, right, color_axis, tree, list("ABCDEFGH"), values, ["size"], "raw"
            )
            artist = right.images[0]
        figure.canvas.draw()
        colors.append(artist.to_rgba(artist.get_array()))
    assert colors[1] == pytest.approx(colors[0], abs=1e-12)
    assert colors[2] == pytest.approx(colors[0], abs=1e-12)


def test_raw_color_map_handles_range_overflow_and_resolves_offset():
    from nwkit.dtt_figure import _original_unit_bar, _raw_color_values

    values = np.array([-1.7e308, 0, 1.7e308])
    normalized, lower, upper = _raw_color_values(values)
    assert normalized == pytest.approx([0, 0.5, 1])
    figure = Figure()
    FigureCanvasAgg(figure)
    axes, color_axis = figure.subplots(1, 2)
    lower = 1e100
    upper = np.nextafter(lower, np.inf)
    artist = axes.imshow([[0, 1]], vmin=0, vmax=1)
    _original_unit_bar(figure, artist, color_axis, lower, upper, "Original units")
    labels = [tick.get_text() for tick in color_axis.get_xticklabels()]
    assert float(labels[0]) == lower
    assert float(labels[-1]) == upper
    assert len(set(labels)) == len(labels)


def test_plot_user_text_is_literal_even_with_tex_enabled(tmp_path, capsys):
    from matplotlib import rc_context

    from nwkit.cli import main
    from tests.test_dtt import command

    with rc_context({"text.usetex": True}):
        # An invalid TeX command must never be interpreted as executable TeX.
        args = command(
            tmp_path, ["--n-sim", "0", "--figure-out", str(tmp_path / "literal.png")]
        )
        table = pd.read_csv(tmp_path / "traits.tsv", sep="\t")
        table = table.rename(columns={"x": r"$\invalidcommand$"})
        table.to_csv(tmp_path / "traits.tsv", sep="\t", index=False)
        args[args.index("--columns") + 1] = r"$\invalidcommand$,y"
        main(args)
    capsys.readouterr()
    assert (tmp_path / "literal.png").read_bytes().startswith(b"\x89PNG")


def test_no_font_failure_mentions_dtt_and_preserves_output(
    tmp_path, monkeypatch, capsys
):
    import nwkit.asr_compare_figure
    from nwkit.cli import main
    from tests.test_dtt import command

    def missing_font(text):
        raise ValueError(
            "The ASR comparison PDF has no single installed font with complete coverage."
        )

    monkeypatch.setattr(nwkit.asr_compare_figure, "_font_family_for_text", missing_font)
    output = tmp_path / "old.png"
    output.write_bytes(b"old")
    with pytest.raises(ValueError, match="DTT figure"):
        main(command(tmp_path, ["--n-sim", "0", "--figure-out", str(output)]))
    assert output.read_bytes() == b"old"
    assert capsys.readouterr().out == ""


def test_frontier_oracle_for_polytomies_zero_edges_and_pruning(tmp_path):
    from scipy.spatial.distance import pdist

    from nwkit.disparity import dtt_curve

    text = "[&R](((A:.25,B:.25):.25,C:.5):.5,((D:.4,E:.4):0,F:.4):.6,G:1);"
    table = pd.DataFrame(
        {
            "leaf_name": list("GFEDCBA"),
            "x": [9, 3, 8, 4, np.nan, 1, 7],
            "y": [1, 8, 5, 6, 2, 3, 4],
        }
    )
    path = tmp_path / "traits.tsv"
    table.to_csv(path, sep="\t", index=False)
    args = parser.parse_args(
        [
            "dtt",
            "-i",
            text,
            "--trait",
            str(path),
            "--columns",
            "x,y",
            "--missing",
            "drop",
            "--n-sim",
            "0",
        ]
    )
    tree, _, names, _, values, depths, _ = _read_data(args)
    design = make_design(tree, names, depths)
    lookup = dict(zip(names, values, strict=True))
    total = pdist(values, metric="sqeuclidean").mean()
    expected, counts = [1.0], [1]
    for time in design.times[1:]:
        frontier = list(tree.children)
        disparities = []
        while frontier:
            node = frontier.pop()
            if node.is_leaf:
                continue
            if depths[node] <= time:
                frontier.extend(node.children)
            else:
                descendants = [lookup[tip.name] for tip in node.leaves()]
                if len(descendants) >= 2:
                    disparities.append(
                        pdist(descendants, metric="sqeuclidean").mean() / total
                    )
        expected.append(np.mean(disparities) if disparities else 0.0)
        counts.append(len(disparities))
    assert dtt_curve(design, values) == pytest.approx(expected)
    assert design.clade_counts.tolist() == counts
    assert names == list("ABDEFG")
