"""Verify plotted quantities, shared geometry, and real ASR figure exports."""

from types import SimpleNamespace

import numpy as np
import pytest
from ete4 import Tree

from nwkit.asr_figure import build_continuous_asr_figure
from nwkit.asr_regimes import RegimeAssignment
from nwkit.cli import main
from nwkit.continuous_asr import compute_bm_marginals
from nwkit.continuous_asr_io import continuous_output_table
from nwkit.rooting_state import set_rooting_info
from nwkit.util import assign_branch_ids


def summary(tree):
    set_rooting_info(tree, True)
    observed = {"A": 1.0, "B": 3.0}
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
    return table, fit


def test_plot_matches_marginals_and_input_depths():
    tree = Tree("((A:1,B:2):0,C:1):90;")
    table, fit = summary(tree)
    figure = build_continuous_asr_figure(tree, table, model="BM", fit=fit)
    tree_ax, trait_ax = figure.axes
    assert tree_ax.get_ylim() == trait_ax.get_ylim()
    assert tree_ax.get_ylim()[1] < 0 < tree_ax.get_ylim()[0] < 3
    ids = assign_branch_ids(tree)
    means = [
        collection.get_offsets()[0].tolist() for collection in trait_ax.collections
    ]
    intervals = [line for line in trait_ax.lines if line.get_linewidth() == 2.4]
    assert len(intervals) == len(table)
    for node in tree.traverse():
        row = table.set_index("branch_id").loc[ids[node]]
        depth = 0.0
        ancestor = node
        while not ancestor.is_root:
            depth += ancestor.dist
            ancestor = ancestor.up
        assert [row["mean"], depth] in means
        assert any(
            np.allclose(line.get_xdata(), [row.ci_lower, row.ci_upper])
            and np.allclose(line.get_ydata(), [depth, depth])
            for line in intervals
        )
    assert any(
        "Imputed tip" == text.get_text() for text in figure.legends[0].get_texts()
    )
    figure.clear()


def test_regimes_and_theta_use_same_colors_and_values():
    tree = Tree("(A:1,B:2,C:1);")
    table, _ = summary(tree)
    regimes = {node: "high" if node.name == "B" else "low" for node in tree.traverse()}
    assignment = RegimeAssignment(("low", "high"), regimes, "test")
    figure = build_continuous_asr_figure(
        tree,
        table,
        model="OUM",
        regime_assignment=assignment,
        fit=SimpleNamespace(theta_by_regime={"low": 1, "high": 4}),
    )
    # Theta references precede node connectors and intervals.
    low, high = figure.axes[1].lines[:2]
    assert list(low.get_xdata()) == [1, 1]
    assert list(high.get_xdata()) == [4, 4]
    assert list(high.get_ydata()) == [0, 2]
    assert high.get_color() in {line.get_color() for line in figure.axes[0].lines}
    assert low.get_color() != high.get_color()
    figure.clear()


def test_long_tip_labels_clear_the_legends():
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    tree = Tree("(A:1,B:2,C:1);")
    table, fit = summary(tree)
    for node in tree.leaves():
        node.name += "_long_species_and_gene_identifier_0123456789"
    figure = build_continuous_asr_figure(tree, table, model="BM", fit=fit)
    canvas = FigureCanvasAgg(figure)
    canvas.draw()
    renderer = canvas.get_renderer()
    legend_top = max(legend.get_window_extent(renderer).y1 for legend in figure.legends)
    assert all(
        text.get_window_extent(renderer).y0 > legend_top
        for text in figure.axes[0].texts
    )
    assert figure.get_figheight() > 7.5
    figure.clear()


def _command(tmp_path, *extra):
    traits = tmp_path / "traits.tsv"
    if not traits.exists():
        traits.write_text("leaf_name\tvalue\tsecond\nA\t1\t2\nB\t3\t4\nC\tNA\t5\n")
    return [
        "asr",
        "-i",
        "[&R](A:1,B:2,C:1);",
        "--trait",
        str(traits),
        "--state-column",
        "value",
        "-o",
        str(tmp_path / "nodes.tsv"),
        *extra,
    ]


@pytest.mark.parametrize(
    "extension,magic", [("pdf", b"%PDF"), ("png", b"\x89PNG"), ("svg", b"<?xml")]
)
def test_cli_exports_all_formats_without_changing_target(
    tmp_path, extension, magic, monkeypatch
):
    import nwkit.asr_figure as plotting

    captured = []
    original = plotting.build_continuous_asr_figure

    def capture(tree, table, **kwargs):
        captured.append(table)
        return original(tree, table, **kwargs)

    monkeypatch.setattr(plotting, "build_continuous_asr_figure", capture)
    output = tmp_path / f"figure.{extension}"
    assert (
        main(
            _command(
                tmp_path,
                "--sigma2",
                "0.5",
                "--target",
                "missing-leaf",
                "--figure-out",
                str(output),
            )
        )
        is None
    )
    assert output.read_bytes().startswith(magic)
    assert len(captured[0]) == 4
    assert len((tmp_path / "nodes.tsv").read_text().splitlines()) == 2


def test_multivariate_cli_has_one_panel_per_trait(tmp_path, monkeypatch):
    import nwkit.asr_figure as plotting

    original = plotting.build_continuous_asr_figure
    panels = []

    def capture(tree, table, **kwargs):
        fig = original(tree, table, **kwargs)
        panels.extend(ax.get_title(loc="left") for ax in fig.axes)
        return fig

    monkeypatch.setattr(plotting, "build_continuous_asr_figure", capture)
    assert (
        main(
            _command(
                tmp_path,
                "--model",
                "MV-BM",
                "--state-column",
                "value,second",
                "--figure-out",
                str(tmp_path / "mv.svg"),
            )
        )
        is None
    )
    assert panels == ["Phylogeny", "value", "second"]


@pytest.mark.parametrize(
    "options,match",
    [
        (["--figure-out", "bad.jpg"], "extension"),
        (["--figure-width", "2"], "requires --figure-out"),
        (["--figure-out", "bad.pdf", "--figure-height", "0"], "positive and finite"),
        (
            ["--figure-out", "bad.pdf", "--trait-type", "discrete", "--rate", "0.2"],
            "continuous ASR only",
        ),
    ],
)
def test_invalid_options_fail_before_output(tmp_path, options, match):
    with pytest.raises(ValueError, match=match):
        main(_command(tmp_path, *options))
    assert not (tmp_path / "nodes.tsv").exists()


def test_figure_cannot_replace_input_or_other_output(tmp_path):
    command = _command(tmp_path)
    protected = tmp_path / "input.svg"
    protected.write_text((tmp_path / "traits.tsv").read_text())
    before = protected.read_bytes()
    with pytest.raises(ValueError, match="overwrite input"):
        main([*command, "--trait", str(protected), "--figure-out", str(protected)])
    assert protected.read_bytes() == before
    with pytest.raises(ValueError, match="same|distinct|different"):
        main(
            [
                *command,
                "-o",
                str(tmp_path / "out.pdf"),
                "--figure-out",
                str(tmp_path / "out.pdf"),
            ]
        )


def test_failed_render_preserves_existing_figure(tmp_path, monkeypatch):
    from matplotlib.figure import Figure

    path = tmp_path / "figure.png"
    path.write_bytes(b"previous figure")

    def fail(self, path, **kwargs):
        with open(path, "wb") as output:
            output.write(b"partial")
        raise OSError("test render failure")

    monkeypatch.setattr(Figure, "savefig", fail)
    with pytest.raises(OSError, match="test render failure"):
        main(_command(tmp_path, "--figure-out", str(path)))
    assert path.read_bytes() == b"previous figure"


@pytest.mark.parametrize("mode", ["unconditional", "conditional"])
@pytest.mark.parametrize("use_colors", [True, False])
def test_simulation_branch_markers_use_sampled_nodes_and_shared_event_colors(
    mode, use_colors
):
    from nwkit.asr_figure import figure_node_types
    from nwkit.asr_paths import simulate_fitted_paths
    from nwkit.draw_helpers import _get_species_overlap_node_types

    tree = Tree(
        "((Homo_sapiens_1:1,Pan_troglodytes_1:1):1,(Homo_sapiens_2:1,Pan_troglodytes_2:1):1);"
    )
    set_rooting_info(tree, True)
    observed = {name: float(index) for index, name in enumerate(tree.leaf_names())}
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
    args = SimpleNamespace(species_overlap_node_plot="auto")
    types = figure_node_types(tree, args)
    assert types == _get_species_overlap_node_types(tree, args, True)[0]
    assert types[tree] == "duplication"
    assert all(types[child] == "speciation" for child in tree.children)
    simulation = simulate_fitted_paths(
        tree,
        observed,
        None,
        posterior,
        fit=fit,
        model="BM",
        count=3,
        steps=10,
        seed=7,
        mode=mode,
    )
    figure = build_continuous_asr_figure(
        tree,
        table,
        model="BM",
        fit=fit,
        simulation=simulation,
        node_types=types if use_colors else {},
    )
    from matplotlib.collections import PathCollection

    markers = [
        item for item in figure.axes[2].collections if isinstance(item, PathCollection)
    ]
    assert len(markers) == 3  # Original branching nodes only, not grid points/tips.
    for node, item in zip([tree, *tree.children], markers, strict=True):
        expected = (
            simulation.root_values[:, 0]
            if node.is_root
            else simulation.branches[node].values[:, -1, 0]
        )
        np.testing.assert_allclose(item.get_offsets()[:, 0], expected)
        np.testing.assert_allclose(item.get_offsets()[:, 1], 0 if node.is_root else 1)
        expected_color = (
            ([1, 0, 0, 1] if node.is_root else [0, 0, 1, 1])
            if use_colors
            else [0.2, 0.2, 0.2, 1]
        )
        np.testing.assert_allclose(item.get_facecolors(), [expected_color])
    legend = [text.get_text() for text in figure.legends[0].get_texts()]
    assert ("Speciation" in legend) == use_colors
    assert ("Duplication" in legend) == use_colors
    figure.clear()


@pytest.mark.parametrize("mode,expected", [("auto", 0), ("yes", 1), ("no", 0)])
def test_species_event_modes_preserve_partial_label_semantics(mode, expected):
    from nwkit.asr_figure import figure_node_types

    tree = Tree("((Homo_sapiens_1:1,Pan_troglodytes_1:1):1,unparseable:2);")
    types = figure_node_types(tree, SimpleNamespace(species_overlap_node_plot=mode))
    assert len(types) == expected


def test_custom_species_regex_is_shared_with_asrcompare(tmp_path, monkeypatch):
    import nwkit.asr_compare_panels as plotting

    traits = tmp_path / "traits.tsv"
    traits.write_text("leaf_name\tx\na_1\t1\nb_1\t2\na_2\t3\nb_2\t4\n")
    original = plotting.build_comparison_panels
    captured = []

    def build(context, table):
        figure = original(context, table)
        types = context.cache["figure_node_types"]
        assert types[context.tree] == "duplication"
        assert set(types.values()) == {"speciation", "duplication"}
        captured.append(True)
        return figure

    monkeypatch.setattr(plotting, "build_comparison_panels", build)
    main(
        [
            "asrcompare",
            "-i",
            "[&R]((a_1:1,b_1:1):1,(a_2:1,b_2:1):1);",
            "--trait",
            str(traits),
            "--state-column",
            "x",
            "--models",
            "BM",
            "--sigma2",
            "0.5",
            "--figure-layout",
            "panels",
            "--figure-out",
            str(tmp_path / "panels.pdf"),
            "--figure-simulations",
            "1",
            "--figure-simulation-steps",
            "4",
            "--species-regex",
            r"^([^_]+)_",
            "-o",
            str(tmp_path / "comparison.tsv"),
        ]
    )
    assert captured


def test_example_is_binary_and_regimes_match_renamed_trait_tips():
    from pathlib import Path

    import pandas as pd

    from nwkit.asr_figure import figure_node_types
    from nwkit.asr_regimes import read_regime_map
    from nwkit.util import read_tree

    example = Path(__file__).resolve().parents[1] / "examples" / "asr_figure"
    tree = read_tree(str(example / "tree.nwk"), "0", True, quiet=True)
    assert all(len(node.children) == 2 for node in tree.traverse() if not node.is_leaf)
    assert set(pd.read_csv(example / "traits.tsv", sep="\t").leaf_name) == set(
        tree.leaf_names()
    )
    assignment = read_regime_map(str(example / "regimes.tsv"), tree)
    assert len(assignment.by_node) == 31
    types = figure_node_types(tree, SimpleNamespace(species_overlap_node_plot="auto"))
    assert types[tree] == "speciation"
    assert list(types.values()).count("duplication") == 3
    assert list(types.values()).count("speciation") == 12
    for node, event in types.items():
        if event == "duplication":
            assert types[node.up] == "speciation"
            assert all(types[child] == "speciation" for child in node.children)
