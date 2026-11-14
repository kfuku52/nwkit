"""Regression coverage for the ASR input, rendering and resource review."""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from ete4 import Tree
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.text import Text

from nwkit.asr import _validate_asr_output_paths
from nwkit.asr_paths import simulate_fitted_paths
from nwkit.cli import main


@pytest.mark.parametrize("output", ["outfile", "model_out", "posterior_samples_out"])
@pytest.mark.parametrize("alias", ["same", "symlink", "hardlink"])
def test_species_map_is_preserved_for_all_output_aliases(tmp_path, output, alias):
    source = tmp_path / "species.tsv"
    content = "leaf_name\tspecies_label\nA\tSpecies A\n"
    source.write_text(content)
    target = source if alias == "same" else tmp_path / "alias.tsv"
    if alias != "same":
        try:
            if alias == "symlink":
                target.symlink_to(source)
            else:
                target.hardlink_to(source)
        except OSError:
            pytest.skip("Filesystem does not support this alias")
    args = SimpleNamespace(species_map_tsv=str(source), **{output: str(target)})
    with pytest.raises(ValueError, match="species-map-tsv"):
        _validate_asr_output_paths(args)
    assert source.read_text() == content


@pytest.mark.parametrize("mode", ["conditional", "unconditional"])
def test_vector_grid_memory_rejected_before_refinement(monkeypatch, mode):
    import nwkit.asr_paths as paths

    tree = Tree("(" + ",".join(f"t{i}:1" for i in range(100)) + ");")

    def forbidden(*args):
        pytest.fail("An excessive vector grid must be rejected before allocation")

    monkeypatch.setattr(paths, "_refined_tree", forbidden)
    with pytest.raises(ValueError, match="matrix memory"):
        simulate_fitted_paths(
            tree,
            {},
            None,
            {},
            fit=SimpleNamespace(trait_names=tuple(range(50))),
            model="MV-BM",
            count=1,
            steps=200,
            mode=mode,
        )


@pytest.mark.parametrize("command", ["asr", "asrcompare"])
def test_single_root_tip_heatmap_cli(tmp_path, command):
    tree = tmp_path / "tree.nwk"
    tree.write_text("A:0;")
    traits = tmp_path / "traits.tsv"
    traits.write_text("leaf_name\tx\nA\t2\n")
    figure = tmp_path / "figure.pdf"
    extra = (
        ["--model", "BM"]
        if command == "asr"
        else ["--models", "BM", "--figure-layout", "panels"]
    )
    main(
        [
            command,
            "-i",
            str(tree),
            "--format",
            "1",
            "--input-rooted",
            "yes",
            "--trait",
            str(traits),
            "--state-column",
            "x",
            "--sigma2",
            "0.5",
            "--figure-tip-heatmap",
            "yes",
            "--figure-trait-tip-labels",
            "yes",
            "--figure-out",
            str(figure),
            "-o",
            str(tmp_path / "nodes.tsv"),
            *extra,
        ]
    )
    assert figure.read_bytes().startswith(b"%PDF")


@pytest.mark.parametrize("command", ["asr", "asrcompare"])
def test_many_tip_labels_do_not_overlap_in_equal_width_columns(
    tmp_path, monkeypatch, command
):
    import nwkit.asr_compare_panels as comparison
    import nwkit.asr_figure as single

    module = single if command == "asr" else comparison
    builder = (
        "build_continuous_asr_figure" if command == "asr" else "build_comparison_panels"
    )
    original = getattr(module, builder)
    checked = []

    def inspect(*args, **kwargs):
        figure = original(*args, **kwargs)
        canvas = FigureCanvasAgg(figure)
        canvas.draw()
        groups = []
        for ax in figure.axes:
            if len(ax.texts) == 64:
                groups.append(ax.texts)
            for band in ax.child_axes:
                if band.get_label() == "trait-tip-labels":
                    groups.append([text for text in band.texts if text.get_gid()])
        assert len(groups) == 7
        for texts in groups:
            boxes = [
                Text.get_window_extent(text, canvas.get_renderer()) for text in texts
            ]
            assert len(boxes) == 64
            assert not any(
                a.overlaps(b) for i, a in enumerate(boxes) for b in boxes[i + 1 :]
            )
        checked.append(True)
        return figure

    monkeypatch.setattr(module, builder, inspect)
    names = [f"Species_{i:03}_gene00000000" for i in range(64)]
    tree = tmp_path / "tree.nwk"
    tree.write_text("[&R](" + ",".join(name + ":1" for name in names) + ");")
    traits = tmp_path / "traits.tsv"
    frame = pd.DataFrame(
        np.random.default_rng(42).normal(size=(64, 3)), columns=["x", "y", "z"]
    )
    frame.insert(0, "leaf_name", names)
    frame.to_csv(traits, sep="\t", index=False)
    extra = (
        ["--model", "MV-BM"]
        if command == "asr"
        else ["--models", "MV-BM", "--figure-layout", "panels"]
    )
    main(
        [
            command,
            "-i",
            str(tree),
            "--trait",
            str(traits),
            "--state-column",
            "x,y,z",
            "--figure-tip-heatmap",
            "yes",
            "--figure-trait-tip-labels",
            "yes",
            "--figure-simulations",
            "1",
            "--figure-simulation-steps",
            "2",
            "--seed",
            "7",
            "--figure-out",
            str(tmp_path / "figure.pdf"),
            "-o",
            str(tmp_path / "nodes.tsv"),
            *extra,
        ]
    )
    assert checked == [True]


def test_single_asr_uses_selected_font_during_build_and_export(tmp_path, monkeypatch):
    from matplotlib import rcParams

    import nwkit.asr_compare_figure as fonts
    import nwkit.asr_figure as plotting

    selected = []
    original = plotting.build_continuous_asr_figure
    previous = list(rcParams["font.family"])

    def select(text):
        assert "発現量" in text and "種A" in text
        selected.append(True)
        return "DejaVu Serif"

    def build(*args, **kwargs):
        assert rcParams["font.family"] == ["DejaVu Serif"]
        figure = original(*args, **kwargs)
        original_save = figure.savefig

        def save(*args, **kwargs):
            assert rcParams["font.family"] == ["DejaVu Serif"]
            return original_save(*args, **kwargs)

        monkeypatch.setattr(figure, "savefig", save)
        return figure

    monkeypatch.setattr(fonts, "_font_family_for_text", select)
    monkeypatch.setattr(plotting, "build_continuous_asr_figure", build)
    traits = tmp_path / "traits.tsv"
    traits.write_text("leaf_name\t発現量\n種A\t1\nB\t2\nC\t3\n", encoding="utf-8")
    # Font selection/context is portable; glyph coverage is checked separately
    # when a system CJK font is available.
    import warnings

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Glyph .* missing from font")
        main(
            [
                "asr",
                "-i",
                "[&R](種A:1,B:1,C:1);",
                "--trait",
                str(traits),
                "--state-column",
                "発現量",
                "--model",
                "BM",
                "--sigma2",
                "0.5",
                "--figure-tip-heatmap",
                "yes",
                "--figure-trait-tip-labels",
                "yes",
                "--figure-out",
                str(tmp_path / "figure.pdf"),
                "-o",
                str(tmp_path / "nodes.tsv"),
            ]
        )
    assert selected == [True]
    assert rcParams["font.family"] == previous


@pytest.mark.parametrize("encoding", ["utf-8", "utf-8-sig"])
def test_shared_input_reader_ignores_legacy_os_encoding(
    tmp_path, monkeypatch, encoding
):
    import builtins

    import nwkit.util as util

    source = tmp_path / "traits.tsv"
    content = "leaf_name\t発現量\n種A\t1\n"
    source.write_bytes(content.encode(encoding))
    original_open = builtins.open

    def legacy_open(path, *args, **kwargs):
        kwargs.setdefault("encoding", "cp1252")
        return original_open(path, *args, **kwargs)

    monkeypatch.setattr(util, "open", legacy_open, raising=False)
    assert util.read_input_text(source) == content
    table = util.read_tsv_preserving_leaf_name(source)
    assert table.columns.tolist() == ["leaf_name", "発現量"]
    assert table.iloc[0]["leaf_name"] == "種A"


@pytest.mark.parametrize("encoding", ["utf-8", "utf-8-sig"])
def test_unicode_tree_file_roundtrip_under_legacy_locale(
    tmp_path, monkeypatch, encoding
):
    import builtins

    import nwkit.util as util

    source = tmp_path / "input.nwk"
    output = tmp_path / "output.nwk"
    source.write_bytes("[&R](種A:1,種B:2);".encode(encoding))
    original_open = builtins.open

    def legacy_open(path, *args, **kwargs):
        kwargs.setdefault("encoding", "cp1252")
        return original_open(path, *args, **kwargs)

    monkeypatch.setattr(util, "open", legacy_open, raising=False)
    tree = util.read_tree(source, "1", True, quiet=True)
    util.write_tree(tree, SimpleNamespace(outfile=str(output)), "1", quiet=True)
    assert "種A" in output.read_text(encoding="utf-8")
    for path in (source, output):
        assert len(list(util.iter_tree_strings(path))) == 1
        restored = util.read_trees(path, "1", True, quiet=True)
        assert list(restored[0].leaf_names()) == ["種A", "種B"]
