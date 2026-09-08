"""Branch-ID labels must agree with ASR input identities after display edits."""

import re

import pandas as pd
import pytest

from nwkit.cli import main


@pytest.mark.parametrize("ladderize", ["yes", "no"])
@pytest.mark.parametrize("layout", ["rectangular", "circular"])
def test_draw_ids_match_table_and_asr(tmp_path, monkeypatch, ladderize, layout):
    import nwkit.draw as drawing

    tree = tmp_path / "tree.nwk"
    # Deliberately unsorted so ladderization changes level-order IDs if done first.
    tree.write_text("[&R](((A:1,B:1):1,C:2):1,D:3);")
    table = tmp_path / "branches.tsv"
    main(["nwk2table", "-i", str(tree), "-o", str(table)])
    expected = pd.read_csv(table, sep="\t").fillna("")
    captured = {}
    original = drawing._draw_tree

    def capture(tree, **kwargs):
        captured.update(
            {
                int(node.props["branch_id"]): str(node.name or "")
                for node in tree.traverse()
            }
        )
        return original(tree, **kwargs)

    monkeypatch.setattr(drawing, "_draw_tree", capture)
    output = tmp_path / "ids.svg"
    main(
        [
            "draw",
            "-i",
            str(tree),
            "-o",
            str(output),
            "--layout",
            layout,
            "--ladderize",
            ladderize,
            "--node-label-property",
            "branch_id",
            "--node-label-target",
            "all",
            "--node-label-prefix",
            "ID=",
            "--support-labels",
            "no",
            "--species-overlap-node-plot",
            "no",
            "--figure-width",
            "7.2",
        ]
    )
    assert captured == dict(zip(expected.branch_id, expected.name, strict=True))
    assert set(re.findall(r">ID=([^<]+)</text>", output.read_text())) == {
        str(i) for i in expected.branch_id
    }
    traits = tmp_path / "traits.tsv"
    traits.write_text("leaf_name\tvalue\nA\t1\nB\t2\nC\t3\nD\t4\n")
    result = tmp_path / "asr.tsv"
    main(
        [
            "asr",
            "-i",
            str(tree),
            "--trait",
            str(traits),
            "--state-column",
            "value",
            "--sigma2",
            "1",
            "-o",
            str(result),
        ]
    )
    asr = pd.read_csv(result, sep="\t").fillna("")
    assert captured == dict(zip(asr.branch_id, asr.name, strict=True))


def test_collapsed_clade_keeps_its_original_id(tmp_path, monkeypatch):
    import nwkit.draw as drawing

    captured = []
    original = drawing._draw_tree

    def capture(tree, **kwargs):
        captured.extend(node.props["branch_id"] for node in tree.traverse())
        return original(tree, **kwargs)

    monkeypatch.setattr(drawing, "_draw_tree", capture)
    output = tmp_path / "collapsed.svg"
    main(
        [
            "draw",
            "-i",
            "[&R](((A:1,B:1):1,C:2):1,D:3);",
            "-o",
            str(output),
            "--node-label-property",
            "branch_id",
            "--node-label-target",
            "all",
            "--node-label-prefix",
            "ID=",
            "--max-visible-tips",
            "2",
            "--collapse-property-aggregation",
            "mean",
            "--ladderize",
            "yes",
            "--species-overlap-node-plot",
            "no",
            "--support-labels",
            "no",
        ]
    )
    assert set(captured) == {0, 1, 2}
    assert set(re.findall(r">ID=([^<]+)</text>", output.read_text())) == {"0", "1", "2"}


def test_input_nhx_ids_are_recomputed_and_filters_use_canonical_ids(tmp_path):
    output = tmp_path / "filtered.svg"
    main(
        [
            "draw",
            "-i",
            "[&R](A:1[&&NHX:branch_id=99],B:1);",
            "-o",
            str(output),
            "--node-label-property",
            "branch_id",
            "--node-label-target",
            "all",
            "--node-label-prefix",
            "ID=",
            "--node-label-filter",
            "branch_id:eq:1",
            "--species-overlap-node-plot",
            "no",
            "--support-labels",
            "no",
        ]
    )
    assert re.findall(r">ID=([^<]+)</text>", output.read_text()) == ["1"]
