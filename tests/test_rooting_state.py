import io
import json
import sys

import pandas as pd
import pytest
from ete4 import Tree

from nwkit.cli import main
from nwkit.consensus import _collect_clade_stats_from_tree_strings
from nwkit.diff import compare_trees
from nwkit.root import midpoint_rooting
from nwkit.rooting_state import ROOTING_PROPERTIES, get_rooting_info, is_rooted
from nwkit.transfer import _validate_property_name
from nwkit.util import (
    copy_tree_iteratively,
    iter_tree_strings,
    read_tree,
    read_trees,
    write_tree,
)
from tests.helpers import make_args


def read(text, mode="auto"):
    return read_tree(text, "auto", True, quiet=True, rooted=mode)


@pytest.mark.parametrize(
    "text,mode,state,source",
    [
        ("(A:1,B:1);", "auto", "rooted", "topology"),
        ("A;", "auto", "rooted", "topology"),
        ("(A:1,B:1,C:1);", "auto", "unknown", "unknown"),
        ("(A:1,B:1,C:1)Root;", "auto", "unknown", "unknown"),
        ("(A:1,B:1,C:1)0;", "auto", "unknown", "unknown"),
        ("(A:1,B:1,C:1)1;", "auto", "unknown", "unknown"),
        ("(A:1,B:1,C:1)100;", "auto", "unknown", "unknown"),
        ("(A:1,B:1,C:1):10;", "auto", "unknown", "unknown"),
        ("[&R] (A:1,B:1,C:1);", "auto", "rooted", "marker"),
        ("[&u] (A:1,(B:1,C:1):1);", "auto", "unrooted", "marker"),
        ("[&U] (A:1,B:1,C:1)Root;", "yes", "rooted", "override"),
        ("[&R] (A:1,B:1);", "no", "unrooted", "override"),
        ("(A:1,B:1,C:1);", "yes", "rooted", "override"),
        ("(A:1,B:1,C:1)[&&NHX:nwkit_rooted=yes];", "auto", "rooted", "nhx"),
        ("(A:1,B:1)[&&NHX:nwkit_rooted=no];", "auto", "unrooted", "nhx"),
    ],
)
def test_input_rooting_precedence(text, mode, state, source):
    tree = read(text, mode)
    info = get_rooting_info(tree)
    assert (info.state, info.source) == (state, source)
    assert is_rooted(tree) is (state == "rooted")
    # Root interpretation is not changed by decorative labels or support.
    tree.name = "new name"
    tree.support = 0
    assert get_rooting_info(tree) == info


def test_forced_conflict_is_reported(capsys):
    tree = read("[&U](A:1,B:1,C:1);", "yes")
    assert "overrides its explicit unrooted" in capsys.readouterr().err
    assert get_rooting_info(tree).declared == "no"


@pytest.mark.parametrize(
    "text",
    [
        "[&R][&U](A,B);",
        "[&U](A,B)[&&NHX:nwkit_rooted=yes];",
        "(A,B)[&&NHX:nwkit_rooted=maybe];",
        "((A,B)[&&NHX:nwkit_rooted=yes],C);",
    ],
)
def test_invalid_rooting_declarations_are_not_bypassed(text):
    with pytest.raises(ValueError):
        read(text, "yes")


def test_tokens_in_quoted_labels_are_not_declarations():
    tree = read("('A[&U]':1,B:1,C:1)'[&R]';")
    assert get_rooting_info(tree).state == "unknown"
    assert "A[&U]" in tree.leaf_names()


@pytest.mark.parametrize("source_kind", ["inline", "file", "stdin"])
def test_nexus_stream_retains_each_rooting_declaration(
    tmp_path, monkeypatch, source_kind
):
    text = "#NEXUS\nBEGIN TREES;\nTREE r = [&R](A:1,B:1,C:1);\nUTREE u = [&U](A:1,(B:1,C:1):1);\nEND;\n"
    if source_kind == "file":
        path = tmp_path / "trees.nex"
        path.write_text(text)
        source = str(path)
    elif source_kind == "stdin":
        monkeypatch.setattr(sys, "stdin", io.StringIO(text))
        source = "-"
    else:
        source = text
    records = list(iter_tree_strings(source))
    assert [get_rooting_info(read(record)).state for record in records] == [
        "rooted",
        "unrooted",
    ]
    assert [is_rooted(tree) for tree in read_trees(text, "auto", True, quiet=True)] == [
        True,
        False,
    ]


@pytest.mark.parametrize("text", ["[&R](A:1,B:1,C:1);", "[&U](A:1,B:1);"])
def test_copy_and_newick_roundtrip_preserve_rooting_without_custom_properties(text):
    tree = copy_tree_iteratively(read(text))
    target = io.StringIO()
    write_tree(tree, make_args(outfile=target), "auto", quiet=True, props=["unrelated"])
    output = target.getvalue()
    assert "nwkit_rooted=" in output
    assert "_nwkit_rooting" not in output
    # Still readable with the existing ETE Newick parser.
    Tree(output, parser=0)
    assert is_rooted(read(output)) is is_rooted(tree)


def test_subtree_write_does_not_attach_tree_metadata_to_original_descendant():
    tree = read("[&R]((A:1,B:1,C:1):1,D:2);")
    subtree = tree.common_ancestor(["A", "B"])
    target = io.StringIO()
    write_tree(subtree, make_args(outfile=target), "auto", quiet=True)
    assert is_rooted(read(target.getvalue()))
    assert not (ROOTING_PROPERTIES & subtree.props.keys())


def test_collapse_and_prune_keep_rooted_polytomy(tmp_path):
    output = tmp_path / "tree.nwk"
    main(["collapse", "-i", "((A:1,B:1):0,C:1);", "--max-dist", "0", "-o", str(output)])
    tree = read(str(output))
    assert is_rooted(tree) and len(tree.children) == 3
    main(["prune", "-i", "((A:1,B:1,C:1):1,D:2);", "--pattern", "D", "-o", str(output)])
    tree = read(str(output))
    assert is_rooted(tree) and len(tree.children) == 3


def test_rooting_operation_replaces_explicit_unrooted_state():
    tree = midpoint_rooting(read("[&U](A:1,B:1,C:2);"))
    assert is_rooted(tree)
    assert not any(
        ROOTING_PROPERTIES & node.props.keys()
        for node in tree.traverse()
        if not node.is_root
    )


@pytest.mark.parametrize("state", ["yes", "no", "unknown"])
def test_tree_table_roundtrip_retains_rooting(tmp_path, state):
    table, output = tmp_path / "nodes.tsv", tmp_path / "tree.nwk"
    main(
        [
            "nwk2table",
            "-i",
            "(A:1,B:1,C:1)[&&NHX:nwkit_rooted=" + state + "];",
            "-o",
            str(table),
        ]
    )
    frame = pd.read_csv(table, sep="\t", keep_default_na=False)
    assert frame.loc[frame.parent == -1, "rooted"].tolist() == [state]
    assert set(frame.loc[frame.parent != -1, "rooted"]) == {""}
    main(["table2nwk", "-i", str(table), "-o", str(output)])
    assert read(str(output)).props["nwkit_rooted"] == state


@pytest.mark.parametrize(
    "root_state,child_state", [("maybe", ""), ("yes", "no"), ("", "yes")]
)
def test_table_rejects_invalid_or_nonroot_declarations(
    tmp_path, root_state, child_state
):
    table = tmp_path / "nodes.tsv"
    table.write_text(
        "branch_id\tparent\tname\trooted\n"
        f"0\t-1\tR\t{root_state}\n1\t0\tA\t{child_state}\n2\t0\tB\t\n"
    )
    with pytest.raises(ValueError, match="nwkit_rooted"):
        main(["table2nwk", "-i", str(table)])


@pytest.mark.parametrize("prop", sorted(ROOTING_PROPERTIES))
def test_annotations_cannot_reassign_rooting_metadata(prop):
    with pytest.raises(ValueError, match="reserved"):
        _validate_property_name(prop)


def test_two_inputs_need_independent_overrides(tmp_path):
    output = tmp_path / "rf.tsv"
    command = [
        "dist",
        "-i",
        "(A:1,B:1,C:1);",
        "--infile2",
        "(A:1,B:1,C:1);",
        "--metric",
        "rf",
        "--input-rooted",
        "yes",
        "-o",
        str(output),
    ]
    with pytest.raises(ValueError, match="--infile2-rooted"):
        main(command)
    main(command + ["--infile2-rooted", "yes"])
    assert pd.read_csv(output, sep="\t").distance.tolist() == [0]


def test_rf_uses_displayed_clades_without_resolving_polytomies(tmp_path):
    output = tmp_path / "rf.tsv"
    main(
        [
            "dist",
            "-i",
            "[&R]((A:1,B:1):1,C:2,D:2);",
            "--infile2",
            "[&R]((A:1,C:1):1,B:2,D:2);",
            "--metric",
            "rf,normalized-rf",
            "-o",
            str(output),
        ]
    )
    frame = pd.read_csv(output, sep="\t")
    assert frame.distance.tolist() == [2, 1]
    assert frame.max_distance.tolist() == [2, 1]


def test_validate_compares_all_root_children_order_independently(tmp_path):
    trees = tmp_path / "trees.nwk"
    trees.write_text(
        "[&R]((A,B),(C,D),E);\n[&R](E,(D,C),(B,A));\n[&R]((A,C),(B,D),E);\n"
    )
    output = tmp_path / "validate.tsv"
    main(
        [
            "validate",
            "-i",
            str(trees),
            "--require-rooted",
            "yes",
            "--require-same-rooting",
            "yes",
            "-o",
            str(output),
        ]
    )
    frame = pd.read_csv(output, sep="\t")
    assert frame.rooting_matches_first.tolist() == [True, True, False]
    assert frame.rooting_state.tolist() == ["rooted"] * 3
    assert frame.rooting_source.tolist() == ["marker"] * 3


def test_validate_distinguishes_unknown_from_explicitly_unrooted(tmp_path):
    trees, output = tmp_path / "trees.nwk", tmp_path / "validation.tsv"
    trees.write_text("(A,B,C);\n[&U](A,B,C);\n")
    main(
        [
            "validate",
            "-i",
            str(trees),
            "--require-same-rooting",
            "yes",
            "-o",
            str(output),
        ]
    )
    frame = pd.read_csv(output, sep="\t")
    assert frame.rooting_state.tolist() == ["unknown", "unrooted"]
    assert frame.rooting_matches_first.tolist() == [True, False]


def test_diff_compares_every_declared_root_child():
    target = read("[&R]((A,B),(C,D),E);")
    for newick, status in (
        ("[&R](E,(B,A),(D,C));", "same"),
        ("[&R]((A,C),(B,D),E);", "different"),
    ):
        rows = compare_trees(target, read(newick))
        root_row = next(row for row in rows if row["comparison"] == "root_split")
        assert root_row["status"] == status
        assert root_row["reason"] == "projected_root_partition"
        assert "target=E|A,B|C,D;" in root_row["shared_taxa_or_split"]


def test_diff_does_not_call_an_explicitly_unrooted_binary_root_a_match():
    rows = compare_trees(read("[&R]((A,B),C);"), read("[&U]((A,B),C);"))
    root_row = next(row for row in rows if row["comparison"] == "root_split")
    assert (root_row["status"], root_row["reason"]) == (
        "unresolved",
        "input_explicitly_unrooted",
    )


@pytest.mark.parametrize("threads", [1, 2])
def test_collection_override_reaches_every_worker(threads):
    trees = ["[&R]((A:1,B:1):1,C:2,D:2);"] + ["((A:1,B:1):1,C:2,D:2);"] * 64
    options = dict(
        tree_strings=trees,
        tree_weights=None,
        format="auto",
        quoted_node_names=True,
        collect_branch_lengths=False,
        require_rooted=True,
        threads=threads,
    )
    with pytest.raises(ValueError, match=r"Input tree \d+ is not rooted"):
        _collect_clade_stats_from_tree_strings(**options)
    leaf_names, _, _, weights, _ = _collect_clade_stats_from_tree_strings(
        **options, input_rooted="yes"
    )
    assert leaf_names == ["A", "B", "C", "D"]
    assert list(weights.values()) == [65]


@pytest.mark.parametrize("command", ["consensus", "cladefreq"])
def test_collection_reference_has_an_independent_override(tmp_path, command):
    output = tmp_path / "result"
    argv = [
        command,
        "-i",
        "[&R]((A:1,B:1):1,C:2,D:2);",
        "--reference",
        "[&U]((A:1,B:1):1,C:2,D:2);",
        "--input-rooted",
        "yes",
        "-o",
        str(output),
    ]
    with pytest.raises(ValueError, match="--reference-rooted"):
        main(argv)
    main(argv + ["--reference-rooted", "yes"])
    assert output.read_text()


@pytest.mark.parametrize("comparison", ["rooted", "unrooted"])
def test_consensus_output_records_its_selected_interpretation(tmp_path, comparison):
    output = tmp_path / "consensus.nwk"
    main(
        [
            "consensus",
            "-i",
            "(A:1,B:1,C:1);",
            "--comparison",
            comparison,
            "-o",
            str(output),
        ]
    )
    tree = read(str(output))
    assert len(tree.children) == 3
    assert get_rooting_info(tree).state == comparison


def test_consensus_cannot_silently_override_an_explicit_unrooted_input(tmp_path):
    output = tmp_path / "consensus.nwk"
    argv = ["consensus", "-i", "[&U](A:1,B:1,C:1);", "-o", str(output)]
    with pytest.raises(ValueError, match="--input-rooted"):
        main(argv)
    main(argv + ["--input-rooted", "yes"])
    assert is_rooted(read(str(output)))


def test_audit_records_declared_and_forced_states(tmp_path):
    audit, output = tmp_path / "audit.jsonl", tmp_path / "validate.tsv"
    main(
        [
            "validate",
            "-i",
            "[&U](A:1,B:1,C:1);",
            "--input-rooted",
            "yes",
            "--audit",
            str(audit),
            "-o",
            str(output),
        ]
    )
    record = json.loads(audit.read_text().strip())
    # The summary is independently read from the original input, with the same mode.
    summary = record["primary_input"]
    assert summary["first_tree_rooting_state"] == "rooted"
    assert summary["first_tree_rooting_source"] == "override"
    assert summary["first_tree_rooting_declared"] == "no"


def test_forcing_rooted_does_not_bypass_binary_requirements(tmp_path):
    trait = tmp_path / "traits.tsv"
    trait.write_text("leaf_name\tx\nA\t1\nB\t2\nC\t3\n")
    with pytest.raises(ValueError, match="strictly bifurcating"):
        main(
            [
                "contrast",
                "-i",
                "(A:1,B:1,C:1);",
                "--input-rooted",
                "yes",
                "--trait",
                str(trait),
                "--columns",
                "x",
            ]
        )
