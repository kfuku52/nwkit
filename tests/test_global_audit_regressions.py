"""Boundary and complexity regressions from the September whole-program audit."""

import itertools

import pandas as pd
import pytest

from nwkit.cli import main
from nwkit.mark import annotate_tree_attr, get_insert_nodes
from nwkit.rf_distance import robinson_foulds
from nwkit.skim import add_group_ids, mark_traits_to_nodes
from nwkit.util import _LeafIntervalSet, read_tree
from tests.helpers import make_args, make_deep_ladder_tree


@pytest.mark.parametrize("force", ["yes", "no"])
def test_label_reserves_existing_non_target_names(tmp_path, force):
    output = tmp_path / "tree.nwk"
    main(
        [
            "label",
            "--infile",
            "(n0:1,:1)n1;",
            "--format",
            "1",
            "--target",
            "leaf",
            "--force",
            force,
            "--outfile",
            str(output),
        ]
    )
    tree = read_tree(str(output), "1", True, quiet=True)
    names = [node.name for node in tree.traverse()]
    assert len(names) == len(set(names))
    assert tree.name == "n1"
    assert set(tree.leaf_names()) == {"n0", "n2"}


def test_skim_mixed_text_is_an_ordinary_category():
    tree = read_tree("((A:1,B:1):1,C:1)R;", "1", True, quiet=True)
    df = pd.DataFrame(
        {"leaf_name": ["A", "B", "C"], "trait": ["_MIXED_", "_MIXED_", "other"]}
    )
    mark_traits_to_nodes(tree, df, make_args(group_by="trait"))
    groups = add_group_ids(df, tree).set_index("leaf_name")["group"]
    assert groups.notna().all()
    assert groups["A"] == groups["B"] != groups["C"]


@pytest.mark.parametrize("target,expected", [("clade", 3), ("mrca", 1)])
def test_mark_all_matches_includes_root_clade(target, expected):
    tree = read_tree("(A:1,B:1)R;", "1", True, quiet=True)
    args = make_args(pattern=".*", target=target, target_only_clade=True)
    annotate_tree_attr(tree, args)
    nodes = get_insert_nodes(tree, args)
    assert len(nodes) == expected
    assert tree in nodes


def test_consensus_deep_tree_is_iterative(tmp_path):
    source = tmp_path / "deep.nwk"
    source.write_text(make_deep_ladder_tree(1100).write(parser=1))
    output = tmp_path / "consensus.nwk"
    main(
        [
            "consensus",
            "--infile",
            str(source),
            "--format",
            "1",
            "--method",
            "strict",
            "--outfile",
            str(output),
        ]
    )
    tree = read_tree(str(output), "0", True, quiet=True)
    original = read_tree(str(source), "1", True, quiet=True)
    assert robinson_foulds(tree, original)[0] == 0
    assert len(list(tree.leaves())) == 1100


def test_leaf_interval_does_not_consume_prefix():
    class IndexOnlyTuple(tuple):
        def __iter__(self):
            raise AssertionError("interval iteration scanned the prefix")

    names = IndexOnlyTuple(str(i) for i in range(10000))
    view = _LeafIntervalSet(names, {"9999": 9999}, 9999, 10000)
    assert list(view) == ["9999"]


@pytest.mark.parametrize("rooted", [True, False])
def test_bitmask_rf_matches_set_oracle(rooted):
    sources = ["((A,B),(C,D))R;", "((A,C),(B,D))R;", "(A,B,C,D)R;", "(D,(C,(B,A)))R;"]
    trees = [read_tree(source, "1", True, quiet=True) for source in sources]

    def oracle(tree):
        taxa = frozenset(tree.leaf_names())
        splits = set()
        for node in tree.traverse():
            if node.is_root:
                continue
            side = frozenset(node.leaf_names())
            if rooted and 1 < len(side) < len(taxa):
                splits.add(side)
            elif not rooted and 1 < len(side) < len(taxa) - 1:
                splits.add(frozenset((side, taxa - side)))
        return splits

    for first, second in itertools.product(trees, repeat=2):
        a, b = oracle(first), oracle(second)
        assert robinson_foulds(first, second, rooted=rooted) == (
            len(a ^ b),
            len(a) + len(b),
        )
