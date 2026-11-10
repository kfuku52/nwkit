"""Regression cases for reusable workflow-facing tree APIs."""

import pytest
from ete4 import Tree

from nwkit.astral import extract_support_labels
from nwkit.parent_index import ParentIndex, validate_parent_links
from nwkit.species_parser import SpeciesParser
from nwkit.tree_features import foreground_components, unrooted_splits


def test_parent_index_handles_deep_tree_without_quadratic_ancestor_storage():
    parents = {node: node - 1 for node in range(4096)}
    index = ParentIndex(parents)
    assert index.common_ancestor(4095, 2048) == 2048
    assert len(index.lca.ancestors) == 12
    assert index[3] == [3, 2, 1, 0]


def test_parent_links_reject_disconnected_cycle():
    with pytest.raises(ValueError, match="disconnected"):
        validate_parent_links({0: -1, 1: 2, 2: 1}, {1: [2], 2: [1]}, 0)


@pytest.mark.parametrize(
    "label, expected",
    [
        ("Oryza_sativa_Indica_Group_gene1", "Oryza_sativa_Indica_Group"),
        (
            "Oryza_sativa_tropical_japonica_subgroup_gene1",
            "Oryza_sativa_tropical_japonica_subgroup",
        ),
        (
            "Cenchrus_americanus_x_Cenchrus_purpureus_gene1",
            "Cenchrus_americanus_x_Cenchrus_purpureus",
        ),
    ],
)
def test_extended_taxonomic_is_opt_in(label, expected):
    assert SpeciesParser("extended-taxonomic").parse(label).species_label == expected
    assert SpeciesParser("taxonomic").parse(label).species_label != expected


def test_astral_leaves_and_missing_annotation_are_preserved():
    text = "(('[q1=9]':1,B:1)'[pp1=.95;q1=.7]':1,C:1);"
    result, missing = extract_support_labels(text, "q1")
    assert "'[q1=9]':1" in result
    assert ").7:1" in result
    assert missing == 0
    result, missing = extract_support_labels(text, "absent")
    assert "'[pp1=.95;q1=.7]'" in result
    assert missing == 1


def test_foreground_clades_and_partial_backbone():
    tree = Tree("((A,B),(C,D));", parser=1)
    assert foreground_components(tree, {"A", "B", "D"}) == [("A", "B"), ("D",)]
    assert unrooted_splits(tree, {"A", "B", "C", "D"}) == {frozenset({"A", "B"})}


@pytest.mark.parametrize(
    "newick",
    ["(((B),D),(C,A,E));", "(((B),D),(D,A,B));", "A;"],
)
def test_clade_rank_ids_preserve_mask_order_and_unary_ties(newick):
    from nwkit.branch_features import assign_clade_rank_ids

    tree = Tree(newick, parser=1)
    nodes = list(tree.traverse())
    bits = {name: 1 << i for i, name in enumerate(sorted(tree.leaf_names()))}
    expected = sorted(
        nodes, key=lambda node: sum(bits[name] for name in set(node.leaf_names()))
    )

    assert assign_clade_rank_ids(tree) is tree
    assert sorted(nodes, key=lambda node: node.props["branch_id"]) == expected
    assert [node.props["branch_id"] for node in expected] == list(range(len(nodes)))


def test_strict_support_normalizes_only_verified_edge_representations():
    from nwkit.branch_features import map_internal_support_by_split

    source = Tree("((A,B),(C,D));", parser=1)
    target = Tree("((A,B),(C,D));", parser=1)
    for index, node in enumerate(target.traverse()):
        node.props["branch_id"] = index
    source.children[0].support = 90
    source.children[1].support = 90 + 1e-9
    before = [node.support for node in source.children]
    values, _ = map_internal_support_by_split(target, source, support_max=100)
    assert set(values.values()) == {90}
    assert [node.support for node in source.children] == before
    assert all(node.support is None for node in target.children)


def test_convert_astral_retains_properties_and_can_drop_comments():
    from nwkit.convert import convert_tree_text
    from nwkit.util import read_tree

    text = "((A:1,B:1)'[pp1=.95;q1=.7]':1[note],C:2);"
    result = convert_tree_text(text, node_label="q1", target="nhx")
    node = next(
        node
        for node in read_tree(result, 1, True, quiet=True).traverse()
        if not node.is_root and not node.is_leaf
    )
    assert node.name == ".7" and node.props["pp1"] == ".95"
    plain = convert_tree_text(
        text, node_label="q1", target="newick", properties="drop", comments="drop"
    )
    assert "[" not in plain


def test_convert_keeps_an_astral_label_when_the_requested_key_is_missing():
    from nwkit.convert import convert_tree_text
    from nwkit.util import read_tree

    text = "((A:1,B:1)'[pp1=.95]':1,C:2);"
    result = convert_tree_text(
        text, node_label="q1", tree_format=0, target="newick", properties="drop"
    )
    node = next(
        node
        for node in read_tree(result, 1, True, quiet=True).traverse()
        if not node.is_root and not node.is_leaf
    )
    assert node.name == "[pp1=.95]"


@pytest.mark.parametrize("alias", ["same", "hardlink", "symlink"])
def test_generax_target_does_not_replace_its_input(tmp_path, alias):
    from nwkit.tool_trees import prepare_target

    source, output = tmp_path / "source.nwk", tmp_path / "output.nwk"
    source.write_text("((A:1,B:1)i:1,(C:1,D:1)j:1)root;")
    original = source.read_bytes()
    if alias == "same":
        output = source
    elif alias == "hardlink":
        output.hardlink_to(source)
    else:
        output.symlink_to(source)
    with pytest.raises(ValueError, match="overwrite"):
        prepare_target(source, output)
    assert source.read_bytes() == original


def test_generax_target_publishes_only_after_roundtrip_validation(
    tmp_path, monkeypatch
):
    from nwkit import tool_trees

    source, output = tmp_path / "source.nwk", tmp_path / "output.nwk"
    source.write_text("((A:1,B:1)i:1,(C:1,D:1)j:1)root;")
    output.write_bytes(b"previous result")
    original, calls = tool_trees.ete4.PhyloTree, []

    def fail_roundtrip(*args, **kwargs):
        calls.append(None)
        if len(calls) == 2:
            raise ValueError("roundtrip failure")
        return original(*args, **kwargs)

    monkeypatch.setattr(tool_trees.ete4, "PhyloTree", fail_roundtrip)
    with pytest.raises(ValueError, match="roundtrip failure"):
        tool_trees.prepare_target(source, output)
    assert output.read_bytes() == b"previous result"


@pytest.mark.parametrize("text", ["((A,B),C);", "((A,B),(C,D),E);", "((A,B),(C,D),D);"])
def test_constraint_rejects_missing_extra_and_duplicate_result_taxa(text):
    from nwkit.tree_features import verify_constraint

    with pytest.raises(ValueError):
        verify_constraint(
            Tree(text, parser=1), Tree("((A,B),(C,D));", parser=1), set("ABCD")
        )


def test_constraint_scope_and_partial_backbone_are_explicit():
    from nwkit.tree_features import verify_constraint

    tree = Tree("((A,B),(C,D),E);", parser=1)
    constraint = Tree("((A,B),(C,D));", parser=1)
    assert verify_constraint(tree, constraint, set("ABCDE"), "partial") == {
        frozenset("AB")
    }
    with pytest.raises(ValueError, match="scope"):
        verify_constraint(tree, constraint, set("ABCDE"), "typo")


def test_partial_constraint_projects_extra_tips_on_both_sides_of_a_split():
    from nwkit.tree_features import verify_constraint

    constraint = Tree("((C,D),(B,A));", parser=1)
    matching = Tree("(E,(A,(F,B)),(C,(G,D)));", parser=1)
    conflicting = Tree("(E,(A,(F,C)),(B,(G,D)));", parser=1)
    taxa = set("ABCDEFG")

    assert verify_constraint(matching, constraint, taxa, "partial") == {frozenset("AB")}
    with pytest.raises(ValueError, match="violates 1 required"):
        verify_constraint(conflicting, constraint, taxa, "partial")


@pytest.mark.parametrize(
    "text", ["(A,(B,(C,(D,(E,F)))));", "((A,B),(C,D),(E,F));", "((A),B,((C,D),E),F);"]
)
def test_foreground_components_match_exhaustive_clade_oracle(text):
    from itertools import combinations

    tree = Tree(text, parser=1)
    leaves = sorted(tree.leaf_names())
    clades = {frozenset(node.leaf_names()) for node in tree.traverse()}
    for size in range(len(leaves) + 1):
        for selected in combinations(leaves, size):
            eligible = {clade for clade in clades if clade <= set(selected)}
            maximal = [
                clade
                for clade in eligible
                if not any(clade < other for other in eligible)
            ]
            assert foreground_components(tree, selected) == sorted(
                tuple(sorted(clade)) for clade in maximal
            )


def test_foreground_components_handle_deep_trees_and_reject_bad_labels():
    tree = Tree()
    current = tree
    for index in range(4096):
        current.add_child(name=str(index))
        current = current.add_child()
    current.name = "last"
    names = set(tree.leaf_names())
    assert foreground_components(tree, names) == [tuple(sorted(names))]
    with pytest.raises(ValueError, match="absent"):
        foreground_components(tree, {"unknown"})
    with pytest.raises(ValueError, match="unique"):
        foreground_components(Tree("(A,A);", parser=1), {"A"})


@pytest.mark.parametrize("profile", ["chain", "calibration"])
@pytest.mark.parametrize("alias", ["same", "hardlink", "symlink"])
def test_posterior_diagnostics_require_distinct_chain_files(
    tmp_path, monkeypatch, profile, alias
):
    from nwkit import posterior_diagnostics

    source, second, output = (
        tmp_path / name for name in ("chain1.tsv", "chain2.tsv", "out.tsv")
    )
    source.write_text("Gen\ta\n1\t1\n")
    output.write_bytes(b"previous diagnostics")
    if alias == "same":
        second = source
    elif alias == "hardlink":
        second.hardlink_to(source)
    else:
        second.symlink_to(source)

    def unexpected(*args, **kwargs):
        raise AssertionError("R must not start for repeated chain inputs")

    monkeypatch.setattr(posterior_diagnostics.subprocess, "run", unexpected)
    with pytest.raises(ValueError, match="distinct"):
        posterior_diagnostics.diagnose_chains([source, second], output, profile=profile)
    assert output.read_bytes() == b"previous diagnostics"


@pytest.mark.parametrize("text", ["((A:1,B:1):1,A:2);", "(A:1,:2);"])
def test_distance_matrix_rejects_labels_that_would_drop_or_alias_rows(text):
    from nwkit.tree_features import patristic_distances

    with pytest.raises(ValueError, match="unique nonempty"):
        patristic_distances(Tree(text, parser=1))
