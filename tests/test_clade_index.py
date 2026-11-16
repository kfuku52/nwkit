import pytest
from ete4 import Tree

from nwkit.clade_index import CladeIndex, LcaIndex


def test_names_for_mask_visits_only_selected_bits_in_name_order():
    index = CladeIndex(Tree("((D,C),(B,A));", parser=1))

    assert index.names_for_mask(0b1010) == ("B", "D")


def _common_ancestor_by_parent_walk(first, second):
    ancestors = set()
    while first is not None:
        ancestors.add(first)
        first = first.up
    while second not in ancestors:
        second = second.up
    return second


def test_names_for_mask_ignores_bits_outside_the_index():
    index = CladeIndex(Tree("(A,B,C);", parser=1))
    assert index.names_for_mask(0b10101) == ("A", "C")
    assert index.names_for_mask(-2) == ("B", "C")


def test_names_for_wide_mask_preserves_order_and_ignores_outside_bits():
    names = tuple(f"T{i:03d}" for i in range(600))
    index = CladeIndex(Tree(f"({','.join(reversed(names))});", parser=1))
    selected = tuple(i for i in range(len(names)) if i % 7 < 3)
    mask = sum(1 << i for i in selected) | (1 << 900)
    assert index.names_for_mask(mask) == tuple(names[i] for i in selected)
    assert index.names_for_mask(-2) == names[1:]


@pytest.mark.parametrize("source", ["A;", "(A,B,C);", "(((A,B),C),(D,(E,F)));"])
def test_lca_matches_parent_walk_for_all_node_pairs(source):
    index = LcaIndex(Tree(source, parser=1))
    for first, first_node in enumerate(index.nodes):
        for second, second_node in enumerate(index.nodes):
            expected = _common_ancestor_by_parent_walk(first_node, second_node)
            assert index.common_ancestor(first_node, second_node) is expected
            assert index.nodes[index.common_ancestor_indices(first, second)] is expected


@pytest.mark.parametrize("difference", [0, 1, 2, 3, 7, 8, 15, 16, 17, 31, 32, 33])
def test_lca_depth_lifting_bit_boundaries(difference):
    tree = Tree()
    branches = []
    for _ in range(2):
        branch = [tree]
        for _ in range(33):
            branch.append(branch[-1].add_child())
        branches.append(branch)
    index = LcaIndex(tree)
    first = branches[0][-1]
    for branch in branches:
        second = branch[33 - difference]
        expected = _common_ancestor_by_parent_walk(first, second)
        for one, two in ((first, second), (second, first)):
            assert index.common_ancestor(one, two) is expected
            assert (
                index.nodes[
                    index.common_ancestor_indices(
                        index.index_by_node[one], index.index_by_node[two]
                    )
                ]
                is expected
            )
