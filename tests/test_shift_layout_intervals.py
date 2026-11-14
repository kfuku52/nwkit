"""Check preorder regime assignment against independent parent propagation."""

import itertools

import numpy as np
import pytest

from nwkit.shift_candidates import set_partitions
from nwkit.shift_native_model import ShiftLayout, ShiftTree
from nwkit.util import read_tree


@pytest.mark.parametrize(
    "newick",
    [
        "(((a:1,b:1):1,(c:1,d:1):1):1,((e:1,f:1):1,(g:1,h:1):1):1);",
        "((((a:1,b:1):1,c:2):1,d:3):1,e:4);",
    ],
)
def test_nested_and_shared_regimes_match_parent_walk(newick):
    tree = ShiftTree.build(read_tree(newick, "auto", True, quiet=True))
    for count in range(4):
        for shifts in itertools.combinations(tree.branch_ids[1:], count):
            for groups in set_partitions((0, *shifts)):
                canonical = tuple(sorted(tuple(sorted(group)) for group in groups))
                aliases = {b: j for j, group in enumerate(canonical) for b in group}
                expected = np.zeros(len(tree.branch_ids), dtype=int)
                for i, branch in enumerate(tree.branch_ids):
                    expected[i] = aliases.get(
                        branch, expected[tree.compiled.parents[i]]
                    )
                invalid = any(
                    expected[tree.indices_by_branch[b]]
                    == expected[tree.compiled.parents[tree.indices_by_branch[b]]]
                    for b in shifts
                )
                if invalid:
                    with pytest.raises(ValueError, match="immediate parent's regime"):
                        ShiftLayout.build(tree, shifts, groups)
                else:
                    layout = ShiftLayout.build(tree, shifts, groups)
                    np.testing.assert_array_equal(layout.node_groups(tree), expected)
