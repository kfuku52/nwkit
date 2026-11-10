"""RF distances from displayed clades/splits, including unresolved roots."""

from nwkit.tree_features import _split_intervals
from nwkit.util import get_subtree_leaf_bitmasks


def _bitmask_splits(tree, rooted, positions):
    count = len(positions)
    all_mask = (1 << count) - 1
    splits = set()
    for node, side in get_subtree_leaf_bitmasks(tree, positions).items():
        if node.is_root:
            continue
        size = side.bit_count()
        if rooted:
            if 1 < size < count:
                splits.add(side)
        elif 1 < size < count - 1:
            splits.add(min(side, all_mask ^ side))
    return splits


def _tree_split_intervals(tree, positions, rooted):
    if len(tree.children) == len(positions):
        # Every child covers exactly one tip, so no split is nontrivial.
        return 0, set()
    return _split_intervals(
        ((node, node.children) for node in tree.traverse(strategy="postorder")),
        positions,
        rooted,
    )


def robinson_foulds(tree1, tree2, *, rooted=True):
    """Return unweighted RF and its sum-of-displayed-splits normalization.

    Polytomies contribute only the clades actually present, not a sampled or
    arbitrary binary resolution.  Root stems, tip edges and unary duplicates
    do not contribute.  Both trees must have the same unique named tips.
    """
    first_tips = list(tree1.leaf_names())
    second_tips = list(tree2.leaf_names())
    positions = {name: index for index, name in enumerate(first_tips)}
    if (
        len(positions) != len(first_tips)
        or len(positions) != len(second_tips)
        or positions.keys() != set(second_tips)
    ):
        # Retain the historical behavior for direct callers with other tip sets.
        first = _bitmask_splits(tree1, rooted, positions)
        second = _bitmask_splits(tree2, rooted, positions)
        return len(first ^ second), len(first) + len(second)

    # Day's DFS ordering makes reference clades contiguous. In unrooted
    # comparisons, use the side excluding the first tip. The complement of a
    # reference clade containing that tip is also contiguous.
    first_count, first = _tree_split_intervals(tree1, positions, rooted)
    second_count, second = _tree_split_intervals(tree2, positions, rooted)
    total = first_count + second_count
    return total - 2 * len(first & second), total
