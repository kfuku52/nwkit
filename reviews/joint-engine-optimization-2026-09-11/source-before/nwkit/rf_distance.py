"""RF distances from displayed clades/splits, including unresolved roots."""

from nwkit.util import get_subtree_leaf_bitmasks


def _splits(tree, rooted, positions):
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


def robinson_foulds(tree1, tree2, *, rooted=True):
    """Return unweighted RF and its sum-of-displayed-splits normalization.

    Polytomies contribute only the clades actually present, not a sampled or
    arbitrary binary resolution.  Root stems, tip edges and unary duplicates
    do not contribute.  Both trees must have the same unique named tips.
    """
    positions = {name: index for index, name in enumerate(tree1.leaf_names())}
    first, second = _splits(tree1, rooted, positions), _splits(tree2, rooted, positions)
    return len(first ^ second), len(first) + len(second)
