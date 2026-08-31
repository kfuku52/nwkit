"""RF distances from displayed clades/splits, including unresolved roots."""

from nwkit.util import get_subtree_leaf_name_sets


def _splits(tree, rooted):
    taxa = frozenset(tree.leaf_names())
    splits = set()
    for node, descendants in get_subtree_leaf_name_sets(tree).items():
        if node.is_root:
            continue
        side = frozenset(descendants)
        if rooted:
            if 1 < len(side) < len(taxa):
                splits.add(side)
        elif 1 < len(side) < len(taxa) - 1:
            splits.add(frozenset((side, taxa - side)))
    return splits


def robinson_foulds(tree1, tree2, *, rooted=True):
    """Return unweighted RF and its sum-of-displayed-splits normalization.

    Polytomies contribute only the clades actually present, not a sampled or
    arbitrary binary resolution.  Root stems, tip edges and unary duplicates
    do not contribute.  Both trees must have the same unique named tips.
    """
    first, second = _splits(tree1, rooted), _splits(tree2, rooted)
    return len(first ^ second), len(first) + len(second)
