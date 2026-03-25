import pandas as pd

from nwkit.consensus import _collect_clade_stats, _read_tree_weights, _scale_support
from nwkit.util import get_subtree_leaf_bitmasks, read_tree, read_trees, support_is_missing


def _mask_to_leaf_set(mask, leaf_names):
    names = list()
    remaining = int(mask)
    index = 0
    while remaining:
        bit = 1 << index
        if remaining & bit:
            names.append(leaf_names[index])
            remaining ^= bit
        index += 1
    return names


def cladefreq_main(args):
    trees = read_trees(args.infile, args.format, args.quoted_node_names)
    if len(trees) == 0:
        raise ValueError('No input trees were found for cladefreq.')
    tree_weights = _read_tree_weights(args.weight_tsv, len(trees))
    total_weight = sum(tree_weights)
    leaf_names, leaf_name_to_bit, _, clade_weights, _ = _collect_clade_stats(
        trees=trees,
        tree_weights=tree_weights,
    )
    reference_mask_to_node = dict()
    if args.reference not in ['', None]:
        reference_tree = read_tree(args.reference, args.reference_format, args.quoted_node_names)
        if set(reference_tree.leaf_names()) != set(leaf_names):
            raise ValueError("Leaf labels in '--reference' must match the input tree collection.")
        subtree_masks = get_subtree_leaf_bitmasks(reference_tree, leaf_name_to_bit)
        for node, mask in subtree_masks.items():
            num_leaves = int(mask).bit_count()
            if node.is_root or node.is_leaf or (num_leaves >= len(leaf_names)):
                continue
            reference_mask_to_node[mask] = node
    rows = list()
    for mask, weight_sum in sorted(
        clade_weights.items(),
        key=lambda item: (-item[1], -int(item[0]).bit_count(), int(item[0])),
    ):
        leaf_set = _mask_to_leaf_set(mask, leaf_names)
        row = {
            'leaf_set': ','.join(leaf_set),
            'num_leaves': len(leaf_set),
            'weight_sum': weight_sum,
            'frequency': _scale_support(weight_sum / total_weight, args.support_scale),
        }
        if args.reference not in ['', None]:
            reference_node = reference_mask_to_node.get(mask)
            row['in_reference'] = reference_node is not None
            if (reference_node is None) or support_is_missing(reference_node.support):
                row['reference_support'] = ''
            else:
                row['reference_support'] = float(reference_node.support)
        rows.append(row)
    out = pd.DataFrame(rows)
    if args.outfile == '-':
        print(out.to_csv(sep='\t', index=False), end='')
    else:
        out.to_csv(args.outfile, sep='\t', index=False)
