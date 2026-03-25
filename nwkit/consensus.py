from collections import defaultdict

import pandas as pd
from ete4 import Tree

from nwkit.util import (
    MISSING_SUPPORT_VALUE,
    TREE_FORMAT_PROP,
    count_set_bits,
    get_subtree_leaf_bitmasks,
    is_rooted,
    read_tree,
    read_trees,
    validate_unique_named_leaves,
    write_tree,
)


def _scale_support(freq, support_scale):
    if support_scale == 'percent':
        return freq * 100.0
    if support_scale == 'proportion':
        return freq
    raise ValueError("Unsupported '--support_scale': {}".format(support_scale))


def _bit_to_order_map(leaf_names):
    return {1 << index: index for index in range(len(leaf_names))}


def _iter_mask_bits(mask):
    remaining = int(mask)
    while remaining:
        bit = remaining & -remaining
        yield bit
        remaining ^= bit


def _masks_are_compatible(mask1, mask2):
    overlap = mask1 & mask2
    return (overlap == 0) or (overlap == mask1) or (overlap == mask2)


def _mask_min_order(mask, bit_to_order):
    return min(bit_to_order[bit] for bit in _iter_mask_bits(mask))


def _read_tree_weights(weight_tsv, num_trees):
    if weight_tsv in ['', None]:
        return [1.0] * num_trees
    weight_df = pd.read_csv(weight_tsv, sep='\t')
    if 'weight' not in weight_df.columns:
        raise ValueError("--weight_tsv must contain a 'weight' column.")
    weights = [None] * num_trees
    if 'tree_id' in weight_df.columns:
        if weight_df['tree_id'].isna().any():
            raise ValueError("--weight_tsv contains missing values in 'tree_id'.")
        for _, row in weight_df.iterrows():
            tree_id = int(row['tree_id'])
            if (tree_id < 1) or (tree_id > num_trees):
                raise ValueError("--weight_tsv tree_id is out of range: {}".format(tree_id))
            if weights[tree_id - 1] is not None:
                raise ValueError("Duplicated 'tree_id' values are not supported in --weight_tsv.")
            weights[tree_id - 1] = float(row['weight'])
    else:
        if len(weight_df.index) != num_trees:
            raise ValueError("--weight_tsv must contain exactly one row per input tree.")
        weights = [float(weight) for weight in weight_df['weight'].tolist()]
    if any(weight is None for weight in weights):
        raise ValueError("--weight_tsv must define weights for every input tree.")
    for weight in weights:
        if weight <= 0:
            raise ValueError("Tree weights must be positive.")
    return weights


def _weighted_median(values, weights):
    sorted_pairs = sorted(zip(values, weights), key=lambda pair: pair[0])
    total_weight = sum(weight for _, weight in sorted_pairs)
    cumulative_weight = 0.0
    for value, weight in sorted_pairs:
        cumulative_weight += weight
        if cumulative_weight >= (total_weight / 2.0):
            return value
    return sorted_pairs[-1][0]


def _aggregate_branch_lengths(observations, method):
    if method == 'none':
        return dict()
    out = dict()
    for mask, records in observations.items():
        values = [value for value, _ in records if value is not None]
        weights = [weight for value, weight in records if value is not None]
        if len(values) == 0:
            out[mask] = None
            continue
        if method == 'mean':
            out[mask] = sum(value * weight for value, weight in zip(values, weights)) / sum(weights)
        elif method == 'median':
            out[mask] = _weighted_median(values, weights)
        else:
            raise ValueError("Unsupported '--branch_length': {}".format(method))
    return out


def _collect_clade_stats(trees, tree_weights):
    first_tree = trees[0]
    validate_unique_named_leaves(first_tree, option_name='--infile', context=" for 'consensus'")
    if not is_rooted(first_tree):
        raise ValueError("Consensus currently requires rooted trees. Root the input trees first.")
    leaf_names = list(first_tree.leaf_names())
    leaf_name_to_bit = {leaf_name: index for index, leaf_name in enumerate(leaf_names)}
    all_mask = (1 << len(leaf_names)) - 1
    clade_weights = defaultdict(float)
    branch_length_observations = defaultdict(list)
    for tree_index, (tree, tree_weight) in enumerate(zip(trees, tree_weights), start=1):
        validate_unique_named_leaves(tree, option_name='--infile', context=" for 'consensus'")
        if not is_rooted(tree):
            raise ValueError("Consensus currently requires rooted trees. Tree {} is unrooted.".format(tree_index))
        if set(tree.leaf_names()) != set(leaf_names):
            raise ValueError('Leaf labels must be identical across all input trees for consensus.')
        subtree_masks = get_subtree_leaf_bitmasks(tree, leaf_name_to_bit)
        for node, mask in subtree_masks.items():
            num_tips = count_set_bits(mask)
            if not node.is_root:
                dist_value = None if (node.dist is None) else float(node.dist)
                branch_length_observations[mask].append((dist_value, tree_weight))
            if node.is_root or (num_tips <= 1) or (num_tips >= len(leaf_names)):
                continue
            clade_weights[mask] += tree_weight
    return leaf_names, leaf_name_to_bit, all_mask, dict(clade_weights), dict(branch_length_observations)


def _get_method_min_freq(method, min_freq):
    if method == 'strict':
        return 1.0
    if method == 'majority':
        return 0.5 + 10 ** -12
    if method == 'greedy':
        return min_freq
    raise ValueError("Unsupported '--method': {}".format(method))


def _select_consensus_masks(clade_weights, total_weight, min_freq):
    selected_masks = list()
    candidates = list()
    for mask, weight in clade_weights.items():
        freq = weight / total_weight
        if freq >= min_freq:
            candidates.append((mask, weight, freq))
    candidates.sort(key=lambda item: (-item[1], -count_set_bits(item[0]), int(item[0])))
    for mask, weight, freq in candidates:
        if all(_masks_are_compatible(mask, selected_mask) for selected_mask, _, _ in selected_masks):
            selected_masks.append((mask, weight, freq))
    return selected_masks


def _get_direct_child_masks(parent_mask, selected_masks, bit_to_order):
    candidates = [
        mask for mask in selected_masks
        if (mask != parent_mask) and ((mask & ~parent_mask) == 0)
    ]
    candidates.sort(key=lambda mask: (-count_set_bits(mask), _mask_min_order(mask, bit_to_order)))
    direct_masks = list()
    for mask in candidates:
        if any((mask & ~direct_mask) == 0 for direct_mask in direct_masks):
            continue
        direct_masks.append(mask)
    direct_masks.sort(key=lambda mask: _mask_min_order(mask, bit_to_order))
    return direct_masks


def _build_consensus_subtree(parent_mask, selected_masks, support_by_mask, dist_by_mask, bit_to_name, bit_to_order, all_mask):
    node = Tree()
    node.dist = None if (parent_mask == all_mask) else dist_by_mask.get(parent_mask)
    if parent_mask == all_mask:
        node.support = None
    else:
        node.support = support_by_mask.get(parent_mask, MISSING_SUPPORT_VALUE)
    direct_masks = _get_direct_child_masks(parent_mask, selected_masks, bit_to_order)
    covered_mask = 0
    child_specs = list()
    for child_mask in direct_masks:
        child_specs.append((_mask_min_order(child_mask, bit_to_order), 'mask', child_mask))
        covered_mask |= child_mask
    remainder_mask = parent_mask & ~covered_mask
    for bit in _iter_mask_bits(remainder_mask):
        child_specs.append((bit_to_order[bit], 'leaf', bit))
    child_specs.sort(key=lambda item: item[0])
    for _, child_type, child_value in child_specs:
        if child_type == 'mask':
            child = _build_consensus_subtree(
                parent_mask=child_value,
                selected_masks=selected_masks,
                support_by_mask=support_by_mask,
                dist_by_mask=dist_by_mask,
                bit_to_name=bit_to_name,
                bit_to_order=bit_to_order,
                all_mask=all_mask,
            )
        else:
            child = Tree()
            child.name = bit_to_name[child_value]
            child.dist = dist_by_mask.get(child_value)
            child.support = MISSING_SUPPORT_VALUE
        node.add_child(child=child)
    return node


def _build_consensus_tree(leaf_names, all_mask, selected_masks, support_by_mask, dist_by_mask):
    bit_to_name = {1 << index: leaf_name for index, leaf_name in enumerate(leaf_names)}
    bit_to_order = _bit_to_order_map(leaf_names)
    consensus_tree = _build_consensus_subtree(
        parent_mask=all_mask,
        selected_masks=selected_masks,
        support_by_mask=support_by_mask,
        dist_by_mask=dist_by_mask,
        bit_to_name=bit_to_name,
        bit_to_order=bit_to_order,
        all_mask=all_mask,
    )
    for node in consensus_tree.traverse():
        node.props[TREE_FORMAT_PROP] = 0
    return consensus_tree


def _annotate_reference_tree(reference_tree, leaf_name_to_bit, num_leaves, clade_weights, total_weight, support_scale):
    validate_unique_named_leaves(reference_tree, option_name='--reference', context=" for 'consensus'")
    if not is_rooted(reference_tree):
        raise ValueError("Consensus currently requires rooted reference trees. Root '--reference' first.")
    subtree_masks = get_subtree_leaf_bitmasks(reference_tree, leaf_name_to_bit)
    for node, mask in subtree_masks.items():
        if node.is_root or node.is_leaf or (count_set_bits(mask) >= num_leaves):
            continue
        freq = clade_weights.get(mask, 0.0) / total_weight
        node.support = _scale_support(freq, support_scale)
    reference_tree.props[TREE_FORMAT_PROP] = 0
    return reference_tree


def consensus_main(args):
    method = getattr(args, 'method', 'greedy')
    branch_length = getattr(args, 'branch_length', 'none')
    weight_tsv = getattr(args, 'weight_tsv', None)
    if (args.min_freq < 0.0) or (args.min_freq > 1.0):
        raise ValueError("'--min_freq' must be between 0 and 1.")
    trees = read_trees(args.infile, args.format, args.quoted_node_names)
    if len(trees) == 0:
        raise ValueError('No input trees were found for consensus.')
    tree_weights = _read_tree_weights(weight_tsv, len(trees))
    total_weight = sum(tree_weights)
    leaf_names, leaf_name_to_bit, all_mask, clade_weights, branch_length_observations = _collect_clade_stats(
        trees=trees,
        tree_weights=tree_weights,
    )
    if args.reference not in ['', None]:
        reference_tree = read_tree(args.reference, args.reference_format, args.quoted_node_names)
        if set(reference_tree.leaf_names()) != set(leaf_names):
            raise ValueError("Leaf labels in '--reference' must match the input tree collection.")
        output_tree = _annotate_reference_tree(
            reference_tree=reference_tree,
            leaf_name_to_bit=leaf_name_to_bit,
            num_leaves=len(leaf_names),
            clade_weights=clade_weights,
            total_weight=total_weight,
            support_scale=args.support_scale,
        )
    else:
        selected_masks = _select_consensus_masks(
            clade_weights=clade_weights,
            total_weight=total_weight,
            min_freq=_get_method_min_freq(method, args.min_freq),
        )
        support_by_mask = {
            mask: _scale_support(freq, args.support_scale)
            for mask, _, freq in selected_masks
        }
        dist_by_mask = _aggregate_branch_lengths(
            observations=branch_length_observations,
            method=branch_length,
        )
        output_tree = _build_consensus_tree(
            leaf_names=leaf_names,
            all_mask=all_mask,
            selected_masks=[mask for mask, _, _ in selected_masks],
            support_by_mask=support_by_mask,
            dist_by_mask=dist_by_mask,
        )
    write_tree(output_tree, args, format=0)
