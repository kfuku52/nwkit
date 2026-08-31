import math
import sys
from itertools import chain

import pandas as pd

from nwkit.consensus import (
    _collect_clade_stats_from_tree_strings,
    _read_tree_weights,
    _scale_support,
)
from nwkit.rooting_state import require_rooted
from nwkit.util import (
    count_set_bits,
    get_subtree_leaf_bitmasks,
    iter_tree_strings,
    read_tree,
    support_is_missing,
    validate_unique_named_leaves,
)

CLADEFREQ_COLUMNS = ("descendant_taxa", "num_taxa", "weight_sum", "frequency")


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
    raw_tree_strings = iter(iter_tree_strings(args.infile))
    try:
        first_tree_string = next(raw_tree_strings)
    except StopIteration:
        raise ValueError("No input trees were found for cladefreq.") from None
    tree_count = [0]

    def counted_tree_strings():
        for tree_string in chain((first_tree_string,), raw_tree_strings):
            tree_count[0] += 1
            yield tree_string

    tree_weights = _read_tree_weights(args.weight_tsv)
    leaf_names, leaf_name_to_bit, _, clade_weights, _ = (
        _collect_clade_stats_from_tree_strings(
            tree_strings=counted_tree_strings(),
            tree_weights=tree_weights,
            format=args.format,
            quoted_node_names=args.quoted_node_names,
            collect_branch_lengths=False,
            threads=getattr(args, "threads", 1),
            require_rooted=True,
            input_rooted=getattr(args, "input_rooted", "auto"),
        )
    )
    num_trees = tree_count[0]
    if tree_weights is not None and len(tree_weights) != num_trees:
        raise ValueError("--weight-tsv must contain exactly one row per input tree.")
    sys.stderr.write("Number of input trees = {:,}\n".format(num_trees))
    try:
        total_weight = (
            float(num_trees) if tree_weights is None else math.fsum(tree_weights)
        )
    except OverflowError as exc:
        raise ValueError(
            "The sum of tree weights is too large for cladefreq output."
        ) from exc
    if not math.isfinite(total_weight):
        raise ValueError("The sum of tree weights must be finite.")
    reference_mask_to_node = dict()
    if args.reference not in ["", None]:
        reference_tree = read_tree(
            args.reference,
            args.reference_format,
            args.quoted_node_names,
            rooted=getattr(args, "reference_rooted", "auto"),
        )
        validate_unique_named_leaves(
            reference_tree,
            option_name="--reference",
            context=" for 'cladefreq'",
        )
        require_rooted(
            reference_tree,
            "'--reference' must be rooted for rooted clade-frequency comparisons.",
            "--reference-rooted",
        )
        if set(reference_tree.leaf_names()) != set(leaf_names):
            raise ValueError(
                "Leaf labels in '--reference' must match the input tree collection."
            )
        subtree_masks = get_subtree_leaf_bitmasks(reference_tree, leaf_name_to_bit)
        for node, mask in subtree_masks.items():
            num_leaves = count_set_bits(mask)
            if node.is_root or node.is_leaf or (num_leaves >= len(leaf_names)):
                continue
            reference_mask_to_node[mask] = node
    rows = list()
    for mask, weight_sum in sorted(
        clade_weights.items(),
        key=lambda item: (-item[1], -count_set_bits(item[0]), int(item[0])),
    ):
        leaf_set = _mask_to_leaf_set(mask, leaf_names)
        row = {
            "descendant_taxa": ",".join(leaf_set),
            "num_taxa": len(leaf_set),
            "weight_sum": weight_sum,
            "frequency": _scale_support(weight_sum / total_weight, args.support_scale),
        }
        if args.reference not in ["", None]:
            reference_node = reference_mask_to_node.get(mask)
            row["in_reference"] = reference_node is not None
            if (reference_node is None) or support_is_missing(reference_node.support):
                row["reference_support"] = ""
            else:
                row["reference_support"] = float(reference_node.support)
        rows.append(row)
    columns = list(CLADEFREQ_COLUMNS)
    if args.reference not in ["", None]:
        columns.extend(("in_reference", "reference_support"))
    out = pd.DataFrame(rows, columns=columns)
    if args.outfile == "-":
        print(out.to_csv(sep="\t", index=False), end="")
    else:
        out.to_csv(args.outfile, sep="\t", index=False)
