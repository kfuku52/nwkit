import re
import sys
from typing import Any

from nwkit.util import read_tree, write_tree


def _propagate_clade_flag(tree, source_property, clade_property):
    """Mark descendants of flagged nodes, visiting parents before children."""
    clade_flags: dict[Any, bool] = {}
    for node in tree.traverse(strategy="preorder"):
        if node.is_root:
            in_clade = bool(node.props.get(source_property))
        else:
            in_clade = clade_flags[node.up] or bool(node.props.get(source_property))
        clade_flags[node] = in_clade
        if in_clade:
            node.props[clade_property] = True


def annotate_tree_attr(tree, args):
    pattern = re.compile(args.pattern)
    leaf_counts = dict()
    target_leaf_counts = dict()
    for node in tree.traverse():
        node.add_props(
            is_target_leaf=False,
            is_descendant_all_target=False,
            is_target_only_mrca=False,
            is_target_only_mrca_clade=False,
            is_all_mrca=False,
            is_all_mrca_clade=False,
        )
    for node in tree.traverse(strategy="postorder"):
        if node.is_leaf:
            is_target_leaf = pattern.fullmatch(node.name or "") is not None
            node.props["is_target_leaf"] = is_target_leaf
            leaf_counts[node] = 1
            target_leaf_counts[node] = int(is_target_leaf)
        else:
            leaf_count = 0
            target_leaf_count = 0
            for child in node.get_children():
                leaf_count += leaf_counts[child]
                target_leaf_count += target_leaf_counts[child]
            leaf_counts[node] = leaf_count
            target_leaf_counts[node] = target_leaf_count
        if (leaf_counts[node] > 0) and (target_leaf_counts[node] == leaf_counts[node]):
            node.props["is_descendant_all_target"] = True
    for node in tree.traverse():
        if node.is_root:
            node.props["is_target_only_mrca"] = node.props["is_descendant_all_target"]
            continue
        if (not node.up.props.get("is_descendant_all_target")) and node.props.get(
            "is_descendant_all_target"
        ):
            node.props["is_target_only_mrca"] = True
    _propagate_clade_flag(tree, "is_target_only_mrca", "is_target_only_mrca_clade")
    target_leaves = [leaf for leaf in tree.leaves() if leaf.props.get("is_target_leaf")]
    num_target_leaves = len(target_leaves)
    if num_target_leaves > 0:
        if num_target_leaves == 1:
            all_mrca_node = target_leaves[0]
        else:
            all_mrca_node = tree.common_ancestor(target_leaves)
        all_mrca_node.props["is_all_mrca"] = True
    _propagate_clade_flag(tree, "is_all_mrca", "is_all_mrca_clade")
    return tree


def get_insert_nodes(tree, args):
    if args.target == "mrca":
        if args.target_only_clade:
            target_attr = "is_target_only_mrca"
        else:
            target_attr = "is_all_mrca"
    elif args.target == "clade":
        if args.target_only_clade:
            target_attr = "is_target_only_mrca_clade"
        else:
            target_attr = "is_all_mrca_clade"
    elif args.target == "leaf":
        target_attr = "is_target_leaf"
    else:
        raise ValueError("Unknown target: {}".format(args.target))
    insert_nodes = [node for node in tree.traverse() if node.props.get(target_attr)]
    return insert_nodes


def label_insert_nodes(tree, args):
    insert_nodes = get_insert_nodes(tree, args)
    sys.stderr.write("{:,} node(s) will be marked.\n".format(len(insert_nodes)))
    for node in insert_nodes:
        if args.insert_pos == "prefix":
            node.name = args.insert_txt + args.insert_sep + (node.name or "")
        elif args.insert_pos == "suffix":
            node.name = (node.name or "") + args.insert_sep + args.insert_txt
    return tree


def mark_main(args):
    tree = read_tree(
        args.infile,
        args.format,
        args.quoted_node_names,
        rooted=getattr(args, "input_rooted", "auto"),
    )
    tree = annotate_tree_attr(tree, args)
    tree = label_insert_nodes(tree, args)
    outformat = 1 if args.outformat == "auto" else args.outformat
    write_tree(tree, args, format=outformat)
