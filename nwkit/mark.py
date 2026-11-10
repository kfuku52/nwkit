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
    first_target = last_target = None
    for node in tree.traverse(strategy="postorder"):
        is_target_leaf = node.is_leaf and pattern.fullmatch(node.name or "") is not None
        if is_target_leaf:
            if first_target is None:
                first_target = node
            last_target = node
        all_target = (
            is_target_leaf
            if node.is_leaf
            else all(child.props["is_descendant_all_target"] for child in node.children)
        )
        node.add_props(
            is_target_leaf=is_target_leaf,
            is_descendant_all_target=all_target,
            is_target_only_mrca=False,
            is_target_only_mrca_clade=False,
            is_all_mrca=False,
            is_all_mrca_clade=False,
        )
    for node in tree.traverse():
        if node.is_root:
            node.props["is_target_only_mrca"] = node.props["is_descendant_all_target"]
            continue
        if (not node.up.props.get("is_descendant_all_target")) and node.props.get(
            "is_descendant_all_target"
        ):
            node.props["is_target_only_mrca"] = True
    _propagate_clade_flag(tree, "is_target_only_mrca", "is_target_only_mrca_clade")
    if first_target is not None and last_target is not None:
        # A clade's leaves form a contiguous span in depth-first order.
        all_mrca_node = (
            first_target
            if first_target is last_target
            else tree.common_ancestor([first_target, last_target])
        )
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
