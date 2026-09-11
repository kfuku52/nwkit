import re
import sys

from nwkit.reconciliation_properties import (
    preserve_collapsed_event_boundaries,
    preserve_pruning_event_boundaries,
)
from nwkit.util import get_tree_property_names, read_tree, remove_singleton, write_tree


def _output_properties(tree, args):
    if not getattr(args, "preserve_properties", False):
        return None
    return get_tree_property_names(tree)


def prune_main(args):
    tree = read_tree(
        args.infile,
        args.format,
        args.quoted_node_names,
        rooted=getattr(args, "input_rooted", "auto"),
    )
    compiled_pattern = re.compile(args.pattern)
    prune_flags = dict()
    leaf_nodes = list(tree.leaves())
    pruned_leaf_names = list()
    for leaf in leaf_nodes:
        is_match = compiled_pattern.fullmatch(leaf.name or "") is not None
        if args.invert_match:
            flag_prune = not is_match
        else:
            flag_prune = is_match
        prune_flags[leaf] = flag_prune
        if flag_prune:
            pruned_leaf_names.append(leaf.name)
    if pruned_leaf_names:
        sys.stderr.write("".join(f"Pruning {name}\n" for name in pruned_leaf_names))
    if len(pruned_leaf_names) == len(leaf_nodes):
        raise ValueError(
            "All leaves would be pruned. Update --pattern/--invert-match to retain at least one leaf."
        )
    if len(pruned_leaf_names) == 0:
        write_tree(
            tree,
            args,
            format=args.outformat,
            props=_output_properties(tree, args),
        )
        return
    for node in tree.traverse(strategy="postorder"):
        if node.is_leaf:
            continue
        prune_flags[node] = all(
            prune_flags.get(child, False) for child in node.get_children()
        )
    if getattr(args, "preserve_properties", False):
        preserve_pruning_event_boundaries(tree, prune_flags)
    for node in list(tree.traverse(strategy="postorder")):
        if node.is_root:
            continue
        if prune_flags.get(node, False):
            node.delete(preserve_branch_length=True)
    if getattr(args, "preserve_properties", False):
        preserve_collapsed_event_boundaries(tree)
    tree = remove_singleton(tree, verbose=False, preserve_branch_length=True)
    write_tree(
        tree,
        args,
        format=args.outformat,
        props=_output_properties(tree, args),
    )
