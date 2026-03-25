import sys

from nwkit.util import (
    read_tree,
    remove_singleton,
    support_is_missing,
    write_tree,
)


def _get_support_value(node):
    support = getattr(node, 'support', None)
    if support_is_missing(support):
        return None
    return float(support)


def _should_collapse(node, args):
    if node.is_root or node.is_leaf:
        return False
    if (args.min_support is not None):
        support = _get_support_value(node)
        if (support is not None) and (support < args.min_support):
            return True
    if (args.max_dist is not None) and (node.dist is not None) and (float(node.dist) <= args.max_dist):
        return True
    return False


def _preserve_deleted_branch_length(node):
    if node.dist is None:
        return
    deleted_dist = float(node.dist)
    for child in node.get_children():
        child_dist = 0.0 if (child.dist is None) else float(child.dist)
        child.dist = child_dist + deleted_dist


def collapse_main(args):
    if (args.min_support is None) and (args.max_dist is None):
        raise ValueError("Specify at least one of '--min_support' or '--max_dist'.")
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    if args.min_support is not None:
        support_values = [
            _get_support_value(node)
            for node in tree.traverse()
            if (not node.is_root) and (not node.is_leaf)
        ]
        support_values = [support for support in support_values if support is not None]
        if len(support_values) == 0:
            raise ValueError("No meaningful support values were found for '--min_support'.")
    collapsed_count = 0
    for node in list(tree.traverse(strategy='postorder')):
        if not _should_collapse(node, args):
            continue
        if args.preserve_branch_length:
            _preserve_deleted_branch_length(node)
        node.delete(
            prevent_nondicotomic=False,
            preserve_branch_length=False,
        )
        collapsed_count += 1
    remove_singleton(tree, preserve_branch_length=args.preserve_branch_length)
    sys.stderr.write('Collapsed {} internal node(s).\n'.format(collapsed_count))
    write_tree(tree, args, format=args.outformat)
