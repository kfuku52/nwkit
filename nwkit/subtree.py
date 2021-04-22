from nwkit.util import *

def subtree_main(args):
    tree = read_tree(args.infile, args.format)
    leaf_names = tree.get_leaf_names()
    if not args.left_leaf in leaf_names:
        raise Exception('--left_leaf not found in the tree: {}'.format(args.left_leaf))
    if not args.right_leaf in leaf_names:
        raise Exception('--right_leaf not found in the tree: {}'.format(args.right_leaf))
    subtree = tree.get_common_ancestor(args.left_leaf, args.right_leaf)
    write_tree(subtree, args, format=args.format)