from nwkit.util import *

def nhx2nwk_main(args):
    tree = read_tree(args.infile, args.format)
    if (args.node_label!=''):
        for node in tree.traverse():
            if node.is_leaf():
                continue
            if args.node_label in dir(node):
                node.name = getattr(node, args.node_label)
    write_tree(tree, args, format=args.format)
