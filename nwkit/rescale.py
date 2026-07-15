import math

from nwkit.util import get_target_nodes, read_tree, write_tree

def rescale_main(args):
    if not math.isfinite(args.factor):
        raise ValueError("'--factor' must be a finite number.")
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    nodes = get_target_nodes(tree=tree, target=args.target)
    for node in nodes:
        if node.dist is not None:
            node.dist = node.dist * args.factor
    write_tree(tree, args, format=args.outformat)
