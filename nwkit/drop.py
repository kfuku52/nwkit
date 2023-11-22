from nwkit.util import *

def drop_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    nodes = get_target_nodes(tree=tree, target=args.target)
    if args.fill is None:
        placeholder = -999999
    else:
        placeholder = args.fill
    for node in nodes:
        if (args.name):
            node.name = placeholder
        if (args.support):
            node.support = placeholder
        if (args.length):
            node.dist = placeholder
    write_tree(tree, args, format=args.outformat)

