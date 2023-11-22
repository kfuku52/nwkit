from nwkit.util import *

def rescale_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    if (args.target=='all'):
        nodes = list(tree.traverse())
    elif (args.target=='root'):
        nodes = [ node for node in tree.traverse() if node.is_root() ]
    elif (args.target=='leaf'):
        nodes = [ node for node in tree.traverse() if node.is_leaf() ]
    elif (args.target=='intnode'):
        nodes = [ node for node in tree.traverse() if not node.is_leaf() ]
    for node in nodes:
        node.dist = node.dist * args.factor
    write_tree(tree, args, format=args.outformat)

