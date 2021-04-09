from nwkit.util import *

def drop_main(args):
    tree = read_tree(args.infile, args.format)
    if (args.target=='all'):
        nodes = list(tree.traverse())
    elif (args.target=='root'):
        nodes = [ node for node in tree.traverse() if node.is_root() ]
    elif (args.target=='leaf'):
        nodes = [ node for node in tree.traverse() if node.is_leaf() ]
    elif (args.target=='intnode'):
        nodes = [ node for node in tree.traverse() if not node.is_leaf() ]
    if args.fill is None:
        placeholder = 123456789
    else:
        placeholder = args.fill
    for node in nodes:
        if (args.name):
            node.name = placeholder
        if (args.support):
            node.support = placeholder
        if (args.length):
            node.dist = placeholder
    write_tree(tree, args, format=args.format)

