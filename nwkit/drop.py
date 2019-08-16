from nwkit.util import *

def drop_main(args):
    tree = read_tree(args)
    if (args.target=='all'):
        nodes = list(tree.traverse())
    elif (args.target=='root'):
        nodes = [ node for node in tree.traverse() if node.is_root() ]
    elif (args.target=='leaf'):
        nodes = [ node for node in tree.traverse() if node.is_leaf() ]
    elif (args.target=='intnode'):
        nodes = [ node for node in tree.traverse() if not node.is_leaf() ]
    for node in nodes:
        if (args.name=='yes'):
            node.name = 123456789 # placeholder
        if (args.support=='yes'):
            node.support = 123456789 # placeholder
        if (args.length=='yes'):
            node.dist = 123456789 # placeholder
    write_tree(tree, args)

