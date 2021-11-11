import re
from nwkit.util import *

def printlabel_main(args):
    tree = read_tree(args.infile, args.format)
    if (args.target=='all'):
        nodes = list(tree.traverse())
    elif (args.target=='root'):
        nodes = [ node for node in tree.traverse() if node.is_root() ]
    elif (args.target=='leaf'):
        nodes = [ node for node in tree.traverse() if node.is_leaf() ]
    elif (args.target=='intnode'):
        nodes = [ node for node in tree.traverse() if not node.is_leaf() ]
    for node in nodes:
        if re.fullmatch(args.pattern, node.name):
            if args.sister:
                print(' '.join( s.name for s in node.get_sisters() ))
            else:
                print(node.name)