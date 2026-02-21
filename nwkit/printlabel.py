import re
from nwkit.util import *

def printlabel_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    nodes = get_target_nodes(tree=tree, target=args.target)
    for node in nodes:
        if re.fullmatch(args.pattern, node.name or ''):
            if args.sister:
                sister_names = [s.name for s in node.get_sisters() if s.name]
                print(' '.join(sister_names))
            else:
                print(node.name or '')
