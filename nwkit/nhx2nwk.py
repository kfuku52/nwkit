from nwkit.util import *

def nhx2nwk_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    if (args.node_label!=''):
        for node in tree.traverse():
            if node.is_leaf:
                continue
            if args.node_label in node.props:
                node.name = node.props.get(args.node_label)
    outformat = args.outformat
    if (outformat == 'auto') and (args.node_label != ''):
        outformat = 1
    write_tree(tree, args, format=outformat)
