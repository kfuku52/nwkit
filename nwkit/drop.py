from nwkit.util import *

def drop_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    nodes = get_target_nodes(tree=tree, target=args.target)
    if args.fill is None:
        placeholder = -999999
    else:
        placeholder = args.fill
    numeric_placeholder = None
    if args.fill is not None and (args.support or args.length):
        try:
            numeric_placeholder = float(args.fill)
        except ValueError as exc:
            raise ValueError("'--fill' must be numeric when '--support' or '--length' is enabled.") from exc
    for node in nodes:
        if (args.name):
            node.name = placeholder
        if (args.support):
            node.support = numeric_placeholder if (numeric_placeholder is not None) else placeholder
        if (args.length):
            node.dist = numeric_placeholder if (numeric_placeholder is not None) else placeholder
    write_tree(tree, args, format=args.outformat)
