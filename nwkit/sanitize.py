from nwkit.util import *

def add_quote(tree, quote_char):
    for node in tree.traverse():
        if not node.name:
            continue
        node.name = quote_char+node.name+quote_char
    return tree

def sanitize_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    if args.remove_singleton:
        tree = remove_singleton(tree, verbose=True)
    if args.resolve_polytomy:
        tree.resolve_polytomy()
    if (args.name_quote=='none'):
        quote_char = ''
    elif (args.name_quote=='single'):
        quote_char = '\''
    elif (args.name_quote=='double'):
        quote_char = '\"'
    else:
        raise ValueError("Unsupported '--name_quote': {}. Choose from none/single/double.".format(args.name_quote))
    tree = add_quote(tree, quote_char)
    write_tree(tree, args, format=args.outformat)
