from nwkit.util import *

def add_quote(tree, quote_char):
    for node in tree.traverse():
        if node.name=='':
            continue
        node.name = quote_char+node.name+quote_char
    return tree

def sanitize_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    if args.remove_singleton:
        tree = remove_singleton(tree, verbose=True)
    if (args.name_quote=='none'):
        quote_char = ''
    if (args.name_quote=='single'):
        quote_char = '\''
    if (args.name_quote=='double'):
        quote_char = '\"'
    tree = add_quote(tree, quote_char)
    write_tree(tree, args, format=args.outformat)