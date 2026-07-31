from nwkit.util import read_tree, remove_singleton, write_tree

def add_quote(tree, quote_char):
    """Deprecated compatibility helper that mutates logical node names."""
    for node in tree.traverse():
        if not node.name:
            continue
        node.name = '{}{}{}'.format(quote_char, node.name, quote_char)
    return tree

def sanitize_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    if args.remove_singleton:
        tree = remove_singleton(tree, verbose=True)
    if args.resolve_polytomy:
        tree.resolve_polytomy()
    if args.name_quote not in ('none', 'single', 'double'):
        raise ValueError("Unsupported '--name-quote': {}. Choose from none/single/double.".format(args.name_quote))
    write_tree(tree, args, format=args.outformat, name_quote=args.name_quote)
