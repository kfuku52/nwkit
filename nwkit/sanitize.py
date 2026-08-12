from nwkit.reconciliation_properties import preserve_collapsed_event_boundaries
from nwkit.util import get_tree_property_names, read_tree, remove_singleton, write_tree


def add_quote(tree, quote_char):
    """Deprecated compatibility helper that mutates logical node names."""
    for node in tree.traverse():
        if not node.name:
            continue
        node.name = "{}{}{}".format(quote_char, node.name, quote_char)
    return tree


def sanitize_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    if args.remove_singleton:
        if getattr(args, "preserve_properties", False):
            preserve_collapsed_event_boundaries(tree)
        tree = remove_singleton(tree, verbose=True)
    if args.resolve_polytomy:
        tree.resolve_polytomy()
    if args.name_quote not in ("none", "single", "double"):
        raise ValueError(
            "Unsupported '--name-quote': {}. Choose from none/single/double.".format(
                args.name_quote
            )
        )
    output_properties = (
        get_tree_property_names(tree)
        if getattr(args, "preserve_properties", False)
        else None
    )
    write_tree(
        tree,
        args,
        format=args.outformat,
        name_quote=args.name_quote,
        props=output_properties,
    )
