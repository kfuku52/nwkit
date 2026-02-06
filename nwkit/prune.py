import re
from nwkit.util import *

def prune_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    for node in tree.traverse():
        node.add_props(flag_prune=False)
    for leaf in tree.leaves():
        if args.invert_match:
            if not re.fullmatch(args.pattern, leaf.name):
                leaf.props['flag_prune'] = True
        else:
            if re.fullmatch(args.pattern, leaf.name):
                leaf.props['flag_prune'] = True
    for node in tree.traverse(strategy='preorder'):
        if node.is_leaf:
            if node.props.get('flag_prune'):
                sys.stderr.write(f'Pruning {node.name}\n')
                node.props['flag_prune'] = True
        else:
            if all([ leaf.props.get('flag_prune') for leaf in node.leaves() ]):
                node.props['flag_prune'] = True
    for node in tree.traverse(strategy='preorder'):
        if node.props.get('flag_prune'):
            node.delete(preserve_branch_length=True)
    tree = remove_singleton(tree, verbose=False, preserve_branch_length=True)
    write_tree(tree, args, format=args.outformat)
