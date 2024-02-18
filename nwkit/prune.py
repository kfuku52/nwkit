import re
from nwkit.util import *

def prune_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    for node in tree.traverse():
        node.flag_prune = False
    for leaf in tree.iter_leaves():
        if args.invert_match:
            if not re.fullmatch(args.pattern, leaf.name):
                leaf.flag_prune = True
        else:
            if re.fullmatch(args.pattern, leaf.name):
                leaf.flag_prune = True
    for node in tree.traverse(strategy='preorder'):
        if node.is_leaf():
            if node.flag_prune:
                sys.stderr.write(f'Pruning {node.name}\n')
                node.flag_prune = True
        else:
            if all([ leaf.flag_prune for leaf in node.iter_leaves() ]):
                node.flag_prune = True
    for node in tree.traverse(strategy='preorder'):
        if node.flag_prune:
            node.delete(preserve_branch_length=True)
    tree = remove_singleton(tree, verbose=False, preserve_branch_length=True)
    write_tree(tree, args, format=args.outformat)
