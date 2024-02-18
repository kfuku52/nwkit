import re
from nwkit.util import *

def prune_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    for leaf in tree.iter_leaves():
        flag_prune = False
        if args.invert_match:
            if not re.fullmatch(args.pattern, leaf.name):
                flag_prune = True
        else:
            if re.fullmatch(args.pattern, leaf.name):
                flag_prune = True
        if flag_prune:
            sys.stderr.write(f'Pruning {leaf.name}\n')
            leaf.delete()
    tree = remove_singleton(tree, verbose=False, preserve_branch_length=True)
    write_tree(tree, args, format=args.outformat)
