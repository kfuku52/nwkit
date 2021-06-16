import sys

from nwkit.util import *

def transfer_root(tree_to, tree_from, verbose=False):
    same_leaf_set = len(set(tree_to.get_leaf_names()) - set(tree_from.get_leaf_names())) == 0
    assert same_leaf_set, 'Input tree and iqtree\'s treefile did not have identical leaves.'
    subroot_leaves = [ n.get_leaf_names() for n in tree_from.get_children() ]
    is_n0_bigger_than_n1 = (len(subroot_leaves[0]) > len(subroot_leaves[1]))
    ingroups = subroot_leaves[0] if is_n0_bigger_than_n1 else subroot_leaves[1]
    outgroups = subroot_leaves[0] if not is_n0_bigger_than_n1 else subroot_leaves[1]
    if verbose:
        sys.stderr.write('outgroups: {}\n'.format(outgroups))
    tree_to.set_outgroup(ingroups[0])
    if (len(outgroups) == 1):
        outgroup_ancestor = [n for n in tree_to.iter_leaves() if n.name == outgroups[0]][0]
    else:
        outgroup_ancestor = tree_to.get_common_ancestor(outgroups)
    tree_to.set_outgroup(outgroup_ancestor)
    subroot_to = tree_to.get_children()
    subroot_from = tree_from.get_children()
    total_subroot_length_to = sum([n.dist for n in subroot_to])
    total_subroot_length_from = sum([n.dist for n in subroot_from])
    for n_to in subroot_to:
        for n_from in subroot_from:
            if (set(n_to.get_leaf_names()) == set(n_from.get_leaf_names())):
                n_to.dist = total_subroot_length_to * (n_from.dist / total_subroot_length_from)
    for n_to in tree_to.traverse():
        if n_to.name == '':
            n_to.name = tree_to.name
            tree_to.name = 'Root'
            break
    return tree_to

def root_main(args):
    tree = read_tree(args.infile, args.format)
    if (args.method=='transfer'):
        tree2 = read_tree(args.infile2, args.format2)
        tree = transfer_root(tree_to=tree, tree_from=tree2, verbose=True)
    write_tree(tree, args, format=args.format)