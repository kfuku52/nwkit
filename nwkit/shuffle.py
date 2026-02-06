import random
import ete4

from nwkit.util import *

def get_shuffled_branch_lengths(nodes):
    branch_lengths = [n.dist for n in nodes]
    random.shuffle(branch_lengths)
    return branch_lengths

def print_rf_dist(tree1, tree2):
    out = tree1.robinson_foulds(t2=tree2, unrooted_trees=False)
    rf, rf_max, common_attrs, edges_t1, edges_t2, discarded_edges_t1, discarded_edges_t2 = out
    sys.stderr.write('Robinson-Foulds distance = {:,} (max = {:,})\n'.format(rf, rf_max))

def shuffle_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    tree_original = tree
    if args.topology:
        num_leaf = len(list(tree.leaves()))
        new_tree = ete4.Tree()
        new_tree.populate(size=num_leaf)
        num_new_branch = len([ n for n in new_tree.traverse() ])
        branch_lengths = get_shuffled_branch_lengths(nodes=tree)
        population = branch_lengths+branch_lengths # Random tree may have more branches
        new_branch_lengths = random.sample(population=population, k=num_new_branch)
        new_branch_lengths[0] = 0 # Root node
        for n,new_branch_length in zip(new_tree.traverse(), new_branch_lengths):
            n.dist = new_branch_length
        leaf_names = list(tree.leaf_names())
        for leaf, leaf_name in zip(new_tree.leaves(), leaf_names):
            leaf.name = leaf_name
        tree = new_tree
    if args.topology|args.branch_length:
        nodes = [n for n in tree.traverse() if not n.is_root]
        branch_lengths = get_shuffled_branch_lengths(nodes)
        for node,new_bl in zip(nodes,branch_lengths):
            node.dist = new_bl
    if args.label:
        leaf_names = list(tree.leaf_names())
        random.shuffle(leaf_names)
        for leaf,new_name in zip(tree.leaves(), leaf_names):
            leaf.name = new_name
    print_rf_dist(tree1=tree_original, tree2=tree)
    write_tree(tree, args, format=args.outformat)
