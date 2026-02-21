import random
import ete4

from nwkit.util import *

def get_shuffled_branch_lengths(nodes):
    branch_lengths = [n.dist for n in nodes]
    random.shuffle(branch_lengths)
    return branch_lengths

def print_rf_dist(tree1, tree2):
    leaf_names1 = list(tree1.leaf_names())
    leaf_names2 = list(tree2.leaf_names())
    all_leaf_names = leaf_names1 + leaf_names2
    if any(name is None for name in all_leaf_names):
        sys.stderr.write('Skipping RF distance: leaf names must be non-empty for RF calculation.\n')
        return
    if (len(leaf_names1) != len(set(leaf_names1))) or (len(leaf_names2) != len(set(leaf_names2))):
        sys.stderr.write('Skipping RF distance: duplicated leaf names are not supported for RF calculation.\n')
        return
    use_unrooted_rf = (not is_rooted(tree1)) or (not is_rooted(tree2))
    out = tree1.robinson_foulds(t2=tree2, unrooted_trees=use_unrooted_rf)
    rf, rf_max, common_attrs, edges_t1, edges_t2, discarded_edges_t1, discarded_edges_t2 = out
    sys.stderr.write('Robinson-Foulds distance = {:,} (max = {:,})\n'.format(rf, rf_max))

def shuffle_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    tree_original = tree.copy(method='deepcopy')
    if args.topology:
        num_leaf = len(list(tree.leaves()))
        new_tree = ete4.Tree()
        new_tree.populate(size=num_leaf)
        num_new_branch = sum(1 for _ in new_tree.traverse())
        branch_lengths = get_shuffled_branch_lengths(nodes=tree.traverse())
        population = branch_lengths+branch_lengths # Random tree may have more branches
        new_branch_lengths = random.sample(population=population, k=num_new_branch)
        new_branch_lengths[0] = 0 # Root node
        for n,new_branch_length in zip(new_tree.traverse(), new_branch_lengths):
            n.dist = new_branch_length
        leaf_names = list(tree.leaf_names())
        for leaf, leaf_name in zip(new_tree.leaves(), leaf_names):
            leaf.name = leaf_name
        tree = new_tree
    if args.topology or args.branch_length:
        nodes = [n for n in tree.traverse() if not n.is_root]
        branch_lengths = get_shuffled_branch_lengths(nodes)
        for node,new_bl in zip(nodes,branch_lengths):
            node.dist = new_bl
    if args.label:
        leaf_nodes = list(tree.leaves())
        leaf_names = [leaf.name for leaf in leaf_nodes]
        random.shuffle(leaf_names)
        for leaf,new_name in zip(leaf_nodes, leaf_names):
            leaf.name = new_name
    print_rf_dist(tree1=tree_original, tree2=tree)
    write_tree(tree, args, format=args.outformat)
