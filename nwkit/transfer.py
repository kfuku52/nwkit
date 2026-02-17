from nwkit.util import *
from nwkit.root import transfer_root

def transfer_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    tree2 = read_tree(args.infile2, args.format2, args.quoted_node_names)
    if is_rooted(tree):
        tree2 = transfer_root(tree_to=tree2, tree_from=tree, verbose=False)
    elif is_rooted(tree2):
        tree = transfer_root(tree_to=tree, tree_from=tree2, verbose=False)
    else:
        tree2 = transfer_root(tree_to=tree2, tree_from=tree, verbose=False)
    if not is_all_leaf_names_identical(tree, tree2, verbose=True):
        raise Exception('Leaf labels in the two trees should be completely matched.')
    target_nodes = get_target_nodes(tree=tree, target=args.target)
    num_target_nodes = len(target_nodes)
    tree2_leaf_index = None
    tree2_clade_index = None
    target_clade_keys = None
    if args.target == 'leaf':
        tree2_leaf_index = {leaf.name: leaf for leaf in tree2.leaves()}
    elif args.target != 'root':
        leaf_name_to_bit = {leaf.name: i for i, leaf in enumerate(tree.leaves())}
        tree2_subtree_leaf_masks = get_subtree_leaf_bitmasks(tree2, leaf_name_to_bit)
        tree2_clade_index = {
            tree2_subtree_leaf_masks[node]: node
            for node in tree2.traverse()
        }
        target_subtree_leaf_masks = get_subtree_leaf_bitmasks(tree, leaf_name_to_bit)
        target_clade_keys = [target_subtree_leaf_masks[node] for node in target_nodes]
    transferred_node_count = 0
    for i, node in enumerate(target_nodes):
        if args.target == 'root':
            tree2_node = tree2
        elif args.target == 'leaf':
            tree2_node = tree2_leaf_index.pop(node.name, None)
        else:
            tree2_node = tree2_clade_index.pop(target_clade_keys[i], None)
        if tree2_node is None:
            txt = 'Skipping. No matching node found in --infile2. Node name in --infile = {}, Leaf names = {}\n'
            sys.stderr.write(txt.format(node.name, ' '.join(list(node.leaf_names()))))
            if args.fill is not None:
                if (args.name):
                    node.name = str(args.fill)
                if (args.support):
                    node.support = float(args.fill)
                if (args.length):
                    node.dist = float(args.fill)
            continue
        if (args.name):
            node.name = tree2_node.name
        if (args.support):
            support = tree2_node.support
            is_missing_support = (support is None)
            if (not is_missing_support):
                is_missing_support = (abs(float(support) + 999999) < 10**-9)
            if is_missing_support:
                if args.fill is not None:
                    node.support = float(args.fill)
            else:
                node.support = support
        if (args.length):
            node.dist = tree2_node.dist
        transferred_node_count += 1
    txt = 'Transferred {} nodes out of {} target nodes.\n'
    sys.stderr.write(txt.format(transferred_node_count, num_target_nodes))
    write_tree(tree, args, format=args.outformat)
