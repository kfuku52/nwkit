from nwkit.util import *
from nwkit.root import transfer_root

def transfer_main(args):
    if args.infile2 in ['', None]:
        raise ValueError("'--infile2' is required for 'transfer'.")
    numeric_fill = None
    if (args.fill is not None) and (args.support or args.length):
        try:
            numeric_fill = float(args.fill)
        except ValueError as exc:
            raise ValueError("'--fill' must be numeric when '--support' or '--length' is enabled.") from exc
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    tree2 = read_tree(args.infile2, args.format2, args.quoted_node_names)
    validate_unique_named_leaves(tree, option_name='--infile', context=" for 'transfer'")
    validate_unique_named_leaves(tree2, option_name='--infile2', context=" for 'transfer'")
    num_leaves_tree = len(list(tree.leaves()))
    num_leaves_tree2 = len(list(tree2.leaves()))
    def transfer_root_or_keep(tree_to, tree_from):
        tree_to_candidate = tree_to.copy(method='deepcopy')
        try:
            return transfer_root(tree_to=tree_to_candidate, tree_from=tree_from, verbose=False)
        except (SystemExit, ValueError):
            sys.stderr.write('Skipping root transfer because no matching root bipartition was found.\n')
            return tree_to
    def has_bifurcating_root(tree_obj):
        return len(tree_obj.get_children()) == 2
    if (args.target != 'leaf') and (num_leaves_tree > 1) and (num_leaves_tree2 > 1):
        if has_bifurcating_root(tree):
            tree2 = transfer_root_or_keep(tree_to=tree2, tree_from=tree)
        elif has_bifurcating_root(tree2):
            tree = transfer_root_or_keep(tree_to=tree, tree_from=tree2)
        else:
            sys.stderr.write('Skipping root transfer because neither input tree has a bifurcating root.\n')
    if not is_all_leaf_names_identical(tree, tree2, verbose=True):
        raise Exception('Leaf labels in the two trees should be completely matched.')
    target_nodes = get_target_nodes(tree=tree, target=args.target)
    num_target_nodes = len(target_nodes)
    tree2_leaf_index = None
    tree2_intnode_index = None
    target_intnode_clade_keys = None
    if args.target in ('leaf', 'all'):
        tree2_leaf_index = {leaf.name: leaf for leaf in tree2.leaves()}
    if args.target in ('intnode', 'all'):
        leaf_name_to_bit = {leaf.name: i for i, leaf in enumerate(tree.leaves())}
        tree2_subtree_leaf_masks = get_subtree_leaf_bitmasks(tree2, leaf_name_to_bit)
        tree2_intnode_index = dict()
        for node in tree2.traverse():
            if node.is_leaf:
                continue
            key = tree2_subtree_leaf_masks[node]
            if key not in tree2_intnode_index:
                tree2_intnode_index[key] = []
            tree2_intnode_index[key].append(node)
        target_subtree_leaf_masks = get_subtree_leaf_bitmasks(tree, leaf_name_to_bit)
        target_intnode_clade_keys = {
            node: target_subtree_leaf_masks[node]
            for node in target_nodes
            if not node.is_leaf
        }
    transferred_node_count = 0
    for node in target_nodes:
        if args.target == 'root':
            tree2_node = tree2
        elif node.is_leaf:
            tree2_node = tree2_leaf_index.pop(node.name, None)
        else:
            matching_nodes = tree2_intnode_index.get(target_intnode_clade_keys[node], [])
            tree2_node = matching_nodes.pop(0) if (len(matching_nodes) > 0) else None
        if tree2_node is None:
            txt = 'Skipping. No matching node found in --infile2. Node name in --infile = {}, Leaf names = {}\n'
            sys.stderr.write(txt.format(node.name, ' '.join(list(node.leaf_names()))))
            if args.fill is not None:
                if (args.name):
                    node.name = str(args.fill)
                if (args.support):
                    node.support = numeric_fill
                if (args.length):
                    node.dist = numeric_fill
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
                    node.support = numeric_fill
            else:
                node.support = support
        if (args.length):
            node.dist = tree2_node.dist
        transferred_node_count += 1
    txt = 'Transferred {} nodes out of {} target nodes.\n'
    sys.stderr.write(txt.format(transferred_node_count, num_target_nodes))
    outformat = args.outformat
    if (outformat == 'auto') and args.name:
        outformat = 1
    write_tree(tree, args, format=outformat)
