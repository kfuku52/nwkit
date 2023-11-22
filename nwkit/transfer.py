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
    tree2_nodes = list(tree2.traverse())
    transferred_node_count = 0
    for node in target_nodes:
        tree2_node = 'undefined'
        for i in range(len(tree2_nodes)):
            if is_all_leaf_names_identical(node, tree2_nodes[i], verbose=False):
                tree2_node = tree2_nodes.pop(i)
                break
        if tree2_node == 'undefined':
            txt = 'Skipping. No matching node found in --infile2. Node name in --infile = {}, Leaf names = {}\n'
            sys.stderr.write(txt.format(node.name, ' '.join(node.get_leaf_names())))
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
            if (tree2_node.support + 999999 < 10**-9):
                if args.fill is not None:
                    node.support = float(args.fill)
            else:
                node.support = tree2_node.support
        if (args.length):
            node.dist = tree2_node.dist
        transferred_node_count += 1
    txt = 'Transferred {} nodes out of {} target nodes.\n'
    sys.stderr.write(txt.format(transferred_node_count, num_target_nodes))
    write_tree(tree, args, format=args.outformat)
