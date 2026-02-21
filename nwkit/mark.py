import re
from nwkit.util import *

def annotate_tree_attr(tree, args):
    pattern = re.compile(args.pattern)
    leaf_counts = dict()
    target_leaf_counts = dict()
    for node in tree.traverse():
        node.add_props(is_target_leaf=False, is_descendant_all_target=False,
                       is_target_only_mrca=False, is_target_only_mrca_clade=False,
                       is_all_mrca=False, is_all_mrca_clade=False)
    for node in tree.traverse(strategy='postorder'):
        if node.is_leaf:
            is_target_leaf = (pattern.fullmatch(node.name or '') is not None)
            node.props['is_target_leaf'] = is_target_leaf
            leaf_counts[node] = 1
            target_leaf_counts[node] = int(is_target_leaf)
        else:
            leaf_count = 0
            target_leaf_count = 0
            for child in node.get_children():
                leaf_count += leaf_counts[child]
                target_leaf_count += target_leaf_counts[child]
            leaf_counts[node] = leaf_count
            target_leaf_counts[node] = target_leaf_count
        if (leaf_counts[node] > 0) and (target_leaf_counts[node] == leaf_counts[node]):
            node.props['is_descendant_all_target'] = True
    for node in tree.traverse():
        if node.is_root:
            continue
        if (not node.up.props.get('is_descendant_all_target')) and node.props.get('is_descendant_all_target'):
            node.props['is_target_only_mrca'] = True
    target_only_clade_flags = dict()
    for node in tree.traverse(strategy='preorder'):
        if node.is_root:
            is_target_only_clade = bool(node.props.get('is_target_only_mrca'))
        else:
            is_target_only_clade = target_only_clade_flags[node.up] or bool(node.props.get('is_target_only_mrca'))
        target_only_clade_flags[node] = is_target_only_clade
        if is_target_only_clade:
            node.props['is_target_only_mrca_clade'] = True
    target_leaves = [ leaf for leaf in tree.leaves() if leaf.props.get('is_target_leaf') ]
    num_target_leaves = len(target_leaves)
    if num_target_leaves > 0:
        if num_target_leaves == 1:
            all_mrca_node = target_leaves[0]
        else:
            all_mrca_node = tree.common_ancestor(target_leaves)
        all_mrca_node.props['is_all_mrca'] = True
    all_mrca_clade_flags = dict()
    for node in tree.traverse(strategy='preorder'):
        if node.is_root:
            is_all_mrca_clade = bool(node.props.get('is_all_mrca'))
        else:
            is_all_mrca_clade = all_mrca_clade_flags[node.up] or bool(node.props.get('is_all_mrca'))
        all_mrca_clade_flags[node] = is_all_mrca_clade
        if is_all_mrca_clade:
            node.props['is_all_mrca_clade'] = True
    return tree

def get_insert_nodes(tree, args):
    if args.target == 'mrca':
        if args.target_only_clade:
            target_attr = 'is_target_only_mrca'
        else:
            target_attr = 'is_all_mrca'
    elif args.target == 'clade':
        if args.target_only_clade:
            target_attr = 'is_target_only_mrca_clade'
        else:
            target_attr = 'is_all_mrca_clade'
    elif args.target == 'leaf':
        target_attr = 'is_target_leaf'
    else:
        raise ValueError("Unknown target: {}".format(args.target))
    insert_nodes = [ node for node in tree.traverse() if node.props.get(target_attr) ]
    return insert_nodes

def label_insert_nodes(tree, args):
    insert_nodes = get_insert_nodes(tree, args)
    sys.stderr.write('{:,} node(s) will be marked.\n'.format(len(insert_nodes)))
    for node in insert_nodes:
        if args.insert_pos=='prefix':
            node.name = args.insert_txt + args.insert_sep + (node.name or '')
        elif args.insert_pos=='suffix':
            node.name = (node.name or '') + args.insert_sep + args.insert_txt
    return tree

def mark_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    tree = annotate_tree_attr(tree, args)
    tree = label_insert_nodes(tree, args)
    write_tree(tree, args, format=args.outformat)
