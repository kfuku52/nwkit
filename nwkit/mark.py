import re
from nwkit.util import *

def annotate_tree_attr(tree, args):
    for node in tree.traverse():
        node.add_props(is_target_leaf=False, is_descendant_all_target=False,
                       is_target_only_mrca=False, is_target_only_mrca_clade=False,
                       is_all_mrca=False, is_all_mrca_clade=False)
    for node in tree.leaves():
        if re.fullmatch(args.pattern, node.name):
            node.props['is_target_leaf'] = True
    for node in tree.traverse(strategy='postorder'):
        if all([l.props.get('is_target_leaf') for l in node.leaves()]):
            node.props['is_descendant_all_target'] = True
    for node in tree.traverse():
        if node.is_root:
            continue
        if (not node.up.props.get('is_descendant_all_target')) and node.props.get('is_descendant_all_target'):
            node.props['is_target_only_mrca'] = True
            node.props['is_target_only_mrca_clade'] = True
            for clade_node in node.descendants():
                clade_node.props['is_target_only_mrca_clade'] = True
    target_leaves = [ leaf for leaf in tree.leaves() if leaf.props.get('is_target_leaf') ]
    num_target_leaves = len(target_leaves)
    if len(target_leaves) > 0:
        for ancestor in target_leaves[0].ancestors():
            num_anc_target_leaves = len([ leaf for leaf in ancestor.leaves() if leaf.props.get('is_target_leaf') ])
            if num_target_leaves == num_anc_target_leaves:
                ancestor.props['is_all_mrca'] = True
                ancestor.props['is_all_mrca_clade'] = True
                for clade_node in ancestor.descendants():
                    clade_node.props['is_all_mrca_clade'] = True
                break
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
