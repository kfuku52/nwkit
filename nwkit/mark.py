import re
from nwkit.util import *

def initialize_tree_attr(tree):
    for node in tree.traverse():
        node.is_target = False
        node.is_descendant_all_target = False
        node.is_mrca = False
    return tree

def get_target_nodes(tree, args):
    if (args.target=='all'):
        nodes = list(tree.traverse())
    elif (args.target=='root'):
        nodes = [ node for node in tree.traverse() if node.is_root() ]
    elif (args.target=='leaf'):
        nodes = [ node for node in tree.traverse() if node.is_leaf() ]
    elif (args.target=='intnode'):
        nodes = [ node for node in tree.traverse() if not node.is_leaf() ]
    return nodes

def get_insert_nodes(tree, args):
    nodes = get_target_nodes(tree, args)
    insert_nodes = list()
    if (args.mrca & ~args.target_only_clade):
        insert_nodes = node[0].get_common_ancestor(nodes[1:len(nodes)])
    elif (args.mrca & args.target_only_clade):
        for node in nodes:
            if re.fullmatch(args.pattern, node.name):
                node.is_target = True
        for node in tree.traverse(strategy='postorder'):
            if all([ l.is_target for l in node.iter_leaves() ]):
                node.is_descendant_all_target = True
        for node in tree.traverse():
            if node.is_root():
                continue
            if (~node.up.is_descendant_all_target & node.is_descendant_all_target):
                insert_nodes.append(node)
    else:
        for node in nodes:
            if re.fullmatch(args.pattern, node.name):
                insert_nodes.append(node)
    return insert_nodes

def label_insert_nodes(tree, args):
    insert_nodes = get_insert_nodes(tree, args)
    sys.stderr.write('{:,} node(s) will be marked.\n'.format(len(insert_nodes)))
    for node in insert_nodes:
        if args.insert_pos=='prefix':
            node.name = args.insert_txt + args.insert_sep + node.name
        elif args.insert_pos=='suffix':
            node.name = node.name + args.insert_sep + args.insert_txt
    return tree

def mark_main(args):
    tree = read_tree(args.infile, args.format)
    tree = initialize_tree_attr(tree)
    tree = label_insert_nodes(tree, args)
    write_tree(tree, args, format=args.format)
