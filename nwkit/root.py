import sys

from nwkit.util import *

def transfer_root(tree_to, tree_from, verbose=False):
    subroot_leaves = [ n.get_leaf_names() for n in tree_from.get_children() ]
    is_n0_bigger_than_n1 = (len(subroot_leaves[0]) > len(subroot_leaves[1]))
    ingroups = subroot_leaves[0] if is_n0_bigger_than_n1 else subroot_leaves[1]
    outgroups = subroot_leaves[0] if not is_n0_bigger_than_n1 else subroot_leaves[1]
    if verbose:
        sys.stderr.write('Outgroups: {}\n'.format(outgroups))
    tree_to.set_outgroup(ingroups[0])
    if (len(outgroups) == 1):
        outgroup_ancestor = [n for n in tree_to.iter_leaves() if n.name == outgroups[0]][0]
    else:
        outgroup_ancestor = tree_to.get_common_ancestor(outgroups)
    if not set(outgroups) == set(outgroup_ancestor.get_leaf_names()):
        sys.stderr.write('No root bipartition found in --infile. Exiting.\n')
        sys.exit(1)
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

def midpoint_rooting(tree):
    outgroup_node = tree.get_midpoint_outgroup()
    tree.set_outgroup(outgroup_node)
    return tree

def mad_rooting(tree):
    """MAD (Minimal Ancestor Deviation) rooting. Tria et al. 2017, DOI:10.1038/s41559-017-0193"""
    import os, subprocess, tempfile
    mad_script = os.path.join(os.path.dirname(__file__), '_mad.py')
    with tempfile.NamedTemporaryFile(suffix='.nwk', mode='w', delete=False) as f:
        f.write(tree.write(format=5, dist_formatter='%0.8f'))
        tmpfile = f.name
    try:
        result = subprocess.run(
            [sys.executable, mad_script, tmpfile],
            check=True, capture_output=True, text=True,
        )
        sys.stderr.write(result.stdout)
        with open(tmpfile + '.rooted') as fh:
            lines = [l.strip() for l in fh if l.strip() and l.strip().endswith(';')]
        rooted_nwk = lines[0]
        rooted_tree = TreeNode(newick=rooted_nwk, format=1)
    finally:
        for p in [tmpfile, tmpfile + '.rooted']:
            if os.path.exists(p):
                os.unlink(p)
    return rooted_tree

def outgroup_rooting(tree, outgroup_str):
    outgroup_list = outgroup_str.split(',')
    sys.stderr.write('Specified outgroup labels: {}\n'.format(' '.join(outgroup_list)))
    outgroup_nodes = [ node for node in tree.traverse() if node.name in outgroup_list ]
    if len(outgroup_nodes)==0:
        sys.stderr.write('Outgroup node not found. Exiting.\n')
        sys.exit(1)
    elif len(outgroup_nodes)==1:
        outgroup_node = outgroup_nodes[0]
    else:
        outgroup_node = outgroup_nodes[0].get_common_ancestor(outgroup_nodes)
    if outgroup_node is tree: # Reroot if the outgroup clade represents the whole tree
        non_outgroup_leaves = [ node for node in tree.iter_leaves() if node not in outgroup_nodes ]
        tree.set_outgroup(non_outgroup_leaves[0])
        outgroup_node = outgroup_nodes[0].get_common_ancestor(outgroup_nodes)
    if outgroup_node is tree:
        sys.stderr.write('Outgroup clade should not represent the whole tree. Please check --outgroup carefully. Exiting.\n')
        sys.exit(1)
    outgroup_leaf_names = outgroup_node.get_leaf_names()
    sys.stderr.write('All leaf labels in the outgroup clade: {}\n'.format(' '.join(outgroup_leaf_names)))
    tree.set_outgroup(outgroup_node)
    return tree

def root_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    if (args.method=='transfer'):
        tree2 = read_tree(args.infile2, args.format2, args.quoted_node_names)
        if not is_all_leaf_names_identical(tree, tree2, verbose=True):
            raise Exception('Leaf labels in the two trees should be completely matched.')
        tree = transfer_root(tree_to=tree, tree_from=tree2, verbose=True)
    elif (args.method=='midpoint'):
        tree = midpoint_rooting(tree=tree)
    elif (args.method=='outgroup'):
        tree = outgroup_rooting(tree=tree, outgroup_str=args.outgroup)
    elif (args.method=='mad'):
        tree = mad_rooting(tree=tree)
    write_tree(tree, args, format=args.outformat)