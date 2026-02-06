import sys

from ete4 import Tree
from ete4.parser.newick import make_parser

from nwkit.util import *

def transfer_root(tree_to, tree_from, verbose=False):
    # Ensure all None dists are 0 (ete4 set_outgroup doesn't handle None well)
    for node in tree_to.traverse():
        if node.dist is None:
            node.dist = 0.0
    subroot_leaves = [ list(n.leaf_names()) for n in tree_from.get_children() ]
    is_n0_bigger_than_n1 = (len(subroot_leaves[0]) > len(subroot_leaves[1]))
    ingroups = subroot_leaves[0] if is_n0_bigger_than_n1 else subroot_leaves[1]
    outgroups = subroot_leaves[0] if not is_n0_bigger_than_n1 else subroot_leaves[1]
    if verbose:
        sys.stderr.write('Outgroups: {}\n'.format(outgroups))
    # Save original root name before set_outgroup (ete4 loses it)
    original_root_name = tree_to.name
    tree_to.set_outgroup(ingroups[0])
    if (len(outgroups) == 1):
        outgroup_ancestor = [n for n in tree_to.leaves() if n.name == outgroups[0]][0]
    else:
        outgroup_ancestor = tree_to.common_ancestor(outgroups)
    if not set(outgroups) == set(outgroup_ancestor.leaf_names()):
        sys.stderr.write('No root bipartition found in --infile. Exiting.\n')
        sys.exit(1)
    tree_to.set_outgroup(outgroup_ancestor)
    subroot_to = tree_to.get_children()
    subroot_from = tree_from.get_children()
    total_subroot_length_to = sum([(n.dist or 0) for n in subroot_to])
    total_subroot_length_from = sum([(n.dist or 0) for n in subroot_from])
    for n_to in subroot_to:
        for n_from in subroot_from:
            if (set(n_to.leaf_names()) == set(n_from.leaf_names())):
                n_to.dist = total_subroot_length_to * (n_from.dist / total_subroot_length_from)
    # Restore root name or assign 'Root' if it was unnamed
    if original_root_name:
        tree_to.name = original_root_name
    else:
        tree_to.name = 'Root'
    return tree_to

def midpoint_rooting(tree):
    # Ensure all None dists are 0 (ete4 set_outgroup doesn't handle None well)
    for node in tree.traverse():
        if node.dist is None:
            node.dist = 0.0
    outgroup_node = tree.get_midpoint_outgroup()
    # If the outgroup is the root itself, tree is already optimally rooted
    if outgroup_node.is_root:
        return tree
    tree.set_outgroup(outgroup_node)
    return tree

def mad_rooting(tree):
    """MAD (Minimal Ancestor Deviation) rooting. Tria et al. 2017, DOI:10.1038/s41559-017-0193"""
    import os, subprocess, tempfile
    mad_script = os.path.join(os.path.dirname(__file__), '_mad.py')
    parser = make_parser(5, dist='%0.8f')
    with tempfile.NamedTemporaryFile(suffix='.nwk', mode='w', delete=False) as f:
        f.write(tree.write(parser=parser))
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
        rooted_tree = Tree(rooted_nwk, parser=1)
    finally:
        for p in [tmpfile, tmpfile + '.rooted']:
            if os.path.exists(p):
                os.unlink(p)
    return rooted_tree

def mv_rooting(tree):
    """Minimum Variance rooting. Mai, Saeedian & Mirarab 2017, DOI:10.1371/journal.pone.0182238"""
    import numpy as np
    # Unroot bifurcating root so each edge is a proper edge in the unrooted tree.
    # Manual unroot because unroot() can drop the dissolved node's branch length.
    children = tree.get_children()
    if len(children) == 2:
        c0, c1 = children
        to_dissolve = c0 if not c0.is_leaf else (c1 if not c1.is_leaf else None)
        if to_dissolve is not None:
            to_keep = c1 if to_dissolve is c0 else c0
            to_keep.dist = (to_keep.dist or 0) + (to_dissolve.dist or 0)
            for gc in list(to_dissolve.get_children()):
                tree.add_child(gc)
            tree.remove_child(to_dissolve)
    leaves = list(tree.leaves())
    all_leaf_names = set(tree.leaf_names())
    best_var = float('inf')
    best_node = None
    best_x = 0.0
    best_L = 0.0
    for node in tree.traverse():
        if node.is_root:
            continue
        L = node.dist if node.dist is not None else 0.0
        subtree_leaf_names = set(node.leaf_names())
        other_leaf_names = all_leaf_names - subtree_leaf_names
        if not subtree_leaf_names or not other_leaf_names:
            continue
        # Distances from node to all leaves
        dists_from_node = {leaf.name: tree.get_distance(node, leaf) for leaf in leaves}
        a = [dists_from_node[name] for name in subtree_leaf_names]
        b = [dists_from_node[name] - L for name in other_leaf_names]
        mean_a = np.mean(a)
        mean_b = np.mean(b)
        # Optimal root position x from node toward parent, constrained to [0, L]
        x = (L + mean_b - mean_a) / 2.0
        x = max(0.0, min(float(L), x))
        # Root-to-tip distances at this position
        all_dists = [ai + x for ai in a] + [bj + (L - x) for bj in b]
        var = np.var(all_dists)
        if var < best_var:
            best_var = var
            best_node = node
            best_x = x
            best_L = L
    sys.stderr.write('MV rooting variance: {:.6g}\n'.format(best_var))
    tree.set_outgroup(best_node)
    # Adjust branch lengths at root: best_x from root to best_node, (L - best_x) to sibling
    best_subtree_leaves = set(best_node.leaf_names())
    for child in tree.get_children():
        if set(child.leaf_names()) == best_subtree_leaves:
            child.dist = best_x
        else:
            child.dist = best_L - best_x
    return tree

def outgroup_rooting(tree, outgroup_str):
    # Ensure all None dists are 0 (ete4 set_outgroup doesn't handle None well)
    for node in tree.traverse():
        if node.dist is None:
            node.dist = 0.0
    outgroup_list = outgroup_str.split(',')
    sys.stderr.write('Specified outgroup labels: {}\n'.format(' '.join(outgroup_list)))
    outgroup_nodes = [ node for node in tree.traverse() if node.name in outgroup_list ]
    if len(outgroup_nodes)==0:
        sys.stderr.write('Outgroup node not found. Exiting.\n')
        sys.exit(1)
    elif len(outgroup_nodes)==1:
        outgroup_node = outgroup_nodes[0]
    else:
        outgroup_node = tree.common_ancestor(outgroup_nodes)
    if outgroup_node is tree: # Reroot if the outgroup clade represents the whole tree
        non_outgroup_leaves = [ node for node in tree.leaves() if node not in outgroup_nodes ]
        tree.set_outgroup(non_outgroup_leaves[0])
        outgroup_node = tree.common_ancestor(outgroup_nodes)
    if outgroup_node is tree:
        sys.stderr.write('Outgroup clade should not represent the whole tree. Please check --outgroup carefully. Exiting.\n')
        sys.exit(1)
    outgroup_leaf_names = list(outgroup_node.leaf_names())
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
    elif (args.method=='mv'):
        tree = mv_rooting(tree=tree)
    write_tree(tree, args, format=args.outformat)
