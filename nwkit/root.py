import sys

from ete4 import Tree
from ete4.parser.newick import make_parser

from nwkit.util import *

def transfer_root(tree_to, tree_from, verbose=False):
    subroot_from = tree_from.get_children()
    tree_from_leaf_sets = tree_from.get_cached_content(prop='name')
    is_n0_bigger_than_n1 = (len(tree_from_leaf_sets[subroot_from[0]]) > len(tree_from_leaf_sets[subroot_from[1]]))
    ingroup_child = subroot_from[0] if is_n0_bigger_than_n1 else subroot_from[1]
    outgroup_child = subroot_from[1] if is_n0_bigger_than_n1 else subroot_from[0]
    ingroup_set = tree_from_leaf_sets[ingroup_child]
    outgroup_set = tree_from_leaf_sets[outgroup_child]
    if verbose:
        sys.stderr.write('Outgroups: {}\n'.format(' '.join(sorted(outgroup_set))))
    # Save original root name before set_outgroup (ete4 loses it)
    original_root_name = tree_to.name
    root_children = tree_to.get_children()
    tree_to_leaf_sets = tree_to.get_cached_content(prop='name')
    is_root_bipartition_already_matching = (
        (len(root_children) == 2) and
        any(outgroup_set == tree_to_leaf_sets[child] for child in root_children)
    )
    support_backup = None
    if not is_root_bipartition_already_matching:
        support_backup = list()
        # Ensure all None dists are 0 and clear support before rerooting.
        for node in tree_to.traverse():
            if node.dist is None:
                node.dist = 0.0
            support_backup.append((node, node.support))
            node.support = None
        tree_to.set_outgroup(next(iter(ingroup_set)))
        if len(outgroup_set) == 1:
            outgroup_name = next(iter(outgroup_set))
            outgroup_ancestor = None
            for leaf in tree_to.leaves():
                if leaf.name == outgroup_name:
                    outgroup_ancestor = leaf
                    break
            if outgroup_ancestor is None:
                sys.stderr.write('No root bipartition found in --infile. Exiting.\n')
                sys.exit(1)
        else:
            outgroup_ancestor = tree_to.common_ancestor(outgroup_set)
        reroot_leaf_sets = tree_to.get_cached_content(prop='name')
        if not outgroup_set == reroot_leaf_sets[outgroup_ancestor]:
            sys.stderr.write('No root bipartition found in --infile. Exiting.\n')
            sys.exit(1)
        tree_to.set_outgroup(outgroup_ancestor)
        tree_to_leaf_sets = tree_to.get_cached_content(prop='name')
    subroot_to = tree_to.get_children()
    total_subroot_length_to = sum((n.dist or 0) for n in subroot_to)
    total_subroot_length_from = sum((n.dist or 0) for n in subroot_from)
    for n_to in subroot_to:
        n_to_leaf_set = tree_to_leaf_sets[n_to]
        for n_from in subroot_from:
            if (n_to_leaf_set == tree_from_leaf_sets[n_from]):
                n_to.dist = total_subroot_length_to * (n_from.dist / total_subroot_length_from)
    # Restore root name or assign 'Root' if it was unnamed
    if original_root_name:
        tree_to.name = original_root_name
    else:
        tree_to.name = 'Root'
    if support_backup is not None:
        for node, support in support_backup:
            node.support = support
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

def _collect_leaf_distance_stats(tree):
    # For each node, store (leaf_count, sum_dist_to_leaves, sumsq_dist_to_leaves).
    subtree_stats = dict()
    for node in tree.traverse(strategy='postorder'):
        if node.is_leaf:
            subtree_stats[node] = (1, 0.0, 0.0)
            continue
        leaf_count = 0
        sum_dist = 0.0
        sumsq_dist = 0.0
        for child in node.get_children():
            child_count, child_sum, child_sumsq = subtree_stats[child]
            edge_length = float(child.dist) if (child.dist is not None) else 0.0
            leaf_count += child_count
            sum_dist += child_sum + (child_count * edge_length)
            sumsq_dist += child_sumsq + (2.0 * edge_length * child_sum) + (child_count * edge_length * edge_length)
        subtree_stats[node] = (leaf_count, sum_dist, sumsq_dist)
    all_stats = {tree: subtree_stats[tree]}
    for node in tree.traverse(strategy='preorder'):
        parent_count, parent_sum, parent_sumsq = all_stats[node]
        for child in node.get_children():
            child_count, child_sum, child_sumsq = subtree_stats[child]
            edge_length = float(child.dist) if (child.dist is not None) else 0.0
            child_sum_from_parent = child_sum + (child_count * edge_length)
            child_sumsq_from_parent = child_sumsq + (2.0 * edge_length * child_sum) + (child_count * edge_length * edge_length)
            outside_count = parent_count - child_count
            outside_sum_from_parent = parent_sum - child_sum_from_parent
            outside_sumsq_from_parent = parent_sumsq - child_sumsq_from_parent
            outside_sum_from_child = outside_sum_from_parent + (outside_count * edge_length)
            outside_sumsq_from_child = (
                outside_sumsq_from_parent +
                (2.0 * edge_length * outside_sum_from_parent) +
                (outside_count * edge_length * edge_length)
            )
            all_stats[child] = (
                parent_count,
                child_sum + outside_sum_from_child,
                child_sumsq + outside_sumsq_from_child,
            )
    return subtree_stats, all_stats

def mv_rooting(tree):
    """Minimum Variance rooting. Mai, Saeedian & Mirarab 2017, DOI:10.1371/journal.pone.0182238"""
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
    subtree_stats, all_stats = _collect_leaf_distance_stats(tree)
    total_leaf_count = subtree_stats[tree][0]
    best_var = float('inf')
    best_node = None
    best_x = 0.0
    best_L = 0.0
    for node in tree.traverse():
        if node.is_root:
            continue
        L = float(node.dist) if (node.dist is not None) else 0.0
        subtree_leaf_count, subtree_sum, subtree_sumsq = subtree_stats[node]
        all_leaf_count, all_sum, all_sumsq = all_stats[node]
        other_leaf_count = all_leaf_count - subtree_leaf_count
        if (subtree_leaf_count == 0) or (other_leaf_count == 0):
            continue
        other_sum = all_sum - subtree_sum
        other_sumsq = all_sumsq - subtree_sumsq
        mean_a = subtree_sum / subtree_leaf_count
        mean_b = (other_sum / other_leaf_count) - L
        # Optimal root position x from node toward parent, constrained to [0, L]
        x = (L + mean_b - mean_a) / 2.0
        x = max(0.0, min(float(L), x))
        # Root-to-tip distance variance at this x (computed from sums and sums of squares).
        total_sum = all_sum + ((subtree_leaf_count - other_leaf_count) * x)
        total_sumsq = (
            subtree_sumsq +
            (2.0 * x * subtree_sum) +
            (subtree_leaf_count * x * x) +
            other_sumsq -
            (2.0 * x * other_sum) +
            (other_leaf_count * x * x)
        )
        mean = total_sum / total_leaf_count
        var = (total_sumsq / total_leaf_count) - (mean * mean)
        if (var < 0) and (abs(var) < 10**-12):
            var = 0.0
        if var < best_var:
            best_var = var
            best_node = node
            best_x = x
            best_L = L
    if best_node is None:
        return tree
    sys.stderr.write('MV rooting variance: {:.6g}\n'.format(best_var))
    tree.set_outgroup(best_node)
    # Adjust branch lengths at root: best_x from root to best_node, (L - best_x) to sibling
    best_subtree_leaves = set(best_node.leaf_names())
    root_leaf_sets = tree.get_cached_content(prop='name')
    for child in tree.get_children():
        if root_leaf_sets[child] == best_subtree_leaves:
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
    outgroup_name_set = set(outgroup_list)
    outgroup_nodes = [node for node in tree.traverse() if node.name in outgroup_name_set]
    if len(outgroup_nodes)==0:
        sys.stderr.write('Outgroup node not found. Exiting.\n')
        sys.exit(1)
    elif len(outgroup_nodes)==1:
        outgroup_node = outgroup_nodes[0]
    else:
        outgroup_node = tree.common_ancestor(outgroup_nodes)
    if outgroup_node is tree: # Reroot if the outgroup clade represents the whole tree
        outgroup_node_set = set(outgroup_nodes)
        non_outgroup_leaf = None
        for node in tree.leaves():
            if node not in outgroup_node_set:
                non_outgroup_leaf = node
                break
        if non_outgroup_leaf is None:
            sys.stderr.write('Outgroup clade should not represent the whole tree. Please check --outgroup carefully. Exiting.\n')
            sys.exit(1)
        tree.set_outgroup(non_outgroup_leaf)
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
