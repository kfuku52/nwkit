import sys
import re

import requests
from ete4 import Tree
from ete4.parser.newick import make_parser

from nwkit.constrain import (
    get_lineages,
    get_lineages_from_taxid,
    get_taxid_counts,
    limit_lineage_to_rank,
    name_to_taxid,
    read_taxid_tsv,
    taxid2tree,
)
from nwkit.util import *

SUPPORTED_TAXONOMY_SOURCES = ('ncbi', 'timetree', 'opentree')
DEFAULT_TAXONOMY_SOURCE_CHAIN = 'ncbi,opentree,timetree'
NCBI_PLACEHOLDER_NAME_PATTERNS = (
    re.compile(r'\bunknown\b', flags=re.I),
    re.compile(r'\bunidentified\b', flags=re.I),
    re.compile(r'\bunclassified\b', flags=re.I),
    re.compile(r'\buncultured\b', flags=re.I),
    re.compile(r'\benvironmental\b', flags=re.I),
    re.compile(r'\bmetagenom(?:e|es|ic)\b', flags=re.I),
    re.compile(r'\bartificial sequences?\b', flags=re.I),
    re.compile(r'\bother sequences?\b', flags=re.I),
)

def _close_ncbi_handle(ncbi):
    if ncbi is None:
        return
    db = getattr(ncbi, 'db', None)
    if db is not None:
        try:
            db.close()
        except Exception:
            pass

def _is_placeholder_ncbi_name(name):
    normalized = re.sub(r'[_\s]+', ' ', str(name or '').strip()).lower()
    if normalized == '':
        return False
    return any(pattern.search(normalized) for pattern in NCBI_PLACEHOLDER_NAME_PATTERNS)

def _is_placeholder_ncbi_resolution(matched_name, lineage_names):
    names_to_check = list()
    if matched_name not in ['', None]:
        names_to_check.append(matched_name)
    names_to_check.extend([name for name in lineage_names if name not in ['', None]])
    return any(_is_placeholder_ncbi_name(name) for name in names_to_check)

def _get_ncbi_lineage_record(label, ncbi, rank, taxid=None):
    if taxid is None:
        query_name = label2sciname(labels=label, in_delim='_', out_delim=' ') if '_' in str(label) else str(label)
        name2taxid = name_to_taxid(query_name, ncbi)
        if len(name2taxid) == 0:
            return None, 'Genus-level match was not found in the NCBI database'
        matched_query_name = query_name if query_name in name2taxid.keys() else re.sub(' .*', '', query_name)
        taxid = int(name2taxid[matched_query_name][0])
    else:
        matched_query_name = None
        taxid = int(taxid)
    lineage = [int(t) for t in ncbi.get_lineage(taxid)]
    lineage_name_map = ncbi.get_taxid_translator(lineage + [taxid])
    lineage_names = [lineage_name_map.get(t) for t in lineage]
    matched_name = lineage_name_map.get(taxid, matched_query_name or str(taxid))
    if _is_placeholder_ncbi_resolution(matched_name, lineage_names):
        return None, 'matched placeholder NCBI taxon: {}'.format(matched_name)
    return {
        'lineage': limit_lineage_to_rank(lineage, ncbi, rank),
        'matched_taxid': taxid,
        'matched_name': matched_name,
    }, None

def _resolve_ncbi_lineages(tree, taxid_tsv=None, rank='no', args=None, verbose=False):
    ncbi = get_ete_ncbitaxa(args=args)
    try:
        lineages = dict()
        unresolved_details = dict()
        if taxid_tsv not in ['', None]:
            taxid_df = read_taxid_tsv(taxid_tsv)
            taxid_df = _order_taxid_tsv_to_match_tree(tree, taxid_df)
            records = zip(taxid_df['leaf_name'], taxid_df['taxid'])
            for label, taxid in records:
                record, reason = _get_ncbi_lineage_record(label=label, taxid=taxid, ncbi=ncbi, rank=rank)
                if record is None:
                    unresolved_details[label] = reason
                    continue
                lineages[label] = record['lineage']
        else:
            for label in tree.leaf_names():
                record, reason = _get_ncbi_lineage_record(label=label, ncbi=ncbi, rank=rank)
                if record is None:
                    unresolved_details[label] = reason
                    continue
                lineages[label] = record['lineage']
        if verbose and unresolved_details:
            details = [
                '{} ({})'.format(label, unresolved_details[label])
                for label in sorted(unresolved_details.keys())
            ]
            sys.stderr.write(
                'Excluding NCBI-unresolved leaf label(s) from taxonomy rooting: {}\n'.format(
                    '; '.join(details)
                )
            )
        return lineages, unresolved_details
    finally:
        _close_ncbi_handle(ncbi)

def _normalize_root_distance_for_reroot(tree):
    if tree.dist is not None:
        tree.dist = 0.0
    return tree

def _collapse_singleton_root(tree):
    while (len(tree.get_children()) == 1) and (len(list(tree.leaves())) > 1):
        child = tree.get_children()[0]
        tree = Tree(child.write(parser=0, format_root_node=True), parser=0)
    return tree

def transfer_root(tree_to, tree_from, verbose=False):
    tree_to = _collapse_singleton_root(tree_to)
    tree_from = _collapse_singleton_root(tree_from)
    validate_unique_named_leaves(tree_to, option_name='--infile', context=' for root transfer')
    validate_unique_named_leaves(tree_from, option_name='--infile2', context=' for root transfer')
    subroot_from = tree_from.get_children()
    if len(subroot_from) != 2:
        raise ValueError('Root transfer requires the source tree root to have exactly two children.')
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
        _normalize_root_distance_for_reroot(tree_to)
        tree_to.set_outgroup(next(iter(ingroup_set)))
        if len(outgroup_set) == 1:
            outgroup_name = next(iter(outgroup_set))
            outgroup_ancestor = None
            for leaf in tree_to.leaves():
                if leaf.name == outgroup_name:
                    outgroup_ancestor = leaf
                    break
            if outgroup_ancestor is None:
                raise ValueError('No root bipartition found in --infile.')
        else:
            outgroup_ancestor = tree_to.common_ancestor(outgroup_set)
        reroot_leaf_sets = tree_to.get_cached_content(prop='name')
        if not outgroup_set == reroot_leaf_sets[outgroup_ancestor]:
            raise ValueError('No root bipartition found in --infile.')
        _normalize_root_distance_for_reroot(tree_to)
        tree_to.set_outgroup(outgroup_ancestor)
        tree_to_leaf_sets = tree_to.get_cached_content(prop='name')
    subroot_to = tree_to.get_children()
    total_subroot_length_to = sum((n.dist or 0) for n in subroot_to)
    total_subroot_length_from = sum((n.dist or 0) for n in subroot_from)
    if abs(total_subroot_length_from) > 10**-15:
        for n_to in subroot_to:
            n_to_leaf_set = tree_to_leaf_sets[n_to]
            for n_from in subroot_from:
                if (n_to_leaf_set == tree_from_leaf_sets[n_from]):
                    n_to.dist = total_subroot_length_to * ((n_from.dist or 0) / total_subroot_length_from)
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
    _normalize_root_distance_for_reroot(tree)
    outgroup_node = tree.get_midpoint_outgroup()
    # If the outgroup is the root itself, tree is already optimally rooted
    if outgroup_node.is_root:
        return tree
    tree.set_outgroup(outgroup_node)
    return tree

def mad_rooting(tree):
    """MAD (Minimal Ancestor Deviation) rooting. Tria et al. 2017, DOI:10.1038/s41559-017-0193"""
    if len(list(tree.leaves())) < 3:
        raise ValueError('MAD rooting requires at least 3 leaves.')
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
    # ete4.set_outgroup requires root.dist to be 0/None.
    if tree.dist is not None:
        tree.dist = 0.0
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
    if outgroup_str is None:
        raise ValueError("Specify at least one outgroup label with '--outgroup'.")
    # Ensure all None dists are 0 (ete4 set_outgroup doesn't handle None well)
    for node in tree.traverse():
        if node.dist is None:
            node.dist = 0.0
    outgroup_list = [label.strip() for label in outgroup_str.split(',') if label.strip()]
    if len(outgroup_list) == 0:
        raise ValueError("Specify at least one outgroup label with '--outgroup'.")
    sys.stderr.write('Specified outgroup labels: {}\n'.format(' '.join(outgroup_list)))
    outgroup_name_set = set(outgroup_list)
    leaf_name_set = set(tree.leaf_names())
    missing_outgroup_names = [name for name in outgroup_list if name not in leaf_name_set]
    if missing_outgroup_names:
        raise ValueError('Outgroup label(s) not found in leaf names: {}'.format(', '.join(missing_outgroup_names)))
    outgroup_nodes = [node for node in tree.leaves() if node.name in outgroup_name_set]
    if len(outgroup_nodes)==0:
        raise ValueError('Outgroup node not found.')
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
            raise ValueError('Outgroup clade should not represent the whole tree. Please check --outgroup carefully.')
        _normalize_root_distance_for_reroot(tree)
        tree.set_outgroup(non_outgroup_leaf)
        outgroup_node = tree.common_ancestor(outgroup_nodes)
    if outgroup_node is tree:
        raise ValueError('Outgroup clade should not represent the whole tree. Please check --outgroup carefully.')
    outgroup_leaf_names = list(outgroup_node.leaf_names())
    sys.stderr.write('All leaf labels in the outgroup clade: {}\n'.format(' '.join(outgroup_leaf_names)))
    _normalize_root_distance_for_reroot(tree)
    tree.set_outgroup(outgroup_node)
    return tree

def _order_taxid_tsv_to_match_tree(tree, taxid_df):
    tree_leaf_names = list(tree.leaf_names())
    tree_leaf_set = set(tree_leaf_names)
    taxid_leaf_set = set(taxid_df['leaf_name'])
    missing_leaf_names = sorted(tree_leaf_set - taxid_leaf_set)
    extra_leaf_names = sorted(taxid_leaf_set - tree_leaf_set)
    if missing_leaf_names or extra_leaf_names:
        messages = list()
        if missing_leaf_names:
            messages.append('missing leaf_name entries for: {}'.format(', '.join(missing_leaf_names)))
        if extra_leaf_names:
            messages.append('unexpected leaf_name entries: {}'.format(', '.join(extra_leaf_names)))
        raise ValueError('--taxid_tsv must match the leaf labels in --infile exactly ({})'.format('; '.join(messages)))
    return taxid_df.set_index('leaf_name').loc[tree_leaf_names].reset_index()

def _parse_taxonomy_sources(taxonomy_source):
    if isinstance(taxonomy_source, (list, tuple)):
        raw_sources = taxonomy_source
    else:
        raw_sources = str(taxonomy_source or DEFAULT_TAXONOMY_SOURCE_CHAIN).split(',')
    sources = list()
    seen = set()
    for raw_source in raw_sources:
        source = str(raw_source).strip().lower()
        if source == '':
            continue
        if source not in SUPPORTED_TAXONOMY_SOURCES:
            raise ValueError(
                "Unknown taxonomy source: {}. Supported sources are: {}.".format(
                    source,
                    ', '.join(SUPPORTED_TAXONOMY_SOURCES),
                )
            )
        if source in seen:
            continue
        seen.add(source)
        sources.append(source)
    if len(sources) == 0:
        raise ValueError('Specify at least one taxonomy source.')
    return sources

def _build_ncbi_reference_tree(tree, taxid_tsv=None, rank='no', args=None, verbose=False):
    lineages, unresolved_details = _resolve_ncbi_lineages(
        tree=tree,
        taxid_tsv=taxid_tsv,
        rank=rank,
        args=args,
        verbose=verbose,
    )
    unresolved_labels = set(unresolved_details.keys())
    if len(lineages) == 0:
        raise ValueError('Failed to resolve usable NCBI lineage for any leaf label.')
    if (len(list(tree.leaves())) > 1) and (len(lineages) < 2):
        raise ValueError(
            'At least two usable NCBI-resolved leaf labels are required after excluding unresolved labels.'
        )
    taxonomy_tree = taxid2tree(lineages, get_taxid_counts(lineages), args=args)
    if set(taxonomy_tree.leaf_names()) != set(lineages.keys()):
        raise ValueError('Leaf labels in the NCBI-derived tree should match the resolved leaf labels in --infile.')
    taxonomy_tree = _collapse_singleton_root(taxonomy_tree)
    if (len(list(taxonomy_tree.leaves())) > 1) and (len(taxonomy_tree.get_children()) != 2):
        raise ValueError(
            'NCBI-derived root is ambiguous: expected exactly 2 root children, found {}.'.format(
                len(taxonomy_tree.get_children())
            )
        )
    return taxonomy_tree, set(lineages.keys()), unresolved_labels

def _resolve_full_outgroup_set_from_resolved_split(tree, resolved_outgroup_set, resolved_leaf_set, source_name):
    analysis_tree = _collapse_singleton_root(tree.copy())
    validate_unique_named_leaves(analysis_tree, option_name='--infile', context=' for taxonomy rooting')
    resolved_leaf_set = set(resolved_leaf_set)
    resolved_outgroup_set = set(resolved_outgroup_set)
    leaf_name_set = set(analysis_tree.leaf_names())
    if len(resolved_outgroup_set) == 0:
        raise ValueError('No root bipartition found in --infile.')
    if not resolved_outgroup_set < resolved_leaf_set:
        raise ValueError('No root bipartition found in --infile.')
    if not resolved_leaf_set <= leaf_name_set:
        raise ValueError('No root bipartition found in --infile.')
    root_children = analysis_tree.get_children()
    subtree_leaf_sets = analysis_tree.get_cached_content(prop='name')
    matching_edges = dict()
    resolved_ingroup_set = resolved_leaf_set - resolved_outgroup_set
    for node in analysis_tree.traverse():
        if node.is_root:
            continue
        node_leaf_set = set(subtree_leaf_sets[node])
        node_resolved_set = node_leaf_set.intersection(resolved_leaf_set)
        if (len(node_resolved_set) == 0) or (node_resolved_set == resolved_leaf_set):
            continue
        if node_resolved_set == resolved_outgroup_set:
            full_outgroup_set = node_leaf_set
        elif node_resolved_set == resolved_ingroup_set:
            full_outgroup_set = leaf_name_set - node_leaf_set
        else:
            continue
        edge_key = ('root_edge',) if ((node.up is analysis_tree) and (len(root_children) == 2)) else node
        previous = matching_edges.get(edge_key)
        if previous is None:
            matching_edges[edge_key] = set(full_outgroup_set)
        elif previous != set(full_outgroup_set):
            raise ValueError('Internal error while deduplicating candidate root edges.')
    if len(matching_edges) == 0:
        raise ValueError('No root bipartition found in --infile.')
    if len(matching_edges) > 1:
        raise ValueError(
            'Unresolved leaf clade(s) interfere with the {}-derived root position in --infile.'.format(
                source_name
            )
        )
    return next(iter(matching_edges.values()))

def _transfer_root_from_reference_with_unresolved(tree_to, tree_from, resolved_leaf_set, verbose=False, source_name='ncbi'):
    resolved_outgroup_set = _get_reference_root_outgroup_set(tree_from, source_name=source_name.upper())
    full_outgroup_set = _resolve_full_outgroup_set_from_resolved_split(
        tree=tree_to,
        resolved_outgroup_set=resolved_outgroup_set,
        resolved_leaf_set=resolved_leaf_set,
        source_name=source_name,
    )
    return _root_by_outgroup_set(tree=tree_to, outgroup_set=full_outgroup_set, verbose=verbose)

def _get_reference_root_outgroup_set(reference_tree, source_name):
    reference_tree = _collapse_singleton_root(reference_tree)
    if (len(list(reference_tree.leaves())) > 1) and (len(reference_tree.get_children()) != 2):
        raise ValueError(
            '{}-derived root is ambiguous: expected exactly 2 root children, found {}.'.format(
                source_name,
                len(reference_tree.get_children()),
            )
        )
    if len(reference_tree.get_children()) != 2:
        return set(reference_tree.leaf_names())
    child_leaf_sets = reference_tree.get_cached_content(prop='name')
    root_children = reference_tree.get_children()
    is_n0_bigger_than_n1 = (len(child_leaf_sets[root_children[0]]) > len(child_leaf_sets[root_children[1]]))
    outgroup_child = root_children[1] if is_n0_bigger_than_n1 else root_children[0]
    return set(child_leaf_sets[outgroup_child])

def _expand_species_outgroup_set(species_outgroup_set, species_to_leaf_labels):
    outgroup_leaf_set = set()
    for species_name in species_outgroup_set:
        if species_name not in species_to_leaf_labels:
            raise ValueError('Unexpected species label returned by the reference tree: {}'.format(species_name))
        outgroup_leaf_set.update(species_to_leaf_labels[species_name])
    return outgroup_leaf_set

def _root_by_outgroup_set(tree, outgroup_set, verbose=False):
    validate_unique_named_leaves(tree, option_name='--infile', context=' for taxonomy rooting')
    outgroup_set = set(outgroup_set)
    if len(outgroup_set) == 0:
        raise ValueError('No root bipartition found in --infile.')
    leaf_name_set = set(tree.leaf_names())
    if not outgroup_set <= leaf_name_set:
        raise ValueError('No root bipartition found in --infile.')
    if outgroup_set == leaf_name_set:
        raise ValueError('No root bipartition found in --infile.')
    if verbose:
        sys.stderr.write('Outgroups: {}\n'.format(' '.join(sorted(outgroup_set))))
    original_root_name = tree.name
    root_children = tree.get_children()
    tree_leaf_sets = tree.get_cached_content(prop='name')
    is_root_bipartition_already_matching = (
        (len(root_children) == 2) and
        any(outgroup_set == tree_leaf_sets[child] for child in root_children)
    )
    support_backup = None
    if not is_root_bipartition_already_matching:
        ingroup_set = leaf_name_set - outgroup_set
        support_backup = list()
        for node in tree.traverse():
            if node.dist is None:
                node.dist = 0.0
            support_backup.append((node, node.support))
            node.support = None
        _normalize_root_distance_for_reroot(tree)
        tree.set_outgroup(next(iter(ingroup_set)))
        if len(outgroup_set) == 1:
            outgroup_name = next(iter(outgroup_set))
            outgroup_ancestor = None
            for leaf in tree.leaves():
                if leaf.name == outgroup_name:
                    outgroup_ancestor = leaf
                    break
            if outgroup_ancestor is None:
                raise ValueError('No root bipartition found in --infile.')
        else:
            outgroup_ancestor = tree.common_ancestor(outgroup_set)
        reroot_leaf_sets = tree.get_cached_content(prop='name')
        if outgroup_set != reroot_leaf_sets[outgroup_ancestor]:
            raise ValueError('No root bipartition found in --infile.')
        _normalize_root_distance_for_reroot(tree)
        tree.set_outgroup(outgroup_ancestor)
    if original_root_name:
        tree.name = original_root_name
    else:
        tree.name = 'Root'
    if support_backup is not None:
        for node, support in support_backup:
            node.support = support
    return tree

def _get_timetree_name_mapping(tree):
    _, species_to_leaf_labels = get_monophyletic_species_groups(
        tree,
        option_name='--infile',
        context=' for TimeTree taxonomy rooting',
    )
    query_names = [species_name.replace('_', ' ') for species_name in species_to_leaf_labels.keys()]
    return query_names, species_to_leaf_labels

def _extract_timetree_unresolved_names(upload_response_text):
    unresolved_count_match = re.search(r'Unresolved Names \((\d+)\)', upload_response_text)
    if unresolved_count_match is None:
        return list()
    unresolved_count = int(unresolved_count_match.group(1))
    unresolved_lines = re.findall(r'<div id="unresolved-names">(.*?)</div>', upload_response_text, flags=re.S)
    unresolved_names = list()
    for unresolved_line in unresolved_lines:
        unresolved_names.extend([
            re.sub(r'<[^>]+>', '', name).strip()
            for name in unresolved_line.split('<br/>')
            if re.sub(r'<[^>]+>', '', name).strip() != ''
        ])
    if unresolved_count > 0 and len(unresolved_names) == 0:
        unresolved_names = ['{} unresolved name(s) reported by TimeTree'.format(unresolved_count)]
    return unresolved_names

def _build_timetree_reference_tree(tree):
    query_names, species_to_leaf_labels = _get_timetree_name_mapping(tree)
    session = requests.Session()
    try:
        try:
            upload_response = session.post(
                'https://timetree.org/ajax/prune/load_names/',
                files={'file': ('species.txt', '\n'.join(query_names) + '\n')},
                timeout=60,
            )
            if upload_response.status_code != 200:
                raise ValueError('Failed to retrieve a pruned guide tree from TimeTree.')
            unresolved_names = _extract_timetree_unresolved_names(upload_response.text)
            if unresolved_names:
                raise ValueError(
                    'TimeTree reported unresolved names: {}'.format('; '.join(unresolved_names))
                )
            newick_response = session.post(
                'https://timetree.org/ajax/newick/prunetree/download',
                data={'export': 'newick', 'rank': ''},
                timeout=60,
            )
        except requests.RequestException as exc:
            raise ValueError('Failed to contact TimeTree.') from exc
        if newick_response.status_code != 200:
            raise ValueError('Failed to download a guide tree from TimeTree.')
        timetree_newick = newick_response.text.strip()
        if not timetree_newick.endswith(';'):
            raise ValueError('Unexpected response format from TimeTree when downloading the guide tree.')
        timetree_tree = Tree(timetree_newick, parser=1)
        expected_species_names = set(species_to_leaf_labels.keys())
        if set(timetree_tree.leaf_names()) != expected_species_names:
            raise ValueError('Leaf labels in the TimeTree-derived tree should match the species names derived from --infile.')
        species_outgroup_set = _get_reference_root_outgroup_set(timetree_tree, source_name='TimeTree')
        if len(species_outgroup_set) == 0:
            raise ValueError(
                'TimeTree-derived root is ambiguous: failed to identify a root bipartition.'
            )
        return _expand_species_outgroup_set(species_outgroup_set, species_to_leaf_labels)
    finally:
        close = getattr(session, 'close', None)
        if callable(close):
            close()

def _get_opentree_name_mapping(tree):
    _, species_to_leaf_labels = get_monophyletic_species_groups(
        tree,
        option_name='--infile',
        context=' for OpenTree taxonomy rooting',
    )
    query_names = [species_name.replace('_', ' ') for species_name in species_to_leaf_labels.keys()]
    return query_names, species_to_leaf_labels

def _extract_opentree_ott_ids(tree):
    query_names, species_to_leaf_labels = _get_opentree_name_mapping(tree)
    session = requests.Session()
    try:
        try:
            response = session.post(
                'https://api.opentreeoflife.org/v3/tnrs/match_names',
                json={
                    'names': query_names,
                    'do_approximate_matching': False,
                },
                timeout=60,
            )
        except requests.RequestException as exc:
            raise ValueError('Failed to contact Open Tree of Life TNRS.') from exc
        if response.status_code != 200:
            raise ValueError('Open Tree of Life TNRS lookup failed.')
        data = response.json()
        results = data.get('results', [])
        if len(results) != len(query_names):
            raise ValueError('Unexpected response format from Open Tree of Life TNRS.')
        ott_ids = list()
        for query_name, result in zip(query_names, results):
            matches = result.get('matches', [])
            valid_matches = list()
            for match in matches:
                taxon = match.get('taxon', {})
                if match.get('is_approximate_match'):
                    continue
                if taxon.get('is_suppressed_from_synth'):
                    continue
                if taxon.get('ott_id') is None:
                    continue
                valid_matches.append(match)
            if len(valid_matches) == 0:
                raise ValueError('Failed to resolve an OpenTree taxon for leaf label: {}'.format(
                    query_name.replace(' ', '_')
                ))
            if len(valid_matches) != 1:
                raise ValueError('OpenTree TNRS returned multiple exact matches for leaf label: {}'.format(
                    query_name.replace(' ', '_')
                ))
            ott_ids.append(int(valid_matches[0]['taxon']['ott_id']))
        return ott_ids, species_to_leaf_labels
    finally:
        close = getattr(session, 'close', None)
        if callable(close):
            close()

def _build_opentree_reference_tree(tree):
    ott_ids, species_to_leaf_labels = _extract_opentree_ott_ids(tree)
    session = requests.Session()
    try:
        try:
            response = session.post(
                'https://api.opentreeoflife.org/v3/tree_of_life/induced_subtree',
                json={
                    'ott_ids': ott_ids,
                    'label_format': 'name',
                },
                timeout=60,
            )
        except requests.RequestException as exc:
            raise ValueError('Failed to contact Open Tree of Life synthetic tree API.') from exc
        if response.status_code != 200:
            raise ValueError('Open Tree of Life induced subtree lookup failed.')
        data = response.json()
        newick = data.get('newick', '').strip()
        if not newick.endswith(';'):
            raise ValueError('Unexpected response format from Open Tree of Life induced subtree.')
        opentree_tree = Tree(newick, parser=1)
        opentree_tree = remove_singleton(opentree_tree, verbose=False, preserve_branch_length=True)
        expected_species_names = set(species_to_leaf_labels.keys())
        if set(opentree_tree.leaf_names()) != expected_species_names:
            raise ValueError('Leaf labels in the OpenTree-derived tree should match the species names derived from --infile.')
        species_outgroup_set = _get_reference_root_outgroup_set(opentree_tree, source_name='OpenTree')
        if len(species_outgroup_set) == 0:
            raise ValueError(
                'OpenTree-derived root is ambiguous: failed to identify a root bipartition.'
            )
        return _expand_species_outgroup_set(species_outgroup_set, species_to_leaf_labels)
    finally:
        close = getattr(session, 'close', None)
        if callable(close):
            close()

def _build_taxonomy_reference_tree(tree, taxonomy_source, taxid_tsv=None, rank='no', args=None):
    validate_unique_named_leaves(tree, option_name='--infile', context=' for taxonomy rooting')
    if taxonomy_source == 'ncbi':
        taxonomy_tree, _, _ = _build_ncbi_reference_tree(tree=tree, taxid_tsv=taxid_tsv, rank=rank, args=args)
        return taxonomy_tree, None
    if taxonomy_source == 'timetree':
        return None, _build_timetree_reference_tree(tree=tree)
    if taxonomy_source == 'opentree':
        return None, _build_opentree_reference_tree(tree=tree)
    raise ValueError("Unknown taxonomy source: {}".format(taxonomy_source))

def taxonomy_rooting(tree, taxonomy_source=DEFAULT_TAXONOMY_SOURCE_CHAIN, taxid_tsv=None, rank='no', verbose=False, args=None):
    if len(list(tree.leaves())) <= 1:
        return tree
    errors = list()
    taxonomy_sources = _parse_taxonomy_sources(taxonomy_source)
    for source in taxonomy_sources:
        if verbose:
            sys.stderr.write('Attempting taxonomy rooting with source: {}\n'.format(source))
        try:
            if source == 'ncbi':
                taxonomy_tree, resolved_leaf_set, unresolved_leaf_set = _build_ncbi_reference_tree(
                    tree=tree,
                    taxid_tsv=taxid_tsv,
                    rank=rank,
                    args=args,
                    verbose=verbose,
                )
                if len(unresolved_leaf_set) == 0:
                    return transfer_root(tree_to=tree, tree_from=taxonomy_tree, verbose=verbose)
                return _transfer_root_from_reference_with_unresolved(
                    tree_to=tree,
                    tree_from=taxonomy_tree,
                    resolved_leaf_set=resolved_leaf_set,
                    verbose=verbose,
                    source_name=source,
                )
            taxonomy_tree, outgroup_set = _build_taxonomy_reference_tree(
                tree=tree,
                taxonomy_source=source,
                taxid_tsv=taxid_tsv,
                rank=rank,
                args=args,
            )
            if taxonomy_tree is not None:
                return transfer_root(tree_to=tree, tree_from=taxonomy_tree, verbose=verbose)
            return _root_by_outgroup_set(tree=tree, outgroup_set=outgroup_set, verbose=verbose)
        except ValueError as exc:
            if str(exc) == 'No root bipartition found in --infile.':
                exc = ValueError('The {}-derived root bipartition was not found in --infile.'.format(source))
            errors.append('{}: {}'.format(source, exc))
            if verbose:
                sys.stderr.write('Taxonomy rooting with source {} failed: {}\n'.format(source, exc))
    raise ValueError('All taxonomy sources failed: {}'.format(' | '.join(errors)))

def root_main(args):
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    if (args.method=='transfer'):
        if args.infile2 in ['', None]:
            raise ValueError("'--infile2' is required when '--method transfer' is used.")
        tree2 = read_tree(args.infile2, args.format2, args.quoted_node_names)
        if not is_rooted(tree2):
            raise ValueError("'--infile2' must be rooted when '--method transfer' is used.")
        if (len(list(tree2.leaves())) > 1) and (len(tree2.get_children()) != 2):
            raise ValueError("'--infile2' root must have exactly two children for '--method transfer'.")
        validate_unique_named_leaves(tree, option_name='--infile', context=' for root transfer')
        validate_unique_named_leaves(tree2, option_name='--infile2', context=' for root transfer')
        if not is_all_leaf_names_identical(tree, tree2, verbose=True):
            raise Exception('Leaf labels in the two trees should be completely matched.')
        if (len(list(tree.leaves())) > 1) and (len(list(tree2.leaves())) > 1):
            tree = transfer_root(tree_to=tree, tree_from=tree2, verbose=True)
    elif (args.method=='midpoint'):
        tree = midpoint_rooting(tree=tree)
    elif (args.method=='outgroup'):
        tree = outgroup_rooting(tree=tree, outgroup_str=args.outgroup)
    elif (args.method=='mad'):
        tree = mad_rooting(tree=tree)
    elif (args.method=='mv'):
        tree = mv_rooting(tree=tree)
    elif (args.method=='taxonomy'):
        tree = taxonomy_rooting(
            tree=tree,
            taxonomy_source=getattr(args, 'taxonomy_source', DEFAULT_TAXONOMY_SOURCE_CHAIN),
            taxid_tsv=getattr(args, 'taxid_tsv', None),
            rank=getattr(args, 'rank', 'no'),
            verbose=True,
            args=args,
        )
    else:
        raise ValueError("Unknown rooting method: {}".format(args.method))
    write_tree(tree, args, format=args.outformat)
