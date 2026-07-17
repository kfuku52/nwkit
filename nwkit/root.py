import sys
import re
from collections import defaultdict

import requests
from ete4 import Tree
from ete4.parser.newick import make_parser

from nwkit.constrain import (
    get_taxid_counts,
    limit_lineage_to_rank,
    name_to_taxid,
    read_taxid_tsv,
    taxid2tree,
)
from nwkit.util import (
    TREE_FORMAT_PROP,
    extract_taxonomy_query,
    get_ete_ncbitaxa,
    get_species_group_records,
    is_all_leaf_names_identical,
    is_rooted,
    read_tree,
    remove_singleton,
    support_is_missing,
    validate_unique_named_leaves,
    warn_cleanup_failure,
    write_tree,
)
from nwkit.clade_mapping import (
    canonical_split,
    find_root_split_candidates,
    projected_root_split,
)

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
        except Exception as exc:
            warn_cleanup_failure('NCBI taxonomy database handle', exc)

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

def _get_ncbi_lineage_record(label, ncbi, rank, args=None, taxid=None):
    if taxid is None:
        query_name = extract_taxonomy_query(label, args=args, out_delim=' ')
        if query_name is None:
            query_name = str(label)
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
                record, reason = _get_ncbi_lineage_record(
                    label=label,
                    taxid=taxid,
                    ncbi=ncbi,
                    rank=rank,
                    args=args,
                )
                if record is None:
                    unresolved_details[label] = reason
                    continue
                lineages[label] = record['lineage']
        else:
            for label in tree.leaf_names():
                record, reason = _get_ncbi_lineage_record(
                    label=label,
                    ncbi=ncbi,
                    rank=rank,
                    args=args,
                )
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


_RESERVED_NODE_PROPERTIES = frozenset((
    'name',
    'dist',
    'support',
    TREE_FORMAT_PROP,
))


def _internal_branch_key(node, all_taxa):
    side = frozenset(str(name) for name in node.leaf_names())
    return canonical_split(side, all_taxa - side)


def _merge_branch_annotation_values(values):
    if not values:
        return None, False
    first = values[0]
    if all(value == first for value in values[1:]):
        return first, True
    return None, False


def _snapshot_internal_branch_annotations(tree):
    all_taxa = frozenset(str(name) for name in tree.leaf_names())
    grouped = defaultdict(list)
    for node in tree.traverse():
        if node.is_root or node.is_leaf:
            continue
        custom_properties = {
            str(prop): value
            for prop, value in node.props.items()
            if str(prop) not in _RESERVED_NODE_PROPERTIES and value is not None
        }
        grouped[_internal_branch_key(node, all_taxa)].append({
            'name': node.name if node.name not in (None, '') else None,
            'support': None if support_is_missing(node.support) else node.support,
            'properties': custom_properties,
        })

    resolved = dict()
    conflicts = set()
    for split, records in grouped.items():
        annotation = {'properties': {}}
        for field in ('name', 'support'):
            values = [record[field] for record in records if record[field] is not None]
            value, usable = _merge_branch_annotation_values(values)
            if usable:
                annotation[field] = value
            elif values:
                conflicts.add((split, field))
        property_names = {
            prop
            for record in records
            for prop in record['properties']
        }
        for prop in property_names:
            values = [
                record['properties'][prop]
                for record in records
                if prop in record['properties']
            ]
            value, usable = _merge_branch_annotation_values(values)
            if usable:
                annotation['properties'][prop] = value
            elif values:
                conflicts.add((split, prop))
        resolved[split] = annotation
    return {
        'all_taxa': all_taxa,
        'by_split': resolved,
        'conflicts': conflicts,
    }


def _restore_internal_branch_annotations(tree, snapshot):
    all_taxa = snapshot['all_taxa']
    for node in tree.traverse():
        if node.is_root or node.is_leaf:
            continue
        node.name = None
        node.support = None
        for prop in list(node.props):
            if str(prop) not in _RESERVED_NODE_PROPERTIES:
                node.props.pop(prop, None)
        annotation = snapshot['by_split'].get(_internal_branch_key(node, all_taxa))
        if annotation is None:
            continue
        if 'name' in annotation:
            node.name = annotation['name']
        if 'support' in annotation:
            node.support = annotation['support']
        for prop, value in annotation['properties'].items():
            node.props[prop] = value


def _prepare_annotations_for_reroot(tree):
    snapshot = _snapshot_internal_branch_annotations(tree)
    root_support = tree.support
    leaf_support = {
        str(leaf.name): leaf.support
        for leaf in tree.leaves()
    }
    for node in tree.traverse():
        node.support = None
    return snapshot, root_support, leaf_support


def _finish_annotations_after_reroot(tree, snapshot, root_support, leaf_support,
                                     verbose=False):
    _restore_internal_branch_annotations(tree, snapshot)
    tree.support = root_support
    for leaf in tree.leaves():
        leaf.support = leaf_support[str(leaf.name)]
    if verbose and snapshot['conflicts']:
        sys.stderr.write(
            'Dropped {} conflicting root-edge annotation(s) while rerooting.\n'.format(
                len(snapshot['conflicts'])
            )
        )


def transfer_root(tree_to, tree_from, verbose=False, redistribute_root_length=True):
    tree_to = tree_to.copy(method='deepcopy')
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
    annotation_backup = None
    if not is_root_bipartition_already_matching:
        # Ensure all None dists are 0 and clear support before rerooting.
        for node in tree_to.traverse():
            if node.dist is None:
                node.dist = 0.0
        annotation_backup = _prepare_annotations_for_reroot(tree_to)
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
    if redistribute_root_length and abs(total_subroot_length_from) > 10**-15:
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
    if annotation_backup is not None:
        _finish_annotations_after_reroot(
            tree_to,
            *annotation_backup,
            verbose=verbose,
        )
    return tree_to


def transfer_root_with_taxon_mode(tree_to, tree_from, taxon_mode='exact', verbose=False,
                                  redistribute_root_length=True):
    if taxon_mode == 'exact':
        return transfer_root(
            tree_to=tree_to,
            tree_from=tree_from,
            verbose=verbose,
            redistribute_root_length=redistribute_root_length,
        )
    if taxon_mode != 'intersection':
        raise ValueError("Unsupported taxon mode for root transfer: {}".format(taxon_mode))
    tree_to = tree_to.copy(method='deepcopy')
    tree_to = _collapse_singleton_root(tree_to)
    tree_from = _collapse_singleton_root(tree_from)
    validate_unique_named_leaves(tree_to, option_name='--infile', context=' for root transfer')
    validate_unique_named_leaves(tree_from, option_name='--infile2', context=' for root transfer')
    shared_taxa, source_split, candidates = find_root_split_candidates(
        target=tree_to,
        source=tree_from,
        taxon_mode=taxon_mode,
    )
    if verbose:
        sys.stderr.write('Shared tips used for root transfer: {}\n'.format(len(shared_taxa)))
    if projected_root_split(tree_to, shared_taxa) == source_split:
        return tree_to
    unique_candidates = list()
    seen_candidate_ids = set()
    for candidate in candidates:
        if id(candidate) in seen_candidate_ids:
            continue
        seen_candidate_ids.add(id(candidate))
        unique_candidates.append(candidate)
    if len(unique_candidates) == 0:
        raise ValueError('No root bipartition matching the shared tips was found in --infile.')
    if len(unique_candidates) > 1:
        raise ValueError(
            'Root transfer is ambiguous after projecting onto shared tips ({} candidate edges).'.format(
                len(unique_candidates)
            )
        )
    candidate = unique_candidates[0]
    original_root_name = tree_to.name
    annotation_backup = _prepare_annotations_for_reroot(tree_to)
    for node in tree_to.traverse():
        if node.dist is None:
            node.dist = 0.0
    _normalize_root_distance_for_reroot(tree_to)
    tree_to.set_outgroup(candidate)
    _normalize_root_distance_for_reroot(tree_to)
    tree_to.name = original_root_name if original_root_name else 'Root'
    _finish_annotations_after_reroot(
        tree_to,
        *annotation_backup,
        verbose=verbose,
    )
    if projected_root_split(tree_to, shared_taxa) != source_split:
        raise ValueError('Root transfer failed to reproduce the source split on shared tips.')
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

def _rename_mad_leaf_names(tree):
    placeholder_to_name = dict()
    for index, leaf in enumerate(tree.leaves()):
        if leaf.name in (None, ''):
            continue
        placeholder = 'NWKITMADLEAF{:010d}'.format(index)
        placeholder_to_name[placeholder] = leaf.name
        leaf.name = placeholder
    return placeholder_to_name

def _restore_mad_leaf_names(tree, placeholder_to_name):
    if not placeholder_to_name:
        return tree
    for leaf in tree.leaves():
        restored_name = placeholder_to_name.get(leaf.name)
        if restored_name is not None:
            leaf.name = restored_name
    return tree

def mad_rooting(tree):
    """MAD (Minimal Ancestor Deviation) rooting. Tria et al. 2017, DOI:10.1038/s41559-017-0193"""
    if len(list(tree.leaves())) < 3:
        raise ValueError('MAD rooting requires at least 3 leaves.')
    import os
    import subprocess
    import tempfile
    mad_script = os.path.join(os.path.dirname(__file__), '_mad.py')
    parser = make_parser(5, dist='%0.8f')
    working_tree = tree.copy(method='deepcopy')
    placeholder_to_name = _rename_mad_leaf_names(working_tree)
    with tempfile.NamedTemporaryFile(suffix='.nwk', mode='w', delete=False) as f:
        f.write(working_tree.write(parser=parser))
        tmpfile = f.name
    try:
        result = subprocess.run(
            [sys.executable, mad_script, tmpfile],
            check=True, capture_output=True, text=True,
        )
        sys.stderr.write(result.stdout)
        with open(tmpfile + '.rooted') as fh:
            lines = [line.strip() for line in fh if line.strip() and line.strip().endswith(';')]
        if not lines:
            raise RuntimeError('MAD did not produce a rooted tree.')
        rooted_nwk = lines[0]
        rooted_tree = Tree(rooted_nwk, parser=1)
        rooted_tree = _restore_mad_leaf_names(rooted_tree, placeholder_to_name)
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
    missing_leaf_names = sorted((str(name) for name in (tree_leaf_set - taxid_leaf_set)))
    extra_leaf_names = sorted((str(name) for name in (taxid_leaf_set - tree_leaf_set)))
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

def _root_by_partial_outgroup_set(tree, resolved_outgroup_set, resolved_leaf_set, verbose=False, source_name='taxonomy'):
    full_outgroup_set = _resolve_full_outgroup_set_from_resolved_split(
        tree=tree,
        resolved_outgroup_set=resolved_outgroup_set,
        resolved_leaf_set=resolved_leaf_set,
        source_name=source_name,
    )
    return _root_by_outgroup_set(tree=tree, outgroup_set=full_outgroup_set, verbose=verbose)

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

def _expand_species_label_set(species_name_set, species_to_leaf_labels):
    leaf_label_set = set()
    for species_name in species_name_set:
        if species_name not in species_to_leaf_labels:
            raise ValueError('Unexpected species label returned by the reference tree: {}'.format(species_name))
        leaf_label_set.update(species_to_leaf_labels[species_name])
    return leaf_label_set

def _build_taxonomy_query_label_mapping(species_label_to_taxonomy_query):
    query_names = list()
    query_label_to_species_labels = defaultdict(list)
    for species_label, taxonomy_query in species_label_to_taxonomy_query.items():
        query_name = str(taxonomy_query).strip()
        query_label = query_name.replace(' ', '_')
        if query_label == '':
            raise ValueError('Taxonomy query must not be empty: {}'.format(species_label))
        if query_label not in query_label_to_species_labels:
            query_names.append(query_name)
        query_label_to_species_labels[query_label].append(species_label)
    return query_names, dict(query_label_to_species_labels)

def _expand_query_label_set(query_label_set, query_label_to_species_labels, species_to_leaf_labels):
    species_label_set = set()
    for query_label in query_label_set:
        if query_label not in query_label_to_species_labels:
            raise ValueError('Unexpected taxonomy query label returned by the reference tree: {}'.format(query_label))
        species_label_set.update(query_label_to_species_labels[query_label])
    return _expand_species_label_set(species_label_set, species_to_leaf_labels)

def _resolve_reference_query_sets(reference_tree, query_label_to_species_labels, species_to_leaf_labels, source_name):
    reference_query_labels = set(reference_tree.leaf_names())
    expected_query_labels = set(query_label_to_species_labels.keys())
    if not reference_query_labels <= expected_query_labels:
        raise ValueError(
            'Leaf labels in the {}-derived tree should match the taxonomy queries derived from --infile.'.format(
                source_name
            )
        )
    if len(reference_query_labels) == 0:
        raise ValueError('Failed to resolve usable {} taxon for any species derived from --infile.'.format(source_name))
    if (len(expected_query_labels) > 1) and (len(reference_query_labels) < 2):
        raise ValueError(
            'At least two usable {}-resolved species are required after excluding unresolved species.'.format(
                source_name
            )
        )
    resolved_species_names = set()
    for query_label in reference_query_labels:
        resolved_species_names.update(query_label_to_species_labels[query_label])
    expected_species_names = set(species_to_leaf_labels.keys())
    resolved_leaf_set = _expand_species_label_set(resolved_species_names, species_to_leaf_labels)
    unresolved_species_names = expected_species_names - resolved_species_names
    unresolved_leaf_set = _expand_species_label_set(unresolved_species_names, species_to_leaf_labels)
    return resolved_leaf_set, unresolved_leaf_set

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

def _get_timetree_name_mapping(tree, args=None):
    _, species_to_leaf_labels, species_label_to_taxonomy_query = get_species_group_records(
        tree,
        option_name='--infile',
        context=' for TimeTree taxonomy rooting',
        args=args,
    )
    query_names, query_label_to_species_labels = _build_taxonomy_query_label_mapping(species_label_to_taxonomy_query)
    return query_names, query_label_to_species_labels, species_to_leaf_labels

def _build_timetree_reference_tree(tree, args=None):
    query_names, query_label_to_species_labels, species_to_leaf_labels = _get_timetree_name_mapping(tree, args=args)
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
        resolved_leaf_set, unresolved_leaf_set = _resolve_reference_query_sets(
            timetree_tree,
            query_label_to_species_labels,
            species_to_leaf_labels,
            source_name='TimeTree',
        )
        query_outgroup_set = _get_reference_root_outgroup_set(timetree_tree, source_name='TimeTree')
        if len(query_outgroup_set) == 0:
            raise ValueError(
                'TimeTree-derived root is ambiguous: failed to identify a root bipartition.'
            )
        return (
            _expand_query_label_set(query_outgroup_set, query_label_to_species_labels, species_to_leaf_labels),
            resolved_leaf_set,
            unresolved_leaf_set,
        )
    finally:
        close = getattr(session, 'close', None)
        if callable(close):
            close()

def _get_opentree_name_mapping(tree, args=None):
    _, species_to_leaf_labels, species_label_to_taxonomy_query = get_species_group_records(
        tree,
        option_name='--infile',
        context=' for OpenTree taxonomy rooting',
        args=args,
    )
    query_names, query_label_to_species_labels = _build_taxonomy_query_label_mapping(species_label_to_taxonomy_query)
    return query_names, query_label_to_species_labels, species_to_leaf_labels

def _extract_opentree_ott_ids(tree, args=None):
    query_names, query_label_to_species_labels, species_to_leaf_labels = _get_opentree_name_mapping(tree, args=args)
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
        unresolved_labels = list()
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
                unresolved_labels.append(query_name.replace(' ', '_'))
                continue
            if len(valid_matches) != 1:
                unresolved_labels.append(query_name.replace(' ', '_'))
                continue
            ott_ids.append(int(valid_matches[0]['taxon']['ott_id']))
        if len(ott_ids) == 0:
            raise ValueError(
                'Failed to resolve a usable OpenTree taxon for any leaf label: {}'.format(
                    '; '.join(unresolved_labels)
                )
            )
        return ott_ids, query_label_to_species_labels, species_to_leaf_labels
    finally:
        close = getattr(session, 'close', None)
        if callable(close):
            close()

def _build_opentree_reference_tree(tree, args=None):
    ott_ids, query_label_to_species_labels, species_to_leaf_labels = _extract_opentree_ott_ids(tree, args=args)
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
        resolved_leaf_set, unresolved_leaf_set = _resolve_reference_query_sets(
            opentree_tree,
            query_label_to_species_labels,
            species_to_leaf_labels,
            source_name='OpenTree',
        )
        query_outgroup_set = _get_reference_root_outgroup_set(opentree_tree, source_name='OpenTree')
        if len(query_outgroup_set) == 0:
            raise ValueError(
                'OpenTree-derived root is ambiguous: failed to identify a root bipartition.'
            )
        return (
            _expand_query_label_set(query_outgroup_set, query_label_to_species_labels, species_to_leaf_labels),
            resolved_leaf_set,
            unresolved_leaf_set,
        )
    finally:
        close = getattr(session, 'close', None)
        if callable(close):
            close()

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
            if source == 'timetree':
                outgroup_set, resolved_leaf_set, unresolved_leaf_set = _build_timetree_reference_tree(tree=tree, args=args)
            elif source == 'opentree':
                outgroup_set, resolved_leaf_set, unresolved_leaf_set = _build_opentree_reference_tree(tree=tree, args=args)
            else:
                raise ValueError("Unknown taxonomy source: {}".format(source))
            if verbose and unresolved_leaf_set:
                sys.stderr.write(
                    'Excluding {}-unresolved leaf label(s) from taxonomy rooting: {}\n'.format(
                        source,
                        '; '.join(sorted(unresolved_leaf_set)),
                    )
                )
            if len(unresolved_leaf_set) == 0:
                return _root_by_outgroup_set(tree=tree, outgroup_set=outgroup_set, verbose=verbose)
            return _root_by_partial_outgroup_set(
                tree=tree,
                resolved_outgroup_set=outgroup_set,
                resolved_leaf_set=resolved_leaf_set,
                verbose=verbose,
                source_name=source,
            )
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
        taxon_mode = getattr(args, 'taxon_mode', 'exact')
        if taxon_mode == 'exact' and not is_all_leaf_names_identical(tree, tree2, verbose=True):
            raise ValueError('Leaf labels must match exactly when --taxon-mode exact.')
        if (len(list(tree.leaves())) > 1) and (len(list(tree2.leaves())) > 1):
            tree = transfer_root_with_taxon_mode(
                tree_to=tree,
                tree_from=tree2,
                taxon_mode=taxon_mode,
                verbose=True,
            )
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
