import json
import math
import re
import sys
from itertools import combinations
from collections import defaultdict
from fractions import Fraction

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from ete4 import Tree

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
    get_subtree_leaf_name_sets,
    get_tree_property_names,
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
TAXONOMY_HTTP_MAX_BYTES = 10 * 1024 * 1024
TAXONOMY_HTTP_CHUNK_BYTES = 64 * 1024
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
MAD_FLOAT_MAX_PATH_DYNAMIC_RANGE = 10 ** 6


def _new_taxonomy_http_session():
    session = requests.Session()
    mount = getattr(session, 'mount', None)
    if callable(mount):
        retry = Retry(
            total=3,
            connect=3,
            read=3,
            status=3,
            backoff_factor=0.5,
            status_forcelist=(429, 500, 502, 503, 504),
            allowed_methods=frozenset(('GET', 'POST')),
        )
        adapter = HTTPAdapter(max_retries=retry)
        mount('https://', adapter)
        mount('http://', adapter)
    return session


def _close_http_response(response):
    close = getattr(response, 'close', None)
    if callable(close):
        close()


def _bounded_response_text(response, source_name):
    headers = getattr(response, 'headers', {}) or {}
    content_length = headers.get('Content-Length')
    if content_length is not None:
        try:
            declared_size = int(content_length)
        except (TypeError, ValueError):
            declared_size = None
        if declared_size is not None and declared_size > TAXONOMY_HTTP_MAX_BYTES:
            _close_http_response(response)
            raise ValueError('{} response exceeds the size limit.'.format(source_name))

    iter_content = getattr(response, 'iter_content', None)
    try:
        if callable(iter_content):
            chunks = list()
            total_size = 0
            try:
                content_iterator = iter_content(chunk_size=TAXONOMY_HTTP_CHUNK_BYTES)
                for chunk in content_iterator:
                    if not chunk:
                        continue
                    if isinstance(chunk, str):
                        chunk = chunk.encode('utf-8')
                    else:
                        chunk = bytes(chunk)
                    total_size += len(chunk)
                    if total_size > TAXONOMY_HTTP_MAX_BYTES:
                        raise ValueError(
                            '{} response exceeds the size limit.'.format(source_name)
                        )
                    chunks.append(chunk)
            except requests.RequestException as exc:
                raise ValueError(
                    'Failed while reading the {} response.'.format(source_name)
                ) from exc
            payload = b''.join(chunks)
            encoding = getattr(response, 'encoding', None) or 'utf-8'
            try:
                return payload.decode(encoding)
            except (LookupError, UnicodeError) as exc:
                raise ValueError(
                    'Unexpected text encoding in the {} response.'.format(source_name)
                ) from exc

        text = str(getattr(response, 'text', ''))
        if len(text.encode('utf-8')) > TAXONOMY_HTTP_MAX_BYTES:
            raise ValueError('{} response exceeds the size limit.'.format(source_name))
        return text
    finally:
        _close_http_response(response)


def _response_json_object(response, source_name):
    text = _bounded_response_text(response, source_name)
    try:
        data = json.loads(text)
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        raise ValueError('Unexpected JSON response from {}.'.format(source_name)) from exc
    if not isinstance(data, dict):
        raise ValueError('Unexpected JSON response from {}.'.format(source_name))
    return data

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
        root_dist = tree.dist
        child_dist = child.dist
        tree = child.copy(method='deepcopy')
        if root_dist is None or child_dist is None:
            tree.dist = None
        else:
            collapsed_root_dist = float(root_dist) + float(child_dist)
            if not math.isfinite(collapsed_root_dist):
                raise ValueError(
                    'Collapsing a singleton root would produce a '
                    'non-finite root branch length.'
                )
            tree.dist = collapsed_root_dist
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


def _snapshot_missing_branch_lengths(tree):
    all_taxa = frozenset(str(name) for name in tree.leaf_names())
    missing_splits = {
        _internal_branch_key(node, all_taxa)
        for node in tree.traverse()
        if (not node.is_root) and node.dist is None
    }
    root_split = None
    missing_root_sides = set()
    root_children = tree.get_children()
    if len(root_children) == 2:
        root_sides = [
            frozenset(str(name) for name in child.leaf_names())
            for child in root_children
        ]
        if root_sides[0] and root_sides[1]:
            root_split = canonical_split(root_sides[0], root_sides[1])
            missing_root_sides = {
                side
                for child, side in zip(root_children, root_sides)
                if child.dist is None
            }
    return {
        'splits': missing_splits,
        'root_split': root_split,
        'missing_root_sides': missing_root_sides,
    }


def _restore_missing_branch_lengths(tree, snapshot):
    missing_splits = snapshot['splits']
    if not missing_splits:
        return
    all_taxa = frozenset(str(name) for name in tree.leaf_names())
    same_root_split = False
    root_children = tree.get_children()
    if snapshot['root_split'] is not None and len(root_children) == 2:
        final_root_sides = [
            frozenset(str(name) for name in child.leaf_names())
            for child in root_children
        ]
        same_root_split = (
            canonical_split(final_root_sides[0], final_root_sides[1])
            == snapshot['root_split']
        )
    for node in tree.traverse():
        if node.is_root:
            continue
        split = _internal_branch_key(node, all_taxa)
        if (
            same_root_split
            and node.up is tree
            and split == snapshot['root_split']
        ):
            side = frozenset(str(name) for name in node.leaf_names())
            if side in snapshot['missing_root_sides']:
                node.dist = None
            continue
        if split in missing_splits:
            node.dist = None


def _first_leaf_name(tree, names):
    names = set(names)
    for leaf in tree.leaves():
        if leaf.name in names:
            return leaf.name
    raise ValueError('No requested leaf was found in the tree.')


def _custom_node_properties(node):
    return {
        str(prop): value
        for prop, value in node.props.items()
        if str(prop) not in _RESERVED_NODE_PROPERTIES and value is not None
    }


def _snapshot_node_annotation(node, include_name=True):
    annotation = {
        'support': node.support,
        'properties': _custom_node_properties(node),
    }
    if include_name:
        annotation['name'] = node.name
    return annotation


def _restore_node_annotation(node, annotation, include_name=True):
    if include_name:
        node.name = annotation['name']
    node.support = annotation['support']
    for prop in list(node.props):
        if str(prop) not in _RESERVED_NODE_PROPERTIES:
            node.props.pop(prop, None)
    for prop, value in annotation['properties'].items():
        node.props[prop] = value


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
        custom_properties = _custom_node_properties(node)
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
    snapshot = {
        'internal': _snapshot_internal_branch_annotations(tree),
        'root': _snapshot_node_annotation(tree),
        'root_dist': tree.dist,
        'missing_branch_lengths': _snapshot_missing_branch_lengths(tree),
        'leaves_by_identity': {
            id(leaf): _snapshot_node_annotation(leaf, include_name=False)
            for leaf in tree.leaves()
        },
        'leaves_by_name': defaultdict(list),
    }
    for leaf in tree.leaves():
        snapshot['leaves_by_name'][str(leaf.name)].append(
            _snapshot_node_annotation(leaf, include_name=False)
        )
    for node in tree.traverse():
        node.support = None
    return snapshot


def _finish_annotations_after_reroot(tree, snapshot, verbose=False):
    _restore_internal_branch_annotations(tree, snapshot['internal'])
    _restore_node_annotation(tree, snapshot['root'])
    leaf_name_offsets = defaultdict(int)
    for leaf in tree.leaves():
        annotation = snapshot['leaves_by_identity'].get(id(leaf))
        if annotation is None:
            leaf_name = str(leaf.name)
            offset = leaf_name_offsets[leaf_name]
            candidates = snapshot['leaves_by_name'].get(leaf_name, ())
            if offset >= len(candidates):
                continue
            annotation = candidates[offset]
            leaf_name_offsets[leaf_name] += 1
        _restore_node_annotation(leaf, annotation, include_name=False)
    _restore_missing_branch_lengths(tree, snapshot['missing_branch_lengths'])
    tree.dist = snapshot['root_dist']
    if verbose and snapshot['internal']['conflicts']:
        sys.stderr.write(
            'Dropped {} conflicting root-edge annotation(s) while rerooting.\n'.format(
                len(snapshot['internal']['conflicts'])
            )
        )


def _redistribute_root_child_lengths(target, source, shared_taxa):
    target_children = target.get_children()
    source_children = source.get_children()
    if len(target_children) != 2 or len(source_children) != 2:
        return
    # Missing, non-finite, or negative halves do not define a usable physical
    # root-edge total or ratio. Leave the target lengths unchanged.
    if (
        any(child.dist is None for child in target_children)
        or any(child.dist is None for child in source_children)
    ):
        return
    try:
        target_lengths = [float(child.dist) for child in target_children]
        source_lengths = [float(child.dist) for child in source_children]
    except (TypeError, ValueError):
        return
    if any(
        (not math.isfinite(length)) or length < 0.0
        for length in target_lengths + source_lengths
    ):
        return

    shared_taxa = frozenset(str(taxon) for taxon in shared_taxa)
    source_by_side = dict()
    for child in source_children:
        side = frozenset(
            str(name)
            for name in child.leaf_names()
            if str(name) in shared_taxa
        )
        if not side or side == shared_taxa or side in source_by_side:
            return
        source_by_side[side] = child

    target_by_side = dict()
    for child in target_children:
        side = frozenset(
            str(name)
            for name in child.leaf_names()
            if str(name) in shared_taxa
        )
        if not side or side == shared_taxa or side in target_by_side:
            return
        target_by_side[side] = child
    if set(target_by_side) != set(source_by_side):
        return

    source_length_by_id = {
        id(child): length
        for child, length in zip(source_children, source_lengths)
    }
    source_total = sum(source_lengths)
    target_total = sum(target_lengths)
    if (
        source_total == 0.0
        or not math.isfinite(source_total)
        or not math.isfinite(target_total)
    ):
        return
    for side, target_child in target_by_side.items():
        source_child = source_by_side[side]
        target_child.dist = (
            target_total
            * (source_length_by_id[id(source_child)] / source_total)
        )


def transfer_root(tree_to, tree_from, verbose=False, redistribute_root_length=True):
    tree_to = tree_to.copy(method='deepcopy')
    tree_to = _collapse_singleton_root(tree_to)
    tree_from = _collapse_singleton_root(tree_from)
    validate_unique_named_leaves(tree_to, option_name='--infile', context=' for root transfer')
    validate_unique_named_leaves(tree_from, option_name='--infile2', context=' for root transfer')
    if not is_all_leaf_names_identical(tree_to, tree_from):
        raise ValueError(
            "Exact root transfer requires identical leaf labels in '--infile' and '--infile2'."
        )
    subroot_from = tree_from.get_children()
    if len(subroot_from) != 2:
        raise ValueError('Root transfer requires the source tree root to have exactly two children.')
    tree_from_leaf_sets = get_subtree_leaf_name_sets(tree_from)
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
    tree_to_leaf_sets = get_subtree_leaf_name_sets(tree_to)
    is_root_bipartition_already_matching = (
        (len(root_children) == 2) and
        any(outgroup_set == tree_to_leaf_sets[child] for child in root_children)
    )
    annotation_backup = None
    if not is_root_bipartition_already_matching:
        annotation_backup = _prepare_annotations_for_reroot(tree_to)
        # Ensure all None dists are 0 and clear support before rerooting.
        for node in tree_to.traverse():
            if node.dist is None:
                node.dist = 0.0
        _normalize_root_distance_for_reroot(tree_to)
        tree_to.set_outgroup(_first_leaf_name(tree_to, ingroup_set))
        if len(outgroup_set) == 1:
            outgroup_name = min(outgroup_set, key=str)
            outgroup_ancestor = None
            for leaf in tree_to.leaves():
                if leaf.name == outgroup_name:
                    outgroup_ancestor = leaf
                    break
            if outgroup_ancestor is None:
                raise ValueError('No root bipartition found in --infile.')
        else:
            outgroup_ancestor = tree_to.common_ancestor(sorted(outgroup_set, key=str))
        reroot_leaf_sets = get_subtree_leaf_name_sets(tree_to)
        if not outgroup_set == reroot_leaf_sets[outgroup_ancestor]:
            raise ValueError('No root bipartition found in --infile.')
        _normalize_root_distance_for_reroot(tree_to)
        tree_to.set_outgroup(outgroup_ancestor)
        tree_to_leaf_sets = get_subtree_leaf_name_sets(tree_to)
    if redistribute_root_length:
        _redistribute_root_child_lengths(
            target=tree_to,
            source=tree_from,
            shared_taxa=tree_to.leaf_names(),
        )
    # Restore root name or assign 'Root' if it was unnamed
    if original_root_name:
        tree_to.name = original_root_name
    else:
        tree_to.name = 'Root'
    if annotation_backup is not None:
        _finish_annotations_after_reroot(
            tree_to,
            annotation_backup,
            verbose=verbose,
        )
        if not original_root_name:
            tree_to.name = 'Root'
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
        if redistribute_root_length:
            _redistribute_root_child_lengths(
                target=tree_to,
                source=tree_from,
                shared_taxa=shared_taxa,
            )
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
        annotation_backup,
        verbose=verbose,
    )
    if redistribute_root_length:
        _redistribute_root_child_lengths(
            target=tree_to,
            source=tree_from,
            shared_taxa=shared_taxa,
        )
        _restore_missing_branch_lengths(
            tree_to,
            annotation_backup['missing_branch_lengths'],
        )
    if not original_root_name:
        tree_to.name = 'Root'
    if projected_root_split(tree_to, shared_taxa) != source_split:
        raise ValueError('Root transfer failed to reproduce the source split on shared tips.')
    return tree_to

def _normalize_reroot_branch_lengths(tree, method_name):
    """Validate and scale non-root branches for numerically stable rooting."""
    branch_lengths = []
    branch_length_scale = 0.0
    for node in tree.traverse():
        if node.is_root:
            continue
        if node.dist is None:
            length = 0.0
        else:
            try:
                length = float(node.dist)
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    '{} rooting requires finite, non-negative branch lengths.'.format(
                        method_name
                    )
                ) from exc
            if not math.isfinite(length) or length < 0.0:
                raise ValueError(
                    '{} rooting requires finite, non-negative branch lengths.'.format(
                        method_name
                    )
                )
        branch_lengths.append((node, length))
        if length > branch_length_scale:
            branch_length_scale = length

    if branch_length_scale == 0.0:
        branch_length_scale = 1.0
    for node, length in branch_lengths:
        node.dist = length / branch_length_scale
    return branch_length_scale


def _restore_reroot_branch_length_scale(tree, branch_length_scale):
    for node in tree.traverse():
        if node.is_root or node.dist is None:
            continue
        node.dist = float(node.dist) * branch_length_scale


def midpoint_rooting(tree):
    tree = tree.copy(method='deepcopy')
    validate_unique_named_leaves(
        tree,
        option_name='--infile',
        context=' for midpoint rooting',
    )
    annotation_backup = _prepare_annotations_for_reroot(tree)
    branch_length_scale = _normalize_reroot_branch_lengths(tree, 'Midpoint')
    _normalize_root_distance_for_reroot(tree)
    outgroup_node = tree.get_midpoint_outgroup()
    # If the outgroup is the root itself, tree is already optimally rooted
    if outgroup_node.is_root:
        _restore_reroot_branch_length_scale(tree, branch_length_scale)
        _finish_annotations_after_reroot(tree, annotation_backup)
        return tree
    tree.set_outgroup(outgroup_node)
    _restore_reroot_branch_length_scale(tree, branch_length_scale)
    _finish_annotations_after_reroot(tree, annotation_backup)
    return tree

def mad_rooting(tree):
    """Root a tree by minimal ancestor deviation (Tria et al. 2017)."""
    if len(list(tree.leaves())) < 3:
        raise ValueError('MAD rooting requires at least 3 leaves.')
    tree = tree.copy(method='deepcopy')
    validate_unique_named_leaves(tree, option_name='--infile', context=' for MAD rooting')
    positive_branch_count = 0
    branch_length_scale = 0.0
    minimum_positive_branch_length = math.inf
    for node in tree.traverse():
        if node.is_root:
            continue
        if node.dist is None:
            raise ValueError('MAD rooting requires every branch length to be present.')
        length = float(node.dist)
        if not math.isfinite(length) or length < 0:
            raise ValueError('MAD rooting requires finite, non-negative branch lengths.')
        if length > 0:
            positive_branch_count += 1
            branch_length_scale = max(branch_length_scale, length)
            minimum_positive_branch_length = min(
                minimum_positive_branch_length,
                length,
            )
    if positive_branch_count < 3:
        raise ValueError('MAD rooting requires at least 3 positive branch lengths.')

    annotation_backup = _prepare_annotations_for_reroot(tree)
    _normalize_root_distance_for_reroot(tree)
    # Work on a dimensionless copy.  Besides improving the conditioning of
    # the quadratic objective, this avoids overflowing when the two finite
    # halves of an arbitrary input root have to be joined into one physical
    # edge.
    for node in tree.traverse():
        if not node.is_root:
            original_length = float(node.dist)
            normalized_length = original_length / branch_length_scale
            if original_length > 0.0 and normalized_length == 0.0:
                raise ValueError(
                    'MAD rooting cannot preserve this branch-length dynamic '
                    'range with finite floating-point values.'
                )
            node.dist = normalized_length

    # Suppress an arbitrary bifurcating root so every candidate is a physical
    # edge of the unrooted tree. The two former root lengths form one edge.
    children = tree.get_children()
    if len(children) == 2:
        first, second = children
        dissolved = first if not first.is_leaf else (second if not second.is_leaf else None)
        if dissolved is not None:
            retained = second if dissolved is first else first
            retained.dist = float(retained.dist) + float(dissolved.dist)
            for grandchild in list(dissolved.get_children()):
                tree.add_child(grandchild)
            tree.remove_child(dissolved)

    # MAD treats zero-distance tips as one effective OTU. Choose a
    # deterministic representative while retaining every tip in the output.
    all_nodes = list(tree.traverse(strategy='preorder'))
    node_index_by_id = {
        id(node): index
        for index, node in enumerate(all_nodes)
    }
    component_parents = list(range(len(all_nodes)))

    def find(index):
        while component_parents[index] != index:
            component_parents[index] = component_parents[
                component_parents[index]
            ]
            index = component_parents[index]
        return index

    def union(first_index, second_index):
        first_root = find(first_index)
        second_root = find(second_index)
        if first_root != second_root:
            component_parents[max(first_root, second_root)] = min(
                first_root,
                second_root,
            )

    for node in all_nodes:
        if not node.is_root and float(node.dist) == 0.0:
            union(
                node_index_by_id[id(node)],
                node_index_by_id[id(node.up)],
            )
    all_tips = list(tree.leaves())
    effective_groups = defaultdict(list)
    for index, tip in enumerate(all_tips):
        component = find(node_index_by_id[id(tip)])
        effective_groups[component].append((index, tip))
    tips = [
        min(group, key=lambda item: (str(item[1].name), item[0]))[1]
        for group in effective_groups.values()
    ]
    tips.sort(key=lambda tip: str(tip.name))
    if len(tips) < 3:
        raise ValueError('MAD rooting requires at least 3 effective leaves.')

    num_tips = len(tips)
    effective_tip_index_by_id = {
        id(tip): index
        for index, tip in enumerate(tips)
    }
    all_tips_by_name = sorted(all_tips, key=lambda tip: str(tip.name))
    all_tip_index_by_id = {
        id(tip): index
        for index, tip in enumerate(all_tips_by_name)
    }
    effective_subtree_masks = dict()
    all_subtree_masks = dict()
    for node in tree.traverse(strategy='postorder'):
        if node.is_leaf:
            effective_index = effective_tip_index_by_id.get(id(node))
            effective_subtree_masks[node] = (
                0 if effective_index is None else (1 << effective_index)
            )
            all_subtree_masks[node] = 1 << all_tip_index_by_id[id(node)]
            continue
        effective_mask = 0
        all_mask = 0
        for child in node.get_children():
            effective_mask |= effective_subtree_masks[child]
            all_mask |= all_subtree_masks[child]
        effective_subtree_masks[node] = effective_mask
        all_subtree_masks[node] = all_mask

    # MAD's objective is a weighted quadratic in the root-to-tip distances.
    # Each tip pair affects exactly the edges on the path between the pair.
    # Accumulate those path contributions once, then obtain every edge's
    # fitted root position and score with two linear tree traversals.
    node_count = len(all_nodes)
    normalized_minimum_branch_length = (
        minimum_positive_branch_length / branch_length_scale
    )
    use_float_arithmetic = (
        positive_branch_count
        <= (
            MAD_FLOAT_MAX_PATH_DYNAMIC_RANGE
            * normalized_minimum_branch_length
        )
    )
    parent_indices = [0] * node_count
    topological_depths = [0] * node_count
    # Ordinary trees use fast float arithmetic. For a path-length dynamic
    # range large enough to lose short terms, dyadic rational arithmetic
    # preserves every input binary float exactly.
    zero = 0.0 if use_float_arithmetic else Fraction(0)
    root_distances_by_index = [zero] * node_count
    for node in all_nodes:
        node_index = node_index_by_id[id(node)]
        if node.is_root:
            parent_indices[node_index] = node_index
            continue
        parent_index = node_index_by_id[id(node.up)]
        parent_indices[node_index] = parent_index
        topological_depths[node_index] = topological_depths[parent_index] + 1
        edge_length = float(node.dist)
        if not use_float_arithmetic:
            edge_length = Fraction.from_float(edge_length)
        root_distances_by_index[node_index] = (
            root_distances_by_index[parent_index] + edge_length
        )

    ancestor_levels = [parent_indices]
    while (1 << len(ancestor_levels)) <= max(topological_depths):
        previous = ancestor_levels[-1]
        ancestor_levels.append([
            previous[previous[index]]
            for index in range(node_count)
        ])

    def lift_node_index(node_index, levels):
        bit = 0
        while levels:
            if levels & 1:
                node_index = ancestor_levels[bit][node_index]
            levels >>= 1
            bit += 1
        return node_index

    def lowest_common_ancestor_index(first_index, second_index):
        first_depth = topological_depths[first_index]
        second_depth = topological_depths[second_index]
        if first_depth > second_depth:
            first_index = lift_node_index(
                first_index,
                first_depth - second_depth,
            )
        elif second_depth > first_depth:
            second_index = lift_node_index(
                second_index,
                second_depth - first_depth,
            )
        if first_index == second_index:
            return first_index
        for ancestors in reversed(ancestor_levels):
            if ancestors[first_index] != ancestors[second_index]:
                first_index = ancestors[first_index]
                second_index = ancestors[second_index]
        return parent_indices[first_index]

    effective_node_indices = [
        node_index_by_id[id(tip)]
        for tip in tips
    ]

    def iter_pair_records():
        for first_index, second_index in combinations(range(num_tips), 2):
            first_node_index = effective_node_indices[first_index]
            second_node_index = effective_node_indices[second_index]
            ancestor_index = lowest_common_ancestor_index(
                first_node_index,
                second_node_index,
            )
            pair_distance = (
                root_distances_by_index[first_node_index]
                + root_distances_by_index[second_node_index]
                - (2 * root_distances_by_index[ancestor_index])
            )
            if pair_distance <= 0:
                raise ValueError(
                    'MAD rooting requires at least 3 effective leaves.'
                )
            yield (
                first_node_index,
                second_node_index,
                ancestor_index,
                pair_distance,
                float(pair_distance),
            )

    minimum_pair_distance = math.inf
    maximum_pair_distance = 0.0
    for _, _, _, _, pair_distance in iter_pair_records():
        minimum_pair_distance = min(minimum_pair_distance, pair_distance)
        maximum_pair_distance = max(maximum_pair_distance, pair_distance)
    # The weight scale is algebraically arbitrary. A geometric midpoint keeps
    # both ends representable; the extreme-range fallback squares it exactly.
    weight_distance_scale = (
        math.sqrt(minimum_pair_distance)
        * math.sqrt(maximum_pair_distance)
    )
    weight_distance_scale_number = (
        weight_distance_scale
        if use_float_arithmetic
        else Fraction.from_float(weight_distance_scale)
    )

    cut_weights = [zero] * node_count
    numerator_constants = [zero] * node_count
    numerator_depths = [zero] * node_count
    root_objective = zero
    for (
        first_node_index,
        second_node_index,
        ancestor_index,
        pair_distance_fraction,
        pair_distance,
    ) in iter_pair_records():
        weight_ratio = weight_distance_scale / pair_distance
        weight_number = (
            weight_ratio
            if use_float_arithmetic
            else Fraction.from_float(weight_ratio)
        )
        weight = weight_number * weight_number

        cut_weights[first_node_index] += weight
        cut_weights[second_node_index] += weight
        cut_weights[ancestor_index] -= 2 * weight

        first_constant = weight * (
            pair_distance_fraction
            - (2 * root_distances_by_index[first_node_index])
        )
        second_constant = weight * (
            pair_distance_fraction
            - (2 * root_distances_by_index[second_node_index])
        )
        numerator_constants[first_node_index] += first_constant
        numerator_constants[second_node_index] += second_constant
        numerator_constants[ancestor_index] -= (
            first_constant + second_constant
        )

        depth_coefficient = 2 * weight
        numerator_depths[first_node_index] += depth_coefficient
        numerator_depths[second_node_index] += depth_coefficient
        numerator_depths[ancestor_index] -= 2 * depth_coefficient

        root_distance_difference = (
            root_distances_by_index[first_node_index]
            - root_distances_by_index[second_node_index]
        )
        root_objective += (
            weight * root_distance_difference * root_distance_difference
        )

    for node in tree.traverse(strategy='postorder'):
        if node.is_root:
            continue
        node_index = node_index_by_id[id(node)]
        parent_index = parent_indices[node_index]
        cut_weights[parent_index] += cut_weights[node_index]
        numerator_constants[parent_index] += numerator_constants[node_index]
        numerator_depths[parent_index] += numerator_depths[node_index]

    objectives_by_index = [zero] * node_count
    root_index = node_index_by_id[id(tree)]
    objectives_by_index[root_index] = root_objective
    numerators_by_index = [zero] * node_count
    for node in tree.traverse(strategy='preorder'):
        if node.is_root:
            continue
        node_index = node_index_by_id[id(node)]
        parent_index = parent_indices[node_index]
        edge_length = float(node.dist)
        if not use_float_arithmetic:
            edge_length = Fraction.from_float(edge_length)
        numerator = (
            numerator_constants[node_index]
            + (
                numerator_depths[node_index]
                * root_distances_by_index[node_index]
            )
        )
        numerators_by_index[node_index] = numerator
        objectives_by_index[node_index] = (
            objectives_by_index[parent_index]
            + (4 * edge_length * numerator)
            - (4 * edge_length * edge_length * cut_weights[node_index])
        )

    best_score = math.inf
    best_node = None
    best_root_distance = 0.0
    best_edge_length = 0.0
    best_side_names = None
    best_candidate_key = None
    pair_count = num_tips * (num_tips - 1) // 2
    all_leaf_mask = (1 << len(all_tips_by_name)) - 1

    for node in tree.traverse():
        if node.is_root:
            continue
        edge_length = float(node.dist)
        if not use_float_arithmetic:
            edge_length = Fraction.from_float(edge_length)
        if edge_length == 0:
            continue
        side_mask = effective_subtree_masks[node]
        if side_mask == 0 or side_mask == ((1 << num_tips) - 1):
            continue

        node_index = node_index_by_id[id(node)]
        cut_weight = cut_weights[node_index]
        numerator = numerators_by_index[node_index]
        denominator = 2 * cut_weight
        root_distance = (
            zero
            if denominator == 0
            else numerator / denominator
        )
        root_distance = min(max(root_distance, zero), edge_length)

        objective = (
            objectives_by_index[node_index]
            - (4 * root_distance * numerator)
            + (4 * root_distance * root_distance * cut_weight)
        )
        normalized_objective = (
            objective
            / (
                pair_count
                * weight_distance_scale_number
                * weight_distance_scale_number
            )
        )
        score = math.sqrt(max(0.0, float(normalized_objective)))
        is_tied = math.isclose(
            score,
            best_score,
            rel_tol=10 ** -12,
            abs_tol=10 ** -15,
        )
        is_improvement = score < best_score and not is_tied
        if not is_improvement and not is_tied:
            continue
        side_name_mask = all_subtree_masks[node]
        other_name_mask = all_leaf_mask ^ side_name_mask
        side_name_indices = tuple(
            index
            for index in range(len(all_tips_by_name))
            if side_name_mask & (1 << index)
        )
        other_name_indices = tuple(
            index
            for index in range(len(all_tips_by_name))
            if other_name_mask & (1 << index)
        )
        candidate_key = (
            (side_name_indices, other_name_indices)
            if (len(side_name_indices), side_name_indices)
            <= (len(other_name_indices), other_name_indices)
            else (other_name_indices, side_name_indices)
        )
        if (
            is_improvement
            or (
                is_tied
                and (best_candidate_key is None or candidate_key < best_candidate_key)
            )
        ):
            best_score = score
            best_node = node
            best_root_distance = float(root_distance)
            best_edge_length = float(edge_length)
            best_side_names = frozenset(
                str(all_tips_by_name[index].name)
                for index in side_name_indices
            )
            best_candidate_key = candidate_key

    if best_node is None:
        raise ValueError('MAD rooting could not identify a candidate root edge.')

    sys.stderr.write('MAD rooting score: {:.6g}\n'.format(best_score))
    tree.set_outgroup(best_node)
    for child in tree.get_children():
        side_names = frozenset(str(name) for name in child.leaf_names())
        if side_names == best_side_names:
            child.dist = best_root_distance
        else:
            child.dist = best_edge_length - best_root_distance
    for node in tree.traverse():
        if node.is_root:
            continue
        restored_length = float(node.dist) * branch_length_scale
        if not math.isfinite(restored_length):
            raise ValueError(
                'MAD rooting produced a branch length that exceeds the '
                'finite floating-point range.'
            )
        node.dist = restored_length
    _finish_annotations_after_reroot(tree, annotation_backup)
    return tree

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
    tree = tree.copy(method='deepcopy')
    validate_unique_named_leaves(
        tree,
        option_name='--infile',
        context=' for MV rooting',
    )
    annotation_backup = _prepare_annotations_for_reroot(tree)
    branch_length_scale = _normalize_reroot_branch_lengths(tree, 'MV')
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
        _restore_reroot_branch_length_scale(tree, branch_length_scale)
        _finish_annotations_after_reroot(tree, annotation_backup)
        return tree
    reported_var = best_var * branch_length_scale
    reported_var *= branch_length_scale
    sys.stderr.write('MV rooting variance: {:.6g}\n'.format(reported_var))
    tree.set_outgroup(best_node)
    # Adjust branch lengths at root: best_x from root to best_node, (L - best_x) to sibling
    best_subtree_leaves = set(best_node.leaf_names())
    root_leaf_sets = get_subtree_leaf_name_sets(tree)
    for child in tree.get_children():
        if root_leaf_sets[child] == best_subtree_leaves:
            child.dist = best_x
        else:
            child.dist = best_L - best_x
    _restore_reroot_branch_length_scale(tree, branch_length_scale)
    _finish_annotations_after_reroot(tree, annotation_backup)
    return tree

def outgroup_rooting(tree, outgroup_str):
    if outgroup_str is None:
        raise ValueError("Specify at least one outgroup label with '--outgroup'.")
    tree = tree.copy(method='deepcopy')
    validate_unique_named_leaves(
        tree,
        option_name='--infile',
        context=' for outgroup rooting',
    )
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
    annotation_backup = _prepare_annotations_for_reroot(tree)
    # ETE cannot reroot through missing distances; restore their missingness
    # from the split-keyed snapshot after the topology operation.
    for node in tree.traverse():
        if node.dist is None:
            node.dist = 0.0
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
    _finish_annotations_after_reroot(tree, annotation_backup)
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
        raise ValueError('--taxid-tsv must match the leaf labels in --infile exactly ({})'.format('; '.join(messages)))
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
    subtree_leaf_sets = get_subtree_leaf_name_sets(analysis_tree)
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
    child_leaf_sets = get_subtree_leaf_name_sets(reference_tree)
    root_children = reference_tree.get_children()
    is_n0_bigger_than_n1 = (len(child_leaf_sets[root_children[0]]) > len(child_leaf_sets[root_children[1]]))
    outgroup_child = root_children[1] if is_n0_bigger_than_n1 else root_children[0]
    return set(child_leaf_sets[outgroup_child])

def _expand_species_label_set(species_name_set, species_to_leaf_labels):
    leaf_label_set = set()
    for species_name in sorted(species_name_set, key=str):
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
    for query_label in sorted(query_label_set, key=str):
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
    tree = tree.copy(method='deepcopy')
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
    tree_leaf_sets = get_subtree_leaf_name_sets(tree)
    is_root_bipartition_already_matching = (
        (len(root_children) == 2) and
        any(outgroup_set == tree_leaf_sets[child] for child in root_children)
    )
    annotation_backup = None
    if not is_root_bipartition_already_matching:
        ingroup_set = leaf_name_set - outgroup_set
        annotation_backup = _prepare_annotations_for_reroot(tree)
        for node in tree.traverse():
            if node.dist is None:
                node.dist = 0.0
        _normalize_root_distance_for_reroot(tree)
        tree.set_outgroup(_first_leaf_name(tree, ingroup_set))
        if len(outgroup_set) == 1:
            outgroup_name = min(outgroup_set, key=str)
            outgroup_ancestor = None
            for leaf in tree.leaves():
                if leaf.name == outgroup_name:
                    outgroup_ancestor = leaf
                    break
            if outgroup_ancestor is None:
                raise ValueError('No root bipartition found in --infile.')
        else:
            outgroup_ancestor = tree.common_ancestor(sorted(outgroup_set, key=str))
        reroot_leaf_sets = get_subtree_leaf_name_sets(tree)
        if outgroup_set != reroot_leaf_sets[outgroup_ancestor]:
            raise ValueError('No root bipartition found in --infile.')
        _normalize_root_distance_for_reroot(tree)
        tree.set_outgroup(outgroup_ancestor)
    if original_root_name:
        tree.name = original_root_name
    else:
        tree.name = 'Root'
    if annotation_backup is not None:
        _finish_annotations_after_reroot(
            tree,
            annotation_backup,
            verbose=verbose,
        )
        if not original_root_name:
            tree.name = 'Root'
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
    session = _new_taxonomy_http_session()
    try:
        try:
            upload_response = session.post(
                'https://timetree.org/ajax/prune/load_names/',
                files={'file': ('species.txt', '\n'.join(query_names) + '\n')},
                timeout=60,
                stream=True,
            )
            _bounded_response_text(upload_response, 'TimeTree upload')
            if upload_response.status_code != 200:
                raise ValueError('Failed to retrieve a pruned guide tree from TimeTree.')
            newick_response = session.post(
                'https://timetree.org/ajax/newick/prunetree/download',
                data={'export': 'newick', 'rank': ''},
                timeout=60,
                stream=True,
            )
        except requests.RequestException as exc:
            raise ValueError('Failed to contact TimeTree.') from exc
        if newick_response.status_code != 200:
            _close_http_response(newick_response)
            raise ValueError('Failed to download a guide tree from TimeTree.')
        timetree_newick = _bounded_response_text(newick_response, 'TimeTree').strip()
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
    session = _new_taxonomy_http_session()
    try:
        try:
            response = session.post(
                'https://api.opentreeoflife.org/v3/tnrs/match_names',
                json={
                    'names': query_names,
                    'do_approximate_matching': False,
                },
                timeout=60,
                stream=True,
            )
        except requests.RequestException as exc:
            raise ValueError('Failed to contact Open Tree of Life TNRS.') from exc
        if response.status_code != 200:
            _close_http_response(response)
            raise ValueError('Open Tree of Life TNRS lookup failed.')
        data = _response_json_object(response, 'Open Tree of Life TNRS')
        results = data.get('results', [])
        if not isinstance(results, list):
            raise ValueError('Unexpected response format from Open Tree of Life TNRS.')
        if len(results) != len(query_names):
            raise ValueError('Unexpected response format from Open Tree of Life TNRS.')
        ott_ids = list()
        unresolved_labels = list()
        for query_name, result in zip(query_names, results):
            if not isinstance(result, dict):
                raise ValueError('Unexpected response format from Open Tree of Life TNRS.')
            matches = result.get('matches', [])
            if not isinstance(matches, list):
                raise ValueError('Unexpected response format from Open Tree of Life TNRS.')
            valid_matches = list()
            for match in matches:
                if not isinstance(match, dict):
                    continue
                taxon = match.get('taxon', {})
                if not isinstance(taxon, dict):
                    continue
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
    session = _new_taxonomy_http_session()
    try:
        try:
            response = session.post(
                'https://api.opentreeoflife.org/v3/tree_of_life/induced_subtree',
                json={
                    'ott_ids': ott_ids,
                    'label_format': 'name',
                },
                timeout=60,
                stream=True,
            )
        except requests.RequestException as exc:
            raise ValueError('Failed to contact Open Tree of Life synthetic tree API.') from exc
        if response.status_code != 200:
            _close_http_response(response)
            raise ValueError('Open Tree of Life induced subtree lookup failed.')
        data = _response_json_object(response, 'Open Tree of Life induced subtree')
        newick_value = data.get('newick', '')
        if not isinstance(newick_value, str):
            raise ValueError('Unexpected response format from Open Tree of Life induced subtree.')
        newick = newick_value.strip()
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
    output_properties = set(get_tree_property_names(tree))
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
    output_properties.update(get_tree_property_names(tree))
    write_tree(tree, args, format=args.outformat, props=output_properties)
