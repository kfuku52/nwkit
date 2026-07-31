import json
import math

import pandas as pd

from nwkit.clade_mapping import (
    build_clade_mapping,
    node_projected_split,
    projected_root_split,
)
from nwkit.transfer import _get_property, parse_property_specs
from nwkit.util import (
    assign_branch_ids,
    get_node_class,
    get_subtree_leaf_name_sets,
    get_target_nodes,
    read_tree,
    support_is_missing,
    validate_unique_named_leaves,
)


DIFF_COLUMNS = (
    'record_type',
    'comparison',
    'status',
    'reason',
    'transferable',
    'target_branch_id',
    'source_branch_id',
    'shared_taxa_or_split',
    'target_taxa',
    'source_taxa',
    'target_name',
    'source_name',
    'target_support',
    'source_support',
    'support_delta',
    'target_length',
    'source_length',
    'length_delta',
    'target_properties',
    'source_properties',
)

NUMERIC_REL_TOLERANCE = 10 ** -9
NUMERIC_ABS_TOLERANCE = 10 ** -12

DIFFERENCE_STATUSES = frozenset((
    'different',
    'unresolved',
    'target_only',
    'source_only',
    'ambiguous',
    'unmatched',
))

MATCHED_STATUSES = frozenset(('exact_match', 'projected_match'))


def _format_taxa(taxa):
    if taxa is None:
        return ''
    return ','.join(sorted(str(taxon) for taxon in taxa))


def _format_split(split):
    if split is None:
        return ''
    return '{}|{}'.format(_format_taxa(split[0]), _format_taxa(split[1]))


def _taxon_sets(tree):
    return get_subtree_leaf_name_sets(tree)


def _delta(target_value, source_value):
    try:
        target_number = float(target_value)
        source_number = float(source_value)
    except (TypeError, ValueError):
        return ''
    if not math.isfinite(target_number) or not math.isfinite(source_number):
        raise ValueError('Numeric values must be finite to calculate a delta.')
    scale = max(abs(target_number), abs(source_number))
    if scale == 0.0:
        return 0.0
    delta = scale * math.fsum((
        target_number / scale,
        -source_number / scale,
    ))
    if not math.isfinite(delta):
        raise ValueError('The numeric delta is too large to represent.')
    return delta


def _support_value(node):
    if node is None or support_is_missing(node.support):
        return ''
    return node.support


def _name_value(node):
    if node is None or node.name is None:
        return ''
    return node.name


def _property_json(node, properties):
    if node is None or not properties:
        return ''
    values = dict()
    for source_prop, output_prop in properties:
        value, present = _get_property(node, source_prop)
        if present:
            if hasattr(value, 'item'):
                value = value.item()
            if not isinstance(value, (str, int, float, bool, type(None))):
                value = str(value)
            values[output_prop] = value
    return json.dumps(values, sort_keys=True, ensure_ascii=False)


def _values_equal(target_value, source_value):
    if isinstance(target_value, dict) and isinstance(source_value, dict):
        return (
            target_value.keys() == source_value.keys()
            and all(
                _values_equal(target_value[key], source_value[key])
                for key in target_value
            )
        )
    if (
        isinstance(target_value, (int, float))
        and not isinstance(target_value, bool)
        and isinstance(source_value, (int, float))
        and not isinstance(source_value, bool)
    ):
        return math.isclose(
            float(target_value),
            float(source_value),
            rel_tol=NUMERIC_REL_TOLERANCE,
            abs_tol=NUMERIC_ABS_TOLERANCE,
        )
    return target_value == source_value


def _property_values(serialized):
    if serialized in ('', None):
        return {}
    return json.loads(serialized)


_MISSING_EDGE_VALUE = object()


def _node_report_values(node, properties):
    if node is None:
        return {
            'name': '',
            'support': '',
            'length': '',
            'properties': '',
        }
    return {
        'name': _name_value(node),
        'support': _support_value(node),
        'length': node.dist if node.dist is not None else '',
        'properties': _property_json(node, properties),
    }


def _coalesce_root_edge_annotation(values):
    """Return the unique/equal non-missing value, or flag a conflict."""
    present = [value for value in values if value is not _MISSING_EDGE_VALUE]
    if not present:
        return _MISSING_EDGE_VALUE, False
    first = present[0]
    if all(_values_equal(first, value) for value in present[1:]):
        return first, False
    return _MISSING_EDGE_VALUE, True


def _root_edge_report_values(nodes, properties):
    """Aggregate the two halves introduced by rooting a physical edge."""
    conflicts = []

    names = [
        node.name if node.name not in (None, '') else _MISSING_EDGE_VALUE
        for node in nodes
    ]
    name, conflict = _coalesce_root_edge_annotation(names)
    if conflict:
        conflicts.append('name')

    supports = [
        _MISSING_EDGE_VALUE if support_is_missing(node.support) else node.support
        for node in nodes
    ]
    support, conflict = _coalesce_root_edge_annotation(supports)
    if conflict:
        conflicts.append('support')

    if any(node.dist is None for node in nodes):
        length = ''
    else:
        try:
            length = math.fsum(float(node.dist) for node in nodes)
        except OverflowError as exc:
            raise ValueError(
                'The physical root-edge length total is too large to report.'
            ) from exc
        if not math.isfinite(length):
            raise ValueError(
                'The physical root-edge length total must be finite.'
            )

    node_properties = [
        _property_values(_property_json(node, properties))
        for node in nodes
    ]
    property_values = {}
    for _, output_prop in properties:
        values = [
            values.get(output_prop, _MISSING_EDGE_VALUE)
            for values in node_properties
        ]
        value, conflict = _coalesce_root_edge_annotation(values)
        if conflict:
            conflicts.append(output_prop)
        elif value is not _MISSING_EDGE_VALUE:
            property_values[output_prop] = value

    return {
        'name': '' if name is _MISSING_EDGE_VALUE else name,
        'support': '' if support is _MISSING_EDGE_VALUE else support,
        'length': length,
        'properties': (
            json.dumps(property_values, sort_keys=True, ensure_ascii=False)
            if property_values else ''
        ),
    }, tuple(conflicts)


def _row_has_value_difference(row):
    if row['record_type'] != 'node' or row['status'] not in MATCHED_STATUSES:
        return False
    return any((
        not _values_equal(row['target_name'], row['source_name']),
        not _values_equal(row['target_support'], row['source_support']),
        not _values_equal(row['target_length'], row['source_length']),
        not _values_equal(
            _property_values(row['target_properties']),
            _property_values(row['source_properties']),
        ),
    ))


def _rows_have_differences(rows, comparison='rooted'):
    return any(
        (
            row['status'] in DIFFERENCE_STATUSES
            or _row_has_value_difference(row)
        )
        and not (
            comparison == 'unrooted'
            and row['record_type'] == 'summary'
            and row['comparison'] == 'root_split'
        )
        for row in rows
    )


def _row_for_nodes(status, reason, comparison, target_node, source_node,
                   shared_key, target_taxa, source_taxa, target_ids, source_ids,
                   properties, target_values=None, source_values=None):
    if target_values is None:
        target_values = _node_report_values(target_node, properties)
    if source_values is None:
        source_values = _node_report_values(source_node, properties)
    target_support = target_values['support']
    source_support = source_values['support']
    target_length = target_values['length']
    source_length = source_values['length']
    return {
        'record_type': 'node',
        'comparison': comparison,
        'status': status,
        'reason': reason,
        'transferable': status == 'exact_match',
        'target_branch_id': target_ids.get(id(target_node), '') if target_node is not None else '',
        'source_branch_id': source_ids.get(id(source_node), '') if source_node is not None else '',
        'shared_taxa_or_split': shared_key,
        'target_taxa': _format_taxa(target_taxa),
        'source_taxa': _format_taxa(source_taxa),
        'target_name': target_values['name'],
        'source_name': source_values['name'],
        'target_support': target_support,
        'source_support': source_support,
        'support_delta': _delta(target_support, source_support),
        'target_length': target_length,
        'source_length': source_length,
        'length_delta': _delta(target_length, source_length),
        'target_properties': target_values['properties'],
        'source_properties': source_values['properties'],
    }


def _summary_row(comparison, status, reason, shared_key='', target_taxa='', source_taxa=''):
    row = {column: '' for column in DIFF_COLUMNS}
    row.update({
        'record_type': 'summary',
        'comparison': comparison,
        'status': status,
        'reason': reason,
        'transferable': '',
        'shared_taxa_or_split': shared_key,
        'target_taxa': target_taxa,
        'source_taxa': source_taxa,
    })
    return row


def _selected(node, target_class):
    if target_class == 'all':
        return True
    return get_node_class(node) == target_class


def _rooted_rows(target, source, mapping, target_class, properties, target_ids, source_ids):
    rows = list()
    selected_target_ids = {id(node) for node in get_target_nodes(target, target_class)}
    matched_source_ids = set()
    for match in mapping.matches:
        if id(match.target) not in selected_target_ids:
            continue
        if match.source is not None:
            matched_source_ids.add(id(match.source))
        status = match.status
        rows.append(_row_for_nodes(
            status=status,
            reason=match.reason,
            comparison='rooted_clade',
            target_node=match.target,
            source_node=match.source,
            shared_key=_format_taxa(match.projected_taxa),
            target_taxa=match.target_taxa,
            source_taxa=match.source_taxa,
            target_ids=target_ids,
            source_ids=source_ids,
            properties=properties,
        ))
    source_taxa_by_node = _taxon_sets(source)
    for source_node in source.traverse():
        if not _selected(source_node, target_class):
            continue
        if id(source_node) in matched_source_ids:
            continue
        projected = source_taxa_by_node[source_node] & mapping.shared_taxa
        rows.append(_row_for_nodes(
            status='source_only',
            reason='clade_absent_from_target',
            comparison='rooted_clade',
            target_node=None,
            source_node=source_node,
            shared_key=_format_taxa(projected),
            target_taxa=None,
            source_taxa=source_taxa_by_node[source_node],
            target_ids=target_ids,
            source_ids=source_ids,
            properties=properties,
        ))
    return rows


def _unrooted_rows(target, source, shared_taxa, target_class, properties,
                   target_ids, source_ids, same_leaf_set):
    target_taxa_by_node = _taxon_sets(target)
    source_taxa_by_node = _taxon_sets(source)

    def split_selected(split):
        if target_class == 'all':
            return True
        if target_class == 'root':
            return False
        is_terminal = min(len(split[0]), len(split[1])) == 1
        if target_class == 'leaf':
            return is_terminal
        if target_class == 'intnode':
            return not is_terminal
        return False

    def eligible_nodes(tree, taxa_by_node):
        result = list()
        for node in tree.traverse():
            if node.is_root:
                continue
            split = node_projected_split(node, shared_taxa, taxon_sets=taxa_by_node)
            if split is not None and split_selected(split):
                result.append((node, split))
        return result

    def is_root_pair(candidates):
        if len(candidates) != 2:
            return False
        parent = candidates[0].up
        return (
            parent is not None
            and parent.is_root
            and candidates[1].up is parent
            and len(parent.get_children()) == 2
            and {id(node) for node in candidates}
            == {id(node) for node in parent.get_children()}
        )

    def resolved_candidate(candidates):
        if len(candidates) == 1:
            node = candidates[0]
            return node, True, _node_report_values(node, properties), ()
        if is_root_pair(candidates):
            values, conflicts = _root_edge_report_values(candidates, properties)
            return None, True, values, conflicts
        return None, False, _node_report_values(None, properties), ()

    def root_edge_conflict_reason(target_conflicts, source_conflicts):
        details = []
        if target_conflicts:
            details.append('target={}'.format(','.join(target_conflicts)))
        if source_conflicts:
            details.append('source={}'.format(','.join(source_conflicts)))
        return 'conflicting_root_edge_annotations:{}'.format(';'.join(details))

    target_nodes = eligible_nodes(target, target_taxa_by_node)
    source_nodes = eligible_nodes(source, source_taxa_by_node)
    target_groups = dict()
    source_groups = dict()
    for node, split in target_nodes:
        target_groups.setdefault(split, list()).append(node)
    for node, split in source_nodes:
        source_groups.setdefault(split, list()).append(node)
    rows = list()
    all_splits = sorted(
        set(target_groups) | set(source_groups),
        key=lambda split: _format_split(split),
    )
    for split in all_splits:
        target_candidates = target_groups.get(split, [])
        source_candidates = source_groups.get(split, [])
        (
            target_node,
            target_resolved,
            target_values,
            target_conflicts,
        ) = resolved_candidate(target_candidates)
        (
            source_node,
            source_resolved,
            source_values,
            source_conflicts,
        ) = resolved_candidate(source_candidates)
        if target_candidates and source_candidates and target_resolved and source_resolved:
            if target_conflicts or source_conflicts:
                status = 'ambiguous'
                reason = root_edge_conflict_reason(
                    target_conflicts,
                    source_conflicts,
                )
            else:
                status = 'exact_match' if same_leaf_set else 'projected_match'
                reason = 'matching_unrooted_split'
        elif not target_candidates:
            status = 'source_only'
            reason = 'split_absent_from_target'
        elif not source_candidates:
            status = 'target_only'
            reason = 'split_absent_from_source'
        else:
            status = 'ambiguous'
            reason = 'projected_split_not_unique'
        rows.append(_row_for_nodes(
            status=status,
            reason=reason,
            comparison='unrooted_split',
            target_node=target_node,
            source_node=source_node,
            shared_key=_format_split(split),
            target_taxa=(
                target_taxa_by_node[target_candidates[0]]
                if target_candidates else None
            ),
            source_taxa=(
                source_taxa_by_node[source_candidates[0]]
                if source_candidates else None
            ),
            target_ids=target_ids,
            source_ids=source_ids,
            properties=properties,
            target_values=target_values,
            source_values=source_values,
        ))
    return rows


def compare_trees(target, source, taxon_mode='exact', comparison='rooted',
                  target_class='intnode', properties=None):
    validate_unique_named_leaves(target, option_name='--infile', context=" for 'diff'")
    validate_unique_named_leaves(source, option_name='--infile2', context=" for 'diff'")
    mapping = build_clade_mapping(target=target, source=source, taxon_mode=taxon_mode)
    target_ids = {id(node): branch_id for node, branch_id in assign_branch_ids(target).items()}
    source_ids = {id(node): branch_id for node, branch_id in assign_branch_ids(source).items()}
    rows = [
        _summary_row(
            comparison='leaf_set',
            status='same' if not mapping.target_only_taxa and not mapping.source_only_taxa else 'different',
            reason='shared={},target_only={},source_only={}'.format(
                len(mapping.shared_taxa),
                len(mapping.target_only_taxa),
                len(mapping.source_only_taxa),
            ),
            shared_key=_format_taxa(mapping.shared_taxa),
            target_taxa=_format_taxa(mapping.target_only_taxa),
            source_taxa=_format_taxa(mapping.source_only_taxa),
        )
    ]
    target_root_split = projected_root_split(target, mapping.shared_taxa)
    source_root_split = projected_root_split(source, mapping.shared_taxa)
    root_status = 'same' if target_root_split == source_root_split and target_root_split is not None else 'different'
    if target_root_split is None or source_root_split is None:
        root_status = 'unresolved'
    rows.append(_summary_row(
        comparison='root_split',
        status=root_status,
        reason='projected_root_bipartition',
        shared_key='target={};source={}'.format(
            _format_split(target_root_split),
            _format_split(source_root_split),
        ),
    ))
    if comparison == 'rooted':
        rows.extend(_rooted_rows(
            target=target,
            source=source,
            mapping=mapping,
            target_class=target_class,
            properties=properties or [],
            target_ids=target_ids,
            source_ids=source_ids,
        ))
    elif comparison == 'unrooted':
        rows.extend(_unrooted_rows(
            target=target,
            source=source,
            shared_taxa=mapping.shared_taxa,
            target_class=target_class,
            properties=properties or [],
            target_ids=target_ids,
            source_ids=source_ids,
            same_leaf_set=not mapping.target_only_taxa and not mapping.source_only_taxa,
        ))
    else:
        raise ValueError("Unsupported comparison mode: {}".format(comparison))
    return rows


def diff_main(args):
    if args.infile2 in ('', None):
        raise ValueError("'--infile2' is required for 'diff'.")
    target = read_tree(args.infile, args.format, args.quoted_node_names)
    source = read_tree(args.infile2, args.format2, args.quoted_node_names)
    properties = parse_property_specs(properties=getattr(args, 'property', None))
    rows = compare_trees(
        target=target,
        source=source,
        taxon_mode=getattr(args, 'taxon_mode', 'exact'),
        comparison=getattr(args, 'comparison', 'rooted'),
        target_class=getattr(args, 'target', 'intnode'),
        properties=properties,
    )
    dataframe = pd.DataFrame(rows, columns=DIFF_COLUMNS)
    if args.outfile == '-':
        print(dataframe.to_csv(sep='\t', index=False), end='')
    else:
        dataframe.to_csv(args.outfile, sep='\t', index=False)
    if (
        getattr(args, 'fail_on_difference', False)
        and _rows_have_differences(rows, comparison=getattr(args, 'comparison', 'rooted'))
    ):
        raise ValueError('Tree differences were detected.')
