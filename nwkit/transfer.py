import json
import math
import re
import sys

import pandas as pd

from nwkit.clade_mapping import (
    build_clade_mapping,
    node_projected_split,
    projected_root_split,
)
from nwkit.root import transfer_root_with_taxon_mode
from nwkit.util import (
    TREE_FORMAT_PROP,
    get_tree_property_names,
    get_target_nodes,
    is_all_leaf_names_identical,
    read_tree,
    support_is_missing,
    validate_unique_named_leaves,
    write_tree,
)


SPECIAL_PROPERTIES = frozenset(('name', 'support', 'length'))
PROPERTY_NAME_PATTERN = re.compile(r'^[A-Za-z_][A-Za-z0-9_.-]*$')
ROOT_EDGE_POLICIES = frozenset((
    'auto',
    'skip',
    'equal-only',
    'matching-side',
    'mean',
    'min',
    'max',
    'edge-total',
))
REPORT_COLUMNS = (
    'source_file',
    'target_node_id',
    'target_node_class',
    'target_taxa',
    'shared_descendant_taxa',
    'shared_split',
    'source_taxa',
    'match_status',
    'match_basis',
    'projection_only',
    'source_property',
    'target_property',
    'source_value',
    'previous_target_value',
    'output_value',
    'status',
    'reason',
    'ambiguity_policy',
    'ambiguity_resolution',
    'target_candidate_count',
    'source_candidate_count',
    'source_candidate_values',
    'projected_value_allowed',
    'taxon_mode',
    'shared_taxon_count',
    'target_only_taxon_count',
    'source_only_taxon_count',
)


def _node_class(node):
    if node.is_root:
        return 'root'
    if node.is_leaf:
        return 'leaf'
    return 'intnode'


def _format_taxa(taxa):
    if taxa is None:
        return ''
    return ','.join(sorted(str(taxon) for taxon in taxa))


def _format_split(split):
    if split is None:
        return ''
    return '{}|{}'.format(_format_taxa(split[0]), _format_taxa(split[1]))


def _validate_property_name(prop):
    if prop == TREE_FORMAT_PROP:
        raise ValueError("Property '{}' is reserved for internal NWKIT use.".format(prop))
    if prop in SPECIAL_PROPERTIES:
        return prop
    if not PROPERTY_NAME_PATTERN.fullmatch(str(prop)):
        raise ValueError(
            "Invalid property name '{}'. Use letters, numbers, underscores, dots, or hyphens.".format(
                prop
            )
        )
    return str(prop)


def parse_property_specs(properties=None, property_maps=None, include_name=False,
                         include_support=False, include_length=False):
    specs = list()
    if include_name:
        specs.append(('name', 'name'))
    if include_support:
        specs.append(('support', 'support'))
    if include_length:
        specs.append(('length', 'length'))
    for raw_property in properties or []:
        for prop in str(raw_property).split(','):
            prop = prop.strip()
            if prop:
                prop = _validate_property_name(prop)
                specs.append((prop, prop))
    for raw_mapping in property_maps or []:
        if '=' not in str(raw_mapping):
            raise ValueError("--property-map must use SOURCE=TARGET syntax: {}".format(raw_mapping))
        source_prop, target_prop = str(raw_mapping).split('=', 1)
        source_prop = _validate_property_name(source_prop.strip())
        target_prop = _validate_property_name(target_prop.strip())
        specs.append((source_prop, target_prop))
    deduplicated = list()
    seen_targets = set()
    for source_prop, target_prop in specs:
        if target_prop in seen_targets:
            raise ValueError("Target property is specified more than once: {}".format(target_prop))
        seen_targets.add(target_prop)
        deduplicated.append((source_prop, target_prop))
    return deduplicated


def parse_root_edge_policies(specs=None):
    if specs is None:
        return {}
    if isinstance(specs, dict):
        raw_items = list(specs.items())
    else:
        raw_items = list()
        for raw in specs:
            if '=' not in str(raw):
                raise ValueError(
                    "--root-edge-policy must use TARGET_PROPERTY=POLICY syntax: {}".format(raw)
                )
            target_prop, policy = str(raw).split('=', 1)
            raw_items.append((target_prop, policy))
    policies = dict()
    for raw_target_prop, raw_policy in raw_items:
        target_prop = str(raw_target_prop).strip()
        if target_prop != '*':
            target_prop = _validate_property_name(target_prop)
        policy = str(raw_policy).strip().lower()
        if policy not in ROOT_EDGE_POLICIES:
            raise ValueError(
                "Unknown root-edge policy '{}'. Choose from: {}.".format(
                    policy,
                    ', '.join(sorted(ROOT_EDGE_POLICIES)),
                )
            )
        if target_prop in policies:
            raise ValueError(
                'Root-edge policy is specified more than once for target property: {}'.format(
                    target_prop
                )
            )
        if target_prop == 'name' and policy in ('mean', 'min', 'max', 'edge-total'):
            raise ValueError("Root-edge policy '{}' is not valid for node names.".format(policy))
        if target_prop == 'length' and policy in ('equal-only', 'mean', 'min', 'max'):
            raise ValueError("Root-edge policy '{}' is not valid for branch lengths.".format(policy))
        if target_prop == 'support' and policy == 'edge-total':
            raise ValueError("Root-edge policy 'edge-total' is only valid for branch lengths.")
        policies[target_prop] = policy
    return policies


def _root_edge_policy_for(target_prop, policies):
    if target_prop in policies:
        return policies[target_prop]
    if '*' in policies:
        return policies['*']
    if target_prop == 'length':
        return 'edge-total'
    return 'auto'


def _validate_root_edge_policy_for_property(target_prop, policy):
    if target_prop == 'name' and policy in ('mean', 'min', 'max', 'edge-total'):
        raise ValueError("Root-edge policy '{}' is not valid for node names.".format(policy))
    if target_prop == 'length' and policy in ('equal-only', 'mean', 'min', 'max'):
        raise ValueError("Root-edge policy '{}' is not valid for branch lengths.".format(policy))
    if target_prop != 'length' and policy == 'edge-total':
        raise ValueError("Root-edge policy 'edge-total' is only valid for branch lengths.")


def _get_property(node, prop):
    if prop == 'name':
        value = node.name
        return value, value not in (None, '')
    if prop == 'support':
        value = node.support
        return value, not support_is_missing(value)
    if prop == 'length':
        value = node.dist
        return value, value is not None
    if prop not in node.props:
        return None, False
    value = node.props.get(prop)
    return value, value is not None


def _coerce_fill(fill, target_prop):
    if fill is None:
        return None
    if target_prop in ('support', 'length'):
        try:
            return float(fill)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "'--fill' must be numeric when support or length is transferred."
            ) from exc
    return str(fill)


def _set_property(node, prop, value):
    if prop == 'name':
        node.name = value
    elif prop == 'support':
        node.support = float(value)
    elif prop == 'length':
        node.dist = float(value)
    else:
        node.props[prop] = value


def _is_bifurcating_root_pair(nodes):
    if len(nodes) != 2:
        return False
    parent = nodes[0].up
    if parent is None or not parent.is_root or len(parent.get_children()) != 2:
        return False
    return {id(node) for node in nodes} == {
        id(node) for node in parent.get_children()
    }


def _physical_split_candidates(tree, shared_taxa):
    candidates = dict()
    taxon_sets = {
        node: frozenset(str(name) for name in taxa)
        for node, taxa in tree.get_cached_content(prop='name').items()
    }
    for node in tree.traverse():
        if node.is_root:
            continue
        split = node_projected_split(
            node,
            shared_taxa=shared_taxa,
            taxon_sets=taxon_sets,
        )
        if split is not None:
            candidates.setdefault(split, list()).append(node)
    return {
        split: tuple(nodes)
        for split, nodes in candidates.items()
    }


def _candidate_projected_taxa(node, shared_taxa):
    return frozenset(str(name) for name in node.leaf_names()).intersection(shared_taxa)


def _candidate_sort_key(node, shared_taxa):
    projected = _candidate_projected_taxa(node, shared_taxa)
    full = frozenset(str(name) for name in node.leaf_names())
    return (
        len(projected),
        tuple(sorted(projected)),
        len(full),
        tuple(sorted(full)),
    )


def _json_safe_value(value):
    if hasattr(value, 'item'):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return str(value)
    if isinstance(value, (str, int, float, bool, type(None))):
        return value
    return str(value)


def _source_candidate_records(candidates, source_prop, shared_taxa):
    records = list()
    for node in sorted(candidates, key=lambda candidate: _candidate_sort_key(candidate, shared_taxa)):
        value, present = _get_property(node, source_prop)
        records.append({
            'node': node,
            'taxa': frozenset(str(name) for name in node.leaf_names()),
            'projected_taxa': _candidate_projected_taxa(node, shared_taxa),
            'value': value,
            'present': present,
        })
    return records


def _candidate_values_json(records):
    return json.dumps(
        [
            {
                'present': record['present'],
                'projected_taxa': sorted(record['projected_taxa']),
                'taxa': sorted(record['taxa']),
                'value': _json_safe_value(record['value']) if record['present'] else None,
            }
            for record in records
        ],
        ensure_ascii=False,
        sort_keys=True,
        separators=(',', ':'),
    )


def _root_edge_resolution_metadata(policy, target_candidates, source_records):
    return {
        'ambiguity_policy': policy,
        'target_candidate_count': len(target_candidates),
        'source_candidate_count': len(source_records),
        'source_candidate_values': _candidate_values_json(source_records),
    }


def _resolved_root_edge_value(mapping, metadata, source_record, value, has_value,
                              reason, resolution):
    return {
        **metadata,
        'resolved': True,
        'status': (
            'projected_match'
            if mapping.target_only_taxa or mapping.source_only_taxa
            else 'exact_match'
        ),
        'reason': reason,
        'ambiguity_resolution': resolution,
        'source': source_record['node'],
        'source_taxa': source_record['taxa'],
        'source_value': value,
        'source_has_value': has_value,
    }


def _unresolved_root_edge_value(metadata, reason, resolution):
    return {
        **metadata,
        'resolved': False,
        'reason': reason,
        'ambiguity_resolution': resolution,
    }


def _resolve_equal_root_edge_value(mapping, metadata, source_records):
    if len(source_records) == 1:
        record = source_records[0]
        return _resolved_root_edge_value(
            mapping=mapping,
            metadata=metadata,
            source_record=record,
            value=record['value'],
            has_value=record['present'],
            reason='matching_unique_canonical_root_edge',
            resolution='unique_edge_value',
        )
    if not all(record['present'] for record in source_records):
        return None
    first_value = source_records[0]['value']
    if not all(record['value'] == first_value for record in source_records[1:]):
        return None
    return _resolved_root_edge_value(
        mapping=mapping,
        metadata=metadata,
        source_record=source_records[0],
        value=first_value,
        has_value=True,
        reason='matching_equal_root_edge_values',
        resolution='equal_edge_values',
    )


def _resolve_matching_root_edge_side(match, mapping, metadata, source_records):
    matching_records = [
        record
        for record in source_records
        if record['projected_taxa'] == match.projected_taxa
    ]
    if len(matching_records) != 1:
        return _unresolved_root_edge_value(
            metadata,
            reason='matching_root_edge_side_not_unique',
            resolution='matching_side_unresolved',
        )
    record = matching_records[0]
    return _resolved_root_edge_value(
        mapping=mapping,
        metadata=metadata,
        source_record=record,
        value=record['value'],
        has_value=record['present'],
        reason='matching_root_edge_side',
        resolution='matching_side',
    )


def _resolve_reduced_root_edge_value(mapping, metadata, source_records, policy):
    if not source_records or not all(record['present'] for record in source_records):
        return _unresolved_root_edge_value(
            metadata,
            reason='incomplete_root_edge_values',
            resolution='numeric_reducer_unresolved',
        )
    try:
        numeric_values = [float(record['value']) for record in source_records]
    except (TypeError, ValueError):
        return _unresolved_root_edge_value(
            metadata,
            reason='non_numeric_root_edge_values',
            resolution='numeric_reducer_unresolved',
        )
    if not all(math.isfinite(value) for value in numeric_values):
        return _unresolved_root_edge_value(
            metadata,
            reason='non_finite_root_edge_values',
            resolution='numeric_reducer_unresolved',
        )
    if policy == 'mean':
        value = sum(numeric_values) / len(numeric_values)
    elif policy == 'min':
        value = min(numeric_values)
    else:
        value = max(numeric_values)
    return _resolved_root_edge_value(
        mapping=mapping,
        metadata=metadata,
        source_record=source_records[0],
        value=value,
        has_value=True,
        reason='reduced_root_edge_values',
        resolution=policy,
    )


def _resolve_root_edge_total(mapping, metadata, source_records, target_candidates,
                             target_node, target_root_ratios):
    if not source_records or not all(record['present'] for record in source_records):
        return _unresolved_root_edge_value(
            metadata,
            reason='incomplete_root_edge_lengths',
            resolution='edge_total_unresolved',
        )
    try:
        numeric_values = [float(record['value']) for record in source_records]
    except (TypeError, ValueError):
        return _unresolved_root_edge_value(
            metadata,
            reason='non_numeric_root_edge_lengths',
            resolution='edge_total_unresolved',
        )
    if not all(math.isfinite(value) for value in numeric_values):
        return _unresolved_root_edge_value(
            metadata,
            reason='non_finite_root_edge_lengths',
            resolution='edge_total_unresolved',
        )
    source_total = sum(numeric_values)
    if _is_bifurcating_root_pair(target_candidates):
        ratio_record = target_root_ratios.get(id(target_node))
        if ratio_record is not None:
            ratio, ratio_was_defined = ratio_record
        else:
            target_total = sum(float(candidate.dist or 0.0) for candidate in target_candidates)
            ratio_was_defined = abs(target_total) > 10 ** -15
            ratio = (
                float(target_node.dist or 0.0) / target_total
                if ratio_was_defined
                else 1.0 / len(target_candidates)
            )
        if ratio_was_defined:
            resolution = 'edge_total_target_ratio'
        else:
            resolution = 'edge_total_equal_zero_ratio'
        output_value = source_total * ratio
    else:
        output_value = source_total
        resolution = 'edge_total_single_branch'
    return _resolved_root_edge_value(
        mapping=mapping,
        metadata=metadata,
        source_record=source_records[0],
        value=output_value,
        has_value=True,
        reason='matching_root_edge_total',
        resolution=resolution,
    )


def _resolve_split_root_edge_property(match, mapping, source_prop, target_prop, policy,
                                      target_root_ratios, target_candidates,
                                      source_candidates):
    if match.match_basis != 'split':
        return None
    if not target_candidates or not source_candidates:
        return None
    if len(target_candidates) > 1 and not _is_bifurcating_root_pair(target_candidates):
        return None
    if len(source_candidates) > 1 and not _is_bifurcating_root_pair(source_candidates):
        return None
    if len(target_candidates) == 1 and len(source_candidates) == 1:
        return None
    source_records = _source_candidate_records(
        source_candidates,
        source_prop=source_prop,
        shared_taxa=mapping.shared_taxa,
    )
    metadata = _root_edge_resolution_metadata(
        policy=policy,
        target_candidates=target_candidates,
        source_records=source_records,
    )
    if policy == 'skip':
        return _unresolved_root_edge_value(
            metadata,
            reason='root_edge_policy_skip',
            resolution='skipped',
        )
    if policy in ('auto', 'equal-only') and target_prop != 'length':
        equal_resolution = _resolve_equal_root_edge_value(
            mapping=mapping,
            metadata=metadata,
            source_records=source_records,
        )
        if equal_resolution is not None:
            return equal_resolution
        if policy == 'equal-only':
            reason = (
                'incomplete_root_edge_values'
                if not all(record['present'] for record in source_records)
                else 'conflicting_root_edge_values'
            )
            return _unresolved_root_edge_value(
                metadata,
                reason=reason,
                resolution='equal_values_required',
            )
    if policy == 'matching-side':
        return _resolve_matching_root_edge_side(
            match=match,
            mapping=mapping,
            metadata=metadata,
            source_records=source_records,
        )
    if policy in ('mean', 'min', 'max'):
        return _resolve_reduced_root_edge_value(
            mapping=mapping,
            metadata=metadata,
            source_records=source_records,
            policy=policy,
        )
    if policy in ('auto', 'edge-total') and target_prop == 'length':
        return _resolve_root_edge_total(
            mapping=mapping,
            metadata=metadata,
            source_records=source_records,
            target_candidates=target_candidates,
            target_node=match.target,
            target_root_ratios=target_root_ratios,
        )
    if policy == 'auto':
        return _resolve_matching_root_edge_side(
            match=match,
            mapping=mapping,
            metadata=metadata,
            source_records=source_records,
        )
    return _unresolved_root_edge_value(
        metadata,
        reason='root_edge_policy_not_applicable',
        resolution='unresolved',
    )


def _resolve_aligned_root_lengths(target, source, mapping, source_prop, policy):
    target_candidates = tuple(target.get_children())
    source_candidates = tuple(source.get_children())
    if len(target_candidates) != 2 or len(source_candidates) != 2:
        return {}
    source_records = _source_candidate_records(
        source_candidates,
        source_prop=source_prop,
        shared_taxa=mapping.shared_taxa,
    )
    metadata = _root_edge_resolution_metadata(
        policy=policy,
        target_candidates=target_candidates,
        source_records=source_records,
    )
    resolutions = dict()
    if policy == 'skip':
        for target_node in target_candidates:
            resolutions[id(target_node)] = _unresolved_root_edge_value(
                metadata,
                reason='root_edge_policy_skip',
                resolution='skipped_aligned_root_edge',
            )
        return resolutions
    if policy == 'matching-side':
        for target_node in target_candidates:
            projected_taxa = _candidate_projected_taxa(target_node, mapping.shared_taxa)
            matching_records = [
                record
                for record in source_records
                if record['projected_taxa'] == projected_taxa
            ]
            if len(matching_records) != 1:
                resolutions[id(target_node)] = _unresolved_root_edge_value(
                    metadata,
                    reason='matching_aligned_root_edge_side_not_unique',
                    resolution='matching_side_unresolved',
                )
                continue
            record = matching_records[0]
            resolutions[id(target_node)] = _resolved_root_edge_value(
                mapping=mapping,
                metadata=metadata,
                source_record=record,
                value=record['value'],
                has_value=record['present'],
                reason='matching_aligned_root_edge_side',
                resolution='matching_side',
            )
        return resolutions
    if policy not in ('auto', 'edge-total'):
        for target_node in target_candidates:
            resolutions[id(target_node)] = _unresolved_root_edge_value(
                metadata,
                reason='root_edge_policy_not_applicable',
                resolution='unresolved_aligned_root_edge',
            )
        return resolutions
    if not all(record['present'] for record in source_records):
        for target_node in target_candidates:
            resolutions[id(target_node)] = _unresolved_root_edge_value(
                metadata,
                reason='incomplete_root_edge_lengths',
                resolution='edge_total_unresolved',
            )
        return resolutions
    try:
        source_values = [float(record['value']) for record in source_records]
    except (TypeError, ValueError):
        source_values = []
    if not source_values or not all(math.isfinite(value) for value in source_values):
        reason = (
            'non_finite_root_edge_lengths'
            if source_values
            else 'non_numeric_root_edge_lengths'
        )
        for target_node in target_candidates:
            resolutions[id(target_node)] = _unresolved_root_edge_value(
                metadata,
                reason=reason,
                resolution='edge_total_unresolved',
            )
        return resolutions
    source_total = sum(source_values)
    target_total = sum(float(candidate.dist or 0.0) for candidate in target_candidates)
    for target_node in target_candidates:
        if abs(target_total) > 10 ** -15:
            ratio = float(target_node.dist or 0.0) / target_total
            resolution = 'edge_total_target_ratio'
        else:
            ratio = 1.0 / len(target_candidates)
            resolution = 'edge_total_equal_zero_ratio'
        matching_records = [
            record
            for record in source_records
            if record['projected_taxa']
            == _candidate_projected_taxa(target_node, mapping.shared_taxa)
        ]
        source_record = matching_records[0] if len(matching_records) == 1 else source_records[0]
        resolutions[id(target_node)] = _resolved_root_edge_value(
            mapping=mapping,
            metadata=metadata,
            source_record=source_record,
            value=source_total * ratio,
            has_value=True,
            reason='matching_aligned_root_edge_total',
            resolution=resolution,
        )
    return resolutions


def _try_align_roots(target, source, taxon_mode, allow_target_reroot=True):
    if (len(list(target.leaves())) <= 1) or (len(list(source.leaves())) <= 1):
        return target, source
    target_bifurcating = len(target.get_children()) == 2
    source_bifurcating = len(source.get_children()) == 2
    try:
        if target_bifurcating:
            source = transfer_root_with_taxon_mode(
                tree_to=source,
                tree_from=target,
                taxon_mode=taxon_mode,
                verbose=False,
                redistribute_root_length=False,
            )
        elif source_bifurcating and allow_target_reroot:
            target = transfer_root_with_taxon_mode(
                tree_to=target,
                tree_from=source,
                taxon_mode=taxon_mode,
                verbose=False,
                redistribute_root_length=False,
            )
        elif source_bifurcating:
            sys.stderr.write(
                'Skipping root alignment because changing the target root is disabled.\n'
            )
        else:
            sys.stderr.write(
                'Skipping root alignment because neither input tree has a bifurcating root.\n'
            )
    except (SystemExit, ValueError) as exc:
        sys.stderr.write('Skipping root alignment: {}\n'.format(exc))
    return target, source


def _write_report(rows, report_path):
    if report_path in (None, ''):
        return
    if report_path == '-':
        raise ValueError("'--report -' cannot be combined with Newick output on stdout.")
    pd.DataFrame(rows, columns=REPORT_COLUMNS).to_csv(report_path, sep='\t', index=False)


def transfer_properties(target, source, property_specs, target_class='all',
                        taxon_mode='exact', fill=None, policy='compatible-only',
                        align_roots=True, source_label='', exclude_root=False,
                        match_basis='clade', allow_projected_values=False,
                        allow_target_reroot=True, root_edge_policies=None):
    if policy not in ('compatible-only', 'strict'):
        raise ValueError("Unsupported transfer policy: {}".format(policy))
    root_edge_policies = parse_root_edge_policies(root_edge_policies)
    validate_unique_named_leaves(target, option_name='--infile', context=" for 'transfer'")
    validate_unique_named_leaves(source, option_name='--infile2', context=" for 'transfer'")
    if taxon_mode == 'exact' and not is_all_leaf_names_identical(target, source, verbose=True):
        raise ValueError('Leaf labels must match exactly when --taxon-mode exact.')
    shared_taxa_before_alignment = frozenset(
        str(name) for name in target.leaf_names()
    ).intersection(str(name) for name in source.leaf_names())
    target_root_split_before_alignment = projected_root_split(
        target,
        shared_taxa_before_alignment,
    )
    source_root_split_before_alignment = projected_root_split(
        source,
        shared_taxa_before_alignment,
    )
    target_was_bifurcating = len(target.get_children()) == 2
    if align_roots and target_class != 'leaf' and match_basis == 'clade':
        target, source = _try_align_roots(
            target=target,
            source=source,
            taxon_mode=taxon_mode,
            allow_target_reroot=allow_target_reroot,
        )
    source_was_aligned_to_target_root = (
        align_roots
        and target_class != 'leaf'
        and match_basis == 'clade'
        and target_was_bifurcating
        and source_root_split_before_alignment != target_root_split_before_alignment
        and projected_root_split(source, shared_taxa_before_alignment)
        == projected_root_split(target, shared_taxa_before_alignment)
    )
    mapping = build_clade_mapping(
        target=target,
        source=source,
        taxon_mode=taxon_mode,
        match_basis=match_basis,
    )
    match_by_target_id = {id(match.target): match for match in mapping.matches}
    target_physical_split_candidates = _physical_split_candidates(
        target,
        mapping.shared_taxa,
    )
    source_physical_split_candidates = _physical_split_candidates(
        source,
        mapping.shared_taxa,
    )
    target_root_ratios = dict()
    target_root_children = target.get_children()
    if len(target_root_children) == 2:
        target_root_total = sum(float(child.dist or 0.0) for child in target_root_children)
        ratio_was_defined = abs(target_root_total) > 10 ** -15
        for child in target_root_children:
            target_root_ratios[id(child)] = (
                (
                    float(child.dist or 0.0) / target_root_total
                    if ratio_was_defined
                    else 1.0 / len(target_root_children)
                ),
                ratio_was_defined,
            )
    aligned_root_length_resolutions = dict()
    source_length_prop = next(
        (
            source_prop
            for source_prop, target_prop in property_specs
            if target_prop == 'length'
        ),
        None,
    )
    if source_was_aligned_to_target_root and source_length_prop is not None:
        length_root_edge_policy = _root_edge_policy_for('length', root_edge_policies)
        _validate_root_edge_policy_for_property('length', length_root_edge_policy)
        aligned_root_length_resolutions = _resolve_aligned_root_lengths(
            target=target,
            source=source,
            mapping=mapping,
            source_prop=source_length_prop,
            policy=length_root_edge_policy,
        )
    target_nodes = get_target_nodes(tree=target, target=target_class)
    if exclude_root:
        target_nodes = [node for node in target_nodes if not node.is_root]
    traversal_ids = {id(node): index for index, node in enumerate(target.traverse(), start=1)}
    rows = list()
    failed = False
    changed_node_ids = set()
    output_properties = set(get_tree_property_names(target))
    for source_prop, target_prop in property_specs:
        _validate_property_name(source_prop)
        _validate_property_name(target_prop)
        _validate_root_edge_policy_for_property(
            target_prop,
            _root_edge_policy_for(target_prop, root_edge_policies),
        )
        if target_prop not in SPECIAL_PROPERTIES or target_prop == 'support':
            output_properties.add(target_prop)
    for target_node in target_nodes:
        match = match_by_target_id[id(target_node)]
        for source_prop, target_prop in property_specs:
            previous_value, _ = _get_property(target_node, target_prop)
            source_value = None
            output_value = previous_value
            effective_match_status = match.status
            effective_reason = match.reason
            effective_source = match.source
            effective_source_taxa = match.source_taxa
            source_has_value = False
            projected_value_allowed = ''
            ambiguity_policy = ''
            ambiguity_resolution = ''
            target_candidate_count = ''
            source_candidate_count = ''
            source_candidate_values = ''
            root_edge_policy = _root_edge_policy_for(target_prop, root_edge_policies)
            if target_prop == 'length' and match.projected_split is not None:
                target_root_edge_candidates = target_physical_split_candidates.get(
                    match.projected_split,
                    (),
                )
                source_root_edge_candidates = source_physical_split_candidates.get(
                    match.projected_split,
                    (),
                )
            else:
                target_root_edge_candidates = match.target_candidates
                source_root_edge_candidates = match.source_candidates
            root_edge_resolution = _resolve_split_root_edge_property(
                match=match,
                mapping=mapping,
                source_prop=source_prop,
                target_prop=target_prop,
                policy=root_edge_policy,
                target_root_ratios=target_root_ratios,
                target_candidates=target_root_edge_candidates,
                source_candidates=source_root_edge_candidates,
            )
            if (
                root_edge_resolution is None
                and target_prop == 'length'
            ):
                root_edge_resolution = aligned_root_length_resolutions.get(id(target_node))
            if root_edge_resolution is not None:
                ambiguity_policy = root_edge_resolution['ambiguity_policy']
                ambiguity_resolution = root_edge_resolution['ambiguity_resolution']
                target_candidate_count = root_edge_resolution['target_candidate_count']
                source_candidate_count = root_edge_resolution['source_candidate_count']
                source_candidate_values = root_edge_resolution['source_candidate_values']
                effective_reason = root_edge_resolution['reason']
                if root_edge_resolution['resolved']:
                    effective_match_status = root_edge_resolution['status']
                    effective_source = root_edge_resolution['source']
                    effective_source_taxa = root_edge_resolution['source_taxa']
                    source_value = root_edge_resolution['source_value']
                    source_has_value = root_edge_resolution['source_has_value']
                else:
                    effective_match_status = 'ambiguous'
            status = effective_match_status
            reason = effective_reason
            projection_only = effective_match_status == 'projected_match'
            if (
                effective_match_status in ('exact_match', 'projected_match')
                and root_edge_resolution is None
            ):
                source_value, source_has_value = _get_property(effective_source, source_prop)
            if effective_match_status in ('exact_match', 'projected_match'):
                projected_value_allowed = (
                    not projection_only
                    or (
                        source_prop not in ('support', 'length')
                        and target_prop not in ('support', 'length')
                    )
                    or bool(allow_projected_values)
                )
            if projection_only and policy == 'strict':
                status = 'projected_match_rejected'
                reason = 'strict_policy_requires_exact_match'
                failed = True
            elif projection_only and not projected_value_allowed:
                status = 'projected_value_rejected'
                reason = 'projected_support_or_length_requires_opt_in'
                failed = True
            elif effective_match_status in ('exact_match', 'projected_match'):
                if source_has_value:
                    _set_property(target_node, target_prop, source_value)
                    output_value = source_value
                    status = 'transferred'
                    reason = effective_reason
                    changed_node_ids.add(id(target_node))
                elif fill is not None:
                    output_value = _coerce_fill(fill, target_prop)
                    _set_property(target_node, target_prop, output_value)
                    status = 'filled'
                    reason = 'source_property_missing'
                    changed_node_ids.add(id(target_node))
                    failed = True
                else:
                    status = 'missing_source_value'
                    reason = 'source_property_missing'
                    failed = True
            elif fill is not None:
                output_value = _coerce_fill(fill, target_prop)
                _set_property(target_node, target_prop, output_value)
                status = 'filled'
                changed_node_ids.add(id(target_node))
                failed = True
            else:
                failed = True
            rows.append({
                'source_file': source_label,
                'target_node_id': traversal_ids[id(target_node)],
                'target_node_class': _node_class(target_node),
                'target_taxa': _format_taxa(match.target_taxa),
                'shared_descendant_taxa': _format_taxa(match.projected_taxa),
                'shared_split': _format_split(match.projected_split),
                'source_taxa': _format_taxa(effective_source_taxa),
                'match_status': effective_match_status,
                'match_basis': match.match_basis,
                'projection_only': projection_only,
                'source_property': source_prop,
                'target_property': target_prop,
                'source_value': source_value,
                'previous_target_value': previous_value,
                'output_value': output_value,
                'status': status,
                'reason': reason,
                'ambiguity_policy': ambiguity_policy,
                'ambiguity_resolution': ambiguity_resolution,
                'target_candidate_count': target_candidate_count,
                'source_candidate_count': source_candidate_count,
                'source_candidate_values': source_candidate_values,
                'projected_value_allowed': projected_value_allowed,
                'taxon_mode': taxon_mode,
                'shared_taxon_count': len(mapping.shared_taxa),
                'target_only_taxon_count': len(mapping.target_only_taxa),
                'source_only_taxon_count': len(mapping.source_only_taxa),
            })
    if policy == 'strict' and failed:
        error = ValueError('Strict transfer failed because at least one requested value was not transferable.')
    else:
        error = None
    return {
        'tree': target,
        'rows': rows,
        'output_properties': output_properties,
        'changed_node_count': len(changed_node_ids),
        'target_node_count': len(target_nodes),
        'mapping': mapping,
        'error': error,
    }


def transfer_main(args):
    if args.infile2 in ('', None):
        raise ValueError("'--infile2' is required for 'transfer'.")
    property_specs = parse_property_specs(
        properties=getattr(args, 'property', None),
        property_maps=getattr(args, 'property_map', None),
        include_name=bool(getattr(args, 'name', False)),
        include_support=bool(getattr(args, 'support', False)),
        include_length=bool(getattr(args, 'length', False)),
    )
    target = read_tree(args.infile, args.format, args.quoted_node_names)
    source = read_tree(args.infile2, args.format2, args.quoted_node_names)
    result = transfer_properties(
        target=target,
        source=source,
        property_specs=property_specs,
        target_class=args.target,
        taxon_mode=getattr(args, 'taxon_mode', 'exact'),
        fill=getattr(args, 'fill', None),
        policy=getattr(args, 'policy', 'compatible-only'),
        align_roots=True,
        source_label=str(args.infile2),
        match_basis=getattr(args, 'match_basis', 'clade'),
        allow_projected_values=bool(getattr(args, 'allow_projected_values', False)),
        root_edge_policies=parse_root_edge_policies(
            getattr(args, 'root_edge_policy', None)
        ),
    )
    _write_report(result['rows'], getattr(args, 'report', None))
    sys.stderr.write(
        'Transferred {} nodes out of {} target nodes.\n'.format(
            result['changed_node_count'],
            result['target_node_count'],
        )
    )
    if result['error'] is not None:
        raise result['error']
    outformat = args.outformat
    if outformat == 'auto' and any(target_prop == 'name' for _, target_prop in property_specs):
        outformat = 1
    write_tree(
        result['tree'],
        args,
        format=outformat,
        props=result['output_properties'],
    )
