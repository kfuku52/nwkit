import re
import sys

import pandas as pd

from nwkit.clade_mapping import build_clade_mapping, projected_root_split
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


def _resolve_split_root_edge_property(match, mapping, source_prop, target_prop):
    if match.match_basis != 'split' or match.status != 'ambiguous':
        return None
    if source_prop == 'length' or target_prop == 'length':
        return None
    target_candidates = match.target_candidates
    source_candidates = match.source_candidates
    if not target_candidates or not source_candidates:
        return None
    if len(target_candidates) > 1 and not _is_bifurcating_root_pair(target_candidates):
        return None
    if len(source_candidates) > 1 and not _is_bifurcating_root_pair(source_candidates):
        return None
    if len(target_candidates) == 1 and len(source_candidates) == 1:
        return None

    values = list()
    for source_node in source_candidates:
        value, has_value = _get_property(source_node, source_prop)
        if has_value:
            values.append(value)
    if values and any(value != values[0] for value in values[1:]):
        return {
            'resolved': False,
            'reason': 'conflicting_root_edge_values',
        }
    preferred_source = next(
        (
            node for node in source_candidates
            if frozenset(str(name) for name in node.leaf_names()) == match.target_taxa
        ),
        source_candidates[0],
    )
    return {
        'resolved': True,
        'status': (
            'projected_match'
            if mapping.target_only_taxa or mapping.source_only_taxa
            else 'exact_match'
        ),
        'reason': 'matching_canonical_root_edge',
        'source': preferred_source,
        'source_taxa': frozenset(str(name) for name in preferred_source.leaf_names()),
        'source_value': values[0] if values else None,
        'source_has_value': bool(values),
    }


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
                        allow_target_reroot=True):
    if policy not in ('compatible-only', 'strict'):
        raise ValueError("Unsupported transfer policy: {}".format(policy))
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
    aligned_root_length_overrides = dict()
    if source_was_aligned_to_target_root:
        target_root_children = target.get_children()
        source_root_children = source.get_children()
        if len(target_root_children) == 2 and len(source_root_children) == 2:
            target_root_total = sum(float(child.dist or 0.0) for child in target_root_children)
            source_root_total = sum(float(child.dist or 0.0) for child in source_root_children)
            if abs(target_root_total) > 10 ** -15:
                for child in target_root_children:
                    aligned_root_length_overrides[id(child)] = (
                        source_root_total * float(child.dist or 0.0) / target_root_total
                    )
    mapping = build_clade_mapping(
        target=target,
        source=source,
        taxon_mode=taxon_mode,
        match_basis=match_basis,
    )
    match_by_target_id = {id(match.target): match for match in mapping.matches}
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
            root_edge_resolution = _resolve_split_root_edge_property(
                match=match,
                mapping=mapping,
                source_prop=source_prop,
                target_prop=target_prop,
            )
            if root_edge_resolution is not None:
                effective_reason = root_edge_resolution['reason']
                if root_edge_resolution['resolved']:
                    effective_match_status = root_edge_resolution['status']
                    effective_source = root_edge_resolution['source']
                    effective_source_taxa = root_edge_resolution['source_taxa']
                    source_value = root_edge_resolution['source_value']
                    source_has_value = root_edge_resolution['source_has_value']
            status = effective_match_status
            reason = effective_reason
            projection_only = effective_match_status == 'projected_match'
            if (
                effective_match_status in ('exact_match', 'projected_match')
                and root_edge_resolution is None
            ):
                source_value, source_has_value = _get_property(effective_source, source_prop)
            if (
                source_prop == 'length'
                and target_prop == 'length'
                and id(target_node) in aligned_root_length_overrides
                and effective_match_status in ('exact_match', 'projected_match')
            ):
                source_value = aligned_root_length_overrides[id(target_node)]
                source_has_value = True
                effective_reason = 'matching_aligned_root_edge_total'
                reason = effective_reason
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
