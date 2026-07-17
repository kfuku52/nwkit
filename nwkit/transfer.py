import re
import sys

import pandas as pd

from nwkit.clade_mapping import build_clade_mapping
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
    'source_taxa',
    'source_property',
    'target_property',
    'source_value',
    'previous_target_value',
    'output_value',
    'status',
    'reason',
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


def _try_align_roots(target, source, taxon_mode):
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
        elif source_bifurcating:
            target = transfer_root_with_taxon_mode(
                tree_to=target,
                tree_from=source,
                taxon_mode=taxon_mode,
                verbose=False,
                redistribute_root_length=False,
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
                        align_roots=True, source_label='', exclude_root=False):
    if policy not in ('compatible-only', 'strict'):
        raise ValueError("Unsupported transfer policy: {}".format(policy))
    validate_unique_named_leaves(target, option_name='--infile', context=" for 'transfer'")
    validate_unique_named_leaves(source, option_name='--infile2', context=" for 'transfer'")
    if taxon_mode == 'exact' and not is_all_leaf_names_identical(target, source, verbose=True):
        raise ValueError('Leaf labels must match exactly when --taxon-mode exact.')
    if align_roots and target_class != 'leaf':
        target, source = _try_align_roots(target=target, source=source, taxon_mode=taxon_mode)
    mapping = build_clade_mapping(target=target, source=source, taxon_mode=taxon_mode)
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
            status = match.status
            reason = match.reason
            if match.status == 'matched':
                source_value, source_has_value = _get_property(match.source, source_prop)
                if source_has_value:
                    _set_property(target_node, target_prop, source_value)
                    output_value = source_value
                    status = 'transferred'
                    reason = 'matching_descendant_taxa'
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
                'source_taxa': _format_taxa(match.source_taxa),
                'source_property': source_prop,
                'target_property': target_prop,
                'source_value': source_value,
                'previous_target_value': previous_value,
                'output_value': output_value,
                'status': status,
                'reason': reason,
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
