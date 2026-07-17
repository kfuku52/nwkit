import json
import os
import sys

import pandas as pd

from nwkit.clade_mapping import build_clade_mapping, projected_root_split
from nwkit.root import transfer_root_with_taxon_mode
from nwkit.transfer import (
    REPORT_COLUMNS,
    parse_property_specs,
    parse_root_edge_policies,
    transfer_properties,
)
from nwkit.util import get_tree_property_names, read_tree, write_tree


def _resolve_manifest_path(value, base_dir):
    if value in (None, '', '-'):
        return value
    if os.path.isabs(str(value)):
        return str(value)
    return os.path.join(base_dir, str(value))


def _load_manifest(path):
    if path in (None, ''):
        return {}
    with open(path) as handle:
        manifest = json.load(handle)
    if not isinstance(manifest, dict):
        raise ValueError('Composition manifest must contain a JSON object.')
    base_dir = os.path.dirname(os.path.realpath(path))
    for key in ('root', 'name', 'support', 'length'):
        if key in manifest:
            manifest[key] = _resolve_manifest_path(manifest[key], base_dir)
    normalized_properties = list()
    for entry in manifest.get('properties', []):
        if not isinstance(entry, dict):
            raise ValueError("Each manifest 'properties' entry must be an object.")
        if 'path' not in entry or 'source' not in entry:
            raise ValueError("Manifest property entries require 'path' and 'source'.")
        normalized = dict(entry)
        normalized['path'] = _resolve_manifest_path(entry['path'], base_dir)
        normalized['target'] = entry.get('target', entry['source'])
        normalized_properties.append(normalized)
    manifest['properties'] = normalized_properties
    if 'root_edge_policies' in manifest and not isinstance(
        manifest['root_edge_policies'], dict
    ):
        raise ValueError("Manifest 'root_edge_policies' must contain a JSON object.")
    return manifest


def _parse_property_source(raw):
    if '@' not in str(raw):
        raise ValueError(
            "--property-source must use SOURCE[@TARGET]@PATH or SOURCE=TARGET@PATH syntax."
        )
    property_spec, path = str(raw).rsplit('@', 1)
    if '=' in property_spec:
        source_prop, target_prop = property_spec.split('=', 1)
    else:
        source_prop = property_spec
        target_prop = property_spec
    specs = parse_property_specs(property_maps=['{}={}'.format(source_prop, target_prop)])
    return {
        'path': path,
        'source': specs[0][0],
        'target': specs[0][1],
    }


def _root_report_row(source_path, target, source, mapping, status, reason, taxon_mode):
    root_match = next(match for match in mapping.matches if match.target.is_root)
    split = projected_root_split(source, mapping.shared_taxa)
    split_text = ''
    if split is not None:
        split_text = '{}|{}'.format(
            ','.join(sorted(split[0])),
            ','.join(sorted(split[1])),
        )
    row = {column: '' for column in REPORT_COLUMNS}
    row.update({
        'source_file': source_path,
        'target_node_id': 1,
        'target_node_class': 'root',
        'target_taxa': ','.join(sorted(str(name) for name in target.leaf_names())),
        'shared_descendant_taxa': ','.join(sorted(mapping.shared_taxa)),
        'shared_split': split_text,
        'source_taxa': ','.join(sorted(str(name) for name in source.leaf_names())),
        'match_status': root_match.status,
        'match_basis': 'split',
        'projection_only': root_match.status == 'projected_match',
        'source_property': 'root',
        'target_property': 'root',
        'status': status,
        'reason': reason,
        'projected_value_allowed': '',
        'taxon_mode': taxon_mode,
        'shared_taxon_count': len(mapping.shared_taxa),
        'target_only_taxon_count': len(mapping.target_only_taxa),
        'source_only_taxon_count': len(mapping.source_only_taxa),
    })
    return row


def _write_report(rows, path):
    if path in (None, ''):
        return
    if path == '-':
        raise ValueError("'--report -' cannot be combined with Newick output on stdout.")
    pd.DataFrame(rows, columns=REPORT_COLUMNS).to_csv(path, sep='\t', index=False)


def _configured_sources(args):
    manifest = _load_manifest(getattr(args, 'manifest', None))
    manifest_root_edge_policies = parse_root_edge_policies(
        manifest.get('root_edge_policies', {})
    )
    cli_root_edge_policies = parse_root_edge_policies(
        getattr(args, 'root_edge_policy', None)
    )
    root_edge_policies = dict(manifest_root_edge_policies)
    root_edge_policies.update(cli_root_edge_policies)
    sources = {
        'root': getattr(args, 'root_source', None) or manifest.get('root'),
        'name': getattr(args, 'name_source', None) or manifest.get('name'),
        'support': getattr(args, 'support_source', None) or manifest.get('support'),
        'length': getattr(args, 'length_source', None) or manifest.get('length'),
        'properties': list(manifest.get('properties', [])),
        'root_edge_policies': root_edge_policies,
        'manifest_root_edge_policies': manifest_root_edge_policies,
        'cli_root_edge_policies': cli_root_edge_policies,
    }
    sources['properties'].extend(
        _parse_property_source(raw)
        for raw in (getattr(args, 'property_source', None) or [])
    )
    return sources


def compose_main(args):
    sources = _configured_sources(args)
    if not any(sources[key] for key in ('root', 'name', 'support', 'length')) and not sources['properties']:
        raise ValueError('At least one composition source must be specified.')
    taxon_mode = getattr(args, 'taxon_mode', 'exact')
    policy = getattr(args, 'policy', 'compatible-only')
    match_basis = getattr(args, 'match_basis', 'clade')
    allow_projected_values = bool(getattr(args, 'allow_projected_values', False))
    source_format = getattr(args, 'source_format', 'auto')
    target = read_tree(args.infile, args.format, args.quoted_node_names)
    report_rows = list()
    output_properties = set(get_tree_property_names(target))
    failures = list()

    root_source_path = sources['root']
    if root_source_path:
        root_source = read_tree(root_source_path, source_format, args.quoted_node_names)
        mapping = build_clade_mapping(target=target, source=root_source, taxon_mode=taxon_mode)
        root_match = next(match for match in mapping.matches if match.target.is_root)
        if policy == 'strict' and root_match.status == 'projected_match':
            reason = 'strict_policy_requires_exact_root_match'
            report_rows.append(_root_report_row(
                source_path=root_source_path,
                target=target,
                source=root_source,
                mapping=mapping,
                status='projected_match_rejected',
                reason=reason,
                taxon_mode=taxon_mode,
            ))
            failures.append(reason)
        else:
            try:
                target = transfer_root_with_taxon_mode(
                    tree_to=target,
                    tree_from=root_source,
                    taxon_mode=taxon_mode,
                    verbose=True,
                )
                report_rows.append(_root_report_row(
                    source_path=root_source_path,
                    target=target,
                    source=root_source,
                    mapping=mapping,
                    status='transferred',
                    reason='matching_root_split',
                    taxon_mode=taxon_mode,
                ))
            except ValueError as exc:
                report_rows.append(_root_report_row(
                    source_path=root_source_path,
                    target=target,
                    source=root_source,
                    mapping=mapping,
                    status='unmatched',
                    reason=str(exc),
                    taxon_mode=taxon_mode,
                ))
                failures.append(str(exc))

    built_in_sources = (
        ('name', sources['name'], ('name', 'name'), 'all', False),
        ('support', sources['support'], ('support', 'support'), 'intnode', True),
        ('length', sources['length'], ('length', 'length'), 'all', True),
    )
    for _, source_path, property_spec, target_class, exclude_root in built_in_sources:
        if not source_path:
            continue
        source_tree = read_tree(source_path, source_format, args.quoted_node_names)
        result = transfer_properties(
            target=target,
            source=source_tree,
            property_specs=[property_spec],
            target_class=target_class,
            taxon_mode=taxon_mode,
            policy=policy,
            align_roots=True,
            source_label=source_path,
            exclude_root=exclude_root,
            match_basis=match_basis,
            allow_projected_values=allow_projected_values,
            allow_target_reroot=False,
            root_edge_policies=sources['root_edge_policies'],
        )
        target = result['tree']
        report_rows.extend(result['rows'])
        output_properties.update(result['output_properties'])
        if result['error'] is not None:
            failures.append(str(result['error']))

    for property_source in sources['properties']:
        source_path = property_source['path']
        source_tree = read_tree(
            source_path,
            property_source.get('format', source_format),
            args.quoted_node_names,
        )
        property_spec = parse_property_specs(property_maps=[
            '{}={}'.format(property_source['source'], property_source['target'])
        ])
        property_root_edge_policies = dict(sources['manifest_root_edge_policies'])
        if 'root_edge_policy' in property_source:
            property_root_edge_policies.update(parse_root_edge_policies({
                property_source['target']: property_source['root_edge_policy'],
            }))
        property_root_edge_policies.update(sources['cli_root_edge_policies'])
        result = transfer_properties(
            target=target,
            source=source_tree,
            property_specs=property_spec,
            target_class=property_source.get('target_class', 'all'),
            taxon_mode=taxon_mode,
            policy=policy,
            align_roots=True,
            source_label=source_path,
            match_basis=match_basis,
            allow_projected_values=allow_projected_values,
            allow_target_reroot=False,
            root_edge_policies=property_root_edge_policies,
        )
        target = result['tree']
        report_rows.extend(result['rows'])
        output_properties.update(result['output_properties'])
        if result['error'] is not None:
            failures.append(str(result['error']))

    _write_report(report_rows, getattr(args, 'report', None))
    if failures and policy == 'strict':
        raise ValueError('Strict composition failed: {}'.format(' | '.join(failures)))
    transferred = sum(row.get('status') == 'transferred' for row in report_rows)
    skipped = sum(row.get('status') not in ('transferred', 'filled') for row in report_rows)
    sys.stderr.write(
        'Composition report: transferred={}, skipped={}\n'.format(transferred, skipped)
    )
    outformat = args.outformat
    if outformat == 'auto' and sources['name']:
        outformat = 1
    write_tree(target, args, format=outformat, props=output_properties)
