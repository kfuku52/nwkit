import sys

import pandas as pd
from ete4 import Tree

from nwkit.util import MISSING_SUPPORT_VALUE, TREE_FORMAT_PROP, is_missing_table_value, write_tree


def _read_table(path):
    if path == '-':
        return pd.read_csv(sys.stdin, sep='\t', dtype=object, keep_default_na=False)
    return pd.read_csv(path, sep='\t', dtype=object, keep_default_na=False)


def _is_empty(value):
    if pd.isna(value):
        return True
    return str(value).strip() == ''


def _parse_int(value, column_name):
    try:
        return int(str(value).strip())
    except ValueError as exc:
        raise ValueError("Column '{}' must contain integers.".format(column_name)) from exc


def _parse_optional_float(value, column_name):
    if is_missing_table_value(value):
        return None
    try:
        return float(str(value).strip())
    except ValueError as exc:
        raise ValueError("Column '{}' must contain numeric values when present.".format(column_name)) from exc


def _normalize_name(value):
    if _is_empty(value):
        return ''
    return str(value)


def _validate_parent_map(parent_by_id, root_id):
    for branch_id in parent_by_id:
        seen = set()
        current = branch_id
        while current != root_id:
            if current in seen:
                raise ValueError('The parent-child table contains a cycle involving branch_id {}.'.format(branch_id))
            seen.add(current)
            parent_id = parent_by_id.get(current)
            if parent_id is None:
                raise ValueError('Parent branch_id not found for branch_id {}.'.format(branch_id))
            if parent_id == -1:
                raise ValueError('branch_id {} is disconnected from the root.'.format(branch_id))
            if parent_id not in parent_by_id:
                raise ValueError('Parent branch_id {} was not found in the table.'.format(parent_id))
            current = parent_id


def _resolve_output_format(args, child_ids, records_by_id):
    internal_ids = {branch_id for branch_id, children in child_ids.items() if len(children) > 0}
    has_internal_names = any(records_by_id[branch_id]['name'] != '' for branch_id in internal_ids)
    has_internal_support = any(records_by_id[branch_id]['support'] is not None for branch_id in internal_ids)
    if args.outformat != 'auto':
        return args.outformat
    if has_internal_support and has_internal_names:
        sys.stderr.write(
            'Warning: Auto output format resolved to 0 because the table contains support values; '
            'internal node names may not be retained in standard Newick.\n'
        )
        return 0
    if has_internal_support:
        return 0
    return 1


def table2nwk_main(args):
    table = _read_table(args.infile)
    if table.empty:
        raise ValueError('Input table is empty.')
    required_columns = {'branch_id', 'parent'}
    missing_columns = sorted(required_columns - set(table.columns))
    if missing_columns:
        raise ValueError("Missing required column(s): {}".format(', '.join(missing_columns)))
    records = list()
    seen_branch_ids = set()
    root_ids = list()
    for row_idx, (_, row) in enumerate(table.iterrows()):
        if _is_empty(row['branch_id']):
            raise ValueError("Column 'branch_id' must not contain missing values.")
        branch_id = _parse_int(row['branch_id'], 'branch_id')
        if branch_id in seen_branch_ids:
            raise ValueError("Duplicated 'branch_id' values are not supported.")
        seen_branch_ids.add(branch_id)
        if _is_empty(row['parent']):
            parent_id = -1
        else:
            parent_id = _parse_int(row['parent'], 'parent')
        if parent_id == -1:
            root_ids.append(branch_id)
        record = {
            'branch_id': branch_id,
            'parent': parent_id,
            'name': _normalize_name(row['name']) if 'name' in row else '',
            'dist': _parse_optional_float(row['dist'], 'dist') if 'dist' in row else None,
            'support': _parse_optional_float(row['support'], 'support') if 'support' in row else None,
            'row_order': row_idx,
        }
        records.append(record)
    if len(root_ids) != 1:
        raise ValueError("The table must contain exactly one root row with parent = -1.")
    root_id = root_ids[0]
    records_by_id = {record['branch_id']: record for record in records}
    parent_by_id = {record['branch_id']: record['parent'] for record in records}
    for record in records:
        parent_id = record['parent']
        if (parent_id != -1) and (parent_id not in records_by_id):
            raise ValueError('Parent branch_id {} was not found in the table.'.format(parent_id))
    _validate_parent_map(parent_by_id, root_id)
    child_ids = dict()
    for record in records:
        child_ids.setdefault(record['branch_id'], list())
    for record in records:
        parent_id = record['parent']
        if parent_id == -1:
            continue
        child_ids[parent_id].append(record['branch_id'])
    for branch_id, children in child_ids.items():
        children.sort(key=lambda child_id: records_by_id[child_id]['row_order'])
    nodes = dict()
    for record in sorted(records, key=lambda rec: rec['row_order']):
        node = Tree()
        node.name = record['name']
        node.dist = record['dist']
        if record['support'] is None:
            node.support = None if record['branch_id'] == root_id else MISSING_SUPPORT_VALUE
        else:
            node.support = record['support']
        nodes[record['branch_id']] = node
    for record in sorted(records, key=lambda rec: rec['row_order']):
        parent_id = record['parent']
        if parent_id == -1:
            continue
        nodes[parent_id].add_child(child=nodes[record['branch_id']])
    output_tree = nodes[root_id]
    output_format = _resolve_output_format(args, child_ids, records_by_id)
    output_tree.props[TREE_FORMAT_PROP] = int(output_format)
    for node in output_tree.traverse():
        node.props[TREE_FORMAT_PROP] = int(output_format)
    write_tree(output_tree, args, format=output_format)
