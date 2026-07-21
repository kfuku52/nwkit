import sys
from collections import Counter

import pandas as pd

from nwkit.transfer import _validate_property_name
from nwkit.util import (
    assign_branch_ids,
    get_node_class,
    get_tree_property_names,
    is_missing_table_value,
    parse_table_missing_values,
    read_tree,
    read_tip_table,
    validate_unique_named_leaves,
    write_tree,
)


ANNOTATION_REPORT_COLUMNS = (
    'branch_id',
    'node_class',
    'descendant_taxa',
    'input_column',
    'property',
    'aggregation',
    'value',
    'status',
    'reason',
)
SUPPORTED_AGGREGATIONS = ('unique', 'mode', 'count', 'mean', 'sum', 'min', 'max', 'list')


def _format_taxa(node):
    return ','.join(sorted(str(name) for name in node.leaf_names()))


def _python_value(value, missing_values):
    if is_missing_table_value(value, missing_values):
        return None
    if hasattr(value, 'item'):
        value = value.item()
    if isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def _parse_column_specs(columns, property_maps, table_columns):
    specs = list()
    if columns:
        for column in str(columns).split(','):
            column = column.strip()
            if column:
                specs.append((column, _validate_property_name(column)))
    for raw_mapping in property_maps or []:
        if '=' not in str(raw_mapping):
            raise ValueError("--property-map must use COLUMN=PROPERTY syntax: {}".format(raw_mapping))
        column, prop = str(raw_mapping).split('=', 1)
        specs.append((column.strip(), _validate_property_name(prop.strip())))
    if not specs:
        specs = [
            (str(column), _validate_property_name(str(column)))
            for column in table_columns
            if str(column) != 'leaf_name'
        ]
    seen_properties = set()
    for column, prop in specs:
        if column not in table_columns:
            raise ValueError("Column '{}' was not found in --table.".format(column))
        if column == 'leaf_name':
            raise ValueError("The 'leaf_name' key column cannot be attached as a property.")
        if prop in seen_properties:
            raise ValueError("Output property is specified more than once: {}".format(prop))
        seen_properties.add(prop)
    return specs


def _parse_aggregation(raw, table_columns):
    parts = str(raw).split(':')
    if len(parts) not in (2, 3):
        raise ValueError("--aggregate must use COLUMN:METHOD[:PROPERTY] syntax: {}".format(raw))
    column, method = parts[0].strip(), parts[1].strip().lower()
    if column not in table_columns:
        raise ValueError("Aggregation column '{}' was not found in --table.".format(column))
    if method not in SUPPORTED_AGGREGATIONS:
        raise ValueError(
            "Unsupported aggregation '{}'. Choose from: {}.".format(
                method,
                ', '.join(SUPPORTED_AGGREGATIONS),
            )
        )
    prop = parts[2].strip() if len(parts) == 3 else '{}_{}'.format(column, method)
    return column, method, _validate_property_name(prop)


def _aggregate(values, method, column):
    if method == 'count':
        return len(values), True, 'count_non_missing'
    if not values:
        return None, False, 'no_non_missing_values'
    unique_by_text = dict()
    for value in values:
        unique_by_text.setdefault(str(value), value)
    if method == 'unique':
        if len(unique_by_text) == 1:
            return next(iter(unique_by_text.values())), True, 'single_unique_value'
        return None, False, 'multiple_values'
    if method == 'mode':
        counts = Counter(str(value) for value in values)
        winning_text = sorted(counts, key=lambda text: (-counts[text], text))[0]
        return unique_by_text[winning_text], True, 'deterministic_mode'
    if method == 'list':
        return '|'.join(sorted(unique_by_text)), True, 'sorted_unique_values'
    try:
        numeric = pd.to_numeric(pd.Series(values), errors='raise')
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "Aggregation '{}' requires numeric values in column '{}'.".format(method, column)
        ) from exc
    if method == 'mean':
        return float(numeric.mean()), True, 'numeric_mean'
    if method == 'sum':
        return float(numeric.sum()), True, 'numeric_sum'
    if method == 'min':
        return float(numeric.min()), True, 'numeric_minimum'
    if method == 'max':
        return float(numeric.max()), True, 'numeric_maximum'
    raise AssertionError('Unhandled aggregation method: {}'.format(method))


def _write_report(rows, path):
    if path in (None, ''):
        return
    if path == '-':
        raise ValueError("'--report' requires a file path, not '-'.")
    pd.DataFrame(rows, columns=ANNOTATION_REPORT_COLUMNS).to_csv(path, sep='\t', index=False)


def annotate_main(args):
    if getattr(args, 'table', None) in (None, ''):
        raise ValueError("'--table' is required for 'annotate'.")
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    validate_unique_named_leaves(tree, option_name='--infile', context=" for 'annotate'")
    table, table_only, tree_only = read_tip_table(
        args.table,
        option_name='--table',
        tree_leaf_names=tree.leaf_names(),
        unmatched=getattr(args, 'unmatched', 'warn'),
        missing_values=getattr(args, 'missing_values', None),
    )
    missing_values = parse_table_missing_values(getattr(args, 'missing_values', None))
    column_specs = _parse_column_specs(
        columns=getattr(args, 'columns', None),
        property_maps=getattr(args, 'property_map', None),
        table_columns=list(table.columns),
    )
    aggregations = [
        _parse_aggregation(raw, list(table.columns))
        for raw in (getattr(args, 'aggregate', None) or [])
    ]
    output_properties = set(get_tree_property_names(tree))
    output_properties.update(prop for _, prop in column_specs)
    output_properties.update(prop for _, _, prop in aggregations)
    rows_by_leaf = {
        str(row['leaf_name']): row
        for _, row in table.iterrows()
    }
    node_ids = assign_branch_ids(tree)
    report_rows = list()
    value_by_leaf_and_column = dict()
    for leaf in tree.leaves():
        leaf_name = str(leaf.name)
        table_row = rows_by_leaf.get(leaf_name)
        for column, prop in column_specs:
            value = None
            if table_row is None:
                status = 'unmatched'
                reason = 'tree_tip_absent_from_table'
            else:
                value = _python_value(table_row[column], missing_values)
                value_by_leaf_and_column[(leaf_name, column)] = value
                if value is None:
                    status = 'missing'
                    reason = 'missing_table_value'
                else:
                    leaf.props[prop] = value
                    status = 'annotated'
                    reason = 'matching_leaf_name'
            report_rows.append({
                'branch_id': node_ids[leaf],
                'node_class': 'leaf',
                'descendant_taxa': leaf_name,
                'input_column': column,
                'property': prop,
                'aggregation': '',
                'value': value,
                'status': status,
                'reason': reason,
            })
    # Populate values needed only for aggregation columns.
    for leaf_name, table_row in rows_by_leaf.items():
        for column, _, _ in aggregations:
            value_by_leaf_and_column[(leaf_name, column)] = _python_value(
                table_row[column],
                missing_values,
            )
    for node in tree.traverse(strategy='postorder'):
        if node.is_leaf:
            continue
        descendant_names = [str(name) for name in node.leaf_names()]
        for column, method, prop in aggregations:
            values = [
                value_by_leaf_and_column.get((leaf_name, column))
                for leaf_name in descendant_names
            ]
            values = [value for value in values if value is not None]
            value, assign, reason = _aggregate(values, method, column)
            if assign:
                node.props[prop] = value
                status = 'aggregated'
            else:
                status = 'not_aggregated'
            report_rows.append({
                'branch_id': node_ids[node],
                'node_class': get_node_class(node),
                'descendant_taxa': _format_taxa(node),
                'input_column': column,
                'property': prop,
                'aggregation': method,
                'value': value,
                'status': status,
                'reason': reason,
            })
    for leaf_name in table_only:
        for column, prop in column_specs:
            report_rows.append({
                'branch_id': '',
                'node_class': '',
                'descendant_taxa': leaf_name,
                'input_column': column,
                'property': prop,
                'aggregation': '',
                'value': _python_value(rows_by_leaf[leaf_name][column], missing_values),
                'status': 'unmatched',
                'reason': 'table_row_absent_from_tree',
            })
    _write_report(report_rows, getattr(args, 'report', None))
    annotated = sum(row['status'] in ('annotated', 'aggregated') for row in report_rows)
    sys.stderr.write('Attached or aggregated {} property values.\n'.format(annotated))
    write_tree(tree, args, format=args.outformat, props=output_properties)
