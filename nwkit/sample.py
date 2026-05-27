import math
import os
import sys

import pandas as pd

from nwkit.util import read_tree, read_tsv_preserving_leaf_name, validate_unique_named_leaves, write_tree


NUMERIC_OPERATORS = {'ge', 'gt', 'le', 'lt'}
COMPARISON_OPERATORS = NUMERIC_OPERATORS | {'eq', 'ne'}
RANK_DIRECTIONS = {'asc', 'desc'}


def read_sample_trait(args, tree):
    if args.trait is None:
        return pd.DataFrame({'leaf_name': list(tree.leaf_names())})
    trait_df = read_tsv_preserving_leaf_name(args.trait)
    if 'leaf_name' not in trait_df.columns:
        raise ValueError("Column 'leaf_name' is required in '--trait'.")
    trait_df = trait_df.copy()
    trait_df['leaf_name'] = [str(leaf_name) for leaf_name in trait_df['leaf_name'].tolist()]
    if any(leaf_name.strip() == '' for leaf_name in trait_df['leaf_name'].tolist()):
        raise ValueError("Column 'leaf_name' in '--trait' must not contain empty values.")
    leaf_names_list = list(tree.leaf_names())
    leaf_name_set = set(leaf_names_list)
    unknown_leaf_names = sorted(
        leaf_name
        for leaf_name in set(trait_df['leaf_name'].tolist())
        if leaf_name not in leaf_name_set
    )
    if unknown_leaf_names:
        raise ValueError(
            "The following 'leaf_name' values in '--trait' were not found in the input tree: {}".format(
                ', '.join(unknown_leaf_names)
            )
        )
    trait_df = trait_df[trait_df['leaf_name'].isin(leaf_name_set)]
    duplicated_leaf_names = trait_df.loc[
        trait_df['leaf_name'].duplicated(keep=False), 'leaf_name'
    ].unique().tolist()
    duplicated_leaf_names = sorted(str(name) for name in duplicated_leaf_names)
    if duplicated_leaf_names:
        raise ValueError("Duplicated 'leaf_name' entries in '--trait': {}".format(', '.join(duplicated_leaf_names)))
    observed_leaf_names = set(trait_df['leaf_name'].tolist())
    missing_leaf_names = [leaf_name for leaf_name in leaf_names_list if leaf_name not in observed_leaf_names]
    if len(missing_leaf_names) > 0:
        log_txt = ''.join(
            f"'{leaf_name}' not found in '{args.trait}'. Treating its metadata as missing.\n"
            for leaf_name in missing_leaf_names
        )
        sys.stderr.write(log_txt)
        trait_df = pd.concat(
            [trait_df, pd.DataFrame({'leaf_name': missing_leaf_names})],
            ignore_index=True,
        )
    return trait_df


def parse_filter_spec(spec):
    parts = str(spec).split(':', 2)
    if len(parts) != 3:
        raise ValueError(
            "Invalid --filter specification: '{}'. Expected COLUMN:OP:VALUE.".format(spec)
        )
    column, operator, value = parts
    column = column.strip()
    operator = operator.strip().lower()
    if column == '':
        raise ValueError("Invalid --filter specification: column name is empty.")
    if operator not in COMPARISON_OPERATORS:
        raise ValueError(
            "Invalid --filter operator '{}'. Supported operators: {}.".format(
                operator, ', '.join(sorted(COMPARISON_OPERATORS))
            )
        )
    return column, operator, value


def parse_rank_spec(spec):
    parts = str(spec).split(':', 1)
    if len(parts) != 2:
        raise ValueError(
            "Invalid --rank specification: '{}'. Expected COLUMN:asc or COLUMN:desc.".format(spec)
        )
    column, direction = parts
    column = column.strip()
    direction = direction.strip().lower()
    if column == '':
        raise ValueError("Invalid --rank specification: column name is empty.")
    if direction not in RANK_DIRECTIONS:
        raise ValueError(
            "Invalid --rank direction '{}'. Supported directions: asc, desc.".format(direction)
        )
    return column, direction


def _require_columns(dataframe, columns, option_name):
    missing = [column for column in columns if column not in dataframe.columns]
    if missing:
        raise ValueError(
            "{} references missing column(s): {}".format(option_name, ', '.join(missing))
        )


def _numeric_series(dataframe, column, option_name):
    numeric = pd.to_numeric(dataframe[column], errors='coerce')
    bad_values = dataframe.loc[dataframe[column].notna() & numeric.isna(), column].astype(str).unique().tolist()
    if bad_values:
        raise ValueError(
            "{} requires numeric values in column '{}'. Non-numeric value(s): {}".format(
                option_name, column, ', '.join(sorted(bad_values))
            )
        )
    return numeric


def apply_filters(dataframe, filter_specs):
    parsed_filters = [parse_filter_spec(spec) for spec in (filter_specs or [])]
    _require_columns(dataframe, [column for column, _operator, _value in parsed_filters], '--filter')
    keep = pd.Series(True, index=dataframe.index)
    for column, operator, raw_value in parsed_filters:
        if operator in NUMERIC_OPERATORS:
            try:
                threshold = float(raw_value)
            except ValueError as exc:
                raise ValueError(
                    "Numeric --filter operator '{}' requires a numeric threshold: {}".format(
                        operator, raw_value
                    )
                ) from exc
            values = _numeric_series(dataframe, column, '--filter')
            present = values.notna()
            if operator == 'ge':
                current = present & (values >= threshold)
            elif operator == 'gt':
                current = present & (values > threshold)
            elif operator == 'le':
                current = present & (values <= threshold)
            else:
                current = present & (values < threshold)
        else:
            values = dataframe[column]
            present = values.notna()
            values_as_text = values.astype(str)
            if operator == 'eq':
                current = present & (values_as_text == str(raw_value))
            else:
                current = present & (values_as_text != str(raw_value))
        keep &= current
    return dataframe.loc[keep].copy()


def sort_candidates(dataframe, rank_specs):
    parsed_ranks = [parse_rank_spec(spec) for spec in (rank_specs or [])]
    _require_columns(dataframe, [column for column, _direction in parsed_ranks], '--rank')
    if dataframe.empty:
        return dataframe.copy()

    sorted_df = dataframe.copy()
    helper_columns = []
    sort_columns = []
    ascending = []
    for index, (column, direction) in enumerate(parsed_ranks):
        helper_column = '__nwkit_rank_{}_{}'.format(index, column)
        numeric_values = pd.to_numeric(sorted_df[column], errors='coerce')
        non_missing = sorted_df[column].notna()
        if non_missing.any() and numeric_values[non_missing].notna().all():
            sorted_df[helper_column] = numeric_values
        else:
            sorted_df[helper_column] = sorted_df[column].astype(str)
        helper_columns.append(helper_column)
        sort_columns.append(helper_column)
        ascending.append(direction == 'asc')

    sort_columns.append('leaf_name')
    ascending.append(True)
    sorted_df = sorted_df.sort_values(
        by=sort_columns,
        ascending=ascending,
        kind='mergesort',
        na_position='last',
    )
    return sorted_df.drop(columns=helper_columns)


def _edge_length(node):
    if node.is_root:
        return 0.0
    if node.dist is None:
        return 1.0
    try:
        return float(node.dist)
    except (TypeError, ValueError):
        return 1.0


def _leaf_path_edges(leaf):
    path_edges = []
    node = leaf
    while not node.is_root:
        path_edges.append((node, _edge_length(node)))
        node = node.up
    return path_edges


def _path_gain(path_edges, covered_edges):
    return sum(length for edge, length in path_edges if edge not in covered_edges)


def select_ranked(candidate_order, path_edges_by_leaf, n):
    covered_edges = set()
    selected = []
    pd_total = 0.0
    for leaf_name in candidate_order[:n]:
        path_edges = path_edges_by_leaf[leaf_name]
        gain = _path_gain(path_edges, covered_edges)
        for edge, _length in path_edges:
            covered_edges.add(edge)
        pd_total += gain
        selected.append(
            {
                'leaf_name': leaf_name,
                'pd_gain': gain,
                'pd_total': pd_total,
            }
        )
    return selected


def select_max_pd(candidate_order, path_edges_by_leaf, n):
    remaining = set(candidate_order)
    rank_index = {leaf_name: index for index, leaf_name in enumerate(candidate_order)}
    covered_edges = set()
    selected = []
    pd_total = 0.0
    for _index in range(n):
        best_leaf = None
        best_gain = None
        for leaf_name in candidate_order:
            if leaf_name not in remaining:
                continue
            gain = _path_gain(path_edges_by_leaf[leaf_name], covered_edges)
            if (
                best_leaf is None
                or gain > best_gain
                or (math.isclose(gain, best_gain) and rank_index[leaf_name] < rank_index[best_leaf])
            ):
                best_leaf = leaf_name
                best_gain = gain
        if best_leaf is None:
            break
        path_edges = path_edges_by_leaf[best_leaf]
        for edge, _length in path_edges:
            covered_edges.add(edge)
        pd_total += best_gain
        remaining.remove(best_leaf)
        selected.append(
            {
                'leaf_name': best_leaf,
                'pd_gain': best_gain,
                'pd_total': pd_total,
            }
        )
    return selected


def _build_output_table(selected_rows, candidate_df):
    selected_df = pd.DataFrame(selected_rows)
    selected_df.insert(0, 'sample_order', range(1, len(selected_df.index) + 1))
    metadata_columns = [column for column in candidate_df.columns if column != 'leaf_name']
    if metadata_columns:
        selected_df = selected_df.merge(
            candidate_df[['leaf_name', *metadata_columns]],
            on='leaf_name',
            how='left',
            sort=False,
        )
    return selected_df


def sample_main(args):
    if args.n < 1:
        raise ValueError("'--n' must be a positive integer.")
    if args.output_table in ['', '-']:
        raise ValueError("'--output-table' must be a file path, not '-' or an empty string.")

    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    validate_unique_named_leaves(tree, option_name='--infile', context=" for 'sample'")
    trait_df = read_sample_trait(args, tree)
    filtered_df = apply_filters(trait_df, args.filter)
    if filtered_df.empty:
        raise ValueError('No leaves passed the requested --filter condition(s).')

    ranked_df = sort_candidates(filtered_df, args.rank)
    candidate_order = ranked_df['leaf_name'].tolist()
    if len(candidate_order) < args.n and not args.allow_fewer:
        raise ValueError(
            "Only {} leaves passed the requested filters, fewer than --n {}. "
            "Adjust --filter/--n or set --allow-fewer yes.".format(len(candidate_order), args.n)
        )
    n_select = min(args.n, len(candidate_order))

    leaf_by_name = {leaf.name: leaf for leaf in tree.leaves()}
    path_edges_by_leaf = {
        leaf_name: _leaf_path_edges(leaf_by_name[leaf_name])
        for leaf_name in candidate_order
    }

    if args.method == 'max-pd':
        selected_rows = select_max_pd(candidate_order, path_edges_by_leaf, n_select)
    elif args.method == 'ranked':
        selected_rows = select_ranked(candidate_order, path_edges_by_leaf, n_select)
    else:
        raise ValueError("Unsupported --method: {}".format(args.method))

    selected_names = [row['leaf_name'] for row in selected_rows]
    if not selected_names:
        raise ValueError('No leaves were selected for output.')
    sys.stderr.write(''.join('Selected leaf: {}\n'.format(name) for name in selected_names))

    if args.output_table is not None:
        output_df = _build_output_table(selected_rows, ranked_df)
        output_table_dir = os.path.dirname(os.path.realpath(args.output_table))
        if output_table_dir:
            os.makedirs(output_table_dir, exist_ok=True)
        output_df.to_csv(args.output_table, sep='\t', index=False)

    tree.prune(selected_names, preserve_branch_length=True)
    write_tree(tree, args, format=args.outformat)
