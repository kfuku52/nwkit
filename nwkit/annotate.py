import math
import sys
from collections import Counter
from typing import Any

import pandas as pd

from nwkit.transfer import _validate_property_name
from nwkit.util import (
    assign_branch_ids,
    get_node_class,
    get_tree_property_names,
    is_missing_table_value,
    parse_table_missing_values,
    read_tip_table,
    read_tree,
    validate_distinct_output_paths,
    validate_unique_named_leaves,
    write_tree,
)

ANNOTATION_REPORT_COLUMNS = (
    "branch_id",
    "node_class",
    "descendant_taxa",
    "input_column",
    "property",
    "aggregation",
    "value",
    "status",
    "reason",
)
SUPPORTED_AGGREGATIONS = (
    "unique",
    "mode",
    "count",
    "mean",
    "sum",
    "min",
    "max",
    "list",
)


def _format_taxa(node):
    return ",".join(sorted(str(name) for name in node.leaf_names()))


def _python_value(value, missing_values):
    if is_missing_table_value(value, missing_values):
        return None
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, (str, int, float, bool)):
        return value
    return str(value)


def _parse_column_specs(columns, property_maps, table_columns):
    specs = list()
    if columns:
        for column in str(columns).split(","):
            column = column.strip()
            if column:
                specs.append((column, _validate_property_name(column)))
    for raw_mapping in property_maps or []:
        if "=" not in str(raw_mapping):
            raise ValueError(
                "--property-map must use COLUMN=PROPERTY syntax: {}".format(raw_mapping)
            )
        column, prop = str(raw_mapping).split("=", 1)
        specs.append((column.strip(), _validate_property_name(prop.strip())))
    if not specs:
        specs = [
            (str(column), _validate_property_name(str(column)))
            for column in table_columns
            if str(column) != "leaf_name"
        ]
    seen_properties = set()
    for column, prop in specs:
        if column not in table_columns:
            raise ValueError("Column '{}' was not found in --table.".format(column))
        if column == "leaf_name":
            raise ValueError(
                "The 'leaf_name' key column cannot be attached as a property."
            )
        if prop in seen_properties:
            raise ValueError(
                "Output property is specified more than once: {}".format(prop)
            )
        seen_properties.add(prop)
    return specs


def _parse_aggregation(raw, table_columns):
    parts = str(raw).split(":")
    if len(parts) not in (2, 3):
        raise ValueError(
            "--aggregate must use COLUMN:METHOD[:PROPERTY] syntax: {}".format(raw)
        )
    column, method = parts[0].strip(), parts[1].strip().lower()
    if column not in table_columns:
        raise ValueError(
            "Aggregation column '{}' was not found in --table.".format(column)
        )
    if method not in SUPPORTED_AGGREGATIONS:
        raise ValueError(
            "Unsupported aggregation '{}'. Choose from: {}.".format(
                method,
                ", ".join(SUPPORTED_AGGREGATIONS),
            )
        )
    prop = parts[2].strip() if len(parts) == 3 else "{}_{}".format(column, method)
    return column, method, _validate_property_name(prop)


def _aggregate(values, method, column):
    if method == "count":
        return len(values), True, "count_non_missing"
    if not values:
        return None, False, "no_non_missing_values"
    unique_by_text: dict[str, Any] = {}
    for value in values:
        unique_by_text.setdefault(str(value), value)
    if method == "unique":
        if len(unique_by_text) == 1:
            return next(iter(unique_by_text.values())), True, "single_unique_value"
        return None, False, "multiple_values"
    if method == "mode":
        counts = Counter(str(value) for value in values)
        winning_text = sorted(counts, key=lambda text: (-counts[text], text))[0]
        return unique_by_text[winning_text], True, "deterministic_mode"
    if method == "list":
        return "|".join(sorted(unique_by_text)), True, "sorted_unique_values"
    try:
        numeric = pd.to_numeric(pd.Series(values), errors="raise")
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "Aggregation '{}' requires numeric values in column '{}'.".format(
                method, column
            )
        ) from exc
    if method == "mean":
        return float(numeric.mean()), True, "numeric_mean"
    if method == "sum":
        return float(numeric.sum()), True, "numeric_sum"
    if method == "min":
        return float(numeric.min()), True, "numeric_minimum"
    if method == "max":
        return float(numeric.max()), True, "numeric_maximum"
    raise AssertionError("Unhandled aggregation method: {}".format(method))


def _aggregation_requirements(aggregations):
    methods_by_column: dict[str, set[str]] = {}
    for column, method, _prop in aggregations:
        methods_by_column.setdefault(column, set()).add(method)
    return {
        column: {
            "numeric": bool(methods.intersection({"mean", "sum", "min", "max"})),
            "unique": bool(methods.intersection({"unique", "mode", "list"})),
        }
        for column, methods in methods_by_column.items()
    }


def _stable_numeric_sum(values):
    # Keep integer columns exact until the final public float conversion.
    if all(isinstance(value, int) and not isinstance(value, bool) for value in values):
        return sum(values)
    float_values = [float(value) for value in values]
    if not all(math.isfinite(value) for value in float_values):
        return sum(float_values)
    scale = max((abs(value) for value in float_values), default=0.0)
    if scale == 0.0:
        return 0.0
    return scale * math.fsum(value / scale for value in float_values)


def _merge_aggregation_summaries(child_summaries, requirements):
    count = sum(summary["count"] for summary in child_summaries)
    summary: dict[str, Any] = {"count": count}
    if requirements["numeric"]:
        summary["sum"] = _stable_numeric_sum(
            [child_summary["sum"] for child_summary in child_summaries]
        )
        populated = [
            child_summary
            for child_summary in child_summaries
            if child_summary["count"] > 0
        ]
        summary["min"] = (
            min(child_summary["min"] for child_summary in populated)
            if populated
            else None
        )
        summary["max"] = (
            max(child_summary["max"] for child_summary in populated)
            if populated
            else None
        )
    if requirements["unique"]:
        child_maps = [child_summary["unique"] for child_summary in child_summaries]
        unique: dict[str, tuple[int, Any, int]] = max(child_maps, key=len, default={})
        for child_map in child_maps:
            if child_map is unique:
                continue
            for text, (first_index, value, value_count) in child_map.items():
                previous = unique.get(text)
                if previous is None:
                    unique[text] = (first_index, value, value_count)
                else:
                    previous_index, previous_value, previous_count = previous
                    if first_index < previous_index:
                        previous_index = first_index
                        previous_value = value
                    unique[text] = (
                        previous_index,
                        previous_value,
                        previous_count + value_count,
                    )
        summary["unique"] = unique
    return summary


def _aggregate_summary(summary, method, column):
    if method == "count":
        return summary["count"], True, "count_non_missing"
    if summary["count"] == 0:
        return None, False, "no_non_missing_values"
    if method == "unique":
        if len(summary["unique"]) == 1:
            return (
                next(iter(summary["unique"].values()))[1],
                True,
                "single_unique_value",
            )
        return None, False, "multiple_values"
    if method == "mode":
        winning_text = min(
            summary["unique"],
            key=lambda text: (-summary["unique"][text][2], text),
        )
        return summary["unique"][winning_text][1], True, "deterministic_mode"
    if method == "list":
        return "|".join(sorted(summary["unique"])), True, "sorted_unique_values"
    if method == "mean":
        return float(summary["sum"] / summary["count"]), True, "numeric_mean"
    if method == "sum":
        return float(summary["sum"]), True, "numeric_sum"
    if method == "min":
        return float(summary["min"]), True, "numeric_minimum"
    if method == "max":
        return float(summary["max"]), True, "numeric_maximum"
    raise AssertionError(
        "Unhandled aggregation method for '{}': {}".format(column, method)
    )


def _write_report(rows, path):
    if path in (None, ""):
        return
    if path == "-":
        raise ValueError("'--report' requires a file path, not '-'.")
    pd.DataFrame(rows, columns=ANNOTATION_REPORT_COLUMNS).to_csv(
        path, sep="\t", index=False
    )


def _annotate_leaf_properties(
    tree,
    rows_by_leaf,
    column_specs,
    missing_values,
    collect_report,
    node_ids,
    report_rows,
):
    annotated = 0
    value_by_leaf_and_column = {}
    for leaf in tree.leaves():
        leaf_name = str(leaf.name)
        table_row = rows_by_leaf.get(leaf_name)
        for column, prop in column_specs:
            value = None
            if table_row is None:
                status, reason = "unmatched", "tree_tip_absent_from_table"
            else:
                value = _python_value(table_row[column], missing_values)
                value_by_leaf_and_column[(leaf_name, column)] = value
                if value is None:
                    status, reason = "missing", "missing_table_value"
                else:
                    leaf.props[prop] = value
                    annotated += 1
                    status, reason = "annotated", "matching_leaf_name"
            if collect_report:
                report_rows.append(
                    {
                        "branch_id": node_ids[leaf],
                        "node_class": "leaf",
                        "descendant_taxa": leaf_name,
                        "input_column": column,
                        "property": prop,
                        "aggregation": "",
                        "value": value,
                        "status": status,
                        "reason": reason,
                    }
                )
    return annotated, value_by_leaf_and_column


def _add_aggregation_input_values(
    value_by_leaf_and_column, rows_by_leaf, aggregations, missing_values
):
    for leaf_name, table_row in rows_by_leaf.items():
        for column, _method, _prop in aggregations:
            value_by_leaf_and_column[(leaf_name, column)] = _python_value(
                table_row[column], missing_values
            )


def _numeric_aggregation_values(
    value_by_leaf_and_column,
    requirements_by_column,
    aggregations,
    leaf_order,
    has_aggregation_target,
):
    numeric_values = {}
    for column, requirements in requirements_by_column.items():
        if not requirements["numeric"] or not has_aggregation_target:
            continue
        leaf_names_and_values = [
            (leaf_name, value)
            for (leaf_name, value_column), value in value_by_leaf_and_column.items()
            if value_column == column and leaf_name in leaf_order and value is not None
        ]
        try:
            converted = pd.to_numeric(
                pd.Series([value for _, value in leaf_names_and_values]),
                errors="raise",
            )
        except (TypeError, ValueError) as exc:
            method = next(
                method
                for aggregate_column, method, _prop in aggregations
                if aggregate_column == column
                and method in {"mean", "sum", "min", "max"}
            )
            raise ValueError(
                "Aggregation '{}' requires numeric values in column '{}'.".format(
                    method, column
                )
            ) from exc
        for (leaf_name, _value), numeric in zip(
            leaf_names_and_values, converted, strict=True
        ):
            numeric_values[(leaf_name, column)] = (
                numeric.item() if hasattr(numeric, "item") else numeric
            )
    return numeric_values


def _leaf_aggregation_summary(
    node, column, requirements, value_by_leaf_and_column, numeric_values, leaf_order
):
    leaf_name = str(node.name)
    value = value_by_leaf_and_column.get((leaf_name, column))
    summary: dict[str, Any] = {"count": 0 if value is None else 1}
    if requirements["numeric"]:
        numeric = numeric_values.get((leaf_name, column))
        summary.update(
            {
                "sum": 0 if numeric is None else numeric,
                "min": numeric,
                "max": numeric,
            }
        )
    if requirements["unique"]:
        summary["unique"] = (
            {} if value is None else {str(value): (leaf_order[leaf_name], value, 1)}
        )
    return summary


def _aggregate_tree_nodes(
    tree,
    aggregations,
    requirements_by_column,
    value_by_leaf_and_column,
    numeric_values,
    leaf_order,
    collect_report,
    node_ids,
    report_rows,
):
    summaries_by_column: dict[str, dict[Any, dict[str, Any]]] = {
        column: {} for column in requirements_by_column
    }
    annotated = 0
    for node in tree.traverse(strategy="postorder"):
        for column, requirements in requirements_by_column.items():
            if node.is_leaf:
                summary = _leaf_aggregation_summary(
                    node,
                    column,
                    requirements,
                    value_by_leaf_and_column,
                    numeric_values,
                    leaf_order,
                )
            else:
                summary = _merge_aggregation_summaries(
                    [
                        summaries_by_column[column][child]
                        for child in node.get_children()
                    ],
                    requirements,
                )
            summaries_by_column[column][node] = summary
        if node.is_leaf:
            continue
        for column, method, prop in aggregations:
            value, assign, reason = _aggregate_summary(
                summaries_by_column[column][node], method, column
            )
            status = "aggregated" if assign else "not_aggregated"
            if assign:
                node.props[prop] = value
                annotated += 1
            if collect_report:
                report_rows.append(
                    {
                        "branch_id": node_ids[node],
                        "node_class": get_node_class(node),
                        "descendant_taxa": _format_taxa(node),
                        "input_column": column,
                        "property": prop,
                        "aggregation": method,
                        "value": value,
                        "status": status,
                        "reason": reason,
                    }
                )
    return annotated


def _append_unmatched_table_rows(
    report_rows, table_only, column_specs, rows_by_leaf, missing_values
):
    for leaf_name in table_only:
        for column, prop in column_specs:
            report_rows.append(
                {
                    "branch_id": "",
                    "node_class": "",
                    "descendant_taxa": leaf_name,
                    "input_column": column,
                    "property": prop,
                    "aggregation": "",
                    "value": _python_value(
                        rows_by_leaf[leaf_name][column], missing_values
                    ),
                    "status": "unmatched",
                    "reason": "table_row_absent_from_tree",
                }
            )


def annotate_main(args):
    if getattr(args, "table", None) in (None, ""):
        raise ValueError("'--table' is required for 'annotate'.")
    validate_distinct_output_paths(
        [
            ("--outfile", getattr(args, "outfile", None)),
            ("--report", getattr(args, "report", None)),
        ]
    )
    tree = read_tree(
        args.infile,
        args.format,
        args.quoted_node_names,
        rooted=getattr(args, "input_rooted", "auto"),
    )
    validate_unique_named_leaves(
        tree, option_name="--infile", context=" for 'annotate'"
    )
    table, table_only, tree_only = read_tip_table(
        args.table,
        option_name="--table",
        tree_leaf_names=tree.leaf_names(),
        unmatched=getattr(args, "unmatched", "warn"),
        missing_values=getattr(args, "missing_values", None),
    )
    missing_values = parse_table_missing_values(getattr(args, "missing_values", None))
    column_specs = _parse_column_specs(
        columns=getattr(args, "columns", None),
        property_maps=getattr(args, "property_map", None),
        table_columns=list(table.columns),
    )
    aggregations = [
        _parse_aggregation(raw, list(table.columns))
        for raw in (getattr(args, "aggregate", None) or [])
    ]
    output_properties = set(get_tree_property_names(tree))
    output_properties.update(prop for _, prop in column_specs)
    output_properties.update(prop for _, _, prop in aggregations)
    collect_report = getattr(args, "report", None) not in (None, "")
    rows_by_leaf = {str(row["leaf_name"]): row for _, row in table.iterrows()}
    node_ids = assign_branch_ids(tree) if collect_report else {}
    report_rows: list[dict[str, Any]] = []
    annotated, value_by_leaf_and_column = _annotate_leaf_properties(
        tree,
        rows_by_leaf,
        column_specs,
        missing_values,
        collect_report,
        node_ids,
        report_rows,
    )
    _add_aggregation_input_values(
        value_by_leaf_and_column, rows_by_leaf, aggregations, missing_values
    )
    requirements_by_column = _aggregation_requirements(aggregations)
    has_aggregation_target = any(not node.is_leaf for node in tree.traverse())
    leaf_order = {str(leaf.name): index for index, leaf in enumerate(tree.leaves())}
    numeric_values = _numeric_aggregation_values(
        value_by_leaf_and_column,
        requirements_by_column,
        aggregations,
        leaf_order,
        has_aggregation_target,
    )
    annotated += _aggregate_tree_nodes(
        tree,
        aggregations,
        requirements_by_column,
        value_by_leaf_and_column,
        numeric_values,
        leaf_order,
        collect_report,
        node_ids,
        report_rows,
    )
    if collect_report:
        _append_unmatched_table_rows(
            report_rows, table_only, column_specs, rows_by_leaf, missing_values
        )
    _write_report(report_rows, getattr(args, "report", None))
    sys.stderr.write("Attached or aggregated {} property values.\n".format(annotated))
    write_tree(tree, args, format=args.outformat, props=output_properties)
