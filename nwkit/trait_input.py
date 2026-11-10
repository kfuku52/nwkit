"""Shared continuous trait-column parsing for comparative summaries."""

import math

import numpy as np

from nwkit.conventions import DEFAULT_TABLE_MISSING_VALUES
from nwkit.util import is_missing_table_value, parse_table_missing_values


def parse_trait_columns(value):
    columns = [part.strip() for part in value.split(",")]
    if not all(columns) or len(set(columns)) != len(columns) or "leaf_name" in columns:
        raise ValueError(
            "Trait/SE columns must be unique, nonempty and different from leaf_name."
        )
    return columns


def numeric_trait_value(value, column, missing):
    if is_missing_table_value(value, missing):
        return np.nan
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"Column '{column}' must contain continuous numeric values."
        ) from exc
    if not math.isfinite(result):
        raise ValueError(f"Column '{column}' contains a non-finite value.")
    return result


def numeric_trait_matrix(table, columns, missing_values):
    """Parse selected columns in the table's current tip order."""
    markers = parse_table_missing_values(missing_values)
    if markers == set(DEFAULT_TABLE_MISSING_VALUES) and all(
        table[column].dtype.kind in "biuf" for column in columns
    ):
        values = np.ascontiguousarray(
            table[columns].to_numpy(dtype=float, na_value=np.nan)
        )
        infinite = np.isinf(values)
        if infinite.any():
            column = columns[int(infinite.argmax()) % len(columns)]
            raise ValueError(f"Column '{column}' contains a non-finite value.")
        return values

    rows = zip(*(table[column].array for column in columns), strict=True)
    return np.asarray(
        [
            [
                numeric_trait_value(value, column, markers)
                for column, value in zip(columns, row, strict=True)
            ]
            for row in rows
        ]
    )
