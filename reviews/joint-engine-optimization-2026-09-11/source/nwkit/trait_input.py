"""Shared continuous trait-column parsing for comparative summaries."""

import numpy as np

from nwkit.util import is_missing_table_value


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
    if not np.isfinite(result):
        raise ValueError(f"Column '{column}' contains a non-finite value.")
    return result
