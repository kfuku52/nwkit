"""External species-age intervals are display evidence, never implicit bounds."""

import numpy as np
import pandas as pd

INPUT_INTERVAL_COLUMNS = (
    "input_interval_lower",
    "input_interval_upper",
    "input_interval_level",
    "input_interval_kind",
    "input_interval_source",
)


def attach_species_intervals(chronology, path):
    table = chronology.species_table
    for key in INPUT_INTERVAL_COLUMNS:
        table[key] = np.nan if key.endswith(("lower", "upper", "level")) else ""
    if not path:
        return
    records = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {"lower", "upper", "level", "kind", "source"}
    if not required.issubset(records.columns) or not (
        {"node", "species_event_id"} & set(records.columns)
    ):
        raise ValueError(
            "Species intervals require node or species_event_id, plus lower, upper, level, kind, source."
        )
    seen = set()
    for row in records.to_dict("records"):
        by_id = bool(row.get("species_event_id"))
        key, value = (
            ("species_event_id", row["species_event_id"])
            if by_id
            else ("node", row.get("node", ""))
        )
        matches = table.index[table[key] == value]
        if not value or len(matches) != 1:
            raise ValueError(f"Unknown or ambiguous species interval node: {value}")
        index = matches[0]
        if index in seen:
            raise ValueError("Duplicate species interval node.")
        seen.add(index)
        if by_id and row.get("node") and table.loc[index, "node"] != row["node"]:
            raise ValueError("Species interval node and clade ID disagree.")
        lo, hi, level = float(row["lower"]), float(row["upper"]), float(row["level"])
        if (
            not np.isfinite([lo, hi, level]).all()
            or not 0 <= lo <= hi
            or not 0 < level < 1
        ):
            raise ValueError(
                "Species intervals require finite 0 <= lower <= upper and 0 < level < 1."
            )
        if (
            row["kind"] not in {"confidence", "credible", "hpd", "percentile"}
            or not row["source"].strip()
        ):
            raise ValueError(
                "Species interval kind must be confidence, credible, hpd, or percentile; source is required."
            )
        if table.loc[index, "age"] == 0 and (lo != 0 or hi != 0):
            raise ValueError("Contemporaneous species tips require zero-age intervals.")
        table.loc[index, list(INPUT_INTERVAL_COLUMNS)] = [
            lo,
            hi,
            level,
            row["kind"],
            row["source"],
        ]
