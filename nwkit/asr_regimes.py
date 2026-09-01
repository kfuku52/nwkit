"""Validated branch-regime inputs shared by discrete and continuous ASR."""

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from nwkit.util import assign_branch_ids


@dataclass(frozen=True)
class RegimeAssignment:
    regimes: tuple[str, ...]
    by_node: dict[Any, str]
    source: str

    @property
    def root_regime(self):
        return next(regime for node, regime in self.by_node.items() if node.is_root)


def _read_tsv(path, option_name):
    source = Path(str(path))
    try:
        return source, pd.read_csv(source, sep="\t", dtype=str, keep_default_na=False)
    except (OSError, pd.errors.ParserError) as exc:
        raise ValueError(f"Failed to read '{option_name}': {source}") from exc


def read_regime_map(path, tree):
    """Read one regime for every node/incoming branch, including the root."""

    if path in (None, ""):
        raise ValueError("The selected model requires --regime-map.")
    source, table = _read_tsv(path, "--regime-map")
    required = {"branch_id", "regime"}
    if not required.issubset(table.columns):
        raise ValueError(
            "'--regime-map' must be a TSV with branch_id and regime columns."
        )
    if len(table.columns) != len(set(table.columns)):
        raise ValueError("'--regime-map' contains duplicated column names.")
    extras = set(table.columns) - required
    if extras:
        raise ValueError(
            "Unsupported '--regime-map' column(s): " + ", ".join(sorted(extras))
        )
    node_ids = assign_branch_ids(tree)
    nodes_by_id = {branch_id: node for node, branch_id in node_ids.items()}
    by_node = {}
    regimes = []
    seen_regimes = set()
    for row_number, row in enumerate(table.itertuples(index=False), start=2):
        raw_branch_id = str(row.branch_id).strip()
        try:
            branch_id = int(raw_branch_id)
        except ValueError as exc:
            raise ValueError(
                f"Invalid branch_id in '--regime-map' row {row_number}: {raw_branch_id}"
            ) from exc
        if str(branch_id) != raw_branch_id and raw_branch_id not in {f"+{branch_id}"}:
            raise ValueError(
                f"Invalid branch_id in '--regime-map' row {row_number}: {raw_branch_id}"
            )
        if branch_id not in nodes_by_id:
            raise ValueError(
                f"Unknown branch_id in '--regime-map' row {row_number}: {branch_id}"
            )
        node = nodes_by_id[branch_id]
        if node in by_node:
            raise ValueError(
                f"Duplicated branch_id in '--regime-map': {branch_id}"
            )
        regime = str(row.regime).strip()
        if regime == "":
            raise ValueError(
                f"Empty regime in '--regime-map' row {row_number}."
            )
        by_node[node] = regime
        if regime not in seen_regimes:
            regimes.append(regime)
            seen_regimes.add(regime)
    missing = sorted(set(nodes_by_id) - {node_ids[node] for node in by_node})
    if missing:
        raise ValueError(
            "'--regime-map' must assign every branch_id, including root 0; missing: "
            + ",".join(str(value) for value in missing)
        )
    return RegimeAssignment(tuple(regimes), by_node, str(source))


def read_regime_parameters(path, regimes, columns):
    """Read complete fixed regime parameters, preserving regime order."""

    if path in (None, ""):
        return None
    source, table = _read_tsv(path, "--regime-parameters")
    required = {"regime", *columns}
    if not required.issubset(table.columns):
        raise ValueError(
            "'--regime-parameters' must contain columns: "
            + ", ".join(("regime", *columns))
        )
    extras = set(table.columns) - required
    if extras:
        raise ValueError(
            "Unsupported '--regime-parameters' column(s): "
            + ", ".join(sorted(extras))
        )
    records = {}
    for row_number, row in enumerate(table.itertuples(index=False), start=2):
        regime = str(row.regime).strip()
        if regime == "" or regime in records:
            problem = "empty" if regime == "" else f"duplicated '{regime}'"
            raise ValueError(
                f"'--regime-parameters' has an {problem} regime at row {row_number}."
            )
        values = {}
        for column in columns:
            raw = str(getattr(row, column)).strip()
            try:
                value = float(raw)
            except ValueError as exc:
                raise ValueError(
                    f"Invalid {column} for regime '{regime}': {raw}"
                ) from exc
            if not math.isfinite(value):
                raise ValueError(
                    f"Invalid {column} for regime '{regime}': values must be finite."
                )
            values[column] = value
        records[regime] = values
    expected, observed = set(regimes), set(records)
    if expected != observed:
        missing = sorted(expected - observed)
        unknown = sorted(observed - expected)
        details = []
        if missing:
            details.append("missing " + ",".join(missing))
        if unknown:
            details.append("unknown " + ",".join(unknown))
        raise ValueError(
            "'--regime-parameters' regimes do not match '--regime-map': "
            + "; ".join(details)
        )
    return {
        regime: {column: records[regime][column] for column in columns}
        for regime in regimes
    }, str(source)
