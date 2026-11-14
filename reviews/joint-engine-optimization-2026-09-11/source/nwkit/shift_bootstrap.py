"""Validate external bootstrap configurations and summarize selection support."""

import csv
from collections import Counter
from typing import Any


def _table(directory, name, columns):
    with (directory / name).open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t", strict=True)
        if reader.fieldnames != columns:
            raise ValueError(f"Invalid kfl1ou {name} columns.")
        try:
            rows = list(reader)
        except csv.Error as error:
            raise ValueError(f"Invalid kfl1ou {name} quoting.") from error
    if any(None in row or any(v is None for v in row.values()) for row in rows):
        raise ValueError(f"Invalid kfl1ou {name} row.")
    return rows


def _count(value):
    if not value.isascii() or not value.isdecimal():
        raise ValueError("Invalid kfl1ou bootstrap count.")
    return int(value)


def _partition(tree, ids, selected):
    regimes: dict[Any, int] = {}
    groups: dict[int, list[str]] = {}
    for node in tree.traverse("preorder"):
        regimes[node] = (
            ids[node] if node.is_root or ids[node] in selected else regimes[node.up]
        )
        if node.is_leaf:
            groups.setdefault(regimes[node], []).append(node.name)
    return tuple(sorted(tuple(sorted(names)) for names in groups.values()))


def collect_bootstrap(directory, args, tree, ids, mapping):
    """Frequencies are conditional on successful refits; no posterior claims."""
    if args.bootstrap == 0:
        return None
    rows = _table(directory, "bootstrap.tsv", ["attempted", "successful", "failed"])
    if len(rows) != 1:
        raise ValueError("Expected one kfl1ou bootstrap summary.")
    counts = {key: _count(value) for key, value in rows[0].items()}
    if (
        counts["attempted"] != args.bootstrap
        or counts["successful"] < 1
        or counts["successful"] + counts["failed"] != args.bootstrap
    ):
        raise ValueError("Inconsistent kfl1ou bootstrap counts.")
    configurations = []
    for index, row in enumerate(
        _table(directory, "bootstrap-configurations.tsv", ["success_index", "clades"]),
        1,
    ):
        keys = row["clades"].split(";") if row["clades"] else []
        if (
            _count(row["success_index"]) != index
            or len(keys) != len(set(keys))
            or len(keys) > args.max_shifts
            or any(key not in mapping or mapping[key] == 0 for key in keys)
        ):
            raise ValueError("Invalid kfl1ou bootstrap configuration.")
        configurations.append(tuple(sorted(mapping[key] for key in keys)))
    successful = counts["successful"]
    if len(configurations) != successful:
        raise ValueError("Bootstrap configurations disagree with successful count.")
    failures = _table(directory, "bootstrap-failures.tsv", ["message", "count"])
    failure_rows = [
        {"message": row["message"], "count": _count(row["count"])} for row in failures
    ]
    if (
        any(not row["message"] or row["count"] < 1 for row in failure_rows)
        or len({row["message"] for row in failure_rows}) != len(failure_rows)
        or sum(row["count"] for row in failure_rows) > counts["failed"]
    ):
        raise ValueError("Inconsistent kfl1ou bootstrap failure messages.")
    inclusion = Counter(branch for config in configurations for branch in config)
    exact = Counter(configurations)
    partitions = Counter(
        _partition(tree, ids, set(config)) for config in configurations
    )
    result = {
        "type": "parametric",
        "seed": args.bootstrap_seed,
        **counts,
        "frequency_denominator": "successful_refits",
        "failure_messages": failure_rows,
        "failures_without_message": counts["failed"]
        - sum(row["count"] for row in failure_rows),
        "successful_configurations": [list(config) for config in configurations],
        "edge_inclusion": [
            {
                "branch_id": branch,
                "count": inclusion[branch],
                "frequency": inclusion[branch] / successful,
            }
            for branch in sorted(ids.values())
            if branch != 0
        ],
        "configuration_frequencies": [
            {
                "shift_branch_ids": list(config),
                "count": count,
                "frequency": count / successful,
            }
            for config, count in sorted(
                exact.items(), key=lambda item: (-item[1], item[0])
            )
        ],
        "tip_partition_frequencies": [
            {
                "groups": [list(group) for group in partition],
                "count": count,
                "frequency": count / successful,
            }
            for partition, count in sorted(
                partitions.items(), key=lambda item: (-item[1], item[0])
            )
        ],
    }
    if args.convergence:
        from nwkit.shift_convergence import convergence_bootstrap

        result.update(
            convergence_bootstrap(directory, configurations, mapping, tree, ids)
        )
    return result
