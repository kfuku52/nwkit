#!/usr/bin/env python3
"""Pair before/after calibration results only when generated input hashes match."""

import argparse
import hashlib
from collections import defaultdict
from pathlib import Path

import numpy as np
from validate_regression_calibration import encode, finite, rate, read_records


def index(records, method):
    indexed = {}
    for task in records:
        matches = [row for row in task["results"] if row["method"] == method]
        if not matches:
            continue
        key = (task["case"]["name"], task["replicate"])
        if key in indexed or len(matches) != 1:
            raise ValueError(f"Duplicate record or method: {key}")
        indexed[key] = (task["input_sha256"], matches[0])
    return indexed


def available(row):
    return (
        row["status"] == "completed"
        and row.get("inference_status") == "ok"
        and finite(row.get("p_value"))
        and 0 <= float(row["p_value"]) <= 1
    )


def compare(before, after, method="wald"):
    before, after = index(before, method), index(after, method)
    grouped = defaultdict(list)
    for key, (digest, updated) in after.items():
        if key not in before or before[key][0] != digest:
            raise ValueError(f"Before/after input mismatch: {key}")
        grouped[key[0]].append((before[key][1], updated))
    report = []
    for case, pairs in sorted(grouped.items()):
        both = [(old, new) for old, new in pairs if available(old) and available(new)]
        changes = [
            int(float(new["p_value"]) < 0.05) - int(float(old["p_value"]) < 0.05)
            for old, new in both
        ]
        report.append(
            {
                "case": case,
                "method": method,
                "paired_datasets": len(pairs),
                "both_available": len(both),
                "became_available": sum(
                    not available(old) and available(new) for old, new in pairs
                ),
                "became_unavailable": sum(
                    available(old) and not available(new) for old, new in pairs
                ),
                "before_reported_rejections_all": rate(
                    sum(
                        available(old) and float(old["p_value"]) < 0.05
                        for old, _ in pairs
                    ),
                    len(pairs),
                ),
                "after_reported_rejections_all": rate(
                    sum(
                        available(new) and float(new["p_value"]) < 0.05
                        for _, new in pairs
                    ),
                    len(pairs),
                ),
                "rejection_difference_given_both_available": float(np.mean(changes))
                if changes
                else None,
                "paired_difference_mcse": float(
                    np.std(changes, ddof=1) / np.sqrt(len(changes))
                )
                if len(changes) > 1
                else None,
                "median_absolute_p_change_given_both_available": float(
                    np.median(
                        [
                            abs(float(new["p_value"]) - float(old["p_value"]))
                            for old, new in both
                        ]
                    )
                )
                if both
                else None,
            }
        )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", required=True, type=Path)
    parser.add_argument("--after", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    report = {
        "scope": "paired numerical-change assessment, not a claim of improved calibration",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "before_protocol_sha256": hashlib.sha256(
            (args.before / "protocol.json").read_bytes()
        ).hexdigest(),
        "after_protocol_sha256": hashlib.sha256(
            (args.after / "protocol.json").read_bytes()
        ).hexdigest(),
        "comparisons": compare(read_records(args.before), read_records(args.after)),
    }
    with args.output.open("x") as handle:
        handle.write(encode(report) + "\n")


if __name__ == "__main__":
    main()
