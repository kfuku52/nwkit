#!/usr/bin/env python3
"""Audit completed Pagel-lambda simulation records and their saved summary."""

import argparse
import hashlib
import json
import math
import tarfile
from collections import Counter, defaultdict
from pathlib import Path

from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parents[1]


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _source_bytes(path, archive):
    if archive is None:
        return (ROOT / path).read_bytes()
    with tarfile.open(archive, "r:gz") as saved:
        member = saved.extractfile(path)
        if member is None:
            raise ValueError(f"Missing archived source: {path}")
        return member.read()


def _rate(value, total, expected):
    _require(
        expected["count"] == value and expected["total"] == total,
        "Rate numerator or denominator mismatch",
    )
    if total == 0:
        _require(
            expected["rate"] is None and expected["mc95"] is None,
            "Undefined rate must have null estimate and interval",
        )
        return
    _require(
        math.isclose(expected["rate"], value / total, abs_tol=1e-12),
        "Rate estimate mismatch",
    )
    interval = binomtest(value, total).proportion_ci(method="wilson")
    _require(
        math.isclose(expected["mc95"][0], interval.low, abs_tol=1e-12)
        and math.isclose(expected["mc95"][1], interval.high, abs_tol=1e-12),
        "Wilson interval mismatch",
    )


def verify(directory, archive=None):
    protocol = json.loads((directory / "protocol.json").read_text())
    _require(protocol["completed"] is True, "Study is incomplete")
    scenarios = protocol["scenarios"]
    _require(len(scenarios) == len(set(scenarios)), "Duplicate scenario")
    for path, digest in protocol["source_sha256"].items():
        actual = hashlib.sha256(_source_bytes(path, archive)).hexdigest()
        _require(actual == digest, f"Source hash mismatch: {path}")
    rows = [
        json.loads(line)
        for line in (directory / "records.jsonl").read_text().splitlines()
    ]
    _require(
        len(rows) == protocol["outer"] * len(scenarios),
        "Missing or extra dataset record",
    )
    grouped = defaultdict(list)
    seen = set()
    for row in rows:
        name, replicate = row["scenario"], row["replicate"]
        _require(name in scenarios and type(replicate) is int, "Unexpected case")
        _require(0 <= replicate < protocol["outer"], "Replicate outside protocol")
        _require((name, replicate) not in seen, "Duplicate dataset")
        seen.add((name, replicate))
        expected_tips = (
            32
            if name == "balanced-32-null"
            else 6
            if name == "balanced-8-missing-null"
            else 8
        )
        _require(row["tips"] == expected_tips, "Tip count mismatch")
        _require(
            row["true_lambda"] == (0.6 if name.endswith("alternative") else 0.0),
            "Generating lambda mismatch",
        )
        _require(len(row["input_sha256"]) == 64, "Invalid input hash")
        if "chi2_p_value" in row:
            _require(0 <= row["chi2_p_value"] <= 1, "Invalid chi-square P-value")
        if "bootstrap_p_value" in row:
            p = row["bootstrap_p_value"]
            _require(
                1 / (protocol["inner"] + 1) <= p <= 1
                and math.isclose(
                    p * (protocol["inner"] + 1),
                    round(p * (protocol["inner"] + 1)),
                    abs_tol=1e-9,
                ),
                "Bootstrap P-value is outside its Monte Carlo grid",
            )
            _require("chi2_p_value" in row, "Bootstrap P-value without fit")
        if "ci_lower" in row:
            _require(
                0 <= row["ci_lower"] <= row["ci_upper"] <= 1,
                "Invalid profile interval",
            )
        grouped[name].append(row)
    summary = json.loads((directory / "summary.json").read_text())
    _require(
        [item["scenario"] for item in summary] == scenarios,
        "Summary scenario order or set mismatch",
    )
    for item in summary:
        rows = grouped[item["scenario"]]
        _require(
            len(rows) == item["generated"] == protocol["outer"],
            "Summary generated count mismatch",
        )
        _require(
            dict(Counter(row["status"] for row in rows)) == item["status_counts"],
            "Fit status counts mismatch",
        )
        errors = dict(
            Counter(row["bootstrap_error"] for row in rows if "bootstrap_error" in row)
        )
        _require(errors == item["bootstrap_failure_reasons"], "Failure counts mismatch")
        for method in ("chi2", "bootstrap"):
            available = [row for row in rows if f"{method}_p_value" in row]
            rejected = sum(
                row[f"{method}_p_value"] <= protocol["level"] for row in available
            )
            _rate(len(available), len(rows), item[f"{method}_p_available"])
            _rate(rejected, len(rows), item[f"{method}_rejection_all"])
            _rate(rejected, len(available), item[f"{method}_rejection_available"])
        intervals = [row for row in rows if "ci_lower" in row]
        if intervals:
            covered = sum(
                row["ci_lower"] <= row["true_lambda"] <= row["ci_upper"]
                for row in intervals
            )
            _rate(len(intervals), len(rows), item["profile_interval_available"])
            _rate(covered, len(rows), item["profile_interval_coverage_all"])
            _rate(covered, len(intervals), item["profile_interval_coverage_available"])
        else:
            _require(
                "profile_interval_available" not in item,
                "Unexpected profile summary",
            )
    return len(seen)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--source-archive", type=Path)
    args = parser.parse_args()
    count = verify(args.directory, args.source_archive)
    print(f"Verified {count} signal datasets in {args.directory}")


if __name__ == "__main__":
    main()
