#!/usr/bin/env python3
"""Audit retained paired RADTE sequence-interval comparison families."""

import argparse
import hashlib
import json
import math
import statistics
import tarfile
from collections import Counter
from pathlib import Path

from scipy.stats import binomtest

INPUTS = {
    "gene_tree": "gene.nwk",
    "species_tree": "species.nwk",
    "species_map_tsv": "mapping.tsv",
    "alignment": "alignment.fasta",
}
SOURCE_FIELDS = {
    "tools/validate_radte_intervals.py": "runner_sha256",
    "tools/benchmark_radte.py": "simulator_sha256",
    "tools/radte_benchmark_cases.py": "scenario_generator_sha256",
    "tools/radte_interval_simulation.py": "independent_generator_sha256",
}


def _require(condition, message):
    if not condition:
        raise ValueError(message)


def _digest(data):
    return hashlib.sha256(data).hexdigest()


def _wilson(value, total, recorded):
    interval = binomtest(value, total).proportion_ci(method="wilson")
    _require(
        math.isclose(recorded[0], interval.low, abs_tol=1e-12)
        and math.isclose(recorded[1], interval.high, abs_tol=1e-12),
        "Wilson interval mismatch",
    )


def verify(root):
    protocol = json.loads((root / "protocol.json").read_text())
    protocol_hash = _digest((root / "protocol.json").read_bytes())
    with tarfile.open(root / "source.tar.gz", "r:gz") as archive:
        source = {
            member.name: _digest(archive.extractfile(member).read())
            for member in archive.getmembers()
            if member.isfile()
        }
    total = 0
    for case in protocol["cases"]:
        name = case["name"]
        directory = root / name
        meta = json.loads((directory / "metadata.json").read_text())
        args = meta["arguments"]
        _require(meta["protocol_sha256"] == protocol_hash, "Protocol hash mismatch")
        _require(
            all(
                source.get(path) == digest
                for path, digest in meta["source_sha256"].items()
            ),
            "Archived NWKIT source mismatch",
        )
        _require(
            all(
                source.get(path) == meta[field] for path, field in SOURCE_FIELDS.items()
            ),
            "Archived study source mismatch",
        )
        _require(
            args["families"] == protocol["families_per_case"]
            and args["seed"] == case["seed"]
            and args["rate_sd"] == case["rate_sd"]
            and args["sites"] == protocol["sites"]
            and args["methods"] == protocol["methods"]
            and args["inference"] == "joint-map"
            and args["likelihood"] == "exact"
            and args["generator"] == "independent"
            and args["study_role"] == "validation",
            "Case arguments differ from the frozen protocol",
        )
        rows = [
            json.loads(line)
            for line in (directory / "cases.jsonl").read_text().splitlines()
        ]
        expected_count = args["families"] * len(args["methods"])
        _require(len(rows) == expected_count, "Missing or extra interval row")
        indexed = {(row["family"], row["method"]): row for row in rows}
        _require(len(indexed) == expected_count, "Duplicate interval row")
        with tarfile.open(directory / "inputs.tar.gz", "r:gz") as inputs:
            for family in range(args["families"]):
                paired = [indexed[family, method] for method in args["methods"]]
                _require(
                    len({row["age"] for row in paired}) == 1
                    and len({row["fitted_rate_sd"] for row in paired}) == 1,
                    "Methods used different point fits",
                )
                for row in paired:
                    _require(
                        row["status"] == "completed"
                        and row["error"] is None
                        and row["actual_estimator"] == "joint-map"
                        and row["truth"] == protocol["true_duplication_age"]
                        and row["target"] == "D",
                        "Unexpected fit or truth status",
                    )
                    available = row["interval_available"]
                    _require(
                        row["covered"]
                        == (available and row["lower"] <= row["truth"] <= row["upper"]),
                        "Coverage indicator mismatch",
                    )
                for key, filename in INPUTS.items():
                    member = inputs.extractfile(f"inputs/f{family:03d}/{filename}")
                    _require(member is not None, "Missing archived input")
                    _require(
                        _digest(member.read()) == paired[0]["input_hashes"][key]
                        and all(
                            row["input_hashes"][key] == paired[0]["input_hashes"][key]
                            for row in paired
                        ),
                        "Input hash mismatch",
                    )
        summary = json.loads((directory / "summary.json").read_text())
        for method in args["methods"]:
            subset = [indexed[family, method] for family in range(args["families"])]
            available = [row for row in subset if row["interval_available"]]
            covered = sum(row["covered"] for row in available)
            item = summary[method]
            _require(
                item["families"] == args["families"]
                and item["completed"] == args["families"]
                and item["intervals_available"] == len(available)
                and item["truth_covered"] == covered,
                "Summary count mismatch",
            )
            _require(
                math.isclose(
                    item["median_width"],
                    statistics.median(row["width"] for row in available),
                    abs_tol=1e-12,
                ),
                "Median interval width mismatch",
            )
            reasons = dict(
                Counter(
                    row["interval_status"]
                    for row in subset
                    if not row["interval_available"]
                )
            )
            _require(
                item["unavailable_reasons"] == reasons, "Unavailable reason mismatch"
            )
            _wilson(len(available), len(subset), item["availability_wilson_95"])
            _wilson(covered, len(subset), item["correct_return_wilson_95"])
            _wilson(covered, len(available), item["coverage_wilson_95"])
        total += args["families"]
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    print(f"Verified {verify(args.directory)} RADTE families in {args.directory}")


if __name__ == "__main__":
    main()
