#!/usr/bin/env python3
"""Verify completed calibration evidence, source archive and regenerated summaries."""

import argparse
import hashlib
import json
import tarfile
from pathlib import Path

from regression_calibration_design import Case, seed_for
from regression_calibration_engine import is_applicable
from validate_regression_calibration import encode, read_records, summarize_records


def verify(directory):
    protocol = json.loads((directory / "protocol.json").read_text())
    if protocol["status"] != "complete":
        raise ValueError(f"Experiment is not complete: {directory}")
    cases = {row["name"]: Case(**row) for row in protocol["cases"]}
    expected_replicates = set(
        range(
            protocol["replicate_start"],
            protocol["replicate_start"] + protocol["outer_replicates"],
        )
    )
    counts = {name: set() for name in cases}
    for row in read_records(directory):
        name, replicate = row["case"]["name"], row["replicate"]
        if (
            name not in cases
            or replicate not in expected_replicates
            or replicate in counts[name]
        ):
            raise ValueError("Unexpected or duplicate dataset")
        counts[name].add(replicate)
        if row["data_seed"] != seed_for(protocol["seed"], name, replicate, "data"):
            raise ValueError("Data seed mismatch")
        if row["fit_seed"] != seed_for(protocol["seed"], name, replicate, "fit"):
            raise ValueError("Fit seed mismatch")
        if (
            hashlib.sha256(encode(row["data"]).encode()).hexdigest()
            != row["input_sha256"]
        ):
            raise ValueError("Generated input hash mismatch")
        expected_methods = {
            method
            for method in protocol["methods"]
            if is_applicable(cases[name], method)
        }
        if (
            len(row["results"]) != len(expected_methods)
            or {item["method"] for item in row["results"]} != expected_methods
        ):
            raise ValueError("Missing or duplicate method result")
    if any(replicates != expected_replicates for replicates in counts.values()):
        raise ValueError("Missing simulation datasets")
    source_hashes = {}
    with tarfile.open(directory / "source.tar.gz") as archive:
        for member in archive.getmembers():
            if not member.isfile() or member.name in source_hashes:
                raise ValueError("Unexpected source archive entry")
            source_hashes[member.name] = hashlib.sha256(
                archive.extractfile(member).read()
            ).hexdigest()
    if source_hashes != protocol["source_sha256"]:
        raise ValueError("Source archive/hash mismatch")
    regenerated = summarize_records(read_records(directory))
    if encode(regenerated) != encode(
        json.loads((directory / "summary.json").read_text())
    ):
        raise ValueError("Summary does not reproduce from records")
    return {
        "directory": str(directory),
        "datasets": sum(map(len, counts.values())),
        "summary_rows": len(regenerated),
        "source_files": len(source_hashes),
        "files_sha256": {
            name: hashlib.sha256((directory / name).read_bytes()).hexdigest()
            for name in (
                "protocol.json",
                "records.jsonl.gz",
                "summary.json",
                "source.tar.gz",
            )
        },
        "status": "verified",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directories", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    results = []
    for directory in args.directories:
        results.append(verify(directory))
        print(
            f"Verified {directory.name}: {results[-1]['datasets']} datasets", flush=True
        )
    with args.output.open("x") as handle:
        handle.write(
            encode(
                {
                    "verifier_source_sha256": hashlib.sha256(
                        Path(__file__).read_bytes()
                    ).hexdigest(),
                    "experiments": results,
                }
            )
            + "\n"
        )


if __name__ == "__main__":
    main()
