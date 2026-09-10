"""Audit saved calibration evidence and check default (no-convergence) null fits."""

import argparse
import gzip
import hashlib
import json
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shift_alpha_design import cases, generate_case  # noqa: E402
from shift_calibration_audit import audit_records  # noqa: E402
from shift_simulation_cases import rate  # noqa: E402

from nwkit.shift_calibration import CalibratedSearch  # noqa: E402
from nwkit.util import read_tree  # noqa: E402

ENGINES = {}


def default_null(row):
    case = row["case"]
    key = row["tree"], case["standard_error"]
    if key not in ENGINES:
        ENGINES[key] = CalibratedSearch(
            read_tree(row["tree"], "auto", True, quiet=True),
            convergence=False,
            variances=np.full(case["tips"], case["standard_error"] ** 2),
        )
    result = ENGINES[key].fit(
        row["truth"]["observations"], seed=case["seed"] + 9000000000, replicates=199
    )
    return {
        "case_id": case["case_id"],
        "family": case["family"],
        "root_model": case["root_model"],
        "any_shift": bool(result["model"]["shift_branch_ids"]),
        "fit": result,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    directory = args.directory
    design = json.loads((directory / "protocol.json").read_text())
    with gzip.open(directory / "records.jsonl.gz", "rt") as source:
        rows = [json.loads(line) for line in source]
    expected = list(cases(design))
    if [r["case"] for r in rows] != expected:
        raise ValueError("Case identities/order do not match the frozen design")
    for path, sha in design["source_sha256"].items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest() != sha:
            raise ValueError(f"Source hash mismatch: {path}")
    for row in rows:
        newick, truth = generate_case(row["case"])
        if (
            row["tree"] != newick
            or row["truth"] != truth
            or row["status"] != "completed"
        ):
            raise ValueError(f"Truth or completion mismatch: {row['case']['case_id']}")
        json.dumps(row, allow_nan=False)
    summary = json.loads((directory / "summary.json").read_text())
    audit_records(rows, design, summary)
    nulls = [row for row in rows if row["case"]["scenario"] == "null"]
    default = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        for start in range(0, len(nulls), 16):
            default.extend(executor.map(default_null, nulls[start : start + 16]))
            print(f"Default-mode null audit: {len(default)}/{len(nulls)}", flush=True)
    grouped = defaultdict(list)
    for row in default:
        grouped[row["family"], row["root_model"]].append(row)
    summary = [
        {
            "family": family,
            "root_model": root,
            "any_shift": rate(sum(row["any_shift"] for row in members), len(members)),
        }
        for (family, root), members in grouped.items()
    ]
    with gzip.open(directory / "default-null-records.jsonl.gz", "wt") as stream:
        for row in default:
            stream.write(json.dumps(row, allow_nan=False) + "\n")
    report = {
        "verified_cases": len(rows),
        "truth_regeneration": "exact_match",
        "source_hashes": "match",
        "metrics_and_summary": "reconstructed_and_matched",
        "observed_likelihoods": "all_candidates_refitted",
        "complete_search_null_rate": rate(
            sum(r["any_shift"] for r in nulls), len(nulls)
        ),
        "default_no_convergence_null_rate": rate(
            sum(r["any_shift"] for r in default), len(default)
        ),
        "default_no_convergence_cells": summary,
        "verifier_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "record_sha256": hashlib.sha256(
            (directory / "records.jsonl.gz").read_bytes()
        ).hexdigest(),
    }
    (directory / "audit.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
