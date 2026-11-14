"""Read-only evidence audit; optional refits write to a new output directory."""

import argparse
import gzip
import hashlib
import json
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stdout
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
        row["truth"]["observations"],
        seed=case["seed"] + 9000000000,
        replicates=row["fit"]["calibration_replicates"],
        level=row["fit"]["calibration_level"],
    )
    return {
        "case_id": case["case_id"],
        "family": case["family"],
        "root_model": case["root_model"],
        "any_shift": bool(result["model"]["shift_branch_ids"]),
        "fit": result,
    }


def audit(directory, *, replay_bootstrap=False, allow_source_revision=False, workers=4):
    if allow_source_revision and not replay_bootstrap:
        raise ValueError("Source revisions require full bootstrap replay")
    design = json.loads((directory / "protocol.json").read_text())
    with gzip.open(directory / "records.jsonl.gz", "rt") as source:
        rows = [json.loads(line) for line in source]
    expected = list(cases(design))
    if [r["case"] for r in rows] != expected:
        raise ValueError("Case identities/order do not match the frozen design")
    source_changes = {}
    for path, sha in design["source_sha256"].items():
        snapshot = directory / "source-snapshot" / Path(path).name
        if hashlib.sha256(snapshot.read_bytes()).hexdigest() != sha:
            raise ValueError(f"Source snapshot hash mismatch: {path}")
        if (
            path.startswith("nwkit/")
            or Path(path).name in {"shift_alpha_design.py", "shift_simulation_cases.py"}
        ) and hashlib.sha256(Path(path).read_bytes()).hexdigest() != sha:
            if not allow_source_revision:
                raise ValueError(
                    f"Active fitting/generator source hash mismatch: {path}"
                )
            source_changes[path] = {
                "frozen": sha,
                "current": hashlib.sha256(Path(path).read_bytes()).hexdigest(),
            }
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
    with redirect_stdout(sys.stderr):
        audit_records(
            rows, design, summary, replay_bootstrap=replay_bootstrap, workers=workers
        )
    nulls = [row for row in rows if row["case"]["scenario"] == "null"]
    report = {
        "verified_cases": len(rows),
        "truth_regeneration": "exact_match",
        "source_hashes": "snapshots_match; revised_sources_replayed"
        if source_changes
        else "snapshots_and_active_fitting_sources_match",
        "source_changes": source_changes,
        "bootstrap_probabilities_replayed": len(rows) if replay_bootstrap else 0,
        "audit_module_sha256": hashlib.sha256(
            Path(__file__).with_name("shift_calibration_audit.py").read_bytes()
        ).hexdigest(),
        "metrics_and_summary": "reconstructed_and_matched",
        "observed_likelihoods": "all_candidates_refitted",
        "complete_search_null_rate": rate(
            sum(r["any_shift"] for r in nulls), len(nulls)
        ),
        "default_null_refitted": False,
        "verifier_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "record_sha256": hashlib.sha256(
            (directory / "records.jsonl.gz").read_bytes()
        ).hexdigest(),
    }
    return report, nulls


def refit_default(nulls, workers):
    default = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
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
    return default, {
        "default_null_refitted": True,
        "default_no_convergence_null_rate": rate(
            sum(r["any_shift"] for r in default), len(default)
        ),
        "default_no_convergence_cells": summary,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument(
        "--output",
        type=Path,
        help="New directory for audit outputs; input evidence is never overwritten",
    )
    parser.add_argument(
        "--refit-default-null",
        action="store_true",
        help="Refit null data without convergence; requires --output",
    )
    parser.add_argument(
        "--replay-bootstrap",
        action="store_true",
        help="Replay every saved bootstrap probability with its original seed",
    )
    parser.add_argument(
        "--allow-source-revision",
        action="store_true",
        help="Require full bootstrap replay when current fitting sources differ from the frozen snapshot",
    )
    args = parser.parse_args()
    if args.allow_source_revision and not args.replay_bootstrap:
        parser.error("--allow-source-revision requires --replay-bootstrap")
    if args.workers < 1:
        parser.error("--workers must be positive")
    if args.refit_default_null and args.output is None:
        parser.error("--refit-default-null requires a new --output directory")
    if args.output is not None:
        source, destination = args.directory.resolve(), args.output.resolve()
        if destination == source or source in destination.parents:
            parser.error("--output must be outside the input evidence directory")
        args.output.mkdir(parents=True, exist_ok=False)
    report, nulls = audit(
        args.directory,
        replay_bootstrap=args.replay_bootstrap,
        allow_source_revision=args.allow_source_revision,
        workers=args.workers,
    )
    default = []
    if args.refit_default_null:
        default, metrics = refit_default(nulls, args.workers)
        report.update(metrics)
    if args.output is not None:
        (args.output / "audit.json").write_text(json.dumps(report, indent=2) + "\n")
        if args.refit_default_null:
            with gzip.open(
                args.output / "default-null-records.jsonl.gz", "wt"
            ) as stream:
                for row in default:
                    stream.write(json.dumps(row, allow_nan=False) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
