#!/usr/bin/env python3
"""Reproducible fixed-model RSC/GLMM calibration with explicit failure denominators."""

import argparse
import fnmatch
import gzip
import hashlib
import importlib.metadata
import json
import os
import platform
import resource
import sys
import tarfile
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.stats import norm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from regression_calibration_design import (  # noqa: E402
    Case,
    case_dict,
    cases,
    generate,
    seed_for,
)
from regression_calibration_engine import evaluate, is_applicable  # noqa: E402

import nwkit  # noqa: E402

METHODS = {
    "wald",
    "parametric-bootstrap",
    "oracle",
    "null-bootstrap",
    "profile-likelihood",
    "likelihood-ratio",
}


def clean(value):
    if isinstance(value, dict):
        return {str(key): clean(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [clean(item) for item in value]
    if isinstance(value, np.generic):
        return clean(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def encode(value):
    return json.dumps(clean(value), sort_keys=True, allow_nan=False)


def finite(value):
    try:
        return np.isfinite(float(value))
    except (ValueError, TypeError):
        return False


def rate(successes, total):
    """Wilson 95% Monte Carlo interval (not a biological uncertainty interval)."""
    if total == 0:
        return {
            "numerator": successes,
            "denominator": total,
            "estimate": None,
            "mcse": None,
            "mc95_lower": None,
            "mc95_upper": None,
        }
    p = successes / total
    z = norm.ppf(0.975)
    denominator = 1 + z * z / total
    center = (p + z * z / (2 * total)) / denominator
    half = z * np.sqrt(p * (1 - p) / total + z * z / (4 * total * total)) / denominator
    return {
        "numerator": successes,
        "denominator": total,
        "estimate": p,
        "mcse": np.sqrt(p * (1 - p) / total),
        "mc95_lower": max(0.0, center - half),
        "mc95_upper": min(1.0, center + half),
    }


def within(metric, lower, upper):
    if metric["estimate"] is None:
        return "unavailable"
    if metric["mc95_lower"] >= lower and metric["mc95_upper"] <= upper:
        return "within_prespecified_band"
    if metric["mc95_upper"] < lower or metric["mc95_lower"] > upper:
        return "outside_prespecified_band"
    return "insufficient_precision_or_overlapping_boundary"


def summarize_records(records):
    grouped = {}
    seen = set()
    for task in records:
        identity = (task["case"]["name"], task["replicate"])
        if identity in seen:
            raise ValueError(f"Duplicate simulation record {identity}")
        seen.add(identity)
        for row in task["results"]:
            if row["status"] == "not_applicable":
                continue
            key = (identity[0], row["method"])
            grouped.setdefault(key, []).append(row)
    summary = []
    for (case, method), rows in sorted(grouped.items()):
        available = [
            row
            for row in rows
            if row["status"] == "completed"
            and row.get("inference_status") == "ok"
            and finite(row.get("p_value"))
            and 0 <= float(row["p_value"]) <= 1
        ]
        ci_rows = [
            row
            for row in rows
            if row["status"] == "completed"
            and row.get("inference_status") == "ok"
            and finite(row.get("confidence_interval_lower"))
            and finite(row.get("confidence_interval_upper"))
            and row["confidence_interval_lower"] <= row["confidence_interval_upper"]
        ]
        rejected = sum(float(row["p_value"]) < 0.05 for row in available)
        covered = sum(
            row["confidence_interval_lower"]
            <= row["target_beta"]
            <= row["confidence_interval_upper"]
            for row in ci_rows
        )
        estimates = [
            float(row["coefficient"]) - row["target_beta"]
            for row in rows
            if row["status"] == "completed" and finite(row.get("coefficient"))
        ]
        durations = [row.get("seconds", 0) for row in rows]
        eligible = sum(row["status"] != "ineligible" for row in rows)
        failed = sum(row["status"] == "fit_failed" for row in rows)
        attempts = sum(row.get("bootstrap_attempts", 0) for row in rows)
        successes = sum(row.get("bootstrap_successes", 0) for row in rows)
        rejection = rate(rejected, len(available))
        coverage = rate(covered, len(ci_rows))
        entry = {
            "case": case,
            "method": method,
            "generated": len(rows),
            "target_beta": rows[0]["target_beta"],
            "status_counts": dict(Counter(row["status"] for row in rows)),
            "error_counts": dict(
                Counter(row["error"] for row in rows if "error" in row)
            ),
            "inference_status_counts": dict(
                Counter(row.get("inference_status", "not-fitted") for row in rows)
            ),
            "rejection_all_generated": rate(rejected, len(rows)),
            "rejection_given_available": rejection,
            "p_available": rate(len(available), len(rows)),
            "coverage_given_interval": coverage,
            "interval_available": rate(len(ci_rows), len(rows))
            if method != "null-bootstrap"
            else None,
            "interval_delivered_and_covers": rate(covered, len(rows))
            if method != "null-bootstrap"
            else None,
            "fit_failure_given_eligible": rate(failed, eligible),
            "bootstrap_refit_failure": rate(attempts - successes, attempts),
            "boundary_count": sum(
                row.get("boundary_warning") in (True, "yes") for row in rows
            ),
            "separation_count": sum(
                row.get("separation_warning") in (True, "yes") for row in rows
            ),
            "bias_given_estimate": float(np.mean(estimates)) if estimates else None,
            "rmse_given_estimate": float(np.sqrt(np.mean(np.square(estimates))))
            if estimates
            else None,
            "estimate_count": len(estimates),
            "mean_interval_width": float(
                np.mean(
                    [
                        row["confidence_interval_upper"]
                        - row["confidence_interval_lower"]
                        for row in ci_rows
                    ]
                )
            )
            if ci_rows
            else None,
            "fit_seconds_total": sum(durations),
            "fit_seconds_median": float(np.median(durations)),
            "fit_seconds_p95": float(np.quantile(durations, 0.95)),
            "null_calibration_band": within(rejection, 0.04, 0.06)
            if rows[0]["target_beta"] == 0
            else "alternative",
            "coverage_band": within(coverage, 0.93, 0.97),
        }
        summary.append(entry)
    return summary


def execute(task):
    case, replicate, seed, methods, bootstrap, timeout = task
    data_seed = seed_for(seed, case.name, replicate, "data")
    fit_seed = seed_for(seed, case.name, replicate, "fit")
    data = generate(case, data_seed)
    encoded = encode(data)
    result = {
        "case": case_dict(case),
        "replicate": replicate,
        "data_seed": data_seed,
        "fit_seed": fit_seed,
        "input_sha256": hashlib.sha256(encoded.encode()).hexdigest(),
        "data": data,
        "results": [
            evaluate(case, data, method, bootstrap, fit_seed, timeout)
            for method in methods
            if is_applicable(case, method)
        ],
        "worker_peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
    }
    return clean(result)


def source_files():
    root = Path(nwkit.__file__).resolve().parent.parent
    files = list((root / "nwkit").rglob("*.py"))
    files.extend(Path(__file__).resolve().parent.glob("*regression_calibration*.py"))
    return root, sorted(set(files))


def snapshot(output):
    root, paths = source_files()
    hashes = {}
    with tarfile.open(output / "source.tar.gz", "w:gz") as archive:
        for path in paths:
            relative = str(path.relative_to(root))
            hashes[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
            archive.add(path, arcname=relative, recursive=False)
    return hashes


def verify_source(hashes):
    root, paths = source_files()
    observed = {
        str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in paths
    }
    if hashes != observed:
        raise RuntimeError(
            "Source changed during calibration; results must not be published as a single implementation"
        )


def read_records(directory):
    with gzip.open(directory / "records.jsonl.gz", "rt") as handle:
        for line in handle:
            yield json.loads(line)


def write_summary(output):
    summary = summarize_records(read_records(output))
    (output / "summary.json").write_text(encode(summary) + "\n")
    lines = [
        "case\tmethod\tgenerated\tp_available\trejection_given_available\tcoverage_given_interval\tfit_failures\tseconds"
    ]
    for row in summary:
        lines.append(
            "\t".join(
                str(value)
                for value in (
                    row["case"],
                    row["method"],
                    row["generated"],
                    row["p_available"]["estimate"],
                    row["rejection_given_available"]["estimate"],
                    row["coverage_given_interval"]["estimate"],
                    row["status_counts"].get("fit_failed", 0),
                    row["fit_seconds_total"],
                )
            )
        )
    (output / "summary.tsv").write_text("\n".join(lines) + "\n")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--summarize", type=Path)
    parser.add_argument("--list-cases", action="store_true")
    parser.add_argument("--cases", default="*", help="Comma-separated shell patterns")
    parser.add_argument(
        "--case-file",
        type=Path,
        help="JSON array of Case fields instead of built-in cases",
    )
    parser.add_argument("--methods", default="wald,oracle")
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--replicate-start", type=int, default=0)
    parser.add_argument("--bootstrap-replicates", type=int, default=199)
    parser.add_argument("--seed", type=int, default=20260910)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--fit-timeout-seconds", type=float, default=120)
    parser.add_argument("--runtime-label", default="unspecified")
    args = parser.parse_args()
    try:
        available_cases = (
            cases()
            if args.case_file is None
            else [Case(**row) for row in json.loads(args.case_file.read_text())]
        )
        if len({case.name for case in available_cases}) != len(available_cases):
            raise ValueError("Case names must be unique")
    except (TypeError, ValueError, OSError) as exc:
        parser.error(str(exc))
    if args.list_cases:
        for case in available_cases:
            print(encode(case_dict(case)))
        return
    if args.summarize:
        write_summary(args.summarize)
        return
    if args.output is None:
        parser.error("--output is required")
    if (
        args.replicates < 1
        or args.workers < 1
        or args.replicate_start < 0
        or args.seed < 0
    ):
        parser.error(
            "Replicates/workers must be positive; start/seed must be nonnegative"
        )
    if args.bootstrap_replicates < 2 or args.fit_timeout_seconds <= 0:
        parser.error(
            "Bootstrap needs at least 2 replicates and timeout must be positive"
        )
    methods = list(dict.fromkeys(args.methods.split(",")))
    if set(methods) - METHODS:
        parser.error("Unknown method")
    selected = [
        case
        for case in available_cases
        if any(
            fnmatch.fnmatchcase(case.name, pattern) for pattern in args.cases.split(",")
        )
    ]
    if not selected:
        parser.error("No cases matched")
    if any(
        not any(is_applicable(case, method) for method in methods) for case in selected
    ):
        parser.error("Every selected case needs at least one applicable method")
    args.output.mkdir(parents=True, exist_ok=False)
    hashes = snapshot(args.output)
    protocol = {
        "status": "running",
        "schema": 1,
        "cases": [case_dict(case) for case in selected],
        "methods": methods,
        "outer_replicates": args.replicates,
        "replicate_start": args.replicate_start,
        "bootstrap_replicates": args.bootstrap_replicates,
        "seed": args.seed,
        "workers": args.workers,
        "timeout_seconds": args.fit_timeout_seconds,
        "runtime_label": args.runtime_label,
        "python": sys.version,
        "platform": platform.platform(),
        "nwkit_file": nwkit.__file__,
        "source_sha256": hashes,
        "package_versions": {
            name: importlib.metadata.version(name)
            for name in ("numpy", "scipy", "pandas", "ete4")
        },
        "thread_environment": {
            name: os.environ.get(name)
            for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
        },
        "alpha": 0.05,
        "confidence_level": 0.95,
        "bands": {"type_i_error": [0.04, 0.06], "coverage": [0.93, 0.97]},
        "scope": "fixed-model engine and raw-tip RSC calibration; selection inference excluded",
        "rsc_objective": "event-balanced composite, not an ordinary Gaussian likelihood",
        "bootstrap_policy": "production retries observed; reference null test uses fixed attempts and failure bounds",
        "missingness": "MAR/MNAR rates are intercept settings, not exact realized fractions",
        "binary_baseline": "probability at zero linear predictor before random effects, not marginal prevalence",
        "negative_binomial_parameterization": "case.dispersion is size r; NWKIT response_dispersion is alpha=1/r",
        "rss_units": "KiB on Linux; cumulative worker high-water mark, not per-fit peak",
    }
    protocol_path = args.output / "protocol.json"
    protocol_path.write_text(encode(protocol) + "\n")
    tasks = [
        (
            case,
            replicate,
            args.seed,
            methods,
            args.bootstrap_replicates,
            args.fit_timeout_seconds,
        )
        for replicate in range(
            args.replicate_start, args.replicate_start + args.replicates
        )
        for case in selected
    ]
    started = time.monotonic()
    last_report = started
    with (
        gzip.open(args.output / "records.jsonl.gz", "wt") as handle,
        ProcessPoolExecutor(max_workers=args.workers) as pool,
    ):
        for completed, row in enumerate(pool.map(execute, tasks, chunksize=1), start=1):
            handle.write(encode(row) + "\n")
            handle.flush()
            if time.monotonic() - last_report >= 10 or completed == len(tasks):
                print(
                    f"{completed}/{len(tasks)} datasets completed; {time.monotonic() - started:.1f}s",
                    flush=True,
                )
                last_report = time.monotonic()
    verify_source(hashes)
    write_summary(args.output)
    protocol.update(
        status="complete",
        datasets_completed=len(tasks),
        elapsed_seconds=time.monotonic() - started,
    )
    protocol_path.write_text(encode(protocol) + "\n")


if __name__ == "__main__":
    main()
