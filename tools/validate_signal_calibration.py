#!/usr/bin/env python3
"""Independent-family calibration study for Pagel's lambda inference."""

import argparse
import hashlib
import json
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from contextlib import nullcontext
from pathlib import Path

import numpy as np
import scipy
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from nwkit.signal_stats import lambda_fit, lambda_parametric_bootstrap  # noqa: E402

SCENARIOS = (
    "balanced-8-null",
    "pectinate-8-null",
    "balanced-8-missing-null",
    "balanced-8-known-se-null",
    "balanced-32-null",
    "balanced-8-alternative",
)


def _covariance(tips, shape):
    """Build a tree-realizable covariance without using NWKIT's tree engine."""
    covariance = np.zeros((tips, tips), dtype=float)
    if shape == "balanced":
        block = 2
        while block < tips:
            for start in range(0, tips, block):
                member = np.zeros(tips)
                member[start : start + block] = 1
                covariance += 0.25 * np.outer(member, member)
            block *= 2
    elif shape == "pectinate":
        weight = 1.2 / (tips - 2)
        for start in range(1, tips - 1):
            member = np.zeros(tips)
            member[start:] = 1
            covariance += weight * np.outer(member, member)
    else:
        raise ValueError(f"Unsupported tree shape: {shape}")
    np.fill_diagonal(covariance, 2.0)
    np.linalg.cholesky(covariance)
    return covariance


def _scenario(name):
    if name not in SCENARIOS:
        raise ValueError(f"Unknown signal scenario: {name}")
    tips = 32 if name == "balanced-32-null" else 8
    shape = "pectinate" if name == "pectinate-8-null" else "balanced"
    covariance = _covariance(tips, shape)
    errors = np.zeros(tips)
    if name == "balanced-8-known-se-null":
        errors = np.array([0.1, 0.2, 0.5, 0.8, 0.15, 0.4, 0.7, 1.0])
    if name == "balanced-8-missing-null":
        observed = np.array([0, 1, 2, 4, 5, 7])
        covariance = covariance[np.ix_(observed, observed)]
        errors = errors[observed]
    return covariance, errors, (0.6 if name.endswith("alternative") else 0.0)


def _seed(master, name, replicate, stream):
    digest = int.from_bytes(hashlib.sha256(name.encode()).digest()[:8], "little")
    return np.random.SeedSequence([master, digest, replicate, stream])


def _one(task):
    name, replicate, master, inner, profile_ci = task
    covariance, errors, true_lambda = _scenario(name)
    diagonal = np.diag(np.diag(covariance))
    generating = 0.7 * (diagonal + true_lambda * (covariance - diagonal)) + np.diag(
        errors**2
    )
    outer_rng = np.random.default_rng(_seed(master, name, replicate, 0))
    values = 1.2 + np.linalg.cholesky(generating) @ outer_rng.standard_normal(
        len(errors)
    )
    result = {
        "scenario": name,
        "replicate": replicate,
        "tips": len(errors),
        "true_lambda": true_lambda,
        "input_sha256": hashlib.sha256(values.tobytes()).hexdigest(),
    }
    try:
        fit = lambda_fit(
            covariance, values, errors, ci_level=0.95 if profile_ci else None
        )
    except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
        result.update(status="fit_failed", error=str(exc))
        return result
    result["status"] = fit["status"]
    if "likelihood_ratio" not in fit:
        return result
    result.update(
        estimate=fit["estimate"],
        likelihood_ratio=fit["likelihood_ratio"],
        chi2_p_value=fit["p_value"],
    )
    if profile_ci:
        result.update(ci_lower=fit["ci_lower"], ci_upper=fit["ci_upper"])
    if inner:
        try:
            result["bootstrap_p_value"] = lambda_parametric_bootstrap(
                covariance,
                values,
                errors,
                fit["likelihood_ratio"],
                inner,
                np.random.default_rng(_seed(master, name, replicate, 1)),
            )
        except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
            result.update(bootstrap_status="failed", bootstrap_error=str(exc))
    return result


def _rate(numerator, denominator):
    if not denominator:
        return {"count": numerator, "total": denominator, "rate": None, "mc95": None}
    fraction = numerator / denominator
    z = norm.ppf(0.975)
    scale = 1 + z**2 / denominator
    center = (fraction + z**2 / (2 * denominator)) / scale
    half = (
        z
        * np.sqrt(fraction * (1 - fraction) / denominator + z**2 / (4 * denominator**2))
        / scale
    )
    return {
        "count": numerator,
        "total": denominator,
        "rate": fraction,
        "mc95": [max(0.0, center - half), min(1.0, center + half)],
    }


def _summarize(records, scenarios, level):
    summary = []
    for name in scenarios:
        rows = [row for row in records if row["scenario"] == name]
        entry = {
            "scenario": name,
            "generated": len(rows),
            "status_counts": {
                status: sum(row["status"] == status for row in rows)
                for status in sorted({row["status"] for row in rows})
            },
            "bootstrap_failure_reasons": {
                error: sum(row.get("bootstrap_error") == error for row in rows)
                for error in sorted(
                    {row["bootstrap_error"] for row in rows if "bootstrap_error" in row}
                )
            },
        }
        for key in ("chi2", "bootstrap"):
            available = [row for row in rows if f"{key}_p_value" in row]
            rejected = sum(row[f"{key}_p_value"] <= level for row in available)
            entry[f"{key}_p_available"] = _rate(len(available), len(rows))
            entry[f"{key}_rejection_all"] = _rate(rejected, len(rows))
            entry[f"{key}_rejection_available"] = _rate(rejected, len(available))
        intervals = [row for row in rows if "ci_lower" in row]
        if intervals:
            covered = sum(
                row["ci_lower"] <= row["true_lambda"] <= row["ci_upper"]
                for row in intervals
            )
            entry["profile_interval_available"] = _rate(len(intervals), len(rows))
            entry["profile_interval_coverage_all"] = _rate(covered, len(rows))
            entry["profile_interval_coverage_available"] = _rate(
                covered, len(intervals)
            )
        summary.append(entry)
    return summary


def _source_hashes():
    paths = (
        Path("nwkit/signal.py"),
        Path("nwkit/signal_stats.py"),
        Path("tools/validate_signal_calibration.py"),
    )
    return {
        str(path): hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
        for path in paths
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--scenarios", default=",".join(SCENARIOS))
    parser.add_argument("--outer", type=int, default=200)
    parser.add_argument(
        "--inner", type=int, default=199, help="Use zero for profile-CI-only studies."
    )
    parser.add_argument("--level", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=20260928)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--profile-ci", action="store_true")
    args = parser.parse_args(argv)
    scenarios = args.scenarios.split(",")
    if not scenarios or len(scenarios) != len(set(scenarios)):
        parser.error("Select distinct scenarios.")
    if any(name not in SCENARIOS for name in scenarios):
        parser.error(f"Scenarios must be chosen from {', '.join(SCENARIOS)}.")
    if (
        args.outer < 1
        or args.inner < 0
        or args.workers < 1
        or args.seed < 0
        or not np.isfinite(args.level)
        or not 0 < args.level < 1
    ):
        parser.error(
            "Positive outer/workers, nonnegative inner/seed, and an interior level are required."
        )
    if args.output.exists():
        parser.error("Output directory must be new.")
    hashes = _source_hashes()
    args.output.mkdir(parents=True)
    protocol = {
        "completed": False,
        "scenarios": scenarios,
        "outer": args.outer,
        "inner": args.inner,
        "level": args.level,
        "seed": args.seed,
        "workers": args.workers,
        "profile_ci": args.profile_ci,
        "generator": "direct Gaussian covariance Cholesky; root mean 1.2, rate 0.7",
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "source_sha256": hashes,
    }
    protocol_path = args.output / "protocol.json"
    protocol_path.write_text(json.dumps(protocol, indent=2, sort_keys=True) + "\n")
    tasks = [
        (name, replicate, args.seed, args.inner, args.profile_ci)
        for name in scenarios
        for replicate in range(args.outer)
    ]
    started = time.perf_counter()
    records = []
    manager = nullcontext() if args.workers == 1 else ProcessPoolExecutor(args.workers)
    with manager as pool, (args.output / "records.jsonl").open("w") as handle:
        iterator = map(_one, tasks) if pool is None else pool.map(_one, tasks)
        for row in iterator:
            records.append(row)
            handle.write(json.dumps(row, sort_keys=True) + "\n")
            handle.flush()
    if _source_hashes() != hashes:
        raise RuntimeError("Signal source changed during calibration.")
    summary = _summarize(records, scenarios, args.level)
    (args.output / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    protocol["elapsed_seconds"] = time.perf_counter() - started
    protocol["completed"] = True
    protocol_path.write_text(json.dumps(protocol, indent=2, sort_keys=True) + "\n")
    for row in summary:
        print(
            row["scenario"],
            "chi2",
            row["chi2_rejection_all"]["count"],
            "bootstrap",
            row["bootstrap_rejection_all"]["count"],
            "of",
            row["generated"],
        )


if __name__ == "__main__":
    main()
