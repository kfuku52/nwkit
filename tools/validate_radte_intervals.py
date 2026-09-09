"""Paired interval coverage on reproducible independent gene families.

Run from a checkout with PYTHONPATH=. and numerical-library threads fixed to 1.
This evaluates statistical coverage, not CLI runtime. Both interval methods use
the same fitted model. --input-root can reuse f000/, f001/, ... input directories.
"""

import argparse
import copy
import hashlib
import json
import os
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from benchmark_radte import simulate
from scipy.stats import binomtest

from nwkit import __version__
from nwkit.cli import parser
from nwkit.radte import run_dating
from nwkit.radte_inputs import read_inputs
from nwkit.radte_model import laplace_intervals
from nwkit.radte_studentized import studentized_intervals

ROOT = Path(__file__).resolve().parents[1]


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def interval_summary(rows):
    available = [row for row in rows if row["interval_available"]]
    covered = sum(row["covered"] for row in available)
    result = {
        "families": len(rows),
        "completed": sum(row["status"] == "completed" for row in rows),
        "intervals_available": len(available),
        "truth_covered": covered,
        "coverage_among_available": covered / len(available) if available else None,
        "correct_interval_returned_fraction": covered / len(rows),
        "median_width": float(np.median([row["width"] for row in available]))
        if available
        else None,
    }
    if available:
        ci = binomtest(covered, len(available)).proportion_ci(method="wilson")
        result["coverage_wilson_95"] = [ci.low, ci.high]
    return result


def dating_args(directory, options, family):
    inputs = (
        ["--generax-nhx", str(directory / "ml.nhx")]
        if (directory / "ml.nhx").exists()
        else ["--gene-tree", str(directory / "gene.nwk"), "--reconcile", "lca"]
    )
    arguments = [
        "radte",
        *inputs,
        "--species-tree",
        str(directory / "species.nwk"),
        "--species-map-tsv",
        str(directory / "mapping.tsv"),
        "--alignment",
        str(directory / "alignment.fasta"),
        "--substitution-model",
        "jc69",
        "--gamma-categories",
        "1",
        "--max-age",
        "100",
        "--starts",
        str(options.starts),
        "--seed",
        str(options.seed + family),
        "--out-prefix",
        str(options.output / "unused"),
    ]
    if options.fit_rate_sd is not None:
        arguments.extend(["--rate-sd", str(options.fit_rate_sd)])
    return parser.parse_args(arguments)


def evaluate_family(directory, options, family):
    truth = json.loads((directory / "truth.json").read_text())["duplication_age"]
    args = dating_args(directory, options, family)
    hashes = {
        name: file_hash(Path(getattr(args, name)))
        for name in (
            "gene_tree",
            "generax_nhx",
            "species_tree",
            "species_map_tsv",
            "alignment",
        )
        if getattr(args, name, None)
    }
    chronology = read_inputs(args)
    fit, problem, _, _ = run_dating(chronology, args)
    target = next(i for i, node in enumerate(chronology.nodes) if node.name == "D")
    group = chronology.group_by_node[target]
    actual = "marginal" if getattr(problem, "marginal", False) else "joint-map"
    rows = []
    for method, evaluate in [
        ("laplace", laplace_intervals),
        ("studentized", studentized_intervals),
    ]:
        result = copy.deepcopy(fit)
        evaluate(result, problem, options.level)
        available = (
            result.interval_lower is not None and result.interval_upper is not None
        )
        lower = (
            float(result.interval_lower[group] * chronology.scale)
            if available
            else None
        )
        upper = (
            float(result.interval_upper[group] * chronology.scale)
            if available
            else None
        )
        assert np.array_equal(result.ages, fit.ages)
        assert np.array_equal(result.parameters, fit.parameters)
        assert problem.feasible(result.parameters)
        if available:
            assert np.isfinite([lower, upper]).all() and lower <= upper
            assert lower >= chronology.lower[group] * chronology.scale
            assert upper <= chronology.upper[group] * chronology.scale
        rows.append(
            {
                "family": family,
                "seed": options.seed + family,
                "method": method,
                "actual_estimator": actual,
                "status": "completed",
                "age": float(fit.ages[group] * chronology.scale),
                "truth": truth,
                "fitted_rate_sd": fit.log_rate_sd,
                "rate_variance_estimated": problem.rate_variance_estimated,
                "interval_status": result.interval_status,
                "interval_available": available,
                "lower": lower,
                "upper": upper,
                "width": upper - lower if available else None,
                "covered": available and lower <= truth <= upper,
                "diagnostics": result.diagnostics,
                "input_hashes": hashes,
            }
        )
    return rows


def options():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--input-root", type=Path)
    ap.add_argument("--families", type=int, default=200)
    ap.add_argument("--species", type=int, default=2)
    ap.add_argument("--sites", type=int, default=2000)
    ap.add_argument("--rate-sd", type=float, default=0.3, help="Simulated log-rate SD.")
    ap.add_argument(
        "--fit-rate-sd",
        type=float,
        help="Supply an SD to inference (oracle/sensitivity check).",
    )
    ap.add_argument("--scenario", choices=["root", "nested", "loss"], default="root")
    ap.add_argument("--seed", type=int, default=91300000)
    ap.add_argument("--starts", type=int, default=3)
    ap.add_argument("--level", type=float, default=0.95)
    a = ap.parse_args()
    if min(a.families, a.species, a.sites, a.starts) < 1:
        ap.error("families, species, sites and starts must be positive")
    if not np.isfinite(a.rate_sd) or a.rate_sd < 0 or not 0 < a.level < 1:
        ap.error(
            "rate SD must be finite/nonnegative and interval level between 0 and 1"
        )
    if a.fit_rate_sd is not None and (
        not np.isfinite(a.fit_rate_sd) or a.fit_rate_sd < 0
    ):
        ap.error("fitted rate SD must be finite/nonnegative")
    return a


def main():
    a = options()
    a.output.mkdir(parents=True, exist_ok=False)
    metadata = {
        "arguments": {
            key: str(value) if isinstance(value, Path) else value
            for key, value in vars(a).items()
        },
        "command": sys.argv,
        "version": __version__,
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "threads": {
            key: value
            for key, value in os.environ.items()
            if key.endswith("NUM_THREADS") or key == "VECLIB_MAXIMUM_THREADS"
        },
        "source_sha256": {
            str(path.relative_to(ROOT)): file_hash(path)
            for path in sorted((ROOT / "nwkit").glob("*.py"))
        },
        "runner_sha256": file_hash(Path(__file__)),
        "simulator_sha256": file_hash(ROOT / "tools/benchmark_radte.py"),
        "scenario_generator_sha256": file_hash(ROOT / "tools/radte_benchmark_cases.py"),
    }
    (a.output / "metadata.json").write_text(json.dumps(metadata, indent=2))
    rows = []
    started = time.perf_counter()
    with (a.output / "cases.jsonl").open("w") as handle:
        for family in range(a.families):
            directory = (a.input_root or a.output / "inputs") / f"f{family:03d}"
            if a.input_root is None:
                simulate(
                    directory,
                    a.species,
                    a.sites,
                    a.seed + family,
                    a.rate_sd,
                    scenario=a.scenario,
                )
            try:
                result = evaluate_family(directory, a, family)
            except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
                result = [
                    {
                        "family": family,
                        "method": method,
                        "status": "failed",
                        "interval_available": False,
                        "covered": False,
                        "error": str(exc),
                    }
                    for method in ["laplace", "studentized"]
                ]
            for row in result:
                rows.append(row)
                handle.write(json.dumps(row) + "\n")
            handle.flush()
            print(f"{family + 1}/{a.families}", flush=True)
    pd.DataFrame(rows).drop(
        columns=["diagnostics", "input_hashes"], errors="ignore"
    ).to_csv(a.output / "cases.csv", index=False)
    summary = {
        method: interval_summary([row for row in rows if row["method"] == method])
        for method in ["laplace", "studentized"]
    }
    summary["elapsed_seconds"] = time.perf_counter() - started
    (a.output / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
