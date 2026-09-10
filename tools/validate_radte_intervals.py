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
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from benchmark_radte import simulate
from scipy.stats import binomtest

from nwkit import __version__
from nwkit.cli import parser
from nwkit.radte import run_dating
from nwkit.radte_calibrated import calibrated_profile_intervals
from nwkit.radte_exact_interval import exact_log_duration_intervals
from nwkit.radte_inputs import read_inputs
from nwkit.radte_model import laplace_intervals
from nwkit.radte_studentized import chronology_domain, studentized_intervals

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
        "interval_return_fraction": len(available) / len(rows) if rows else None,
        "truth_covered": covered,
        "coverage_among_available": covered / len(available) if available else None,
        "correct_interval_returned_fraction": covered / len(rows) if rows else None,
        "full_domain_return_fraction": sum(
            row.get("full_domain", False) for row in available
        )
        / len(rows)
        if rows
        else None,
        "unavailable_reasons": dict(
            Counter(
                row.get("interval_status", row["status"])
                for row in rows
                if not row["interval_available"]
            )
        ),
        "median_width": float(np.median([row["width"] for row in available]))
        if available
        else None,
    }
    if available:
        ci = binomtest(covered, len(available)).proportion_ci(method="wilson")
        result["coverage_wilson_95"] = [ci.low, ci.high]
    if rows:
        for name, count in [
            ("availability", len(available)),
            ("correct_return", covered),
        ]:
            ci = binomtest(count, len(rows)).proportion_ci(method="wilson")
            result[name + "_wilson_95"] = [ci.low, ci.high]
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
        "--substitution-model",
        options.fit_model,
        "--gamma-categories",
        "1",
        "--max-age",
        str(options.max_age),
        "--rate-correlation",
        str(options.fit_rho),
        "--inference",
        options.inference,
        "--likelihood",
        options.likelihood,
        "--starts",
        str(options.starts),
        "--seed",
        str(
            (
                options.optimizer_seed
                if options.optimizer_seed is not None
                else options.seed
            )
            + family
        ),
        "--out-prefix",
        str(options.output / "unused"),
    ]
    if not options.branch_only:
        arguments.extend(["--alignment", str(directory / "alignment.fasta")])
    if (directory / "bounds.tsv").exists():
        arguments.extend(["--species-node-bounds-tsv", str(directory / "bounds.tsv")])
    if options.fit_rate_sd is not None:
        arguments.extend(["--rate-sd", str(options.fit_rate_sd)])
    return parser.parse_args(arguments)


def evaluate_family(directory, options, family):
    truth_data = json.loads((directory / "truth.json").read_text())
    truths = truth_data.get("target_ages", {"D": truth_data["duplication_age"]})
    args = dating_args(directory, options, family)
    hashes = {
        name: file_hash(Path(getattr(args, name)))
        for name in (
            "gene_tree",
            "generax_nhx",
            "species_tree",
            "species_map_tsv",
            "alignment",
            "species_node_bounds_tsv",
        )
        if getattr(args, name, None)
    }
    chronology = read_inputs(args)
    fit, problem, _, _ = run_dating(chronology, args)
    actual = (
        "branch-only"
        if problem.likelihood is None
        else ("marginal" if getattr(problem, "marginal", False) else "joint-map")
    )
    bounds, _ = problem.bounds_and_constraint()
    x = fit.parameters
    boundary = {
        "active_parameter_lower": np.flatnonzero(x - bounds.lb < 1e-7).tolist(),
        "active_parameter_upper": np.flatnonzero(bounds.ub - x < 1e-7).tolist(),
        "gradient": problem.value_gradient(x)[1].tolist(),
    }
    if getattr(problem, "marginal", False) and problem.fixed_sd is None:
        zero = x.copy()
        zero[-1] = 0
        boundary["variance_score_at_zero_fixed_other_parameters"] = float(
            problem.value_gradient(zero)[1][-1]
        )
    rows = []
    evaluators = {"laplace": laplace_intervals, "studentized": studentized_intervals}
    evaluators["exact-log-duration"] = exact_log_duration_intervals
    for method in options.methods:
        result = copy.deepcopy(fit)
        error = None
        try:
            if method == "calibrated-profile":
                calibrated_profile_intervals(
                    result,
                    problem,
                    level=options.level,
                    replicates=options.replicates,
                    grid_points=options.grid_points,
                    starts=options.starts,
                    seed=args.seed,
                )
            else:
                evaluated = evaluators[method](result, problem, options.level)
                if evaluated is False:
                    result.interval_status = (
                        "unavailable-unsupported-exact-log-duration-model"
                    )
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
            error = str(exc)
            result.interval_lower = result.interval_upper = None
            result.interval_status = "unavailable-interval-exception"
        for target in options.targets:
            matches = [
                i for i, node in enumerate(chronology.nodes) if node.name == target
            ]
            if len(matches) != 1 or target not in truths:
                rows.append(
                    dict(
                        family=family,
                        target=target,
                        method=method,
                        status="target-unmatched",
                        interval_available=False,
                        covered=False,
                    )
                )
                continue
            group = chronology.group_by_node[matches[0]]
            truth = truths[target]
            row = interval_row(result, fit, problem, group, truth)
            row.update(
                family=family,
                seed=options.seed + family,
                optimizer_seed=args.seed,
                target=target,
                method=method,
                actual_estimator=actual,
                error=error,
                boundary_diagnostics=boundary,
                input_hashes=hashes,
                calibration_profile=result.calibration_profile,
            )
            rows.append(row)
    return rows


def interval_row(result, fit, problem, group, truth):
    chronology = problem.chronology
    domain_lower, domain_upper = chronology_domain(chronology)
    available = result.interval_lower is not None and result.interval_upper is not None
    lower = (
        float(result.interval_lower[group] * chronology.scale) if available else None
    )
    upper = (
        float(result.interval_upper[group] * chronology.scale) if available else None
    )
    assert np.array_equal(result.ages, fit.ages)
    assert np.array_equal(result.parameters, fit.parameters)
    assert problem.feasible(result.parameters)
    if available:
        assert np.isfinite([lower, upper]).all() and lower <= upper
        assert lower >= chronology.lower[group] * chronology.scale
        assert upper <= chronology.upper[group] * chronology.scale
    return {
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
        "fixed_target": group not in problem.free,
        "truth_within_domain": bool(
            domain_lower[group] * chronology.scale
            <= truth
            <= domain_upper[group] * chronology.scale
        ),
        "full_domain": bool(
            available
            and np.isclose(lower, domain_lower[group] * chronology.scale)
            and np.isclose(upper, domain_upper[group] * chronology.scale)
        ),
    }


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
    ap.add_argument("--optimizer-seed", type=int)
    ap.add_argument(
        "--study-role", choices=["development", "validation"], default="development"
    )
    ap.add_argument(
        "--protocol",
        type=Path,
        help="Frozen study protocol; required for validation runs.",
    )
    ap.add_argument("--generator", choices=["legacy", "independent"], default="legacy")
    ap.add_argument("--shape", choices=["balanced", "pectinate"], default="balanced")
    ap.add_argument("--generate-rho", type=float, default=0.0)
    ap.add_argument("--fit-rho", type=float, default=0.0)
    ap.add_argument("--generate-model", choices=["jc69", "hky"], default="jc69")
    ap.add_argument("--fit-model", choices=["jc69", "hky", "gtr"], default="jc69")
    ap.add_argument("--calibration-width", type=float, default=0.0)
    ap.add_argument("--copy-ratio", type=float, default=1.0)
    ap.add_argument("--branch-only", action="store_true")
    ap.add_argument("--max-age", type=float, default=100.0)
    ap.add_argument(
        "--inference", choices=["auto", "marginal", "joint-map"], default="auto"
    )
    ap.add_argument(
        "--likelihood", choices=["auto", "quadratic", "exact"], default="auto"
    )
    ap.add_argument("--targets", nargs="+", default=["D"])
    ap.add_argument(
        "--methods",
        nargs="+",
        choices=["laplace", "studentized", "exact-log-duration", "calibrated-profile"],
        default=["laplace", "studentized"],
    )
    ap.add_argument("--replicates", type=int, default=199)
    ap.add_argument("--grid-points", type=int, default=17)
    ap.add_argument("--level", type=float, default=0.95)
    a = ap.parse_args()
    if a.study_role == "validation" and (
        a.protocol is None or not a.protocol.is_file()
    ):
        ap.error("Validation requires an existing --protocol frozen before the run")
    if (
        not all(
            np.isfinite(v)
            for v in (
                a.generate_rho,
                a.fit_rho,
                a.calibration_width,
                a.copy_ratio,
                a.max_age,
            )
        )
        or not 0 <= a.generate_rho < 1
        or not 0 <= a.fit_rho < 1
        or not 0 <= a.calibration_width < 1
        or a.copy_ratio <= 0
        or a.max_age <= 0
    ):
        ap.error("Invalid correlation, calibration width, copy ratio or maximum age")
    if a.generator == "legacy" and (
        a.shape != "balanced" or a.generate_rho != 0 or a.generate_model != "jc69"
    ):
        ap.error(
            "Nondefault shape, generating correlation/model require --generator independent"
        )
    if len(set(a.targets)) != len(a.targets) or len(set(a.methods)) != len(a.methods):
        ap.error("Targets and methods must be unique")
    if a.species < 2 or (a.scenario == "nested" and a.species < 4):
        ap.error("At least two species required (four for nested duplication)")
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
        "protocol_sha256": file_hash(a.protocol) if a.protocol is not None else None,
        "independent_generator_sha256": file_hash(
            ROOT / "tools/radte_interval_simulation.py"
        ),
    }
    (a.output / "metadata.json").write_text(json.dumps(metadata, indent=2))
    rows = []
    started = time.perf_counter()
    with (a.output / "cases.jsonl").open("w") as handle:
        for family in range(a.families):
            directory = (a.input_root or a.output / "inputs") / f"f{family:03d}"
            try:
                if a.input_root is None:
                    generator = simulate
                    extra = {}
                    if a.generator == "independent":
                        from radte_interval_simulation import (
                            simulate as independent_simulate,
                        )

                        generator = independent_simulate
                        extra = dict(
                            shape=a.shape, rho=a.generate_rho, model=a.generate_model
                        )
                    generator(
                        directory,
                        a.species,
                        a.sites,
                        a.seed + family,
                        a.rate_sd,
                        scenario=a.scenario,
                        width=a.calibration_width,
                        copy_ratio=a.copy_ratio,
                        **extra,
                    )
                result = evaluate_family(directory, a, family)
            except (ValueError, RuntimeError, np.linalg.LinAlgError) as exc:
                result = [
                    {
                        "family": family,
                        "method": method,
                        "target": target,
                        "seed": a.seed + family,
                        "status": "failed",
                        "interval_available": False,
                        "covered": False,
                        "error": str(exc),
                    }
                    for method in a.methods
                    for target in a.targets
                ]
            for row in result:
                rows.append(row)
                handle.write(json.dumps(row) + "\n")
            handle.flush()
            print(f"{family + 1}/{a.families}", flush=True)
    pd.DataFrame(rows).drop(
        columns=[
            "diagnostics",
            "input_hashes",
            "boundary_diagnostics",
            "calibration_profile",
        ],
        errors="ignore",
    ).to_csv(a.output / "cases.csv", index=False)
    summary = {
        method + (":" + target if len(a.targets) > 1 else ""): interval_summary(
            [
                row
                for row in rows
                if row["method"] == method
                and row["target"] == target
                and not row.get("fixed_target", False)
            ]
        )
        for method in a.methods
        for target in a.targets
    }
    summary["elapsed_seconds"] = time.perf_counter() - started
    (a.output / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
