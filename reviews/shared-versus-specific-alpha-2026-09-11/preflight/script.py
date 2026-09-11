"""Direct estimated shared/trait-specific-alpha comparison in native full OU.

Run inside the frozen GeneGalleon development runtime. Controller processes are
isolated from each fit/search, with explicit time limits and durable outcomes.
"""

import argparse
import hashlib
import json
import math
import os
import platform
import resource
import signal
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from benchmark_shift_covariance import balanced_tree

from nwkit.shift_native_fit import NativeFitOptions, fit_native_layout
from nwkit.shift_native_heuristic import NativeSearchOptions, heuristic_native_search
from nwkit.shift_native_ic import native_information_criterion
from nwkit.shift_native_model import ShiftData, ShiftLayout
from nwkit.shift_native_provenance import (
    native_implementation,
    native_implementation_sha256,
)
from nwkit.shift_simulation import simulate_shift
from nwkit.shift_simulation_cli import explicit_simulation


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, np.ndarray):
        return clean(value.tolist())
    if isinstance(value, (float, np.floating)) and not math.isfinite(value):
        return str(value)
    if isinstance(value, np.generic):
        return value.item()
    return value


def fixture(job):
    p = job["traits"]
    tree = balanced_tree(100)
    indices = [
        i
        for i in range(1, len(tree.branch_ids))
        if 10 <= tree.tip_intervals[i][1] - tree.tip_intervals[i][0] <= 30
    ]
    index = indices[len(indices) // 2]
    branch = tree.branch_ids[index]
    alpha = np.ones(p) if job["truth"] == "shared" else np.geomspace(0.25, 4.0, p)
    # Correlation applies to diffusion, not an arbitrary stationary covariance.
    covariance = 0.2 * np.eye(p) + 0.8 * np.ones((p, p))
    scenario = job["scenario"]
    optima, shifts = np.zeros((1, p)), []
    if scenario != "null":
        direction = np.ones(p) / np.sqrt(p)
        if scenario == "opposed":
            direction[:] = 0
            direction[:2] = [1 / np.sqrt(2), -1 / np.sqrt(2)]
        age = tree.remaining_times[tree.compiled.parents[index]]
        optima = np.stack([np.zeros(p), 2 * direction / -np.expm1(-alpha * age)])
        shifts = [branch]
    seed = (
        710000
        + 10000 * (job["truth"] == "different")
        + 1000 * p
        + 100 * ["null", "aligned", "opposed"].index(scenario)
        + job["replicate"]
    )
    parameters = {
        "trait_names": [f"x{j}" for j in range(p)],
        "alpha": alpha.tolist(),
        "process_tip_covariance": covariance.tolist(),
        "regime_optima": optima.tolist(),
        "shift_branch_ids": shifts,
    }
    spec = explicit_simulation(tree, parameters)
    values, _ = simulate_shift(spec, seed=seed)
    data = ShiftData.build(tree, values[0], spec.trait_names)
    return (
        data,
        branch if shifts else None,
        {
            "parameters": parameters,
            "seed": seed,
            "data_sha256": hashlib.sha256(values.tobytes()).hexdigest(),
            "tip_names": list(tree.leaf_names),
            "true_alpha": alpha.tolist(),
        },
    )


def describe_fit(data, fitted):
    joint = fitted["joint_covariance"]
    return {
        "log_likelihood": fitted["log_likelihood"],
        "shift_branch_ids": list(fitted["layout"].shifts),
        "alpha": [r["alpha"] for r in fitted["traits"]],
        "alpha_status": [r["alpha_status"] for r in fitted["traits"]],
        "process_tip_covariance": joint["process_tip_covariance"],
        "optimizer": joint["optimizer"],
        "engine": joint["engine"],
        "criteria": {
            c: native_information_criterion(data, fitted, c) for c in ["AIC", "BIC"]
        },
    }


def alarm_handler(signum, frame):
    raise TimeoutError("Prespecified per-fit/search wall-time budget exceeded.")


def worker(job, result_path):
    signal.signal(signal.SIGALRM, alarm_handler)
    outcome = {"job": job, "implementation_sha256": native_implementation_sha256()}
    start = time.perf_counter()
    stage = "fixture"
    try:
        data, true, truth = fixture(job)
        outcome["generating"] = truth
        options = NativeFitOptions(trait_covariance="full", alpha_model=job["mode"])
        if job["kind"] == "timing":
            layout = ShiftLayout.build(data.tree, [] if true is None else [true])
            seconds, results = [], []
            for run in range(job["repeats"] + 1):
                stage = "warmup" if run == 0 else f"repeat_{run}"
                signal.alarm(job["timeout"])
                t = time.perf_counter()
                fitted = fit_native_layout(data, layout, options=options)
                elapsed = time.perf_counter() - t
                signal.alarm(0)
                results.append(describe_fit(data, fitted))
                if run:
                    seconds.append(elapsed)
            outcome.update(
                seconds=seconds, median_seconds=float(np.median(seconds)), fits=results
            )
        else:
            stage = "search"
            signal.alarm(job["timeout"])
            search = heuristic_native_search(
                data,
                options=NativeSearchOptions(
                    max_shifts=1,
                    candidate_pool=4,
                    refit_budget=5,
                    beam_width=4,
                    screening_budget=10000,
                    lasso_iterations=200,
                ),
                fit_arguments={"options": options},
                criterion="AIC",
            )
            signal.alarm(0)
            retained = [
                describe_fit(data, fit) for fit in search.best_by_complexity.values()
            ]
            selected = {}
            for criterion in ["AIC", "BIC"]:
                best = min(retained, key=lambda f: f["criteria"][criterion]["score"])
                branches = best["shift_branch_ids"]
                selected[criterion] = {
                    "branches": branches,
                    "any_shift": bool(branches),
                    "exact_branch": true is not None and branches == [true],
                    "false_branch": any(b != true for b in branches),
                    "fit": best,
                }
            outcome.update(
                selected=selected,
                retained=retained,
                records=search.records,
                search_metadata=search.metadata,
            )
        outcome["status"] = "complete"
    except TimeoutError as exc:
        outcome.update(status="timeout", error=str(exc), stage=stage)
    except Exception as exc:
        outcome.update(
            status="failed", error=f"{type(exc).__name__}: {exc}", stage=stage
        )
    finally:
        signal.alarm(0)
    outcome["elapsed_seconds_not_benchmark"] = time.perf_counter() - start
    outcome["worker_peak_rss_kib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    Path(result_path).write_text(
        json.dumps(clean(outcome), indent=2, allow_nan=False) + "\n"
    )


def run_job(pair):
    job, path = pair
    log = path.with_suffix(".log")
    with log.open("w") as stream:
        completed = subprocess.run(
            [
                sys.executable,
                __file__,
                "--worker",
                json.dumps(job),
                "--output",
                str(path),
            ],
            stdout=stream,
            stderr=stream,
        )
    if completed.returncode or not path.exists():
        return {
            "job": job,
            "status": "worker_failed",
            "returncode": completed.returncode,
        }
    return json.loads(path.read_text())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--part", choices=["timing", "pilot"], default="timing")
    parser.add_argument("--traits", default="2,5,10")
    parser.add_argument("--replicates", type=int, default=10)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=int, default=120)
    args = parser.parse_args()
    if args.worker:
        worker(json.loads(args.worker), args.output)
        return
    args.output.mkdir(parents=True, exist_ok=True)
    jobs = []
    traits = [int(p) for p in args.traits.split(",")]
    for p in traits:
        for truth in ["shared", "different"]:
            scenarios = (
                ["opposed"] if args.part == "timing" else ["null", "aligned", "opposed"]
            )
            for scenario in scenarios:
                count = 1 if args.part == "timing" else args.replicates
                for rep in range(count):
                    # Alternate execution order in the serial timing experiment.
                    modes = (
                        ["shared", "trait-specific"]
                        if (p + rep + (truth == "different")) % 2 == 0
                        else ["trait-specific", "shared"]
                    )
                    for mode in modes:
                        job = {
                            "kind": args.part,
                            "traits": p,
                            "truth": truth,
                            "scenario": scenario,
                            "replicate": rep,
                            "mode": mode,
                            "repeats": args.repeats,
                            "timeout": args.timeout,
                        }
                        name = f"{args.part}-p{p}-{truth}-{scenario}-{rep}-{mode}.json"
                        jobs.append((job, args.output / name))
    manifest = {
        "arguments": vars(args) | {"output": str(args.output)},
        "implementation": native_implementation(),
        "environment": {
            "machine": platform.machine(),
            "threads": {
                k: os.environ.get(k)
                for k in ["OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"]
            },
        },
        "jobs": [job for job, path in jobs],
    }
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    rows = []
    workers = 1 if args.part == "timing" else args.workers
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [executor.submit(run_job, pair) for pair in jobs]
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            (args.output / "results.json").write_text(
                json.dumps(clean(rows), indent=2) + "\n"
            )
            print(
                json.dumps(
                    {
                        "completed": len(rows),
                        "total": len(jobs),
                        "status": row["status"],
                        "job": row["job"],
                    }
                ),
                flush=True,
            )


if __name__ == "__main__":
    main()
