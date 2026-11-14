"""Paired 100-tip benchmarks of the actual native joint-covariance implementation.

Timing runs are sequential with a warmup and repeated measurements. Calibration
jobs may run in separate processes, but their times are not speed measurements.
The calibration pilot fixes alpha at its generating value: covariance and means
are estimated and the entire search is replayed from each fitted null model.
"""

import argparse
import json
import math
import os
import platform
import resource
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import replace
from functools import partial
from pathlib import Path

import numpy as np
import scipy
from ete4 import Tree

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from nwkit.shift_joint_model import evaluate_joint  # noqa: E402
from nwkit.shift_native_bootstrap import calibrate_native_search  # noqa: E402
from nwkit.shift_native_fit import NativeFitOptions, fit_native_layout  # noqa: E402
from nwkit.shift_native_model import ShiftData, ShiftLayout, ShiftTree  # noqa: E402
from nwkit.shift_native_provenance import native_implementation_sha256  # noqa: E402
from nwkit.shift_native_search import exhaustive_native_search  # noqa: E402
from nwkit.shift_simulation import simulate_shift  # noqa: E402
from nwkit.shift_simulation_cli import explicit_simulation  # noqa: E402


def balanced_tree(n):
    step = 1 / (2 ** math.ceil(math.log2(math.ceil(math.log2(n)))))

    def node(ids, depth):
        if len(ids) == 1:
            return f"t{ids[0]}"
        mid = len(ids) // 2
        children = (ids[:mid], ids[mid:])
        return (
            "("
            + ",".join(
                node(group, depth + step)
                + ":"
                + repr(1 - depth if len(group) == 1 else step)
                for group in children
            )
            + ")"
        )

    return ShiftTree.build(Tree(node(list(range(n)), 0.0) + ";"))


def generated_data(n, p, rho, scenario, seed):
    tree = balanced_tree(n)
    possible = [
        i
        for i in range(1, len(tree.branch_ids))
        if 10 <= tree.tip_intervals[i][1] - tree.tip_intervals[i][0] <= 30
    ]
    index = possible[len(possible) // 2]
    branch = tree.branch_ids[index]
    covariance = (1 - rho) * np.eye(p) + rho * np.ones((p, p))
    optimum = np.zeros((1, p))
    shifts = []
    if scenario != "null":
        shifts = [branch]
        effect = np.ones(p) / math.sqrt(p)
        if scenario == "opposed":
            effect[:] = 0
            effect[:2] = [1 / math.sqrt(2), -1 / math.sqrt(2)]
        age = tree.remaining_times[tree.compiled.parents[index]]
        optimum = np.stack((np.zeros(p), 2 * effect / -np.expm1(-age)))
    spec = explicit_simulation(
        tree,
        {
            "trait_names": [f"x{j}" for j in range(p)],
            "alpha": 1.0,
            "process_tip_covariance": covariance.tolist(),
            "shift_branch_ids": shifts,
            "regime_optima": optimum.tolist(),
        },
    )
    values, _ = simulate_shift(spec, seed=seed)
    data = ShiftData.build(tree, values[0], spec.trait_names)
    return data, branch if shifts else None


def timed(function, repeats):
    function()
    seconds = []
    result = None
    for _ in range(repeats):
        start = time.perf_counter()
        result = function()
        seconds.append(time.perf_counter() - start)
    return result, seconds


def timing(args):
    rows = []
    for p in args.traits:
        data, true = generated_data(args.tips, p, 0.8, "opposed", 100 + p)
        layout = ShiftLayout.build(data.tree, [true])
        options = NativeFitOptions(trait_covariance="full", alpha_model="shared")
        fits = {}
        engines = ("auto", "pruning") if p <= args.pruning_max_traits else ("auto",)
        for engine in engines:
            fits[engine], seconds = timed(
                partial(
                    fit_native_layout,
                    data,
                    layout,
                    options=replace(options, joint_engine=engine),
                    alpha_height=1.0,
                ),
                args.repeats,
            )
            row = {
                "kind": "fixed_layout_covariance_ML",
                "traits": p,
                "engine": engine,
                "seconds": seconds,
                "median_seconds": float(np.median(seconds)),
                "log_likelihood": fits[engine]["log_likelihood"],
            }
            rows.append(row)
            print(json.dumps(row), flush=True)
        if "pruning" in fits:
            likelihood_error = abs(
                fits["auto"]["log_likelihood"] - fits["pruning"]["log_likelihood"]
            )
            covariance_error = float(
                np.max(
                    np.abs(
                        np.asarray(
                            fits["auto"]["joint_covariance"]["process_tip_covariance"]
                        )
                        - np.asarray(
                            fits["pruning"]["joint_covariance"][
                                "process_tip_covariance"
                            ]
                        )
                    )
                )
            )
            if likelihood_error > 1e-5 or covariance_error > 1e-3:
                raise AssertionError(
                    f"Numerical fit disagreement: {likelihood_error}, {covariance_error}"
                )
            rows.append(
                {
                    "kind": "equivalence",
                    "traits": p,
                    "absolute_loglik_error": likelihood_error,
                    "max_covariance_error": covariance_error,
                }
            )
        else:
            joint = fits["auto"]["joint_fit"]
            evaluated, seconds = timed(
                partial(
                    evaluate_joint,
                    data,
                    layout,
                    joint.alpha_height,
                    joint.covariance_coordinate,
                    joint.measurement_variance,
                ),
                args.repeats,
            )
            error = abs(evaluated.log_likelihood - joint.log_likelihood)
            if error > 1e-5:
                raise AssertionError(
                    f"Profile/pruning likelihood disagreement: {error}"
                )
            rows.append(
                {
                    "kind": "pruning_evaluation_at_profile_MLE",
                    "traits": p,
                    "seconds": seconds,
                    "median_seconds": float(np.median(seconds)),
                    "absolute_loglik_error": error,
                    "numerical_optimization_timed": False,
                }
            )
        for mode in ("diagonal", "full"):
            local = replace(options, trait_covariance=mode)
            search, seconds = timed(
                partial(
                    exhaustive_native_search,
                    data,
                    max_shifts=1,
                    fit_arguments={"options": local, "alpha_height": 1.0},
                    criterion="AIC",
                ),
                args.repeats,
            )
            row = {
                "kind": "exhaustive_zero_or_one_shift",
                "traits": p,
                "mode": mode,
                "alpha": "fixed_1",
                "seconds": seconds,
                "median_seconds": float(np.median(seconds)),
                "candidates": len(search.records),
            }
            rows.append(row)
            print(json.dumps(row), flush=True)
        # Estimate alpha too for a fixed layout. Keep this distinct from search
        # timings and from the conditional-alpha calibration experiment.
        fit, seconds = timed(
            partial(fit_native_layout, data, layout, options=options), args.repeats
        )
        rows.append(
            {
                "kind": "fixed_layout_estimated_alpha",
                "traits": p,
                "seconds": seconds,
                "median_seconds": float(np.median(seconds)),
                "alpha_mode": [r["alpha_status"] for r in fit["traits"]],
                "alpha": [r["alpha"] for r in fit["traits"]],
            }
        )
    return rows


def calibration_job(job):
    n, p, scenario, rep, draws = job
    seed = (
        2026091100
        + p * 10000
        + ("null", "aligned", "opposed").index(scenario) * 1000
        + rep
    )
    data, true = generated_data(n, p, 0.8, scenario, seed)
    rows = []
    for mode in ("diagonal", "full"):
        options = NativeFitOptions(trait_covariance=mode, alpha_model="shared")
        run = partial(
            exhaustive_native_search,
            max_shifts=1,
            fit_arguments={"options": options, "alpha_height": 1.0},
        )
        try:
            searched = run(data)
            selected, calibration = calibrate_native_search(
                data, searched, run, replicates=draws, seed=seed + 500000, level=0.05
            )
            shifts = list(selected["layout"].shifts)
            rows.append(
                {
                    "traits": p,
                    "scenario": scenario,
                    "replicate": rep,
                    "mode": mode,
                    "true_branch": true,
                    "selected_branches": shifts,
                    "detected": bool(shifts),
                    "correct_branch": true is not None and shifts == [true],
                    "false_branch": bool(shifts) and shifts != [true],
                    "calibration": calibration,
                    "status": "complete",
                }
            )
        except (ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
            rows.append(
                {
                    "traits": p,
                    "scenario": scenario,
                    "replicate": rep,
                    "mode": mode,
                    "status": "failed",
                    "message": str(exc),
                }
            )
    return rows


def calibration(args):
    jobs = [
        (args.tips, p, scenario, rep, args.calibration_replicates)
        for p in args.traits
        for scenario in ("null", "aligned", "opposed")
        for rep in range(args.replicates)
    ]
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(calibration_job, job) for job in jobs]
        for i, future in enumerate(as_completed(futures), 1):
            rows.extend(future.result())
            if i % 5 == 0 or i == len(jobs):
                print(
                    json.dumps(
                        {
                            "completed": i,
                            "jobs": len(jobs),
                            "failures": sum(r["status"] != "complete" for r in rows),
                        }
                    ),
                    flush=True,
                )
    return sorted(
        rows, key=lambda r: (r["traits"], r["scenario"], r["replicate"], r["mode"])
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--part", choices=["timing", "calibration"], required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--tips", type=int, default=100)
    parser.add_argument(
        "--traits",
        type=lambda value: list(map(int, value.split(","))),
        default=[2, 5, 10],
    )
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument(
        "--pruning-max-traits",
        type=int,
        default=2,
        help="Largest trait count for repeated general numerical ML; larger models verify pruning likelihood at the exact profile MLE.",
    )
    parser.add_argument("--replicates", type=int, default=10)
    parser.add_argument("--calibration-replicates", type=int, default=19)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if (
        args.tips < 32
        or min(args.traits) < 2
        or min(args.repeats, args.replicates, args.workers) < 1
    ):
        parser.error("Need >=32 tips, >=2 traits, and positive counts.")
    implementation = native_implementation_sha256()
    rows = timing(args) if args.part == "timing" else calibration(args)
    if native_implementation_sha256() != implementation:
        raise RuntimeError(
            "Source changed during benchmark; results are not published."
        )
    result = {
        "configuration": {**vars(args), "output": str(args.output)},
        "implementation_sha256": implementation,
        "environment": {
            "python": sys.version,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "platform": platform.platform(),
            "machine": platform.machine(),
            "threads": {
                k: os.environ.get(k)
                for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")
            },
            "parent_peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        },
        "rows": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
