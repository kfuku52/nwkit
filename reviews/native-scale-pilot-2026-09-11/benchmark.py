"""Frozen-seed scale and nested-bootstrap pilot; not an adoption study."""

import argparse
import cProfile
import hashlib
import json
import os
import platform
import resource
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np
import scipy
from benchmark_native_shift import balanced_tree, pectinate_tree

from nwkit.shift_native_bootstrap import (
    calibrate_native_search,
    native_selection_support,
)
from nwkit.shift_native_fit import NativeFitOptions, fit_native_layout
from nwkit.shift_native_heuristic import NativeSearchOptions, heuristic_native_search
from nwkit.shift_native_model import ShiftData, ShiftLayout, covariance_geometry
from nwkit.util import read_tree


def fixture(args):
    newick = (balanced_tree if args.shape == "balanced" else pectinate_tree)(args.tips)
    rng = np.random.default_rng(args.seed)
    tree = read_tree(newick, "auto", True, quiet=True)
    placeholder = ShiftData.build(
        tree,
        rng.normal(size=(args.tips, args.traits)),
        [f"trait{i}" for i in range(args.traits)],
    )
    tree = placeholder.tree
    leaves = np.asarray(tree.compiled.leaf_indices)
    chosen = leaves[np.linspace(0, args.tips - 1, args.shifts, dtype=int)]
    shifts = [tree.branch_ids[i] for i in chosen]
    groups = None
    if args.shared:
        groups = [[0]] + [
            shifts[j :: min(5, args.shifts)] for j in range(min(5, args.shifts))
        ]
    layout = ShiftLayout.build(tree, shifts, groups)
    values = np.empty((args.tips, args.traits))
    errors = np.zeros_like(values)
    if args.error != "none":
        errors[:] = np.linspace(0.01, 0.09, args.tips)[:, None]
    for trait in range(args.traits):
        alpha = 1.0 + trait
        design = layout.design(tree, alpha)
        beta = np.zeros(design.shape[1])
        # Strong terminal shifts isolate scale cost; not representative power evidence.
        for j in range(1, len(beta)):
            beta[j] = (6.0 if j % 2 else -6.0) / max(
                np.max(np.abs(design[:, j])), 1e-12
            )
        slopes, innovations, root = covariance_geometry(tree, alpha, 1.0, "OUfixedRoot")
        state = np.zeros(len(tree.branch_ids))
        state[0] = rng.normal() * np.sqrt(root)
        noise = rng.normal(size=len(state))
        parents = np.asarray(tree.compiled.parents)
        for indices in tree.levels:
            state[indices] = (
                slopes[indices] * state[parents[indices]]
                + np.sqrt(innovations[indices]) * noise[indices]
            )
        extra = 0.04 if args.error == "estimated" else 0.0
        values[:, trait] = (
            design @ beta
            + state[leaves]
            + rng.normal(size=args.tips) * np.sqrt(errors[:, trait] + extra)
        )
    data = ShiftData.build(tree, values, placeholder.trait_names, errors)
    digest = hashlib.sha256(
        newick.encode() + values.tobytes() + errors.tobytes()
    ).hexdigest()
    return data, layout, digest


def fit_summary(result):
    return {
        "shifts": list(result["layout"].shifts),
        "groups": [list(g) for g in result["layout"].groups],
        "log_likelihood": result["log_likelihood"],
        "traits": result["traits"],
    }


def execute(args, output):
    data, layout, digest = fixture(args)
    output["input_sha256"] = digest
    output["truth"] = {
        "shifts": list(layout.shifts),
        "groups": [list(g) for g in layout.groups],
    }
    fit_options = NativeFitOptions(estimate_measurement_error=args.error == "estimated")
    search_options = NativeSearchOptions(
        max_shifts=args.shifts,
        convergence=args.shared,
        candidate_pool=args.pool,
        refit_budget=args.refits,
        screening_budget=args.screening,
        beam_width=2,
    )
    output["fit_options"] = asdict(fit_options)
    output["search_options"] = asdict(search_options)
    calls = []
    output["search_calls"] = calls

    def search(sample):
        begin = time.perf_counter()
        record = {"status": "running"}
        calls.append(record)
        try:
            result = heuristic_native_search(
                sample, options=search_options, fit_arguments={"options": fit_options}
            )
            record.update(
                status="ok",
                metadata=result.metadata,
                largest_fitted_shift_count=max(
                    len(r["shift_branch_ids"]) for r in result.records
                ),
                evaluations=sum(r["evaluations"] for r in result.records),
            )
            return result
        except Exception as exc:
            record.update(status="failed", error=str(exc))
            raise
        finally:
            record["seconds"] = time.perf_counter() - begin

    begin = time.perf_counter()
    if args.mode == "fit":
        output["fit"] = fit_summary(
            fit_native_layout(data, layout, options=fit_options)
        )
        output["fit_seconds"] = time.perf_counter() - begin
        return
    result = search(data)
    output["search_seconds"] = time.perf_counter() - begin
    output["candidates"] = result.records
    output["best_unpenalized"] = fit_summary(result.families()[-1][1])
    if args.mode == "search":
        return
    begin = time.perf_counter()
    selected, calibration = calibrate_native_search(
        data, result, search, replicates=args.draws, seed=args.seed + 1
    )
    output["calibration_seconds"] = time.perf_counter() - begin
    output["calibration"] = calibration
    output["selected"] = fit_summary(selected)

    def select(sample, seed):
        candidate = search(sample)
        fit, record = calibrate_native_search(
            sample, candidate, search, replicates=args.draws, seed=seed
        )
        output.setdefault("nested_calibrations", []).append(record)
        return fit

    begin = time.perf_counter()
    output["support"] = native_selection_support(
        data, selected, select, replicates=args.support, seed=args.seed + 2
    )
    output["support_seconds"] = time.perf_counter() - begin


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode", choices=["fit", "search", "workflow"], default="search"
    )
    parser.add_argument("--tips", type=int, default=1000)
    parser.add_argument("--traits", type=int, default=1)
    parser.add_argument("--shifts", type=int, default=10)
    parser.add_argument(
        "--shape", choices=["balanced", "pectinate"], default="balanced"
    )
    parser.add_argument(
        "--error", choices=["none", "known", "estimated"], default="none"
    )
    parser.add_argument("--shared", action="store_true")
    parser.add_argument("--pool", type=int, default=128)
    parser.add_argument("--refits", type=int, default=256)
    parser.add_argument("--screening", type=int, default=100000)
    parser.add_argument("--draws", type=int, default=19)
    parser.add_argument("--support", type=int, default=2)
    parser.add_argument("--seed", type=int, default=2026091107)
    parser.add_argument("--profile", action="store_true")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    path = Path(args.output)
    if path.exists() or path.with_suffix(".prof").exists():
        parser.error("Use fresh output paths.")
    if args.shifts < 1 or args.tips <= args.shifts + 2 or args.traits < 1:
        parser.error("Require positive shifts/traits and tips > shifts + 2.")
    output = {
        "configuration": vars(args),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
        "threads": {
            k: os.environ.get(k)
            for k in (
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
                "OMP_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS",
            )
        },
    }
    profiler = cProfile.Profile() if args.profile else None
    begin = time.perf_counter()
    try:
        if profiler:
            profiler.enable()
        execute(args, output)
        output["status"] = "ok"
    except Exception as exc:
        output.update(status="failed", error=f"{type(exc).__name__}: {exc}")
    finally:
        if profiler:
            profiler.disable()
            profiler.dump_stats(path.with_suffix(".prof"))
        output["total_seconds"] = time.perf_counter() - begin
        output["peak_rss_bytes"] = resource.getrusage(
            resource.RUSAGE_SELF
        ).ru_maxrss * (1 if platform.system() == "Darwin" else 1024)
        path.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    print(
        json.dumps(
            {k: output[k] for k in ("status", "total_seconds", "peak_rss_bytes")}
        ),
        flush=True,
    )
    if output["status"] != "ok":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
