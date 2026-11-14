"""Measure native search separately from fixed GLS and calibration costs."""

import argparse
import hashlib
import json
import platform
import resource
import time
from pathlib import Path

import numpy as np
from benchmark_native_shift import balanced_tree, pectinate_tree

from nwkit.shift_native_heuristic import NativeSearchOptions, heuristic_native_search
from nwkit.shift_native_model import ShiftData
from nwkit.shift_native_search import exhaustive_native_search
from nwkit.util import read_tree


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tips", type=int, default=32)
    parser.add_argument("--traits", type=int, default=2)
    parser.add_argument("--shifts", type=int, default=2)
    parser.add_argument("--candidate-pool", type=int, default=24)
    parser.add_argument("--refit-budget", type=int, default=48)
    parser.add_argument("--screening-budget", type=int, default=2000)
    parser.add_argument("--beam-width", type=int, default=2)
    parser.add_argument("--search-memory-mb", type=int, default=512)
    parser.add_argument(
        "--convergence", action=argparse.BooleanOptionalAction, default=True
    )
    parser.add_argument(
        "--shape", choices=["balanced", "pectinate"], default="balanced"
    )
    parser.add_argument("--fit-covariance", action="store_true")
    parser.add_argument("--compare-exhaustive", action="store_true")
    parser.add_argument("--seed", type=int, default=2026091001)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    path = Path(args.output)
    if path.exists():
        parser.error("Use a new output file.")
    if args.compare_exhaustive and (args.tips > 8 or args.shifts > 2):
        parser.error(
            "This benchmark's exhaustive reference is limited to 8 tips/2 shifts."
        )
    newick = (balanced_tree if args.shape == "balanced" else pectinate_tree)(args.tips)
    values = np.random.default_rng(args.seed).normal(size=(args.tips, args.traits))
    values[args.tips // 2 :] += 2
    data = ShiftData.build(
        read_tree(newick, "auto", True, quiet=True),
        values,
        [f"trait{i}" for i in range(args.traits)],
    )
    arguments = (
        {}
        if args.fit_covariance
        else {
            "alpha_height": np.geomspace(0.3, 3, args.traits),
            "process_variance": 1.0,
        }
    )
    options = NativeSearchOptions(
        max_shifts=args.shifts,
        convergence=args.convergence,
        candidate_pool=args.candidate_pool,
        refit_budget=args.refit_budget,
        screening_budget=args.screening_budget,
        beam_width=args.beam_width,
        memory_limit=args.search_memory_mb * 1024**2,
    )
    started = time.perf_counter()
    result = heuristic_native_search(data, options=options, fit_arguments=arguments)
    elapsed = time.perf_counter() - started
    best = result.families()[-1][1]
    output = {
        "scope": "single search; excludes calibration/support; not a speedup claim",
        "configuration": vars(args),
        "input_sha256": hashlib.sha256(newick.encode() + values.tobytes()).hexdigest(),
        "wall_seconds": elapsed,
        "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        * (1 if platform.system() == "Darwin" else 1024),
        "best_log_likelihood": best["log_likelihood"],
        "selected_shifts_without_penalty_or_calibration": list(best["layout"].shifts),
        "largest_fitted_shift_count": max(
            len(row["shift_branch_ids"]) for row in result.records
        ),
        "metadata": result.metadata,
        "candidates": result.records,
    }
    if args.compare_exhaustive:
        started = time.perf_counter()
        exact = exhaustive_native_search(
            data,
            max_shifts=args.shifts,
            convergence=args.convergence,
            fit_arguments=arguments,
        )
        output["exhaustive_seconds"] = time.perf_counter() - started
        output["exhaustive_log_likelihood"] = exact.families()[-1][1]["log_likelihood"]
        output["heuristic_log_likelihood_gap"] = (
            output["exhaustive_log_likelihood"] - best["log_likelihood"]
        )
    path.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    print(
        json.dumps(
            {
                key: value
                for key, value in output.items()
                if key not in {"candidates", "metadata"}
            }
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
