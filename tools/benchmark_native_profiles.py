"""Measure repeated 100-shift candidate profiles, excluding search/calibration.

Run the same script against each checkout via PYTHONPATH. A warmup and three
repeats report all scores for equivalence checking, time, and process peak RSS.
"""

import argparse
import hashlib
import json
import platform
import resource
import sys
import time
from pathlib import Path

import numpy as np
from benchmark_native_shift import balanced_tree, pectinate_tree

from nwkit.shift_native_fit import fit_native_layout
from nwkit.shift_native_model import ShiftData, ShiftLayout
from nwkit.shift_native_quick import NativeQuickProfile
from nwkit.util import read_tree


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--shape", choices=["balanced", "pectinate"], default="balanced"
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--traits", type=int, choices=[1, 2, 4], default=2)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Use a new output file.")
    rng = np.random.default_rng(20260910)
    newick = (balanced_tree if args.shape == "balanced" else pectinate_tree)(1000)
    values = rng.normal(size=(1000, args.traits))
    values[100:200] += 3
    if args.traits > 1:
        values[20:50, 1] = np.nan
    started = time.perf_counter()
    data = ShiftData.build(
        read_tree(newick, "auto", True, quiet=True),
        values,
        [f"trait{i}" for i in range(args.traits)],
    )
    branches = tuple(
        data.tree.branch_ids[i]
        for tip, i in enumerate(data.tree.compiled.leaf_indices)
        if np.isfinite(data.values[tip]).all()
    )[::7][:128]
    layouts = [
        ShiftLayout.build(data.tree, sorted(rng.choice(branches, 100, replace=False)))
        for _ in range(20)
    ]
    null = fit_native_layout(
        data, ShiftLayout.build(data.tree), alpha_height=0.7, process_variance=1.0
    )
    quick = NativeQuickProfile(data, null, branches, 0.7)
    setup = time.perf_counter() - started
    measurements = []
    for repeat in range(4):
        started, cpu = time.perf_counter(), time.process_time()
        scores = [quick.score(layout) for layout in layouts]
        if not np.isfinite(scores).all():
            raise ValueError(
                "The benchmark requires finite, identifiable candidate scores."
            )
        measurements.append(
            {
                "wall_seconds": time.perf_counter() - started,
                "cpu_seconds": time.process_time() - cpu,
                "scores": scores,
                "warmup": repeat == 0,
            }
        )
    result = {
        "scope": "candidate ranking only; excludes full search, covariance estimation and calibration",
        "shape": args.shape,
        "tips": 1000,
        "traits": args.traits,
        "shifts": 100,
        "pool": 128,
        "layouts": 20,
        "setup_seconds": setup,
        "measurements": measurements,
        "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        * (1 if sys.platform == "darwin" else 1024),
        "python": sys.version,
        "platform": platform.platform(),
        "input_sha256": hashlib.sha256(newick.encode() + values.tobytes()).hexdigest(),
    }
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
