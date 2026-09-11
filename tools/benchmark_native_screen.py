"""Measure screening matrix setup in a fresh process, with reproducible outputs."""

import argparse
import hashlib
import importlib.util
import json
import os
import platform
import resource
import time
from pathlib import Path

import numpy as np
import scipy
from benchmark_native_shift import balanced_tree, pectinate_tree

from nwkit import shift_native_screen
from nwkit.shift_native_fit import fit_native_layout
from nwkit.shift_native_model import ShiftData, ShiftLayout
from nwkit.util import read_tree


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tips", type=int, default=1000)
    parser.add_argument("--traits", type=int, default=4)
    parser.add_argument(
        "--shape", choices=["balanced", "pectinate"], default="balanced"
    )
    parser.add_argument("--optimum-increments", action="store_true")
    parser.add_argument("--source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    arrays_path = args.output.with_suffix(".npz")
    if args.tips < 4 or args.traits < 1:
        parser.error(
            "Use at least four tips and one trait for the missing-data fixture."
        )
    if arrays_path == args.output:
        parser.error("JSON output and the .npz matrix output must be distinct.")
    if args.output.exists() or arrays_path.exists():
        parser.error("Use new output paths.")
    module = shift_native_screen
    if args.source:
        spec = importlib.util.spec_from_file_location("screen_reference", args.source)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    newick = (balanced_tree if args.shape == "balanced" else pectinate_tree)(args.tips)
    values = np.random.default_rng(20260911).normal(size=(args.tips, args.traits))
    for trait in range(args.traits):
        values[trait::11, trait] = np.nan
    data = ShiftData.build(
        read_tree(newick, "auto", True, quiet=True),
        values,
        [f"trait{i}" for i in range(args.traits)],
        np.full(values.shape, 0.04),
    )
    fit = fit_native_layout(
        data,
        ShiftLayout.build(data.tree),
        alpha_height=np.geomspace(0.3, 3, args.traits),
        process_variance=1.0,
    )
    measurements = []
    for repetition in range(4):
        started = time.perf_counter()
        matrices, responses, branches = module._whitened_matrices(
            data, fit, 8 * 1024**3, optimum_increments=args.optimum_increments
        )
        measurements.append(
            {"warmup": repetition == 0, "wall_seconds": time.perf_counter() - started}
        )
        if repetition < 3:
            del matrices, responses
    result = {
        "scope": "screening matrix setup only; excludes lasso/search/calibration",
        "configuration": {
            key: str(v) if isinstance(v, Path) else v for key, v in vars(args).items()
        },
        "measurements": measurements,
        "peak_rss_bytes": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        * (1 if platform.system() == "Darwin" else 1024),
        "input_sha256": hashlib.sha256(newick.encode() + values.tobytes()).hexdigest(),
        "source_sha256": hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "thread_limits": {
            key: os.environ.get(key)
            for key in (
                "OPENBLAS_NUM_THREADS",
                "MKL_NUM_THREADS",
                "OMP_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS",
            )
        },
        "platform": platform.platform(),
    }
    np.savez(
        arrays_path,
        branches=branches,
        **{f"x{i}": x for i, x in enumerate(matrices)},
        **{f"y{i}": y for i, y in enumerate(responses)},
    )
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
