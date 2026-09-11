"""Diagnose screening convergence without claiming selection accuracy."""

import argparse
import json
import time
from pathlib import Path
from types import SimpleNamespace

from benchmark_scale_pilot import fixture

from nwkit.shift_native_fit import fit_native_layout
from nwkit.shift_native_model import ShiftLayout
from nwkit.shift_native_screen import group_lasso_screen

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--iterations", type=int, required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()
path = Path(args.output)
if path.exists():
    parser.error("Use a fresh output.")
data, _, digest = fixture(
    SimpleNamespace(
        tips=1000,
        traits=1,
        shifts=100,
        shape="balanced",
        shared=False,
        error="none",
        seed=2026091107,
    )
)
fit = fit_native_layout(data, ShiftLayout.build(data.tree))
started = time.perf_counter()
pool, metadata = group_lasso_screen(
    data, fit, pool_size=128, iterations=args.iterations, memory_limit=512 * 1024**2
)
record = {
    "scope": "screening only; same frozen 100-shift distinct-regime fixture",
    "input_sha256": digest,
    "iterations": args.iterations,
    "seconds": time.perf_counter() - started,
    "pool": list(pool),
    "metadata": metadata,
}
path.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
print(record["seconds"], metadata["all_paths_converged"])
