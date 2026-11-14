"""Verify nested workflows and the largest/slowest pilot on the final source."""

import json
import os
import subprocess
import sys
import time
from pathlib import Path

report = Path(__file__).resolve().parent
sources = json.loads((report / "final-variants.json").read_text())
env = dict(os.environ)
for key in [
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OMP_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
]:
    env[key] = "1"
workflow = [
    "--mode",
    "workflow",
    "--shifts",
    "1",
    "--shared",
    "--pool",
    "8",
    "--refits",
    "8",
    "--screening",
    "200",
]
cases = [
    ("workflow128-before", "before", [*workflow, "--tips", "128"]),
    ("workflow128-after", "after", [*workflow, "--tips", "128"]),
    ("workflow32-after", "after", [*workflow, "--tips", "32"]),
    (
        "fit100-estimated-after",
        "after",
        ["--mode", "fit", "--shifts", "100", "--error", "estimated"],
    ),
    ("search100-after", "after", ["--shifts", "100"]),
    ("search10-known-after", "after", ["--shifts", "10", "--error", "known"]),
]
for name, variant, arguments in cases:
    source = Path(sources[variant]["path"])
    env["PYTHONPATH"] = str(source)
    path = report / "results" / (name + ".json")
    if any(
        path.with_suffix(suffix).exists()
        for suffix in (".json", ".log", ".process.json")
    ):
        raise ValueError(f"Refusing to overwrite {path}")
    command = [
        sys.executable,
        str(source / "tools/benchmark_scale_pilot.py"),
        *arguments,
        "--output",
        str(path),
    ]
    print("START", name, flush=True)
    started = time.perf_counter()
    with path.with_suffix(".log").open("w") as log:
        try:
            result = subprocess.run(
                command, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=600
            )
            status = {"returncode": result.returncode}
        except subprocess.TimeoutExpired:
            status = {"status": "censored", "completed": False, "timeout_seconds": 600}
    status.update(command=command, process_seconds=time.perf_counter() - started)
    path.with_suffix(".process.json").write_text(json.dumps(status, indent=2) + "\n")
    print("END", name, status, flush=True)
