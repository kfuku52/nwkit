"""Run predeclared pilot processes sequentially with explicit censoring."""

import json
import os
import subprocess
import sys
import time
from pathlib import Path

baseline = Path(sys.argv[1]).resolve()
output = Path(sys.argv[2]).resolve()
env = dict(os.environ, PYTHONPATH=str(baseline))
for key in [
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OMP_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
]:
    env[key] = "1"
cases = [
    ("fit10-none", ["--mode", "fit", "--shifts", "10"]),
    ("fit100-none", ["--mode", "fit", "--shifts", "100"]),
    ("fit100-known", ["--mode", "fit", "--shifts", "100", "--error", "known"]),
    ("fit100-estimated", ["--mode", "fit", "--shifts", "100", "--error", "estimated"]),
    ("fit100-four", ["--mode", "fit", "--shifts", "100", "--traits", "4"]),
    ("fit10-pectinate", ["--mode", "fit", "--shape", "pectinate", "--shifts", "10"]),
    ("search10", ["--shifts", "10"]),
    ("search100", ["--shifts", "100"]),
    ("search10-shared", ["--shifts", "10", "--shared"]),
    ("search100-shared", ["--shifts", "100", "--shared"]),
    ("search10-known", ["--shifts", "10", "--error", "known"]),
    ("search10-estimated", ["--shifts", "10", "--error", "estimated"]),
    (
        "workflow",
        [
            "--mode",
            "workflow",
            "--tips",
            "32",
            "--shifts",
            "1",
            "--shared",
            "--pool",
            "8",
            "--refits",
            "8",
            "--screening",
            "200",
        ],
    ),
]
for name, arguments in cases:
    path = output / (name + ".json")
    if any(
        path.with_suffix(suffix).exists()
        for suffix in (".json", ".log", ".process.json", ".prof")
    ):
        raise RuntimeError(f"Refusing to reuse artifacts for {path}")
    command = [
        sys.executable,
        str(baseline / "tools/benchmark_scale_pilot.py"),
        *arguments,
        "--output",
        str(path),
    ]
    started = time.perf_counter()
    print("START", name, flush=True)
    with path.with_suffix(".log").open("w") as log:
        try:
            result = subprocess.run(
                command, env=env, stdout=log, stderr=subprocess.STDOUT, timeout=600
            )
            status = {"returncode": result.returncode}
        except subprocess.TimeoutExpired:
            status = {"timeout_seconds": 600, "status": "censored", "completed": False}
    status.update(command=command, process_seconds=time.perf_counter() - started)
    path.with_suffix(".process.json").write_text(json.dumps(status, indent=2) + "\n")
    print("END", name, status, flush=True)
