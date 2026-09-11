"""Alternating fresh-process trials after the sequential pilot completes."""

import json
import os
import subprocess
import sys
from pathlib import Path

output = Path(sys.argv[1]).resolve()
sources = [Path(value).resolve() for value in sys.argv[2:]]
env = dict(os.environ)
for key in [
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OMP_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
]:
    env[key] = "1"
for case, args in [
    ("fit100", ["--mode", "fit", "--shifts", "100"]),
    ("fit10-pectinate", ["--mode", "fit", "--shifts", "10", "--shape", "pectinate"]),
    ("search10", ["--shifts", "10"]),
    ("fit10-known", ["--mode", "fit", "--shifts", "10", "--error", "known"]),
]:
    for repeat in range(3):
        order = list(enumerate(sources))
        if repeat % 2:
            order.reverse()
        for index, source in order:
            variant = ("before", "after")[index]
            path = output / f"{case}-{variant}-{repeat}.json"
            if path.exists() or path.with_suffix(".log").exists():
                raise RuntimeError(f"Refusing to reuse artifacts for {path}")
            env["PYTHONPATH"] = str(source)
            command = [
                sys.executable,
                str(source / "tools/benchmark_scale_pilot.py"),
                *args,
                "--output",
                str(path),
            ]
            print("START", path.name, flush=True)
            with path.with_suffix(".log").open("w") as log:
                subprocess.run(
                    command,
                    env=env,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    timeout=600,
                    check=True,
                )
            result = json.loads(path.read_text())
            print("END", path.name, result["total_seconds"], flush=True)
