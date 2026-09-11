"""Reuse verified baseline trials and measure the final cutoff adjustment."""

import hashlib
import json
import os
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np
import scipy
from check_equivalence import check

report = Path(__file__).resolve().parent
manifest = json.loads((report / "final-variants.json").read_text())
for variant in manifest.values():
    root = Path(variant["path"])
    actual_sources = {
        str(path.relative_to(root)) for path in (root / "nwkit").glob("*.py")
    }
    if actual_sources != variant["sources"].keys():
        raise ValueError("Measured source inventory differs")
    for relative, expected in variant["sources"].items():
        if hashlib.sha256((root / relative).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Changed measured source: {relative}")
    if (
        hashlib.sha256(
            (root / "tools/benchmark_scale_pilot.py").read_bytes()
        ).hexdigest()
        != variant["harness_sha256"]
    ):
        raise ValueError("Changed harness")
output = report / "final-comparisons"
expected_trials = {
    f"{case}-before-{repeat}.json"
    for case in ("fit100", "fit10-pectinate", "search10", "fit10-known")
    for repeat in range(3)
}
if {path.name for path in output.glob("*-before-*.json")} != expected_trials:
    raise ValueError("Require all twelve baseline trials")
source = Path(manifest["after"]["path"])
env = dict(os.environ, PYTHONPATH=str(source))
for key in [
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OMP_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
]:
    env[key] = "1"
versions = {
    "python": platform.python_version(),
    "numpy": np.__version__,
    "scipy": scipy.__version__,
    "platform": platform.platform(),
}
for path in sorted(output.glob("*-before-*.json")):
    baseline = json.loads(path.read_text())
    if any(baseline[key] != value for key, value in versions.items()):
        raise ValueError("Baseline environment differs")
    if any(baseline["threads"][key] != env[key] for key in baseline["threads"]):
        raise ValueError("Baseline thread configuration differs")
    target = path.with_name(path.name.replace("-before-", "-after-"))
    if target.exists() or target.with_suffix(".log").exists():
        raise ValueError(f"Refusing to overwrite {target}")
    args = baseline["configuration"].copy()
    args["output"] = str(target)
    command = [sys.executable, str(source / "tools/benchmark_scale_pilot.py")]
    for key, value in args.items():
        option = "--" + key.replace("_", "-")
        if isinstance(value, bool):
            if value:
                command.append(option)
        else:
            command.extend([option, str(value)])
    print("START", target.name, flush=True)
    with target.with_suffix(".log").open("w") as log:
        subprocess.run(
            command,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            timeout=600,
            check=True,
        )
    after = json.loads(target.read_text())
    check(baseline, after)
    for key in ("fit", "candidates", "best_unpenalized"):
        if key in baseline and baseline[key] != after[key]:
            raise ValueError(f"Scientific output is not exactly equal: {key}")
    print("PASS", target.name, after["total_seconds"], flush=True)
