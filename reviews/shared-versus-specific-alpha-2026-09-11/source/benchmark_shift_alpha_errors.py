"""Serial estimated-alpha timing supplement with known and estimated errors."""

import argparse
import hashlib
import json
import resource
import signal
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from benchmark_shift_alpha_models import alarm_handler, clean, describe_fit, fixture

from nwkit.shift_native_fit import NativeFitOptions, fit_native_layout
from nwkit.shift_native_model import ShiftData, ShiftLayout
from nwkit.shift_native_provenance import (
    native_implementation,
    native_implementation_sha256,
)
from nwkit.shift_simulation import simulate_shift
from nwkit.shift_simulation_cli import explicit_simulation


def worker(job, path):
    signal.signal(signal.SIGALRM, alarm_handler)
    outcome = {"job": job, "implementation_sha256": native_implementation_sha256()}
    started = time.perf_counter()
    stage = "fixture"
    try:
        original, branch, truth = fixture(job)
        parameters = dict(
            truth["parameters"],
            sampling_standard_errors=[0.1, 0.1],
            measurement_covariance=[[0.04, 0], [0, 0.04]],
        )
        spec = explicit_simulation(original.tree, parameters)
        values, _ = simulate_shift(spec, seed=truth["seed"])
        data = ShiftData.build(
            original.tree, values[0], spec.trait_names, spec.sampling_variances
        )
        outcome["generating"] = dict(
            truth,
            parameters=parameters,
            data_sha256=hashlib.sha256(values.tobytes()).hexdigest(),
        )
        options = NativeFitOptions(
            trait_covariance="full",
            alpha_model=job["mode"],
            estimate_measurement_error=True,
        )
        layout = ShiftLayout.build(data.tree, [branch])
        seconds, fits = [], []
        for run in range(4):
            stage = "warmup" if run == 0 else f"repeat_{run}"
            signal.alarm(120)
            before = time.perf_counter()
            fitted = fit_native_layout(data, layout, options=options)
            elapsed = time.perf_counter() - before
            signal.alarm(0)
            described = describe_fit(data, fitted)
            described["measurement_covariance"] = fitted["joint_covariance"][
                "measurement_covariance"
            ]
            described["identifiability"] = fitted["joint_covariance"]["identifiability"]
            described["measurement_parameter_count"] = fitted["joint_covariance"][
                "measurement_parameter_count"
            ]
            if described["measurement_parameter_count"] != 2:
                raise AssertionError(
                    "Additional measurement variances were not estimated."
                )
            fits.append(described)
            if run:
                seconds.append(elapsed)
        outcome.update(
            status="complete",
            seconds=seconds,
            median_seconds=float(np.median(seconds)),
            fits=fits,
            all_evaluated_modes_succeeded=all(
                f["optimizer"]["evaluated_alpha_modes_succeeded"] for f in fits
            ),
        )
    except TimeoutError as exc:
        outcome.update(status="timeout", stage=stage, error=str(exc))
    except Exception as exc:
        outcome.update(
            status="failed", stage=stage, error=f"{type(exc).__name__}: {exc}"
        )
    finally:
        signal.alarm(0)
    outcome["elapsed_seconds_not_benchmark"] = time.perf_counter() - started
    outcome["worker_peak_rss_kib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    path.write_text(json.dumps(clean(outcome), indent=2, allow_nan=False) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--worker")
    args = parser.parse_args()
    if args.worker:
        worker(json.loads(args.worker), args.output)
        return
    args.output.mkdir(parents=True, exist_ok=True)
    jobs = []
    for truth in ["shared", "different"]:
        modes = (
            ["shared", "trait-specific"]
            if truth == "shared"
            else ["trait-specific", "shared"]
        )
        for mode in modes:
            jobs.append(
                {
                    "kind": "timing_errors",
                    "traits": 2,
                    "truth": truth,
                    "scenario": "opposed",
                    "replicate": 0,
                    "mode": mode,
                    "repeats": 3,
                    "timeout": 120,
                    "known_sampling_se": 0.1,
                    "generating_measurement_variance": 0.04,
                    "estimate_measurement_error": True,
                }
            )
    manifest = {
        "jobs": jobs,
        "implementation": native_implementation(),
        "benchmark_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "base_benchmark_sha256": hashlib.sha256(
            Path(__file__).with_name("benchmark_shift_alpha_models.py").read_bytes()
        ).hexdigest(),
    }
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    rows = []
    for job in jobs:
        path = args.output / f"{job['truth']}-{job['mode']}.json"
        with path.with_suffix(".log").open("w") as log:
            subprocess.run(
                [
                    sys.executable,
                    __file__,
                    "--worker",
                    json.dumps(job),
                    "--output",
                    str(path),
                ],
                stdout=log,
                stderr=log,
                check=True,
            )
        result = json.loads(path.read_text())
        rows.append(result)
        (args.output / "results.json").write_text(json.dumps(rows, indent=2) + "\n")
        print(
            json.dumps(
                {
                    "completed": len(rows),
                    "total": len(jobs),
                    "job": job,
                    "status": result["status"],
                }
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
