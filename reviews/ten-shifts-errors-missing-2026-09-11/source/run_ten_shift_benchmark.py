"""Durable subprocess controller for the ten-shift benchmark."""

import argparse
import hashlib
import json
import os
import platform
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from nwkit.shift_native_provenance import native_implementation


def run_job(pair):
    job, path = pair
    with path.with_suffix(".log").open("w") as log:
        process = subprocess.run(
            [
                sys.executable,
                str(Path(__file__).with_name("benchmark_ten_shifts_missing.py")),
                "--job",
                json.dumps(job),
                "--output",
                str(path),
            ],
            stdout=log,
            stderr=log,
        )
    if path.exists():
        return json.loads(path.read_text())
    result = dict(job=job, status="worker_failed", returncode=process.returncode)
    path.write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=10)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=int, default=1800)
    parser.add_argument("--traits", type=int, default=2)
    parser.add_argument("--part", choices=["accuracy", "timing"], default="accuracy")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    jobs = []
    for r in range(args.replicates):
        for truth in ("shared", "different"):
            modes = (
                ["shared", "trait-specific"]
                if (r + (truth == "different")) % 2 == 0
                else ["trait-specific", "shared"]
            )
            for mode in modes:
                job = dict(
                    truth=truth,
                    mode=mode,
                    replicate=r,
                    traits=args.traits,
                    missing_rate=0.2,
                    timeout=args.timeout,
                    part=args.part,
                )
                jobs.append((job, args.output / f"{truth}-{r}-{mode}.json"))
    manifest = dict(
        jobs=[j for j, p in jobs],
        implementation=native_implementation(),
        scripts={
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [
                Path(__file__),
                Path(__file__).with_name("benchmark_ten_shifts_missing.py"),
                Path(__file__).with_name("benchmark_shift_alpha_models.py"),
                Path(__file__).with_name("benchmark_shift_covariance.py"),
            ]
        },
        workers=args.workers,
        platform=platform.platform(),
        python=sys.version,
        threads={
            k: os.environ.get(k)
            for k in ["OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"]
        },
    )
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    rows = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        pending = [executor.submit(run_job, j) for j in jobs]
        for future in as_completed(pending):
            row = future.result()
            rows.append(row)
            (args.output / "results.json").write_text(json.dumps(rows, indent=2) + "\n")
            print(
                json.dumps(
                    dict(
                        completed=len(rows),
                        total=len(jobs),
                        job=row["job"],
                        status=row["status"],
                    )
                ),
                flush=True,
            )


if __name__ == "__main__":
    main()
