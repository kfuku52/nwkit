"""Host driver for before/after and derivative-order ablation in Docker runtimes.

The timed work runs in fresh containers. Use images built from the snapshots
being compared; no source checkout is overlaid onto their installed packages.
"""

import argparse
import hashlib
import json
import os
import subprocess
from pathlib import Path

import numpy as np


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-image", required=True)
    parser.add_argument("--candidate-image", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workloads", default="4:150,16:1500,64:1500")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--evaluations", type=int, default=10)
    args = parser.parse_args()
    if args.repeats < 2 or args.evaluations < 1 or args.output.exists():
        parser.error("Use a new output, >=2 repetitions and positive evaluations.")
    try:
        workloads = [
            tuple(map(int, pair.split(":"))) for pair in args.workloads.split(",")
        ]
    except ValueError:
        parser.error("Workloads must be integer tips:codons pairs.")
    if any(
        len(pair) != 2 or pair[0] < 4 or pair[0] % 2 or pair[1] < 1
        for pair in workloads
    ):
        parser.error("Use even tip counts >=4 and positive codon counts.")
    trial = Path(__file__).with_name("benchmark_radte_iqtree.py").resolve()
    variants = [
        ("before", args.baseline_image, 2),
        ("cached-full", args.candidate_image, 2),
        ("cached-score", args.candidate_image, 1),
    ]
    images = {
        name: subprocess.check_output(
            ["docker", "image", "inspect", name, "--format", "{{.Id}}"], text=True
        ).strip()
        for name in {args.baseline_image, args.candidate_image}
    }
    report = {
        "images": images,
        "trial_script_sha256": hashlib.sha256(trial.read_bytes()).hexdigest(),
        "evaluations_per_trial": args.evaluations,
        "threads": 1,
        "trials": [],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for tips, sites in workloads:
        reference = None
        for repeat in range(args.repeats):
            # Rotate which variant runs first; retain every observation.
            order = variants[repeat % 3 :] + variants[: repeat % 3]
            for label, image, derivative_order in order:
                command = [
                    "docker",
                    "run",
                    "--rm",
                    "-e",
                    "OPENBLAS_NUM_THREADS=1",
                    "-e",
                    "OMP_NUM_THREADS=1",
                    "-v",
                    f"{trial}:/benchmark.py:ro",
                    "-w",
                    "/tmp",
                    image,
                    "python",
                    "/benchmark.py",
                    "--trial",
                    "--tips",
                    str(tips),
                    "--sites",
                    str(sites),
                    "--mode",
                    "persistent",
                    "--evaluations",
                    str(args.evaluations),
                    "--derivative-order",
                    str(derivative_order),
                ]
                if label != "before":
                    command.append("--worker-stats")
                row = json.loads(
                    subprocess.check_output(command, text=True, env=os.environ)
                )
                row.update(variant=label, image=images[image], repeat=repeat)
                if reference is None:
                    reference = row
                np.testing.assert_allclose(
                    row["nll"], reference["nll"], rtol=0, atol=2e-6
                )
                np.testing.assert_allclose(
                    row["gradient"], reference["gradient"], rtol=2e-7, atol=2e-5
                )
                row["max_nll_difference"] = float(
                    np.max(np.abs(np.array(row["nll"]) - reference["nll"]))
                )
                report["trials"].append(row)
                args.output.write_text(json.dumps(report, indent=2) + "\n")
                print(
                    json.dumps(
                        {k: v for k, v in row.items() if k not in {"nll", "gradient"}}
                    ),
                    flush=True,
                )


if __name__ == "__main__":
    main()
