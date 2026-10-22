"""Re-evaluate the two conditional prototypes against the integrated library."""

import fcntl
import json
import os
import tempfile
from pathlib import Path

for key in (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "MKL_NUM_THREADS",
):
    os.environ[key] = "1"

from benchmark import fixture_module  # noqa: E402


def main():
    with Path(tempfile.gettempdir(), "nwkit-performance-adoption.lock").open(
        "a"
    ) as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        scalar = fixture_module("06-gaussian-batched-solves")
        multiple = fixture_module("10-observed-covariance")
        result = {
            "scope": "current integrated library vs unchanged isolated conditional prototype",
            "contention": "other user-owned CPU jobs active outside our lock",
            "batched_solves": [
                scalar.benchmark_case(*case)
                for case in (
                    ("dense", 128, 3, False),
                    ("dense-fixed", 256, 5, True),
                    ("dense", 512, 8, True),
                    ("diagonal", 2000, 5, True),
                )
            ],
            "observed_covariance": [
                multiple.full_case(n, 4, fraction)
                for n, fraction in ((64, 0.5), (128, 0.65), (128, 0.9), (256, 0.9))
            ],
        }
        Path(__file__).with_name("conditional-results.json").write_text(
            json.dumps(result, indent=2) + "\n"
        )
        print(
            json.dumps(
                {
                    "solve_ratios": [x["speedup"] for x in result["batched_solves"]],
                    "observed_ratios": [
                        x["speedup"] for x in result["observed_covariance"]
                    ],
                }
            )
        )


if __name__ == "__main__":
    main()
