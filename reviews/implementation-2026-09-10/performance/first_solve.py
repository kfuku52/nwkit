"""Time the first and subsequent identical Cholesky solves in a fresh process."""

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.linalg import cho_solve


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("factor", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    with np.load(args.factor) as data:
        factor, rhs = data["factor"], data["values"]
    records = []
    reference = None
    for index in range(6):
        cpu, wall = time.process_time(), time.perf_counter()
        value = cho_solve((factor, True), rhs)
        records.append(
            dict(
                call=index,
                cpu_seconds=time.process_time() - cpu,
                wall_seconds=time.perf_counter() - wall,
            )
        )
        if reference is None:
            reference = value
        else:
            np.testing.assert_array_equal(value, reference)
    args.output.write_text(json.dumps(records, indent=2) + "\n")


if __name__ == "__main__":
    main()
