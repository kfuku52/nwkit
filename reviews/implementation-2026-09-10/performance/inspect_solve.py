"""Record actual slow Cholesky solves; use a disposable working directory."""

import argparse
import contextlib
import json
import sys
import time
from pathlib import Path

import numpy as np
from threadpoolctl import threadpool_info


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.root / "after"))
    import nwkit.gaussian as gaussian
    from nwkit.cli import main as cli

    original = gaussian.cho_solve
    slow = []

    def observed(factor, values):
        start = time.process_time()
        result = original(factor, values)
        elapsed = time.process_time() - start
        if elapsed > 0.001 and len(slow) < 5:
            path = Path("factor" + str(len(slow)) + ".npz").resolve()
            np.savez(path, factor=factor[0], values=values)
            finite = np.abs(factor[0][factor[0] != 0])
            slow.append(
                dict(
                    cpu=elapsed,
                    shape=list(factor[0].shape),
                    rhs_shape=list(values.shape),
                    min_abs=float(finite.min()),
                    max_abs=float(finite.max()),
                    file=str(path),
                )
            )
        return result

    gaussian.cho_solve = observed
    manifest = json.loads((args.root / "manifest.json").read_text())
    case = next(c for c in manifest["cases"] if c["name"] == "regress-lambda-128")
    with open("stdout.txt", "w") as output, contextlib.redirect_stdout(output):
        cli(case["args"])
    args.output.write_text(
        json.dumps(dict(slow=slow, threadpools=threadpool_info()), indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
