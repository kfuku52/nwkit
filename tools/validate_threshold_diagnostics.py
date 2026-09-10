"""Export independent R posterior references (requires R package posterior).

Run from the repository root with PYTHONPATH=. and, if needed, R_LIBS pointing
at the reference-only R library. This never invokes NWKIT diagnostics.
"""

import argparse
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from tests.threshold_diagnostic_support import generate_diagnostic_cases

R_REFERENCE = """
library(posterior)
args <- commandArgs(trailingOnly=TRUE)
x <- as.matrix(read.table(args[1], header=FALSE))
cat(as.character(packageVersion("posterior")), "\\n")
cat(format(c(rhat(x), ess_bulk(x), ess_tail(x), ess_mean(x), mcse_mean(x)),
           digits=17), sep="\\n")
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cases = {}
    draws = generate_diagnostic_cases()
    version = None
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "draws.tsv"
        for name, traces in draws.items():
            np.savetxt(path, traces.T, fmt="%.17g", delimiter="\t")
            completed = subprocess.run(
                ["Rscript", "-e", R_REFERENCE, str(path)],
                check=True,
                capture_output=True,
                text=True,
            )
            lines = completed.stdout.split()
            version = lines[0]
            cases[name] = dict(
                zip(
                    ("rhat", "ess_bulk", "ess_tail", "ess_mean", "mcse_mean"),
                    map(float, lines[1:]),
                    strict=True,
                )
            )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    draws_path = args.output.with_name("threshold_diagnostic_draws.npz")
    np.savez_compressed(draws_path, **draws)
    args.output.write_text(
        json.dumps(
            {
                "reference": f"R posterior {version}",
                "draws_sha256": hashlib.sha256(draws_path.read_bytes()).hexdigest(),
                "cases": cases,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
