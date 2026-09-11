"""Check the numeric batching crossover without changing statistical fitting."""

import argparse
import json
import time
from pathlib import Path

import numpy as np
from benchmark_native_shift import balanced_tree, pectinate_tree

from nwkit.compiled_tree import CompiledTree
from nwkit.gaussian_whitening import TreeWhitening
from nwkit.util import read_tree

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", required=True)
args = parser.parse_args()
path = Path(args.output)
if path.exists():
    parser.error("Use a fresh output.")
records = []
for shape, sizes in [
    ("balanced", [16, 32, 64, 128, 256, 512, 1000]),
    ("pectinate", [128, 1000]),
]:
    for size in sizes:
        newick = (balanced_tree if shape == "balanced" else pectinate_tree)(size)
        compiled = CompiledTree.from_tree(read_tree(newick, "auto", True, quiet=True))
        slopes = np.full(len(compiled.nodes), 0.95)
        innovations = np.full(len(compiled.nodes), 0.2)
        errors = np.full(size, 0.03)

        # This callback is invoked synchronously within the current loop iteration.
        def build():
            return TreeWhitening.build(
                compiled,  # noqa: B023
                compiled.leaf_indices,  # noqa: B023
                slopes,  # noqa: B023
                innovations,  # noqa: B023
                errors,  # noqa: B023
                root_variance=0.4,
            )

        build()
        times = []
        for _ in range(3):
            start = time.perf_counter()
            for _ in range(50):
                factor = build()
            times.append((time.perf_counter() - start) / 50)
        values = np.random.default_rng(109).normal(size=(size, 3))
        white = factor.apply(values)
        records.append(
            {
                "shape": shape,
                "tips": size,
                "seconds_per_build": times,
                "log_determinant": factor.log_determinant,
                "whitened_gram": (white.T @ white).tolist(),
            }
        )
path.write_text(json.dumps(records, indent=2, allow_nan=False) + "\n")
