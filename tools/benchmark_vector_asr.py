"""Reproducible missing/noisy vector ASR fit benchmark (run from repository root)."""

import argparse
import json
import platform
import resource
import statistics
import time

import numpy as np

from nwkit.multivariate_gaussian_asr import fit_dense_mvbm
from nwkit.util import read_tree


def workload(tips):
    names = [f"t{i}" for i in range(tips)]
    parts = [name + ":1" for name in names]
    while len(parts) > 1:
        parts = [
            "(" + ",".join(parts[i : i + 2]) + "):1" if i + 1 < len(parts) else parts[i]
            for i in range(0, len(parts), 2)
        ]
    tree = read_tree(parts[0] + ";", "1", True, quiet=True, rooted="yes")
    rng = np.random.default_rng(311)
    values = rng.multivariate_normal([0, 0], [[1, 0.3], [0.3, 2]], size=tips)
    observed = {
        name: [float(values[i, 0]), None if i % 5 == 0 else float(values[i, 1])]
        for i, name in enumerate(names)
    }
    errors = {name: [0.1, 0.2] for name in names}
    return tree, observed, errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backend", choices=("dense", "pruning"), default="dense")
    parser.add_argument("--tips", type=int, default=64)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--model", choices=("BM", "OU"), default="BM")
    args = parser.parse_args()
    tree, observed, errors = workload(args.tips)
    fitter = fit_dense_mvbm
    options = {}
    if args.model == "OU":
        from nwkit.multivariate_gaussian_asr import fit_dense_mvou

        fitter = fit_dense_mvou
        options["alpha"] = 0.4
    if args.backend == "pruning":
        from nwkit.vector_fit import fit_pruning_mvbm

        fitter = fit_pruning_mvbm
        if args.model == "OU":
            from nwkit.vector_ou_fit import fit_pruning_mvou

            fitter = fit_pruning_mvou
    durations = []
    fit = None
    for iteration in range(args.repeats + 1):
        start = time.perf_counter()
        _, fit = fitter(tree, observed, ("x", "y"), standard_errors=errors, **options)
        if iteration:
            durations.append(time.perf_counter() - start)
    print(
        json.dumps(
            {
                "backend": args.backend,
                "model": args.model,
                "tips": args.tips,
                "seed": 311,
                "python": platform.python_version(),
                "platform": platform.platform(),
                "seconds": durations,
                "median_seconds": statistics.median(durations),
                "peak_rss_native_units": resource.getrusage(
                    resource.RUSAGE_SELF
                ).ru_maxrss,
                "sigma": fit.sigma.tolist(),
                "restricted_log_likelihood": fit.restricted_log_likelihood,
                "log_likelihood": fit.log_likelihood,
            }
        )
    )


if __name__ == "__main__":
    main()
