"""Compare equivalent fixed-layout dense and tree GLS in separate processes.

Run each engine with the same size/traits/seed/repeats. These are covariance
and mean-fit measurements, not an end-to-end search or kfl1ou speed claim.
"""

import argparse
import hashlib
import json
import platform
import resource
import statistics
import time
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import solve_triangular

from nwkit.shift_native_model import ShiftData, ShiftLayout, evaluate_trait
from nwkit.util import read_tree


def balanced_tree(size):
    def subtree(start, count):
        if count == 1:
            return f"t{start}", 0
        left, lh = subtree(start, count // 2)
        right, rh = subtree(start + count // 2, count - count // 2)
        height = max(lh, rh) + 1
        return f"({left}:{height - lh},{right}:{height - rh})", height

    return subtree(0, size)[0] + ";"


def pectinate_tree(size):
    subtree = "t0"
    for index in range(1, size):
        subtree = f"({subtree}:1,t{index}:{index})"
    return subtree + ";"


def shared_heights(tree):
    n = len(tree.leaf_names)
    matrix = np.eye(n)
    ranges = {node: (i, i + 1) for i, node in enumerate(tree.compiled.leaf_indices)}
    depths = np.zeros(len(tree.branch_ids))
    for i in range(1, len(depths)):
        depths[i] = depths[tree.compiled.parents[i]] + tree.times[i]
    for i in tree.compiled.postorder:
        children = tree.compiled.children[i]
        if not children:
            continue
        left, right = [ranges[child] for child in children]
        matrix[left[0] : left[1], right[0] : right[1]] = depths[i]
        matrix[right[0] : right[1], left[0] : left[1]] = depths[i]
        ranges[i] = (left[0], right[1])
    return matrix


def dense_fit(data, layout, trait, alpha, shared):
    covariance = (
        np.exp(-2 * alpha * (1 - shared))
        * (-np.expm1(-2 * alpha * shared))
        / -np.expm1(-2 * alpha)
    )
    covariance.flat[:: len(shared) + 1] += data.variances[:, trait] + 0.03
    matrix = np.column_stack((data.values[:, trait], layout.design(data.tree, alpha)))
    factor = np.linalg.cholesky(covariance)
    white = solve_triangular(factor, matrix, lower=True)
    q, r = np.linalg.qr(white[:, 1:])
    beta = np.linalg.solve(r, q.T @ white[:, 0])
    residual = white[:, 0] - white[:, 1:] @ beta
    likelihood = -0.5 * (
        len(shared) * np.log(2 * np.pi)
        + 2 * np.log(np.diag(factor)).sum()
        + residual @ residual
    )
    return [float(likelihood), *beta.tolist()]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", choices=["dense", "tree"], required=True)
    parser.add_argument(
        "--shape", choices=["balanced", "pectinate"], default="balanced"
    )
    parser.add_argument("--tips", type=int, default=128)
    parser.add_argument("--traits", type=int, default=4)
    parser.add_argument("--shifts", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20261001)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)
    if args.tips < args.shifts + 3 or args.traits < 1 or args.repeats < 1:
        parser.error("Require tips > shifts + 2, and positive traits/repeats.")
    destination = Path(args.output)
    if destination.exists():
        parser.error("Use a new output file.")
    start = time.perf_counter()
    newick = (balanced_tree if args.shape == "balanced" else pectinate_tree)(args.tips)
    tree = read_tree(newick, "auto", True, quiet=True)
    y = np.random.default_rng(args.seed).normal(size=(args.tips, args.traits))
    data = ShiftData.build(tree, y, [f"trait{i}" for i in range(args.traits)])
    selected = [
        data.tree.branch_ids[i] for i in data.tree.compiled.leaf_indices[: args.shifts]
    ]
    layout = ShiftLayout.build(data.tree, selected)
    shared = shared_heights(data.tree) if args.engine == "dense" else None
    setup = time.perf_counter() - start
    alphas = np.geomspace(0.2, 20, args.traits)

    def evaluate():
        values = []
        for j, alpha in enumerate(alphas):
            if args.engine == "dense":
                values.append(dense_fit(data, layout, j, alpha, shared))
            else:
                fit = evaluate_trait(data, layout, j, alpha, 1, 0.03)
                values.append([fit.log_likelihood, *fit.coefficients.tolist()])
        return values

    evaluate()  # One warm-up, identical to each timed batch.
    times = []
    for _ in range(args.repeats):
        start = time.perf_counter()
        values = evaluate()
        times.append(time.perf_counter() - start)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_bytes = int(rss if platform.system() == "Darwin" else rss * 1024)
    record = {
        "engine": args.engine,
        "shape": args.shape,
        "tips": args.tips,
        "traits": args.traits,
        "shifts": args.shifts,
        "seed": args.seed,
        "setup_seconds": setup,
        "wall_seconds": times,
        "median_seconds": statistics.median(times),
        "peak_rss_bytes": peak_bytes,
        "values": values,
        "input_sha256": hashlib.sha256(newick.encode() + y.tobytes()).hexdigest(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "scope": "fixed covariance, profiled mean, no search or calibration",
    }
    destination.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    print(
        json.dumps(
            {
                key: record[key]
                for key in (
                    "engine",
                    "tips",
                    "traits",
                    "median_seconds",
                    "peak_rss_bytes",
                )
            }
        )
    )


if __name__ == "__main__":
    main()
