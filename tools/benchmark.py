"""Compare real CLI startup or Gaussian fits in isolated, single-threaded processes."""

import argparse
import contextlib
import hashlib
import io
import json
import math
import os
import runpy
import subprocess
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CLI_ARGUMENTS = {
    "version": ["--version"],
    "help": ["--help"],
    "regress-help": ["regress", "--help"],
}


def gaussian_workload(tips):
    import numpy as np
    from ete4 import Tree

    from nwkit.ordinary_regression import _fit_ordinary_gaussian

    nodes = []
    for index in range(tips):
        leaf = Tree()
        leaf.name, leaf.dist = f"T{index}", 1.0
        nodes.append(leaf)
    while len(nodes) > 1:
        parents = []
        for index in range(0, len(nodes), 2):
            if index + 1 == len(nodes):
                parents.append(nodes[index])
                continue
            parent = Tree()
            parent.dist = 1.0
            parent.add_child(nodes[index])
            parent.add_child(nodes[index + 1])
            parents.append(parent)
        nodes = parents
    tree = nodes[0]
    tree.dist = 0.0
    names = list(tree.leaf_names())
    rng = np.random.default_rng(1234)
    predictor = rng.normal(size=tips)
    design = np.column_stack([np.ones(tips), predictor])
    response = 0.7 + 1.3 * predictor + rng.normal(size=tips)
    start = time.perf_counter()
    fit = _fit_ordinary_gaussian(
        response,
        design,
        np.zeros(tips),
        tree,
        names,
        evolution_model="lambda",
        evolution_parameter=None,
        branch_length="original",
        custom_covariance=None,
        reml=True,
    )
    elapsed = time.perf_counter() - start
    return elapsed, {
        "beta": fit["beta"].tolist(),
        "beta_covariance": fit["beta_covariance"].tolist(),
        "evolution_parameter": fit["evolution_parameter"],
        "objective": fit["objective"],
        "component_variances": fit["component_variances"],
        "optimizer_converged": bool(fit["optimizer_converged"]),
    }


def peak_rss_bytes():
    try:
        import resource
    except ImportError:
        return None  # resource is unavailable on Windows
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak if sys.platform == "darwin" else peak * 1024


def worker(args):
    sys.path.insert(0, str(args.checkout))
    if args.case == "gaussian":
        elapsed, result = gaussian_workload(args.tips)
    else:
        stream = io.StringIO()
        sys.argv = ["nwkit", *CLI_ARGUMENTS[args.case]]
        start = time.perf_counter()
        with contextlib.redirect_stdout(stream):
            try:
                runpy.run_module("nwkit.cli", run_name="__main__")
            except SystemExit as exc:
                if exc.code not in (None, 0):
                    raise
        elapsed = time.perf_counter() - start
        from nwkit import __version__

        output = stream.getvalue()
        # A required patch-version bump is not a change to command behavior.
        normalized = output.replace(__version__, "<version>")
        result = {"stdout_sha256": hashlib.sha256(normalized.encode()).hexdigest()}
    print(
        json.dumps(
            {
                "wall_seconds": elapsed,
                "peak_rss_bytes": peak_rss_bytes(),
                "result": result,
            }
        )
    )


def equivalent(first, second):
    if isinstance(first, dict) and isinstance(second, dict):
        return first.keys() == second.keys() and all(
            equivalent(first[key], second[key]) for key in first
        )
    if isinstance(first, list) and isinstance(second, list):
        return len(first) == len(second) and all(
            equivalent(a, b) for a, b in zip(first, second, strict=True)
        )
    if isinstance(first, (int, float)) and not isinstance(first, bool):
        return math.isclose(first, second, rel_tol=1e-4, abs_tol=1e-6)
    return first == second


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--case", choices=[*CLI_ARGUMENTS, "gaussian"], default="version"
    )
    parser.add_argument("--checkout", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--repeat", type=int, default=3)
    parser.add_argument("--tips", type=int, default=256)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    args.checkout = args.checkout.resolve()
    if args.tips < 3 or args.repeat < 1:
        parser.error("--tips must be at least 3 and --repeat must be positive")
    if args.worker:
        worker(args)
        return
    environment = os.environ.copy()
    environment.update(
        OPENBLAS_NUM_THREADS="1",
        OMP_NUM_THREADS="1",
        VECLIB_MAXIMUM_THREADS="1",
        MKL_NUM_THREADS="1",
    )
    sources = [("current", args.checkout)]
    if args.baseline is not None:
        sources.insert(0, ("baseline", args.baseline.resolve()))
    records = []
    for repetition in range(args.repeat):
        for label, checkout in sources:
            if not (checkout / "nwkit" / "__init__.py").is_file():
                parser.error(f"Checkout does not contain NWKIT: {checkout}")
            completed = subprocess.run(
                [
                    sys.executable,
                    str(Path(__file__).resolve()),
                    "--worker",
                    "--case",
                    args.case,
                    "--tips",
                    str(args.tips),
                    "--checkout",
                    str(checkout),
                ],
                cwd=checkout,
                env=environment,
                capture_output=True,
                text=True,
            )
            if completed.returncode:
                raise RuntimeError(
                    f"Benchmark failed for {checkout}:\n{completed.stderr}"
                )
            record = json.loads(completed.stdout)
            record.update(label=label, checkout=str(checkout), repetition=repetition)
            records.append(record)
    matched = all(
        equivalent(records[0]["result"], record["result"]) for record in records
    )
    print(
        json.dumps(
            {
                "case": args.case,
                "tips": args.tips if args.case == "gaussian" else None,
                "equivalent": matched,
                "records": records,
            },
            indent=2,
        )
    )
    if not matched:
        raise RuntimeError("Outputs differ beyond the declared numerical tolerances.")


if __name__ == "__main__":
    main()
