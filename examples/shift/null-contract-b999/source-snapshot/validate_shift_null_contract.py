"""Frozen, cellwise null validation of the actual calibrated search.

Primary cells use independent branch-recursion covariances on balanced and
pectinate trees. Known errors and exact alpha limits form a separate extension.
The reference generator imports no production inference helpers.
"""

import argparse
import gzip
import hashlib
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.stats import beta

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shift_continuous_reference import DenseOUReference  # noqa: E402

from nwkit.shift_calibration import CalibratedSearch  # noqa: E402
from nwkit.util import read_tree  # noqa: E402

ENGINES = {}


def tree_text(n, shape):
    def subtree(labels, height):
        if len(labels) == 1:
            return labels[0], 0.0
        size = len(labels)
        left_size = size // 2 if shape == "balanced" else 1
        child_height = (
            height * (size - 2) / (size - 1)
            if shape == "pectinate"
            else height - 1 / np.log2(n)
        )
        children = []
        for group in (labels[:left_size], labels[left_size:]):
            name, age = (
                subtree(group, child_height) if len(group) > 1 else (group[0], 0.0)
            )
            children.append(f"{name}:{height - age:.17g}")
        return "(" + ",".join(children) + ")", height

    return subtree([f"t{i}" for i in range(n)], 1.0)[0] + ";"


def cells(suite):
    result = []
    for n in (4, 8, 16):
        for shape in ("balanced", "pectinate"):
            settings = (
                [(0.2, False), (2.1, False)]
                if suite == "primary"
                else [(0.0, False), (None, False), (0.2, True), (2.1, True)]
            )
            for alpha, known_error in settings:
                result.append(
                    {
                        "cell_id": len(result),
                        "tips": n,
                        "shape": shape,
                        "alpha_height": alpha,
                        "known_error": known_error,
                        "tree": tree_text(n, shape),
                        "process_tip_variance": 1.0,
                    }
                )
    return result


def execute(case):
    cell, replicate, specification = case
    seed = np.random.SeedSequence([specification["seed"], cell["cell_id"], replicate])
    data_seed, bootstrap_seed = seed.spawn(2)
    boot = int(bootstrap_seed.generate_state(1)[0] % (2**31))
    tree = read_tree(cell["tree"], "auto", True, quiet=True)
    variances = (
        np.geomspace(0.01, 0.25, cell["tips"])
        if cell["known_error"]
        else np.zeros(cell["tips"])
    )
    null = {"groups": [[0]], "shift_branch_ids": []}
    reference = DenseOUReference(tree, null, variances)
    a = np.inf if cell["alpha_height"] is None else cell["alpha_height"]
    _, K = reference.branch_moments(a)
    covariance = cell["process_tip_variance"] * K + np.diag(variances)
    y = np.linalg.cholesky(covariance) @ np.random.default_rng(data_seed).normal(
        size=cell["tips"]
    )
    row = {
        "cell_id": cell["cell_id"],
        "replicate": replicate,
        "observations": y.tolist(),
        "known_variances": variances.tolist(),
        "generating_covariance": covariance.tolist(),
        "bootstrap_seed": boot,
        "lanes": [],
    }
    for convergence in specification["convergence_lanes"]:
        key = cell["tree"], cell["known_error"], convergence
        if key not in ENGINES:
            ENGINES[key] = CalibratedSearch(
                tree, convergence=convergence, variances=variances
            )
        start = time.monotonic()
        lane = {"convergence": convergence}
        try:
            fit = ENGINES[key].fit(
                y, seed=boot, replicates=specification["bootstrap_replicates"]
            )
            lane.update(
                status="completed",
                any_shift=bool(fit["model"]["shift_branch_ids"]),
                fit=fit,
            )
        except (ValueError, np.linalg.LinAlgError) as exc:
            lane.update(status="failed", error=str(exc))
        lane["seconds"] = time.monotonic() - start
        row["lanes"].append(lane)
    return row


def summarize(rows, specification):
    result = []
    comparisons = len(specification["cells"]) * len(specification["convergence_lanes"])
    tail = 0.05 / comparisons
    for cell in specification["cells"]:
        members = [row for row in rows if row["cell_id"] == cell["cell_id"]]
        for convergence in specification["convergence_lanes"]:
            lanes = [
                next(
                    lane for lane in row["lanes"] if lane["convergence"] == convergence
                )
                for row in members
            ]
            n = len(lanes)
            failed = sum(lane["status"] != "completed" for lane in lanes)
            count = sum(lane.get("any_shift", False) for lane in lanes)
            worst = count + failed
            upper = (
                1.0 if worst == n else float(beta.ppf(1 - tail, worst + 1, n - worst))
            )
            result.append(
                {
                    **cell,
                    "convergence": convergence,
                    "attempted": n,
                    "completed": n - failed,
                    "failed": failed,
                    "any_shift": count,
                    "completed_rate": count / (n - failed) if n > failed else None,
                    "worst_case_rate": worst / n,
                    "simultaneous_one_sided_95_upper": upper,
                    "acceptance_upper": specification["acceptance_upper"],
                    "passed": upper <= specification["acceptance_upper"],
                }
            )
    return {
        "phase": specification["phase"],
        "independent_datasets": len(rows),
        "fits": sum(len(row["lanes"]) for row in rows),
        "comparisons": comparisons,
        "cell_results": result,
        "all_cells_pass": all(row["passed"] for row in result),
        "interpretation": "Cellwise finite-design validation; no claim of uniform continuous-parameter error control. Pilot rates must not be used as validation evidence.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--suite", choices=["primary", "extension"], default="primary")
    parser.add_argument(
        "--phase", choices=["pilot", "validation"], default="validation"
    )
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--bootstrap-replicates", type=int, default=199)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument(
        "--cell-ids",
        help="Prespecified comma-separated subset; omitted means all suite cells",
    )
    args = parser.parse_args()
    if args.replicates < 1 or args.workers < 1 or args.bootstrap_replicates < 19:
        parser.error("Positive replicate/worker counts and B >= 19 are required")
    selected = cells(args.suite)
    if args.cell_ids:
        indexes = set(map(int, args.cell_ids.split(",")))
        if not indexes <= {c["cell_id"] for c in selected}:
            parser.error("Unknown cell id")
        selected = [c for c in selected if c["cell_id"] in indexes]
    specification = {
        "schema_version": 1,
        "suite": args.suite,
        "phase": args.phase,
        "seed": 20260921 if args.phase == "pilot" else 20260923,
        "replicates": args.replicates,
        "bootstrap_replicates": args.bootstrap_replicates,
        "workers": args.workers,
        "convergence_lanes": [False, True],
        "acceptance_upper": 0.075,
        "cells": selected,
        "source_sha256": {
            name: hashlib.sha256(Path(name).read_bytes()).hexdigest()
            for name in (
                "nwkit/shift_calibration.py",
                "nwkit/shift_candidates.py",
                "tools/shift_continuous_reference.py",
                "tools/validate_shift_null_contract.py",
            )
        },
        "thread_environment": {
            key: os.environ.get(key)
            for key in (
                "OPENBLAS_NUM_THREADS",
                "OMP_NUM_THREADS",
                "MKL_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS",
            )
        },
    }
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "protocol.json").write_text(
        json.dumps(specification, indent=2) + "\n"
    )
    snapshot = args.output / "source-snapshot"
    snapshot.mkdir()
    for name in specification["source_sha256"]:
        (snapshot / Path(name).name).write_bytes(Path(name).read_bytes())
    cases = [
        (cell, r, specification) for cell in selected for r in range(args.replicates)
    ]
    rows = []
    with (
        ProcessPoolExecutor(max_workers=args.workers) as executor,
        gzip.open(args.output / "records.jsonl.gz", "wt") as stream,
    ):
        for row in executor.map(execute, cases, chunksize=4):
            stream.write(json.dumps(row, allow_nan=False) + "\n")
            rows.append(row)
            if len(rows) % 20 == 0:
                stream.flush()
                print(
                    f"Completed {len(rows)}/{len(cases)} independent datasets",
                    flush=True,
                )
    summary = summarize(rows, specification)
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(
        json.dumps(
            {key: value for key, value in summary.items() if key != "cell_results"}
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
