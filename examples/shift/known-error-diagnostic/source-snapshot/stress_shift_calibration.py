"""Prespecified null stress checks on random ultrametric trees and unequal errors."""

import argparse
import gzip
import hashlib
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shift_simulation_cases import rate  # noqa: E402

from nwkit.shift_calibration import CalibratedSearch  # noqa: E402
from nwkit.util import read_tree  # noqa: E402


def generate(case):
    rng = np.random.default_rng(case["seed"])
    clusters = [(f"t{i}", 0.0) for i in range(case["tips"])]
    age = 0.0
    while len(clusters) > 1:
        k = len(clusters)
        age += rng.exponential(2 / (k * (k - 1)))
        i, j = sorted(rng.choice(k, 2, replace=False), reverse=True)
        left, right = clusters.pop(i), clusters.pop(j)
        clusters.append(
            (f"({left[0]}:{age - left[1]:.17g},{right[0]}:{age - right[1]:.17g})", age)
        )
    tree = read_tree(clusters[0][0] + ";", "auto", True, quiet=True)
    for node in tree.traverse():
        if not node.is_root:
            node.dist /= age

    # Store exact normalized branch lengths before both generation and fitting.
    def newick(node):
        name = (
            node.name
            if node.is_leaf
            else "(" + ",".join(newick(c) for c in node.children) + ")"
        )
        return name + (";" if node.is_root else f":{node.dist:.17g}")

    text = newick(tree)
    tree = read_tree(text, "auto", True, quiet=True)
    nodes = list(tree.traverse("preorder"))
    indexes = {node: i for i, node in enumerate(nodes)}
    loading = {tree: np.zeros(len(nodes))}
    a = case["alpha_height"]
    for node in nodes[1:]:
        if a is None:
            decay, variance = 0.0, 1.0
        elif a == 0:
            decay, variance = 1.0, node.dist
        else:
            decay = np.exp(-a * node.dist)
            variance = -np.expm1(-2 * a * node.dist) / -np.expm1(-2 * a)
        loading[node] = decay * loading[node.up]
        loading[node][indexes[node]] += np.sqrt(variance)
    leaves = list(tree.leaves())
    matrix = np.asarray([loading[t] for t in leaves])
    variances = (
        np.geomspace(0.01, 0.25, len(leaves))
        if case["known_error"]
        else np.zeros(len(leaves))
    )
    covariance = matrix @ matrix.T + np.diag(variances)
    y = matrix @ rng.normal(size=len(nodes)) + np.sqrt(variances) * rng.normal(
        size=len(leaves)
    )
    return text, tree, y, variances, covariance


def execute(case):
    text, tree, y, variances, covariance = generate(case)
    row = {
        "case": case,
        "tree": text,
        "values": y.tolist(),
        "variances": variances.tolist(),
        "generating_covariance": covariance.tolist(),
    }
    try:
        search = CalibratedSearch(tree, convergence=False, variances=variances)
        fit = search.fit(y, seed=case["seed"] + 7000000000, replicates=199)
        row.update(
            status="completed",
            any_shift=bool(fit["model"]["shift_branch_ids"]),
            fit=fit,
        )
    except (ValueError, np.linalg.LinAlgError) as exc:
        row.update(status="failed", error=str(exc))
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--workers", default=2, type=int)
    parser.add_argument("--seed", default=20260916, type=int)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    cases = []
    for n in (4, 8, 16):
        for a in (0, 0.01, 2.0, 100.0, None):
            for rep in range(8):
                index = len(cases)
                seed = int(
                    np.random.SeedSequence([args.seed, index]).generate_state(1)[0]
                )
                cases.append(
                    dict(
                        case_id=index,
                        tips=n,
                        alpha_height=a,
                        known_error=bool(rep % 2),
                        seed=seed,
                    )
                )
    protocol = {
        "master_seed": args.seed,
        "calibration_replicates": 199,
        "convergence": False,
        "cases": cases,
        "source_sha256": {
            p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
            for p in (
                "nwkit/shift_calibration.py",
                "nwkit/shift_candidates.py",
                "tools/stress_shift_calibration.py",
            )
        },
    }
    (args.output / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    snapshot = args.output / "source-snapshot"
    snapshot.mkdir()
    for source in protocol["source_sha256"]:
        (snapshot / Path(source).name).write_bytes(Path(source).read_bytes())
    rows = []
    with (
        ProcessPoolExecutor(max_workers=args.workers) as executor,
        gzip.open(args.output / "records.jsonl.gz", "wt") as stream,
    ):
        for start in range(0, len(cases), 8):
            for row in executor.map(execute, cases[start : start + 8]):
                rows.append(row)
                stream.write(json.dumps(row, allow_nan=False) + "\n")
                stream.flush()
            print(f"{len(rows)}/{len(cases)} stress cases completed", flush=True)
    good = [r for r in rows if r["status"] == "completed"]
    summary = {
        "attempted": len(rows),
        "completed": len(good),
        "false_selection": rate(sum(r["any_shift"] for r in good), len(good)),
        "by_tips": [
            {
                "tips": n,
                "false_selection": rate(
                    sum(r["any_shift"] for r in good if r["case"]["tips"] == n),
                    sum(r["case"]["tips"] == n for r in good),
                ),
            }
            for n in (4, 8, 16)
        ],
        "by_error": [
            {
                "known_error": error,
                "false_selection": rate(
                    sum(
                        r["any_shift"]
                        for r in good
                        if r["case"]["known_error"] == error
                    ),
                    sum(r["case"]["known_error"] == error for r in good),
                ),
            }
            for error in (False, True)
        ],
    }
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
