"""Audit diagnostic inputs/fits and replay prespecified bootstrap pairs."""

import argparse
import gzip
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from diagnose_shift_known_error import execute, generate, summarize

from nwkit.shift_calibration import CalibratedSearch
from nwkit.util import assign_branch_ids


def check_row(row):
    text, tree, y, variances, covariance = generate(row["case"])
    if row["tree"] != text:
        raise ValueError("Tree regeneration mismatch")
    for name, expected in (
        ("values", y),
        ("variances", variances),
        ("generating_covariance", covariance),
    ):
        np.testing.assert_array_equal(row[name], expected)
    if row["status"] != "completed":
        raise ValueError("Incomplete diagnostic pair")
    branches = {i: n for n, i in assign_branch_ids(tree).items()}
    eligible = [
        i
        for i, node in branches.items()
        if not node.is_root and 2 <= len(node) <= len(tree) - 2
    ]
    expected_branch = min(
        eligible, key=lambda i: (abs(len(branches[i]) - len(tree) / 2), i)
    )
    if row["shift_branch_id"] != expected_branch:
        raise ValueError("Generating shift location differs from protocol")
    descendants = set(branches[expected_branch].leaves())
    expected_mean = np.array([4.0 if t in descendants else 0.0 for t in tree.leaves()])
    np.testing.assert_array_equal(row["one_shift_mean"], expected_mean)
    search = CalibratedSearch(tree, convergence=False, variances=variances)
    for stage, family, values in (
        ("null", search.families[0], y),
        ("one_shift", search.families[1], y + expected_mean),
    ):
        z = search.q @ (values - values.mean())
        best, at = search.profile(z)
        best, at = best[:, 0], at[:, 0]
        winner = int(family[np.argmax(best[family])])
        item = search.cache[at[winner]]
        _, beta, _ = search._at(z[:, None], item)
        result = row[stage]
        expected_seed = row["case"]["seed"] + (
            8000000000 if stage == "null" else 9000000000
        )
        if result["bootstrap_seed"] != expected_seed:
            raise ValueError("Bootstrap seed differs from protocol")
        np.testing.assert_allclose(
            result["statistic"], 2 * (best.max() - best[winner]), atol=1e-10, rtol=1e-12
        )
        np.testing.assert_allclose(
            result["fitted_contrast_mean"],
            item[2] @ item[3][winner] @ beta[winner, :, 0],
            atol=1e-10,
            rtol=1e-12,
        )
        alpha = search.grid[item[0]]
        if (
            result["fitted_null_model"] != search.models[winner]
            or result["best_alternative_model"] != search.models[int(np.argmax(best))]
            or result["fitted_process_variance"] != float(item[1])
            or result["fitted_alpha_height"]
            != (None if np.isinf(alpha) else float(alpha))
        ):
            raise ValueError("Fitted null/alternative mismatch")
        for method in ("plugin", "oracle"):
            p = result[method + "_p"]
            count = p * (row["case"]["replicates"] + 1)
            if (
                not 0 < p <= 1
                or abs(count - round(count)) > 1e-10
                or result[method + "_reject"] != (p <= 0.05)
            ):
                raise ValueError("Invalid probability/count/decision")
    return row["case"]["case_id"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    root = args.directory
    protocol = json.loads((root / "protocol.json").read_text())
    for source, expected in protocol["source_sha256"].items():
        for path in (Path(source), root / "source-snapshot" / Path(source).name):
            if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
                raise ValueError(f"Source hash mismatch: {path}")
    with gzip.open(root / "records.jsonl.gz", "rt") as stream:
        rows = [json.loads(line) for line in stream]
    if [r["case"] for r in rows] != protocol["cases"]:
        raise ValueError("Case identity/order mismatch")
    if summarize(rows) != json.loads((root / "summary.json").read_text()):
        raise ValueError("Summary mismatch")
    first_by_alpha: dict = {}
    for row in rows:
        first_by_alpha.setdefault(row["case"]["alpha_height"], row)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for start in range(0, len(rows), 8):
            list(pool.map(check_row, rows[start : start + 8]))
            print(
                f"Observed-pair audit: {min(start + 8, len(rows))}/{len(rows)}",
                flush=True,
            )
        selected = list(first_by_alpha.values())
        for original, replay in zip(
            selected, pool.map(execute, [r["case"] for r in selected]), strict=True
        ):
            if original != replay:
                raise ValueError("Bootstrap replay mismatch")
    report = {
        "source_hashes": "active_and_frozen_match",
        "pairs_regenerated_and_refitted": len(rows),
        "bootstrap_pairs_replayed": [r["case"]["case_id"] for r in selected],
        "replay_selection": "First prespecified replicate of each alpha cell; remaining probabilities not independently replayed.",
        "summary": "reconstructed",
        "records_sha256": hashlib.sha256(
            (root / "records.jsonl.gz").read_bytes()
        ).hexdigest(),
        "verifier_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (root / "audit.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
