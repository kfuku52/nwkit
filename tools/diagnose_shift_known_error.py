"""Paired plug-in/oracle diagnostics for known-error OU selection (not inference)."""

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
from stress_shift_calibration import generate  # noqa: E402

from nwkit.shift_calibration import CalibratedSearch  # noqa: E402
from nwkit.util import assign_branch_ids  # noqa: E402

SOURCES = (
    "nwkit/shift_calibration.py",
    "nwkit/shift_candidates.py",
    "tools/stress_shift_calibration.py",
    "tools/diagnose_shift_known_error.py",
)


def diagnose(search, y, truth_mean, covariance, family, seed, replicates):
    """Test one family unconditionally, with identical noise for both generators."""
    z = search.q @ (y - y.mean())
    best, at = search.profile(z)
    best, at = best[:, 0], at[:, 0]
    winner = int(family[np.argmax(best[family])])
    item = search.cache[at[winner]]
    _, beta, _ = search._at(z[:, None], item)
    fitted_mean = item[2] @ item[3][winner] @ beta[winner, :, 0]
    observed = float(2 * (best.max() - best[winner]))
    noise = np.random.default_rng(seed).normal(size=(search.d, replicates))
    plugin = search._probability_from_noise(
        fitted_mean, item[2], 1, family, observed, noise
    )
    oracle_L = np.linalg.cholesky(search.q @ covariance @ search.q.T)
    oracle = search._probability_from_noise(
        search.q @ truth_mean, oracle_L, 1, family, observed, noise
    )
    alpha = search.grid[item[0]]
    return {
        "statistic": observed,
        "plugin_p": plugin,
        "oracle_p": oracle,
        "plugin_reject": plugin <= 0.05,
        "oracle_reject": oracle <= 0.05,
        "fitted_alpha_height": None if np.isinf(alpha) else float(alpha),
        "fitted_process_variance": float(item[1]),
        "fitted_null_model": search.models[winner],
        "best_alternative_model": search.models[int(np.argmax(best))],
        "fitted_contrast_mean": fitted_mean.tolist(),
        "bootstrap_seed": seed,
    }


def execute(case):
    text, tree, y, variances, covariance = generate(case)
    leaves = list(tree.leaves())
    ids = assign_branch_ids(tree)
    branches = [
        n for n in tree.traverse() if not n.is_root and 2 <= len(n) <= len(leaves) - 2
    ]
    branch = min(branches, key=lambda n: (abs(len(n) - len(leaves) / 2), ids[n]))
    descendants = set(branch.leaves())
    effect = np.array([4.0 if leaf in descendants else 0.0 for leaf in leaves])
    row = {
        "case": case,
        "tree": text,
        "values": y.tolist(),
        "variances": variances.tolist(),
        "generating_covariance": covariance.tolist(),
        "shift_branch_id": ids[branch],
        "one_shift_mean": effect.tolist(),
    }
    try:
        search = CalibratedSearch(tree, convergence=False, variances=variances)
        row["null"] = diagnose(
            search,
            y,
            np.zeros(len(leaves)),
            covariance,
            search.families[0],
            case["seed"] + 8000000000,
            case["replicates"],
        )
        row["one_shift"] = diagnose(
            search,
            y + effect,
            effect,
            covariance,
            search.families[1],
            case["seed"] + 9000000000,
            case["replicates"],
        )
        row["status"] = "completed"
    except (ValueError, np.linalg.LinAlgError) as exc:
        row.update(status="failed", error=str(exc))
    return row


def summarize(rows):
    cells = []
    for alpha in (0, 0.01, 2.0, 100.0, None):
        members = [r for r in rows if r["case"]["alpha_height"] == alpha]
        good = [r for r in members if r["status"] == "completed"]
        for stage in ("null", "one_shift"):
            cells.append(
                {
                    "alpha_height": alpha,
                    "stage": stage,
                    "attempted": len(members),
                    "completed": len(good),
                    "plugin": rate(
                        sum(r[stage]["plugin_reject"] for r in good), len(good)
                    ),
                    "oracle": rate(
                        sum(r[stage]["oracle_reject"] for r in good), len(good)
                    ),
                    "plugin_only": sum(
                        r[stage]["plugin_reject"] and not r[stage]["oracle_reject"]
                        for r in good
                    ),
                    "oracle_only": sum(
                        r[stage]["oracle_reject"] and not r[stage]["plugin_reject"]
                        for r in good
                    ),
                }
            )
    return {
        "attempted_pairs": len(rows),
        "completed_pairs": sum(r["status"] == "completed" for r in rows),
        "cells": cells,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--workers", default=4, type=int)
    parser.add_argument("--per-alpha", default=20, type=int)
    parser.add_argument("--replicates", default=199, type=int)
    parser.add_argument("--seed", default=20260921, type=int)
    args = parser.parse_args()
    if args.per_alpha < 1 or args.workers < 1 or args.replicates < 19:
        parser.error(
            "Positive workers/cases and at least 19 bootstrap replicates required"
        )
    cases = []
    for alpha in (0, 0.01, 2.0, 100.0, None):
        for _ in range(args.per_alpha):
            index = len(cases)
            seed = int(np.random.SeedSequence([args.seed, index]).generate_state(1)[0])
            cases.append(
                dict(
                    case_id=index,
                    tips=16,
                    alpha_height=alpha,
                    known_error=True,
                    seed=seed,
                    replicates=args.replicates,
                )
            )
    args.output.mkdir(parents=True, exist_ok=False)
    protocol = {
        "master_seed": args.seed,
        "cases": cases,
        "scope": "Paired unconditional tests of no-shift and <=1-shift families; oracle uses generating mean/covariance, unavailable with real data. Not gated full-selection error rates.",
        "process_tip_variance": 1.0,
        "one_shift_effective_tip_mean": 4.0,
        "level": 0.05,
        "source_sha256": {
            s: hashlib.sha256(Path(s).read_bytes()).hexdigest() for s in SOURCES
        },
    }
    (args.output / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    snapshot = args.output / "source-snapshot"
    snapshot.mkdir()
    for source in SOURCES:
        (snapshot / Path(source).name).write_bytes(Path(source).read_bytes())
    rows = []
    with (
        ProcessPoolExecutor(max_workers=args.workers) as executor,
        gzip.open(args.output / "records.jsonl.gz", "wt") as stream,
    ):
        for start in range(0, len(cases), 4):
            for row in executor.map(execute, cases[start : start + 4]):
                rows.append(row)
                stream.write(json.dumps(row, allow_nan=False) + "\n")
                stream.flush()
            print(f"Paired diagnostics: {len(rows)}/{len(cases)}", flush=True)
    (args.output / "summary.json").write_text(
        json.dumps(summarize(rows), indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
