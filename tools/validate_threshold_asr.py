"""Bounded prior-predictive threshold computation audit.

Example: PYTHONPATH=. python tools/validate_threshold_asr.py --replicates 8
--samples 400 --output /tmp/threshold-audit.json. Increase repetitions only for
an explicitly budgeted calibration run. This is not a hypothesis-test audit.
"""

import argparse
import hashlib
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
from scipy.stats import beta

from nwkit.threshold_asr import compute_threshold_marginals
from nwkit.threshold_diagnostics import diagnose_threshold_draws
from nwkit.util import read_tree

CASES = {
    "balanced": ("((A:1,B:1)I:1,(C:1,D:1)J:1)R;", (0.0,), False),
    "pectinate": ("(((A:0.2,B:2)I:0.3,C:1)J:0.7,D:3)R;", (0.0,), False),
    "short_tips": ("((A:0.01,B:0.01)I:0.5,(C:0.01,D:0.01)J:0.5)R;", (0.0,), False),
    "polytomy": ("(A:0.2,B:1,C:3,D:5)R;", (0.0,), False),
    "rare_ordinal": ("((A:0.2,B:0.2)I:0.5,(C:0.2,D:0.2)J:0.5)R;", (0.0, 2.0), False),
    "missing_tip": ("((A:1,B:1)I:1,(C:1,D:1)J:1)R;", (0.0,), True),
}


def binomial_interval(success, total):
    if total == 0:
        return None
    return [
        0.0 if success == 0 else float(beta.ppf(0.025, success, total - success + 1)),
        1.0
        if success == total
        else float(beta.ppf(0.975, success + 1, total - success)),
    ]


def replicate(case, index, samples, burnin, seed):
    source, thresholds, missing = CASES[case]
    tree = read_tree(source, "1", True, quiet=True, rooted="yes")
    # Independent ancestral generator: no NWKIT process or sampler helpers.
    rng = np.random.default_rng(
        np.random.SeedSequence([seed, list(CASES).index(case), index])
    )
    truth = {}
    for node in tree.traverse("preorder"):
        truth[node] = float(
            rng.normal()
            if node.is_root
            else rng.normal(truth[node.up], np.sqrt(node.dist))
        )
    states = tuple(str(i) for i in range(len(thresholds) + 1))
    observed = {}
    likelihoods = {}
    for node in tree.leaves():
        category = int(np.searchsorted(thresholds, truth[node], side="left"))
        observed[node.name] = None if missing and node.name == "D" else states[category]
        likelihoods[node.name] = np.eye(len(states))[category]
    captures = []

    def capture(tree, nodes, states, constraints, traces, *, estimated):
        # Read-only observation of retained draws before their temporary storage
        # is released; sampling and production diagnostics are not changed.
        for node_index, node in enumerate(nodes):
            if not node.is_leaf:
                values = traces[:, :, node_index]
                lo, hi = np.quantile(values, [0.025, 0.975])
                captures.append(
                    {
                        "name": node.name,
                        "truth": truth[node],
                        "lower": float(lo),
                        "upper": float(hi),
                        "covered": bool(lo <= truth[node] <= hi),
                        "posterior_cdf_at_truth": float(np.mean(values <= truth[node])),
                    }
                )
        return diagnose_threshold_draws(
            tree, nodes, states, constraints, traces, estimated=estimated
        )

    result = {"case": case, "replicate": index, "observed": observed}
    try:
        with patch("nwkit.threshold_asr.diagnose_threshold_draws", capture):
            _, fit = compute_threshold_marginals(
                tree,
                states,
                observed,
                likelihoods,
                thresholds=",".join(map(str, thresholds)),
                num_samples=samples,
                burnin=burnin,
                chains=4,
                seed=int(rng.integers(2**32)),
            )
        result.update(
            status=fit.fit_status,
            nodes=captures,
            rhat_max=float(fit.rhat_max),
            ess_min=float(fit.ess_min),
        )
    except (ValueError, FloatingPointError) as exc:
        result.update(status="failed", error=str(exc), nodes=[])
    return result


def summarize(rows):
    output = {}
    for case in CASES:
        group = [row for row in rows if row["case"] == case]
        available = [row for row in group if row["nodes"]]
        passing = [row for row in available if row["status"] == "ok"]

        # Root outcomes are independent between replicates. Ancestors within a
        # tree are not independent binomial replicates and are not pooled here.
        def coverage(subset):
            covered = sum(
                next(node for node in row["nodes"] if node["name"] == "R")["covered"]
                for row in subset
            )
            return {
                "covered": covered,
                "available": len(subset),
                "binomial_95_interval": binomial_interval(covered, len(subset)),
            }

        output[case] = {
            "runs": len(group),
            "failed": len(group) - len(available),
            "diagnostics_passed": len(passing),
            "root_coverage_all_returned": coverage(available),
            "root_coverage_diagnostics_passed": coverage(passing),
        }
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=8)
    parser.add_argument("--samples", type=int, default=400)
    parser.add_argument("--burnin", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260910)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.replicates < 1 or args.samples < 1 or args.burnin < 0:
        parser.error("positive replicates/samples and nonnegative burnin required")
    paths = [
        Path("nwkit") / name
        for name in (
            "threshold_asr.py",
            "mcmc_diagnostics.py",
            "threshold_diagnostics.py",
        )
    ] + [Path(__file__)]
    source_hashes = {
        str(
            path.relative_to(Path.cwd()) if path.is_absolute() else path
        ): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in paths
    }
    rows = []
    for case in CASES:
        for index in range(args.replicates):
            rows.append(replicate(case, index, args.samples, args.burnin, args.seed))
        print(f"{case}: {args.replicates} replicates completed", flush=True)
    current_hashes = {
        str(
            path.relative_to(Path.cwd()) if path.is_absolute() else path
        ): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in paths
    }
    if current_hashes != source_hashes:
        raise RuntimeError(
            "Source changed during the audit; rerun on a stable snapshot."
        )
    payload = {
        "scope": "Small fixed-threshold prior-predictive pilot; not a precise coverage or convergence guarantee",
        "settings": {
            key: value for key, value in vars(args).items() if key != "output"
        },
        "source_sha256": source_hashes,
        "summary": summarize(rows),
        "replicates": rows,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
