"""Enumerate <=2 shifts and optimum equalities using public kfl1ou fits.

Research reference only: at most 16 tips; no bootstrap or production CLI change.
Run with PYTHONPATH=. from a checkout. Input case needs tree.nwk and traits.tsv
(leaf_name, value, se); optional truth.json adds known-truth recovery metrics.
"""

import argparse
import csv
import hashlib
import json
import math
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from shift_joint_backend import R_SCRIPT
from shift_joint_candidates import canonical_groups, enumerate_candidates, tip_groups

from nwkit.cli import parser
from nwkit.shift_backend_probe import collect_pbic_validation
from nwkit.shift_reference import evaluate_shift_model
from nwkit.util import assign_branch_ids, read_tree


def write_json(path, value):
    path.write_text(
        json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )


def _baseline(case, output, options):
    arguments = [
        "shift",
        "--selection",
        "ic",
        "-i",
        str(case / "tree.nwk"),
        "--trait",
        str(case / "traits.tsv"),
        "--state-column",
        "value",
        "--standard-error-column",
        "se",
        "--criterion",
        options.criterion,
        "--root-model",
        options.root_model,
        "--max-shifts",
        str(options.max_shifts),
        "--search-strategy",
        "exhaustive",
        "--convergence",
        "--rscript",
        options.rscript,
        "--model-out",
        str(output / "two-stage.json"),
        "--fit-out",
        str(output / "two-stage.rds"),
        "-o",
        str(output / "two-stage-regimes.tsv"),
    ]
    args = parser.parse_args(arguments)
    args.handler(args)
    return json.loads((output / "two-stage.json").read_text())


def _write_backend_inputs(output, tree, baseline, candidates):
    ids = assign_branch_ids(tree)
    tokens = {name: token for token, name in baseline["tip_tokens"].items()}
    pd.DataFrame(
        [
            {
                "branch_id": ids[node],
                "clade": "/".join(sorted(tokens[name] for name in node.leaf_names())),
            }
            for node in ids
            if not node.is_root
        ]
    ).to_csv(output / "branches.tsv", sep="\t", index=False)
    pd.DataFrame(
        [
            {
                "candidate_id": row["candidate_id"],
                "shifts": ";".join(map(str, row["shift_branch_ids"])),
                "groups": "|".join(
                    ";".join(map(str, group)) for group in row["groups"]
                ),
            }
            for row in candidates
        ]
    ).to_csv(output / "candidates.tsv", sep="\t", index=False)


def _collect(output, candidates):
    with (output / "candidate-results.tsv").open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if [int(row["candidate_id"]) for row in rows] != list(range(len(candidates))):
        raise ValueError("Joint candidate records are missing, duplicated or reordered")
    results = []
    for row, candidate in zip(rows, candidates, strict=True):
        result = {
            **candidate,
            "status": row["status"],
            "error": row["error"],
            "warnings": row["warnings"],
        }
        if row["status"] == "completed":
            result.update(
                {
                    key: float(row[key])
                    for key in ("score", "log_likelihood", "alpha", "sigma2")
                }
            )
            if not all(
                math.isfinite(result[key])
                for key in ("score", "log_likelihood", "alpha", "sigma2")
            ):
                raise ValueError("Nonfinite successful candidate")
        elif row["status"] != "failed" or not row["error"]:
            raise ValueError("Invalid failed candidate record")
        for key in ("free_score", "free_log_likelihood"):
            result[key] = None if row[key] == "NA" else float(row[key])
        results.append(result)
    write_json(output / "candidate-results.json", results)
    return results


def _verify_joint(output, tree, baseline, candidate):
    metrics = pd.read_csv(output / "joint.tsv", sep="\t").iloc[0].to_dict()
    tips = pd.read_csv(output / "joint-tips.tsv", sep="\t")
    effects = pd.read_csv(output / "joint-effects.tsv", sep="\t")
    if int(metrics["candidate_id"]) != candidate["candidate_id"]:
        raise ValueError("Saved joint fit does not match the best candidate")
    parameters = {key: float(metrics[key]) for key in ("alpha", "sigma2", "intercept")}
    by_name = {
        baseline["tip_tokens"][row["token"]]: row for row in tips.to_dict("records")
    }
    observations = {row["leaf_name"]: row for row in baseline["tip_predictions"]}
    ids = assign_branch_ids(tree)
    nodes = {branch: node for node, branch in ids.items()}
    mean, _, log_likelihood = evaluate_shift_model(
        tree,
        observations=[observations[node.name]["observed"] for node in tree.leaves()],
        standard_errors=[
            observations[node.name]["standard_error"] for node in tree.leaves()
        ],
        root_model=baseline["root_model"],
        mean_effects={
            nodes[int(row.branch_id)]: row.mean_effect for row in effects.itertuples()
        },
        **parameters,
    )
    predictions = [float(by_name[node.name]["predicted"]) for node in tree.leaves()]
    np.testing.assert_allclose(mean, predictions, rtol=1e-6, atol=1e-7)
    np.testing.assert_allclose(
        log_likelihood, candidate["log_likelihood"], rtol=1e-6, atol=1e-6
    )
    for names in tip_groups(tree, candidate["shift_branch_ids"], candidate["groups"]):
        values = [float(by_name[name]["optimum"]) for name in names]
        np.testing.assert_allclose(values, values[0], rtol=1e-7, atol=1e-7)
    return {
        **candidate,
        "parameters": parameters,
        "tip_predictions": {
            name: float(row["predicted"]) for name, row in by_name.items()
        },
        "reference_log_likelihood": log_likelihood,
        "backend_version": str(metrics["backend_version"]),
    }


def _score_audit(output, results, baseline):
    anchor = pd.read_csv(
        output / "baseline-refit.tsv", sep="\t", keep_default_na=False
    ).iloc[0]
    anchor_score = None if anchor["score"] == "NA" else float(anchor["score"])
    comparable = [
        row
        for row in results
        if row["status"] == "completed" and row["free_score"] is not None
    ]
    gaps = [
        row["score"]
        - row["free_score"]
        + 2 * (row["log_likelihood"] - row["free_log_likelihood"])
        for row in comparable
    ]
    max_gap = max(map(abs, gaps)) if gaps else None
    anchor_matches = anchor_score is not None and math.isclose(
        anchor_score, baseline["parameters"]["score"], rel_tol=1e-7, abs_tol=1e-6
    )
    return {
        "baseline_refit_score": anchor_score,
        "baseline_refit_error": str(anchor["error"]),
        "baseline_refit_matches": anchor_matches,
        "unconstrained_equivalence_comparisons": len(gaps),
        "maximum_information_penalty_gap": max_gap,
        "score_equivalence_passed": max_gap is not None and max_gap < 1e-5,
        "comparison_validated": anchor_matches
        and max_gap is not None
        and max_gap < 1e-5,
    }


def _truth_metrics(tree, selected, groups, predictions, truth):
    return {
        "exact_edges": sorted(selected) == truth["shift_branch_ids"],
        "ancestry_recovered": [list(g) for g in tip_groups(tree, selected)]
        == truth["ancestry_partition"],
        "shared_recovered": [list(g) for g in tip_groups(tree, selected, groups)]
        == truth["shared_partition"],
        "tip_mean_rmse": float(
            np.sqrt(
                np.mean(
                    [
                        (predictions[name] - value) ** 2
                        for name, value in truth["tip_mean"].items()
                    ]
                )
            )
        ),
    }


def run(options):
    case = options.case.resolve()
    tree = read_tree((case / "tree.nwk").read_text(), "auto", True, quiet=True)
    candidates, counts = enumerate_candidates(
        tree, options.max_shifts, options.candidate_limit
    )
    executable = shutil.which(options.rscript)
    if executable is None:
        raise ValueError("Rscript executable was not found")
    output = options.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    shutil.copy2(case / "tree.nwk", output / "tree.nwk")
    shutil.copy2(case / "traits.tsv", output / "traits.tsv")
    for name in (
        "validate_shift_joint.py",
        "shift_joint_backend.py",
        "shift_joint_candidates.py",
    ):
        target = output / "source" / name
        target.parent.mkdir(exist_ok=True)
        target.write_bytes(Path(__file__).with_name(name).read_bytes())
    write_json(
        output / "manifest.json",
        {
            "schema_version": 1,
            "options": {**vars(options), "case": str(case), "output": str(output)},
            "input_sha256": {
                name: hashlib.sha256((case / name).read_bytes()).hexdigest()
                for name in ("tree.nwk", "traits.tsv")
            },
            "candidate_generation": counts,
            "candidates": len(candidates),
            "continuous_global_optimum_certified": False,
        },
    )
    baseline = _baseline(output, output, options)
    _write_backend_inputs(output, tree, baseline, candidates)
    (output / "joint.R").write_text(R_SCRIPT)
    with (output / "joint.log").open("w") as log:
        completed = subprocess.run(
            [executable, "--vanilla", str(output / "joint.R")],
            cwd=output,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    results = (
        _collect(output, candidates)
        if (output / "candidate-results.tsv").exists()
        else []
    )
    if completed.returncode:
        raise ValueError(f"Joint backend failed; inspect {output / 'joint.log'}")
    good = [row for row in results if row["status"] == "completed"]
    winner = min(good, key=lambda row: (row["score"], row["candidate_id"]))
    joint = _verify_joint(output, tree, baseline, winner)
    base_groups = [row["branch_ids"] for row in baseline["convergence"]["groups"]]
    matching = [
        row
        for row in good
        if row["shift_branch_ids"] == sorted(baseline["shift_branch_ids"])
        and canonical_groups(row["groups"]) == canonical_groups(base_groups)
    ]
    summary = {
        "criterion": options.criterion,
        "pbic_validation": collect_pbic_validation(output)
        if options.criterion == "pBIC"
        else {"status": "not_applicable", "criterion": options.criterion},
        "root_model": options.root_model,
        "attempted": len(results),
        "successful": len(good),
        "failed": len(results) - len(good),
        "discrete_coverage_complete": len(good) == len(candidates),
        "continuous_global_optimum_certified": False,
        "two_stage_score": baseline["parameters"]["score"],
        "joint_score": joint["score"],
        "score_improvement": baseline["parameters"]["score"] - joint["score"],
        "same_candidate_refit_score": matching[0]["score"] if matching else None,
        "joint": joint,
    }
    summary["score_audit"] = _score_audit(output, results, baseline)
    if not summary["score_audit"]["comparison_validated"]:
        summary["score_improvement"] = None
    if (case / "truth.json").exists():
        truth = json.loads((case / "truth.json").read_text())
        shutil.copy2(case / "truth.json", output / "truth.json")
        summary["joint_truth"] = _truth_metrics(
            tree,
            joint["shift_branch_ids"],
            joint["groups"],
            joint["tip_predictions"],
            truth,
        )
        summary["two_stage_truth"] = _truth_metrics(
            tree,
            baseline["shift_branch_ids"],
            base_groups,
            {row["leaf_name"]: row["predicted"] for row in baseline["tip_predictions"]},
            truth,
        )
        truth_scores = [
            row["score"]
            for row in good
            if [
                list(g)
                for g in tip_groups(tree, row["shift_branch_ids"], row["groups"])
            ]
            == truth["shared_partition"]
        ]
        summary["best_truth_partition_score"] = (
            min(truth_scores) if truth_scores else None
        )
    write_json(output / "summary.json", summary)
    print(
        json.dumps(
            {key: value for key, value in summary.items() if key != "joint"}, indent=2
        )
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--case", type=Path, required=True)
    cli.add_argument("--output", type=Path, required=True)
    cli.add_argument("--rscript", default="Rscript")
    cli.add_argument("--criterion", choices=["BIC", "pBIC", "AICc"], default="BIC")
    cli.add_argument(
        "--root-model", choices=["OUfixedRoot", "OUrandomRoot"], default="OUfixedRoot"
    )
    cli.add_argument("--max-shifts", type=int, choices=[0, 1, 2], default=2)
    cli.add_argument("--candidate-limit", type=int, default=5000)
    run(cli.parse_args())
