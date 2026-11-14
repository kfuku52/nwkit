"""Audit and summarize the frozen paired OU simulation without dropping failures."""

import argparse
import gzip
import json
import shutil
from collections import defaultdict
from pathlib import Path

import numpy as np
from shift_alpha_audit import audit_case, read_table
from shift_alpha_design import cases
from shift_simulation_cases import rate


def summarize_group(rows):
    good = [r for r in rows if r["status"] == "completed"]
    result = {
        "attempted": len(rows),
        "completed": len(good),
        "failed": len(rows) - len(good),
    }
    for metric in (
        "any_shift",
        "shared_recovered",
        "exact_edges",
        "lower_boundary",
        "upper_boundary",
    ):
        count = sum(r[metric] for r in good)
        result[metric] = rate(count, len(good))
        result[metric]["failure_bounds"] = [
            count / len(rows),
            (count + len(rows) - len(good)) / len(rows),
        ]
    result["mean_rmse"] = (
        float(np.mean([r["tip_mean_rmse"] for r in good])) if good else None
    )
    result["audit_failed"] = sum(not r["audit_passed"] for r in good)
    result["candidate_failed"] = sum(r["candidate_failed"] for r in rows)
    return result


def paired(rows, dimension, left, right):
    identity = ["case_id", "floor_id", "criterion", "method"]
    identity.remove(dimension)
    lookup = {tuple(r[k] for k in identity) + (r[dimension],): r for r in rows}
    pairs = [
        (r, lookup[tuple(r[k] for k in identity) + (right,)])
        for r in rows
        if r[dimension] == left
    ]
    good = [(a, b) for a, b in pairs if a["status"] == b["status"] == "completed"]
    result = {
        "dimension": dimension,
        "left": left,
        "right": right,
        "attempted": len(pairs),
        "completed": len(good),
    }
    result["partition_changed"] = rate(
        sum(a["shared_partition"] != b["shared_partition"] for a, b in good), len(good)
    )
    for key in ("any_shift", "shared_recovered", "exact_edges"):
        result[key] = {
            "left_only": sum(a[key] and not b[key] for a, b in good),
            "right_only": sum(b[key] and not a[key] for a, b in good),
        }
    result["mean_rmse_right_minus_left"] = (
        float(np.mean([b["tip_mean_rmse"] - a["tip_mean_rmse"] for a, b in good]))
        if good
        else None
    )
    if dimension == "method":
        eligible = [
            (a, b)
            for a, b in good
            if a["audit_passed"]
            and b["audit_passed"]
            and not a["candidate_failed"]
            and a["baseline_in_candidate_set"]
            and a["same_model_refit_score_gap"] is not None
            and abs(a["same_model_refit_score_gap"]) <= 1e-4
        ]
        result["search_score_comparison"] = {
            "eligible": len(eligible),
            "excluded": len(good) - len(eligible),
            "joint_lower": sum(b["score"] < a["score"] - 1e-4 for a, b in eligible),
            "tied": sum(abs(b["score"] - a["score"]) <= 1e-4 for a, b in eligible),
            "joint_higher": sum(b["score"] > a["score"] + 1e-4 for a, b in eligible),
        }
    return result


def collect(root, output):
    specification = json.loads((root / "protocol.json").read_text())
    design = list(cases(specification))
    if len({r["seed"] for r in design}) != len(design):
        raise ValueError("Duplicate case seed")
    jobs = [json.loads(line) for line in (root / "jobs.jsonl").read_text().splitlines()]
    if len(jobs) != len(design) or {r["case_id"] for r in jobs} != {
        r["case_id"] for r in design
    }:
        raise ValueError("Not all planned cases have finished")
    rows, ledgers, inputs = [], [], []
    for case in design:
        folder = root / case["case_id"]
        inputs.append(
            {
                "case": case,
                "tree": (folder / "tree.nwk").read_text().strip(),
                "truth": json.loads((folder / "truth.json").read_text()),
            }
        )
        if next(j for j in jobs if j["case_id"] == case["case_id"])["returncode"]:
            rows.extend(
                {
                    **case,
                    "floor_id": floor["floor_id"],
                    "method": method,
                    "criterion": criterion,
                    "status": "failed",
                    "error": "Backend process failed; see raw backend.log",
                    "candidate_failed": 0,
                }
                for floor in specification["floors"]
                for method in specification["methods"]
                for criterion in specification["criteria"]
            )
            continue
        rows.extend(audit_case(folder))
        ledgers.extend(
            {"case_id": case["case_id"], **r}
            for r in read_table(folder / "candidates-results.tsv")
        )
    keys = ("family", "scenario", "root_model", "floor_id", "criterion", "method")
    groups = defaultdict(list)
    for row in rows:
        groups[tuple(row[k] for k in keys)].append(row)
    summaries = [
        {**dict(zip(keys, key, strict=True)), **summarize_group(group)}
        for key, group in sorted(groups.items())
    ]
    comparisons = []
    for dimension, left, right in (
        ("floor_id", "small", "raised"),
        ("criterion", "BIC", "pBIC"),
        ("method", "two_stage", "joint"),
    ):
        paired_keys = [k for k in keys if k != dimension]
        paired_groups = defaultdict(list)
        for row in rows:
            paired_groups[tuple(row[k] for k in paired_keys)].append(row)
        comparisons.extend(
            {
                **dict(zip(paired_keys, key, strict=True)),
                **paired(group, dimension, left, right),
            }
            for key, group in sorted(paired_groups.items())
        )
    good = [r for r in rows if r["status"] == "completed"]
    audit = {
        "datasets_planned": len(design),
        "datasets_finished": len(jobs),
        "point_models_attempted": len(rows),
        "point_models_completed": len(good),
        "point_models_failed": len(rows) - len(good),
        "point_models_audit_failed": sum(not r["audit_passed"] for r in good),
        "candidate_attempted": len(ledgers),
        "candidate_planned": sum(2 * (195 if c["tips"] == 8 else 899) for c in design),
        "candidate_failed": sum(r["status"] != "completed" for r in ledgers),
        "max_mean_error": max(r["mean_audit_error"] for r in good),
        "max_likelihood_error": max(r["likelihood_audit_error"] for r in good),
        "max_bic_error": max(r.get("bic_audit_error", 0) for r in good),
        "max_optimum_equality_error": max(r["optimum_equality_error"] for r in good),
    }
    output.mkdir(parents=True, exist_ok=True)
    for filename, value in (
        ("summary.json", summaries),
        ("paired.json", comparisons),
        ("audit.json", audit),
    ):
        (output / filename).write_text(
            json.dumps(value, indent=2, allow_nan=False) + "\n"
        )
    for filename, values in (
        ("records.jsonl.gz", rows),
        ("candidate-ledger.jsonl.gz", ledgers),
        ("inputs.jsonl.gz", inputs),
    ):
        with gzip.open(output / filename, "wt") as handle:
            for row in values:
                handle.write(json.dumps(row, allow_nan=False) + "\n")
    for filename in (
        "protocol.json",
        "backend-probe.tsv",
        "installed-backend-sha256.json",
        "jobs.jsonl",
    ):
        shutil.copy2(root / filename, output / filename)
    shutil.copytree(root / "source", output / "source", dirs_exist_ok=True)
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    collect(arguments.input, arguments.output)
