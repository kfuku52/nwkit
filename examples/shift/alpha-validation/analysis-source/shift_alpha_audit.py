"""Independent fixed-parameter audits and truth metrics for paired OU fits."""

import csv
import json
import math
from pathlib import Path

import numpy as np
from shift_joint_candidates import tip_groups

from nwkit.shift_reference import evaluate_shift_model
from nwkit.util import assign_branch_ids, read_tree


def read_table(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def ids(value):
    return [int(x) for x in value.split(";") if x]


def audit_case(folder):
    case = json.loads((folder / "case.json").read_text())
    truth = json.loads((folder / "truth.json").read_text())
    tree = read_tree(str(folder / "tree.nwk"), "auto", True, quiet=True)
    branches = assign_branch_ids(tree)
    by_id = {branch: node for node, branch in branches.items()}
    height = max(tree.get_distance(tree, node) for node in tree.leaves())
    settings = {row["floor_id"]: row for row in read_table(folder / "settings.tsv")}
    candidates = read_table(folder / "candidates.tsv")
    ledger = read_table(folder / "candidates-results.tsv")
    if len(ledger) != 2 * len(candidates):
        raise ValueError(f"Incomplete candidate ledger: {folder}")
    models = read_table(folder / "models.tsv")
    expected = {
        (floor, method, criterion)
        for floor in settings
        for method in ("two_stage", "joint")
        for criterion in ("BIC", "pBIC")
    }
    if (
        len(models) != 8
        or {(r["floor_id"], r["method"], r["criterion"]) for r in models} != expected
    ):
        raise ValueError("Missing or duplicate model record")
    for floor in settings:
        rows = [r for r in ledger if r["floor_id"] == floor]
        if sorted(int(r["candidate_id"]) for r in rows) != list(range(len(candidates))):
            raise ValueError("Missing or duplicate candidate")
    output = []
    for raw in models:
        row = {**case, **raw}
        local_ledger = [r for r in ledger if r["floor_id"] == raw["floor_id"]]
        good_candidates = [r for r in local_ledger if r["status"] == "completed"]
        row["candidate_attempted"] = len(local_ledger)
        row["candidate_failed"] = len(local_ledger) - len(good_candidates)
        if raw["status"] != "completed":
            output.append(row)
            continue
        for name in (
            "score",
            "reported_score",
            "alpha",
            "sigma2",
            "intercept",
            "log_likelihood",
        ):
            row[name] = float(raw[name])
            if not math.isfinite(row[name]):
                raise ValueError(f"Nonfinite successful fit: {name}")
        selected = ids(raw["shift_branch_ids"])
        groups = [ids(g) for g in raw["groups"].split("|")]
        if sorted(b for g in groups for b in g) != sorted([0, *selected]):
            raise ValueError("Partition does not cover selected shifts exactly once")
        aliases = {branch: min(group) for group in groups for branch in group}
        labels, effective = {}, []
        for node in tree.traverse("preorder"):
            labels[node] = (
                aliases[branches[node]]
                if branches[node] in aliases
                else labels[node.up]
            )
            if (
                not node.is_root
                and branches[node] in selected
                and labels[node] != labels[node.up]
            ):
                effective.append(branches[node])
        label = "-".join(raw[k] for k in ("floor_id", "method", "criterion"))
        tips = read_table(folder / "models" / f"{label}-tips.tsv")
        if sorted(r["leaf_name"] for r in tips) != sorted(truth["tip_names"]):
            raise ValueError("Tip coverage mismatch")
        predictions = {r["leaf_name"]: float(r["predicted"]) for r in tips}
        optima = {r["leaf_name"]: float(r["optimum"]) for r in tips}
        optimum_groups = {}
        for node in tree.leaves():
            optimum_groups.setdefault(labels[node], []).append(optima[node.name])
        row["optimum_equality_error"] = max(
            max(v) - min(v) for v in optimum_groups.values()
        )
        optimum_tolerance = max(1e-6, max(abs(v) for v in optima.values()) * 1e-10)
        effects = read_table(folder / "models" / f"{label}-effects.tsv")
        row["mean_effects"] = {r["branch_id"]: float(r["mean_effect"]) for r in effects}
        row["tip_predictions"] = [
            {
                "leaf_name": r["leaf_name"],
                "predicted": float(r["predicted"]),
                "optimum": float(r["optimum"]),
            }
            for r in tips
        ]
        mean, _, ll = evaluate_shift_model(
            tree,
            observations=truth["observations"],
            standard_errors=[truth["standard_error"]] * len(tips),
            alpha=row["alpha"],
            sigma2=row["sigma2"],
            intercept=row["intercept"],
            root_model=case["root_model"],
            mean_effects={
                by_id[int(r["branch_id"])]: float(r["mean_effect"]) for r in effects
            },
        )
        row["mean_audit_error"] = float(
            np.max(np.abs(mean - [predictions[n] for n in truth["tip_names"]]))
        )
        row["likelihood_audit_error"] = abs(ll - row["log_likelihood"])
        row["audit_passed"] = (
            row["mean_audit_error"] < 1e-6 and row["likelihood_audit_error"] < 1e-6
        )
        row["audit_passed"] &= row["optimum_equality_error"] <= optimum_tolerance
        if raw["criterion"] == "BIC":
            expected_bic = -2 * ll + (len(groups) + 2 + len(selected)) * math.log(
                len(tips)
            )
            row["bic_audit_error"] = abs(row["score"] - expected_bic)
            row["audit_passed"] &= row["bic_audit_error"] < 1e-6
        row["selected_shifts"] = len(selected)
        row["effective_shift_count"] = len(effective)
        row["any_shift"] = bool(effective)
        row["shared_partition"] = [list(g) for g in tip_groups(tree, selected, groups)]
        row["shared_recovered"] = row["shared_partition"] == truth["shared_partition"]
        row["exact_edges"] = sorted(effective) == truth["shift_branch_ids"]
        row["tip_mean_rmse"] = float(
            np.sqrt(
                np.mean(
                    [
                        (predictions[n] - truth["tip_mean"][n]) ** 2
                        for n in truth["tip_names"]
                    ]
                )
            )
        )
        row["alpha_height_estimate"] = row["alpha"] * height
        setting = settings[raw["floor_id"]]
        row["lower_boundary"] = row["alpha"] <= float(setting["lower"]) * (1 + 1e-4)
        row["upper_boundary"] = row["alpha"] >= float(setting["upper"]) * (1 - 1e-4)
        minimum = min(float(r[raw["criterion"]]) for r in good_candidates)
        row["score_minus_joint_minimum"] = row["score"] - minimum
        signature = (
            tuple(sorted(selected)),
            tuple(sorted(tuple(sorted(g)) for g in groups)),
        )
        matching = [
            c
            for c in candidates
            if (
                tuple(sorted(ids(c["shifts"]))),
                tuple(sorted(tuple(sorted(ids(g))) for g in c["groups"].split("|"))),
            )
            == signature
        ]
        row["baseline_in_candidate_set"] = bool(matching)
        if matching:
            match = next(
                r
                for r in local_ledger
                if r["candidate_id"] == matching[0]["candidate_id"]
            )
            row["same_model_refit_score_gap"] = (
                row["score"] - float(match[raw["criterion"]])
                if match["status"] == "completed"
                else None
            )
        else:
            row["same_model_refit_score_gap"] = None
        if raw["method"] == "joint" and abs(row["score_minus_joint_minimum"]) > 1e-7:
            raise ValueError("Joint selection is not the ledger minimum")
        output.append(row)
    return output
