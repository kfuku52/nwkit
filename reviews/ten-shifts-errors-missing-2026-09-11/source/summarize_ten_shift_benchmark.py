"""Summarize paired outcomes without silently dropping failures."""

import argparse
import json
from collections import Counter
from pathlib import Path
from statistics import mean, median

METRICS = (
    "tp",
    "fp",
    "fn",
    "precision",
    "recall",
    "f1",
    "estimated_shifts",
    "exact_set",
)


def criterion_summary(records):
    return {key: mean(r[key] for r in records) if records else None for key in METRICS}


def group_summary(all_rows, valid, truth, mode):
    completed = [r for r in all_rows if r["status"] == "complete"]
    group = dict(
        truth=truth,
        mode=mode,
        statuses=dict(Counter(r["status"] for r in all_rows)),
        paired_complete=len(valid),
        criteria={},
        completed_candidate_pool_truth_counts=[
            len(
                set(r["metadata"]["screening"]["pool"])
                & set(r["generating"]["true_branches"])
            )
            for r in completed
        ],
        completed_lasso_paths_converged=[
            r["metadata"]["screening"]["all_paths_converged"] for r in completed
        ],
        operational_recall_with_no_result_as_zero={},
        all_completed_metrics={},
    )
    for criterion in ("AIC", "BIC"):
        group["operational_recall_with_no_result_as_zero"][criterion] = (
            sum(r["selected"][criterion]["recall"] for r in completed) / len(all_rows)
            if all_rows
            else None
        )
        group["all_completed_metrics"][criterion] = [
            dict(
                replicate=r["job"]["replicate"],
                **{k: v for k, v in r["selected"][criterion].items() if k != "fit"},
            )
            for r in completed
        ]
        group["criteria"][criterion] = criterion_summary(
            [p[mode]["selected"][criterion] for p in valid]
        )
    seconds = [p[mode]["search_seconds"] for p in valid]
    group["paired_median_search_seconds"] = median(seconds) if seconds else None
    group["search_seconds"] = seconds
    return group


def collect_pairs(rows):
    paired = {}
    for row in rows:
        j = row["job"]
        pair = paired.setdefault((j["truth"], j["replicate"]), {})
        assert j["mode"] not in pair, j
        pair[j["mode"]] = row
    hash_checks = []
    for key, pair in paired.items():
        if len(pair) == 2 and all("generating" in r for r in pair.values()):
            hashes = {r["generating"]["data_sha256"] for r in pair.values()}
            assert len(hashes) == 1, key
            hash_checks.append(key)
    return paired, hash_checks


def nesting_checks(paired):
    checks = []
    for (truth, replicate), pair in paired.items():
        if len(pair) != 2 or any(r["status"] != "complete" for r in pair.values()):
            continue
        tables = {}
        for mode, row in pair.items():
            tables[mode] = {
                json.dumps([r["shift_branch_ids"], r["groups"]]): r["log_likelihood"]
                for r in row["records"]
                if r["log_likelihood"] is not None
            }
        for key in sorted(tables["shared"].keys() & tables["trait-specific"].keys()):
            delta = tables["trait-specific"][key] - tables["shared"][key]
            checks.append(
                dict(
                    truth=truth,
                    replicate=replicate,
                    layout=json.loads(key),
                    specific_minus_shared_loglik=delta,
                    nesting_violation=delta < -1e-5,
                )
            )
    return checks


def summarize(rows):
    groups = []
    paired, hash_checks = collect_pairs(rows)
    for truth in ("shared", "different"):
        valid = [
            p
            for (t, r), p in paired.items()
            if t == truth
            and len(p) == 2
            and all(v["status"] == "complete" for v in p.values())
        ]
        for mode in ("shared", "trait-specific"):
            all_rows = [
                r
                for r in rows
                if r["job"]["truth"] == truth and r["job"]["mode"] == mode
            ]
            groups.append(group_summary(all_rows, valid, truth, mode))
    return dict(
        groups=groups,
        paired_hash_checks=len(hash_checks),
        planned_observed=len(rows),
        nesting_checks=nesting_checks(paired),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    parts = {}
    for part in ("accuracy", "timing"):
        path = args.directory / part / "results.json"
        if path.exists():
            parts[part] = summarize(json.loads(path.read_text()))
    (args.directory / "summary.json").write_text(json.dumps(parts, indent=2) + "\n")
    print(json.dumps(parts, indent=2))


if __name__ == "__main__":
    main()
