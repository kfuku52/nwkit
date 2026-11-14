"""Verify the portable evidence itself, without R or the raw model directory."""

import argparse
import gzip
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from shift_alpha_design import cases, generate_case

from nwkit.shift_reference import evaluate_shift_model
from nwkit.util import assign_branch_ids, read_tree


def json_lines(path):
    with gzip.open(path, "rt") as handle:
        yield from (json.loads(line) for line in handle)


def verify(folder):
    specification = json.loads((folder / "protocol.json").read_text())
    design = {r["case_id"]: r for r in cases(specification)}
    inputs = list(json_lines(folder / "inputs.jsonl.gz"))
    assert len(inputs) == len(design)
    assert len({r["case"]["case_id"] for r in inputs}) == len(design)
    trees, truths = {}, {}
    for item in inputs:
        case = item["case"]
        assert case == design[case["case_id"]]
        newick, truth = generate_case(case)
        assert newick == item["tree"] and truth == item["truth"]
        trees[case["case_id"]] = read_tree(newick, "auto", True, quiet=True)
        truths[case["case_id"]] = truth
    records = list(json_lines(folder / "records.jsonl.gz"))
    expected = {
        (key, floor["floor_id"], method, criterion)
        for key in design
        for floor in specification["floors"]
        for method in specification["methods"]
        for criterion in specification["criteria"]
    }
    actual = [
        (r["case_id"], r["floor_id"], r["method"], r["criterion"]) for r in records
    ]
    assert len(actual) == len(expected) and set(actual) == expected
    maximum_ll_error = 0.0
    for record in records:
        if record["status"] != "completed":
            continue
        tree = trees[record["case_id"]]
        truth = truths[record["case_id"]]
        branches = {branch: node for node, branch in assign_branch_ids(tree).items()}
        mean, _, ll = evaluate_shift_model(
            tree,
            observations=truth["observations"],
            standard_errors=[truth["standard_error"]] * truth["tips"],
            alpha=record["alpha"],
            sigma2=record["sigma2"],
            intercept=record["intercept"],
            root_model=record["root_model"],
            mean_effects={
                branches[int(key)]: value
                for key, value in record["mean_effects"].items()
            },
        )
        prediction = {r["leaf_name"]: r["predicted"] for r in record["tip_predictions"]}
        np.testing.assert_allclose(
            mean, [prediction[name] for name in truth["tip_names"]], atol=1e-6, rtol=0
        )
        maximum_ll_error = max(maximum_ll_error, abs(ll - record["log_likelihood"]))
        rmse = float(
            np.sqrt(
                np.mean(
                    [
                        (prediction[name] - truth["tip_mean"][name]) ** 2
                        for name in truth["tip_names"]
                    ]
                )
            )
        )
        assert math.isclose(rmse, record["tip_mean_rmse"], abs_tol=1e-12)
    assert maximum_ll_error < 1e-6
    groups = defaultdict(list)
    fields = ("family", "scenario", "root_model", "floor_id", "criterion", "method")
    for row in records:
        groups[tuple(row[field] for field in fields)].append(row)
    summaries = json.loads((folder / "summary.json").read_text())
    assert len(summaries) == len(groups)
    for row in summaries:
        members = groups[tuple(row[field] for field in fields)]
        good = [r for r in members if r["status"] == "completed"]
        assert row["attempted"] == len(members) and row["completed"] == len(good)
        for metric in (
            "any_shift",
            "shared_recovered",
            "exact_edges",
            "lower_boundary",
            "upper_boundary",
        ):
            count, n = sum(r[metric] for r in good), len(good)
            observed = row[metric]
            assert observed["count"] == count and observed["denominator"] == n
            np.testing.assert_allclose(
                observed["failure_bounds"],
                [count / len(members), (count + len(members) - n) / len(members)],
                atol=1e-12,
            )
            if n:
                p, z = count / n, 1.959963984540054
                assert math.isclose(observed["rate"], p, abs_tol=1e-12)
                center = (p + z * z / (2 * n)) / (1 + z * z / n)
                half = (
                    z
                    * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
                    / (1 + z * z / n)
                )
                np.testing.assert_allclose(
                    observed["wilson_95"], [center - half, center + half], atol=1e-12
                )
        if good:
            assert math.isclose(
                row["mean_rmse"],
                sum(r["tip_mean_rmse"] for r in good) / len(good),
                abs_tol=1e-12,
            )
    comparisons = json.loads((folder / "paired.json").read_text())
    for comparison in comparisons:
        dimension = comparison["dimension"]
        conditions = {
            field: comparison[field] for field in fields if field != dimension
        }
        members = [
            r
            for r in records
            if all(r[field] == value for field, value in conditions.items())
        ]
        left = {r["case_id"]: r for r in members if r[dimension] == comparison["left"]}
        right = {
            r["case_id"]: r for r in members if r[dimension] == comparison["right"]
        }
        assert left.keys() == right.keys() and comparison["attempted"] == len(left)
        good = [
            (left[key], right[key])
            for key in left
            if left[key]["status"] == right[key]["status"] == "completed"
        ]
        assert comparison["completed"] == len(good)
        assert comparison["partition_changed"]["count"] == sum(
            a["shared_partition"] != b["shared_partition"] for a, b in good
        )
        for metric in ("any_shift", "shared_recovered", "exact_edges"):
            assert comparison[metric]["left_only"] == sum(
                a[metric] and not b[metric] for a, b in good
            )
            assert comparison[metric]["right_only"] == sum(
                b[metric] and not a[metric] for a, b in good
            )
        if good:
            assert math.isclose(
                comparison["mean_rmse_right_minus_left"],
                sum(b["tip_mean_rmse"] - a["tip_mean_rmse"] for a, b in good)
                / len(good),
                abs_tol=1e-12,
            )
    candidate_counts = Counter()
    minima = {}
    candidate_ids = defaultdict(set)
    for row in json_lines(folder / "candidate-ledger.jsonl.gz"):
        key = (row["case_id"], row["floor_id"])
        candidate_counts[key] += 1
        assert row["candidate_id"] not in candidate_ids[key]
        candidate_ids[key].add(row["candidate_id"])
        if row["status"] == "completed":
            for criterion in ("BIC", "pBIC"):
                score_key = (*key, criterion)
                minima[score_key] = min(
                    minima.get(score_key, math.inf), float(row[criterion])
                )
    for key, count in candidate_counts.items():
        assert count == (195 if design[key[0]]["tips"] == 8 else 899)
    for record in records:
        if record["method"] == "joint" and record["status"] == "completed":
            assert math.isclose(
                record["score"],
                minima[(record["case_id"], record["floor_id"], record["criterion"])],
                abs_tol=1e-7,
            )
    return {
        "status": "passed",
        "regenerated_datasets": len(inputs),
        "selected_records": len(records),
        "summary_cells": len(summaries),
        "paired_comparisons": len(comparisons),
        "joint_candidate_rows": sum(candidate_counts.values()),
        "max_archived_likelihood_error": maximum_ll_error,
        "checks": [
            "exact regenerated trees/truth/observations",
            "unique complete selected-fit grid",
            "archived parameters reproduce means/likelihood/RMSE",
            "independent counts/means/Wilson formula",
            "unique candidate IDs and selected joint minima",
        ],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folder", type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.folder), indent=2))
