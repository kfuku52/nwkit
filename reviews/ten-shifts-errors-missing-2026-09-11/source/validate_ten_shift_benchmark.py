"""Audit benchmark inputs, complete outcomes, metrics and source provenance."""

import argparse
import hashlib
import json
import math
import tarfile
from pathlib import Path

import numpy as np


def sha(data):
    return hashlib.sha256(data).hexdigest()


def check_sources(root, part):
    manifest = json.loads((root / part / "manifest.json").read_text())
    with tarfile.open(root / "source/nwkit-source.tar.gz") as archive:
        for name, expected in manifest["implementation"]["sources"].items():
            assert sha(archive.extractfile("nwkit/" + name).read()) == expected, name
    for name, expected in manifest["scripts"].items():
        assert sha((root / "source" / name).read_bytes()) == expected, name
    return manifest


def check_fit(fit, mode, observed_tips):
    k = len(fit["shift_branch_ids"])
    assert len(set(fit["shift_branch_ids"])) == k <= 12
    assert fit["engine"] == "vector_tree_pruning"
    assert fit["optimizer"]["evaluated_alpha_modes_succeeded"]
    parameters = 3 * k + (8 if mode == "shared" else 9)
    for criterion in ("AIC", "BIC"):
        ic = fit["criteria"][criterion]
        assert ic["parameter_count"] == parameters
        assert ic["sample_size"] == observed_tips
        penalty = parameters * (2 if criterion == "AIC" else math.log(observed_tips))
        expected = -2 * fit["log_likelihood"] + penalty
        assert abs(ic["score"] - expected) < 1e-8


def check_selection(row):
    true = set(row["generating"]["true_branches"])
    for criterion, selected in row["selected"].items():
        found = set(selected["fit"]["shift_branch_ids"])
        tp, fp, fn = len(true & found), len(found - true), len(true - found)
        assert (selected["tp"], selected["fp"], selected["fn"]) == (tp, fp, fn)
        assert abs(selected["f1"] - 2 * tp / (len(found) + len(true))) < 1e-12
        assert abs(selected["recall"] - tp / len(true)) < 1e-12
        expected_precision = tp / len(found) if found else 0.0
        assert abs(selected["precision"] - expected_precision) < 1e-12
        assert selected["estimated_shifts"] == len(found)
        assert selected["exact_set"] == (true == found)
        best_score = min(f["criteria"][criterion]["score"] for f in row["retained"])
        assert abs(selected["fit"]["criteria"][criterion]["score"] - best_score) < 1e-8


def check_part(root, part, expected_count):
    manifest = check_sources(root, part)
    rows = json.loads((root / part / "results.json").read_text())
    jobs = [json.dumps(j, sort_keys=True) for j in manifest["jobs"]]
    observed = [json.dumps(r["job"], sort_keys=True) for r in rows]
    assert len(rows) == expected_count
    assert len(set(observed)) == expected_count
    assert set(jobs) == set(observed)
    checks = 0
    for row in rows:
        job = row["job"]
        input_path = root / "inputs" / f"{job['truth']}-{job['replicate']}.npz"
        with np.load(input_path, allow_pickle=False) as data:
            values = data["values"]
            assert values.shape == (100, 2)
            assert sha(values.tobytes()) == row["generating"]["data_sha256"]
            assert np.allclose(data["sampling_variances"], 0.01)
            assert np.array_equal(
                np.isnan(values), row["generating"]["parameters"]["missing"]
            )
            observed_tips = int(np.any(np.isfinite(values), axis=1).sum())
        assert len(set(row["generating"]["true_branches"])) == 10
        if row["status"] == "complete":
            assert row["metadata"]["refit_budget"] == 26
            assert len(row["records"]) <= 26
            for fit in row["retained"]:
                check_fit(fit, job["mode"], observed_tips)
                checks += 1
            check_selection(row)
        else:
            assert row["status"] in ("failed", "timeout"), row
    return dict(outcomes=len(rows), retained_fits_checked=checks)


def check_kernel(root):
    manifest = json.loads((root / "kernel-manifest.json").read_text())
    for name, expected in manifest["scripts"].items():
        assert sha((root / "source" / name).read_bytes()) == expected, name
    rows = json.loads((root / "kernel.json").read_text())
    assert len(rows) == 16
    maximum_error = 0.0
    for row in rows:
        assert row["observed_coordinates"] == 162
        assert len(row["pruning_seconds"]) == len(row["dense_seconds"]) == 7
        assert min(row["pruning_seconds"] + row["dense_seconds"]) > 0
        maximum_error = max(maximum_error, *row["errors"].values())
        assert max(row["errors"].values()) < 1e-8
        covariance = np.asarray(row["covariance_coordinate"])
        assert np.allclose(
            np.asarray(row["measurement_variance"]) / np.diag(covariance), 0.04
        )
    return dict(cases=16, maximum_output_error=maximum_error)


def compare_replays(root):
    accuracy = json.loads((root / "accuracy/results.json").read_text())
    timing = json.loads((root / "timing/results.json").read_text())
    index = {
        (r["job"]["truth"], r["job"]["replicate"], r["job"]["mode"]): r
        for r in accuracy
    }
    compared = []
    for row in timing:
        job = row["job"]
        original = index[(job["truth"], job["replicate"], job["mode"])]
        item = dict(
            truth=job["truth"],
            mode=job["mode"],
            accuracy_status=original["status"],
            timing_status=row["status"],
        )
        if original["status"] == row["status"] == "complete":
            item["criteria"] = {
                c: dict(
                    same_branches=original["selected"][c]["fit"]["shift_branch_ids"]
                    == row["selected"][c]["fit"]["shift_branch_ids"],
                    loglik_difference=row["selected"][c]["fit"]["log_likelihood"]
                    - original["selected"][c]["fit"]["log_likelihood"],
                )
                for c in ("AIC", "BIC")
            }
        compared.append(item)
    return compared


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    report = {
        p: check_part(args.directory, p, n)
        for p, n in [("accuracy", 20), ("timing", 4)]
    }
    report["kernel"] = check_kernel(args.directory)
    report["replays"] = compare_replays(args.directory)
    (args.directory / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
