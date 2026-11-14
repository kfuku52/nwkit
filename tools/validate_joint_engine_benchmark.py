"""Verify frozen-source engine benchmarks and compare detection outcomes."""

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


def read(path):
    return json.loads(path.read_text())


def check_sources(root, manifest, folder="source"):
    for name, digest in manifest["implementation"]["sources"].items():
        assert (
            hashlib.sha256((root / folder / "nwkit" / name).read_bytes()).hexdigest()
            == digest
        ), name


def check_outcome(row, before, engine="dense_observed_gls"):
    assert row["generating"] == before["generating"]
    if row["status"] != "complete":
        assert row["status"] in {"failed", "timeout"}
        return 0
    missing = np.array(row["generating"]["parameters"]["missing"])
    observations = int(np.any(~missing, axis=1).sum())
    assert row["metadata"]["refit_budget"] == 26 and len(row["records"]) <= 26
    true = set(row["generating"]["true_branches"])
    for fit in row["retained"]:
        assert fit["engine"] == engine
        assert fit["optimizer"]["evaluated_alpha_modes_succeeded"]
        k = len(fit["shift_branch_ids"])
        parameters = 3 * k + (8 if row["job"]["mode"] == "shared" else 9)
        for criterion in ("AIC", "BIC"):
            ic = fit["criteria"][criterion]
            assert ic["parameter_count"] == parameters
            assert ic["sample_size"] == observations
            penalty = parameters * (2 if criterion == "AIC" else math.log(observations))
            assert abs(ic["score"] - (-2 * fit["log_likelihood"] + penalty)) < 1e-8
    for criterion, selected in row["selected"].items():
        branches = set(selected["fit"]["shift_branch_ids"])
        tp, fp, fn = len(true & branches), len(branches - true), len(true - branches)
        assert [selected[k] for k in ("tp", "fp", "fn")] == [tp, fp, fn]
        assert abs(selected["f1"] - 2 * tp / (len(branches) + len(true))) < 1e-12
        assert (
            abs(
                selected["fit"]["criteria"][criterion]["score"]
                - min(f["criteria"][criterion]["score"] for f in row["retained"])
            )
            < 1e-8
        )
    return len(row["retained"])


def compare_part(root, baseline, part, old_part, count):
    manifest = read(root / part / "manifest.json")
    check_sources(root, manifest)
    rows = read(root / part / "results.json")
    assert len(rows) == count
    assert sorted(json.dumps(x["job"], sort_keys=True) for x in rows) == sorted(
        json.dumps(x, sort_keys=True) for x in manifest["jobs"]
    )
    output = []
    fits = 0
    implementation_hash = hashlib.sha256(
        json.dumps(manifest["implementation"], sort_keys=True, allow_nan=False).encode()
    ).hexdigest()
    for row in rows:
        assert row["implementation_sha256"] == implementation_hash
        job = row["job"]
        before = read(
            baseline
            / old_part
            / f"{job['truth']}-{job['replicate']}-{job['mode']}.json"
        )
        fits += check_outcome(row, before)
        record = dict(
            job=job,
            before_status=before["status"],
            after_status=row["status"],
            before_seconds=before.get("search_seconds"),
            after_seconds=row.get("search_seconds"),
            before_peak_rss_kib=before["peak_rss_kib"],
            after_peak_rss_kib=row["peak_rss_kib"],
        )
        if row["status"] == "complete":
            record["after_f1"] = {c: s["f1"] for c, s in row["selected"].items()}
        if before["status"] == row["status"] == "complete":
            record["before_f1"] = {c: s["f1"] for c, s in before["selected"].items()}
            record["same_selected_branches"] = {
                c: row["selected"][c]["fit"]["shift_branch_ids"]
                == before["selected"][c]["fit"]["shift_branch_ids"]
                for c in ("AIC", "BIC")
            }
            record["selected_likelihood_delta"] = {
                c: row["selected"][c]["fit"]["log_likelihood"]
                - before["selected"][c]["fit"]["log_likelihood"]
                for c in ("AIC", "BIC")
            }
        output.append(record)
    return dict(retained_fits_verified=fits, comparisons=output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("baseline", type=Path)
    args = parser.parse_args()
    result = dict(
        accuracy=compare_part(
            args.root, args.baseline, "accuracy-final", "accuracy", 20
        ),
        timing=compare_part(args.root, args.baseline, "timing-final-v2", "timing", 4),
    )
    old = read(args.root / "kernel-baseline.json")
    new = read(args.root / "kernel-optimized.json")
    check_sources(args.root, new)
    check_sources(args.root, old, folder="source-before")
    assert len(old["records"]) == len(new["records"]) == 8
    kernels = []
    for before, after in zip(old["records"], new["records"], strict=True):
        for name in ("tips", "traits", "root", "observed"):
            assert before[name] == after[name]
        errors = {
            name: float(
                np.max(np.abs(np.asarray(before[name]) - np.asarray(after[name])))
            )
            for name in (
                "log_likelihood",
                "coefficients",
                "coefficient_covariance",
                "predicted",
            )
        }
        assert max(errors.values()) < 1e-8
        kernels.append(
            dict(
                tips=after["tips"],
                traits=after["traits"],
                root=after["root"],
                engine=after["engine"],
                errors=errors,
                before_median_seconds=float(np.median(before["seconds"])),
                after_median_seconds=float(np.median(after["seconds"])),
                median_ratio=float(
                    np.median(before["seconds"]) / np.median(after["seconds"])
                ),
            )
        )
    result["kernels"] = kernels
    before = read(args.root / "baseline-current-shared-0.json")
    after = read(args.root / "candidate-current-shared-0.json")
    assert before["status"] == after["status"] == "complete"
    for row, manifest in ((before, old), (after, new)):
        digest = hashlib.sha256(
            json.dumps(
                manifest["implementation"], sort_keys=True, allow_nan=False
            ).encode()
        ).hexdigest()
        assert row["implementation_sha256"] == digest
    check_outcome(before, before, engine="vector_tree_pruning")
    check_outcome(after, before)
    for criterion in ("AIC", "BIC"):
        assert (
            before["selected"][criterion]["fit"]["shift_branch_ids"]
            == after["selected"][criterion]["fit"]["shift_branch_ids"]
        )
    result["immediate_baseline_pair"] = dict(
        before_seconds=before["search_seconds"],
        after_seconds=after["search_seconds"],
        ratio=before["search_seconds"] / after["search_seconds"],
        before_peak_rss_kib=before["peak_rss_kib"],
        after_peak_rss_kib=after["peak_rss_kib"],
        same_selected_branches=True,
        likelihood_delta={
            c: after["selected"][c]["fit"]["log_likelihood"]
            - before["selected"][c]["fit"]["log_likelihood"]
            for c in ("AIC", "BIC")
        },
    )
    (args.root / "validation.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
    print(
        json.dumps(
            {
                "accuracy_fits": result["accuracy"]["retained_fits_verified"],
                "timing_fits": result["timing"]["retained_fits_verified"],
                "kernel_cases": len(kernels),
            }
        )
    )


if __name__ == "__main__":
    main()
