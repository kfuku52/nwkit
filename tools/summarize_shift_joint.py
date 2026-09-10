"""Verify and export a small joint-search run without binary R objects."""

import argparse
import hashlib
import json
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from validate_shift_joint import (
    _collect,
    _score_audit,
    _truth_metrics,
    _verify_joint,
    canonical_groups,
    enumerate_candidates,
    read_tree,
    tip_groups,
    write_json,
)


def verify_fit(path, summary):
    tree = read_tree((path / "tree.nwk").read_text(), "auto", True, quiet=True)
    manifest = json.loads((path / "manifest.json").read_text())
    for name, digest in manifest["input_sha256"].items():
        if hashlib.sha256((path / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f"Changed fit input: {path / name}")
    options = manifest["options"]
    candidates, _ = enumerate_candidates(
        tree, options["max_shifts"], options["candidate_limit"]
    )
    # _collect serializes its parsed ledger; use a temporary copy to keep the run read-only.
    with tempfile.TemporaryDirectory() as temporary:
        scratch = Path(temporary)
        shutil.copy2(path / "candidate-results.tsv", scratch)
        ledger = _collect(scratch, candidates)
    good = [row for row in ledger if row["status"] == "completed"]
    winner = min(good, key=lambda row: (row["score"], row["candidate_id"]))
    baseline = json.loads((path / "two-stage.json").read_text())
    verified = _verify_joint(path, tree, baseline, winner)
    if (
        verified != summary["joint"]
        or _score_audit(path, ledger, baseline) != summary["score_audit"]
    ):
        raise ValueError("Joint fit or score audit disagrees with summary")
    if [summary[key] for key in ("attempted", "successful", "failed")] != [
        len(ledger),
        len(good),
        len(ledger) - len(good),
    ]:
        raise ValueError("Incorrect candidate denominators")
    expected_scalars = {
        "criterion": options["criterion"],
        "root_model": options["root_model"],
        "two_stage_score": baseline["parameters"]["score"],
        "joint_score": winner["score"],
        "discrete_coverage_complete": len(good) == len(candidates),
        "continuous_global_optimum_certified": False,
        "score_improvement": baseline["parameters"]["score"] - winner["score"]
        if summary["score_audit"]["comparison_validated"]
        else None,
    }
    if any(summary[key] != value for key, value in expected_scalars.items()):
        raise ValueError("Incorrect comparison summary")
    truth = json.loads((path / "truth.json").read_text())
    groups = [row["branch_ids"] for row in baseline["convergence"]["groups"]]
    fits = [
        (
            "joint",
            winner["shift_branch_ids"],
            winner["groups"],
            verified["tip_predictions"],
        ),
        (
            "two_stage",
            baseline["shift_branch_ids"],
            groups,
            {row["leaf_name"]: row["predicted"] for row in baseline["tip_predictions"]},
        ),
    ]
    for label, selected, aliases, predictions in fits:
        if (
            _truth_metrics(tree, selected, aliases, predictions, truth)
            != summary[f"{label}_truth"]
        ):
            raise ValueError("Incorrect truth recovery metrics")
    truth_scores = [
        row["score"]
        for row in good
        if canonical_groups(tip_groups(tree, row["shift_branch_ids"], row["groups"]))
        == canonical_groups(truth["shared_partition"])
    ]
    np.testing.assert_allclose(min(truth_scores), summary["best_truth_partition_score"])
    return ledger, baseline


def verify_job_coverage(rows, manifest):
    expected = []
    for index, (tips, scenario, _, se) in enumerate(manifest["grid"]):
        criteria = (
            ["BIC", "pBIC"]
            if tips == 8 and scenario == "convergent" and se == 0
            else ["BIC"]
        )
        expected.extend((f"case-{index:02d}", criterion) for criterion in criteria)
    actual = [(row["case"], row["criterion"]) for row in rows]
    if sorted(actual) != sorted(expected):
        raise ValueError("Missing, duplicate or unexpected pilot jobs")


def export(root, output):
    if output.exists():
        raise ValueError("Evidence output already exists")
    rows = json.loads((root / "results.json").read_text())
    manifest = json.loads((root / "manifest.json").read_text())
    verify_job_coverage(rows, manifest)
    ledgers, models, inputs, manifests = {}, {}, {}, {}
    for row in rows:
        key = row["case"].replace("case-", "fit-") + "-" + row["criterion"]
        if row["status"] != "completed":
            raise ValueError(
                "This exporter requires completed runs; inspect the retained failure logs"
            )
        path = root / key
        summary = json.loads((path / "summary.json").read_text())
        if any(row[name] != value for name, value in summary.items()):
            raise ValueError("Parent record disagrees with saved fit")
        ledgers[key], baseline = verify_fit(path, summary)
        models[key] = {
            "two_stage": baseline,
            **{
                name: pd.read_csv(path / f"{name}.tsv", sep="\t").to_dict("records")
                for name in ("joint", "joint-tips", "joint-effects")
            },
        }
        manifests[key] = json.loads((path / "manifest.json").read_text())
        case = root / row["case"]
        if any(
            (case / name).read_bytes() != (path / name).read_bytes()
            for name in ("tree.nwk", "traits.tsv", "truth.json")
        ):
            raise ValueError("Case inputs disagree with fitted inputs")
        inputs[row["case"]] = {
            "tree": (case / "tree.nwk").read_text(),
            "traits_tsv": (case / "traits.tsv").read_text(),
            "truth": json.loads((case / "truth.json").read_text()),
        }
    manifest = json.loads((root / "manifest.json").read_text())
    sources = {}
    for name, digest in manifest["source_sha256"].items():
        source = (root / "source" / name).read_bytes()
        if hashlib.sha256(source).hexdigest() != digest:
            raise ValueError("Execution source snapshot changed")
        sources[name] = source.decode()
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output.parent) as temporary:
        stage = Path(temporary) / "bundle"
        stage.mkdir()
        for name, value in {
            "records": rows,
            "candidates": ledgers,
            "models": models,
            "inputs": inputs,
            "fit-manifests": manifests,
            "manifest": manifest,
            "source-snapshot": sources,
            "export": {
                "source": Path(__file__).read_text(),
                "verification": "Input hashes, candidate enumeration, winner density, score audit, truth metrics and source snapshots; continuous global optima are not certified.",
            },
        }.items():
            write_json(stage / f"{name}.json", value)
        stage.rename(output)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--run", type=Path, required=True)
    cli.add_argument("--output", type=Path, required=True)
    args = cli.parse_args()
    export(args.run.resolve(), args.output.resolve())
