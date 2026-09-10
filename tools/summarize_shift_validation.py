"""Export compact, traceable evidence from completed simulation runs."""

import argparse
import hashlib
import json
import tempfile
from collections import defaultdict
from pathlib import Path

import pandas as pd
from shift_simulation_audit import audit_run
from shift_simulation_cases import SCENARIOS, summarize


def write_json(path, value):
    path.write_text(
        json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )


def validate_records(manifest, rows):
    options = manifest["options"]
    cells = (
        len(options["tips"])
        * len(SCENARIOS)
        * 2
        * len(options["standard_errors"])
        * len(options["effects"])
    )
    modes = (
        ("shift", "convergence", "bootstrap")
        if options["bootstrap"]
        else ("shift", "convergence")
    )
    expected = {
        (f"c{cell:03d}-r{replicate:04d}", mode)
        for cell in range(cells)
        for replicate in range(options["replicates"])
        for mode in modes
    }
    actual = [(row["case"], row["mode"]) for row in rows]
    if len(actual) != len(expected) or set(actual) != expected:
        raise ValueError(
            "Run records are incomplete or duplicated relative to the manifest"
        )
    if any(row["status"] not in ("completed", "failed") for row in rows):
        raise ValueError("Unknown fit status")


def _export(runs, output):
    manifests, records, inputs, cells = [], [], {}, []
    sources = {}
    export_hashes = {}
    for path in (
        Path(__file__),
        Path(__file__).with_name("shift_simulation_cases.py"),
        Path(__file__).with_name("shift_simulation_audit.py"),
    ):
        content = path.read_bytes()
        digest = hashlib.sha256(content).hexdigest()
        relative = "tools/" + path.name
        sources[digest] = {"path": relative, "text": content.decode("utf-8")}
        export_hashes[relative] = digest
    for run in runs:
        manifest = json.loads((run / "manifest.json").read_text())
        criterion = manifest["options"]["criterion"]
        manifests.append(manifest)
        for relative, expected in manifest["source_sha256"].items():
            source_path = Path(relative)
            if source_path.is_absolute() or ".." in source_path.parts:
                raise ValueError("Invalid source snapshot path")
            content = (run / "source-snapshot" / relative).read_bytes()
            if hashlib.sha256(content).hexdigest() != expected:
                raise ValueError(f"Source snapshot checksum mismatch: {relative}")
            sources[expected] = {"path": relative, "text": content.decode("utf-8")}
        raw = [
            json.loads(line)
            for line in (run / "records.jsonl").read_text().splitlines()
        ]
        validate_records(manifest, raw)
        audited_inputs, audited_cells = audit_run(run, manifest, raw)
        cells.extend({"criterion": criterion, **row} for row in audited_cells)
        records.extend({**row, "criterion": criterion} for row in raw)
        for case, value in audited_inputs.items():
            if case in inputs and inputs[case] != value:
                raise ValueError(f"Inputs differ between paired runs: {case}")
            inputs[case] = value
    keys = [(row["criterion"], row["case"], row["mode"]) for row in records]
    if len(set(keys)) != len(keys):
        raise ValueError("Duplicate criterion/case/mode records")
    pooled = defaultdict(list)
    for row in records:
        pooled[
            row["criterion"], row["tips"], row["effect"], row["scenario"], row["mode"]
        ].append(row)
    scenarios = []
    for (criterion, tips, effect, scenario, mode), rows in pooled.items():
        summary = summarize(rows)
        # An additional conditional diagnostic: a merge under the exact true
        # shift placement. For null/single/distinct, no merge is warranted.
        correct = [row for row in rows if row.get("exact_edges")]
        summary["merges_given_correct_edges"] = {
            "count": sum(row["any_merge"] for row in correct),
            "denominator": len(correct),
        }
        scenarios.append(
            {
                "criterion": criterion,
                "tips": tips,
                "effect": effect,
                "scenario": scenario,
                "mode": mode,
                **summary,
            }
        )
    output.mkdir(parents=True, exist_ok=False)
    write_json(
        output / "export.json", {"schema_version": 1, "source_sha256": export_hashes}
    )
    write_json(output / "manifests.json", manifests)
    write_json(output / "source-snapshot.json", sources)
    write_json(output / "inputs.json", inputs)
    write_json(output / "records.json", records)
    write_json(output / "cell-summary.json", cells)
    write_json(output / "scenario-summary.json", scenarios)
    pd.DataFrame(records).to_csv(output / "records.csv", index=False)
    lines = [
        "| Criterion | Tips | Effect | Truth | Mode | Completed / attempted | Any shift | Any merge | Ancestry recovered | Shared groups recovered | Mean tip-mean RMSE |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for row in scenarios:
        if row["mode"] == "bootstrap":
            continue
        metrics = [
            f"{row[key]['count']}/{row[key]['denominator']}"
            for key in (
                "any_shift",
                "any_merge",
                "ancestry_recovered",
                "shared_recovered",
            )
        ]
        rmse = row["mean_tip_mean_rmse"]
        value = "NA" if rmse is None else f"{rmse:.3f}"
        lines.append(
            f"| {row['criterion']} | {row['tips']} | {row['effect']} | {row['scenario']} | {row['mode']} | {row['completed']}/{row['attempted']} | "
            + " | ".join(metrics)
            + f" | {value} |"
        )
    (output / "table.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def export(runs, output):
    output = Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Evidence output already exists: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix=".shift-evidence-", dir=output.parent
    ) as temporary:
        staged = Path(temporary) / "bundle"
        _export(runs, staged)
        if output.exists() or output.is_symlink():
            raise FileExistsError(f"Evidence output already exists: {output}")
        staged.rename(output)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("runs", nargs="+", type=Path)
    cli.add_argument("--output", required=True, type=Path)
    options = cli.parse_args()
    export(options.runs, options.output)
