#!/usr/bin/env python3
"""Build the Markdown supplement from analysis outputs and source inventory."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from inventory import build_report


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS_DIR = PROJECT_ROOT / "paper" / "results"
DEFAULT_OUTPUT = PROJECT_ROOT / "paper" / "supplement.md"


def cell(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).replace("|", "\\|").replace("\n", " ")


def markdown_table(dataframe: pd.DataFrame) -> str:
    columns = list(dataframe.columns)
    lines = [
        "| " + " | ".join(cell(column) for column in columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in dataframe.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(cell(value) for value in row) + " |")
    return "\n".join(lines)


def status_only(value: object) -> str:
    text = str(value)
    if text.startswith("Native"):
        return "Native"
    if text.startswith("Partial"):
        return "Partial"
    if text.startswith("Not identified"):
        return "—"
    return text


def build_text(results_dir: Path, inventory: dict[str, object]) -> str:
    software = pd.read_csv(PROJECT_ROOT / "paper" / "analysis" / "software_versions.tsv", sep="\t")
    capabilities = pd.read_csv(PROJECT_ROOT / "paper" / "analysis" / "capability_matrix.tsv", sep="\t")
    correctness = pd.read_csv(results_dir / "correctness_checks.tsv", sep="\t")
    pipeline = pd.read_csv(results_dir / "pipeline_smoke.tsv", sep="\t")
    input_cases = pd.read_csv(results_dir / "input_case_checks.tsv", sep="\t")
    benchmark = pd.read_csv(results_dir / "benchmark_summary.tsv", sep="\t")
    hardware = pd.read_csv(results_dir / "hardware.tsv", sep="\t")

    status_matrix = capabilities.copy()
    tool_columns = [column for column in status_matrix.columns if column not in ("task", "decision_rule")]
    for column in tool_columns:
        status_matrix[column] = status_matrix[column].map(status_only)
    status_matrix = status_matrix[["task", *tool_columns]]

    evidence_rows: list[dict[str, str]] = []
    for _, row in capabilities.iterrows():
        for column in tool_columns:
            evidence_rows.append({
                "task": row["task"],
                "software": column.replace("_", " "),
                "classification and evidence note": row[column],
                "decision rule": row["decision_rule"],
            })
    evidence = pd.DataFrame(evidence_rows)

    commands = pd.DataFrame(inventory["functional_subcommands"])
    commands = commands.rename(columns={
        "name": "command",
        "reads_stdin_by_default": "stdin default",
        "writes_stdout_by_default": "stdout default",
        "has_seed_option": "seed option",
    })

    benchmark_display = benchmark.copy()
    numeric_columns = [column for column in benchmark_display.columns if column.endswith(("_median", "_min", "_max"))]
    for column in numeric_columns:
        benchmark_display[column] = benchmark_display[column].map(lambda value: f"{float(value):.4g}")

    sections = [
        "---\ntitle: \"Supplementary Material for NWKIT\"\n---",
        "# Supplementary Methods",
        (
            "The files in `paper/analysis` generate the inventory, independent-output checks, "
            "predefined input corpus, worked PEPC example, benchmarks, and figures. Randomized "
            "analyses use seed 20260716. Raw TSV and JSON results are retained in `paper/results`; "
            "this document contains human-readable summaries."
        ),
        "## Table S1. Software versions and source provenance",
        markdown_table(software),
        "## Table S2. Task-level command-line capability classifications",
        (
            "Native denotes a directly documented CLI workflow; Partial denotes a narrower or "
            "distributed implementation; — means that an equivalent was not identified in the "
            "reviewed version."
        ),
        markdown_table(status_matrix),
        "## Table S3. Capability evidence notes and decision rules",
        markdown_table(evidence),
        "## Table S4. Independent-output checks",
        markdown_table(correctness),
        "## Table S5. Standard-stream pipeline smoke test",
        markdown_table(pipeline),
        "## Table S6. Predefined input corpus",
        markdown_table(input_cases),
        "## Table S7. Benchmark hardware",
        markdown_table(hardware),
        "## Table S8. Benchmark summaries",
        (
            "Values are medians, minima, and maxima across three process-level runs. Peak RSS "
            "includes descendant processes and was sampled every 5 ms."
        ),
        markdown_table(benchmark_display),
        "## Table S9. NWKIT command and stream inventory",
        (
            f"Inventory for NWKIT {inventory['nwkit_version']} at commit "
            f"`{inventory['git_commit_short']}`; {inventory['pytest_test_count']} tests were collected."
        ),
        markdown_table(commands),
    ]
    return "\n\n".join(sections) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    inventory = build_report()
    args.results_dir.mkdir(parents=True, exist_ok=True)
    (args.results_dir / "inventory.json").write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n")
    args.output.write_text(build_text(args.results_dir, inventory))


if __name__ == "__main__":
    main()
