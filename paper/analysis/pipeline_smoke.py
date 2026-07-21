#!/usr/bin/env python3
"""Exercise a representative three-command standard-stream pipeline."""

from __future__ import annotations

import json
from pathlib import Path
import subprocess
import tempfile

from ete4 import Tree
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESULTS = PROJECT_ROOT / "paper" / "results" / "pipeline_smoke.tsv"


def pipe_step(arguments: list[str], input_text: str, audit_path: Path) -> str:
    result = subprocess.run(
        ["nwkit", *arguments, "--audit", str(audit_path)],
        cwd=PROJECT_ROOT,
        input=input_text,
        text=True,
        capture_output=True,
        check=True,
    )
    return result.stdout


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="nwkit-pipeline-") as temporary_directory:
        audit_path = Path(temporary_directory) / "pipeline.audit.jsonl"
        tree_text = "(((T1:1):1,T2:1):1,(T3:1,T4:1):1);"
        tree_text = pipe_step(["sanitize", "--remove-singleton", "yes"], tree_text, audit_path)
        tree_text = pipe_step(
            ["rename", "--pattern", "^T", "--replacement", "Taxon_", "-t", "leaf"],
            tree_text,
            audit_path,
        )
        tree_text = pipe_step(["rescale", "--factor", "2"], tree_text, audit_path)
        audit_records = [json.loads(line) for line in audit_path.read_text().splitlines()]
    tree = Tree(tree_text, parser=1)
    leaf_names = set(tree.leaf_names())
    singleton_count = sum(
        len(list(node.get_children())) == 1
        for node in tree.traverse()
        if not node.is_leaf
    )
    distances = sorted(float(node.dist) for node in tree.traverse() if not node.is_root)
    passed = (
        leaf_names == {"Taxon_1", "Taxon_2", "Taxon_3", "Taxon_4"}
        and singleton_count == 0
        and distances == [2.0, 2.0, 2.0, 2.0, 2.0, 4.0]
        and [record["command"] for record in audit_records] == ["sanitize", "rename", "rescale"]
        and all(record["status"] == "ok" for record in audit_records)
        and all(record["stdin"]["sha256"] for record in audit_records)
    )
    RESULTS.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([
        {
            "commands": "sanitize|rename|rescale",
            "passed": passed,
            "leaf_names": ",".join(sorted(leaf_names)),
            "singleton_nodes": singleton_count,
            "nonroot_branch_lengths": ",".join(str(value) for value in distances),
            "audit_records": len(audit_records),
            "audit_commands": ",".join(record["command"] for record in audit_records),
            "audit_statuses": ",".join(record["status"] for record in audit_records),
            "stdin_hashes_recorded": all(record["stdin"]["sha256"] for record in audit_records),
            "output_newick": tree_text.strip(),
        }
    ]).to_csv(RESULTS, sep="\t", index=False)
    if not passed:
        raise SystemExit("Standard-stream pipeline check failed")


if __name__ == "__main__":
    main()
