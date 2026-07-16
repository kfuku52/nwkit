#!/usr/bin/env python3
"""Reproduce a C4-trait workflow on the public CSUBST PEPC example tree."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import re
import subprocess

from ete4 import Tree
import pandas as pd
import requests


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DATA_DIR = PROJECT_ROOT / "paper" / "data" / "pepc"
DEFAULT_RESULTS_DIR = PROJECT_ROOT / "paper" / "results" / "pepc"
CSUBST_COMMIT = "1b0d98eec71702538e2c09fef60ff1d50dc7e2f4"
FILES = {
    "PEPC.tree.nwk": (
        f"https://raw.githubusercontent.com/kfuku52/csubst/{CSUBST_COMMIT}/csubst/dataset/PEPC.tree.nwk",
        "fe287d9ce9802fa35e6707cb33f187d96f06e0e55d54aadde27348993df7c8f4",
    ),
    "PEPC.foreground.tsv": (
        f"https://raw.githubusercontent.com/kfuku52/csubst/{CSUBST_COMMIT}/csubst/dataset/PEPC.foreground.tsv",
        "51a46eb5634cf8077303a6514536a38b6c1e1d4f0237a6d739449e8074528227",
    ),
}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def fetch_inputs(data_dir: Path) -> dict[str, Path]:
    data_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    for filename, (url, expected_hash) in FILES.items():
        path = data_dir / filename
        if path.exists() and sha256_bytes(path.read_bytes()) == expected_hash:
            paths[filename] = path
            continue
        response = requests.get(url, timeout=60)
        response.raise_for_status()
        if sha256_bytes(response.content) != expected_hash:
            raise ValueError(f"Checksum mismatch for {filename}")
        path.write_bytes(response.content)
        paths[filename] = path
    return paths


def expand_traits(tree_path: Path, foreground_path: Path, output_path: Path) -> pd.DataFrame:
    tree = Tree(tree_path.read_text(), parser=0)
    foreground = pd.read_csv(foreground_path, sep="\t")
    trait_columns = [column for column in foreground.columns if column != "name"]
    rows: list[dict[str, object]] = []
    for leaf_name in tree.leaf_names():
        matches = foreground.loc[
            [re.fullmatch(pattern, leaf_name) is not None for pattern in foreground["name"].astype(str)]
        ]
        if len(matches.index) != 1:
            raise ValueError(f"Expected one foreground pattern for {leaf_name}; observed {len(matches.index)}")
        match = matches.iloc[0]
        row: dict[str, object] = {"leaf_name": leaf_name}
        row.update({column: int(match[column]) for column in trait_columns})
        rows.append(row)
    traits = pd.DataFrame(rows)
    traits.to_csv(output_path, sep="\t", index=False)
    return traits


def run_cli(arguments: list[str]) -> None:
    subprocess.run(["nwkit", *arguments], cwd=PROJECT_ROOT, check=True)


def run_workflow(tree_path: Path, trait_path: Path, results_dir: Path) -> None:
    results_dir.mkdir(parents=True, exist_ok=True)
    run_cli([
        "validate",
        "-i",
        str(tree_path),
        "-o",
        str(results_dir / "pepc_validation.tsv"),
        "--require-rooted",
        "yes",
        "--require-binary",
        "yes",
        "--require-unambiguous-format",
        "yes",
    ])
    run_cli([
        "monophyly",
        "-i",
        str(tree_path),
        "--trait",
        str(trait_path),
        "--group-by",
        "C4",
        "-o",
        str(results_dir / "pepc_c4_monophyly.tsv"),
    ])
    run_cli([
        "asr",
        "-i",
        str(tree_path),
        "--trait",
        str(trait_path),
        "--state-column",
        "C4",
        "--states",
        "0,1",
        "--model",
        "ER",
        "--root-prior",
        "equal",
        "--target",
        "all",
        "--output",
        "probabilities",
        "--model-out",
        str(results_dir / "pepc_c4_asr_model.tsv"),
        "--tree-out",
        str(results_dir / "pepc_c4_asr_tree.nwk"),
        "--tree-annotation",
        "all",
        "-o",
        str(results_dir / "pepc_c4_asr.tsv"),
    ])

    traits = pd.read_csv(trait_path, sep="\t")
    num_c4 = int((traits["C4"] == 1).sum())
    sample_size = min(8, num_c4)
    sample_arguments = [
        "sample",
        "-i",
        str(tree_path),
        "--trait",
        str(trait_path),
        "--method",
        "max-pd",
        "--filter",
        "C4:eq:1",
    ]
    run_cli([
        *sample_arguments,
        "-n",
        str(sample_size),
        "--output-table",
        str(results_dir / "pepc_c4_representatives.tsv"),
        "-o",
        str(results_dir / "pepc_c4_representatives.nwk"),
    ])
    run_cli([
        *sample_arguments,
        "-n",
        str(num_c4),
        "--output-table",
        str(results_dir / "pepc_c4_all_candidates.tsv"),
        "-o",
        str(results_dir / "pepc_c4_all_candidates.nwk"),
    ])
    run_cli([
        "draw",
        "-i",
        str(tree_path),
        "--trait",
        str(trait_path),
        "--group-by",
        "C4",
        "--image_format",
        "svg",
        "-o",
        str(results_dir / "pepc_c4_tree.svg"),
    ])


def summarize(tree_path: Path, trait_path: Path, results_dir: Path) -> pd.DataFrame:
    tree = Tree(tree_path.read_text(), parser=0)
    traits = pd.read_csv(trait_path, sep="\t")
    monophyly = pd.read_csv(results_dir / "pepc_c4_monophyly.tsv", sep="\t")
    asr = pd.read_csv(results_dir / "pepc_c4_asr.tsv", sep="\t")
    model = pd.read_csv(results_dir / "pepc_c4_asr_model.tsv", sep="\t")
    selected = pd.read_csv(results_dir / "pepc_c4_representatives.tsv", sep="\t")
    all_candidates = pd.read_csv(results_dir / "pepc_c4_all_candidates.tsv", sep="\t")

    state_by_branch = dict(zip(asr["branch_id"], asr["map_state"].astype(str)))
    map_shift_count = 0
    for row in asr.itertuples():
        if row.parent == -1:
            continue
        if state_by_branch.get(row.parent) != str(row.map_state):
            map_shift_count += 1
    root_row = asr.loc[asr["node_type"] == "root"].iloc[0]
    pd_coverage = float(selected["pd_total"].iloc[-1] / all_candidates["pd_total"].iloc[-1])
    c4_monophyly = monophyly.loc[monophyly["group"].astype(str) == "1"].iloc[0]
    summary = pd.DataFrame([
        {
            "num_tips": len(list(tree.leaves())),
            "num_c4_tips": int((traits["C4"] == 1).sum()),
            "c4_is_monophyletic": bool(c4_monophyly["is_monophyletic"]),
            "c4_mrca_intruder_tips": int(c4_monophyly["num_intruder_leaves"]),
            "estimated_er_rate": float(model["rate"].iloc[0]),
            "root_p_c4": float(root_row["p_1"]),
            "map_state_shift_edges": map_shift_count,
            "representative_count": len(selected.index),
            "representative_pd_coverage": pd_coverage,
            "selected_c4_tips": ",".join(selected["leaf_name"].astype(str)),
        }
    ])
    summary.to_csv(results_dir / "pepc_summary.tsv", sep="\t", index=False)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    args = parser.parse_args()
    inputs = fetch_inputs(args.data_dir)
    args.results_dir.mkdir(parents=True, exist_ok=True)
    trait_path = args.results_dir / "pepc_traits.tsv"
    expand_traits(inputs["PEPC.tree.nwk"], inputs["PEPC.foreground.tsv"], trait_path)
    run_workflow(inputs["PEPC.tree.nwk"], trait_path, args.results_dir)
    summarize(inputs["PEPC.tree.nwk"], trait_path, args.results_dir)


if __name__ == "__main__":
    main()
