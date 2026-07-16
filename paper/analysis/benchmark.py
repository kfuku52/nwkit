#!/usr/bin/env python3
"""Measure NWKIT scaling and consensus operating costs across CLI toolkits."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import platform
import random
import subprocess
import tempfile
import time

import pandas as pd
import psutil

from evaluate_correctness import random_rooted_tree


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS_DIR = PROJECT_ROOT / "paper" / "results"


def balanced_tree(start: int, end: int) -> str:
    if end - start == 1:
        return f"T{start}"
    midpoint = (start + end) // 2
    return f"({balanced_tree(start, midpoint)}:1,{balanced_tree(midpoint, end)}:1)"


def caterpillar_tree(num_tips: int) -> str:
    current = "T0"
    for index in range(1, num_tips):
        current = f"({current}:1,T{index}:1)"
    return current


def peak_process_rss(process: psutil.Process) -> int:
    total = 0
    try:
        total += process.memory_info().rss
        for child in process.children(recursive=True):
            try:
                total += child.memory_info().rss
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        pass
    return total


def measure_command(command: list[str], cwd: Path) -> dict[str, object]:
    with tempfile.TemporaryFile(mode="w+") as error_handle:
        start = time.perf_counter()
        child = subprocess.Popen(
            command,
            cwd=cwd,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=error_handle,
        )
        process = psutil.Process(child.pid)
        peak_rss = 0
        while child.poll() is None:
            peak_rss = max(peak_rss, peak_process_rss(process))
            time.sleep(0.005)
        peak_rss = max(peak_rss, peak_process_rss(process))
        elapsed = time.perf_counter() - start
        error_handle.seek(0)
        error_text = error_handle.read().strip().replace("\n", " | ")
    return {
        "elapsed_s": elapsed,
        "peak_rss_mb": peak_rss / (1024 * 1024),
        "returncode": child.returncode,
        "error": error_text[:500],
    }


def hardware_metadata() -> pd.DataFrame:
    values = {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "python": platform.python_version(),
        "logical_cpus": psutil.cpu_count(logical=True),
        "physical_cpus": psutil.cpu_count(logical=False),
        "total_memory_gb": psutil.virtual_memory().total / (1024 ** 3),
    }
    return pd.DataFrame([{"key": key, "value": value} for key, value in values.items()])


def benchmark_validation(directory: Path, repetitions: int) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for num_tips in (128, 512, 2048, 8192, 32768):
        for shape, tree_string in (
            ("balanced", balanced_tree(0, num_tips)),
            ("caterpillar", caterpillar_tree(num_tips)),
        ):
            input_path = directory / f"{shape}-{num_tips}.nwk"
            input_path.write_text(tree_string + ";\n")
            command = [
                "nwkit",
                "validate",
                "-i",
                str(input_path),
                "-o",
                os.devnull,
                "--require-binary",
                "yes",
            ]
            for replicate in range(1, repetitions + 1):
                measured = measure_command(command, PROJECT_ROOT)
                rows.append({
                    "benchmark": "validation_tip_scaling",
                    "software": "NWKIT",
                    "shape": shape,
                    "num_tips": num_tips,
                    "num_trees": 1,
                    "replicate": replicate,
                    **measured,
                })
    return rows


def write_tree_collection(path: Path, pool: list[str], num_trees: int) -> None:
    path.write_text("\n".join(pool[index % len(pool)] for index in range(num_trees)) + "\n")


def benchmark_consensus(
    directory: Path,
    repetitions: int,
    gotree: Path | None,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    rng = random.Random(20260716)
    labels = [f"T{index}" for index in range(256)]
    pool = [random_rooted_tree(labels, rng) for _ in range(7)]
    for num_trees in (10, 50, 250):
        input_path = directory / f"consensus-{num_trees}.nwk"
        write_tree_collection(input_path, pool, num_trees)
        commands: list[tuple[str, list[str]]] = [
            (
                "NWKIT",
                [
                    "nwkit",
                    "consensus",
                    "-i",
                    str(input_path),
                    "-o",
                    os.devnull,
                    "--method",
                    "majority",
                    "--branch-length",
                    "none",
                ],
            ),
            ("PhyKIT", ["phykit", "consensus_tree", "-t", str(input_path), "-m", "majority"]),
        ]
        if gotree is not None:
            commands.append(
                (
                    "Gotree",
                    [str(gotree), "compute", "consensus", "-i", str(input_path), "-o", os.devnull, "-f", "0.5"],
                )
            )
        for software, command in commands:
            for replicate in range(1, repetitions + 1):
                measured = measure_command(command, PROJECT_ROOT)
                rows.append({
                    "benchmark": "consensus_tree_count_scaling",
                    "software": software,
                    "shape": "mixed_binary",
                    "num_tips": len(labels),
                    "num_trees": num_trees,
                    "replicate": replicate,
                    **measured,
                })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--repetitions", type=int, default=3)
    parser.add_argument("--gotree", type=Path)
    args = parser.parse_args()
    args.results_dir.mkdir(parents=True, exist_ok=True)
    hardware_metadata().to_csv(args.results_dir / "hardware.tsv", sep="\t", index=False)

    with tempfile.TemporaryDirectory(prefix="nwkit-benchmark-") as temporary_directory:
        directory = Path(temporary_directory)
        rows = benchmark_validation(directory, args.repetitions)
        rows.extend(benchmark_consensus(directory, args.repetitions, args.gotree))
    results = pd.DataFrame(rows)
    results.to_csv(args.results_dir / "benchmark.tsv", sep="\t", index=False)
    failures = results.loc[results["returncode"] != 0]
    if not failures.empty:
        print(failures.to_string(index=False))


if __name__ == "__main__":
    main()
