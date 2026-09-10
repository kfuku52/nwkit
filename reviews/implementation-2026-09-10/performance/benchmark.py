"""Paired fresh-process CLI wall/CPU time and peak RSS, with full TSV comparison.

Run --manifest MANIFEST --output RESULTS; each manifest case declares two
checkouts, arguments and output paths relative to its disposable run directory.
One warmup precedes three measured pairs, alternating order. No tracemalloc or
coverage instrumentation is used. Timeout and output differences remain visible.
"""

import argparse
import contextlib
import hashlib
import json
import os
import platform
import resource
import statistics
import subprocess
import sys
import time
from pathlib import Path


def worker(options):
    checkout = Path(options.checkout).resolve()
    sys.path.insert(0, str(checkout))
    command = json.loads(options.command)
    started = time.perf_counter()
    cpu_started = time.process_time()
    with open("stdout.txt", "w", encoding="utf-8") as stdout:
        with contextlib.redirect_stdout(stdout):
            from nwkit.cli import main

            try:
                main(command)
            except SystemExit as error:
                if error.code not in (None, 0):
                    raise
    elapsed, cpu = time.perf_counter() - started, time.process_time() - cpu_started
    import nwkit

    assert Path(nwkit.__file__).resolve().parent.parent == checkout
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print(
        json.dumps(
            dict(
                wall_seconds=elapsed,
                cpu_seconds=cpu,
                peak_rss_bytes=rss if sys.platform == "darwin" else rss * 1024,
                version=nwkit.__version__,
                module=str(nwkit.__file__),
            )
        )
    )


def compare_tables(first, second):
    import numpy as np
    import pandas as pd

    a, b = pd.read_csv(first, sep="\t"), pd.read_csv(second, sep="\t")
    result = dict(
        rows=[len(a), len(b)],
        added_columns=sorted(set(b) - set(a)),
        removed_columns=sorted(set(a) - set(b)),
        changed_columns=[],
    )
    if len(a) != len(b):
        result["equivalent_shared_columns"] = False
        return result
    for column in a.columns.intersection(b.columns):
        x, y = a[column], b[column]
        if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
            same = np.allclose(x, y, rtol=1e-5, atol=1e-7, equal_nan=True)
        else:
            same = x.fillna("<NA>").astype(str).equals(y.fillna("<NA>").astype(str))
        if not same:
            result["changed_columns"].append(column)
    result["equivalent_shared_columns"] = not result["changed_columns"]
    return result


def execute(case, label, repeat, directory, timeout):
    destination = directory / case["name"] / str(repeat) / label
    destination.mkdir(parents=True, exist_ok=True)
    environment = os.environ.copy()
    environment.update(
        OPENBLAS_NUM_THREADS="1",
        OMP_NUM_THREADS="1",
        VECLIB_MAXIMUM_THREADS="1",
        MKL_NUM_THREADS="1",
        NUMEXPR_NUM_THREADS="1",
        PYTHONHASHSEED="0",
        MPLBACKEND="Agg",
    )
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--worker",
        "--checkout",
        case[label],
        "--command",
        json.dumps(case["args"]),
    ]
    start = time.perf_counter()
    load = os.getloadavg()
    try:
        completed = subprocess.run(
            command,
            cwd=destination,
            env=environment,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        (destination / "stderr.txt").write_text(completed.stderr)
        record = (
            json.loads(completed.stdout)
            if completed.returncode == 0
            else dict(error=completed.stderr[-4000:])
        )
        record["exit_code"] = completed.returncode
    except subprocess.TimeoutExpired:
        record = dict(error="timeout", exit_code=None)
    record.update(
        label=label,
        repeat=repeat,
        process_wall_seconds=time.perf_counter() - start,
        load_average=load,
        directory=str(destination),
    )
    return record


def summarize(case, records):
    summary = {}
    for label in ["before", "after"]:
        rows = [r for r in records if r["label"] == label and r["repeat"] >= 0]
        successful = [r for r in rows if r["exit_code"] == 0]
        summary[label] = dict(successful=len(successful), attempted=len(rows))
        if successful:
            for key in [
                "wall_seconds",
                "cpu_seconds",
                "peak_rss_bytes",
                "process_wall_seconds",
            ]:
                values = [r[key] for r in successful]
                summary[label][key] = dict(
                    median=statistics.median(values), min=min(values), max=max(values)
                )
    good = {
        label: next(
            (
                r
                for r in records
                if r["label"] == label and r["repeat"] >= 0 and r["exit_code"] == 0
            ),
            None,
        )
        for label in ["before", "after"]
    }
    if all(good.values()):
        summary["outputs"] = {
            name: compare_tables(
                Path(good["before"]["directory"]) / name,
                Path(good["after"]["directory"]) / name,
            )
            for name in case["tables"]
        }
        summary["ratios"] = {
            key: summary["after"][key]["median"] / summary["before"][key]["median"]
            for key in ["wall_seconds", "cpu_seconds", "peak_rss_bytes"]
        }
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--checkout")
    parser.add_argument("--command")
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--only", nargs="*")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--timeout", type=int, default=900)
    args = parser.parse_args()
    if args.worker:
        worker(args)
        return
    manifest = json.loads(args.manifest.read_text())
    directory = Path(manifest["runs"])
    results = dict(
        environment=dict(
            python=sys.version,
            executable=sys.executable,
            platform=platform.platform(),
            cpu_count=os.cpu_count(),
            thread_limit=1,
            warmups=1,
            repeats=args.repeats,
        ),
        manifest=manifest,
        cases={},
    )
    for case in manifest["cases"]:
        if args.only and case["name"] not in args.only:
            continue
        records = []
        for repeat in range(-1, args.repeats):
            for label in ["before", "after"] if repeat % 2 else ["after", "before"]:
                record = execute(case, label, repeat, directory, args.timeout)
                records.append(record)
                print(
                    case["name"],
                    repeat,
                    label,
                    record.get("wall_seconds", record.get("error")),
                    flush=True,
                )
            results["cases"][case["name"]] = dict(
                records=records, summary=summarize(case, records)
            )
            args.output.write_text(json.dumps(results, indent=2) + "\n")
    for case in results["cases"].values():
        for record in case["records"]:
            if record["exit_code"] == 0:
                path = Path(record["directory"]) / "stdout.txt"
                record["stdout_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    args.output.write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
