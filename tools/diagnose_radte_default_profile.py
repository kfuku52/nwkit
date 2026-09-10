"""Replay the first N families with generating rate SD supplied (oracle only).

These are paired diagnostic reanalyses, not new independent validation families.
No successful-fit selection, model tuning, or change to the public estimator.
"""

import argparse
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from validate_radte_default_profile import (
    digest,
    finite_number,
    read_target,
    save_json,
    summarize,
)


def replay(source, output, timeout):
    original = json.loads((source / "result.json").read_text())
    truth = json.loads((source / "truth.json").read_text())
    row = dict(
        case=original["case"],
        family=original["family"],
        truth=truth["age"],
        point_success=False,
        interval_available=False,
        covered=False,
        original_result_sha256=digest(source / "result.json"),
        supplied_rate_sd=truth["case"]["rate_sd"],
    )
    output.mkdir(parents=True, exist_ok=False)
    command_path = source / "profile.command.json"
    if not command_path.exists():
        command_path = source / "none.command.json"
    if not command_path.exists():
        row["status"] = "unavailable-source-inference-command"
        save_json(output / "result.json", row)
        return row
    command = json.loads(command_path.read_text())
    command[0] = sys.executable
    command[command.index("--uncertainty") + 1] = "profile"
    command[command.index("--out-prefix") + 1] = str(output / "oracle")
    command.extend(["--rate-sd", str(row["supplied_rate_sd"])])
    save_json(output / "command.json", command)
    with (output / "run.log").open("w") as log:
        try:
            process = subprocess.run(
                command,
                cwd=source,
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=timeout,
            )
            row["status"] = (
                "complete" if process.returncode == 0 else f"exit-{process.returncode}"
            )
        except subprocess.TimeoutExpired:
            row["status"] = "timeout"
    if row["status"] == "complete":
        manifest, target = read_target(output / "oracle")
        if target is None:
            row["status"] = "target-unmatched"
        else:
            row["profile_estimator"] = manifest["method"]
            row["diagnostics"] = manifest["diagnostics"]
            row["estimate"] = finite_number(target["estimated_age"])
            row["point_success"] = row["estimate"] is not None
            if not row["point_success"]:
                row["status"] = "nonfinite-estimate"
                save_json(output / "result.json", row)
                return row
            row["bias"] = row["estimate"] - row["truth"]
            lower, upper = (
                finite_number(target[k]) for k in ["interval_lower", "interval_upper"]
            )
            row.update(
                lower=lower, upper=upper, interval_status=target["interval_status"]
            )
            row["interval_available"] = (
                lower is not None and upper is not None and lower <= upper
            )
            row["covered"] = (
                row["interval_available"] and lower <= row["truth"] <= upper
            )
            row["width"] = upper - lower if row["interval_available"] else None
    save_json(output / "result.json", row)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--case", required=True)
    parser.add_argument("--families", type=int, default=50)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--timeout", type=float, default=180)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.input, args.output = args.input.resolve(), args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=False)
    # Preserve the exact metadata bytes even if the parent study subsequently
    # appends its end-of-run source-integrity checks.
    (args.output / "source-metadata.json").write_bytes(
        (args.input / "metadata.json").read_bytes()
    )
    sources = [args.input / args.case / f"f{i:04d}" for i in range(args.families)]
    originals = [json.loads((p / "result.json").read_text()) for p in sources]
    save_json(
        args.output / "protocol.json",
        dict(
            role="diagnostic-oracle",
            selection="First N family indices, irrespective of outcomes",
            arguments={
                k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()
            },
            script_sha256=digest(__file__),
            original_metadata_sha256=digest(args.output / "source-metadata.json"),
        ),
    )
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = [
            pool.submit(replay, source, args.output / source.name, args.timeout)
            for source in sources
        ]
        rows = []
        for future in as_completed(futures):
            rows.append(future.result())
            print(f"{len(rows)}/{len(futures)}", flush=True)
    save_json(
        args.output / "summary.json",
        dict(original=summarize(originals), known_sd=summarize(rows)),
    )
    save_json(args.output / "cases.json", rows)


if __name__ == "__main__":
    main()
