"""Prespecified small held-out joint-search diagnostic grid, not error calibration."""

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from shift_simulation_cases import SCENARIOS, simulate


def execute(case, output, criterion, root_model, rscript):
    command = [
        sys.executable,
        str(Path(__file__).with_name("validate_shift_joint.py")),
        "--case",
        str(case),
        "--output",
        str(output),
        "--criterion",
        criterion,
        "--root-model",
        root_model,
        "--rscript",
        rscript,
    ]
    with output.with_suffix(".log").open("w") as log:
        completed = subprocess.run(
            command,
            stdout=log,
            stderr=subprocess.STDOUT,
            env={**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"},
            check=False,
        )
    row = {
        "case": case.name,
        "criterion": criterion,
        "status": "completed" if completed.returncode == 0 else "failed",
    }
    if completed.returncode == 0:
        row.update(json.loads((output / "summary.json").read_text()))
    else:
        row["error_log"] = str(output.with_suffix(".log"))
    return row


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, required=True)
    cli.add_argument("--rscript", default="Rscript")
    cli.add_argument("--seed", type=int, default=20260910)
    options = cli.parse_args()
    root = options.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    grid = [
        (8, scenario, model, 0.0)
        for scenario in SCENARIOS
        for model in ("OUfixedRoot", "OUrandomRoot")
    ]
    grid += [(8, "convergent", model, 0.2) for model in ("OUfixedRoot", "OUrandomRoot")]
    grid += [
        (16, "convergent", model, 0.0) for model in ("OUfixedRoot", "OUrandomRoot")
    ]
    sources = [
        Path(__file__),
        *[
            Path(__file__).with_name(name)
            for name in (
                "validate_shift_joint.py",
                "shift_joint_candidates.py",
                "shift_joint_backend.py",
                "shift_simulation_cases.py",
            )
        ],
    ]
    source_root = Path(__file__).resolve().parents[1]
    sources += sorted((source_root / "nwkit").glob("shift*.py"))
    hashes = {}
    for source in sources:
        relative = source.relative_to(source_root)
        target = root / "source" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        hashes[str(relative)] = hashlib.sha256(source.read_bytes()).hexdigest()
    (root / "manifest.json").write_text(
        json.dumps(
            {
                "seed": options.seed,
                "grid": grid,
                "effect": 2,
                "criteria": "BIC for all cases; pBIC also for eight-tip convergent/no-SE cases",
                "source_sha256": hashes,
            },
            indent=2,
        )
        + "\n"
    )
    jobs = []
    for index, (tips, scenario, model, se) in enumerate(grid):
        case = root / f"case-{index:02d}"
        case.mkdir()
        seed = int(np.random.SeedSequence([options.seed, index]).generate_state(1)[0])
        newick, truth = simulate(
            tips=tips, scenario=scenario, root_model=model, se=se, effect=2, seed=seed
        )
        (case / "tree.nwk").write_text(newick + "\n")
        (case / "truth.json").write_text(json.dumps(truth, indent=2) + "\n")
        pd.DataFrame(
            {"leaf_name": truth["tip_names"], "value": truth["observations"], "se": se}
        ).to_csv(case / "traits.tsv", sep="\t", index=False)
        criteria = (
            ["BIC", "pBIC"]
            if tips == 8 and scenario == "convergent" and se == 0
            else ["BIC"]
        )
        for criterion in criteria:
            jobs.append(
                (
                    case,
                    root / f"fit-{index:02d}-{criterion}",
                    criterion,
                    model,
                    options.rscript,
                )
            )
    rows = []
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(execute, *job) for job in jobs]
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            with (root / "results.jsonl").open("a") as stream:
                stream.write(json.dumps(row, allow_nan=False) + "\n")
            print(row["case"], row["criterion"], row["status"], flush=True)
    rows.sort(key=lambda row: (row["case"], row["criterion"]))
    (root / "results.json").write_text(
        json.dumps(rows, indent=2, allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    main()
