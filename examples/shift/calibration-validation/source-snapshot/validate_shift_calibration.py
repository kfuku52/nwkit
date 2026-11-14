#!/usr/bin/env python3
"""Frozen independent validation of production complete-search calibration."""

import argparse
import gzip
import hashlib
import json
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shift_alpha_design import cases, generate_case, protocol  # noqa: E402
from shift_simulation_cases import rate  # noqa: E402

from nwkit.shift_calibration import CalibratedSearch  # noqa: E402
from nwkit.shift_candidates import tip_groups  # noqa: E402
from nwkit.util import read_tree  # noqa: E402

ENGINES = {}


def execute(case):
    newick, truth = generate_case(case)
    key = (newick, case["standard_error"])
    if key not in ENGINES:
        ENGINES[key] = CalibratedSearch(
            read_tree(newick, "auto", True, quiet=True),
            variances=np.full(case["tips"], case["standard_error"] ** 2),
        )
    search = ENGINES[key]
    row = {"case": case, "tree": newick, "truth": truth}
    try:
        fit = search.fit(
            truth["observations"], seed=case["seed"] + 9000000000, replicates=199
        )
        model = fit["model"]
        groups = [
            list(group)
            for group in tip_groups(
                search.tree, model["shift_branch_ids"], model["groups"]
            )
        ]
        row.update(
            status="completed",
            fit=fit,
            any_shift=bool(model["shift_branch_ids"]),
            partition_recovered=sorted(groups) == truth["shared_partition"],
            mean_rmse=float(
                np.sqrt(
                    np.mean(
                        (
                            np.array(fit["predicted"])
                            - np.array(list(truth["tip_mean"].values()))
                        )
                        ** 2
                    )
                )
            ),
        )
    except (ValueError, np.linalg.LinAlgError) as exc:
        row.update(status="failed", error=str(exc))
    return row


def bounded_results(executor, design):
    all_cases = list(cases(design))
    for start in range(0, len(all_cases), 16):
        yield from executor.map(execute, all_cases[start : start + 16])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=50)
    parser.add_argument("--extension-replicates", type=int, default=25)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    design = protocol(args.replicates, args.extension_replicates, seed=20260912)
    for key in (
        "floors",
        "criteria",
        "methods",
        "alpha_upper_height",
        "alpha_start_height",
        "boundary_relative_tolerance",
        "comparison_score_tolerance",
    ):
        design.pop(key)
    design.update(
        method="complete_search_parametric_bootstrap",
        calibration_replicates=199,
        calibration_level=0.05,
        source_sha256={
            path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
            for path in (
                "nwkit/shift_calibration.py",
                "nwkit/shift_candidates.py",
                "tools/validate_shift_calibration.py",
                "tools/shift_simulation_cases.py",
                "tools/shift_alpha_design.py",
            )
        },
    )
    (args.output / "protocol.json").write_text(json.dumps(design, indent=2) + "\n")
    rows = []
    with (
        ProcessPoolExecutor(max_workers=args.workers) as executor,
        gzip.open(args.output / "records.jsonl.gz", "wt") as stream,
    ):
        for index, row in enumerate(bounded_results(executor, design)):
            rows.append(row)
            stream.write(json.dumps(row, allow_nan=False) + "\n")
            stream.flush()
            if index % 25 == 0:
                print(
                    f"{index + 1} completed: {row['case']['family']} / {row['case']['scenario']}",
                    flush=True,
                )
    groups = defaultdict(list)
    for row in rows:
        c = row["case"]
        groups[c["family"], c["scenario"], c["root_model"]].append(row)
    summary = []
    for (family, scenario, root), members in groups.items():
        good = [r for r in members if r["status"] == "completed"]
        summary.append(
            dict(
                family=family,
                scenario=scenario,
                root_model=root,
                attempted=len(members),
                completed=len(good),
                any_shift=rate(sum(r["any_shift"] for r in good), len(good)),
                partition_recovered=rate(
                    sum(r["partition_recovered"] for r in good), len(good)
                ),
                alpha_limit_supported=rate(
                    sum(r["fit"]["alpha_limit_supported"] for r in good), len(good)
                ),
            )
        )
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    lines = [
        "# Independent calibration validation",
        "",
        "Seed 20260912; B=199; stagewise nominal level 0.05. Plug-in tests and finite nuisance grids do not guarantee uniform 5% error.",
        "",
        "| Family | Truth | Root | Completed | Any shift | True partition |",
        "|---|---|---|---:|---:|---:|",
    ]
    for row in summary:
        lines.append(
            f"| {row['family']} | {row['scenario']} | {row['root_model']} | {row['completed']}/{row['attempted']} | {row['any_shift']['count']}/{row['completed']} | {row['partition_recovered']['count']}/{row['completed']} |"
        )
    (args.output / "README.md").write_text("\n".join(lines) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
