"""Frozen bounded paired experiment for confidence-set restricted calibration."""

import argparse
import gzip
import hashlib
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.stats import beta as beta_distribution

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from shift_alpha_design import generate_case  # noqa: E402
from shift_confidence_calibration import NullConfidenceBank  # noqa: E402

from nwkit.shift_calibration import CalibratedSearch  # noqa: E402
from nwkit.util import read_tree  # noqa: E402

ENGINE = None
SCENARIOS = ("null", "single", "distinct", "convergent")
BETAS = (0.005, 0.01)


def generating_case(seed, scenario):
    return dict(
        tips=8,
        scenario=scenario,
        root_model="OUfixedRoot",
        alpha_height=2.1,
        sigma2_height=0.25,
        standard_error=0.0,
        effect=2.0,
        seed=seed,
    )


def execute(case):
    global ENGINE
    index, specification = case
    data_seed, bootstrap_seed = np.random.SeedSequence(
        [specification["master_seed"], index]
    ).spawn(2)
    data_seed = int(data_seed.generate_state(1)[0])
    bootstrap_seed = int(bootstrap_seed.generate_state(1)[0])
    started = time.perf_counter()
    row = dict(block=index, data_seed=data_seed, bootstrap_seed=bootstrap_seed)
    try:
        text, _ = generate_case(generating_case(data_seed, "null"))
        if ENGINE is None:
            ENGINE = CalibratedSearch(
                read_tree(text, "auto", True, quiet=True), convergence=False
            )
        bank = NullConfidenceBank(ENGINE)
        row["bank"] = bank.simulate(
            bootstrap_seed, specification["bootstrap_replicates"]
        )
        row["tree"] = text
        row["conditions"] = []
        for scenario in SCENARIOS:
            _, truth = generate_case(generating_case(data_seed, scenario))
            row["conditions"].append(
                dict(
                    scenario=scenario,
                    truth=truth,
                    result=bank.evaluate(truth["observations"], BETAS),
                )
            )
        row["status"] = "completed"
    except (ValueError, np.linalg.LinAlgError) as exc:
        row.update(status="failed", error=str(exc))
    row["seconds"] = time.perf_counter() - started
    return row


def summarize(rows, phase):
    completed = [r for r in rows if r["status"] == "completed"]
    failures = len(rows) - len(completed)
    cells = []
    for scenario in SCENARIOS:
        results = [
            next(c["result"] for c in r["conditions"] if c["scenario"] == scenario)
            for r in completed
        ]
        for b in BETAS:
            candidates = [
                next(c for c in r["candidates"] if c["beta"] == b) for r in results
            ]
            baseline = [r["full_envelope_reject"] for r in results]
            chosen = [c["reject"] for c in candidates]
            gains = sum(c and not f for c, f in zip(chosen, baseline, strict=True))
            losses = sum(f and not c for c, f in zip(chosen, baseline, strict=True))
            worst_gains = gains + failures
            upper = (
                1.0
                if worst_gains == len(rows)
                else float(
                    beta_distribution.ppf(
                        0.975,
                        worst_gains + 1,
                        len(rows) - worst_gains,
                    )
                )
            )
            cells.append(
                dict(
                    scenario=scenario,
                    beta=b,
                    attempted=len(rows),
                    completed=len(completed),
                    failed=failures,
                    baseline_rejections=sum(baseline),
                    candidate_rejections=sum(chosen),
                    paired_gains=gains,
                    paired_losses=losses,
                    net_gain=(gains - losses) / len(rows),
                    gross_gain_simultaneous_95_upper=upper,
                    mean_retained_grid_points=float(
                        np.mean([c["retained_count"] for c in candidates])
                    )
                    if candidates
                    else None,
                    empty_confidence_sets=sum(
                        c["empty_confidence_set"] for c in candidates
                    ),
                )
            )
    advance = []
    for b in BETAS:
        selected = {r["scenario"]: r for r in cells if r["beta"] == b}
        if (
            phase == "development"
            and not failures
            and selected["distinct"]["net_gain"] >= 0.05
            and all(selected[s]["net_gain"] >= -0.05 for s in ("single", "convergent"))
            and selected["null"]["candidate_rejections"] / len(rows) <= 0.075
        ):
            advance.append(b)
    return dict(
        phase=phase,
        independent_blocks=len(rows),
        completed_blocks=len(completed),
        failed_blocks=failures,
        paired_observations=4 * len(rows),
        cells=cells,
        candidates_passing_screen=advance,
        production_adoption=False,
        interpretation="Development screening only. Stop if no candidate passes; positive screening still requires new independent validation.",
        gain_bound_scope="Bonferroni over the two beta candidates for the prespecified distinct-shift futility assessment; no joint claim across all scenarios",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--phase", choices=("pilot", "development"), default="development"
    )
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("Positive workers required")
    sources = [
        "nwkit/shift_calibration.py",
        "nwkit/shift_candidates.py",
        "tools/shift_confidence_calibration.py",
        "tools/try_shift_confidence_calibration.py",
        "tools/shift_alpha_design.py",
        "tools/shift_simulation_cases.py",
        "reviews/shift-calibration-final-protocol.md",
    ]
    specification = dict(
        phase=args.phase,
        blocks=2 if args.phase == "pilot" else 100,
        master_seed=20260928 if args.phase == "pilot" else 20260930,
        bootstrap_replicates=999,
        betas=BETAS,
        level=0.05,
        convergence=False,
        candidate_grid="unchanged production 27-point alpha-height grid",
        workers=args.workers,
        source_sha256={
            s: hashlib.sha256(Path(s).read_bytes()).hexdigest() for s in sources
        },
    )
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "protocol.json").write_text(
        json.dumps(specification, indent=2) + "\n"
    )
    snapshot = args.output / "source-snapshot"
    snapshot.mkdir()
    for s in sources:
        (snapshot / Path(s).name).write_bytes(Path(s).read_bytes())
    rows = []
    with (
        ProcessPoolExecutor(max_workers=args.workers) as executor,
        gzip.open(args.output / "records.jsonl.gz", "wt") as stream,
    ):
        cases = [(i, specification) for i in range(specification["blocks"])]
        for row in executor.map(execute, cases):
            rows.append(row)
            stream.write(json.dumps(row, allow_nan=False) + "\n")
            stream.flush()
            print(
                f"Completed {len(rows)}/{len(cases)} paired blocks ({row['status']})",
                flush=True,
            )
    result = summarize(rows, args.phase)
    (args.output / "summary.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps({k: v for k, v in result.items() if k != "cells"}), flush=True)


if __name__ == "__main__":
    main()
