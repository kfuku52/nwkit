"""Independent finite-null-grid calibration study with all failures retained."""

import argparse
import hashlib
import json
import logging
import platform
import sys
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from pathlib import Path

import Bio
import msprime
import numpy as np
import pandas as pd
import scipy
from pilot import PARSER, SPECIES, model, reference_family, tree

from nwkit.mul_locus import LocusParameters
from nwkit.mul_locus_cli import bank_record, finite_json
from nwkit.mul_locus_mc import (
    build_bank,
    calibrate,
    integration_alpha,
    make_tasks,
    probability_interval,
    validate_model,
)
from nwkit.mul_msc import dated_text
from nwkit.mul_msc_fit import species_topology_signature

SCENARIOS = {
    "baseline": ((0.03, 0.08), 0.03, (0.3, 1.0), 0.9),
    "turnover": ((0.08, 0.16), 0.12, (0.3, 1.0), 0.9),
    "missing-ils": ((0.03, 0.08), 0.03, (1.0, 3.0), 0.6),
}
_BANKS = None
_MODEL = None


def scenario_model(name, samples, seed):
    duplication, loss, ne, detection = SCENARIOS[name]
    config = model(samples)
    config.update(
        integration="ancestral-stratified",
        seed=seed,
        detection={s: detection for s in ("A", "X", "B")},
        parameter_grid=[
            {"duplication": d, "loss": loss, "ne": n, "hybridization_age": a}
            for d in duplication
            for n in ne
            for a in (0.3, 0.7)
        ],
    )
    return config


def scenario_cases(name):
    duplication, loss, ne, _ = SCENARIOS[name]
    cases = [
        (f"null-d{i}-ne{j}", 0, LocusParameters(d, loss, n, 0.5), True)
        for i, d in enumerate(duplication)
        for j, n in enumerate(ne)
    ]
    midpoint = sum(duplication) / 2
    cases.append(
        ("null-off-grid", 0, LocusParameters(midpoint, loss, max(ne), 0.5), False)
    )
    cases.extend(
        (
            f"allop-{parent}",
            candidate,
            LocusParameters(midpoint, loss, max(ne), 0.5),
            False,
        )
        for parent, candidate in (("A", 1), ("B", 2))
    )
    return cases


def worker_init(banks, config):
    global _BANKS, _MODEL
    _BANKS, _MODEL = banks, config
    logging.getLogger("msprime").setLevel(logging.WARNING)


def integration_bank(task, config):
    try:
        return build_bank(task, config)
    except Exception as error:
        return {
            "candidate": task[0],
            "grid": task[2],
            "status": "bank-failed",
            "error": f"{type(error).__name__}: {error}",
        }


def calibration_seed(args, scenario_index, case_index, replicate):
    return int(
        np.random.SeedSequence(
            [args.calibration_seed, scenario_index, case_index, replicate]
        ).generate_state(1)[0]
    )


def trial_jobs(args, name, index):
    return [
        (name, case, truth, point, grid, index, i, r, args)
        for i, (case, truth, point, grid) in enumerate(scenario_cases(name))
        for r in range(args.replicates)
    ]


def trial_record(job):
    name, case, truth, point, on_grid, index, case_index, replicate, args = job
    return {
        "scenario": name,
        "case": case,
        "replicate": replicate + 1,
        "families": args.families,
        "true_candidate": truth,
        "null_on_grid": truth == 0 and on_grid,
        "true_parameters": json.dumps(asdict(point), sort_keys=True),
        "calibration_seed": calibration_seed(args, index, case_index, replicate),
    }


def validate_calibration_seeds(args):
    seen = {}
    for name in args.scenarios:
        index = tuple(SCENARIOS).index(name)
        for job in trial_jobs(args, name, index):
            row = trial_record(job)
            identity = name, row["case"], row["replicate"]
            seed = row["calibration_seed"]
            if seed in seen:
                raise ValueError(
                    f"Calibration seed collision between {seen[seed]} and {identity}: {seed}."
                )
            seen[seed] = identity


def failed_trial(job, status, error):
    row = trial_record(job)
    directory = job[-1].output / row["scenario"] / f"{row['case']}-r{row['replicate']}"
    directory.mkdir(exist_ok=True)
    if status == "worker-failed":
        for name in ("failure.json", "summary.json"):
            path = directory / name
            if path.exists():
                try:
                    saved = json.loads(path.read_text())
                except (OSError, ValueError):
                    continue
                if (
                    isinstance(saved, dict)
                    and all(saved.get(key) == value for key, value in row.items())
                    and saved.get("status")
                    in ("completed", "generation-failed", "calibration-failed")
                    and (
                        saved["status"] != "completed"
                        or all(
                            type(saved.get(key)) is bool
                            for key in (
                                "reported_event",
                                "point_event",
                                "fitted_null_event_same_draws",
                                "parent_correct",
                            )
                        )
                    )
                ):
                    return saved
    row.update(status=status, error=error)
    name = (
        "worker-failure.json"
        if (directory / "failure.json").exists()
        else "failure.json"
    )
    (directory / name).write_text(json.dumps(row, indent=2) + "\n")
    return row


def evaluate(job):
    name, case, truth, point, on_grid, scenario_index, case_index, replicate, args = job
    directory = args.output / name / f"{case}-r{replicate + 1}"
    directory.mkdir()
    row = trial_record(job)
    phase = "generation"
    try:
        population = next(
            t[3]
            for t in make_tasks(tree(SPECIES), "X", "A B", [point], 100)[0]
            if t[0] == truth
        )
        rng = np.random.default_rng(
            np.random.SeedSequence([args.seed, scenario_index, case_index, replicate])
        )
        families = [
            reference_family(population, point, _MODEL, rng)
            for _ in range(args.families)
        ]
        (directory / "families.jsonl").write_text(
            "".join(json.dumps(f[3]) + "\n" for f in families)
        )
        for view, index in (("true", 1), ("estimated", 2)):
            (directory / f"{view}.nwk").write_text(
                "".join(dated_text(f[index]) + "\n" for f in families)
            )
        phase = "calibration"
        observations = [species_topology_signature(f[1], PARSER) for f in families]
        config = {**_MODEL, "seed": row["calibration_seed"]}
        fitted, calibration = calibrate(
            _BANKS,
            observations,
            config,
            integration_alpha(config, len(_BANKS)),
            args.bootstrap,
            null_calibration="grid-supremum",
        )
        (directory / "results.json").write_text(
            json.dumps(
                finite_json(
                    {
                        "observations": observations,
                        "fit": fitted,
                        "calibration": calibration,
                    }
                ),
                indent=2,
                allow_nan=False,
            )
            + "\n"
        )
        plug_in = next(
            s
            for s in calibration["null_grid_calibrations"]
            if s["generating_null_grid"] == fitted["null"]["grid"]
        )
        row.update(
            status="completed",
            selected_alternative=fitted["alternative"]["mul.tree"],
            parent_correct=truth != 0 and fitted["alternative"]["mul.tree"] == truth,
            p_value=calibration["p_value"],
            mc_p_lower=calibration["mc_p_lower"],
            mc_p_upper=calibration["mc_p_upper"],
            reported_event=calibration["mc_p_upper"] <= args.event_alpha,
            point_event=calibration["p_value"] <= args.event_alpha,
            fitted_null_p_value_same_draws=plug_in["p_value"],
            fitted_null_mc_p_upper_same_draws=plug_in["mc_p_upper"],
            fitted_null_event_same_draws=plug_in["mc_p_upper"] <= args.event_alpha,
            least_favorable_grid=calibration["least_favorable_grid"],
            least_favorable_mc_upper_grid=calibration["least_favorable_mc_upper_grid"],
            tree_accuracy=sum(
                species_topology_signature(f[1], PARSER)
                == species_topology_signature(f[2], PARSER)
                for f in families
            )
            / len(families),
        )
        (directory / "summary.json").write_text(json.dumps(row, indent=2) + "\n")
    except Exception as error:
        row.update(status=f"{phase}-failed", error=f"{type(error).__name__}: {error}")
        (directory / "failure.json").write_text(json.dumps(row, indent=2) + "\n")
    return row


def summarize(rows):
    groups = {}
    for row in rows:
        key = row["scenario"], row["case"], row["true_candidate"], row["null_on_grid"]
        groups.setdefault(key, []).append(row)
    results = []
    for (scenario, case, truth, on_grid), group in groups.items():
        completed = [r for r in group if r["status"] == "completed"]
        failures = len(group) - len(completed)
        reports = sum(r["reported_event"] for r in completed)
        low = probability_interval(reports, len(group), 0.05)[0]
        high = probability_interval(reports + failures, len(group), 0.05)[1]
        results.append(
            {
                "scenario": scenario,
                "case": case,
                "true_candidate": truth,
                "null_on_grid": on_grid,
                "planned": len(group),
                "completed": len(completed),
                "failed": failures,
                "reported": reports,
                "point_reported": sum(r["point_event"] for r in completed),
                "fitted_null_reported_same_draws": sum(
                    r["fitted_null_event_same_draws"] for r in completed
                ),
                "parent_correct": sum(r["parent_correct"] for r in completed),
                "operational_report_rate": reports / len(group),
                "failure_aware_rate_lower": reports / len(group),
                "failure_aware_rate_upper": (reports + failures) / len(group),
                "failure_aware_95_lower": low,
                "failure_aware_95_upper": high,
            }
        )
    return results


def write_protocol(args):
    validate_calibration_seeds(args)
    source = Path(__file__).resolve().parents[2]
    files = [
        Path(__file__),
        Path(__file__).with_name("pilot.py"),
        *(
            source / "nwkit" / name
            for name in (
                "mul_locus.py",
                "mul_locus_mc.py",
                "mul_locus_cli.py",
                "mul_msc_model.py",
                "mul_coalescent.py",
                "mul_msc_fit.py",
                "species_parser.py",
                "util.py",
            )
        ),
    ]
    protocol = {
        "arguments": {
            k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()
        },
        "species_tree": SPECIES,
        "h1": "X",
        "h2": "A B",
        "scenario_indices": {
            name: tuple(SCENARIOS).index(name) for name in args.scenarios
        },
        "models": {
            name: scenario_model(
                name, args.samples, args.bank_seed + tuple(SCENARIOS).index(name)
            )
            for name in args.scenarios
        },
        "cases": {
            name: [
                {
                    "name": n,
                    "candidate": t,
                    "parameters": asdict(p),
                    "null_on_grid": t == 0 and g,
                }
                for n, t, p, g in scenario_cases(name)
            ]
            for name in args.scenarios
        },
        "python": sys.version,
        "platform": platform.platform(),
        "versions": {
            "msprime": msprime.__version__,
            "biopython": Bio.__version__,
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "sources": {
            str(p.relative_to(source)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in files
        },
        "method": "independent-SSA-msprime-evaluation; true-rooted-topology finite-grid supremum calibration",
        "decision": "reported_event iff mc_p_upper <= event_alpha; point decisions and fitted-null comparisons separate",
        "limitations": "True-tree event calibration only. Estimated trees retained for later pipeline-aware work, not calibrated here. Off-grid nulls are misspecification diagnostics without a finite-grid validity guarantee. All failures count in planned denominators; failure-aware bounds allow every failed decision either outcome. Intervals are pointwise descriptive, conditional on shared frozen banks; not simultaneous or general biological guarantees.",
    }
    (args.output / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    return protocol


def collect_banks(tasks, config, cpus):
    banks = []
    try:
        if cpus == 1:
            banks = [integration_bank(task, config) for task in tasks]
        else:
            with ProcessPoolExecutor(max_workers=cpus) as pool:
                for bank in pool.map(integration_bank, tasks, [config] * len(tasks)):
                    banks.append(bank)
    except Exception as error:
        banks.extend(
            {
                "candidate": task[0],
                "grid": task[2],
                "status": "bank-failed",
                "error": f"{type(error).__name__}: {error}",
            }
            for task in tasks[len(banks) :]
        )
    return banks


def run_scenario(args, name, index, config):
    directory = args.output / name
    directory.mkdir()
    tasks, excluded = make_tasks(
        tree(SPECIES), "X", "A B", validate_model(config, tree(SPECIES)), 100
    )
    jobs = trial_jobs(args, name, index)
    (directory / "excluded.json").write_text(json.dumps(excluded, indent=2) + "\n")
    banks = collect_banks(tasks, config, args.cpus)
    failures = [b for b in banks if isinstance(b, dict)]
    if failures:
        (directory / "bank-failures.json").write_text(
            json.dumps(failures, indent=2) + "\n"
        )
        (directory / "partial-banks.json").write_text(
            json.dumps(
                [bank_record(b) for b in banks if not isinstance(b, dict)],
                allow_nan=False,
            )
            + "\n"
        )
        rows = [failed_trial(job, "bank-failed", json.dumps(failures)) for job in jobs]
        pd.DataFrame(rows).to_csv(directory / "summary.tsv", sep="\t", index=False)
        return rows
    (directory / "banks.json").write_text(
        json.dumps([bank_record(b) for b in banks], allow_nan=False) + "\n"
    )
    print(f"{name}: {len(banks)} independent banks completed", flush=True)
    rows = []
    if args.cpus == 1:
        worker_init(banks, config)
        iterator = map(evaluate, jobs)
        for row in iterator:
            rows.append(row)
            save_progress(directory, rows, row)
    else:
        try:
            with ProcessPoolExecutor(
                max_workers=args.cpus, initializer=worker_init, initargs=(banks, config)
            ) as pool:
                for row in pool.map(evaluate, jobs, chunksize=1):
                    rows.append(row)
                    save_progress(directory, rows, row)
        except Exception as error:
            for job in jobs[len(rows) :]:
                row = failed_trial(
                    job, "worker-failed", f"{type(error).__name__}: {error}"
                )
                rows.append(row)
                save_progress(directory, rows, row)
    return rows


def save_progress(directory, rows, row):
    pd.DataFrame(rows).to_csv(directory / "summary.tsv", sep="\t", index=False)
    print(
        json.dumps({k: row[k] for k in ("scenario", "case", "replicate", "status")}),
        flush=True,
    )


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, required=True)
    cli.add_argument(
        "--scenarios", nargs="+", choices=tuple(SCENARIOS), default=list(SCENARIOS)
    )
    for name, default in (
        ("samples", 100000),
        ("replicates", 5),
        ("families", 30),
        ("bootstrap", 99),
        ("cpus", 1),
        ("seed", 20261033),
        ("bank-seed", 20261030),
        ("calibration-seed", 20261034),
    ):
        cli.add_argument("--" + name, type=int, default=default)
    cli.add_argument("--event-alpha", type=float, default=0.05)
    args = cli.parse_args()
    if min(
        args.samples, args.replicates, args.families, args.bootstrap, args.cpus
    ) < 1 or (
        min(args.seed, args.bank_seed, args.calibration_seed) < 0
        or not 0 < args.event_alpha < 1
        or len(set(args.scenarios)) != len(args.scenarios)
    ):
        cli.error(
            "Positive budgets/CPUs, nonnegative seeds, alpha in (0,1), and distinct scenarios required."
        )
    args.output.mkdir(parents=True, exist_ok=False)
    logging.getLogger("msprime").setLevel(logging.WARNING)
    protocol = write_protocol(args)
    rows = []
    for name in args.scenarios:
        rows.extend(
            run_scenario(
                args, name, protocol["scenario_indices"][name], protocol["models"][name]
            )
        )
        pd.DataFrame(rows).to_csv(args.output / "summary.tsv", sep="\t", index=False)
        pd.DataFrame(summarize(rows)).to_csv(
            args.output / "rates.tsv", sep="\t", index=False
        )
    statuses = Counter(r["status"] for r in rows)
    print(json.dumps({"analyses": len(rows), "statuses": statuses}), flush=True)
    return int(any(r["status"] != "completed" for r in rows))


if __name__ == "__main__":
    raise SystemExit(main())
