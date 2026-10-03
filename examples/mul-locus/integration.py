"""Frozen-protocol support diagnosis and independent paired integration probe."""

import argparse
import hashlib
import json
import logging
import platform
import re
import sys
import time
from collections import Counter
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd
from calibration import SCENARIOS, scenario_cases, scenario_model
from pilot import SPECIES, tree
from pilot import reference_true_selected as reference_selected

from nwkit.mul_locus import sample_selected
from nwkit.mul_locus_cli import bank_record, finite_json
from nwkit.mul_locus_integral import (
    build_paired_banks,
    integrated_probability,
    score_integrated_bank,
)
from nwkit.mul_locus_mc import (
    LocusBank,
    calibrate,
    integration_alpha,
    make_tasks,
    pattern_probability,
    score_bank,
    validate_model,
)

METHODS = ("histogram", "detection", "hybrid")
SOURCES = (
    "examples/mul-locus/integration.py",
    "examples/mul-locus/calibration.py",
    "examples/mul-locus/pilot.py",
    "examples/mul-locus/reference.py",
    "nwkit/mul_locus_integral.py",
    "nwkit/mul_locus.py",
    "nwkit/mul_locus_mc.py",
    "nwkit/mul_coalescent.py",
    "nwkit/mul_msc_model.py",
    "nwkit/mul_msc_fit.py",
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    path.write_text(json.dumps(finite_json(value), indent=2, allow_nan=False) + "\n")


def freeze(value):
    return tuple(freeze(v) for v in value) if isinstance(value, list) else value


def copy_vector(signature):
    if signature is None:
        return ()
    if signature[0] == "tip":
        return (signature[1],)
    return tuple(sorted(s for child in signature[1:] for s in copy_vector(child)))


def read_banks(path, config):
    tasks = make_tasks(
        tree(SPECIES), "X", "A B", validate_model(config, tree(SPECIES)), 100
    )[0]
    by_key = {(t[0], t[2]): t for t in tasks}
    banks = []
    for record in json.loads(path.read_text()):
        task = by_key[(record["candidate"], record["grid"])]
        counts = Counter({freeze(p["signature"]): p["hits"] for p in record["counts"]})
        strata = record.get("strata")
        if strata is not None:
            strata = [
                {
                    **s,
                    "counts": Counter(
                        {freeze(p["signature"]): p["hits"] for p in s["counts"]}
                    ),
                }
                for s in strata
            ]
        banks.append(
            LocusBank(*task, counts, record["samples"], record["attempts"], strata)
        )
    return banks


def support_rows(banks, observations):
    rows = []
    vectors = [
        {
            vector: sum(
                hits for key, hits in bank.counts.items() if copy_vector(key) == vector
            )
            for vector in set(copy_vector(key) for key in observations)
        }
        for bank in banks
    ]
    for signature, multiplicity in Counter(observations).items():
        for bank, counts in zip(banks, vectors, strict=True):
            hits = bank.counts[signature]
            count_hits = counts[copy_vector(signature)]
            rows.append(
                {
                    "signature": signature,
                    "multiplicity": multiplicity,
                    "candidate": bank.candidate,
                    "grid": bank.grid,
                    "pattern_hits": hits,
                    "copy_vector": copy_vector(signature),
                    "copy_vector_hits": count_hits,
                    "absence": "none"
                    if hits
                    else "copy-count-or-detection"
                    if not count_hits
                    else "topology",
                }
            )
    return rows


def diagnose(source, output):
    protocol = json.loads((source / "protocol.json").read_text())
    rows = []
    for name, config in protocol["models"].items():
        banks = read_banks(source / name / "banks.json", config)
        for path in sorted((source / name).glob("*/failure.json")):
            failure = json.loads(path.read_text())
            if failure["status"] != "calibration-failed":
                raise ValueError(
                    "Support diagnosis requires a calibration failure, not a worker failure."
                )
            match = re.search(
                r"generating grid (\d+) replicate (\d+)", failure["error"]
            )
            if match is None:
                raise ValueError(
                    "Failure has no reproducible generating grid/replicate."
                )
            grid, replicate = map(int, match.groups())
            null = next(b for b in banks if b.candidate == 0 and b.grid == grid)
            rng = np.random.default_rng(
                np.random.SeedSequence(
                    [failure["calibration_seed"], 2, grid, replicate - 1]
                )
            )
            observations, attempts = sample_selected(
                null.population, null.parameters, config, rng, count=failure["families"]
            )
            patterns = support_rows(banks, observations)
            finite = {
                kind: sum(
                    all(b.counts[s] for s in observations)
                    for b in banks
                    if bool(b.candidate) == kind
                )
                for kind in (False, True)
            }
            if finite[False] and finite[True]:
                raise ArithmeticError(
                    "Recorded failure did not reproduce missing support."
                )
            row = {
                **failure,
                "dataset": path.parent.name,
                "generating_grid": grid,
                "null_replicate": replicate,
                "attempts": attempts,
                "finite_null_banks": finite[False],
                "finite_alternative_banks": finite[True],
                "observations": observations,
                "patterns": patterns,
                "source_failure_sha256": digest(path),
                "source_bank_sha256": digest(source / name / "banks.json"),
            }
            rows.append(row)
    save(output / "diagnosis.json", rows)
    return rows


def integrated_record(bank):
    return {
        "candidate": bank.candidate,
        "h2": bank.h2,
        "grid": bank.grid,
        "parameters": asdict(bank.parameters),
        "method": bank.method,
        "interval_method": bank.interval_method,
        "strata": [
            {
                **{
                    key: value
                    for key, value in stratum.items()
                    if key not in ("selected", "patterns")
                },
                "selected": asdict(stratum["selected"]),
                "patterns": [
                    {"signature": key, "moments": asdict(moment)}
                    for key, moment in stratum["patterns"].items()
                ],
            }
            for stratum in bank.strata
        ],
    }


def integration_tables(banks, observations, config):
    rows = []
    alpha = integration_alpha(config, len(banks["histogram"]))
    for method, candidates in banks.items():
        probability = (
            pattern_probability if method == "histogram" else integrated_probability
        )
        for bank in candidates:
            for signature in dict.fromkeys(observations):
                estimate, lo, hi = probability(bank, signature, alpha)
                rows.append(
                    {
                        "method": method,
                        "candidate": bank.candidate,
                        "grid": bank.grid,
                        "signature": signature,
                        "estimate": estimate,
                        "lower": lo,
                        "upper": hi,
                        "width": hi - lo,
                        "zero_estimate": estimate == 0,
                    }
                )
    return rows


def case_records(args, name, case_index, replicate, case):
    label, truth, point, on_grid = case
    index = tuple(SCENARIOS).index(name)
    words = np.random.SeedSequence(
        [args.calibration_seed, 7, index, case_index, replicate]
    ).generate_state(4)
    seed = sum(int(word) << (32 * i) for i, word in enumerate(words))
    return [
        {
            "scenario": name,
            "case": label,
            "replicate": replicate + 1,
            "method": method,
            "true_candidate": truth,
            "true_parameters": asdict(point),
            "null_on_grid": truth == 0 and on_grid,
            "families": args.families,
            "calibration_seed": str(seed),
            "data_seed_namespace": [args.seed, 6, index, case_index, replicate],
        }
        for method in METHODS
    ]


def evaluate_case(args, name, case_index, replicate, case, banks, config):
    sampler = reference_selected
    if getattr(args, "reference", "msprime") == "conditional":
        from reference import reference_selected as sampler
    label, truth, point, on_grid = case
    index = tuple(SCENARIOS).index(name)
    directory = args.output / name / f"{label}-r{replicate + 1}"
    directory.mkdir()
    population = next(
        t[3]
        for t in make_tasks(tree(SPECIES), "X", "A B", [point], 100)[0]
        if t[0] == truth
    )
    rng = np.random.default_rng(
        np.random.SeedSequence([args.seed, 6, index, case_index, replicate])
    )
    planned = case_records(args, name, case_index, replicate, case)
    try:
        observations, attempts = sampler(
            population, point, config, rng, count=args.families
        )
        if len(observations) != args.families:
            raise ValueError("Independent generator returned the wrong family count.")
    except Exception as error:
        for row in planned:
            row.update(
                status="generation-failed", error=f"{type(error).__name__}: {error}"
            )
            save(directory / f"{row['method']}-summary.json", row)
        return planned
    save(
        directory / "observations.json",
        {"observations": observations, "attempts": attempts, "truth": asdict(point)},
    )
    save(
        directory / "probabilities.json",
        integration_tables(banks, observations, config),
    )
    rows = []
    for row in planned:
        method = row["method"]
        candidates = banks[method]
        started = time.perf_counter()
        try:
            fitted, calibration = calibrate(
                candidates,
                observations,
                {**config, "seed": int(row["calibration_seed"])},
                integration_alpha(config, len(candidates)),
                args.bootstrap,
                sampler=sampler,
                null_calibration="grid-supremum",
                scorer=score_bank if method == "histogram" else score_integrated_bank,
            )
            save(
                directory / f"{method}.json",
                {"fit": fitted, "calibration": calibration},
            )
            row.update(
                status="completed",
                reported_event=calibration["mc_p_upper"] <= 0.05,
                point_event=calibration["p_value"] <= 0.05,
                parent_correct=fitted["alternative"]["mul.tree"] == truth,
                mc_p_upper=calibration["mc_p_upper"],
                p_value=calibration["p_value"],
                contrast=fitted["contrast"],
                contrast_width=fitted["contrast_upper"] - fitted["contrast_lower"],
            )
        except Exception as error:
            row.update(
                status="calibration-failed", error=f"{type(error).__name__}: {error}"
            )
        row["seconds"] = time.perf_counter() - started
        save(directory / f"{method}-summary.json", row)
        rows.append(row)
    return rows


def run_probe(args, name):
    index = tuple(SCENARIOS).index(name)
    config = scenario_model(name, args.samples, args.bank_seed + index)
    directory = args.output / name
    directory.mkdir()
    banks = {method: [] for method in METHODS}
    audit = []
    for task in make_tasks(
        tree(SPECIES), "X", "A B", validate_model(config, tree(SPECIES)), 100
    )[0]:
        started = time.perf_counter()
        try:
            paired, history_counts = build_paired_banks(
                task,
                config,
                exact_tip_limit=args.exact_tip_limit,
                max_states=args.max_states,
            )
        except Exception as error:
            save(
                directory / "bank-failure.json",
                {
                    "candidate": task[0],
                    "grid": task[2],
                    "error": f"{type(error).__name__}: {error}",
                },
            )
            save(directory / "partial-bank-audit.json", audit)
            rows = []
            for case_index, case in enumerate(scenario_cases(name)):
                for r in range(args.replicates):
                    planned = case_records(args, name, case_index, r, case)
                    for row in planned:
                        row.update(
                            status="bank-failed",
                            error=f"{type(error).__name__}: {error}",
                        )
                    rows.extend(planned)
            save(directory / "failed-trials.json", rows)
            return rows
        for method in METHODS:
            if method != "histogram":
                paired[method].interval_method = getattr(
                    args, "interval", "empirical-bernstein"
                )
            banks[method].append(paired[method])
        audit.append(
            {
                "candidate": task[0],
                "grid": task[2],
                **history_counts,
                "paired_build_seconds": time.perf_counter() - started,
            }
        )
        print(f"{name}: paired bank {len(audit)}/20", flush=True)
    for method, candidates in banks.items():
        record = bank_record if method == "histogram" else integrated_record
        save(directory / f"{method}-banks.json", [record(bank) for bank in candidates])
    save(directory / "bank-audit.json", audit)
    rows = []
    for case_index, case in enumerate(scenario_cases(name)):
        for replicate in range(args.replicates):
            rows.extend(
                evaluate_case(args, name, case_index, replicate, case, banks, config)
            )
            pd.DataFrame(rows).to_csv(directory / "summary.tsv", sep="\t", index=False)
            print(
                f"{name}: {case[0]} r{replicate + 1} {[r['status'] for r in rows[-3:]]}",
                flush=True,
            )
    return rows


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, required=True)
    cli.add_argument("--source", type=Path)
    cli.add_argument(
        "--scenarios", nargs="+", choices=tuple(SCENARIOS), default=list(SCENARIOS)
    )
    cli.add_argument("--samples", type=int, default=1000)
    cli.add_argument("--families", type=int, default=10)
    cli.add_argument("--bootstrap", type=int, default=19)
    cli.add_argument("--replicates", type=int, default=1)
    cli.add_argument("--seed", type=int, default=20261101)
    cli.add_argument("--bank-seed", type=int, default=20261102)
    cli.add_argument("--calibration-seed", type=int, default=20261103)
    cli.add_argument("--exact-tip-limit", type=int, default=4)
    cli.add_argument("--max-states", type=int, default=100000)
    cli.add_argument(
        "--reference", choices=("msprime", "conditional"), default="msprime"
    )
    cli.add_argument(
        "--interval",
        choices=("empirical-bernstein", "chernoff-kl"),
        default="empirical-bernstein",
    )
    args = cli.parse_args()
    if (
        min(
            args.samples,
            args.families,
            args.bootstrap,
            args.replicates,
            args.max_states,
        )
        < 1
        or min(args.seed, args.bank_seed, args.calibration_seed, args.exact_tip_limit)
        < 0
        or len(set(args.scenarios)) != len(args.scenarios)
    ):
        cli.error(
            "Positive budgets, nonnegative seeds/hidden-tip limit and unique scenarios required."
        )
    logging.getLogger("msprime").setLevel(logging.WARNING)
    root = Path(__file__).resolve().parents[2]
    args.output.mkdir(parents=True, exist_ok=False)
    save(
        args.output / "protocol.json",
        {
            "arguments": {
                k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()
            },
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "sources": {name: digest(root / name) for name in SOURCES},
            "models": {
                name: scenario_model(
                    name, args.samples, args.bank_seed + tuple(SCENARIOS).index(name)
                )
                for name in args.scenarios
            },
            "exact_tip_limit": args.exact_tip_limit,
            "decision": "max upper score-MC P <= 0.05; failed analyses retained",
            "purpose": "Independent exploratory paired probe, not a power/5% error guarantee or a retuning of historical results",
            "selection": "2 <= observed tips <= 4; full hidden histories retained",
            "interval": args.interval
            + "; union budget covers two strata, every possible pattern and selection denominator",
            "reference": args.reference
            + "; topology only, not empirical inference-error calibration",
            "source_protocol_sha256": digest(args.source / "protocol.json")
            if args.source
            else None,
        },
    )
    if args.source:
        diagnose(args.source, args.output)
    rows = []
    for name in args.scenarios:
        rows.extend(run_probe(args, name))
        pd.DataFrame(rows).to_csv(args.output / "summary.tsv", sep="\t", index=False)
    return int(any(row["status"] != "completed" for row in rows))


if __name__ == "__main__":
    raise SystemExit(main())
