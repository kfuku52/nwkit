"""Run reproducible paired shift/convergence trials against independently drawn OU data.

Run with PYTHONPATH=. from the checkout. Outputs require a new directory.
Bootstrap is a separate run so failures cannot erase point-estimation results.
"""

import argparse
import contextlib
import hashlib
import itertools
import json
import platform
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from shift_simulation_cases import SCENARIOS, selection_metrics, simulate, summarize

from nwkit import __version__
from nwkit.cli import parser


def write_json(path, value):
    path.write_text(
        json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )


def fit(directory, truth, options, mode):
    destination = directory / mode
    destination.mkdir()
    args = [
        "shift",
        "--selection",
        "ic",
        "-i",
        str(directory / "tree.nwk"),
        "--trait",
        str(directory / "traits.tsv"),
        "--state-column",
        "value",
        "--standard-error-column",
        "se",
        "--model-out",
        str(destination / "model.json"),
        "-o",
        str(destination / "regimes.tsv"),
        "--rscript",
        options.rscript,
        "--max-shifts",
        "2",
        "--criterion",
        options.criterion,
        "--search-strategy",
        "exhaustive",
        "--root-model",
        truth["root_model"],
    ]
    if mode != "shift":
        args.append("--convergence")
    if mode == "bootstrap":
        args.extend(
            [
                "--bootstrap",
                str(options.bootstrap),
                "--bootstrap-seed",
                str(truth["seed"] % 2147483647),
            ]
        )
    row = {
        "case": directory.name,
        "mode": mode,
        **{
            key: truth[key]
            for key in (
                "seed",
                "scenario",
                "tips",
                "root_model",
                "standard_error",
                "effect",
            )
        },
    }
    with (
        (destination / "run.log").open("w", encoding="utf-8") as log,
        contextlib.redirect_stdout(log),
        contextlib.redirect_stderr(log),
    ):
        try:
            parsed = parser.parse_args(args)
            parsed.handler(parsed)
            model = json.loads((destination / "model.json").read_text())
            row.update(status="completed", **selection_metrics(model, truth))
        except Exception as error:
            row.update(status="failed", error=f"{type(error).__name__}: {error}")
    write_json(destination / "result.json", row)
    return row


def validate_options(options):
    if (
        options.replicates < 1
        or not 0 <= options.bootstrap <= 2147483647
        or options.seed < 0
    ):
        raise ValueError(
            "replicates must be positive; seed nonnegative; bootstrap between 0 and 2147483647"
        )
    for name in ("tips", "standard_errors", "effects"):
        values = getattr(options, name)
        if not values or len(values) != len(set(values)):
            raise ValueError(f"{name} must contain unique values")
    for tips, se, effect in itertools.product(
        options.tips, options.standard_errors, options.effects
    ):
        simulate(
            tips=tips,
            scenario="null",
            root_model="OUfixedRoot",
            se=se,
            effect=effect,
            seed=options.seed,
        )


def run(options):
    validate_options(options)
    options.output.mkdir(parents=True, exist_ok=False)
    source_root = Path(__file__).resolve().parents[1]
    sources = [
        Path(__file__),
        Path(__file__).with_name("shift_simulation_cases.py"),
        *sorted((source_root / "nwkit").glob("shift*.py")),
    ]
    manifest = {
        "schema_version": 1,
        "nwkit_version": __version__,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "platform": platform.platform(),
        "options": {**vars(options), "output": str(options.output)},
        "source_sha256": {
            str(path.relative_to(source_root)): hashlib.sha256(
                path.read_bytes()
            ).hexdigest()
            for path in sources
        },
    }
    write_json(options.output / "manifest.json", manifest)
    for path in sources:
        snapshot = options.output / "source-snapshot" / path.relative_to(source_root)
        snapshot.parent.mkdir(parents=True, exist_ok=True)
        snapshot.write_bytes(path.read_bytes())
    records = []
    cells = itertools.product(
        options.tips,
        SCENARIOS,
        ("OUfixedRoot", "OUrandomRoot"),
        options.standard_errors,
        options.effects,
    )
    for cell_index, (tips, scenario, root_model, se, effect) in enumerate(cells):
        for replicate in range(options.replicates):
            seed = int(
                np.random.SeedSequence(
                    [options.seed, cell_index, replicate]
                ).generate_state(1)[0]
            )
            newick, truth = simulate(
                tips=tips,
                scenario=scenario,
                root_model=root_model,
                se=se,
                effect=effect,
                seed=seed,
            )
            directory = options.output / f"c{cell_index:03d}-r{replicate:04d}"
            directory.mkdir()
            (directory / "tree.nwk").write_text(newick + "\n")
            write_json(directory / "truth.json", truth)
            pd.DataFrame(
                {
                    "leaf_name": truth["tip_names"],
                    "value": truth["observations"],
                    "se": se,
                }
            ).to_csv(directory / "traits.tsv", sep="\t", index=False)
            modes = (
                ("shift", "convergence", "bootstrap")
                if options.bootstrap
                else ("shift", "convergence")
            )
            for mode in modes:
                row = fit(directory, truth, options, mode)
                records.append(row)
                with (options.output / "records.jsonl").open("a") as stream:
                    stream.write(json.dumps(row, allow_nan=False) + "\n")
                print(
                    f"{directory.name} {scenario} {root_model} SE={se} {mode}: {row['status']}",
                    flush=True,
                )
    groups = defaultdict(list)
    keys = ("tips", "scenario", "root_model", "standard_error", "effect", "mode")
    for row in records:
        groups[tuple(row[key] for key in keys)].append(row)
    summary = [
        {**dict(zip(keys, cell, strict=True)), **summarize(rows)}
        for cell, rows in groups.items()
    ]
    write_json(options.output / "summary.json", summary)
    pd.DataFrame(records).to_csv(options.output / "records.csv", index=False)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", required=True, type=Path)
    cli.add_argument("--rscript", default="Rscript")
    cli.add_argument("--replicates", type=int, default=5)
    cli.add_argument("--bootstrap", type=int, default=3)
    cli.add_argument("--seed", type=int, default=20260909)
    cli.add_argument("--tips", type=int, nargs="+", default=[8])
    cli.add_argument("--standard-errors", type=float, nargs="+", default=[0.0, 0.2])
    cli.add_argument("--effects", type=float, nargs="+", default=[2.0])
    cli.add_argument("--criterion", choices=["BIC", "pBIC", "AICc"], default="BIC")
    options = cli.parse_args()
    try:
        validate_options(options)
    except ValueError as error:
        cli.error(str(error))
    run(options)


if __name__ == "__main__":
    main()
