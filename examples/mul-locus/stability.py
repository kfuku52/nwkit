"""Separate site-bootstrap parent stability for a completed locus pilot."""

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from pilot import PARSER, SPECIES, estimate_tree, tree

from nwkit.mul_locus_cli import finite_json
from nwkit.mul_locus_mc import (
    LocusBank,
    integration_alpha,
    make_tasks,
    search_banks,
    validate_model,
)
from nwkit.mul_msc_fit import species_topology_signature


def signature(value):
    return tuple(signature(v) if isinstance(v, list) else v for v in value)


def histogram(values):
    return Counter({signature(v["signature"]): v["hits"] for v in values})


def load_banks(study, model):
    tasks, _ = make_tasks(
        tree(SPECIES), "X", "A B", validate_model(model, tree(SPECIES)), 100
    )
    lookup = {(t[0], t[2]): t for t in tasks}
    banks = []
    for record in json.loads((study / "banks.json").read_text()):
        candidate, h2, grid, population, point = lookup[
            (record["candidate"], record["grid"])
        ]
        strata = record.get("strata")
        if strata is not None:
            strata = [{**s, "counts": histogram(s["counts"])} for s in strata]
        banks.append(
            LocusBank(
                candidate,
                h2,
                grid,
                population,
                point,
                histogram(record["counts"]),
                record["samples"],
                record["attempts"],
                strata,
            )
        )
    if len(banks) != len(tasks) or len({(b.candidate, b.grid) for b in banks}) != len(
        tasks
    ):
        raise ValueError("Study must contain exactly every declared integration bank.")
    return banks


def resample_alignment(sequences, rng):
    lengths = {len(s) for s in sequences.values()}
    if len(lengths) != 1 or not lengths or next(iter(lengths)) < 1:
        raise ValueError("Site bootstrap requires a nonempty equal-length alignment.")
    indices = rng.integers(next(iter(lengths)), size=next(iter(lengths)))
    return {name: "".join(seq[i] for i in indices) for name, seq in sequences.items()}


def dataset(study, name, banks, model, replicates, seed):
    path = study / name / "families.jsonl"
    families = [json.loads(line) for line in path.read_text().splitlines()]
    rows, audit = [], []
    alpha = integration_alpha(model, len(banks))
    for replicate in range(replicates):
        rng = np.random.default_rng(np.random.SeedSequence([*seed, replicate]))
        try:
            genes = [
                estimate_tree(resample_alignment(f["sequences"], rng)) for f in families
            ]
            observed = [species_topology_signature(g, PARSER) for g in genes]
            fitted = search_banks(banks, observed, alpha)
            alternatives = [r for r in fitted["rows"] if r["mul.tree"]]
            best_lower = max(r["mc_lower"] for r in alternatives)
            overlap = sorted(
                {r["mul.tree"] for r in alternatives if r["mc_upper"] >= best_lower}
            )
            row = {
                "dataset": name,
                "replicate": replicate + 1,
                "families": len(families),
                "status": "completed",
                "selected_alternative": fitted["alternative"]["mul.tree"],
                "null_grid": fitted["null"]["grid"],
                "alternative_grid": fitted["alternative"]["grid"],
                "contrast": fitted["contrast"],
                "contrast_lower": fitted["contrast_lower"],
                "contrast_upper": fitted["contrast_upper"],
                "mc_overlapping_parents": json.dumps(overlap),
            }
            audit.append(
                {"replicate": replicate + 1, "observations": observed, "fit": fitted}
            )
        except Exception as error:
            row = {
                "dataset": name,
                "replicate": replicate + 1,
                "status": "failed",
                "error": str(error),
            }
            audit.append(row)
        rows.append(row)
    return rows, audit


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--study", type=Path, required=True)
    cli.add_argument("--output", type=Path, required=True)
    cli.add_argument("--replicates", type=int, default=19)
    cli.add_argument("--seed", type=int, default=20261029)
    args = cli.parse_args()
    if args.replicates < 1 or args.seed < 0:
        cli.error("Positive replicates and nonnegative seed required.")
    protocol = json.loads((args.study / "protocol.json").read_text())
    summary = pd.read_csv(args.study / "summary.tsv", sep="\t")
    expected = 8 * protocol["arguments"]["replicates"]
    if len(summary) != expected:
        raise ValueError("Original paired study must finish, retaining all failures.")
    names = sorted(p.parent.name for p in args.study.glob("*/families.jsonl"))
    if len(names) != expected // 2:
        raise ValueError(
            "Every original dataset must retain its complete family alignments."
        )
    if any(
        len((args.study / n / "families.jsonl").read_text().splitlines())
        != protocol["arguments"]["families"]
        for n in names
    ):
        raise ValueError("Original family counts must be preserved.")
    banks = load_banks(args.study, protocol["model"])
    args.output.mkdir(parents=True, exist_ok=False)
    files = [
        args.study / "protocol.json",
        args.study / "banks.json",
        args.study / "summary.tsv",
        *[args.study / n / "families.jsonl" for n in names],
    ]
    inputs = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    source = Path(__file__).resolve().parents[2]
    sources = {
        str(p.relative_to(source)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in [
            Path(__file__),
            Path(__file__).with_name("pilot.py"),
            source / "nwkit/mul_locus_mc.py",
            source / "nwkit/mul_locus.py",
            source / "nwkit/mul_locus_cli.py",
            source / "nwkit/mul_msc_model.py",
            source / "nwkit/mul_msc_fit.py",
        ]
    }
    (args.output / "protocol.json").write_text(
        json.dumps(
            {
                "seed": args.seed,
                "replicates": args.replicates,
                "inputs": inputs,
                "sources": sources,
                "method": "within-family-column-bootstrap-JC69-NJ-midpoint",
                "limitations": "Descriptive parent stability conditional on observed alignments and finite integration banks; not event P-values, posterior support or biological confidence intervals.",
            },
            indent=2,
        )
        + "\n"
    )
    rows = []
    for i, name in enumerate(names):
        values, audit = dataset(
            args.study, name, banks, protocol["model"], args.replicates, [args.seed, i]
        )
        rows.extend(values)
        (args.output / (name + ".json")).write_text(
            json.dumps(finite_json(audit)) + "\n"
        )
        pd.DataFrame(finite_json(rows)).to_csv(
            args.output / "summary.tsv", sep="\t", index=False
        )
        print(
            f"{name}: {sum(r['status'] == 'completed' for r in values)}/{len(values)} completed",
            flush=True,
        )
    failures = sum(r["status"] != "completed" for r in rows)
    if failures:
        raise SystemExit(f"{failures} bootstrap analyses failed; all retained")


if __name__ == "__main__":
    main()
