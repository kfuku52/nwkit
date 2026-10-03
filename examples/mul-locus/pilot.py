"""Independent SSA/msprime study of the finite-grid DL+ILS event contrast."""

import argparse
import hashlib
import json
import logging
import math
import platform
import sys
from io import StringIO
from pathlib import Path

import Bio
import msprime
import numpy as np
import pandas as pd
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from ete4 import Tree

from nwkit.mul_locus import Locus, LocusParameters
from nwkit.mul_locus_cli import bank_record, finite_json
from nwkit.mul_locus_mc import (
    build_bank,
    calibrate,
    integration_alpha,
    make_tasks,
    validate_model,
)
from nwkit.mul_msc import dated_text
from nwkit.mul_msc_fit import species_topology_signature
from nwkit.species_parser import SpeciesParser
from nwkit.util import compute_node_ages, read_tree

SPECIES = "((A:1,X:1):1,B:2);"
PARSER = SpeciesParser(species_regex=r".*_([^_]+)$")


def tree(text):
    return read_tree("[&R]" + text, "auto", True, quiet=True)


def model(samples):
    return {
        "schema": "nwkit-mul-locus-mc-model-v1",
        "copy_role": "distinct-loci",
        "root_locus_count": 1,
        "species_time_unit": "generations",
        "ancestral_stem": 0.5,
        "detection": {"A": 0.9, "X": 0.9, "B": 0.9},
        "max_observed_tips": 4,
        "parameter_grid": [
            {"duplication": d, "loss": 0.03, "ne": n, "hybridization_age": a}
            for d in (0.03, 0.08)
            for n in (0.3, 1)
            for a in (0.3, 0.7)
        ],
        "samples": samples,
        "seed": 20261026,
        "confidence": 0.99,
        "max_attempts": 1000000,
        "max_locus_nodes": 10000,
        "max_coalescent_states": 100000,
    }


def reference_locus(population, point, config, rng):
    ages = compute_node_ages(population)
    age = ages[population] + config["ancestral_stem"]
    root = Locus(age, "origin")
    alive = [(root, population, False)]
    nodes = 1
    while alive:
        boundary = max(ages[pop] for _, pop, _ in alive)
        rate = len(alive) * (point.duplication + point.loss)
        wait = rng.exponential(1 / rate) if rate else math.inf
        if age - wait > boundary:
            age -= wait
            i = int(rng.integers(len(alive)))
            parent, pop, daughter = alive.pop(i)
            dup = rng.random() < point.duplication / (point.duplication + point.loss)
            node = Locus(age, "duplication" if dup else "loss", daughter=daughter)
            parent.children.append(node)
            if dup:
                alive.extend([(node, pop, False), (node, pop, True)])
            nodes += 1
        else:
            age = boundary
            following = []
            for parent, pop, daughter in alive:
                if ages[pop] != age:
                    following.append((parent, pop, daughter))
                    continue
                node = Locus(
                    age,
                    "tip" if pop.is_leaf else "speciation",
                    pop.props.get("mul_species", pop.name) if pop.is_leaf else "",
                    daughter=daughter,
                )
                parent.children.append(node)
                following.extend((node, child, False) for child in pop.children)
                nodes += 1
            alive = following
        if nodes > config["max_locus_nodes"]:
            raise ValueError("Independent locus reference exceeded node cap.")
    return root


def population_from_locus(locus):
    source = locus.children[0]
    root = Tree()
    mapping = {source: root}
    tips = {}
    bounds = []
    histories = []
    pending = [source]
    while pending:
        node = pending.pop()
        target = mapping[node]
        target.name = f"p{len(histories)}"
        histories.append(
            {
                "name": target.name,
                "age": node.age,
                "kind": node.kind,
                "species": node.species,
                "daughter": node.daughter,
            }
        )
        if node.kind == "tip":
            tips[target.name] = f"g{len(tips)}_{node.species}"
        for child in node.children:
            converted = target.add_child()
            converted.dist = node.age - (child.age if child.children else 0)
            mapping[child] = converted
            pending.append(child)
        if node.daughter:
            descendants = []
            todo = [node]
            while todo:
                child = todo.pop()
                if child.kind == "tip":
                    descendants.append(child)
                todo.extend(child.children)
            bounds.append((node, descendants))
    return root, mapping, tips, bounds, histories


def reference_family(
    population, point, config, rng, estimated=False, *, true_only=False
):
    if true_only and estimated:
        raise ValueError("True-only reference sampling cannot estimate a gene tree.")
    for attempt in range(1, config["max_attempts"] + 1):
        locus = reference_locus(population, point, config, rng)
        root, mapping, labels, bounds, histories = population_from_locus(locus)
        detected = {
            name
            for name in labels.values()
            if rng.random() < config["detection"][PARSER.parse(name).species_label]
        }
        if not 2 <= len(detected) <= config["max_observed_tips"]:
            continue
        demography = msprime.Demography.from_species_tree(
            dated_text(root)[4:], initial_size=point.ne, time_units="gen"
        )
        samples = [msprime.SampleSet(1, population=pop, ploidy=1) for pop in labels]
        for _draw in range(1, 50001):
            seed = int(rng.integers(1, 2**31))
            ancestry = msprime.sim_ancestry(
                samples=samples,
                demography=demography,
                ploidy=2,
                sequence_length=600,
                recombination_rate=0,
                random_seed=seed,
            )
            sample_by_pop = {
                demography.populations[ancestry.node(s).population].name: int(s)
                for s in ancestry.samples()
            }
            genealogy = ancestry.first()
            valid = True
            for daughter, descendants in bounds:
                members = [sample_by_pop[mapping[n].name] for n in descendants]
                if len(members) > 1:
                    mrca = members[0]
                    for member in members[1:]:
                        mrca = genealogy.mrca(mrca, member)
                    parent_age = next(
                        n["age"]
                        for n in histories
                        if n["name"] == mapping[daughter].up.name
                    )
                    if genealogy.time(mrca) > parent_age:
                        valid = False
                        break
            if valid:
                break
        else:
            raise ValueError("Independent bounded-coalescent rejection cap reached.")
        sample_labels = {s: labels[pop] for pop, s in sample_by_pop.items()}
        true = tree(genealogy.as_newick(node_labels=sample_labels, precision=17))
        true.prune(list(detected), preserve_branch_length=True)
        if true_only:
            return (
                true,
                true,
                None,
                {
                    "selection_attempts": attempt,
                    "coalescent_attempts": _draw,
                    "ancestry_seed": seed,
                    "locus_history": histories,
                },
            )
        mutation_seed = int(rng.integers(1, 2**31))
        mutated = msprime.sim_mutations(
            ancestry, rate=0.01, model=msprime.JC69(), random_seed=mutation_seed
        )
        reference = "".join(rng.choice(list("ACGT"), size=600))
        sequences = {
            name: seq
            for name, seq in zip(
                sample_labels.values(),
                mutated.alignments(reference_sequence=reference),
                strict=True,
            )
            if name in detected
        }
        inferred = estimate_tree(sequences)
        audit = {
            "selection_attempts": attempt,
            "coalescent_attempts": _draw,
            "ancestry_seed": seed,
            "mutation_seed": mutation_seed,
            "locus_history": histories,
            "sequences": sequences,
        }
        return inferred if estimated else true, true, inferred, audit
    raise ValueError("Independent family selection cap reached.")


def estimate_tree(sequences):
    names = list(sequences)
    matrix = []
    for i, name in enumerate(names):
        row = []
        for other in names[:i]:
            p = sum(
                a != b for a, b in zip(sequences[name], sequences[other], strict=True)
            ) / len(sequences[name])
            if p >= 0.75:
                raise ValueError("JC69 saturation, no family discarded.")
            row.append(-0.75 * math.log1p(-4 * p / 3))
        matrix.append([*row, 0])
    inferred = DistanceTreeConstructor().nj(DistanceMatrix(names, matrix))
    inferred.root_at_midpoint()
    output = StringIO()
    Phylo.write(inferred, output, "newick")
    return tree(output.getvalue())


def reference_selected(
    population, parameters, model, rng, *, count, estimated=False, true_only=False
):
    values = []
    attempts = 0
    for _ in range(count):
        gene, _, _, audit = reference_family(
            population, parameters, model, rng, estimated, true_only=true_only
        )
        attempts += audit["selection_attempts"]
        values.append(species_topology_signature(gene, PARSER))
    return values, attempts


def reference_true_selected(population, parameters, model, rng, *, count):
    return reference_selected(
        population, parameters, model, rng, count=count, true_only=True
    )


def main():
    logging.getLogger("msprime").setLevel(logging.WARNING)
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, required=True)
    cli.add_argument("--seed", type=int, default=20261025)
    cli.add_argument("--replicates", type=int, default=5)
    cli.add_argument("--families", type=int, default=30)
    cli.add_argument("--samples", type=int, default=10000)
    cli.add_argument("--bank-seed", type=int, default=20261026)
    cli.add_argument(
        "--integration",
        choices=("selected-histogram", "ancestral-stratified"),
        default="selected-histogram",
    )
    cli.add_argument("--true-bootstrap", type=int, default=39)
    cli.add_argument("--estimated-bootstrap", type=int, default=19)
    args = cli.parse_args()
    if (
        min(
            args.replicates,
            args.families,
            args.samples,
            args.true_bootstrap,
            args.estimated_bootstrap,
        )
        < 1
    ):
        cli.error("positive counts required")
    args.output.mkdir(parents=True, exist_ok=False)
    config = model(args.samples)
    config["seed"] = args.bank_seed
    config["integration"] = args.integration
    parameters = validate_model(config, tree(SPECIES))
    tasks, excluded = make_tasks(tree(SPECIES), "X", "A B", parameters, 100)
    source = Path(__file__).resolve().parents[2]
    files = [
        Path(__file__),
        *(
            source / "nwkit" / name
            for name in (
                "mul_locus.py",
                "mul_locus_mc.py",
                "mul_msc_model.py",
                "mul_coalescent.py",
            )
        ),
    ]
    protocol = {
        "arguments": {
            k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()
        },
        "model": config,
        "excluded": excluded,
        "python": sys.version,
        "platform": platform.platform(),
        "versions": {
            "msprime": msprime.__version__,
            "biopython": Bio.__version__,
            "numpy": np.__version__,
        },
        "sources": {
            str(p.relative_to(source)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in files
        },
    }
    (args.output / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    banks = []
    for task in tasks:
        banks.append(build_bank(task, config))
        print(f"bank {len(banks)}/{len(tasks)}", flush=True)
    (args.output / "banks.json").write_text(
        json.dumps([bank_record(bank) for bank in banks]) + "\n"
    )
    alpha = integration_alpha(config, len(banks))
    cases = [
        ("null-low", 0, LocusParameters(0.03, 0.03, 0.3, 0.5)),
        ("null-high", 0, LocusParameters(0.08, 0.03, 1, 0.5)),
        ("allop-A", 1, LocusParameters(0.05, 0.03, 1, 0.5)),
        ("allop-B", 2, LocusParameters(0.05, 0.03, 1, 0.5)),
    ]
    summaries = []
    failures = 0
    for case_index, (name, truth, point) in enumerate(cases):
        population = next(
            task[3]
            for task in make_tasks(tree(SPECIES), "X", "A B", [point], 100)[0]
            if task[0] == truth
        )
        for replicate in range(args.replicates):
            directory = args.output / f"{name}-r{replicate + 1}"
            directory.mkdir()
            rng = np.random.default_rng(
                np.random.SeedSequence([args.seed, case_index, replicate])
            )
            try:
                families = [
                    reference_family(population, point, config, rng)
                    for _ in range(args.families)
                ]
            except Exception as error:
                failure = {
                    "case": name,
                    "replicate": replicate + 1,
                    "status": "generation-failed",
                    "error": str(error),
                }
                (directory / "generation-failure.json").write_text(
                    json.dumps(failure) + "\n"
                )
                summaries.extend(
                    {**failure, "view": view} for view in ("true", "estimated")
                )
                failures += 2
                pd.DataFrame(summaries).to_csv(
                    args.output / "summary.tsv", sep="\t", index=False
                )
                print(json.dumps(failure), flush=True)
                continue
            (directory / "families.jsonl").write_text(
                "".join(json.dumps(f[3]) + "\n" for f in families)
            )
            for view, index in (("true", 1), ("estimated", 2)):
                genes = [f[index] for f in families]
                (directory / (view + ".nwk")).write_text(
                    "".join(dated_text(g) + "\n" for g in genes)
                )
                try:
                    observations = [
                        species_topology_signature(g, PARSER) for g in genes
                    ]
                    if view == "true":
                        fitted, calibration = calibrate(
                            banks, observations, config, alpha, args.true_bootstrap
                        )
                    else:

                        def sampler(p, q, m, r, *, count):
                            return reference_selected(
                                p, q, m, r, count=count, estimated=True
                            )

                        fitted, calibration = calibrate(
                            banks,
                            observations,
                            config,
                            alpha,
                            args.estimated_bootstrap,
                            sampler=sampler,
                        )
                    (directory / (view + ".json")).write_text(
                        json.dumps(
                            finite_json({"fit": fitted, "calibration": calibration}),
                            indent=2,
                        )
                        + "\n"
                    )
                    alt = fitted["alternative"]
                    row = {
                        "case": name,
                        "replicate": replicate + 1,
                        "view": view,
                        "families": args.families,
                        "status": "completed",
                        "true_candidate": truth,
                        "selected_alternative": alt["mul.tree"],
                        "parent_correct": truth != 0 and alt["mul.tree"] == truth,
                        "p_value": calibration["p_value"],
                        "mc_p_lower": calibration["mc_p_lower"],
                        "mc_p_upper": calibration["mc_p_upper"],
                        "reported_event": calibration["mc_p_upper"] <= 0.05,
                        "tree_accuracy": sum(
                            species_topology_signature(f[1], PARSER)
                            == species_topology_signature(f[2], PARSER)
                            for f in families
                        )
                        / len(families),
                    }
                except Exception as error:
                    failures += 1
                    row = {
                        "case": name,
                        "replicate": replicate + 1,
                        "view": view,
                        "families": args.families,
                        "status": "failed",
                        "error": str(error),
                    }
                    (directory / (view + "-failure.json")).write_text(
                        json.dumps(row) + "\n"
                    )
                summaries.append(row)
                pd.DataFrame(summaries).to_csv(
                    args.output / "summary.tsv", sep="\t", index=False
                )
                print(json.dumps(row), flush=True)
    if failures:
        raise SystemExit(f"{failures} analyses failed; retained in summary")


if __name__ == "__main__":
    main()
