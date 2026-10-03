"""Independent ancestry/sequence simulation and paired conditional parent fits."""

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
import scipy
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

from nwkit.mul_msc import dated_text
from nwkit.mul_msc_fit import FitSettings, fit_candidate, species_topology_signature
from nwkit.mul_msc_model import dated_candidates, validate_sampling
from nwkit.mul_reconcile_model import Reconciliation
from nwkit.species_parser import SpeciesParser
from nwkit.util import read_tree

SPECIES = "(((A:2,X:2):1,B:3):2,C:5);"
AGE_ERROR_SPECIES = "(((A:2.4,X:2.4):0.6,B:3):2,C:5);"


def tree(text):
    return read_tree("[&R]" + text, "auto", True, quiet=True)


def parser():
    return SpeciesParser(species_regex=r".*_([^_]+)$")


def cases():
    return [
        {
            "name": f"standard-{parent}-{ne}",
            "parent": parent,
            "ne": ne,
            "kind": "standard",
        }
        for ne in (0.5, 5.0)
        for parent in ("A", "B", "C")
    ] + [
        {"name": kind, "parent": "B", "ne": 5.0, "kind": kind}
        for kind in ("missing", "age-error", "ghost-proxy", "omitted-parent")
    ]


def simulate_families(case, count, rng):
    donor_age = 1.3 if case["kind"] == "ghost-proxy" else 0.7
    truth = next(
        c
        for c in dated_candidates(tree(SPECIES), "X", case["parent"], donor_age)[0]
        if c.status == "evaluated"
    )
    population_tree = truth.tree.copy()
    labels = {}
    for i, node in enumerate(population_tree.traverse("postorder")):
        if node.is_leaf:
            labels[f"p{i}"] = {
                "A": "a_A",
                "B": "b_B",
                "C": "c_C",
                "X+": "x1_X",
                "X*": "x2_X",
            }[node.name]
        node.name = f"p{i}"
    demography = msprime.Demography.from_species_tree(
        dated_text(population_tree)[4:], initial_size=case["ne"], time_units="gen"
    )
    true_genes, estimated_genes, alignments, audits = [], [], [], []
    for family in range(count):
        names = dict(labels)
        if rng.random() < 0.5:
            names = {
                pop: ("x2_X" if name == "x1_X" else "x1_X" if name == "x2_X" else name)
                for pop, name in names.items()
            }
        if case["kind"] == "missing" and rng.random() < 0.25:
            missing = rng.choice(
                [pop for pop, name in names.items() if name.endswith("_X")]
            )
            names.pop(missing)
        samples = [msprime.SampleSet(1, population=pop, ploidy=1) for pop in names]
        ancestry_seed, mutation_seed = (int(rng.integers(1, 2**31)) for _ in range(2))
        ancestry = msprime.sim_ancestry(
            samples=samples,
            demography=demography,
            sequence_length=600,
            recombination_rate=0,
            ploidy=2,
            model="hudson",
            random_seed=ancestry_seed,
        )
        sample_labels = {
            int(sample): names[
                demography.populations[ancestry.node(sample).population].name
            ]
            for sample in ancestry.samples()
        }
        gene = tree(ancestry.first().as_newick(node_labels=sample_labels, precision=17))
        mutated = msprime.sim_mutations(
            ancestry, rate=0.003, model=msprime.JC69(), random_seed=mutation_seed
        )
        reference = "".join(rng.choice(list("ACGT"), size=600))
        sequences = dict(
            zip(
                sample_labels.values(),
                mutated.alignments(reference_sequence=reference),
                strict=True,
            )
        )
        estimated = estimate_tree(sequences)
        true_genes.append(gene)
        estimated_genes.append(estimated)
        alignments.append({"family": family + 1, "sequences": sequences})
        audits.append(
            {
                "family": family + 1,
                "ancestry_seed": ancestry_seed,
                "mutation_seed": mutation_seed,
                "num_tips": len(names),
                "topology_correct": labeled_topology(gene)
                == labeled_topology(estimated),
                "alignment_sha256": hashlib.sha256(
                    json.dumps(sequences, sort_keys=True).encode()
                ).hexdigest(),
            }
        )
    return true_genes, estimated_genes, alignments, audits


def labeled_topology(gene):
    values = {}
    for node in gene.traverse("postorder"):
        values[node] = (
            ("tip", node.name)
            if node.is_leaf
            else ("node", *sorted(values[c] for c in node.children))
        )
    return values[gene]


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
                raise ArithmeticError(
                    "JC69 distance saturated; dataset is not silently pruned."
                )
            row.append(-0.75 * math.log1p(-4 * p / 3))
        matrix.append([*row, 0.0])
    estimated = DistanceTreeConstructor().nj(DistanceMatrix(names, matrix))
    estimated.root_at_midpoint()
    output = StringIO()
    Phylo.write(estimated, output, "newick")
    return tree(output.getvalue())


def compare(case, genes):
    species = tree(AGE_ERROR_SPECIES if case["kind"] == "age-error" else SPECIES)
    h2 = "A C" if case["kind"] == "omitted-parent" else None
    candidates, polyploid = dated_candidates(
        species, "X", h2, None, age_bounds=(0.1, 1.9)
    )
    validate_sampling(genes, parser(), polyploid, species.leaf_names())
    valid = [c for c in candidates if c.status == "evaluated"]
    settings = FitSettings("joint", (0.1, 1.9), (0.1, 10), grid_points=3, starts=3)
    fits = [fit_candidate(c, genes, parser(), settings) for c in valid]
    fits.sort(key=lambda fit: (-fit["log_likelihood"], fit["id"]))
    parent = {c.id: c.h2 for c in valid}
    dl = [
        (math.fsum(Reconciliation(g, c.tree, parser()).score for g in genes), c.id)
        for c in valid
    ]
    dl.sort()
    best = fits[0]
    msc_ties = [
        parent[f["id"]]
        for f in fits
        if best["log_likelihood"] - f["log_likelihood"] <= 1e-7
    ]
    dl_ties = [parent[id_] for score, id_ in dl if score == dl[0][0]]
    estimates = best["estimates"]
    summary = {
        "msc_parent": parent[best["id"]],
        "dl_parent": parent[dl[0][1]],
        "msc_correct": len(msc_ties) == 1 and msc_ties[0] == case["parent"],
        "dl_correct": len(dl_ties) == 1 and dl_ties[0] == case["parent"],
        "msc_ties": msc_ties,
        "dl_ties": dl_ties,
        "parameter_status": best["diagnostics"]["status"],
        "ne_estimate": estimates.get("effective_population_size"),
        "attachment_age_estimate": estimates.get("hybridization_age"),
        "ne_error": None
        if estimates.get("effective_population_size") is None
        else estimates["effective_population_size"] - case["ne"],
        "attachment_error": None
        if estimates.get("hybridization_age") is None
        else estimates["hybridization_age"]
        - (1.3 if case["kind"] == "ghost-proxy" else 0.7),
        "unique_species_topologies": len(
            {species_topology_signature(g, parser()) for g in genes}
        ),
    }
    evidence = {
        "summary": summary,
        "fits": [{k: v for k, v in fit.items() if k != "candidate"} for fit in fits],
        "dl_scores": [{"parent": parent[id_], "score": score} for score, id_ in dl],
        "settings": vars(settings),
    }
    return summary, evidence


def main():
    logging.getLogger("msprime").setLevel(logging.WARNING)
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output", type=Path, required=True)
    cli.add_argument("--seed", type=int, default=20261018)
    cli.add_argument("--replicates", type=int, default=2)
    cli.add_argument("--families", type=int, default=100)
    cli.add_argument("--cases", nargs="+", default=None)
    args = cli.parse_args()
    if min(args.replicates, args.families) < 1:
        cli.error("replicates/families must be positive")
    args.output.mkdir(parents=True, exist_ok=False)
    selected = [
        case for case in cases() if args.cases is None or case["name"] in args.cases
    ]
    if not selected or (
        args.cases is not None
        and set(args.cases) != {case["name"] for case in selected}
    ):
        cli.error("unknown/empty case selection")
    source = Path(__file__).resolve().parents[2]
    paths = [
        Path(__file__),
        *(
            source / "nwkit" / name
            for name in ("mul_coalescent.py", "mul_msc_model.py", "mul_msc_fit.py")
        ),
    ]
    protocol = {
        "seed": args.seed,
        "replicates": args.replicates,
        "families": args.families,
        "cases": selected,
        "python": sys.version,
        "platform": platform.platform(),
        "versions": {
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "msprime": msprime.__version__,
            "biopython": Bio.__version__,
        },
        "source_sha256": {
            str(path.relative_to(source)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in paths
        },
        "interpretation": "Conditional known-species-age parent pilot; not WGD detection, global identifiability or calibrated support; no Ne/attachment truth supplied to fitting.",
    }
    (args.output / "protocol.json").write_text(
        json.dumps(protocol, indent=2, allow_nan=False) + "\n"
    )
    summaries = []
    for case in selected:
        # Stable case-position seeds keep a selected case identical to a full run.
        case_index = next(
            i for i, value in enumerate(cases()) if value["name"] == case["name"]
        )
        for replicate in range(args.replicates):
            directory = args.output / f"{case['name']}-r{replicate + 1}"
            directory.mkdir()
            rng = np.random.default_rng(
                np.random.SeedSequence([args.seed, case_index, replicate])
            )
            try:
                true, estimated, alignments, audits = simulate_families(
                    case, args.families, rng
                )
                (directory / "alignments.jsonl").write_text(
                    "".join(json.dumps(row) + "\n" for row in alignments)
                )
                pd.DataFrame(audits).to_csv(
                    directory / "families.tsv", sep="\t", index=False
                )
                for label, genes in (("true", true), ("estimated", estimated)):
                    (directory / f"{label}.nwk").write_text(
                        "".join(dated_text(g) + "\n" for g in genes)
                    )
                    summary, evidence = compare(case, genes)
                    (directory / f"{label}.json").write_text(
                        json.dumps(evidence, indent=2, allow_nan=False) + "\n"
                    )
                    row = {
                        "case": case["name"],
                        "replicate": replicate + 1,
                        "trees": label,
                        "true_parent": case["parent"],
                        "true_ne": case["ne"],
                        "families": args.families,
                        "gene_topology_accuracy": sum(
                            a["topology_correct"] for a in audits
                        )
                        / len(audits),
                        **summary,
                    }
                    summaries.append(row)
                    pd.DataFrame(summaries).to_csv(
                        args.output / "summary.tsv", sep="\t", index=False
                    )
                    print(json.dumps(row, allow_nan=False), flush=True)
            except Exception as error:
                (directory / "failure.json").write_text(
                    json.dumps({"error": str(error), "type": type(error).__name__})
                    + "\n"
                )
                raise


if __name__ == "__main__":
    main()
