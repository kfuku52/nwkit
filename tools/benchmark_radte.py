"""Reproducible RADTE simulation/CLI benchmark; Unix /usr/bin/time required.

All runs include interpreter startup, reconciliation, sequence prefit, dating,
and outputs. Peak RSS is the largest process, not simultaneous process-tree RSS.
PAML is a different estimator with soft priors; timings are not speedup claims.
"""

import argparse
import json
import platform
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from ete4 import Tree
from radte_benchmark_cases import calibration_intervals, gene_chronogram

ROOT = Path(__file__).resolve().parents[1]


def species_newick(names, age=10.0):
    if len(names) == 1:
        return names[0], 0.0
    middle = len(names) // 2
    children = [
        species_newick(part, age / 2) for part in [names[:middle], names[middle:]]
    ]
    return "(" + ",".join(
        f"{tree}:{age - child_age}" for tree, child_age in children
    ) + ")", age


def simulate(
    directory,
    count,
    sites,
    seed,
    sd,
    *,
    scenario="root",
    width=0.0,
    copy_ratio=1.0,
    mapping_errors=0,
):
    directory.mkdir(parents=True, exist_ok=True)
    names = [f"S{i}" for i in range(count)]
    text, _ = species_newick(names)
    species = Tree(text + ";", parser=1)
    for i, node in enumerate(species.traverse()):
        if not node.is_leaf:
            node.name = f"N{i}"
    species.write(
        outfile=str(directory / "species.nwk"), parser=1, format_root_node=True
    )
    gene, true_age = gene_chronogram(species, scenario)
    (directory / "truth.json").write_text(
        json.dumps(dict(duplication_age=true_age, scenario=scenario))
    )
    bounds = calibration_intervals(species, width)
    if bounds is not None:
        (directory / "bounds.tsv").write_text(bounds)
    rng = np.random.default_rng(seed)
    for node in gene.traverse():
        if node is not gene:
            shift = (
                copy_ratio if all(n.name.endswith("_2") for n in node.leaves()) else 1.0
            )
            node.dist *= 0.01 * shift * np.exp(rng.normal(0, sd))
    gene.write(outfile=str(directory / "gene.nwk"), parser=1, format_root_node=True)
    mapped = {n.name: n.name.rsplit("_", 1)[0] for n in gene.leaves()}
    for name in list(mapped)[:mapping_errors]:
        mapped[name] = names[(names.index(mapped[name]) + 1) % len(names)]
    mapping = "leaf_name\tspecies_label\n" + "".join(
        f"{name}\t{label}\n" for name, label in mapped.items()
    )
    (directory / "mapping.tsv").write_text(mapping)
    sequences = {gene: rng.integers(0, 4, sites)}
    for node in gene.traverse():
        if node is gene:
            continue
        state = sequences[node.up].copy()
        redraw = rng.random(sites) > np.exp(-4 * node.dist / 3)
        state[redraw] = rng.integers(0, 4, redraw.sum())
        sequences[node] = state
    alphabet = np.array(list("ACGT"))
    (directory / "alignment.fasta").write_text(
        "".join(
            f">{n.name}\n{''.join(alphabet[sequences[n]])}\n" for n in gene.leaves()
        )
    )


def command(directory, method, prefix, seed, args):
    cmd = [
        sys.executable,
        "-m",
        "nwkit",
        "radte",
        "--gene-tree",
        str(directory / "gene.nwk"),
        "--species-tree",
        str(directory / "species.nwk"),
        "--species-map-tsv",
        str(directory / "mapping.tsv"),
        "--reconcile",
        "lca",
        "--max-age",
        "100",
        "--seed",
        str(seed),
        "--out-prefix",
        str(prefix),
    ]
    if (directory / "bounds.tsv").is_file():
        cmd += ["--species-node-bounds-tsv", str(directory / "bounds.tsv")]
    if method != "branch":
        cmd += [
            "--alignment",
            str(directory / "alignment.fasta"),
            "--substitution-model",
            "jc69",
        ]
    if method.startswith("paml"):
        cmd += [
            "--mcmctree-likelihood",
            "approximate" if method == "paml-approximate" else "exact",
            "--backend",
            "mcmctree",
            "--mcmctree-samples",
            str(args.paml_samples),
            "--mcmctree-burnin",
            "2000",
        ]
    else:
        if method == "sequence":
            cmd += ["--gamma-categories", "1"]
        cmd += ["--uncertainty", "laplace"]
    return cmd


def timed_run(cmd, prefix):
    log = Path(str(prefix) + ".process.log")
    timing = Path(str(prefix) + ".time.txt")
    flag = "-l" if sys.platform == "darwin" else "-v"
    start = time.perf_counter()
    with log.open("w") as handle:
        result = subprocess.run(
            ["/usr/bin/time", flag, "-o", str(timing), *cmd],
            cwd=ROOT,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    seconds = time.perf_counter() - start
    text = timing.read_text()
    if sys.platform == "darwin":
        match = re.search(r"(\d+)\s+maximum resident set size", text)
        divisor = 1024**2
    else:
        match = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", text)
        divisor = 1024
    rss = int(match[1]) / divisor if match else None
    return dict(
        seconds=seconds,
        peak_process_rss_mib=rss,
        exit_code=result.returncode,
        command=cmd,
    )


def inspect_result(prefix, row):
    if row["exit_code"]:
        return row
    nodes = pd.read_csv(str(prefix) + ".nodes.tsv", sep="\t")
    shared = pd.read_csv(str(prefix) + ".shared-ages.tsv", sep="\t")
    manifest = json.loads(Path(str(prefix) + ".manifest.json").read_text())
    root = nodes.loc[nodes.gene_name == "D"].iloc[0]
    truth = json.loads((prefix.parent / "truth.json").read_text())["duplication_age"]
    row.update(
        age=float(root.estimated_age),
        true_age=truth,
        error=float(root.estimated_age - truth),
        interval_available=bool(pd.notna(root.interval_lower)),
        covered=bool(root.interval_lower <= truth <= root.interval_upper),
        method=manifest["method"],
        diagnostics=manifest["diagnostics"],
        shared_max_difference=float(shared.max_member_age_difference.max()),
        bounds_satisfied=bool(nodes.within_original_bounds.all()),
    )
    if row["shared_max_difference"] != 0:
        raise RuntimeError("Shared-age contract failed")
    if (
        manifest["calibration_policy"] == "hard-all-events"
        and not row["bounds_satisfied"]
    ):
        raise RuntimeError("Hard-bound contract failed")
    return row


def run_case(args, directory, method, seed):
    rows = []
    for repeat in range(-args.warmups, args.repeats):
        prefix = directory / f"{method}-{repeat}"
        row = inspect_result(
            prefix, timed_run(command(directory, method, prefix, seed, args), prefix)
        )
        row.update(repeat=repeat, seed=seed, requested_method=method)
        if repeat >= 0:
            rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--species", type=int, nargs="+", default=[2, 8, 32])
    parser.add_argument("--families", type=int, default=10)
    parser.add_argument("--sites", type=int, default=2000)
    parser.add_argument("--rate-sd", type=float, default=0.3)
    parser.add_argument(
        "--scenario", choices=["root", "nested", "loss"], default="root"
    )
    parser.add_argument("--calibration-width", type=float, default=0.0)
    parser.add_argument("--copy-rate-ratio", type=float, default=1.0)
    parser.add_argument("--mapping-errors", type=int, default=0)
    parser.add_argument(
        "--methods",
        nargs="+",
        choices=["branch", "sequence", "paml", "paml-approximate"],
        default=["branch", "sequence"],
    )
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--warmups", type=int, default=1)
    parser.add_argument("--paml-samples", type=int, default=20000)
    args = parser.parse_args()
    if (
        min(args.species) < 2
        or args.families < 1
        or args.repeats < 1
        or args.warmups < 0
        or args.sites < 1
        or not 0 <= args.calibration_width < 1
        or args.copy_rate_ratio <= 0
        or args.rate_sd < 0
        or args.mapping_errors < 0
    ):
        parser.error(
            "Need >=2 species, positive families/repeats, and nonnegative warmups"
        )
    args.outdir = args.outdir.resolve()
    if args.outdir.exists():
        parser.error("Use a new output directory to preserve previous measurements")
    args.outdir.mkdir(parents=True)
    rows = []
    for count in args.species:
        for family in range(args.families):
            seed = 1701 + family
            directory = args.outdir / f"species-{count}-seed-{seed}"
            simulate(
                directory,
                count,
                args.sites,
                seed,
                args.rate_sd,
                scenario=args.scenario,
                width=args.calibration_width,
                copy_ratio=args.copy_rate_ratio,
                mapping_errors=args.mapping_errors,
            )
            for method in args.methods:
                case = run_case(args, directory, method, seed)
                for row in case:
                    row["species"] = count
                rows.extend(case)
                print(count, seed, method, [r["exit_code"] for r in case], flush=True)
                (args.outdir / "results.json").write_text(
                    json.dumps(
                        dict(
                            environment=dict(
                                python=sys.version,
                                platform=platform.platform(),
                                numpy=np.__version__,
                            ),
                            settings={
                                key: str(value) if isinstance(value, Path) else value
                                for key, value in vars(args).items()
                            },
                            rows=rows,
                        ),
                        indent=2,
                    )
                )


if __name__ == "__main__":
    main()
