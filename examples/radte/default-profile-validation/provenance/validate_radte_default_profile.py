"""External AliSim validation of the GeneGalleon native GY94/profile route.

Generation never calls NWKIT likelihood/rate functions. Each family has a
separate rate realization and alignment. Truth is one named duplication, so
the coverage denominator does not treat within-family nodes as independent.
Point-only and profile runs are recorded separately: a failed profile never
erases a completed point fit. No unavailable interval is replaced by a point.
"""

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import shutil
import subprocess
import sys
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from ete4 import Tree
from scipy.stats import binomtest

from radte_benchmark_cases import gene_chronogram
from radte_interval_simulation import species_text

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def run(command, directory, name, timeout):
    save_json(directory / (name + ".command.json"), command)
    started = time.monotonic()
    with (directory / (name + ".log")).open("w") as log:
        try:
            result = subprocess.run(command, cwd=directory, stdout=log,
                                    stderr=subprocess.STDOUT, timeout=timeout)
            status = "complete" if result.returncode == 0 else "exit-" + str(result.returncode)
        except subprocess.TimeoutExpired:
            status = "timeout"
    return status, time.monotonic() - started


def generate_trees(directory, case, seed):
    """Shared chronology fixtures only; independent edge rates, no native Q."""
    text, _ = species_text([f"S{i}" for i in range(case["species"])], "balanced")
    species = Tree(text + ";", parser=1)
    for i, node in enumerate(species.traverse()):
        if not node.is_leaf:
            node.name = f"N{i}"
    species.write(outfile=str(directory / "species.nwk"), parser=1, format_root_node=True)
    gene, truth = gene_chronogram(species, case.get("scenario", "root"))
    # Save time tree before substitution lengths are assigned.
    gene.write(outfile=str(directory / "truth.nwk"), parser=1, format_root_node=True)
    rng = np.random.default_rng(seed)
    rates = []
    for node in gene.traverse():
        if node is gene:
            continue
        shift = case.get("copy_ratio", 1.0) if all(
            tip.name.endswith("_2") for tip in node.leaves()
        ) else 1.0
        rate = 0.01 * shift * math.exp(rng.normal(0, case["rate_sd"]))
        rates.append(dict(clade=sorted(t.name for t in node.leaves()),
                          duration=node.dist, rate=rate))
        node.dist *= rate
    gene.write(outfile=str(directory / "gene.nwk"), parser=1, format_root_node=True)
    mapping = {tip.name: tip.name.rsplit("_", 1)[0] for tip in gene.leaves()}
    for tip in sorted(mapping)[:case.get("mapping_errors", 0)]:
        mapping[tip] = f"S{(int(mapping[tip][1:]) + 1) % case['species']}"
    (directory / "mapping.tsv").write_text(
        "leaf_name\tspecies_label\n" + "".join(f"{k}\t{v}\n" for k, v in mapping.items())
    )
    save_json(directory / "truth.json", dict(age=truth, target="D", rates=rates,
                                             seed=seed, case=case))
    return truth


def read_target(prefix):
    manifest = json.loads(Path(str(prefix) + ".manifest.json").read_text())
    with Path(str(prefix) + ".nodes.tsv").open() as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t")
                if row["gene_name"] == "D"]
    if len(rows) != 1 or rows[0]["event_type"] != "duplication":
        return manifest, None
    return manifest, rows[0]


def finite_number(value):
    try:
        value = float(value)
        return value if math.isfinite(value) else None
    except (TypeError, ValueError):
        return None


def evaluate(case, family, options):
    directory = options.output / case["name"] / f"f{family:04d}"
    directory.mkdir(parents=True, exist_ok=False)
    seed = case["seed"] + family
    row = dict(case=case["name"], family=family, seed=seed, point_success=False,
               interval_available=False, covered=False, status="not-started")
    try:
        truth = generate_trees(directory, case, seed)
        row["truth"] = truth
        command = [options.iqtree, "--alisim", "alignment", "-t", "gene.nwk",
                   "-m", "GY{0.5,2}+FQ+G4{0.7}", "--length", str(3 * case["codons"]),
                   "--out-format", "fasta", "-seed", str(seed + 10000000), "-nt", "1"]
        row["generation_status"], row["generation_seconds"] = run(
            command, directory, "generate", options.timeout)
        if row["generation_status"] != "complete":
            row["status"] = "generation-" + row["generation_status"]
            return row
        # Independently check AliSim length semantics and standard-code stops.
        sequences = {}
        for line in (directory / "alignment.fa").read_text().splitlines():
            if line.startswith(">"):
                name = line[1:].split()[0]
                sequences[name] = ""
            else:
                sequences[name] += line.strip()
        if not sequences or any(len(seq) != 3 * case["codons"] for seq in sequences.values()):
            raise ValueError("AliSim did not return the requested number of codons")
        if any(seq[i:i+3] in {"TAA", "TAG", "TGA"}
               for seq in sequences.values() for i in range(0, len(seq), 3)):
            raise ValueError("AliSim generated a standard-code stop")
        row["input_sha256"] = {name: digest(directory / name) for name in
                               ["gene.nwk", "species.nwk", "mapping.tsv", "alignment.fa", "truth.json"]}
        common = [sys.executable, "-m", "nwkit", "radte", "--gene-tree", "gene.nwk",
                  "--species-tree", "species.nwk", "--species-map-tsv", "mapping.tsv",
                  "--reconcile", "lca", "--alignment", "alignment.fa",
                  "--backend", "native", "--sequence-engine", "native",
                  "--substitution-model", "gy94", "--gamma-categories", "4",
                  "--inference", "auto", "--likelihood", "auto", "--maxiter", "1000",
                  "--max-age", "1000", "--interval-level", "0.95",
                  "--seed", str(seed + 20000000)]
        for method in ["none", "profile"]:
            row[method + "_status"], row[method + "_seconds"] = run(
                common + ["--uncertainty", method, "--out-prefix", method],
                directory, method, options.timeout)
            if row[method + "_status"] != "complete":
                row["status"] = method + "-" + row[method + "_status"]
                break
            manifest, target = read_target(directory / method)
            row[method + "_estimator"] = manifest["method"]
            row[method + "_diagnostics"] = manifest["diagnostics"]
            if target is None:
                row["status"] = "target-unmatched"
                break
            row["estimate"] = finite_number(target["estimated_age"])
            row[method + "_estimate"] = row["estimate"]
            row["point_success"] = row["estimate"] is not None
            row["bias"] = row["estimate"] - truth
            if method == "profile":
                lower, upper = (finite_number(target[k]) for k in ["interval_lower", "interval_upper"])
                row.update(lower=lower, upper=upper, interval_status=target["interval_status"])
                row["interval_available"] = lower is not None and upper is not None and lower <= upper
                row["covered"] = row["interval_available"] and lower <= truth <= upper
                row["miss_side"] = ("below" if truth < lower else "above" if truth > upper else "covered") if row["interval_available"] else "unavailable"
                row["width"] = upper - lower if row["interval_available"] else None
                row["status"] = "complete"
    except (ValueError, OSError, KeyError) as exc:
        row["status"] = "evaluation-error"
        row["error"] = str(exc)
    finally:
        save_json(directory / "result.json", row)
    return row


def summarize(rows):
    n = len(rows)
    points = [r for r in rows if r["point_success"]]
    available = [r for r in rows if r["interval_available"]]
    covered = sum(r["covered"] for r in rows)
    result = dict(families=n, point_successes=len(points), intervals_returned=len(available),
                  truth_covered=covered, coverage_among_returned=covered / len(available) if available else None,
                  correct_return_fraction=covered / n, availability=len(available) / n,
                  bias=float(np.mean([r["bias"] for r in points])) if points else None,
                  rmse=float(np.sqrt(np.mean([r["bias"]**2 for r in points]))) if points else None,
                  median_width=float(np.median([r["width"] for r in available])) if available else None,
                  statuses=dict(Counter(r["status"] for r in rows)),
                  interval_statuses=dict(Counter(r.get("interval_status", r["status"]) for r in rows)),
                  estimators=dict(Counter(r.get("profile_estimator", r.get("none_estimator", "unavailable")) for r in rows)))
    for label, numerator, denominator in [("coverage", covered, len(available)),
                                         ("availability", len(available), n), ("correct_return", covered, n)]:
        ci = binomtest(numerator, denominator).proportion_ci(method="wilson") if denominator else None
        result[label + "_wilson95"] = [ci.low, ci.high] if ci else None
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--iqtree", default="iqtree3")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--timeout", type=float, default=300)
    options = parser.parse_args()
    options.output = options.output.resolve()
    options.iqtree = shutil.which(options.iqtree)
    if not options.iqtree or options.workers < 1 or options.timeout <= 0:
        parser.error("Require IQ-TREE, positive workers and timeout")
    protocol = json.loads(options.protocol.read_text())
    cases = protocol["cases"]
    if len({c["name"] for c in cases}) != len(cases) or any(
        c["families"] < 1 or c["species"] < 2 or c["codons"] < 1 for c in cases
    ):
        parser.error("Invalid or duplicate protocol cases")
    options.output.mkdir(parents=True, exist_ok=False)
    save_json(options.output / "protocol.json", protocol)
    metadata = dict(
        command=sys.argv, python=sys.version, platform=platform.platform(),
        protocol_sha256=digest(options.protocol), runner_sha256=digest(__file__),
        generator_sha256={p.name: digest(p) for p in [ROOT / "tools/radte_benchmark_cases.py",
                                                    ROOT / "tools/radte_interval_simulation.py"]},
        source_sha256={p.name: digest(p) for p in sorted((ROOT / "nwkit").glob("*.py"))},
        iqtree_sha256=digest(options.iqtree),
        iqtree_version=subprocess.check_output([options.iqtree, "--version"], text=True),
        threads={k: v for k, v in os.environ.items() if k.endswith("NUM_THREADS")},
        workers=options.workers, timeout_per_stage=options.timeout)
    save_json(options.output / "metadata.json", metadata)
    rows = []
    with ThreadPoolExecutor(max_workers=options.workers) as pool, (options.output / "cases.jsonl").open("w") as handle:
        futures = [pool.submit(evaluate, c, i, options) for c in cases for i in range(c["families"])]
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            handle.write(json.dumps(row, allow_nan=False) + "\n")
            handle.flush()
            print(f"{len(rows)}/{len(futures)} {row['case']} {row['status']}", flush=True)
    save_json(options.output / "summary.json", {
        c["name"]: summarize([r for r in rows if r["case"] == c["name"]]) for c in cases})
    metadata["changed_source_files_during_run"] = [name for name, value in metadata["source_sha256"].items()
                                                   if digest(ROOT / "nwkit" / name) != value]
    metadata["runner_changed_during_run"] = digest(__file__) != metadata["runner_sha256"]
    save_json(options.output / "metadata.json", metadata)


if __name__ == "__main__":
    main()
