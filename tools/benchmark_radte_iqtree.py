"""Measure IQ-TREE 3 setup and uncached branch evaluations for one interface."""

import argparse
import hashlib
import json
import os
import platform
import resource
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
from ete4 import Tree
from scipy.linalg import expm

from nwkit.radte_codon import CODONS, codon_matrix
from nwkit.radte_inputs import build_chronology
from nwkit.radte_iqtree import IQTreeLikelihood
from nwkit.reconcile import build_reconciliation_table


def fixture(directory, tips, sites):
    names = [f"{species}_{i}" for i in range(tips // 2) for species in ["A", "B"]]

    def balanced(labels):
        if len(labels) == 1:
            return labels[0] + ":0.1"
        middle = len(labels) // 2
        return (
            "(" + balanced(labels[:middle]) + "," + balanced(labels[middle:]) + "):0.1"
        )

    gene = Tree(balanced(names).rsplit(":", 1)[0] + ";", parser=1)
    species = Tree("(A:10,B:10)AB;", parser=1)
    table = build_reconciliation_table(
        gene, species, {n: n.split("_")[0] for n in names}
    )
    chronology = build_chronology(gene, species, table, None, 100)
    q, pi = codon_matrix(
        "gy94", np.ones((tips, 1), dtype=np.uint64), np.ones(1), 2, 0.5, "fq"
    )
    transition = np.maximum(expm(q * 0.1), 0)
    transition /= transition.sum(axis=1, keepdims=True)
    rng = np.random.default_rng(918)
    sequences = {gene: rng.choice(len(CODONS), sites, p=pi)}
    for node in gene.traverse("preorder"):
        if node is not gene:
            sequences[node] = np.array(
                [rng.choice(len(CODONS), p=transition[i]) for i in sequences[node.up]]
            )
    path = directory / "alignment.fa"
    path.write_text(
        "".join(
            ">" + n.name + "\n" + "".join(CODONS[i] for i in sequences[n]) + "\n"
            for n in gene.leaves()
        )
    )
    return chronology, path


def trial(args):
    with tempfile.TemporaryDirectory() as temporary:
        c, path = fixture(Path(temporary), args.tips, args.sites)
        started = time.perf_counter()
        exact = IQTreeLikelihood(
            c,
            path,
            "GY{0.5,2}+FQ+G4{1}",
            executable=args.executable,
            interface=args.interface,
            worker=args.worker,
            threads=args.threads,
        )
        setup = time.perf_counter() - started
        try:

            def evaluate(lengths):
                if args.derivative_order == 1:
                    return exact.value_gradient(lengths)
                _, nll, gradient, _, mapping, _ = exact.evaluate(lengths)
                return nll, gradient[mapping]

            base = np.array([n.dist for n in c.edges])
            evaluate(base)
            probes = [
                base * np.exp(np.sin(np.arange(len(base)) + i) * 0.08)
                for i in range(1, args.evaluations + 1)
            ]
            started = time.perf_counter()
            values = [evaluate(lengths) for lengths in probes]
            elapsed = time.perf_counter() - started
            result = dict(
                interface=exact.interface,
                library=exact.worker_info,
                tips=args.tips,
                codons=args.sites,
                setup_s=setup,
                evaluations=args.evaluations,
                derivative_order_requested=args.derivative_order,
                threads=args.threads,
                evaluation_s=elapsed,
                nll=[v[0] for v in values],
                gradient=[v[1].tolist() for v in values],
            )
        finally:
            exact.close()
        # macOS reports bytes; Linux reports KiB.
        rss = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
        result["peak_iqtree_rss_kib"] = rss / 1024 if sys.platform == "darwin" else rss
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", default="iqtree3")
    parser.add_argument(
        "--interface", choices=["cli", "library", "auto"], default="cli"
    )
    parser.add_argument("--worker")
    parser.add_argument("--output", type=Path, default=Path("iqtree-benchmark.json"))
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--evaluations", type=int, default=20)
    parser.add_argument("--derivative-order", type=int, choices=[1, 2], default=1)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--trial", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--tips", type=int, default=4)
    parser.add_argument("--sites", type=int, default=150)
    args = parser.parse_args()
    if args.evaluations < 1 or args.threads < 1:
        parser.error("Evaluation and thread counts must be positive.")
    if args.trial:
        print(json.dumps(trial(args)))
        return
    executable = shutil.which(args.executable)
    if executable is None or args.repeats < 2:
        parser.error("An IQ-TREE executable and at least two repetitions are required.")
    report = dict(
        environment=dict(
            platform=platform.platform(),
            python=sys.version,
            iqtree=subprocess.check_output([executable, "--version"], text=True),
            executable_sha256=hashlib.sha256(Path(executable).read_bytes()).hexdigest(),
            threads=args.threads,
        ),
        trials=[],
    )
    environment = {**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"}
    for tips, sites in [(4, 150), (16, 150), (16, 1500), (64, 1500)]:
        reference = None
        for repeat in range(args.repeats):
            command = [
                sys.executable,
                str(Path(__file__).resolve()),
                "--trial",
                "--executable",
                executable,
                "--interface",
                args.interface,
                "--tips",
                str(tips),
                "--sites",
                str(sites),
                "--evaluations",
                str(args.evaluations),
                "--derivative-order",
                str(args.derivative_order),
                "--threads",
                str(args.threads),
            ]
            if args.worker:
                command += ["--worker", args.worker]
            data = json.loads(
                subprocess.check_output(command, text=True, env=environment)
            )
            data["repeat"] = repeat
            if reference is None:
                reference = data
            np.testing.assert_allclose(data["nll"], reference["nll"], rtol=0, atol=2e-5)
            np.testing.assert_allclose(
                data["gradient"], reference["gradient"], rtol=1e-5, atol=2e-4
            )
            data["max_nll_difference"] = float(
                np.max(np.abs(np.array(data["nll"]) - reference["nll"]))
            )
            report["trials"].append(data)
            args.output.write_text(json.dumps(report, indent=2) + "\n")
            print(
                json.dumps(
                    {k: v for k, v in data.items() if k not in {"nll", "gradient"}}
                ),
                flush=True,
            )


if __name__ == "__main__":
    main()
