"""Measure constrained marginal-problem construction in an installed runtime.

Mount this script and benchmark_radte_iqtree.py into each baseline/candidate
runtime. The synthetic quadratic isolates profile setup from sequence fitting.
Parent construction and evaluation are outside the timed loop.
"""

import argparse
import inspect
import json
import tempfile
import time
from dataclasses import replace
from pathlib import Path

import numpy as np
from benchmark_radte_iqtree import fixture

from nwkit.radte_marginal import MarginalDatingProblem
from nwkit.radte_sequence import QuadraticLikelihood


def quadratic(chronology):
    roots = {i for i, node in enumerate(chronology.edges) if node.up is chronology.gene}
    mapping, next_row = [], 1
    for i in range(len(chronology.edges)):
        mapping.append(0 if i in roots else next_row)
        next_row += i not in roots
    mapping = np.array(mapping)
    center = np.log(
        np.bincount(mapping, weights=[node.dist for node in chronology.edges])
    )
    return QuadraticLikelihood(
        center, np.zeros(next_row), np.eye(next_row) * 40 + 0.5, 20, mapping
    )


def trial(chronology, likelihood, rho, constructors):
    parent = MarginalDatingProblem(
        chronology, likelihood=likelihood, rho=rho, rate_sd=0.4
    )
    extra = (
        {"shared_structure": parent.shared_structure}
        if "shared_structure" in inspect.signature(MarginalDatingProblem).parameters
        else {}
    )
    group = int(np.flatnonzero(chronology.upper > chronology.lower)[0])
    lower, upper = chronology.lower.copy(), chronology.upper.copy()
    lower[group] = upper[group] = chronology.initial[group]
    constrained = replace(chronology, lower=lower, upper=upper)
    started = time.perf_counter()
    for _ in range(constructors):
        problem = MarginalDatingProblem(
            constrained, likelihood=likelihood, rho=rho, rate_sd=0.4, **extra
        )
    elapsed = time.perf_counter() - started
    value, gradient = problem.value_gradient(problem.initial_parameters())
    return dict(
        rho=rho,
        constructors=constructors,
        shared=bool(extra),
        seconds=elapsed,
        nll=value,
        gradient=gradient.tolist(),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--constructors", type=int, default=40)
    args = parser.parse_args()
    if args.repeats < 1 or args.constructors < 1:
        parser.error("Counts must be positive.")
    rows = []
    for tips in [4, 64, 256]:
        with tempfile.TemporaryDirectory() as temporary:
            chronology, _ = fixture(Path(temporary), tips, 1)
            likelihood = quadratic(chronology)
            for rho in [0.0, 0.5]:
                for repeat in range(args.repeats):
                    row = trial(chronology, likelihood, rho, args.constructors)
                    row.update(tips=tips, repeat=repeat)
                    rows.append(row)
                    print(json.dumps({k: v for k, v in row.items() if k != "gradient"}))
    args.output.write_text(json.dumps(rows, indent=2) + "\n")


if __name__ == "__main__":
    main()
