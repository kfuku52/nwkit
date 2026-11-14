"""Frozen independent branch-only validation of exact log-duration intervals.

This does not validate sequence intervals or general internal duplications.
The explicit six-edge covariance and dense QR reference are independent of
the native precision builder and interval implementation.
"""

import argparse
import copy
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from ete4 import Tree
from radte_interval_reference import gaussian_contrast_interval
from scipy.stats import binomtest

from nwkit.radte_exact_interval import exact_log_duration_intervals
from nwkit.radte_inputs import build_chronology
from nwkit.radte_model import fit_dates, laplace_intervals
from nwkit.radte_studentized import studentized_intervals
from nwkit.reconcile import build_reconciliation_table

DISTANCE = np.array(
    [
        [0, 2, 1, 1, 3, 3],
        [2, 0, 3, 3, 1, 1],
        [1, 3, 0, 2, 4, 4],
        [1, 3, 2, 0, 4, 4],
        [3, 1, 4, 4, 0, 2],
        [3, 1, 4, 4, 2, 0],
    ]
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    args = parser.parse_args()
    protocol = json.loads(args.protocol.read_text())
    if any(
        protocol.get(key) != expected
        for key, expected in [
            ("level", 0.95),
            ("truth_age", 20),
            ("species_age", 10),
            ("max_age", 100),
        ]
    ):
        parser.error(
            "This exact-reference design requires level=.95, truth_age=20, species_age=10, max_age=100."
        )
    if not protocol.get("cells") or any(
        not isinstance(cell.get("families"), int)
        or cell["families"] < 1
        or not 0 <= cell["rho"] < 1
        or not np.isfinite(cell["sd"])
        or cell["sd"] <= 0
        for cell in protocol["cells"]
    ):
        parser.error("Each cell requires positive families/SD and 0 <= rho < 1.")
    args.output.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    sources = list((root / "nwkit").glob("radte*.py")) + [
        Path(__file__),
        root / "tools/radte_interval_reference.py",
    ]
    metadata = dict(
        protocol=protocol,
        protocol_sha256=hashlib.sha256(args.protocol.read_bytes()).hexdigest(),
        source_sha256={
            str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sources
        },
    )
    (args.output / "metadata.json").write_text(json.dumps(metadata, indent=2))
    rows = []
    for cell in protocol["cells"]:
        gene = Tree("((A_1:.1,B_1:.1)S1:.1,(A_2:.1,B_2:.1)S2:.1)D;", parser=1)
        species = Tree("(A:10,B:10)AB;", parser=1)
        table = build_reconciliation_table(
            gene, species, {n.name: n.name.split("_")[0] for n in gene.leaves()}
        )
        c = build_chronology(gene, species, table, max_age=100)
        covariance = cell["rho"] ** DISTANCE
        factor = np.linalg.cholesky(covariance)
        design = np.column_stack([np.ones(6), [1, 1, 0, 0, 0, 0]])
        rng = np.random.default_rng(cell["seed"])
        for family in range(cell["families"]):
            y = np.log(0.1) + cell["sd"] * (factor @ rng.normal(size=6))
            for node, length in zip(c.edges, np.exp(y), strict=True):
                node.dist = float(length)
            reference = gaussian_contrast_interval(y, design, covariance, [0, 1])
            expected = np.array(
                [
                    max(1 + c.min_duration, 1 + np.exp(reference["lower"])),
                    min(10, 1 + np.exp(reference["upper"])),
                ]
            )
            try:
                fit, problem = fit_dates(
                    c, rho=cell["rho"], starts=1, seed=family + cell["seed"]
                )
                for method, evaluate in [
                    ("exact", exact_log_duration_intervals),
                    ("laplace", laplace_intervals),
                    ("studentized", studentized_intervals),
                ]:
                    result = copy.deepcopy(fit)
                    evaluate(result, problem)
                    available = result.interval_lower is not None
                    group = problem.free[0]
                    endpoints = (
                        np.array(
                            [result.interval_lower[group], result.interval_upper[group]]
                        )
                        if available
                        else None
                    )
                    error = (
                        float(np.max(abs(endpoints - expected)))
                        if method == "exact" and available
                        else None
                    )
                    rows.append(
                        dict(
                            cell=cell["name"],
                            family=family,
                            method=method,
                            available=available,
                            covered=bool(
                                available and endpoints[0] <= 2 <= endpoints[1]
                            ),
                            lower=float(endpoints[0] * c.scale) if available else None,
                            upper=float(endpoints[1] * c.scale) if available else None,
                            status=result.interval_status,
                            reference_error=error,
                        )
                    )
            except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
                for method in ["exact", "laplace", "studentized"]:
                    rows.append(
                        dict(
                            cell=cell["name"],
                            family=family,
                            method=method,
                            available=False,
                            covered=False,
                            lower=None,
                            upper=None,
                            status=str(exc),
                            reference_error=None,
                        )
                    )
        print(cell["name"], flush=True)
    with (args.output / "cases.csv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary = {}
    for cell in protocol["cells"]:
        for method in ["exact", "laplace", "studentized"]:
            selected = [
                r for r in rows if r["cell"] == cell["name"] and r["method"] == method
            ]
            n = len(selected)
            returned = sum(r["available"] for r in selected)
            covered = sum(r["covered"] for r in selected)
            errors = [
                r["reference_error"]
                for r in selected
                if r["reference_error"] is not None
            ]
            lower = (
                binomtest(covered, n, alternative="greater")
                .proportion_ci(
                    confidence_level=1 - 0.05 / len(protocol["cells"]), method="exact"
                )
                .low
            )
            summary[cell["name"] + ":" + method] = dict(
                families=n,
                returned=returned,
                covered=covered,
                availability=returned / n,
                conditional_coverage=covered / returned if returned else None,
                correct_return=covered / n,
                simultaneous_correct_return_lower=float(lower),
                max_reference_error=max(errors) if errors else None,
            )
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
