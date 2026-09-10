#!/usr/bin/env python3
"""Plot recorded calibration estimates and Monte Carlo intervals, without fitting."""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def load(path):
    return {(row["case"], row["method"]): row for row in json.loads(path.read_text())}


def points(axis, source, cases, method, metric, label, color, offset):
    xs, ys, lower, upper = [], [], [], []
    for index, case in enumerate(cases):
        row = source[(case, method)][metric]
        if row is None or row["estimate"] is None:
            continue
        xs.append(index + offset)
        ys.append(row["estimate"])
        lower.append(max(0, row["estimate"] - row["mc95_lower"]))
        upper.append(max(0, row["mc95_upper"] - row["estimate"]))
    axis.errorbar(
        xs,
        ys,
        yerr=[lower, upper],
        fmt="o",
        ms=5,
        capsize=3,
        lw=1.3,
        color=color,
        label=label,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    bootstrap = load(args.results / "bootstrap-pilot-20260910" / "summary.json")
    before = load(args.results / "wald-pilot-20260910" / "summary.json")
    after = load(args.results / "glmm-fixed-pilot-20260910" / "summary.json")
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    rsc_cases = ["rsc-e2", "rsc-e5", "rsc-e20"]
    methods = [
        ("wald", "Wald t", "#0072B2", -0.2),
        ("parametric-bootstrap", "Coefficient bootstrap", "#D55E00", 0),
        ("null-bootstrap", "Null bootstrap (test only)", "#009E73", 0.2),
    ]
    for axis, metric, target, title in [
        (axes[0, 0], "rejection_given_available", 0.05, "RSC: null rejection rate"),
        (axes[0, 1], "coverage_given_interval", 0.95, "RSC: 95% interval coverage"),
    ]:
        for method, label, color, offset in methods:
            points(axis, bootstrap, rsc_cases, method, metric, label, color, offset)
        axis.axhline(target, ls="--", lw=1, color="0.4")
        axis.set(
            xticks=range(3),
            xticklabels=["2 events", "5 events", "20 events"],
            title=title,
        )
        axis.set_ylim((0, 0.38) if target == 0.05 else (0.6, 1.01))
    axes[0, 0].legend(fontsize=8, loc="upper right")
    glmm_cases = [
        "binomial-n8-p0.05",
        "binomial-n30-p0.5",
        "poisson-n30",
        "negative-binomial-n30",
    ]
    labels = [
        "Binary n=8\nrare baseline",
        "Binary n=30\nbalanced baseline",
        "Poisson n=30",
        "NB2 n=30",
    ]
    for axis, metric, title in [
        (axes[1, 0], "p_available", "GLMM: fraction with available P value"),
        (
            axes[1, 1],
            "rejection_given_available",
            "GLMM: null rejection, available fits only",
        ),
    ]:
        points(
            axis,
            before,
            glmm_cases,
            "wald",
            metric,
            "Before multistart fix",
            "#777777",
            -0.09,
        )
        points(
            axis,
            after,
            glmm_cases,
            "wald",
            metric,
            "After multistart fix",
            "#0072B2",
            0.09,
        )
        axis.set(xticks=range(4), xticklabels=labels, title=title)
        axis.set_ylim((0, 1.02) if metric == "p_available" else (0, 0.22))
    axes[1, 1].axhline(0.05, ls="--", lw=1, color="0.4")
    axes[1, 0].legend(fontsize=8, loc="lower right")
    for axis in axes.flat:
        axis.grid(axis="y", alpha=0.2)
        axis.spines[["top", "right"]].set_visible(False)
        axis.tick_params(labelsize=9)
    fig.suptitle(
        "Fixed-model calibration pilot · 200 generated datasets per condition\n"
        "Error bars: Wilson 95% Monte Carlo intervals; RSC bootstrap B=999",
        fontsize=13,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for suffix in ("png", "pdf"):
        target = args.output.with_suffix("." + suffix)
        if target.exists():
            raise FileExistsError(target)
        fig.savefig(target, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
