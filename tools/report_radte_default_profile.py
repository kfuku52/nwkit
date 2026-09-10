"""Audit and summarize completed external default-profile validation runs."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ete4 import Tree
from scipy.stats import binomtest

LABELS = {
    "root": "Root duplication",
    "high_sd": "High rate variation",
    "internal": "Internal duplication",
    "long_alignment": "Longer alignment",
    "larger_family": "Larger family",
    "loss_rate_shift": "Loss + copy rate shift",
    "wrong_mapping": "Wrong species mapping",
    "low_sd": "Low rate variation",
}


def audit(root):
    protocol = json.loads((root / "protocol.json").read_text())
    summary = json.loads((root / "summary.json").read_text())
    rows = [
        json.loads(line) for line in (root / "cases.jsonl").read_text().splitlines()
    ]
    keys = [(r["case"], r["family"]) for r in rows]
    expected = {(c["name"], i) for c in protocol["cases"] for i in range(c["families"])}
    if len(set(keys)) != len(keys) or set(keys) != expected:
        raise ValueError("Missing/duplicated families in " + str(root))
    for row in rows:
        directory = root / row["case"] / f"f{row['family']:04d}"
        if json.loads((directory / "result.json").read_text()) != row:
            raise ValueError(
                "Per-family result disagrees with JSONL: " + str(directory)
            )
        if (directory / "truth.json").exists():
            truth = json.loads((directory / "truth.json").read_text())
            np.testing.assert_allclose(row["truth"], truth["age"])
            tree = Tree(str(directory / "truth.nwk"), parser=1)
            target = next(n for n in tree.traverse() if n.name == truth["target"])
            np.testing.assert_allclose(
                [tree.get_distance(target, tip) for tip in target.leaves()],
                row["truth"],
            )
        for name, value in row.get("input_sha256", {}).items():
            if hashlib.sha256((directory / name).read_bytes()).hexdigest() != value:
                raise ValueError("Input changed: " + str(directory / name))
    for case in protocol["cases"]:
        subset = [r for r in rows if r["case"] == case["name"]]
        result = summary[case["name"]]
        if (
            result["families"],
            result["intervals_returned"],
            result["truth_covered"],
        ) != (
            len(subset),
            sum(r["interval_available"] for r in subset),
            sum(r["covered"] for r in subset),
        ):
            raise ValueError("Summary denominators disagree: " + case["name"])
        for row in subset:
            if row["interval_available"]:
                if row["covered"] != (row["lower"] <= row["truth"] <= row["upper"]):
                    raise ValueError("Incorrect coverage indicator")
        points = [r for r in subset if r["point_success"]]
        available = [r for r in subset if r["interval_available"]]
        covered = sum(r["covered"] for r in subset)
        if result["point_successes"] != len(points):
            raise ValueError("Point denominator disagrees: " + case["name"])
        for row in subset:
            if row["covered"] and not row["interval_available"]:
                raise ValueError("Missing interval counted as covered")
        for row in points:
            np.testing.assert_allclose(row["bias"], row["estimate"] - row["truth"])
        errors = [r["estimate"] - r["truth"] for r in points]
        widths = [r["upper"] - r["lower"] for r in available]
        for field, value in {
            "bias": float(np.mean(errors)) if errors else None,
            "rmse": float(np.sqrt(np.mean(np.square(errors)))) if errors else None,
            "median_width": float(np.median(widths)) if widths else None,
        }.items():
            if value is None:
                if result[field] is not None:
                    raise ValueError("Missing metric reported: " + field)
            else:
                np.testing.assert_allclose(result[field], value)
        # Independent recomputation of all fractions and Wilson limits.
        for prefix, field, numerator, denominator in [
            ("coverage", "coverage_among_returned", covered, len(available)),
            ("availability", "availability", len(available), len(subset)),
            ("correct_return", "correct_return_fraction", covered, len(subset)),
        ]:
            if denominator:
                np.testing.assert_allclose(result[field], numerator / denominator)
                ci = binomtest(numerator, denominator).proportion_ci(method="wilson")
                np.testing.assert_allclose(
                    [ci.low, ci.high], result[prefix + "_wilson95"]
                )
            elif result[field] is not None or result[prefix + "_wilson95"] is not None:
                raise ValueError("Empty denominator has reported coverage")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    summaries = {}
    for root in args.results:
        current = audit(root)
        if summaries.keys() & current.keys():
            raise ValueError("Duplicate case names across studies")
        summaries.update(current)
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "summary.json").write_text(json.dumps(summaries, indent=2) + "\n")
    with (args.output / "summary.tsv").open("w") as handle:
        columns = [
            "case",
            "families",
            "point_successes",
            "intervals_returned",
            "truth_covered",
            "coverage_among_returned",
            "correct_return_fraction",
            "bias",
            "rmse",
            "median_width",
        ]
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        for name, result in summaries.items():
            writer.writerow(dict(case=name, **result))
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True, layout="constrained")
    positions = np.arange(len(summaries))
    labels = []
    for i, (name, result) in enumerate(summaries.items()):
        labels.append(f"{LABELS.get(name, name)} (n={result['families']})")
        for axis, value, ci in [
            (axes[0], result["coverage_among_returned"], result["coverage_wilson95"]),
            (
                axes[1],
                result["correct_return_fraction"],
                result["correct_return_wilson95"],
            ),
        ]:
            if value is not None:
                axis.errorbar(
                    100 * value,
                    i,
                    xerr=[[100 * (value - ci[0])], [100 * (ci[1] - value)]],
                    fmt="o",
                    color="#0072B2",
                    capsize=3,
                )
        axes[1].plot(
            100 * result["availability"], i, "|", color="#D55E00", markersize=13
        )
    axes[0].set_yticks(positions, labels)
    axes[0].invert_yaxis()
    for axis in axes:
        axis.set_xlim(0, 102)
        axis.axvline(95, color="#777777", linestyle="--", linewidth=1)
        axis.grid(axis="x", alpha=0.2)
        axis.set_xlabel("Percent (Wilson 95% intervals)")
        axis.spines[["right", "top"]].set_visible(False)
    axes[0].set_title("Truth covered / returned intervals")
    axes[1].set_title(
        "Correct interval returned / all trials\nOrange tick: interval availability"
    )
    fig.suptitle(
        "Native GY94 / nominal 95% profile: external simulation\nFixed species ages; one duplication target per independent family",
        fontsize=12,
    )
    fig.savefig(args.output / "coverage.png", dpi=180)
    fig.savefig(args.output / "coverage.svg")
    plt.close(fig)
    print(f"Audited {sum(s['families'] for s in summaries.values())} families")


if __name__ == "__main__":
    main()
