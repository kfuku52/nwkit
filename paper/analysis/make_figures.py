#!/usr/bin/env python3
"""Generate the three main manuscript figures from versioned analysis outputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd
from ete4 import Tree

from nwkit.util import assign_branch_ids


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS_DIR = PROJECT_ROOT / "paper" / "results"
DEFAULT_FIGURES_DIR = PROJECT_ROOT / "paper" / "figures"

C3 = "#0072B2"
C4 = "#E69F00"
GREEN = "#009E73"
PURPLE = "#7B61A8"
GRAY = "#666666"
LIGHT_GRAY = "#E7E7E7"
INK = "#222222"


def save_figure(figure: plt.Figure, output_stem: Path) -> None:
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(output_stem.with_suffix(".svg"), bbox_inches="tight")
    figure.savefig(output_stem.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(figure)


def rounded_box(axis, x, y, width, height, text, facecolor, fontsize=8, edgecolor="none"):
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        facecolor=facecolor,
        edgecolor=edgecolor,
        linewidth=0.8,
    )
    axis.add_patch(patch)
    axis.text(x + width / 2, y + height / 2, text, ha="center", va="center", fontsize=fontsize)
    return patch


def figure1_architecture(figures_dir: Path) -> None:
    figure, axis = plt.subplots(figsize=(7.2, 6.4))
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    rounded_box(
        axis,
        0.18,
        0.81,
        0.64,
        0.065,
        "Newick tree(s)  +  taxon/trait tables  +  reference trees",
        "#DCEAF3",
        fontsize=9,
    )
    axis.text(
        0.5,
        0.91,
        "33 stdin defaults  ·  31 stdout defaults  ·  612 tests",
        ha="center",
        va="center",
        fontsize=7.5,
        color=GRAY,
    )
    axis.add_patch(FancyArrowPatch((0.5, 0.81), (0.5, 0.775), arrowstyle="-|>", mutation_scale=12, color=GRAY))
    rounded_box(
        axis,
        0.07,
        0.71,
        0.86,
        0.065,
        "Shared inputs: stdin or path · explicit/automatic Newick formats · quoted-name control",
        "#F2F2F2",
        fontsize=8,
        edgecolor="#C8C8C8",
    )

    families = [
        (
            0.06,
            0.48,
            0.27,
            0.19,
            "Inspect & normalize  (7)\n\nvalidate · info · sanitize\nnwk2table · table2nwk\nnhx2nwk · printlabel",
            "#E1F3EE",
        ),
        (
            0.365,
            0.48,
            0.27,
            0.19,
            "Transform  (13)\n\ncollapse · compose · drop\nintersection · label · mark\nprune · rename · rescale\nroot · shuffle · subtree\ntransfer",
            "#FFF0DB",
        ),
        (
            0.67,
            0.48,
            0.27,
            0.19,
            "Compare & summarize  (4)\n\ndiff · dist\ncladefreq · consensus",
            "#E8E2F3",
        ),
        (
            0.205,
            0.25,
            0.27,
            0.17,
            "Use taxonomy & traits  (5)\n\nannotate · constrain · monophyly\nsample · skim",
            "#F8E1E3",
        ),
        (
            0.525,
            0.25,
            0.27,
            0.17,
            "Analyze & present  (4)\n\nasr · mcmctree\ndraw · image",
            "#E5EDF8",
        ),
    ]
    for x, y, width, height, text, color in families:
        rounded_box(axis, x, y, width, height, text, color, fontsize=7.5, edgecolor="#B8B8B8")
    axis.add_patch(FancyArrowPatch((0.5, 0.71), (0.5, 0.675), arrowstyle="-|>", mutation_scale=9, color="#999999", linewidth=0.7))

    rounded_box(
        axis,
        0.07,
        0.12,
        0.86,
        0.065,
        "Shared outputs: stdout or path · documented tables · seeds · JSONL provenance audit",
        "#F2F2F2",
        fontsize=7.7,
        edgecolor="#C8C8C8",
    )
    axis.add_patch(FancyArrowPatch((0.5, 0.25), (0.5, 0.185), arrowstyle="-|>", mutation_scale=9, color="#999999", linewidth=0.7))
    rounded_box(
        axis,
        0.11,
        0.025,
        0.78,
        0.055,
        "Newick · TSV · annotated trees · audit records · figures · image manifests",
        "#DCEAF3",
        fontsize=8,
    )
    axis.add_patch(FancyArrowPatch((0.5, 0.12), (0.5, 0.08), arrowstyle="-|>", mutation_scale=10, color=GRAY))

    save_figure(figure, figures_dir / "figure1_architecture")


def node_coordinates(tree: Tree) -> tuple[dict[Tree, float], dict[Tree, float]]:
    x_by_node: dict[Tree, float] = {tree: 0.0}
    for node in tree.traverse(strategy="preorder"):
        for child in node.get_children():
            x_by_node[child] = x_by_node[node] + float(child.dist or 0.0)
    leaves = list(tree.leaves())
    y_by_node: dict[Tree, float] = {leaf: float(len(leaves) - index - 1) for index, leaf in enumerate(leaves)}
    for node in tree.traverse(strategy="postorder"):
        if node.is_leaf:
            continue
        children = list(node.get_children())
        y_by_node[node] = sum(y_by_node[child] for child in children) / len(children)
    return x_by_node, y_by_node


def probability_color(probability: float) -> tuple[float, float, float]:
    left = np.asarray(mcolors.to_rgb(C3))
    right = np.asarray(mcolors.to_rgb(C4))
    return tuple(left * (1.0 - probability) + right * probability)


def selected_tip_label(label: str) -> str:
    parts = label.replace("'", "").split("_")
    if len(parts) < 3:
        return label
    abbreviated_species = rf"{parts[0][0]}.\ {parts[1]}"
    return rf"$\mathit{{{abbreviated_species}}}$ {parts[-1]}"


def figure2_pepc(results_dir: Path, figures_dir: Path) -> None:
    example_dir = results_dir / "pepc"
    tree_path = PROJECT_ROOT / "paper" / "data" / "pepc" / "PEPC.tree.nwk"
    tree = Tree(tree_path.read_text(), parser=0)
    node_to_branch = assign_branch_ids(tree)
    asr = pd.read_csv(example_dir / "pepc_c4_asr.tsv", sep="\t").set_index("branch_id")
    traits = pd.read_csv(example_dir / "pepc_traits.tsv", sep="\t").set_index("leaf_name")
    representatives = set(pd.read_csv(example_dir / "pepc_c4_representatives.tsv", sep="\t")["leaf_name"])
    summary = pd.read_csv(example_dir / "pepc_summary.tsv", sep="\t").iloc[0]
    x_by_node, y_by_node = node_coordinates(tree)

    figure = plt.figure(figsize=(7.2, 5.2))
    figure.subplots_adjust(top=0.96, bottom=0.08)
    grid = figure.add_gridspec(1, 2, width_ratios=(4.9, 1.6), wspace=0.08)
    tree_axis = figure.add_subplot(grid[0, 0])
    summary_axis = figure.add_subplot(grid[0, 1])
    tree_axis.set_title(
        "a)  Marginal P(C4) and observed tip states",
        loc="left",
        fontsize=9,
        pad=8,
        fontweight="bold",
    )

    for node in tree.traverse(strategy="preorder"):
        children = list(node.get_children())
        if children:
            tree_axis.plot(
                [x_by_node[node], x_by_node[node]],
                [min(y_by_node[child] for child in children), max(y_by_node[child] for child in children)],
                color="#A7A7A7",
                linewidth=0.5,
                zorder=1,
            )
        for child in children:
            probability = float(asr.loc[node_to_branch[child], "p_1"])
            tree_axis.plot(
                [x_by_node[node], x_by_node[child]],
                [y_by_node[child], y_by_node[child]],
                color=probability_color(probability),
                linewidth=1.05,
                solid_capstyle="round",
                zorder=2,
            )

    maximum_x = max(x_by_node.values())
    label_x = maximum_x * 1.025
    for leaf in tree.leaves():
        state = int(traits.loc[leaf.name, "C4"])
        selected = leaf.name in representatives
        tree_axis.scatter(
            x_by_node[leaf],
            y_by_node[leaf],
            s=20 if selected else 10,
            marker="^" if state == 1 else "o",
            color=C4 if state == 1 else C3,
            edgecolor="black" if selected else "white",
            linewidth=0.65 if selected else 0.2,
            zorder=4,
        )
        if selected:
            tree_axis.plot([x_by_node[leaf], label_x], [y_by_node[leaf], y_by_node[leaf]], color="#BDBDBD", linewidth=0.35, zorder=0)
            tree_axis.text(
                label_x,
                y_by_node[leaf],
                selected_tip_label(leaf.name),
                va="center",
                fontsize=6.5,
            )

    tree_axis.set_xlim(0, maximum_x * 1.55)
    tree_axis.set_ylim(-1.5, len(list(tree.leaves())) + 0.5)
    tree_axis.set_yticks([])
    tree_axis.set_xlabel("Expected substitutions per site", fontsize=8.2)
    tree_axis.tick_params(axis="x", labelsize=7.5)
    for spine in ("left", "right", "top"):
        tree_axis.spines[spine].set_visible(False)

    cmap = mcolors.LinearSegmentedColormap.from_list("c4_probability", [C3, C4])
    colorbar_axis = tree_axis.inset_axes([0.02, 0.045, 0.28, 0.013])
    matplotlib.colorbar.ColorbarBase(colorbar_axis, cmap=cmap, norm=mcolors.Normalize(0, 1), orientation="horizontal")
    colorbar_axis.tick_params(labelsize=6, length=2)
    colorbar_axis.set_title("marginal P(C4)", fontsize=6.8, pad=2)
    tree_axis.scatter([], [], color=C3, marker="o", s=16, label="observed C3")
    tree_axis.scatter([], [], color=C4, marker="^", s=20, label="observed C4")
    tree_axis.scatter(
        [],
        [],
        facecolors="none",
        edgecolor="black",
        marker="^",
        s=22,
        label="max-PD representative",
    )
    tree_axis.legend(loc="lower right", frameon=False, fontsize=6.7, handletextpad=0.4, borderpad=0.2)

    summary_axis.axis("off")
    summary_axis.set_title("b)  Workflow readout", loc="left", fontsize=9, pad=8, fontweight="bold")
    steps = [
        ("Preflight", "71 unique tips\nrooted and binary", "#E1F3EE"),
        ("Trait diagnosis", "21 C4 tips\nC4 polyphyletic", "#FFF0DB"),
        ("Mk (ER)", f"root P(C4) = {summary.root_p_c4:.3f}\n{int(summary.map_state_shift_edges)} MAP-shift edges", "#E8E2F3"),
        (
            "Max-PD sample",
            f"8 representatives\n{summary.representative_pd_coverage * 100:.1f}% of C4-tip PD",
            "#F8E1E3",
        ),
    ]
    y_positions = [0.83, 0.63, 0.43, 0.23]
    for index, ((heading, body, color), y_position) in enumerate(zip(steps, y_positions)):
        rounded_box(summary_axis, 0.05, y_position, 0.9, 0.13, "", color, fontsize=7, edgecolor="#B8B8B8")
        summary_axis.text(0.11, y_position + 0.095, heading, fontsize=8, fontweight="bold", va="center")
        summary_axis.text(0.11, y_position + 0.045, body, fontsize=7.5, va="center", linespacing=1.25)
        if index < len(steps) - 1:
            summary_axis.add_patch(FancyArrowPatch((0.5, y_position), (0.5, y_positions[index + 1] + 0.13), arrowstyle="-|>", mutation_scale=8, color="#888888"))
    summary_axis.add_patch(Rectangle((0.1, 0.12), 0.8, 0.025, facecolor=LIGHT_GRAY, edgecolor="none"))
    summary_axis.add_patch(Rectangle((0.1, 0.12), 0.8 * float(summary.representative_pd_coverage), 0.025, facecolor=GREEN, edgecolor="none"))
    summary_axis.text(0.1, 0.092, "union of root-to-tip branch length retained", fontsize=6.4, color=GRAY)
    summary_axis.text(
        0.05,
        0.015,
        "Data: CSUBST PEPC example\ncommit 1b0d98eec717",
        fontsize=6.4,
        color=GRAY,
        va="bottom",
    )
    save_figure(figure, figures_dir / "figure2_pepc_workflow")


def median_range(group: pd.DataFrame, column: str) -> tuple[float, float, float]:
    values = group[column].astype(float).to_numpy()
    return float(np.median(values)), float(np.min(values)), float(np.max(values))


def figure3_evaluation(results_dir: Path, figures_dir: Path) -> None:
    correctness = pd.read_csv(results_dir / "correctness_checks.tsv", sep="\t")
    input_cases = pd.read_csv(results_dir / "input_case_checks.tsv", sep="\t")
    benchmark = pd.read_csv(results_dir / "benchmark.tsv", sep="\t")
    benchmark = benchmark.loc[benchmark["returncode"] == 0].copy()

    figure, axes = plt.subplots(2, 2, figsize=(7.2, 6.5))
    figure.subplots_adjust(left=0.10, right=0.98, bottom=0.10, top=0.96, hspace=0.47, wspace=0.34)

    axis = axes[0, 0]
    check_labels = [
        "Rooted RF\nvs DendroPy",
        "Clade frequency\nvs direct count",
        "Majority consensus\nvs direct count",
        "Mk marginals\nvs enumeration",
        "Seeded shuffle\nreproduction",
    ]
    y_positions = np.arange(len(correctness.index))[::-1]
    values = correctness["agreement"].to_numpy() * 100
    axis.scatter(values, y_positions, s=34, color=GREEN, zorder=3)
    axis.hlines(y_positions, 0, values, color="#B8DCD2", linewidth=2)
    for y, row in zip(y_positions, correctness.itertuples()):
        axis.text(96, y + 0.17, f"{int(row.passed)}/{int(row.n)}", ha="right", va="bottom", fontsize=7, color=INK)
    axis.set_yticks(y_positions, check_labels, fontsize=7.2)
    axis.set_xlim(0, 105)
    axis.set_xlabel("Agreement (%)", fontsize=8.2)
    axis.set_title("a)  Predefined output checks", loc="left", fontsize=9, fontweight="bold")
    axis.grid(axis="x", color=LIGHT_GRAY, linewidth=0.6)

    axis = axes[0, 1]
    compact_labels = [
        "valid",
        "duplicate tip",
        "empty tip",
        "negative length",
        "singleton",
        "multifurcation",
        "missing support",
        "non-ultrametric",
        "unrooted",
        "ambiguous label",
        "quoted name",
        "leaf-set mismatch",
    ]
    for index, (label, matched) in enumerate(
        zip(compact_labels, input_cases["outcome_matched"])
    ):
        row = 3 - index // 3
        column = index % 3
        color = GREEN if bool(matched) else C4
        axis.scatter(column, row, s=46, color=color)
        axis.text(column, row - 0.24, label, ha="center", va="top", fontsize=7)
    axis.set_xlim(-0.5, 2.5)
    axis.set_ylim(-0.6, 3.55)
    axis.axis("off")
    axis.set_title("b)  Predefined input cases", loc="left", fontsize=9, fontweight="bold")
    matched_count = int(input_cases["outcome_matched"].sum())
    axis.text(
        2.45,
        3.48,
        f"{matched_count}/{len(input_cases.index)} outcomes matched",
        ha="right",
        va="top",
        fontsize=7.5,
        color=GREEN,
    )

    axis = axes[1, 0]
    validation = benchmark.loc[benchmark["benchmark"] == "validation_tip_scaling"]
    for shape, color, marker in (("balanced", C3, "o"), ("caterpillar", C4, "s")):
        subset = validation.loc[validation["shape"] == shape]
        x_values = []
        medians = []
        lower = []
        upper = []
        for num_tips, group in subset.groupby("num_tips"):
            median, minimum, maximum = median_range(group, "elapsed_s")
            x_values.append(num_tips)
            medians.append(median)
            lower.append(median - minimum)
            upper.append(maximum - median)
        axis.errorbar(x_values, medians, yerr=[lower, upper], marker=marker, color=color, linewidth=1.3, markersize=4, capsize=2, label=shape)
    axis.set_xscale("log", base=2)
    axis.set_xlabel("Tips in one tree", fontsize=8.2)
    axis.set_ylabel("Wall time (s)", fontsize=8.2)
    axis.set_title("c)  NWKIT preflight-report scaling", loc="left", fontsize=9, fontweight="bold")
    axis.legend(frameon=False, fontsize=7.2)
    axis.grid(color=LIGHT_GRAY, linewidth=0.6)
    axis.tick_params(labelsize=7.2)
    axis.set_ylim(0, max(validation["elapsed_s"]) * 1.08)

    axis = axes[1, 1]
    consensus = benchmark.loc[benchmark["benchmark"] == "consensus_tree_count_scaling"]
    style = {
        "NWKIT": (C3, "o"),
        "Gotree": (GREEN, "s"),
        "PhyKIT": (PURPLE, "^"),
    }
    for software in ("NWKIT", "Gotree", "PhyKIT"):
        subset = consensus.loc[consensus["software"] == software]
        if subset.empty:
            continue
        x_values = []
        medians = []
        lower = []
        upper = []
        for num_trees, group in subset.groupby("num_trees"):
            median, minimum, maximum = median_range(group, "elapsed_s")
            x_values.append(num_trees)
            medians.append(median)
            lower.append(median - minimum)
            upper.append(maximum - median)
        color, marker = style[software]
        axis.errorbar(x_values, medians, yerr=[lower, upper], marker=marker, color=color, linewidth=1.3, markersize=4, capsize=2, label=software)
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("Trees (256 tips each)", fontsize=8.2)
    axis.set_ylabel("Wall time (s; log scale)", fontsize=8.2)
    axis.set_title("d)  Majority-consensus CLI runs", loc="left", fontsize=9, fontweight="bold")
    axis.legend(frameon=False, fontsize=7.2)
    axis.grid(color=LIGHT_GRAY, linewidth=0.6, which="both")
    axis.tick_params(labelsize=7.2)

    for axis in axes.flat:
        for spine in ("top", "right"):
            axis.spines[spine].set_visible(False)
    save_figure(figure, figures_dir / "figure3_evaluation")


def benchmark_summary(results_dir: Path) -> None:
    benchmark = pd.read_csv(results_dir / "benchmark.tsv", sep="\t")
    benchmark = benchmark.loc[benchmark["returncode"] == 0]
    summary_rows: list[dict[str, object]] = []
    group_columns = ["benchmark", "software", "shape", "num_tips", "num_trees"]
    for keys, group in benchmark.groupby(group_columns, dropna=False):
        row = dict(zip(group_columns, keys))
        for column in ("elapsed_s", "peak_rss_mb"):
            median, minimum, maximum = median_range(group, column)
            row[f"{column}_median"] = median
            row[f"{column}_min"] = minimum
            row[f"{column}_max"] = maximum
        row["repetitions"] = len(group.index)
        summary_rows.append(row)
    pd.DataFrame(summary_rows).to_csv(results_dir / "benchmark_summary.tsv", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--figures-dir", type=Path, default=DEFAULT_FIGURES_DIR)
    args = parser.parse_args()
    benchmark_summary(args.results_dir)
    figure1_architecture(args.figures_dir)
    figure2_pepc(args.results_dir, args.figures_dir)
    figure3_evaluation(args.results_dir, args.figures_dir)


if __name__ == "__main__":
    main()
