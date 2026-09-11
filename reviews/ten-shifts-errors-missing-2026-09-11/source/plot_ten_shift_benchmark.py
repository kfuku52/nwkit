"""Plot every completed dataset and paired exact-branch F1."""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    rows = json.loads((args.directory / "accuracy/results.json").read_text())
    assert len(rows) == 20, "Wait for all prespecified outcomes."
    index = {
        (r["job"]["truth"], r["job"]["replicate"], r["job"]["mode"]): r for r in rows
    }
    colors = ["#176d9c", "#db7a24"]
    modes = ["shared", "trait-specific"]
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharey=True)
    for rowidx, criterion in enumerate(("BIC", "AIC")):
        for colidx, truth in enumerate(("shared", "different")):
            ax = axes[rowidx, colidx]
            for rep in range(5):
                offset = 0.045 * (rep - 2)
                pair = [index[(truth, rep, mode)] for mode in modes]
                complete = [r["status"] == "complete" for r in pair]
                if all(complete):
                    ax.plot(
                        [offset, 1 + offset],
                        [r["selected"][criterion]["f1"] for r in pair],
                        color="#bac2c9",
                        lw=1.2,
                        zorder=1,
                    )
                for x, r in enumerate(pair):
                    if complete[x]:
                        ax.scatter(
                            x + offset,
                            r["selected"][criterion]["f1"],
                            s=54,
                            edgecolors=colors[x],
                            facecolors=colors[x] if all(complete) else "white",
                            zorder=2,
                        )
            for x, mode in enumerate(modes):
                count = sum(
                    index[(truth, rep, mode)]["status"] == "complete"
                    for rep in range(5)
                )
                ax.text(
                    x,
                    -0.19,
                    f"{count}/5 completed",
                    ha="center",
                    transform=ax.get_xaxis_transform(),
                    fontsize=10,
                )
            true_alpha = "[1, 1]" if truth == "shared" else "[0.25, 4]"
            ax.set_title(
                f"{criterion} | True alpha = {true_alpha}", fontsize=12, pad=12
            )
            ax.set_xticks([0, 1], ["Shared alpha", "Trait-specific alpha"])
            ax.set_xlim(-0.25, 1.25)
            ax.set_ylim(-0.03, 1.03)
            ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
            ax.grid(axis="y", color="#e5e9ed", zorder=0)
            ax.spines[["top", "right"]].set_visible(False)
            if colidx == 0:
                ax.set_ylabel("Exact-branch F1")
    fig.suptitle(
        "100 tips, 10 true shifts: errors and missing values", fontsize=16, y=0.98
    )
    fig.text(
        0.5,
        0.925,
        "2 traits | 20% MCAR missingness | known SE 0.1 + extra error variance 0.04",
        ha="center",
        fontsize=10,
    )
    handles = [
        Line2D(
            [0],
            [0],
            color="#bac2c9",
            lw=1.2,
            label="Same dataset; both models completed",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            markerfacecolor="white",
            color="#333333",
            linestyle="none",
            label="Only this model completed",
        ),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.01),
        ncol=2,
        frameon=False,
        fontsize=9,
    )
    fig.subplots_adjust(top=0.86, bottom=0.15, hspace=0.6, wspace=0.2)
    fig.savefig(args.directory / "accuracy.png", dpi=180, facecolor="white")
    fig.savefig(args.directory / "accuracy.svg", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
