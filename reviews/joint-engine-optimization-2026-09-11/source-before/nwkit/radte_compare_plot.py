"""Tree-based comparison panels with common ages and matched node positions."""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from nwkit.result_plot import (
    _age_extent,
    _draw_tree_panel,
    _layout,
    _species_tip_labels,
)

MODES = ("fixed", "bounded", "ensemble")
TITLES = (
    "Fixed species ages",
    "Species ages fitted within bounds",
    "Sampled species chronograms",
)


def _grid(reference, *, width=18, top=0.79):
    gene_rows = _layout(reference.gene)[-1]
    species_rows = _layout(reference.species, _species_tip_labels(reference.species))[
        -1
    ]
    heights = [max(2, 0.30 * gene_rows), max(1.7, 0.30 * species_rows)]
    fig, axes = plt.subplots(
        2,
        3,
        figsize=(width, 3.5 + sum(heights)),
        gridspec_kw={"height_ratios": heights},
    )
    fig.subplots_adjust(
        left=0.045, right=0.98, top=top, bottom=0.12, wspace=0.18, hspace=0.5
    )
    return fig, axes


def _interval_caption(data):
    rows = [r for r in data.events.values() if r["event_type"] != "leaf"]
    available = sum(r.get("interval_lower") is not None for r in rows)
    method = data.manifest.get("options", {}).get("uncertainty", "none")
    method = "input percentiles" if method == "input-ensemble" else method
    return f"{method}; gene intervals available: {available}/{len(rows)}"


def comparison_figure(runs, level):
    fig, axes = _grid(runs["fixed"], top=0.75)
    extent = max(_age_extent(data) for data in runs.values())
    fig.suptitle(
        "Compare ages on the same trees",
        x=0.045,
        ha="left",
        fontsize=19,
        weight="bold",
        y=0.98,
    )
    fig.text(
        0.045,
        0.928,
        "Each column shows the same gene tree (top) and species tree (bottom). All six panels share the same age axis.",
        fontsize=10,
    )
    fig.text(
        0.045,
        0.897,
        f"Dots/squares mark node ages; horizontal bars show {100 * level:g}% RADTE intervals. Older ages are to the left.",
        fontsize=10,
    )
    fig.legend(
        handles=[
            Line2D(
                [], [], marker="o", linestyle="", color="#0072B2", label="Speciation"
            ),
            Line2D(
                [], [], marker="s", linestyle="", color="#D55E00", label="Duplication"
            ),
            Line2D([], [], linewidth=4, color="#7153aa", label="Gene age interval"),
            Line2D(
                [], [], linewidth=4, color="#94a3b8", label="External species interval"
            ),
            Line2D([], [], color="#475569", linestyle=":", label="Hard range"),
        ],
        loc="upper left",
        bbox_to_anchor=(0.04, 0.875),
        ncol=5,
        frameon=False,
        fontsize=9,
    )
    for col, (mode, title) in enumerate(zip(MODES, TITLES, strict=True)):
        for row, gene in enumerate((True, False)):
            ax = axes[row, col]
            _draw_tree_panel(ax, runs[mode], gene=gene, font_size=8, extent=extent)
            ax.set_title(
                title + "\n" + _interval_caption(runs[mode])
                if gene
                else "Species tree",
                fontsize=10,
                loc="left",
                pad=12,
            )
            ax.set_xlabel("Age (input time units)", fontsize=9)
            ax.set_gid(f"comparison:{mode}:{'gene' if gene else 'species'}")
    fig.text(
        0.045,
        0.035,
        "Blue species bars: fitted intervals (middle) or input-sample percentiles (right). Gray bars: external input evidence.\n"
        "Ensemble branches retain the reference point fit; their bars describe variation across chronograms. Conditional intervals and ensemble percentiles are not a combined CI.",
        fontsize=9,
    )
    return fig


def components_figure(components, reference):
    fig, axes = _grid(reference)
    extent = _age_extent(reference)
    fig.suptitle(
        "Uncertainty summaries aligned with tree nodes",
        x=0.045,
        ha="left",
        fontsize=19,
        weight="bold",
        y=0.98,
    )
    fig.text(
        0.045,
        0.927,
        "Read horizontally from a node in the reference tree to its two uncertainty summaries. Dashed rows identify the matching node.",
        fontsize=10,
    )
    fig.text(
        0.045,
        0.89,
        "The two summary columns measure different quantities; their horizontal scales are not age axes and must not be added.",
        fontsize=10,
    )
    indexed = components.set_index("shared_age_id")
    for row, gene in enumerate((True, False)):
        tree_ax = axes[row, 0]
        _draw_tree_panel(tree_ax, reference, gene=gene, font_size=8, extent=extent)
        tree_ax.set_title(
            "Reference gene tree" if gene else "Reference species tree",
            loc="left",
            fontsize=11,
        )
        tree_ax.set_xlabel("Age (input time units)", fontsize=9)
        tree = reference.gene if gene else reference.species
        index = reference.gene_index if gene else reference.species_index
        y = _layout(tree, None if gene else _species_tip_labels(tree))[1]
        for col, (field, title) in enumerate(
            (
                ("input_refit_sd", "Across chronograms: SD of point ages"),
                (
                    "mean_conditional_interval_width",
                    "Within chronograms: mean interval width",
                ),
            ),
            1,
        ):
            ax = axes[row, col]
            for node in tree.traverse():
                if node.is_leaf:
                    continue
                sid = index.clade_id_for_node(node)
                key = reference.events[sid]["shared_age_id"] if gene else "S:" + sid
                value = (
                    float(indexed.loc[key, field]) if key in indexed.index else np.nan
                )
                ax.axhline(y[node], color="#d6dde5", linestyle="--", linewidth=0.7)
                if np.isfinite(value):
                    artist = ax.scatter(value, y[node], color="#0072B2", s=25)
                    artist.set_gid("component:" + field + ":" + sid)
                else:
                    ax.text(0, y[node], "Unavailable", va="center", fontsize=8)
            ax.set_ylim(tree_ax.get_ylim())
            ax.set_yticks([])
            ax.set_xlim(left=-0.05 * max(1, ax.get_xlim()[1]))
            ax.set_title(title, fontsize=10, loc="left", pad=12)
            ax.set_xlabel("Input time units", fontsize=9)
            for side in ("top", "right", "left"):
                ax.spines[side].set_visible(False)
    for col in (1, 2):
        upper = max(ax.get_xlim()[1] for ax in axes[:, col])
        for ax in axes[:, col]:
            ax.set_xlim(-0.03 * upper, upper)
    fig.text(
        0.045,
        0.035,
        "Species ages fixed within a chronogram have zero conditional width but may vary between chronograms.\n"
        "Repeated gene speciation nodes share the same species-age summary. Coverage, failures and bootstrap variances are retained in the TSV.",
        fontsize=9,
    )
    return fig
