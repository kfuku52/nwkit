"""One-page continuous model panels, reusing comparison fits without refitting."""

import hashlib
import textwrap

import numpy as np
import pandas as pd

from nwkit.asr_figure import (
    _coordinates,
    _draw_simulation,
    _draw_trait,
    _draw_tree,
    _label_overhang,
    _regime_styles,
    _theta_values,
    continuous_figure_table,
    event_legend,
    figure_node_types,
    validate_figure_options,
)
from nwkit.asr_heatmap import draw_tip_heatmap, heatmap_extra_space
from nwkit.asr_paths import PATH_MODELS, simulate_fitted_paths, simulation_options
from nwkit.util import assign_branch_ids


def panels_requested(args):
    return getattr(args, "figure_layout", "table") == "panels"


def validate_panel_options(args, trait_type):
    settings = (
        "figure_width",
        "figure_height",
        "figure_simulation_mode",
        "figure_simulation_steps",
    )
    if not panels_requested(args):
        if (
            getattr(args, "figure_simulations", 0)
            or getattr(args, "figure_tip_heatmap", "no") == "yes"
            or getattr(args, "figure_trait_tip_labels", "no") == "yes"
            or any(getattr(args, key, None) is not None for key in settings)
        ):
            raise ValueError(
                "Figure dimensions, heatmaps and simulations require --figure-layout panels."
            )
        return
    if not getattr(args, "figure_out", None):
        raise ValueError("--figure-layout panels requires --figure-out.")
    if trait_type != "continuous":
        raise ValueError(
            "--figure-layout panels currently supports continuous traits only."
        )
    validate_figure_options(args)


def _candidate_seed(seed, model_id):
    """Stable across sorting, exclusions and Python processes."""
    if seed is None:
        return None
    digest = hashlib.sha256(model_id.encode("utf-8")).digest()
    return np.random.SeedSequence([seed, *np.frombuffer(digest, dtype="<u4").tolist()])


def _simulation(context, model_id, cached):
    posterior, fit, settings, assignment = cached
    count, steps, mode = simulation_options(context.args)
    if not count or settings.model not in PATH_MODELS:
        return None
    observed, errors = context.cache["continuous_data"]
    return simulate_fitted_paths(
        context.tree,
        observed,
        errors,
        posterior,
        fit=fit,
        model=settings.model,
        root_prior=settings.root_prior,
        regime_assignment=assignment,
        count=count,
        steps=steps,
        mode=mode,
        seed=_candidate_seed(getattr(context.args, "seed", None), model_id),
    )


def _row_title(row, criterion):
    value = row[criterion]
    score = (
        f"{criterion.upper()}={float(value):.3f}"
        if np.isfinite(value)
        else f"{criterion.upper()} unavailable"
    )
    rank = row["criterion_rank"]
    ranking = f" | rank {int(rank)} within set" if pd.notna(rank) else ""
    group = row["comparison_group"] or "Not assigned"
    return f"{row['model_id']} | {row['status']} | {score}{ranking}\nSet: {group}"


def _draw_model(context, row, cached, figure, slot, trait_axes):
    posterior, fit, settings, assignment = cached
    observed, errors = context.cache["continuous_data"]
    table = continuous_figure_table(
        context.tree, observed, errors, posterior, context.trait_columns, settings
    )
    simulation = _simulation(context, row["model_id"], cached)
    include_simulation = bool(getattr(context.args, "figure_simulations", 0))
    stride = 2 if include_simulation else 1
    nodes, tips, positions, depths = _coordinates(context.tree)
    by_node, styles = _regime_styles(context.tree, assignment)
    node_types = context.cache["figure_node_types"]
    heatmap = getattr(context.args, "figure_tip_heatmap", "no") == "yes"
    ids = assign_branch_ids(context.tree)
    grid = slot.subgridspec(1, 1 + stride * len(context.trait_columns), wspace=0.24)
    tree_ax = figure.add_subplot(grid[0])
    _draw_tree(
        tree_ax,
        nodes,
        tips,
        positions,
        depths,
        by_node,
        styles,
        node_types,
        heatmap_traits=len(context.trait_columns) if heatmap else 0,
    )
    if heatmap:
        draw_tip_heatmap(tree_ax, table, tips, positions, _label_overhang(tips) + 0.65)
    from nwkit.asr_tip_labels import draw_trait_tip_labels

    tip_labels = getattr(context.args, "figure_trait_tip_labels", "no") == "yes"
    tip_colors = {tip.name: styles[by_node[tip]][0] for tip in tips}
    for index, trait in enumerate(context.trait_columns):
        rows = table[table.trait == trait].set_index("branch_id")
        theta = _theta_values(fit, index, styles)
        ax = figure.add_subplot(grid[1 + index * stride], sharey=tree_ax)
        _draw_trait(ax, rows, nodes, depths, ids, by_node, styles, theta, node_types)
        if tip_labels:
            draw_trait_tip_labels(
                ax,
                {tip.name: float(rows.loc[ids[tip], "mean"]) for tip in tips},
                tip_colors,
            )
        ax.set_title(f"{trait} | ASR", loc="left", fontsize=10, fontweight="bold")
        trait_axes[trait].append(ax)
        if include_simulation:
            sim_ax = figure.add_subplot(grid[2 + index * stride], sharey=tree_ax)
            if simulation is None:
                sim_ax.text(
                    0.5,
                    0.5,
                    "Branch simulation unavailable\nfor this transformed model",
                    ha="center",
                    va="center",
                    transform=sim_ax.transAxes,
                    fontsize=9,
                )
                sim_ax.tick_params(axis="y", labelleft=False)
                sim_ax.set_xticks([])
            else:
                _draw_simulation(
                    sim_ax,
                    simulation,
                    index,
                    trait,
                    nodes,
                    depths,
                    by_node,
                    styles,
                    theta,
                    node_types,
                )
                sim_ax.set_title(sim_ax.get_title(loc="left"), loc="left", fontsize=10)
                trait_axes[trait].append(sim_ax)
                if tip_labels:
                    draw_trait_tip_labels(
                        sim_ax,
                        {
                            tip.name: float(
                                simulation.root_values[0, index]
                                if tip.is_root
                                else simulation.branches[tip].values[0, -1, index]
                            )
                            for tip in tips
                        },
                        tip_colors,
                        first_history=simulation.count > 1,
                    )
    extent = max(depths.values()) or 1.0
    tree_ax.set_ylim(extent * 1.04, -extent * 0.05)
    tree_ax.set_title("Phylogeny", loc="left", fontsize=10, fontweight="bold")
    tree_ax.set_ylabel("Distance from root")
    return simulation, styles


def _row_notes(ax, simulation, styles, node_types, *, legend_left):
    from matplotlib.lines import Line2D

    ax.set_axis_off()
    note = simulation.root_description if simulation is not None else ""
    if simulation is not None and simulation.mode == "unconditional":
        note = "New histories, not conditioned on tips. " + note
    ax.text(
        0, 0.06, note, fontsize=8, color="#555555", va="bottom", transform=ax.transAxes
    )
    if any(styles) or node_types:
        handles = [
            Line2D([], [], color=color, linestyle=style, label=regime)
            for regime, (color, style) in styles.items()
            if regime
        ] + event_legend(node_types)
        ax.legend(
            handles=handles,
            loc="upper left",
            bbox_to_anchor=(legend_left, 1),
            ncol=min(6, len(handles)),
            frameon=False,
            fontsize=8,
            borderaxespad=0,
        )


def build_comparison_panels(context, table):
    from matplotlib.figure import Figure

    records = table.sort_values(
        ["comparison_group", "criterion_rank"], kind="stable", na_position="last"
    ).to_dict("records")
    fits = context.cache.get("figure_fits", {})
    context.cache["figure_node_types"] = figure_node_types(context.tree, context.args)
    overhang = _label_overhang(list(context.tree.leaves()))
    if getattr(context.args, "figure_tip_heatmap", "no") == "yes":
        overhang += heatmap_extra_space(len(context.trait_columns))
    if getattr(context.args, "figure_trait_tip_labels", "no") == "yes":
        overhang += 0.45
    heights = [
        5.3 + overhang
        if row["model_id"] in fits
        else 1.1 + 0.14 * max(0, len(str(row["message"])) // 120)
        for row in records
    ]
    count = len(context.trait_columns)
    panels = count * (2 if getattr(context.args, "figure_simulations", 0) else 1)
    width = (
        getattr(context.args, "figure_width", None)
        or max(3.6, len(list(context.tree.leaves())) * 0.36) + 3.8 * panels
    )
    height = getattr(context.args, "figure_height", None) or sum(heights) + 1.8
    figure = Figure(figsize=(width, height))
    grid = figure.add_gridspec(
        len(records),
        1,
        height_ratios=heights,
        left=0.08,
        right=0.97,
        bottom=0.9 / height,
        top=1 - 0.9 / height,
        hspace=0.1,
    )
    trait_axes: dict[str, list] = {trait: [] for trait in context.trait_columns}
    for index, row in enumerate(records):
        cached = fits.get(row["model_id"])
        row_grid = (
            grid[index].subgridspec(
                3,
                1,
                height_ratios=[0.65, 3.65, 1.0 + overhang],
                hspace=0.3,
            )
            if cached is not None
            else grid[index].subgridspec(1, 1)
        )
        title_ax = figure.add_subplot(row_grid[0])
        title_ax.set_axis_off()
        title_ax.text(
            0,
            1,
            _row_title(row, getattr(context.args, "criterion", "aic")),
            ha="left",
            va="top",
            fontsize=10,
            fontweight="bold",
            transform=title_ax.transAxes,
        )
        if cached is None:
            title_ax.text(
                0,
                0,
                textwrap.fill(str(row["message"]), width=max(60, int(width * 13))),
                fontsize=8,
                va="bottom",
                transform=title_ax.transAxes,
            )
            continue
        simulation, styles = _draw_model(
            context, row, cached, figure, row_grid[1], trait_axes
        )
        notes_grid = row_grid[2].subgridspec(2, 1, height_ratios=[0.1 + overhang, 0.9])
        _row_notes(
            figure.add_subplot(notes_grid[1]),
            simulation,
            styles,
            context.cache["figure_node_types"],
            legend_left=0.03 + 1 / (1 + panels),
        )
    for axes in trait_axes.values():
        limits = (
            (min(ax.get_xlim()[0] for ax in axes), max(ax.get_xlim()[1] for ax in axes))
            if axes
            else None
        )
        for ax in axes:
            ax.set_xlim(limits)
    for ax in figure.axes:
        if ax.axison:
            ax.spines[:].set_visible(True)
            ax.spines[:].set_color("#BBBBBB")
            ax.tick_params(labelsize=8)
            ax.grid(axis="y", color="#EEEEEE", linewidth=0.6)
    figure.text(
        0.08,
        1 - 0.25 / height,
        "ASR model comparison | ancestral traits and histories",
        va="top",
        fontsize=16,
        fontweight="bold",
    )
    figure.text(
        0.08,
        0.15 / height,
        "Open circles: node means; bars: 95% marginal intervals; squares: observed tips; diamonds: imputed tips; dashed lines: OU optima.\nASR segments connect node means, not branch paths. Intervals condition on fitted parameters and the input tree.\nDepth uses input branch-length units. Trait scales are shared across models. IC ranks apply only within the stated comparison set.",
        fontsize=8,
        va="bottom",
        color="#555555",
    )
    return figure


def draw_comparison_panels(context, table, path):
    from matplotlib import rc_context

    from nwkit.asr_compare_figure import _font_family_for_text

    text = " ".join(
        [
            *context.trait_columns,
            *context.tree.leaf_names(),
            *table.model_id.astype(str),
        ]
    )
    family = _font_family_for_text(text)
    with rc_context({"font.family": family} if family else {}):
        figure = build_comparison_panels(context, table)
        try:
            figure.savefig(path, format="pdf", facecolor="white")
        finally:
            figure.clear()
