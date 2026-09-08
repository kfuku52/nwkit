"""Continuous ancestral traits beside a tree on a shared branch-depth axis.

Connectors join node summaries only: they are not conditional branch paths.
Intervals are node marginals conditional on the fitted model parameters.
"""

import math
import sys
from pathlib import Path

import numpy as np

from nwkit.output_transaction import output_transaction
from nwkit.util import assign_branch_ids, validate_outputs_do_not_replace_inputs

_COLORS = ("#333333", "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#B58800")
_STYLES = ("-", "--", "-.", ":")
_NODE_MARKER_SIZE = 24
_NODE_MARKER_LINEWIDTH = 0.8
_TIP_LABEL_ROTATION = 90


def validate_figure_trait(args, trait_type):
    if trait_type != "continuous" and getattr(args, "figure_out", None) not in (
        None,
        "",
    ):
        raise ValueError("--figure-out currently supports continuous ASR only.")


def validate_figure_options(args):
    """Reject invalid destinations/settings before fitting or writing outputs."""
    from nwkit.asr_paths import simulation_options

    simulation_options(args)
    path = getattr(args, "figure_out", None)
    if getattr(args, "figure_trait_tip_labels", "no") == "yes" and path in (None, ""):
        raise ValueError("--figure-trait-tip-labels yes requires --figure-out.")
    if getattr(args, "figure_tip_heatmap", "no") == "yes" and path in (None, ""):
        raise ValueError("--figure-tip-heatmap yes requires --figure-out.")
    for name in ("figure_width", "figure_height"):
        value = getattr(args, name, None)
        if value is not None:
            if not math.isfinite(value) or value <= 0:
                raise ValueError(
                    f"--{name.replace('_', '-')} must be positive and finite."
                )
            if path in (None, ""):
                raise ValueError(f"--{name.replace('_', '-')} requires --figure-out.")
    if path in (None, ""):
        return
    if Path(path).suffix.lower() not in {".pdf", ".svg", ".png"}:
        raise ValueError("--figure-out must use a .pdf, .svg, or .png extension.")
    inputs = [
        ("--" + name.replace("_", "-"), getattr(args, name, None))
        for name in (
            "infile",
            "trait",
            "regime_map",
            "regime_parameters",
            "rate_matrix",
            "rate_design",
            "transition_graph",
            "species_map_tsv",
        )
    ]
    validate_outputs_do_not_replace_inputs(
        inputs, [("--figure-out", path)], label="ASR figure output"
    )


def _coordinates(tree):
    nodes = list(tree.traverse("preorder"))
    tips = list(tree.leaves())
    positions = {node: float(index) for index, node in enumerate(tips)}
    for node in tree.traverse("postorder"):
        if not node.is_leaf:
            positions[node] = (
                positions[node.children[0]] + positions[node.children[-1]]
            ) / 2
    depths = {tree: 0.0}
    for node in nodes[1:]:
        depths[node] = depths[node.up] + float(node.dist)
        if not math.isfinite(depths[node]):
            raise ValueError(
                "ASR figure branch depths overflow; rescale branch lengths."
            )
    return nodes, tips, positions, depths


def _regime_styles(tree, assignment):
    regimes = assignment.regimes if assignment is not None else ("",)
    by_node = (
        assignment.by_node
        if assignment is not None
        else dict.fromkeys(tree.traverse(), "")
    )
    styles = {
        regime: (
            _COLORS[index % len(_COLORS)],
            _STYLES[(index // len(_COLORS)) % len(_STYLES)],
        )
        for index, regime in enumerate(regimes)
    }
    return by_node, styles


def figure_node_types(tree, args):
    """Use the same label parsing and overlap decisions as `draw`."""
    from nwkit.draw_helpers import _get_species_overlap_node_types

    mode = getattr(args, "species_overlap_node_plot", "auto")
    if mode not in {"auto", "yes", "no"}:
        raise ValueError("--species-overlap-node-plot must be yes, no, or auto.")
    if mode == "no":
        return {}
    types, parsed = _get_species_overlap_node_types(
        tree, args, require_all_tip_labels=mode == "auto"
    )
    if mode == "auto" and not parsed:
        sys.stderr.write(
            "Skipping speciation/duplication figure colors because some leaf labels "
            "did not match the configured species parser.\n"
        )
    return types


def _event_color(node, node_types, default):
    from nwkit.draw_helpers import DUPLICATION_COLOR, SPECIATION_COLOR

    event = (node_types or {}).get(node, "")
    return {"speciation": SPECIATION_COLOR, "duplication": DUPLICATION_COLOR}.get(
        event, default
    )


def event_legend(node_types):
    from matplotlib.lines import Line2D

    return [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="none",
            markerfacecolor=_event_color(node, node_types, "#333333"),
            markeredgecolor="white",
            label=event.capitalize(),
        )
        for event, node in {
            event: node for node, event in (node_types or {}).items()
        }.items()
    ]


def _draw_tree(
    ax,
    nodes,
    tips,
    positions,
    depths,
    by_node,
    styles,
    node_types=None,
    heatmap_traits=0,
):
    from nwkit.asr_heatmap import heatmap_height

    for node in nodes:
        color, linestyle = styles[by_node[node]]
        if not node.is_leaf:
            ax.plot(
                [positions[node.children[0]], positions[node.children[-1]]],
                [depths[node], depths[node]],
                color="#999999",
                linewidth=0.8,
            )
        if not node.is_root:
            ax.plot(
                [positions[node], positions[node]],
                [depths[node.up], depths[node]],
                color=color,
                linestyle=linestyle,
                linewidth=1.5,
            )
        if (
            node.is_root
            or by_node[node] != by_node[node.up]
            or node in (node_types or {})
        ):
            ax.scatter(
                positions[node],
                depths[node],
                s=_NODE_MARKER_SIZE,
                color=_event_color(node, node_types, color),
                linewidths=_NODE_MARKER_LINEWIDTH,
                zorder=4,
            )
    for node in tips:
        ax.annotate(
            str(node.name),
            (positions[node], 0 if heatmap_traits else depths[node]),
            xycoords=ax.get_xaxis_transform() if heatmap_traits else "data",
            xytext=(
                0,
                -8 - (72 * heatmap_height(heatmap_traits) if heatmap_traits else 0),
            ),
            textcoords="offset points",
            rotation=_TIP_LABEL_ROTATION,
            ha="center",
            va="top",
            fontsize=9,
            color=styles[by_node[node]][0],
            annotation_clip=False,
        )
    ax.set_xlim(-0.8, max(0.8, len(tips) - 0.2))
    ax.set_xticks([])
    ax.set_ylabel("Distance from root (input branch-length units)")
    ax.set_title("Phylogeny", loc="left", fontweight="bold", pad=14)


def _theta_values(fit, trait_index, styles):
    values = getattr(fit, "theta_by_regime", None)
    if values is not None:
        return values
    theta = getattr(fit, "theta", None)
    if theta is None:
        return {}
    vector = np.asarray(theta)
    value = float(vector) if vector.ndim == 0 else float(vector[trait_index])
    return dict.fromkeys(styles, value)


def _label_overhang(tips):
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import TextToPath

    measure = TextToPath()
    font = FontProperties(size=9)
    angle = math.radians(_TIP_LABEL_ROTATION)
    extents = []
    for node in tips:
        width, height, _ = measure.get_text_width_height_descent(
            str(node.name), font, False
        )
        extents.append(width * math.sin(angle) + height * math.cos(angle))
    return max(0.0, (max(extents, default=0.0) + 8) / 72 - 0.65)


def _draw_theta(ax, nodes, depths, by_node, styles, theta):
    for regime, value in theta.items():
        members = [node for node in nodes if by_node[node] == regime]
        low = min(depths[node] if node.is_root else depths[node.up] for node in members)
        high = max(depths[node] for node in members)
        ax.plot(
            [value, value],
            [low, high],
            color=styles[regime][0],
            linestyle=(0, (4, 4)),
            linewidth=1,
            alpha=0.65,
        )


def _draw_trait(ax, rows, nodes, depths, ids, by_node, styles, theta, node_types=None):
    _draw_theta(ax, nodes, depths, by_node, styles, theta)
    for node in nodes:
        row = rows.loc[ids[node]]
        color, linestyle = styles[by_node[node]]
        mean = float(row["mean"])
        depth = depths[node]
        if not node.is_root:
            ax.plot(
                [float(rows.loc[ids[node.up], "mean"]), mean],
                [depths[node.up], depth],
                color=color,
                linestyle=linestyle,
                linewidth=1,
                alpha=0.65,
                zorder=2,
            )
        ax.plot(
            [row["ci_lower"], row["ci_upper"]],
            [depth, depth],
            color=_event_color(node, node_types, color),
            linewidth=2.4,
            zorder=3,
        )
        ax.scatter(
            mean,
            depth,
            s=_NODE_MARKER_SIZE,
            marker="D" if row["is_imputed"] else "o",
            facecolor=_event_color(node, node_types, "white"),
            edgecolor=_event_color(node, node_types, color),
            linewidth=_NODE_MARKER_LINEWIDTH,
            zorder=4,
        )
        if row["observed_value"] != "":
            ax.scatter(
                float(row["observed_value"]),
                depth,
                s=18,
                marker="s",
                color=color,
                zorder=5,
            )
    ax.set_xlabel(str(rows.iloc[0]["trait"]))
    ax.set_title(str(rows.iloc[0]["trait"]), loc="left", fontweight="bold", pad=14)
    ax.tick_params(axis="y", labelleft=False)
    ax.margins(x=0.12)


def _draw_simulation(
    ax,
    simulation,
    trait_index,
    trait,
    nodes,
    depths,
    by_node,
    styles,
    theta,
    node_types=None,
):
    from matplotlib.collections import LineCollection

    _draw_theta(ax, nodes, depths, by_node, styles, theta)
    segments: dict[str, list[np.ndarray]] = {regime: [] for regime in styles}
    for node, path in simulation.branches.items():
        vertices = np.empty((simulation.count, len(path.elapsed), 2))
        vertices[:, :, 0] = path.values[:, :, trait_index]
        vertices[:, :, 1] = depths[node.up] + path.elapsed
        segments[by_node[node]].extend(vertices)
    for regime, paths in segments.items():
        color, linestyle = styles[regime]
        ax.add_collection(
            LineCollection(
                paths,
                colors=color,
                linestyles=linestyle,
                linewidths=0.8,
                alpha=min(0.8, 1 / math.sqrt(simulation.count)),
                zorder=2,
            )
        )
    ax.autoscale_view()
    for node in nodes:
        if not node.is_root and len(node.children) < 2:
            continue
        values = (
            simulation.root_values[:, trait_index]
            if node.is_root
            else simulation.branches[node].values[:, -1, trait_index]
        )
        ax.scatter(
            values,
            np.full(simulation.count, depths[node]),
            s=_NODE_MARKER_SIZE,
            color=_event_color(node, node_types, "#333333"),
            linewidths=_NODE_MARKER_LINEWIDTH,
            zorder=4,
        )
    label = "simulation" if simulation.mode == "unconditional" else "posterior paths"
    ax.set_title(f"{trait} | {label}", loc="left", fontweight="bold", pad=14)
    ax.set_xlabel(trait)
    ax.tick_params(axis="y", labelleft=False)
    ax.margins(x=0.12)


def build_continuous_asr_figure(
    tree,
    table,
    *,
    model,
    fit,
    regime_assignment=None,
    width=None,
    height=None,
    simulation=None,
    node_types=None,
    tip_heatmap=False,
    trait_tip_labels=False,
):
    """Build a headless Matplotlib Figure using a complete continuous summary table."""
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D

    from nwkit.asr_heatmap import draw_tip_heatmap, heatmap_extra_space

    nodes, tips, positions, depths = _coordinates(tree)
    ids = assign_branch_ids(tree)
    traits = list(dict.fromkeys(table["trait"]))
    by_node, styles = _regime_styles(tree, regime_assignment)
    tree_width = max(3.4, len(tips) * 0.36)
    overhang = _label_overhang(tips) + (0.45 if trait_tip_labels else 0)
    extra_footer = 0.2 if simulation is not None else 0.0
    heatmap_space = heatmap_extra_space(len(traits)) if tip_heatmap else 0.0
    panels_per_trait = 2 if simulation is not None else 1
    num_panels = len(traits) * panels_per_trait
    figure_height = height or 7.5 + overhang + extra_footer + heatmap_space
    figure = Figure(figsize=(width or tree_width + 3.8 * num_panels, figure_height))
    axes = figure.subplots(
        1,
        1 + num_panels,
        sharey=True,
    )
    _draw_tree(
        axes[0],
        nodes,
        tips,
        positions,
        depths,
        by_node,
        styles,
        node_types,
        heatmap_traits=len(traits) if tip_heatmap else 0,
    )
    if tip_heatmap:
        draw_tip_heatmap(axes[0], table, tips, positions, overhang + 0.65)
    from nwkit.asr_tip_labels import draw_trait_tip_labels

    tip_colors = {tip.name: styles[by_node[tip]][0] for tip in tips}
    has_theta = False
    for index, trait in enumerate(traits):
        rows = table[table.trait == trait].set_index("branch_id")
        if set(rows.index) != set(ids.values()) or not rows.index.is_unique:
            raise ValueError(
                "ASR figures require exactly one summary per node and trait."
            )
        theta = _theta_values(fit, index, styles)
        has_theta = has_theta or bool(theta)
        ax = axes[1 + index * panels_per_trait]
        _draw_trait(ax, rows, nodes, depths, ids, by_node, styles, theta, node_types)
        if trait_tip_labels:
            draw_trait_tip_labels(
                ax,
                {tip.name: float(rows.loc[ids[tip], "mean"]) for tip in tips},
                tip_colors,
            )
        if simulation is not None:
            ax.set_title(f"{trait} | ASR", loc="left", fontweight="bold", pad=14)
            simulated_ax = axes[2 + index * panels_per_trait]
            _draw_simulation(
                simulated_ax,
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
            limits = (
                min(ax.get_xlim()[0], simulated_ax.get_xlim()[0]),
                max(ax.get_xlim()[1], simulated_ax.get_xlim()[1]),
            )
            ax.set_xlim(limits)
            simulated_ax.set_xlim(limits)
            if trait_tip_labels:
                draw_trait_tip_labels(
                    simulated_ax,
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
    axes[0].set_ylim(extent * 1.04, -extent * 0.05)
    for ax in axes:
        ax.spines[:].set_visible(True)
        ax.spines[:].set_color("#BBBBBB")
        ax.tick_params(labelsize=9, color="#BBBBBB")
        ax.grid(axis="y", color="#EEEEEE", linewidth=0.7)
        ax.set_axisbelow(True)
    figure.suptitle(
        f"Ancestral trait reconstruction | {model}",
        x=0.08,
        ha="left",
        fontsize=16,
        fontweight="bold",
        y=0.97,
    )
    level = float(table.iloc[0]["ci_level"])
    legend = [
        Line2D(
            [],
            [],
            marker="o",
            color="#333333",
            markerfacecolor="white",
            linestyle="none",
            label="Node mean",
        ),
        Line2D(
            [], [], marker="s", color="#333333", linestyle="none", label="Observed tip"
        ),
        Line2D(
            [],
            [],
            color="#999999",
            linewidth=3,
            label=f"{100 * level:g}% node interval",
        ),
    ]
    if bool(table["is_imputed"].any()):
        legend.append(
            Line2D(
                [],
                [],
                marker="D",
                color="#333333",
                markerfacecolor="white",
                linestyle="none",
                label="Imputed tip",
            )
        )
    if has_theta:
        legend.append(
            Line2D([], [], color="#555555", linestyle="--", label="OU optimum (theta)")
        )
    if simulation is not None:
        legend.append(
            Line2D(
                [],
                [],
                color="#555555",
                linewidth=0.8,
                label=f"{simulation.count} sampled {'history' if simulation.count == 1 else 'histories'}",
            )
        )
    figure.legend(
        handles=legend + event_legend(node_types),
        loc="lower left",
        bbox_to_anchor=(0.07, (0.7875 + extra_footer) / figure_height),
        frameon=False,
        ncol=min(3, len(legend)),
        fontsize=9,
    )
    if regime_assignment is not None:
        figure.legend(
            handles=[
                Line2D([], [], color=color, linestyle=style, label=regime)
                for regime, (color, style) in styles.items()
            ],
            loc="lower left",
            bbox_to_anchor=(0.07, (0.3375 + extra_footer) / figure_height),
            frameon=False,
            ncol=min(6, len(styles)),
            fontsize=9,
            title="Branch regimes",
        )
    notes = (
        "ASR lines connect node means; they are not reconstructed branch paths.\n"
        "Intervals condition on model parameters and the input tree."
    )
    if simulation is not None:
        conditioning = (
            "New histories, not conditioned on tips. "
            if simulation.mode == "unconditional"
            else ""
        )
        notes += f"\n{conditioning}{simulation.root_description}"
    figure.text(
        0.08,
        0.1125 / figure_height,
        notes,
        fontsize=8,
        color="#555555",
        va="bottom",
    )
    figure.subplots_adjust(
        left=0.09,
        right=0.97,
        top=0.86,
        bottom=min(
            0.8, (2.175 + overhang + extra_footer + heatmap_space) / figure_height
        ),
        wspace=0.2,
    )
    return figure


def write_continuous_asr_figure(
    tree, observed, errors, posterior, traits, args, settings, fit, assignment
):
    """Render all nodes, independently of the TSV --target selection."""
    if getattr(args, "figure_out", None) in (None, ""):
        return
    from nwkit.asr_paths import simulate_fitted_paths, simulation_options

    count, steps, mode = simulation_options(args)
    simulation = None
    if count:
        simulation = simulate_fitted_paths(
            tree,
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
            seed=getattr(args, "seed", None),
        )
    table = continuous_figure_table(tree, observed, errors, posterior, traits, settings)
    figure = build_continuous_asr_figure(
        tree,
        table,
        model=settings.model,
        fit=fit,
        regime_assignment=assignment,
        width=getattr(args, "figure_width", None),
        height=getattr(args, "figure_height", None),
        simulation=simulation,
        node_types=figure_node_types(tree, args),
        tip_heatmap=getattr(args, "figure_tip_heatmap", "no") == "yes",
        trait_tip_labels=getattr(args, "figure_trait_tip_labels", "no") == "yes",
    )
    try:
        with output_transaction([args.figure_out]) as staged:
            figure.savefig(
                staged[args.figure_out],
                format=Path(args.figure_out).suffix[1:].lower(),
                dpi=180,
                bbox_inches="tight",
                facecolor="white",
            )
    finally:
        figure.clear()


def continuous_figure_table(tree, observed, errors, posterior, traits, settings):
    """All-node summaries shared by individual and comparison figures."""
    if settings.model in {"MV-BM", "MV-OU", "MV-OU-DIAG", "MV-OU-FULL"}:
        from nwkit.multivariate_asr import multivariate_output_table

        return multivariate_output_table(
            tree,
            list(tree.traverse()),
            observed,
            posterior,
            traits,
            settings.ci_level,
            errors=errors,
        )
    else:
        from nwkit.continuous_asr_io import continuous_output_table

        return continuous_output_table(
            tree,
            list(tree.traverse()),
            observed,
            errors,
            posterior,
            trait=traits[0],
            ci_level=settings.ci_level,
        )
