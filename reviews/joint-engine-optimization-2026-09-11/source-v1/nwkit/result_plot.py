"""Source-backed reconciliation and RADTE reports using NWKIT tree geometry."""

import math
import textwrap
from collections import Counter
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

from nwkit.draw_layouts import get_rectangular_coordinates
from nwkit.output_transaction import output_transaction
from nwkit.result_plot_data import (
    load_dating_plot_data,
    radte_plot_protected_paths,
    read_result_table,
    reconciliation_plot_data,
)
from nwkit.util import read_tree, validate_outputs_do_not_replace_inputs

EVENT_STYLE = {
    "speciation": ("#0072B2", "o", "Speciation"),
    "duplication": ("#D55E00", "s", "Duplication"),
    "transfer": ("#CC79A7", ">", "Transfer"),
    "unresolved": ("#777777", "X", "Unresolved"),
    "leaf": ("#334155", "o", "Tip"),
}
INK = "#213247"
MUTED = "#596b80"
INTERVAL = "#7153aa"


def figure_format(path, requested="auto"):
    fmt = Path(path).suffix.lower().lstrip(".") if requested == "auto" else requested
    if fmt not in {"pdf", "svg", "png"}:
        raise ValueError(
            "Result figures require a .pdf, .svg, or .png output (or --image-format)."
        )
    if str(path) == "-":
        raise ValueError("Result figures require a filesystem output path.")
    return fmt


def _wrap(value, width=30):
    return textwrap.fill(
        str(value), width=width, break_long_words=True, break_on_hyphens=False
    )


def _wrap_identifier(value, width=30):
    remaining = str(value)
    lines = []
    while len(remaining) > width:
        split = remaining.rfind("_", 0, width) + 1
        if split < width // 3:
            split = width
        lines.append(remaining[:split])
        remaining = remaining[split:]
    return "\n".join([*lines, remaining])


def _layout(tree, tip_names=None):
    labels = {
        n: _wrap_identifier((tip_names or {}).get(n, n.name)) for n in tree.leaves()
    }
    weights = {n: max(1, label.count("\n") + 1) for n, label in labels.items()}
    x, y, leaves = get_rectangular_coordinates(
        tree, use_topology_depth=True, leaf_weight_by_leaf=weights
    )
    return x, y, leaves, labels, sum(weights.values())


def _internal_labels(tree, prefix):
    nodes = [n for n in tree.traverse() if not n.is_leaf]
    names = Counter(str(n.name or "") for n in nodes)
    labels = {
        n: str(n.name)
        for n in nodes
        if n.name and len(str(n.name)) <= 16 and names[str(n.name)] == 1
    }
    used = set(labels.values())
    counter = 1
    for node in nodes:
        if node in labels:
            continue
        while f"{prefix}{counter}" in used:
            counter += 1
        labels[node] = f"{prefix}{counter}"
        used.add(labels[node])
    return labels


def _species_labels(tree):
    labels = _internal_labels(tree, "s")
    leaves = sorted(tree.leaves(), key=lambda n: str(n.name))
    labels.update({n: str(n.name) for n in leaves if len(str(n.name)) <= 12})
    used = set(labels.values())
    counter = 1
    for node in leaves:
        if node in labels:
            continue
        while f"t{counter}" in used:
            counter += 1
        labels[node] = f"t{counter}"
        used.add(labels[node])
    return labels


def _species_tip_labels(tree):
    labels = _species_labels(tree)
    return {
        n: str(n.name) if labels[n] == str(n.name) else labels[n] + " | " + str(n.name)
        for n in tree.leaves()
    }


def _fixed(row, data):
    return (
        not str(data.manifest.get("calibration_policy", "hard")).startswith("PAML")
        and row["age_min"] == row["age_max"]
    )


def _age_extent(data, *, include_bounds=False):
    values = [0.0]
    for row in [*data.events.values(), *data.species_rows.values()]:
        if row["estimated_age"] is not None:
            values.append(row["estimated_age"])
        if "age" in row:
            values.append(float(row["age"]))
        if include_bounds:
            values.append(row["age_max"])
        if row.get("interval_upper") is not None:
            values.append(row["interval_upper"])
        if row.get("input_interval_upper") is not None:
            values.append(row["input_interval_upper"])
    values.extend(r["age_max"] for r in data.species_rows.values() if data.dated)
    return max(values) or 1.0


def _uses_reference_species(data):
    return data.dated and any(
        row["estimated_age"] is None for row in data.species_rows.values()
    )


def _shared_ages(data):
    return sorted(
        {
            r["estimated_age"]
            for r in data.events.values()
            if r["event_type"] == "speciation"
        }
    )


def _draw_species_evidence(ax, row, y, sid, *, soft=False):
    if row.get("input_interval_lower") is not None:
        line = ax.hlines(
            y + 0.16,
            row["input_interval_lower"],
            row["input_interval_upper"],
            color="#94a3b8",
            linewidth=4,
            zorder=1,
        )
        line.set_gid("species-input-interval:" + sid)
        ax.scatter([float(row["age"])], [y + 0.16], color="#64748b", s=12, zorder=3)
    if row["age_min"] != row["age_max"]:
        line = ax.hlines(
            y + 0.3,
            row["age_min"],
            row["age_max"],
            color="#475569",
            linewidth=1,
            linestyles=":",
            zorder=1,
        )
        line.set_gid(("species-prior-range:" if soft else "species-hard-range:") + sid)
    if (
        row.get("interval_lower") is not None
        and row["interval_upper"] > row["interval_lower"]
    ):
        line = ax.hlines(
            y,
            row["interval_lower"],
            row["interval_upper"],
            color=EVENT_STYLE["speciation"][0],
            alpha=0.45,
            linewidth=4,
            zorder=1,
        )
        line.set_gid("species-estimated-interval:" + sid)


def _draw_node_intervals(ax, row, y, sid, gene, is_leaf, data):
    if not gene:
        if not is_leaf:
            _draw_species_evidence(
                ax,
                row,
                y,
                sid,
                soft=str(data.manifest.get("calibration_policy", "")).startswith(
                    "PAML"
                ),
            )
        return
    lower, upper = row["interval_lower"], row["interval_upper"]
    if lower is not None and upper > lower:
        interval = ax.hlines(
            y, lower, upper, color=INTERVAL, alpha=0.45, linewidth=4.0, zorder=1
        )
        interval.set_gid("age-interval:" + sid)


def _result_node_label(label, row, species_names, age, gene, dated):
    if gene:
        label += " / " + species_names.get(row.get("species_event_id"), "unmapped")
    if dated:
        label += f"  {age:.4g}"
        if not gene and row.get("estimated_age") is None:
            label += " [input only]"
        if (
            not gene
            and row.get("estimation_status") == "calibration-only"
            and row["age_min"] != row["age_max"]
        ):
            label += " [constraint only]"
    return label


def _species_point_style(ax, row, age, y, color, fill, marker, *, reference_species):
    if not reference_species:
        return age, color, fill, marker
    if row["estimated_age"] is None:
        return age, "#64748b", "white", "D"
    point_age = row["estimated_age"]
    ax.plot([age, point_age], [y, y], color="#94a3b8", linestyle=":", linewidth=0.7)
    return point_age, color, fill, marker


def _draw_tree_panel(ax, data, *, gene, font_size, extent=None):
    tree = data.gene if gene else data.species
    reference_species = not gene and _uses_reference_species(data)
    index = data.gene_index if gene else data.species_index
    x, y, leaves, tip_labels, _ = _layout(
        tree, None if gene else _species_tip_labels(tree)
    )
    labels = _internal_labels(tree, "g" if gene else "s")
    species_labels = _species_labels(data.species)
    species_names = {
        data.species_index.clade_id_for_node(n): label
        for n, label in species_labels.items()
    }
    rows = data.events if gene else data.species_rows
    if data.dated:
        x = {
            n: float(rows[index.clade_id_for_node(n)]["age"])
            if reference_species
            else rows[index.clade_id_for_node(n)]["estimated_age"]
            for n in tree.traverse()
        }
        for age in _shared_ages(data):
            line = ax.axvline(
                age, color="#a8bfd3", linestyle="--", linewidth=0.8, zorder=0
            )
            line.set_gid("shared-speciation-age")
    for node in tree.traverse():
        if node is not tree:
            (line,) = ax.plot(
                [x[node.up], x[node.up], x[node]],
                [y[node.up], y[node], y[node]],
                color="#94a3b8" if reference_species else INK,
                linewidth=1.05,
                solid_capstyle="round",
                zorder=2,
            )
            line.set_gid("branch:" + index.clade_id_for_node(node))
        if data.dated:
            sid = index.clade_id_for_node(node)
            _draw_node_intervals(ax, rows[sid], y[node], sid, gene, node.is_leaf, data)
    span = extent if data.dated else max(x.values(), default=1.0) or 1.0
    label_x = -0.035 * span if data.dated else max(x.values()) + 0.05 * span
    for node in tree.traverse():
        row = rows[index.clade_id_for_node(node)] if rows else {}
        kind = row.get("event_type", "leaf" if node.is_leaf else "speciation")
        if gene and not row.get("species_event_id"):
            kind = "unresolved"
        color, marker, _ = EVENT_STYLE[kind]
        size = 15 if node.is_leaf else 38
        if data.dated and not gene and not node.is_leaf and _fixed(row, data):
            marker = "D"
            fill = "white"
        else:
            fill = color
        point_age, color, fill, marker = _species_point_style(
            ax,
            row,
            x[node],
            y[node],
            color,
            fill,
            marker,
            reference_species=reference_species,
        )
        point = ax.scatter(
            [point_age],
            [y[node]],
            s=size,
            marker=marker,
            facecolors=fill,
            edgecolors=color,
            linewidths=1.1,
            zorder=4,
        )
        point.set_gid("event:" + kind + ":" + index.clade_id_for_node(node))
        if node.is_leaf:
            if not data.dated and x[node] < max(x.values()):
                ax.plot(
                    [x[node], label_x],
                    [y[node], y[node]],
                    color="#d5dce5",
                    linestyle=":",
                    linewidth=0.7,
                )
            text = ax.text(
                label_x,
                y[node],
                tip_labels[node],
                va="center",
                ha="left",
                fontsize=font_size,
                color=INK,
                linespacing=1.35,
            )
            text.set_gid("tip-label:" + str(node.name))
            continue
        label = _result_node_label(
            labels[node], row, species_names, point_age, gene, data.dated
        )
        annotation = ax.annotate(
            _wrap(label, 32),
            (point_age, y[node]),
            xytext=(-5, 5) if data.dated else (5, 5),
            ha="right" if data.dated else "left",
            textcoords="offset points",
            fontsize=max(6, font_size - 1),
            color=color,
            va="bottom",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=0.2),
        )
        annotation.set_gid("node-label:" + index.clade_id_for_node(node))
    ymax = max(y.values(), default=0.0)
    ax.set_ylim(ymax + 0.65, -0.85)
    if data.dated:
        ax.set_xlim(span * 1.05, -span * 0.38)
        ticks = MaxNLocator(nbins=6).tick_values(0, span)
        ax.set_xticks([t for t in ticks if 0 <= t <= span])
        ax.tick_params(axis="x", labelsize=font_size, colors=MUTED, length=3)
        ax.spines["bottom"].set_color("#ccd5df")
    else:
        ax.set_xlim(-0.06 * span, span * 1.8)
        ax.set_xticks([])
        ax.spines["bottom"].set_visible(False)
    ax.set_yticks([])
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    return labels


def _legend(data):
    kinds = {r["event_type"] for r in data.events.values()} - {"leaf"}
    if any(not r.get("species_event_id") for r in data.events.values()):
        kinds.add("unresolved")
    handles = [
        Line2D(
            [], [], linestyle="", marker=marker, color=color, markersize=5, label=label
        )
        for kind, (color, marker, label) in EVENT_STYLE.items()
        if kind in kinds
    ]
    if data.dated:
        if any(
            _fixed(data.species_rows[data.species_index.clade_id_for_node(n)], data)
            for n in data.species.traverse()
            if not n.is_leaf
        ):
            handles.append(
                Line2D(
                    [],
                    [],
                    linestyle="",
                    marker="D",
                    markerfacecolor="white",
                    color=EVENT_STYLE["speciation"][0],
                    markersize=5,
                    label="Fixed calibration",
                )
            )
        if any(
            r["interval_lower"] is not None
            and r["interval_upper"] > r["interval_lower"]
            for r in data.events.values()
        ):
            handles.append(
                Line2D(
                    [],
                    [],
                    color=INTERVAL,
                    linewidth=4,
                    alpha=0.55,
                    label="Age interval",
                )
            )
        if any(
            r.get("input_interval_lower") is not None
            for r in data.species_rows.values()
        ):
            handles.append(
                Line2D(
                    [],
                    [],
                    color="#94a3b8",
                    linewidth=4,
                    label="External species interval",
                )
            )
        if any(r["age_min"] != r["age_max"] for r in data.species_rows.values()):
            handles.append(
                Line2D(
                    [], [], color="#475569", linestyle=":", label="Species hard range"
                )
            )
        if any(
            r.get("interval_lower") is not None
            and r["interval_upper"] > r["interval_lower"]
            for r in data.species_rows.values()
        ):
            handles.append(
                Line2D(
                    [],
                    [],
                    color=EVENT_STYLE["speciation"][0],
                    linewidth=4,
                    alpha=0.45,
                    label="Species ensemble percentile"
                    if data.manifest.get("input_ensemble")
                    else "Species estimated interval",
                )
            )
        handles.append(
            Line2D([], [], color="#a8bfd3", linestyle="--", label="Shared species age")
        )
    return handles


def _diagnostic_lines(data):
    manifest = data.manifest
    interval = str(manifest.get("uncertainty", "none"))
    if interval in {"none", "not-requested"}:
        interval_text = "Not requested; points only"
    else:
        level = manifest.get("interval_level")
        interval_text = (
            f"{100 * float(level):g}% " if level is not None else ""
        ) + interval.replace("-", " ")
    methods = {
        "marginal-lognormal": "Marginal lognormal clock",
        "sequence-marginal-quadratic": "Marginal sequence likelihood (quadratic approximation)",
        "sequence-empirical-bayes-map": "Conditional sequence MAP",
        "mcmctree-posterior-mean": "MCMCTree posterior mean",
    }
    method = manifest.get("method", "unreported")
    lines = [f"Method: {methods.get(method, method)}", f"Intervals: {interval_text}"]
    policy = manifest.get("calibration_policy", "unreported")
    lines.append(
        "Calibration: "
        + ("Hard bounds on all events" if policy == "hard-all-events" else policy)
    )
    groups = Counter(
        r.get("shared_age_id")
        for r in data.events.values()
        if r["event_type"] == "speciation"
    )
    lines.append(
        f"Shared species events: {len(groups)} ({sum(v > 1 for v in groups.values())} represented by multiple gene nodes)"
    )
    if manifest.get("experimental_native_estimator"):
        lines.extend(
            [
                "Exploratory native RADTE estimate; general interval coverage is not established.",
                "Native intervals are conditional estimates, not MCMC posterior intervals.",
                "Conditioned on the species calibration and reconciliation assumptions.",
            ]
        )
        if "profile" in interval and manifest.get("sequence_model"):
            lines.append(
                "Profile intervals also condition on the fitted substitution model."
            )
    descriptions = {
        "rate_sd_changes_uncertainty_only_in_tree_mode": "Fixed rate SD affects interval width, not point ages (tree-only input).",
        "conditional_on_input_branch_lengths_and_root_split": "Conditional on the supplied branch lengths and root split.",
    }
    species_names = {
        "S:" + data.species_index.clade_id_for_node(node): label
        for node, label in _species_labels(data.species).items()
    }
    for diagnostic in manifest.get("diagnostics", []):
        for prefix, message in (
            (
                "profile_multiple_local_optima:",
                "Profile fits found multiple local optima at ",
            ),
            ("profile_nonunique_age_optimum:", "Profile fits found nonunique ages at "),
        ):
            if str(diagnostic).startswith(prefix):
                key = str(diagnostic)[len(prefix) :]
                descriptions[str(diagnostic)] = (
                    message + species_names.get(key, key[:18]) + "."
                )

        lines.append(
            descriptions.get(
                str(diagnostic), "Diagnostic: " + str(diagnostic).replace("_", " ")
            )
        )
    external = sorted(
        {
            (
                r.get("input_interval_kind"),
                r.get("input_interval_level"),
                r.get("input_interval_source"),
            )
            for r in data.species_rows.values()
            if r.get("input_interval_lower") is not None
        }
    )
    for kind, level, source in external:
        lines.append(
            f"External species interval: {100 * level:g}% {kind}; {source} (display only)."
        )
    return lines


def _draw_age_summary(ax, data, labels, font_size):
    rows = [
        (n, data.events[data.gene_index.clade_id_for_node(n)])
        for n in data.gene.traverse()
        if data.events[data.gene_index.clade_id_for_node(n)]["event_type"]
        == "duplication"
    ]
    rows.sort(key=lambda item: item[1]["estimated_age"], reverse=True)
    soft = str(data.manifest.get("calibration_policy", "")).startswith("PAML")
    ax.set_title(
        "Duplication ages and " + ("original" if soft else "allowed") + " ranges",
        loc="left",
        fontsize=font_size + 1,
        color=INK,
        pad=12,
    )
    if not rows:
        ax.text(
            0,
            0.6,
            "No duplication events in this result.",
            transform=ax.transAxes,
            color=MUTED,
            fontsize=font_size,
        )
        ax.axis("off")
        return
    for y, (node, row) in enumerate(rows):
        ax.hlines(
            y, row["age_min"], row["age_max"], color="#d4dce6", linewidth=6, zorder=1
        )
        if row["interval_lower"] is not None:
            ax.hlines(
                y,
                row["interval_lower"],
                row["interval_upper"],
                color=INTERVAL,
                linewidth=3,
                zorder=2,
            )
        age = row["estimated_age"]
        ax.scatter(
            [age], [y], marker="s", color=EVENT_STYLE["duplication"][0], s=24, zorder=3
        )
        tolerance = max(1.0, row["age_max"]) * 1e-6
        at_bound = (
            min(abs(age - row["age_min"]), abs(age - row["age_max"])) <= tolerance
        )
        label = labels[node] + (" *" if at_bound else "")
        if at_bound:
            ax.annotate(
                "bound",
                (age, y),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=6,
                color=MUTED,
            )
        ax.text(
            -0.03,
            y,
            label,
            transform=ax.get_yaxis_transform(),
            ha="right",
            va="center",
            fontsize=font_size,
        )
    extent = _age_extent(data, include_bounds=True)
    ax.set_xlim(extent * 1.04, -extent * 0.04)
    ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.set_yticks([])
    ax.xaxis.set_major_locator(MaxNLocator(4))
    ax.tick_params(labelsize=font_size - 1, colors=MUTED)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color("#ccd5df")
    ax.set_xlabel(
        "Gray: original ranges (not all used by PAML)"
        if soft
        else "Gray: allowed ranges; * estimate at a boundary",
        fontsize=font_size - 1,
        color=MUTED,
    )


def build_result_figure(data, *, width=None, height=None, font_size=8.0, time_unit=""):
    for label, value in (
        ("figure width", width),
        ("figure height", height),
        ("font size", font_size),
    ):
        if value is not None and (not math.isfinite(float(value)) or float(value) <= 0):
            raise ValueError(f"{label} must be finite and positive.")
    gene_rows, species_rows = (
        _layout(data.gene)[-1],
        _layout(data.species, _species_tip_labels(data.species))[-1],
    )
    font_scale = max(1.0, font_size / 8)
    width = 13.0 * font_scale if width is None else width
    if width < 8:
        raise ValueError(
            "Reconciliation/RADTE reports require --figure-width of at least 8 inches."
        )
    counts = Counter(r["event_type"] for r in data.events.values())
    if data.dated:
        gene_height = max(1.8, 0.26 * gene_rows) * font_scale
        species_height = max(1.35, 0.26 * species_rows) * font_scale
        diagnostic_lines = [_wrap(line, 78) for line in _diagnostic_lines(data)]
        footer_height = (
            max(
                1.4,
                0.24 * counts["duplication"],
                0.17 * sum(x.count("\n") + 1 for x in diagnostic_lines),
            )
            * font_scale
        )
        natural_height = 3.7 + gene_height + species_height + footer_height
    else:
        natural_height = (
            2.4 + max(1.8, 0.26 * max(gene_rows, species_rows)) * font_scale
        )
    legend_extra = max(0, math.ceil(len(_legend(data)) / 5) - 1) * 0.3 * font_scale
    natural_height += legend_extra
    height = natural_height if height is None else height
    if height < natural_height:
        raise ValueError(
            f"Figure height is too small for the labels; use at least {natural_height:.1f} inches."
        )
    with matplotlib.rc_context(
        {
            "font.family": "DejaVu Sans",
            "font.size": font_size,
            "svg.fonttype": "none",
            "svg.hashsalt": "nwkit-results",
        }
    ):
        fig = plt.figure(figsize=(width, height), facecolor="white")
        try:
            fig.text(
                0.065,
                1 - 0.25 / height,
                "NWKIT  /  " + ("RADTE" if data.dated else "RECONCILIATION"),
                fontsize=8,
                color=MUTED,
                va="top",
                weight="bold",
            )
            fig.text(
                0.065,
                1 - 0.48 / height,
                "Shared speciation ages"
                if data.dated
                else "Gene events and species mapping",
                fontsize=18,
                weight="bold",
                color=INK,
                va="top",
            )
            sources = ", ".join(
                sorted(
                    {
                        str(r.get("event_source", "unreported"))
                        for r in data.events.values()
                    }
                )
            )
            subtitle = f"{counts['leaf']} gene tips  |  {len(list(data.species.leaves()))} species  |  {counts['duplication']} duplications  |  source: {sources}"
            fig.text(
                0.065,
                1 - 0.83 / height,
                subtitle,
                fontsize=font_size,
                color=MUTED,
                va="top",
            )
            fig.legend(
                handles=_legend(data),
                loc="upper left",
                bbox_to_anchor=(0.06, 1 - 1.02 / height),
                ncol=5,
                frameon=False,
                fontsize=font_size,
                handlelength=1.7,
                columnspacing=1.7,
            )
            if data.dated:
                extent = _age_extent(data)
                footer_bottom = 0.6 / height
                footer_top = footer_bottom + footer_height / natural_height
                species_bottom = footer_top + 0.7 / height
                species_top = species_bottom + species_height / natural_height
                gene_bottom = species_top + 0.45 / height
                gene_top = 1 - (1.6 + legend_extra) / height
                gene_ax = fig.add_axes(
                    [0.07, gene_bottom, 0.87, gene_top - gene_bottom]
                )
                species_ax = fig.add_axes(
                    [0.07, species_bottom, 0.87, species_top - species_bottom]
                )
                labels = _draw_tree_panel(
                    gene_ax, data, gene=True, font_size=font_size, extent=extent
                )
                _draw_tree_panel(
                    species_ax, data, gene=False, font_size=font_size, extent=extent
                )
                gene_ax.set_title(
                    "Dated gene tree",
                    loc="left",
                    fontsize=font_size + 2,
                    color=INK,
                    pad=12,
                )
                gene_ax.tick_params(labelbottom=False)
                species_ax.set_title(
                    "Input species tree (gray); estimated ages on matching rows"
                    if _uses_reference_species(data)
                    else "Species tree",
                    loc="left",
                    fontsize=font_size + 2,
                    color=INK,
                    pad=12,
                )
                species_ax.set_xlabel(
                    "Age (" + (time_unit or "input time units") + ")",
                    color=MUTED,
                    fontsize=font_size,
                )
                age_ax = fig.add_axes(
                    [0.085, footer_bottom, 0.32, footer_top - footer_bottom]
                )
                _draw_age_summary(age_ax, data, labels, font_size)
                diag_ax = fig.add_axes(
                    [0.46, footer_bottom, 0.48, footer_top - footer_bottom]
                )
                diag_ax.axis("off")
                diag_ax.set_title(
                    "Inference and interval diagnostics",
                    loc="left",
                    fontsize=font_size + 1,
                    color=INK,
                    pad=12,
                )
                diag_ax.text(
                    0,
                    1,
                    "\n".join(diagnostic_lines),
                    va="top",
                    fontsize=font_size,
                    color=MUTED,
                    linespacing=1.6,
                )
            else:
                bottom, top = 0.9 / height, 1 - 1.6 / height
                gene_ax = fig.add_axes([0.065, bottom, 0.55, top - bottom])
                species_ax = fig.add_axes([0.70, bottom, 0.27, top - bottom])
                _draw_tree_panel(gene_ax, data, gene=True, font_size=font_size)
                _draw_tree_panel(species_ax, data, gene=False, font_size=font_size)
                gene_ax.set_title(
                    "Gene tree", loc="left", fontsize=font_size + 2, color=INK, pad=12
                )
                species_ax.set_title(
                    "Input species tree (gray); estimated ages on matching rows"
                    if _uses_reference_species(data)
                    else "Species tree",
                    loc="left",
                    fontsize=font_size + 2,
                    color=INK,
                    pad=12,
                )
                unmapped = sum(
                    not r.get("species_event_id") for r in data.events.values()
                )
                note = f"Node labels: gene node / mapped species node.  Unmapped: {unmapped}; unresolved events: {counts['unresolved']}."
                fig.text(0.065, 0.48 / height, note, color=MUTED, fontsize=font_size)
                fig.text(
                    0.065,
                    0.25 / height,
                    "Topology view: branch lengths are omitted. Events are read from the supplied result, never re-inferred.",
                    color=MUTED,
                    fontsize=font_size,
                )
            return fig
        except BaseException:
            plt.close(fig)
            raise


def save_result_figure(data, path, *, image_format="auto", **options):
    fmt = figure_format(path, image_format)
    fig = build_result_figure(data, **options)
    try:
        with matplotlib.rc_context(
            {"svg.fonttype": "none", "svg.hashsalt": "nwkit-results"}
        ):
            fig.savefig(path, format=fmt, dpi=180, bbox_inches="tight", pad_inches=0.18)
    finally:
        plt.close(fig)


def draw_saved_results(args):
    if not args.species_tree or args.species_tree == "-":
        raise ValueError("Result plotting requires --species-tree PATH.")
    _validate_report_draw_options(args)
    fmt = figure_format(args.outfile, args.image_format)
    inputs = [("--species-tree", args.species_tree)]
    if args.radte_prefix:
        if args.infile != "-":
            raise ValueError("--radte-prefix supplies the dated tree; omit --infile.")
        inputs += list(radte_plot_protected_paths(args.radte_prefix).items())
    else:
        if args.infile == "-" or args.reconciliation == "-":
            raise ValueError(
                "Saved reconciliation plotting requires file paths for --infile and --reconciliation."
            )
        inputs += [("--infile", args.infile), ("--reconciliation", args.reconciliation)]
    validate_outputs_do_not_replace_inputs(
        inputs, [("--outfile", args.outfile)], label="Result figure"
    )
    if args.radte_prefix:
        data = load_dating_plot_data(
            args.radte_prefix,
            args.species_tree,
            species_rooted=args.species_tree_rooted,
            species_format=args.species_tree_format,
        )
    else:
        data = reconciliation_plot_data(
            read_tree(
                args.infile,
                args.format,
                args.quoted_node_names,
                rooted=args.input_rooted,
            ),
            read_tree(
                args.species_tree,
                args.species_tree_format,
                args.quoted_node_names,
                rooted=args.species_tree_rooted,
            ),
            read_result_table(args.reconciliation),
        )
    with output_transaction([args.outfile], follow_symlinks=False) as staged:
        save_result_figure(
            data,
            staged[args.outfile],
            image_format=fmt,
            width=args.figure_width,
            height=args.figure_height,
            font_size=args.font_size,
            time_unit=args.branch_length_unit
            or data.manifest.get("options", {}).get("branch_length_unit", ""),
        )


def _validate_report_draw_options(args):
    # Report layouts have two trees and a diagnostic panel. Do not silently
    # accept single-tree layout/trait options that cannot affect these reports.
    from nwkit.cli import pdraw

    defaults = pdraw.parse_args([])
    supported = {
        "command",
        "handler",
        "audit",
        "debug",
        "infile",
        "outfile",
        "format",
        "quoted_node_names",
        "input_rooted",
        "species_tree",
        "species_tree_rooted",
        "species_tree_format",
        "reconciliation",
        "radte_prefix",
        "figure_width",
        "figure_height",
        "font_size",
        "branch_length_unit",
        "image_format",
    }
    for action in pdraw._actions:
        if action.dest not in supported and hasattr(args, action.dest):
            if getattr(args, action.dest) != getattr(defaults, action.dest, None):
                raise ValueError(
                    f"--{action.dest.replace('_', '-')} is a single-tree option, unavailable for result reports."
                )
