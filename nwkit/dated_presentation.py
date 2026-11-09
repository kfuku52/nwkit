"""Publication panels for saved dated trees, parsed and validated by NWKIT.

Geological boundaries: ICS International Chronostratigraphic Chart 2026/06,
https://stratigraphy.org/ICSchart/ChronostratChart2026-06.pdf . Precambrian is
shown as a single broad interval. Backgrounds alternate between two light greys.
"""

import csv
import hashlib
import json
import math
import re
from contextlib import ExitStack
from pathlib import Path
from typing import Any

STATUS = ("single", "duplicated", "fragmented", "missing")
STATUS_LABELS = ("Single-copy", "Duplicated", "Fragmented", "Missing")
STATUS_COLOURS = ("#000000", "#B22222", "#666666", "#CCCCCC")
GEOLOGICAL_DATA = Path(__file__).with_name("data_geology") / "geological_periods.tsv"
GEOLOGICAL_SOURCE = "https://stratigraphy.org/ICSchart/ChronostratChart2026-06.pdf"


def species_key(value):
    return str(value).strip().replace(" ", "_")


def read_metadata_table(path, required=()):
    """Reject ambiguous TSVs before DictReader can silently discard columns."""
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if (
            not fields
            or not all(fields)
            or len(fields) != len(set(fields))
            or set(required) - set(fields)
        ):
            raise ValueError(
                f"Invalid metadata header in {path}: columns must be unique and include required fields."
            )
        rows = list(reader)
    if any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError(f"Metadata row width does not match its header in {path}.")
    return rows


def geological_intervals(max_age, data_path=None):
    if not math.isfinite(max_age) or not 0 < max_age <= 4567:
        raise ValueError("Geological background requires ages between 0 and 4567 Ma.")
    with Path(data_path or GEOLOGICAL_DATA).open(
        encoding="utf-8", newline=""
    ) as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    return [
        {
            **row,
            "young_Ma": float(row["young_Ma"]),
            "old_Ma": min(max_age, float(row["old_Ma"])),
        }
        for row in records
        if float(row["young_Ma"]) < max_age
    ]


def interval_label(nodes):
    intervals = {
        (str(node.props.get("age_ci_kind", "")).upper(), node.props.get("age_ci_level"))
        for node in nodes
        if node.props.get("age_ci_low") is not None
    }
    if not intervals:
        return None
    labels = []
    for kind, level in sorted(intervals, key=str):
        percentage = f"{100 * float(level):g}% " if level is not None else ""
        description = {
            "HPD": "highest posterior density intervals",
            "ETI": "equal-tailed credible intervals",
            "EQUAL-TAIL": "equal-tailed credible intervals",
        }.get(kind, "credible intervals")
        labels.append(percentage + description)
    return "; ".join(labels)


def read_branch_annotations(path, nodes, ys):
    """Map presentation symbols to exact stem branches, never estimate dates."""
    from matplotlib.colors import is_color_like

    clades = {frozenset(node.leaf_names()): node for node in nodes}
    required = {"descendant_species", "label", "symbol"}
    optional = {"event_id", "colour", "branch_fraction"}
    events, identifiers, positions = [], set(), set()
    legend_styles: dict[str, Any] = {}
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if (
            len(fields) != len(set(fields))
            or required - set(fields)
            or set(fields) - required - optional
        ):
            raise ValueError(
                "Branch annotations require descendant_species,label,symbol and only documented optional columns."
            )
        for index, row in enumerate(reader, 1):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(
                    "Branch annotation row width does not match its header."
                )
            row = {key: value.strip() for key, value in row.items()}
            names = [name.strip() for name in row["descendant_species"].split(",")]
            node = clades.get(frozenset(names))
            if (
                not all(names)
                or len(names) != len(set(names))
                or node is None
                or node.is_root
            ):
                raise ValueError(
                    "Branch annotation must select an exact non-root clade or tip."
                )
            identifier = row.get("event_id") or f"event-{index}"
            if (
                not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", identifier)
                or identifier in identifiers
            ):
                raise ValueError("Branch event_id must be a unique safe identifier.")
            label, symbol, colour = (
                row["label"],
                row["symbol"],
                row.get("colour") or "#202020",
            )
            if (
                not label
                or symbol not in {"^", "v", "o", "s", "D", "x", "+", "*"}
                or not is_color_like(colour)
            ):
                raise ValueError("Invalid branch annotation label, symbol or colour.")
            fraction = float(row.get("branch_fraction") or 0.5)
            child, parent = float(node.props["age"]), float(node.up.props["age"])
            if not math.isfinite(fraction) or not 0 < fraction < 1 or parent <= child:
                raise ValueError(
                    "Branch fraction must be between 0 and 1 on a positive-length branch."
                )
            position = (frozenset(names), fraction)
            if position in positions:
                raise ValueError("Branch annotations share the same display position.")
            if label in legend_styles and legend_styles[label] != (symbol, colour):
                raise ValueError(
                    "A branch legend label has conflicting symbols or colours."
                )
            identifiers.add(identifier)
            positions.add(position)
            legend_styles[label] = (symbol, colour)
            events.append(
                {
                    "event_id": identifier,
                    "label": label,
                    "symbol": symbol,
                    "colour": colour,
                    "descendant_species": sorted(names),
                    "branch_fraction": fraction,
                    "branch_child_age_Ma": child,
                    "branch_parent_age_Ma": parent,
                    "display_position_Ma": child + fraction * (parent - child),
                    "y": ys[node],
                    "position_interpretation": "Graphical position along the branch, not an estimated event date.",
                }
            )
    return events


def _prepare_tree_view(
    age_clades,
    branch_annotations,
    busco_summary,
    figure_width,
    font_size,
    geological_background,
    infile,
    layout_report,
    node_ages,
    outfile,
    tip_annotations,
    tip_order,
):
    """Prepare tree view."""
    from nwkit.file_paths import (
        validate_distinct_output_paths,
        validate_outputs_do_not_replace_inputs,
    )
    from nwkit.time_tree import (
        infer_node_ages_from_branch_lengths,
        prepare_time_tree_annotations,
    )
    from nwkit.util import read_tree

    outfile = Path(outfile)
    if geological_background not in {"none", "period"} or node_ages not in {
        "none",
        "root",
        "all",
    }:
        raise ValueError("Invalid geological background or node age label mode.")
    inputs = [("tree", infile)] + [
        (name, value)
        for name, value in [
            ("busco_summary", busco_summary),
            ("tip_order", tip_order),
            ("tip_annotations", tip_annotations),
            ("branch_annotations", branch_annotations),
            ("age_clades", age_clades),
        ]
        if value is not None
    ]
    outputs = [("plot", outfile)] + (
        [("layout_report", layout_report)] if layout_report else []
    )
    validate_distinct_output_paths(outputs)
    validate_outputs_do_not_replace_inputs(inputs, outputs)
    if outfile.suffix.lower() not in {".pdf", ".svg", ".png"}:
        raise ValueError("Presentation plots support PDF, SVG and PNG.")
    if min(figure_width, font_size) <= 0 or not all(
        math.isfinite(value) for value in [figure_width, font_size]
    ):
        raise ValueError("Figure width and font size must be finite and positive.")
    tree = read_tree(str(infile), "auto", True, rooted="yes")
    annotations = prepare_time_tree_annotations(tree)
    if not annotations["credible_interval_count"]:
        infer_node_ages_from_branch_lengths(tree)
    leaves = list(tree.leaves())
    species = [leaf.name for leaf in leaves]
    if not all(species) or len({species_key(name) for name in species}) != len(species):
        raise ValueError(
            "Tree tip names must be nonempty and unique, including normalized species IDs."
        )
    order = species
    if tip_order:
        order = [
            row["species_id"] for row in read_metadata_table(tip_order, ["species_id"])
        ]
        if len(order) != len(species) or set(order) != set(species):
            raise ValueError("Tip order must contain each tree species exactly once.")
        ranks = {name: rank for rank, name in enumerate(order)}
        for node in tree.traverse(strategy="postorder"):
            if not node.is_leaf:
                node.children.sort(
                    key=lambda child: min(ranks[name] for name in child.leaf_names())
                )
        if list(tree.leaf_names()) != order:
            raise ValueError("Requested tip order is not compatible with the topology.")
    return inputs, leaves, order, outfile, outputs, species, tree


def _prepare_node_annotations(
    age_clades,
    branch_annotations,
    geological_background,
    geological_data,
    leaves,
    order,
    species,
    tip_annotations,
    tree,
):
    """Prepare node annotations."""
    from matplotlib.font_manager import FontProperties

    styles = {}
    if tip_annotations:
        from matplotlib.colors import is_color_like

        for row in read_metadata_table(tip_annotations, ["species_id"]):
            name = row["species_id"]
            if name not in species or name in styles:
                raise ValueError("Unknown or duplicate tip annotation.")
            row["font_weight"] = row.get("font_weight") or "normal"
            row["colour"] = row.get("colour") or "#202020"
            FontProperties(weight=row["font_weight"])
            if not is_color_like(row["colour"]):
                raise ValueError("Invalid tip annotation colour.")
            styles[name] = row
    selected = set()
    if age_clades:
        for row in read_metadata_table(age_clades, ["descendant_species"]):
            names = [name.strip() for name in row["descendant_species"].split(",")]
            clade = frozenset(names)
            if not all(names) or len(names) != len(clade) or clade in selected:
                raise ValueError(
                    "Age clades require unique nonempty species and no repeated clade."
                )
            selected.add(clade)
        actual = {
            frozenset(node.leaf_names()) for node in tree.traverse() if not node.is_leaf
        }
        if not selected <= actual:
            raise ValueError(
                "A requested age label is not an internal clade in the tree."
            )
    nodes = list(tree.traverse())
    count = len(species)
    ranks = {name: index for index, name in enumerate(order)}
    ys = {leaf: float(count - 1 - ranks[leaf.name]) for leaf in leaves}
    for node in tree.traverse(strategy="postorder"):
        if not node.is_leaf:
            ys[node] = (ys[node.children[0]] + ys[node.children[-1]]) / 2
    events = (
        read_branch_annotations(branch_annotations, nodes, ys)
        if branch_annotations
        else []
    )
    max_age = max(
        float(node.props.get("age_ci_high", node.props["age"])) for node in nodes
    )
    if geological_background == "period" and max_age > 4567:
        raise ValueError("Geological background requires ages between 0 and 4567 Ma.")
    max_age = (
        min(4567, math.ceil(max_age / 10) * 10)
        if geological_background == "period"
        else math.ceil(max_age / 10) * 10
    )
    if max_age <= 0:
        raise ValueError("A dated tree must span positive time.")
    periods = (
        geological_intervals(max_age, geological_data)
        if geological_background == "period"
        else []
    )
    return count, events, max_age, nodes, periods, selected, styles, ys


def _prepare_side_panels(
    busco_prefix,
    busco_results,
    busco_summary,
    figure_height,
    figure_width,
    font_size,
    geological_data,
    inputs,
    outputs,
    panel_loader,
    periods,
    row_spacing_points,
    species,
):
    """Prepare side panels."""
    from nwkit.file_paths import validate_outputs_do_not_replace_inputs

    metadata_inputs: list[Path] = []
    counts, dataset, dataset_sources = (
        panel_loader(
            busco_summary,
            species,
            busco_results,
            busco_prefix,
            input_paths=metadata_inputs,
        )
        if busco_summary
        else ({}, None, {})
    )
    inputs.extend(("BUSCO metadata", source) for source in metadata_inputs)
    if periods:
        inputs.append(
            ("geological boundaries", Path(geological_data or GEOLOGICAL_DATA))
        )
    validate_outputs_do_not_replace_inputs(inputs, outputs)
    spacing = (
        row_spacing_points
        if row_spacing_points is not None
        else max(font_size + 1, font_size * 1.12)
    )
    if not math.isfinite(spacing) or spacing < font_size:
        raise ValueError(
            "Row spacing must be finite and at least the font size, in points."
        )
    if figure_width < 3.6 or (
        figure_height is not None
        and (not math.isfinite(figure_height) or figure_height < 2.5)
    ):
        raise ValueError(
            "Presentation plots need a width of at least 3.6 and a height of at least 2.5 inches."
        )
    return counts, dataset, dataset_sources, spacing


def _measure_tree_labels(
    cleanup,
    counts,
    figure_height,
    figure_width,
    font_family,
    font_size,
    nodes,
    species,
    styles,
):
    """Measure tree labels."""
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties

    figure = plt.figure(figsize=(figure_width, figure_height or 3))
    cleanup.callback(plt.close, figure)
    renderer = figure.canvas.get_renderer()

    def text_width(text, style="normal", weight="normal"):
        properties = FontProperties(
            family=font_family, size=font_size, style=style, weight=weight
        )
        return (
            max(
                renderer.get_text_width_height_descent(line, properties, False)[0]
                for line in text.splitlines()
            )
            * 72
            / figure.dpi
        )

    width_points = figure_width * 72
    left_points, right_points, gap_points = 12, 10, 6
    label_points = (
        max(
            text_width(
                name.replace("_", " "),
                "italic",
                styles.get(name, {}).get("font_weight", "normal"),
            )
            for name in species
        )
        + 2
    )
    available = (
        width_points
        - left_points
        - right_points
        - label_points
        - gap_points * (2 if counts else 1)
    )
    bar_points = max(58, min(110, available * 0.27)) if counts else 0
    tree_points = available - bar_points
    if tree_points < 48:
        raise ValueError("Figure is too narrow for the tip labels and panels.")
    credible_label = interval_label(nodes)
    credible_labels = [credible_label] if credible_label else []
    if (
        credible_label
        and text_width(credible_label) + font_size * 2.15 > width_points - 24
    ):
        credible_labels = credible_label.split("; ")
    return (
        bar_points,
        credible_label,
        credible_labels,
        figure,
        gap_points,
        label_points,
        left_points,
        text_width,
        tree_points,
        width_points,
    )


def _measure_legend_layout(
    STATUS_LABELS,
    count,
    counts,
    credible_labels,
    events,
    figure,
    figure_height,
    figure_width,
    font_size,
    periods,
    show_geological_source,
    spacing,
    text_width,
    width_points,
):
    """Measure legend layout."""

    def pack_legend_labels(items):
        groups: list[list[str]] = []
        group: list[str] = []
        group_width = 0.0
        for label in items:
            size = text_width(label) + font_size * 1.35
            if size + font_size * 0.8 > width_points - 24:
                raise ValueError(
                    "Legend label is too wide; shorten it or increase figure width."
                )
            if group and group_width + size + font_size * 1.7 > width_points - 24:
                groups.append(group)
                group, group_width = [], 0
            group_width += size + (font_size * 0.9 if group else 0)
            group.append(label)
        return groups + ([group] if group else [])

    status_labels = list(STATUS_LABELS) if counts else []
    base_groups = pack_legend_labels(credible_labels + status_labels)
    if len(base_groups) > 1 and credible_labels and status_labels:
        base_groups = pack_legend_labels(credible_labels) + pack_legend_labels(
            status_labels
        )
    base_groups.reverse()  # BUSCO below age intervals, as in the compact publication panel.
    event_labels: dict[str, Any] = {}
    for event in events:
        event_labels.setdefault(event["label"], event)
    event_groups = pack_legend_labels(event_labels)
    legend_rows = len(base_groups) + len(event_groups)
    row_heights = [
        font_size
        + 8
        + (max(len(label.splitlines()) for label in group) - 1) * font_size * 1.5
        for group in base_groups + event_groups
    ]
    source_text = "Geological periods: ICS 2026/06."
    if any(p["name"] == "Precambrian" for p in periods):
        source_text += " Precambrian shown as one interval."
    words = source_text.split()
    source_lines: list[str] = []
    for word in words:
        if (
            not source_lines
            or text_width(source_lines[-1] + " " + word) > width_points - 24
        ):
            source_lines.append(word)
        else:
            source_lines[-1] += " " + word
    source_text = "\n".join(source_lines)
    source_points = 20 if periods and show_geological_source else 0
    if source_points:
        source_points += (len(source_lines) - 1) * font_size * 1.5
    footer_points = (
        (48 if counts else font_size * 3 + 7) + source_points + sum(row_heights)
    )
    period_offset = 4
    header_points = max(
        30,
        (max(text_width(p["name"]) for p in periods) + period_offset + 8)
        if periods
        else 30,
    )
    height = (
        figure_height
        if figure_height is not None
        else max(2.5, (header_points + footer_points + (count + 0.2) * spacing) / 72)
    )
    figure.set_size_inches(figure_width, height)
    bottom, top = footer_points / (height * 72), 1 - header_points / (height * 72)
    if (top - bottom) * height < 1:
        if count > 6:
            raise ValueError(
                "Figure is too short for the geological names and tree; increase its height."
            )
    return (
        base_groups,
        bottom,
        event_groups,
        event_labels,
        height,
        legend_rows,
        period_offset,
        row_heights,
        source_points,
        source_text,
        top,
    )


def _draw_geological_axes(
    bottom,
    count,
    figure,
    font_size,
    gap_points,
    label_points,
    left_points,
    max_age,
    period_offset,
    periods,
    top,
    tree_points,
    width_points,
):
    """Draw geological axes."""
    from matplotlib.ticker import MaxNLocator

    axis = figure.add_axes(
        [
            left_points / width_points,
            bottom,
            tree_points / width_points,
            top - bottom,
        ]
    )
    label_left = left_points + tree_points + gap_points
    labels = figure.add_axes(
        [
            label_left / width_points,
            bottom,
            label_points / width_points,
            top - bottom,
        ],
        sharey=axis,
    )
    labels.set_xlim(0, 1)
    labels.axis("off")
    axis.set_ylim(-0.6, count - 0.4)
    axis.set_xlim(max_age, 0)
    axis.set_xticks(
        [
            value
            for value in MaxNLocator(
                nbins=max(2, min(5, int(tree_points / 38)))
            ).tick_values(0, max_age)
            if 0 <= value <= max_age
        ]
    )
    axis.set_yticks([])
    axis.spines[["left", "right", "top"]].set_visible(False)
    axis.spines["bottom"].set_linewidth(0.6)
    axis.xaxis.tick_bottom()
    axis.xaxis.set_label_position("bottom")
    axis.set_xlabel("Divergence time (Ma)")
    axis.tick_params(axis="x", labelsize=font_size, length=2.5, width=0.6)
    period_texts = []
    for period in periods:
        axis.axvspan(
            period["young_Ma"], period["old_Ma"], color=period["colour"], zorder=-5
        )
        text = axis.annotate(
            period["name"],
            xy=((period["young_Ma"] + period["old_Ma"]) / 2, 1),
            xycoords=("data", "axes fraction"),
            xytext=(0, period_offset),
            textcoords="offset points",
            rotation=90,
            ha="center",
            va="bottom",
            fontsize=font_size,
            annotation_clip=False,
        )
        text.set_gid("geological-period-" + period["name"])
        period_texts.append(text)
    if period_texts:
        # Spread crowded names within the header while retaining their
        # original time anchors. A short leader connects a shifted name.
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        ordered = list(reversed(period_texts))
        widths = [text.get_window_extent(renderer).width for text in ordered]
        gap = 2 * figure.dpi / 72
        left, right = axis.bbox.x0 - max(widths) / 2, axis.bbox.x1 + max(widths) / 2
        if sum(widths) + gap * (len(widths) - 1) > right - left:
            raise ValueError(
                "Geological period names need more space; increase figure width."
            )
        anchors = [axis.transData.transform((text.xy[0], 0))[0] for text in ordered]
        centres = []
        cursor = left
        for anchor, width in zip(anchors, widths, strict=True):
            centre = max(anchor, cursor + width / 2)
            centres.append(centre)
            cursor = centre + width / 2 + gap
        cursor = right
        for index in range(len(ordered) - 1, -1, -1):
            centres[index] = min(centres[index], cursor - widths[index] / 2)
            cursor = centres[index] - widths[index] / 2 - gap
        for text, centre, anchor in zip(ordered, centres, anchors, strict=True):
            offset = (centre - anchor) * 72 / figure.dpi
            text.set_position((offset, period_offset))
            if abs(offset) > 0.25:
                axis.annotate(
                    "",
                    xy=text.xy,
                    xycoords=("data", "axes fraction"),
                    xytext=(offset, period_offset - 1),
                    textcoords="offset points",
                    annotation_clip=False,
                    arrowprops={
                        "arrowstyle": "-",
                        "color": "#777777",
                        "lw": 0.4,
                        "shrinkA": 0,
                        "shrinkB": 0,
                    },
                )
    return axis, label_left, labels, period_texts


def _draw_branches_and_intervals(axis, events, font_size, nodes, ys):
    """Draw branches and intervals."""
    from matplotlib.lines import Line2D

    branches = []
    for node in nodes:
        x, y = float(node.props["age"]), ys[node]
        if not node.is_leaf:
            points = [(x, ys[child]) for child in node.children]
            branches.append(
                (
                    (x, min(point[1] for point in points)),
                    (x, max(point[1] for point in points)),
                )
            )
        if not node.is_root:
            branches.append(((float(node.up.props["age"]), y), (x, y)))
    branch_lines = []
    for index, (first, second) in enumerate(branches):
        line = axis.plot(
            [first[0], second[0]],
            [first[1], second[1]],
            color="#202020",
            lw=0.7,
            zorder=2,
        )[0]
        line.set_gid(f"tree-branch-{index}")
        branch_lines.append(line)
    interval_count = 0
    interval_lines = []
    for node in nodes:
        if node.props.get("age_ci_low") is not None:
            low, high = (
                float(node.props["age_ci_low"]),
                float(node.props["age_ci_high"]),
            )
            # Stop the bar at its exact bounds. Projecting line ends can
            # otherwise extend beyond a thinner endpoint marker.
            line = axis.plot(
                [low, high],
                [ys[node]] * 2,
                color="#D55E00",
                lw=1.2,
                solid_capstyle="butt",
                zorder=1,
            )[0]
            caps = axis.plot(
                [low, high],
                [ys[node]] * 2,
                linestyle="none",
                marker="|",
                markersize=3,
                markeredgewidth=1.2,
                color="#D55E00",
                zorder=1,
            )[0]
            line.set_gid(f"age-interval-{interval_count}")
            caps.set_gid(f"age-interval-cap-{interval_count}")
            interval_lines.extend([line, caps])
            interval_count += 1

    def event_marker(event, x=(), y=(), **kwargs):
        import matplotlib.patheffects as effects

        filled = event["symbol"] not in {"x", "+"}
        marker = Line2D(
            x,
            y,
            linestyle="none",
            marker=event["symbol"],
            markersize=font_size * 0.625,
            color=event["colour"],
            markerfacecolor=event["colour"],
            markeredgecolor="white" if filled else event["colour"],
            markeredgewidth=0.6 if filled else 1,
            **kwargs,
        )
        if not filled:
            marker.set_path_effects(
                [
                    effects.Stroke(linewidth=2.3, foreground="white"),
                    effects.Normal(),
                ]
            )
        return marker

    event_markers = []
    for event in events:
        marker = event_marker(
            event,
            [event["display_position_Ma"]],
            [event["y"]],
            zorder=5,
            clip_on=False,
        )
        marker.set_gid("branch-event-" + event["event_id"])
        axis.add_line(marker)
        event_markers.append(marker)
    return (
        branch_lines,
        branches,
        event_marker,
        event_markers,
        interval_count,
        interval_lines,
    )


def _draw_tip_and_age_labels(
    axis, font_size, labels, leaves, node_ages, nodes, selected, styles, ys
):
    """Draw tip and age labels."""
    tip_texts = []
    for leaf in leaves:
        style = styles.get(leaf.name, {})
        tip_texts.append(
            labels.text(
                0,
                ys[leaf],
                leaf.name.replace("_", " "),
                va="center",
                style="italic",
                weight=style.get("font_weight", "normal"),
                color=style.get("colour", "#202020"),
                fontsize=font_size,
            )
        )
    age_texts = []
    for node in nodes:
        if node.is_leaf or not (
            node_ages == "all"
            or (node_ages == "root" and node.is_root)
            or frozenset(node.leaf_names()) in selected
        ):
            continue
        text = axis.annotate(
            f"{float(node.props['age']):.1f}",
            xy=(float(node.props["age"]), ys[node]),
            xytext=(2, 3),
            textcoords="offset points",
            fontsize=font_size,
            zorder=4,
        )
        age_texts.append(text)
    return age_texts, tip_texts


def _place_age_labels(age_texts, axis, branches, event_markers, figure):
    """Place age labels."""
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    branch_pixels = [
        (axis.transData.transform(a), axis.transData.transform(b)) for a, b in branches
    ]
    event_boxes = [
        marker.get_window_extent(renderer).padded(1.2 * figure.dpi / 72)
        for marker in event_markers
    ]
    if any(
        box.overlaps(other)
        for index, box in enumerate(event_boxes)
        for other in event_boxes[index + 1 :]
    ):
        raise ValueError(
            "Branch symbols overlap; adjust branch_fraction or increase figure size."
        )
    occupied = list(event_boxes)
    for text in age_texts:
        placed = False
        positions = [
            (2, 3, "left"),
            (-2, 3, "right"),
            (2, -10, "left"),
            (-2, -10, "right"),
            (2, 12, "left"),
            (-2, 12, "right"),
            (2, -19, "left"),
            (-2, -19, "right"),
        ]
        # Compact rows may leave space only between branch levels. Search
        # nearby point offsets, retaining a leader for displaced labels.
        positions.extend(
            (dx, dy, align)
            for distance in range(1, 49)
            for dy in (distance, -distance)
            for dx, align in ((2, "left"), (-2, "right"))
        )
        for dx, dy, align in positions:
            text.set_position((dx, dy))
            text.set_ha(align)
            box = text.get_window_extent(renderer).padded(0.5 * figure.dpi / 72)
            crosses = any(
                (
                    min(a[0], b[0]) < box.x1
                    and max(a[0], b[0]) > box.x0
                    and min(a[1], b[1]) < box.y1
                    and max(a[1], b[1]) > box.y0
                )
                for a, b in branch_pixels
            )
            if (
                not crosses
                and not any(box.overlaps(other) for other in occupied)
                and box.x0 >= axis.bbox.x0
                and box.y0 >= axis.bbox.y0 - 12 * figure.dpi / 72
                and axis.bbox.contains(box.x1, box.y1)
            ):
                occupied.append(box)
                placed = True
                if abs(dy) > 10:
                    axis.annotate(
                        "",
                        xy=text.xy,
                        xytext=(dx, dy),
                        textcoords="offset points",
                        arrowprops={
                            "arrowstyle": "-",
                            "color": "#888888",
                            "lw": 0.35,
                            "shrinkA": 4,
                            "shrinkB": 1,
                        },
                    )
                break
        if not placed:
            raise ValueError(
                f"Age label {text.get_text()} Ma cannot fit without overlap; "
                "increase figure size or select fewer labelled clades."
            )
    return event_boxes


def _draw_count_panel(
    STATUS_COLOURS,
    axis,
    bar_points,
    bottom,
    counts,
    dataset,
    figure,
    font_size,
    gap_points,
    label_left,
    label_points,
    leaves,
    text_width,
    top,
    width_points,
    ys,
):
    """Draw count panel."""
    from matplotlib.ticker import MaxNLocator

    busco_axis_label = None
    bars = percent = None
    if counts:
        bar_left = label_left + label_points + gap_points
        bars = figure.add_axes(
            [
                bar_left / width_points,
                bottom,
                bar_points / width_points,
                top - bottom,
            ],
            sharey=axis,
        )
        total = sum(next(iter(counts.values())))
        for leaf in leaves:
            left = 0
            for value, colour in zip(counts[leaf.name], STATUS_COLOURS, strict=True):
                bars.barh(
                    ys[leaf],
                    value,
                    left=left,
                    height=0.62,
                    color=colour,
                    edgecolor="none",
                )
                left += value
        bars.set_xlim(0, total)
        bars.set_xticks(
            [
                value
                for value in MaxNLocator(
                    nbins=3, integer=True, steps=[1, 2, 5, 10]
                ).tick_values(0, total)
                if 0 <= value <= total
            ]
        )
        bars.set_yticks([])
        bars.spines[["left", "right", "top"]].set_visible(False)
        bars.tick_params(axis="x", labelsize=font_size, length=2.5, width=0.6)
        busco_axis_label = "Number of BUSCO genes" + (
            f"\n({dataset})" if dataset else ""
        )
        bars.set_xlabel(busco_axis_label)
        if max(text_width(line) for line in busco_axis_label.splitlines()) > bar_points:
            bars.xaxis.label.set_x(1)
            bars.xaxis.label.set_ha("right")
        percent = bars.twiny()
        percent.set_xlim(0, 100)
        percent.set_xticks([0, 50, 100])
        percent.set_xlabel("BUSCO genes (%)")
        percent.tick_params(axis="x", labelsize=font_size, length=2.5, width=0.6)
        percent.spines[["left", "right", "bottom"]].set_visible(False)
        percent.spines["top"].set_linewidth(0.6)
    return bars, busco_axis_label, percent


def _draw_legends(
    STATUS_COLOURS,
    STATUS_LABELS,
    base_groups,
    counts,
    credible_labels,
    event_groups,
    event_labels,
    event_marker,
    figure,
    font_size,
    height,
    left_points,
    periods,
    row_heights,
    show_geological_source,
    source_points,
    source_text,
    width_points,
):
    """Draw legends."""
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    handles = {
        label: Line2D([0], [0], color="#D55E00", lw=1.2, label=label)
        for label in credible_labels
    }
    if counts:
        handles.update(
            {
                name: Patch(facecolor=colour, label=name)
                for colour, name in zip(STATUS_COLOURS, STATUS_LABELS, strict=True)
            }
        )
    legend_groups = [[handles[label] for label in group] for group in base_groups] + [
        [event_marker(event_labels[label], label=label) for label in group]
        for group in event_groups
    ]
    legends = []
    legend_bottom = source_points
    for group, row_height in zip(legend_groups, row_heights, strict=True):
        legends.append(
            figure.legend(
                handles=group,
                loc="lower center",
                bbox_to_anchor=(0.5, legend_bottom / (height * 72)),
                ncol=len(group),
                frameon=False,
                fontsize=font_size,
                handlelength=1,
                columnspacing=0.9,
                handletextpad=0.35,
            )
        )
        legend_bottom += row_height
    if periods and show_geological_source:
        figure.text(
            left_points / width_points,
            5 / (height * 72),
            source_text,
            fontsize=font_size,
        )
    return legends


def _validate_presentation_bounds(
    bars, counts, figure, legends, period_texts, tip_texts
):
    """Validate presentation bounds."""
    import matplotlib

    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    period_boxes = [text.get_window_extent(renderer) for text in period_texts]
    if any(
        a.overlaps(b) for i, a in enumerate(period_boxes) for b in period_boxes[i + 1 :]
    ):
        raise ValueError("Geological period names overlap; increase figure width.")
    legend_boxes = [legend.get_window_extent(renderer) for legend in legends]
    if any(
        box.overlaps(other)
        for index, box in enumerate(legend_boxes)
        for other in legend_boxes[index + 1 :]
    ):
        raise ValueError(
            "Legend rows overlap; increase figure height or shorten multiline labels."
        )
    for text in figure.findobj(matplotlib.text.Text):
        if text.get_visible() and text.get_text():
            box = text.get_window_extent(renderer)
            if (
                box.x0 < 0
                or box.y0 < 0
                or box.x1 > figure.bbox.x1 + 0.5
                or box.y1 > figure.bbox.y1 + 0.5
            ):
                raise ValueError(
                    f"Text extends outside the figure: {text.get_text()!r}; increase its width or height."
                )
    if (
        counts
        and bars is not None
        and max(text.get_window_extent(renderer).x1 for text in tip_texts)
        + 4 * figure.dpi / 72
        > bars.bbox.x0
    ):
        raise ValueError(
            "Species labels extend into the BUSCO panel; increase figure width."
        )
    return renderer


def render_dated_tree(
    infile,
    outfile,
    *,
    busco_summary=None,
    panel_loader=None,
    category_names=STATUS,
    category_labels=STATUS_LABELS,
    category_colours=STATUS_COLOURS,
    busco_results=None,
    busco_prefix="busco_cds",
    geological_background="period",
    show_geological_source=False,
    figure_width=4.8,
    figure_height=None,
    row_spacing_points=None,
    font_family="Helvetica",
    font_size=8,
    tip_order=None,
    tip_annotations=None,
    branch_annotations=None,
    node_ages="none",
    age_clades=None,
    layout_report=None,
    geological_data=None,
):
    if len(category_names) != len(category_labels) or len(category_names) != len(
        category_colours
    ):
        raise ValueError("Side-panel categories, labels and colours must match")
    STATUS_LABELS = category_labels
    STATUS_COLOURS = category_colours
    if busco_summary is not None and panel_loader is None:
        raise ValueError("A side-panel source requires a panel_loader")

    import matplotlib

    matplotlib.use("Agg")

    import nwkit
    from nwkit.output_transaction import output_transaction

    inputs, leaves, order, outfile, outputs, species, tree = _prepare_tree_view(
        age_clades,
        branch_annotations,
        busco_summary,
        figure_width,
        font_size,
        geological_background,
        infile,
        layout_report,
        node_ages,
        outfile,
        tip_annotations,
        tip_order,
    )
    count, events, max_age, nodes, periods, selected, styles, ys = (
        _prepare_node_annotations(
            age_clades,
            branch_annotations,
            geological_background,
            geological_data,
            leaves,
            order,
            species,
            tip_annotations,
            tree,
        )
    )
    counts, dataset, dataset_sources, spacing = _prepare_side_panels(
        busco_prefix,
        busco_results,
        busco_summary,
        figure_height,
        figure_width,
        font_size,
        geological_data,
        inputs,
        outputs,
        panel_loader,
        periods,
        row_spacing_points,
        species,
    )
    with (
        matplotlib.rc_context(
            {
                "font.family": font_family,
                "font.size": font_size,
                "pdf.fonttype": 42,
                "svg.fonttype": "none",
            }
        ),
        ExitStack() as cleanup,
    ):
        (
            bar_points,
            credible_label,
            credible_labels,
            figure,
            gap_points,
            label_points,
            left_points,
            text_width,
            tree_points,
            width_points,
        ) = _measure_tree_labels(
            cleanup,
            counts,
            figure_height,
            figure_width,
            font_family,
            font_size,
            nodes,
            species,
            styles,
        )

        (
            base_groups,
            bottom,
            event_groups,
            event_labels,
            height,
            legend_rows,
            period_offset,
            row_heights,
            source_points,
            source_text,
            top,
        ) = _measure_legend_layout(
            STATUS_LABELS,
            count,
            counts,
            credible_labels,
            events,
            figure,
            figure_height,
            figure_width,
            font_size,
            periods,
            show_geological_source,
            spacing,
            text_width,
            width_points,
        )
        axis, label_left, labels, period_texts = _draw_geological_axes(
            bottom,
            count,
            figure,
            font_size,
            gap_points,
            label_points,
            left_points,
            max_age,
            period_offset,
            periods,
            top,
            tree_points,
            width_points,
        )
        (
            branch_lines,
            branches,
            event_marker,
            event_markers,
            interval_count,
            interval_lines,
        ) = _draw_branches_and_intervals(axis, events, font_size, nodes, ys)
        age_texts, tip_texts = _draw_tip_and_age_labels(
            axis, font_size, labels, leaves, node_ages, nodes, selected, styles, ys
        )
        # Choose a nearby label position that clears the branches and earlier
        # age labels. Keep an explicit leader when extra separation is needed.
        event_boxes = _place_age_labels(
            age_texts, axis, branches, event_markers, figure
        )
        bars, busco_axis_label, percent = _draw_count_panel(
            STATUS_COLOURS,
            axis,
            bar_points,
            bottom,
            counts,
            dataset,
            figure,
            font_size,
            gap_points,
            label_left,
            label_points,
            leaves,
            text_width,
            top,
            width_points,
            ys,
        )
        legends = _draw_legends(
            STATUS_COLOURS,
            STATUS_LABELS,
            base_groups,
            counts,
            credible_labels,
            event_groups,
            event_labels,
            event_marker,
            figure,
            font_size,
            height,
            left_points,
            periods,
            row_heights,
            show_geological_source,
            source_points,
            source_text,
            width_points,
        )
        renderer = _validate_presentation_bounds(
            bars, counts, figure, legends, period_texts, tip_texts
        )
        report = {
            "renderer": "genegalleon_dated_tree_presentation",
            "nwkit_version": nwkit.__version__,
            "input_sha256": hashlib.sha256(Path(infile).read_bytes()).hexdigest(),
            "species_order": order,
            "credible_interval_count": interval_count,
            "credible_interval_label": credible_label,
            "credible_interval_style": {
                "alpha": interval_lines[0].get_alpha() if interval_lines else None,
                "zorder": interval_lines[0].get_zorder() if interval_lines else None,
                "tree_zorder": branch_lines[0].get_zorder() if branch_lines else None,
                "bar_capstyle": interval_lines[0].get_solid_capstyle()
                if interval_lines
                else None,
                "bar_linewidth_points": interval_lines[0].get_linewidth()
                if interval_lines
                else None,
                "cap_linewidth_points": interval_lines[1].get_markeredgewidth()
                if interval_lines
                else None,
            },
            "tree_x_axis_position": "bottom",
            "tree_x_axis_y_points": float(axis.bbox.y0) * 72 / figure.dpi,
            "busco_percentage_axis_y_points": float(percent.bbox.y1) * 72 / figure.dpi
            if percent is not None
            else None,
            "busco_plot_bbox_points": [
                float(value) * 72 / figure.dpi for value in bars.bbox.extents
            ]
            if bars is not None
            else None,
            "row_spacing_points": axis.bbox.height * 72 / figure.dpi / (count + 0.2),
            "legend_rows": legend_rows,
            "legend_bbox_points": [
                [
                    float(value) * 72 / figure.dpi
                    for value in legend.get_window_extent(renderer).extents
                ]
                for legend in legends
            ],
            "mean_age_label_count": len(age_texts),
            "geological_background": geological_background,
            "geological_source": GEOLOGICAL_SOURCE if periods else None,
            "geological_source_credit_visible": bool(
                periods and show_geological_source
            ),
            "branch_annotations": events,
            "branch_symbol_bbox_points": [
                [float(value) * 72 / figure.dpi for value in box.extents]
                for box in event_boxes
            ],
            "geological_intervals": periods,
            "geological_label_placement": "above_tree" if periods else None,
            "tree_plot_bbox_points": [
                float(value) * 72 / figure.dpi for value in axis.bbox.extents
            ],
            "geological_labels": [
                {
                    "name": text.get_text(),
                    "anchor_Ma": text.xy[0],
                    "rotation_degrees": 90,
                    "horizontal_offset_points": text.get_position()[0],
                    "bbox_points": [
                        float(value) * 72 / figure.dpi
                        for value in text.get_window_extent(renderer).extents
                    ],
                }
                for text in period_texts
            ],
            "busco_dataset": dataset,
            "busco_dataset_sources": dataset_sources,
            "busco_axis_label": busco_axis_label,
            "busco_counts": counts,
            "tip_y_coordinates": {leaf.name: ys[leaf] for leaf in leaves},
            "all_ages_Ma": [
                {
                    "descendant_species": sorted(node.leaf_names()),
                    "mean": node.props["age"],
                    "low": float(node.props["age_ci_low"])
                    if node.props.get("age_ci_low") is not None
                    else None,
                    "high": float(node.props["age_ci_high"])
                    if node.props.get("age_ci_high") is not None
                    else None,
                }
                for node in nodes
                if not node.is_leaf
            ],
            "figure_size_inches": [figure_width, height],
            "font_size_points": font_size,
            "font_family": font_family,
            "inference": "Saved tree and age intervals reused without inference.",
        }
        paths = [Path(target) for _, target in outputs]
        with output_transaction(paths, create_parents=True) as staged:
            image_metadata: dict[str, str | None] = {
                "Creator": "GeneGalleon / NWKIT " + nwkit.__version__
            }
            if outfile.suffix.lower() == ".pdf":
                image_metadata.update(CreationDate=None, ModDate=None)
            elif outfile.suffix.lower() == ".svg":
                image_metadata["Date"] = None
            else:
                image_metadata = {}
            figure.savefig(
                staged[outfile],
                format=outfile.suffix[1:].lower(),
                dpi=300,
                metadata=image_metadata or None,
            )
            if layout_report:
                Path(staged[Path(layout_report)]).write_text(
                    json.dumps(report, indent=2) + "\n", encoding="utf-8"
                )
        return report
