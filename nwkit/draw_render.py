"""Render ordered artist layers and evaluate a completed tree drawing."""

import hashlib
import math
import sys
from typing import Any

import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, PathPatch

from nwkit import __version__
from nwkit.draw_helpers import (
    DUPLICATION_COLOR,
    LABEL_COLOR,
    LEGEND_EDGE_COLOR,
    SPECIATION_COLOR,
    SUPPORT_LABEL_OFFSET_PT,
    TREE_LINE_CAPSTYLE,
    _add_bottom_guide_title,
    _add_radial_depth_guide,
    _add_scale_bar,
    _add_slanted_depth_guide,
    _add_spiral_depth_key,
    _age_interval_path,
    _branch_marker_artists,
    _compound_polygon_path,
    _densitree_branch_envelope,
    _densitree_topology_groups,
    _densitree_topology_signature,
    _distance_label,
    _draw_probability_pie,
    _draw_tip_images,
    _fit_artists_within_figure,
    _format_property_value,
    _format_support_value,
    _has_meaningful_support,
    _matches_property_filters,
    _node_matches_target,
    _parse_legend_columns,
    _property_color,
    _readable_text_angle,
    _spatial_node_label_placement,
    _spatial_text_alignment,
    _tip_font_style,
    _tip_track_artist,
    _tree_drawing_fingerprint,
)
from nwkit.draw_quality import DrawingArtist, evaluate_drawing
from nwkit.draw_types import (
    AnnotationOptions,
    AxesPlacement,
    DrawingPanels,
    DrawingPlacement,
    DrawingSpacing,
    FontStyle,
    PropertyStyle,
    RenderFrame,
    SpatialOptions,
    TimeDisplay,
    TimeLayout,
    TipLabelGeometry,
    TreeDepths,
)


def _draw_initial_depth_guides(
    *,
    depths: TreeDepths,
    font_style: FontStyle,
    frame: RenderFrame,
    geometry: DrawingPlacement,
    time_layout: TimeLayout,
    tree,
):
    depth_guide_color = "#8c8c8c"
    radial_depth_labels_drawn = None
    if (
        depths.requested_depth_guide is not None
        and time_layout.layout_name == "slanted"
    ):
        for artist in _add_slanted_depth_guide(
            axes=frame.ax,
            root_x=geometry.xcoord[tree],
            ticks=depths.depth_ticks,
            color=depth_guide_color,
            font_family=font_style.font_family,
            font_size=font_style.font_size,
        ):
            frame.drawing_artists.append(
                DrawingArtist(
                    artist,
                    kind="depth_guide",
                    priority=100,
                )
            )
    elif (
        depths.requested_depth_guide is not None and time_layout.layout_name == "radial"
    ):
        radial_artists, radial_depth_labels_drawn = _add_radial_depth_guide(
            axes=frame.ax,
            drawing_layout=geometry.drawing_layout,
            tree=tree,
            ticks=depths.depth_ticks,
            color=depth_guide_color,
            font_family=font_style.font_family,
            font_size=font_style.font_size,
        )
        for artist in radial_artists:
            frame.drawing_artists.append(
                DrawingArtist(
                    artist,
                    kind="depth_guide",
                    priority=100,
                )
            )
    return depth_guide_color, radial_depth_labels_drawn


def _draw_densitree_layers(
    *,
    frame: RenderFrame,
    geometry: DrawingPlacement,
    time_display: TimeDisplay,
    branch_width,
    densitree_ci_color,
    densitree_color,
):
    densitree_topology_groups = _densitree_topology_groups(
        geometry.posterior_layouts,
        time_display.densitree_trees,
    )
    densitree_ci_group_report = []
    if time_display.densitree_mode in {"ci", "both"}:
        envelope_padding = max(min(geometry.x_span, geometry.y_span), 1e-9) * 0.0015
        maximum_group_count = max(
            (len(group) for _, group in densitree_topology_groups),
            default=1,
        )
        for topology_rank, (signature, group) in enumerate(
            densitree_topology_groups,
            start=1,
        ):
            path_groups: dict[Any, list[Any]] = {}
            for _, path_by_clade in group:
                for clade, path in path_by_clade.items():
                    path_groups.setdefault(clade, []).append(path)
            topology_frequency = len(group) / len(geometry.posterior_layouts)
            relative_frequency = len(group) / maximum_group_count
            group_alpha = time_display.densitree_ci_alpha * math.sqrt(
                relative_frequency
            )
            envelope_polygons = []
            retained_counts = []
            for paths in path_groups.values():
                polygons, retained_indices = _densitree_branch_envelope(
                    paths,
                    level=time_display.densitree_ci_level,
                    padding=envelope_padding,
                )
                envelope_polygons.extend(polygons)
                retained_counts.append(len(retained_indices))
            compound_path = _compound_polygon_path(envelope_polygons)
            if compound_path is not None:
                frame.ax.add_patch(
                    PathPatch(
                        compound_path,
                        facecolor=densitree_ci_color,
                        edgecolor="none",
                        alpha=group_alpha,
                        zorder=0.6,
                    )
                )
            densitree_ci_group_report.append(
                {
                    "topology_rank": topology_rank,
                    "sample_count": len(group),
                    "sample_fraction": topology_frequency,
                    "opacity": group_alpha,
                    "branch_count": len(path_groups),
                    "retained_path_count_min": min(retained_counts, default=0),
                    "retained_path_count_max": max(retained_counts, default=0),
                    "topology_sha256": hashlib.sha256(
                        repr(signature).encode("utf-8")
                    ).hexdigest(),
                }
            )

    if time_display.densitree_mode in {"all", "both"}:
        posterior_line_width = max(0.25, float(branch_width) * 0.65)
        posterior_segments = [
            np.asarray(path, dtype=float)
            for _, path_by_clade in geometry.posterior_layouts
            for path in path_by_clade.values()
        ]
        frame.ax.add_collection(
            LineCollection(
                posterior_segments,
                colors=[densitree_color],
                linewidths=posterior_line_width,
                alpha=time_display.densitree_alpha,
                zorder=1,
                capstyle=TREE_LINE_CAPSTYLE,
                joinstyle="round",
            )
        )
    return densitree_ci_group_report


def _draw_tree_branches(
    *,
    frame: RenderFrame,
    geometry: DrawingPlacement,
    branch_color,
    branch_color_property,
    branch_width,
    color_by_node,
    terminal_branch_color,
    tree,
    width_by_node,
):
    for node, path in geometry.drawing_layout.edge_paths.items():
        if branch_color_property not in (None, "") and node.props.get(
            branch_color_property
        ) not in (None, ""):
            edge_color = color_by_node[node]
        elif node.is_leaf and terminal_branch_color:
            edge_color = terminal_branch_color
        else:
            edge_color = color_by_node[node]
        path_x, path_y = zip(*path, strict=True)
        line = frame.ax.plot(
            path_x,
            path_y,
            color=edge_color,
            linewidth=width_by_node[node],
            zorder=2,
            solid_capstyle=TREE_LINE_CAPSTYLE,
            solid_joinstyle="round",
        )[0]
        frame.branch_lines.append((node, line))

    if geometry.drawing_layout.root_path:
        root_x, root_y = zip(*geometry.drawing_layout.root_path, strict=True)
        root_line = frame.ax.plot(
            root_x,
            root_y,
            color=branch_color,
            linewidth=float(branch_width),
            zorder=2,
            solid_capstyle=TREE_LINE_CAPSTYLE,
            solid_joinstyle="round",
        )[0]
        frame.branch_lines.append((tree, root_line))


def _draw_time_intervals(
    *,
    font_style: FontStyle,
    frame: RenderFrame,
    geometry: DrawingPlacement,
    panels: DrawingPanels,
    time_display: TimeDisplay,
    time_layout: TimeLayout,
    branch_width,
    tree,
):
    credible_interval_count = 0
    show_credible_intervals = time_display.credible_interval_mode == "yes" or (
        time_display.credible_interval_mode == "auto"
        and time_layout.has_credible_intervals
        and not time_layout.topology_only_time_layout
    )
    if show_credible_intervals:
        for node in tree.traverse():
            interval_path = _age_interval_path(node, geometry.drawing_layout)
            if interval_path is None:
                if (
                    node.is_root
                    and "age_ci_low" in node.props
                    and "age_ci_high" in node.props
                ):
                    credible_interval_count += 1
                    kind = str(node.props.get("age_ci_kind", "CI"))
                    level = float(node.props.get("age_ci_level", 0.95))
                    label = "{}% {} {:g}–{:g}".format(
                        int(round(level * 100.0)),
                        kind,
                        float(node.props["age_ci_low"]),
                        float(node.props["age_ci_high"]),
                    )
                    artist = frame.ax.annotate(
                        label,
                        xy=(geometry.xcoord[node], geometry.ycoord[node]),
                        xytext=(-4.0, 6.0),
                        textcoords="offset points",
                        ha="right",
                        va="bottom",
                        fontsize=font_style.font_size,
                        fontfamily=font_style.font_family,
                        color="#D55E00",
                        bbox={
                            "boxstyle": "round,pad=0.16",
                            "facecolor": "#ffffff",
                            "edgecolor": "#D55E00",
                            "linewidth": 0.6,
                            "alpha": 0.94,
                        },
                        zorder=10,
                    )
                    frame.drawing_artists.append(
                        DrawingArtist(
                            artist,
                            kind="time_credible_interval",
                            priority=52,
                            movable=True,
                            shift_direction=(0.0, -1.0),
                            owner=node,
                        )
                    )
                continue
            credible_interval_count += 1
            frame.ax.plot(
                interval_path[:, 0],
                interval_path[:, 1],
                color="#D55E00",
                linewidth=max(float(branch_width) * 1.5, 1.0),
                alpha=0.9,
                zorder=4,
                solid_capstyle="round",
            )
            interval_vector = interval_path[1] - interval_path[0]
            display_vector = np.asarray(
                [
                    interval_vector[0]
                    * panels.tree_panel_width_in
                    / max(geometry.x_span, 1e-15),
                    interval_vector[1]
                    * geometry.tree_panel_height_in
                    / max(geometry.y_span, 1e-15),
                ]
            )
            display_length = float(np.linalg.norm(display_vector))
            if display_length > 1e-15:
                display_normal = (
                    np.asarray(
                        [-display_vector[1], display_vector[0]],
                        dtype=float,
                    )
                    / display_length
                )
                cap_half_inches = max(float(font_style.font_size) * 0.16, 1.1) / 72.0
                cap_vector = np.asarray(
                    [
                        display_normal[0]
                        * cap_half_inches
                        * geometry.x_span
                        / max(panels.tree_panel_width_in, 1e-15),
                        display_normal[1]
                        * cap_half_inches
                        * geometry.y_span
                        / max(geometry.tree_panel_height_in, 1e-15),
                    ]
                )
                for endpoint in interval_path:
                    cap = np.vstack(
                        (
                            endpoint - cap_vector,
                            endpoint + cap_vector,
                        )
                    )
                    frame.ax.plot(
                        cap[:, 0],
                        cap[:, 1],
                        color="#D55E00",
                        linewidth=0.75,
                        zorder=4.1,
                    )
    return credible_interval_count


def _draw_time_constraints(
    *,
    font_style: FontStyle,
    frame: RenderFrame,
    geometry: DrawingPlacement,
    spatial_options: SpatialOptions,
    time_display: TimeDisplay,
    tree,
):
    calibration_nodes = [
        node
        for node in tree.traverse()
        if node.props.get("calibration_type") not in (None, "")
    ]
    show_constraints = time_display.constraint_mode == "yes" or (
        time_display.constraint_mode == "auto" and bool(calibration_nodes)
    )
    constraint_color_by_type = {
        "point": "#CC79A7",
        "bounded": "#009E73",
        "lower": "#0072B2",
        "upper": "#E69F00",
    }
    constraint_artist_count = 0
    if show_constraints:
        for node in calibration_nodes:
            calibration_type = str(node.props["calibration_type"])
            lower = node.props.get("calibration_lower")
            upper = node.props.get("calibration_upper")
            if calibration_type == "point":
                label = "@{:g}".format(float(lower))
            elif lower is not None and upper is not None:
                label = "{:g}–{:g}".format(float(lower), float(upper))
            elif lower is not None:
                label = "≥{:g}".format(float(lower))
            elif upper is not None:
                label = "≤{:g}".format(float(upper))
            else:
                label = str(node.props.get("calibration_raw", "calibration"))
            color = constraint_color_by_type.get(calibration_type, "#666666")
            if spatial_options.spatial_layout:
                placement = _spatial_node_label_placement(
                    node,
                    geometry.drawing_layout,
                    clearance_points=max(float(font_style.font_size) * 1.35, 7.0),
                )
            else:
                placement = {
                    "offset": (0.0, max(float(font_style.font_size) * 1.25, 8.0)),
                    "horizontal_alignment": "center",
                    "vertical_alignment": "bottom",
                    "rotation": 0.0,
                    "shift_direction": (0.0, 1.0),
                }
            frame.ax.scatter(
                [geometry.xcoord[node]],
                [geometry.ycoord[node]],
                s=max(float(font_style.font_size) * 0.42, 2.8) ** 2,
                marker="o",
                facecolor="#ffffff",
                edgecolor=color,
                linewidth=0.75,
                zorder=9.5,
            )
            artist = frame.ax.annotate(
                label,
                xy=(geometry.xcoord[node], geometry.ycoord[node]),
                xytext=placement["offset"],
                textcoords="offset points",
                ha=placement["horizontal_alignment"],
                va=placement["vertical_alignment"],
                rotation=placement["rotation"],
                rotation_mode="anchor",
                fontsize=font_style.font_size,
                fontfamily=font_style.font_family,
                color=color,
                bbox={
                    "boxstyle": "round,pad=0.18",
                    "facecolor": "#ffffff",
                    "edgecolor": color,
                    "linewidth": 0.6,
                    "alpha": 0.94,
                },
                arrowprops={
                    "arrowstyle": "-",
                    "color": color,
                    "linewidth": 0.65,
                    "shrinkA": 1.5,
                    "shrinkB": max(float(font_style.font_size) * 0.22, 1.5),
                },
                zorder=10,
            )
            frame.drawing_artists.append(
                DrawingArtist(
                    artist,
                    kind="time_constraint",
                    priority=50,
                    movable=True,
                    shift_direction=placement["shift_direction"],
                    owner=node,
                )
            )
            constraint_artist_count += 1
    return constraint_artist_count


def _draw_node_type_markers(
    *,
    frame: RenderFrame,
    geometry: DrawingPlacement,
    branch_markers,
    group_color_by_name,
    node_type_by_node,
    tree,
):
    marker_area = 18.0 * (0.5**2)
    marker_size_pt = 4.8 * 0.5
    marker_drawing_artists, marker_legend_handles = _branch_marker_artists(
        frame.ax,
        frame.fig,
        tree,
        geometry.drawing_layout,
        branch_markers,
    )
    frame.drawing_artists.extend(marker_drawing_artists)
    legend_handles = list(marker_legend_handles)
    if len(node_type_by_node) > 0:
        for node, node_type in node_type_by_node.items():
            marker_color = (
                DUPLICATION_COLOR if (node_type == "duplication") else SPECIATION_COLOR
            )
            frame.ax.scatter(
                [geometry.xcoord[node]],
                [geometry.ycoord[node]],
                s=marker_area,
                marker="o",
                facecolor=marker_color,
                edgecolor=LEGEND_EDGE_COLOR,
                linewidth=0.4,
                zorder=5,
            )
        legend_handles.extend(
            [
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    linestyle="None",
                    markerfacecolor=SPECIATION_COLOR,
                    markeredgecolor=LEGEND_EDGE_COLOR,
                    markeredgewidth=0.4,
                    markersize=marker_size_pt,
                    label="Speciation node",
                ),
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    linestyle="None",
                    markerfacecolor=DUPLICATION_COLOR,
                    markeredgecolor=LEGEND_EDGE_COLOR,
                    markeredgewidth=0.4,
                    markersize=marker_size_pt,
                    label="Duplication node",
                ),
            ]
        )
    if len(group_color_by_name) > 0:
        legend_handles.extend(
            Patch(facecolor=color, edgecolor="none", label=group_name)
            for group_name, color in sorted(group_color_by_name.items())
        )
    return legend_handles


def _property_legend_handles(
    *,
    annotations: AnnotationOptions,
    panels: DrawingPanels,
    property_style: PropertyStyle,
    branch_color,
    branch_color_property,
    branch_width,
    branch_width_property,
    legend_handles,
    tip_badge_property,
    tree,
    width_by_node,
):
    if annotations.node_pie_properties:
        legend_handles.extend(
            Patch(
                facecolor=_property_color(
                    prop=prop,
                    value=prop,
                    property_colors=property_style.property_colors,
                    fallback_index=index,
                    palette=property_style.trait_palette,
                ),
                edgecolor="none",
                label=prop,
            )
            for index, prop in enumerate(annotations.node_pie_properties)
        )
    if panels.badge_values:
        legend_handles.extend(
            Patch(
                facecolor=_property_color(
                    prop=tip_badge_property,
                    value=value,
                    property_colors=property_style.property_colors,
                    fallback_index=panels.badge_value_index[value],
                    palette=property_style.trait_palette,
                ),
                edgecolor="none",
                label="{}: {}".format(tip_badge_property, value),
            )
            for value in panels.badge_values
        )
    legend_handles.extend(
        Patch(facecolor=color, edgecolor="none", label=label)
        for label, color in panels.tip_track_legend_entries
    )
    if branch_color_property not in (None, ""):
        branch_values = sorted(
            {
                str(node.props[branch_color_property])
                for node in tree.traverse()
                if not node.is_root
                and node.props.get(branch_color_property) not in (None, "")
            }
        )
        branch_value_index = {value: index for index, value in enumerate(branch_values)}
        legend_handles.extend(
            Line2D(
                [0],
                [0],
                color=_property_color(
                    prop=branch_color_property,
                    value=value,
                    property_colors=property_style.property_colors,
                    fallback_index=branch_value_index[value],
                    palette=property_style.trait_palette,
                ),
                linewidth=float(branch_width),
                label="{}: {}".format(branch_color_property, value),
            )
            for value in branch_values
        )
        if any(
            node.props.get(branch_color_property) in (None, "")
            for node in tree.traverse()
            if not node.is_root
        ):
            legend_handles.append(
                Line2D(
                    [0],
                    [0],
                    color=branch_color,
                    linewidth=float(branch_width),
                    label="{}: missing".format(branch_color_property),
                )
            )
    if branch_width_property not in (None, ""):
        numeric_width_nodes = []
        for node in tree.traverse():
            if node.is_root:
                continue
            try:
                value = float(node.props.get(branch_width_property))
            except (TypeError, ValueError):
                continue
            if math.isfinite(value):
                numeric_width_nodes.append((value, node))
        if numeric_width_nodes:
            numeric_width_nodes.sort(key=lambda item: item[0])
            width_examples = [numeric_width_nodes[0]]
            if numeric_width_nodes[-1][0] != numeric_width_nodes[0][0]:
                width_examples.append(numeric_width_nodes[-1])
            legend_handles.extend(
                Line2D(
                    [0],
                    [0],
                    color=branch_color,
                    linewidth=width_by_node[node],
                    label="{}: {:g}".format(branch_width_property, value),
                )
                for value, node in width_examples
            )
    unique_legend_handles = []
    seen_legend_labels = set()
    for handle in legend_handles:
        label = handle.get_label()
        if label in seen_legend_labels:
            continue
        seen_legend_labels.add(label)
        unique_legend_handles.append(handle)
    legend_handles = unique_legend_handles
    return legend_handles


def _draw_tree_legend(
    *,
    font_style: FontStyle,
    frame: RenderFrame,
    legend,
    legend_columns,
    legend_handles,
    legend_position,
):
    if legend and len(legend_handles) > 0:
        resolved_legend_position = str(legend_position).strip().lower()
        if resolved_legend_position == "auto":
            resolved_legend_position = "right" if len(legend_handles) > 12 else "top"
        if resolved_legend_position not in {"top", "right"}:
            raise ValueError("'--legend-position' must be auto, top, or right.")
        if resolved_legend_position == "right":
            legend_location = "upper left"
            legend_anchor = (1.005, 1.0)
            legend_ncol = 1
        else:
            legend_location = "lower right"
            legend_anchor = (1.0, 1.005)
            legend_ncol = _parse_legend_columns(
                legend_columns,
                len(legend_handles),
            )
        legend_artist = frame.ax.legend(
            handles=legend_handles,
            loc=legend_location,
            bbox_to_anchor=legend_anchor,
            frameon=False,
            borderaxespad=0.1,
            handletextpad=0.3,
            labelspacing=0.2,
            ncol=legend_ncol,
            prop={"family": font_style.font_family, "size": font_style.font_size},
        )
        frame.drawing_artists.append(
            DrawingArtist(
                legend_artist,
                kind="legend",
                priority=100,
            )
        )


def _draw_support_labels(
    *,
    font_style: FontStyle,
    frame: RenderFrame,
    geometry: DrawingPlacement,
    spatial_options: SpatialOptions,
    support_labels,
    support_min,
    tree,
):
    for node in tree.traverse():
        if not support_labels:
            continue
        if not _has_meaningful_support(node):
            continue
        support_value = float(node.support)
        if not math.isfinite(support_value):
            raise ValueError("Displayed support values must be finite numbers.")
        if (support_min is not None) and (support_value < support_min):
            continue
        label_x, label_y = geometry.drawing_layout.support_anchors[node]
        edge_angle = geometry.drawing_layout.support_angles.get(node, 0.0)
        if spatial_options.spatial_layout:
            radians = math.radians(edge_angle)
            support_offset = (
                -math.sin(radians) * 1.5,
                math.cos(radians) * 1.5,
            )
        else:
            support_offset = (0.0, SUPPORT_LABEL_OFFSET_PT)
        support_artist = frame.ax.annotate(
            _format_support_value(node.support),
            xy=(label_x, label_y),
            xytext=support_offset,
            textcoords="offset points",
            va="bottom",
            ha="center",
            rotation=_readable_text_angle(edge_angle)
            if spatial_options.spatial_layout
            else 0.0,
            rotation_mode="anchor",
            fontsize=font_style.font_size,
            fontfamily=font_style.font_family,
            color=LABEL_COLOR,
        )
        frame.drawing_artists.append(
            DrawingArtist(
                support_artist,
                kind="support_label",
                priority=55,
                movable=True,
                shift_direction=(
                    -math.sin(math.radians(edge_angle)),
                    math.cos(math.radians(edge_angle)),
                )
                if spatial_options.spatial_layout
                else (0.0, 1.0),
                owner=node,
            )
        )


def _draw_tip_labels(
    *,
    annotations: AnnotationOptions,
    font_style: FontStyle,
    frame: RenderFrame,
    geometry: DrawingPlacement,
    panels: DrawingPanels,
    property_style: PropertyStyle,
    spatial_options: SpatialOptions,
    branch_color,
    leaf_label_color_by_leaf,
    tip_badge_missing_label,
    tip_badge_property,
    tip_label_font_style,
    tip_label_position,
    tip_labels,
    tip_track_properties,
) -> TipLabelGeometry:
    data_per_inch = geometry.x_span / max(panels.tree_panel_width_in, 0.2)
    leaf_pies_drawn = any(
        annotations.node_pie_properties
        and _node_matches_target(leaf, annotations.node_pie_target)
        and (
            (not annotations.node_pie_leaf_filters)
            or _matches_property_filters(
                leaf,
                annotations.node_pie_leaf_filters,
                option_name="--node-pie-leaf-filter",
            )
        )
        and all(prop in leaf.props for prop in annotations.node_pie_properties)
        for leaf in geometry.leaf_order
    )
    tip_pie_clearance = 0.0
    if leaf_pies_drawn:
        tip_pie_radius_pt = max(float(font_style.font_size), 4.5) / 2.0
        tip_pie_clearance = ((tip_pie_radius_pt + 2.0) / 72.0) * data_per_inch
    label_offset = max(geometry.x_span * 0.02, 0.06, tip_pie_clearance)
    tip_label_artists = []
    track_span_data = (spatial_options.track_span_pt / 72.0) * data_per_inch
    for leaf in geometry.leaf_order:
        label_angle = geometry.drawing_layout.label_angles.get(leaf, 0.0)
        if spatial_options.spatial_layout:
            text_rotation, text_alignment = _spatial_text_alignment(label_angle)
            radians = math.radians(label_angle)
            for track_index, property_name in enumerate(tip_track_properties):
                track_distance = (
                    3.0
                    + (spatial_options.track_size_pt / 2.0)
                    + track_index * spatial_options.track_stride_pt
                )
                track_artist = _tip_track_artist(
                    ax=frame.ax,
                    xy=(geometry.xcoord[leaf], geometry.ycoord[leaf]),
                    offset_points=(
                        math.cos(radians) * track_distance,
                        math.sin(radians) * track_distance,
                    ),
                    color=panels.tip_track_color_by_leaf_property[
                        (leaf, property_name)
                    ],
                    size_points=spatial_options.track_size_pt,
                )
                frame.drawing_artists.append(
                    DrawingArtist(
                        track_artist,
                        kind="tip_track",
                        priority=85,
                        owner=leaf,
                    )
                )
            label_offset_points = (
                math.cos(radians) * (4.0 + spatial_options.track_span_pt),
                math.sin(radians) * (4.0 + spatial_options.track_span_pt),
            )
            label_x = geometry.xcoord[leaf]
            label_y = geometry.ycoord[leaf]
            if tip_labels:
                tip_artist = frame.ax.annotate(
                    panels.tip_label_text_by_leaf[leaf],
                    xy=(label_x, label_y),
                    xytext=label_offset_points,
                    textcoords="offset points",
                    va="center",
                    ha=text_alignment,
                    rotation=text_rotation,
                    rotation_mode="anchor",
                    fontsize=font_style.font_size,
                    fontfamily=font_style.font_family,
                    fontstyle=_tip_font_style(leaf.name or "", tip_label_font_style),
                    linespacing=1.15,
                    color=leaf_label_color_by_leaf.get(leaf, LABEL_COLOR),
                    annotation_clip=False,
                )
                tip_label_artists.append(tip_artist)
                frame.drawing_artists.append(
                    DrawingArtist(
                        tip_artist,
                        kind="tip_label",
                        priority=95,
                        owner=leaf,
                    )
                )
        elif spatial_options.resolved_tip_label_position == "branch-end":
            track_origin_x = geometry.xcoord[leaf] + label_offset
            for track_index, property_name in enumerate(tip_track_properties):
                track_artist = _tip_track_artist(
                    ax=frame.ax,
                    xy=(track_origin_x, geometry.ycoord[leaf]),
                    offset_points=(
                        (spatial_options.track_size_pt / 2.0)
                        + track_index * spatial_options.track_stride_pt,
                        0.0,
                    ),
                    color=panels.tip_track_color_by_leaf_property[
                        (leaf, property_name)
                    ],
                    size_points=spatial_options.track_size_pt,
                )
                frame.drawing_artists.append(
                    DrawingArtist(
                        track_artist,
                        kind="tip_track",
                        priority=85,
                        owner=leaf,
                    )
                )
            label_x = track_origin_x + track_span_data + ((1.0 / 72.0) * data_per_inch)
            label_y = geometry.ycoord[leaf]
            if tip_labels:
                tip_artist = frame.ax.text(
                    label_x,
                    label_y,
                    panels.tip_label_text_by_leaf[leaf],
                    va="center",
                    ha="left",
                    fontsize=font_style.font_size,
                    fontfamily=font_style.font_family,
                    fontstyle=_tip_font_style(leaf.name or "", tip_label_font_style),
                    linespacing=1.15,
                    color=leaf_label_color_by_leaf.get(leaf, LABEL_COLOR),
                )
                tip_label_artists.append(tip_artist)
                frame.drawing_artists.append(
                    DrawingArtist(
                        tip_artist,
                        kind="tip_label",
                        priority=95,
                        owner=leaf,
                    )
                )
        elif spatial_options.resolved_tip_label_position == "aligned":
            track_origin_x = geometry.x_max + label_offset
            for track_index, property_name in enumerate(tip_track_properties):
                track_artist = _tip_track_artist(
                    ax=frame.ax,
                    xy=(track_origin_x, geometry.ycoord[leaf]),
                    offset_points=(
                        (spatial_options.track_size_pt / 2.0)
                        + track_index * spatial_options.track_stride_pt,
                        0.0,
                    ),
                    color=panels.tip_track_color_by_leaf_property[
                        (leaf, property_name)
                    ],
                    size_points=spatial_options.track_size_pt,
                )
                frame.drawing_artists.append(
                    DrawingArtist(
                        track_artist,
                        kind="tip_track",
                        priority=85,
                        owner=leaf,
                    )
                )
            label_x = track_origin_x + track_span_data + ((1.0 / 72.0) * data_per_inch)
            label_y = geometry.ycoord[leaf]
            if tip_labels:
                tip_artist = frame.ax.text(
                    label_x,
                    label_y,
                    panels.tip_label_text_by_leaf[leaf],
                    va="center",
                    ha="left",
                    fontsize=font_style.font_size,
                    fontfamily=font_style.font_family,
                    fontstyle=_tip_font_style(leaf.name or "", tip_label_font_style),
                    linespacing=1.15,
                    color=leaf_label_color_by_leaf.get(leaf, LABEL_COLOR),
                )
                tip_label_artists.append(tip_artist)
                frame.drawing_artists.append(
                    DrawingArtist(
                        tip_artist,
                        kind="tip_label",
                        priority=95,
                        owner=leaf,
                    )
                )
        else:
            raise ValueError(
                "Unsupported '--tip-label-position': {}".format(tip_label_position)
            )
        if str(leaf.props.get("nwkit_collapsed", "")).lower() == "true":
            collapsed_marker = (
                (3, 0, label_angle - 90.0) if spatial_options.spatial_layout else ">"
            )
            frame.ax.scatter(
                [geometry.xcoord[leaf]],
                [geometry.ycoord[leaf]],
                s=max(float(font_style.font_size), 5.0) ** 2,
                marker=collapsed_marker,
                facecolor=branch_color,
                edgecolor=LEGEND_EDGE_COLOR,
                linewidth=0.4,
                zorder=7,
            )
        if tip_badge_property not in (None, ""):
            raw_value = leaf.props.get(tip_badge_property)
            is_missing = raw_value in (None, "")
            if is_missing and tip_badge_missing_label in (None, ""):
                continue
            badge_value = str(tip_badge_missing_label if is_missing else raw_value)
            badge_color = _property_color(
                prop=tip_badge_property,
                value=badge_value,
                property_colors=property_style.property_colors,
                fallback_index=panels.badge_value_index[badge_value],
                palette=property_style.trait_palette,
            )
            badge_style = {
                "boxstyle": "round,pad=0.16",
                "facecolor": badge_color,
                "edgecolor": LEGEND_EDGE_COLOR,
                "linewidth": 0.4,
                "alpha": 0.92,
            }
            if spatial_options.spatial_layout:
                text_rotation, _ = _spatial_text_alignment(label_angle)
                radians = math.radians(label_angle)
                badge_offset_points = (
                    spatial_options.track_span_pt
                    + (
                        panels.tip_label_size_by_leaf[leaf][0] * 72.0
                        if tip_labels
                        else 0.0
                    )
                    + 4.0
                    + max(len(badge_value), 1) * font_style.font_size * 0.29
                    + (font_style.font_size * 0.16)
                )
                badge_artist = frame.ax.annotate(
                    badge_value,
                    xy=(label_x, label_y),
                    xytext=(
                        math.cos(radians) * badge_offset_points,
                        math.sin(radians) * badge_offset_points,
                    ),
                    textcoords="offset points",
                    va="center",
                    ha="center",
                    rotation=text_rotation,
                    rotation_mode="anchor",
                    fontsize=font_style.font_size,
                    fontfamily=font_style.font_family,
                    fontweight="bold",
                    color=LABEL_COLOR,
                    bbox=badge_style,
                    zorder=8,
                    annotation_clip=False,
                )
                frame.drawing_artists.append(
                    DrawingArtist(
                        badge_artist,
                        kind="tip_badge",
                        priority=70,
                        movable=True,
                        shift_direction=(-math.sin(radians), math.cos(radians)),
                        owner=leaf,
                    )
                )
            else:
                approximate_label_width = (
                    panels.tip_label_size_by_leaf[leaf][0] * data_per_inch
                    if tip_labels
                    else 0.0
                )
                approximate_badge_width = (
                    max(len(badge_value), 1)
                    * (font_style.font_size * 0.58 / 72.0)
                    * data_per_inch
                )
                badge_padding = (font_style.font_size * 0.16 / 72.0) * data_per_inch
                badge_gap = (4.0 / 72.0) * data_per_inch
                badge_offset = (
                    approximate_label_width
                    + badge_gap
                    + (approximate_badge_width / 2.0)
                    + badge_padding
                )
                badge_artist = frame.ax.text(
                    label_x + badge_offset,
                    label_y,
                    badge_value,
                    va="center",
                    ha="center",
                    fontsize=font_style.font_size,
                    fontfamily=font_style.font_family,
                    fontweight="bold",
                    color=LABEL_COLOR,
                    bbox=badge_style,
                    zorder=8,
                )
                frame.drawing_artists.append(
                    DrawingArtist(
                        badge_artist,
                        kind="tip_badge",
                        priority=70,
                        owner=leaf,
                    )
                )
    return TipLabelGeometry(
        data_per_inch=data_per_inch,
        label_offset=label_offset,
        tip_label_artists=tip_label_artists,
    )


def _draw_node_annotations(
    *,
    annotations: AnnotationOptions,
    font_style: FontStyle,
    frame: RenderFrame,
    geometry: DrawingPlacement,
    property_style: PropertyStyle,
    spatial_options: SpatialOptions,
    tree,
):
    for node in tree.traverse():
        leaf_passes_pie_filters = (
            (not node.is_leaf)
            or (not annotations.node_pie_leaf_filters)
            or _matches_property_filters(
                node,
                annotations.node_pie_leaf_filters,
                option_name="--node-pie-leaf-filter",
            )
        )
        if (
            annotations.node_pie_properties
            and _node_matches_target(node, annotations.node_pie_target)
            and leaf_passes_pie_filters
        ):
            if all(prop in node.props for prop in annotations.node_pie_properties):
                probabilities = [
                    node.props[prop] for prop in annotations.node_pie_properties
                ]
                colors = [
                    _property_color(
                        prop=prop,
                        value=prop,
                        property_colors=property_style.property_colors,
                        fallback_index=index,
                        palette=property_style.trait_palette,
                    )
                    for index, prop in enumerate(annotations.node_pie_properties)
                ]
                node_pie_artist = _draw_probability_pie(
                    frame.ax,
                    geometry.xcoord[node],
                    geometry.ycoord[node],
                    probabilities=probabilities,
                    colors=colors,
                    marker_size_pt=max(font_style.font_size, 4.5),
                )
                if node_pie_artist is not None:
                    frame.drawing_artists.append(
                        DrawingArtist(
                            node_pie_artist,
                            kind="node_pie",
                            priority=80,
                            owner=node,
                        )
                    )
        if (
            annotations.node_label_property not in (None, "")
            and annotations.node_label_property in node.props
            and _node_matches_target(node, annotations.node_label_target)
            and _matches_property_filters(node, annotations.node_label_filters)
        ):
            node_has_pie = (
                bool(annotations.node_pie_properties)
                and _node_matches_target(node, annotations.node_pie_target)
                and leaf_passes_pie_filters
                and all(prop in node.props for prop in annotations.node_pie_properties)
            )
            pie_radius = (
                max(float(font_style.font_size), 4.5) / 2.0 if node_has_pie else 0.0
            )
            if spatial_options.spatial_layout:
                placement = _spatial_node_label_placement(
                    node,
                    geometry.drawing_layout,
                    clearance_points=(
                        pie_radius + max(float(font_style.font_size) * 0.7, 4.0)
                    ),
                )
            else:
                clearance = pie_radius + 2.0
                placement = {
                    "offset": (
                        (-clearance, clearance)
                        if node.is_leaf
                        else (clearance, clearance)
                    ),
                    "rotation": 0.0,
                    "horizontal_alignment": "right" if node.is_leaf else "left",
                    "vertical_alignment": "bottom",
                    "shift_direction": (0.0, 1.0),
                }
            node_label_artist = frame.ax.annotate(
                "{}{}".format(
                    annotations.node_label_prefix,
                    _format_property_value(
                        node.props[annotations.node_label_property],
                        decimals=annotations.node_label_decimals,
                    ),
                ),
                xy=(geometry.xcoord[node], geometry.ycoord[node]),
                xytext=placement["offset"],
                textcoords="offset points",
                va=placement["vertical_alignment"],
                ha=placement["horizontal_alignment"],
                rotation=placement["rotation"],
                rotation_mode="anchor",
                fontsize=font_style.font_size,
                fontfamily=font_style.font_family,
                color=LABEL_COLOR,
                zorder=9,
            )
            frame.drawing_artists.append(
                DrawingArtist(
                    node_label_artist,
                    kind="node_label",
                    priority=45,
                    movable=True,
                    shift_direction=placement["shift_direction"],
                    owner=node,
                )
            )


def _draw_root_marker(
    *,
    font_style: FontStyle,
    frame: RenderFrame,
    geometry: DrawingPlacement,
    root_marker,
    root_marker_color,
    root_marker_size,
    tree,
):
    root_marker_size_pt = max(float(font_style.font_size), 5.0)
    if root_marker_size is not None:
        root_marker_size_pt = float(root_marker_size)
        if root_marker_size_pt <= 0.0:
            raise ValueError("--root-marker-size must be greater than zero.")
    if root_marker != "none":
        marker_by_name = {"circle": "o", "diamond": "D"}
        if root_marker not in marker_by_name:
            raise ValueError("Unsupported '--root-marker': {}".format(root_marker))
        frame.ax.scatter(
            [geometry.xcoord[tree]],
            [geometry.ycoord[tree]],
            s=root_marker_size_pt**2,
            marker=marker_by_name[root_marker],
            facecolor=root_marker_color,
            edgecolor=LEGEND_EDGE_COLOR,
            linewidth=0.45,
            zorder=8,
        )
    return root_marker_size_pt


def _configure_drawing_axes(
    *,
    annotations: AnnotationOptions,
    font_style: FontStyle,
    frame: RenderFrame,
    geometry: DrawingPlacement,
    label_geometry: TipLabelGeometry,
    panels: DrawingPanels,
    spacing: DrawingSpacing,
    spatial_options: SpatialOptions,
    root_marker,
    root_marker_size_pt,
    tip_badge_property,
    tip_labels,
    tip_track_properties,
    tree,
) -> AxesPlacement:
    label_panel_span = geometry.x_span * (
        panels.label_panel_width_in / max(panels.tree_panel_width_in, 0.2)
    )
    root_has_pie = (
        bool(annotations.node_pie_properties)
        and _node_matches_target(tree, annotations.node_pie_target)
        and all(prop in tree.props for prop in annotations.node_pie_properties)
    )
    root_has_symbol = (root_marker != "none") or root_has_pie
    root_symbol_diameter_pt = 0.0
    if root_marker != "none":
        root_symbol_diameter_pt = root_marker_size_pt
    if root_has_pie:
        root_symbol_diameter_pt = max(
            root_symbol_diameter_pt, max(float(font_style.font_size), 4.5)
        )
    root_symbol_margin = 0.0
    if root_has_symbol:
        root_symbol_radius_pt = root_symbol_diameter_pt / 2.0
        root_symbol_margin = (
            (root_symbol_radius_pt + 1.5) / 72.0
        ) * label_geometry.data_per_inch
    if spatial_options.spatial_layout:
        horizontal_label_extent_in = 0.0
        vertical_label_extent_in = 0.0
        if tip_labels:
            for leaf in geometry.leaf_order:
                width_in, height_in = panels.tip_label_size_by_leaf[leaf]
                angle = math.radians(
                    geometry.drawing_layout.label_angles.get(leaf, 0.0)
                )
                horizontal_label_extent_in = max(
                    horizontal_label_extent_in,
                    abs(math.cos(angle)) * width_in + abs(math.sin(angle)) * height_in,
                )
                vertical_label_extent_in = max(
                    vertical_label_extent_in,
                    abs(math.sin(angle)) * width_in + abs(math.cos(angle)) * height_in,
                )
            horizontal_label_extent_in += 5.0 / 72.0
            vertical_label_extent_in += 5.0 / 72.0
        if tip_badge_property not in (None, ""):
            badge_extent = max(0.3, float(font_style.font_size) * 2.0 / 72.0)
            horizontal_label_extent_in += badge_extent
            vertical_label_extent_in += badge_extent
        if tip_track_properties:
            horizontal_label_extent_in += spatial_options.track_span_pt / 72.0
            vertical_label_extent_in += spatial_options.track_span_pt / 72.0
        horizontal_label_margin = horizontal_label_extent_in * (
            geometry.x_span / max(panels.tree_panel_width_in, 0.2)
        )
        vertical_label_margin = vertical_label_extent_in * (
            geometry.y_span / max(geometry.tree_panel_height_in, 0.2)
        )
        plot_pad_x = max(geometry.x_span * 0.035, horizontal_label_margin)
        plot_pad_y = max(geometry.y_span * 0.035, vertical_label_margin)
        frame.ax.set_xlim(geometry.x_min - plot_pad_x, geometry.x_max + plot_pad_x)
        frame.ax.set_ylim(geometry.y_min - plot_pad_y, geometry.y_max + plot_pad_y)
        if geometry.drawing_layout.equal_aspect:
            frame.ax.set_aspect("equal", adjustable="box")
    else:
        root_path_span = 0.0
        if geometry.drawing_layout.root_path:
            root_path_span = abs(
                geometry.drawing_layout.root_path[-1][0]
                - geometry.drawing_layout.root_path[0][0]
            )
        left_plot_margin = max(root_path_span * 1.2, root_symbol_margin)
        frame.ax.set_xlim(
            geometry.x_min - left_plot_margin,
            geometry.x_max + label_geometry.label_offset + label_panel_span,
        )
        if len(geometry.leaf_order) == 0:
            frame.ax.set_ylim(0.5, -0.5)
        else:
            frame.ax.set_ylim(geometry.y_max + 0.5, geometry.y_min - 0.5)
    frame.ax.axis("off")
    axes_left = spacing.left_margin_in / panels.fig_width
    axes_right = 1.0 - (
        (spacing.right_margin_in + panels.image_panel_width_in) / panels.fig_width
    )
    axes_top = 1.0 - (spacing.top_margin_in / geometry.fig_height)
    axes_bottom = spacing.bottom_margin_in / geometry.fig_height
    frame.fig.subplots_adjust(
        left=axes_left,
        right=axes_right,
        top=axes_top,
        bottom=axes_bottom,
    )
    return AxesPlacement(
        axes_left=axes_left,
        axes_right=axes_right,
        axes_top=axes_top,
        axes_bottom=axes_bottom,
    )


def _draw_distance_guides(
    *,
    axes_placement: AxesPlacement,
    depths: TreeDepths,
    font_style: FontStyle,
    frame: RenderFrame,
    geometry: DrawingPlacement,
    spacing: DrawingSpacing,
    time_layout: TimeLayout,
    branch_color,
    branch_length_unit,
    depth_guide_color,
):
    if depths.requested_scale_bar is not None:
        scale_label = "{:g}".format(depths.requested_scale_bar)
        if str(branch_length_unit).strip():
            scale_label = "{} {}".format(
                scale_label,
                str(branch_length_unit).strip(),
            )
        scale_artist = _add_scale_bar(
            figure=frame.fig,
            axes=frame.ax,
            size=depths.requested_scale_bar,
            label=scale_label,
            color=branch_color,
            font_family=font_style.font_family,
            font_size=font_style.font_size,
            anchor_x=axes_placement.axes_left,
            anchor_y=(2.0 / 72.0) / geometry.fig_height,
        )
        frame.drawing_artists.append(
            DrawingArtist(
                scale_artist,
                kind="scale_bar",
                priority=100,
            )
        )
    if depths.requested_depth_guide is not None:
        radial_guide_name = (
            "Concentric arcs"
            if (
                time_layout.layout_name == "radial"
                and geometry.drawing_layout.metadata.get(
                    "angular_span_degrees",
                    360.0,
                )
                < 360.0
            )
            else "Concentric rings"
        )
        guide_prefix = {
            "slanted": "Root-to-node distance",
            "radial": ("{} every {:g}: root-to-node distance; root = 0").format(
                radial_guide_name, depths.requested_depth_guide
            ),
            "spiral": "Root-to-node distance encoded across spiral band",
        }[time_layout.layout_name]
        title_artist = _add_bottom_guide_title(
            figure=frame.fig,
            text=_distance_label(guide_prefix, branch_length_unit),
            x=axes_placement.axes_left,
            y=(2.0 / 72.0) / geometry.fig_height,
            font_family=font_style.font_family,
            font_size=font_style.font_size,
        )
        frame.drawing_artists.append(
            DrawingArtist(
                title_artist,
                kind="depth_guide_title",
                priority=100,
            )
        )
        if time_layout.layout_name == "spiral":
            for artist in _add_spiral_depth_key(
                figure=frame.fig,
                ticks=depths.depth_ticks,
                tree_span=depths.branch_depth_span,
                axes_left=axes_placement.axes_left,
                axes_right=axes_placement.axes_right,
                strip_height_points=spacing.depth_guide_strip_pt,
                font_family=font_style.font_family,
                font_size=font_style.font_size,
                color=depth_guide_color,
            ):
                frame.drawing_artists.append(
                    DrawingArtist(
                        artist,
                        kind="depth_guide",
                        priority=100,
                    )
                )


def _draw_tip_image_panel(
    *,
    axes_placement: AxesPlacement,
    frame: RenderFrame,
    geometry: DrawingPlacement,
    panels: DrawingPanels,
    tip_image_by_leaf,
):
    if tip_image_by_leaf:
        image_ax = frame.fig.add_axes(
            [
                axes_placement.axes_right,
                axes_placement.axes_bottom,
                panels.image_panel_width_in / panels.fig_width,
                axes_placement.axes_top - axes_placement.axes_bottom,
            ],
            sharey=frame.ax,
        )
        image_ax.set_xlim(0.0, 1.0)
        image_ax.set_ylim(frame.ax.get_ylim())
        image_ax.set_facecolor("none")
        image_ax.axis("off")
        _draw_tip_images(
            ax=image_ax,
            leaf_order=geometry.leaf_order,
            ycoord=geometry.ycoord,
            image_by_leaf=tip_image_by_leaf,
            image_size_pt=panels.tip_image_size_pt,
        )


def _finalize_drawing_quality(
    *,
    frame: RenderFrame,
    collision_policy,
):
    fit_report = _fit_artists_within_figure(
        figure=frame.fig,
        axes=frame.ax,
        artists=[item.artist for item in frame.drawing_artists],
    )
    requested_collision_policy = str(collision_policy).strip().lower()
    if requested_collision_policy not in {"resolve", "warn", "error", "ignore"}:
        raise ValueError(
            "'--collision-policy' must be resolve, warn, error, or ignore."
        )
    quality_report = evaluate_drawing(
        figure=frame.fig,
        artists=frame.drawing_artists,
        branch_lines=frame.branch_lines,
        policy=("resolve" if requested_collision_policy == "resolve" else "ignore"),
        emit_warning=False,
    )
    final_fit_report = _fit_artists_within_figure(
        figure=frame.fig,
        axes=frame.ax,
        artists=[item.artist for item in frame.drawing_artists],
    )
    if not final_fit_report["fits_within_figure"]:
        overflow_message = (
            "Drawing annotations exceed the fixed figure boundary by up to "
            "{:.2f} point(s). Increase the figure size or enable label wrapping."
        ).format(final_fit_report["maximum_overflow_points"])
        if str(collision_policy).strip().lower() == "error":
            raise ValueError(overflow_message)
        if str(collision_policy).strip().lower() in {"resolve", "warn"}:
            sys.stderr.write(overflow_message + "\n")
    final_quality = evaluate_drawing(
        figure=frame.fig,
        artists=frame.drawing_artists,
        branch_lines=frame.branch_lines,
        policy="ignore",
        max_iterations=0,
    )
    final_collision_count = final_quality["final_collision_count"]
    if final_collision_count:
        collision_message = (
            "Drawing contains {} unresolved collision(s) after final fitting. "
            "Increase --figure-width (and --figure-height when fixed), use "
            "--tip-label-wrap auto with --tip-spacing label-aware, or "
            "reduce annotation density."
        ).format(final_collision_count)
        if requested_collision_policy == "error":
            raise ValueError(collision_message)
        if requested_collision_policy in {"resolve", "warn"}:
            sys.stderr.write(collision_message + "\n")
    branch_crossing_count = final_quality["branch_crossing_count"]
    if branch_crossing_count:
        crossing_message = (
            "Drawing contains {} crossing branch pair(s); the selected "
            "geometry is not planar for this input."
        ).format(branch_crossing_count)
        if requested_collision_policy == "error":
            raise ValueError(crossing_message)
        if requested_collision_policy in {"resolve", "warn"}:
            sys.stderr.write(crossing_message + "\n")
    quality_report["collision_policy"] = requested_collision_policy
    quality_report["final_collision_count"] = final_quality["final_collision_count"]
    quality_report["final_collisions_by_kind"] = final_quality[
        "final_collisions_by_kind"
    ]
    quality_report["annotation_occupied_fraction"] = final_quality[
        "annotation_occupied_fraction"
    ]
    quality_report["branch_collision_check_complete"] = final_quality[
        "branch_collision_check_complete"
    ]
    quality_report["branch_crossing_count"] = final_quality["branch_crossing_count"]
    quality_report["branch_crossing_check_complete"] = final_quality[
        "branch_crossing_check_complete"
    ]
    quality_report.update(final_fit_report)
    return quality_report, fit_report


def _add_drawing_report_metadata(
    *,
    depths: TreeDepths,
    font_style: FontStyle,
    geometry: DrawingPlacement,
    panels: DrawingPanels,
    spatial_options: SpatialOptions,
    time_display: TimeDisplay,
    time_layout: TimeLayout,
    branch_length_unit,
    collapsed_clades,
    constraint_artist_count,
    credible_interval_count,
    densitree_ci_group_report,
    fit_report,
    image_format,
    mcmctree_posterior,
    quality_report,
    radial_depth_labels_drawn,
    tip_track_palette,
    tip_track_properties,
    tree,
    unrooted_method,
):
    quality_report.update(
        {
            "nwkit_version": __version__,
            "output_format": str(image_format),
            "layout_requested": time_layout.layout_name,
            "layout": geometry.drawing_layout.name,
            "subtree_packing": spatial_options.resolved_subtree_packing,
            "unrooted_method": str(unrooted_method),
            "daylight_iterations_requested": spatial_options.daylight_iterations,
            **geometry.drawing_layout.metadata,
            "figure_width_inches": panels.fig_width,
            "figure_height_inches": geometry.fig_height,
            "font_family": str(font_style.font_family),
            "font_size_points": float(font_style.font_size),
            "input_tip_count": sum(
                int(item.get("tip_count", 1)) for item in collapsed_clades
            )
            + len(geometry.leaf_order)
            - len(collapsed_clades),
            "visible_tip_count": len(geometry.leaf_order),
            "display_tree_sha256": _tree_drawing_fingerprint(tree),
            "collapsed_clades": collapsed_clades,
            "wrapped_tip_count": sum(
                "\n" in panels.tip_label_text_by_leaf[leaf]
                for leaf in geometry.leaf_order
            ),
            "branch_lengths_encoded": (
                time_layout.layout_name in depths.branch_length_layouts
                and not depths.use_topology_depth
            ),
            "branch_length_encoding": (
                "segment"
                if time_layout.layout_name in depths.segment_length_layouts
                and not depths.use_topology_depth
                else (
                    "depth-projection"
                    if time_layout.layout_name in depths.depth_projection_layouts
                    and not depths.use_topology_depth
                    else (
                        "warped-depth"
                        if time_layout.layout_name in depths.warped_depth_layouts
                        and not depths.use_topology_depth
                        else "none"
                    )
                )
            ),
            "scale_bar": depths.requested_scale_bar,
            "scale_bar_position": (
                "bottom-reserved" if depths.requested_scale_bar is not None else None
            ),
            "scale_bar_label_position": (
                "above" if depths.requested_scale_bar is not None else None
            ),
            "depth_guide_interval": depths.requested_depth_guide,
            "depth_guide_type": (
                {
                    "slanted": "axis-grid",
                    "radial": (
                        "concentric-arcs"
                        if geometry.drawing_layout.metadata.get(
                            "angular_span_degrees",
                            360.0,
                        )
                        < 360.0
                        else "concentric-rings"
                    ),
                    "spiral": "spiral-depth-key",
                }[time_layout.layout_name]
                if depths.requested_depth_guide is not None
                else None
            ),
            "depth_guide_in_panel_labels": (
                radial_depth_labels_drawn
                if depths.requested_depth_guide is not None
                and time_layout.layout_name == "radial"
                else None
            ),
            "branch_length_unit": str(branch_length_unit),
            "tip_track_properties": list(tip_track_properties),
            "tip_track_palette": str(tip_track_palette),
            "tip_spacing": spatial_options.resolved_tip_spacing,
            "time_constraints": time_display.constraint_mode,
            "time_constraint_count": constraint_artist_count,
            "time_credible_intervals": time_display.credible_interval_mode,
            "time_credible_interval_count": credible_interval_count,
            "densitree": time_display.densitree_mode,
            "densitree_source": (
                "tree-collection"
                if time_display.densitree_trees
                else ("mcmctree-ages" if mcmctree_posterior is not None else None)
            ),
            "densitree_sample_count": len(geometry.posterior_layouts),
            "densitree_topology_count": (
                len(
                    {
                        _densitree_topology_signature(sample)
                        for sample in time_display.densitree_trees
                    }
                )
                if time_display.densitree_trees
                else (1 if geometry.posterior_layouts else 0)
            ),
            "densitree_root_alignment": (
                "reference-coordinate"
                if time_display.densitree_trees
                else ("shared-topology" if geometry.posterior_layouts else None)
            ),
            "densitree_ci_level": (
                time_display.densitree_ci_level
                if time_display.densitree_mode in {"ci", "both"}
                else None
            ),
            "densitree_ci_interpretation": (
                (
                    "topology-stratified branchwise empirical central-path "
                    "envelope; opacity encodes relative topology frequency"
                    if time_display.densitree_trees
                    else "branchwise empirical central-path envelope"
                )
                if time_display.densitree_mode in {"ci", "both"}
                else None
            ),
            "densitree_ci_path_coverage": (
                time_display.densitree_ci_level
                if time_display.densitree_mode in {"ci", "both"}
                else None
            ),
            "densitree_ci_topology_group_count": (
                len(densitree_ci_group_report)
                if time_display.densitree_mode in {"ci", "both"}
                else 0
            ),
            "densitree_ci_topology_frequency_encoding": (
                "square-root opacity relative to the most frequent topology"
                if time_display.densitree_trees
                and time_display.densitree_mode in {"ci", "both"}
                else None
            ),
            "densitree_ci_topology_groups": densitree_ci_group_report,
            "initial_fit": fit_report,
        }
    )
