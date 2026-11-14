"""Validate drawing options and measure and place tree layouts."""

import math
import sys

from nwkit.draw_helpers import (
    TIP_LABEL_GAP_PT,
    _angular_sector_bounds,
    _depth_guide_ticks,
    _get_tree_plot_coordinates,
    _has_positive_branch_length,
    _legend_is_needed,
    _make_posterior_drawing_layouts,
    _make_tree_collection_drawing_layouts,
    _matches_property_filters,
    _node_matches_target,
    _parse_depth_guide,
    _parse_scale_bar,
    _parse_tip_label_wrap,
    _polar_auto_figure_height,
    _prepare_tip_label_text,
    _tip_font_style,
    _tip_track_colors,
    _validated_tip_image_size,
)
from nwkit.draw_layouts import make_tree_layout
from nwkit.draw_types import (
    AnnotationOptions,
    DrawingPanels,
    DrawingPlacement,
    DrawingSpacing,
    FontStyle,
    PropertyStyle,
    SpatialOptions,
    TimeDisplay,
    TimeLayout,
    TreeDepths,
)


def _resolve_time_display(
    *,
    densitree,
    densitree_alpha,
    densitree_ci_alpha,
    densitree_ci_level,
    densitree_trees,
    mcmctree_posterior,
    time_constraints,
    time_credible_intervals,
) -> TimeDisplay:
    densitree_mode = str(densitree).strip().lower()
    if densitree_mode not in {"none", "all", "ci", "both"}:
        raise ValueError("'--densitree' must be none, all, ci, or both.")
    densitree_trees = list(densitree_trees or [])
    if densitree_mode != "none" and mcmctree_posterior is None and not densitree_trees:
        raise ValueError(
            "'--densitree' requires '--mcmctree-posterior' or '--densitree-trees'."
        )
    if mcmctree_posterior is not None and densitree_trees:
        raise ValueError(
            "'--mcmctree-posterior' and '--densitree-trees' are mutually exclusive."
        )
    densitree_alpha = float(densitree_alpha)
    densitree_ci_alpha = float(densitree_ci_alpha)
    densitree_ci_level = float(densitree_ci_level)
    if not 0.0 < densitree_alpha <= 1.0:
        raise ValueError(
            "'--densitree-alpha' must be greater than 0 and no greater than 1."
        )
    if not 0.0 < densitree_ci_alpha <= 1.0:
        raise ValueError(
            "'--densitree-ci-alpha' must be greater than 0 and no greater than 1."
        )
    if not 0.0 < densitree_ci_level < 1.0:
        raise ValueError("'--densitree-ci-level' must be between 0 and 1.")
    constraint_mode = str(time_constraints).strip().lower()
    credible_interval_mode = str(time_credible_intervals).strip().lower()
    if constraint_mode not in {"auto", "yes", "no"}:
        raise ValueError("'--time-constraints' must be auto, yes, or no.")
    if credible_interval_mode not in {"auto", "yes", "no"}:
        raise ValueError("'--time-credible-intervals' must be auto, yes, or no.")
    return TimeDisplay(
        densitree_mode=densitree_mode,
        densitree_trees=densitree_trees,
        densitree_alpha=densitree_alpha,
        densitree_ci_alpha=densitree_ci_alpha,
        densitree_ci_level=densitree_ci_level,
        constraint_mode=constraint_mode,
        credible_interval_mode=credible_interval_mode,
    )


def _validate_node_annotation_options(
    *,
    annotations: AnnotationOptions,
    tree,
):
    if annotations.node_pie_properties:
        _node_matches_target(tree, annotations.node_pie_target)
    if annotations.node_label_property not in (None, ""):
        _node_matches_target(tree, annotations.node_label_target)
        if int(annotations.node_label_decimals) < 0:
            raise ValueError("--node-label-decimals must be zero or greater.")
        # Parse every filter before rendering so malformed later filters are
        # not masked when an earlier property happens to be absent.
        _matches_property_filters(
            tree,
            annotations.node_label_filters,
            option_name="--node-label-filter",
        )
    if annotations.node_pie_leaf_filters:
        _matches_property_filters(
            tree,
            annotations.node_pie_leaf_filters,
            option_name="--node-pie-leaf-filter",
        )


def _validate_time_layout(
    *,
    time_display: TimeDisplay,
    layout,
    tree,
) -> TimeLayout:
    layout_name = str(layout).strip().lower()
    supported_layouts = {
        "rectangular",
        "slanted",
        "cladogram",
        "circular",
        "radial",
        "unrooted",
        "spiral",
        "fractal",
    }
    if layout_name not in supported_layouts:
        raise ValueError("Unsupported '--layout': {}".format(layout))
    if (
        time_display.densitree_mode != "none"
        and layout_name in {"cladogram", "fractal"}
        and not time_display.densitree_trees
    ):
        raise ValueError(
            "'--densitree' requires a layout that encodes branch-length or "
            "root-to-node time variation; cladogram and fractal do not."
        )
    if time_display.densitree_trees and layout_name in {"unrooted", "fractal"}:
        raise ValueError(
            "'--densitree-trees' currently supports rectangular, slanted, "
            "cladogram, circular, radial, and spiral layouts."
        )
    topology_only_time_layout = layout_name in {"cladogram", "fractal"}
    has_credible_intervals = any("age_ci_low" in node.props for node in tree.traverse())
    if (
        topology_only_time_layout
        and time_display.credible_interval_mode == "yes"
        and has_credible_intervals
    ):
        raise ValueError(
            "'--time-credible-intervals yes' requires a layout that encodes "
            "branch-length or root-to-node time variation."
        )
    if (
        topology_only_time_layout
        and time_display.credible_interval_mode == "auto"
        and has_credible_intervals
    ):
        sys.stderr.write(
            "Skipping time credible intervals because the selected layout "
            "encodes topology rather than time.\n"
        )
    return TimeLayout(
        layout_name=layout_name,
        topology_only_time_layout=topology_only_time_layout,
        has_credible_intervals=has_credible_intervals,
    )


def _resolve_spatial_options(
    *,
    time_layout: TimeLayout,
    angular_center,
    angular_span,
    daylight_iterations,
    label_panel_width,
    subtree_packing,
    tip_image_by_leaf,
    tip_label_font_style,
    tip_label_position,
    tip_spacing,
    tip_track_properties,
    tip_track_size,
) -> SpatialOptions:
    angular_span = float(angular_span)
    angular_center = float(angular_center)
    if not math.isfinite(angular_span) or angular_span <= 0.0 or angular_span > 360.0:
        raise ValueError(
            "--angular-span must be greater than zero and no greater than 360."
        )
    if not math.isfinite(angular_center):
        raise ValueError("--angular-center must be a finite number.")
    if time_layout.layout_name not in {"circular", "radial"}:
        if angular_span != 360.0:
            raise ValueError(
                "'--angular-span' is supported only with circular and radial layouts."
            )
        if angular_center % 360.0 != 90.0:
            raise ValueError(
                "'--angular-center' is supported only with circular and radial layouts."
            )
    resolved_subtree_packing = str(subtree_packing).strip().lower()
    if resolved_subtree_packing not in {"standard", "tidy"}:
        raise ValueError("'--subtree-packing' must be standard or tidy.")
    if resolved_subtree_packing == "tidy" and time_layout.layout_name not in {
        "rectangular",
        "circular",
        "spiral",
    }:
        raise ValueError(
            "'--subtree-packing tidy' is supported only with rectangular, "
            "circular, and spiral layouts."
        )
    spatial_layout = time_layout.layout_name in {
        "circular",
        "radial",
        "unrooted",
        "spiral",
        "fractal",
    }
    if spatial_layout and tip_image_by_leaf:
        raise ValueError(
            "'--tip-image-manifest' is currently supported only with "
            "rectangular, slanted, and cladogram layouts."
        )
    if spatial_layout and label_panel_width is not None:
        raise ValueError(
            "'--label-panel-width' is not used by two-dimensional layouts."
        )
    track_size_pt = float(tip_track_size)
    if not math.isfinite(track_size_pt) or track_size_pt <= 0.0:
        raise ValueError("'--tip-track-size' must be greater than zero.")
    track_stride_pt = track_size_pt + 2.0
    track_span_pt = len(tip_track_properties) * track_stride_pt
    _tip_font_style("", tip_label_font_style)
    resolved_tip_spacing = str(tip_spacing).strip().lower()
    if resolved_tip_spacing not in {"uniform", "label-aware"}:
        raise ValueError("'--tip-spacing' must be uniform or label-aware.")
    daylight_iterations = int(daylight_iterations)
    if daylight_iterations <= 0:
        raise ValueError("'--daylight-iterations' must be greater than zero.")
    resolved_tip_label_position = str(tip_label_position).strip().lower()
    if resolved_tip_label_position == "auto":
        resolved_tip_label_position = (
            "aligned"
            if time_layout.layout_name == "rectangular"
            and resolved_subtree_packing == "standard"
            else "branch-end"
        )
    if spatial_layout and resolved_tip_label_position != "branch-end":
        raise ValueError(
            "Two-dimensional layouts require '--tip-label-position branch-end'."
        )
    return SpatialOptions(
        angular_span=angular_span,
        angular_center=angular_center,
        resolved_subtree_packing=resolved_subtree_packing,
        spatial_layout=spatial_layout,
        track_size_pt=track_size_pt,
        track_stride_pt=track_stride_pt,
        track_span_pt=track_span_pt,
        resolved_tip_spacing=resolved_tip_spacing,
        daylight_iterations=daylight_iterations,
        resolved_tip_label_position=resolved_tip_label_position,
    )


def _prepare_tree_depths(
    *,
    time_layout: TimeLayout,
    depth_guide,
    scale_bar,
    tree,
) -> TreeDepths:
    use_topology_depth = not _has_positive_branch_length(tree)
    if use_topology_depth:
        sys.stderr.write(
            "Tree has no positive branch lengths; drawing positions by topology depth.\n"
        )
    base_xcoord, _, base_leaf_order = _get_tree_plot_coordinates(
        tree=tree,
        use_topology_depth=use_topology_depth,
    )
    base_x_values = list(base_xcoord.values())
    base_x_min = min(base_x_values) if base_x_values else 0.0
    base_x_max = max(base_x_values) if base_x_values else 1.0
    base_x_span = max(base_x_max - base_x_min, 1.0)
    segment_length_layouts = {
        "rectangular",
        "circular",
        "unrooted",
    }
    depth_projection_layouts = {"slanted", "radial"}
    warped_depth_layouts = {"spiral"}
    depth_guide_layouts = depth_projection_layouts | warped_depth_layouts
    branch_length_layouts = (
        segment_length_layouts | depth_projection_layouts | warped_depth_layouts
    )
    branch_depth_span = max(base_x_max - base_x_min, 0.0)
    requested_scale_bar = _parse_scale_bar(scale_bar, branch_depth_span)
    requested_depth_guide = _parse_depth_guide(
        depth_guide,
        branch_depth_span,
    )
    if (
        requested_scale_bar is not None
        and time_layout.layout_name not in segment_length_layouts
    ):
        raise ValueError(
            "'--scale-bar' requires a branch-length-preserving layout with "
            "directly measurable branch segments."
        )
    if requested_scale_bar is not None and use_topology_depth:
        raise ValueError("'--scale-bar' requires positive input branch lengths.")
    if (
        requested_scale_bar is not None
        and requested_scale_bar > branch_depth_span + 1e-12
    ):
        raise ValueError(
            "'--scale-bar' must not exceed the displayed tree-depth span "
            "({:g}).".format(branch_depth_span)
        )
    if (
        requested_depth_guide is not None
        and time_layout.layout_name not in depth_guide_layouts
    ):
        raise ValueError(
            "'--depth-guide' is supported by slanted, radial, and spiral layouts."
        )
    if requested_depth_guide is not None and use_topology_depth:
        raise ValueError("'--depth-guide' requires positive input branch lengths.")
    depth_ticks = _depth_guide_ticks(
        branch_depth_span,
        requested_depth_guide,
    )
    return TreeDepths(
        use_topology_depth=use_topology_depth,
        base_xcoord=base_xcoord,
        base_leaf_order=base_leaf_order,
        base_x_min=base_x_min,
        base_x_max=base_x_max,
        base_x_span=base_x_span,
        segment_length_layouts=segment_length_layouts,
        depth_projection_layouts=depth_projection_layouts,
        warped_depth_layouts=warped_depth_layouts,
        branch_length_layouts=branch_length_layouts,
        branch_depth_span=branch_depth_span,
        requested_scale_bar=requested_scale_bar,
        requested_depth_guide=requested_depth_guide,
        depth_ticks=depth_ticks,
    )


def _prepare_drawing_panels(
    *,
    depths: TreeDepths,
    font_style: FontStyle,
    property_style: PropertyStyle,
    spatial_options: SpatialOptions,
    time_layout: TimeLayout,
    figure_height,
    figure_width,
    label_panel_width,
    tip_badge_missing_label,
    tip_badge_property,
    tip_image_by_leaf,
    tip_image_gap,
    tip_image_size,
    tip_label_wrap,
    tip_labels,
    tip_track_palette,
    tip_track_properties,
    tip_track_type,
    tree,
) -> DrawingPanels:
    leaf_order = depths.base_leaf_order
    fig_width = float(figure_width)
    if fig_width <= 0.0:
        raise ValueError("--figure-width must be greater than zero.")
    tip_image_size_pt = _validated_tip_image_size(tip_image_size)
    tip_image_gap_pt = float(tip_image_gap)
    if tip_image_gap_pt < 0.0:
        raise ValueError("--tip-image-gap must be zero or greater.")
    image_panel_width_in = 0.0
    if tip_image_by_leaf:
        image_panel_width_in = (tip_image_size_pt + (2.0 * tip_image_gap_pt)) / 72.0
    main_panel_width_in = fig_width - image_panel_width_in
    if main_panel_width_in <= 0.2:
        raise ValueError(
            "--figure-width is too small for the requested tip-image column."
        )
    label_height_guess_in = (
        float(figure_height)
        if figure_height is not None
        else (fig_width * 0.72 if time_layout.layout_name == "fractal" else fig_width)
    )
    if tip_labels:
        tip_label_text_by_leaf, tip_label_size_by_leaf = _prepare_tip_label_text(
            leaf_order=leaf_order,
            wrap=tip_label_wrap,
            font_size=font_style.font_size,
            font_family=font_style.font_family,
            layout_name=time_layout.layout_name,
            panel_width_in=main_panel_width_in,
            panel_height_in=max(label_height_guess_in, 0.2),
        )
    else:
        _parse_tip_label_wrap(tip_label_wrap)
        tip_label_text_by_leaf = {leaf: str(leaf.name or "") for leaf in leaf_order}
        tip_label_size_by_leaf = {leaf: (0.0, 0.0) for leaf in leaf_order}
    max_tip_label_width_in = max(
        (size[0] for size in tip_label_size_by_leaf.values()),
        default=0.0,
    )
    max_tip_label_height_in = max(
        (size[1] for size in tip_label_size_by_leaf.values()),
        default=(float(font_style.font_size) / 72.0 if tip_labels else 0.0),
    )
    tip_track_color_by_leaf_property, tip_track_legend_entries = _tip_track_colors(
        tree=tree,
        properties=tip_track_properties,
        mode=tip_track_type,
        property_colors=property_style.property_colors,
        categorical_palette=property_style.trait_palette,
        continuous_palette=tip_track_palette,
    )
    badge_values = []
    if tip_badge_property not in (None, ""):
        badge_values = sorted(
            {
                str(
                    tip_badge_missing_label
                    if leaf.props.get(tip_badge_property) in (None, "")
                    else leaf.props.get(tip_badge_property)
                )
                for leaf in tree.leaves()
                if (
                    leaf.props.get(tip_badge_property) not in (None, "")
                    or tip_badge_missing_label not in (None, "")
                )
            }
        )
    badge_value_index = {value: index for index, value in enumerate(badge_values)}
    if spatial_options.spatial_layout:
        label_panel_width_in = 0.0
    elif (
        not tip_labels
        and tip_badge_property in (None, "")
        and not tip_track_properties
        and label_panel_width is None
    ):
        label_panel_width_in = 0.0
    elif label_panel_width is None:
        label_panel_width_in = min(
            max(0.8, max_tip_label_width_in + (5.0 / 72.0)),
            main_panel_width_in * 0.58,
        )
    else:
        label_panel_width_in = float(label_panel_width)
        if label_panel_width_in <= 0.0 or label_panel_width_in >= main_panel_width_in:
            raise ValueError(
                "--label-panel-width must be greater than zero and leave room "
                "for the tree and tip-image column within --figure-width."
            )
    tree_panel_width_in = main_panel_width_in - label_panel_width_in
    if tree_panel_width_in < 0.2:
        raise ValueError(
            "--figure-width is too small for the requested label and tip-image panels."
        )
    return DrawingPanels(
        leaf_order=leaf_order,
        fig_width=fig_width,
        tip_image_size_pt=tip_image_size_pt,
        image_panel_width_in=image_panel_width_in,
        main_panel_width_in=main_panel_width_in,
        tip_label_text_by_leaf=tip_label_text_by_leaf,
        tip_label_size_by_leaf=tip_label_size_by_leaf,
        max_tip_label_width_in=max_tip_label_width_in,
        max_tip_label_height_in=max_tip_label_height_in,
        tip_track_color_by_leaf_property=tip_track_color_by_leaf_property,
        tip_track_legend_entries=tip_track_legend_entries,
        badge_values=badge_values,
        badge_value_index=badge_value_index,
        label_panel_width_in=label_panel_width_in,
        tree_panel_width_in=tree_panel_width_in,
    )


def _prepare_drawing_spacing(
    *,
    annotations: AnnotationOptions,
    depths: TreeDepths,
    font_style: FontStyle,
    panels: DrawingPanels,
    spatial_options: SpatialOptions,
    time_layout: TimeLayout,
    branch_color_property,
    branch_markers,
    branch_width_property,
    figure_height,
    group_color_by_name,
    legend,
    node_type_by_node,
    tip_badge_property,
    tip_image_by_leaf,
    tip_labels,
    tip_track_properties,
) -> DrawingSpacing:
    spacing_height_by_leaf = {}
    for leaf in panels.leaf_order:
        spacing_height = panels.tip_label_size_by_leaf[leaf][1] if tip_labels else 0.0
        if tip_track_properties:
            spacing_height = max(spacing_height, spatial_options.track_size_pt / 72.0)
        if tip_badge_property not in (None, ""):
            spacing_height = max(
                spacing_height,
                float(font_style.font_size) * 1.3 / 72.0,
            )
        if str(leaf.name) in tip_image_by_leaf:
            spacing_height = max(spacing_height, panels.tip_image_size_pt / 72.0)
        spacing_height_by_leaf[leaf] = max(
            spacing_height,
            float(font_style.font_size) / 72.0,
        )
    if spatial_options.resolved_tip_spacing == "label-aware" and spacing_height_by_leaf:
        row_pitch_pt = (
            sum(spacing_height_by_leaf.values()) / len(spacing_height_by_leaf)
        ) * 72.0 + TIP_LABEL_GAP_PT
    else:
        row_pitch_pt = (
            max(
                float(font_style.font_size),
                panels.max_tip_label_height_in * 72.0,
            )
            + TIP_LABEL_GAP_PT
        )
        if tip_image_by_leaf:
            row_pitch_pt = max(
                row_pitch_pt,
                panels.tip_image_size_pt + TIP_LABEL_GAP_PT,
            )
    row_pitch_in = row_pitch_pt / 72.0
    property_legend_needed = bool(
        panels.badge_values or annotations.node_pie_properties
    )
    legend_needed = _legend_is_needed(
        legend,
        node_type_by_node,
        group_color_by_name,
        property_legend_needed,
        panels.tip_track_legend_entries,
        branch_color_property not in (None, ""),
        branch_width_property not in (None, ""),
        branch_markers,
    )
    top_margin_in = (44.0 / 72.0) if legend_needed else (4.0 / 72.0)
    scale_bar_strip_pt = (
        max(float(font_style.font_size) * 1.4 + 8.0, 18.0)
        if depths.requested_scale_bar is not None
        else 0.0
    )
    depth_guide_strip_pt = (
        max(float(font_style.font_size) * 2.6 + 14.0, 34.0)
        if depths.requested_depth_guide is not None
        else 0.0
    )
    bottom_margin_in = (3.0 + max(scale_bar_strip_pt, depth_guide_strip_pt)) / 72.0
    left_margin_in = 2.0 / 72.0
    right_margin_in = 2.0 / 72.0
    if figure_height is None and spatial_options.spatial_layout:
        if time_layout.layout_name in {"circular", "radial"}:
            radial_label_extent = (
                panels.max_tip_label_width_in if tip_labels else 0.0
            ) + (spatial_options.track_span_pt / 72.0)
            if tip_badge_property not in (None, ""):
                radial_label_extent += (float(font_style.font_size) * 1.3) / 72.0
            fig_height = _polar_auto_figure_height(
                panel_width=panels.tree_panel_width_in,
                angular_span=spatial_options.angular_span,
                angular_center=spatial_options.angular_center,
                radial_label_extent=radial_label_extent,
                tangential_label_extent=panels.max_tip_label_height_in,
                top_margin=top_margin_in,
                bottom_margin=bottom_margin_in,
            )
        else:
            fig_height = (
                panels.fig_width * 0.72
                if time_layout.layout_name == "fractal"
                else panels.fig_width
            )
    elif figure_height is None:
        # Refined after tidy coordinates have been calculated.
        fig_height = (
            max(len(panels.leaf_order), 1) * row_pitch_in
            + top_margin_in
            + bottom_margin_in
        )
    else:
        fig_height = float(figure_height)
        if fig_height <= top_margin_in + bottom_margin_in:
            raise ValueError("--figure-height is too small for the required margins.")
    return DrawingSpacing(
        spacing_height_by_leaf=spacing_height_by_leaf,
        row_pitch_pt=row_pitch_pt,
        row_pitch_in=row_pitch_in,
        legend_needed=legend_needed,
        top_margin_in=top_margin_in,
        depth_guide_strip_pt=depth_guide_strip_pt,
        bottom_margin_in=bottom_margin_in,
        left_margin_in=left_margin_in,
        right_margin_in=right_margin_in,
        fig_height=fig_height,
    )


def _place_drawing(
    *,
    depths: TreeDepths,
    panels: DrawingPanels,
    spacing: DrawingSpacing,
    spatial_options: SpatialOptions,
    time_display: TimeDisplay,
    time_layout: TimeLayout,
    figure_height,
    mcmctree_posterior,
    mcmctree_posterior_topology,
    posterior_ladderize,
    spiral_turns,
    tip_image_by_leaf,
    tip_labels,
    tree,
    unrooted_method,
) -> DrawingPlacement:
    fig_height = spacing.fig_height
    leaf_order = panels.leaf_order
    available_height_in = max(
        fig_height - spacing.top_margin_in - spacing.bottom_margin_in, 0.2
    )
    terminal_extent_by_leaf = {}
    label_data_per_inch = depths.base_x_span / max(panels.tree_panel_width_in, 0.2)
    for leaf in leaf_order:
        label_width_in = panels.tip_label_size_by_leaf.get(leaf, (0.0, 0.0))[0]
        label_extent = (
            label_width_in + (spatial_options.track_span_pt / 72.0)
        ) * label_data_per_inch
        if tip_labels and spatial_options.resolved_tip_label_position == "aligned":
            label_extent += max(depths.base_x_max - depths.base_xcoord[leaf], 0.0)
        terminal_extent_by_leaf[leaf] = label_extent

    layout_aspect_ratio = panels.tree_panel_width_in / available_height_in
    if (
        time_layout.layout_name == "fractal"
        and spatial_options.resolved_tip_spacing == "label-aware"
        and tip_labels
    ):
        normal_label_margin = panels.max_tip_label_width_in + (5.0 / 72.0)
        label_aware_width = max(
            panels.tree_panel_width_in - (2.0 * normal_label_margin),
            0.2,
        )
        label_aware_height = max(
            available_height_in - (2.0 * normal_label_margin),
            0.2,
        )
        layout_aspect_ratio = label_aware_width / label_aware_height
    spacing_size_by_leaf = {
        leaf: (
            panels.tip_label_size_by_leaf[leaf][0],
            spacing.spacing_height_by_leaf[leaf],
        )
        for leaf in leaf_order
    }
    drawing_layout = make_tree_layout(
        tree,
        layout=time_layout.layout_name,
        use_topology_depth=depths.use_topology_depth,
        aspect_ratio=layout_aspect_ratio,
        spiral_turns=spiral_turns,
        angular_span=spatial_options.angular_span,
        angular_center=spatial_options.angular_center,
        terminal_extent_by_leaf=terminal_extent_by_leaf,
        label_size_by_leaf=spacing_size_by_leaf,
        tip_spacing=spatial_options.resolved_tip_spacing,
        subtree_packing=spatial_options.resolved_subtree_packing,
        unrooted_method=unrooted_method,
        daylight_iterations=spatial_options.daylight_iterations,
    )
    if time_display.densitree_mode == "none":
        posterior_layouts = []
    elif time_display.densitree_trees:
        posterior_layouts = _make_tree_collection_drawing_layouts(
            tree,
            time_display.densitree_trees,
            layout=time_layout.layout_name,
            use_topology_depth=depths.use_topology_depth,
            aspect_ratio=layout_aspect_ratio,
            spiral_turns=spiral_turns,
            angular_span=spatial_options.angular_span,
            angular_center=spatial_options.angular_center,
            terminal_extent_by_leaf=terminal_extent_by_leaf,
            label_size_by_leaf=spacing_size_by_leaf,
            tip_spacing=spatial_options.resolved_tip_spacing,
            subtree_packing=spatial_options.resolved_subtree_packing,
            unrooted_method=unrooted_method,
            daylight_iterations=spatial_options.daylight_iterations,
        )
    else:
        posterior_layouts = _make_posterior_drawing_layouts(
            tree,
            mcmctree_posterior,
            mcmctree_posterior_topology,
            bool(posterior_ladderize),
            layout=time_layout.layout_name,
            use_topology_depth=depths.use_topology_depth,
            aspect_ratio=layout_aspect_ratio,
            spiral_turns=spiral_turns,
            angular_span=spatial_options.angular_span,
            angular_center=spatial_options.angular_center,
            terminal_extent_by_leaf=terminal_extent_by_leaf,
            label_size_by_leaf=spacing_size_by_leaf,
            tip_spacing=spatial_options.resolved_tip_spacing,
            subtree_packing=spatial_options.resolved_subtree_packing,
            unrooted_method=unrooted_method,
            daylight_iterations=spatial_options.daylight_iterations,
        )
    xcoord = drawing_layout.xcoord
    ycoord = drawing_layout.ycoord
    leaf_order = drawing_layout.leaf_order
    x_min, x_max, y_min, y_max = drawing_layout.bounds
    if posterior_layouts:
        posterior_bounds = [bounds for bounds, _ in posterior_layouts]
        x_min = min([x_min] + [bounds[0] for bounds in posterior_bounds])
        x_max = max([x_max] + [bounds[1] for bounds in posterior_bounds])
        y_min = min([y_min] + [bounds[2] for bounds in posterior_bounds])
        y_max = max([y_max] + [bounds[3] for bounds in posterior_bounds])
    if (
        depths.requested_depth_guide is not None
        and time_layout.layout_name == "radial"
        and depths.depth_ticks
    ):
        guide_radius = max(depths.depth_ticks)
        root_x = drawing_layout.xcoord[tree]
        root_y = drawing_layout.ycoord[tree]
        guide_x_min, guide_x_max, guide_y_min, guide_y_max = _angular_sector_bounds(
            guide_radius,
            drawing_layout.metadata.get("angular_span_degrees", 360.0),
            drawing_layout.metadata.get("angular_center_degrees", 90.0),
        )
        x_min = min(x_min, root_x + guide_x_min)
        x_max = max(x_max, root_x + guide_x_max)
        y_min = min(y_min, root_y + guide_y_min)
        y_max = max(y_max, root_y + guide_y_max)
    x_span = max(x_max - x_min, 1e-9)
    y_span = max(y_max - y_min, 1e-9)

    if not spatial_options.spatial_layout:
        tree_panel_height_in = max(y_span + 1.0, 1.0) * spacing.row_pitch_in
        if figure_height is None:
            fig_height = (
                tree_panel_height_in + spacing.top_margin_in + spacing.bottom_margin_in
            )
        if (
            tip_image_by_leaf
            and fig_height + 10**-9
            < tree_panel_height_in + spacing.top_margin_in + spacing.bottom_margin_in
        ):
            raise ValueError(
                "--figure-height is too small for non-overlapping tip images. "
                "Increase it or reduce --tip-image-size."
            )
    else:
        tree_panel_height_in = max(
            fig_height - spacing.top_margin_in - spacing.bottom_margin_in, 0.2
        )

    if (
        time_layout.layout_name in {"cladogram", "fractal"}
        and not depths.use_topology_depth
    ):
        sys.stderr.write(
            "{} layout encodes topology; input branch lengths do not determine geometry.\n".format(
                time_layout.layout_name.capitalize(),
            )
        )
    return DrawingPlacement(
        drawing_layout=drawing_layout,
        posterior_layouts=posterior_layouts,
        xcoord=xcoord,
        ycoord=ycoord,
        leaf_order=leaf_order,
        x_min=x_min,
        x_max=x_max,
        y_min=y_min,
        y_max=y_max,
        x_span=x_span,
        y_span=y_span,
        tree_panel_height_in=tree_panel_height_in,
        fig_height=fig_height,
    )
