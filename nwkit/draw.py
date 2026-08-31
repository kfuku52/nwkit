"""Drawing command orchestration; setup, rendering, and saving have separate stages."""

import os
import sys
from typing import Any

import matplotlib
import matplotlib.pyplot as plt

from nwkit.draw_helpers import (
    BRANCH_COLOR as BRANCH_COLOR,
)
from nwkit.draw_helpers import (
    DUPLICATION_COLOR as DUPLICATION_COLOR,
)
from nwkit.draw_helpers import (
    FIGURE_WIDTH_IN as FIGURE_WIDTH_IN,
)
from nwkit.draw_helpers import (
    FONT_FAMILY as FONT_FAMILY,
)
from nwkit.draw_helpers import (
    FONT_SIZE_PT as FONT_SIZE_PT,
)
from nwkit.draw_helpers import (
    LABEL_COLOR as LABEL_COLOR,
)
from nwkit.draw_helpers import (
    LEGEND_EDGE_COLOR as LEGEND_EDGE_COLOR,
)
from nwkit.draw_helpers import (
    MAX_TIP_IMAGE_EDGE_PX as MAX_TIP_IMAGE_EDGE_PX,
)
from nwkit.draw_helpers import (
    MAX_TIP_IMAGE_SIZE_PT as MAX_TIP_IMAGE_SIZE_PT,
)
from nwkit.draw_helpers import (
    SPECIATION_COLOR as SPECIATION_COLOR,
)
from nwkit.draw_helpers import (
    SUPPORT_LABEL_OFFSET_PT as SUPPORT_LABEL_OFFSET_PT,
)
from nwkit.draw_helpers import (
    TIP_IMAGE_GAP_PT as TIP_IMAGE_GAP_PT,
)
from nwkit.draw_helpers import (
    TIP_IMAGE_SIZE_PT as TIP_IMAGE_SIZE_PT,
)
from nwkit.draw_helpers import (
    TIP_LABEL_GAP_PT as TIP_LABEL_GAP_PT,
)
from nwkit.draw_helpers import (
    TIP_ROW_PITCH_PT as TIP_ROW_PITCH_PT,
)
from nwkit.draw_helpers import (
    TREE_LINE_CAPSTYLE as TREE_LINE_CAPSTYLE,
)
from nwkit.draw_helpers import (
    _add_bottom_guide_title as _add_bottom_guide_title,
)
from nwkit.draw_helpers import (
    _add_radial_depth_guide as _add_radial_depth_guide,
)
from nwkit.draw_helpers import (
    _add_scale_bar as _add_scale_bar,
)
from nwkit.draw_helpers import (
    _add_slanted_depth_guide as _add_slanted_depth_guide,
)
from nwkit.draw_helpers import (
    _add_spiral_depth_key as _add_spiral_depth_key,
)
from nwkit.draw_helpers import (
    _age_interval_path as _age_interval_path,
)
from nwkit.draw_helpers import (
    _angular_sector_bounds as _angular_sector_bounds,
)
from nwkit.draw_helpers import (
    _apply_font_style as _apply_font_style,
)
from nwkit.draw_helpers import (
    _assign_group_colors as _assign_group_colors,
)
from nwkit.draw_helpers import (
    _auto_wrap_candidates as _auto_wrap_candidates,
)
from nwkit.draw_helpers import (
    _branch_marker_artists as _branch_marker_artists,
)
from nwkit.draw_helpers import (
    _branch_style_maps as _branch_style_maps,
)
from nwkit.draw_helpers import (
    _coerce_filter_value as _coerce_filter_value,
)
from nwkit.draw_helpers import (
    _compound_polygon_path as _compound_polygon_path,
)
from nwkit.draw_helpers import (
    _convex_hull as _convex_hull,
)
from nwkit.draw_helpers import (
    _densitree_branch_envelope as _densitree_branch_envelope,
)
from nwkit.draw_helpers import (
    _densitree_edge_key as _densitree_edge_key,
)
from nwkit.draw_helpers import (
    _densitree_topology_groups as _densitree_topology_groups,
)
from nwkit.draw_helpers import (
    _densitree_topology_signature as _densitree_topology_signature,
)
from nwkit.draw_helpers import (
    _depth_guide_ticks as _depth_guide_ticks,
)
from nwkit.draw_helpers import (
    _distance_label as _distance_label,
)
from nwkit.draw_helpers import (
    _draw_branch_marker_overlays as _draw_branch_marker_overlays,
)
from nwkit.draw_helpers import (
    _draw_probability_pie as _draw_probability_pie,
)
from nwkit.draw_helpers import (
    _draw_tip_images as _draw_tip_images,
)
from nwkit.draw_helpers import (
    _first_nonzero_direction as _first_nonzero_direction,
)
from nwkit.draw_helpers import (
    _fit_artists_within_figure as _fit_artists_within_figure,
)
from nwkit.draw_helpers import (
    _format_property_value as _format_property_value,
)
from nwkit.draw_helpers import (
    _format_support_value as _format_support_value,
)
from nwkit.draw_helpers import (
    _get_species_by_leaf as _get_species_by_leaf,
)
from nwkit.draw_helpers import (
    _get_species_overlap_node_types as _get_species_overlap_node_types,
)
from nwkit.draw_helpers import (
    _get_tree_plot_coordinates as _get_tree_plot_coordinates,
)
from nwkit.draw_helpers import (
    _has_meaningful_support as _has_meaningful_support,
)
from nwkit.draw_helpers import (
    _has_positive_branch_length as _has_positive_branch_length,
)
from nwkit.draw_helpers import (
    _isolated_matplotlib_state as _isolated_matplotlib_state,
)
from nwkit.draw_helpers import (
    _legend_is_needed as _legend_is_needed,
)
from nwkit.draw_helpers import (
    _load_tip_image as _load_tip_image,
)
from nwkit.draw_helpers import (
    _load_tip_images as _load_tip_images,
)
from nwkit.draw_helpers import (
    _make_posterior_drawing_layouts as _make_posterior_drawing_layouts,
)
from nwkit.draw_helpers import (
    _make_tree_collection_drawing_layouts as _make_tree_collection_drawing_layouts,
)
from nwkit.draw_helpers import (
    _matches_property_filters as _matches_property_filters,
)
from nwkit.draw_helpers import (
    _measure_texts_in_inches as _measure_texts_in_inches,
)
from nwkit.draw_helpers import (
    _node_matches_target as _node_matches_target,
)
from nwkit.draw_helpers import (
    _normalized_branch_markers as _normalized_branch_markers,
)
from nwkit.draw_helpers import (
    _parse_branch_width_range as _parse_branch_width_range,
)
from nwkit.draw_helpers import (
    _parse_depth_guide as _parse_depth_guide,
)
from nwkit.draw_helpers import (
    _parse_legend_columns as _parse_legend_columns,
)
from nwkit.draw_helpers import (
    _parse_property_colors as _parse_property_colors,
)
from nwkit.draw_helpers import (
    _parse_property_names as _parse_property_names,
)
from nwkit.draw_helpers import (
    _parse_scale_bar as _parse_scale_bar,
)
from nwkit.draw_helpers import (
    _parse_tip_label_wrap as _parse_tip_label_wrap,
)
from nwkit.draw_helpers import (
    _physical_edge_paths as _physical_edge_paths,
)
from nwkit.draw_helpers import (
    _point_on_path as _point_on_path,
)
from nwkit.draw_helpers import (
    _polar_auto_figure_height as _polar_auto_figure_height,
)
from nwkit.draw_helpers import (
    _prepare_tip_label_text as _prepare_tip_label_text,
)
from nwkit.draw_helpers import (
    _property_color as _property_color,
)
from nwkit.draw_helpers import (
    _read_tip_image_manifest as _read_tip_image_manifest,
)
from nwkit.draw_helpers import (
    _read_trait_table as _read_trait_table,
)
from nwkit.draw_helpers import (
    _readable_text_angle as _readable_text_angle,
)
from nwkit.draw_helpers import (
    _remap_leaf_dictionary as _remap_leaf_dictionary,
)
from nwkit.draw_helpers import (
    _resample_polyline as _resample_polyline,
)
from nwkit.draw_helpers import (
    _resolve_image_format as _resolve_image_format,
)
from nwkit.draw_helpers import (
    _sector_contains_direction as _sector_contains_direction,
)
from nwkit.draw_helpers import (
    _spatial_node_label_placement as _spatial_node_label_placement,
)
from nwkit.draw_helpers import (
    _spatial_text_alignment as _spatial_text_alignment,
)
from nwkit.draw_helpers import (
    _taxonomic_prefix as _taxonomic_prefix,
)
from nwkit.draw_helpers import (
    _tip_font_style as _tip_font_style,
)
from nwkit.draw_helpers import (
    _tip_track_artist as _tip_track_artist,
)
from nwkit.draw_helpers import (
    _tip_track_colors as _tip_track_colors,
)
from nwkit.draw_helpers import (
    _tree_drawing_fingerprint as _tree_drawing_fingerprint,
)
from nwkit.draw_helpers import (
    _validated_tip_image_size as _validated_tip_image_size,
)
from nwkit.draw_helpers import (
    _wrap_taxonomic_label as _wrap_taxonomic_label,
)
from nwkit.draw_helpers import (
    _wrap_tip_label as _wrap_tip_label,
)
from nwkit.draw_output import save_drawing
from nwkit.draw_prep import collapse_tree_for_drawing
from nwkit.draw_render import (
    _add_drawing_report_metadata,
    _configure_drawing_axes,
    _draw_densitree_layers,
    _draw_distance_guides,
    _draw_initial_depth_guides,
    _draw_node_annotations,
    _draw_node_type_markers,
    _draw_root_marker,
    _draw_support_labels,
    _draw_time_constraints,
    _draw_time_intervals,
    _draw_tip_image_panel,
    _draw_tip_labels,
    _draw_tree_branches,
    _draw_tree_legend,
    _finalize_drawing_quality,
    _property_legend_handles,
)
from nwkit.draw_setup import (
    _place_drawing,
    _prepare_drawing_panels,
    _prepare_drawing_spacing,
    _prepare_tree_depths,
    _resolve_spatial_options,
    _resolve_time_display,
    _validate_node_annotation_options,
    _validate_time_layout,
)
from nwkit.draw_types import (
    AnnotationOptions,
    FontStyle,
    PropertyStyle,
    RenderFrame,
)
from nwkit.file_paths import (
    validate_distinct_output_paths,
    validate_outputs_do_not_replace_inputs,
)
from nwkit.output_transaction import validate_output_targets
from nwkit.time_tree import (
    prepare_time_tree_annotations,
    read_mcmctree_posterior,
    summarize_mcmctree_posterior,
)
from nwkit.util import (
    is_rooted,
    read_tree,
    read_trees,
    validate_unique_named_leaves,
)


@_isolated_matplotlib_state
def _draw_tree(
    tree,
    outfile,
    image_format,
    node_type_by_node,
    leaf_label_color_by_leaf,
    group_color_by_name,
    support_labels=True,
    support_min=None,
    figure_width=FIGURE_WIDTH_IN,
    figure_height=None,
    label_panel_width=None,
    font_size=FONT_SIZE_PT,
    font_family=FONT_FAMILY,
    branch_color=BRANCH_COLOR,
    terminal_branch_color=None,
    tip_labels=True,
    tip_label_position="aligned",
    tip_label_wrap="none",
    tip_spacing="uniform",
    subtree_packing="standard",
    root_marker="none",
    root_marker_color="#0072B2",
    root_marker_size=None,
    tip_badge_property=None,
    tip_badge_missing_label=None,
    node_pie_properties=None,
    node_pie_target="root,intnode",
    node_pie_leaf_filters=None,
    node_label_property=None,
    node_label_target="intnode",
    node_label_filters=None,
    node_label_decimals=2,
    node_label_prefix="",
    property_colors=None,
    legend=True,
    transparent=False,
    trait_palette="tab10",
    tip_image_by_leaf=None,
    tip_image_size=TIP_IMAGE_SIZE_PT,
    tip_image_gap=TIP_IMAGE_GAP_PT,
    layout="rectangular",
    spiral_turns=None,
    angular_span=360.0,
    angular_center=90.0,
    unrooted_method="equal-angle",
    daylight_iterations=5,
    scale_bar="none",
    depth_guide="none",
    branch_length_unit="",
    branch_width=0.8,
    branch_color_property=None,
    branch_width_property=None,
    branch_width_range="0.4,2.5",
    branch_markers=None,
    tip_label_font_style="plain",
    tip_track_properties=None,
    tip_track_type="auto",
    tip_track_size=5.0,
    tip_track_palette="viridis",
    legend_columns="auto",
    legend_position="auto",
    collision_policy="resolve",
    layout_report=None,
    collapsed_clades=None,
    time_constraints="auto",
    time_credible_intervals="auto",
    mcmctree_posterior=None,
    mcmctree_posterior_topology=None,
    densitree_trees=None,
    posterior_ladderize=False,
    densitree="none",
    densitree_alpha=0.035,
    densitree_color="#0072B2",
    densitree_ci_level=0.95,
    densitree_ci_alpha=0.18,
    densitree_ci_color="#56B4E9",
    transactional_output=True,
):
    _apply_font_style(font_size=font_size, font_family=font_family)
    matplotlib.rcParams["svg.hashsalt"] = "nwkit"
    property_colors = property_colors or {}
    node_pie_properties = node_pie_properties or []
    node_pie_leaf_filters = node_pie_leaf_filters or []
    node_label_filters = node_label_filters or []
    tip_image_by_leaf = tip_image_by_leaf or {}
    tip_track_properties = tip_track_properties or []
    branch_markers = _normalized_branch_markers(branch_markers)
    collapsed_clades = collapsed_clades or []
    font_style = FontStyle(
        font_size=font_size,
        font_family=font_family,
    )
    annotations = AnnotationOptions(
        node_pie_properties=node_pie_properties,
        node_pie_target=node_pie_target,
        node_pie_leaf_filters=node_pie_leaf_filters,
        node_label_property=node_label_property,
        node_label_target=node_label_target,
        node_label_filters=node_label_filters,
        node_label_decimals=node_label_decimals,
        node_label_prefix=node_label_prefix,
    )
    property_style = PropertyStyle(
        property_colors=property_colors,
        trait_palette=trait_palette,
    )
    time_display = _resolve_time_display(
        densitree=densitree,
        densitree_alpha=densitree_alpha,
        densitree_ci_alpha=densitree_ci_alpha,
        densitree_ci_level=densitree_ci_level,
        densitree_trees=densitree_trees,
        mcmctree_posterior=mcmctree_posterior,
        time_constraints=time_constraints,
        time_credible_intervals=time_credible_intervals,
    )

    _validate_node_annotation_options(
        annotations=annotations,
        tree=tree,
    )

    time_layout = _validate_time_layout(
        time_display=time_display,
        layout=layout,
        tree=tree,
    )

    spatial_options = _resolve_spatial_options(
        time_layout=time_layout,
        angular_center=angular_center,
        angular_span=angular_span,
        daylight_iterations=daylight_iterations,
        label_panel_width=label_panel_width,
        subtree_packing=subtree_packing,
        tip_image_by_leaf=tip_image_by_leaf,
        tip_label_font_style=tip_label_font_style,
        tip_label_position=tip_label_position,
        tip_spacing=tip_spacing,
        tip_track_properties=tip_track_properties,
        tip_track_size=tip_track_size,
    )

    depths = _prepare_tree_depths(
        time_layout=time_layout,
        depth_guide=depth_guide,
        scale_bar=scale_bar,
        tree=tree,
    )

    panels = _prepare_drawing_panels(
        depths=depths,
        font_style=font_style,
        property_style=property_style,
        spatial_options=spatial_options,
        time_layout=time_layout,
        figure_height=figure_height,
        figure_width=figure_width,
        label_panel_width=label_panel_width,
        tip_badge_missing_label=tip_badge_missing_label,
        tip_badge_property=tip_badge_property,
        tip_image_by_leaf=tip_image_by_leaf,
        tip_image_gap=tip_image_gap,
        tip_image_size=tip_image_size,
        tip_label_wrap=tip_label_wrap,
        tip_labels=tip_labels,
        tip_track_palette=tip_track_palette,
        tip_track_properties=tip_track_properties,
        tip_track_type=tip_track_type,
        tree=tree,
    )

    spacing = _prepare_drawing_spacing(
        annotations=annotations,
        depths=depths,
        font_style=font_style,
        panels=panels,
        spatial_options=spatial_options,
        time_layout=time_layout,
        branch_color_property=branch_color_property,
        branch_markers=branch_markers,
        branch_width_property=branch_width_property,
        figure_height=figure_height,
        group_color_by_name=group_color_by_name,
        legend=legend,
        node_type_by_node=node_type_by_node,
        tip_badge_property=tip_badge_property,
        tip_image_by_leaf=tip_image_by_leaf,
        tip_labels=tip_labels,
        tip_track_properties=tip_track_properties,
    )

    geometry = _place_drawing(
        depths=depths,
        panels=panels,
        spacing=spacing,
        spatial_options=spatial_options,
        time_display=time_display,
        time_layout=time_layout,
        figure_height=figure_height,
        mcmctree_posterior=mcmctree_posterior,
        mcmctree_posterior_topology=mcmctree_posterior_topology,
        posterior_ladderize=posterior_ladderize,
        spiral_turns=spiral_turns,
        tip_image_by_leaf=tip_image_by_leaf,
        tip_labels=tip_labels,
        tree=tree,
        unrooted_method=unrooted_method,
    )

    color_by_node, width_by_node = _branch_style_maps(
        tree=tree,
        base_color=branch_color,
        base_width=branch_width,
        color_property=branch_color_property,
        width_property=branch_width_property,
        width_range=branch_width_range,
        property_colors=property_style.property_colors,
        palette=property_style.trait_palette,
    )
    fig, ax = plt.subplots(figsize=(panels.fig_width, geometry.fig_height))
    frame = RenderFrame(fig=fig, ax=ax, drawing_artists=[], branch_lines=[])
    depth_guide_color, radial_depth_labels_drawn = _draw_initial_depth_guides(
        depths=depths,
        font_style=font_style,
        frame=frame,
        geometry=geometry,
        time_layout=time_layout,
        tree=tree,
    )

    densitree_ci_group_report = _draw_densitree_layers(
        frame=frame,
        geometry=geometry,
        time_display=time_display,
        branch_width=branch_width,
        densitree_ci_color=densitree_ci_color,
        densitree_color=densitree_color,
    )

    _draw_tree_branches(
        frame=frame,
        geometry=geometry,
        branch_color=branch_color,
        branch_color_property=branch_color_property,
        branch_width=branch_width,
        color_by_node=color_by_node,
        terminal_branch_color=terminal_branch_color,
        tree=tree,
        width_by_node=width_by_node,
    )

    credible_interval_count = _draw_time_intervals(
        font_style=font_style,
        frame=frame,
        geometry=geometry,
        panels=panels,
        time_display=time_display,
        time_layout=time_layout,
        branch_width=branch_width,
        tree=tree,
    )

    constraint_artist_count = _draw_time_constraints(
        font_style=font_style,
        frame=frame,
        geometry=geometry,
        spatial_options=spatial_options,
        time_display=time_display,
        tree=tree,
    )

    legend_handles = _draw_node_type_markers(
        frame=frame,
        geometry=geometry,
        branch_markers=branch_markers,
        group_color_by_name=group_color_by_name,
        node_type_by_node=node_type_by_node,
        tree=tree,
    )

    legend_handles = _property_legend_handles(
        annotations=annotations,
        panels=panels,
        property_style=property_style,
        branch_color=branch_color,
        branch_color_property=branch_color_property,
        branch_width=branch_width,
        branch_width_property=branch_width_property,
        legend_handles=legend_handles,
        tip_badge_property=tip_badge_property,
        tree=tree,
        width_by_node=width_by_node,
    )

    _draw_tree_legend(
        font_style=font_style,
        frame=frame,
        legend=legend,
        legend_columns=legend_columns,
        legend_handles=legend_handles,
        legend_position=legend_position,
    )

    _draw_support_labels(
        font_style=font_style,
        frame=frame,
        geometry=geometry,
        spatial_options=spatial_options,
        support_labels=support_labels,
        support_min=support_min,
        tree=tree,
    )

    label_geometry = _draw_tip_labels(
        annotations=annotations,
        font_style=font_style,
        frame=frame,
        geometry=geometry,
        panels=panels,
        property_style=property_style,
        spatial_options=spatial_options,
        branch_color=branch_color,
        leaf_label_color_by_leaf=leaf_label_color_by_leaf,
        tip_badge_missing_label=tip_badge_missing_label,
        tip_badge_property=tip_badge_property,
        tip_label_font_style=tip_label_font_style,
        tip_label_position=tip_label_position,
        tip_labels=tip_labels,
        tip_track_properties=tip_track_properties,
    )

    _draw_node_annotations(
        annotations=annotations,
        font_style=font_style,
        frame=frame,
        geometry=geometry,
        property_style=property_style,
        spatial_options=spatial_options,
        tree=tree,
    )

    root_marker_size_pt = _draw_root_marker(
        font_style=font_style,
        frame=frame,
        geometry=geometry,
        root_marker=root_marker,
        root_marker_color=root_marker_color,
        root_marker_size=root_marker_size,
        tree=tree,
    )

    axes_placement = _configure_drawing_axes(
        annotations=annotations,
        font_style=font_style,
        frame=frame,
        geometry=geometry,
        label_geometry=label_geometry,
        panels=panels,
        spacing=spacing,
        spatial_options=spatial_options,
        root_marker=root_marker,
        root_marker_size_pt=root_marker_size_pt,
        tip_badge_property=tip_badge_property,
        tip_labels=tip_labels,
        tip_track_properties=tip_track_properties,
        tree=tree,
    )

    _draw_distance_guides(
        axes_placement=axes_placement,
        depths=depths,
        font_style=font_style,
        frame=frame,
        geometry=geometry,
        spacing=spacing,
        time_layout=time_layout,
        branch_color=branch_color,
        branch_length_unit=branch_length_unit,
        depth_guide_color=depth_guide_color,
    )

    _draw_tip_image_panel(
        axes_placement=axes_placement,
        frame=frame,
        geometry=geometry,
        panels=panels,
        tip_image_by_leaf=tip_image_by_leaf,
    )

    quality_report, fit_report = _finalize_drawing_quality(
        frame=frame,
        collision_policy=collision_policy,
    )

    _add_drawing_report_metadata(
        depths=depths,
        font_style=font_style,
        geometry=geometry,
        panels=panels,
        spatial_options=spatial_options,
        time_display=time_display,
        time_layout=time_layout,
        branch_length_unit=branch_length_unit,
        collapsed_clades=collapsed_clades,
        constraint_artist_count=constraint_artist_count,
        credible_interval_count=credible_interval_count,
        densitree_ci_group_report=densitree_ci_group_report,
        fit_report=fit_report,
        image_format=image_format,
        mcmctree_posterior=mcmctree_posterior,
        quality_report=quality_report,
        radial_depth_labels_drawn=radial_depth_labels_drawn,
        tip_track_palette=tip_track_palette,
        tip_track_properties=tip_track_properties,
        tree=tree,
        unrooted_method=unrooted_method,
    )

    save_drawing(
        frame.fig,
        outfile,
        image_format,
        transparent,
        layout_report,
        quality_report,
        transactional=transactional_output,
    )


def _validate_draw_paths(outfile, layout_report, input_paths):
    outputs = [("--outfile", outfile), ("--layout-report", layout_report)]
    validate_distinct_output_paths(outputs)
    validate_outputs_do_not_replace_inputs(
        list(input_paths.items()), outputs, label="Drawing output"
    )
    validate_output_targets(
        [path for _, path in outputs if path not in (None, "", "-")]
    )


def draw_main(args):
    if args.outfile == "-":
        raise ValueError("STDOUT is not supported for 'draw'. Use --outfile PATH.")
    trait_path = getattr(args, "trait", None)
    tip_image_manifest = getattr(args, "tip_image_manifest", None)
    posterior_path = getattr(args, "mcmctree_posterior", None)
    densitree_tree_path = getattr(args, "densitree_trees", None)
    layout_report = getattr(args, "layout_report", None)
    _validate_draw_paths(
        outfile=args.outfile,
        layout_report=layout_report,
        input_paths={
            "--infile": args.infile,
            "--trait": trait_path,
            "--tip-image-manifest": tip_image_manifest,
            "--mcmctree-posterior": posterior_path,
            "--densitree-trees": densitree_tree_path,
            "--species-map-tsv": getattr(args, "species_map_tsv", None),
        },
    )
    tree = read_tree(
        args.infile,
        args.format,
        args.quoted_node_names,
        rooted=getattr(args, "input_rooted", "auto"),
    )
    prepare_time_tree_annotations(tree)
    densitree_mode = str(getattr(args, "densitree", "none")).strip().lower()
    if posterior_path not in (None, "") and densitree_tree_path not in (None, ""):
        raise ValueError(
            "'--mcmctree-posterior' and '--densitree-trees' are mutually exclusive."
        )
    if densitree_tree_path not in (None, "") and densitree_mode == "none":
        raise ValueError("'--densitree-trees' requires '--densitree all, ci, or both'.")
    if (
        densitree_mode != "none"
        and str(getattr(args, "layout", "rectangular")).strip().lower()
        in {"cladogram", "fractal"}
        and densitree_tree_path in (None, "")
    ):
        raise ValueError(
            "'--densitree' requires a layout that encodes branch-length or "
            "root-to-node time variation; cladogram and fractal do not."
        )
    posterior = None
    posterior_topology = None
    densitree_sample_trees = []
    if posterior_path not in (None, ""):
        if posterior_path == "-" and args.infile == "-":
            raise ValueError(
                "'--infile' and '--mcmctree-posterior' cannot both read from STDIN."
            )
        posterior = read_mcmctree_posterior(
            posterior_path,
            tree=tree,
            burnin=getattr(args, "posterior_burnin", 0),
            thin=getattr(args, "posterior_thin", 1),
        )
        summarize_mcmctree_posterior(
            tree,
            posterior=posterior,
            point=getattr(args, "posterior_point", "mean"),
            ci=getattr(args, "posterior_ci", "hpd"),
            level=getattr(args, "posterior_ci_level", 0.95),
        )
        if densitree_mode != "none":
            posterior_topology = tree.copy()
            sys.stderr.write(
                "Rendering {:,} posterior dated trees reconstructed from MCMCtree samples.\n".format(
                    posterior.sample_count,
                )
            )
    elif densitree_mode != "none":
        if densitree_tree_path in (None, ""):
            raise ValueError(
                "'--densitree' requires '--mcmctree-posterior' or '--densitree-trees'."
            )
    if densitree_tree_path not in (None, ""):
        if densitree_tree_path == "-" and args.infile == "-":
            raise ValueError(
                "'--infile' and '--densitree-trees' cannot both read from STDIN."
            )
        burnin = int(getattr(args, "posterior_burnin", 0))
        thin = int(getattr(args, "posterior_thin", 1))
        if burnin < 0:
            raise ValueError("'--posterior-burnin' must be zero or greater.")
        if thin < 1:
            raise ValueError("'--posterior-thin' must be at least one.")
        all_sample_trees = read_trees(
            densitree_tree_path,
            args.format,
            args.quoted_node_names,
            quiet=True,
            rooted=getattr(args, "densitree_trees_rooted", "auto"),
        )
        densitree_sample_trees = all_sample_trees[burnin::thin]
        if not densitree_sample_trees:
            raise ValueError(
                "'--posterior-burnin' and '--posterior-thin' retained no "
                "DensiTree sample trees."
            )
        validate_unique_named_leaves(
            tree,
            option_name="--infile",
            context=" for '--densitree-trees'",
        )
        reference_tip_set = set(str(name) for name in tree.leaf_names())
        for index, sample_tree in enumerate(densitree_sample_trees, start=1):
            validate_unique_named_leaves(
                sample_tree,
                option_name="--densitree-trees",
                context=" in retained tree {}".format(index),
            )
            sample_tip_set = set(str(name) for name in sample_tree.leaf_names())
            if sample_tip_set != reference_tip_set:
                missing = sorted(reference_tip_set - sample_tip_set)
                extra = sorted(sample_tip_set - reference_tip_set)
                raise ValueError(
                    "'--densitree-trees' retained tree {} has a different "
                    "tip set (missing: {}; extra: {}).".format(
                        index,
                        ",".join(missing) or "none",
                        ",".join(extra) or "none",
                    )
                )
            prepare_time_tree_annotations(sample_tree)
            if bool(args.ladderize):
                sample_tree.ladderize()
        sys.stderr.write(
            "Rendering {:,} retained trees representing {:,} rooted topologies.\n".format(
                len(densitree_sample_trees),
                len(
                    {
                        _densitree_topology_signature(sample_tree)
                        for sample_tree in densitree_sample_trees
                    }
                ),
            )
        )
    if bool(args.ladderize):
        tree.ladderize()
    max_visible_tips = getattr(args, "max_visible_tips", None)
    if max_visible_tips not in (None, "") and int(max_visible_tips) < 2:
        raise ValueError("'--max-visible-tips' must be at least 2.")
    collapse_required = max_visible_tips not in (None, "") and len(
        list(tree.leaves())
    ) > int(max_visible_tips)
    if collapse_required and densitree_mode != "none":
        raise ValueError(
            "'--max-visible-tips' cannot currently be combined with '--densitree'."
        )
    if collapse_required and (
        trait_path not in (None, "") or tip_image_manifest not in (None, "")
    ):
        raise ValueError(
            "'--max-visible-tips' cannot currently be combined with "
            "'--trait' or '--tip-image-manifest'."
        )
    tree, collapsed_clades = collapse_tree_for_drawing(
        tree,
        max_visible_tips=max_visible_tips,
        label_template=getattr(args, "collapse_label", None),
        property_aggregation=getattr(
            args,
            "collapse_property_aggregation",
            "none",
        ),
    )
    group_by = getattr(args, "group_by", None)
    trait_palette = getattr(args, "trait_palette", "tab10")
    support_labels = getattr(args, "support_labels", True)
    support_min = getattr(args, "support_min", None)
    property_colors = _parse_property_colors(getattr(args, "property_color", None))
    node_pie_properties = _parse_property_names(
        getattr(args, "node_pie_properties", None)
    )
    tip_image_path_by_leaf = _read_tip_image_manifest(
        path=tip_image_manifest,
        tree=tree,
        image_root=getattr(args, "tip_image_root", None),
        unmatched=getattr(args, "unmatched", "warn"),
        missing_values=getattr(args, "missing_values", None),
    )
    args._nwkit_tip_image_paths = sorted(set(tip_image_path_by_leaf.values()))
    _validate_draw_paths(
        outfile=args.outfile,
        layout_report=layout_report,
        input_paths={
            "--infile": args.infile,
            "--trait": trait_path,
            "--tip-image-manifest": tip_image_manifest,
            "--mcmctree-posterior": posterior_path,
            "--densitree-trees": densitree_tree_path,
            "--species-map-tsv": getattr(args, "species_map_tsv", None),
            **{
                "--tip-image '{}'".format(leaf_name): path
                for leaf_name, path in tip_image_path_by_leaf.items()
            },
        },
    )
    tip_image_by_leaf = _load_tip_images(
        tip_image_path_by_leaf,
        image_size_pt=getattr(args, "tip_image_size", 18.0),
    )
    if tip_image_by_leaf:
        sys.stderr.write(
            "Loaded tip images for {} tree tip(s) from: {}\n".format(
                len(tip_image_by_leaf),
                tip_image_manifest,
            )
        )
    if (trait_path not in ["", None]) and (group_by in ["", None]):
        raise ValueError("'--group-by' is required when '--trait' is specified.")
    image_format = _resolve_image_format(
        outfile=args.outfile, image_format=str(args.image_format).lower()
    )
    node_plot_mode = str(args.species_overlap_node_plot).strip().lower()
    if node_plot_mode == "no":
        node_type_by_node: dict[Any, str] = {}
    else:
        if not is_rooted(tree):
            if node_plot_mode == "yes":
                raise ValueError(
                    "Speciation/duplication node plotting requires a rooted tree. "
                    "Use '--species-overlap-node-plot no' for unrooted trees."
                )
            node_type_by_node = dict()
            sys.stderr.write(
                "Skipping speciation/duplication node markers because the input tree is unrooted.\n"
            )
        else:
            node_type_by_node, all_tip_labels_parsed = _get_species_overlap_node_types(
                tree=tree,
                args=args,
                require_all_tip_labels=(node_plot_mode == "auto"),
            )
            if (node_plot_mode == "auto") and (not all_tip_labels_parsed):
                sys.stderr.write(
                    "Skipping speciation/duplication node markers because some leaf labels did not match the configured species parser.\n"
                )
    leaf_label_color_by_leaf = dict()
    group_color_by_name = dict()
    if trait_path not in ["", None]:
        leaf_to_group = _read_trait_table(
            path=trait_path,
            group_by=group_by,
            tree=tree,
            unmatched=getattr(args, "unmatched", "warn"),
            missing_values=getattr(args, "missing_values", None),
        )
        group_color_by_name = _assign_group_colors(
            group_names=set(leaf_to_group.values()),
            palette=trait_palette,
        )
        for leaf in tree.leaves():
            group_name = leaf_to_group.get(str(leaf.name))
            if group_name is None:
                continue
            leaf_label_color_by_leaf[leaf] = group_color_by_name[group_name]
    outdir = os.path.dirname(os.path.realpath(args.outfile))
    os.makedirs(outdir, exist_ok=True)
    _draw_tree(
        tree=tree,
        outfile=args.outfile,
        image_format=image_format,
        node_type_by_node=node_type_by_node,
        leaf_label_color_by_leaf=leaf_label_color_by_leaf,
        group_color_by_name=group_color_by_name,
        support_labels=support_labels,
        support_min=support_min,
        figure_width=getattr(args, "figure_width", FIGURE_WIDTH_IN),
        figure_height=getattr(args, "figure_height", None),
        label_panel_width=getattr(args, "label_panel_width", None),
        font_size=getattr(args, "font_size", FONT_SIZE_PT),
        font_family=getattr(args, "font_family", FONT_FAMILY),
        branch_color=getattr(args, "branch_color", BRANCH_COLOR),
        terminal_branch_color=getattr(args, "terminal_branch_color", None),
        tip_labels=getattr(args, "tip_labels", True),
        tip_label_position=getattr(args, "tip_label_position", "aligned"),
        tip_label_wrap=getattr(args, "tip_label_wrap", "none"),
        tip_spacing=getattr(args, "tip_spacing", "uniform"),
        subtree_packing=getattr(args, "subtree_packing", "standard"),
        root_marker=getattr(args, "root_marker", "none"),
        root_marker_color=getattr(args, "root_marker_color", "#0072B2"),
        root_marker_size=getattr(args, "root_marker_size", None),
        tip_badge_property=getattr(args, "tip_badge_property", None),
        tip_badge_missing_label=getattr(args, "tip_badge_missing_label", None),
        node_pie_properties=node_pie_properties,
        node_pie_target=getattr(args, "node_pie_target", "root,intnode"),
        node_pie_leaf_filters=getattr(args, "node_pie_leaf_filter", None),
        node_label_property=getattr(args, "node_label_property", None),
        node_label_target=getattr(args, "node_label_target", "intnode"),
        node_label_filters=getattr(args, "node_label_filter", None),
        node_label_decimals=getattr(args, "node_label_decimals", 2),
        node_label_prefix=getattr(args, "node_label_prefix", ""),
        property_colors=property_colors,
        legend=getattr(args, "legend", True),
        transparent=getattr(args, "transparent", False),
        trait_palette=trait_palette,
        tip_image_by_leaf=tip_image_by_leaf,
        tip_image_size=getattr(args, "tip_image_size", TIP_IMAGE_SIZE_PT),
        tip_image_gap=getattr(args, "tip_image_gap", TIP_IMAGE_GAP_PT),
        layout=getattr(args, "layout", "rectangular"),
        spiral_turns=getattr(args, "spiral_turns", None),
        angular_span=getattr(args, "angular_span", 360.0),
        angular_center=getattr(args, "angular_center", 90.0),
        unrooted_method=getattr(args, "unrooted_method", "equal-angle"),
        daylight_iterations=getattr(args, "daylight_iterations", 5),
        scale_bar=getattr(args, "scale_bar", "none"),
        depth_guide=getattr(args, "depth_guide", "none"),
        branch_length_unit=getattr(args, "branch_length_unit", ""),
        branch_width=getattr(args, "branch_width", 0.8),
        branch_color_property=getattr(args, "branch_color_property", None),
        branch_width_property=getattr(args, "branch_width_property", None),
        branch_width_range=getattr(args, "branch_width_range", "0.4,2.5"),
        tip_label_font_style=getattr(args, "tip_label_font_style", "plain"),
        tip_track_properties=getattr(args, "tip_track", None),
        tip_track_type=getattr(args, "tip_track_type", "auto"),
        tip_track_size=getattr(args, "tip_track_size", 5.0),
        tip_track_palette=getattr(args, "tip_track_palette", "viridis"),
        legend_columns=getattr(args, "legend_columns", "auto"),
        legend_position=getattr(args, "legend_position", "auto"),
        collision_policy=getattr(args, "collision_policy", "resolve"),
        layout_report=layout_report,
        collapsed_clades=collapsed_clades,
        time_constraints=getattr(args, "time_constraints", "auto"),
        time_credible_intervals=getattr(args, "time_credible_intervals", "auto"),
        mcmctree_posterior=posterior,
        mcmctree_posterior_topology=posterior_topology,
        densitree_trees=densitree_sample_trees,
        posterior_ladderize=bool(args.ladderize),
        densitree=densitree_mode,
        densitree_alpha=getattr(args, "densitree_alpha", 0.035),
        densitree_color=getattr(args, "densitree_color", "#0072B2"),
        densitree_ci_level=getattr(args, "densitree_ci_level", 0.95),
        densitree_ci_alpha=getattr(args, "densitree_ci_alpha", 0.18),
        densitree_ci_color=getattr(args, "densitree_ci_color", "#56B4E9"),
    )
    sys.stderr.write("Wrote tree image: {}\n".format(os.path.realpath(args.outfile)))
