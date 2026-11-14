"""Typed data exchanged between drawing setup, rendering, and output stages."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D
    from matplotlib.text import Text

    from nwkit.draw_layouts import TreeDrawingLayout
    from nwkit.draw_quality import DrawingArtist


@dataclass(frozen=True)
class FontStyle:
    """Text metrics shared by measurement and rendering."""

    font_size: float
    font_family: str


@dataclass(frozen=True)
class AnnotationOptions:
    """Normalized node annotation options."""

    node_pie_properties: list[str]
    node_pie_target: str
    node_pie_leaf_filters: list[str]
    node_label_property: str | None
    node_label_target: str
    node_label_filters: list[str]
    node_label_decimals: int
    node_label_prefix: str


@dataclass(frozen=True)
class PropertyStyle:
    """Shared categorical color configuration."""

    property_colors: dict
    trait_palette: str


@dataclass(frozen=True)
class TimeDisplay:
    """Validated posterior and time annotation options."""

    densitree_mode: str
    densitree_trees: list
    densitree_alpha: float
    densitree_ci_alpha: float
    densitree_ci_level: float
    constraint_mode: str
    credible_interval_mode: str


@dataclass(frozen=True)
class TimeLayout:
    """Layout compatibility with time annotations."""

    layout_name: str
    topology_only_time_layout: bool
    has_credible_intervals: bool


@dataclass(frozen=True)
class SpatialOptions:
    """Validated layout, spacing, and track dimensions."""

    angular_span: float
    angular_center: float
    resolved_subtree_packing: str
    spatial_layout: bool
    track_size_pt: float
    track_stride_pt: float
    track_span_pt: float
    resolved_tip_spacing: str
    daylight_iterations: int
    resolved_tip_label_position: str


@dataclass(frozen=True)
class TreeDepths:
    """Unscaled tree depths and requested distance guides."""

    use_topology_depth: bool
    base_xcoord: dict
    base_leaf_order: list
    base_x_min: float
    base_x_max: float
    base_x_span: float
    segment_length_layouts: set[str]
    depth_projection_layouts: set[str]
    warped_depth_layouts: set[str]
    branch_length_layouts: set[str]
    branch_depth_span: float
    requested_scale_bar: float | None
    requested_depth_guide: float | None
    depth_ticks: list[float]


@dataclass(frozen=True)
class DrawingPanels:
    """Measured labels, tracks, badges, and panel widths."""

    leaf_order: list
    fig_width: float
    tip_image_size_pt: float
    image_panel_width_in: float
    main_panel_width_in: float
    tip_label_text_by_leaf: dict
    tip_label_size_by_leaf: dict
    max_tip_label_width_in: float
    max_tip_label_height_in: float
    tip_track_color_by_leaf_property: dict
    tip_track_legend_entries: list
    badge_values: list[str]
    badge_value_index: dict[str, int]
    label_panel_width_in: float
    tree_panel_width_in: float


@dataclass(frozen=True)
class DrawingSpacing:
    """Physical row spacing and figure margins."""

    spacing_height_by_leaf: dict
    row_pitch_pt: float
    row_pitch_in: float
    legend_needed: bool
    top_margin_in: float
    depth_guide_strip_pt: float
    bottom_margin_in: float
    left_margin_in: float
    right_margin_in: float
    fig_height: float


@dataclass(frozen=True)
class DrawingPlacement:
    """Final tree coordinates and physical drawing bounds."""

    drawing_layout: TreeDrawingLayout
    posterior_layouts: list
    xcoord: dict
    ycoord: dict
    leaf_order: list
    x_min: float
    x_max: float
    y_min: float
    y_max: float
    x_span: float
    y_span: float
    tree_panel_height_in: float
    fig_height: float


@dataclass(frozen=True)
class RenderFrame:
    """The figure and mutable artist collections for one render."""

    fig: Figure
    ax: Axes
    drawing_artists: list[DrawingArtist]
    branch_lines: list[tuple[Any, Line2D]]


@dataclass(frozen=True)
class TipLabelGeometry:
    """Label clearance and rendered tip text."""

    data_per_inch: float
    label_offset: float
    tip_label_artists: list[Text]


@dataclass(frozen=True)
class AxesPlacement:
    """Figure-relative bounds of the main drawing axes."""

    axes_left: float
    axes_right: float
    axes_top: float
    axes_bottom: float
