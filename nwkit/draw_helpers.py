"""Reusable drawing primitives, annotation helpers, and input readers."""

import bisect
import hashlib
import math
import os
import sys
from typing import Any

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colors as mcolors
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnnotationBbox, DrawingArea, OffsetImage
from matplotlib.patches import Arc, Circle, Rectangle, Wedge
from matplotlib.path import Path as MatplotlibPath
from matplotlib.text import Text
from matplotlib.transforms import ScaledTranslation
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

from nwkit.clade_mapping import canonical_split
from nwkit.draw_layouts import get_rectangular_coordinates, make_tree_layout
from nwkit.draw_quality import DrawingArtist
from nwkit.time_tree import (
    posterior_sample_tree,
)
from nwkit.util import (
    extract_species_label,
    get_subtree_leaf_name_sets,
    read_tip_table,
    validate_unique_named_leaves,
)

TREE_LINE_CAPSTYLE = "round"
FONT_FAMILY = "Helvetica"
FONT_SIZE_PT = 8.0
TIP_LABEL_GAP_PT = 0.5
TIP_ROW_PITCH_PT = FONT_SIZE_PT + TIP_LABEL_GAP_PT
TIP_IMAGE_SIZE_PT = 18.0
TIP_IMAGE_GAP_PT = 4.0
MAX_TIP_IMAGE_SIZE_PT = 512.0
MAX_TIP_IMAGE_EDGE_PX = 2048
SUPPORT_LABEL_OFFSET_PT = 0.0
FIGURE_WIDTH_IN = 3.6
SPECIATION_COLOR = (0.0, 0.0, 1.0)
DUPLICATION_COLOR = (1.0, 0.0, 0.0)
BRANCH_COLOR = "#000000"
LABEL_COLOR = "#2d2d2d"
LEGEND_EDGE_COLOR = "white"


def _resolve_image_format(outfile, image_format):
    if image_format != "auto":
        return image_format
    ext = os.path.splitext(str(outfile))[1].lower()
    if ext in (".pdf", ".png", ".svg"):
        return ext[1:]
    return "pdf"


def _tree_drawing_fingerprint(tree):
    """Hash ordered topology and node properties without recursive Newick I/O."""

    digest = hashlib.sha256()

    def add(value):
        encoded = str(value).encode("utf-8", errors="backslashreplace")
        digest.update(str(len(encoded)).encode("ascii"))
        digest.update(b":")
        digest.update(encoded)

    stack = [("node", tree)]
    while stack:
        event, value = stack.pop()
        if event == "close":
            digest.update(b")")
            continue
        node = value
        digest.update(b"(")
        for key in sorted(node.props):
            add(key)
            add(node.props[key])
        stack.append(("close", None))
        stack.extend(("node", child) for child in reversed(node.get_children()))
    return digest.hexdigest()


def _has_positive_branch_length(tree):
    for node in tree.traverse():
        if node.is_root:
            continue
        if (node.dist is not None) and (float(node.dist) > 0):
            return True
    return False


def _get_tree_plot_coordinates(tree, use_topology_depth=False):
    return get_rectangular_coordinates(
        tree,
        use_topology_depth=use_topology_depth,
    )


def _get_species_by_leaf(tree, args):
    species_by_leaf = dict()
    is_all_parsed = True
    num_leaf = 0
    for leaf in tree.leaves():
        num_leaf += 1
        species_label = extract_species_label(leaf.name or "", args=args)
        if species_label is None:
            is_all_parsed = False
            continue
        species_by_leaf[leaf] = species_label
    return species_by_leaf, (is_all_parsed and (num_leaf > 0))


def _get_species_overlap_node_types(tree, args, require_all_tip_labels=False):
    species_by_leaf, all_tip_labels_parsed = _get_species_by_leaf(tree=tree, args=args)
    if bool(require_all_tip_labels) and (not all_tip_labels_parsed):
        return dict(), all_tip_labels_parsed
    species_to_bit: dict[str, int] = {}
    species_mask_by_node = dict()
    has_missing_species_by_node = dict()
    node_type_by_node = dict()
    for node in tree.traverse(strategy="postorder"):
        if node.is_leaf:
            species_label = species_by_leaf.get(node)
            if species_label is None:
                species_mask_by_node[node] = 0
                has_missing_species_by_node[node] = True
                continue
            bit = species_to_bit.get(species_label)
            if bit is None:
                bit = 1 << len(species_to_bit)
                species_to_bit[species_label] = bit
            species_mask_by_node[node] = bit
            has_missing_species_by_node[node] = False
            continue
        children = node.get_children()
        missing_species = False
        union_mask = 0
        child_species_count_sum = 0
        for child in children:
            child_mask = species_mask_by_node[child]
            union_mask |= child_mask
            child_species_count_sum += child_mask.bit_count()
            if has_missing_species_by_node[child]:
                missing_species = True
        species_mask_by_node[node] = union_mask
        has_missing_species_by_node[node] = missing_species
        if missing_species:
            continue
        if len(children) < 2:
            continue
        is_duplication = child_species_count_sum > union_mask.bit_count()
        node_type_by_node[node] = "duplication" if is_duplication else "speciation"
    return node_type_by_node, all_tip_labels_parsed


def _apply_font_style(font_size=FONT_SIZE_PT, font_family=FONT_FAMILY):
    matplotlib.rcParams["font.size"] = font_size
    matplotlib.rc("xtick", labelsize=font_size)
    matplotlib.rc("ytick", labelsize=font_size)
    matplotlib.rc("font", size=font_size, family=font_family)
    matplotlib.rc("axes", titlesize=font_size, labelsize=font_size)
    matplotlib.rc("legend", fontsize=font_size)
    matplotlib.rc("figure", titlesize=font_size)


def _has_meaningful_support(node):
    if node.is_root or node.is_leaf:
        return False
    support = node.support
    if support is None:
        return False
    support_value = float(support)
    if abs(support_value - (-999999.0)) < 10**-9:
        return False
    return True


def _format_support_value(support):
    support_value = float(support)
    if abs(support_value - round(support_value)) < 10**-9:
        return str(int(round(support_value)))
    return "{:g}".format(support_value)


def _read_trait_table(path, group_by, tree, unmatched="warn", missing_values=None):
    trait_df, _, _ = read_tip_table(
        path,
        option_name="--trait",
        tree_leaf_names=tree.leaf_names(),
        required_columns=(group_by,),
        unmatched=unmatched,
        missing_values=missing_values,
    )
    tree_leaf_name_set = set(tree.leaf_names())
    trait_df = trait_df[trait_df["leaf_name"].isin(tree_leaf_name_set)].copy()
    leaf_to_group = dict()
    for _, row in trait_df.iterrows():
        if pd.isna(row[group_by]):
            continue
        leaf_to_group[str(row["leaf_name"])] = str(row[group_by])
    return leaf_to_group


def _read_tip_image_manifest(
    path,
    tree,
    image_root=None,
    unmatched="warn",
    missing_values=None,
):
    if path in (None, ""):
        return dict()
    if path == "-" and image_root in (None, ""):
        raise ValueError(
            "'--tip-image-root' is required when '--tip-image-manifest' reads from STDIN."
        )
    validate_unique_named_leaves(
        tree,
        option_name="--infile",
        context=" for '--tip-image-manifest'",
    )
    manifest_df, _, _ = read_tip_table(
        path,
        option_name="--tip-image-manifest",
        tree_leaf_names=tree.leaf_names(),
        required_columns=("local_path",),
        unmatched=unmatched,
        missing_values=missing_values,
        duplicate_leaf_names="first",
    )
    if image_root in (None, ""):
        base_dir = os.path.dirname(os.path.realpath(path))
    else:
        base_dir = os.path.realpath(image_root)
    if not os.path.isdir(base_dir):
        raise ValueError("'--tip-image-root' is not a directory: {}".format(base_dir))

    tree_leaf_names = set(str(name) for name in tree.leaf_names())
    path_by_leaf = dict()
    missing_path_leaf_names = list()
    for _, row in manifest_df.iterrows():
        leaf_name = str(row["leaf_name"])
        if leaf_name not in tree_leaf_names:
            continue
        raw_local_path = row["local_path"]
        if pd.isna(raw_local_path) or str(raw_local_path).strip() == "":
            missing_path_leaf_names.append(leaf_name)
            continue
        local_path = os.path.expanduser(str(raw_local_path).strip())
        resolved_path = (
            os.path.realpath(local_path)
            if os.path.isabs(local_path)
            else os.path.realpath(os.path.join(base_dir, local_path))
        )
        if not os.path.isfile(resolved_path):
            raise FileNotFoundError(
                "Image for tree tip '{}' in '--tip-image-manifest' was not found: {}".format(
                    leaf_name,
                    resolved_path,
                )
            )
        path_by_leaf[leaf_name] = resolved_path

    if missing_path_leaf_names:
        message = "Tree tips have missing local_path values in --tip-image-manifest: {}".format(
            " ".join(sorted(missing_path_leaf_names))
        )
        if unmatched == "error":
            raise ValueError(message)
        if unmatched == "warn":
            sys.stderr.write(message + "\n")
    return path_by_leaf


def _load_tip_image(path, max_edge_px=None):
    from nwkit.image import (
        load_pillow_modules,
        rasterize_svg_to_image,
        validate_image_dimensions,
    )

    Image, _, ImageOps = load_pillow_modules()
    try:
        if os.path.splitext(path)[1].lower() == ".svg":
            image = rasterize_svg_to_image(path, max_edge=max_edge_px)
        else:
            with Image.open(path) as source_image:
                width, height = source_image.size
                validate_image_dimensions(width, height, label="Tip image")
                if max_edge_px not in (None, 0):
                    source_image.draft(
                        source_image.mode,
                        (int(max_edge_px), int(max_edge_px)),
                    )
                if max_edge_px not in (None, 0):
                    resampling = getattr(Image, "Resampling", Image)
                    source_image.thumbnail(
                        (int(max_edge_px), int(max_edge_px)),
                        resampling.LANCZOS,
                    )
                image = ImageOps.exif_transpose(source_image)
                image = image.convert("RGBA")
                image.load()
    except (OSError, ValueError) as exc:
        raise ValueError("Failed to read tip image '{}': {}".format(path, exc)) from exc
    return np.asarray(image)


def _load_tip_images(path_by_leaf, image_size_pt=18.0):
    image_size_pt = _validated_tip_image_size(image_size_pt)
    max_edge_px = max(
        1,
        min(
            MAX_TIP_IMAGE_EDGE_PX,
            int(math.ceil(image_size_pt * 4.0)),
        ),
    )
    image_by_leaf = dict()
    image_by_path: dict[str, np.ndarray] = {}
    for leaf_name, path in path_by_leaf.items():
        image = image_by_path.get(path)
        if image is None:
            image = _load_tip_image(path, max_edge_px=max_edge_px)
            image_by_path[path] = image
        image_by_leaf[leaf_name] = image
    return image_by_leaf


def _draw_tip_images(ax, leaf_order, ycoord, image_by_leaf, image_size_pt):
    for leaf in leaf_order:
        image = image_by_leaf.get(str(leaf.name))
        if image is None:
            continue
        height_px, width_px = image.shape[:2]
        max_edge_px = max(height_px, width_px)
        if max_edge_px <= 0:
            continue
        image_artist = OffsetImage(
            image,
            zoom=float(image_size_pt) / float(max_edge_px),
            interpolation="hanning",
        )
        annotation = AnnotationBbox(
            image_artist,
            (0.5, ycoord[leaf]),
            xycoords=("axes fraction", "data"),
            box_alignment=(0.5, 0.5),
            frameon=False,
            pad=0.0,
            annotation_clip=False,
        )
        annotation.set_clip_on(False)
        ax.add_artist(annotation)


def _assign_group_colors(group_names, palette="tab10"):
    if len(group_names) == 0:
        return dict()
    cmap = plt.get_cmap(palette)
    colors = dict()
    for index, group_name in enumerate(sorted(group_names)):
        colors[group_name] = mcolors.to_hex(
            cmap(index % max(getattr(cmap, "N", len(group_names)), 1))
        )
    return colors


def _parse_property_colors(raw_colors):
    colors = dict()
    for raw in raw_colors or []:
        if "=" not in str(raw):
            raise ValueError(
                "--property-color must use VALUE=COLOR syntax: {}".format(raw)
            )
        value, color = str(raw).split("=", 1)
        value = value.strip()
        color = color.strip()
        if value == "" or color == "":
            raise ValueError(
                "--property-color must use VALUE=COLOR syntax: {}".format(raw)
            )
        try:
            mcolors.to_rgba(color)
        except ValueError as exc:
            raise ValueError(
                "Invalid color for --property-color '{}': {}".format(value, color)
            ) from exc
        colors[value] = color
    return colors


def _parse_property_names(raw):
    if raw in (None, ""):
        return list()
    names = [name.strip() for name in str(raw).split(",") if name.strip()]
    if len(names) != len(set(names)):
        raise ValueError("Property names must not be repeated: {}".format(raw))
    return names


def _node_matches_target(node, raw_target):
    targets = {
        value.strip().lower()
        for value in str(raw_target or "all").split(",")
        if value.strip()
    }
    allowed = {"all", "root", "intnode", "leaf"}
    unknown = targets - allowed
    if unknown:
        raise ValueError(
            "Unsupported node target(s): {}".format(", ".join(sorted(unknown)))
        )
    if "all" in targets:
        return True
    if node.is_root:
        return "root" in targets
    if node.is_leaf:
        return "leaf" in targets
    return "intnode" in targets


def _coerce_filter_value(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return str(value)


def _matches_property_filters(node, raw_filters, option_name="--node-label-filter"):
    operators = {
        "ge": lambda left, right: left >= right,
        "gt": lambda left, right: left > right,
        "le": lambda left, right: left <= right,
        "lt": lambda left, right: left < right,
        "eq": lambda left, right: left == right,
        "ne": lambda left, right: left != right,
    }
    parsed_filters = []
    for raw in raw_filters or []:
        parts = str(raw).split(":", 2)
        if len(parts) != 3 or parts[1].lower() not in operators:
            raise ValueError(
                "{} must use PROPERTY:OP:VALUE syntax with "
                "OP in ge,gt,le,lt,eq,ne: {}".format(option_name, raw)
            )
        prop, operator, expected = parts[0].strip(), parts[1].lower(), parts[2].strip()
        if not prop:
            raise ValueError(
                "{} requires a non-empty PROPERTY: {}".format(option_name, raw)
            )
        parsed_filters.append((prop, operator, expected))
    for prop, operator, expected in parsed_filters:
        if prop not in node.props:
            return False
        observed = _coerce_filter_value(node.props[prop])
        expected = _coerce_filter_value(expected)
        if isinstance(observed, float) != isinstance(expected, float):
            observed = str(node.props[prop])
            expected = str(parts[2].strip())
        if not operators[operator](observed, expected):
            return False
    return True


def _property_color(prop, value, property_colors, fallback_index=0, palette="tab10"):
    candidates = ["{}:{}".format(prop, value), str(prop), str(value)]
    if str(prop).startswith("asr_p_"):
        candidates.append(str(prop)[len("asr_p_") :])
    for candidate in candidates:
        if candidate in property_colors:
            return property_colors[candidate]
    cmap = plt.get_cmap(palette)
    return mcolors.to_hex(cmap(fallback_index % max(getattr(cmap, "N", 1), 1)))


def _draw_probability_pie(ax, x, y, probabilities, colors, marker_size_pt):
    try:
        values = [float(probability) for probability in probabilities]
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "Node-pie properties must be finite, non-negative numbers."
        ) from exc
    if any((not math.isfinite(value)) or value < 0.0 for value in values):
        raise ValueError("Node-pie properties must be finite, non-negative numbers.")
    total = sum(values)
    if total <= 0.0:
        return None
    size = float(marker_size_pt)
    drawing = DrawingArea(size, size, 0.0, 0.0)
    center = size / 2.0
    radius = max((size / 2.0) - 0.25, 0.1)
    normalized = [probability / total for probability in values]
    cumulative = 0.0
    for probability, color in zip(normalized, colors, strict=True):
        if probability <= 0.0:
            continue
        start = cumulative
        cumulative += probability
        drawing.add_artist(
            Wedge(
                (center, center),
                radius,
                90.0 - (360.0 * cumulative),
                90.0 - (360.0 * start),
                facecolor=color,
                edgecolor="none",
            )
        )
    drawing.add_artist(
        Circle(
            (center, center),
            radius,
            facecolor="none",
            edgecolor=LEGEND_EDGE_COLOR,
            linewidth=0.45,
        )
    )
    artist = AnnotationBbox(
        drawing,
        (x, y),
        xybox=(0.0, 0.0),
        xycoords="data",
        boxcoords="offset points",
        frameon=False,
        box_alignment=(0.5, 0.5),
        annotation_clip=False,
        zorder=7,
    )
    ax.add_artist(artist)
    return artist


def _format_property_value(value, decimals=2):
    try:
        return ("{:.%df}" % int(decimals)).format(float(value))
    except (TypeError, ValueError):
        return str(value)


def _isolated_matplotlib_state(function):
    def wrapped(*args, **kwargs):
        existing_figures = set(plt.get_fignums())
        try:
            with matplotlib.rc_context(
                {
                    "font.family": [FONT_FAMILY],
                    "font.sans-serif": [FONT_FAMILY],
                    "svg.fonttype": "none",
                }
            ):
                return function(*args, **kwargs)
        finally:
            for figure_number in set(plt.get_fignums()) - existing_figures:
                plt.close(figure_number)

    return wrapped


def _validated_tip_image_size(value):
    try:
        image_size_pt = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("--tip-image-size must be a finite number.") from exc
    if not math.isfinite(image_size_pt):
        raise ValueError("--tip-image-size must be a finite number.")
    if image_size_pt <= 0.0:
        raise ValueError("--tip-image-size must be greater than zero.")
    if image_size_pt > MAX_TIP_IMAGE_SIZE_PT:
        raise ValueError(
            "--tip-image-size must be no greater than {} points.".format(
                int(MAX_TIP_IMAGE_SIZE_PT)
            )
        )
    return image_size_pt


def _readable_text_angle(angle):
    """Return a left-to-right equivalent of an edge direction in degrees."""

    angle = ((float(angle) + 180.0) % 360.0) - 180.0
    if angle > 90.0:
        angle -= 180.0
    elif angle < -90.0:
        angle += 180.0
    return angle


def _spatial_text_alignment(angle):
    """Choose label alignment while retaining the outward offset direction."""

    normalized = ((float(angle) + 180.0) % 360.0) - 180.0
    points_right = -90.0 <= normalized <= 90.0
    return _readable_text_angle(normalized), ("left" if points_right else "right")


def _first_nonzero_direction(points, from_end=False):
    """Return a unit direction along the first or last nonzero path segment."""

    pairs = list(zip(points, points[1:], strict=False))
    if from_end:
        pairs = list(reversed(pairs))
    for start, end in pairs:
        dx = float(end[0]) - float(start[0])
        dy = float(end[1]) - float(start[1])
        norm = math.hypot(dx, dy)
        if norm > 1e-12:
            return (dx / norm, dy / norm)
    return (1.0, 0.0)


def _spatial_node_label_placement(node, drawing_layout, clearance_points):
    """Place a node label in locally open space beside its incident edges."""

    if node.is_root:
        if drawing_layout.root_path:
            # Root paths run from the empty-side stub into the root. Reverse
            # that direction so labels extend into the reserved empty sector.
            inward = _first_nonzero_direction(
                drawing_layout.root_path,
                from_end=True,
            )
            direction = (-inward[0], -inward[1])
        else:
            child_angles = sorted(
                math.atan2(direction[1], direction[0]) % (2.0 * math.pi)
                for direction in (
                    _first_nonzero_direction(
                        drawing_layout.edge_paths.get(child, ()),
                    )
                    for child in node.get_children()
                )
            )
            if child_angles:
                cyclic = child_angles + [child_angles[0] + (2.0 * math.pi)]
                gap_start, gap_end = max(
                    zip(cyclic, cyclic[1:], strict=False),
                    key=lambda gap: (gap[1] - gap[0], -gap[0]),
                )
                empty_angle = (gap_start + gap_end) / 2.0
                direction = (math.cos(empty_angle), math.sin(empty_angle))
            else:
                direction = (-1.0, 0.0)
    else:
        incoming = _first_nonzero_direction(
            drawing_layout.edge_paths.get(node, ()),
            from_end=True,
        )
        incident_angles = [math.atan2(-incoming[1], -incoming[0])]
        incident_angles.extend(
            math.atan2(direction[1], direction[0])
            for direction in (
                _first_nonzero_direction(
                    drawing_layout.edge_paths.get(child, ()),
                )
                for child in node.get_children()
            )
        )
        incident_angles = sorted(angle % (2.0 * math.pi) for angle in incident_angles)
        cyclic = incident_angles + [incident_angles[0] + (2.0 * math.pi)]
        gap_start, gap_end = max(
            zip(cyclic, cyclic[1:], strict=False),
            key=lambda gap: (gap[1] - gap[0], -gap[0]),
        )
        empty_angle = (gap_start + gap_end) / 2.0
        direction = (math.cos(empty_angle), math.sin(empty_angle))

    if direction[0] > 0.2:
        horizontal_alignment = "left"
    elif direction[0] < -0.2:
        horizontal_alignment = "right"
    else:
        horizontal_alignment = "center"
    if direction[1] > 0.4:
        vertical_alignment = "bottom"
    elif direction[1] < -0.4:
        vertical_alignment = "top"
    else:
        vertical_alignment = "center"
    return {
        "offset": (
            direction[0] * clearance_points,
            direction[1] * clearance_points,
        ),
        "rotation": 0.0,
        "horizontal_alignment": horizontal_alignment,
        "vertical_alignment": vertical_alignment,
        # Horizontal node labels separate most efficiently along their short
        # (vertical) axis. Preserve the locally empty half-plane for branch
        # avoidance while still allowing artist--artist resolution to reverse.
        "shift_direction": (0.0, 1.0 if direction[1] >= 0.0 else -1.0),
    }


def _parse_tip_label_wrap(value):
    text = str(value).strip().lower()
    if text in {"none", "auto", "taxonomy"}:
        return text
    try:
        width = int(text)
    except ValueError as error:
        raise ValueError(
            "'--tip-label-wrap' must be none, auto, taxonomy, or a positive integer."
        ) from error
    if width <= 0:
        raise ValueError(
            "'--tip-label-wrap' must be none, auto, taxonomy, or a positive integer."
        )
    return width


def _wrap_tip_label(text, width):
    """Insert display-only line breaks while preserving visible delimiters."""

    text = str(text)
    width = int(width)
    if width <= 0 or len(text) <= width:
        return text
    lines = []
    remaining = text
    while len(remaining) > width:
        whitespace_breaks = [
            index
            for index, char in enumerate(remaining[: width + 1])
            if char.isspace() and index > 0
        ]
        delimiter_breaks = [
            index + 1
            for index, char in enumerate(remaining[:width])
            if char in "_-/|" and index > 0
        ]
        if whitespace_breaks:
            cut = whitespace_breaks[-1]
            line = remaining[:cut].rstrip()
            remaining = remaining[cut:].lstrip()
        elif delimiter_breaks:
            cut = delimiter_breaks[-1]
            line = remaining[:cut]
            remaining = remaining[cut:]
        else:
            cut = width
            line = remaining[:cut]
            remaining = remaining[cut:]
        lines.append(line)
    lines.append(remaining)
    return "\n".join(lines)


def _measure_texts_in_inches(texts, font_size, font_family):
    """Measure multiline strings with the same Matplotlib text renderer."""

    figure = Figure(figsize=(1.0, 1.0), dpi=72.0)
    canvas = FigureCanvasAgg(figure)
    renderer = canvas.get_renderer()
    measured = {}
    for text in dict.fromkeys(texts):
        artist = Text(
            x=0.0,
            y=0.0,
            text=text,
            fontsize=float(font_size),
            fontfamily=font_family,
            linespacing=1.15,
        )
        artist.set_figure(figure)
        bounds = artist.get_window_extent(renderer=renderer)
        measured[text] = (
            max(float(bounds.width) / 72.0, 0.0),
            max(float(bounds.height) / 72.0, float(font_size) / 72.0),
        )
    return measured


def _auto_wrap_candidates(text):
    length = len(str(text))
    if length <= 1:
        return [length]
    widths = {length}
    for line_count in range(2, min(5, length) + 1):
        widths.add(max(1, int(math.ceil(length / line_count))))
    return sorted(widths, reverse=True)


def _taxonomic_prefix(text):
    """Return a conservative underscore-delimited binomial prefix."""

    parts = str(text).split("_")
    if len(parts) < 2:
        return None
    genus, species = parts[:2]
    if not genus or not species:
        return None
    if not genus[0].isupper() or not species[0].islower():
        return None
    return "{}_{}".format(genus, species)


def _wrap_taxonomic_label(text, width):
    prefix = _taxonomic_prefix(text)
    if prefix is None:
        return _wrap_tip_label(text, width)
    suffix = str(text)[len(prefix) :].lstrip("_")
    if not suffix:
        return str(text)
    wrapped_suffix = _wrap_tip_label(suffix, width)
    return "{}\n{}".format(prefix, wrapped_suffix)


def _prepare_tip_label_text(
    leaf_order,
    wrap,
    font_size,
    font_family,
    layout_name,
    panel_width_in,
    panel_height_in,
):
    """Resolve wrapping and return displayed labels with exact dimensions."""

    wrap = _parse_tip_label_wrap(wrap)
    raw_text_by_leaf = {leaf: str(leaf.name or "") for leaf in leaf_order}
    if wrap == "none":
        text_by_leaf = dict(raw_text_by_leaf)
    elif isinstance(wrap, int):
        text_by_leaf = {
            leaf: _wrap_tip_label(text, wrap) for leaf, text in raw_text_by_leaf.items()
        }
    else:
        candidates_by_leaf = {
            leaf: [
                _wrap_taxonomic_label(text, width)
                for width in _auto_wrap_candidates(text)
            ]
            for leaf, text in raw_text_by_leaf.items()
        }
        all_candidates = [
            candidate
            for candidates in candidates_by_leaf.values()
            for candidate in candidates
        ]
        candidate_sizes = _measure_texts_in_inches(
            all_candidates,
            font_size=font_size,
            font_family=font_family,
        )
        spatial = layout_name in {
            "circular",
            "radial",
            "unrooted",
            "spiral",
            "fractal",
        }
        normal_budget = max(min(panel_width_in, panel_height_in) * 0.20, 0.55)
        if spatial:
            perimeter_per_tip = (2.0 * (panel_width_in + panel_height_in)) / max(
                len(leaf_order), 1
            )
            tangential_budget = max(perimeter_per_tip * 0.72, float(font_size) / 72.0)
        else:
            tangential_budget = max(float(font_size) * 1.25 / 72.0, 1e-6)
        text_by_leaf = {}
        for leaf, candidates in candidates_by_leaf.items():

            def score(candidate):
                width, height = candidate_sizes[candidate]
                congestion = max(
                    width / normal_budget,
                    height / tangential_budget,
                )
                line_penalty = 0.22 * candidate.count("\n")
                return congestion + line_penalty, width * height, candidate.count("\n")

            unwrapped = raw_text_by_leaf[leaf]
            if unwrapped not in candidate_sizes:
                candidate_sizes.update(
                    _measure_texts_in_inches(
                        [unwrapped],
                        font_size=font_size,
                        font_family=font_family,
                    )
                )
                candidates = [unwrapped] + candidates
            if wrap == "taxonomy" and _taxonomic_prefix(unwrapped) is not None:
                semantic_candidates = [
                    candidate
                    for candidate in candidates
                    if "\n" in candidate or candidate == unwrapped
                ]
                candidates = semantic_candidates or candidates
            unwrapped_width, unwrapped_height = candidate_sizes[unwrapped]
            if (
                unwrapped_width <= normal_budget
                and unwrapped_height <= tangential_budget
            ):
                text_by_leaf[leaf] = unwrapped
                continue
            best = min(candidates, key=score)
            if score(best)[0] < score(unwrapped)[0] * 0.95:
                text_by_leaf[leaf] = best
            else:
                text_by_leaf[leaf] = unwrapped
    sizes = _measure_texts_in_inches(
        text_by_leaf.values(),
        font_size=font_size,
        font_family=font_family,
    )
    return (
        text_by_leaf,
        {leaf: sizes[text] for leaf, text in text_by_leaf.items()},
    )


def _fit_artists_within_figure(figure, axes, artists, padding_points=2.0):
    """Shrink an axes until rendered artists fit the fixed-size figure."""

    if not artists:
        return {
            "fits_within_figure": True,
            "overflow_left_points": 0.0,
            "overflow_right_points": 0.0,
            "overflow_bottom_points": 0.0,
            "overflow_top_points": 0.0,
            "maximum_overflow_points": 0.0,
        }
    padding_pixels = float(padding_points) * figure.dpi / 72.0
    for _ in range(4):
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        bounds = [artist.get_window_extent(renderer=renderer) for artist in artists]
        left = min(bound.x0 for bound in bounds)
        right = max(bound.x1 for bound in bounds)
        bottom = min(bound.y0 for bound in bounds)
        top = max(bound.y1 for bound in bounds)
        figure_bounds = figure.bbox
        overflow_left = max((figure_bounds.x0 + padding_pixels) - left, 0.0)
        overflow_right = max(right - (figure_bounds.x1 - padding_pixels), 0.0)
        overflow_bottom = max((figure_bounds.y0 + padding_pixels) - bottom, 0.0)
        overflow_top = max(top - (figure_bounds.y1 - padding_pixels), 0.0)
        if max(overflow_left, overflow_right, overflow_bottom, overflow_top) < 0.25:
            break
        position = axes.get_position(original=True)
        new_left = position.x0 + (overflow_left / figure_bounds.width)
        new_right = position.x1 - (overflow_right / figure_bounds.width)
        new_bottom = position.y0 + (overflow_bottom / figure_bounds.height)
        new_top = position.y1 - (overflow_top / figure_bounds.height)
        if new_right - new_left <= 0.1 or new_top - new_bottom <= 0.1:
            break
        axes.set_position(
            [new_left, new_bottom, new_right - new_left, new_top - new_bottom]
        )
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    bounds = [artist.get_window_extent(renderer=renderer) for artist in artists]
    figure_bounds = figure.bbox
    overflow = {
        "overflow_left_points": max(
            (figure_bounds.x0 + padding_pixels) - min(bound.x0 for bound in bounds),
            0.0,
        )
        * 72.0
        / figure.dpi,
        "overflow_right_points": max(
            max(bound.x1 for bound in bounds) - (figure_bounds.x1 - padding_pixels),
            0.0,
        )
        * 72.0
        / figure.dpi,
        "overflow_bottom_points": max(
            (figure_bounds.y0 + padding_pixels) - min(bound.y0 for bound in bounds),
            0.0,
        )
        * 72.0
        / figure.dpi,
        "overflow_top_points": max(
            max(bound.y1 for bound in bounds) - (figure_bounds.y1 - padding_pixels),
            0.0,
        )
        * 72.0
        / figure.dpi,
    }
    overflow["maximum_overflow_points"] = max(overflow.values(), default=0.0)
    overflow["fits_within_figure"] = bool(overflow["maximum_overflow_points"] < 0.25)
    return overflow


def _parse_scale_bar(value, tree_span):
    text = str(value).strip().lower()
    if text in {"none", ""}:
        return None
    if text == "auto":
        target = max(float(tree_span) / 5.0, 1e-12)
        exponent = math.floor(math.log10(target))
        scaled = target / (10.0**exponent)
        base = 1.0 if scaled < 1.5 else (2.0 if scaled < 3.5 else 5.0)
        return base * (10.0**exponent)
    try:
        result = float(text)
    except ValueError as error:
        raise ValueError(
            "'--scale-bar' must be none, auto, or a positive number."
        ) from error
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError("'--scale-bar' must be none, auto, or a positive number.")
    return result


def _parse_depth_guide(value, tree_span):
    text = str(value).strip().lower()
    if text in {"none", ""}:
        return None
    if text == "auto":
        target = max(float(tree_span) / 4.0, 1e-12)
        exponent = math.floor(math.log10(target))
        scaled = target / (10.0**exponent)
        base = 1.0 if scaled < 1.5 else (2.0 if scaled < 3.5 else 5.0)
        return base * (10.0**exponent)
    try:
        interval = float(text)
    except ValueError as error:
        raise ValueError(
            "'--depth-guide' must be none, auto, or a positive interval."
        ) from error
    if not math.isfinite(interval) or interval <= 0.0:
        raise ValueError("'--depth-guide' must be none, auto, or a positive interval.")
    if interval > float(tree_span) + 1e-12:
        raise ValueError(
            "'--depth-guide' interval must not exceed the displayed tree-depth "
            "span ({:g}).".format(tree_span)
        )
    return interval


def _depth_guide_ticks(tree_span, interval):
    if interval is None:
        return []
    count = int(math.floor((float(tree_span) + 1e-12) / float(interval)))
    if count > 50:
        raise ValueError(
            "'--depth-guide' would draw {} intervals; use an interval of at "
            "least {:g}.".format(count, float(tree_span) / 50.0)
        )
    return [float(interval) * index for index in range(count + 1)]


def _distance_label(prefix, unit):
    unit = str(unit).strip()
    return "{} ({})".format(prefix, unit) if unit else prefix


def _add_scale_bar(
    figure,
    axes,
    size,
    label,
    color,
    font_family,
    font_size,
    anchor_x,
    anchor_y,
):
    """Add a scale bar to the dedicated strip below the tree axes."""

    artist = AnchoredSizeBar(
        axes.transData,
        size,
        label,
        loc="lower left",
        bbox_to_anchor=(float(anchor_x), float(anchor_y)),
        bbox_transform=figure.transFigure,
        pad=0.0,
        borderpad=0.0,
        sep=2.0,
        frameon=False,
        size_vertical=0.0,
        color=color,
        label_top=True,
        fontproperties=FontProperties(
            family=font_family,
            size=font_size,
        ),
    )
    artist.set_zorder(10)
    axes.add_artist(artist)
    return artist


def _add_slanted_depth_guide(
    axes,
    root_x,
    ticks,
    color,
    font_family,
    font_size,
):
    """Draw a projected-depth axis and vertical guides for a slanted tree."""

    artists = []
    transform = axes.get_xaxis_transform()
    maximum = max(ticks, default=0.0)
    axes.plot(
        [root_x, root_x + maximum],
        [0.0, 0.0],
        transform=transform,
        color=color,
        linewidth=0.55,
        clip_on=False,
        zorder=0.4,
    )
    for depth in ticks:
        x = root_x + depth
        axes.axvline(
            x,
            color=color,
            linewidth=0.45,
            linestyle=(0, (1.5, 2.5)),
            alpha=0.55,
            zorder=0.25,
        )
        axes.plot(
            [x, x],
            [0.0, -0.014],
            transform=transform,
            color=color,
            linewidth=0.55,
            clip_on=False,
            zorder=0.4,
        )
        label = axes.annotate(
            "{:g}".format(depth),
            xy=(x, 0.0),
            xycoords=("data", "axes fraction"),
            xytext=(0.0, -3.0),
            textcoords="offset points",
            ha="center",
            va="top",
            fontsize=font_size,
            fontfamily=font_family,
            color=LABEL_COLOR,
            annotation_clip=False,
            zorder=9,
        )
        artists.append(label)
    return artists


def _angular_sector_bounds(radius, angular_span, angular_center):
    """Return Cartesian bounds for an origin-anchored circular sector."""

    radius = max(float(radius), 0.0)
    span = min(max(float(angular_span), 0.0), 360.0)
    center = float(angular_center)
    if span >= 360.0:
        return (-radius, radius, -radius, radius)
    start = math.radians(center - (span / 2.0))
    end = math.radians(center + (span / 2.0))
    angles = [start, end]
    first_quadrant = math.ceil(start / (math.pi / 2.0))
    last_quadrant = math.floor(end / (math.pi / 2.0))
    angles.extend(
        index * (math.pi / 2.0) for index in range(first_quadrant, last_quadrant + 1)
    )
    points = [(0.0, 0.0)] + [
        (radius * math.cos(angle), radius * math.sin(angle)) for angle in angles
    ]
    x_values, y_values = zip(*points, strict=True)
    return min(x_values), max(x_values), min(y_values), max(y_values)


def _sector_contains_direction(angular_span, angular_center, direction):
    """Return whether a sector contains a cardinal direction in degrees."""

    span = float(angular_span)
    if span >= 360.0:
        return True
    delta = ((float(direction) - float(angular_center) + 180.0) % 360.0) - 180.0
    return abs(delta) <= (span / 2.0) + 1e-12


def _polar_auto_figure_height(
    panel_width,
    angular_span,
    angular_center,
    radial_label_extent,
    tangential_label_extent,
    top_margin,
    bottom_margin,
):
    """Estimate a sector-aware height including radially oriented labels."""

    x_min, x_max, y_min, y_max = _angular_sector_bounds(
        1.0,
        angular_span,
        angular_center,
    )
    x_span = max(x_max - x_min, 1e-6)
    y_span = max(y_max - y_min, 1e-6)
    radial = max(float(radial_label_extent), 0.0)
    tangential = max(float(tangential_label_extent), 0.0)
    horizontal_padding = 2.0 * tangential
    vertical_padding = 2.0 * tangential
    if _sector_contains_direction(angular_span, angular_center, 0.0):
        horizontal_padding += radial
    if _sector_contains_direction(angular_span, angular_center, 180.0):
        horizontal_padding += radial
    if _sector_contains_direction(angular_span, angular_center, 90.0):
        vertical_padding += radial
    if _sector_contains_direction(angular_span, angular_center, 270.0):
        vertical_padding += radial
    available_width = max(float(panel_width) - horizontal_padding, 0.2)
    geometry_scale = available_width / x_span
    panel_height = (geometry_scale * y_span) + vertical_padding
    return max(
        panel_height + float(top_margin) + float(bottom_margin),
        float(top_margin) + float(bottom_margin) + 0.2,
    )


def _add_radial_depth_guide(
    axes,
    drawing_layout,
    tree,
    ticks,
    color,
    font_family,
    font_size,
):
    """Draw root-depth arcs and label them within the occupied sector."""

    root_x = drawing_layout.xcoord[tree]
    root_y = drawing_layout.ycoord[tree]
    leaf_angles = sorted(
        math.radians(drawing_layout.label_angles.get(leaf, 0.0)) % (2.0 * math.pi)
        for leaf in tree.leaves()
    )
    segments = [
        (start, end)
        for path in drawing_layout.edge_paths.values()
        for start, end in zip(path, path[1:], strict=False)
    ]

    def point_segment_distance(point, start, end):
        delta_x = end[0] - start[0]
        delta_y = end[1] - start[1]
        denominator = (delta_x * delta_x) + (delta_y * delta_y)
        if denominator <= 1e-18:
            return math.hypot(point[0] - start[0], point[1] - start[1])
        fraction = (
            ((point[0] - start[0]) * delta_x) + ((point[1] - start[1]) * delta_y)
        ) / denominator
        fraction = min(max(fraction, 0.0), 1.0)
        closest = (
            start[0] + (fraction * delta_x),
            start[1] + (fraction * delta_y),
        )
        return math.hypot(point[0] - closest[0], point[1] - closest[1])

    positive_ticks = [depth for depth in ticks if depth > 1e-12]
    outer_radius = max(positive_ticks, default=1.0)
    if len(positive_ticks) <= 12:
        score_ticks = positive_ticks
    else:
        score_ticks = [
            positive_ticks[round(index * (len(positive_ticks) - 1) / 11.0)]
            for index in range(12)
        ]
    if len(segments) <= 512:
        score_segments = segments
    else:
        score_segments = [
            segments[round(index * (len(segments) - 1) / 511.0)] for index in range(512)
        ]

    def nearest_leaf_angle_distance(direction):
        if not leaf_angles:
            return math.pi
        position = bisect.bisect_left(leaf_angles, direction)
        neighbors = (
            leaf_angles[position % len(leaf_angles)],
            leaf_angles[(position - 1) % len(leaf_angles)],
        )
        return min(
            abs(
                math.atan2(
                    math.sin(direction - angle),
                    math.cos(direction - angle),
                )
            )
            for angle in neighbors
        )

    def direction_score(direction):
        cosine = math.cos(direction)
        sine = math.sin(direction)
        branch_clearance = min(
            (
                point_segment_distance(
                    (root_x + (depth * cosine), root_y + (depth * sine)),
                    start,
                    end,
                )
                for depth in score_ticks
                for start, end in score_segments
            ),
            default=outer_radius,
        )
        leaf_clearance = nearest_leaf_angle_distance(direction) * outer_radius
        return min(branch_clearance, leaf_clearance)

    angular_span = float(drawing_layout.metadata.get("angular_span_degrees", 360.0))
    angular_center = float(drawing_layout.metadata.get("angular_center_degrees", 90.0))
    available = math.radians(angular_span)
    start = math.radians(angular_center) - (available / 2.0)
    candidate_count = 360 if len(segments) <= 512 else 180
    if angular_span < 360.0:
        candidates = [
            start + (available * index / max(candidate_count - 1, 1))
            for index in range(candidate_count)
        ]
    else:
        candidates = [
            2.0 * math.pi * index / candidate_count for index in range(candidate_count)
        ]
    direction = max(candidates, key=direction_score)
    best_score = direction_score(direction)
    labels_have_clear_sector = (
        len(segments) <= 2000 and best_score >= outer_radius * 0.012
    )
    artists = []
    for depth in positive_ticks:
        guide_style = {
            "facecolor": "none",
            "edgecolor": color,
            "linewidth": 0.45,
            "linestyle": (0, (1.5, 2.5)),
            "alpha": 0.60,
            "zorder": 0.25,
        }
        if angular_span < 360.0:
            axes.add_patch(
                Arc(
                    (root_x, root_y),
                    2.0 * depth,
                    2.0 * depth,
                    theta1=math.degrees(start),
                    theta2=math.degrees(start + available),
                    **guide_style,
                )
            )
        else:
            axes.add_patch(
                Circle(
                    (root_x, root_y),
                    depth,
                    **guide_style,
                )
            )
        if not labels_have_clear_sector:
            continue
        x = root_x + (depth * math.cos(direction))
        y = root_y + (depth * math.sin(direction))
        label = axes.annotate(
            "{:g}".format(depth),
            xy=(x, y),
            xytext=(3.0 * -math.sin(direction), 3.0 * math.cos(direction)),
            textcoords="offset points",
            ha="left" if math.cos(direction) >= 0.0 else "right",
            va="center",
            fontsize=font_size,
            fontfamily=font_family,
            color=LABEL_COLOR,
            bbox={
                "boxstyle": "square,pad=0.08",
                "facecolor": "#ffffff",
                "edgecolor": "none",
                "alpha": 0.88,
            },
            annotation_clip=False,
            zorder=9,
        )
        artists.append(label)
    return artists, labels_have_clear_sector


def _add_bottom_guide_title(
    figure,
    text,
    x,
    y,
    font_family,
    font_size,
):
    return figure.text(
        float(x),
        float(y),
        text,
        ha="left",
        va="bottom",
        fontsize=font_size,
        fontfamily=font_family,
        color=LABEL_COLOR,
        zorder=10,
    )


def _add_spiral_depth_key(
    figure,
    ticks,
    tree_span,
    axes_left,
    axes_right,
    strip_height_points,
    font_family,
    font_size,
    color,
):
    """Add a linear key for depth encoded across the spiral track."""

    if not ticks or tree_span <= 0.0:
        return []
    figure_height_points = figure.get_figheight() * 72.0
    figure_width_points = figure.get_figwidth() * 72.0
    start_x = float(axes_left) + (
        max(float(font_size) * 0.45, 3.0) / figure_width_points
    )
    end_x = start_x + min((float(axes_right) - start_x) * 0.46, 0.42)
    line_y = (
        float(strip_height_points) - float(font_size) - 5.0
    ) / figure_height_points
    tick_height = 3.0 / figure_height_points
    line = Line2D(
        [start_x, end_x],
        [line_y, line_y],
        transform=figure.transFigure,
        color=color,
        linewidth=0.65,
        zorder=10,
        clip_on=False,
    )
    figure.add_artist(line)
    artists = []
    for depth in ticks:
        fraction = depth / float(tree_span)
        x = start_x + ((end_x - start_x) * fraction)
        tick = Line2D(
            [x, x],
            [line_y - tick_height, line_y + tick_height],
            transform=figure.transFigure,
            color=color,
            linewidth=0.55,
            zorder=10,
            clip_on=False,
        )
        figure.add_artist(tick)
        label = figure.text(
            x,
            line_y + tick_height + (1.0 / figure_height_points),
            "{:g}".format(depth),
            ha="center",
            va="bottom",
            fontsize=font_size,
            fontfamily=font_family,
            color=LABEL_COLOR,
            zorder=10,
        )
        artists.append(label)
    return artists


def _parse_legend_columns(value, item_count):
    text = str(value).strip().lower()
    if text == "auto":
        return max(1, min(4, int(math.ceil(max(item_count, 1) / 3.0))))
    try:
        columns = int(text)
    except ValueError as error:
        raise ValueError(
            "'--legend-columns' must be auto or a positive integer."
        ) from error
    if columns <= 0:
        raise ValueError("'--legend-columns' must be auto or a positive integer.")
    return columns


def _tip_font_style(text, mode):
    mode = str(mode).strip().lower()
    if mode == "plain":
        return "normal"
    if mode == "italic":
        return "italic"
    if mode == "taxonomy":
        prefix = _taxonomic_prefix(str(text).replace("\n", ""))
        return "italic" if prefix == str(text).replace("\n", "") else "normal"
    raise ValueError("'--tip-label-font-style' must be plain, italic, or taxonomy.")


def _parse_branch_width_range(value):
    parts = [part.strip() for part in str(value).split(",")]
    if len(parts) != 2:
        raise ValueError("'--branch-width-range' must contain MIN,MAX.")
    try:
        minimum, maximum = (float(part) for part in parts)
    except ValueError as error:
        raise ValueError("'--branch-width-range' must contain MIN,MAX.") from error
    if (
        minimum <= 0.0
        or maximum < minimum
        or not all(map(math.isfinite, (minimum, maximum)))
    ):
        raise ValueError("'--branch-width-range' must satisfy 0 < MIN <= MAX.")
    return minimum, maximum


def _branch_style_maps(
    tree,
    base_color,
    base_width,
    color_property,
    width_property,
    width_range,
    property_colors,
    palette,
):
    width = float(base_width)
    if not math.isfinite(width) or width <= 0.0:
        raise ValueError("'--branch-width' must be greater than zero.")
    color_by_node = {}
    width_by_node = {}
    nodes = [node for node in tree.traverse() if not node.is_root]
    if color_property not in (None, ""):
        values = sorted(
            {
                str(node.props[color_property])
                for node in nodes
                if node.props.get(color_property) not in (None, "")
            }
        )
        value_index = {value: index for index, value in enumerate(values)}
        for node in nodes:
            value = node.props.get(color_property)
            if value not in (None, ""):
                color_by_node[node] = _property_color(
                    prop=color_property,
                    value=value,
                    property_colors=property_colors,
                    fallback_index=value_index[str(value)],
                    palette=palette,
                )
    if width_property not in (None, ""):
        numeric = {}
        for node in nodes:
            try:
                value = float(node.props.get(width_property))
            except (TypeError, ValueError):
                continue
            if math.isfinite(value):
                numeric[node] = value
        if numeric:
            output_min, output_max = _parse_branch_width_range(width_range)
            input_min = min(numeric.values())
            input_max = max(numeric.values())
            for node, value in numeric.items():
                fraction = (
                    0.5
                    if input_max <= input_min
                    else (value - input_min) / (input_max - input_min)
                )
                width_by_node[node] = output_min + (
                    (output_max - output_min) * fraction
                )
    return (
        {node: color_by_node.get(node, base_color) for node in nodes},
        {node: width_by_node.get(node, width) for node in nodes},
    )


def _physical_edge_paths(tree, drawing_layout):
    all_taxa = frozenset(str(name) for name in tree.leaf_names())
    taxon_sets = get_subtree_leaf_name_sets(tree)
    nodes_by_split: dict[Any, list[Any]] = {}
    for node in tree.traverse():
        if node.is_root or node not in drawing_layout.edge_paths:
            continue
        side = frozenset(str(name) for name in taxon_sets[node])
        split = canonical_split(side, all_taxa - side)
        nodes_by_split.setdefault(split, []).append(node)

    paths = {}
    for split, nodes in nodes_by_split.items():
        if len(nodes) == 1:
            node = nodes[0]
            path = list(drawing_layout.edge_paths[node])
            side = frozenset(str(name) for name in taxon_sets[node])
            if side == split[0]:
                path.reverse()
            paths[split] = path
            continue
        root_children = tree.get_children()
        if len(nodes) == 2 and {id(node) for node in nodes} == {
            id(child) for child in root_children
        }:
            side_a = next(
                node
                for node in nodes
                if frozenset(str(name) for name in taxon_sets[node]) == split[0]
            )
            side_b = nodes[0] if nodes[1] is side_a else nodes[1]
            path_a = list(reversed(drawing_layout.edge_paths[side_a]))
            path_b = list(drawing_layout.edge_paths[side_b])
            paths[split] = path_a + path_b[1:]
    return paths


def _point_on_path(path, fraction):
    points = [np.asarray(point, dtype=float) for point in path]
    if not points:
        return None
    if len(points) == 1:
        return tuple(points[0])
    segment_lengths = [
        float(np.linalg.norm(second - first))
        for first, second in zip(points, points[1:], strict=False)
    ]
    total_length = math.fsum(segment_lengths)
    if total_length <= 1e-15:
        return tuple(points[len(points) // 2])
    target = min(max(float(fraction), 0.0), 1.0) * total_length
    traversed = 0.0
    for index, (first, second, segment_length) in enumerate(
        zip(
            points,
            points[1:],
            segment_lengths,
            strict=True,
        )
    ):
        if target <= traversed + segment_length or index == len(segment_lengths) - 1:
            local_fraction = (target - traversed) / max(segment_length, 1e-15)
            return tuple(first + ((second - first) * local_fraction))
        traversed += segment_length
    return tuple(points[-1])


def _draw_branch_marker_overlays(ax, figure, tree, drawing_layout, branch_markers):
    edge_paths = _physical_edge_paths(tree, drawing_layout)
    resolved = []
    for marker in branch_markers:
        split = marker["split"]
        path = edge_paths.get(split)
        if path is None:
            raise ValueError(
                "A rooting marker split was not found in the displayed tree."
            )
        fraction = marker.get("position_fraction_from_side_a")
        plotted_fraction = 0.5 if fraction is None else float(fraction)
        point = _point_on_path(path, plotted_fraction)
        if point is None:
            raise ValueError("A rooting marker branch has no drawable path.")
        display_point = np.asarray(ax.transData.transform(point), dtype=float)
        point_units = display_point * (72.0 / figure.dpi)
        resolved.append((marker, point, point_units))

    parents = list(range(len(resolved)))

    def find(index):
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(first, second):
        first_root = find(first)
        second_root = find(second)
        if first_root != second_root:
            parents[max(first_root, second_root)] = min(first_root, second_root)

    largest_marker = max(
        (float(marker.get("size", 6.0)) for marker, _, _ in resolved),
        default=6.0,
    )
    cell_size = largest_marker + 2.0
    collision_headroom = 1.75
    neighbor_span = math.ceil(collision_headroom)
    marker_grid: dict[tuple[int, int], list[int]] = {}
    for index, (marker, _, point_units) in enumerate(resolved):
        cell = (
            math.floor(point_units[0] / cell_size),
            math.floor(point_units[1] / cell_size),
        )
        neighbor_indices = {
            previous
            for x_offset in range(-neighbor_span, neighbor_span + 1)
            for y_offset in range(-neighbor_span, neighbor_span + 1)
            for previous in marker_grid.get(
                (cell[0] + x_offset, cell[1] + y_offset),
                (),
            )
        }
        marker_size = float(marker.get("size", 6.0))
        for previous in neighbor_indices:
            previous_marker, _, previous_point = resolved[previous]
            separation = float(np.linalg.norm(point_units - previous_point))
            collision_distance = (
                marker_size + float(previous_marker.get("size", 6.0))
            ) / 2.0 + 2.0
            # Leave headroom for the axes shrink applied during final fitting.
            if separation <= collision_distance * collision_headroom:
                union(index, previous)
        marker_grid.setdefault(cell, []).append(index)

    grouped: dict[int, list[Any]] = {}
    for index, record in enumerate(resolved):
        grouped.setdefault(find(index), []).append(record)

    artists = []
    legend_handles = []
    for group in grouped.values():
        count = len(group)
        largest_marker = max(float(marker.get("size", 6.0)) for marker, _, _ in group)
        point_xs = [point_units[0] for _, _, point_units in group]
        point_ys = [point_units[1] for _, _, point_units in group]
        anchor_diameter = math.hypot(
            max(point_xs) - min(point_xs),
            max(point_ys) - min(point_ys),
        )
        fan_radius = (
            0.0
            if count == 1
            else anchor_diameter
            + max(
                5.0,
                (largest_marker + 2.0) / (2.0 * math.sin(math.pi / count)),
            )
        )
        for index, (marker, point, _) in enumerate(group):
            if count == 1:
                offset_x = 0.0
                offset_y = 0.0
            else:
                angle = (2.0 * math.pi * index / count) - (math.pi / 2.0)
                offset_x = fan_radius * math.cos(angle)
                offset_y = fan_radius * math.sin(angle)
                ax.annotate(
                    "",
                    xy=point,
                    xytext=(offset_x, offset_y),
                    textcoords="offset points",
                    arrowprops={
                        "arrowstyle": "-",
                        "color": marker["color"],
                        "linewidth": 0.45,
                    },
                    zorder=7,
                )
            transform = ax.transData + ScaledTranslation(
                offset_x / 72.0,
                offset_y / 72.0,
                figure.dpi_scale_trans,
            )
            filled = bool(marker.get("filled", True))
            artist = ax.plot(
                [point[0]],
                [point[1]],
                linestyle="None",
                marker=marker["marker"],
                markersize=float(marker.get("size", 6.0)),
                markerfacecolor=marker["color"] if filled else "white",
                markeredgecolor=marker["color"],
                markeredgewidth=1.0,
                transform=transform,
                zorder=8,
            )[0]
            artists.append(artist)
            legend_handles.append(
                Line2D(
                    [0],
                    [0],
                    linestyle="None",
                    marker=marker["marker"],
                    markersize=float(marker.get("size", 6.0)),
                    markerfacecolor=marker["color"] if filled else "white",
                    markeredgecolor=marker["color"],
                    markeredgewidth=1.0,
                    label=marker["label"],
                )
            )
    return artists, legend_handles


def _normalized_branch_markers(branch_markers):
    return [] if branch_markers is None else list(branch_markers)


def _branch_marker_artists(ax, figure, tree, drawing_layout, branch_markers):
    if not branch_markers:
        return [], []
    artists, legend_handles = _draw_branch_marker_overlays(
        ax,
        figure,
        tree,
        drawing_layout,
        branch_markers,
    )
    drawing_artists = [
        DrawingArtist(
            artist,
            kind="branch_marker",
            priority=90,
        )
        for artist in artists
    ]
    return drawing_artists, legend_handles


def _legend_is_needed(legend, *content_groups):
    return bool(legend) and any(bool(group) for group in content_groups)


def _tip_track_colors(
    tree,
    properties,
    mode,
    property_colors,
    categorical_palette,
    continuous_palette,
):
    mode = str(mode).strip().lower()
    if mode not in {"auto", "categorical", "continuous"}:
        raise ValueError("'--tip-track-type' must be auto, categorical, or continuous.")
    colors = {}
    legend_entries = []
    for property_name in properties:
        raw = {leaf: leaf.props.get(property_name) for leaf in tree.leaves()}
        present = [value for value in raw.values() if value not in (None, "")]
        numeric = {}
        for leaf, value in raw.items():
            try:
                number = float(value)
            except (TypeError, ValueError):
                continue
            if math.isfinite(number):
                numeric[leaf] = number
        continuous = mode == "continuous" or (
            mode == "auto" and present and len(numeric) == len(present)
        )
        if mode == "continuous" and len(numeric) != len(present):
            raise ValueError(
                "Tip-track property '{}' contains non-numeric values but "
                "--tip-track-type continuous was requested.".format(property_name)
            )
        if continuous and numeric:
            minimum = min(numeric.values())
            maximum = max(numeric.values())
            cmap = plt.get_cmap(continuous_palette)
            for leaf in tree.leaves():
                if leaf not in numeric:
                    colors[(leaf, property_name)] = "#d9d9d9"
                    continue
                fraction = (
                    0.5
                    if maximum <= minimum
                    else (numeric[leaf] - minimum) / (maximum - minimum)
                )
                colors[(leaf, property_name)] = mcolors.to_hex(cmap(fraction))
            if maximum <= minimum:
                legend_entries.append(
                    (
                        "{}: {:g}".format(property_name, minimum),
                        mcolors.to_hex(cmap(0.5)),
                    )
                )
            else:
                legend_entries.extend(
                    [
                        (
                            "{}: {:g}".format(property_name, minimum),
                            mcolors.to_hex(cmap(0.0)),
                        ),
                        (
                            "{}: {:g}".format(property_name, maximum),
                            mcolors.to_hex(cmap(1.0)),
                        ),
                    ]
                )
        else:
            values = sorted({str(value) for value in present})
            for leaf, value in raw.items():
                if value in (None, ""):
                    colors[(leaf, property_name)] = "#d9d9d9"
                else:
                    colors[(leaf, property_name)] = _property_color(
                        prop=property_name,
                        value=value,
                        property_colors=property_colors,
                        fallback_index=values.index(str(value)),
                        palette=categorical_palette,
                    )
            legend_entries.extend(
                (
                    "{}: {}".format(property_name, value),
                    _property_color(
                        prop=property_name,
                        value=value,
                        property_colors=property_colors,
                        fallback_index=index,
                        palette=categorical_palette,
                    ),
                )
                for index, value in enumerate(values)
            )
        if len(present) != len(raw):
            legend_entries.append(
                (
                    "{}: missing".format(property_name),
                    "#d9d9d9",
                )
            )
    return colors, legend_entries


def _tip_track_artist(ax, xy, offset_points, color, size_points):
    drawing = DrawingArea(size_points, size_points, 0.0, 0.0)
    drawing.add_artist(
        Rectangle(
            (0.0, 0.0),
            size_points,
            size_points,
            facecolor=color,
            edgecolor="#ffffff",
            linewidth=0.35,
        )
    )
    artist = AnnotationBbox(
        drawing,
        xy,
        xybox=offset_points,
        xycoords="data",
        boxcoords="offset points",
        frameon=False,
        box_alignment=(0.5, 0.5),
        annotation_clip=False,
        zorder=8,
    )
    ax.add_artist(artist)
    return artist


def _remap_leaf_dictionary(source_tree, target_tree, values):
    source_leaves = list(source_tree.leaves())
    target_leaves = list(target_tree.leaves())
    source_by_name = {str(leaf.name): leaf for leaf in source_leaves}
    target_by_name = {str(leaf.name): leaf for leaf in target_leaves}
    if set(source_by_name) != set(target_by_name):
        raise ValueError("Posterior sample tree has a different tip set.")
    return {
        target: values[source_by_name[name]]
        for name, target in target_by_name.items()
        if source_by_name[name] in values
    }


def _densitree_edge_key(node):
    return frozenset(str(name) for name in node.leaf_names())


def _densitree_topology_signature(tree):
    all_tips = frozenset(str(name) for name in tree.leaf_names())
    return tuple(
        sorted(
            (
                tuple(sorted(_densitree_edge_key(node)))
                for node in tree.traverse()
                if (
                    not node.is_root
                    and 1 < len(_densitree_edge_key(node)) < len(all_tips)
                )
            ),
            key=lambda clade: (len(clade), clade),
        )
    )


def _make_posterior_drawing_layouts(
    tree,
    posterior,
    posterior_topology,
    ladderize,
    *,
    layout,
    use_topology_depth,
    aspect_ratio,
    spiral_turns,
    angular_span,
    angular_center,
    terminal_extent_by_leaf,
    label_size_by_leaf,
    tip_spacing,
    subtree_packing,
    unrooted_method,
    daylight_iterations,
):
    """Render corresponding posterior trees with exactly the main-tree settings."""

    reference_nodes = list(tree.traverse(strategy="preorder"))
    layouts = []
    for sample_index in range(posterior.sample_count):
        posterior_tree = posterior_sample_tree(
            posterior_topology,
            posterior,
            sample_index,
        )
        if ladderize:
            posterior_tree.ladderize()
        sample_nodes = list(posterior_tree.traverse(strategy="preorder"))
        if len(sample_nodes) != len(reference_nodes):
            raise ValueError(
                "Posterior sample topology differs from the displayed topology."
            )
        if [node.name for node in posterior_tree.leaves()] != [
            node.name for node in tree.leaves()
        ]:
            raise ValueError(
                "Posterior sample topology has a different ordered tip set."
            )
        sample_layout = make_tree_layout(
            posterior_tree,
            layout=layout,
            use_topology_depth=use_topology_depth,
            aspect_ratio=aspect_ratio,
            spiral_turns=spiral_turns,
            angular_span=angular_span,
            angular_center=angular_center,
            terminal_extent_by_leaf=_remap_leaf_dictionary(
                tree,
                posterior_tree,
                terminal_extent_by_leaf,
            ),
            label_size_by_leaf=_remap_leaf_dictionary(
                tree,
                posterior_tree,
                label_size_by_leaf,
            ),
            tip_spacing=tip_spacing,
            subtree_packing=subtree_packing,
            unrooted_method=unrooted_method,
            daylight_iterations=daylight_iterations,
        )
        path_by_clade = {
            _densitree_edge_key(reference): sample_layout.edge_paths[sample]
            for reference, sample in zip(reference_nodes, sample_nodes, strict=True)
            if not reference.is_root
        }
        layouts.append((sample_layout.bounds, path_by_clade))
    return layouts


def _make_tree_collection_drawing_layouts(
    tree,
    sample_trees,
    *,
    layout,
    use_topology_depth,
    aspect_ratio,
    spiral_turns,
    angular_span,
    angular_center,
    terminal_extent_by_leaf,
    label_size_by_leaf,
    tip_spacing,
    subtree_packing,
    unrooted_method,
    daylight_iterations,
):
    """Render topology-varying trees against one fixed reference tip order."""

    if layout in {"unrooted", "fractal"}:
        raise ValueError(
            "'--densitree-trees' currently supports rectangular, slanted, "
            "cladogram, circular, radial, and spiral layouts; fixed tip "
            "positions are not defined for unrooted or fractal geometry."
        )
    alignment_packing = (
        subtree_packing
        if layout in {"rectangular", "circular", "spiral"}
        else "standard"
    )
    alignment_layout = make_tree_layout(
        tree,
        layout="rectangular",
        use_topology_depth=use_topology_depth,
        terminal_extent_by_leaf=terminal_extent_by_leaf,
        label_size_by_leaf=label_size_by_leaf,
        tip_spacing=tip_spacing,
        subtree_packing=alignment_packing,
    )
    fixed_leaf_position_by_name = {
        str(leaf.name): float(alignment_layout.ycoord[leaf])
        for leaf in alignment_layout.leaf_order
    }
    fixed_root_position = float(alignment_layout.ycoord[tree])
    layouts = []
    for sample_tree in sample_trees:
        sample_layout = make_tree_layout(
            sample_tree,
            layout=layout,
            use_topology_depth=use_topology_depth,
            aspect_ratio=aspect_ratio,
            spiral_turns=spiral_turns,
            angular_span=angular_span,
            angular_center=angular_center,
            terminal_extent_by_leaf=_remap_leaf_dictionary(
                tree,
                sample_tree,
                terminal_extent_by_leaf,
            ),
            label_size_by_leaf=_remap_leaf_dictionary(
                tree,
                sample_tree,
                label_size_by_leaf,
            ),
            tip_spacing=tip_spacing,
            subtree_packing=subtree_packing,
            unrooted_method=unrooted_method,
            daylight_iterations=daylight_iterations,
            fixed_leaf_position_by_name=fixed_leaf_position_by_name,
            fixed_root_position=fixed_root_position,
        )
        path_by_clade = {
            _densitree_edge_key(node): path
            for node, path in sample_layout.edge_paths.items()
        }
        layouts.append((sample_layout.bounds, path_by_clade))
    return layouts


def _resample_polyline(path, sample_count=32):
    points = np.asarray(path, dtype=float)
    if len(points) == 0:
        return np.zeros((sample_count, 2), dtype=float)
    if len(points) == 1:
        return np.repeat(points, sample_count, axis=0)
    segment_lengths = np.sqrt(np.sum(np.diff(points, axis=0) ** 2, axis=1))
    cumulative = np.concatenate(([0.0], np.cumsum(segment_lengths)))
    if cumulative[-1] <= 1e-15:
        return np.repeat(points[:1], sample_count, axis=0)
    targets = np.linspace(0.0, cumulative[-1], sample_count)
    return np.column_stack(
        [np.interp(targets, cumulative, points[:, dimension]) for dimension in range(2)]
    )


def _convex_hull(points):
    """Return the counter-clockwise convex hull of two-dimensional points."""

    unique = sorted({(float(x), float(y)) for x, y in points})
    if len(unique) <= 1:
        return np.asarray(unique, dtype=float)

    def cross(origin, first, second):
        return (first[0] - origin[0]) * (second[1] - origin[1]) - (
            first[1] - origin[1]
        ) * (second[0] - origin[0])

    lower: list[tuple[float, float]] = []
    for point in unique:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0.0:
            lower.pop()
        lower.append(point)
    upper: list[tuple[float, float]] = []
    for point in reversed(unique):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0.0:
            upper.pop()
        upper.append(point)
    return np.asarray(lower[:-1] + upper[:-1], dtype=float)


def _densitree_branch_envelope(
    paths,
    level,
    padding,
    sample_count=33,
    segments_per_polygon=4,
):
    """Return a tube covering a central fraction of whole sampled paths.

    Paths are selected as complete objects by their maximum displacement from
    the coordinate-wise median path.  Consecutive points of every retained
    resampled path are then enclosed by a padded convex hull.  The resulting
    strip therefore does not claim pointwise or joint posterior probability,
    but it does contain at least the requested empirical fraction of whole
    resampled paths.
    """

    sampled = np.stack(
        [_resample_polyline(path, sample_count=sample_count) for path in paths], axis=0
    )
    center_path = np.median(sampled, axis=0)
    maximum_displacement = np.max(
        np.linalg.norm(sampled - center_path[None, :, :], axis=2),
        axis=1,
    )
    retained_count = max(
        1,
        min(
            len(paths),
            int(math.ceil(float(level) * len(paths) - 1e-12)),
        ),
    )
    retained_indices = np.argsort(
        maximum_displacement,
        kind="stable",
    )[:retained_count]
    retained = sampled[retained_indices]
    diagonal = float(padding) / math.sqrt(2.0)
    offsets = np.asarray(
        [
            (padding, 0.0),
            (diagonal, diagonal),
            (0.0, padding),
            (-diagonal, diagonal),
            (-padding, 0.0),
            (-diagonal, -diagonal),
            (0.0, -padding),
            (diagonal, -diagonal),
        ]
    )
    polygons = []
    for start in range(0, sample_count - 1, segments_per_polygon):
        stop = min(start + segments_per_polygon, sample_count - 1)
        path_piece = retained[:, start : stop + 1, :].reshape(-1, 2)
        padded = (path_piece[:, None, :] + offsets[None, :, :]).reshape(-1, 2)
        hull = _convex_hull(padded)
        if len(hull) >= 3:
            polygons.append(hull)
    return polygons, retained_indices


def _densitree_topology_groups(posterior_layouts, sample_trees):
    """Group rendered samples without interpolating between tree topologies."""

    if not posterior_layouts:
        return []
    if not sample_trees:
        return [((), list(posterior_layouts))]
    grouped: dict[tuple, list[Any]] = {}
    for sample_tree, rendered in zip(sample_trees, posterior_layouts, strict=True):
        signature = _densitree_topology_signature(sample_tree)
        grouped.setdefault(signature, []).append(rendered)
    return sorted(
        grouped.items(),
        key=lambda item: (-len(item[1]), item[0]),
    )


def _compound_polygon_path(polygons):
    """Combine overlapping envelope pieces into one consistently filled path."""

    vertices = []
    codes = []
    for polygon in polygons:
        points = np.asarray(polygon, dtype=float)
        if len(points) < 3:
            continue
        vertices.extend(points)
        codes.extend(
            [
                MatplotlibPath.MOVETO,
                *([MatplotlibPath.LINETO] * (len(points) - 1)),
            ]
        )
        vertices.append(points[0])
        codes.append(MatplotlibPath.CLOSEPOLY)
    if not vertices:
        return None
    return MatplotlibPath(np.asarray(vertices), np.asarray(codes, dtype=np.uint8))


def _age_interval_path(node, drawing_layout):
    try:
        age = float(node.props["age"])
        low = float(node.props["age_ci_low"])
        high = float(node.props["age_ci_high"])
    except (KeyError, TypeError, ValueError):
        return None
    if not all(math.isfinite(value) for value in (age, low, high)):
        return None
    center = np.asarray(
        [drawing_layout.xcoord[node], drawing_layout.ycoord[node]],
        dtype=float,
    )
    direction = None
    scale = None
    if node.is_root:
        if drawing_layout.spatial or not drawing_layout.root_path:
            return None
        vector = center - np.asarray(drawing_layout.root_path[0], dtype=float)
        length = float(np.linalg.norm(vector))
        root_age = max(age, 1e-15)
        display_span = max(
            abs(float(drawing_layout.xcoord[other]) - float(center[0]))
            for other in drawing_layout.xcoord
        )
        if length > 1e-15 and display_span > 1e-15:
            direction = vector / length
            scale = display_span / root_age
    elif drawing_layout.name in {"rectangular", "slanted"}:
        direction = np.asarray([1.0, 0.0])
        scale = 1.0
    elif drawing_layout.name in {"circular", "radial"}:
        root_center = np.asarray(
            [
                drawing_layout.xcoord[node.root],
                drawing_layout.ycoord[node.root],
            ]
        )
        vector = center - root_center
        length = float(np.linalg.norm(vector))
        if length > 1e-15:
            direction = vector / length
            scale = 1.0
    elif node.dist not in (None, 0):
        path = drawing_layout.edge_paths.get(node, ())
        for start, end in reversed(list(zip(path, path[1:], strict=False))):
            vector = np.asarray(end, dtype=float) - np.asarray(start, dtype=float)
            length = float(np.linalg.norm(vector))
            if length > 1e-15:
                direction = vector / length
                scale = length / float(node.dist)
                break
    if direction is None:
        for child in node.get_children():
            vector = (
                np.asarray(
                    [drawing_layout.xcoord[child], drawing_layout.ycoord[child]],
                    dtype=float,
                )
                - center
            )
            length = float(np.linalg.norm(vector))
            child_dist = 0.0 if child.dist is None else float(child.dist)
            if length > 1e-15 and child_dist > 0.0:
                direction = vector / length
                scale = length / child_dist
                break
    if direction is None or scale is None:
        return None
    younger = center + direction * max(age - low, 0.0) * scale
    older = center - direction * max(high - age, 0.0) * scale
    return np.vstack((older, younger))
