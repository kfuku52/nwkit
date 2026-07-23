import math
import os
import sys

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
from matplotlib.path import Path as MplPath
from matplotlib.patches import Patch

from nwkit import __version__
from nwkit.util import (
    extract_species_label,
    is_rooted,
    read_tip_table,
    read_tree,
    validate_unique_named_leaves,
)


TREE_LINE_CAPSTYLE = 'round'
FONT_FAMILY = 'Helvetica'
FONT_SIZE_PT = 8.0
TIP_LABEL_GAP_PT = 0.5
TIP_ROW_PITCH_PT = FONT_SIZE_PT + TIP_LABEL_GAP_PT
TIP_IMAGE_SIZE_PT = 18.0
TIP_IMAGE_GAP_PT = 4.0
SUPPORT_LABEL_OFFSET_PT = 0.0
FIGURE_WIDTH_IN = 3.6
SPECIATION_COLOR = (0.0, 0.0, 1.0)
DUPLICATION_COLOR = (1.0, 0.0, 0.0)
BRANCH_COLOR = '#000000'
LABEL_COLOR = '#2d2d2d'
LEGEND_EDGE_COLOR = 'white'


matplotlib.rcParams['font.family'] = [FONT_FAMILY]
matplotlib.rcParams['font.sans-serif'] = [FONT_FAMILY]
matplotlib.rcParams['svg.fonttype'] = 'none'


def _resolve_image_format(outfile, image_format):
    if image_format != 'auto':
        return image_format
    ext = os.path.splitext(str(outfile))[1].lower()
    if ext in ('.pdf', '.png', '.svg'):
        return ext[1:]
    return 'pdf'


def _has_positive_branch_length(tree):
    for node in tree.traverse():
        if node.is_root:
            continue
        if (node.dist is not None) and (float(node.dist) > 0):
            return True
    return False


def _get_tree_plot_coordinates(tree, use_topology_depth=False):
    xcoord = dict()
    ycoord = dict()
    leaf_order = list()

    def assign_x(node, parent_x):
        xcoord[node] = float(parent_x)
        for child in node.get_children():
            if use_topology_depth:
                child_dist = 1.0
            else:
                child_dist = 0.0 if (child.dist is None) else float(child.dist)
            assign_x(node=child, parent_x=parent_x + child_dist)

    def assign_y(node, current_y):
        if node.is_leaf:
            ycoord[node] = float(current_y)
            leaf_order.append(node)
            return current_y + 1
        child_ys = list()
        for child in node.get_children():
            current_y = assign_y(node=child, current_y=current_y)
            child_ys.append(ycoord[child])
        if len(child_ys) == 0:
            ycoord[node] = float(current_y)
            return current_y + 1
        ycoord[node] = float(sum(child_ys) / len(child_ys))
        return current_y

    assign_x(node=tree, parent_x=0.0)
    assign_y(node=tree, current_y=0)
    return xcoord, ycoord, leaf_order


def _get_species_by_leaf(tree, args):
    species_by_leaf = dict()
    is_all_parsed = True
    num_leaf = 0
    for leaf in tree.leaves():
        num_leaf += 1
        species_label = extract_species_label(leaf.name or '', args=args)
        if species_label is None:
            is_all_parsed = False
            continue
        species_by_leaf[leaf] = species_label
    return species_by_leaf, (is_all_parsed and (num_leaf > 0))

def _get_species_overlap_node_types(tree, args, require_all_tip_labels=False):
    species_by_leaf, all_tip_labels_parsed = _get_species_by_leaf(tree=tree, args=args)
    if bool(require_all_tip_labels) and (not all_tip_labels_parsed):
        return dict(), all_tip_labels_parsed
    species_to_bit = dict()
    species_mask_by_node = dict()
    has_missing_species_by_node = dict()
    node_type_by_node = dict()
    for node in tree.traverse(strategy='postorder'):
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
        node_type_by_node[node] = 'duplication' if is_duplication else 'speciation'
    return node_type_by_node, all_tip_labels_parsed


def _apply_font_style(font_size=FONT_SIZE_PT, font_family=FONT_FAMILY):
    matplotlib.rcParams['font.size'] = font_size
    matplotlib.rc('xtick', labelsize=font_size)
    matplotlib.rc('ytick', labelsize=font_size)
    matplotlib.rc('font', size=font_size, family=font_family)
    matplotlib.rc('axes', titlesize=font_size, labelsize=font_size)
    matplotlib.rc('legend', fontsize=font_size)
    matplotlib.rc('figure', titlesize=font_size)


def _has_meaningful_support(node):
    if node.is_root or node.is_leaf:
        return False
    support = node.support
    if support is None:
        return False
    support_value = float(support)
    if abs(support_value - (-999999.0)) < 10 ** -9:
        return False
    return True


def _format_support_value(support):
    support_value = float(support)
    if abs(support_value - round(support_value)) < 10 ** -9:
        return str(int(round(support_value)))
    return '{:g}'.format(support_value)


def _read_trait_table(path, group_by, tree, unmatched='warn', missing_values=None):
    trait_df, _, _ = read_tip_table(
        path,
        option_name='--trait',
        tree_leaf_names=tree.leaf_names(),
        required_columns=(group_by,),
        unmatched=unmatched,
        missing_values=missing_values,
    )
    tree_leaf_name_set = set(tree.leaf_names())
    trait_df = trait_df[trait_df['leaf_name'].isin(tree_leaf_name_set)].copy()
    leaf_to_group = dict()
    for _, row in trait_df.iterrows():
        if pd.isna(row[group_by]):
            continue
        leaf_to_group[str(row['leaf_name'])] = str(row[group_by])
    return leaf_to_group


def _read_tip_image_manifest(
    path,
    tree,
    image_root=None,
    unmatched='warn',
    missing_values=None,
):
    if path in (None, ''):
        return dict()
    if path == '-' and image_root in (None, ''):
        raise ValueError(
            "'--tip-image-root' is required when '--tip-image-manifest' reads from STDIN."
        )
    validate_unique_named_leaves(
        tree,
        option_name='--infile',
        context=" for '--tip-image-manifest'",
    )
    manifest_df, _, _ = read_tip_table(
        path,
        option_name='--tip-image-manifest',
        tree_leaf_names=tree.leaf_names(),
        required_columns=('local_path',),
        unmatched=unmatched,
        missing_values=missing_values,
    )
    if image_root in (None, ''):
        base_dir = os.path.dirname(os.path.realpath(path))
    else:
        base_dir = os.path.realpath(image_root)
    if not os.path.isdir(base_dir):
        raise ValueError(
            "'--tip-image-root' is not a directory: {}".format(base_dir)
        )

    tree_leaf_names = set(str(name) for name in tree.leaf_names())
    path_by_leaf = dict()
    missing_path_leaf_names = list()
    for _, row in manifest_df.iterrows():
        leaf_name = str(row['leaf_name'])
        if leaf_name not in tree_leaf_names:
            continue
        raw_local_path = row['local_path']
        if pd.isna(raw_local_path) or str(raw_local_path).strip() == '':
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
        message = (
            "Tree tips have missing local_path values in --tip-image-manifest: {}".format(
                ' '.join(sorted(missing_path_leaf_names))
            )
        )
        if unmatched == 'error':
            raise ValueError(message)
        if unmatched == 'warn':
            sys.stderr.write(message + '\n')
    return path_by_leaf


def _load_tip_image(path):
    from nwkit.image import load_pillow_modules, rasterize_svg_to_image

    Image, _, ImageOps = load_pillow_modules()
    try:
        if os.path.splitext(path)[1].lower() == '.svg':
            image = rasterize_svg_to_image(path)
        else:
            with Image.open(path) as source_image:
                image = ImageOps.exif_transpose(source_image).convert('RGBA')
                image.load()
    except (OSError, ValueError) as exc:
        raise ValueError("Failed to read tip image '{}': {}".format(path, exc)) from exc
    return np.asarray(image)


def _load_tip_images(path_by_leaf):
    image_by_leaf = dict()
    image_by_path = dict()
    for leaf_name, path in path_by_leaf.items():
        image = image_by_path.get(path)
        if image is None:
            image = _load_tip_image(path)
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
            interpolation='hanning',
        )
        annotation = AnnotationBbox(
            image_artist,
            (0.5, ycoord[leaf]),
            xycoords=('axes fraction', 'data'),
            box_alignment=(0.5, 0.5),
            frameon=False,
            pad=0.0,
            annotation_clip=False,
        )
        annotation.set_clip_on(False)
        ax.add_artist(annotation)


def _assign_group_colors(group_names, palette='tab10'):
    if len(group_names) == 0:
        return dict()
    cmap = plt.get_cmap(palette)
    colors = dict()
    for index, group_name in enumerate(sorted(group_names)):
        colors[group_name] = mcolors.to_hex(cmap(index % max(getattr(cmap, 'N', len(group_names)), 1)))
    return colors


def _parse_property_colors(raw_colors):
    colors = dict()
    for raw in raw_colors or []:
        if '=' not in str(raw):
            raise ValueError("--property-color must use VALUE=COLOR syntax: {}".format(raw))
        value, color = str(raw).split('=', 1)
        value = value.strip()
        color = color.strip()
        if value == '' or color == '':
            raise ValueError("--property-color must use VALUE=COLOR syntax: {}".format(raw))
        try:
            mcolors.to_rgba(color)
        except ValueError as exc:
            raise ValueError("Invalid color for --property-color '{}': {}".format(value, color)) from exc
        colors[value] = color
    return colors


def _parse_property_names(raw):
    if raw in (None, ''):
        return list()
    names = [name.strip() for name in str(raw).split(',') if name.strip()]
    if len(names) != len(set(names)):
        raise ValueError('Property names must not be repeated: {}'.format(raw))
    return names


def _node_matches_target(node, raw_target):
    targets = {value.strip().lower() for value in str(raw_target or 'all').split(',') if value.strip()}
    allowed = {'all', 'root', 'intnode', 'leaf'}
    unknown = targets - allowed
    if unknown:
        raise ValueError('Unsupported node target(s): {}'.format(', '.join(sorted(unknown))))
    if 'all' in targets:
        return True
    if node.is_root:
        return 'root' in targets
    if node.is_leaf:
        return 'leaf' in targets
    return 'intnode' in targets


def _coerce_filter_value(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return str(value)


def _matches_property_filters(node, raw_filters, option_name='--node-label-filter'):
    operators = {
        'ge': lambda left, right: left >= right,
        'gt': lambda left, right: left > right,
        'le': lambda left, right: left <= right,
        'lt': lambda left, right: left < right,
        'eq': lambda left, right: left == right,
        'ne': lambda left, right: left != right,
    }
    for raw in raw_filters or []:
        parts = str(raw).split(':', 2)
        if len(parts) != 3 or parts[1].lower() not in operators:
            raise ValueError(
                '{} must use PROPERTY:OP:VALUE syntax with '
                'OP in ge,gt,le,lt,eq,ne: {}'.format(option_name, raw)
            )
        prop, operator, expected = parts[0].strip(), parts[1].lower(), parts[2].strip()
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


def _property_color(prop, value, property_colors, fallback_index=0, palette='tab10'):
    candidates = ["{}:{}".format(prop, value), str(prop), str(value)]
    if str(prop).startswith('asr_p_'):
        candidates.append(str(prop)[len('asr_p_'):])
    for candidate in candidates:
        if candidate in property_colors:
            return property_colors[candidate]
    cmap = plt.get_cmap(palette)
    return mcolors.to_hex(cmap(fallback_index % max(getattr(cmap, 'N', 1), 1)))


def _draw_probability_pie(ax, x, y, probabilities, colors, marker_size_pt):
    total = sum(probabilities)
    if total <= 0.0:
        return
    normalized = [probability / total for probability in probabilities]
    cumulative = 0.0
    for probability, color in zip(normalized, colors):
        if probability <= 0.0:
            continue
        start = cumulative
        cumulative += probability
        angles = [
            (math.pi / 2.0) - (2.0 * math.pi * (start + probability * index / 24.0))
            for index in range(25)
        ]
        vertices = [(0.0, 0.0)] + [(math.cos(angle), math.sin(angle)) for angle in angles] + [(0.0, 0.0)]
        codes = [MplPath.MOVETO] + ([MplPath.LINETO] * len(angles)) + [MplPath.CLOSEPOLY]
        ax.scatter(
            [x],
            [y],
            s=marker_size_pt ** 2,
            marker=MplPath(vertices, codes),
            facecolor=color,
            edgecolor='none',
            linewidth=0.0,
            zorder=6,
        )
    ax.scatter(
        [x],
        [y],
        s=marker_size_pt ** 2,
        marker='o',
        facecolor='none',
        edgecolor=LEGEND_EDGE_COLOR,
        linewidth=0.45,
        zorder=7,
    )


def _format_property_value(value, decimals=2):
    try:
        return ('{:.%df}' % int(decimals)).format(float(value))
    except (TypeError, ValueError):
        return str(value)


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
    tip_label_position='aligned',
    root_marker='none',
    root_marker_color='#0072B2',
    root_marker_size=None,
    tip_badge_property=None,
    tip_badge_missing_label=None,
    node_pie_properties=None,
    node_pie_target='root,intnode',
    node_pie_leaf_filters=None,
    node_label_property=None,
    node_label_target='intnode',
    node_label_filters=None,
    node_label_decimals=2,
    node_label_prefix='',
    property_colors=None,
    legend=True,
    transparent=False,
    trait_palette='tab10',
    tip_image_by_leaf=None,
    tip_image_size=TIP_IMAGE_SIZE_PT,
    tip_image_gap=TIP_IMAGE_GAP_PT,
):
    _apply_font_style(font_size=font_size, font_family=font_family)
    matplotlib.rcParams['svg.hashsalt'] = 'nwkit'
    property_colors = property_colors or {}
    node_pie_properties = node_pie_properties or []
    node_pie_leaf_filters = node_pie_leaf_filters or []
    node_label_filters = node_label_filters or []
    tip_image_by_leaf = tip_image_by_leaf or {}
    use_topology_depth = (not _has_positive_branch_length(tree))
    if use_topology_depth:
        sys.stderr.write('Tree has no positive branch lengths; drawing x positions by topology depth.\n')
    xcoord, ycoord, leaf_order = _get_tree_plot_coordinates(tree=tree, use_topology_depth=use_topology_depth)
    x_values = list(xcoord.values())
    x_min = min(x_values) if len(x_values) > 0 else 0.0
    x_max = max(x_values) if len(x_values) > 0 else 1.0
    x_span = max(x_max - x_min, 1.0)
    max_tip_label_chars = max([len(leaf.name or '') for leaf in leaf_order], default=1)

    fig_width = float(figure_width)
    if fig_width <= 0.0:
        raise ValueError('--figure-width must be greater than zero.')
    tip_image_size_pt = float(tip_image_size)
    tip_image_gap_pt = float(tip_image_gap)
    if tip_image_size_pt <= 0.0:
        raise ValueError('--tip-image-size must be greater than zero.')
    if tip_image_gap_pt < 0.0:
        raise ValueError('--tip-image-gap must be zero or greater.')
    image_panel_width_in = 0.0
    if tip_image_by_leaf:
        image_panel_width_in = (
            tip_image_size_pt + (2.0 * tip_image_gap_pt)
        ) / 72.0
    main_panel_width_in = fig_width - image_panel_width_in
    if main_panel_width_in <= 0.2:
        raise ValueError(
            '--figure-width is too small for the requested tip-image column.'
        )
    if label_panel_width is None:
        label_panel_width_in = min(
            max(0.8, max_tip_label_chars * 0.07),
            main_panel_width_in * 0.58,
        )
    else:
        label_panel_width_in = float(label_panel_width)
        if label_panel_width_in <= 0.0 or label_panel_width_in >= main_panel_width_in:
            raise ValueError(
                '--label-panel-width must be greater than zero and leave room '
                'for the tree and tip-image column within --figure-width.'
            )
    tree_panel_width_in = main_panel_width_in - label_panel_width_in
    if tree_panel_width_in < 0.2:
        raise ValueError(
            '--figure-width is too small for the requested label and tip-image panels.'
        )
    row_pitch_pt = float(font_size) + TIP_LABEL_GAP_PT
    if tip_image_by_leaf:
        row_pitch_pt = max(
            row_pitch_pt,
            tip_image_size_pt + TIP_LABEL_GAP_PT,
        )
    row_pitch_in = row_pitch_pt / 72.0
    tree_panel_height_in = max(len(leaf_order), 1) * row_pitch_in
    property_legend_needed = bool(tip_badge_property or node_pie_properties) and bool(property_colors)
    legend_needed = bool(legend) and bool(node_type_by_node or group_color_by_name or property_legend_needed)
    top_margin_in = (28.0 / 72.0) if legend_needed else (4.0 / 72.0)
    bottom_margin_in = 3.0 / 72.0
    left_margin_in = 2.0 / 72.0
    right_margin_in = 2.0 / 72.0
    if figure_height is None:
        fig_height = tree_panel_height_in + top_margin_in + bottom_margin_in
    else:
        fig_height = float(figure_height)
        if fig_height <= top_margin_in + bottom_margin_in:
            raise ValueError('--figure-height is too small for the required margins.')
        if (
            tip_image_by_leaf
            and fig_height + 10 ** -9
            < tree_panel_height_in + top_margin_in + bottom_margin_in
        ):
            raise ValueError(
                '--figure-height is too small for non-overlapping tip images. '
                'Increase it or reduce --tip-image-size.'
            )

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    for node in tree.traverse():
        if node.is_leaf:
            continue
        children = node.get_children()
        if len(children) <= 1:
            continue
        for child in children:
            ax.plot(
                [xcoord[node], xcoord[node]],
                [ycoord[node], ycoord[child]],
                color=branch_color,
                linewidth=0.8,
                zorder=1,
                solid_capstyle=TREE_LINE_CAPSTYLE,
            )

    for node in tree.traverse():
        if node.is_root:
            continue
        edge_color = terminal_branch_color if node.is_leaf and terminal_branch_color else branch_color
        ax.plot(
            [xcoord[node.up], xcoord[node]],
            [ycoord[node], ycoord[node]],
            color=edge_color,
            linewidth=0.8,
            zorder=2,
            solid_capstyle=TREE_LINE_CAPSTYLE,
        )

    root_stub = max(x_span * 0.03, 0.03)
    ax.plot(
        [xcoord[tree] - root_stub, xcoord[tree]],
        [ycoord[tree], ycoord[tree]],
        color=branch_color,
        linewidth=0.8,
        zorder=2,
        solid_capstyle=TREE_LINE_CAPSTYLE,
    )

    marker_area = 18.0 * (0.5 ** 2)
    marker_size_pt = 4.8 * 0.5
    legend_handles = list()
    if len(node_type_by_node) > 0:
        for node, node_type in node_type_by_node.items():
            marker_color = DUPLICATION_COLOR if (node_type == 'duplication') else SPECIATION_COLOR
            ax.scatter(
                [xcoord[node]],
                [ycoord[node]],
                s=marker_area,
                marker='o',
                facecolor=marker_color,
                edgecolor=LEGEND_EDGE_COLOR,
                linewidth=0.4,
                zorder=5,
            )
        legend_handles.extend([
            Line2D(
                [0],
                [0],
                marker='o',
                linestyle='None',
                markerfacecolor=SPECIATION_COLOR,
                markeredgecolor=LEGEND_EDGE_COLOR,
                markeredgewidth=0.4,
                markersize=marker_size_pt,
                label='Speciation node',
            ),
            Line2D(
                [0],
                [0],
                marker='o',
                linestyle='None',
                markerfacecolor=DUPLICATION_COLOR,
                markeredgecolor=LEGEND_EDGE_COLOR,
                markeredgewidth=0.4,
                markersize=marker_size_pt,
                label='Duplication node',
            ),
        ])
    if len(group_color_by_name) > 0:
        legend_handles.extend(
            Patch(facecolor=color, edgecolor='none', label=group_name)
            for group_name, color in sorted(group_color_by_name.items())
        )
    if property_legend_needed:
        legend_handles.extend(
            Patch(facecolor=color, edgecolor='none', label=value)
            for value, color in property_colors.items()
        )
    if legend and len(legend_handles) > 0:
        ax.legend(
            handles=legend_handles,
            loc='lower left',
            bbox_to_anchor=(0.0, 1.005),
            frameon=False,
            borderaxespad=0.1,
            handletextpad=0.3,
            labelspacing=0.2,
            ncol=1,
            prop={'family': font_family, 'size': font_size},
        )

    for node in tree.traverse():
        if not support_labels:
            continue
        if not _has_meaningful_support(node):
            continue
        support_value = float(node.support)
        if (support_min is not None) and (support_value < support_min):
            continue
        parent_x = xcoord[node.up]
        node_x = xcoord[node]
        label_x = parent_x + (node_x - parent_x) * 0.5
        ax.annotate(
            _format_support_value(node.support),
            xy=(label_x, ycoord[node]),
            xytext=(0.0, SUPPORT_LABEL_OFFSET_PT),
            textcoords='offset points',
            va='bottom',
            ha='center',
            fontsize=font_size,
            fontfamily=font_family,
            color=LABEL_COLOR,
        )

    data_per_inch = x_span / tree_panel_width_in
    leaf_pies_drawn = any(
        node_pie_properties
        and _node_matches_target(leaf, node_pie_target)
        and (
            (not node_pie_leaf_filters)
            or _matches_property_filters(
                leaf,
                node_pie_leaf_filters,
                option_name='--node-pie-leaf-filter',
            )
        )
        and all(prop in leaf.props for prop in node_pie_properties)
        for leaf in leaf_order
    )
    tip_pie_clearance = 0.0
    if leaf_pies_drawn:
        tip_pie_radius_pt = max(float(font_size), 4.5) / 2.0
        tip_pie_clearance = (
            (tip_pie_radius_pt + 2.0) / 72.0
        ) * data_per_inch
    label_offset = max(x_span * 0.02, 0.06, tip_pie_clearance)
    for leaf in leaf_order:
        if tip_label_position == 'branch-end':
            label_x = xcoord[leaf] + label_offset
        elif tip_label_position == 'aligned':
            label_x = x_max + label_offset
        else:
            raise ValueError("Unsupported '--tip-label-position': {}".format(tip_label_position))
        ax.text(
            label_x,
            ycoord[leaf],
            leaf.name or '',
            va='center',
            ha='left',
            fontsize=font_size,
            fontfamily=font_family,
            color=leaf_label_color_by_leaf.get(leaf, LABEL_COLOR),
        )
        if tip_badge_property not in (None, ''):
            raw_value = leaf.props.get(tip_badge_property)
            is_missing = raw_value in (None, '')
            if is_missing and tip_badge_missing_label in (None, ''):
                continue
            value = str(tip_badge_missing_label if is_missing else raw_value)
            badge_color = _property_color(
                prop=tip_badge_property,
                value=value,
                property_colors=property_colors,
                palette=trait_palette,
            )
            approximate_label_width = max(len(leaf.name or ''), 1) * (font_size * 0.58 / 72.0) * data_per_inch
            approximate_badge_width = max(len(value), 1) * (font_size * 0.58 / 72.0) * data_per_inch
            badge_padding = (font_size * 0.16 / 72.0) * data_per_inch
            badge_gap = (4.0 / 72.0) * data_per_inch
            badge_offset = (
                approximate_label_width
                + badge_gap
                + (approximate_badge_width / 2.0)
                + badge_padding
            )
            ax.text(
                label_x + badge_offset,
                ycoord[leaf],
                value,
                va='center',
                ha='center',
                fontsize=font_size,
                fontfamily=font_family,
                fontweight='bold',
                color=LABEL_COLOR,
                bbox={
                    'boxstyle': 'round,pad=0.16',
                    'facecolor': badge_color,
                    'edgecolor': LEGEND_EDGE_COLOR,
                    'linewidth': 0.4,
                    'alpha': 0.92,
                },
                zorder=8,
            )

    for node in tree.traverse():
        leaf_passes_pie_filters = (
            (not node.is_leaf)
            or (not node_pie_leaf_filters)
            or _matches_property_filters(
                node,
                node_pie_leaf_filters,
                option_name='--node-pie-leaf-filter',
            )
        )
        if (
            node_pie_properties
            and _node_matches_target(node, node_pie_target)
            and leaf_passes_pie_filters
        ):
            if all(prop in node.props for prop in node_pie_properties):
                probabilities = [float(node.props[prop]) for prop in node_pie_properties]
                colors = [
                    _property_color(
                        prop=prop,
                        value=prop,
                        property_colors=property_colors,
                        fallback_index=index,
                        palette=trait_palette,
                    )
                    for index, prop in enumerate(node_pie_properties)
                ]
                _draw_probability_pie(
                    ax,
                    xcoord[node],
                    ycoord[node],
                    probabilities=probabilities,
                    colors=colors,
                    marker_size_pt=max(font_size, 4.5),
                )
        if (
            node_label_property not in (None, '')
            and node_label_property in node.props
            and _node_matches_target(node, node_label_target)
            and _matches_property_filters(node, node_label_filters)
        ):
            label_offset_points = (-3.0, 3.0) if node.is_leaf else (2.0, 2.0)
            ax.annotate(
                '{}{}'.format(
                    node_label_prefix,
                    _format_property_value(node.props[node_label_property], decimals=node_label_decimals),
                ),
                xy=(xcoord[node], ycoord[node]),
                xytext=label_offset_points,
                textcoords='offset points',
                va='bottom',
                ha='right' if node.is_leaf else 'left',
                fontsize=font_size,
                fontfamily=font_family,
                color=LABEL_COLOR,
                zorder=9,
            )

    root_marker_size_pt = max(float(font_size), 5.0)
    if root_marker_size is not None:
        root_marker_size_pt = float(root_marker_size)
        if root_marker_size_pt <= 0.0:
            raise ValueError('--root-marker-size must be greater than zero.')
    if root_marker != 'none':
        marker_by_name = {'circle': 'o', 'diamond': 'D'}
        if root_marker not in marker_by_name:
            raise ValueError("Unsupported '--root-marker': {}".format(root_marker))
        ax.scatter(
            [xcoord[tree]],
            [ycoord[tree]],
            s=root_marker_size_pt ** 2,
            marker=marker_by_name[root_marker],
            facecolor=root_marker_color,
            edgecolor=LEGEND_EDGE_COLOR,
            linewidth=0.45,
            zorder=8,
        )

    label_panel_span = x_span * (label_panel_width_in / tree_panel_width_in)
    root_has_pie = (
        bool(node_pie_properties)
        and _node_matches_target(tree, node_pie_target)
        and all(prop in tree.props for prop in node_pie_properties)
    )
    root_has_symbol = (root_marker != 'none') or root_has_pie
    root_symbol_diameter_pt = 0.0
    if root_marker != 'none':
        root_symbol_diameter_pt = root_marker_size_pt
    if root_has_pie:
        root_symbol_diameter_pt = max(root_symbol_diameter_pt, max(float(font_size), 4.5))
    root_symbol_margin = 0.0
    if root_has_symbol:
        root_symbol_radius_pt = root_symbol_diameter_pt / 2.0
        root_symbol_margin = (
            (root_symbol_radius_pt + 1.5) / 72.0
        ) * data_per_inch
    left_plot_margin = max(root_stub * 1.2, root_symbol_margin)
    ax.set_xlim(x_min - left_plot_margin, x_max + label_offset + label_panel_span)
    if len(leaf_order) == 0:
        ax.set_ylim(0.5, -0.5)
    else:
        ax.set_ylim(float(len(leaf_order)) - 0.5, -0.5)
    ax.axis('off')
    axes_left = left_margin_in / fig_width
    axes_right = 1.0 - ((right_margin_in + image_panel_width_in) / fig_width)
    axes_top = 1.0 - (top_margin_in / fig_height)
    axes_bottom = bottom_margin_in / fig_height
    fig.subplots_adjust(
        left=axes_left,
        right=axes_right,
        top=axes_top,
        bottom=axes_bottom,
    )
    if tip_image_by_leaf:
        image_ax = fig.add_axes(
            [
                axes_right,
                axes_bottom,
                image_panel_width_in / fig_width,
                axes_top - axes_bottom,
            ],
            sharey=ax,
        )
        image_ax.set_xlim(0.0, 1.0)
        image_ax.set_ylim(ax.get_ylim())
        image_ax.set_facecolor('none')
        image_ax.axis('off')
        _draw_tip_images(
            ax=image_ax,
            leaf_order=leaf_order,
            ycoord=ycoord,
            image_by_leaf=tip_image_by_leaf,
            image_size_pt=tip_image_size_pt,
        )
    metadata = {'Creator': 'NWKIT {}'.format(__version__)}
    if image_format == 'svg':
        metadata['Date'] = None
    elif image_format == 'pdf':
        metadata['CreationDate'] = None
        metadata['ModDate'] = None
    fig.savefig(
        outfile,
        format=image_format,
        dpi=300,
        transparent=bool(transparent),
        metadata=metadata,
    )
    plt.close(fig)


def draw_main(args):
    if args.outfile == '-':
        raise ValueError("STDOUT is not supported for 'draw'. Use --outfile PATH.")
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    if bool(args.ladderize):
        tree.ladderize()
    trait_path = getattr(args, 'trait', None)
    group_by = getattr(args, 'group_by', None)
    trait_palette = getattr(args, 'trait_palette', 'tab10')
    support_labels = getattr(args, 'support_labels', True)
    support_min = getattr(args, 'support_min', None)
    property_colors = _parse_property_colors(getattr(args, 'property_color', None))
    node_pie_properties = _parse_property_names(getattr(args, 'node_pie_properties', None))
    tip_image_manifest = getattr(args, 'tip_image_manifest', None)
    tip_image_path_by_leaf = _read_tip_image_manifest(
        path=tip_image_manifest,
        tree=tree,
        image_root=getattr(args, 'tip_image_root', None),
        unmatched=getattr(args, 'unmatched', 'warn'),
        missing_values=getattr(args, 'missing_values', None),
    )
    setattr(
        args,
        '_nwkit_tip_image_paths',
        sorted(set(tip_image_path_by_leaf.values())),
    )
    tip_image_by_leaf = _load_tip_images(tip_image_path_by_leaf)
    if tip_image_by_leaf:
        sys.stderr.write(
            'Loaded tip images for {} tree tip(s) from: {}\n'.format(
                len(tip_image_by_leaf),
                tip_image_manifest,
            )
        )
    if (trait_path not in ['', None]) and (group_by in ['', None]):
        raise ValueError("'--group-by' is required when '--trait' is specified.")
    image_format = _resolve_image_format(outfile=args.outfile, image_format=str(args.image_format).lower())
    node_plot_mode = str(args.species_overlap_node_plot).strip().lower()
    if node_plot_mode == 'no':
        node_type_by_node = dict()
    else:
        if not is_rooted(tree):
            if node_plot_mode == 'yes':
                raise ValueError(
                    "Speciation/duplication node plotting requires a rooted tree. "
                    "Use '--species-overlap-node-plot no' for unrooted trees."
                )
            node_type_by_node = dict()
            sys.stderr.write(
                'Skipping speciation/duplication node markers because the input tree is unrooted.\n'
            )
        else:
            node_type_by_node, all_tip_labels_parsed = _get_species_overlap_node_types(
                tree=tree,
                args=args,
                require_all_tip_labels=(node_plot_mode == 'auto'),
            )
            if (node_plot_mode == 'auto') and (not all_tip_labels_parsed):
                sys.stderr.write(
                    'Skipping speciation/duplication node markers because some leaf labels did not match the configured species parser.\n'
                )
    leaf_label_color_by_leaf = dict()
    group_color_by_name = dict()
    if trait_path not in ['', None]:
        leaf_to_group = _read_trait_table(
            path=trait_path,
            group_by=group_by,
            tree=tree,
            unmatched=getattr(args, 'unmatched', 'warn'),
            missing_values=getattr(args, 'missing_values', None),
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
        figure_width=getattr(args, 'figure_width', FIGURE_WIDTH_IN),
        figure_height=getattr(args, 'figure_height', None),
        label_panel_width=getattr(args, 'label_panel_width', None),
        font_size=getattr(args, 'font_size', FONT_SIZE_PT),
        font_family=getattr(args, 'font_family', FONT_FAMILY),
        branch_color=getattr(args, 'branch_color', BRANCH_COLOR),
        terminal_branch_color=getattr(args, 'terminal_branch_color', None),
        tip_label_position=getattr(args, 'tip_label_position', 'aligned'),
        root_marker=getattr(args, 'root_marker', 'none'),
        root_marker_color=getattr(args, 'root_marker_color', '#0072B2'),
        root_marker_size=getattr(args, 'root_marker_size', None),
        tip_badge_property=getattr(args, 'tip_badge_property', None),
        tip_badge_missing_label=getattr(args, 'tip_badge_missing_label', None),
        node_pie_properties=node_pie_properties,
        node_pie_target=getattr(args, 'node_pie_target', 'root,intnode'),
        node_pie_leaf_filters=getattr(args, 'node_pie_leaf_filter', None),
        node_label_property=getattr(args, 'node_label_property', None),
        node_label_target=getattr(args, 'node_label_target', 'intnode'),
        node_label_filters=getattr(args, 'node_label_filter', None),
        node_label_decimals=getattr(args, 'node_label_decimals', 2),
        node_label_prefix=getattr(args, 'node_label_prefix', ''),
        property_colors=property_colors,
        legend=getattr(args, 'legend', True),
        transparent=getattr(args, 'transparent', False),
        trait_palette=trait_palette,
        tip_image_by_leaf=tip_image_by_leaf,
        tip_image_size=getattr(args, 'tip_image_size', TIP_IMAGE_SIZE_PT),
        tip_image_gap=getattr(args, 'tip_image_gap', TIP_IMAGE_GAP_PT),
    )
    sys.stderr.write('Wrote tree image: {}\n'.format(os.path.realpath(args.outfile)))
