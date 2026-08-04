import bisect
import hashlib
import math
import os
import sys

import matplotlib
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnnotationBbox, DrawingArea, OffsetImage
from matplotlib.path import Path as MplPath
from matplotlib.patches import Circle, Patch, Rectangle
from matplotlib.text import Text
from matplotlib.font_manager import FontProperties
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

from nwkit import __version__
from nwkit.draw_layouts import get_rectangular_coordinates, make_tree_layout
from nwkit.draw_prep import collapse_tree_for_drawing
from nwkit.draw_quality import DrawingArtist, evaluate_drawing, write_layout_report
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
MAX_TIP_IMAGE_SIZE_PT = 512.0
MAX_TIP_IMAGE_EDGE_PX = 2048
SUPPORT_LABEL_OFFSET_PT = 0.0
FIGURE_WIDTH_IN = 3.6
SPECIATION_COLOR = (0.0, 0.0, 1.0)
DUPLICATION_COLOR = (1.0, 0.0, 0.0)
BRANCH_COLOR = '#000000'
LABEL_COLOR = '#2d2d2d'
LEGEND_EDGE_COLOR = 'white'


def _resolve_image_format(outfile, image_format):
    if image_format != 'auto':
        return image_format
    ext = os.path.splitext(str(outfile))[1].lower()
    if ext in ('.pdf', '.png', '.svg'):
        return ext[1:]
    return 'pdf'


def _tree_drawing_fingerprint(tree):
    """Hash ordered topology and node properties without recursive Newick I/O."""

    digest = hashlib.sha256()

    def add(value):
        encoded = str(value).encode('utf-8', errors='backslashreplace')
        digest.update(str(len(encoded)).encode('ascii'))
        digest.update(b':')
        digest.update(encoded)

    stack = [('node', tree)]
    while stack:
        event, value = stack.pop()
        if event == 'close':
            digest.update(b')')
            continue
        node = value
        digest.update(b'(')
        for key in sorted(node.props):
            add(key)
            add(node.props[key])
        stack.append(('close', None))
        stack.extend(('node', child) for child in reversed(node.get_children()))
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
        duplicate_leaf_names='first',
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


def _load_tip_image(path, max_edge_px=None):
    from nwkit.image import (
        load_pillow_modules,
        rasterize_svg_to_image,
        validate_image_dimensions,
    )

    Image, _, ImageOps = load_pillow_modules()
    try:
        if os.path.splitext(path)[1].lower() == '.svg':
            image = rasterize_svg_to_image(path, max_edge=max_edge_px)
        else:
            with Image.open(path) as source_image:
                validate_image_dimensions(*source_image.size, label='Tip image')
                if max_edge_px not in (None, 0):
                    source_image.draft(
                        source_image.mode,
                        (int(max_edge_px), int(max_edge_px)),
                    )
                if max_edge_px not in (None, 0):
                    resampling = getattr(Image, 'Resampling', Image)
                    source_image.thumbnail(
                        (int(max_edge_px), int(max_edge_px)),
                        getattr(resampling, 'LANCZOS'),
                    )
                image = ImageOps.exif_transpose(source_image)
                image = image.convert('RGBA')
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
    image_by_path = dict()
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
        return []
    artists = []
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
        artists.append(ax.scatter(
            [x],
            [y],
            s=marker_size_pt ** 2,
            marker=MplPath(vertices, codes),
            facecolor=color,
            edgecolor='none',
            linewidth=0.0,
            zorder=6,
        ))
    artists.append(ax.scatter(
        [x],
        [y],
        s=marker_size_pt ** 2,
        marker='o',
        facecolor='none',
        edgecolor=LEGEND_EDGE_COLOR,
        linewidth=0.45,
        zorder=7,
    ))
    return artists


def _format_property_value(value, decimals=2):
    try:
        return ('{:.%df}' % int(decimals)).format(float(value))
    except (TypeError, ValueError):
        return str(value)


def _isolated_matplotlib_state(function):
    def wrapped(*args, **kwargs):
        existing_figures = set(plt.get_fignums())
        try:
            with matplotlib.rc_context({
                'font.family': [FONT_FAMILY],
                'font.sans-serif': [FONT_FAMILY],
                'svg.fonttype': 'none',
            }):
                return function(*args, **kwargs)
        finally:
            for figure_number in set(plt.get_fignums()) - existing_figures:
                plt.close(figure_number)
    return wrapped


def _validated_tip_image_size(value):
    try:
        image_size_pt = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError('--tip-image-size must be a finite number.') from exc
    if not math.isfinite(image_size_pt):
        raise ValueError('--tip-image-size must be a finite number.')
    if image_size_pt <= 0.0:
        raise ValueError('--tip-image-size must be greater than zero.')
    if image_size_pt > MAX_TIP_IMAGE_SIZE_PT:
        raise ValueError(
            '--tip-image-size must be no greater than {} points.'.format(
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
    return _readable_text_angle(normalized), ('left' if points_right else 'right')


def _parse_tip_label_wrap(value):
    text = str(value).strip().lower()
    if text in {'none', 'auto', 'taxonomy'}:
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
            for index, char in enumerate(remaining[:width + 1])
            if char.isspace() and index > 0
        ]
        delimiter_breaks = [
            index + 1
            for index, char in enumerate(remaining[:width])
            if char in '_-/|' and index > 0
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
    return '\n'.join(lines)


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

    parts = str(text).split('_')
    if len(parts) < 2:
        return None
    genus, species = parts[:2]
    if not genus or not species:
        return None
    if not genus[0].isupper() or not species[0].islower():
        return None
    return '{}_{}'.format(genus, species)


def _wrap_taxonomic_label(text, width):
    prefix = _taxonomic_prefix(text)
    if prefix is None:
        return _wrap_tip_label(text, width)
    suffix = str(text)[len(prefix):].lstrip('_')
    if not suffix:
        return str(text)
    wrapped_suffix = _wrap_tip_label(suffix, width)
    return '{}\n{}'.format(prefix, wrapped_suffix)


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
    raw_text_by_leaf = {leaf: str(leaf.name or '') for leaf in leaf_order}
    if wrap == 'none':
        text_by_leaf = dict(raw_text_by_leaf)
    elif isinstance(wrap, int):
        text_by_leaf = {
            leaf: _wrap_tip_label(text, wrap)
            for leaf, text in raw_text_by_leaf.items()
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
            'circular', 'fan', 'radial', 'unrooted', 'spiral', 'fractal',
        }
        normal_budget = max(min(panel_width_in, panel_height_in) * 0.20, 0.55)
        if spatial:
            perimeter_per_tip = (
                (2.0 * (panel_width_in + panel_height_in))
                / max(len(leaf_order), 1)
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
                line_penalty = 0.22 * candidate.count('\n')
                return congestion + line_penalty, width * height, candidate.count('\n')

            unwrapped = raw_text_by_leaf[leaf]
            if unwrapped not in candidate_sizes:
                candidate_sizes.update(_measure_texts_in_inches(
                    [unwrapped],
                    font_size=font_size,
                    font_family=font_family,
                ))
                candidates = [unwrapped] + candidates
            if wrap == 'taxonomy' and _taxonomic_prefix(unwrapped) is not None:
                semantic_candidates = [
                    candidate
                    for candidate in candidates
                    if '\n' in candidate or candidate == unwrapped
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
            'fits_within_figure': True,
            'overflow_left_points': 0.0,
            'overflow_right_points': 0.0,
            'overflow_bottom_points': 0.0,
            'overflow_top_points': 0.0,
            'maximum_overflow_points': 0.0,
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
        axes.set_position([new_left, new_bottom, new_right - new_left, new_top - new_bottom])
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    bounds = [artist.get_window_extent(renderer=renderer) for artist in artists]
    figure_bounds = figure.bbox
    overflow = {
        'overflow_left_points': max(
            (figure_bounds.x0 + padding_pixels) - min(bound.x0 for bound in bounds),
            0.0,
        ) * 72.0 / figure.dpi,
        'overflow_right_points': max(
            max(bound.x1 for bound in bounds) - (figure_bounds.x1 - padding_pixels),
            0.0,
        ) * 72.0 / figure.dpi,
        'overflow_bottom_points': max(
            (figure_bounds.y0 + padding_pixels) - min(bound.y0 for bound in bounds),
            0.0,
        ) * 72.0 / figure.dpi,
        'overflow_top_points': max(
            max(bound.y1 for bound in bounds) - (figure_bounds.y1 - padding_pixels),
            0.0,
        ) * 72.0 / figure.dpi,
    }
    overflow['maximum_overflow_points'] = max(overflow.values(), default=0.0)
    overflow['fits_within_figure'] = bool(
        overflow['maximum_overflow_points'] < 0.25
    )
    return overflow


def _parse_scale_bar(value, tree_span):
    text = str(value).strip().lower()
    if text in {'none', ''}:
        return None
    if text == 'auto':
        target = max(float(tree_span) / 5.0, 1e-12)
        exponent = math.floor(math.log10(target))
        scaled = target / (10.0 ** exponent)
        base = 1.0 if scaled < 1.5 else (2.0 if scaled < 3.5 else 5.0)
        return base * (10.0 ** exponent)
    try:
        result = float(text)
    except ValueError as error:
        raise ValueError("'--scale-bar' must be none, auto, or a positive number.") from error
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError("'--scale-bar' must be none, auto, or a positive number.")
    return result


def _parse_depth_guide(value, tree_span):
    text = str(value).strip().lower()
    if text in {'none', ''}:
        return None
    if text == 'auto':
        target = max(float(tree_span) / 4.0, 1e-12)
        exponent = math.floor(math.log10(target))
        scaled = target / (10.0 ** exponent)
        base = 1.0 if scaled < 1.5 else (2.0 if scaled < 3.5 else 5.0)
        return base * (10.0 ** exponent)
    try:
        interval = float(text)
    except ValueError as error:
        raise ValueError(
            "'--depth-guide' must be none, auto, or a positive interval."
        ) from error
    if not math.isfinite(interval) or interval <= 0.0:
        raise ValueError(
            "'--depth-guide' must be none, auto, or a positive interval."
        )
    if interval > float(tree_span) + 1e-12:
        raise ValueError(
            "'--depth-guide' interval must not exceed the displayed tree-depth "
            'span ({:g}).'.format(tree_span)
        )
    return interval


def _depth_guide_ticks(tree_span, interval):
    if interval is None:
        return []
    count = int(math.floor((float(tree_span) + 1e-12) / float(interval)))
    if count > 50:
        raise ValueError(
            "'--depth-guide' would draw {} intervals; use an interval of at "
            'least {:g}.'.format(count, float(tree_span) / 50.0)
        )
    return [float(interval) * index for index in range(count + 1)]


def _distance_label(prefix, unit):
    unit = str(unit).strip()
    return '{} ({})'.format(prefix, unit) if unit else prefix


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
        loc='lower left',
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
            '{:g}'.format(depth),
            xy=(x, 0.0),
            xycoords=('data', 'axes fraction'),
            xytext=(0.0, -3.0),
            textcoords='offset points',
            ha='center',
            va='top',
            fontsize=font_size,
            fontfamily=font_family,
            color=LABEL_COLOR,
            annotation_clip=False,
            zorder=9,
        )
        artists.append(label)
    return artists


def _add_radial_depth_guide(
    axes,
    drawing_layout,
    tree,
    ticks,
    color,
    font_family,
    font_size,
):
    """Draw concentric root-depth rings and label them in an empty sector."""

    root_x = drawing_layout.xcoord[tree]
    root_y = drawing_layout.ycoord[tree]
    leaf_angles = sorted(
        math.radians(drawing_layout.label_angles.get(leaf, 0.0))
        % (2.0 * math.pi)
        for leaf in tree.leaves()
    )
    segments = [
        (start, end)
        for path in drawing_layout.edge_paths.values()
        for start, end in zip(path, path[1:])
    ]

    def point_segment_distance(point, start, end):
        delta_x = end[0] - start[0]
        delta_y = end[1] - start[1]
        denominator = (delta_x * delta_x) + (delta_y * delta_y)
        if denominator <= 1e-18:
            return math.hypot(point[0] - start[0], point[1] - start[1])
        fraction = (
            ((point[0] - start[0]) * delta_x)
            + ((point[1] - start[1]) * delta_y)
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
            positive_ticks[round(
                index * (len(positive_ticks) - 1) / 11.0
            )]
            for index in range(12)
        ]
    if len(segments) <= 512:
        score_segments = segments
    else:
        score_segments = [
            segments[round(index * (len(segments) - 1) / 511.0)]
            for index in range(512)
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
            abs(math.atan2(
                math.sin(direction - angle),
                math.cos(direction - angle),
            ))
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

    candidate_count = 360 if len(segments) <= 512 else 180
    candidates = [
        2.0 * math.pi * index / candidate_count
        for index in range(candidate_count)
    ]
    direction = max(candidates, key=direction_score)
    best_score = direction_score(direction)
    labels_have_clear_sector = (
        len(segments) <= 2000
        and best_score >= outer_radius * 0.012
    )
    artists = []
    for depth in positive_ticks:
        axes.add_patch(Circle(
            (root_x, root_y),
            depth,
            facecolor='none',
            edgecolor=color,
            linewidth=0.45,
            linestyle=(0, (1.5, 2.5)),
            alpha=0.60,
            zorder=0.25,
        ))
        if not labels_have_clear_sector:
            continue
        x = root_x + (depth * math.cos(direction))
        y = root_y + (depth * math.sin(direction))
        label = axes.annotate(
            '{:g}'.format(depth),
            xy=(x, y),
            xytext=(3.0 * -math.sin(direction), 3.0 * math.cos(direction)),
            textcoords='offset points',
            ha='left' if math.cos(direction) >= 0.0 else 'right',
            va='center',
            fontsize=font_size,
            fontfamily=font_family,
            color=LABEL_COLOR,
            bbox={
                'boxstyle': 'square,pad=0.08',
                'facecolor': '#ffffff',
                'edgecolor': 'none',
                'alpha': 0.88,
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
        ha='left',
        va='bottom',
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
    line_y = (float(strip_height_points) - float(font_size) - 5.0) / figure_height_points
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
            '{:g}'.format(depth),
            ha='center',
            va='bottom',
            fontsize=font_size,
            fontfamily=font_family,
            color=LABEL_COLOR,
            zorder=10,
        )
        artists.append(label)
    return artists


def _parse_legend_columns(value, item_count):
    text = str(value).strip().lower()
    if text == 'auto':
        return max(1, min(4, int(math.ceil(max(item_count, 1) / 3.0))))
    try:
        columns = int(text)
    except ValueError as error:
        raise ValueError("'--legend-columns' must be auto or a positive integer.") from error
    if columns <= 0:
        raise ValueError("'--legend-columns' must be auto or a positive integer.")
    return columns


def _tip_font_style(text, mode):
    mode = str(mode).strip().lower()
    if mode == 'plain':
        return 'normal'
    if mode == 'italic':
        return 'italic'
    if mode == 'taxonomy':
        prefix = _taxonomic_prefix(str(text).replace('\n', ''))
        return 'italic' if prefix == str(text).replace('\n', '') else 'normal'
    raise ValueError("'--tip-label-font-style' must be plain, italic, or taxonomy.")


def _parse_branch_width_range(value):
    parts = [part.strip() for part in str(value).split(',')]
    if len(parts) != 2:
        raise ValueError("'--branch-width-range' must contain MIN,MAX.")
    try:
        minimum, maximum = (float(part) for part in parts)
    except ValueError as error:
        raise ValueError("'--branch-width-range' must contain MIN,MAX.") from error
    if minimum <= 0.0 or maximum < minimum or not all(map(math.isfinite, (minimum, maximum))):
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
    if color_property not in (None, ''):
        values = sorted({
            str(node.props[color_property])
            for node in nodes
            if node.props.get(color_property) not in (None, '')
        })
        value_index = {value: index for index, value in enumerate(values)}
        for node in nodes:
            value = node.props.get(color_property)
            if value not in (None, ''):
                color_by_node[node] = _property_color(
                    prop=color_property,
                    value=value,
                    property_colors=property_colors,
                    fallback_index=value_index[str(value)],
                    palette=palette,
                )
    if width_property not in (None, ''):
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
                width_by_node[node] = output_min + ((output_max - output_min) * fraction)
    return (
        {node: color_by_node.get(node, base_color) for node in nodes},
        {node: width_by_node.get(node, width) for node in nodes},
    )


def _tip_track_colors(
    tree,
    properties,
    mode,
    property_colors,
    categorical_palette,
    continuous_palette,
):
    mode = str(mode).strip().lower()
    if mode not in {'auto', 'categorical', 'continuous'}:
        raise ValueError("'--tip-track-type' must be auto, categorical, or continuous.")
    colors = {}
    legend_entries = []
    for property_name in properties:
        raw = {
            leaf: leaf.props.get(property_name)
            for leaf in tree.leaves()
        }
        present = [value for value in raw.values() if value not in (None, '')]
        numeric = {}
        for leaf, value in raw.items():
            try:
                number = float(value)
            except (TypeError, ValueError):
                continue
            if math.isfinite(number):
                numeric[leaf] = number
        continuous = mode == 'continuous' or (
            mode == 'auto' and present and len(numeric) == len(present)
        )
        if mode == 'continuous' and len(numeric) != len(present):
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
                    colors[(leaf, property_name)] = '#d9d9d9'
                    continue
                fraction = (
                    0.5
                    if maximum <= minimum
                    else (numeric[leaf] - minimum) / (maximum - minimum)
                )
                colors[(leaf, property_name)] = mcolors.to_hex(cmap(fraction))
            if maximum <= minimum:
                legend_entries.append((
                    '{}: {:g}'.format(property_name, minimum),
                    mcolors.to_hex(cmap(0.5)),
                ))
            else:
                legend_entries.extend([
                    ('{}: {:g}'.format(property_name, minimum), mcolors.to_hex(cmap(0.0))),
                    ('{}: {:g}'.format(property_name, maximum), mcolors.to_hex(cmap(1.0))),
                ])
        else:
            values = sorted({str(value) for value in present})
            for leaf, value in raw.items():
                if value in (None, ''):
                    colors[(leaf, property_name)] = '#d9d9d9'
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
                    '{}: {}'.format(property_name, value),
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
            legend_entries.append((
                '{}: missing'.format(property_name),
                '#d9d9d9',
            ))
    return colors, legend_entries


def _tip_track_artist(ax, xy, offset_points, color, size_points):
    drawing = DrawingArea(size_points, size_points, 0.0, 0.0)
    drawing.add_artist(Rectangle(
        (0.0, 0.0),
        size_points,
        size_points,
        facecolor=color,
        edgecolor='#ffffff',
        linewidth=0.35,
    ))
    artist = AnnotationBbox(
        drawing,
        xy,
        xybox=offset_points,
        xycoords='data',
        boxcoords='offset points',
        frameon=False,
        box_alignment=(0.5, 0.5),
        annotation_clip=False,
        zorder=8,
    )
    ax.add_artist(artist)
    return artist


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
    tip_label_position='aligned',
    tip_label_wrap='none',
    tip_spacing='uniform',
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
    layout='rectangular',
    spiral_turns=None,
    fan_span=180.0,
    unrooted_method='equal-angle',
    daylight_iterations=5,
    scale_bar='none',
    depth_guide='none',
    branch_length_unit='',
    branch_width=0.8,
    branch_color_property=None,
    branch_width_property=None,
    branch_width_range='0.4,2.5',
    tip_label_font_style='plain',
    tip_track_properties=None,
    tip_track_type='auto',
    tip_track_size=5.0,
    tip_track_palette='viridis',
    legend_columns='auto',
    legend_position='auto',
    collision_policy='resolve',
    layout_report=None,
    collapsed_clades=None,
):
    _apply_font_style(font_size=font_size, font_family=font_family)
    matplotlib.rcParams['svg.hashsalt'] = 'nwkit'
    property_colors = property_colors or {}
    node_pie_properties = node_pie_properties or []
    node_pie_leaf_filters = node_pie_leaf_filters or []
    node_label_filters = node_label_filters or []
    tip_image_by_leaf = tip_image_by_leaf or {}
    tip_track_properties = tip_track_properties or []
    collapsed_clades = collapsed_clades or []
    layout_name = str(layout).strip().lower()
    supported_layouts = {
        'rectangular',
        'slanted',
        'cladogram',
        'tidy',
        'circular',
        'fan',
        'radial',
        'unrooted',
        'spiral',
        'fractal',
    }
    if layout_name not in supported_layouts:
        raise ValueError("Unsupported '--layout': {}".format(layout))
    spatial_layout = layout_name in {
        'circular',
        'fan',
        'radial',
        'unrooted',
        'spiral',
        'fractal',
    }
    if spatial_layout and tip_image_by_leaf:
        raise ValueError(
            "'--tip-image-manifest' is currently supported only with "
            "rectangular, slanted, cladogram, and tidy layouts."
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
    _tip_font_style('', tip_label_font_style)
    resolved_tip_spacing = str(tip_spacing).strip().lower()
    if resolved_tip_spacing not in {'uniform', 'label-aware'}:
        raise ValueError("'--tip-spacing' must be uniform or label-aware.")
    daylight_iterations = int(daylight_iterations)
    if daylight_iterations <= 0:
        raise ValueError("'--daylight-iterations' must be greater than zero.")
    resolved_tip_label_position = str(tip_label_position).strip().lower()
    if resolved_tip_label_position == 'auto':
        resolved_tip_label_position = (
            'aligned' if layout_name == 'rectangular' else 'branch-end'
        )
    if spatial_layout and resolved_tip_label_position != 'branch-end':
        raise ValueError(
            "Two-dimensional layouts require "
            "'--tip-label-position branch-end'."
        )
    use_topology_depth = (not _has_positive_branch_length(tree))
    if use_topology_depth:
        sys.stderr.write(
            'Tree has no positive branch lengths; drawing positions by topology depth.\n'
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
        'rectangular', 'tidy', 'circular', 'fan', 'unrooted',
    }
    depth_projection_layouts = {'slanted', 'radial'}
    warped_depth_layouts = {'spiral'}
    depth_guide_layouts = depth_projection_layouts | warped_depth_layouts
    branch_length_layouts = (
        segment_length_layouts
        | depth_projection_layouts
        | warped_depth_layouts
    )
    branch_depth_span = max(base_x_max - base_x_min, 0.0)
    requested_scale_bar = _parse_scale_bar(scale_bar, branch_depth_span)
    requested_depth_guide = _parse_depth_guide(
        depth_guide,
        branch_depth_span,
    )
    if requested_scale_bar is not None and layout_name not in segment_length_layouts:
        raise ValueError(
            "'--scale-bar' requires a branch-length-preserving layout with "
            'directly measurable branch segments.'
        )
    if requested_scale_bar is not None and use_topology_depth:
        raise ValueError(
            "'--scale-bar' requires positive input branch lengths."
        )
    if (
        requested_scale_bar is not None
        and requested_scale_bar > branch_depth_span + 1e-12
    ):
        raise ValueError(
            "'--scale-bar' must not exceed the displayed tree-depth span "
            '({:g}).'.format(branch_depth_span)
        )
    if (
        requested_depth_guide is not None
        and layout_name not in depth_guide_layouts
    ):
        raise ValueError(
            "'--depth-guide' is supported by slanted, radial, and spiral "
            'layouts.'
        )
    if requested_depth_guide is not None and use_topology_depth:
        raise ValueError(
            "'--depth-guide' requires positive input branch lengths."
        )
    depth_ticks = _depth_guide_ticks(
        branch_depth_span,
        requested_depth_guide,
    )
    leaf_order = base_leaf_order
    fig_width = float(figure_width)
    if fig_width <= 0.0:
        raise ValueError('--figure-width must be greater than zero.')
    tip_image_size_pt = _validated_tip_image_size(tip_image_size)
    tip_image_gap_pt = float(tip_image_gap)
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
    label_height_guess_in = (
        float(figure_height)
        if figure_height is not None
        else (fig_width * 0.72 if layout_name == 'fractal' else fig_width)
    )
    if tip_labels:
        tip_label_text_by_leaf, tip_label_size_by_leaf = _prepare_tip_label_text(
            leaf_order=leaf_order,
            wrap=tip_label_wrap,
            font_size=font_size,
            font_family=font_family,
            layout_name=layout_name,
            panel_width_in=main_panel_width_in,
            panel_height_in=max(label_height_guess_in, 0.2),
        )
    else:
        _parse_tip_label_wrap(tip_label_wrap)
        tip_label_text_by_leaf = {leaf: str(leaf.name or '') for leaf in leaf_order}
        tip_label_size_by_leaf = {leaf: (0.0, 0.0) for leaf in leaf_order}
    max_tip_label_width_in = max(
        (size[0] for size in tip_label_size_by_leaf.values()),
        default=0.0,
    )
    max_tip_label_height_in = max(
        (size[1] for size in tip_label_size_by_leaf.values()),
        default=(float(font_size) / 72.0 if tip_labels else 0.0),
    )
    tip_track_color_by_leaf_property, tip_track_legend_entries = _tip_track_colors(
        tree=tree,
        properties=tip_track_properties,
        mode=tip_track_type,
        property_colors=property_colors,
        categorical_palette=trait_palette,
        continuous_palette=tip_track_palette,
    )
    badge_values = []
    if tip_badge_property not in (None, ''):
        badge_values = sorted({
            str(
                tip_badge_missing_label
                if leaf.props.get(tip_badge_property) in (None, '')
                else leaf.props.get(tip_badge_property)
            )
            for leaf in tree.leaves()
            if (
                leaf.props.get(tip_badge_property) not in (None, '')
                or tip_badge_missing_label not in (None, '')
            )
        })
    badge_value_index = {
        value: index for index, value in enumerate(badge_values)
    }
    if spatial_layout:
        label_panel_width_in = 0.0
    elif (
        not tip_labels
        and tip_badge_property in (None, '')
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
                '--label-panel-width must be greater than zero and leave room '
                'for the tree and tip-image column within --figure-width.'
            )
    tree_panel_width_in = main_panel_width_in - label_panel_width_in
    if tree_panel_width_in < 0.2:
        raise ValueError(
            '--figure-width is too small for the requested label and tip-image panels.'
        )
    spacing_height_by_leaf = {}
    for leaf in leaf_order:
        spacing_height = (
            tip_label_size_by_leaf[leaf][1]
            if tip_labels
            else 0.0
        )
        if tip_track_properties:
            spacing_height = max(spacing_height, track_size_pt / 72.0)
        if tip_badge_property not in (None, ''):
            spacing_height = max(
                spacing_height,
                float(font_size) * 1.3 / 72.0,
            )
        if str(leaf.name) in tip_image_by_leaf:
            spacing_height = max(spacing_height, tip_image_size_pt / 72.0)
        spacing_height_by_leaf[leaf] = max(
            spacing_height,
            float(font_size) / 72.0,
        )
    if resolved_tip_spacing == 'label-aware' and spacing_height_by_leaf:
        row_pitch_pt = (
            sum(spacing_height_by_leaf.values())
            / len(spacing_height_by_leaf)
        ) * 72.0 + TIP_LABEL_GAP_PT
    else:
        row_pitch_pt = max(
            float(font_size),
            max_tip_label_height_in * 72.0,
        ) + TIP_LABEL_GAP_PT
        if tip_image_by_leaf:
            row_pitch_pt = max(
                row_pitch_pt,
                tip_image_size_pt + TIP_LABEL_GAP_PT,
            )
    row_pitch_in = row_pitch_pt / 72.0
    property_legend_needed = bool(badge_values or node_pie_properties)
    legend_needed = bool(legend) and bool(
        node_type_by_node
        or group_color_by_name
        or property_legend_needed
        or tip_track_legend_entries
        or branch_color_property not in (None, '')
        or branch_width_property not in (None, '')
    )
    top_margin_in = (44.0 / 72.0) if legend_needed else (4.0 / 72.0)
    scale_bar_strip_pt = (
        max(float(font_size) * 1.4 + 8.0, 18.0)
        if requested_scale_bar is not None
        else 0.0
    )
    depth_guide_strip_pt = (
        max(float(font_size) * 2.6 + 14.0, 34.0)
        if requested_depth_guide is not None
        else 0.0
    )
    bottom_margin_in = (
        3.0 + max(scale_bar_strip_pt, depth_guide_strip_pt)
    ) / 72.0
    left_margin_in = 2.0 / 72.0
    right_margin_in = 2.0 / 72.0
    if figure_height is None and spatial_layout:
        fig_height = fig_width * 0.72 if layout_name == 'fractal' else fig_width
    elif figure_height is None:
        # Refined after tidy coordinates have been calculated.
        fig_height = max(len(leaf_order), 1) * row_pitch_in + top_margin_in + bottom_margin_in
    else:
        fig_height = float(figure_height)
        if fig_height <= top_margin_in + bottom_margin_in:
            raise ValueError('--figure-height is too small for the required margins.')

    available_height_in = max(fig_height - top_margin_in - bottom_margin_in, 0.2)
    terminal_extent_by_leaf = {}
    label_data_per_inch = base_x_span / max(tree_panel_width_in, 0.2)
    for leaf in leaf_order:
        label_width_in = tip_label_size_by_leaf.get(leaf, (0.0, 0.0))[0]
        label_extent = (
            label_width_in + (track_span_pt / 72.0)
        ) * label_data_per_inch
        if tip_labels and resolved_tip_label_position == 'aligned':
            label_extent += max(base_x_max - base_xcoord[leaf], 0.0)
        terminal_extent_by_leaf[leaf] = label_extent

    layout_aspect_ratio = tree_panel_width_in / available_height_in
    if (
        layout_name == 'fractal'
        and resolved_tip_spacing == 'label-aware'
        and tip_labels
    ):
        normal_label_margin = max_tip_label_width_in + (5.0 / 72.0)
        label_aware_width = max(
            tree_panel_width_in - (2.0 * normal_label_margin),
            0.2,
        )
        label_aware_height = max(
            available_height_in - (2.0 * normal_label_margin),
            0.2,
        )
        layout_aspect_ratio = label_aware_width / label_aware_height
    spacing_size_by_leaf = {
        leaf: (
            tip_label_size_by_leaf[leaf][0],
            spacing_height_by_leaf[leaf],
        )
        for leaf in leaf_order
    }
    drawing_layout = make_tree_layout(
        tree,
        layout=layout_name,
        use_topology_depth=use_topology_depth,
        aspect_ratio=layout_aspect_ratio,
        spiral_turns=spiral_turns,
        fan_span=fan_span,
        terminal_extent_by_leaf=terminal_extent_by_leaf,
        label_size_by_leaf=spacing_size_by_leaf,
        tip_spacing=resolved_tip_spacing,
        unrooted_method=unrooted_method,
        daylight_iterations=daylight_iterations,
    )
    xcoord = drawing_layout.xcoord
    ycoord = drawing_layout.ycoord
    leaf_order = drawing_layout.leaf_order
    x_min, x_max, y_min, y_max = drawing_layout.bounds
    if (
        requested_depth_guide is not None
        and layout_name == 'radial'
        and depth_ticks
    ):
        guide_radius = max(depth_ticks)
        root_x = drawing_layout.xcoord[tree]
        root_y = drawing_layout.ycoord[tree]
        x_min = min(x_min, root_x - guide_radius)
        x_max = max(x_max, root_x + guide_radius)
        y_min = min(y_min, root_y - guide_radius)
        y_max = max(y_max, root_y + guide_radius)
    x_span = max(x_max - x_min, 1e-9)
    y_span = max(y_max - y_min, 1e-9)

    if not spatial_layout:
        tree_panel_height_in = max(y_span + 1.0, 1.0) * row_pitch_in
        if figure_height is None:
            fig_height = tree_panel_height_in + top_margin_in + bottom_margin_in
        if (
            tip_image_by_leaf
            and fig_height + 10 ** -9
            < tree_panel_height_in + top_margin_in + bottom_margin_in
        ):
            raise ValueError(
                '--figure-height is too small for non-overlapping tip images. '
                'Increase it or reduce --tip-image-size.'
            )
    else:
        tree_panel_height_in = max(fig_height - top_margin_in - bottom_margin_in, 0.2)

    if layout_name in {'cladogram', 'fractal'} and not use_topology_depth:
        sys.stderr.write(
            '{} layout encodes topology; input branch lengths do not determine geometry.\n'.format(
                layout_name.capitalize(),
            )
        )

    color_by_node, width_by_node = _branch_style_maps(
        tree=tree,
        base_color=branch_color,
        base_width=branch_width,
        color_property=branch_color_property,
        width_property=branch_width_property,
        width_range=branch_width_range,
        property_colors=property_colors,
        palette=trait_palette,
    )
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    drawing_artists = []
    branch_lines = []
    depth_guide_color = '#8c8c8c'
    radial_depth_labels_drawn = None
    if requested_depth_guide is not None and layout_name == 'slanted':
        for artist in _add_slanted_depth_guide(
            axes=ax,
            root_x=xcoord[tree],
            ticks=depth_ticks,
            color=depth_guide_color,
            font_family=font_family,
            font_size=font_size,
        ):
            drawing_artists.append(DrawingArtist(
                artist,
                kind='depth_guide',
                priority=100,
            ))
    elif requested_depth_guide is not None and layout_name == 'radial':
        radial_artists, radial_depth_labels_drawn = _add_radial_depth_guide(
            axes=ax,
            drawing_layout=drawing_layout,
            tree=tree,
            ticks=depth_ticks,
            color=depth_guide_color,
            font_family=font_family,
            font_size=font_size,
        )
        for artist in radial_artists:
            drawing_artists.append(DrawingArtist(
                artist,
                kind='depth_guide',
                priority=100,
            ))

    for node, path in drawing_layout.edge_paths.items():
        if branch_color_property not in (None, '') and node.props.get(branch_color_property) not in (None, ''):
            edge_color = color_by_node[node]
        elif node.is_leaf and terminal_branch_color:
            edge_color = terminal_branch_color
        else:
            edge_color = color_by_node[node]
        path_x, path_y = zip(*path)
        line = ax.plot(
            path_x,
            path_y,
            color=edge_color,
            linewidth=width_by_node[node],
            zorder=2,
            solid_capstyle=TREE_LINE_CAPSTYLE,
            solid_joinstyle='round',
        )[0]
        branch_lines.append((node, line))

    if drawing_layout.root_path:
        root_x, root_y = zip(*drawing_layout.root_path)
        root_line = ax.plot(
            root_x,
            root_y,
            color=branch_color,
            linewidth=float(branch_width),
            zorder=2,
            solid_capstyle=TREE_LINE_CAPSTYLE,
            solid_joinstyle='round',
        )[0]
        branch_lines.append((tree, root_line))

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
    if node_pie_properties:
        legend_handles.extend(
            Patch(
                facecolor=_property_color(
                    prop=prop,
                    value=prop,
                    property_colors=property_colors,
                    fallback_index=index,
                    palette=trait_palette,
                ),
                edgecolor='none',
                label=prop,
            )
            for index, prop in enumerate(node_pie_properties)
        )
    if badge_values:
        legend_handles.extend(
            Patch(
                facecolor=_property_color(
                    prop=tip_badge_property,
                    value=value,
                    property_colors=property_colors,
                    fallback_index=badge_value_index[value],
                    palette=trait_palette,
                ),
                edgecolor='none',
                label='{}: {}'.format(tip_badge_property, value),
            )
            for value in badge_values
        )
    legend_handles.extend(
        Patch(facecolor=color, edgecolor='none', label=label)
        for label, color in tip_track_legend_entries
    )
    if branch_color_property not in (None, ''):
        branch_values = sorted({
            str(node.props[branch_color_property])
            for node in tree.traverse()
            if not node.is_root and node.props.get(branch_color_property) not in (None, '')
        })
        branch_value_index = {
            value: index for index, value in enumerate(branch_values)
        }
        legend_handles.extend(
            Line2D(
                [0], [0],
                color=_property_color(
                    prop=branch_color_property,
                    value=value,
                    property_colors=property_colors,
                    fallback_index=branch_value_index[value],
                    palette=trait_palette,
                ),
                linewidth=float(branch_width),
                label='{}: {}'.format(branch_color_property, value),
            )
            for value in branch_values
        )
        if any(
            node.props.get(branch_color_property) in (None, '')
            for node in tree.traverse()
            if not node.is_root
        ):
            legend_handles.append(Line2D(
                [0], [0],
                color=branch_color,
                linewidth=float(branch_width),
                label='{}: missing'.format(branch_color_property),
            ))
    if branch_width_property not in (None, ''):
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
                    [0], [0],
                    color=branch_color,
                    linewidth=width_by_node[node],
                    label='{}: {:g}'.format(branch_width_property, value),
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
    if legend and len(legend_handles) > 0:
        resolved_legend_position = str(legend_position).strip().lower()
        if resolved_legend_position == 'auto':
            resolved_legend_position = (
                'right'
                if len(legend_handles) > 12
                else 'top'
            )
        if resolved_legend_position not in {'top', 'right'}:
            raise ValueError("'--legend-position' must be auto, top, or right.")
        if resolved_legend_position == 'right':
            legend_location = 'upper left'
            legend_anchor = (1.005, 1.0)
            legend_ncol = 1
        else:
            legend_location = 'lower right'
            legend_anchor = (1.0, 1.005)
            legend_ncol = _parse_legend_columns(
                legend_columns,
                len(legend_handles),
            )
        legend_artist = ax.legend(
            handles=legend_handles,
            loc=legend_location,
            bbox_to_anchor=legend_anchor,
            frameon=False,
            borderaxespad=0.1,
            handletextpad=0.3,
            labelspacing=0.2,
            ncol=legend_ncol,
            prop={'family': font_family, 'size': font_size},
        )
        drawing_artists.append(DrawingArtist(
            legend_artist,
            kind='legend',
            priority=100,
        ))

    for node in tree.traverse():
        if not support_labels:
            continue
        if not _has_meaningful_support(node):
            continue
        support_value = float(node.support)
        if (support_min is not None) and (support_value < support_min):
            continue
        label_x, label_y = drawing_layout.support_anchors[node]
        edge_angle = drawing_layout.support_angles.get(node, 0.0)
        if spatial_layout:
            radians = math.radians(edge_angle)
            support_offset = (
                -math.sin(radians) * 1.5,
                math.cos(radians) * 1.5,
            )
        else:
            support_offset = (0.0, SUPPORT_LABEL_OFFSET_PT)
        support_artist = ax.annotate(
            _format_support_value(node.support),
            xy=(label_x, label_y),
            xytext=support_offset,
            textcoords='offset points',
            va='bottom',
            ha='center',
            rotation=_readable_text_angle(edge_angle) if spatial_layout else 0.0,
            rotation_mode='anchor',
            fontsize=font_size,
            fontfamily=font_family,
            color=LABEL_COLOR,
        )
        drawing_artists.append(DrawingArtist(
            support_artist,
            kind='support_label',
            priority=55,
            movable=True,
            shift_direction=(
                -math.sin(math.radians(edge_angle)),
                math.cos(math.radians(edge_angle)),
            ) if spatial_layout else (0.0, 1.0),
            owner=node,
        ))

    data_per_inch = x_span / max(tree_panel_width_in, 0.2)
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
    tip_label_artists = []
    track_span_data = (track_span_pt / 72.0) * data_per_inch
    for leaf in leaf_order:
        label_angle = drawing_layout.label_angles.get(leaf, 0.0)
        if spatial_layout:
            text_rotation, text_alignment = _spatial_text_alignment(label_angle)
            radians = math.radians(label_angle)
            for track_index, property_name in enumerate(tip_track_properties):
                track_distance = (
                    3.0
                    + (track_size_pt / 2.0)
                    + track_index * track_stride_pt
                )
                track_artist = _tip_track_artist(
                    ax=ax,
                    xy=(xcoord[leaf], ycoord[leaf]),
                    offset_points=(
                        math.cos(radians) * track_distance,
                        math.sin(radians) * track_distance,
                    ),
                    color=tip_track_color_by_leaf_property[(leaf, property_name)],
                    size_points=track_size_pt,
                )
                drawing_artists.append(DrawingArtist(
                    track_artist,
                    kind='tip_track',
                    priority=85,
                    owner=leaf,
                ))
            label_offset_points = (
                math.cos(radians) * (4.0 + track_span_pt),
                math.sin(radians) * (4.0 + track_span_pt),
            )
            label_x = xcoord[leaf]
            label_y = ycoord[leaf]
            if tip_labels:
                tip_artist = ax.annotate(
                    tip_label_text_by_leaf[leaf],
                    xy=(label_x, label_y),
                    xytext=label_offset_points,
                    textcoords='offset points',
                    va='center',
                    ha=text_alignment,
                    rotation=text_rotation,
                    rotation_mode='anchor',
                    fontsize=font_size,
                    fontfamily=font_family,
                    fontstyle=_tip_font_style(leaf.name or '', tip_label_font_style),
                    linespacing=1.15,
                    color=leaf_label_color_by_leaf.get(leaf, LABEL_COLOR),
                    annotation_clip=False,
                )
                tip_label_artists.append(tip_artist)
                drawing_artists.append(DrawingArtist(
                    tip_artist,
                    kind='tip_label',
                    priority=95,
                    owner=leaf,
                ))
        elif resolved_tip_label_position == 'branch-end':
            track_origin_x = xcoord[leaf] + label_offset
            for track_index, property_name in enumerate(tip_track_properties):
                track_artist = _tip_track_artist(
                    ax=ax,
                    xy=(track_origin_x, ycoord[leaf]),
                    offset_points=(
                        (track_size_pt / 2.0)
                        + track_index * track_stride_pt,
                        0.0,
                    ),
                    color=tip_track_color_by_leaf_property[(leaf, property_name)],
                    size_points=track_size_pt,
                )
                drawing_artists.append(DrawingArtist(
                    track_artist,
                    kind='tip_track',
                    priority=85,
                    owner=leaf,
                ))
            label_x = track_origin_x + track_span_data + ((1.0 / 72.0) * data_per_inch)
            label_y = ycoord[leaf]
            if tip_labels:
                tip_artist = ax.text(
                    label_x,
                    label_y,
                    tip_label_text_by_leaf[leaf],
                    va='center',
                    ha='left',
                    fontsize=font_size,
                    fontfamily=font_family,
                    fontstyle=_tip_font_style(leaf.name or '', tip_label_font_style),
                    linespacing=1.15,
                    color=leaf_label_color_by_leaf.get(leaf, LABEL_COLOR),
                )
                tip_label_artists.append(tip_artist)
                drawing_artists.append(DrawingArtist(
                    tip_artist,
                    kind='tip_label',
                    priority=95,
                    owner=leaf,
                ))
        elif resolved_tip_label_position == 'aligned':
            track_origin_x = x_max + label_offset
            for track_index, property_name in enumerate(tip_track_properties):
                track_artist = _tip_track_artist(
                    ax=ax,
                    xy=(track_origin_x, ycoord[leaf]),
                    offset_points=(
                        (track_size_pt / 2.0)
                        + track_index * track_stride_pt,
                        0.0,
                    ),
                    color=tip_track_color_by_leaf_property[(leaf, property_name)],
                    size_points=track_size_pt,
                )
                drawing_artists.append(DrawingArtist(
                    track_artist,
                    kind='tip_track',
                    priority=85,
                    owner=leaf,
                ))
            label_x = track_origin_x + track_span_data + ((1.0 / 72.0) * data_per_inch)
            label_y = ycoord[leaf]
            if tip_labels:
                tip_artist = ax.text(
                    label_x,
                    label_y,
                    tip_label_text_by_leaf[leaf],
                    va='center',
                    ha='left',
                    fontsize=font_size,
                    fontfamily=font_family,
                    fontstyle=_tip_font_style(leaf.name or '', tip_label_font_style),
                    linespacing=1.15,
                    color=leaf_label_color_by_leaf.get(leaf, LABEL_COLOR),
                )
                tip_label_artists.append(tip_artist)
                drawing_artists.append(DrawingArtist(
                    tip_artist,
                    kind='tip_label',
                    priority=95,
                    owner=leaf,
                ))
        else:
            raise ValueError(
                "Unsupported '--tip-label-position': {}".format(tip_label_position)
            )
        if str(leaf.props.get('nwkit_collapsed', '')).lower() == 'true':
            collapsed_marker = (
                (3, 0, label_angle - 90.0)
                if spatial_layout
                else '>'
            )
            ax.scatter(
                [xcoord[leaf]],
                [ycoord[leaf]],
                s=max(float(font_size), 5.0) ** 2,
                marker=collapsed_marker,
                facecolor=branch_color,
                edgecolor=LEGEND_EDGE_COLOR,
                linewidth=0.4,
                zorder=7,
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
                fallback_index=badge_value_index[value],
                palette=trait_palette,
            )
            badge_style = {
                'boxstyle': 'round,pad=0.16',
                'facecolor': badge_color,
                'edgecolor': LEGEND_EDGE_COLOR,
                'linewidth': 0.4,
                'alpha': 0.92,
            }
            if spatial_layout:
                text_rotation, _ = _spatial_text_alignment(label_angle)
                radians = math.radians(label_angle)
                badge_offset_points = (
                    track_span_pt
                    + (tip_label_size_by_leaf[leaf][0] * 72.0 if tip_labels else 0.0)
                    + 4.0
                    + max(len(value), 1) * font_size * 0.29
                    + (font_size * 0.16)
                )
                badge_artist = ax.annotate(
                    value,
                    xy=(label_x, label_y),
                    xytext=(
                        math.cos(radians) * badge_offset_points,
                        math.sin(radians) * badge_offset_points,
                    ),
                    textcoords='offset points',
                    va='center',
                    ha='center',
                    rotation=text_rotation,
                    rotation_mode='anchor',
                    fontsize=font_size,
                    fontfamily=font_family,
                    fontweight='bold',
                    color=LABEL_COLOR,
                    bbox=badge_style,
                    zorder=8,
                    annotation_clip=False,
                )
                drawing_artists.append(DrawingArtist(
                    badge_artist,
                    kind='tip_badge',
                    priority=70,
                    movable=True,
                    shift_direction=(-math.sin(radians), math.cos(radians)),
                    owner=leaf,
                ))
            else:
                approximate_label_width = (
                    tip_label_size_by_leaf[leaf][0] * data_per_inch
                    if tip_labels
                    else 0.0
                )
                approximate_badge_width = max(len(value), 1) * (font_size * 0.58 / 72.0) * data_per_inch
                badge_padding = (font_size * 0.16 / 72.0) * data_per_inch
                badge_gap = (4.0 / 72.0) * data_per_inch
                badge_offset = (
                    approximate_label_width
                    + badge_gap
                    + (approximate_badge_width / 2.0)
                    + badge_padding
                )
                badge_artist = ax.text(
                    label_x + badge_offset,
                    label_y,
                    value,
                    va='center',
                    ha='center',
                    fontsize=font_size,
                    fontfamily=font_family,
                    fontweight='bold',
                    color=LABEL_COLOR,
                    bbox=badge_style,
                    zorder=8,
                )
                drawing_artists.append(DrawingArtist(
                    badge_artist,
                    kind='tip_badge',
                    priority=70,
                    owner=leaf,
                ))

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
            node_label_artist = ax.annotate(
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
            drawing_artists.append(DrawingArtist(
                node_label_artist,
                kind='node_label',
                priority=45,
                movable=True,
                shift_direction=(0.0, 1.0),
                owner=node,
            ))

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

    label_panel_span = x_span * (
        label_panel_width_in / max(tree_panel_width_in, 0.2)
    )
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
    if spatial_layout:
        horizontal_label_extent_in = 0.0
        vertical_label_extent_in = 0.0
        if tip_labels:
            for leaf in leaf_order:
                width_in, height_in = tip_label_size_by_leaf[leaf]
                angle = math.radians(drawing_layout.label_angles.get(leaf, 0.0))
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
        if tip_badge_property not in (None, ''):
            badge_extent = max(0.3, float(font_size) * 2.0 / 72.0)
            horizontal_label_extent_in += badge_extent
            vertical_label_extent_in += badge_extent
        if tip_track_properties:
            horizontal_label_extent_in += track_span_pt / 72.0
            vertical_label_extent_in += track_span_pt / 72.0
        horizontal_label_margin = horizontal_label_extent_in * (
            x_span / max(tree_panel_width_in, 0.2)
        )
        vertical_label_margin = vertical_label_extent_in * (
            y_span / max(tree_panel_height_in, 0.2)
        )
        plot_pad_x = max(x_span * 0.035, horizontal_label_margin)
        plot_pad_y = max(y_span * 0.035, vertical_label_margin)
        ax.set_xlim(x_min - plot_pad_x, x_max + plot_pad_x)
        ax.set_ylim(y_min - plot_pad_y, y_max + plot_pad_y)
        if drawing_layout.equal_aspect:
            ax.set_aspect('equal', adjustable='box')
    else:
        root_path_span = 0.0
        if drawing_layout.root_path:
            root_path_span = abs(
                drawing_layout.root_path[-1][0] - drawing_layout.root_path[0][0]
            )
        left_plot_margin = max(root_path_span * 1.2, root_symbol_margin)
        ax.set_xlim(
            x_min - left_plot_margin,
            x_max + label_offset + label_panel_span,
        )
        if len(leaf_order) == 0:
            ax.set_ylim(0.5, -0.5)
        else:
            ax.set_ylim(y_max + 0.5, y_min - 0.5)
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
    if requested_scale_bar is not None:
        scale_label = '{:g}'.format(requested_scale_bar)
        if str(branch_length_unit).strip():
            scale_label = '{} {}'.format(
                scale_label,
                str(branch_length_unit).strip(),
            )
        scale_artist = _add_scale_bar(
            figure=fig,
            axes=ax,
            size=requested_scale_bar,
            label=scale_label,
            color=branch_color,
            font_family=font_family,
            font_size=font_size,
            anchor_x=axes_left,
            anchor_y=(2.0 / 72.0) / fig_height,
        )
        drawing_artists.append(DrawingArtist(
            scale_artist,
            kind='scale_bar',
            priority=100,
        ))
    if requested_depth_guide is not None:
        guide_prefix = {
            'slanted': 'Root-to-node distance',
            'radial': (
                'Concentric rings every {:g}: root-to-node distance; root = 0'
            ).format(requested_depth_guide),
            'spiral': 'Root-to-node distance encoded across spiral band',
        }[layout_name]
        title_artist = _add_bottom_guide_title(
            figure=fig,
            text=_distance_label(guide_prefix, branch_length_unit),
            x=axes_left,
            y=(2.0 / 72.0) / fig_height,
            font_family=font_family,
            font_size=font_size,
        )
        drawing_artists.append(DrawingArtist(
            title_artist,
            kind='depth_guide_title',
            priority=100,
        ))
        if layout_name == 'spiral':
            for artist in _add_spiral_depth_key(
                figure=fig,
                ticks=depth_ticks,
                tree_span=branch_depth_span,
                axes_left=axes_left,
                axes_right=axes_right,
                strip_height_points=depth_guide_strip_pt,
                font_family=font_family,
                font_size=font_size,
                color=depth_guide_color,
            ):
                drawing_artists.append(DrawingArtist(
                    artist,
                    kind='depth_guide',
                    priority=100,
                ))
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
    fit_report = _fit_artists_within_figure(
        figure=fig,
        axes=ax,
        artists=[item.artist for item in drawing_artists],
    )
    quality_report = evaluate_drawing(
        figure=fig,
        artists=drawing_artists,
        branch_lines=branch_lines,
        policy=collision_policy,
    )
    final_fit_report = _fit_artists_within_figure(
        figure=fig,
        axes=ax,
        artists=[item.artist for item in drawing_artists],
    )
    if not final_fit_report['fits_within_figure']:
        overflow_message = (
            'Drawing annotations exceed the fixed figure boundary by up to '
            '{:.2f} point(s). Increase the figure size or enable label wrapping.'
        ).format(final_fit_report['maximum_overflow_points'])
        if str(collision_policy).strip().lower() == 'error':
            raise ValueError(overflow_message)
        if str(collision_policy).strip().lower() in {'resolve', 'warn'}:
            sys.stderr.write(overflow_message + '\n')
    final_quality = evaluate_drawing(
        figure=fig,
        artists=drawing_artists,
        branch_lines=branch_lines,
        policy='ignore',
        max_iterations=0,
    )
    quality_report['final_collision_count'] = final_quality['final_collision_count']
    quality_report['final_collisions_by_kind'] = final_quality['final_collisions_by_kind']
    quality_report['annotation_occupied_fraction'] = final_quality[
        'annotation_occupied_fraction'
    ]
    quality_report['branch_collision_check_complete'] = final_quality[
        'branch_collision_check_complete'
    ]
    quality_report.update(final_fit_report)
    quality_report.update({
        'nwkit_version': __version__,
        'output_format': str(image_format),
        'layout_requested': layout_name,
        'layout': drawing_layout.name,
        'unrooted_method': str(unrooted_method),
        'daylight_iterations_requested': daylight_iterations,
        **drawing_layout.metadata,
        'figure_width_inches': fig_width,
        'figure_height_inches': fig_height,
        'font_family': str(font_family),
        'font_size_points': float(font_size),
        'input_tip_count': sum(
            int(item.get('tip_count', 1)) for item in collapsed_clades
        ) + len(leaf_order) - len(collapsed_clades),
        'visible_tip_count': len(leaf_order),
        'display_tree_sha256': _tree_drawing_fingerprint(tree),
        'collapsed_clades': collapsed_clades,
        'wrapped_tip_count': sum(
            '\n' in tip_label_text_by_leaf[leaf] for leaf in leaf_order
        ),
        'branch_lengths_encoded': (
            layout_name in branch_length_layouts and not use_topology_depth
        ),
        'branch_length_encoding': (
            'segment'
            if layout_name in segment_length_layouts and not use_topology_depth
            else (
                'depth-projection'
                if layout_name in depth_projection_layouts and not use_topology_depth
                else (
                    'warped-depth'
                    if layout_name in warped_depth_layouts and not use_topology_depth
                    else 'none'
                )
            )
        ),
        'scale_bar': requested_scale_bar,
        'scale_bar_position': (
            'bottom-reserved' if requested_scale_bar is not None else None
        ),
        'scale_bar_label_position': (
            'above' if requested_scale_bar is not None else None
        ),
        'depth_guide_interval': requested_depth_guide,
        'depth_guide_type': (
            {
                'slanted': 'axis-grid',
                'radial': 'concentric-rings',
                'spiral': 'spiral-depth-key',
            }[layout_name]
            if requested_depth_guide is not None
            else None
        ),
        'depth_guide_in_panel_labels': (
            radial_depth_labels_drawn
            if requested_depth_guide is not None and layout_name == 'radial'
            else None
        ),
        'branch_length_unit': str(branch_length_unit),
        'tip_track_properties': list(tip_track_properties),
        'tip_track_palette': str(tip_track_palette),
        'tip_spacing': resolved_tip_spacing,
        'initial_fit': fit_report,
    })
    metadata = {'Creator': 'NWKIT {}'.format(__version__)}
    if image_format == 'svg':
        metadata['Date'] = None
    elif image_format == 'pdf':
        metadata['CreationDate'] = None
        metadata['ModDate'] = None
    try:
        fig.savefig(
            outfile,
            format=image_format,
            dpi=300,
            transparent=bool(transparent),
            metadata=metadata,
        )
        write_layout_report(layout_report, quality_report)
    finally:
        plt.close(fig)


def _canonical_path(path):
    if path in (None, '', '-'):
        return None
    return os.path.normcase(os.path.realpath(os.path.expanduser(str(path))))


def _validate_draw_paths(outfile, layout_report, input_paths):
    outputs = {
        '--outfile': _canonical_path(outfile),
        '--layout-report': _canonical_path(layout_report),
    }
    if outputs['--outfile'] == outputs['--layout-report'] and outputs['--outfile']:
        raise ValueError("'--outfile' and '--layout-report' must use different paths.")
    inputs = {
        option: _canonical_path(path)
        for option, path in input_paths.items()
        if _canonical_path(path) is not None
    }
    for output_option, output_path in outputs.items():
        if output_path is None:
            continue
        for input_option, input_path in inputs.items():
            if output_path == input_path:
                raise ValueError(
                    '{} must not overwrite the {} input.'.format(
                        output_option,
                        input_option,
                    )
                )


def draw_main(args):
    if args.outfile == '-':
        raise ValueError("STDOUT is not supported for 'draw'. Use --outfile PATH.")
    trait_path = getattr(args, 'trait', None)
    tip_image_manifest = getattr(args, 'tip_image_manifest', None)
    layout_report = getattr(args, 'layout_report', None)
    _validate_draw_paths(
        outfile=args.outfile,
        layout_report=layout_report,
        input_paths={
            '--infile': args.infile,
            '--trait': trait_path,
            '--tip-image-manifest': tip_image_manifest,
        },
    )
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    if bool(args.ladderize):
        tree.ladderize()
    max_visible_tips = getattr(args, 'max_visible_tips', None)
    if max_visible_tips not in (None, '') and int(max_visible_tips) < 2:
        raise ValueError("'--max-visible-tips' must be at least 2.")
    collapse_required = (
        max_visible_tips not in (None, '')
        and len(list(tree.leaves())) > int(max_visible_tips)
    )
    if collapse_required and (
        trait_path not in (None, '') or tip_image_manifest not in (None, '')
    ):
        raise ValueError(
            "'--max-visible-tips' cannot currently be combined with "
            "'--trait' or '--tip-image-manifest'."
        )
    tree, collapsed_clades = collapse_tree_for_drawing(
        tree,
        max_visible_tips=max_visible_tips,
        label_template=getattr(args, 'collapse_label', None),
        property_aggregation=getattr(
            args,
            'collapse_property_aggregation',
            'none',
        ),
    )
    group_by = getattr(args, 'group_by', None)
    trait_palette = getattr(args, 'trait_palette', 'tab10')
    support_labels = getattr(args, 'support_labels', True)
    support_min = getattr(args, 'support_min', None)
    property_colors = _parse_property_colors(getattr(args, 'property_color', None))
    node_pie_properties = _parse_property_names(getattr(args, 'node_pie_properties', None))
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
    _validate_draw_paths(
        outfile=args.outfile,
        layout_report=layout_report,
        input_paths={
            '--infile': args.infile,
            '--trait': trait_path,
            '--tip-image-manifest': tip_image_manifest,
            **{
                "--tip-image '{}'".format(leaf_name): path
                for leaf_name, path in tip_image_path_by_leaf.items()
            },
        },
    )
    tip_image_by_leaf = _load_tip_images(
        tip_image_path_by_leaf,
        image_size_pt=getattr(args, 'tip_image_size', 18.0),
    )
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
        tip_labels=getattr(args, 'tip_labels', True),
        tip_label_position=getattr(args, 'tip_label_position', 'aligned'),
        tip_label_wrap=getattr(args, 'tip_label_wrap', 'none'),
        tip_spacing=getattr(args, 'tip_spacing', 'uniform'),
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
        layout=getattr(args, 'layout', 'rectangular'),
        spiral_turns=getattr(args, 'spiral_turns', None),
        fan_span=getattr(args, 'fan_span', 180.0),
        unrooted_method=getattr(args, 'unrooted_method', 'equal-angle'),
        daylight_iterations=getattr(args, 'daylight_iterations', 5),
        scale_bar=getattr(args, 'scale_bar', 'none'),
        depth_guide=getattr(args, 'depth_guide', 'none'),
        branch_length_unit=getattr(args, 'branch_length_unit', ''),
        branch_width=getattr(args, 'branch_width', 0.8),
        branch_color_property=getattr(args, 'branch_color_property', None),
        branch_width_property=getattr(args, 'branch_width_property', None),
        branch_width_range=getattr(args, 'branch_width_range', '0.4,2.5'),
        tip_label_font_style=getattr(args, 'tip_label_font_style', 'plain'),
        tip_track_properties=getattr(args, 'tip_track', None),
        tip_track_type=getattr(args, 'tip_track_type', 'auto'),
        tip_track_size=getattr(args, 'tip_track_size', 5.0),
        tip_track_palette=getattr(args, 'tip_track_palette', 'viridis'),
        legend_columns=getattr(args, 'legend_columns', 'auto'),
        legend_position=getattr(args, 'legend_position', 'auto'),
        collision_policy=getattr(args, 'collision_policy', 'resolve'),
        layout_report=layout_report,
        collapsed_clades=collapsed_clades,
    )
    sys.stderr.write('Wrote tree image: {}\n'.format(os.path.realpath(args.outfile)))
