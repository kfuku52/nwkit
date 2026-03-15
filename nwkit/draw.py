import os
import re
import sys

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from nwkit.util import compile_species_regex, extract_species_label, is_rooted, read_tree


TREE_LINE_CAPSTYLE = 'round'
FONT_FAMILY = 'Helvetica'
FONT_SIZE_PT = 8.0
TIP_LABEL_GAP_PT = 0.5
TIP_ROW_PITCH_PT = FONT_SIZE_PT + TIP_LABEL_GAP_PT
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


def _get_species_by_leaf(tree, species_pattern):
    species_by_leaf = dict()
    is_all_parsed = True
    num_leaf = 0
    for leaf in tree.leaves():
        num_leaf += 1
        species_label = extract_species_label(leaf.name or '', species_regex=species_pattern.pattern)
        if species_label is None:
            is_all_parsed = False
            continue
        species_by_leaf[leaf] = species_label
    return species_by_leaf, (is_all_parsed and (num_leaf > 0))


def _get_species_overlap_node_types(tree, species_regex, require_all_tip_labels=False):
    species_pattern = compile_species_regex(species_regex=species_regex)
    if species_pattern is None:
        return dict(), False
    species_by_leaf, all_tip_labels_parsed = _get_species_by_leaf(tree=tree, species_pattern=species_pattern)
    if bool(require_all_tip_labels) and (not all_tip_labels_parsed):
        return dict(), all_tip_labels_parsed
    node_type_by_node = dict()
    for node in tree.traverse():
        if node.is_leaf:
            continue
        children = node.get_children()
        if len(children) < 2:
            continue
        child_species_sets = list()
        missing_species = False
        for child in children:
            species_set = set()
            for leaf in child.leaves():
                species_label = species_by_leaf.get(leaf, None)
                if species_label is None:
                    missing_species = True
                    break
                species_set.add(species_label)
            if missing_species:
                break
            child_species_sets.append(species_set)
        if missing_species:
            continue
        is_duplication = False
        for i in range(len(child_species_sets)):
            for j in range(i + 1, len(child_species_sets)):
                if len(child_species_sets[i].intersection(child_species_sets[j])) > 0:
                    is_duplication = True
                    break
            if is_duplication:
                break
        node_type_by_node[node] = 'duplication' if is_duplication else 'speciation'
    return node_type_by_node, all_tip_labels_parsed


def _apply_font_style():
    matplotlib.rcParams['font.size'] = FONT_SIZE_PT
    matplotlib.rc('xtick', labelsize=FONT_SIZE_PT)
    matplotlib.rc('ytick', labelsize=FONT_SIZE_PT)
    matplotlib.rc('font', size=FONT_SIZE_PT, family=FONT_FAMILY)
    matplotlib.rc('axes', titlesize=FONT_SIZE_PT, labelsize=FONT_SIZE_PT)
    matplotlib.rc('legend', fontsize=FONT_SIZE_PT)
    matplotlib.rc('figure', titlesize=FONT_SIZE_PT)


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


def _draw_tree(tree, outfile, image_format, node_type_by_node):
    _apply_font_style()
    use_topology_depth = (not _has_positive_branch_length(tree))
    if use_topology_depth:
        sys.stderr.write('Tree has no positive branch lengths; drawing x positions by topology depth.\n')
    xcoord, ycoord, leaf_order = _get_tree_plot_coordinates(tree=tree, use_topology_depth=use_topology_depth)
    x_values = list(xcoord.values())
    x_min = min(x_values) if len(x_values) > 0 else 0.0
    x_max = max(x_values) if len(x_values) > 0 else 1.0
    x_span = max(x_max - x_min, 1.0)
    max_tip_label_chars = max([len(leaf.name or '') for leaf in leaf_order], default=1)

    fig_width = FIGURE_WIDTH_IN
    label_panel_width_in = min(max(0.8, max_tip_label_chars * 0.07), fig_width * 0.58)
    tree_panel_width_in = max(fig_width - label_panel_width_in, 0.8)
    row_pitch_in = TIP_ROW_PITCH_PT / 72.0
    tree_panel_height_in = max(len(leaf_order), 1) * row_pitch_in
    top_margin_in = (28.0 / 72.0) if (len(node_type_by_node) > 0) else (4.0 / 72.0)
    bottom_margin_in = 3.0 / 72.0
    left_margin_in = 2.0 / 72.0
    right_margin_in = 2.0 / 72.0
    fig_height = tree_panel_height_in + top_margin_in + bottom_margin_in

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
                color=BRANCH_COLOR,
                linewidth=0.8,
                zorder=1,
                solid_capstyle=TREE_LINE_CAPSTYLE,
            )

    for node in tree.traverse():
        if node.is_root:
            continue
        ax.plot(
            [xcoord[node.up], xcoord[node]],
            [ycoord[node], ycoord[node]],
            color=BRANCH_COLOR,
            linewidth=0.8,
            zorder=2,
            solid_capstyle=TREE_LINE_CAPSTYLE,
        )

    root_stub = max(x_span * 0.03, 0.03)
    ax.plot(
        [xcoord[tree] - root_stub, xcoord[tree]],
        [ycoord[tree], ycoord[tree]],
        color=BRANCH_COLOR,
        linewidth=0.8,
        zorder=2,
        solid_capstyle=TREE_LINE_CAPSTYLE,
    )

    marker_area = 18.0 * (0.5 ** 2)
    marker_size_pt = 4.8 * 0.5
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
        legend_handles = [
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
        ]
        ax.legend(
            handles=legend_handles,
            loc='lower left',
            bbox_to_anchor=(0.0, 1.005),
            frameon=False,
            borderaxespad=0.1,
            handletextpad=0.3,
            labelspacing=0.2,
            ncol=1,
            prop={'family': FONT_FAMILY, 'size': FONT_SIZE_PT},
        )

    for node in tree.traverse():
        if not _has_meaningful_support(node):
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
            fontsize=FONT_SIZE_PT,
            fontfamily=FONT_FAMILY,
            color=LABEL_COLOR,
        )

    label_offset = max(x_span * 0.02, 0.06)
    for leaf in leaf_order:
        ax.text(
            x_max + label_offset,
            ycoord[leaf],
            leaf.name or '',
            va='center',
            ha='left',
            fontsize=FONT_SIZE_PT,
            fontfamily=FONT_FAMILY,
            color=LABEL_COLOR,
        )

    label_panel_span = x_span * (label_panel_width_in / tree_panel_width_in)
    ax.set_xlim(x_min - root_stub * 1.2, x_max + label_offset + label_panel_span)
    if len(leaf_order) == 0:
        ax.set_ylim(0.5, -0.5)
    else:
        ax.set_ylim(float(len(leaf_order)) - 0.5, -0.5)
    ax.axis('off')
    fig.subplots_adjust(
        left=left_margin_in / fig_width,
        right=1.0 - (right_margin_in / fig_width),
        top=1.0 - (top_margin_in / fig_height),
        bottom=bottom_margin_in / fig_height,
    )
    fig.savefig(outfile, format=image_format, dpi=300)
    plt.close(fig)


def draw_main(args):
    if args.outfile == '-':
        raise ValueError("STDOUT is not supported for 'draw'. Use --outfile PATH.")
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    if bool(args.ladderize):
        tree.ladderize()
    image_format = _resolve_image_format(outfile=args.outfile, image_format=str(args.image_format).lower())
    node_plot_mode = str(args.species_overlap_node_plot).strip().lower()
    if node_plot_mode == 'no':
        node_type_by_node = dict()
    else:
        if not is_rooted(tree):
            if node_plot_mode == 'yes':
                raise ValueError(
                    "Speciation/duplication node plotting requires a rooted tree. "
                    "Use '--species_overlap_node_plot no' for unrooted trees."
                )
            node_type_by_node = dict()
            sys.stderr.write(
                'Skipping speciation/duplication node markers because the input tree is unrooted.\n'
            )
        else:
            node_type_by_node, all_tip_labels_parsed = _get_species_overlap_node_types(
                tree=tree,
                species_regex=args.species_regex,
                require_all_tip_labels=(node_plot_mode == 'auto'),
            )
            if (node_plot_mode == 'auto') and (not all_tip_labels_parsed):
                sys.stderr.write(
                    'Skipping speciation/duplication node markers because some leaf labels did not match --species_regex.\n'
                )
    outdir = os.path.dirname(os.path.realpath(args.outfile))
    os.makedirs(outdir, exist_ok=True)
    _draw_tree(
        tree=tree,
        outfile=args.outfile,
        image_format=image_format,
        node_type_by_node=node_type_by_node,
    )
    sys.stderr.write('Wrote tree image: {}\n'.format(os.path.realpath(args.outfile)))
