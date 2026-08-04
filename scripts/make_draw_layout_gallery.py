#!/usr/bin/env python3
"""Generate the reproducible 10-versus-100-tip ``nwkit draw`` gallery."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib.pyplot as plt
from ete4 import Tree

from nwkit.draw import _measure_texts_in_inches
from nwkit.draw_layouts import make_tree_layout


FONT_FAMILY = 'Helvetica'
FONT_SIZE = 8.0
SMALL_NEWICK = (
    '((((A:0.5,J:3.2):0.4,C:1.0):1.0,(D:0.5,E:2.3):1.8):0.6,'
    '((F:3.0,(G:0.4,H:1.7):0.5):0.7,(I:0.6,B:2.5):1.2):1.0);'
)
SMALL_LABELS = {
    'A': 'A',
    'B': 'B\nexceptionally\nlong taxon\ncollection\nlabel',
    'C': 'C',
    'D': 'D',
    'E': 'E',
    'F': 'F',
    'G': 'G\nmultiline\nlabel',
    'H': 'H',
    'I': 'I',
    'J': 'J',
}
PANELS = (
    ('rectangular', 'rectangular', 'standard', 360.0, 'orthogonal phylogram'),
    ('rectangular + tidy', 'rectangular', 'tidy', 360.0, '2D compact packing'),
    ('slanted', 'slanted', 'standard', 360.0, 'straight rooted phylogram'),
    ('cladogram', 'cladogram', 'standard', 360.0, 'aligned topology'),
    ('circular · 360°', 'circular', 'standard', 360.0, 'complete circle'),
    ('circular · 180°', 'circular', 'standard', 180.0, 'upper-half fan'),
    ('radial', 'radial', 'standard', 360.0, 'straight rooted radial'),
    ('unrooted', 'unrooted', 'standard', 360.0, 'equal-angle unrooted'),
    ('spiral', 'spiral', 'standard', 360.0, 'injective spiral band'),
    ('fractal', 'fractal', 'standard', 360.0, 'rectangle-fitted topology'),
)


def _hundred_tip_newick():
    ratios = (0.44, 0.58, 0.36, 0.62)

    def crown(start, count, depth, attached=True):
        if count == 1:
            body = 'T{:03d}'.format(start)
            length = 0.10 + (((start * 37) % 17) / 80.0)
        else:
            left_count = round(count * ratios[depth % len(ratios)])
            left_count = max(1, min(count - 1, left_count))
            body = '({},{})'.format(
                crown(start, left_count, depth + 1),
                crown(start + left_count, count - left_count, depth + 1),
            )
            length = 0.055 + (((start + count + depth * 7) % 11) / 100.0)
        return '{}:{:.3f}'.format(body, length) if attached else body

    clades = []
    start = 1
    for count, stem in zip((18, 22, 27, 33), (0.45, 2.45, 4.45, 6.45)):
        clades.append('{}:{:.3f}'.format(
            crown(start, count, 0, attached=False),
            stem,
        ))
        start += count
    return '(({},{}):0.45,({},{}):0.45);'.format(*clades)


def _layout_for(panel, dense):
    _, layout_name, packing, span, _ = panel
    newick = _hundred_tip_newick() if dense else SMALL_NEWICK
    labels = (
        {'T{:03d}'.format(index): 'T{:03d}'.format(index) for index in range(1, 101)}
        if dense
        else SMALL_LABELS
    )
    tree = Tree(newick, parser=1)
    if dense:
        sizes = {leaf: (0.0, 0.0) for leaf in tree.leaves()}
        extents = {leaf: 0.0 for leaf in tree.leaves()}
    else:
        measured = _measure_texts_in_inches(
            labels.values(),
            font_size=FONT_SIZE,
            font_family=FONT_FAMILY,
        )
        sizes = {leaf: measured[labels[leaf.name]] for leaf in tree.leaves()}
        extents = {
            leaf: 1.0 + (0.3 * len(labels[leaf.name].splitlines()))
            for leaf in tree.leaves()
        }
    drawing = make_tree_layout(
        tree,
        layout=layout_name,
        subtree_packing=packing,
        angular_span=span,
        angular_center=90.0,
        spiral_turns=1.65,
        aspect_ratio=1.6,
        tip_spacing='uniform' if dense else 'label-aware',
        label_size_by_leaf=sizes,
        terminal_extent_by_leaf=extents,
    )
    return tree, drawing, labels


def _spatial_label_alignment(angle):
    degrees = float(angle) % 360.0
    if 90.0 < degrees < 270.0:
        return degrees + 180.0, 'right'
    return degrees, 'left'


def _draw_panel(axes, panel, dense, show_title):
    title, layout_name, packing, span, descriptor = panel
    tree, drawing, labels = _layout_for(panel, dense)
    terminal_width = 0.35 if dense else 0.68
    internal_width = 0.27 if dense else 0.56
    for node, path in drawing.edge_paths.items():
        x_values, y_values = zip(*path)
        axes.plot(
            x_values,
            y_values,
            color='#315A72' if node.is_leaf else '#222222',
            linewidth=terminal_width if node.is_leaf else internal_width,
            solid_capstyle='round',
            solid_joinstyle='round',
            zorder=2,
        )
    if drawing.root_path and layout_name != 'unrooted':
        x_values, y_values = zip(*drawing.root_path)
        axes.plot(
            x_values,
            y_values,
            color='#222222',
            linewidth=internal_width,
            solid_capstyle='round',
        )
    if not dense:
        for leaf in drawing.leaf_order:
            if drawing.spatial:
                angle = drawing.label_angles.get(leaf, 0.0)
                radians = math.radians(angle)
                rotation, alignment = _spatial_label_alignment(angle)
                offset = (4.0 * math.cos(radians), 4.0 * math.sin(radians))
                axes.annotate(
                    labels[leaf.name],
                    (drawing.xcoord[leaf], drawing.ycoord[leaf]),
                    xytext=offset,
                    textcoords='offset points',
                    ha=alignment,
                    va='center',
                    rotation=rotation,
                    rotation_mode='anchor',
                    fontsize=FONT_SIZE,
                    fontfamily=FONT_FAMILY,
                    linespacing=1.12,
                    color='#202020',
                    clip_on=False,
                )
            else:
                axes.annotate(
                    labels[leaf.name],
                    (drawing.xcoord[leaf], drawing.ycoord[leaf]),
                    xytext=(3.0, 0.0),
                    textcoords='offset points',
                    ha='left',
                    va='center',
                    fontsize=FONT_SIZE,
                    fontfamily=FONT_FAMILY,
                    linespacing=1.12,
                    color='#202020',
                    clip_on=False,
                )
    x_min, x_max, y_min, y_max = drawing.bounds
    width = max(x_max - x_min, 1e-6)
    height = max(y_max - y_min, 1e-6)
    if drawing.equal_aspect:
        axes.set_aspect('equal', adjustable='datalim')
    x_padding = width * (0.10 if dense else (0.46 if drawing.spatial else 0.16))
    y_padding = height * (0.10 if dense else (0.48 if drawing.spatial else 0.16))
    axes.set_xlim(x_min - x_padding, x_max + x_padding)
    axes.set_ylim(y_min - y_padding, y_max + y_padding)
    axes.axis('off')
    if show_title:
        axes.text(
            0.0,
            1.08,
            title,
            transform=axes.transAxes,
            ha='left',
            va='bottom',
            fontsize=FONT_SIZE,
            fontfamily=FONT_FAMILY,
            fontweight='bold',
            color='#111111',
        )
        axes.text(
            1.0,
            1.08,
            descriptor,
            transform=axes.transAxes,
            ha='right',
            va='bottom',
            fontsize=FONT_SIZE,
            fontfamily=FONT_FAMILY,
            color='#315A72',
        )
    if dense and layout_name == 'rectangular':
        axes.text(
            0.0,
            0.0,
            'standard rows' if packing == 'standard' else 'tidy contour packing',
            transform=axes.transAxes,
            ha='left',
            va='bottom',
            fontsize=FONT_SIZE,
            fontfamily=FONT_FAMILY,
            color='#666666',
        )


def make_gallery(output_prefix):
    plt.rcParams.update({
        'font.family': FONT_FAMILY,
        'font.size': FONT_SIZE,
        'svg.fonttype': 'none',
    })
    figure = plt.figure(figsize=(14.4, 15.0), dpi=220)
    grid = figure.add_gridspec(
        5,
        5,
        left=0.035,
        right=0.965,
        top=0.90,
        bottom=0.045,
        width_ratios=(1.0, 1.0, 0.18, 1.0, 1.0),
        hspace=0.54,
        wspace=0.12,
    )
    axes = [
        [figure.add_subplot(grid[row, column]) for column in (0, 1, 3, 4)]
        for row in range(5)
    ]
    for index, panel in enumerate(PANELS):
        row = index // 2
        column = (index % 2) * 2
        _draw_panel(axes[row][column], panel, dense=False, show_title=True)
        _draw_panel(axes[row][column + 1], panel, dense=True, show_title=False)
    figure.suptitle(
        'NWKIT tree drawing configurations — 10 versus 100 tips',
        x=0.035,
        y=0.992,
        ha='left',
        va='top',
        fontsize=FONT_SIZE,
        fontfamily=FONT_FAMILY,
        fontweight='bold',
    )
    figure.text(
        0.035,
        0.970,
        'Eight geometries plus tidy-packed and 180° circular configurations; every label and caption is 8 pt Helvetica.',
        ha='left',
        va='top',
        fontsize=FONT_SIZE,
        fontfamily=FONT_FAMILY,
        color='#444444',
    )
    for axis, label in zip(
        axes[0],
        ('10 tips · labels shown', '100 tips · labels suppressed') * 2,
    ):
        bounds = axis.get_position()
        figure.text(
            (bounds.x0 + bounds.x1) / 2.0,
            0.922,
            label,
            ha='center',
            va='bottom',
            fontsize=FONT_SIZE,
            fontfamily=FONT_FAMILY,
            fontweight='bold',
        )
    output = Path(output_prefix).expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output.with_suffix('.png'), dpi=300, facecolor='white')
    figure.savefig(output.with_suffix('.svg'), facecolor='white')
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--output-prefix',
        default='img/nwkit_draw_layout_gallery',
        help='Output path without extension; both PNG and SVG are written.',
    )
    args = parser.parse_args()
    make_gallery(args.output_prefix)


if __name__ == '__main__':
    main()
