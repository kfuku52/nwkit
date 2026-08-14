#!/usr/bin/env python3
"""Generate the reproducible 7.2-inch ``nwkit draw`` gallery."""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
from ete4 import Tree

from nwkit.draw import _measure_texts_in_inches
from nwkit.draw_layouts import make_tree_layout

FONT_FAMILY = "Helvetica"
FONT_SIZE = 8.0
SMALL_NEWICK = (
    "((((A:0.5,J:3.2):0.4,C:1.0):1.0,(D:0.5,E:2.3):1.8):0.6,"
    "((F:3.0,(G:0.4,H:1.7):0.5):0.7,(I:0.6,B:2.5):1.2):1.0);"
)
SMALL_LABELS = {
    "A": "A",
    "B": "B\nexceptionally\nlong taxon\ncollection\nlabel",
    "C": "C",
    "D": "D",
    "E": "E",
    "F": "F",
    "G": "G\nmultiline\nlabel",
    "H": "H",
    "I": "I",
    "J": "J",
}
PANELS = (
    ("rectangular", "rectangular", "standard", 360.0, "orthogonal phylogram"),
    ("rectangular + tidy", "rectangular", "tidy", 360.0, "2D compact packing"),
    ("slanted", "slanted", "standard", 360.0, "straight rooted phylogram"),
    ("cladogram", "cladogram", "standard", 360.0, "aligned topology"),
    ("circular · 360°", "circular", "standard", 360.0, "complete circle"),
    ("circular · 180°", "circular", "standard", 180.0, "upper-half fan"),
    ("radial", "radial", "standard", 360.0, "straight rooted radial"),
    ("unrooted", "unrooted", "standard", 360.0, "equal-angle unrooted"),
    ("spiral", "spiral", "standard", 360.0, "injective spiral band"),
    ("fractal", "fractal", "standard", 360.0, "rectangle-fitted topology"),
)

TIME_PANELS = (
    (
        "calibration constraints",
        "bounded · min · max · fixed",
        "constraints",
    ),
    (
        "dated tree + 95% HPD",
        "posterior ages + uncertainty",
        "dated",
    ),
    (
        "DensiTree samples",
        "160 trees · topology + age",
        "densitree-all",
    ),
    (
        "topology-stratified envelopes",
        "95% paths · 70/20/10%",
        "densitree-ci",
    ),
)


def _hundred_tip_newick():
    ratios = (0.44, 0.58, 0.36, 0.62)

    def crown(start, count, depth, attached=True):
        if count == 1:
            body = "T{:03d}".format(start)
            length = 0.10 + (((start * 37) % 17) / 80.0)
        else:
            left_count = round(count * ratios[depth % len(ratios)])
            left_count = max(1, min(count - 1, left_count))
            body = "({},{})".format(
                crown(start, left_count, depth + 1),
                crown(start + left_count, count - left_count, depth + 1),
            )
            length = 0.055 + (((start + count + depth * 7) % 11) / 100.0)
        return "{}:{:.3f}".format(body, length) if attached else body

    clades = []
    start = 1
    for count, stem in zip((18, 22, 27, 33), (0.45, 2.45, 4.45, 6.45), strict=True):
        clades.append(
            "{}:{:.3f}".format(
                crown(start, count, 0, attached=False),
                stem,
            )
        )
        start += count
    return "(({},{}):0.45,({},{}):0.45);".format(*clades)


def _layout_for(panel, dense):
    _, layout_name, packing, span, _ = panel
    newick = _hundred_tip_newick() if dense else SMALL_NEWICK
    labels = (
        {"T{:03d}".format(index): "T{:03d}".format(index) for index in range(1, 101)}
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
        tip_spacing="uniform" if dense else "label-aware",
        label_size_by_leaf=sizes,
        terminal_extent_by_leaf=extents,
    )
    return tree, drawing, labels


def _spatial_label_alignment(angle):
    degrees = float(angle) % 360.0
    if 90.0 < degrees < 270.0:
        return degrees + 180.0, "right"
    return degrees, "left"


def _draw_panel(axes, panel, dense, show_title):
    title, layout_name, packing, span, descriptor = panel
    tree, drawing, labels = _layout_for(panel, dense)
    terminal_width = 0.35 if dense else 0.68
    internal_width = 0.27 if dense else 0.56
    for node, path in drawing.edge_paths.items():
        x_values, y_values = zip(*path, strict=True)
        axes.plot(
            x_values,
            y_values,
            color="#315A72" if node.is_leaf else "#222222",
            linewidth=terminal_width if node.is_leaf else internal_width,
            solid_capstyle="round",
            solid_joinstyle="round",
            zorder=2,
        )
    if drawing.root_path and layout_name != "unrooted":
        x_values, y_values = zip(*drawing.root_path, strict=True)
        axes.plot(
            x_values,
            y_values,
            color="#222222",
            linewidth=internal_width,
            solid_capstyle="round",
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
                    textcoords="offset points",
                    ha=alignment,
                    va="center",
                    rotation=rotation,
                    rotation_mode="anchor",
                    fontsize=FONT_SIZE,
                    fontfamily=FONT_FAMILY,
                    linespacing=1.12,
                    color="#202020",
                    clip_on=False,
                )
            else:
                axes.annotate(
                    labels[leaf.name],
                    (drawing.xcoord[leaf], drawing.ycoord[leaf]),
                    xytext=(3.0, 0.0),
                    textcoords="offset points",
                    ha="left",
                    va="center",
                    fontsize=FONT_SIZE,
                    fontfamily=FONT_FAMILY,
                    linespacing=1.12,
                    color="#202020",
                    clip_on=False,
                )
    x_min, x_max, y_min, y_max = drawing.bounds
    width = max(x_max - x_min, 1e-6)
    height = max(y_max - y_min, 1e-6)
    if drawing.equal_aspect:
        axes.set_aspect("equal", adjustable="datalim")
    x_padding = width * (0.10 if dense else (0.46 if drawing.spatial else 0.16))
    y_padding = height * (0.10 if dense else (0.64 if drawing.spatial else 0.16))
    axes.set_xlim(x_min - x_padding, x_max + x_padding)
    axes.set_ylim(y_min - y_padding, y_max + y_padding)
    axes.axis("off")
    if show_title:
        axes.text(
            0.0,
            1.08,
            title,
            transform=axes.transAxes,
            ha="left",
            va="bottom",
            fontsize=FONT_SIZE,
            fontfamily=FONT_FAMILY,
            fontweight="bold",
            color="#111111",
        )
        axes.text(
            1.0,
            1.08,
            descriptor,
            transform=axes.transAxes,
            ha="right",
            va="bottom",
            fontsize=FONT_SIZE,
            fontfamily=FONT_FAMILY,
            color="#315A72",
        )
    if dense and layout_name == "rectangular":
        axes.text(
            0.0,
            0.0,
            "standard rows" if packing == "standard" else "tidy contour packing",
            transform=axes.transAxes,
            ha="left",
            va="bottom",
            fontsize=FONT_SIZE,
            fontfamily=FONT_FAMILY,
            color="#666666",
        )


def _time_example_files(directory):
    constraint_tree = directory / "constraints.nwk"
    constraint_tree.write_text(
        "(((A:18,B:25)'B(18,28,0.025,0.025)':40,"
        "(C:27,D:34)'>22':31)'<80':38,"
        "(E:20,F:31)'@30':70);\n"
    )
    topology = directory / "topology.nwk"
    topology.write_text("(((A:22,B:22):40,(C:31,D:31):31):38,(E:37,F:37):63);\n")
    posterior = directory / "mcmc.txt"
    rows = ["Gen t_n7 t_n8 t_n9 t_n10 t_n11"]
    for index in range(1, 161):
        phase = 2.0 * math.pi * index / 41.0
        root = 100.0 + 5.0 * math.sin(phase)
        left = 62.0 + 4.0 * math.sin(phase * 1.37 + 0.2)
        ab = 22.0 + 3.0 * math.sin(phase * 2.11 + 0.7)
        cd = 31.0 + 4.0 * math.sin(phase * 1.73 + 1.3)
        ef = 37.0 + 3.5 * math.sin(phase * 1.91 + 2.8)
        rows.append(
            "{} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f}".format(
                index,
                root,
                left,
                ab,
                cd,
                ef,
            )
        )
    posterior.write_text("\n".join(rows) + "\n")
    tree_collection = directory / "posterior-trees.nwk"
    sample_newicks = []
    for index in range(1, 161):
        phase = 2.0 * math.pi * index / 41.0
        root = 100.0 + 5.0 * math.sin(phase)
        left = 62.0 + 4.0 * math.sin(phase * 1.37 + 0.2)
        first_left = 22.0 + 3.0 * math.sin(phase * 2.11 + 0.7)
        second_left = 31.0 + 4.0 * math.sin(phase * 1.73 + 1.3)
        ef = 37.0 + 3.5 * math.sin(phase * 1.91 + 2.8)
        mixture_position = (index - 1) % 20

        def pair_newick(pair, age, parent_age):
            return "({}:{:.6f},{}:{:.6f}):{:.6f}".format(
                pair[0],
                age,
                pair[1],
                age,
                parent_age - age,
            )

        if mixture_position < 14:
            left_pairs = (("A", "B"), ("C", "D"))
            left_newick = "({},{}):{:.6f}".format(
                pair_newick(left_pairs[0], first_left, left),
                pair_newick(left_pairs[1], second_left, left),
                root - left,
            )
            right_newick = pair_newick(("E", "F"), ef, root)
        elif mixture_position < 18:
            left_pairs = (("A", "C"), ("B", "D"))
            left_newick = "({},{}):{:.6f}".format(
                pair_newick(left_pairs[0], first_left, left),
                pair_newick(left_pairs[1], second_left, left),
                root - left,
            )
            right_newick = pair_newick(("E", "F"), ef, root)
        else:
            mixed_age = 70.0 + 3.0 * math.sin(phase * 1.19 + 1.8)
            left_newick = "({},{}):{:.6f}".format(
                pair_newick(("A", "B"), first_left, mixed_age),
                pair_newick(("E", "F"), ef, mixed_age),
                root - mixed_age,
            )
            right_newick = pair_newick(("C", "D"), second_left, root)
        sample_newicks.append("({},{});".format(left_newick, right_newick))
    tree_collection.write_text("\n".join(sample_newicks) + "\n")
    return constraint_tree, topology, posterior, tree_collection


def _render_time_examples(directory, figure_width, figure_height):
    constraint_tree, topology, posterior, tree_collection = _time_example_files(
        directory
    )
    common = [
        "--image-format",
        "png",
        "--figure-width",
        "{:.6f}".format(figure_width),
        "--figure-height",
        "{:.6f}".format(figure_height),
        "--font-size",
        str(FONT_SIZE),
        "--font-family",
        FONT_FAMILY,
        "--support-labels",
        "no",
        "--species-overlap-node-plot",
        "no",
        "--collision-policy",
        "resolve",
        "--branch-width",
        "0.7",
    ]
    specifications = {
        "constraints": [
            "-i",
            str(constraint_tree),
            "--layout",
            "rectangular",
            "--tip-label-position",
            "aligned",
            "--time-constraints",
            "yes",
            "--scale-bar",
            "auto",
            "--branch-length-unit",
            "Ma",
        ],
        "dated": [
            "-i",
            str(topology),
            "--layout",
            "rectangular",
            "--tip-label-position",
            "aligned",
            "--mcmctree-posterior",
            str(posterior),
            "--time-credible-intervals",
            "yes",
            "--scale-bar",
            "auto",
            "--branch-length-unit",
            "Ma",
        ],
        "densitree-all": [
            "-i",
            str(topology),
            "--layout",
            "slanted",
            "--tip-label-position",
            "aligned",
            "--densitree-trees",
            str(tree_collection),
            "--time-credible-intervals",
            "no",
            "--densitree",
            "all",
            "--densitree-alpha",
            "0.065",
        ],
        "densitree-ci": [
            "-i",
            str(topology),
            "--layout",
            "slanted",
            "--tip-label-position",
            "aligned",
            "--densitree-trees",
            str(tree_collection),
            "--time-credible-intervals",
            "no",
            "--densitree",
            "ci",
            "--densitree-ci-alpha",
            "0.34",
        ],
    }
    images = {}
    for key, arguments in specifications.items():
        output = directory / "{}.png".format(key)
        report = directory / "{}.json".format(key)
        command = [
            sys.executable,
            "-m",
            "nwkit",
            "draw",
            *arguments,
            *common,
            "--layout-report",
            str(report),
            "-o",
            str(output),
        ]
        subprocess.run(
            command,
            cwd=Path(__file__).resolve().parents[1],
            check=True,
            capture_output=True,
            text=True,
        )
        quality = json.loads(report.read_text())
        if quality["final_collision_count"] != 0:
            raise RuntimeError(
                "{} panel retained {} annotation collision(s).".format(
                    key,
                    quality["final_collision_count"],
                )
            )
        if quality["maximum_overflow_points"] > 0.25:
            raise RuntimeError(
                "{} panel exceeds its figure bounds by {:.3f} pt.".format(
                    key,
                    quality["maximum_overflow_points"],
                )
            )
        images[key] = plt.imread(output)
    return images


def _draw_time_panel(axes, image):
    axes.imshow(image)
    axes.axis("off")


def _draw_card_header(axes, title, descriptor):
    axes.axis("off")
    axes.text(
        0.0,
        0.98,
        title,
        transform=axes.transAxes,
        ha="left",
        va="top",
        fontsize=FONT_SIZE,
        fontfamily=FONT_FAMILY,
        fontweight="bold",
        color="#111111",
    )
    axes.text(
        0.0,
        0.08,
        descriptor,
        transform=axes.transAxes,
        ha="left",
        va="bottom",
        fontsize=FONT_SIZE,
        fontfamily=FONT_FAMILY,
        color="#315A72",
    )


def make_gallery(output_prefix):
    plt.rcParams.update(
        {
            "font.family": FONT_FAMILY,
            "font.size": FONT_SIZE,
            "svg.fonttype": "none",
        }
    )
    figure = plt.figure(figsize=(7.2, 14.7), dpi=220)
    grid = figure.add_gridspec(
        5,
        2,
        left=0.055,
        right=0.955,
        top=0.94,
        bottom=0.225,
        hspace=0.12,
        wspace=0.14,
    )
    layout_axes = []
    for index, panel in enumerate(PANELS):
        row = index // 2
        column = index % 2
        card = grid[row, column].subgridspec(
            2,
            2,
            height_ratios=(0.12, 0.88),
            width_ratios=(1.95, 0.85),
            hspace=0.0,
            wspace=0.02,
        )
        header = figure.add_subplot(card[0, :])
        small_axis = figure.add_subplot(card[1, 0])
        dense_axis = figure.add_subplot(card[1, 1])
        _draw_card_header(header, panel[0], panel[4])
        _draw_panel(small_axis, panel, dense=False, show_title=False)
        _draw_panel(dense_axis, panel, dense=True, show_title=False)
        layout_axes.append((small_axis, dense_axis))
    time_grid = figure.add_gridspec(
        1,
        4,
        left=0.055,
        right=0.955,
        top=0.172,
        bottom=0.012,
        wspace=0.10,
    )
    time_image_axes = []
    for index, (title, descriptor, key) in enumerate(TIME_PANELS):
        row = 0
        column = index
        card = time_grid[row, column].subgridspec(
            2,
            1,
            height_ratios=(0.15, 0.85),
            hspace=0.0,
        )
        header = figure.add_subplot(card[0, 0])
        image_axis = figure.add_subplot(card[1, 0])
        _draw_card_header(header, title, descriptor)
        time_image_axes.append((image_axis, key))
    figure.canvas.draw()
    image_bounds = time_image_axes[0][0].get_position()
    time_figure_width = image_bounds.width * figure.get_figwidth()
    time_figure_height = image_bounds.height * figure.get_figheight()
    with tempfile.TemporaryDirectory(prefix="nwkit-draw-gallery-") as temporary:
        time_images = _render_time_examples(
            Path(temporary),
            figure_width=time_figure_width,
            figure_height=time_figure_height,
        )
    for image_axis, key in time_image_axes:
        _draw_time_panel(image_axis, time_images[key])
    figure.suptitle(
        "NWKIT tree drawing — layouts and time-aware overlays",
        x=0.055,
        y=0.993,
        ha="left",
        va="top",
        fontsize=FONT_SIZE,
        fontfamily=FONT_FAMILY,
        fontweight="bold",
    )
    figure.text(
        0.055,
        0.976,
        "Eight geometries plus tidy-packed and 180° circular configurations.",
        ha="left",
        va="top",
        fontsize=FONT_SIZE,
        fontfamily=FONT_FAMILY,
        color="#444444",
    )
    figure.text(
        0.055,
        0.959,
        "Each card: 10 tips with labels (left) · 100 tips without labels (right) · all text 8 pt Helvetica",
        ha="left",
        va="top",
        fontsize=FONT_SIZE,
        fontfamily=FONT_FAMILY,
        color="#444444",
    )
    figure.text(
        0.055,
        0.212,
        "Time-aware drawing modes",
        ha="left",
        va="top",
        fontsize=FONT_SIZE,
        fontfamily=FONT_FAMILY,
        fontweight="bold",
        color="#111111",
    )
    figure.text(
        0.055,
        0.195,
        "One row: calibration forms · dated-node HPD · topology/age samples · frequency-scaled within-topology envelopes.",
        ha="left",
        va="top",
        fontsize=FONT_SIZE,
        fontfamily=FONT_FAMILY,
        color="#444444",
    )
    output = Path(output_prefix).expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output.with_suffix(".png"), dpi=300, facecolor="white")
    figure.savefig(output.with_suffix(".svg"), facecolor="white")
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-prefix",
        default="img/nwkit_draw_layout_gallery",
        help="Output path without extension; both PNG and SVG are written.",
    )
    args = parser.parse_args()
    make_gallery(args.output_prefix)


if __name__ == "__main__":
    main()
