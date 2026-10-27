"""Vector SVG/PDF output of a neutral scene, with explicit PDF pagination."""

from __future__ import annotations

import logging
import math
import textwrap
from contextlib import contextmanager

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch

from nwkit import __version__
from nwkit.reconciliation_bounds import fitted_scene
from nwkit.reconciliation_scene import Primitive

PAGE_SIZES = {"A4": (841.89, 595.28), "A3": (1190.55, 841.89)}


@contextmanager
def _quiet_font_subsetting():
    logger = logging.getLogger("fontTools.subset")
    level = logger.level
    logger.setLevel(logging.WARNING)
    try:
        yield
    finally:
        logger.setLevel(level)


def _append_notes(scene, font_size):
    primitives = list(scene.primitives)
    y = scene.height + 20
    line_width = max(20, int((scene.width - 60) / (font_size * 0.6)))
    for index, notice in enumerate(scene.notices):
        for line_index, text in enumerate(textwrap.wrap(notice, width=line_width)):
            primitives.append(
                Primitive(
                    "text",
                    f"notice:{index}:{line_index}",
                    [(30, y)],
                    text=text,
                    color="#596b80",
                )
            )
            y += font_size * 1.5
    return primitives, max(scene.height, y + 10)


def _primitive(ax, primitive, font_size, scale=1):
    points = primitive.points
    if primitive.kind in {"line", "upper"}:
        if primitive.kind == "upper" and primitive.outline_color:
            outline = ax.plot(
                [x for x, _ in points],
                [y for _, y in points],
                color=primitive.outline_color,
                linewidth=(primitive.width + 1.8) * scale,
                solid_capstyle="butt",
                solid_joinstyle="miter",
                zorder=0.5,
            )[0]
            outline.set_gid(f"outline:{primitive.id}")
        artist = ax.plot(
            [x for x, _ in points],
            [y for _, y in points],
            color=primitive.color,
            linewidth=primitive.width * scale,
            solid_capstyle="butt" if primitive.kind == "upper" else "round",
            solid_joinstyle="miter" if primitive.kind == "upper" else "round",
            linestyle=(0, (3, 2)) if primitive.dashed else "-",
            zorder=1 if primitive.kind == "upper" else 2 if primitive.dashed else 3,
        )[0]
    elif primitive.kind == "transfer":
        artist = FancyArrowPatch(
            points[0],
            points[-1],
            connectionstyle="arc3,rad=0.25",
            arrowstyle="-|>",
            mutation_scale=8 * scale,
            linewidth=primitive.width * scale,
            color=primitive.color,
            zorder=4,
        )
        ax.add_patch(artist)
    elif primitive.kind == "event":
        artist = ax.plot(
            *points[0],
            marker="o" if primitive.marker == "O" else primitive.marker,
            markersize=(5.5 if primitive.marker in {"x", "O"} else 4) * scale,
            color=primitive.color,
            markerfacecolor="white" if primitive.marker == "O" else primitive.color,
            linestyle="none",
            zorder=5,
        )[0]
    elif primitive.kind == "text":
        artist = ax.text(
            *points[0],
            textwrap.fill(primitive.text, 24)
            if not primitive.id.startswith("notice:")
            else primitive.text,
            fontsize=font_size,
            color=primitive.color,
            va="center",
            zorder=6,
            clip_on=True,
        )
    else:
        raise ValueError(f"Unknown scene primitive: {primitive.kind}.")
    artist.set_gid(primitive.id)


def _boundary_crossings(primitive, boundary):
    """Intersect the drawn polyline/Arc3 curve with a scene-space page edge."""
    if primitive.kind == "transfer":
        (x0, y0), (x1, y1) = primitive.points
        # Arc3 uses display coordinates. Our axes invert y and scale x/y equally.
        cx, cy = (x0 + x1) / 2 - (y1 - y0) * 0.25, (y0 + y1) / 2 + (x1 - x0) * 0.25
        a, b, c = y0 - 2 * cy + y1, 2 * (cy - y0), y0 - boundary
        if abs(a) < 1e-12:
            roots = [-c / b] if abs(b) >= 1e-12 else []
        else:
            discriminant = b * b - 4 * a * c
            roots = (
                [(-b + sign * math.sqrt(discriminant)) / (2 * a) for sign in (-1, 1)]
                if discriminant > 0
                else []
            )
        return [
            (1 - t) ** 2 * x0 + 2 * (1 - t) * t * cx + t**2 * x1
            for t in sorted(roots)
            if 0 < t < 1
        ]
    ys = [y for _, y in primitive.points]
    if not min(ys) < boundary < max(ys):
        return []
    crossings = []
    for (x0, y0), (x1, y1) in zip(
        primitive.points[:-1], primitive.points[1:], strict=True
    ):
        if y0 != y1 and min(y0, y1) <= boundary <= max(y0, y1):
            x = x0 + (x1 - x0) * (boundary - y0) / (y1 - y0)
            if x not in crossings:
                crossings.append(x)
    return crossings


def _continuations(ax, primitives, bottom, top, page_index, page_count, font_size):
    for primitive in primitives:
        if primitive.kind not in {"line", "transfer"}:
            continue
        for boundary in (bottom, top):
            destination = page_index + 2 if boundary == top else page_index
            if not 1 <= destination <= page_count:
                continue
            for crossing, x in enumerate(_boundary_crossings(primitive, boundary)):
                label = (
                    "Transfer" if primitive.kind == "transfer" else "Branch"
                ) + f" continues: p{destination}"
                artist = ax.text(
                    x,
                    boundary + (3 if boundary == bottom else -3),
                    label,
                    fontsize=font_size,
                    color=primitive.color,
                    va="top" if boundary == bottom else "bottom",
                    ha="right" if x > ax.get_xlim()[1] / 2 else "left",
                    clip_on=True,
                )
                artist.set_gid(
                    f"continuation:{primitive.id}:p{page_index + 1}:{boundary}:{crossing}"
                )


def _visible(primitive, bottom, top, font_size):
    ys = [y for _, y in primitive.points]
    if primitive.kind == "transfer":
        (x0, y0), (x1, y1) = primitive.points
        ys.append((y0 + y1) / 2 + (x1 - x0) * 0.25)
        ys.append((y0 + y1) / 2 - (x1 - x0) * 0.25)
    margin = (
        font_size * max(1, len(primitive.text) // 24 + 1)
        if primitive.kind == "text"
        else primitive.width / 2 + 4
    )
    return min(ys) - margin <= top and max(ys) + margin >= bottom


def _figure(
    scene,
    primitives,
    title,
    font_size,
    width,
    height,
    bottom,
    top,
    page_index,
    page_count,
    capacity,
    legend,
):
    fig = plt.figure(figsize=(width / 72, height / 72))
    scale = (height - 140) / (top - bottom)
    ax = fig.add_axes(
        [0, 40 / height, scene.width * scale / width, (height - 140) / height]
    )
    ax.set_xlim(0, scene.width)
    ax.set_ylim(top, bottom)
    ax.axis("off")
    for primitive in primitives:
        if _visible(primitive, bottom, top, font_size):
            _primitive(ax, primitive, font_size * scale, scale)
    if page_count > 1:
        _continuations(ax, primitives, bottom, top, page_index, page_count, font_size)
    fig.text(
        0.04,
        1 - 20 / height,
        title or "Reconciliation",
        fontsize=font_size,
        weight="bold",
        va="top",
    )
    fig.text(
        0.04,
        1 - 41 / height,
        scene.caption,
        fontsize=font_size,
        va="top",
        color="#596b80",
    )
    if legend:
        styles = {p.text: p for p in primitives if p.kind == "event"}
        handles = [
            Line2D(
                [],
                [],
                color=p.color,
                marker="o" if p.marker == "O" else p.marker,
                markerfacecolor="white" if p.marker == "O" else p.color,
                linestyle="none",
                markersize=4,
                label=kind.replace("_", " ").capitalize(),
            )
            for kind, p in sorted(styles.items())
        ]
        if handles:
            fig.legend(
                handles=handles,
                loc="upper left",
                bbox_to_anchor=(0.035, 1 - 52 / height),
                ncol=max(1, min(4, int(width / 150))),
                frameon=False,
                fontsize=font_size,
            )
    fig.text(
        0.04,
        15 / height,
        f"NWKIT {__version__}  |  {page_index + 1}/{page_count}",
        fontsize=font_size,
        color="#596b80",
    )
    return fig


def write_scene(
    scene,
    path,
    *,
    title="",
    image_format="svg",
    font_size=8,
    font_family="DejaVu Sans",
    transparent=False,
    pdf_pages="auto",
    pdf_page_size="A4",
    legend=True,
):
    if not (0 < font_size <= 72):
        raise ValueError("Font size must be positive and at most 72 points.")
    if image_format != "pdf" and (pdf_pages != "auto" or pdf_page_size != "A4"):
        raise ValueError("PDF page options require PDF output.")
    scene = fitted_scene(scene, font_size)
    primitives, natural_height = _append_notes(scene, font_size)
    page_width, page_height = PAGE_SIZES[pdf_page_size]
    scale = min(1.0, page_width / scene.width)
    capacity = (page_height - 140) / scale
    tiled = image_format == "pdf" and (
        pdf_pages == "tiled" or (pdf_pages == "auto" and natural_height > capacity)
    )
    count = math.ceil(natural_height / capacity) if tiled else 1
    width = page_width if tiled else scene.width
    height = page_height if tiled else natural_height + 140
    if max(width, height) > 14400:
        raise ValueError(
            "Single-page output exceeds 200 inches; select tiled PDF output."
        )
    metadata = {"Creator": f"NWKIT {__version__}", "Title": title}
    if image_format == "pdf":
        metadata.update(CreationDate=None, ModDate=None)
    elif image_format == "svg":
        metadata["Date"] = None
    with (
        _quiet_font_subsetting(),
        matplotlib.rc_context(
            {
                "pdf.fonttype": 42,
                "svg.fonttype": "none",
                "svg.hashsalt": "nwkit-reconciliation",
                "font.family": font_family,
            }
        ),
    ):
        if image_format == "pdf":
            with PdfPages(path, metadata=metadata) as pdf:
                for index in range(count):
                    bottom, top = (
                        (index * capacity, (index + 1) * capacity)
                        if tiled
                        else (0, natural_height)
                    )
                    fig = _figure(
                        scene,
                        primitives,
                        title,
                        font_size,
                        width,
                        height,
                        bottom,
                        top,
                        index,
                        count,
                        capacity,
                        legend,
                    )
                    try:
                        pdf.savefig(fig, transparent=transparent)
                    finally:
                        plt.close(fig)
        else:
            fig = _figure(
                scene,
                primitives,
                title,
                font_size,
                width,
                height,
                0,
                natural_height,
                0,
                1,
                capacity,
                legend,
            )
            try:
                fig.savefig(
                    path,
                    format=image_format,
                    metadata=metadata,
                    transparent=transparent,
                )
            finally:
                plt.close(fig)
