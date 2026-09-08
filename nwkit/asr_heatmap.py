"""Observed tip traits aligned to phylogeny columns; independent of fitted means."""

import numpy as np

_CELL_HEIGHT = 0.18
_GAP = 0.10
_KEY_HEIGHT = 0.42
_MISSING = "#D0D0D0"


def heatmap_height(num_traits):
    return _CELL_HEIGHT * num_traits + _GAP


def heatmap_extra_space(num_traits):
    return heatmap_height(num_traits) + _KEY_HEIGHT * num_traits + 0.2


def tip_trait_values(table, tips):
    """Use observed values only, retaining missing entries as NaN."""
    traits = list(dict.fromkeys(table["trait"]))
    values = []
    for trait in traits:
        rows = table[(table.trait == trait) & (table.node_class == "leaf")]
        rows = rows.set_index("name")
        row = []
        for tip in tips:
            value = rows.loc[tip.name, "observed_value"]
            row.append(np.nan if value == "" or value is None else float(value))
        values.append(row)
    return traits, np.asarray(values, dtype=float)


def _below_axis(parent, offset, height, *, width=1.0):
    """Inset with physical-inch spacing that follows later subplot layout changes."""
    from matplotlib.transforms import Bbox

    child = parent.inset_axes([0, 0, 1, 1])

    def locator(ax, renderer):
        bounds = parent.get_position()
        figure_height = parent.figure.get_figheight()
        return Bbox.from_bounds(
            bounds.x0,
            bounds.y0 - (offset + height) / figure_height,
            bounds.width * width,
            height / figure_height,
        )

    child.set_axes_locator(locator)
    return child


def _trait_scale(values):
    from matplotlib.colors import Normalize

    finite = values[np.isfinite(values)]
    if not len(finite):
        return Normalize(0, 1), []
    low, high = float(finite.min()), float(finite.max())
    if low == high:
        margin = max(1.0, abs(low)) * 0.01
        return Normalize(low - margin, high + margin), [low]
    return Normalize(low, high), [low, high]


def _draw_key(parent, offset, trait, norm, ticks, cmap, missing, row_number):
    ax = _below_axis(parent, offset, 0.08, width=0.72)
    ax.set_label(f"tip-heatmap-scale:{trait}")
    if ticks:
        ax.imshow(
            np.linspace(0, 1, 256)[None, :],
            aspect="auto",
            cmap=cmap,
            extent=(norm.vmin, norm.vmax, 0, 1),
            vmin=0,
            vmax=1,
        )
    ax.set_xticks(ticks)
    ax.set_yticks([])
    ax.tick_params(axis="x", labelsize=7, length=2, pad=1)
    prefix = f"{row_number}. " if row_number is not None else ""
    observed_label = (
        "Observed trait value" if trait == "Trait value" else f"Observed {trait}"
    )
    ax.set_title(
        prefix + (observed_label if ticks else f"{trait}: all missing"),
        fontsize=8,
        loc="left",
        pad=3,
    )
    ax.spines[:].set_visible(False)
    if missing:
        from matplotlib.patches import Rectangle

        ax.add_patch(
            Rectangle(
                (1.08, 0),
                0.06,
                1,
                transform=ax.transAxes,
                color=_MISSING,
                clip_on=False,
            )
        )
        ax.text(1.17, 0.5, "Missing", transform=ax.transAxes, va="center", fontsize=7)


def draw_tip_heatmap(parent, table, tips, positions, label_depth):
    from matplotlib import colormaps
    from matplotlib.patches import Rectangle

    traits, values = tip_trait_values(table, tips)
    ax = _below_axis(parent, _GAP, _CELL_HEIGHT * len(traits))
    ax.set_label("observed-tip-heatmap")
    ax.set_xlim(parent.get_xlim())
    ax.set_ylim(len(traits), 0)
    ax.set_axis_off()
    cmap = colormaps["viridis"]
    for index, (trait, row) in enumerate(zip(traits, values, strict=True)):
        if len(traits) > 1:
            ax.text(
                parent.get_xlim()[0] + 0.05,
                index + 0.5,
                str(index + 1),
                va="center",
                fontsize=7,
            )
        norm, ticks = _trait_scale(row)
        for tip, value in zip(tips, row, strict=True):
            color = _MISSING if not np.isfinite(value) else cmap(norm(value))
            cell = Rectangle(
                (positions[tip] - 0.5, index),
                1,
                1,
                facecolor=color,
                edgecolor="white",
                linewidth=0.5,
            )
            cell.set_gid(f"tip-heatmap:{trait}:{tip.name}")
            ax.add_patch(cell)
        offset = heatmap_height(len(traits)) + label_depth + 0.22 + index * _KEY_HEIGHT
        _draw_key(
            parent,
            offset,
            trait,
            norm,
            ticks,
            cmap,
            bool(np.isnan(row).any()),
            index + 1 if len(traits) > 1 else None,
        )
    return ax
