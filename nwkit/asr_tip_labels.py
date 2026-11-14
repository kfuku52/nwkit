"""Detached, value-ordered tip labels for continuous trait panels."""

import numpy as np

from nwkit.asr_heatmap import _below_axis


def draw_trait_tip_labels(parent, values, colors, *, first_history=False):
    """Fan endpoints into evenly spaced labels without covering the data panel.

    The top of the detached strip repeats tip values on the parent's x scale.
    Sorting by value keeps the fan connectors from crossing. Full names are
    retained, so abbreviations cannot make two distinct tips ambiguous.
    """
    names = sorted(values, key=lambda name: (values[name], name))
    band = _below_axis(parent, 0.62, 0.32)
    band.set_label("trait-tip-labels")
    band.set_axis_off()
    band.set_ylim(0, 1)
    band.set_xlim(parent.get_xlim())

    def sync_limits(ax):
        band.set_xlim(ax.get_xlim())

    parent.callbacks.connect("xlim_changed", sync_limits)
    for name, slot in zip(names, np.linspace(0.04, 0.96, len(names)), strict=True):
        color = colors[name]
        band.annotate(
            name,
            xy=(values[name], 1),
            xycoords="data",
            xytext=(slot, 0),
            textcoords="axes fraction",
            rotation=90,
            ha="center",
            va="top",
            fontsize=7,
            color=color,
            annotation_clip=False,
            arrowprops={
                "arrowstyle": "-",
                "color": color,
                "alpha": 0.45,
                "linewidth": 0.6,
                "shrinkA": 2,
                "shrinkB": 0,
                "relpos": (0.5, 1),
            },
        ).set_gid(f"trait-tip-label:{name}")
        band.plot(values[name], 1, "o", color=color, markersize=2, clip_on=False)
    band.text(
        0,
        1.2,
        "Tip labels (first history)"
        if first_history
        else "Tip labels · ordered by value",
        transform=band.transAxes,
        fontsize=7,
        color="#666666",
    )
    return band
