"""DTT and conditional Brownian-null plots without pyplot state."""

import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.text import Text

from nwkit.disparity import transform_traits
from nwkit.evolution import tree_depths
from nwkit.trait_input import parse_trait_columns


def _tree_branches(ax, tree):
    """Draw the shared chronology and return tip positions."""
    tips = list(tree.leaves())
    depths = tree_depths(tree, allow_zero=True)
    height = max(depths[tip] for tip in tips)
    y = {tip: float(len(tips) - i - 1) for i, tip in enumerate(tips)}
    for node in tree.traverse("postorder"):
        if not node.is_leaf:
            y[node] = (y[node.children[0]] + y[node.children[-1]]) / 2
            ax.plot(
                [depths[node] / height] * 2,
                [y[node.children[0]], y[node.children[-1]]],
                color="#a3adb6",
                linewidth=1,
            )
        if not node.is_root:
            ax.plot(
                [depths[node.up] / height, depths[node] / height],
                [y[node], y[node]],
                color="#a3adb6",
                linewidth=1,
            )
    ax.set(
        xlim=(0, 1),
        ylim=(-0.7, len(tips) - 0.3),
        yticks=[],
        xticks=np.linspace(0, 1, 6),
        xlabel="Relative crown time",
    )
    ax.spines["left"].set_visible(False)
    return tips, depths, height, y


def _raw_color_values(values):
    """Map to [0, 1] before Matplotlib can expand tiny raw ranges to +/-0.1."""
    lower, upper = float(np.min(values)), float(np.max(values))
    with np.errstate(over="ignore"):
        span = upper - lower
    if np.isfinite(span):
        normalized = (values - lower) / span
    else:
        # Opposite-sign traits may span more than the largest finite float.
        unit = max(abs(lower), abs(upper))
        normalized = (values / unit - lower / unit) / (upper / unit - lower / unit)
    return normalized, lower, upper


def _original_unit_bar(figure, artist, axis, lower, upper, label):
    ticks = np.array([0.0, 0.5, 1.0])
    # The midpoint expression is finite even if upper-lower overflows.
    originals = np.array([lower, lower / 2 + upper / 2, upper])
    if len(np.unique(originals)) < 3:
        ticks, originals = ticks[[0, 2]], originals[[0, 2]]
    labels = []
    for precision in range(6, 18):
        labels = [f"{value:.{precision}g}" for value in originals]
        if len(set(labels)) == len(labels):
            break
    bar = figure.colorbar(artist, cax=axis, orientation="horizontal")
    bar.set_ticks(ticks, labels=labels)
    bar.set_label(label, fontsize=8)
    bar.ax.tick_params(labelsize=8)


def _trait_tree(ax, color_axis, tree, names, values, column):
    """Show observed tip values in original units; branches carry no inferred states."""
    tips, depths, height, y = _tree_branches(ax, tree)
    lookup = dict(zip(names, values, strict=True))
    ordered = np.array([lookup[tip.name] for tip in tips])
    color_values, lower, upper = _raw_color_values(ordered)
    dots = ax.scatter(
        [depths[tip] / height for tip in tips],
        [y[tip] for tip in tips],
        c=color_values,
        cmap="viridis",
        vmin=0,
        vmax=1,
        s=45 if len(tips) <= 40 else 9,
        edgecolors="white",
        linewidths=0.5,
        zorder=3,
        clip_on=False,
    )
    if len(tips) <= 40:
        for tip, value in zip(tips, ordered, strict=True):
            ax.annotate(
                f"{tip.name}  {value:.4g}",
                (depths[tip] / height, y[tip]),
                xytext=(8, 0),
                textcoords="offset points",
                va="center",
                fontsize=9,
                parse_math=False,
            )
    ax.set_title(column, loc="left", fontsize=13, fontweight="bold", parse_math=False)
    ax.set(
        xlim=(0, 1),
        ylim=(-0.7, len(tips) - 0.3),
        yticks=[],
        xticks=np.linspace(0, 1, 6),
        xlabel="Relative crown time",
    )
    ax.spines["left"].set_visible(False)
    _original_unit_bar(
        ax.figure, dots, color_axis, lower, upper, "Observed tip value · original units"
    )


def _trait_heatmap(tree_ax, heat_ax, color_ax, tree, names, values, columns, scale):
    tips, depths, height, y = _tree_branches(tree_ax, tree)
    tree_ax.set_title("Retained tree", loc="left", fontsize=13, fontweight="bold")
    if len(tips) <= 40:
        for tip in tips:
            tree_ax.annotate(
                tip.name,
                (depths[tip] / height, y[tip]),
                xytext=(8, 0),
                textcoords="offset points",
                va="center",
                fontsize=9,
                parse_math=False,
            )
    if scale == "standardize":
        values, _, _ = transform_traits(values, "standardize")
        values = values - values.mean(axis=0)
    lookup = dict(zip(names, values, strict=True))
    ordered = np.array([lookup[tip.name] for tip in tips])
    if scale == "standardize":
        limit = float(np.max(np.abs(ordered)))
        settings = {"cmap": "RdBu_r", "vmin": -limit, "vmax": limit}
    else:
        ordered, lower, upper = _raw_color_values(ordered)
        settings = {"cmap": "viridis", "vmin": 0.0, "vmax": 1.0}
    heat = heat_ax.imshow(
        ordered,
        aspect="auto",
        interpolation="nearest",
        extent=(-0.5, len(columns) - 0.5, -0.5, len(tips) - 0.5),
        origin="upper",
        **settings,
    )
    heat_ax.set(ylim=tree_ax.get_ylim(), yticks=[], xticks=np.arange(len(columns)))
    heat_ax.set_xticklabels(
        columns,
        rotation=90 if len(columns) > 20 else 45,
        ha="center" if len(columns) > 20 else "right",
        fontsize=9,
    )
    for label in heat_ax.get_xticklabels():
        label.set_parse_math(False)
    heat_ax.set_title("Observed traits", loc="left", fontsize=13, fontweight="bold")
    if scale == "standardize":
        bar = heat_ax.figure.colorbar(heat, cax=color_ax, orientation="horizontal")
        bar.set_label("Display z score · per trait", fontsize=9)
    else:
        _original_unit_bar(
            heat_ax.figure,
            heat,
            color_ax,
            lower,
            upper,
            "Observed value · shared original units",
        )


def draw_dtt(
    curve,
    summary,
    null_mdi,
    path,
    format_name,
    columns,
    tree,
    names,
    values,
    *,
    figure_layout="heatmap",
    figure_columns=None,
    figure_scale="standardize",
):
    simulated = len(null_mdi) > 0
    calculated_count = len(columns)
    displayed = (
        parse_trait_columns(figure_columns) if figure_columns is not None else columns
    )
    values = values[:, [columns.index(column) for column in displayed]]
    columns = displayed
    heatmap = figure_layout != "trees"
    tree_height = min(10, max(2.5, len(names) * 0.20 + 0.7))
    heights = (
        [tree_height, 0.4, 5.2] if heatmap else [tree_height] * len(columns) + [5.2]
    )
    right_width = max(3.4, min(12, len(columns) * 0.45)) if heatmap else 3.4
    width, height = 7.8 + right_width, sum(heights)
    if width * height * 170**2 > 32_000_000:
        raise ValueError(
            "DTT figure exceeds 32 million pixels; use --figure-layout heatmap, fewer --figure-columns or omit --figure-out."
        )
    figure = Figure(figsize=(width, height), layout="constrained")
    FigureCanvasAgg(figure)
    grid = figure.add_gridspec(
        len(heights), 2, height_ratios=heights, width_ratios=[7.8, right_width]
    )
    ax = figure.add_subplot(grid[-1, 0])
    if heatmap:
        tree_ax = figure.add_subplot(grid[0, 0], sharex=ax)
        _trait_heatmap(
            tree_ax,
            figure.add_subplot(grid[0, 1], sharey=tree_ax),
            figure.add_subplot(grid[1, 1]),
            tree,
            names,
            values,
            columns,
            figure_scale,
        )
    else:
        for index, column in enumerate(columns):
            color_grid = grid[index, 1].subgridspec(3, 1, height_ratios=[1, 0.10, 1])
            _trait_tree(
                figure.add_subplot(grid[index, 0], sharex=ax),
                figure.add_subplot(color_grid[1]),
                tree,
                names,
                values[:, index],
                column,
            )
    distribution = figure.add_subplot(grid[-1, 1]) if simulated else None
    x, observed = curve.relative_time.to_numpy(), curve.relative_disparity.to_numpy()
    color, null_color = "#087e8b", "#65758b"
    if simulated:
        level = f"{100 * summary.envelope_level:g}%"
        ax.fill_between(
            x,
            curve.bm_lower.to_numpy(),
            curve.bm_upper.to_numpy(),
            color="#c6d3df",
            alpha=0.65,
            linewidth=0,
            label=f"{level} pointwise BM envelope",
        )
        ax.plot(
            x,
            curve.bm_median.to_numpy(),
            color=null_color,
            linewidth=1.7,
            linestyle="--",
            label="BM median",
        )
    ax.plot(x, observed, color=color, linewidth=2.5, label="Observed DTT", zorder=3)
    ax.set(
        xlim=(0, 1),
        xlabel="Relative time  ·  crown root → present",
        ylabel="Mean relative within-clade disparity",
    )
    ax.set_ylim(bottom=0)
    ax.set_title("Disparity through time", loc="left", fontsize=15, fontweight="bold")
    ax.legend(frameon=False, loc="best", fontsize=9)
    ax.grid(axis="y", color="#e7ebef", linewidth=0.7)
    ax.set_axisbelow(True)
    if distribution is not None:
        bins = min(35, max(5, int(np.sqrt(len(null_mdi)))))
        distribution.hist(
            null_mdi, bins=bins, color="#c6d3df", edgecolor="white", linewidth=0.6
        )
        distribution.axvline(
            summary.mdi, color=color, linewidth=2.4, label="Observed MDI"
        )
        distribution.axvline(0, color=null_color, linewidth=0.8, linestyle=":")
        distribution.set(xlabel="Area relative to BM median", ylabel="BM simulations")
        distribution.set_title(
            f"MDI = {summary.mdi:.3f}", loc="left", fontsize=13, fontweight="bold"
        )
        distribution.legend(frameon=False, fontsize=9)
    for panel in figure.axes:
        for side in ("top", "right"):
            panel.spines[side].set_visible(False)
        panel.spines["left"].set_color("#9aa5af")
        panel.spines["bottom"].set_color("#9aa5af")
    description = f"{int(summary.num_tips)} tips · {calculated_count} trait{'s' if calculated_count != 1 else ''} analyzed · {len(columns)} displayed · {int(summary.num_simulations)} BM simulations"
    figure.suptitle(description, x=0.01, ha="left", fontsize=10, color="#52616d")
    figure.supxlabel(
        "BM comparison conditions on fitted rates and tree; the envelope is not a confidence band."
        if simulated
        else "Observed disparity only; singleton lineages are excluded from the clade average.",
        fontsize=9,
        color="#52616d",
    )
    from nwkit.asr_compare_figure import _font_family_for_text

    visible_names = names if len(names) <= 40 else []
    try:
        family = _font_family_for_text(" ".join([*columns, *visible_names]))
    except ValueError as exc:
        raise ValueError(str(exc).replace("ASR comparison PDF", "DTT figure")) from exc
    for label in figure.findobj(Text):
        # User labels are literal text, independent of the caller's TeX settings.
        label.set_usetex(False)
        if family:
            label.set_fontfamily(family)
    figure.savefig(path, format=format_name.lower(), dpi=170)
