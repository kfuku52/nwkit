"""CLI contracts, transactional exports and probability-ribbon figures."""

import io
import sys
from contextlib import redirect_stdout
from copy import copy
from pathlib import Path

import numpy as np

from nwkit.output_transaction import output_transaction
from nwkit.stochastic_history import sample_histories
from nwkit.stochastic_map_tables import (
    duration_table,
    history_table,
    probability_table,
    time_table,
)

MAP_OUTPUTS = (
    "map_history_out",
    "map_summary_out",
    "map_time_out",
    "map_probabilities_out",
    "map_figure_out",
)


def extended_maps_requested(args):
    return any(getattr(args, name, None) not in (None, "") for name in MAP_OUTPUTS)


def any_maps_requested(args):
    return extended_maps_requested(args) or getattr(
        args, "stochastic_map_out", None
    ) not in (None, "")


def validate_map_options(args):
    requested = extended_maps_requested(args)
    if requested and getattr(args, "model", None) in {"MK-MIXTURE", "THRESHOLD"}:
        raise ValueError(
            "Stochastic-map histories are not defined for MK-MIXTURE or THRESHOLD."
        )
    if requested and getattr(args, "tree_ensemble", None):
        raise ValueError(
            "Stochastic-map histories currently require a single input tree."
        )
    for name, output, minimum, maximum in (
        ("map_time_bins", "map_time_out", 1, 10000),
        ("map_grid_points", "map_probabilities_out", 2, 10001),
    ):
        value = getattr(args, name, None)
        if value is None:
            continue
        if not getattr(args, output, None) and not (
            name == "map_grid_points" and getattr(args, "map_figure_out", None)
        ):
            raise ValueError(f"--{name.replace('_', '-')} requires its map output.")
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, np.integer))
            or not minimum <= value <= maximum
        ):
            raise ValueError(
                f"--{name.replace('_', '-')} must be an integer between {minimum} and {maximum}."
            )
    figure = getattr(args, "map_figure_out", None)
    if (
        figure
        and not getattr(args, "_map_figure_format", None)
        and Path(figure).suffix.lower() not in {".png", ".pdf", ".svg"}
    ):
        raise ValueError("--map-figure-out requires PNG, PDF or SVG.")
    if requested and getattr(args, "outfile", None) in (None, ""):
        raise ValueError("--outfile must be a path or '-' for standard output.")


def run_map_transaction(args, handler):
    """One transaction includes ordinary ASR and all optional map outputs."""
    from nwkit.asr import _validate_asr_output_paths

    _validate_asr_output_paths(args)
    paths = {
        name: path
        for name, path in vars(args).items()
        if (name == "outfile" or name.endswith("_out")) and path not in (None, "", "-")
    }
    staged_args = copy(args)
    staged_args._map_output_staged = True
    if getattr(args, "map_figure_out", None):
        staged_args._map_figure_format = Path(args.map_figure_out).suffix.lower()[1:]
    captured = io.StringIO()
    with output_transaction(paths.values()) as staged:
        for name, path in paths.items():
            setattr(staged_args, name, staged[path])
        with redirect_stdout(captured):
            handler(staged_args)
    sys.stdout.write(captured.getvalue())


def _validate_table_sizes(tree, states, args, count):
    branches = sum(not node.is_root for node in tree.traverse())
    points = getattr(args, "map_grid_points", None) or 51
    bins = getattr(args, "map_time_bins", None) or 20
    p = len(states)
    if getattr(args, "map_probabilities_out", None) or getattr(
        args, "map_figure_out", None
    ):
        if branches * points * count > 20_000_000:
            raise ValueError(
                "Map probabilities exceed 20,000,000 draw/grid evaluations; "
                "reduce --n-sim or --map-grid-points."
            )
        if branches * points * p > 500_000:
            raise ValueError(
                "Map probability output exceeds 500,000 rows; reduce --map-grid-points."
            )
    if getattr(args, "map_time_out", None) and bins * p * p > 500_000:
        raise ValueError(
            "Map time output exceeds 500,000 rows; reduce --map-time-bins."
        )
    if branches * p > 500_000:
        raise ValueError(
            "Map duration output exceeds 500,000 rows; reduce tree/state size."
        )
    if getattr(args, "map_figure_out", None) and (
        p > 20 or sum(1 for _ in tree.leaves()) > 200
    ):
        raise ValueError(
            "Map figures support at most 20 states and 200 tips; use the history/probability tables for larger data."
        )
    return points, bins


def write_extended_maps(tree, states, fit, args):
    import pandas as pd

    from nwkit.asr import (
        _merge_simulation_counts,
        _stochastic_map_rows,
        _validated_simulation_threads,
        _write_table,
    )

    count, threads = _validated_simulation_threads(
        100 if getattr(args, "n_sim", None) is None else args.n_sim,
        1 if getattr(args, "threads", None) is None else args.threads,
    )
    points, bins = _validate_table_sizes(tree, states, args, count)
    sample = sample_histories(
        tree,
        states,
        fit,
        count,
        seed=getattr(args, "seed", None),
        threads=threads,
    )
    builders = {
        "map_history_out": lambda: history_table(sample),
        "map_summary_out": lambda: duration_table(sample),
        "map_time_out": lambda: time_table(sample, bins),
    }
    for name, build in builders.items():
        if getattr(args, name, None):
            _write_table(build(), getattr(args, name))
    if getattr(args, "stochastic_map_out", None):
        total, any_count = _merge_simulation_counts(
            draw.counts for draw in sample.draws
        )
        rows = _stochastic_map_rows(
            tree, states, sample.branch_ids, total, any_count, len(sample.draws)
        )
        # Preserve count-only columns, including for a root-only tree.
        columns = (
            "branch_id",
            "parent",
            "node_class",
            "name",
            "from_state",
            "to_state",
            "total_count",
            "mean_count",
            "posterior_frequency",
            "num_simulations",
        )
        _write_table(pd.DataFrame(rows, columns=columns), args.stochastic_map_out)
    if getattr(args, "map_probabilities_out", None) or getattr(
        args, "map_figure_out", None
    ):
        probabilities = probability_table(sample, points)
        if getattr(args, "map_probabilities_out", None):
            _write_table(probabilities, args.map_probabilities_out)
        if getattr(args, "map_figure_out", None):
            draw_map_figure(tree, sample, probabilities, args)
    sys.stderr.write(
        f"Stochastic maps: {len(sample.draws)} conditional histories; time is measured forward from the input root. Rate/root-model parameter uncertainty is excluded.\n"
    )


def draw_map_figure(tree, sample, probabilities, args):
    from matplotlib import colormaps
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from matplotlib.patches import Patch

    from nwkit.asr_figure import _coordinates

    nodes, tips, positions, depths = _coordinates(tree)
    palette = colormaps["tab20"].colors
    colors = (*palette[::2], *palette[1::2])
    figure = Figure(
        figsize=(10, max(3.0, len(tips) * 0.28 + 1.6)), layout="constrained"
    )
    FigureCanvasAgg(figure)
    ax = figure.subplots()
    for node in nodes:
        y = positions[node]
        if not node.is_leaf:
            child_y = [positions[child] for child in node.children]
            ax.plot(
                [depths[node]] * 2,
                [min(child_y), max(child_y)],
                color=".55",
                linewidth=0.7,
                zorder=0,
            )
        if node.is_root:
            continue
        table = probabilities[probabilities.branch_id == sample.branch_ids[node]]
        x = (
            table[table.state == sample.states[0]]
            .sort_values("position")
            .time_from_root.to_numpy()
        )
        if sample.lengths[sample.branch_ids[node]] == 0:
            continue
        lower = np.full(len(x), y - 0.32)
        for index, state in enumerate(sample.states):
            state_rows = table[table.state == state].sort_values("position")
            upper = lower + 0.64 * state_rows.probability.to_numpy()
            ax.fill_between(x, lower, upper, color=colors[index], linewidth=0)
            lower = upper
    ax.set_yticks([positions[tip] for tip in tips], [str(tip.name) for tip in tips])
    ax.yaxis.tick_right()
    for label in ax.get_yticklabels():
        label.set_parse_math(False)
    ax.tick_params(axis="y", length=0, labelsize=9)
    ax.set_ylim(-0.7, max(positions.values()) + 0.7)
    ax.invert_yaxis()
    ax.set_xlabel("Time from input root (branch-length units)")
    ax.set_title(
        f"Conditional state probabilities · {len(sample.draws)} maps",
        loc="left",
        fontsize=12,
    )
    legend = ax.legend(
        handles=[
            Patch(color=colors[i], label=str(state))
            for i, state in enumerate(sample.states)
        ],
        loc="upper left",
        bbox_to_anchor=(0, -0.10),
        ncol=min(5, len(sample.states)),
        frameon=False,
    )
    for label in legend.get_texts():
        label.set_parse_math(False)
    for side in ("top", "left", "right"):
        ax.spines[side].set_visible(False)
    format_name = (
        getattr(args, "_map_figure_format", None)
        or Path(args.map_figure_out).suffix[1:]
    )
    figure.savefig(args.map_figure_out, format=format_name, dpi=160)
