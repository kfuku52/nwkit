"""History segments, duration summaries, time bins and branch probabilities."""

import numpy as np
import pandas as pd

from nwkit.util import get_node_class

HISTORY_COLUMNS = (
    "simulation",
    "branch_id",
    "parent",
    "node_class",
    "name",
    "segment",
    "state",
    "start",
    "end",
    "duration",
    "start_from_root",
    "end_from_root",
)
SUMMARY_COLUMNS = (
    "branch_id",
    "parent",
    "node_class",
    "name",
    "state",
    "branch_length",
    "total_duration",
    "mean_duration",
    "duration_sd",
    "duration_q025",
    "duration_q975",
    "mean_fraction",
    "num_simulations",
)
PROBABILITY_COLUMNS = (
    "branch_id",
    "parent",
    "node_class",
    "name",
    "position",
    "distance",
    "time_from_root",
    "state",
    "probability",
    "mc_se",
    "num_simulations",
)
TIME_COLUMNS = (
    "bin",
    "time_start",
    "time_end",
    "lineage_time",
    "quantity",
    "state",
    "other_state",
    "total",
    "mean",
    "per_lineage_time",
    "num_simulations",
)


def _branches(sample):
    return [node for node in sample.branch_ids if not node.is_root]


def _metadata(sample, node):
    return {
        "branch_id": sample.branch_ids[node],
        "parent": sample.branch_ids[node.up],
        "node_class": get_node_class(node),
        "name": str(node.name or ""),
    }


def _root_times(depth, local_times):
    local = np.asarray(local_times, dtype=float)
    times = depth + local
    if np.any((np.diff(local) > 0) & (np.diff(times) <= 0)):
        raise ValueError(
            "Map times are not distinguishable on the root-time axis; "
            "reduce the root-depth/branch-time scale disparity or use duration/time-bin output."
        )
    return times


def history_table(sample):
    rows = []
    nodes = {identifier: node for node, identifier in sample.branch_ids.items()}
    for draw in sample.draws:
        for branch, segments in draw.branches.items():
            node = nodes[branch]
            depth = sample.depths[node.up]
            metadata = _metadata(sample, node)
            root_times = _root_times(
                depth, [item[0] for item in segments] + [segments[-1][1]]
            )
            for index, (start, end, state) in enumerate(segments):
                rows.append(
                    {
                        "simulation": draw.simulation,
                        **metadata,
                        "segment": index + 1,
                        "state": sample.states[state],
                        "start": start,
                        "end": end,
                        "duration": end - start,
                        "start_from_root": root_times[index],
                        "end_from_root": root_times[index + 1],
                    }
                )
    return pd.DataFrame(rows, columns=HISTORY_COLUMNS)


def duration_table(sample):
    rows = []
    count = len(sample.draws)
    for node in _branches(sample):
        branch = sample.branch_ids[node]
        length = sample.lengths[branch]
        values = np.zeros((count, len(sample.states)))
        for index, draw in enumerate(sample.draws):
            for start, end, state in draw.branches[branch]:
                if length:
                    values[index, state] += (end - start) / length
        for state, label in enumerate(sample.states):
            durations = values[:, state]
            low, high = np.quantile(durations, [0.025, 0.975]) * length
            mean_fraction = float(durations.mean())
            total = float(durations.sum()) * length
            if not np.isfinite(total):
                raise ValueError(
                    "Map total duration overflows; reduce --n-sim or branch lengths."
                )
            mean = mean_fraction * length
            rows.append(
                {
                    **_metadata(sample, node),
                    "state": label,
                    "branch_length": sample.lengths[branch],
                    "total_duration": total,
                    "mean_duration": mean,
                    "duration_sd": float(durations.std(ddof=1)) * length
                    if count > 1
                    else "",
                    "duration_q025": low,
                    "duration_q975": high,
                    "mean_fraction": mean_fraction if length else "",
                    "num_simulations": count,
                }
            )
    return pd.DataFrame(rows, columns=SUMMARY_COLUMNS)


def probability_table(sample, points):
    rows = []
    fractions = np.linspace(0, 1, points)
    count = len(sample.draws)
    for node in _branches(sample):
        branch = sample.branch_ids[node]
        distances = fractions * sample.lengths[branch]
        root_times = _root_times(sample.depths[node.up], distances)
        frequencies = np.zeros((points, len(sample.states)), dtype=int)
        for draw in sample.draws:
            segments = draw.branches[branch]
            starts = np.asarray([item[0] for item in segments])
            states = np.asarray([item[2] for item in segments])
            selected = states[np.searchsorted(starts, distances, side="right") - 1]
            frequencies[np.arange(points), selected] += 1
        for index, position in enumerate(fractions):
            for state, label in enumerate(sample.states):
                probability = frequencies[index, state] / count
                rows.append(
                    {
                        **_metadata(sample, node),
                        "position": position,
                        "distance": distances[index],
                        "time_from_root": root_times[index],
                        "state": label,
                        "probability": probability,
                        "mc_se": np.sqrt(probability * (1 - probability) / count),
                        "num_simulations": count,
                    }
                )
    return pd.DataFrame(rows, columns=PROBABILITY_COLUMNS)


def time_table(sample, bins):
    segments = sum(
        len(branch) for draw in sample.draws for branch in draw.branches.values()
    )
    if segments * bins > 20_000_000:
        raise ValueError(
            "Map time-bin aggregation exceeds 20,000,000 segment/bin intersections; reduce --map-time-bins or --n-sim."
        )
    height = max(sample.depths.values())
    if height <= 0:
        raise ValueError("--map-time-out requires a positive tree time span.")
    edges = np.linspace(0, height, bins + 1)
    if np.any(np.diff(edges) <= 0):
        raise ValueError(
            "Time bins are not numerically distinct; reduce --map-time-bins."
        )
    p, count = len(sample.states), len(sample.draws)
    duration = np.zeros((bins, p))
    transitions = np.zeros((bins, p, p), dtype=int)
    exposure = np.zeros(bins)
    for node in _branches(sample):
        branch, depth = sample.branch_ids[node], sample.depths[node.up]
        local_edges = edges - depth
        exposure += np.maximum(
            0,
            np.minimum(local_edges[1:], sample.lengths[branch])
            - np.maximum(local_edges[:-1], 0),
        )
        for draw in sample.draws:
            segments = draw.branches[branch]
            for index, (start, end, state) in enumerate(segments):
                duration[:, state] += np.maximum(
                    0,
                    np.minimum(local_edges[1:], end)
                    - np.maximum(local_edges[:-1], start),
                )
                if index:
                    bin_index = min(
                        bins - 1,
                        int(np.searchsorted(local_edges, start, side="right") - 1),
                    )
                    transitions[bin_index, segments[index - 1][2], state] += 1
    if not np.isfinite(duration).all() or not np.isfinite(exposure).all():
        raise ValueError(
            "Map time-bin totals overflow; reduce --n-sim or branch lengths."
        )
    rows = []
    for index in range(bins):
        base = {
            "bin": index + 1,
            "time_start": edges[index],
            "time_end": edges[index + 1],
            "lineage_time": exposure[index],
            "num_simulations": count,
        }
        for first, state in enumerate(sample.states):
            rows.append(
                _time_row(base, "duration", state, "", duration[index, first], count)
            )
            for second, other in enumerate(sample.states):
                if first != second:
                    rows.append(
                        _time_row(
                            base,
                            "transitions",
                            state,
                            other,
                            int(transitions[index, first, second]),
                            count,
                        )
                    )
    return pd.DataFrame(rows, columns=TIME_COLUMNS)


def _time_row(base, quantity, state, other, total, count):
    mean = total / count
    exposure = float(base["lineage_time"])
    rate = float(mean) / exposure if exposure else ""
    if rate != "" and not np.isfinite(rate):
        raise ValueError("Map time-bin rate overflows; rescale branch lengths.")
    return {
        **base,
        "quantity": quantity,
        "state": state,
        "other_state": other,
        "total": total,
        "mean": mean,
        "per_lineage_time": rate,
    }
