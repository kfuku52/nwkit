"""Storage and variable-level diagnostics for threshold posterior draws."""

import math
import tempfile
from contextlib import contextmanager

import numpy as np
import pandas as pd

from nwkit.mcmc_diagnostics import diagnose
from nwkit.util import assign_branch_ids


@contextmanager
def trace_storage(shape, *, memory_limit=64 * 1024**2):
    """Bound resident trace storage; disk-backed traces are always temporary."""
    if math.prod(shape) * 8 <= memory_limit:
        yield np.empty(shape, dtype=float)
    else:
        with tempfile.TemporaryFile() as handle:
            handle.truncate(math.prod(shape) * 8)
            traces = np.memmap(handle, dtype=float, mode="r+", shape=shape)
            try:
                yield traces
            finally:
                traces._mmap.close()


def diagnose_threshold_draws(tree, nodes, states, constraints, traces, *, estimated):
    """Inspect every node, including imputed tips, and all free thresholds."""
    ids = assign_branch_ids(tree)
    count = len(nodes)
    thresholds = traces[:, :, count:]
    rows = []

    def append(values, node, kind, label="", structural=False):
        rows.append(
            {
                "branch_id": ids[node] if node is not None else -1,
                "name": str(node.name or "") if node is not None else "",
                "variable": kind,
                "state_or_threshold": label,
                **diagnose(values, indicator=kind == "category", structural=structural),
            }
        )

    for index, node in enumerate(nodes):
        values = traces[:, :, index]
        append(values, node, "liability")
        # Second raw moment determines the reported variance along with mean.
        append(values * values, node, "liability_second_moment")
        categories = np.zeros(values.shape, dtype=int)
        for threshold_index in range(thresholds.shape[2]):
            categories += values > thresholds[:, :, threshold_index]
        allowed = constraints.get(index)
        for category, state in enumerate(states):
            structural = allowed is not None and (
                len(allowed) == 1 or category not in allowed
            )
            append(categories == category, node, "category", state, structural)
    for index in range(len(states) - 1):
        append(
            thresholds[:, :, index],
            None,
            "threshold",
            str(index + 1),
            structural=not estimated or index == 0,
        )
    return pd.DataFrame(rows)


def summarize_diagnostics(table):
    active = table[table.status != "structural_constant"]
    reasons = sorted(
        {
            reason
            for status in active.status
            for reason in status.split("+")
            if reason != "ok"
        }
    )

    def extreme(column, operation):
        values = active[column].dropna()
        return float(getattr(values, operation)()) if len(values) else math.nan

    bulk = extreme("ess_bulk", "min")
    tail = extreme("ess_tail", "min")
    ess_values = [v for v in (bulk, tail) if math.isfinite(v)]
    probabilities = active[active.variable == "category"]
    precision = probabilities.mcse_mean.dropna()
    return {
        "rhat_max": extreme("rhat", "max"),
        "ess_min": min(ess_values) if ess_values else math.nan,
        "ess_bulk_min": bulk,
        "ess_tail_min": tail,
        "probability_mcse_max": float(precision.max()) if len(precision) else math.nan,
        "fit_status": "+".join(reasons) if reasons else "ok",
    }
