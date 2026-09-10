"""Deterministic inputs shared by tests and the independent R reference exporter."""

from pathlib import Path

import numpy as np


def generate_diagnostic_cases():
    rng = np.random.default_rng(72)
    base = rng.normal(size=(4, 1000))
    lag_two = base.copy()
    antithetic = base.copy()
    for index in range(2, base.shape[1]):
        lag_two[:, index] += 0.9 * lag_two[:, index - 2]
        antithetic[:, index] -= 0.7 * antithetic[:, index - 1]
    return {
        "iid": base,
        "drift": base + np.linspace(-2, 2, 1000),
        "location": base + np.arange(4)[:, None],
        "scale": base * np.array([0.1, 1, 2, 4])[:, None],
        "heavy_tail": rng.standard_cauchy(size=base.shape),
        "lag_two": lag_two,
        "antithetic": antithetic,
        "ties": np.round(base),
        "odd": rng.normal(size=(4, 1001)),
    }


def diagnostic_cases():
    """Load exact reference inputs, independent of future RNG implementations."""
    with np.load(
        Path(__file__).parent / "data/threshold_diagnostic_draws.npz",
        allow_pickle=False,
    ) as archive:
        return {name: archive[name] for name in archive.files}
