"""Independent references and adversarial convergence cases."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from nwkit.mcmc_diagnostics import diagnose
from nwkit.threshold_diagnostics import (
    diagnose_threshold_draws,
    summarize_diagnostics,
    trace_storage,
)
from nwkit.util import read_tree
from tests.threshold_diagnostic_support import diagnostic_cases

REFERENCE = json.loads(
    (Path(__file__).parent / "data/threshold_diagnostics.json").read_text()
)


@pytest.mark.parametrize("name,values", diagnostic_cases().items())
def test_matches_independent_r_posterior(name, values):
    result = diagnose(values)
    for metric, expected in REFERENCE["cases"][name].items():
        assert result[metric] == pytest.approx(expected, rel=1e-6, abs=1e-8)


@pytest.mark.parametrize("name", ["drift", "location", "scale"])
def test_nonstationarity_and_scale_mismatch_do_not_pass(name):
    assert "mcmc_rhat" in diagnose(diagnostic_cases()[name])["status"]


def test_multilag_and_antithetic_ess():
    cases = diagnostic_cases()
    assert diagnose(cases["lag_two"])["ess_bulk"] < 400
    # Negative correlation can legitimately provide ESS greater than draw count.
    assert diagnose(cases["antithetic"])["ess_mean"] > cases["antithetic"].size


@pytest.mark.parametrize(
    "shape,status",
    [
        ((1, 100), "mcmc_rhat_unavailable"),
        ((4, 7), "insufficient_draws"),
        ((4, 100), "constant_trace"),
    ],
)
def test_unavailable_is_never_success(shape, status):
    result = diagnose(np.zeros(shape))
    assert result["status"] == status
    assert np.isnan(result["rhat"])
    assert np.isnan(result["ess_bulk"])


def test_constants_are_classified_by_model_not_draws():
    x = np.zeros((4, 100))
    assert diagnose(x, structural=True)["status"] == "structural_constant"
    assert diagnose(x, indicator=True)["status"] == "unresolved_rare_category"
    x[1] = 1
    assert diagnose(x)["status"] == "constant_trace"
    x[0, 0] = np.nan
    assert diagnose(x)["status"] == "nonfinite_trace"


def test_one_stuck_chain_and_stuck_half_are_not_ignored():
    x = diagnostic_cases()["iid"].copy()
    x[0, :500] = 0
    assert diagnose(x)["status"] == "constant_trace"


def test_binary_mcse_uses_probability_ess_not_continuous_tail_quantiles():
    x = (diagnostic_cases()["iid"] > 0).astype(float)
    result = diagnose(x, indicator=True)
    assert result["status"] == "ok"
    assert np.isnan(result["ess_tail"])
    assert result["mcse_mean"] == pytest.approx(
        np.sqrt(np.var(x, ddof=1) / result["ess_mean"])
    )
    assert result["mcse_mean"] < 0.01


def test_memory_and_disk_storage_match_and_cleanup_on_error():
    for limit in (0, 10000):
        with trace_storage((2, 10, 3), memory_limit=limit) as traces:
            traces[:] = 2
            assert traces.sum() == 120
            backing = traces._mmap if isinstance(traces, np.memmap) else None
        if backing is not None:
            assert backing.closed
    with (
        pytest.raises(RuntimeError),
        trace_storage((2, 10, 3), memory_limit=0) as traces,
    ):
        backing = traces._mmap
        raise RuntimeError("interrupted sampling")
    assert backing.closed


def test_nonroot_stall_and_structural_categories_are_visible():
    tree = read_tree("((A:1,B:1)I:1,C:1)R;", "1", True, quiet=True, rooted="yes")
    nodes = list(tree.traverse("preorder"))
    rng = np.random.default_rng(10)
    traces = rng.normal(size=(4, 1000, len(nodes) + 1))
    traces[:, :, -1] = 0  # fixed binary threshold
    traces[:, :, nodes.index(tree["I"])] = 1  # root remains well mixed
    tip_index = nodes.index(tree["A"])
    traces[:, :, tip_index] = -np.abs(traces[:, :, tip_index])
    table = diagnose_threshold_draws(
        tree,
        nodes,
        ("low", "high"),
        {tip_index: np.array([0])},
        traces,
        estimated=False,
    )
    assert (
        table[(table.name == "R") & (table.variable == "liability")].status.iloc[0]
        == "ok"
    )
    assert "constant_trace" in summarize_diagnostics(table)["fit_status"]
    assert set(table[(table.name == "A") & (table.variable == "category")].status) == {
        "structural_constant"
    }
    assert set(table[table.variable == "threshold"].status) == {"structural_constant"}
    assert len(table[table.variable == "liability"]) == len(nodes)


@pytest.mark.parametrize("values", [np.zeros(3), np.empty((0, 2)), np.empty((2, 0))])
def test_invalid_trace_shape(values):
    with pytest.raises(ValueError, match="chains-by-draws"):
        diagnose(values)


def test_folded_constant_with_different_scales_does_not_pass():
    x = np.tile([-1.0, 1.0], (4, 500)) * np.arange(1, 5)[:, None]
    with np.errstate(divide="raise", invalid="raise"):
        result = diagnose(x)
    assert result["rhat"] == np.inf
    assert "mcmc_rhat" in result["status"]


def test_reference_inputs_match_exported_hash():
    path = Path(__file__).parent / "data/threshold_diagnostic_draws.npz"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == REFERENCE["draws_sha256"]
