import math

import pandas as pd
import pytest

from nwkit.asr_regimes import RegimeAssignment
from nwkit.cli import main
from nwkit.nonstationary_continuous_asr import compute_bm_drift_marginals
from nwkit.regime_drift_asr import (
    build_regime_drift_process,
    compute_regime_drift_marginals,
)
from nwkit.util import read_tree


def tree_from(source):
    return read_tree(source, "1", True, quiet=True, rooted="yes")


def assignment(tree, values):
    nodes = list(tree.traverse())
    return RegimeAssignment(
        tuple(dict.fromkeys(values)), dict(zip(nodes, values, strict=True)), "memory"
    )


def test_regime_drift_process_uses_incoming_branch_parameters():
    tree = tree_from("(A:2,B:3)R;")
    regimes = assignment(tree, ["root", "slow", "fast"])
    process = build_regime_drift_process(
        tree,
        regimes,
        sigma2_by_regime={"root": 1.0, "slow": 0.5, "fast": 2.0},
        drift_by_regime={"root": 0.0, "slow": -0.2, "fast": 0.4},
    )
    assert process.transitions[tree["A"]].intercept == pytest.approx(-0.4)
    assert process.transitions[tree["A"]].variance == pytest.approx(1.0)
    assert process.transitions[tree["B"]].intercept == pytest.approx(1.2)
    assert process.transitions[tree["B"]].variance == pytest.approx(6.0)


def test_equal_fixed_regimes_reduce_to_brownian_drift():
    tree = tree_from("((A:0.5,B:1.2):0.4,(C:0.7,D:1.5):0.2)R;")
    regimes = assignment(tree, ["cold", "cold", "warm", "cold", "warm", "warm", "cold"])
    values = {"A": -1.0, "B": 1.0, "C": 2.0, "D": None}
    errors = {"A": 0.1, "B": 0.2, "C": 0.3}
    expected, expected_fit = compute_bm_drift_marginals(
        tree,
        values,
        sigma2=1.2,
        drift=0.3,
        standard_errors=errors,
    )
    fixed = {regime: {"sigma2": 1.2, "drift": 0.3} for regime in regimes.regimes}
    observed, fit = compute_regime_drift_marginals(
        tree,
        values,
        regimes,
        regime_parameters=fixed,
        standard_errors=errors,
    )
    assert fit.restricted_log_likelihood == pytest.approx(
        expected_fit.restricted_log_likelihood, abs=2e-12
    )
    for node in tree.traverse():
        assert observed[node].mean == pytest.approx(expected[node].mean, abs=2e-12)
        assert observed[node].variance == pytest.approx(
            expected[node].variance, abs=2e-12
        )


def test_regime_drifts_reject_rank_deficient_tip_design():
    tree = tree_from("(A:1,B:1,C:1,D:1)R;")
    regimes = assignment(tree, ["root", "tip", "tip", "tip", "tip"])
    with pytest.raises(ValueError, match="drifts are not identifiable"):
        compute_regime_drift_marginals(
            tree,
            {"A": 0.0, "B": 1.0, "C": 2.0, "D": 3.0},
            regimes,
            sigma2=1.0,
        )


@pytest.mark.integration
def test_regime_drift_cli_fixed_parameters(tmp_path):
    traits = tmp_path / "traits.tsv"
    regimes = tmp_path / "regimes.tsv"
    parameters = tmp_path / "parameters.tsv"
    output = tmp_path / "asr.tsv"
    metadata = tmp_path / "model.tsv"
    traits.write_text("leaf_name\tvalue\nA\t0\nB\t1\nC\t3\nD\t4\n")
    regimes.write_text(
        "branch_id\tregime\n0\tlow\n1\tlow\n2\thigh\n3\tlow\n4\tlow\n5\thigh\n6\thigh\n"
    )
    parameters.write_text("regime\tsigma2\tdrift\nlow\t0.8\t-0.1\nhigh\t1.5\t0.3\n")
    main(
        [
            "asr",
            "-i",
            "((A:1,B:1):1,(C:1,D:1):1)R;",
            "--input-rooted",
            "yes",
            "--trait",
            str(traits),
            "--state-column",
            "value",
            "--model",
            "BMS-DRIFT",
            "--regime-map",
            str(regimes),
            "--regime-parameters",
            str(parameters),
            "--model-out",
            str(metadata),
            "-o",
            str(output),
        ]
    )
    assert len(pd.read_csv(output, sep="\t")) == 7
    row = pd.read_csv(metadata, sep="\t").iloc[0]
    assert row["model"] == "BMS-DRIFT"
    assert row["sigma2_low"] == pytest.approx(0.8)
    assert row["drift_high"] == pytest.approx(0.3)
    assert math.isfinite(row["restricted_log_likelihood"])
