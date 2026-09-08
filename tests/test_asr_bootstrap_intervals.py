import numpy as np
import pandas as pd
import pytest

from nwkit.asr_bootstrap_intervals import bootstrap_prediction_intervals
from nwkit.asr_continuous_diagnostics import _simulated_tip_data
from nwkit.cli import main
from nwkit.evolution import build_evolutionary_process
from nwkit.util import read_tree


def test_prediction_error_intervals_reproduce_conditional_gaussian_uncertainty():
    tree = read_tree("(A:1,B:1)R;", "1", True, quiet=True, rooted="yes")
    process = build_evolutionary_process(tree, model="brownian", root_mode="flat")
    observed = {"A": 0.0, "B": 1.0}

    def run():
        return bootstrap_prediction_intervals(
            process,
            observed,
            lambda seed: _simulated_tip_data(
                process, observed, None, seed, include_latent=True
            ),
            lambda values: process,
            num_simulations=1000,
            seed=101,
        )

    result = run()
    pd.testing.assert_frame_equal(result, run())
    root = result.loc[result.name == "R"].iloc[0]
    assert root.error_sd == pytest.approx(np.sqrt(0.5), rel=0.08)
    assert root.lower == pytest.approx(0.5 - 1.96 * np.sqrt(0.5), abs=0.15)
    assert root.upper == pytest.approx(0.5 + 1.96 * np.sqrt(0.5), abs=0.15)
    for tip in ("A", "B"):
        row = result.loc[result.name == tip].iloc[0]
        assert row.lower == pytest.approx(observed[tip], abs=1e-14)
        assert row.upper == pytest.approx(observed[tip], abs=1e-14)


@pytest.mark.integration
def test_bootstrap_intervals_cli(tmp_path):
    traits = tmp_path / "traits.tsv"
    traits.write_text("leaf_name\tvalue\nA\t0\nB\t1\nC\t3\nD\t4\n")
    out = tmp_path / "intervals.tsv"
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
            "--bootstrap-intervals-out",
            str(out),
            "--bootstrap-interval-simulations",
            "5",
            "--seed",
            "9",
            "-o",
            str(tmp_path / "asr.tsv"),
        ]
    )
    table = pd.read_csv(out, sep="\t")
    assert len(table) == 7
    assert set(table.num_simulations) == {5}
    assert (table.lower <= table.upper).all()


@pytest.mark.parametrize("count", [True, 1, 2.5])
def test_bootstrap_interval_counts(count):
    with pytest.raises(ValueError, match="at least two"):
        bootstrap_prediction_intervals(None, None, None, None, num_simulations=count)
