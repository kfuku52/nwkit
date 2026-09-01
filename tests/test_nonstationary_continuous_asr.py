import math

import numpy as np
import pandas as pd
import pytest

from nwkit.cli import main
from nwkit.continuous_asr import compute_bm_marginals
from nwkit.nonstationary_continuous_asr import (
    compute_bm_drift_marginals,
    compute_eb_marginals,
)
from nwkit.util import read_tree


def tree_from(source):
    return read_tree(source, "1", True, quiet=True, rooted="yes")


@pytest.mark.parametrize("model", ["EB", "BM-DRIFT"])
def test_zero_extension_parameter_is_exactly_brownian(model):
    tree = tree_from("((A:0.4,B:1.2)I:0.7,C:1.5,D:0.8)R;")
    values = {"A": -1.0, "B": 2.0, "C": 4.0, "D": None}
    expected, expected_fit = compute_bm_marginals(tree, values, sigma2=0.8)
    if model == "EB":
        posterior, fit = compute_eb_marginals(
            tree, values, sigma2=0.8, eb_rate=0.0
        )
    else:
        posterior, fit = compute_bm_drift_marginals(
            tree, values, sigma2=0.8, drift=0.0
        )
    assert fit.restricted_log_likelihood == pytest.approx(
        expected_fit.restricted_log_likelihood, abs=1e-12
    )
    for node in tree.traverse():
        assert posterior[node] == pytest.approx(expected[node], abs=1e-12)


def test_early_burst_fixed_rate_matches_analytic_star_variances():
    tree = tree_from("(A:1,B:2,C:1.5)R;")
    rate = math.log(2.0)
    sigma2 = 0.5
    posterior, _ = compute_eb_marginals(
        tree,
        {"A": 0.0, "B": 3.0, "C": None},
        sigma2=sigma2,
        eb_rate=rate,
    )
    variances = [sigma2 * math.expm1(rate * length) / rate for length in (1.0, 2.0)]
    precisions = [1 / value for value in variances]
    expected_mean = 3.0 * precisions[1] / sum(precisions)
    expected_variance = 1 / sum(precisions)
    assert posterior[tree].mean == pytest.approx(expected_mean)
    assert posterior[tree].variance == pytest.approx(expected_variance)


def test_early_burst_rate_is_rejected_when_star_profile_is_flat():
    tree = tree_from("(A:1,B:1,C:1,D:1)R;")
    with pytest.raises(ValueError, match="not identifiable"):
        compute_eb_marginals(tree, {"A": 0.0, "B": 1.0, "C": 2.0, "D": 3.0})


def test_drift_is_estimated_from_heterochronous_tips():
    tree = tree_from("(A:1,B:2,C:3,D:4)R;")
    values = {
        name: 1.0 + 2.0 * depth
        for name, depth in zip("ABCD", (1, 2, 3, 4), strict=True)
    }
    posterior, fit = compute_bm_drift_marginals(tree, values, sigma2=1.0)
    assert fit.drift == pytest.approx(2.0, abs=2e-7)
    assert fit.drift_estimated
    assert posterior[tree].mean == pytest.approx(1.0, abs=2e-7)


def test_exact_linear_trend_prefers_explicit_zero_diffusion_boundary():
    tree = tree_from("(A:1,B:2,C:3,D:4)R;")
    values = {
        name: 1.0 + 2.0 * depth
        for name, depth in zip("ABCD", (1, 2, 3, 4), strict=True)
    }
    posterior, fit = compute_bm_drift_marginals(tree, values)
    assert fit.drift == pytest.approx(2.0)
    assert fit.sigma2 == 0.0
    assert fit.fit_status == "singular_zero_boundary"
    assert fit.restricted_log_likelihood is None
    assert posterior[tree].mean == pytest.approx(1.0)


def test_free_drift_rejects_contemporaneous_observations():
    with pytest.raises(ValueError, match="one sampling depth"):
        compute_bm_drift_marginals(
            tree_from("(A:1,B:1,C:1)R;"),
            {"A": 0.0, "B": 1.0, "C": 2.0},
            sigma2=1.0,
        )


@pytest.mark.integration
@pytest.mark.parametrize(
    "model,parameter,value",
    [("EB", "--eb-rate", "-0.2"), ("BM-DRIFT", "--drift", "0.4")],
)
def test_nonstationary_cli_reports_selected_parameter(
    model, parameter, value, tmp_path
):
    trait = tmp_path / "traits.tsv"
    output = tmp_path / "asr.tsv"
    metadata = tmp_path / "model.tsv"
    trait.write_text("leaf_name\tvalue\nA\t0\nB\t1\nC\t3\n")
    main(
        [
            "asr",
            "-i",
            "[&R](A:1,B:2,C:3)R;",
            "--trait",
            str(trait),
            "--state-column",
            "value",
            "--model",
            model,
            "--sigma2",
            "1.2",
            parameter,
            value,
            "--model-out",
            str(metadata),
            "-o",
            str(output),
        ]
    )
    result = pd.read_csv(output, sep="\t")
    fit = pd.read_csv(metadata, sep="\t").iloc[0]
    column = "eb_rate" if model == "EB" else "drift"
    assert np.all(np.isfinite(result[["mean", "variance"]]))
    assert fit["model"] == model
    assert fit[column] == pytest.approx(float(value))
