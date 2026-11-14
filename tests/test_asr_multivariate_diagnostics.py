import numpy as np
import pandas as pd
import pytest

from nwkit.asr_multivariate_diagnostics import fitted_vector_process
from nwkit.cli import main
from nwkit.multivariate_asr import compute_mvbm_marginals
from nwkit.util import read_tree
from nwkit.vector_gaussian import condition_vector_tree


def test_fitted_vector_process_matches_mvbm_marginals():
    tree = read_tree("((A:1,B:1):1,(C:1,D:1):1)R;", "1", True, quiet=True, rooted="yes")
    observed = {"A": [0.0, 1.0], "B": [1.0, 0.0], "C": [3.0, 4.0], "D": [4.0, 2.0]}
    expected, fit = compute_mvbm_marginals(tree, observed, ("x", "y"))
    result = condition_vector_tree(fitted_vector_process(tree, "MV-BM", fit), observed)
    for i, node in enumerate(result.nodes):
        assert result.means[i] == pytest.approx(expected[node].mean)
        assert result.covariances[i] == pytest.approx(expected[node].covariance)


@pytest.mark.integration
def test_multivariate_sample_cli_reproducible(tmp_path):
    traits = tmp_path / "traits.tsv"
    traits.write_text("leaf_name\tx\ty\nA\t0\t1\nB\t1\t0\nC\t3\t4\nD\t4\t2\n")
    outputs = []
    for index in range(2):
        out = tmp_path / f"samples{index}.tsv"
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
                "x,y",
                "--model",
                "MV-BM",
                "--posterior-samples-out",
                str(out),
                "--posterior-samples",
                "10",
                "--seed",
                "2",
                "--posterior-predictive-out",
                str(tmp_path / f"ppc{index}.tsv"),
                "--posterior-predictive-simulations",
                "10",
                "--bootstrap-out",
                str(tmp_path / f"bootstrap{index}.tsv"),
                "--bootstrap-simulations",
                "3",
                "-o",
                str(tmp_path / "asr.tsv"),
            ]
        )
        outputs.append(pd.read_csv(out, sep="\t"))
    pd.testing.assert_frame_equal(*outputs)
    assert len(outputs[0]) == 10 * 7 * 2
    assert np.isfinite(outputs[0].value).all()
    for prefix in ("ppc", "bootstrap"):
        first, second = [
            pd.read_csv(tmp_path / f"{prefix}{i}.tsv", sep="\t") for i in range(2)
        ]
        pd.testing.assert_frame_equal(first, second)
    assert set(first.fit_status) == {"ok"}
    assert "sigma_0_1" in first.columns
    predictive = pd.read_csv(tmp_path / "ppc0.tsv", sep="\t")
    assert len(predictive) == 11
    assert "covariance" in set(predictive.statistic)


@pytest.mark.parametrize("model", ["MV-OU", "MV-OU-DIAG"])
def test_fitted_ou_vector_matches_dense_marginals(model):
    from nwkit.multivariate_gaussian_asr import fit_dense_mvou, fit_dense_mvou_diag

    tree = read_tree("((A:1,B:2):1,(C:1,D:1):2)R;", "1", True, quiet=True, rooted="yes")
    observed = {"A": [0.0, 1.0], "B": [1.0, 0.0], "C": [3.0, 4.0], "D": [4.0, 2.0]}
    fitter = fit_dense_mvou if model == "MV-OU" else fit_dense_mvou_diag
    expected, fit = fitter(tree, observed, ("x", "y"), alpha=0.4)
    result = condition_vector_tree(fitted_vector_process(tree, model, fit), observed)
    assert result.log_likelihood == pytest.approx(fit.log_likelihood, abs=1e-7)
    for i, node in enumerate(result.nodes):
        assert result.means[i] == pytest.approx(expected[node].mean, abs=1e-8)
        assert result.covariances[i] == pytest.approx(
            expected[node].covariance, abs=1e-8
        )
