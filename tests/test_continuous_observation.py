import numpy as np
import pandas as pd
import pytest

from nwkit.cli import main
from nwkit.continuous_observation import (
    read_measurement_covariances,
    validate_measurement_covariances,
)
from nwkit.util import read_tree


def covariance_file(path, names, matrix):
    rows = [
        dict(leaf_name=name, trait=first, other_trait=second, covariance=matrix[i, j])
        for name in names
        for i, first in enumerate(("x", "y"))
        for j, second in enumerate(("x", "y"))
    ]
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def test_covariance_tsv_validation(tmp_path):
    path = tmp_path / "cov.tsv"
    matrix = np.array([[0.1, 0.03], [0.03, 0.2]])
    covariance_file(path, ["A"], matrix)
    assert read_measurement_covariances(path, {"A": [1, None]}, ("x", "y"))[
        "A"
    ] == pytest.approx(matrix)
    with pytest.raises(ValueError, match="required"):
        read_measurement_covariances(path, {"A": [1, 2], "B": [2, 3]}, ("x", "y"))


@pytest.mark.parametrize(
    "matrix", [[[1, 2], [2, 1]], [[0, 1], [1, 1]], [[1, 0.1], [0.2, 1]]]
)
def test_invalid_covariance_rejected(matrix):
    with pytest.raises(ValueError, match="covariance|Covariance"):
        validate_measurement_covariances({"A": [0, 1]}, ("x", "y"), {"A": matrix})


def test_correlated_fit_matches_independent_dense_likelihood():
    from scipy.linalg import block_diag

    from nwkit.multivariate_asr import compute_mvbm_marginals

    tree = read_tree("(A:1,B:2,C:1,D:2)R;", "1", True, quiet=True, rooted="yes")
    observed = {"A": [0, 1], "B": [1, 0], "C": [3, 4], "D": [4, 2]}
    error = np.array([[0.1, 0.04], [0.04, 0.2]])
    posterior, fit = compute_mvbm_marginals(
        tree,
        observed,
        ("x", "y"),
        measurement_covariances={name: error for name in observed},
    )
    covariance = block_diag(*(fit.sigma * t + error for t in (1, 2, 1, 2)))
    design = np.tile(np.eye(2), (4, 1))
    y = np.array(list(observed.values())).ravel()
    inverse = np.linalg.inv(covariance)
    precision = design.T @ inverse @ design
    mean = np.linalg.solve(precision, design.T @ inverse @ y)
    residual = y - design @ mean
    loglike = -0.5 * (
        6 * np.log(2 * np.pi)
        + np.linalg.slogdet(covariance)[1]
        + np.linalg.slogdet(precision)[1]
        + residual @ inverse @ residual
    )
    assert fit.restricted_log_likelihood == pytest.approx(loglike, abs=1e-8)
    assert posterior[tree].mean == pytest.approx(mean, abs=1e-8)
    assert posterior[tree].covariance == pytest.approx(
        np.linalg.inv(precision), abs=1e-8
    )


@pytest.mark.integration
def test_correlated_covariance_cli_and_bootstrap(tmp_path):
    traits = tmp_path / "traits.tsv"
    traits.write_text("leaf_name\tx\ty\nA\t0\t1\nB\t1\t0\nC\t3\t4\nD\t4\t2\n")
    path = tmp_path / "cov.tsv"
    covariance_file(path, "ABCD", np.array([[0.1, 0.04], [0.04, 0.2]]))
    model = tmp_path / "model.tsv"
    samples = tmp_path / "samples.tsv"
    bootstrap = tmp_path / "bootstrap.tsv"
    base = [
        "-i",
        "(A:1,B:2,C:1,D:2)R;",
        "--input-rooted",
        "yes",
        "--trait",
        str(traits),
        "--state-column",
        "x,y",
        "--measurement-covariance",
        str(path),
    ]
    main(
        [
            "asr",
            *base,
            "--model",
            "MV-BM",
            "--model-out",
            str(model),
            "--posterior-samples-out",
            str(samples),
            "--posterior-samples",
            "2",
            "--bootstrap-out",
            str(bootstrap),
            "--bootstrap-simulations",
            "2",
            "--seed",
            "4",
            "-o",
            str(tmp_path / "out.tsv"),
        ]
    )
    assert len(pd.read_csv(samples, sep="\t")) == 20
    assert set(pd.read_csv(bootstrap, sep="\t").fit_status) == {"ok"}
    fitted = pd.read_csv(model, sep="\t")
    assert fitted.measurement_covariance.iloc[0] == str(path)
    comparison = tmp_path / "compare.tsv"
    main(["asrcompare", *base, "--models", "MV-BM", "-o", str(comparison)])
    assert pd.read_csv(comparison, sep="\t").log_likelihood.iloc[0] == pytest.approx(
        fitted.restricted_log_likelihood.iloc[0]
    )


def test_replicate_sufficient_statistic_retains_density():
    from scipy.stats import norm

    from nwkit.continuous_observation import summarize_independent_replicates

    values, errors = [1.0, 2.0, 4.0], [0.2, 0.5, 0.7]
    mean, error, constant = summarize_independent_replicates(values, errors)
    for latent in [-1.0, 0.0, 3.0]:
        expected = float(np.sum(norm.logpdf(values, loc=latent, scale=errors)))
        assert constant + norm.logpdf(mean, loc=latent, scale=error) == pytest.approx(
            expected, abs=1e-10
        )
    with pytest.raises(ValueError, match="positive"):
        summarize_independent_replicates([1, 2], [0, 1])


@pytest.mark.integration
def test_replicate_cli_likelihood_uses_raw_measurements(tmp_path):
    traits = tmp_path / "traits.tsv"
    traits.write_text("leaf_name\tx\nA\t999\nB\t3\n")
    replicates = tmp_path / "replicates.tsv"
    replicates.write_text(
        "leaf_name\ttrait\tvalue\tstandard_error\nA\tx\t0\t0.2\nA\tx\t2\t0.5\nB\tx\t3\t0.4\n"
    )
    model = tmp_path / "model.tsv"
    base = [
        "-i",
        "(A:1,B:2)R;",
        "--input-rooted",
        "yes",
        "--trait",
        str(traits),
        "--state-column",
        "x",
        "--replicate-observations",
        str(replicates),
        "--sigma2",
        "1",
    ]
    main(
        [
            "asr",
            *base,
            "--model",
            "BM",
            "--model-out",
            str(model),
            "-o",
            str(tmp_path / "out.tsv"),
        ]
    )
    covariance = np.array([[1.04, 1, 0], [1, 1.25, 0], [0, 0, 2.16]])
    inverse = np.linalg.inv(covariance)
    y = np.array([0, 2, 3])
    mean = np.sum(inverse @ y) / np.sum(inverse)
    residual = y - mean
    loglike = -0.5 * (
        2 * np.log(2 * np.pi)
        + np.linalg.slogdet(covariance)[1]
        + np.log(np.sum(inverse))
        + residual @ inverse @ residual
    )
    fitted = pd.read_csv(model, sep="\t")
    assert fitted.restricted_log_likelihood.iloc[0] == pytest.approx(loglike, abs=1e-8)
    assert fitted.replicate_observations.iloc[0] == str(replicates)
    comparison = tmp_path / "compare.tsv"
    main(["asrcompare", *base, "--models", "BM", "-o", str(comparison)])
    assert pd.read_csv(comparison, sep="\t").log_likelihood.iloc[0] == pytest.approx(
        loglike, abs=1e-8
    )
