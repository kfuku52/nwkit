import numpy as np
import pytest
from scipy.linalg import expm, solve_continuous_lyapunov
from scipy.stats import multivariate_normal

from nwkit.full_ou import (
    decode_full_ou,
    full_ou_identifiability,
    full_ou_initial,
    observation_moment_design,
)
from nwkit.full_ou_fit import fit_full_mvou
from nwkit.multivariate_gaussian_asr import _prepare_observations
from nwkit.util import read_tree


def _tree(newick="(A:1,B:2,C:3,D:4)R;"):
    return read_tree(newick, "1", True, quiet=True, rooted="yes")


def test_stable_parameterization_allows_negative_diagonal_and_complex_eigenvalues():
    attraction = np.array([[-1.0, 10.0], [-1.0, 3.0]])
    diffusion = np.eye(2)
    covariance = solve_continuous_lyapunov(attraction, diffusion)
    parameters = []
    for matrix in (covariance, diffusion):
        lower = np.linalg.cholesky(matrix)
        parameters += [np.log(lower[0, 0]), lower[1, 0], np.log(lower[1, 1])]
    parameters += [(attraction @ covariance - diffusion / 2)[1, 0]]
    actual, noise, stationary = decode_full_ou(parameters, 2)
    assert actual == pytest.approx(attraction)
    assert noise == pytest.approx(diffusion)
    assert stationary == pytest.approx(covariance)
    assert actual @ stationary + stationary @ actual.T == pytest.approx(noise)


@pytest.mark.parametrize("scale", [np.ones(2), np.array([1e-7, 1e7])])
def test_fixed_full_ou_matches_independent_dense_gls(scale):
    tree = _tree()
    attraction = np.array([[0.7, -0.3], [0.2, 0.5]])
    diffusion = np.array([[1.2, 0.3], [0.3, 0.8]])
    values = np.array([[0.2, 1.0], [1.1, -0.4], [2.2, 0.3], [-1.0, 1.5]])
    observed = {
        node.name: list(values[i] * scale) for i, node in enumerate(tree.leaves())
    }
    observed["B"][1] = None
    errors = {name: np.array([0.1, 0.2]) * scale for name in observed}
    posterior, fit = fit_full_mvou(
        tree,
        observed,
        ("x", "y"),
        attraction=attraction * scale[:, None] / scale[None, :],
        diffusion=diffusion * np.outer(scale, scale),
        standard_errors=errors,
    )
    stationary = solve_continuous_lyapunov(attraction, diffusion)
    covariance = np.zeros((8, 8))
    for i in range(4):
        for j in range(4):
            block = (
                stationary
                if i == j
                else expm(-attraction * (i + 1))
                @ stationary
                @ expm(-attraction.T * (j + 1))
            )
            covariance[2 * i : 2 * i + 2, 2 * j : 2 * j + 2] = block
    covariance += np.diag([0.01, 0.04] * 4)
    keep = np.array([0, 1, 2, 4, 5, 6, 7])
    covariance = covariance[np.ix_(keep, keep)]
    design = np.tile(np.eye(2), (4, 1))[keep]
    values = values.ravel()[keep]
    precision_design = np.linalg.solve(covariance, design)
    theta = np.linalg.solve(design.T @ precision_design, precision_design.T @ values)
    assert fit.theta / scale == pytest.approx(theta, abs=2e-6)
    assert fit.log_likelihood == pytest.approx(
        multivariate_normal.logpdf(values, mean=design @ theta, cov=covariance)
        - np.log(scale)[keep % 2].sum(),
        abs=1e-8,
    )
    assert fit.identifiability_status == "fixed_covariance_parameters"
    assert posterior[tree].covariance.shape == (2, 2)


@pytest.mark.slow
def test_free_full_ou_real_optimizer_returns_stable_finite_process():
    tree = _tree("(((A:1,B:2):1,(C:1,D:3):1):1,((E:2,F:1):2,(G:1,H:4):1):1)R;")
    rng = np.random.default_rng(33)
    observed = {name: rng.normal(size=2) for name in tree.leaf_names()}
    errors = {name: [0.2, 0.2] for name in observed}
    posterior, fit = fit_full_mvou(tree, observed, ("x", "y"), standard_errors=errors)
    assert fit.optimizer_success
    assert fit.optimizer_converged_starts > 0
    assert np.isfinite(fit.log_likelihood)
    assert np.linalg.eigvals(fit.attraction_matrix).real.min() > 0
    assert np.linalg.eigvalsh(fit.diffusion_sigma).min() > 0
    assert (
        fit.attraction_matrix @ fit.sigma + fit.sigma @ fit.attraction_matrix.T
        == pytest.approx(fit.diffusion_sigma, rel=1e-7, abs=1e-8)
    )
    assert all(np.isfinite(marginal.mean).all() for marginal in posterior.values())
    # The unrestricted optimum cannot be worse than a feasible fixed process
    # obtained from its own fitted covariance parameters.
    _, fixed = fit_full_mvou(
        tree,
        observed,
        ("x", "y"),
        standard_errors=errors,
        attraction=fit.attraction_matrix,
        diffusion=fit.diffusion_sigma,
    )
    assert fit.log_likelihood == pytest.approx(fixed.log_likelihood, abs=1e-4)


def test_local_identifiability_detects_ultrametric_rotation():
    initial, _ = full_ou_initial(2)
    values = {name: [float(i), float(i % 2)] for i, name in enumerate("ABCD")}
    data = _prepare_observations(_tree("(A:1,B:1,C:1,D:1)R;"), values, ("x", "y"), None)
    design, complete = observation_moment_design(data)
    status, rank, _ = full_ou_identifiability(initial, 2, 1, design, complete)
    assert status == "local_rank_deficient"
    assert rank < len(initial)
    data = _prepare_observations(_tree(), values, ("x", "y"), None)
    design, complete = observation_moment_design(data)
    assert full_ou_identifiability(initial, 2, 1, design, complete)[:2] == (
        "local_full_rank",
        len(initial),
    )
    subset, complete = observation_moment_design(data, limit=2)
    assert not complete
    assert (
        full_ou_identifiability(initial, 2, 1, subset, complete)[0]
        == "inconclusive_design_subset"
    )


def test_full_ou_rejects_partial_fixed_matrices_and_insufficient_data():
    values = {"A": [1, 2], "B": [3, 4], "C": [2, 1], "D": [0, 3]}
    with pytest.raises(ValueError, match="both fixed"):
        fit_full_mvou(_tree(), values, ("x", "y"), attraction=np.eye(2))
    with pytest.raises(ValueError, match="too few"):
        fit_full_mvou(_tree(), values, ("x", "y"))


@pytest.mark.integration
def test_full_ou_cli_comparison_and_diagnostics(tmp_path):
    import pandas as pd

    from nwkit.cli import main

    traits = tmp_path / "traits.tsv"
    traits.write_text("leaf_name\tx\ty\nA\t0\t1\nB\t1\t0\nC\t3\t4\nD\t4\t2\n")
    base = [
        "-i",
        "(A:1,B:2,C:3,D:4)R;",
        "--input-rooted",
        "yes",
        "--trait",
        str(traits),
        "--state-column",
        "x,y",
        "--attraction-matrix",
        "0.7,-0.3;0.2,0.5",
        "--diffusion-matrix",
        "1.2,0.3;0.3,0.8",
    ]
    model = tmp_path / "model.tsv"
    samples = tmp_path / "samples.tsv"
    bootstrap = tmp_path / "bootstrap.tsv"
    predictive = tmp_path / "predictive.tsv"
    main(
        [
            "asr",
            *base,
            "--model",
            "MV-OU-FULL",
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
            "--posterior-predictive-out",
            str(predictive),
            "--posterior-predictive-simulations",
            "2",
            "--seed",
            "5",
            "-o",
            str(tmp_path / "out.tsv"),
        ]
    )
    fitted = pd.read_csv(model, sep="\t").iloc[0]
    assert fitted.identifiability_status == "fixed_covariance_parameters"
    assert fitted.attraction_78_to_79 == pytest.approx(-0.3)
    assert len(pd.read_csv(samples, sep="\t")) == 20
    assert set(pd.read_csv(bootstrap, sep="\t").fit_status) == {"ok"}
    comparison = tmp_path / "compare.tsv"
    main(["asrcompare", *base, "--models", "MV-OU-FULL", "-o", str(comparison)])
    row = pd.read_csv(comparison, sep="\t").iloc[0]
    assert row.num_parameters == 2
    assert row.log_likelihood == pytest.approx(fitted.log_likelihood, abs=1e-8)
    assert "attraction_matrix" in row.fixed_parameters


def test_full_ou_nonidentifiable_fit_excluded_from_ranking():
    from types import SimpleNamespace

    from nwkit.asr_compare import ComparisonCandidate, _classify_fit
    from nwkit.asr_comparison import _continuous_parameter_count

    fit = SimpleNamespace(
        model="MV-OU-FULL",
        trait_names=("x", "y"),
        theta_estimated=True,
        attraction_estimated=True,
    )
    assert _continuous_parameter_count(fit) == 9
    summary = {"fit_status": "local_rank_deficient", "log_likelihood": -20}
    assert _classify_fit(
        ComparisonCandidate("MV-OU-FULL", "stationary", "default"), fit, summary
    )[:2] == ("nonregular", "no")


def test_reported_optimizer_success_cannot_leave_mean_unoptimized(monkeypatch):
    from nwkit.optimization import MultistartResult

    tree = _tree("(((A:1,B:2):1,(C:1,D:3):1):1,((E:2,F:1):2,(G:1,H:4):1):1)R;")
    rng = np.random.default_rng(33)
    observed = {name: rng.normal(size=2) for name in tree.leaf_names()}
    errors = {name: [0.2, 0.2] for name in observed}

    def premature(objective, initial, bounds, **kwargs):
        parameters = np.asarray(initial, dtype=float)
        parameters[-2:] = [4.0, -3.0]
        return MultistartResult(
            parameters, objective(parameters), True, "premature", 1, 1, 0
        )

    monkeypatch.setattr("nwkit.full_ou_fit.deterministic_multistart", premature)
    posterior, fit = fit_full_mvou(tree, observed, ("x", "y"), standard_errors=errors)
    fixed_posterior, fixed = fit_full_mvou(
        tree,
        observed,
        ("x", "y"),
        standard_errors=errors,
        attraction=fit.attraction_matrix,
        diffusion=fit.diffusion_sigma,
    )
    assert fit.theta == pytest.approx(fixed.theta, abs=1e-10)
    assert fit.log_likelihood == pytest.approx(fixed.log_likelihood, abs=1e-10)
    for node in tree.traverse():
        assert posterior[node].mean == pytest.approx(fixed_posterior[node].mean)
