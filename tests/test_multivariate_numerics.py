"""Independent numerical and resource regressions from the September audit."""

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from nwkit import multivariate_gaussian_asr as dense
from nwkit import multivariate_pgls as pgls
from nwkit.asr_compare import ComparisonCandidate, ComparisonContext, _fit_continuous
from nwkit.cli import main
from nwkit.gaussian import DiagonalLowRankCovariance
from nwkit.multivariate_asr import compute_mvbm_marginals
from nwkit.optimization import FitResourceError
from nwkit.sparse_laplace import SparseCovarianceModel
from nwkit.util import read_tree

VALUES = {
    "A": (0.0, 1.0),
    "B": (0.4, 0.8),
    "C": (1.2, 0.2),
    "D": (1.5, -0.3),
    "E": (1.7, 0.8),
    "F": (2.3, -0.6),
}


def tree(source="((A:.5,B:.7):.4,(C:.6,D:.8):.3,E:1,F:1.2)R;"):
    return read_tree(source, "1", True, quiet=True, rooted="yes")


def missing_design():
    design = np.column_stack([np.ones(8), [0, 0, 0, 0, 1, 1, 1, 1]])
    responses = np.array(
        [
            [0.0, 1.0],
            [0.4, 0.8],
            [1.2, 0.2],
            [1.5, -0.3],
            [np.nan, 0.8],
            [np.nan, -0.6],
            [np.nan, 0.3],
            [np.nan, 1.3],
        ]
    )
    return design, responses


@pytest.mark.parametrize("reml", [True, False])
@pytest.mark.parametrize("sparse_path", [True, False])
def test_missing_response_design_must_be_identifiable(reml, sparse_path, monkeypatch):
    if sparse_path:
        monkeypatch.setattr(pgls, "MAX_DENSE_MULTIVARIATE_DIMENSION", 1)
    design, responses = missing_design()
    assert np.linalg.matrix_rank(design) == 2
    with pytest.raises(ValueError, match="response 1.*rank-deficient"):
        pgls.fit_multivariate_pgls(responses, design, {"tree": np.eye(8)}, reml=reml)


def test_regress_cli_rejects_nonidentified_missing_response(tmp_path):
    source = tmp_path / "tree.nwk"
    source.write_text("((A:1,B:1,C:1,D:1):.2,(E:1,F:1,G:1,H:1):.2)R;")
    design, responses = missing_design()
    data = tmp_path / "data.tsv"
    pd.DataFrame(
        {
            "leaf_name": list("ABCDEFGH"),
            "x": design[:, 1],
            "y1": responses[:, 0],
            "y2": responses[:, 1],
        }
    ).to_csv(data, sep="\t", index=False)
    output = tmp_path / "result.tsv"
    with pytest.raises(ValueError, match="not identifiable"):
        main(
            [
                "regress",
                "--tree",
                str(source),
                "--tree-format",
                "1",
                "--data",
                str(data),
                "--responses",
                "y1,y2",
                "--predictors",
                "x",
                "--multivariate-responses",
                "yes",
                "--allow-missing-responses",
                "yes",
                "--reml",
                "no",
                "--outfile",
                str(output),
            ]
        )
    assert not output.exists()


def test_singular_covariance_is_not_given_artificial_noise():
    with pytest.raises((ValueError, np.linalg.LinAlgError)):
        pgls.fit_multivariate_pgls(
            np.array(list(VALUES.values())), np.ones((6, 1)), {"tree": np.zeros((6, 6))}
        )


def test_materialized_sparse_covariance_preserves_extreme_finite_entries():
    loading = np.array([[1e154, 0.0], [-1e154, 0.0], [0.0, 2.0**-537]])
    precision = sparse.eye(2, format="csc")
    projected = sparse.csr_matrix(loading)
    fitted = pgls.SparseFittedCovariance(precision, projected)
    latent = SparseCovarianceModel(
        precision=precision,
        tip_loading=projected,
        logdet_covariance=0.0,
        sampling_parent=np.full(2, -1),
        sampling_transition=np.zeros(2),
        sampling_variance=np.ones(2),
    )
    expected = loading @ loading.T
    assert expected[-1, -1] == np.nextafter(0.0, 1.0)
    with np.errstate(over="raise", invalid="raise"):
        np.testing.assert_array_equal(latent.materialize(), expected)
        np.testing.assert_array_equal(fitted.materialize(), expected)


def test_pgls_cholesky_factors_finite_covariance_near_maximum():
    correlation = np.array([[1.0, -0.5], [-0.5, 1.0]])
    with np.errstate(over="raise", invalid="raise"):
        factor = pgls._positive_cholesky(1e308 * correlation)
    np.testing.assert_allclose(factor / 1e154, np.linalg.cholesky(correlation))


@pytest.mark.parametrize("scale", [2.0**-1060, 2.0**1023])
def test_pgls_validates_covariance_in_its_own_units(scale):
    valid = scale * np.array([[1.0, 0.5], [0.5, 1.0]])
    asymmetric = scale * np.array([[1.0, 0.25], [0.5, 1.0]])
    indefinite = scale * np.array([[1.0, 1.5], [1.5, 1.0]])
    with np.errstate(over="raise", invalid="raise"):
        np.testing.assert_array_equal(
            pgls._validate_psd_covariance(valid, 2, "Test covariance"), valid
        )
        with pytest.raises(ValueError, match="must be symmetric"):
            pgls._validate_psd_covariance(asymmetric, 2, "Test covariance")
        with pytest.raises(ValueError, match="must be positive semidefinite"):
            pgls._validate_psd_covariance(indefinite, 2, "Test covariance")


def test_pgls_likelihood_matches_reported_covariance():
    responses = np.array(list(VALUES.values()))
    design = np.ones((6, 1))
    fit = pgls.fit_multivariate_pgls(responses, design, {"tree": np.eye(6)}, reml=False)
    residual = (responses - fit.coefficients[:, 0]).ravel(order="F")
    sign, logdet = np.linalg.slogdet(fit.fitted_covariance)
    assert sign == 1
    expected = -0.5 * (
        12 * np.log(2 * np.pi)
        + logdet
        + residual @ np.linalg.solve(fit.fitted_covariance, residual)
    )
    assert fit.log_likelihood == pytest.approx(expected, abs=1e-9)


@pytest.mark.parametrize("reml", [False, True])
@pytest.mark.parametrize("sparse_path", [False, True])
def test_pgls_rescales_traits_with_missing_values_and_sampling_error(
    reml, sparse_path, monkeypatch
):
    phylogeny = tree()
    names = list(VALUES)
    covariance = pgls.EvolutionaryCovarianceFactory(phylogeny, tuple(names))
    responses = np.asarray(list(VALUES.values()))
    responses[2, 0] = np.nan
    design = np.ones((len(names), 1))
    exponents = np.array([-260, 300])
    scales = np.ldexp(np.ones(2), exponents)
    cell_scales = np.repeat(scales, len(names))
    fixed_loading = np.linspace(-0.03, 0.06, responses.size)[:, None]
    fixed_diagonal = np.full(responses.size, 0.005)
    if sparse_path:
        monkeypatch.setattr(pgls, "MAX_DENSE_MULTIVARIATE_DIMENSION", 1)
    before = pgls.fit_multivariate_pgls(
        responses,
        design,
        {"tree": covariance},
        fixed_covariance=DiagonalLowRankCovariance(
            fixed_diagonal, sparse.csr_matrix(fixed_loading)
        ),
        reml=reml,
    )
    after = pgls.fit_multivariate_pgls(
        responses * scales,
        design,
        {"tree": covariance},
        fixed_covariance=DiagonalLowRankCovariance(
            fixed_diagonal * cell_scales**2,
            sparse.csr_matrix(fixed_loading * cell_scales[:, None]),
        ),
        reml=reml,
    )

    counts = np.isfinite(responses).sum(axis=0) - int(reml)
    assert after.log_likelihood + counts @ np.log(scales) == pytest.approx(
        before.log_likelihood, abs=1e-7
    )
    assert after.boundary_warning == before.boundary_warning
    np.testing.assert_allclose(
        after.coefficients[:, 0] / scales, before.coefficients[:, 0], atol=1e-6
    )
    np.testing.assert_allclose(
        after.coefficient_covariance / np.outer(scales, scales),
        before.coefficient_covariance,
        rtol=2e-5,
    )
    before_covariance, after_covariance = (
        before.fitted_covariance,
        after.fitted_covariance,
    )
    if sparse_path:
        before_covariance = before_covariance.materialize()
        after_covariance = after_covariance.materialize()
    np.testing.assert_allclose(
        after_covariance / np.outer(cell_scales, cell_scales),
        before_covariance,
        rtol=2e-5,
    )


@pytest.mark.parametrize("reml", [False, True])
@pytest.mark.parametrize("known_noise", [False, True])
def test_high_missing_multivariate_fit_retains_full_covariance_and_likelihood(
    reml, known_noise
):
    n, traits = 128, 4
    rng = np.random.default_rng(1142)
    x = rng.normal(size=n)
    design = np.column_stack([np.ones(n), x])
    covariance = 0.65 ** np.abs(np.arange(n)[:, None] - np.arange(n))
    trait_factor = np.eye(traits) + 0.15 * np.tril(np.ones((traits, traits)), -1)
    responses = (
        0.3
        + x[:, None] * np.linspace(0.5, 1.2, traits)
        + np.linalg.cholesky(covariance) @ rng.normal(size=(n, traits)) @ trait_factor.T
    )
    missing = rng.uniform(size=(n, traits)) < 0.9
    missing[np.arange(n), np.arange(n) % traits] = False
    responses[missing] = np.nan
    observed = np.flatnonzero(~missing.ravel(order="F"))
    assert 3 * len(observed) <= responses.size
    fixed = np.full(responses.size, 0.05) if known_noise else None
    fit = pgls.fit_multivariate_pgls(
        responses, design, {"tree": covariance}, fixed_covariance=fixed, reml=reml
    )
    assert fit.fitted_covariance.shape == (responses.size, responses.size)
    expected_covariance = np.kron(fit.component_trait_covariances["tree"], covariance)
    if known_noise:
        expected_covariance += np.diag(fixed)
    np.testing.assert_array_equal(fit.fitted_covariance, expected_covariance)
    selected = fit.fitted_covariance[np.ix_(observed, observed)]
    selected_design = np.kron(np.eye(traits), design)[observed]
    residual = (
        responses.ravel(order="F")[observed]
        - selected_design @ fit.coefficients.ravel()
    )
    _, logdet = np.linalg.slogdet(selected)
    quadratic = residual @ np.linalg.solve(selected, residual)
    count = len(observed)
    if reml:
        count -= selected_design.shape[1]
        _, information = np.linalg.slogdet(
            selected_design.T @ np.linalg.solve(selected, selected_design)
        )
        logdet += information
    expected = -0.5 * (count * np.log(2 * np.pi) + logdet + quadratic)
    assert fit.log_likelihood == pytest.approx(expected, abs=1e-8)


def test_pgls_cached_components_preserve_covariance_scale(monkeypatch):
    responses = np.array(list(VALUES.values()))
    design = np.ones((6, 1))
    model = SparseCovarianceModel(
        precision=sparse.eye(6, format="csc"),
        tip_loading=sparse.eye(6, format="csr"),
        logdet_covariance=0.0,
        sampling_parent=np.full(6, -1),
        sampling_transition=np.zeros(6),
        sampling_variance=np.ones(6),
        covariance_scale=3.0,
    )
    first = pgls.fit_multivariate_pgls(responses, design, {"tree": model})
    monkeypatch.setattr(pgls, "MAX_DENSE_MULTIVARIATE_DIMENSION", 1)
    second = pgls.fit_multivariate_pgls(responses, design, {"tree": model})
    assert first.log_likelihood == pytest.approx(second.log_likelihood, abs=1e-7)
    np.testing.assert_allclose(
        first.component_trait_covariances["tree"],
        second.component_trait_covariances["tree"],
        rtol=1e-4,
    )


@pytest.mark.parametrize(
    "fitter", [dense.fit_dense_mvbm, dense.fit_dense_mvou, dense.fit_dense_mvou_diag]
)
@pytest.mark.parametrize("scale", [1e-16, 1e16])
def test_multivariate_time_units_preserve_fits_and_posteriors(fitter, scale):
    original = tree()
    transformed = tree()
    for node in transformed.traverse():
        if not node.is_root:
            node.dist *= scale
    is_bm = fitter is dense.fit_dense_mvbm
    errors = {name: (0.05, 0.08) for name in VALUES}
    kwargs = {} if is_bm else {"alpha": 0.7}
    first_post, first = fitter(
        original, VALUES, ("x", "y"), standard_errors=errors, **kwargs
    )
    kwargs = {} if is_bm else {"alpha": 0.7 / scale}
    second_post, second = fitter(
        transformed, VALUES, ("x", "y"), standard_errors=errors, **kwargs
    )
    assert second.sigma_rank == first.sigma_rank == 2
    np.testing.assert_allclose(
        second.sigma * (scale if is_bm else 1.0), first.sigma, rtol=2e-4, atol=1e-6
    )
    first_likelihood = (
        first.restricted_log_likelihood if is_bm else first.log_likelihood
    )
    second_likelihood = (
        second.restricted_log_likelihood if is_bm else second.log_likelihood
    )
    assert second_likelihood == pytest.approx(first_likelihood, abs=1e-6)
    np.testing.assert_allclose(
        first_post[original].mean, second_post[transformed].mean, rtol=2e-4, atol=1e-6
    )
    np.testing.assert_allclose(
        first_post[original].covariance,
        second_post[transformed].covariance,
        rtol=2e-4,
        atol=1e-6,
    )


@pytest.mark.parametrize("scale", [1e-200, 1e200])
def test_free_mvou_attraction_fit_is_invariant_to_extreme_time_units(scale):
    original = tree()
    scaled = tree()
    for node in scaled.traverse():
        if not node.is_root:
            node.dist *= scale
    errors = {name: (0.05, 0.08) for name in VALUES}
    _, first = dense.fit_dense_mvou(
        original, VALUES, ("x", "y"), standard_errors=errors
    )
    _, second = dense.fit_dense_mvou(scaled, VALUES, ("x", "y"), standard_errors=errors)
    assert second.alpha * scale == pytest.approx(first.alpha, rel=2e-5)
    assert second.log_likelihood == pytest.approx(first.log_likelihood, abs=1e-9)
    assert second.fit_status == first.fit_status


@pytest.mark.parametrize("fitter", [dense.fit_dense_mvou, dense.fit_dense_mvou_diag])
def test_fixed_mvou_attraction_accepts_subnormal_tree_time_units(fitter):
    original = tree()
    scaled = tree()
    scale = 1e-308
    for node in scaled.traverse():
        if not node.is_root:
            node.dist *= scale
    expected, first = fitter(original, VALUES, ("x", "y"), alpha=0.3)
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        posterior, second = fitter(scaled, VALUES, ("x", "y"), alpha=0.3 / scale)
    assert second.log_likelihood == pytest.approx(first.log_likelihood, abs=1e-9)
    np.testing.assert_allclose(second.sigma, first.sigma, rtol=1e-5)
    np.testing.assert_allclose(
        posterior[scaled].mean, expected[original].mean, atol=1e-6
    )


def test_diagonal_ou_with_extreme_fixed_attraction_has_finite_diffusion():
    original = tree()
    values = np.asarray(list(VALUES.values()))
    estimated_mean = values.mean(axis=0)
    centered = values - estimated_mean
    expected_sigma = centered.T @ centered / len(values)

    with np.errstate(over="raise", invalid="raise", divide="raise"):
        _, fit = dense.fit_dense_mvou_diag(original, VALUES, ("x", "y"), alpha=1e308)

    assert fit.optimizer_success
    assert fit.alpha == 1e308
    np.testing.assert_allclose(fit.theta, estimated_mean, atol=1e-7)
    np.testing.assert_allclose(fit.sigma, expected_sigma, rtol=1e-5)
    np.testing.assert_allclose(
        fit.diffusion_sigma / 1e308, 2 * expected_sigma, rtol=1e-5
    )


@pytest.mark.parametrize("scale", [1e-180, 1e16, 1e180])
def test_fast_bm_rank_and_likelihood_are_time_unit_invariant(scale):
    current = tree()
    expected, baseline = compute_mvbm_marginals(current, VALUES, ("x", "y"))
    for node in current.traverse():
        if not node.is_root:
            node.dist *= scale
    with np.errstate(over="raise", invalid="raise"):
        posterior, fit = compute_mvbm_marginals(current, VALUES, ("x", "y"))
    assert fit.sigma_rank == baseline.sigma_rank == 2
    assert fit.restricted_log_likelihood == pytest.approx(
        baseline.restricted_log_likelihood, abs=1e-9
    )
    np.testing.assert_allclose(fit.sigma * scale, baseline.sigma)
    for node in current.traverse():
        np.testing.assert_allclose(
            posterior[node].mean, expected[node].mean, atol=1e-12
        )
        np.testing.assert_allclose(
            posterior[node].covariance, expected[node].covariance
        )


@pytest.mark.parametrize("model", ["MV-BM", "MV-OU", "MV-OU-DIAG"])
def test_identifiability_checks_are_time_unit_invariant(model):
    for scale in (1.0, 1e-16):
        current = tree()
        for node in current.traverse():
            if not node.is_root:
                node.dist *= scale
        data = dense._prepare_observations(current, VALUES, ("x", "y"))
        geometry = dense._geometry(data)
        if model == "MV-BM":
            dense._validate_mvbm_covariance_design(
                data, geometry.observed_shared_depth, np.full(2, np.nan)
            )
        elif model == "MV-OU":
            dense._validate_mvou_alpha_design(data, geometry.observed_distance)
        else:
            dense._validate_mvou_diag_alpha_design(data, geometry)


def test_fast_bm_likelihood_only_matches_full_fit(monkeypatch):
    from nwkit import multivariate_asr

    current = tree()
    _, expected = compute_mvbm_marginals(current, VALUES, ("x", "y"))

    def forbidden(*args, **kwargs):
        raise AssertionError("likelihood-only fit attempted posterior")

    monkeypatch.setattr(multivariate_asr, "_mvbm_posterior", forbidden)
    posterior, actual = compute_mvbm_marginals(
        current, VALUES, ("x", "y"), compute_posterior=False
    )
    assert posterior == {}
    assert actual.restricted_log_likelihood == expected.restricted_log_likelihood
    np.testing.assert_array_equal(actual.sigma, expected.sigma)


@pytest.mark.parametrize("fitter", [dense.fit_dense_mvou, dense.fit_dense_mvou_diag])
def test_stationary_ou_preserves_distance_below_long_common_stem(fitter):
    likelihoods = []
    for stem in (0.0, 1e16):
        current = tree(f"(((A:1,B:1):1,(C:1,D:1):1,E:2,F:2):{stem})R;")
        data = dense._prepare_observations(current, VALUES, ("x", "y"))
        geometry = dense._geometry(data)
        assert geometry.observed_distance[0, 2] == pytest.approx(2.0)
        _, fit = fitter(current, VALUES, ("x", "y"), alpha=0.7)
        likelihoods.append(fit.log_likelihood)
    assert likelihoods[0] == pytest.approx(likelihoods[1], abs=1e-6)


def test_geometry_does_not_allocate_nodes_by_observations():
    current = tree("(" + ",".join(f"t{i}:1" for i in range(2000)) + ")R;")
    data = dense._prepare_observations(
        current, {f"t{i}": (float(i), float(i % 3)) for i in range(10)}, ("x", "y")
    )
    geometry = dense._geometry(data)
    arrays = [
        value for value in vars(geometry).values() if isinstance(value, np.ndarray)
    ]
    assert all(value.shape != (2001, 20) for value in arrays)
    shared, left, right = geometry.cross(0)
    assert shared.shape == left.shape == right.shape == (20,)
    np.testing.assert_array_equal(shared, np.zeros(20))
    np.testing.assert_array_equal(right, np.ones(20))


def test_comparison_skips_posterior_and_shares_geometry(monkeypatch):
    current = tree()
    dataframe = pd.DataFrame(
        {
            "leaf_name": list(VALUES),
            "x": [v[0] for v in VALUES.values()],
            "y": [v[1] for v in VALUES.values()],
        }
    )
    context = ComparisonContext(
        current, dataframe, "continuous", ("x", "y"), None, SimpleNamespace(alpha=0.7)
    )

    def forbidden(*args, **kwargs):
        raise AssertionError("model comparison attempted ancestral reconstruction")

    monkeypatch.setattr(dense, "_posterior", forbidden)
    monkeypatch.setattr(dense, "_diagonal_ou_posterior", forbidden)
    first = _fit_continuous(
        context, ComparisonCandidate("MV-OU", "stationary", "default")
    )
    monkeypatch.setattr(dense, "_geometry", forbidden)
    second = _fit_continuous(
        context, ComparisonCandidate("MV-OU-DIAG", "stationary", "default")
    )
    assert first.log_likelihood == pytest.approx(second.log_likelihood, abs=1e-7)


def test_posterior_storage_limit_does_not_block_likelihood_only(monkeypatch):
    monkeypatch.setattr(dense, "_MAX_POSTERIOR_BYTES", 1)
    with pytest.raises(FitResourceError, match="result-storage"):
        dense.fit_dense_mvou(tree(), VALUES, ("x", "y"), alpha=0.7)
    posterior, fit = dense.fit_dense_mvou(
        tree(), VALUES, ("x", "y"), alpha=0.7, compute_posterior=False
    )
    assert posterior == {}
    assert np.isfinite(fit.log_likelihood)
