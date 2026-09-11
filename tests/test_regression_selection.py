import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.optimize import approx_fprime
from scipy.special import expit

from nwkit.cli import main
from nwkit.penalized_phylogenetic import (
    LaplaceLoss,
    fit_elastic_net,
    prediction_loss,
    validate_inputs,
)
from nwkit.regression_selection import _read_table, nested_selection


@pytest.mark.parametrize("unit", [1e-200, 1e-18, 1e18, 1e200])
def test_selection_is_invariant_to_predictor_units(unit):
    x, y, covariance, _ = example()
    original = fit_elastic_net(x, y, covariance)
    scaled = fit_elastic_net(x * unit, y, covariance)
    np.testing.assert_array_equal(scaled.active, original.active)
    np.testing.assert_allclose(scaled.predict(x * unit), original.predict(x), atol=1e-5)
    np.testing.assert_allclose(
        scaled.standardized_coefficients, original.standardized_coefficients, atol=1e-5
    )


def test_prediction_preserves_accuracy_with_large_predictor_origin():
    x = np.arange(12.0)[:, None]
    y = 3 * x[:, 0] + np.sin(x[:, 0])
    original = fit_elastic_net(x, y, np.eye(12))
    shifted = fit_elastic_net(x + 1e15, y, np.eye(12))
    np.testing.assert_allclose(
        shifted.predict(x + 1e15), original.predict(x), atol=1e-5
    )


def test_nonfinite_squared_prediction_loss_is_rejected():
    with pytest.raises(FloatingPointError):
        prediction_loss(np.array([1e200]), np.array([-1e200]), "gaussian")


def test_tiny_covariance_cannot_hide_material_asymmetry():
    covariance = np.eye(6) * 1e-12
    covariance[0, 1] = 5e-13
    with pytest.raises(ValueError, match="symmetric"):
        validate_inputs(np.arange(6.0)[:, None], np.arange(6.0), covariance, "gaussian")


def test_binomial_variance_boundary_requires_stationarity_after_relative_stop():
    rng = np.random.default_rng(52)
    x = rng.poisson(2, (3, 18)).T[6:12]
    y = np.tile([0, 1], 3)
    fit = fit_elastic_net(x, y, np.eye(6), family="binomial", strength=0.5)
    assert fit.projected_gradient < 2e-4
    assert fit.objective < 0.6802


@pytest.mark.parametrize(
    "text", ["leaf_name\tx\na\t1\t2\n", "leaf_name\tx\na\n", "leaf_name\t\na\t1\n"]
)
def test_malformed_table_cannot_silently_shift_species_or_columns(tmp_path, text):
    path = tmp_path / "bad.tsv"
    path.write_text(text)
    with pytest.raises(ValueError, match="fields|non-empty"):
        _read_table(path)


def example():
    rng = np.random.default_rng(91)
    x = rng.normal(size=(18, 4))
    groups = np.repeat(["a", "b", "c"], 6)
    covariance = np.eye(18) + 0.6 * (groups[:, None] == groups)
    y = 0.5 + 2 * x[:, 0] + rng.normal(scale=0.2, size=18)
    return x, y, covariance, groups


@pytest.mark.parametrize("family", ["binomial", "poisson", "negative-binomial"])
def test_laplace_fixed_linear_gradient_matches_numeric(family):
    y = (
        np.array([0, 1, 0, 1, 1, 0.0])
        if family == "binomial"
        else np.array([0, 1, 2, 1, 5, 3.0])
    )
    covariance = np.eye(6) + 0.3
    linear = np.linspace(-0.5, 0.8, 6)
    nuisance = (
        np.array([-0.5, -0.2]) if family == "negative-binomial" else np.array([-0.5])
    )
    loss = LaplaceLoss(y, covariance, family)
    value, gradient, _ = loss.value_gradient(linear, nuisance)
    numeric = approx_fprime(linear, lambda eta: loss.state(eta, nuisance)[0], 1e-5)
    assert np.isfinite(value)
    np.testing.assert_allclose(gradient, numeric, atol=3e-5, rtol=2e-3)


@pytest.mark.parametrize("ratio", [1.0, 0.5])
def test_gaussian_elastic_net_matches_orthogonal_soft_threshold(ratio):
    rng = np.random.default_rng(2)
    q, _ = np.linalg.qr(np.column_stack([np.ones(20), rng.normal(size=(20, 3))]))
    x = q[:, 1:] * np.sqrt(20)
    y = 0.7 + x @ np.array([2.0, -0.5, 0.1])
    fit = fit_elastic_net(x, y, np.eye(20), strength=0.2, l1_ratio=ratio)
    beta = np.array([2.0, -0.5, 0.1])
    expected = (
        np.sign(beta)
        * np.maximum(np.abs(beta) - 0.2 * ratio, 0)
        / (1 + 0.2 * (1 - ratio))
    )
    np.testing.assert_allclose(fit.coefficients, np.r_[0.7, expected], atol=2e-5)


def test_gaussian_coordinate_refinement_does_not_trust_early_optimizer_stop(
    monkeypatch,
):
    from types import SimpleNamespace

    import nwkit.penalized_phylogenetic as penalized

    def stop_at_start(objective, initial, **kwargs):
        return SimpleNamespace(x=initial, fun=objective(initial)[0], nit=0)

    monkeypatch.setattr(penalized, "minimize", stop_at_start)
    x = np.array([-1, 1, -1, 1, -1, 1, -1, 1], dtype=float)[:, None]
    y = 0.3 + 2 * x[:, 0]
    fit = fit_elastic_net(x, y, np.eye(8), strength=0.2, l1_ratio=0.5)
    np.testing.assert_allclose(fit.coefficients, [0.3, 1.9 / 1.1], atol=1e-8)
    assert fit.projected_gradient < 1e-8


def test_gaussian_matches_whitened_gls_at_small_penalty():
    x, y, covariance, _ = example()
    design = np.column_stack([np.ones(len(y)), x])
    expected = np.linalg.solve(
        design.T @ np.linalg.solve(covariance, design),
        design.T @ np.linalg.solve(covariance, y),
    )
    fitted = fit_elastic_net(x, y, covariance, strength=1e-8, l1_ratio=0.5)
    np.testing.assert_allclose(fitted.coefficients, expected, atol=2e-5)


def test_p_greater_n_correlated_predictors_and_free_covariate():
    x, y, covariance, _ = example()
    design = np.column_stack([x, np.tile(x[:, :1], (1, 20)), np.ones(len(y))])
    fit = fit_elastic_net(
        design, y, covariance, strength=0.1, l1_ratio=0.5, unpenalized=(1,)
    )
    assert fit.coefficients.shape == (26,)
    assert not fit.active[-1]
    assert fit.coefficients[-1] == 0
    assert fit.projected_gradient < 2e-4
    assert fit.standardized_coefficients[1] == pytest.approx(
        fit.standardized_coefficients[5], abs=1e-5
    )


@pytest.mark.parametrize("family", ["binomial", "poisson", "negative-binomial"])
def test_non_gaussian_fit(family):
    x, _, covariance, _ = example()
    rng = np.random.default_rng(4)
    y = (
        rng.binomial(1, expit(x[:, 0]))
        if family == "binomial"
        else rng.poisson(np.exp(0.5 + 0.4 * x[:, 0]))
    )
    fit = fit_elastic_net(x, y, covariance, family=family, strength=0.1)
    assert fit.projected_gradient < 2e-4
    assert np.isfinite(fit.predict(x)).all()
    assert (
        fit.dispersion is not None
        if family == "negative-binomial"
        else fit.dispersion is None
    )


def test_nested_cv_learns_preprocessing_in_training_fold_only():
    x, y, covariance, groups = example()
    kwargs = dict(strengths=(0.2, 0.05), l1_ratios=(0.5,))
    _, _, fits, _, prediction = nested_selection(x, y, covariance, groups, **kwargs)
    shifted = y.copy()
    shifted[groups == "a"] += 100
    _, _, altered_fits, _, altered_prediction = nested_selection(
        x, shifted, covariance, groups, **kwargs
    )
    np.testing.assert_allclose(fits[0][1].center, x[groups != "a"].mean(axis=0))
    np.testing.assert_allclose(fits[0][1].coefficients, altered_fits[0][1].coefficients)
    np.testing.assert_allclose(
        prediction.loc[prediction.fold == "a", "predicted"],
        altered_prediction.loc[altered_prediction.fold == "a", "predicted"],
    )


def test_fail_closed_for_unusable_inner_binary_fold():
    x, _, covariance, groups = example()
    y = (groups == "a").astype(float)
    with pytest.raises(ValueError, match="variable response"):
        nested_selection(x, y, covariance, groups, family="binomial")
    with pytest.raises(ValueError, match="0/1"):
        validate_inputs(x, np.arange(18), covariance, "binomial")
    with pytest.raises(ValueError, match="integer"):
        validate_inputs(x, np.arange(18) + 0.1, covariance, "negative-binomial")


@pytest.mark.parametrize("predictor_file", [False, True])
def test_cli_bundle_and_input_protection(tmp_path, monkeypatch, predictor_file):
    x, y, _, groups = example()
    names = [f"s{i}" for i in range(len(y))]
    tree = tmp_path / "tree.nwk"
    tree.write_text("(" + ",".join(f"{name}:1" for name in names) + ");")
    data = tmp_path / "data.tsv"
    columns = (
        ["遺伝子1", "OG2", "OG3", "OG4"]
        if predictor_file
        else ["OG1", "OG2", "OG3", "OG4"]
    )
    table = pd.DataFrame(x, columns=columns)
    table.insert(0, "leaf_name", names)
    table["trait"] = y
    table.to_csv(data, sep="\t", index=False)
    folds = tmp_path / "folds.tsv"
    pd.DataFrame(dict(leaf_name=names, fold=groups)).to_csv(
        folds, sep="\t", index=False
    )
    prefix = tmp_path / "result"
    args = [
        "regress-select",
        "--input-rooted",
        "yes",
        "--tree",
        str(tree),
        "--data",
        str(data),
        "--folds",
        str(folds),
        "--response",
        "trait",
        "--predictors",
        ",".join(columns),
        "--strengths",
        "0.2,0.05",
        "--l1-ratios",
        "0.5",
        "--out-prefix",
        str(prefix),
    ]
    if predictor_file:
        path = tmp_path / "predictors.txt"
        path.write_text("\n".join(columns), encoding="utf-8")
        original = Path.read_text

        def locale_read(source, encoding=None, errors=None, **kwargs):
            if source == path and encoding is None:
                encoding = "cp1252"
            return original(source, encoding=encoding, errors=errors, **kwargs)

        monkeypatch.setattr(Path, "read_text", locale_read)
        index = args.index("--predictors")
        args[index : index + 2] = ["--predictor-file", str(path)]
    audit = tmp_path / "audit.json"
    assert main([*args, "--audit", str(audit)]) == 0
    record = json.loads(audit.read_text())
    assert len(record["outputs"]) == 6
    assert {row["argument"] for row in record["inputs"]} == (
        {"tree", "data", "folds", "predictor_file"}
        if predictor_file
        else {"tree", "data", "folds"}
    )
    assert record["primary_input"]["kind"] == "newick"
    coefficients = pd.read_csv(str(prefix) + ".coefficients.tsv", sep="\t")
    assert coefficients.loc[coefficients.term == columns[0], "selected"].item()
    assert not any("p_value" in col or "pval" in col for col in coefficients)
    metadata = json.loads((tmp_path / "result.metadata.json").read_text())
    assert metadata["n_folds"] == 3
    assert metadata["nested_cv_mean_loss"] >= 0
    assert metadata["nested_cv_mean_loss"] < metadata["nested_cv_baseline_mean_loss"]
    with pytest.raises(ValueError, match="rooted"):
        unrooted_args = list(args)
        unrooted_args[unrooted_args.index("--input-rooted") + 1] = "no"
        main(unrooted_args)
    saved = (tmp_path / "result.coefficients.tsv").read_bytes()
    with pytest.raises(ValueError):
        main([*args, "--audit", str(tmp_path / "result.coefficients.tsv")])
    assert (tmp_path / "result.coefficients.tsv").read_bytes() == saved
    protected = tmp_path / "result.cv.tsv"
    protected.write_text(data.read_text())
    args[args.index("--data") + 1] = str(protected)
    with pytest.raises(ValueError, match="input"):
        main(args)
    assert protected.read_text() == data.read_text()


@pytest.mark.parametrize("invalid_mode", [False, True])
def test_scalar_mode_verifies_score_instead_of_solver_flag(monkeypatch, invalid_mode):
    import nwkit.phylogenetic_glmm as glmm

    cholesky = glmm._positive_definite_cholesky
    minimize = glmm.minimize
    calls = 0

    def force_fallback(matrix):
        nonlocal calls
        calls += 1
        if calls == 1:
            raise np.linalg.LinAlgError("exercise fallback")
        return cholesky(matrix)

    def replace_flag(*args, **kwargs):
        result = minimize(*args, **kwargs)
        result.success = invalid_mode
        if invalid_mode:
            result.x[:] = np.nan
        result.message = "test optimizer flag"
        return result

    monkeypatch.setattr(glmm, "_positive_definite_cholesky", force_fallback)
    monkeypatch.setattr(glmm, "minimize", replace_flag)
    args = (
        np.array([4.0, 2.0, 1.0, 4.0, 1.0, 1.0]),
        np.linspace(0.5, 0.8, 6),
        np.eye(6) * 12.85,
        np.arange(6),
    )
    kwargs = dict(
        family="negative-binomial",
        dispersion=0.071,
        zero_probability=None,
        offset=np.zeros(6),
        trials=None,
        censor_lower=None,
        censor_upper=None,
    )
    if invalid_mode:
        with (
            np.errstate(invalid="ignore"),
            pytest.raises(RuntimeError, match="mode optimization"),
        ):
            glmm._scalar_random_mode(*args, **kwargs)
    else:
        mode, value, _ = glmm._scalar_random_mode(*args, **kwargs)
        assert np.isfinite(mode).all() and np.isfinite(value)


@pytest.mark.parametrize("family", ["poisson", "negative-binomial", "binomial"])
def test_small_penalty_matches_existing_glmm(family):
    from nwkit.phylogenetic_glmm import fit_phylogenetic_glmm

    x, _, covariance, _ = example()
    x = x[:, :1]
    rng = np.random.default_rng(12)
    y = (
        rng.binomial(1, expit(x[:, 0]))
        if family == "binomial"
        else rng.negative_binomial(2, 0.45, 18)
    )
    design = np.column_stack([np.ones(len(y)), x])
    response = y.astype(str) if family == "binomial" else y
    extra = dict(reference="0", levels=["1", "0"]) if family == "binomial" else {}
    reference = fit_phylogenetic_glmm(
        response,
        design,
        lambda _: covariance,
        family=family,
        coefficient_penalty="none",
        **extra,
    )
    fitted = fit_elastic_net(x, y, covariance, family=family, strength=1e-7)
    np.testing.assert_allclose(
        fitted.coefficients, reference.coefficients.ravel(), atol=0.015
    )


@pytest.mark.parametrize(
    "family", ["gaussian", "binomial", "poisson", "negative-binomial"]
)
def test_simulated_strong_phylogeny_with_more_predictors_than_tips(family):
    rng = np.random.default_rng(73)
    groups = np.repeat(np.arange(3), 8)
    covariance = np.eye(24) + 5 * (groups[:, None] == groups)
    x = rng.normal(size=(24, 32))
    random = rng.multivariate_normal(np.zeros(24), covariance * 0.15)
    linear = 0.5 * x[:, 0] + random
    if family == "gaussian":
        y = linear
    elif family == "binomial":
        y = rng.binomial(1, expit(linear))
    elif family == "poisson":
        y = rng.poisson(np.exp(linear))
    else:
        y = rng.negative_binomial(2, 2 / (2 + np.exp(linear)))
    fitted = fit_elastic_net(
        x, y, covariance, family=family, strength=0.2, l1_ratio=0.5
    )
    assert fitted.projected_gradient < 2e-4
    assert np.isfinite(fitted.coefficients).all()


def test_nested_selection_is_invariant_to_branch_length_units():
    x, y, covariance, groups = example()
    kwargs = dict(strengths=(0.2,), l1_ratios=(0.5,))
    fitted, _, _, _, prediction = nested_selection(x, y, covariance, groups, **kwargs)
    rescaled, _, _, _, rescaled_prediction = nested_selection(
        x, y, covariance * 1e6, groups, **kwargs
    )
    np.testing.assert_allclose(fitted.coefficients, rescaled.coefficients, atol=1e-6)
    np.testing.assert_allclose(
        prediction.predicted, rescaled_prediction.predicted, atol=1e-6
    )


@pytest.mark.parametrize("ratio", [1.0, 0.5])
def test_gaussian_lasso_flags_rank_deficient_equicorrelation_design(ratio):
    x = np.tile(np.array([-1, 1] * 4, dtype=float)[:, None], (1, 2))
    fit = fit_elastic_net(x, 0.3 + 2 * x[:, 0], np.eye(8), l1_ratio=ratio)
    assert (fit.boundary_warning == "lasso_equicorrelation_design_rank_deficient") == (
        ratio == 1
    )


def test_gaussian_iteration_budget_does_not_waive_kkt_check():
    rng = np.random.default_rng(820)
    x = rng.normal(size=(20, 5))
    x[:, 1] += 3 * x[:, 0]
    y = x @ np.array([1.2, -0.8, 0.5, 0.4, 0.3])
    with pytest.raises(RuntimeError, match="coordinates did not converge"):
        fit_elastic_net(x, y, np.eye(20), strength=0.01, maxiter=1)
    fitted = fit_elastic_net(x, y, np.eye(20), strength=0.01)
    assert fitted.projected_gradient <= 1e-6
