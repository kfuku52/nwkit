"""Numerical, calibration and CLI contracts independent of the R backend."""

import json

import numpy as np
import pytest
from scipy.stats import multivariate_normal

from nwkit.cli import main
from nwkit.shift_calibration import CalibratedSearch
from nwkit.util import read_tree

TREE = "(((A:1,B:1):1,(C:1,D:1):1):1,((E:1,F:1):1,(G:1,H:1):1):1);"
Y = np.array([0.2, -0.4, 1.1, 0.9, -0.1, 0.3, 1.5, 0.7])


def engine(**kwargs):
    return CalibratedSearch(read_tree(TREE, "auto", True, quiet=True), **kwargs)


def test_covariance_limits_and_root_invariance():
    s = engine()
    bm, bm_weights = s.geometry(0)
    iid, iid_weights = s.geometry(np.inf)
    small, sw = s.geometry(1e-9)
    large, lw = s.geometry(1e4)
    np.testing.assert_allclose(small, bm, atol=1e-9)
    np.testing.assert_allclose(sw, bm_weights, atol=1e-9)
    np.testing.assert_allclose(large, iid, atol=1e-9)
    np.testing.assert_allclose(lw, iid_weights, atol=1e-9)
    a = 0.7
    fixed, _ = s.geometry(a)
    random = np.exp(-a * s.distance) / -np.expm1(-2 * a)
    np.testing.assert_allclose(s.q @ fixed @ s.q.T, s.q @ random @ s.q.T, atol=1e-12)


def test_known_error_likelihood_matches_direct_gaussian():
    variance = np.linspace(0.01, 0.08, 8)
    s = engine(variances=variance, alpha_grid=[0, 0.7, np.inf])
    z = s.q @ Y
    for item in s.cache[::37]:
        score, beta, _ = s._at(z[:, None], item)
        L, W = item[2:4]
        for m in (0, 5, len(s.models) - 1):
            mean = L @ W[m] @ beta[m, :, 0]
            expected = multivariate_normal.logpdf(z, mean=mean, cov=L @ L.T)
            assert score[m, 0] == pytest.approx(expected, rel=1e-10, abs=1e-9)


def test_reproducibility_scale_and_zero_error_equivalence():
    s = engine()
    result = s.fit(Y, seed=42, replicates=19)
    assert result == engine(variances=np.zeros(8)).fit(Y, seed=42, replicates=19)
    scaled = s.fit(3 * Y + 10, seed=42, replicates=19)
    assert result["model"] == scaled["model"]
    assert result["tests"][0]["p_value"] == scaled["tests"][0]["p_value"]
    np.testing.assert_allclose(
        scaled["predicted"], np.array(result["predicted"]) * 3 + 10
    )
    assert scaled["process_tip_variance"] == pytest.approx(
        result["process_tip_variance"] * 9
    )


@pytest.mark.parametrize(
    "alpha,status", [(0, "brownian_limit"), (np.inf, "independent_limit")]
)
def test_exact_limits_have_predictions_but_no_finite_optima(alpha, status):
    fit = engine(alpha_grid=[alpha]).fit(Y, replicates=19)
    assert fit["alpha_status"] == status
    assert fit["optimum_identifiable"] is False
    assert fit["optimum_effects"] == [None, None]
    assert np.all(np.isfinite(fit["predicted"]))
    json.dumps(fit, allow_nan=False)


def test_bootstrap_repeats_entire_search(monkeypatch):
    s = engine()
    original = s.profile
    calls = []

    def record(z):
        calls.append(np.asarray(z).shape)
        return original(z)

    monkeypatch.setattr(s, "profile", record)
    fit = s.fit(Y, replicates=19)
    assert calls[0] == (7,)
    assert len(calls) == 1 + sum(
        len(test.get("null_alpha_evaluations", [None])) for test in fit["tests"]
    )
    assert all(shape == (7, 19) for shape in calls[1:])
    assert all(test["p_value"] >= 1 / 20 for test in fit["tests"])


def test_degenerate_data_and_invalid_calibration_rejected():
    s = engine()
    with pytest.raises(ValueError, match="Zero residual"):
        s.fit(np.ones(8), replicates=19)
    with pytest.raises(ValueError, match="resolution"):
        s.fit(Y, replicates=19, level=0.01)


def test_cli_default_calibrated_atomic_outputs(tmp_path):
    tree = tmp_path / "tree.nwk"
    tree.write_text(TREE)
    trait = tmp_path / "trait.tsv"
    trait.write_text(
        "leaf_name\tx\n"
        + "".join(f"{name}\t{y}\n" for name, y in zip("ABCDEFGH", Y, strict=True))
    )
    model = tmp_path / "model.json"
    regimes = tmp_path / "regimes.tsv"
    command = [
        "shift",
        "-i",
        str(tree),
        "--trait",
        str(trait),
        "--state-column",
        "x",
        "--model-out",
        str(model),
        "-o",
        str(regimes),
        "--calibration-replicates",
        "19",
        "--convergence",
    ]
    main(command)
    result = json.loads(model.read_text())
    assert result["schema_version"] == 7
    assert result["selection"] == "calibrated"
    assert result["calibration"]["replicates"] == 19
    assert "contrast_log_likelihood" in result["parameters"]
    assert "score" not in result["parameters"]
    prior = model.read_bytes(), regimes.read_bytes()
    with pytest.raises(ValueError, match="IC backend"):
        main([*command, "--fit-out", str(tmp_path / "fit.rds")])
    assert prior == (model.read_bytes(), regimes.read_bytes())
    assert not (tmp_path / "fit.rds").exists()


def test_zero_process_variance_with_exact_observations():
    one_exact = engine(variances=np.r_[0, np.ones(7)], alpha_grid=[0])
    assert any(item[1] == 0 for item in one_exact.cache)
    fit = one_exact.fit(Y * 0.001, replicates=19)
    assert fit["process_tip_variance"] == 0
    assert fit["predicted"][0] == pytest.approx(Y[0] * 0.001)
    with pytest.raises(ValueError, match="unbounded"):
        engine(variances=np.r_[0, 0, np.ones(6)], alpha_grid=[0])


def test_fixed_parameter_mean_design_matches_independent_generator():
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
    from shift_simulation_cases import simulate

    newick, truth = simulate(
        tips=8,
        scenario="convergent",
        root_model="OUfixedRoot",
        se=0.2,
        effect=2,
        seed=482,
        alpha=0.7,
        sigma2=0.25,
    )
    s = CalibratedSearch(read_tree(newick, "auto", True, quiet=True))
    a = truth["alpha"] * s.height
    K, weights = s.geometry(a)
    v = truth["sigma2"] * -np.expm1(-2 * a) / (2 * truth["alpha"])
    np.testing.assert_allclose(
        v * K + 0.04 * np.eye(8), truth["covariance"], atol=1e-12
    )
    m = next(
        i
        for i, model in enumerate(s.models)
        if model["shift_branch_ids"] == truth["shift_branch_ids"]
        and len(model["groups"]) == 2
    )
    coefficients = np.array([2 * -np.expm1(-a), 0])
    mean = s.loads[m] @ (weights[m] * (s.transform[m] @ coefficients))
    np.testing.assert_allclose(mean, list(truth["tip_mean"].values()), atol=1e-12)


def test_limit_counts_only_enabled_candidates():
    s = engine(convergence=False)
    at_limit = engine(convergence=False, limit=len(s.models))
    assert at_limit.models == s.models
    assert [m["candidate_id"] for m in s.models] == list(range(len(s.models)))
    with pytest.raises(ValueError, match="limit"):
        engine(convergence=False, limit=len(s.models) - 1)


def test_short_branches_agree_with_svd_least_squares():
    e = 1e-12
    newick = f"(((A:{e},B:{e}):{1 - e},(C:{e},D:{e}):{1 - e}):1,((E:{e},F:{e}):{1 - e},(G:{e},H:{e}):{1 - e}):1);"
    s = CalibratedSearch(read_tree(newick, "auto", True, quiet=True))
    z = s.q @ np.array([0.2, 0.1, 0.5, 0.6, 1, 2, 3, 4])
    for item in s.cache:
        _, _, rss = s._at(z[:, None], item)
        L, W = item[2:4]
        w = np.linalg.solve(L, z)
        for index, dimension in enumerate(s.dim):
            design = W[index, :, :dimension]
            beta = np.linalg.lstsq(design, w, rcond=None)[0]
            residual = w - design @ beta
            assert rss[index, 0] == pytest.approx(residual @ residual, rel=1e-10)


def test_calibration_batches_preserve_unbatched_rng_and_probability(monkeypatch):
    s = engine()
    family = s.families[0]
    L = s.cache[0][2]
    mean = np.arange(s.d) / 10
    B, scale, observed = 137, 0.8, 7.0
    noise = np.random.default_rng(493).normal(size=(s.d, B))
    sb, _ = s.profile(mean[:, None] + np.sqrt(scale) * L @ noise)
    expected = (
        1
        + np.count_nonzero(
            2 * (sb.max(axis=0) - sb[family].max(axis=0)) >= observed - 1e-10
        )
    ) / (B + 1)
    calls = []
    original = s.profile

    def record(z):
        calls.append(z.shape[1])
        return original(z)

    monkeypatch.setattr(s, "profile", record)
    actual = s._calibrated_probability(
        mean, L, scale, family, observed, np.random.default_rng(493), B
    )
    assert actual == expected
    assert calls == [64, 64, 9]


@pytest.mark.parametrize("replicates", [0, -1, 19.5, True])
def test_invalid_replicate_count_is_rejected(replicates):
    with pytest.raises(ValueError, match="integer"):
        engine().fit(Y, replicates=replicates)


def test_known_error_extreme_unit_rescaling():
    s = engine(variances=np.full(8, 0.04), alpha_grid=[0, 1, np.inf])
    fit = s.fit(Y, replicates=19)
    tiny = engine(variances=np.full(8, 0.04e-200), alpha_grid=[0, 1, np.inf]).fit(
        Y * 1e-100, replicates=19
    )
    assert fit["model"] == tiny["model"]
    assert [r["p_value"] for r in fit["tests"]] == [r["p_value"] for r in tiny["tests"]]
    np.testing.assert_allclose(
        np.array(tiny["predicted"]) * 1e100, fit["predicted"], rtol=1e-10, atol=1e-10
    )


@pytest.mark.parametrize("observed", [0.0, 1e6])
def test_null_envelope_bounds_full_grid_maximum_and_preserves_decision(observed):
    s = engine(alpha_grid=[0, 1, np.inf])
    seed, B, level = 821, 19, 0.05
    noise = np.random.default_rng(seed).normal(size=(s.d, B))
    probabilities = []
    for item in s.cache:
        sb, _ = s.profile(item[2] @ noise)
        statistics = 2 * (sb.max(axis=0) - sb[0])
        probabilities.append(
            (1 + np.count_nonzero(statistics >= observed - 1e-10)) / (B + 1)
        )
    exact = max(probabilities)
    reported, info = s._null_probability(
        observed, np.random.default_rng(seed), B, level
    )
    assert reported >= exact
    assert (reported <= level) == (exact <= level)
    if reported <= level:
        assert info["p_value_kind"] == "grid_supremum"
        assert len(info["null_alpha_evaluations"]) == len(s.grid)
        assert reported == exact
    else:
        assert info["p_value_lower_bound"] > level
        assert info["p_value_kind"] == "conservative_upper_bound"


def test_unsaturated_exact_observations_can_fit_but_not_a_clipped_lower_boundary():
    s = engine(variances=np.r_[np.zeros(4), np.ones(4)], alpha_grid=[0])
    fit = s.fit(Y, replicates=19)
    assert fit["process_tip_variance"] > s.variance_grid[1]
    null_only = engine(max_shifts=0, variances=np.r_[0, 0, np.ones(6)], alpha_grid=[0])
    with pytest.raises(ValueError, match="singular lower grid boundary"):
        null_only.fit(np.r_[0.2, 0.2, Y[2:]], replicates=19)


@pytest.mark.parametrize("exact_count,max_shifts", [(2, 1), (2, 2), (3, 2)])
def test_saturated_mixed_error_candidates_rejected_before_bootstrap(
    exact_count, max_shifts
):
    with pytest.raises(ValueError, match="unbounded"):
        engine(
            max_shifts=max_shifts,
            variances=np.r_[np.zeros(exact_count), np.ones(8 - exact_count)],
        )


@pytest.mark.parametrize("alpha_height", [0.0, 0.7, np.inf])
def test_every_candidate_mean_matches_independent_branch_recursion(alpha_height):
    from nwkit.util import assign_branch_ids

    s = engine()
    ids = assign_branch_ids(s.tree)
    _, weights = s.geometry(alpha_height)
    baseline = 2.3
    for index, model in enumerate(s.models):
        aliases = {
            branch: group
            for group, branches in enumerate(model["groups"])
            for branch in branches
        }
        base = aliases[0]
        group_values = {base: 0.0}
        beta = np.zeros(2)
        for j, group in enumerate(g for g in range(len(model["groups"])) if g != base):
            group_values[group] = (-1.0) ** j * (j + 1.2)
            beta[j] = group_values[group] * (
                1
                if alpha_height == 0 or np.isinf(alpha_height)
                else -np.expm1(-alpha_height)
            )
        model_mean = baseline + s.loads[index] @ (
            weights[index] * (s.transform[index] @ beta)
        )
        means, regime = {s.tree: baseline}, {s.tree: base}
        for node in s.tree.traverse("preorder"):
            if node.is_root:
                continue
            regime[node] = aliases.get(ids[node], regime[node.up])
            offset = group_values[regime[node]]
            if alpha_height == 0:
                means[node] = means[node.up] + node.dist / s.height * offset
            elif np.isinf(alpha_height):
                means[node] = baseline + offset
            else:
                decay = np.exp(-alpha_height / s.height * node.dist)
                means[node] = decay * means[node.up] + (1 - decay) * (baseline + offset)
        np.testing.assert_allclose(
            model_mean, [means[t] for t in s.tree.leaves()], rtol=1e-12, atol=1e-12
        )
