import json

import numpy as np
import pandas as pd
import pytest
from scipy.special import logsumexp
from scipy.stats import norm, poisson

from nwkit.cli import main
from nwkit.continuous_asr import compute_bm_marginals
from nwkit.jump_asr import fit_jump_bm
from nwkit.markov_gaussian_asr import fit_markov_gaussian, history_transition
from nwkit.util import read_tree


def tree():
    return read_tree("(A:1,B:1)R;", "1", True, quiet=True, rooted="yes")


def test_zero_jump_rate_matches_brownian_and_joint_draws():
    phylogeny = tree()
    observed = {"A": 0.0, "B": 2.0}
    fit = fit_jump_bm(
        phylogeny,
        observed,
        sigma2=1,
        jump_rate=0,
        jump_sd=4,
        history_samples=20,
        seed=1,
    )
    posterior, bm = compute_bm_marginals(phylogeny, observed, sigma2=1)
    assert fit.continuous_log_likelihood == pytest.approx(bm.restricted_log_likelihood)
    assert fit.effective_sample_size == pytest.approx(20)
    assert fit.relative_mcse < 1e-8
    assert fit.marginals[phylogeny].variance == pytest.approx(
        posterior[phylogeny].variance
    )
    lower, upper = fit.marginals[phylogeny].interval(0.95)
    assert [lower, upper] == pytest.approx(norm.ppf([0.025, 0.975], 1, np.sqrt(0.5)))
    samples = fit.sample(observed, None, 4000, seed=1)
    assert samples.values[:, 0].var() == pytest.approx(0.5, rel=0.05)
    for index, node in enumerate(samples.nodes):
        if node.is_leaf:
            assert samples.values[:, index] == pytest.approx(observed[node.name])


@pytest.mark.slow
def test_jump_importance_likelihood_matches_independent_poisson_sum():
    # Two unit edges collapse to a difference with Poisson(2*lambda) count.
    fit = fit_jump_bm(
        tree(),
        {"A": 0.0, "B": 4.0},
        sigma2=1,
        jump_rate=0.4,
        jump_sd=2,
        history_samples=12000,
        seed=4,
    )
    counts = np.arange(40)
    exact = logsumexp(
        poisson.logpmf(counts, 0.8) + norm.logpdf(4, scale=np.sqrt(2 + 4 * counts))
    )
    assert abs(fit.log_likelihood - exact) < 4 * fit.relative_mcse
    assert fit.effective_sample_size > 1000
    # Large tip difference shifts posterior jump counts upward.
    posterior_count = sum(
        weight * sum(history.values())
        for weight, history in zip(fit.weights, fit.latent, strict=True)
    )
    assert posterior_count > 0.8


def config(ou=False):
    result = {"states": ["x", "y"], "q": [[-0.4, 0.4], [0.2, -0.2]], "sigma2": [1, 1]}
    if ou:
        result.update(alpha=[0.5, 0.5], theta=[2, 2])
    return result


def test_equal_regime_rates_reduce_to_bm_and_preserve_tip_states():
    phylogeny = tree()
    observed = {"A": 0.0, "B": 2.0}
    fit = fit_markov_gaussian(
        phylogeny, observed, {"A": "x", "B": "y"}, config(), history_samples=40, seed=4
    )
    posterior, bm = compute_bm_marginals(phylogeny, observed, sigma2=1)
    assert fit.continuous_log_likelihood == pytest.approx(bm.restricted_log_likelihood)
    assert fit.log_likelihood == pytest.approx(
        fit.continuous_log_likelihood + fit.discrete_log_likelihood
    )
    assert fit.marginals[phylogeny].variance == pytest.approx(
        posterior[phylogeny].variance
    )
    for history in fit.latent:
        for node in phylogeny.leaves():
            assert history["node_states"][node] == (node.name == "B")
            segments = history["segments"][node]
            assert sum(duration for _, duration in segments) == pytest.approx(node.dist)
            assert segments[-1][0] == history["node_states"][node]


def test_ou_segments_compose_affine_moments():
    parameters = {"sigma2": [1, 3], "alpha": [0.2, 0.7], "theta": [-2, 4]}
    first = history_transition([(0, 1)], parameters)
    second = history_transition([(1, 2)], parameters)
    both = history_transition([(0, 1), (1, 2)], parameters)
    assert both.slope == pytest.approx(second.slope * first.slope)
    assert both.intercept == pytest.approx(
        second.slope * first.intercept + second.intercept
    )
    assert both.variance == pytest.approx(
        second.slope**2 * first.variance + second.variance
    )
    fit = fit_markov_gaussian(
        tree(),
        {"A": 0.0, "B": 2.0},
        {},
        config(True),
        model="MM-OU",
        history_samples=20,
        seed=4,
    )
    assert fit.effective_sample_size == pytest.approx(20)


@pytest.mark.parametrize("model", ["JUMP-BM", "MM-BM", "MM-OU"])
def test_latent_cli_outputs(tmp_path, model):
    traits = tmp_path / "traits.tsv"
    traits.write_text("leaf_name\tvalue\tregime\nA\t0\tx\nB\t2\ty\n")
    source = tmp_path / "config.json"
    source.write_text(json.dumps(config(model == "MM-OU")))
    output, metadata, history, draws = [
        tmp_path / name for name in ("out.tsv", "model.tsv", "history.tsv", "draws.tsv")
    ]
    options = (
        ["--sigma2", "1", "--jump-rate", "0.4", "--jump-sd", "2"]
        if model == "JUMP-BM"
        else ["--latent-regime-config", str(source), "--regime-column", "regime"]
    )
    main(
        [
            "asr",
            "-i",
            "(A:1,B:1)R;",
            "--input-rooted",
            "yes",
            "--trait",
            str(traits),
            "--state-column",
            "value",
            "--model",
            model,
            "--history-samples",
            "30",
            "--seed",
            "2",
            "-o",
            str(output),
            "--model-out",
            str(metadata),
            "--latent-history-out",
            str(history),
            "--posterior-samples-out",
            str(draws),
            "--posterior-samples",
            "3",
            *options,
        ]
    )
    row = pd.read_csv(metadata, sep="\t").iloc[0]
    assert row.model == model
    assert not row.rankable
    assert row.history_uncertainty_included
    assert len(pd.read_csv(draws, sep="\t")) == 9
    history_table = pd.read_csv(history, sep="\t")
    if model != "JUMP-BM":
        assert set(history_table.loc[history_table.branch_id == 1, "regime"]) == {"x"}
        assert set(history_table.loc[history_table.branch_id == 2, "regime"]) == {"y"}
        assert row.discrete_log_likelihood < 0
    assert len(history_table) == 30 * (2 if model == "JUMP-BM" else 3)
    assert not pd.read_csv(output, sep="\t").empty


def test_latent_invalid_parameters_and_resource_limits():
    with pytest.raises(ValueError, match="positive"):
        fit_jump_bm(tree(), {"A": 0}, sigma2=0, jump_rate=1, jump_sd=1)
    with pytest.raises(ValueError, match="history-node"):
        fit_jump_bm(
            tree(), {"A": 0}, sigma2=1, jump_rate=1, jump_sd=1, history_samples=400000
        )
    bad = config()
    bad["q"][0][0] = 1
    with pytest.raises(ValueError, match="row sums"):
        fit_markov_gaussian(tree(), {"A": 0}, {}, bad)


@pytest.mark.slow
def test_continuous_observations_reweight_unknown_regimes():
    # With Q=0 one root regime governs the whole tree. The analytic likelihood
    # is a 50:50 mixture of two Brownian difference densities.
    phylogeny = tree()
    specification = {
        "states": ["slow", "fast"],
        "q": [[0, 0], [0, 0]],
        "sigma2": [1, 9],
    }
    fit = fit_markov_gaussian(
        phylogeny,
        {"A": 0, "B": 5},
        {},
        specification,
        history_samples=2000,
        seed=12,
    )
    component_ll = norm.logpdf(5, scale=np.sqrt([2, 18]))
    expected_ll = logsumexp(component_ll) - np.log(2)
    expected_fast = np.exp(component_ll[1] - logsumexp(component_ll))
    actual_fast = sum(
        weight * (history["node_states"][phylogeny] == 1)
        for weight, history in zip(fit.weights, fit.latent, strict=True)
    )
    assert fit.log_likelihood == pytest.approx(expected_ll, abs=4 * fit.relative_mcse)
    assert actual_fast == pytest.approx(expected_fast, abs=0.015)
    assert actual_fast > 0.95


@pytest.mark.parametrize(
    "extra",
    [
        ["--history-samples", "0"],
        ["--posterior-predictive-out", "unused.tsv"],
        ["--tree-ensemble", "unused.nwk", "--tree-ensemble-out", "unused.tsv"],
    ],
)
def test_latent_cli_rejects_unsupported_options_before_outputs(tmp_path, extra):
    traits = tmp_path / "traits.tsv"
    traits.write_text("leaf_name\tvalue\nA\t0\nB\t1\n")
    output = tmp_path / "out.tsv"
    with pytest.raises(ValueError):
        main(
            [
                "asr",
                "-i",
                "(A:1,B:1)R;",
                "--input-rooted",
                "yes",
                "--trait",
                str(traits),
                "--state-column",
                "value",
                "--model",
                "JUMP-BM",
                "--sigma2",
                "1",
                "--jump-rate",
                "0",
                "--jump-sd",
                "1",
                "-o",
                str(output),
                *extra,
            ]
        )
    assert not output.exists()


def test_latent_history_output_cannot_overwrite_trait_input(tmp_path):
    traits = tmp_path / "traits.tsv"
    original = "leaf_name\tvalue\nA\t0\nB\t1\n"
    traits.write_text(original)
    with pytest.raises(ValueError, match="input"):
        main(
            [
                "asr",
                "-i",
                "(A:1,B:1)R;",
                "--trait",
                str(traits),
                "--state-column",
                "value",
                "--model",
                "JUMP-BM",
                "--sigma2",
                "1",
                "--jump-rate",
                "0",
                "--jump-sd",
                "1",
                "--latent-history-out",
                str(traits),
            ]
        )
    assert traits.read_text() == original
