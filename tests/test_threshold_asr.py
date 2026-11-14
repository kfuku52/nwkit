import math

import numpy as np
import pandas as pd
import pytest
from scipy.stats import truncnorm

from nwkit.cli import main
from nwkit.compiled_tree import CompiledTree
from nwkit.threshold_asr import (
    _allowed_intervals,
    _initialize_values,
    _rhat_and_ess,
    _truncated_draw,
    _truncated_quantile,
    _update_thresholds,
    build_threshold_process,
    compute_threshold_marginals,
    parse_thresholds,
    threshold_liability_table,
)
from nwkit.util import read_tree


def tree_from(source):
    return read_tree(source, "1", True, quiet=True, rooted="yes")


def likelihoods(tree, states, observed):
    result = {}
    for leaf in tree.leaves():
        state = observed[str(leaf.name)]
        result[str(leaf.name)] = (
            np.ones(len(states), dtype=float)
            if state is None
            else np.asarray([float(state == candidate) for candidate in states])
        )
    return result


def test_parse_thresholds_validates_count_and_order():
    assert parse_thresholds("-1,2", 3) == (-1.0, 2.0)
    with pytest.raises(ValueError, match="requires 2"):
        parse_thresholds("0", 3)
    with pytest.raises(ValueError, match="strictly increasing"):
        parse_thresholds("1,1", 3)


def test_threshold_process_is_identified_unit_brownian():
    tree = tree_from("(A:1,B:2)R;")
    process = build_threshold_process(tree)
    assert process.root.mode == "gaussian"
    assert process.root.mean == 0.0
    assert process.root.variance == 1.0
    assert process.transitions[tree["A"]].variance == 1.0
    assert process.transitions[tree["B"]].variance == 2.0


def test_liability_initialization_stays_inside_extreme_finite_thresholds():
    tree = tree_from("(A:1,B:2)R;")
    compiled = CompiledTree.from_tree(tree)
    index = compiled.index_by_node[tree["A"]]
    smallest_positive = np.nextafter(0.0, 1.0)
    for lower, upper in ((1e308, 1.7e308), (0.0, smallest_positive)):
        values = _initialize_values(
            compiled, {index: np.array([1])}, np.array([lower, upper])
        )
        assert lower < values[index] <= upper


def test_disjoint_ordinal_intervals_keep_equal_mass_in_both_far_tails():
    rng = np.random.default_rng(23)
    allowed = np.array([0, 2])
    thresholds = np.array([-40.0, 40.0])
    bounds = _allowed_intervals(allowed, thresholds)
    draws = np.array([_truncated_draw(0.0, 1.0, bounds, rng) for _ in range(300)])
    assert np.all(np.abs(draws) >= 40.0)
    assert 0.4 < np.mean(draws > 0.0) < 0.6


def test_threshold_draw_rounds_to_first_allowed_value_in_extreme_tail():
    rng = np.random.default_rng(5)
    result = _truncated_draw(0.0, 1.0, [(1e200, math.inf)], rng)
    assert result == np.nextafter(1e200, np.inf)


def test_disjoint_ordinal_intervals_compare_tails_beyond_log_probability_range():
    bounds = [(-math.inf, -1e200), (1e200, math.inf)]
    rng = np.random.default_rng(35)
    symmetric = [_truncated_draw(0.0, 1.0, bounds, rng) for _ in range(500)]
    shifted = [_truncated_draw(1e-200, 1.0, bounds, rng) for _ in range(1500)]
    assert np.all(np.abs(symmetric) >= 1e200)
    assert np.mean(np.asarray(symmetric) > 0.0) == pytest.approx(0.5, abs=0.07)
    # The normalized ratio between the tails tends to exp(2 * mean * boundary).
    assert np.mean(np.asarray(shifted) > 0.0) == pytest.approx(
        1.0 / (1.0 + math.exp(-2.0)), abs=0.03
    )


def test_disjoint_ordinal_intervals_compare_extreme_tail_widths_after_centering():
    width = 1e-200
    bounds = [(0.0, width), (2 * width, 4 * width)]
    rng = np.random.default_rng(19)
    draws = np.array([_truncated_draw(-1e16, 1.0, bounds, rng) for _ in range(1200)])
    first = (draws > 0.0) & (draws <= width)
    second = (draws > 2 * width) & (draws <= 4 * width)
    assert np.all(first | second)
    assert np.mean(second) == pytest.approx(2.0 / 3.0, abs=0.04)


def test_threshold_draw_reports_tail_with_no_finite_allowed_value():
    rng = np.random.default_rng(5)
    with pytest.raises(
        ValueError, match="ordinal liability exceeds floating-point range"
    ):
        _truncated_draw(0.0, 1.0, [(float(np.finfo(float).max), math.inf)], rng)


@pytest.mark.parametrize("sign", [-1.0, 1.0])
def test_threshold_draw_resolves_tail_excess_after_centering(sign):
    rng = np.random.default_rng(331)
    bounds = [(0.0, math.inf)] if sign > 0.0 else [(-math.inf, 0.0)]
    draws = (
        np.array([_truncated_draw(-sign * 1e16, 1.0, bounds, rng) for _ in range(1500)])
        * sign
        * 1e16
    )
    assert np.all(draws > 0.0)
    assert np.mean(draws) == pytest.approx(1.0, abs=0.07)
    assert np.median(draws) == pytest.approx(math.log(2.0), abs=0.07)


@pytest.mark.parametrize("mean", [-1.0, 0.0])
def test_threshold_draw_preserves_width_when_cdf_endpoints_coincide(mean):
    width = 1e-200
    rng = np.random.default_rng(3)
    sample = np.array(
        [_truncated_draw(mean, 1.0, [(0.0, width)], rng) for _ in range(400)]
    )
    assert np.all((sample > 0.0) & (sample <= width))
    assert np.mean(sample / width) == pytest.approx(0.5, abs=0.04)


def test_disjoint_narrow_ordinal_intervals_are_selected_in_proportion_to_width():
    width = 1e-200
    thresholds = width * np.array([0, 1, 2, 4])
    allowed = np.array([1, 3])
    rng = np.random.default_rng(7)
    bounds = _allowed_intervals(allowed, thresholds)
    sample = np.array([_truncated_draw(-1.0, 1.0, bounds, rng) for _ in range(800)])
    first = (sample > 0.0) & (sample <= width)
    second = (sample > 2 * width) & (sample <= 4 * width)
    assert np.all(first | second)
    assert np.mean(second) == pytest.approx(2 / 3, abs=0.04)


@pytest.mark.parametrize(
    ("lower", "upper"),
    [(-np.inf, -40.0), (40.0, 40.01), (-2.0, 3.0), (0.0, np.inf)],
)
def test_truncated_quantile_matches_scipy_at_fixed_seed(lower, upper):
    rng = np.random.default_rng(2916)
    expected_rng = np.random.default_rng(2916)
    actual = [_truncated_quantile(lower, upper, rng) for _ in range(50)]
    expected = truncnorm.rvs(lower, upper, size=50, random_state=expected_rng)
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-12)
    assert rng.random() == expected_rng.random()


def test_liability_sampling_accepts_a_positive_subnormal_branch_length():
    tree = tree_from("(A:1e-310,B:1,C:2)R;")
    states = ("low", "high")
    observed = {"A": "low", "B": "high", "C": None}
    with np.errstate(over="raise", divide="raise", invalid="raise"):
        posterior, fit = compute_threshold_marginals(
            tree,
            states,
            observed,
            likelihoods(tree, states, observed),
            num_samples=20,
            burnin=5,
            chains=1,
            seed=7,
        )
    np.testing.assert_array_equal(posterior[tree["A"]], [1.0, 0.0])
    np.testing.assert_array_equal(posterior[tree["B"]], [0.0, 1.0])
    assert all(np.isfinite(entry.mean) for entry in fit.liability_marginals.values())


@pytest.mark.parametrize("threshold", [1e5, 1e8])
def test_large_fixed_threshold_preserves_category_and_liability_variance(threshold):
    tree = tree_from("(A:1,B:1,C:1)R;")
    states = ("low", "high")
    observed = {"A": "low", "B": "high", "C": None}
    posterior, fit = compute_threshold_marginals(
        tree,
        states,
        observed,
        likelihoods(tree, states, observed),
        thresholds=str(threshold),
        num_samples=100,
        burnin=20,
        chains=2,
        seed=5,
    )
    np.testing.assert_array_equal(posterior[tree["B"]], [0.0, 1.0])
    marginal = fit.liability_marginals[tree["B"]]
    assert marginal.mean > threshold
    assert 0.0 < marginal.variance < (10.0 / threshold) ** 2


def test_binary_threshold_samples_probabilities_and_liabilities():
    tree = tree_from("((A:1,B:1):1,(C:1,D:1):1)R;")
    states = ("low", "high")
    observed = {"A": "low", "B": "low", "C": "high", "D": None}
    posterior, fit = compute_threshold_marginals(
        tree,
        states,
        observed,
        likelihoods(tree, states, observed),
        num_samples=120,
        burnin=60,
        chains=2,
        seed=17,
    )
    assert posterior[tree["A"]] == pytest.approx([1.0, 0.0])
    assert posterior[tree["C"]] == pytest.approx([0.0, 1.0])
    assert sum(posterior[tree["D"]]) == pytest.approx(1.0)
    assert 0.0 < posterior[tree][0] < 1.0
    assert fit.thresholds == pytest.approx((0.0,))
    assert fit.thresholds_estimated is False
    assert len(threshold_liability_table(tree, fit)) == 7


def test_fixed_ordinal_thresholds_are_preserved_exactly():
    tree = tree_from("(A:1,B:1,C:1)R;")
    states = ("a", "b", "c", "d", "e", "f")
    observed = {"A": "c", "B": "d", "C": None}
    thresholds = (-1e308, -1.0, 0.3, 1.0, 1e308)
    with np.errstate(over="raise", invalid="raise"):
        posterior, fit = compute_threshold_marginals(
            tree,
            states,
            observed,
            likelihoods(tree, states, observed),
            thresholds=",".join(map(str, thresholds)),
            num_samples=12,
            burnin=4,
            chains=2,
            seed=7,
        )
    assert fit.thresholds == thresholds
    np.testing.assert_array_equal(posterior[tree["A"]], [0, 0, 1, 0, 0, 0])
    np.testing.assert_array_equal(posterior[tree["B"]], [0, 0, 0, 1, 0, 0])


def test_threshold_allows_zero_burnin_and_marks_single_chain_rhat_unavailable():
    tree = tree_from("(A:1,B:1,C:1)R;")
    states = ("low", "high")
    observed = {"A": "low", "B": "high", "C": "low"}
    _, fit = compute_threshold_marginals(
        tree,
        states,
        observed,
        likelihoods(tree, states, observed),
        num_samples=10,
        burnin=0,
        chains=1,
        seed=3,
    )
    assert fit.burnin == 0
    assert np.isnan(fit.rhat_max)
    assert "mcmc_rhat_unavailable" in fit.fit_status


def test_threshold_ess_detects_dependence_beyond_the_first_lag():
    rng = np.random.default_rng(2874)
    traces = rng.normal(size=(4, 4096, 1))
    for sample in range(2, traces.shape[1]):
        traces[:, sample, :] += 0.75 * traces[:, sample - 2, :]
    rhat, ess = _rhat_and_ess(traces)
    # Independent even and odd AR(1) subchains have zero lag-one correlation.
    # Their integrated autocorrelation time is (1 + 0.75) / (1 - 0.75).
    assert rhat[0] == pytest.approx(1.0, abs=0.01)
    assert ess[0] == pytest.approx(traces.size / 7.0, rel=0.2)
    scaled_rhat, scaled_ess = _rhat_and_ess(traces * 1e200)
    np.testing.assert_allclose(scaled_rhat, rhat)
    np.testing.assert_allclose(scaled_ess, ess)


def test_threshold_ess_rejects_apparent_agreement_of_constant_separated_chains():
    traces = np.broadcast_to(np.arange(4.0)[:, None, None], (4, 100, 1))
    rhat, ess = _rhat_and_ess(traces)
    assert np.isinf(rhat[0])
    assert ess[0] == pytest.approx(400.0 / 199.0)


def test_threshold_liability_variance_stays_finite_when_sum_of_squares_overflows():
    variance = 1.7e308
    tree = read_tree(
        f"(A:{variance},B:{variance})R;", "1", True, rooted="yes", quiet=True
    )
    with np.errstate(over="raise", invalid="raise"):
        _, fit = compute_threshold_marginals(
            tree,
            ("low", "high"),
            {"A": "low", "B": "high"},
            {"A": np.array([1.0, 0.0]), "B": np.array([0.0, 1.0])},
            thresholds="0",
            num_samples=500,
            burnin=30,
            chains=2,
            seed=4,
        )

    for name, sign in (("A", -1), ("B", 1)):
        marginal = fit.liability_marginals[tree[name]]
        assert marginal.mean / math.sqrt(variance) == pytest.approx(
            sign * math.sqrt(2.0 / math.pi), abs=0.07
        )
        assert marginal.variance / variance == pytest.approx(
            1.0 - 2.0 / math.pi, abs=0.06
        )


def test_ordinal_cut_update_preserves_the_only_representable_valid_cut():
    first = math.ulp(0.0)
    second = math.nextafter(first, math.inf)
    values = np.array([first, second, 1e100])
    allowed = np.eye(3, dtype=bool)[:, [1, 2, 2]]
    thresholds = np.array([0.0, first])
    _update_thresholds(values, allowed, thresholds, np.random.default_rng(9))
    np.testing.assert_array_equal(thresholds, [0.0, first])
    np.testing.assert_array_equal(
        np.searchsorted(thresholds, values, side="left"), [1, 2, 2]
    )


def test_ordinal_thresholds_are_sampled_in_identified_order():
    tree = tree_from("((A:1,B:1):1,(C:1,D:1):1,E:2,F:2)R;")
    states = ("small", "medium", "large")
    observed = {
        "A": "small",
        "B": "small",
        "C": "medium",
        "D": "medium",
        "E": "large",
        "F": "large",
    }
    posterior, fit = compute_threshold_marginals(
        tree,
        states,
        observed,
        likelihoods(tree, states, observed),
        num_samples=100,
        burnin=80,
        chains=2,
        seed=4,
    )
    assert fit.thresholds_estimated is True
    assert fit.thresholds[0] == pytest.approx(0.0)
    assert fit.thresholds[1] > fit.thresholds[0]
    assert all(
        sum(probabilities) == pytest.approx(1.0) for probabilities in posterior.values()
    )


def test_estimated_thresholds_preserve_noncontiguous_ambiguous_tip_states():
    tree = tree_from("(A:1,B:1,C:1,D:1,E:1,F:1,G:1)R;")
    states = ("small", "medium", "large")
    observed = {
        "A": "small|large",
        "B": "small",
        "C": "small",
        "D": "medium",
        "E": "medium",
        "F": "large",
        "G": "large",
    }
    tip_likelihoods = {
        name: (
            np.asarray([1.0, 0.0, 1.0])
            if value == "small|large"
            else np.asarray([float(value == state) for state in states])
        )
        for name, value in observed.items()
    }
    posterior, _ = compute_threshold_marginals(
        tree,
        states,
        observed,
        tip_likelihoods,
        num_samples=120,
        burnin=60,
        chains=2,
        seed=9,
    )
    assert posterior[tree["A"]][1] == 0.0
    assert posterior[tree["A"]][0] + posterior[tree["A"]][2] == pytest.approx(1.0)


@pytest.mark.integration
def test_threshold_cli_writes_category_and_liability_outputs(tmp_path):
    traits = tmp_path / "traits.tsv"
    output = tmp_path / "asr.tsv"
    model = tmp_path / "model.tsv"
    liability = tmp_path / "liability.tsv"
    traits.write_text("leaf_name\tstate\nA\tlow\nB\tlow\nC\thigh\nD\tNA\n")
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
            "state",
            "--missing-values",
            "NA",
            "--model",
            "THRESHOLD",
            "--states",
            "low,high",
            "--liability-samples",
            "80",
            "--liability-burnin",
            "40",
            "--liability-chains",
            "2",
            "--seed",
            "9",
            "--model-out",
            str(model),
            "--liability-out",
            str(liability),
            "-o",
            str(output),
        ]
    )
    result = pd.read_csv(output, sep="\t")
    metadata = pd.read_csv(model, sep="\t").iloc[0]
    latent = pd.read_csv(liability, sep="\t")
    assert {"p_low", "p_high"}.issubset(result.columns)
    assert metadata["model"] == "THRESHOLD"
    assert metadata["liability_process"] == "standardized_brownian"
    assert metadata["mcmc_chains"] == 2
    assert len(latent) == 7


def test_threshold_cli_requires_explicit_state_order(tmp_path):
    traits = tmp_path / "traits.tsv"
    traits.write_text("leaf_name\tstate\nA\tlow\nB\thigh\n")
    with pytest.raises(ValueError, match="requires ordered categories"):
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
                "state",
                "--model",
                "THRESHOLD",
                "-o",
                str(tmp_path / "out.tsv"),
            ]
        )
