import numpy as np
import pandas as pd
import pytest
from scipy.special import logsumexp

import nwkit.asr as asr
from nwkit.asr import (
    _discrete_gamma_category_rates,
    compute_mk_mixture_marginals,
)
from nwkit.cli import main
from nwkit.util import read_tree


def tree_from(source):
    return read_tree(source, "1", True, quiet=True, rooted="yes")


def character_data(patterns):
    states = ["0", "1"]
    observed_by_character = {}
    likelihoods_by_character = {}
    for character, pattern in patterns.items():
        observed = dict(zip("ABCDEF", pattern, strict=True))
        likelihoods = {}
        for leaf, state in observed.items():
            values = np.zeros(2, dtype=float)
            values[int(state)] = 1.0
            likelihoods[leaf] = values
        observed_by_character[character] = observed
        likelihoods_by_character[character] = likelihoods
    return states, observed_by_character, likelihoods_by_character


def test_discrete_gamma_uses_equal_probability_bin_means():
    rates = _discrete_gamma_category_rates(1.0, 4)
    assert rates == pytest.approx(
        [0.136953782644657, 0.476751856235452, 1.0, 2.386294361119891],
        abs=2e-15,
    )
    assert np.mean(rates) == pytest.approx(1.0)


@pytest.mark.parametrize("mixture", ["gamma", "free"])
def test_joint_mk_rate_mixtures_return_character_posteriors(mixture):
    tree = tree_from("(((A:1,B:1):1,C:1):1,((D:1,E:1):1,F:1):1)R;")
    states, observed, likelihoods = character_data(
        {
            "slow": "000111",
            "fast": "010101",
            "middle": "001101",
            "other": "110010",
        }
    )
    posterior, fit = compute_mk_mixture_marginals(
        tree,
        states,
        observed,
        likelihoods,
        base_model="ER",
        mixture=mixture,
        categories=2,
        rate=0.4,
        rate_bounds=(1e-4, 10.0),
    )
    assert fit["model"] == "MK-MIXTURE"
    assert fit["rate_mixture"] == mixture
    assert np.dot(
        fit["mixture_weights"], fit["mixture_category_rates"]
    ) == pytest.approx(1.0)
    assert set(posterior) == set(observed)
    for by_node in posterior.values():
        assert all(values.sum() == pytest.approx(1.0) for values in by_node.values())


def test_mk_mixture_cli_stacks_multiple_characters(tmp_path):
    trait = tmp_path / "traits.tsv"
    output = tmp_path / "asr.tsv"
    model = tmp_path / "model.tsv"
    trait.write_text(
        "leaf_name\tc1\tc2\tc3\n"
        "A\t0\t0\t1\nB\t0\t1\t1\nC\t0\t0\t0\n"
        "D\t1\t1\t0\nE\t1\t0\t1\nF\t1\t1\t0\n"
    )
    main(
        [
            "asr",
            "-i",
            "[&R](((A:1,B:1):1,C:1):1,((D:1,E:1):1,F:1):1)R;",
            "--trait",
            str(trait),
            "--state-column",
            "c1,c2,c3",
            "--trait-type",
            "discrete",
            "--model",
            "MK-MIXTURE",
            "--rate",
            "0.4",
            "--rate-categories",
            "2",
            "--model-out",
            str(model),
            "-o",
            str(output),
        ]
    )
    result = pd.read_csv(output, sep="\t")
    metadata = pd.read_csv(model, sep="\t").iloc[0]
    assert set(result["trait"]) == {"c1", "c2", "c3"}
    assert metadata["model"] == "MK-MIXTURE"
    assert metadata["rate_mixture"] == "gamma"
    assert metadata["rate_categories"] == 2


def test_mk_mixture_rejects_unidentifiable_zero_fixed_base_rate():
    tree = tree_from("(((A:1,B:1):1,C:1):1,((D:1,E:1):1,F:1):1)R;")
    states, observed, likelihoods = character_data(
        {"first": "000111", "second": "010101"}
    )
    with pytest.raises(ValueError, match="zero base rate makes the mixture"):
        compute_mk_mixture_marginals(
            tree,
            states,
            observed,
            likelihoods,
            base_model="ER",
            categories=2,
            rate=0.0,
        )


def test_mk_mixture_averages_only_categories_with_posterior_mass():
    tree = tree_from("(A:1,B:1)R;")
    states = ["0", "1"]
    observed = {"first": {"A": "0", "B": "1"}, "second": {"A": "1", "B": "0"}}
    likelihoods = {
        character: {name: np.eye(2)[int(value)] for name, value in values.items()}
        for character, values in observed.items()
    }
    posterior, fit = compute_mk_mixture_marginals(
        tree, states, observed, likelihoods, rate=np.nextafter(0.0, 1.0)
    )
    assert np.isfinite(fit["log_likelihood"])
    for character, values in posterior.items():
        probabilities = fit["category_posterior_by_character"][character]
        assert np.any(probabilities == 0.0)
        assert probabilities.sum() == pytest.approx(1.0)
        np.testing.assert_allclose(values[tree], [0.5, 0.5])
        np.testing.assert_allclose(values[tree["A"]], likelihoods[character]["A"])
        np.testing.assert_allclose(values[tree["B"]], likelihoods[character]["B"])


def test_batched_mk_mixture_matches_independent_character_likelihoods(monkeypatch):
    tree = tree_from("((A:0.7,B:0.2,C:0.4):0.9,D:0.3,E:0.1,F:0.8)R;")
    states, observed, likelihoods = character_data(
        {"first": "000111", "second": "010101", "third": "001101"}
    )
    observed["second"]["A"] = None
    likelihoods["second"]["A"] = np.ones(2)
    monkeypatch.setattr(asr, "_MAX_MIXTURE_EDGE_VALUES", 8)
    _, fit = compute_mk_mixture_marginals(
        tree, states, observed, likelihoods, rate=0.6, root_prior_mode="empirical"
    )
    expected_total = 0.0
    for character, tips in likelihoods.items():
        terms = [
            np.log(weight)
            + asr._log_likelihood(
                tree, tips, fit["root_prior"], fit["rate_matrix"] * rate
            )
            for rate, weight in zip(
                fit["mixture_category_rates"], fit["mixture_weights"], strict=True
            )
        ]
        score = logsumexp(terms)
        expected_total += score
        np.testing.assert_allclose(
            fit["category_posterior_by_character"][character],
            np.exp(terms - score),
            atol=1e-13,
        )
    assert fit["log_likelihood"] == pytest.approx(expected_total, abs=1e-12)


@pytest.mark.parametrize("informative_zero_edges", [False, True])
def test_mk_mixture_rejects_data_without_an_informative_positive_branch(
    informative_zero_edges,
):
    tree = tree_from("((A:0)I:0,B:1)R;")
    observed = {
        character: {"A": "0" if informative_zero_edges else None, "B": None}
        for character in ("first", "second")
    }
    likelihoods = {
        character: {
            "A": np.array([1.0, 0.0]) if informative_zero_edges else np.ones(2),
            "B": np.ones(2),
        }
        for character in observed
    }
    with pytest.raises(ValueError, match="positive-length branch ancestral"):
        compute_mk_mixture_marginals(tree, ["0", "1"], observed, likelihoods, rate=0.4)


def test_mk_mixture_rejects_fractional_category_count():
    tree = tree_from("(((A:1,B:1):1,C:1):1,((D:1,E:1):1,F:1):1)R;")
    states, observed, likelihoods = character_data(
        {"first": "000111", "second": "010101"}
    )
    with pytest.raises(ValueError, match="must be an integer"):
        compute_mk_mixture_marginals(
            tree,
            states,
            observed,
            likelihoods,
            categories=2.5,
        )


def test_mk_mixture_cli_does_not_replace_zero_category_count(tmp_path):
    trait = tmp_path / "traits.tsv"
    trait.write_text("leaf_name\tc1\tc2\nA\t0\t1\nB\t0\t1\nC\t1\t0\nD\t1\t0\n")
    with pytest.raises(ValueError, match="between 2 and 8"):
        main(
            [
                "asr",
                "-i",
                "[&R]((A:1,B:1):1,(C:1,D:1):1)R;",
                "--trait",
                str(trait),
                "--state-column",
                "c1,c2",
                "--trait-type",
                "discrete",
                "--model",
                "MK-MIXTURE",
                "--rate-categories",
                "0",
                "-o",
                str(tmp_path / "asr.tsv"),
            ]
        )
