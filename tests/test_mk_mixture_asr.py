import numpy as np
import pandas as pd
import pytest

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
