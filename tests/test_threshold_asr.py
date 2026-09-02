import numpy as np
import pandas as pd
import pytest

from nwkit.cli import main
from nwkit.threshold_asr import (
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
