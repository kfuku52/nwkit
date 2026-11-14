import json

import numpy as np
import pandas as pd
import pytest
from scipy.linalg import expm

from nwkit.asr import _get_root_prior, compute_mk_marginals
from nwkit.asr_discrete_cross_validation import discrete_cross_validation
from nwkit.cli import main
from nwkit.util import read_tree


@pytest.mark.parametrize("mode", ["tip", "clade"])
def test_ctmc_holdout_matches_independent_two_tip_conditional(mode):
    tree = read_tree("(A:1,B:2)R;", "1", True, quiet=True, rooted="yes")
    states = ["x", "y"]
    observed = {"A": "x", "B": "y"}
    likelihoods = {"A": np.array([1.0, 0]), "B": np.array([0.0, 1])}
    q = np.array([[-0.3, 0.3], [0.3, -0.3]])
    calls = []

    def refit(training, training_likelihoods):
        calls.append(training)
        held = next(name for name in training if training[name] is None)
        assert training_likelihoods[held] == pytest.approx([1, 1])
        return compute_mk_marginals(
            tree,
            states,
            training,
            training_likelihoods,
            model="CUSTOM",
            fixed_rate_matrix=q,
        )[1]

    result = discrete_cross_validation(
        tree, states, observed, likelihoods, refit, mode=mode
    )
    probability = expm(q * 3)[0, 1]
    assert result.observation_probability.to_numpy() == pytest.approx([probability] * 2)
    assert result.log_score.to_numpy() == pytest.approx([np.log(probability)] * 2)
    assert result.brier_score.to_numpy() == pytest.approx(
        [2 * (1 - probability) ** 2] * 2
    )
    assert len(calls) == 2
    assert observed == {"A": "x", "B": "y"}


def test_soft_likelihood_scores_observation_not_fabricated_true_state():
    tree = read_tree("(A:1,B:1)R;", "1", True, quiet=True, rooted="yes")
    likelihoods = {"A": np.array([0.8, 0.1]), "B": np.array([0, 1])}

    def refit(*unused):
        return {
            "posterior_by_node": {
                node: np.array([0.25, 0.75]) for node in tree.traverse()
            }
        }

    result = discrete_cross_validation(
        tree, ["x", "y"], {"A": None, "B": "y"}, likelihoods, refit
    ).set_index("name")
    assert result.loc["A", "observation_probability"] == pytest.approx(0.275)
    assert pd.isna(result.loc["A", "brier_score"])
    assert _get_root_prior(
        "empirical", ["x", "y"], {"A": None}, {"A": likelihoods["A"]}
    ) == pytest.approx([8 / 9, 1 / 9])


def test_discrete_cv_cli_refits_without_file_backed_observation_leakage(tmp_path):
    traits = tmp_path / "traits.tsv"
    traits.write_text("leaf_name\tstate\nA\tx\nB\ty\n")
    observations = tmp_path / "observations.tsv"
    observations.write_text("leaf_name\tx\ty\nA\t0.9\t0.1\nB\t0.2\t0.8\n")
    cv = tmp_path / "cv.tsv"
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
            "ER",
            "--states",
            "x,y",
            "--rate",
            "0.3",
            "--tip-likelihoods",
            str(observations),
            "--cross-validation-out",
            str(cv),
            "-o",
            str(tmp_path / "out.tsv"),
        ]
    )
    result = pd.read_csv(cv, sep="\t").set_index("name")
    transition = expm(np.array([[-0.3, 0.3], [0.3, -0.3]]) * 2)
    predicted_a = transition @ np.array([0.2, 0.8])
    assert result.loc["A", "observation_probability"] == pytest.approx(
        predicted_a @ np.array([0.9, 0.1])
    )
    assert list(
        json.loads(result.loc["A", "state_probabilities"]).values()
    ) == pytest.approx(predicted_a)
    assert result.brier_score.isna().all()
