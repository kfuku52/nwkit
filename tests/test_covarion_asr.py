import math

import numpy as np
import pandas as pd
import pytest

from nwkit.asr import _covarion_parameter_setup, compute_covarion_marginals
from nwkit.cli import main
from nwkit.hidden_rate_asr import covarion_rate_matrix
from nwkit.util import read_tree


def tree_from(source):
    return read_tree(source, "1", True, quiet=True, rooted="yes")


def likelihoods(states):
    index = {state: position for position, state in enumerate(states)}
    observed = {"A": "0", "B": "0", "C": "1", "D": "1", "E": "0", "F": "1"}
    result = {}
    for name, state in observed.items():
        values = np.zeros(len(states), dtype=float)
        values[index[state]] = 1.0
        result[name] = values
    return observed, result


def test_structured_covarion_generator_has_ordered_multipliers():
    graph = np.ones((2, 2), dtype=bool)
    np.fill_diagonal(graph, False)
    matrix, multipliers = covarion_rate_matrix(graph, 3, 0.7, 1.2, 0.3)
    assert matrix.shape == (6, 6)
    assert np.allclose(matrix.sum(axis=1), 0.0)
    assert np.all(np.diff(multipliers) > 0.0)
    assert np.prod(multipliers) == pytest.approx(1.0)


def test_structured_covarion_generator_enforces_effective_rate_bounds():
    graph = np.ones((2, 2), dtype=bool)
    np.fill_diagonal(graph, False)
    with pytest.raises(ValueError, match="within the fitted rate bounds"):
        covarion_rate_matrix(
            graph,
            2,
            10.0,
            2.0,
            0.3,
            effective_rate_bounds=(0.1, 20.0),
        )


def test_covarion_rejects_bounds_indistinguishable_on_log_scale():
    lower = 1e100
    upper = float(np.nextafter(lower, math.inf))
    assert lower < upper
    with pytest.raises(ValueError, match="distinguishable on the fitted log scale"):
        _covarion_parameter_setup((lower, upper), lower, None)


def test_covarion_fit_aggregates_hidden_class_probabilities():
    tree = tree_from("(((A:1,B:1):1,C:1):1,((D:1,E:1):1,F:1):1)R;")
    states = ["0", "1"]
    observed, tip_likelihoods = likelihoods(states)
    posterior, fit = compute_covarion_marginals(
        tree,
        states,
        observed,
        tip_likelihoods,
        hidden_categories=2,
        rate=0.4,
        rate_bounds=(1e-4, 10.0),
    )
    assert fit["model"] == "COVARION"
    assert fit["base_rate"] == 0.4
    assert len(fit["hidden_rate_multipliers"]) == 2
    assert np.min(fit["hidden_rate_effective_rates"]) >= 1e-4
    assert np.max(fit["hidden_rate_effective_rates"]) <= 10.0
    assert np.allclose([values.sum() for values in posterior.values()], 1.0)
    assert all(len(values) == 2 for values in posterior.values())


def test_covarion_fixed_rate_must_leave_room_for_hidden_spread():
    tree = tree_from("(((A:1,B:1):1,C:1):1,((D:1,E:1):1,F:1):1)R;")
    states = ["0", "1"]
    observed, tip_likelihoods = likelihoods(states)
    with pytest.raises(ValueError, match="strictly inside"):
        compute_covarion_marginals(
            tree,
            states,
            observed,
            tip_likelihoods,
            rate=10.0,
            rate_bounds=(1e-4, 10.0),
        )


def test_covarion_rejects_fractional_or_excessive_hidden_state_spaces():
    tree = tree_from("(((A:1,B:1):1,C:1):1,((D:1,E:1):1,F:1):1)R;")
    states = ["0", "1"]
    observed, tip_likelihoods = likelihoods(states)
    with pytest.raises(ValueError, match="must be an integer"):
        compute_covarion_marginals(
            tree,
            states,
            observed,
            tip_likelihoods,
            hidden_categories=2.5,
        )
    with pytest.raises(ValueError, match="more than 64 expanded states"):
        compute_covarion_marginals(
            tree,
            states,
            observed,
            tip_likelihoods,
            hidden_categories=33,
        )


def test_covarion_cli_writes_structured_parameters(tmp_path):
    trait = tmp_path / "traits.tsv"
    output = tmp_path / "asr.tsv"
    model = tmp_path / "model.tsv"
    trait.write_text("leaf_name\tstate\nA\t0\nB\t0\nC\t1\nD\t1\nE\t0\nF\t1\n")
    main(
        [
            "asr",
            "-i",
            "[&R](((A:1,B:1):1,C:1):1,((D:1,E:1):1,F:1):1)R;",
            "--trait",
            str(trait),
            "--state-column",
            "state",
            "--model",
            "COVARION",
            "--trait-type",
            "discrete",
            "--rate",
            "0.4",
            "--hidden-categories",
            "2",
            "--model-out",
            str(model),
            "-o",
            str(output),
        ]
    )
    result = pd.read_csv(output, sep="\t")
    metadata = pd.read_csv(model, sep="\t").iloc[0]
    assert {"p_0", "p_1"}.issubset(result.columns)
    assert metadata["model"] == "COVARION"
    assert metadata["base_rate"] == pytest.approx(0.4)
    assert len(str(metadata["hidden_rate_multipliers"]).split(",")) == 2
    effective = [
        float(value)
        for value in str(metadata["hidden_rate_effective_rates"]).split(",")
    ]
    assert min(effective) >= 1e-9
    assert max(effective) <= 1e3
