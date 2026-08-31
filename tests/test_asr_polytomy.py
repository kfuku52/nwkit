import itertools
import math

import numpy as np
import pandas as pd
import pytest
from scipy.linalg import expm

import nwkit.asr as asr
from nwkit.cli import main
from nwkit.util import read_tree


def tree_from(text):
    return read_tree(text, "1", True, quiet=True, rooted="yes")


def test_large_star_matches_analytic_likelihood():
    count, rate = 1024, 0.2
    tree = tree_from("(" + ",".join(f"T{i}:1" for i in range(count)) + ");")
    observed = {f"T{i}": str(i % 2) for i in range(count)}
    likelihood = {name: np.eye(2)[int(state)] for name, state in observed.items()}
    posterior, fit = asr.compute_mk_marginals(
        tree, ["0", "1"], observed, likelihood, rate=rate
    )
    same = (1 + math.exp(-2 * rate)) / 2
    different = (1 - math.exp(-2 * rate)) / 2
    expected = (count / 2) * (math.log(same) + math.log(different))
    assert fit["log_likelihood"] == pytest.approx(expected, abs=1e-9)
    np.testing.assert_allclose(posterior[tree], [0.5, 0.5], atol=1e-12)
    for leaf in tree.leaves():
        np.testing.assert_allclose(posterior[leaf], likelihood[leaf.name])


def test_outside_evidence_can_restore_a_very_unlikely_inside_state():
    count, rate = 600, 0.2
    tree = tree_from("((" + ",".join(f"T{i}:1" for i in range(count)) + ")I:0,Z:0)R;")
    observed = {leaf.name: ("1" if leaf.name == "Z" else "0") for leaf in tree.leaves()}
    likelihood = {name: np.eye(2)[int(state)] for name, state in observed.items()}
    posterior, fit = asr.compute_mk_marginals(
        tree, ["0", "1"], observed, likelihood, rate=rate
    )
    expected = math.log(0.5) + count * math.log((1 - math.exp(-2 * rate)) / 2)
    assert fit["log_likelihood"] == pytest.approx(expected, abs=1e-9)
    np.testing.assert_allclose(posterior[tree["I"]], [0, 1])
    sampled = asr._sample_node_states(tree, ["0", "1"], fit, np.random.default_rng(17))
    assert sampled[tree] == sampled[tree["I"]] == sampled[tree["Z"]] == 1
    mapping = asr._simulate_stochastic_maps(tree, ["0", "1"], fit, 1, seed=17)
    assert mapping.loc[mapping.name.isin(["I", "Z"]), "total_count"].sum() == 0


@pytest.mark.parametrize(
    "model,rates",
    [("ER", [0.2]), ("SYM", [0.1, 0.3, 0.7]), ("ARD", [0.1, 0.2, 0.3, 0.4, 0.6, 0.8])],
)
def test_polytomy_marginals_match_exhaustive_enumeration_and_zero_edge_resolution(
    monkeypatch, model, rates
):
    states = ["x", "y", "z"]
    matrix = asr._build_rate_matrix(model, states, rates)
    observed = {"A": "x", "B": "y", "C": "x|z", "D": None}
    likelihood = {
        "A": np.array([1.0, 0.0, 0.0]),
        "B": np.array([0.0, 1.0, 0.0]),
        "C": np.array([1.0, 0.0, 1.0]),
        "D": np.ones(3),
    }

    def fixed_fit(tree, model, states, likelihood_by_leaf, root_prior, **kwargs):
        return {
            "rate_matrix": matrix,
            "log_likelihood": asr._log_likelihood(
                tree, likelihood_by_leaf, root_prior, matrix
            ),
        }

    monkeypatch.setattr(asr, "_fit_rate_matrix", fixed_fit)
    tree = tree_from("(A:0.3,B:0.7,C:1.1,D:0.5)R;")
    nodes = list(tree.traverse())
    transitions = {
        node: expm(matrix * float(node.dist)) for node in nodes if not node.is_root
    }
    expected = {node: np.zeros(3) for node in nodes}
    total = 0.0
    for assignment in itertools.product(range(3), repeat=len(nodes)):
        state_by_node = dict(zip(nodes, assignment, strict=True))
        weight = 1 / 3
        for node in nodes:
            if node.is_root:
                continue
            weight *= transitions[node][state_by_node[node.up], state_by_node[node]]
            weight *= likelihood[node.name][state_by_node[node]]
        total += weight
        for node, state in state_by_node.items():
            expected[node][state] += weight
    posterior, fit = asr.compute_mk_marginals(
        tree, states, observed, likelihood, model=model
    )
    assert fit["log_likelihood"] == pytest.approx(math.log(total), abs=1e-12)
    for node in nodes:
        np.testing.assert_allclose(posterior[node], expected[node] / total, atol=1e-12)

    resolved = tree_from("(((A:0.3,B:0.7):0,C:1.1):0,D:0.5)R;")
    resolved_posterior, resolved_fit = asr.compute_mk_marginals(
        resolved, states, observed, likelihood, model=model
    )
    assert resolved_fit["log_likelihood"] == pytest.approx(
        fit["log_likelihood"], abs=1e-12
    )
    for node in nodes:
        np.testing.assert_allclose(
            posterior[node], resolved_posterior[resolved[node.name]], atol=1e-12
        )


@pytest.mark.parametrize("model", ["ER", "SYM", "ARD"])
def test_asr_cli_polytomy_with_missing_tip_and_stochastic_maps(tmp_path, model):
    trait = tmp_path / "traits.tsv"
    trait.write_text("leaf_name\tstate\nA\tx\nB\ty\nC\tx|y\nD\tNA\n")
    output, annotated, mapping = (
        tmp_path / "asr.tsv",
        tmp_path / "asr.nwk",
        tmp_path / "maps.tsv",
    )
    main(
        [
            "asr",
            "-i",
            "(A:1,B:1,C:1,D:1);",
            "--input-rooted",
            "yes",
            "--trait",
            str(trait),
            "--state-column",
            "state",
            "--model",
            model,
            "--rate",
            "0.2",
            "--rate-bounds",
            "0.001,1",
            "--tree-out",
            str(annotated),
            "--tree-annotation",
            "all",
            "--stochastic-map-out",
            str(mapping),
            "--n-sim",
            "3",
            "--seed",
            "5",
            "-o",
            str(output),
        ]
    )
    frame = pd.read_csv(output, sep="\t")
    assert len(frame) == 5
    assert frame.loc[frame.name == "D", "is_imputed"].all()
    np.testing.assert_allclose(frame[["p_x", "p_y"]].sum(axis=1), 1, atol=1e-12)
    tree = read_tree(str(annotated), "auto", True, quiet=True)
    assert len(tree.children) == 4
    assert "asr_p_x" in tree.props
    assert set(pd.read_csv(mapping, sep="\t").name) == {"A", "B", "C", "D"}


def test_asr_force_does_not_bypass_branch_length_validation(tmp_path):
    with pytest.raises(ValueError, match="non-negative"):
        main(
            [
                "asr",
                "-i",
                "(A:-1,B:1,C:1);",
                "--input-rooted",
                "yes",
                "--trait",
                str(tmp_path / "not_read.tsv"),
                "--state-column",
                "state",
            ]
        )
