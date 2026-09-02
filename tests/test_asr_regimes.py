import math

import numpy as np
import pandas as pd
import pytest

from nwkit import asr
from nwkit.asr_regimes import read_regime_map, read_regime_parameters
from nwkit.cli import main
from nwkit.util import read_tree


def tree_from(source):
    return read_tree(source, "1", True, quiet=True, rooted="yes")


def write_regime_map(path, regimes):
    rows = ["branch_id\tregime"]
    rows.extend(f"{branch_id}\t{regime}" for branch_id, regime in enumerate(regimes))
    path.write_text("\n".join(rows) + "\n")


def test_regime_map_requires_one_assignment_for_every_branch(tmp_path):
    tree = tree_from("(A:1,B:1,C:1)R;")
    path = tmp_path / "regimes.tsv"
    write_regime_map(path, ["slow", "slow", "fast"])
    with pytest.raises(ValueError, match="including root 0; missing: 3"):
        read_regime_map(path, tree)

    write_regime_map(path, ["root", "slow", "fast", "fast"])
    assignment = read_regime_map(path, tree)
    assert assignment.regimes == ("root", "slow", "fast")
    assert assignment.root_regime == "root"
    assert [assignment.by_node[node] for node in tree.traverse()] == [
        "root",
        "slow",
        "fast",
        "fast",
    ]


def test_regime_parameters_are_complete_finite_and_ordered(tmp_path):
    path = tmp_path / "parameters.tsv"
    path.write_text("regime\tsigma2\nfast\t2\nslow\t0.5\n")
    parameters, source = read_regime_parameters(path, ("slow", "fast"), ("sigma2",))
    assert list(parameters) == ["slow", "fast"]
    assert parameters["slow"]["sigma2"] == 0.5
    assert source == str(path)


def test_edge_specific_q_likelihood_matches_direct_star_calculation(tmp_path):
    tree = tree_from("(A:1,B:2,C:0.5)R;")
    path = tmp_path / "regimes.tsv"
    write_regime_map(path, ["slow", "slow", "fast", "fast"])
    assignment = read_regime_map(path, tree)
    slow = asr._build_rate_matrix("ER", ["0", "1"], [0.1])
    fast = asr._build_rate_matrix("ER", ["0", "1"], [1.2])
    by_node = {
        node: {"slow": slow, "fast": fast}[assignment.by_node[node]]
        for node in tree.traverse()
        if not node.is_root
    }
    observations = {
        "A": np.array([1.0, 0.0]),
        "B": np.array([0.0, 1.0]),
        "C": np.array([0.0, 1.0]),
    }
    observed = [0, 1, 1]
    expected = 0.0
    for root_state in range(2):
        probability = 0.5
        for child, state in zip(tree.children, observed, strict=True):
            probability *= asr._transition_matrix(by_node[child], child.dist)[
                root_state, state
            ]
        expected += probability
    assert asr._log_likelihood(
        tree, observations, np.array([0.5, 0.5]), by_node
    ) == pytest.approx(math.log(expected))


@pytest.mark.parametrize(
    "tree_source,regime_by_branch",
    [
        ("(A:1,B:1,M:1)R;", ["background", "background", "background", "missing"]),
        ("(A:1,B:1,M:0)R;", ["background", "background", "background", "zero"]),
    ],
)
def test_mk_regime_rejects_branches_without_informative_likelihood_effect(
    tree_source, regime_by_branch, tmp_path
):
    tree = tree_from(tree_source)
    path = tmp_path / "regimes.tsv"
    write_regime_map(path, regime_by_branch)
    assignment = read_regime_map(path, tree)
    states = ["0", "1"]
    observed = {"A": "0", "B": "1", "M": None}
    likelihood = {
        "A": np.array([1.0, 0.0]),
        "B": np.array([0.0, 1.0]),
        "M": np.ones(2),
    }
    if tree_source.endswith("M:0)R;"):
        observed["M"] = "0"
        likelihood["M"] = np.array([1.0, 0.0])

    with pytest.raises(ValueError, match="uninformative regime"):
        asr.compute_mk_marginals(
            tree,
            states,
            observed,
            likelihood,
            model="MK-REGIME",
            regime_assignment=assignment,
        )


def test_hrm_does_not_mutate_a_read_only_transition_graph(monkeypatch):
    tree = tree_from("(A:1,B:1)R;")
    graph = np.array([[True, True], [True, True]], dtype=bool)
    original = graph.copy()
    graph.setflags(write=False)
    captured = {}

    def fake_compute_mk_marginals(
        tree, states, observed_state_by_leaf, likelihood_by_leaf, **options
    ):
        captured["graph"] = options["transition_graph"]
        prior = np.full(len(states), 1.0 / len(states))
        posterior = {node: prior.copy() for node in tree.traverse()}
        return posterior, {"root_prior": prior, "posterior_by_node": posterior}

    monkeypatch.setattr(asr, "compute_mk_marginals", fake_compute_mk_marginals)
    asr.compute_hrm_marginals(
        tree,
        ["0", "1"],
        {"A": "0", "B": "1"},
        {"A": np.array([1.0, 0.0]), "B": np.array([0.0, 1.0])},
        transition_graph=graph,
    )

    np.testing.assert_array_equal(graph, original)
    assert not np.any(np.diag(captured["graph"]))


def test_prior_only_hrm_rejects_fitted_hidden_rates():
    tree = tree_from("(A:1,B:1)R;")
    with pytest.raises(ValueError, match="fully fixed transition process"):
        asr.compute_hrm_marginals(
            tree,
            ["0", "1"],
            {"A": None, "B": None},
            {"A": np.ones(2), "B": np.ones(2)},
        )


def test_hrm_validates_hidden_category_count_before_expansion():
    tree = tree_from("(A:1,B:1)R;")
    observed = {"A": "0", "B": "1"}
    likelihoods = {"A": np.array([1.0, 0.0]), "B": np.array([0.0, 1.0])}
    with pytest.raises(ValueError, match="must be an integer"):
        asr.compute_hrm_marginals(
            tree,
            ["0", "1"],
            observed,
            likelihoods,
            hidden_categories=2.5,
        )
    with pytest.raises(ValueError, match="more than 256 free transition rates"):
        asr.compute_hrm_marginals(
            tree,
            ["0", "1"],
            observed,
            likelihoods,
            hidden_categories=1000,
        )

    many_states = [str(index) for index in range(33)]
    sparse_graph = np.zeros((33, 33), dtype=bool)
    sparse_observed = {"A": "0", "B": "1"}
    sparse_likelihoods = {
        "A": np.eye(33)[0],
        "B": np.eye(33)[1],
    }
    with pytest.raises(ValueError, match="more than 64 expanded states"):
        asr.compute_hrm_marginals(
            tree,
            many_states,
            sparse_observed,
            sparse_likelihoods,
            hidden_categories=2,
            transition_graph=sparse_graph,
        )


@pytest.mark.integration
def test_mk_regime_cli_fits_and_reports_each_regime_q(tmp_path):
    tree = tmp_path / "tree.nwk"
    traits = tmp_path / "traits.tsv"
    regimes = tmp_path / "regimes.tsv"
    output = tmp_path / "asr.tsv"
    model = tmp_path / "model.tsv"
    stochastic = tmp_path / "maps.tsv"
    tree.write_text("[&R]((A:1,B:1):1,(C:1,D:1):1)R;\n")
    traits.write_text("leaf_name\tstate\nA\t0\nB\t0\nC\t1\nD\t1\n")
    # level order: root, two internals, then four leaves
    write_regime_map(
        regimes,
        [
            "background",
            "background",
            "foreground",
            "background",
            "background",
            "foreground",
            "foreground",
        ],
    )
    main(
        [
            "asr",
            "-i",
            str(tree),
            "--trait",
            str(traits),
            "--state-column",
            "state",
            "--trait-type",
            "discrete",
            "--model",
            "MK-REGIME",
            "--regime-map",
            str(regimes),
            "--rate-bounds",
            "0.001,10",
            "--model-out",
            str(model),
            "--stochastic-map-out",
            str(stochastic),
            "--n-sim",
            "8",
            "--threads",
            "2",
            "--seed",
            "17",
            "-o",
            str(output),
        ]
    )
    result = pd.read_csv(output, sep="\t")
    metadata = pd.read_csv(model, sep="\t").iloc[0]
    maps = pd.read_csv(stochastic, sep="\t")
    assert np.allclose(result[["p_0", "p_1"]].sum(axis=1), 1.0)
    assert metadata["model"] == "MK-REGIME"
    assert metadata["regime_model"] == "ER"
    assert metadata["regimes"] == "background,foreground"
    assert math.isfinite(metadata["rate_background_0_to_1"])
    assert math.isfinite(metadata["rate_foreground_0_to_1"])
    assert len(maps) == 12
    assert set(maps["num_simulations"]) == {8}


@pytest.mark.integration
def test_hrm_cli_marginalizes_hidden_classes_and_maps_observed_changes(tmp_path):
    trait = tmp_path / "traits.tsv"
    output = tmp_path / "asr.tsv"
    model = tmp_path / "model.tsv"
    stochastic = tmp_path / "maps.tsv"
    trait.write_text("leaf_name\tstate\nA\t0\nB\t0\nC\t1\nD\t1\nE\t0\nF\t1\n")
    main(
        [
            "asr",
            "-i",
            "[&R]((A:1,B:1):1,(C:1,D:1):1,E:2,F:2)R;",
            "--trait",
            str(trait),
            "--state-column",
            "state",
            "--trait-type",
            "discrete",
            "--model",
            "HRM",
            "--hidden-categories",
            "2",
            "--rate-bounds",
            "0.01,3",
            "--model-out",
            str(model),
            "--stochastic-map-out",
            str(stochastic),
            "--n-sim",
            "8",
            "--seed",
            "7",
            "-o",
            str(output),
        ]
    )
    result = pd.read_csv(output, sep="\t")
    metadata = pd.read_csv(model, sep="\t").iloc[0]
    maps = pd.read_csv(stochastic, sep="\t")
    assert np.allclose(result[["p_0", "p_1"]].sum(axis=1), 1.0)
    assert metadata["model"] == "HRM"
    assert metadata["hidden_categories"] == 2
    assert metadata["num_expanded_states"] == 4
    assert set(maps["from_state"]) == {0, 1}
    assert set(maps["to_state"]) == {0, 1}
