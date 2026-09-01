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
    parameters, source = read_regime_parameters(
        path, ("slow", "fast"), ("sigma2",)
    )
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
            probability *= asr._transition_matrix(
                by_node[child], child.dist
            )[root_state, state]
        expected += probability
    assert asr._log_likelihood(
        tree, observations, np.array([0.5, 0.5]), by_node
    ) == pytest.approx(math.log(expected))


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
        regimes, ["background", "background", "foreground", "background", "background", "foreground", "foreground"]
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
    trait.write_text(
        "leaf_name\tstate\nA\t0\nB\t0\nC\t1\nD\t1\nE\t0\nF\t1\n"
    )
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
