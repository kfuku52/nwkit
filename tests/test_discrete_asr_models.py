import math

import numpy as np
import pandas as pd
import pytest

import nwkit.asr as asr
from nwkit.asr_models import default_model
from nwkit.cli import main
from nwkit.discrete_asr_models import (
    build_rate_matrix,
    parameter_labels,
    read_rate_matrix,
    read_transition_graph,
    stationary_distribution,
    validate_rate_matrix,
)
from nwkit.util import read_tree


def tree_from(source):
    return read_tree(source, "1", True, quiet=True, rooted="yes")


def test_f81_uses_target_specific_rates_and_stationary_frequencies():
    states = ["a", "b", "c"]
    matrix = build_rate_matrix("F81", states, [1.0, 2.0, 3.0])
    np.testing.assert_allclose(
        matrix,
        [[-5.0, 2.0, 3.0], [1.0, -4.0, 3.0], [1.0, 2.0, -3.0]],
    )
    np.testing.assert_allclose(stationary_distribution(matrix), [1 / 6, 2 / 6, 3 / 6])


def test_single_state_f81_has_no_unidentifiable_rate_parameter():
    assert parameter_labels("F81", ["only"]) == []
    np.testing.assert_array_equal(build_rate_matrix("F81", ["only"], []), [[0.0]])


def test_gtr_satisfies_detailed_balance_with_fitted_frequency_ratios():
    states = ["a", "b", "c"]
    # Three pair exchangeabilities followed by pi_b/pi_a and pi_c/pi_a.
    matrix = build_rate_matrix("GTR", states, [2.0, 3.0, 5.0, 4.0, 2.0])
    equilibrium = stationary_distribution(matrix)
    np.testing.assert_allclose(equilibrium, [1 / 7, 4 / 7, 2 / 7])
    flux = equilibrium[:, None] * matrix
    np.testing.assert_allclose(flux, flux.T, atol=1e-15)


@pytest.mark.parametrize("model", ["F81", "GTR"])
def test_frequency_models_require_complete_transition_graph(model):
    graph, _ = read_transition_graph("ordered", ["a", "b", "c"])
    with pytest.raises(ValueError, match="requires a complete"):
        build_rate_matrix(model, ["a", "b", "c"], [1.0], graph)


@pytest.mark.integration
@pytest.mark.parametrize("model", ["F81", "GTR"])
def test_frequency_model_cli_fits_stationary_root_and_reports_equilibrium(
    model, tmp_path
):
    trait = tmp_path / "traits.tsv"
    output = tmp_path / "asr.tsv"
    metadata = tmp_path / "model.tsv"
    trait.write_text(
        "leaf_name\tstate\nA\tx\nB\ty\nC\tz\nD\tx\nE\ty\nF\tz\nG\tx\nH\tz\n"
    )
    main(
        [
            "asr",
            "-i",
            "[&R]((A:0.2,B:0.4):0.3,(C:0.7,D:1.1):0.5,E:0.6,F:1.3,G:1.7,H:2)R;",
            "--trait",
            str(trait),
            "--state-column",
            "state",
            "--model",
            model,
            "--rate-bounds",
            "0.001,20",
            "--model-out",
            str(metadata),
            "-o",
            str(output),
        ]
    )
    result = pd.read_csv(output, sep="\t")
    fit = pd.read_csv(metadata, sep="\t").iloc[0]
    assert np.allclose(result[["p_x", "p_y", "p_z"]].sum(axis=1), 1.0)
    assert fit["model"] == model
    assert fit["root_prior"] == "stationary"
    equilibrium = np.asarray([fit[f"equilibrium_{state}"] for state in "xyz"])
    assert np.all(equilibrium > 0.0)
    assert equilibrium.sum() == pytest.approx(1.0)
    assert math.isfinite(fit["log_likelihood"])


def test_ordered_graph_uses_state_order_and_removes_nonadjacent_rates():
    states = ["low", "medium", "high"]
    graph, source = read_transition_graph("ordered", states)
    matrix = build_rate_matrix("ER", states, [0.4], graph)
    assert source == "ordered"
    np.testing.assert_array_equal(
        graph,
        [[False, True, False], [True, False, True], [False, True, False]],
    )
    np.testing.assert_allclose(
        matrix,
        [[-0.4, 0.4, 0.0], [0.4, -0.8, 0.4], [0.0, 0.4, -0.4]],
    )


def test_directed_graph_supports_irreversible_ard_and_rejects_sym(tmp_path):
    path = tmp_path / "graph.tsv"
    path.write_text("from_state\tto_state\n0\t1\n1\t2\n")
    graph, _ = read_transition_graph(path, ["0", "1", "2"])
    matrix = build_rate_matrix("ARD", ["0", "1", "2"], [0.2, 0.7], graph)
    np.testing.assert_allclose(
        matrix, [[-0.2, 0.2, 0.0], [0.0, -0.7, 0.7], [0.0, 0.0, 0.0]]
    )
    with pytest.raises(ValueError, match="SYM requires a symmetric"):
        build_rate_matrix("SYM", ["0", "1", "2"], [], graph)


def test_stationary_distribution_handles_asymmetry_and_unique_absorbing_class():
    asymmetric = np.array([[-0.2, 0.2], [0.4, -0.4]])
    np.testing.assert_allclose(stationary_distribution(asymmetric), [2 / 3, 1 / 3])
    irreversible = np.array([[-0.2, 0.2], [0.0, 0.0]])
    np.testing.assert_allclose(stationary_distribution(irreversible), [0.0, 1.0])
    with pytest.raises(ValueError, match="unique stationary"):
        stationary_distribution(np.zeros((2, 2)))


def test_stationary_prior_rejects_structurally_impossible_tip_states_early():
    tree = tree_from("[&R]((a:0.3,b:0.7):0.5,c:0.2)R;")
    states = ["x", "y", "z"]
    observed = {"a": "x", "b": "y", "c": "z"}
    likelihood = {
        name: np.array([state == item for item in states], dtype=float)
        for name, state in observed.items()
    }
    graph = np.zeros((3, 3), dtype=bool)
    graph[0, 1] = graph[1, 2] = True
    with pytest.raises(ValueError, match="zero likelihood"):
        asr.compute_mk_marginals(
            tree,
            states,
            observed,
            likelihood,
            model="ARD",
            transition_graph=graph,
            root_prior_mode="stationary",
        )


def test_labelled_rate_matrix_derives_or_validates_diagonal(tmp_path):
    derived = tmp_path / "derived.tsv"
    derived.write_text("state\tx\ty\nx\t0\t0.2\ny\t0.4\t0\n")
    states, matrix = read_rate_matrix(derived)
    assert states == ["x", "y"]
    np.testing.assert_allclose(matrix, [[-0.2, 0.2], [0.4, -0.4]])

    invalid = tmp_path / "invalid.tsv"
    invalid.write_text("state\tx\ty\nx\t-1\t0.2\ny\t0.4\t-0.4\n")
    with pytest.raises(ValueError, match="diagonals"):
        read_rate_matrix(invalid)
    with pytest.raises(ValueError, match="rows must sum to zero"):
        validate_rate_matrix([[-0.1, 0.2], [0.4, -0.4]])


def test_generator_validation_is_relative_to_the_rate_scale():
    invalid = np.array([[-2e-15, 1e-15], [3e-15, -3e-15]])
    with pytest.raises(ValueError, match="rows must sum to zero"):
        validate_rate_matrix(invalid)
    valid = np.array([[-1e-15, 1e-15], [3e-15, -3e-15]])
    np.testing.assert_array_equal(validate_rate_matrix(valid), valid)
    np.testing.assert_allclose(stationary_distribution(valid), [0.75, 0.25])


def test_transition_matrix_rejects_material_exponential_repairs(monkeypatch):
    monkeypatch.setattr(
        asr,
        "expm",
        lambda matrix: np.array([[1.001, -0.001], [0.2, 0.8]]),
    )
    with pytest.raises(ValueError, match="materially negative"):
        asr._transition_matrix(np.array([[-1.0, 1.0], [0.5, -0.5]]), 1.0)


def test_transition_matrix_accepts_roundoff_on_valid_long_branches():
    matrix = np.array([[-1.0, 0.6, 0.4], [0.2, -0.5, 0.3], [0.9, 0.1, -1.0]])
    transition = asr._transition_matrix(matrix, 1e6)
    np.testing.assert_allclose(transition.sum(axis=1), 1.0, atol=1e-15)
    np.testing.assert_allclose(
        transition,
        np.tile(stationary_distribution(matrix), (3, 1)),
        atol=2e-10,
    )


def test_transition_matrix_preserves_slow_changes_on_extreme_long_branches():
    rate = 1e-20
    matrix = np.array(
        [[-rate, rate, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, -1.0]]
    )
    transition = asr._transition_matrix(matrix, 1e20)
    np.testing.assert_allclose(
        transition[0], [math.exp(-1.0), -math.expm1(-1.0), 0.0], rtol=2e-13
    )
    np.testing.assert_allclose(transition[1], [0.0, 1.0, 0.0])
    np.testing.assert_allclose(transition[2], [0.0, 1.0, 0.0], atol=1e-15)
    np.testing.assert_allclose(transition.sum(axis=1), 1.0, atol=1e-15)


@pytest.mark.parametrize("model", ["ER", "F81", "GTR"])
def test_prior_only_asr_rejects_fitted_transition_processes(model):
    tree = tree_from("(A:1,B:2,C:3)R;")
    states = ["0", "1"]
    observed = dict.fromkeys(["A", "B", "C"])
    likelihood = {name: np.ones(2) for name in observed}
    with pytest.raises(ValueError, match="fully fixed transition process"):
        asr.compute_mk_marginals(
            tree,
            states,
            observed,
            likelihood,
            model=model,
        )


def test_prior_only_asr_accepts_fixed_er_and_custom_generators():
    tree = tree_from("(A:1,B:2,C:3)R;")
    states = ["0", "1"]
    observed = dict.fromkeys(["A", "B", "C"])
    likelihood = {name: np.ones(2) for name in observed}
    fixed = np.array([[-0.2, 0.2], [0.4, -0.4]])
    for options in (
        {"model": "ER", "rate": 0.2},
        {"model": "CUSTOM", "fixed_rate_matrix": fixed},
    ):
        posterior, fit = asr.compute_mk_marginals(
            tree, states, observed, likelihood, **options
        )
        assert not fit["rate_estimated"]
        for probabilities in posterior.values():
            assert probabilities.sum() == pytest.approx(1.0)


def test_custom_q_and_stationary_root_match_direct_likelihood():
    tree = tree_from("(A:0.3,B:0.7,C:1.1)R;")
    states = ["x", "y"]
    observed = {"A": "x", "B": "y", "C": "x"}
    likelihood = {
        name: np.array([state == "x", state == "y"], dtype=float)
        for name, state in observed.items()
    }
    matrix = np.array([[-0.2, 0.2], [0.4, -0.4]])
    posterior, fit = asr.compute_mk_marginals(
        tree,
        states,
        observed,
        likelihood,
        model="CUSTOM",
        fixed_rate_matrix=matrix,
        root_prior_mode="stationary",
    )
    expected_prior = np.array([2 / 3, 1 / 3])
    np.testing.assert_allclose(fit["root_prior"], expected_prior)
    assert fit["log_likelihood"] == pytest.approx(
        asr._log_likelihood(tree, likelihood, expected_prior, matrix), abs=1e-13
    )
    assert math.isclose(float(posterior[tree].sum()), 1.0)


def test_custom_q_cli_infers_state_order_and_supports_stochastic_maps(tmp_path):
    trait = tmp_path / "trait.tsv"
    q_path = tmp_path / "q.tsv"
    output = tmp_path / "asr.tsv"
    model = tmp_path / "model.tsv"
    mapping = tmp_path / "mapping.tsv"
    trait.write_text("leaf_name\tstate\nA\tx\nB\ty\nC\tx\n")
    q_path.write_text("state\tx\ty\nx\t0\t0.2\ny\t0.4\t0\n")
    main(
        [
            "asr",
            "-i",
            "[&R](A:1,B:1,C:1)R;",
            "--trait",
            str(trait),
            "--state-column",
            "state",
            "--model",
            "CUSTOM",
            "--rate-matrix",
            str(q_path),
            "--root-prior",
            "stationary",
            "--model-out",
            str(model),
            "--stochastic-map-out",
            str(mapping),
            "--n-sim",
            "3",
            "--seed",
            "4",
            "-o",
            str(output),
        ]
    )
    metadata = pd.read_csv(model, sep="\t").iloc[0]
    assert metadata.model == "CUSTOM"
    assert metadata.states == "x,y"
    np.testing.assert_allclose(
        np.fromstring(metadata.root_prior_values, sep=","), [2 / 3, 1 / 3]
    )
    assert metadata.q_source == f"fixed:{q_path}"
    assert not bool(metadata.rate_estimated)
    assert len(pd.read_csv(mapping, sep="\t")) == 6


def test_custom_model_requires_matrix_and_rejects_conflicting_states(tmp_path):
    trait = tmp_path / "trait.tsv"
    q_path = tmp_path / "q.tsv"
    trait.write_text("leaf_name\tstate\nA\tx\nB\ty\n")
    q_path.write_text("state\tx\ty\nx\t0\t0.2\ny\t0.4\t0\n")
    common = [
        "asr",
        "-i",
        "[&R](A:1,B:1)R;",
        "--trait",
        str(trait),
        "--state-column",
        "state",
        "--model",
        "CUSTOM",
    ]
    with pytest.raises(ValueError, match="requires --rate-matrix"):
        main(common)
    with pytest.raises(ValueError, match="--states must exactly match"):
        main(common + ["--rate-matrix", str(q_path), "--states", "y,x"])
    with pytest.raises(ValueError, match="--rate-bounds cannot be combined"):
        main(common + ["--rate-matrix", str(q_path), "--rate-bounds", "0.01,1"])

    trait.write_text("leaf_name\tstate\nA\tx\nB\tz\n")
    with pytest.raises(ValueError, match="defined by '--rate-matrix'.*z"):
        main(common + ["--rate-matrix", str(q_path)])


def test_ordered_graph_requires_explicit_state_order(tmp_path):
    trait = tmp_path / "trait.tsv"
    trait.write_text("leaf_name\tstate\nA\tlow\nB\thigh\n")
    arguments = [
        "asr",
        "-i",
        "[&R](A:1,B:1)R;",
        "--trait",
        str(trait),
        "--state-column",
        "state",
        "--transition-graph",
        "ordered",
    ]
    with pytest.raises(ValueError, match="requires explicit --states"):
        main(arguments)
    main(
        arguments
        + [
            "--states",
            "low,high",
            "--rate",
            "0.2",
            "-o",
            str(tmp_path / "ordered.tsv"),
        ]
    )


def test_directed_ard_multistart_avoids_single_start_local_optimum():
    tree = tree_from("[&R](L0:0.2,L1:0.5,L2:0.8,L3:1.1)R;")
    states = ["0", "1", "2"]
    graph = np.array([[False, True, False], [False, False, True], [True, False, False]])
    observed = {"L0": "1", "L1": "2", "L2": "1", "L3": "0"}
    likelihood = {
        name: np.array([state == item for item in states], dtype=float)
        for name, state in observed.items()
    }
    _, fit = asr.compute_mk_marginals(
        tree,
        states,
        observed,
        likelihood,
        model="ARD",
        transition_graph=graph,
        rate_bounds=(1e-9, 1e3),
    )
    np.testing.assert_allclose(fit["rates"], [18.2267, 9.09335, 18.0971], rtol=2e-4)
    assert fit["log_likelihood"] == pytest.approx(-4.1571426014, abs=2e-8)
    assert fit["optimizer_starts"] > 1
    assert fit["optimizer_converged_starts"] > 1
    assert fit["fit_status"] == "ok"


def test_mk_multistart_handles_penalty_regions_on_long_branches():
    tree = tree_from("[&R]((a:0.3,b:0.7):0.5,(c:50,d:0.9):0.25)R;")
    states = ["x", "y", "z"]
    observed = {"a": "x", "b": "y", "c": "z", "d": "x"}
    likelihood = {
        name: np.array([state == item for item in states], dtype=float)
        for name, state in observed.items()
    }
    _, fit = asr.compute_mk_marginals(
        tree,
        states,
        observed,
        likelihood,
        model="SYM",
    )
    assert fit["log_likelihood"] == pytest.approx(-3.666938522, abs=2e-7)
    assert fit["optimizer_converged_starts"] == fit["optimizer_starts"]
    assert fit["fit_status"] == "ok"


def test_default_models_do_not_depend_on_registry_order():
    assert default_model("discrete") == "ER"
    assert default_model("continuous") == "BM"


def test_transition_graph_error_names_inferred_trait_state_space(tmp_path):
    graph_path = tmp_path / "graph.tsv"
    graph_path.write_text("from_state\tto_state\nx\tz\n")
    with pytest.raises(ValueError, match="defined by '--trait'.*z"):
        read_transition_graph(graph_path, ["x", "y"], state_source="--trait")
