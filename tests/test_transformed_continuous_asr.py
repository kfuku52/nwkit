from types import SimpleNamespace

import pytest

from nwkit.continuous_asr import compute_bm_marginals
from nwkit.continuous_asr_io import continuous_model_table
from nwkit.evolution import build_evolutionary_process
from nwkit.gaussian_inference import simulate_gaussian_process
from nwkit.transformed_continuous_asr import (
    compute_transformed_bm_marginals,
    default_parameter_bounds,
)
from nwkit.util import read_tree


def tree_from(source):
    return read_tree(source, "1", True, quiet=True, rooted="yes")


@pytest.mark.parametrize(
    "model,parameter",
    [("lambda", 1.0), ("kappa", 1.0), ("delta", 1.0), ("eb", 0.0), ("acdc", 0.0)],
)
def test_brownian_boundary_matches_bm(model, parameter):
    tree = tree_from("((A:1,B:1)I:1,(C:1,D:1)J:1)R;")
    values = {"A": 0.0, "B": 1.0, "C": 2.0, "D": 4.0}
    expected, expected_fit = compute_bm_marginals(tree, values, sigma2=0.8)
    posterior, fit = compute_transformed_bm_marginals(
        tree,
        values,
        model=model,
        sigma2=0.8,
        evolution_parameter=parameter,
    )
    assert fit.restricted_log_likelihood == pytest.approx(
        expected_fit.restricted_log_likelihood, abs=1e-12
    )
    for node in tree.traverse():
        assert posterior[node] == pytest.approx(expected[node])


@pytest.mark.parametrize(
    "model,parameter",
    [("lambda", 0.4), ("kappa", 0.7), ("delta", 1.4), ("eb", -0.2), ("acdc", 0.2)],
)
def test_fixed_transformed_models_return_finite_fit(model, parameter):
    tree = tree_from("((A:1,B:1)I:1,(C:1,D:1)J:1)R;")
    values = {"A": 0.0, "B": 1.0, "C": 2.0, "D": 4.0}
    posterior, fit = compute_transformed_bm_marginals(
        tree, values, model=model, evolution_parameter=parameter
    )
    assert fit.evolution_parameter == parameter
    assert fit.sigma2 > 0.0
    assert fit.restricted_log_likelihood is not None
    assert len(posterior) == len(list(tree.traverse()))


@pytest.mark.parametrize("model", ["lambda", "kappa", "delta", "eb", "acdc"])
def test_estimated_transformed_models_respect_bounds(model):
    tree = tree_from("(((A:0.5,B:0.5)I:0.5,C:1)K:1,(D:1,E:1)J:1)R;")
    values = {"A": -0.2, "B": 0.1, "C": 0.8, "D": 2.0, "E": 2.3}
    bounds = default_parameter_bounds(tree, model)
    _, fit = compute_transformed_bm_marginals(tree, values, model=model)
    assert bounds[0] <= fit.evolution_parameter <= bounds[1]
    assert fit.evolution_parameter_estimated
    assert fit.optimizer_grid_evaluations >= 25


def test_delta_requires_ultrametric_tree():
    tree = tree_from("((A:1,B:2)I:1,C:2)R;")
    values = {"A": 0.0, "B": 1.0, "C": 2.0}
    with pytest.raises(ValueError, match="requires an ultrametric tree"):
        compute_transformed_bm_marginals(
            tree,
            values,
            model="delta",
            evolution_parameter=1.2,
        )
    with pytest.raises(ValueError, match="requires an ultrametric tree"):
        compute_transformed_bm_marginals(tree, values, model="delta")


def test_eb_rejects_acceleration_but_acdc_accepts_it():
    tree = tree_from("((A:1,B:1)I:1,C:2)R;")
    values = {"A": 0.0, "B": 1.0, "C": 2.0}
    with pytest.raises(ValueError, match="use ACDC"):
        compute_transformed_bm_marginals(
            tree, values, model="eb", evolution_parameter=0.2
        )
    _, fit = compute_transformed_bm_marginals(
        tree, values, model="acdc", evolution_parameter=0.2
    )
    assert fit.evolution_parameter == 0.2


def test_kappa_zero_keeps_exact_zero_length_tip_constraint():
    tree = tree_from("((A:0,B:1)I:1,C:2)R;")
    posterior, _ = compute_transformed_bm_marginals(
        tree,
        {"A": 2.0, "B": 1.0, "C": 0.0},
        model="kappa",
        sigma2=1.0,
        evolution_parameter=0.0,
    )
    assert posterior[tree["A"]].variance == 0.0
    assert posterior[tree["I"]].mean == 2.0
    assert posterior[tree["I"]].variance == 0.0


def test_profile_interval_is_reported_in_model_table():
    tree = tree_from(
        "((((A:0.3,B:0.7):0.4,(C:0.5,D:0.9):0.6):0.8,"
        "((E:0.2,F:1.0):0.5,(G:0.4,H:0.8):0.7):0.9):0.6,"
        "(((I:0.6,J:0.3):0.7,(K:0.9,L:0.4):0.5):0.8,"
        "((M:0.7,N:0.2):0.6,(O:1.0,P:0.5):0.4):0.9):0.7)R;"
    )
    process = build_evolutionary_process(
        tree,
        model="lambda",
        parameter=0.55,
        root_mode="fixed",
        root_mean=1.0,
        variance_scale=0.8,
        allow_zero=True,
    )
    draw = simulate_gaussian_process(process, num_samples=1, seed=15)
    values = {
        str(node.name): float(draw.values[0, index])
        for index, node in enumerate(draw.nodes)
        if node.is_leaf
    }
    _, fit = compute_transformed_bm_marginals(
        tree, values, model="lambda", profile_ci_level=0.8
    )
    assert fit.evolution_parameter_ci_lower is not None
    assert fit.evolution_parameter_ci_upper is not None
    assert (
        fit.evolution_parameter_ci_lower
        < fit.evolution_parameter
        < fit.evolution_parameter_ci_upper
    )
    table = continuous_model_table(
        fit,
        SimpleNamespace(
            model="LAMBDA",
            trait_type="continuous",
            state_column="x",
            standard_error_column=None,
        ),
        0.95,
    )
    row = table.iloc[0]
    assert row["profile_ci_level"] == pytest.approx(0.8)
    assert row["evolution_parameter_ci_lower"] == pytest.approx(
        fit.evolution_parameter_ci_lower
    )
