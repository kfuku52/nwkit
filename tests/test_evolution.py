import numpy as np
import pandas as pd
import pytest
from ete4 import Tree

from nwkit.contrast import calculate_contrasts
from nwkit.evolution import (
    build_evolutionary_covariance,
    build_sparse_evolutionary_model,
    read_custom_covariance,
    transformed_edge_variances,
    validate_custom_covariance,
    validate_evolution_parameter,
)
from nwkit.ordinary_regression import fit_ordinary_regression
from nwkit.sparse_laplace import (
    combine_sparse_covariance_models,
    condition_sparse_tip_model,
    prepare_sparse_latent_sampler,
)

TREE_TEXT = "(((A:1,B:1):1,C:2):1,(D:1,E:1):2);"
LEAF_NAMES = ["A", "B", "C", "D", "E"]


def _tree(text=TREE_TEXT):
    return Tree(text, parser=1)


def _values(values):
    return dict(zip(LEAF_NAMES, values, strict=True))


@pytest.mark.parametrize(
    "model,parameter",
    [
        ("brownian", None),
        ("lambda", 0.6),
        ("ou", 0.4),
        ("kappa", 1.4),
        ("delta", 1.3),
        ("eb", -0.2),
        ("acdc", 0.1),
        ("independent", None),
    ],
)
def test_sparse_evolutionary_model_matches_normalized_dense_covariance(
    model, parameter
):
    dense = build_evolutionary_covariance(
        _tree(), LEAF_NAMES, model=model, parameter=parameter
    )
    sparse_model = build_sparse_evolutionary_model(
        _tree(), LEAF_NAMES, model=model, parameter=parameter
    )
    precision_factor = sparse_model.precision_factor()
    np.testing.assert_allclose(
        (precision_factor.T @ precision_factor).toarray(),
        sparse_model.precision.toarray(),
        rtol=2e-10,
        atol=2e-10,
    )
    np.testing.assert_allclose(
        sparse_model.materialize(),
        dense / np.mean(np.diag(dense)),
        rtol=2e-10,
        atol=2e-10,
    )


def test_sparse_predictor_conditioning_matches_dense_gaussian_posterior():
    dense = build_evolutionary_covariance(_tree(), LEAF_NAMES)
    sparse_prior = build_sparse_evolutionary_model(_tree(), LEAF_NAMES)
    prior_variance = 1.7 * float(np.mean(np.diag(dense)))
    sampling = np.asarray([0.2, 0.3, 0.4, 0.5, 0.6])
    conditioned = condition_sparse_tip_model(sparse_prior, prior_variance, sampling)
    prior = 1.7 * dense
    expected = prior - prior @ np.linalg.solve(prior + np.diag(sampling), prior)
    np.testing.assert_allclose(
        conditioned.materialize(), expected, rtol=2e-10, atol=2e-10
    )


def test_conditioned_sparse_model_supports_exact_reusable_sampling():
    sparse_prior = build_sparse_evolutionary_model(_tree(), LEAF_NAMES)
    conditioned = condition_sparse_tip_model(
        sparse_prior,
        1.7,
        np.asarray([0.2, 0.3, 0.4, 0.5, 0.6]),
    )
    precision_factor = conditioned.precision_factor()
    np.testing.assert_allclose(
        (precision_factor.T @ precision_factor).toarray(),
        conditioned.precision.toarray(),
        rtol=2e-12,
        atol=2e-12,
    )
    latent = combine_sparse_covariance_models({"posterior": (1.0, conditioned)})
    sampler = prepare_sparse_latent_sampler(latent)
    rng = np.random.default_rng(17)
    draws = np.asarray([sampler.sample(rng) for _ in range(30_000)])
    np.testing.assert_allclose(
        np.cov(draws, rowvar=False),
        conditioned.materialize(),
        rtol=0.06,
        atol=0.01,
    )


def test_kappa_covariance_raises_each_branch_length_to_power():
    covariance = build_evolutionary_covariance(
        _tree(), LEAF_NAMES, model="kappa", parameter=2.0
    )
    expected = np.asarray(
        [
            [3, 2, 1, 0, 0],
            [2, 3, 1, 0, 0],
            [1, 1, 5, 0, 0],
            [0, 0, 0, 5, 4],
            [0, 0, 0, 4, 5],
        ],
        dtype=float,
    )
    np.testing.assert_allclose(covariance, expected)


def test_delta_covariance_raises_ultrametric_node_depths_to_power():
    covariance = build_evolutionary_covariance(
        _tree(), LEAF_NAMES, model="delta", parameter=2.0
    )
    assert covariance[0, 0] == pytest.approx(3.0)
    assert covariance[0, 1] == pytest.approx(4.0 / 3.0)
    assert covariance[0, 2] == pytest.approx(1.0 / 3.0)
    assert covariance[3, 4] == pytest.approx(4.0 / 3.0)


@pytest.mark.parametrize("model,rate_change", [("eb", -0.4), ("acdc", 0.2)])
def test_rate_change_covariance_integrates_exponential_rate(model, rate_change):
    covariance = build_evolutionary_covariance(
        _tree(), LEAF_NAMES, model=model, parameter=rate_change
    )

    def transformed(depth):
        return np.expm1(rate_change * depth) / rate_change

    assert covariance[0, 0] == pytest.approx(transformed(3.0))
    assert covariance[0, 1] == pytest.approx(transformed(2.0))
    assert covariance[0, 2] == pytest.approx(transformed(1.0))
    assert covariance[0, 3] == pytest.approx(0.0)


@pytest.mark.parametrize(
    "model,parameter",
    [
        ("lambda", 1.0),
        ("kappa", 1.0),
        ("delta", 1.0),
        ("eb", 0.0),
        ("acdc", 0.0),
    ],
)
def test_tree_transform_identity_parameters_recover_brownian(model, parameter):
    brownian = build_evolutionary_covariance(_tree(), LEAF_NAMES)
    transformed = build_evolutionary_covariance(
        _tree(), LEAF_NAMES, model=model, parameter=parameter
    )
    np.testing.assert_allclose(transformed, brownian, rtol=1e-12, atol=1e-12)


def test_ou_tends_to_brownian_as_alpha_approaches_zero():
    brownian = build_evolutionary_covariance(_tree(), LEAF_NAMES)
    ou = build_evolutionary_covariance(_tree(), LEAF_NAMES, model="ou", parameter=1e-9)
    np.testing.assert_allclose(ou, brownian, rtol=1e-7, atol=1e-7)


def test_ultrametric_ou_contrast_tree_reproduces_tip_covariance():
    tree = _tree()
    alpha = 0.6
    edge_variances = transformed_edge_variances(
        tree,
        model="ou",
        parameter=alpha,
    )
    depths = {}
    for node in tree.traverse(strategy="preorder"):
        depths[node] = 0.0 if node.is_root else depths[node.up] + edge_variances[node]
    represented = np.zeros((len(LEAF_NAMES), len(LEAF_NAMES)))
    leaf_by_name = {str(leaf.name): leaf for leaf in tree.leaves()}
    leaves = [leaf_by_name[name] for name in LEAF_NAMES]
    for first, first_leaf in enumerate(leaves):
        for second, second_leaf in enumerate(leaves):
            represented[first, second] = depths[
                tree.common_ancestor([first_leaf, second_leaf])
            ]
    direct = build_evolutionary_covariance(
        tree,
        LEAF_NAMES,
        model="ou",
        parameter=alpha,
    )
    np.testing.assert_allclose(represented, direct)


@pytest.mark.parametrize(
    ("model", "parameter"),
    [
        ("brownian", None),
        ("lambda", 0.35),
        ("ou", 0.6),
        ("kappa", 0.7),
        ("delta", 1.4),
        ("eb", -0.2),
        ("acdc", 0.2),
        ("independent", None),
    ],
)
def test_supported_contrast_transform_whitens_its_tip_covariance(model, parameter):
    tree = _tree()
    contrasts, coefficients, leaf_names = calculate_contrasts(
        tree,
        _values([0.0] * len(LEAF_NAMES)),
        evolution_model=model,
        evolution_parameter=parameter,
        return_coefficients=True,
    )
    nodes = list(coefficients)
    transform = np.vstack([coefficients[node] for node in nodes])
    covariance = build_evolutionary_covariance(
        tree,
        leaf_names,
        model=model,
        parameter=parameter,
    )
    expected = np.diag([contrasts[node]["contrast_variance"] for node in nodes])
    np.testing.assert_allclose(
        transform @ covariance @ transform.T,
        expected,
        rtol=1e-10,
        atol=1e-10,
    )


def test_delta_and_contrast_ou_reject_non_ultrametric_trees():
    tree = _tree(TREE_TEXT.replace("C:2", "C:1"))
    with pytest.raises(ValueError, match="ultrametric"):
        build_evolutionary_covariance(tree, LEAF_NAMES, model="delta", parameter=1.2)
    with pytest.raises(ValueError, match="ultrametric"):
        transformed_edge_variances(tree, model="ou", parameter=0.5)


def test_independent_model_does_not_consume_input_branch_lengths():
    tree = _tree(TREE_TEXT.replace("A:1", "A:0"))
    covariance = build_evolutionary_covariance(
        tree,
        LEAF_NAMES,
        model="independent",
    )
    np.testing.assert_allclose(covariance, np.eye(len(LEAF_NAMES)))


@pytest.mark.parametrize(
    ("model", "parameter"),
    [
        ("brownian", 0.5),
        ("lambda", -0.1),
        ("lambda", 1.1),
        ("ou", 0.0),
        ("kappa", -0.1),
        ("delta", 0.0),
        ("eb", 0.1),
        ("acdc", float("inf")),
    ],
)
def test_evolution_parameter_domains_are_enforced(model, parameter):
    with pytest.raises(ValueError):
        validate_evolution_parameter(model, parameter)


def test_extreme_fixed_kappa_reports_transformation_overflow():
    with pytest.raises(ValueError, match="non-finite transformed branches"):
        build_evolutionary_covariance(
            _tree(),
            LEAF_NAMES,
            model="kappa",
            parameter=1e308,
        )


def test_custom_covariance_is_reordered_by_species_name_and_validated():
    order = ["E", "C", "A", "D", "B"]
    base = build_evolutionary_covariance(_tree(), LEAF_NAMES)
    frame = pd.DataFrame(base, index=LEAF_NAMES, columns=LEAF_NAMES).loc[order, order]
    observed = validate_custom_covariance(frame, LEAF_NAMES)
    np.testing.assert_allclose(observed, base)
    invalid = frame.copy()
    invalid.iloc[0, 1] += 1.0
    with pytest.raises(ValueError, match="symmetric"):
        validate_custom_covariance(invalid, LEAF_NAMES)
    with pytest.raises(ValueError, match="exactly match tree tips"):
        validate_custom_covariance(frame.drop(index="A", columns="A"), LEAF_NAMES)
    nonnumeric = frame.astype(object)
    nonnumeric.iloc[0, 0] = "not-a-number"
    with pytest.raises(ValueError, match="numeric"):
        validate_custom_covariance(nonnumeric, LEAF_NAMES)


def test_custom_covariance_tsv_rejects_non_positive_definite_matrix(tmp_path):
    path = tmp_path / "covariance.tsv"
    frame = pd.DataFrame(
        np.ones((len(LEAF_NAMES), len(LEAF_NAMES))),
        index=LEAF_NAMES,
        columns=LEAF_NAMES,
    )
    frame.index.name = "leaf_name"
    frame.to_csv(path, sep="\t")
    with pytest.raises(ValueError, match="positive definite"):
        read_custom_covariance(str(path), LEAF_NAMES)


@pytest.mark.parametrize(
    "model,parameter",
    [
        ("kappa", 0.7),
        ("delta", 1.4),
        ("eb", -0.2),
        ("acdc", 0.2),
    ],
)
def test_new_evolution_models_fit_fixed_parameter_pgls(model, parameter):
    result = fit_ordinary_regression(
        _tree(),
        {"expression": _values([2.0, 5.0, 7.5, 8.0, 12.5])},
        {"body_size": _values([1.0, 2.0, 4.0, 3.0, 7.0])},
        ["expression"],
        ["body_size"],
        evolution_model=model,
        evolution_parameter=parameter,
    )
    assert set(result["evolution_model"]) == {model}
    assert set(result["evolution_parameter_status"]) == {"fixed"}
    assert set(result["evolution_parameter"]) == {parameter}


@pytest.mark.parametrize("model", ["kappa", "delta", "eb", "acdc"])
def test_new_evolution_model_parameters_can_be_estimated(model):
    result = fit_ordinary_regression(
        _tree(),
        {"expression": _values([2.0, 5.0, 7.5, 8.0, 12.5])},
        {"body_size": _values([1.0, 2.0, 4.0, 3.0, 7.0])},
        ["expression"],
        ["body_size"],
        evolution_model=model,
    )
    assert set(result["evolution_parameter_status"]) == {"estimated"}
    assert np.isfinite(result.iloc[0]["evolution_parameter"])


def test_custom_covariance_pgls_matches_direct_gls():
    tree = _tree()
    response = np.asarray([2.0, 5.0, 7.5, 8.0, 12.5])
    predictor = np.asarray([1.0, 2.0, 4.0, 3.0, 7.0])
    covariance = build_evolutionary_covariance(tree, LEAF_NAMES)
    result = fit_ordinary_regression(
        tree,
        {"expression": _values(response)},
        {"body_size": _values(predictor)},
        ["expression"],
        ["body_size"],
        evolution_model="custom",
        custom_covariance=covariance,
    ).set_index("term")
    inverse = np.linalg.inv(covariance)
    design = np.column_stack([np.ones(len(predictor)), predictor])
    expected = np.linalg.solve(
        design.T @ inverse @ design,
        design.T @ inverse @ response,
    )
    assert result.loc["(intercept)", "coefficient"] == pytest.approx(expected[0])
    assert result.loc["body_size", "coefficient"] == pytest.approx(expected[1])
