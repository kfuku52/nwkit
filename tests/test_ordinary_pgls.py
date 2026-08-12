import hashlib
import io
import json
import sys

import numpy as np
import pandas as pd
import pytest
from ete4 import Tree

from nwkit.cli import main
from nwkit.ordinary_pgls import (
    _global_bounded_scalar_minimize,
    build_phylogenetic_covariance,
    estimate_marginal_evolution_parameter,
    fit_ordinary_model_comparison,
    fit_ordinary_pgls,
)

TREE_TEXT = "(((A:1,B:1):1,C:2):1,(D:1,E:1):2);"
LEAF_NAMES = ["A", "B", "C", "D", "E"]


def _tree():
    return Tree(TREE_TEXT, parser=1)


def test_global_scalar_search_checks_multiple_likelihood_basins():
    def objective(value):
        return min((value + 1.5) ** 2, 0.1 + (value - 1.5) ** 2)

    result = _global_bounded_scalar_minimize(objective, (-3.0, 3.0))

    assert result.success
    assert result.x == pytest.approx(-1.5, abs=1e-5)
    assert result.fun == pytest.approx(0.0, abs=1e-10)


def test_marginal_evolution_parameter_estimation_uses_trait_only_ml():
    diagnostics = estimate_marginal_evolution_parameter(
        _tree(),
        _values([1.0, 1.2, 2.0, 4.0, 4.2]),
        "body_size",
        evolution_model="lambda",
    )

    assert 0.0 <= diagnostics["parameter"] <= 1.0
    assert diagnostics["parameter_status"] == "estimated"
    assert np.isfinite(diagnostics["log_likelihood"])
    assert isinstance(diagnostics["optimizer_converged"], bool)


def _values(values):
    return dict(zip(LEAF_NAMES, values, strict=True))


def _write_inputs(tmp_path, *, biological_replicates=False):
    tree_path = tmp_path / "species.nwk"
    data_path = tmp_path / "data.tsv"
    tree_path.write_text(TREE_TEXT)
    responses = [2.0, 5.0, 7.5, 8.0, 12.5]
    predictors = [1.0, 2.0, 4.0, 3.0, 7.0]
    if biological_replicates:
        rows = []
        for leaf_name, response, predictor in zip(
            LEAF_NAMES, responses, predictors, strict=True
        ):
            rows.extend(
                [
                    {
                        "leaf_name": leaf_name,
                        "sample_id": "{}_1".format(leaf_name),
                        "expression": response - 0.5,
                        "body_size": predictor,
                    },
                    {
                        "leaf_name": leaf_name,
                        "sample_id": "{}_2".format(leaf_name),
                        "expression": response + 0.5,
                        "body_size": predictor,
                    },
                ]
            )
    else:
        rows = [
            {
                "leaf_name": leaf_name,
                "expression": response,
                "body_size": predictor,
            }
            for leaf_name, response, predictor in zip(
                LEAF_NAMES, responses, predictors, strict=True
            )
        ]
    pd.DataFrame(rows).to_csv(data_path, sep="\t", index=False)
    return tree_path, data_path


def test_brownian_tip_covariance_matches_shared_root_path_lengths():
    covariance = build_phylogenetic_covariance(_tree(), LEAF_NAMES)
    expected = np.asarray(
        [
            [3, 2, 1, 0, 0],
            [2, 3, 1, 0, 0],
            [1, 1, 3, 0, 0],
            [0, 0, 0, 3, 2],
            [0, 0, 0, 2, 3],
        ],
        dtype=float,
    )
    np.testing.assert_allclose(covariance, expected)


def test_lambda_ou_and_independent_covariance_transformations():
    brownian = build_phylogenetic_covariance(_tree(), LEAF_NAMES)
    lambda_covariance = build_phylogenetic_covariance(
        _tree(), LEAF_NAMES, evolution_model="lambda", parameter=0.25
    )
    np.testing.assert_allclose(np.diag(lambda_covariance), np.diag(brownian))
    off_diagonal = ~np.eye(len(LEAF_NAMES), dtype=bool)
    np.testing.assert_allclose(
        lambda_covariance[off_diagonal],
        0.25 * brownian[off_diagonal],
    )
    np.testing.assert_allclose(
        build_phylogenetic_covariance(
            _tree(), LEAF_NAMES, evolution_model="independent"
        ),
        np.eye(len(LEAF_NAMES)),
    )
    alpha = 0.5
    ou = build_phylogenetic_covariance(
        _tree(), LEAF_NAMES, evolution_model="ou", parameter=alpha
    )
    expected_diagonal = -np.expm1(-2.0 * alpha * 3.0) / (2.0 * alpha)
    assert ou[0, 0] == pytest.approx(expected_diagonal)
    expected_ab = -np.expm1(-2.0 * alpha * 2.0) / (2.0 * alpha) * np.exp(-2.0 * alpha)
    assert ou[0, 1] == pytest.approx(expected_ab)
    assert ou[0, 3] == pytest.approx(0.0)


def test_ordinary_brownian_pgls_matches_direct_gls_with_intercept():
    tree = _tree()
    response = np.asarray([2.0, 5.0, 7.5, 8.0, 12.5])
    predictor = np.asarray([1.0, 2.0, 4.0, 3.0, 7.0])
    result = fit_ordinary_pgls(
        tree,
        {"expression": _values(response)},
        {"body_size": _values(predictor)},
        ["expression"],
        ["body_size"],
    ).set_index("term")
    covariance = build_phylogenetic_covariance(tree, LEAF_NAMES)
    inverse = np.linalg.inv(covariance)
    design = np.column_stack([np.ones(len(predictor)), predictor])
    expected = np.linalg.solve(
        design.T @ inverse @ design,
        design.T @ inverse @ response,
    )
    assert result.loc["(intercept)", "coefficient"] == pytest.approx(expected[0])
    assert result.loc["body_size", "coefficient"] == pytest.approx(expected[1])
    assert set(result["intercept"]) == {"yes"}
    assert set(result["evolution_model"]) == {"brownian"}


def test_ordinary_pgls_encodes_categorical_predictor_and_reports_omnibus_test():
    habitats = dict(
        zip(
            LEAF_NAMES,
            ["aquatic", "terrestrial", "arboreal", "terrestrial", "arboreal"],
            strict=True,
        )
    )
    result = fit_ordinary_pgls(
        _tree(),
        {"expression": _values([2.0, 4.0, 8.0, 5.0, 9.0])},
        {"habitat": habitats},
        ["expression"],
        ["habitat"],
        categorical_predictors=["habitat"],
        factor_references={"habitat": "aquatic"},
    )

    coefficient_rows = result[result["term_test"] == "coefficient"].set_index("term")
    assert set(coefficient_rows.index) == {
        "(intercept)",
        "habitat[arboreal]",
        "habitat[terrestrial]",
    }
    assert coefficient_rows.loc["habitat[arboreal]", "source_term"] == "habitat"
    assert coefficient_rows.loc["habitat[arboreal]", "predictor_type"] == "categorical"
    omnibus = result[result["term_test"] == "omnibus"].iloc[0]
    assert omnibus["term"] == "habitat"
    assert omnibus["degrees_of_freedom"] == 2
    assert 0.0 <= omnibus["p_value"] <= 1.0


def test_ordinary_pgls_auto_detects_string_predictor_as_categorical():
    result = fit_ordinary_pgls(
        _tree(),
        {"expression": _values([2.0, 4.0, 8.0, 5.0, 9.0])},
        {
            "habitat": dict(
                zip(
                    LEAF_NAMES,
                    ["water", "land", "water", "land", "water"],
                    strict=True,
                )
            )
        },
        ["expression"],
        ["habitat"],
    )

    factor = result[result["term"] == "habitat[water]"].iloc[0]
    assert factor["predictor_reference"] == "land"
    assert factor["factor_coding"] == "treatment"


def test_predictor_measurement_error_reduces_attenuation_in_simulation():
    names = ["S{}".format(index) for index in range(32)]
    nodes = ["{}:1".format(name) for name in names]
    while len(nodes) > 1:
        nodes = [
            "({},{}):1".format(nodes[index], nodes[index + 1])
            for index in range(0, len(nodes), 2)
        ]
    tree = Tree(nodes[0] + ";", parser=1)
    rng = np.random.default_rng(7)
    latent = rng.normal(0.0, 1.5, len(names))
    observed = latent + rng.normal(0.0, 1.0, len(names))
    response = 0.5 + 2.0 * latent + rng.normal(0.0, 0.3, len(names))

    def values(array):
        return dict(zip(names, array, strict=True))

    sampling = pd.DataFrame(np.eye(len(names)), index=names, columns=names)

    naive = fit_ordinary_pgls(
        tree,
        {"response": values(response)},
        {"predictor": values(observed)},
        ["response"],
        ["predictor"],
        evolution_model="independent",
    ).set_index("term")
    corrected = fit_ordinary_pgls(
        tree,
        {"response": values(response)},
        {"predictor": values(observed)},
        ["response"],
        ["predictor"],
        evolution_model="independent",
        predictor_evolution_model="independent",
        predictor_sampling_covariance={"predictor": sampling},
        reml=False,
    ).set_index("term")

    naive_error = abs(float(naive.loc["predictor", "coefficient"]) - 2.0)
    corrected_error = abs(float(corrected.loc["predictor", "coefficient"]) - 2.0)
    assert corrected_error < naive_error
    assert corrected.loc["predictor", "measurement_error_model"] == "latent-predictor"


def test_zero_predictor_sampling_covariance_recovers_exact_predictor_pgls():
    arguments = dict(
        tree=_tree(),
        response_values_by_trait={"expression": _values([2.0, 5.0, 7.5, 8.0, 12.5])},
        predictor_values_by_trait={"body_size": _values([1.0, 2.0, 4.0, 3.0, 7.0])},
        responses=["expression"],
        predictors=["body_size"],
        reml=False,
    )
    exact = fit_ordinary_pgls(**arguments)
    zero_sampling = pd.DataFrame(
        np.zeros((len(LEAF_NAMES), len(LEAF_NAMES))),
        index=LEAF_NAMES,
        columns=LEAF_NAMES,
    )
    latent = fit_ordinary_pgls(
        **arguments,
        predictor_sampling_covariance={"body_size": zero_sampling},
    )

    np.testing.assert_allclose(latent["coefficient"], exact["coefficient"], rtol=1e-6)
    np.testing.assert_allclose(
        latent["standard_error"], exact["standard_error"], rtol=5e-4
    )


def test_independent_correlation_matches_ordinary_least_squares():
    response = np.asarray([2.0, 5.0, 7.5, 8.0, 12.5])
    predictor = np.asarray([1.0, 2.0, 4.0, 3.0, 7.0])
    result = fit_ordinary_pgls(
        _tree(),
        {"expression": _values(response)},
        {"body_size": _values(predictor)},
        ["expression"],
        ["body_size"],
        evolution_model="independent",
    ).set_index("term")
    expected = np.linalg.lstsq(
        np.column_stack([np.ones(len(predictor)), predictor]),
        response,
        rcond=None,
    )[0]
    assert result.loc["(intercept)", "coefficient"] == pytest.approx(expected[0])
    assert result.loc["body_size", "coefficient"] == pytest.approx(expected[1])


def test_independent_pgls_ignores_invalid_input_branch_lengths():
    tree = Tree(TREE_TEXT.replace("A:1", "A:0"), parser=1)
    result = fit_ordinary_pgls(
        tree,
        {"expression": _values([2.0, 5.0, 7.5, 8.0, 12.5])},
        {"body_size": _values([1.0, 2.0, 4.0, 3.0, 7.0])},
        ["expression"],
        ["body_size"],
        evolution_model="independent",
    )
    assert set(result["branch_length_mode"]) == {"not-applicable"}


@pytest.mark.parametrize(
    ("evolution_model", "parameter"),
    [("lambda", 0.4), ("ou", 0.7)],
)
def test_fixed_phylogenetic_correlation_parameters_are_reported(
    evolution_model, parameter
):
    result = fit_ordinary_pgls(
        _tree(),
        {"expression": _values([2.0, 5.0, 7.5, 8.0, 12.5])},
        {"body_size": _values([1.0, 2.0, 4.0, 3.0, 7.0])},
        ["expression"],
        ["body_size"],
        evolution_model=evolution_model,
        evolution_parameter=parameter,
    )
    assert set(result["evolution_parameter_status"]) == {"fixed"}
    assert set(result["evolution_parameter"]) == {parameter}


@pytest.mark.parametrize("evolution_model", ["lambda", "ou"])
def test_phylogenetic_correlation_parameter_can_be_estimated(evolution_model):
    result = fit_ordinary_pgls(
        _tree(),
        {"expression": _values([2.0, 5.0, 7.5, 8.0, 12.5])},
        {"body_size": _values([1.0, 2.0, 4.0, 3.0, 7.0])},
        ["expression"],
        ["body_size"],
        evolution_model=evolution_model,
    )
    assert set(result["evolution_parameter_status"]) == {"estimated"}
    assert result.iloc[0]["evolution_parameter"] > 0.0 or (
        evolution_model == "lambda" and result.iloc[0]["evolution_parameter"] == 0.0
    )
    if evolution_model == "lambda":
        assert result.iloc[0]["evolution_parameter"] <= 1.0


def test_ordinary_pgls_parametric_bootstrap_is_reproducible():
    predictor_sampling = pd.DataFrame(
        np.eye(len(LEAF_NAMES)) * 0.04,
        index=LEAF_NAMES,
        columns=LEAF_NAMES,
    )
    arguments = dict(
        tree=_tree(),
        response_values_by_trait={"expression": _values([2.0, 5.0, 7.5, 8.0, 12.5])},
        predictor_values_by_trait={"body_size": _values([1.0, 2.0, 4.0, 3.0, 7.0])},
        responses=["expression"],
        predictors=["body_size"],
        evolution_model="lambda",
        predictor_sampling_covariance={"body_size": predictor_sampling},
        inference="parametric-bootstrap",
        bootstrap_replicates=6,
        seed=19,
    )
    first = fit_ordinary_pgls(**arguments)
    second = fit_ordinary_pgls(**arguments)
    pd.testing.assert_frame_equal(first, second)
    assert set(first["inference_method"]) == {"parametric-bootstrap"}
    assert set(first["measurement_error_model"]) == {"latent-predictor"}
    assert set(first["evolution_parameter_status"]) == {"estimated"}


def test_ordinary_sampling_covariance_requires_exact_unique_tree_tip_names():
    covariance = pd.DataFrame(
        np.eye(len(LEAF_NAMES) + 1),
        index=LEAF_NAMES + ["extra"],
        columns=LEAF_NAMES + ["extra"],
    )
    with pytest.raises(ValueError, match="exactly match tree tips"):
        fit_ordinary_pgls(
            _tree(),
            {"expression": _values([2.0, 5.0, 7.5, 8.0, 12.5])},
            {"body_size": _values([1.0, 2.0, 4.0, 3.0, 7.0])},
            ["expression"],
            ["body_size"],
            response_sampling_covariance={"expression": covariance},
        )


@pytest.mark.integration
def test_conventional_pgls_cli_reads_tree_and_tip_table(tmp_path):
    tree_path, data_path = _write_inputs(tmp_path)
    output_path = tmp_path / "pgls.tsv"

    main(
        [
            "pgls",
            "--tree",
            str(tree_path),
            "--data",
            str(data_path),
            "--responses",
            "expression",
            "--predictors",
            "body_size",
            "--outfile",
            str(output_path),
        ]
    )

    result = pd.read_csv(output_path, sep="\t")
    assert list(result["term"]) == ["(intercept)", "body_size"]
    assert set(result["model"]) == {"ordinary"}
    assert set(result["n_species"]) == {5}


@pytest.mark.integration
def test_conventional_pgls_cli_can_replace_invalid_branch_lengths_with_units(
    tmp_path,
):
    tree_path, data_path = _write_inputs(tmp_path)
    tree_path.write_text(TREE_TEXT.replace("A:1", "A:0"))
    output_path = tmp_path / "pgls.tsv"
    arguments = [
        "pgls",
        "--tree",
        str(tree_path),
        "--data",
        str(data_path),
        "--responses",
        "expression",
        "--predictors",
        "body_size",
        "--outfile",
        str(output_path),
    ]

    with pytest.raises(ValueError, match="positive finite"):
        main(arguments)

    main(arguments + ["--branch-length", "unit"])
    result = pd.read_csv(output_path, sep="\t")
    assert set(result["evolution_model"]) == {"brownian"}


@pytest.mark.integration
def test_conventional_pgls_cli_supports_multiple_traits_and_no_intercept(tmp_path):
    tree_path, data_path = _write_inputs(tmp_path)
    dataframe = pd.read_csv(data_path, sep="\t")
    dataframe["expression_alt"] = [9.0, 7.0, 8.0, 3.0, 2.0]
    dataframe["temperature"] = [5.0, 3.0, 6.0, 2.0, 8.0]
    dataframe.to_csv(data_path, sep="\t", index=False)
    output_path = tmp_path / "pgls.tsv"

    main(
        [
            "pgls",
            "--tree",
            str(tree_path),
            "--data",
            str(data_path),
            "--responses",
            "expression,expression_alt",
            "--predictors",
            "body_size,temperature",
            "--intercept",
            "no",
            "--outfile",
            str(output_path),
        ]
    )

    result = pd.read_csv(output_path, sep="\t")
    assert set(result["response"]) == {"expression", "expression_alt"}
    assert set(result["term"]) == {"body_size", "temperature"}
    assert set(result["intercept"]) == {"no"}
    assert len(result) == 4


@pytest.mark.integration
def test_conventional_pgls_cli_propagates_biological_replicate_uncertainty(
    tmp_path,
):
    tree_path, data_path = _write_inputs(tmp_path, biological_replicates=True)
    output_path = tmp_path / "pgls.tsv"
    covariance_path = tmp_path / "sampling.tsv"
    summary_path = tmp_path / "tip-summary.tsv"

    main(
        [
            "pgls",
            "--tree",
            str(tree_path),
            "--data",
            str(data_path),
            "--responses",
            "expression",
            "--predictors",
            "body_size",
            "--biological-id",
            "sample_id",
            "--outfile",
            str(output_path),
            "--sampling-covariance-out",
            str(covariance_path),
            "--tip-summary-out",
            str(summary_path),
        ]
    )

    result = pd.read_csv(output_path, sep="\t")
    covariance = pd.read_csv(covariance_path, sep="\t")
    summary = pd.read_csv(summary_path, sep="\t")
    assert result.iloc[0]["mean_sampling_variance"] > 0.0
    assert len(covariance) == 15
    assert set(summary["n_biological"]) == {2}
    assert set(summary["variance_method"]) == {"pooled"}


@pytest.mark.integration
def test_conventional_pgls_cli_supports_known_standard_errors(tmp_path):
    tree_path, data_path = _write_inputs(tmp_path)
    dataframe = pd.read_csv(data_path, sep="\t")
    dataframe["expression_se"] = [0.2, 0.3, 0.25, 0.4, 0.35]
    dataframe["expression_n"] = [4, 5, 6, 4, 5]
    dataframe.to_csv(data_path, sep="\t", index=False)
    output_path = tmp_path / "pgls.tsv"
    summary_path = tmp_path / "tip-summary.tsv"

    main(
        [
            "pgls",
            "--tree",
            str(tree_path),
            "--data",
            str(data_path),
            "--responses",
            "expression",
            "--predictors",
            "body_size",
            "--within-variance",
            "known-se",
            "--sample-size-columns",
            "expression_n",
            "--outfile",
            str(output_path),
            "--tip-summary-out",
            str(summary_path),
        ]
    )

    result = pd.read_csv(output_path, sep="\t")
    summary = pd.read_csv(summary_path, sep="\t")
    assert result.iloc[0]["mean_sampling_variance"] > 0.0
    assert set(summary["variance_method"]) == {"known-se"}
    np.testing.assert_allclose(summary["standard_error"], dataframe["expression_se"])


@pytest.mark.integration
def test_conventional_pgls_supports_response_and_predictor_known_se(tmp_path):
    tree_path, data_path = _write_inputs(tmp_path)
    dataframe = pd.read_csv(data_path, sep="\t")
    dataframe["expression_se"] = [0.2, 0.3, 0.25, 0.4, 0.35]
    dataframe["body_size_se"] = [0.1, 0.2, 0.15, 0.25, 0.2]
    dataframe.to_csv(data_path, sep="\t", index=False)
    output_path = tmp_path / "pgls.tsv"
    response_covariance = tmp_path / "response-covariance.tsv"
    predictor_covariance = tmp_path / "predictor-covariance.tsv"
    response_summary = tmp_path / "response-summary.tsv"
    predictor_summary = tmp_path / "predictor-summary.tsv"
    model_comparison = tmp_path / "model-comparison.tsv"

    main(
        [
            "pgls",
            "--tree",
            str(tree_path),
            "--data",
            str(data_path),
            "--responses",
            "expression",
            "--predictors",
            "body_size",
            "--within-variance",
            "known-se",
            "--predictor-within-variance",
            "known-se",
            "--outfile",
            str(output_path),
            "--sampling-covariance-out",
            str(response_covariance),
            "--predictor-sampling-covariance-out",
            str(predictor_covariance),
            "--tip-summary-out",
            str(response_summary),
            "--predictor-tip-summary-out",
            str(predictor_summary),
            "--compare-evolution-models",
            "brownian,independent",
            "--model-comparison-out",
            str(model_comparison),
        ]
    )

    result = pd.read_csv(output_path, sep="\t")
    slope = result[result["term"] == "body_size"].iloc[0]
    assert set(result["measurement_error_model"]) == {"latent-predictor"}
    assert slope["standard_error"] > 0
    assert slope["mean_predictor_sampling_variance"] > 0
    assert len(pd.read_csv(response_covariance, sep="\t")) == 15
    assert len(pd.read_csv(predictor_covariance, sep="\t")) == 15
    assert set(pd.read_csv(response_summary, sep="\t")["trait"]) == {"expression"}
    assert set(pd.read_csv(predictor_summary, sep="\t")["trait"]) == {"body_size"}
    assert set(pd.read_csv(model_comparison, sep="\t")["evolution_model"]) == {
        "brownian",
        "independent",
    }


@pytest.mark.integration
def test_conventional_pgls_supports_unpaired_replicate_rows_for_both_roles(tmp_path):
    tree_path, data_path = _write_inputs(tmp_path)
    original = pd.read_csv(data_path, sep="\t")
    rows = []
    for record in original.to_dict("records"):
        for replicate, offset in [("r1", -0.2), ("r2", 0.2)]:
            rows.append(
                {
                    "leaf_name": record["leaf_name"],
                    "response_sample": replicate,
                    "predictor_sample": "",
                    "expression": record["expression"] + offset,
                    "body_size": "",
                }
            )
        for replicate, offset in [("p1", -0.1), ("p2", 0.1)]:
            rows.append(
                {
                    "leaf_name": record["leaf_name"],
                    "response_sample": "",
                    "predictor_sample": replicate,
                    "expression": "",
                    "body_size": record["body_size"] + offset,
                }
            )
    pd.DataFrame(rows).to_csv(data_path, sep="\t", index=False)
    output_path = tmp_path / "unpaired.tsv"

    main(
        [
            "pgls",
            "--tree",
            str(tree_path),
            "--data",
            str(data_path),
            "--responses",
            "expression",
            "--predictors",
            "body_size",
            "--biological-id",
            "response_sample",
            "--predictor-biological-id",
            "predictor_sample",
            "--outfile",
            str(output_path),
        ]
    )

    result = pd.read_csv(output_path, sep="\t")
    assert set(result["measurement_error_model"]) == {"latent-predictor"}
    assert result[result["term"] == "body_size"].iloc[0]["standard_error"] > 0.0


@pytest.mark.integration
def test_conventional_pgls_supports_known_se_and_unpaired_raw_predictor_rows(
    tmp_path,
):
    tree_path, data_path = _write_inputs(tmp_path)
    original = pd.read_csv(data_path, sep="\t")
    rows = []
    for record in original.to_dict("records"):
        rows.append(
            {
                "leaf_name": record["leaf_name"],
                "predictor_sample": "",
                "expression": record["expression"],
                "expression_se": 0.2,
                "body_size": "",
            }
        )
        for replicate, offset in [("p1", -0.1), ("p2", 0.1)]:
            rows.append(
                {
                    "leaf_name": record["leaf_name"],
                    "predictor_sample": replicate,
                    "expression": "",
                    "expression_se": "",
                    "body_size": record["body_size"] + offset,
                }
            )
    pd.DataFrame(rows).to_csv(data_path, sep="\t", index=False)
    output_path = tmp_path / "mixed-replicates.tsv"

    main(
        [
            "pgls",
            "--tree",
            str(tree_path),
            "--data",
            str(data_path),
            "--responses",
            "expression",
            "--predictors",
            "body_size",
            "--within-variance",
            "known-se",
            "--predictor-biological-id",
            "predictor_sample",
            "--outfile",
            str(output_path),
        ]
    )

    result = pd.read_csv(output_path, sep="\t")
    assert set(result["measurement_error_model"]) == {"latent-predictor"}
    assert result.iloc[0]["mean_sampling_variance"] > 0.0
    assert (
        result[result["term"] == "body_size"].iloc[0][
            "mean_predictor_sampling_variance"
        ]
        > 0.0
    )


@pytest.mark.integration
def test_conventional_pgls_cli_supports_new_models_custom_covariance_and_comparison(
    tmp_path,
):
    tree_path, data_path = _write_inputs(tmp_path)
    covariance_path = tmp_path / "evolution-covariance.tsv"
    covariance = build_phylogenetic_covariance(_tree(), LEAF_NAMES)
    covariance_frame = pd.DataFrame(
        covariance,
        index=LEAF_NAMES,
        columns=LEAF_NAMES,
    )
    covariance_frame.index.name = "leaf_name"
    covariance_frame.to_csv(covariance_path, sep="\t")
    output_path = tmp_path / "pgls.tsv"
    comparison_path = tmp_path / "comparison.tsv"

    main(
        [
            "pgls",
            "--tree",
            str(tree_path),
            "--data",
            str(data_path),
            "--responses",
            "expression",
            "--predictors",
            "body_size",
            "--evolution-model",
            "kappa",
            "--evolution-parameter",
            "0.8",
            "--evolution-covariance",
            str(covariance_path),
            "--compare-evolution-models",
            "brownian,independent,custom",
            "--model-comparison-out",
            str(comparison_path),
            "--outfile",
            str(output_path),
        ]
    )

    result = pd.read_csv(output_path, sep="\t")
    comparison = pd.read_csv(comparison_path, sep="\t")
    assert set(result["evolution_model"]) == {"kappa"}
    assert set(result["evolution_parameter"]) == {0.8}
    assert set(result["branch_length_mode"]) == {"original"}
    assert list(comparison["evolution_model"]) == [
        "brownian",
        "independent",
        "custom",
    ]
    assert (
        comparison.set_index("evolution_model").loc["custom", "branch_length_mode"]
        == "not-applicable"
    )
    assert comparison["akaike_weight"].sum() == pytest.approx(1.0)
    assert comparison["aicc_weight"].sum() == pytest.approx(1.0)
    assert np.isfinite(comparison[["aic", "aicc", "bic"]]).all(axis=None)


def test_model_comparison_estimates_parameterized_models_by_ml():
    result = fit_ordinary_model_comparison(
        _tree(),
        {"expression": _values([2.0, 5.0, 7.5, 8.0, 12.5])},
        {"body_size": _values([1.0, 2.0, 4.0, 3.0, 7.0])},
        ["expression"],
        ["body_size"],
        ["lambda", "kappa", "delta", "eb", "acdc", "ou"],
    )
    assert set(result["evolution_parameter_status"]) == {"estimated"}
    assert np.isfinite(result["evolution_parameter"]).all()
    assert np.isfinite(result[["log_likelihood", "aic", "bic"]]).all(axis=None)
    assert result["aicc"].isna().all()
    assert result["akaike_weight"].sum() == pytest.approx(1.0)
    assert result["aicc_weight"].eq("").all()


def test_model_comparison_rejects_a_single_string_as_a_model_sequence():
    with pytest.raises(ValueError, match="non-empty unique sequence"):
        fit_ordinary_model_comparison(
            _tree(),
            {"expression": _values([2.0, 5.0, 7.5, 8.0, 12.5])},
            {"body_size": _values([1.0, 2.0, 4.0, 3.0, 7.0])},
            ["expression"],
            ["body_size"],
            "brownian",
        )


def test_conventional_pgls_rejects_predictors_that_vary_among_replicates(tmp_path):
    tree_path, data_path = _write_inputs(tmp_path, biological_replicates=True)
    dataframe = pd.read_csv(data_path, sep="\t")
    dataframe.loc[1, "body_size"] += 1.0
    dataframe.to_csv(data_path, sep="\t", index=False)

    with pytest.raises(ValueError, match="differs among rows"):
        main(
            [
                "pgls",
                "--tree",
                str(tree_path),
                "--data",
                str(data_path),
                "--responses",
                "expression",
                "--predictors",
                "body_size",
                "--biological-id",
                "sample_id",
            ]
        )


def test_conventional_pgls_output_cannot_overwrite_an_input(tmp_path):
    tree_path, data_path = _write_inputs(tmp_path)
    original_data = data_path.read_text()

    with pytest.raises(ValueError, match="must not overwrite an input"):
        main(
            [
                "pgls",
                "--tree",
                str(tree_path),
                "--data",
                str(data_path),
                "--responses",
                "expression",
                "--predictors",
                "body_size",
                "--outfile",
                str(data_path),
            ]
        )

    assert data_path.read_text() == original_data


def test_model_comparison_output_cannot_overwrite_custom_covariance(tmp_path):
    tree_path, data_path = _write_inputs(tmp_path)
    covariance_path = tmp_path / "covariance.tsv"
    covariance = pd.DataFrame(
        build_phylogenetic_covariance(_tree(), LEAF_NAMES),
        index=LEAF_NAMES,
        columns=LEAF_NAMES,
    )
    covariance.index.name = "leaf_name"
    covariance.to_csv(covariance_path, sep="\t")
    original_covariance = covariance_path.read_text()

    with pytest.raises(ValueError, match="must not overwrite an input"):
        main(
            [
                "pgls",
                "--tree",
                str(tree_path),
                "--data",
                str(data_path),
                "--responses",
                "expression",
                "--predictors",
                "body_size",
                "--evolution-model",
                "custom",
                "--evolution-covariance",
                str(covariance_path),
                "--compare-evolution-models",
                "custom",
                "--model-comparison-out",
                str(covariance_path),
            ]
        )

    assert covariance_path.read_text() == original_covariance


@pytest.mark.integration
def test_custom_covariance_and_model_comparison_are_recorded_in_audit(tmp_path):
    tree_path, data_path = _write_inputs(tmp_path)
    covariance_path = tmp_path / "covariance.tsv"
    covariance = pd.DataFrame(
        build_phylogenetic_covariance(_tree(), LEAF_NAMES),
        index=LEAF_NAMES,
        columns=LEAF_NAMES,
    )
    covariance.index.name = "leaf_name"
    covariance.to_csv(covariance_path, sep="\t")
    output_path = tmp_path / "pgls.tsv"
    comparison_path = tmp_path / "comparison.tsv"
    audit_path = tmp_path / "audit.jsonl"

    main(
        [
            "pgls",
            "--tree",
            str(tree_path),
            "--data",
            str(data_path),
            "--responses",
            "expression",
            "--predictors",
            "body_size",
            "--evolution-model",
            "custom",
            "--evolution-covariance",
            str(covariance_path),
            "--compare-evolution-models",
            "brownian,custom",
            "--model-comparison-out",
            str(comparison_path),
            "--outfile",
            str(output_path),
            "--audit",
            str(audit_path),
        ]
    )

    record = json.loads(audit_path.read_text())
    input_paths = {item["path"] for item in record["inputs"]}
    output_paths = {item["path"] for item in record["outputs"]}
    assert input_paths == {
        str(tree_path.resolve()),
        str(data_path.resolve()),
        str(covariance_path.resolve()),
    }
    assert output_paths == {str(output_path.resolve()), str(comparison_path.resolve())}


def test_conventional_pgls_mode_rejects_incomplete_mixed_and_invalid_options():
    common = ["pgls", "--responses", "expression", "--predictors", "body_size"]
    with pytest.raises(ValueError, match="Conventional PGLS requires"):
        main(common + ["--tree", "species.nwk"])
    with pytest.raises(ValueError, match="cannot use reconciled/precomputed"):
        main(
            common
            + [
                "--tree",
                "species.nwk",
                "--data",
                "data.tsv",
                "--gene-tree",
                "gene.nwk",
            ]
        )
    with pytest.raises(ValueError, match="gene-evolution-model"):
        main(
            common
            + [
                "--tree",
                "species.nwk",
                "--data",
                "data.tsv",
                "--gene-evolution-model",
                "kappa",
            ]
        )
    with pytest.raises(ValueError, match="event-weighting"):
        main(
            common
            + [
                "--tree",
                "species.nwk",
                "--data",
                "data.tsv",
                "--event-weighting",
                "equal",
            ]
        )
    with pytest.raises(ValueError, match="does not take a parameter"):
        main(
            common
            + [
                "--tree",
                "species.nwk",
                "--data",
                "data.tsv",
                "--evolution-parameter",
                "0.5",
            ]
        )
    with pytest.raises(ValueError, match="must be used together"):
        main(
            common
            + [
                "--tree",
                "species.nwk",
                "--data",
                "data.tsv",
                "--compare-evolution-models",
                "brownian,lambda",
            ]
        )
    with pytest.raises(ValueError, match="requires '--evolution-covariance'"):
        main(
            common
            + [
                "--tree",
                "species.nwk",
                "--data",
                "data.tsv",
                "--evolution-model",
                "custom",
            ]
        )
    with pytest.raises(ValueError, match="requires a custom selected"):
        main(
            common
            + [
                "--tree",
                "species.nwk",
                "--data",
                "data.tsv",
                "--evolution-covariance",
                "covariance.tsv",
            ]
        )
    with pytest.raises(ValueError, match="must not be empty"):
        main(
            common
            + [
                "--tree",
                "species.nwk",
                "--data",
                "data.tsv",
                "--compare-evolution-models",
                "",
                "--model-comparison-out",
                "comparison.tsv",
            ]
        )


@pytest.mark.integration
def test_conventional_pgls_tree_stdin_is_summarized_in_audit(monkeypatch, tmp_path):
    _, data_path = _write_inputs(tmp_path)
    monkeypatch.setattr(sys, "stdin", io.StringIO(TREE_TEXT))
    output_path = tmp_path / "pgls.tsv"
    audit_path = tmp_path / "audit.jsonl"

    main(
        [
            "pgls",
            "--tree",
            "-",
            "--data",
            str(data_path),
            "--responses",
            "expression",
            "--predictors",
            "body_size",
            "--outfile",
            str(output_path),
            "--audit",
            str(audit_path),
        ]
    )

    record = json.loads(audit_path.read_text())
    assert record["stdin"] == {
        "argument": "tree",
        "sha256": hashlib.sha256(TREE_TEXT.encode()).hexdigest(),
        "bytes": len(TREE_TEXT.encode()),
    }
    assert record["primary_input"]["kind"] == "newick"
    assert record["primary_input"]["first_tree_tip_count"] == 5


@pytest.mark.parametrize(
    ("response_values", "options", "family"),
    [
        (
            ["absent", "present", "absent", "present", "present"],
            {},
            "binomial",
        ),
        (
            ["low", "middle", "high", "low", "high"],
            {"categorical_responses": ["state"]},
            "multinomial",
        ),
        (
            ["low", "middle", "high", "low", "high"],
            {"ordered_responses": {"state": ("low", "middle", "high")}},
            "ordinal",
        ),
    ],
)
def test_conventional_phylogenetic_glmm_response_families(
    response_values, options, family
):
    result = fit_ordinary_pgls(
        _tree(),
        {"state": _values(response_values)},
        {"body_size": _values([1.0, 2.0, 4.0, 3.0, 7.0])},
        ["state"],
        ["body_size"],
        **options,
    )

    coefficients = result[result["term_test"] == "coefficient"]
    assert set(coefficients["response_family"]) == {family}
    assert set(coefficients["model"]) == {"ordinary-pglmm"}
    assert np.isfinite(pd.to_numeric(coefficients["coefficient"])).all()
    if family == "ordinal":
        assert "(intercept)" not in set(coefficients["term"])
        assert len(result[result["term_test"] == "threshold"]) == 2


@pytest.mark.integration
def test_conventional_categorical_response_and_predictor_replicates(tmp_path):
    tree_path = tmp_path / "species.nwk"
    data_path = tmp_path / "categorical-replicates.tsv"
    output_path = tmp_path / "categorical-pglmm.tsv"
    response_summary_path = tmp_path / "response-summary.tsv"
    predictor_summary_path = tmp_path / "predictor-summary.tsv"
    tree_path.write_text(TREE_TEXT)
    rows = []
    response_states = [
        ("absent", "present"),
        ("absent", "absent"),
        ("present", "absent"),
        ("present", "present"),
        ("present", "absent"),
    ]
    habitats = [
        ("water", "land"),
        ("land", "land"),
        ("water", "water"),
        ("land", "land"),
        ("water", "water"),
    ]
    for leaf, states, habitat_states in zip(
        LEAF_NAMES, response_states, habitats, strict=True
    ):
        for replicate in range(2):
            rows.append(
                {
                    "leaf_name": leaf,
                    "sample": "{}_{}".format(leaf, replicate + 1),
                    "state": states[replicate],
                    "habitat": habitat_states[replicate],
                }
            )
    pd.DataFrame(rows).to_csv(data_path, sep="\t", index=False)

    main(
        [
            "pgls",
            "--tree",
            str(tree_path),
            "--data",
            str(data_path),
            "--responses",
            "state",
            "--predictors",
            "habitat",
            "--categorical-responses",
            "state",
            "--categorical-predictors",
            "habitat",
            "--biological-id",
            "sample",
            "--predictor-biological-id",
            "sample",
            "--categorical-replicate-policy",
            "latent",
            "--tip-summary-out",
            str(response_summary_path),
            "--predictor-tip-summary-out",
            str(predictor_summary_path),
            "--outfile",
            str(output_path),
        ]
    )

    result = pd.read_csv(output_path, sep="\t")
    coefficients = result[result["term_test"] == "coefficient"]
    assert set(coefficients["response_family"]) == {"binomial"}
    assert set(coefficients["measurement_error_model"]) == {"latent-predictor"}
    response_summary = pd.read_csv(response_summary_path, sep="\t")
    assert set(response_summary["variance_method"]) == {"categorical-counts"}
    predictor_summary = pd.read_csv(predictor_summary_path, sep="\t")
    assert "categorical-latent" in set(predictor_summary["variance_method"])


@pytest.mark.parametrize(
    ("family", "values", "extra"),
    [
        ("poisson", [1, 2, 2, 4, 5], {}),
        ("negative-binomial", [0, 1, 3, 4, 8], {}),
        ("zero-inflated-poisson", [0, 0, 1, 3, 6], {}),
        ("zero-inflated-negative-binomial", [0, 0, 1, 4, 9], {}),
        ("hurdle-poisson", [0, 1, 2, 3, 6], {}),
        ("hurdle-negative-binomial", [0, 1, 2, 3, 7], {}),
        ("gamma", [1.0, 1.5, 2.0, 3.0, 4.5], {}),
        ("lognormal", [1.0, 1.5, 2.0, 3.0, 4.5], {}),
        ("beta", [0.1, 0.2, 0.4, 0.6, 0.8], {}),
        (
            "beta-binomial",
            [1, 2, 4, 6, 8],
            {"response_trials": {"state": _values([10, 10, 10, 10, 10])}},
        ),
    ],
)
def test_conventional_scalar_non_gaussian_response_families(family, values, extra):
    result = fit_ordinary_pgls(
        _tree(),
        {"state": _values(values)},
        {"body_size": _values([1.0, 2.0, 4.0, 3.0, 7.0])},
        ["state"],
        ["body_size"],
        response_families={"state": family},
        **extra,
    )

    assert set(result["response_family"]) == {family}
    assert set(result["term_test"]) == {"coefficient"}
    assert np.isfinite(pd.to_numeric(result["coefficient"])).all()
    assert set(result["coefficient_penalty"]) == {"student-t"}


def test_conventional_censored_gaussian_uses_bounds_for_missing_observations():
    result = fit_ordinary_pgls(
        _tree(),
        {"state": _values([1.0, np.nan, 3.0, np.nan, 5.0])},
        {"body_size": _values([1.0, 2.0, 4.0, 3.0, 7.0])},
        ["state"],
        ["body_size"],
        response_families={"state": "censored-gaussian"},
        response_censor_lower={"state": _values([np.nan, np.nan, np.nan, 3.5, np.nan])},
        response_censor_upper={"state": _values([np.nan, 2.5, np.nan, np.nan, np.nan])},
    )

    assert set(result["response_family"]) == {"censored-gaussian"}
    assert set(result["link_function"]) == {"identity"}
    assert np.isfinite(pd.to_numeric(result["response_dispersion"])).all()
    assert np.isfinite(pd.to_numeric(result["coefficient"])).all()


def test_response_family_configuration_rejects_contradictory_or_ignored_options():
    tree = _tree()
    predictors = {"body_size": _values([1.0, 2.0, 4.0, 3.0, 7.0])}
    with pytest.raises(ValueError, match="unordered categorical"):
        fit_ordinary_pgls(
            tree,
            {"state": _values([0, 1, 0, 1, 1])},
            predictors,
            ["state"],
            ["body_size"],
            categorical_responses=["state"],
            response_families={"state": "gaussian"},
        )
    with pytest.raises(ValueError, match="offset"):
        fit_ordinary_pgls(
            tree,
            {"state": _values([1.0, 1.5, 2.0, 3.0, 4.5])},
            predictors,
            ["state"],
            ["body_size"],
            response_families={"state": "gamma"},
            response_offsets={"state": _values([0.0] * 5)},
        )
    with pytest.raises(ValueError, match="strictly in"):
        fit_ordinary_pgls(
            tree,
            {"state": _values([0, 0, 1, 3, 6])},
            predictors,
            ["state"],
            ["body_size"],
            response_families={"state": "zero-inflated-poisson"},
            response_zero_probabilities={"state": 1.0},
        )
    with pytest.raises(ValueError, match="finite or missing"):
        fit_ordinary_pgls(
            tree,
            {"state": _values([1.0, 1.5, 2.0, 3.0, 4.5])},
            predictors,
            ["state"],
            ["body_size"],
            response_families={"state": "censored-gaussian"},
            response_censor_lower={
                "state": _values([np.inf, np.nan, np.nan, np.nan, np.nan])
            },
        )
    with pytest.raises(ValueError, match="must have a missing response"):
        fit_ordinary_pgls(
            tree,
            {"state": _values([1.0, 1.5, 2.0, 3.0, 4.5])},
            predictors,
            ["state"],
            ["body_size"],
            response_families={"state": "censored-gaussian"},
            response_censor_lower={
                "state": _values([0.5, np.nan, np.nan, np.nan, np.nan])
            },
        )
    with pytest.raises(ValueError, match="positive finite coefficient prior SD"):
        fit_ordinary_pgls(
            tree,
            {"state": _values([1, 1, 2, 3, 5])},
            predictors,
            ["state"],
            ["body_size"],
            response_families={"state": "poisson"},
            coefficient_prior_sd=float("nan"),
        )


def test_categorical_separation_is_regularized_and_likelihood_tested():
    result = fit_ordinary_pgls(
        _tree(),
        {"state": _values(["no", "no", "no", "yes", "yes"])},
        {"body_size": _values([1.0, 2.0, 3.0, 4.0, 5.0])},
        ["state"],
        ["body_size"],
        inference="likelihood-ratio",
    )

    assert set(result["separation_warning"]) == {"yes"}
    assert set(result["inference_method"]) == {"likelihood-ratio"}
    assert np.isfinite(pd.to_numeric(result["coefficient"])).all()
    assert np.isfinite(pd.to_numeric(result["p_value"])).all()


def test_non_gaussian_profile_likelihood_reports_asymmetric_intervals():
    result = fit_ordinary_pgls(
        _tree(),
        {"state": _values([1, 2, 2, 4, 5])},
        {"body_size": _values([1.0, 2.0, 4.0, 3.0, 7.0])},
        ["state"],
        ["body_size"],
        response_families={"state": "poisson"},
        inference="profile-likelihood",
    )

    assert set(result["inference_method"]) == {"profile-likelihood"}
    lower = pd.to_numeric(result["confidence_interval_lower"])
    upper = pd.to_numeric(result["confidence_interval_upper"])
    estimates = pd.to_numeric(result["coefficient"])
    assert np.isfinite(lower).all() and np.isfinite(upper).all()
    assert (lower < estimates).all() and (estimates < upper).all()


@pytest.mark.parametrize(
    ("values", "response_families"),
    [
        (["no", "no", "yes", "yes", "yes"], None),
        ([1, 2, 2, 4, 5], {"state": "poisson"}),
    ],
)
def test_non_gaussian_parametric_bootstrap_refits_the_family(values, response_families):
    result = fit_ordinary_pgls(
        _tree(),
        {"state": _values(values)},
        {"body_size": _values([1.0, 2.0, 4.0, 3.0, 7.0])},
        ["state"],
        ["body_size"],
        response_families=response_families,
        inference="parametric-bootstrap",
        bootstrap_replicates=4,
        seed=17,
    )

    coefficients = result[result["term_test"] == "coefficient"]
    assert set(coefficients["inference_method"]) == {"parametric-bootstrap"}
    assert np.isfinite(pd.to_numeric(coefficients["standard_error"])).all()
    assert pd.to_numeric(coefficients["p_value"]).between(0.0, 1.0).all()
    assert np.isfinite(pd.to_numeric(coefficients["confidence_interval_lower"])).all()


def test_multivariate_gaussian_pgls_retains_partially_observed_tips():
    responses = {
        "first": _values([1.0, 2.0, 3.0, 4.0, 5.0]),
        "second": _values([2.0, 3.5, np.nan, 6.5, 8.0]),
    }
    result = fit_ordinary_pgls(
        _tree(),
        responses,
        {"body_size": _values([1.0, 2.0, 4.0, 3.0, 7.0])},
        ["first", "second"],
        ["body_size"],
        multivariate_responses=True,
        allow_missing_responses=True,
    )

    coefficients = result[result["term_test"] == "coefficient"]
    assert set(coefficients["response"]) == {"first", "second"}
    assert set(coefficients["model"]) == {"ordinary-multivariate-pgls"}
    covariance_rows = result[result["term_test"] == "response-covariance"]
    assert len(covariance_rows) == 3
    assert set(result["reml"]) == {"yes"}


def test_missing_response_flag_requires_multivariate_model():
    with pytest.raises(ValueError, match="requires multivariate_responses"):
        fit_ordinary_pgls(
            _tree(),
            {"state": _values([1.0, 2.0, 3.0, 4.0, 5.0])},
            {"body_size": _values([1.0, 2.0, 4.0, 3.0, 7.0])},
            ["state"],
            ["body_size"],
            allow_missing_responses=True,
        )


def test_multivariate_cli_combines_biological_replicates_with_partial_missingness(
    tmp_path,
):
    tree_path = tmp_path / "species.nwk"
    data_path = tmp_path / "multivariate-replicates.tsv"
    output_path = tmp_path / "multivariate.tsv"
    tree_path.write_text(TREE_TEXT)
    rows = []
    for index, (leaf_name, predictor) in enumerate(
        zip(LEAF_NAMES, [1.0, 2.0, 4.0, 3.0, 7.0], strict=True), start=1
    ):
        for replicate in range(2):
            rows.append(
                {
                    "leaf_name": leaf_name,
                    "sample": "{}{}".format(leaf_name, replicate),
                    "first": index + replicate * 0.2,
                    "second": (
                        np.nan if leaf_name == "C" else 2.0 * index + replicate * 0.2
                    ),
                    "body_size": predictor,
                }
            )
    pd.DataFrame(rows).to_csv(data_path, sep="\t", index=False)

    main(
        [
            "pgls",
            "--tree",
            str(tree_path),
            "--data",
            str(data_path),
            "--responses",
            "first,second",
            "--predictors",
            "body_size",
            "--biological-id",
            "sample",
            "--multivariate-responses",
            "yes",
            "--allow-missing-responses",
            "yes",
            "--outfile",
            str(output_path),
        ]
    )

    result = pd.read_csv(output_path, sep="\t")
    assert set(result["model"]) == {"ordinary-multivariate-pgls"}
    assert set(result.loc[result["term_test"] == "coefficient", "response"]) == {
        "first",
        "second",
    }


def test_censored_gaussian_cli_preserves_censored_biological_replicates(tmp_path):
    tree_path = tmp_path / "species.nwk"
    data_path = tmp_path / "censored-replicates.tsv"
    output_path = tmp_path / "censored.tsv"
    tree_path.write_text(TREE_TEXT)
    rows = []
    for index, (leaf_name, predictor) in enumerate(
        zip(LEAF_NAMES, [1.0, 2.0, 4.0, 3.0, 7.0], strict=True), start=1
    ):
        rows.extend(
            [
                {
                    "leaf_name": leaf_name,
                    "sample": "{}-exact".format(leaf_name),
                    "state": float(index),
                    "upper": np.nan,
                    "body_size": predictor,
                },
                {
                    "leaf_name": leaf_name,
                    "sample": "{}-left".format(leaf_name),
                    "state": np.nan,
                    "upper": float(index) + 0.5,
                    "body_size": predictor,
                },
            ]
        )
    pd.DataFrame(rows).to_csv(data_path, sep="\t", index=False)

    main(
        [
            "pgls",
            "--tree",
            str(tree_path),
            "--data",
            str(data_path),
            "--responses",
            "state",
            "--predictors",
            "body_size",
            "--biological-id",
            "sample",
            "--response-family",
            "state=censored-gaussian",
            "--response-censor-upper",
            "state=upper",
            "--outfile",
            str(output_path),
        ]
    )

    result = pd.read_csv(output_path, sep="\t")
    assert set(result["response_family"]) == {"censored-gaussian"}
    assert np.isfinite(result["coefficient"]).all()
