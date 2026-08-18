import hashlib
import io
import json
import sys
import threading
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from ete4 import Tree
from scipy import sparse

from nwkit import regress as regression_mod
from nwkit import regression_pipeline as regression_pipeline_mod
from nwkit.cli import main
from nwkit.contrast import build_contrast_table
from nwkit.gaussian import (
    NestedLowRankFactor,
    SparseCovarianceFactor,
    factor_diagonal_low_rank_updates,
    factor_logdet,
    solve_factor,
)
from nwkit.model_matrix import PredictorTerm
from nwkit.reconcile import build_reconciliation_table
from nwkit.regress import _profile_covariance_fit, fit_reconciled_pgls
from nwkit.regression_pipeline import (
    RegressionPipelineArtifacts,
    write_regression_bundle,
)
from nwkit.sparse_laplace import (
    JointPredictorUncertainty,
    continuous_predictor_loading,
)


def _write_raw_regression_inputs(tmp_path, *, biological_replicates=False):
    gene_tree = tmp_path / "gene.nwk"
    species_tree = tmp_path / "species.nwk"
    expression = tmp_path / "expression.tsv"
    species_traits = tmp_path / "species-traits.tsv"
    gene_names = ["Genus_{}_g1".format(letter) for letter in "abcde"]
    species_names = ["Genus_{}".format(letter) for letter in "abcde"]
    gene_tree.write_text("((({}:1,{}:1):1,{}:2):1,({}:1,{}:1):2);".format(*gene_names))
    species_tree.write_text(
        "((({}:1,{}:1):1,{}:2):1,({}:1,{}:1):2);".format(*species_names)
    )
    expression_values = [2.0, 5.0, 7.0, 8.0, 12.0]
    if biological_replicates:
        rows = []
        for leaf_name, value in zip(gene_names, expression_values, strict=True):
            rows.extend(
                [
                    {
                        "leaf_name": leaf_name,
                        "sample_id": "{}_1".format(leaf_name),
                        "expression": value - 0.5,
                    },
                    {
                        "leaf_name": leaf_name,
                        "sample_id": "{}_2".format(leaf_name),
                        "expression": value + 0.5,
                    },
                ]
            )
    else:
        rows = [
            {"leaf_name": leaf_name, "expression": value}
            for leaf_name, value in zip(gene_names, expression_values, strict=True)
        ]
    pd.DataFrame(rows).to_csv(expression, sep="\t", index=False)
    pd.DataFrame(
        {
            "leaf_name": species_names,
            "body_size": [1.0, 2.0, 4.0, 3.0, 7.0],
        }
    ).to_csv(species_traits, sep="\t", index=False)
    return gene_tree, species_tree, expression, species_traits


def _predictor_table(values=(1.0, 2.0, 3.0)):
    rows = []
    for index, value in enumerate(values, start=1):
        rows.append(
            {
                "tree_id": "species",
                "branch_clade_id": "event{}".format(index),
                "descendant_taxa": "taxa{}".format(index),
                "numerator_clade_id": "num{}".format(index),
                "denominator_clade_id": "den{}".format(index),
                "trait": "body_size",
                "evolution_model": "brownian",
                "evolution_parameter_name": "",
                "evolution_parameter": "",
                "branch_length_mode": "original",
                "raw_contrast": value,
            }
        )
    return pd.DataFrame(rows)


def _response_row(event_index, value, *, gene_index=1, tree_id="OG1"):
    return {
        "tree_id": tree_id,
        "gene_clade_id": "gene{}_{}".format(event_index, gene_index),
        "lineage_clade_id": "lineage{}".format(gene_index),
        "event_type": "speciation",
        "eligible": "yes",
        "coverage_status": "complete",
        "species_event_id": "event{}".format(event_index),
        "species_event_taxa": "taxa{}".format(event_index),
        "species_numerator_event_id": "num{}".format(event_index),
        "species_denominator_event_id": "den{}".format(event_index),
        "trait": "expression",
        "evolution_model": "brownian",
        "evolution_parameter_name": "",
        "evolution_parameter": "",
        "branch_length_mode": "original",
        "raw_contrast": value,
        "contrast_variance": 1.0,
    }


def _sampling_covariance_table(response, matrix):
    ids = response["gene_clade_id"].tolist()
    rows = []
    for first in range(len(ids)):
        for second in range(first, len(ids)):
            rows.append(
                {
                    "tree_id": response.iloc[first]["tree_id"],
                    "trait": response.iloc[first]["trait"],
                    "contrast_id_1": ids[first],
                    "contrast_id_2": ids[second],
                    "sampling_covariance": matrix[first, second],
                }
            )
    return pd.DataFrame(rows)


def test_hierarchical_gaussian_components_remain_unmaterialized():
    design = np.asarray([[1.0], [1.5], [2.0], [2.5], [3.0], [3.5]])
    event_inverse = np.asarray([0, 0, 1, 1, 2, 2])
    lineage_inverse = np.arange(6)
    output = regression_mod._build_covariance_components(
        design,
        np.ones(6),
        np.zeros(6),
        event_inverse,
        np.asarray([2, 2, 2]),
        lineage_inverse,
        np.ones(6, dtype=int),
        np.asarray(["e1", "e2", "e3"]),
        np.asarray(["l1", "l2", "l3", "l4", "l5", "l6"]),
        ["body_size"],
        n_events=3,
        num_parameters=1,
        event_weighting="none",
        model="hierarchical",
        event_random_effect="auto",
        lineage_random_slope="no",
    )
    fixed, components, factors, raw_components, _designs, use_event, _use_lineage = (
        output
    )
    assert fixed.ndim == 1
    assert components[0][1].ndim == 1
    assert dict(components)["species_event_variance"] is None
    assert factors["species_event_variance"].shape == (6, 3)
    assert raw_components["species_event_variance"].ndim == 1
    assert use_event


def test_collinear_hierarchical_covariance_components_are_not_both_fitted():
    design = np.ones((6, 1))
    grouping = np.asarray([0, 0, 1, 1, 2, 2])
    counts = np.asarray([2, 2, 2])
    labels = np.asarray(["a", "b", "c"])
    output = regression_mod._build_covariance_components(
        design,
        np.ones(6),
        np.zeros(6),
        grouping,
        counts,
        grouping,
        counts,
        labels,
        labels,
        ["intercept"],
        n_events=3,
        num_parameters=1,
        event_weighting="none",
        model="hierarchical",
        event_random_effect="auto",
        lineage_random_slope="auto",
    )
    assert output[-2:] == (True, False)
    with pytest.raises(ValueError, match="lineage random slope.*not identifiable"):
        regression_mod._build_covariance_components(
            design,
            np.ones(6),
            np.zeros(6),
            grouping,
            counts,
            grouping,
            counts,
            labels,
            labels,
            ["intercept"],
            n_events=3,
            num_parameters=1,
            event_weighting="none",
            model="hierarchical",
            event_random_effect="yes",
            lineage_random_slope="yes",
        )


def test_large_dense_gaussian_fit_is_rejected_before_optimization(monkeypatch):
    monkeypatch.setattr(regression_mod, "MAX_DENSE_GAUSSIAN_OBSERVATIONS", 3)
    dense_covariance = np.eye(4)
    dense_covariance[0, 1] = dense_covariance[1, 0] = 0.1
    with pytest.raises(ValueError, match="Dense Gaussian covariance fitting"):
        _profile_covariance_fit(
            np.arange(4.0),
            np.ones((4, 1)),
            dense_covariance,
            [("evolutionary_rate", np.ones(4))],
            reml=False,
        )


def test_large_structured_gaussian_fit_remains_available(monkeypatch):
    monkeypatch.setattr(regression_mod, "MAX_DENSE_GAUSSIAN_OBSERVATIONS", 3)
    result = _profile_covariance_fit(
        np.asarray([1.0, 2.0, 4.0, 8.0]),
        np.ones((4, 1)),
        np.zeros(4),
        [("evolutionary_rate", np.ones(4))],
        reml=False,
    )
    assert np.isfinite(result["objective"])


def test_pgls_matches_standard_pic_regression_for_one_to_one_orthologs():
    species_tree = Tree("(((A:1,B:1):1,C:2):1,(D:1,E:1):2);", parser=1)
    gene_tree = Tree("(((A:1,B:1):1,C:2):1,(D:1,E:1):2);", parser=1)
    mapping = {name: name for name in ["A", "B", "C", "D", "E"]}
    reconciliation = build_reconciliation_table(
        gene_tree, species_tree, mapping, tree_id="OG1"
    )
    reconciliation_by_id = {
        row["gene_clade_id"]: row for row in reconciliation.to_dict("records")
    }
    predictor = build_contrast_table(
        species_tree,
        {"body_size": {"A": 1, "B": 2, "C": 4, "D": 3, "E": 7}},
        tree_id="species",
    )
    response = build_contrast_table(
        gene_tree,
        {"expression": {"A": 2, "B": 5, "C": 7, "D": 8, "E": 12}},
        reconciliation_by_id=reconciliation_by_id,
        event_type="speciation",
        tree_id="OG1",
    )

    result = fit_reconciled_pgls(response, predictor, ["expression"], ["body_size"])

    joined = response.merge(
        predictor[predictor["trait"] == "body_size"],
        left_on="species_event_id",
        right_on="branch_clade_id",
        suffixes=("_response", "_predictor"),
    )
    scale = np.sqrt(joined["contrast_variance_response"].to_numpy(float))
    y = joined["raw_contrast_response"].to_numpy(float) / scale
    x = joined["raw_contrast_predictor"].to_numpy(float) / scale
    expected = float((x @ y) / (x @ x))
    assert result.iloc[0]["coefficient"] == pytest.approx(expected)
    tip_x = np.asarray([1, 2, 4, 3, 7], dtype=float)
    tip_y = np.asarray([2, 5, 7, 8, 12], dtype=float)
    covariance = np.asarray(
        [
            [3, 2, 1, 0, 0],
            [2, 3, 1, 0, 0],
            [1, 1, 3, 0, 0],
            [0, 0, 0, 3, 2],
            [0, 0, 0, 2, 3],
        ],
        dtype=float,
    )
    inverse_covariance = np.linalg.inv(covariance)
    tip_design = np.column_stack([np.ones(len(tip_x)), tip_x])
    direct_gls = np.linalg.solve(
        tip_design.T @ inverse_covariance @ tip_design,
        tip_design.T @ inverse_covariance @ tip_y,
    )
    assert result.iloc[0]["coefficient"] == pytest.approx(direct_gls[1])
    assert result.iloc[0]["n_gene_contrasts"] == 4
    assert result.iloc[0]["n_species_events"] == 4
    assert result.iloc[0]["degrees_of_freedom"] == 3
    assert result.iloc[0]["intercept"] == "no"


@pytest.mark.parametrize("reml", [True, False])
def test_diagonal_profile_fit_matches_dense_reference(reml):
    response = np.asarray([1.0, 2.5, 3.2, 5.0])
    design = np.column_stack([np.ones(4), [0.5, 1.5, 2.0, 4.0]])
    diagonal = np.asarray([0.8, 1.2, 0.7, 1.5])
    optimized = _profile_covariance_fit(
        response,
        design,
        np.zeros((4, 4)),
        [("evolutionary_rate", np.diag(diagonal))],
        reml=reml,
    )
    inverse = np.diag(1.0 / diagonal)
    expected_beta = np.linalg.solve(
        design.T @ inverse @ design,
        design.T @ inverse @ response,
    )
    residual = response - design @ expected_beta
    effective_n = len(response) - design.shape[1] if reml else len(response)
    expected_rate = float(residual @ inverse @ residual) / effective_n
    covariance = expected_rate * np.diag(diagonal)
    cholesky = np.linalg.cholesky(covariance)
    gram = design.T @ np.linalg.solve(covariance, design)
    quadratic = float(residual @ np.linalg.solve(covariance, residual))
    expected_objective = 0.5 * (
        effective_n * np.log(2.0 * np.pi)
        + 2.0 * np.log(np.diag(cholesky)).sum()
        + quadratic
        + (np.linalg.slogdet(gram)[1] if reml else 0.0)
    )

    np.testing.assert_allclose(optimized["beta"], expected_beta)
    assert optimized["component_variances"]["evolutionary_rate"] == pytest.approx(
        expected_rate
    )
    assert optimized["objective"] == pytest.approx(expected_objective)
    assert optimized["cholesky"].ndim == 1


def test_grouped_low_rank_factor_matches_dense_linear_algebra():
    diagonal = np.asarray([0.8, 1.2, 0.7, 1.5, 0.9, 1.1])
    group = np.asarray(
        [
            [1.0, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            [0.0, 0.8, 0.0],
            [0.0, 1.1, 0.0],
            [0.0, 0.0, 0.7],
            [0.0, 0.0, 1.2],
        ]
    )
    low_rank = np.asarray(
        [
            [0.2, -0.1],
            [0.3, 0.4],
            [-0.2, 0.5],
            [0.1, 0.2],
            [0.4, -0.3],
            [-0.1, 0.3],
        ]
    )
    covariance = np.diag(diagonal) + group @ group.T + low_rank @ low_rank.T
    factor = factor_diagonal_low_rank_updates(diagonal, [group, low_rank])
    values = np.column_stack([np.arange(1.0, 7.0), np.linspace(-1.0, 1.0, 6)])

    np.testing.assert_allclose(
        solve_factor(factor, values),
        np.linalg.solve(covariance, values),
        rtol=1e-12,
        atol=1e-12,
    )
    assert factor_logdet(factor) == pytest.approx(
        np.linalg.slogdet(covariance)[1], rel=1e-12, abs=1e-12
    )


def test_wide_sparse_loading_uses_observation_space_factorization():
    diagonal = np.asarray([0.8, 1.2, 0.7, 1.5])
    loading = sparse.csr_matrix(
        np.asarray(
            [
                [0.2, 0.0, 0.1, 0.0, 0.3, 0.0],
                [0.0, 0.4, 0.1, 0.0, 0.0, 0.2],
                [0.3, 0.0, 0.0, 0.5, 0.0, 0.1],
                [0.0, 0.2, 0.0, 0.1, 0.4, 0.0],
            ]
        )
    )
    covariance = np.diag(diagonal) + (loading @ loading.T).toarray()
    factor = factor_diagonal_low_rank_updates(diagonal, [loading])
    values = np.column_stack([np.arange(1.0, 5.0), np.linspace(-1.0, 1.0, 4)])

    assert isinstance(factor, SparseCovarianceFactor)
    np.testing.assert_allclose(
        solve_factor(factor, values), np.linalg.solve(covariance, values)
    )
    assert factor_logdet(factor) == pytest.approx(np.linalg.slogdet(covariance)[1])


def test_sparse_base_and_dense_update_use_nested_factorization():
    size = 40
    diagonal = np.linspace(0.8, 1.2, size)
    groups = sparse.csr_matrix(
        (
            np.ones(size),
            (np.arange(size), np.repeat(np.arange(8), 5)),
        ),
        shape=(size, 8),
    )
    dense_update = sparse.csr_matrix(
        np.column_stack(
            [
                np.sin(np.arange(size) * 0.17),
                np.cos(np.arange(size) * 0.11),
                np.linspace(-0.5, 0.5, size),
            ]
        )
    )
    covariance = (
        np.diag(diagonal)
        + (groups @ groups.T).toarray()
        + (dense_update @ dense_update.T).toarray()
    )
    factor = factor_diagonal_low_rank_updates(diagonal, [groups, dense_update])
    values = np.column_stack([np.arange(1.0, size + 1.0), np.ones(size)])

    assert isinstance(factor, NestedLowRankFactor)
    np.testing.assert_allclose(
        solve_factor(factor, values),
        np.linalg.solve(covariance, values),
        rtol=1e-11,
        atol=1e-11,
    )
    assert factor_logdet(factor) == pytest.approx(np.linalg.slogdet(covariance)[1])


def test_equal_event_weighting_is_invariant_to_identical_paralog_copies():
    base = pd.DataFrame([_response_row(1, 1.0), _response_row(2, 6.0)])
    repeated_rows = [_response_row(1, 1.0, gene_index=index) for index in range(1, 11)]
    repeated = pd.DataFrame(repeated_rows + [_response_row(2, 6.0)])
    predictors = _predictor_table(values=(1.0, 2.0))

    base_result = fit_reconciled_pgls(base, predictors, ["expression"], ["body_size"])
    equal_result = fit_reconciled_pgls(
        repeated, predictors, ["expression"], ["body_size"]
    )
    observation_result = fit_reconciled_pgls(
        repeated,
        predictors,
        ["expression"],
        ["body_size"],
        event_weighting="observation",
    )

    assert equal_result.iloc[0]["coefficient"] == pytest.approx(
        base_result.iloc[0]["coefficient"]
    )
    assert equal_result.iloc[0]["coefficient"] == pytest.approx(2.6)
    assert observation_result.iloc[0]["coefficient"] == pytest.approx(22.0 / 14.0)
    assert equal_result.iloc[0]["n_gene_contrasts"] == 11
    assert equal_result.iloc[0]["n_species_events"] == 2
    assert equal_result.iloc[0]["degrees_of_freedom"] == 1


def test_equal_event_pseudolikelihood_is_copy_invariant():
    base = pd.DataFrame(
        [
            _response_row(index, value)
            for index, value in enumerate([1.0, 6.0, 7.0, 11.0, 8.0, 13.0], start=1)
        ]
    )
    predictors = _predictor_table(values=(1.0, 2.0, 3.0, 4.0, 5.0, 6.0))
    duplicated = []
    for copy_index in range(5):
        copy = base.copy()
        copy["gene_clade_id"] = (
            copy["gene_clade_id"].astype(str) + "_" + str(copy_index)
        )
        duplicated.append(copy)

    baseline = fit_reconciled_pgls(
        base,
        predictors,
        ["expression"],
        ["body_size"],
        model="replicate-reml",
    ).iloc[0]
    repeated = fit_reconciled_pgls(
        pd.concat(duplicated, ignore_index=True),
        predictors,
        ["expression"],
        ["body_size"],
        model="replicate-reml",
    ).iloc[0]

    for column in [
        "coefficient",
        "standard_error",
        "evolutionary_rate",
        "log_likelihood",
    ]:
        assert repeated[column] == pytest.approx(baseline[column], rel=1e-10)


def test_equal_event_eiv_pseudolikelihood_is_uneven_copy_invariant():
    values = (1.0, 6.0, 7.0, 11.0, 8.0, 13.0)
    predictor = _predictor_table(values=tuple(float(i) for i in range(1, 7)))
    predictor["contrast_variance"] = 1.0
    predictor_sampling = _sampling_covariance_table(
        predictor.rename(columns={"branch_clade_id": "gene_clade_id"}),
        np.eye(6) * 0.05,
    )
    base = pd.DataFrame(
        [_response_row(index, value) for index, value in enumerate(values, start=1)]
    )
    repeated = pd.concat(
        [
            pd.DataFrame(
                [_response_row(1, values[0], gene_index=index) for index in range(5)]
            ),
            base.iloc[1:],
        ],
        ignore_index=True,
    )
    common = dict(
        predictor_contrasts=predictor,
        responses=["expression"],
        predictors=["body_size"],
        predictor_sampling_covariance=predictor_sampling,
        model="replicate-reml",
        reml=False,
        event_random_effect="no",
        lineage_random_slope="no",
    )

    baseline = fit_reconciled_pgls(response_contrasts=base, **common).iloc[0]
    duplicated = fit_reconciled_pgls(response_contrasts=repeated, **common).iloc[0]

    for column in [
        "coefficient",
        "standard_error",
        "evolutionary_rate",
        "log_likelihood",
    ]:
        assert duplicated[column] == pytest.approx(baseline[column], rel=2e-6)


def test_categorical_omnibus_row_reports_its_actual_wald_inference():
    response = pd.DataFrame(
        [
            _response_row(index, value)
            for index, value in enumerate([1.0, 4.0, 3.0, 7.0, 6.0, 9.0], start=1)
        ]
    )
    first = _predictor_table(values=(0.0, 1.0, 0.0, 1.0, 0.0, 1.0))
    first["trait"] = "habitat[a]"
    second = _predictor_table(values=(0.0, 0.0, 1.0, 1.0, 0.0, 1.0))
    second["trait"] = "habitat[b]"
    terms = ["habitat[a]", "habitat[b]"]
    metadata = {
        term: PredictorTerm(
            term,
            "habitat",
            "categorical",
            level,
            "reference",
            "treatment",
        )
        for term, level in zip(terms, ["a", "b"], strict=True)
    }

    result = fit_reconciled_pgls(
        response,
        pd.concat([first, second], ignore_index=True),
        ["expression"],
        terms,
        model="replicate-reml",
        inference="parametric-bootstrap",
        bootstrap_replicates=2,
        predictor_metadata=metadata,
        predictor_groups={"habitat": tuple(terms)},
    )

    assert set(
        result.loc[result["term_test"] == "coefficient", "inference_method"]
    ) == {"parametric-bootstrap"}
    assert (
        result.loc[result["term_test"] == "omnibus", "inference_method"].iloc[0]
        == "wald"
    )


def test_species_event_cluster_hc1_standard_error_matches_reference_formula():
    predictor_values = np.asarray([1.0, 2.0, 3.0])
    response_values = np.asarray([2.0, 3.5, 7.0])
    response = pd.DataFrame(
        [
            _response_row(index, value)
            for index, value in enumerate(response_values, start=1)
        ]
    )

    result = fit_reconciled_pgls(
        response,
        _predictor_table(values=tuple(predictor_values)),
        ["expression"],
        ["body_size"],
        model="legacy",
    ).iloc[0]

    coefficient = float(
        (predictor_values @ response_values) / (predictor_values @ predictor_values)
    )
    residuals = response_values - coefficient * predictor_values
    event_count = len(predictor_values)
    degrees_of_freedom = event_count - 1
    variance = (
        (event_count / degrees_of_freedom)
        * np.sum((predictor_values * residuals) ** 2)
        / np.sum(predictor_values**2) ** 2
    )
    assert result["coefficient"] == pytest.approx(coefficient)
    assert result["standard_error"] == pytest.approx(np.sqrt(variance))
    assert result["statistic"] == pytest.approx(coefficient / np.sqrt(variance))
    assert result["degrees_of_freedom"] == degrees_of_freedom


def test_pgls_fits_each_tree_and_response_separately():
    rows = []
    for tree_id, multiplier in [("OG1", 2.0), ("OG2", -1.0)]:
        for event_index, predictor in enumerate((1.0, 2.0, 3.0), start=1):
            rows.append(
                _response_row(
                    event_index,
                    multiplier * predictor,
                    tree_id=tree_id,
                )
            )
    result = fit_reconciled_pgls(
        pd.DataFrame(rows),
        _predictor_table(),
        ["expression"],
        ["body_size"],
    ).set_index("tree_id")

    assert result.loc["OG1", "coefficient"] == pytest.approx(2.0)
    assert result.loc["OG2", "coefficient"] == pytest.approx(-1.0)
    assert set(result["model_id"]) == {"OG1:expression", "OG2:expression"}


def test_pgls_supports_multiple_predictors():
    predictor = _predictor_table(values=(1.0, 2.0, 3.0, 4.0))
    second = predictor.copy()
    second["trait"] = "temperature"
    second["raw_contrast"] = [0.0, 1.0, 0.0, 1.0]
    predictor = pd.concat([predictor, second], ignore_index=True)
    rows = []
    for event_index, (x1, x2) in enumerate(
        zip((1.0, 2.0, 3.0, 4.0), (0.0, 1.0, 0.0, 1.0), strict=True),
        start=1,
    ):
        rows.append(_response_row(event_index, 2.0 * x1 - 3.0 * x2))

    result = fit_reconciled_pgls(
        pd.DataFrame(rows),
        predictor,
        ["expression"],
        ["body_size", "temperature"],
    ).set_index("term")

    assert result.loc["body_size", "coefficient"] == pytest.approx(2.0)
    assert result.loc["temperature", "coefficient"] == pytest.approx(-3.0)
    assert set(result["matrix_rank"]) == {2}


def test_precomputed_pgls_validates_and_reports_both_evolutionary_transforms():
    response = pd.DataFrame(
        [_response_row(1, 2.0), _response_row(2, 4.0), _response_row(3, 6.0)]
    )
    response["evolution_model"] = "kappa"
    response["evolution_parameter_name"] = "kappa"
    response["evolution_parameter"] = 0.7
    predictor = _predictor_table()
    predictor["tree_id"] = "species-tree-1"
    predictor["evolution_model"] = "delta"
    predictor["evolution_parameter_name"] = "delta"
    predictor["evolution_parameter"] = 1.4

    result = fit_reconciled_pgls(
        response,
        predictor,
        ["expression"],
        ["body_size"],
    ).iloc[0]

    assert result["predictor_tree_id"] == "species-tree-1"
    assert result["response_evolution_model"] == "kappa"
    assert result["response_evolution_parameter"] == pytest.approx(0.7)
    assert result["response_evolution_parameter_status"] == "recorded"
    assert result["response_evolution_parameter_bootstrap_refit"] == "no"
    assert result["predictor_evolution_model"] == "delta"
    assert result["predictor_evolution_parameter"] == pytest.approx(1.4)
    assert result["predictor_evolution_parameter_status"] == "recorded"


def test_precomputed_pgls_rejects_mixed_or_invalid_evolutionary_metadata():
    response = pd.DataFrame(
        [_response_row(1, 2.0), _response_row(2, 4.0), _response_row(3, 6.0)]
    )
    response.loc[0, "evolution_model"] = "kappa"
    response.loc[0, "evolution_parameter_name"] = "kappa"
    response.loc[0, "evolution_parameter"] = "0.7"
    with pytest.raises(ValueError, match="mixes evolutionary transforms"):
        fit_reconciled_pgls(
            response,
            _predictor_table(),
            ["expression"],
            ["body_size"],
        )

    predictor = _predictor_table()
    predictor.loc[0, "branch_length_mode"] = "not-applicable"
    with pytest.raises(ValueError, match="branch_length_mode is inconsistent"):
        fit_reconciled_pgls(
            pd.DataFrame(
                [
                    _response_row(1, 2.0),
                    _response_row(2, 4.0),
                    _response_row(3, 6.0),
                ]
            ),
            predictor,
            ["expression"],
            ["body_size"],
        )


def test_precomputed_pgls_rejects_predictors_combined_from_multiple_trees():
    predictor = _predictor_table()
    predictor.loc[0, "tree_id"] = "another-species-tree"
    with pytest.raises(ValueError, match="exactly one tree_id"):
        fit_reconciled_pgls(
            pd.DataFrame(
                [
                    _response_row(1, 2.0),
                    _response_row(2, 4.0),
                    _response_row(3, 6.0),
                ]
            ),
            predictor,
            ["expression"],
            ["body_size"],
        )


def test_pgls_filters_ineligible_and_partial_rows_with_reported_counts():
    rows = [
        _response_row(1, 2.0),
        _response_row(2, 3.0),
        _response_row(3, 4.0),
        _response_row(3, 5.0, gene_index=2),
    ]
    rows[2]["coverage_status"] = "partial"
    rows[3].update(
        {
            "eligible": "no",
            "coverage_status": "not-applicable",
            "species_event_id": "",
            "species_event_taxa": "",
            "species_numerator_event_id": "",
            "species_denominator_event_id": "",
        }
    )

    result = fit_reconciled_pgls(
        pd.DataFrame(rows),
        _predictor_table(),
        ["expression"],
        ["body_size"],
    )

    assert result.iloc[0]["n_species_events"] == 2
    assert result.iloc[0]["n_excluded_coverage"] == 1
    assert result.iloc[0]["n_excluded_ineligible"] == 1


def test_pgls_rejects_species_event_orientation_mismatch():
    response = pd.DataFrame([_response_row(1, 1.0), _response_row(2, 2.0)])
    response.loc[0, "species_numerator_event_id"] = "wrong"
    with pytest.raises(ValueError, match="orientation disagrees"):
        fit_reconciled_pgls(
            response,
            _predictor_table(values=(1.0, 2.0)),
            ["expression"],
            ["body_size"],
        )


def test_pgls_requires_more_species_events_than_predictors():
    response = pd.DataFrame([_response_row(1, 1.0)])
    with pytest.raises(ValueError, match="more unique species events"):
        fit_reconciled_pgls(
            response,
            _predictor_table(values=(1.0,)),
            ["expression"],
            ["body_size"],
        )


def test_pgls_rejects_inconsistent_eligibility_and_coverage():
    response = pd.DataFrame([_response_row(1, 1.0), _response_row(2, 2.0)])
    response.loc[0, "eligible"] = "no"
    with pytest.raises(ValueError, match="eligible and coverage_status"):
        fit_reconciled_pgls(
            response,
            _predictor_table(values=(1.0, 2.0)),
            ["expression"],
            ["body_size"],
        )


def test_pgls_cli_writes_coefficient_table(tmp_path):
    response_path = tmp_path / "gene-contrasts.tsv"
    predictor_path = tmp_path / "species-contrasts.tsv"
    output_path = tmp_path / "pgls.tsv"
    pd.DataFrame(
        [
            _response_row(1, 2.0),
            _response_row(2, 3.5),
            _response_row(3, 7.0),
        ]
    ).to_csv(response_path, sep="\t", index=False)
    _predictor_table().to_csv(predictor_path, sep="\t", index=False)

    assert (
        main(
            [
                "regress",
                "--infile",
                str(response_path),
                "--predictor-contrasts",
                str(predictor_path),
                "--responses",
                "expression",
                "--predictors",
                "body_size",
                "--outfile",
                str(output_path),
            ]
        )
        is None
    )
    output = pd.read_csv(output_path, sep="\t")
    assert output.iloc[0]["term"] == "body_size"
    assert output.iloc[0]["model"] == "hierarchical"
    assert output.iloc[0]["covariance_estimator"] == "gaussian-REML"
    assert output.iloc[0]["n_species_events"] == 3


def test_precomputed_pgls_output_cannot_overwrite_an_input(tmp_path):
    response_path = tmp_path / "gene-contrasts.tsv"
    predictor_path = tmp_path / "species-contrasts.tsv"
    pd.DataFrame(
        [
            _response_row(1, 2.0),
            _response_row(2, 4.0),
            _response_row(3, 6.0),
        ]
    ).to_csv(response_path, sep="\t", index=False)
    _predictor_table().to_csv(predictor_path, sep="\t", index=False)
    original = response_path.read_text()

    with pytest.raises(ValueError, match="must not overwrite an input"):
        main(
            [
                "regress",
                "--infile",
                str(response_path),
                "--predictor-contrasts",
                str(predictor_path),
                "--responses",
                "expression",
                "--predictors",
                "body_size",
                "--outfile",
                str(response_path),
            ]
        )
    assert response_path.read_text() == original


def test_replicate_reml_matches_gls_with_evolutionary_and_sampling_covariance():
    response = pd.DataFrame(
        [
            _response_row(1, 1.8),
            _response_row(2, 4.5),
            _response_row(3, 5.7),
            _response_row(4, 9.0),
        ]
    )
    response["contrast_variance"] = [1.0, 2.0, 1.5, 0.8]
    sampling = np.asarray(
        [
            [0.20, 0.05, 0.00, 0.00],
            [0.05, 0.30, 0.02, 0.00],
            [0.00, 0.02, 0.10, -0.01],
            [0.00, 0.00, -0.01, 0.25],
        ]
    )
    sidecar = _sampling_covariance_table(response, sampling)

    result = fit_reconciled_pgls(
        response,
        _predictor_table(values=(1.0, 2.0, 3.0, 4.0)),
        ["expression"],
        ["body_size"],
        model="replicate-reml",
        response_sampling_covariance=sidecar,
    ).iloc[0]

    covariance = (
        float(result["evolutionary_rate"])
        * np.diag(response["contrast_variance"].to_numpy(float))
        + sampling
    )
    inverse = np.linalg.inv(covariance)
    x = np.arange(1.0, 5.0)
    y = response["raw_contrast"].to_numpy(float)
    expected = float((x @ inverse @ y) / (x @ inverse @ x))
    assert result["coefficient"] == pytest.approx(expected)
    assert result["mean_sampling_variance"] == pytest.approx(np.diag(sampling).mean())
    assert result["event_random_effect"] == "no"
    assert result["lineage_random_slope"] == "no"


def test_replicate_response_columns_require_full_covariance_sidecar():
    response = pd.DataFrame(
        [_response_row(1, 2.0), _response_row(2, 4.0), _response_row(3, 6.0)]
    )
    response["sampling_variance"] = 0.2
    with pytest.raises(ValueError, match="full.*sampling-covariance"):
        fit_reconciled_pgls(
            response,
            _predictor_table(),
            ["expression"],
            ["body_size"],
        )


def test_sampling_covariance_must_be_complete_and_positive_semidefinite():
    response = pd.DataFrame([_response_row(1, 2.0), _response_row(2, 4.0)])
    incomplete = pd.DataFrame(
        [
            {
                "tree_id": "OG1",
                "trait": "expression",
                "contrast_id_1": "gene1_1",
                "contrast_id_2": "gene1_1",
                "sampling_covariance": 1.0,
            },
            {
                "tree_id": "OG1",
                "trait": "expression",
                "contrast_id_1": "gene2_1",
                "contrast_id_2": "gene2_1",
                "sampling_covariance": 1.0,
            },
        ]
    )
    with pytest.raises(ValueError, match="incomplete"):
        fit_reconciled_pgls(
            response,
            _predictor_table(values=(1.0, 2.0)),
            ["expression"],
            ["body_size"],
            model="replicate-reml",
            response_sampling_covariance=incomplete,
        )

    invalid = _sampling_covariance_table(response, np.asarray([[1.0, 2.0], [2.0, 1.0]]))
    with pytest.raises(ValueError, match="positive semidefinite"):
        fit_reconciled_pgls(
            response,
            _predictor_table(values=(1.0, 2.0)),
            ["expression"],
            ["body_size"],
            model="replicate-reml",
            response_sampling_covariance=invalid,
        )


def test_hierarchical_model_partially_pools_lineage_slopes():
    rows = []
    residuals = [0.1, -0.1, 0.2, -0.2, 0.1, -0.1]
    for event_index in range(1, 7):
        for gene_index, slope, sign in [(1, 1.5, 1.0), (2, 2.5, -1.0)]:
            row = _response_row(
                event_index,
                slope * event_index + sign * residuals[event_index - 1],
                gene_index=gene_index,
            )
            rows.append(row)
    result, random_effects = fit_reconciled_pgls(
        pd.DataFrame(rows),
        _predictor_table(values=tuple(float(value) for value in range(1, 7))),
        ["expression"],
        ["body_size"],
        event_random_effect="no",
        lineage_random_slope="yes",
        return_random_effects=True,
    )

    assert result.iloc[0]["coefficient"] == pytest.approx(2.0, abs=0.02)
    assert result.iloc[0]["lineage_slope_variance"] > 0.1
    assert result.iloc[0]["lineage_random_slope"] == "yes"
    assert set(random_effects["effect_type"]) == {"lineage_slope"}
    assert set(random_effects["group_id"]) == {"lineage1", "lineage2"}
    assert np.isfinite(random_effects["conditional_standard_error"]).all()
    assert np.isfinite(random_effects["total_standard_error"]).all()
    assert (
        random_effects["total_interval_lower"] <= random_effects["total_coefficient"]
    ).all()
    assert (
        random_effects["total_coefficient"] <= random_effects["total_interval_upper"]
    ).all()
    assert random_effects["reliability"].between(0.0, 1.0).all()


def test_lineage_slope_groups_are_invariant_to_predictor_units_and_basis():
    event_values = np.arange(1.0, 9.0)
    first_values = event_values
    second_values = np.asarray([1.0, -1.0, 2.0, -2.0, 3.0, -3.0, 4.0, -4.0])
    rows = []
    for event_index, (first, second) in enumerate(
        zip(first_values, second_values, strict=True), start=1
    ):
        for gene_index, first_slope, second_slope in [
            (1, 1.2, 0.4),
            (2, 2.0, -0.3),
        ]:
            rows.append(
                _response_row(
                    event_index,
                    first_slope * first
                    + second_slope * second
                    + 0.05 * (-1.0) ** (event_index + gene_index),
                    gene_index=gene_index,
                )
            )

    def predictor_table(first, second, first_name, second_name):
        first_table = _predictor_table(values=tuple(first))
        first_table["trait"] = first_name
        second_table = _predictor_table(values=tuple(second))
        second_table["trait"] = second_name
        return pd.concat([first_table, second_table], ignore_index=True)

    common = dict(
        response_contrasts=pd.DataFrame(rows),
        responses=["expression"],
        event_random_effect="no",
        lineage_random_slope="yes",
        reml=False,
    )
    original = fit_reconciled_pgls(
        predictor_contrasts=predictor_table(
            first_values, second_values, "first", "second"
        ),
        predictors=["first", "second"],
        predictor_groups={"basis": ("first", "second")},
        **common,
    )
    rescaled = fit_reconciled_pgls(
        predictor_contrasts=predictor_table(
            first_values, 100.0 * second_values, "first", "second"
        ),
        predictors=["first", "second"],
        predictor_groups={"basis": ("first", "second")},
        **common,
    )
    transformed = fit_reconciled_pgls(
        predictor_contrasts=predictor_table(
            first_values + second_values,
            first_values - second_values,
            "sum",
            "difference",
        ),
        predictors=["sum", "difference"],
        predictor_groups={"basis": ("sum", "difference")},
        **common,
    )
    separate_sources = fit_reconciled_pgls(
        predictor_contrasts=predictor_table(
            first_values, second_values, "first", "second"
        ),
        predictors=["first", "second"],
        predictor_groups={"first": ("first",), "second": ("second",)},
        **common,
    )

    by_term = original.set_index("term")
    scaled_by_term = rescaled.set_index("term")
    assert scaled_by_term.loc["first", "coefficient"] == pytest.approx(
        by_term.loc["first", "coefficient"], rel=2e-5
    )
    assert scaled_by_term.loc["first", "standard_error"] == pytest.approx(
        by_term.loc["first", "standard_error"], rel=2e-4
    )
    assert scaled_by_term.loc["second", "coefficient"] * 100.0 == pytest.approx(
        by_term.loc["second", "coefficient"], rel=2e-5
    )
    assert scaled_by_term.loc["second", "standard_error"] * 100.0 == pytest.approx(
        by_term.loc["second", "standard_error"], rel=2e-4
    )
    assert rescaled.iloc[0]["log_likelihood"] == pytest.approx(
        original.iloc[0]["log_likelihood"], rel=2e-5
    )
    assert transformed.iloc[0]["log_likelihood"] == pytest.approx(
        original.iloc[0]["log_likelihood"], rel=2e-5
    )
    separate_by_term = separate_sources.set_index("term")
    assert separate_by_term.loc["first", "lineage_slope_variance"] != pytest.approx(
        separate_by_term.loc["second", "lineage_slope_variance"], rel=1e-3
    )


def test_grouped_lineage_joint_test_uses_fixed_terms_plus_its_variance_df():
    rows = []
    first_values = np.arange(1.0, 9.0)
    second_values = np.asarray([1.0, -1.0, 2.0, -2.0, 3.0, -3.0, 4.0, -4.0])
    for event_index, (first, second) in enumerate(
        zip(first_values, second_values, strict=True), start=1
    ):
        for gene_index, shift in [(1, -0.4), (2, 0.4)]:
            rows.append(
                _response_row(
                    event_index,
                    (1.5 + shift) * first + (0.2 - shift) * second,
                    gene_index=gene_index,
                )
            )
    first = _predictor_table(values=tuple(first_values))
    first["trait"] = "state[a]"
    second = _predictor_table(values=tuple(second_values))
    second["trait"] = "state[b]"
    terms = ["state[a]", "state[b]"]

    result = fit_reconciled_pgls(
        pd.DataFrame(rows),
        pd.concat([first, second], ignore_index=True),
        ["expression"],
        terms,
        predictor_groups={"state": tuple(terms)},
        event_random_effect="no",
        lineage_random_slope="yes",
        lineage_inference="likelihood-ratio",
    )

    tests = result.set_index("term_test")
    assert tests.loc["lineage-heterogeneity", "degrees_of_freedom"] == 1
    assert tests.loc["average-and-lineage-joint", "degrees_of_freedom"] == 3


def test_lineage_random_intervals_include_predictor_dependent_eiv_covariance():
    rows = []
    for event_index in range(1, 7):
        for gene_index, slope in [(1, 1.4), (2, 2.6)]:
            rows.append(
                _response_row(
                    event_index,
                    slope * event_index + 0.03 * (-1.0) ** event_index,
                    gene_index=gene_index,
                )
            )
    predictor = _predictor_table(values=tuple(float(value) for value in range(1, 7)))
    predictor["contrast_variance"] = 1.0
    predictor_sampling = _sampling_covariance_table(
        predictor.rename(columns={"branch_clade_id": "gene_clade_id"}),
        np.eye(6) * 0.05,
    )

    _, random_effects = fit_reconciled_pgls(
        pd.DataFrame(rows),
        predictor,
        ["expression"],
        ["body_size"],
        predictor_sampling_covariance=predictor_sampling,
        event_random_effect="no",
        lineage_random_slope="yes",
        return_random_effects=True,
    )

    assert set(random_effects["inference_status"]) == {
        "empirical-bayes-delta-eiv-conditional-on-variance"
    }
    assert np.isfinite(random_effects["conditional_standard_error"]).all()
    assert np.isfinite(random_effects["total_standard_error"]).all()


def test_lineage_inference_and_leave_one_out_separate_average_and_subset_effects():
    rows = []
    for event_index in range(1, 7):
        for gene_index, slope in [(1, 1.5), (2, 2.5)]:
            rows.append(
                _response_row(
                    event_index,
                    slope * event_index + (0.05 if gene_index == 1 else -0.05),
                    gene_index=gene_index,
                )
            )
    result, sensitivity = fit_reconciled_pgls(
        pd.DataFrame(rows),
        _predictor_table(values=tuple(float(value) for value in range(1, 7))),
        ["expression"],
        ["body_size"],
        event_random_effect="no",
        lineage_random_slope="yes",
        lineage_inference="likelihood-ratio",
        lineage_leave_one_out=True,
        return_sensitivity=True,
    )

    by_test = result.set_index("term_test")
    assert by_test.loc["lineage-heterogeneity", "statistic"] > 0.0
    assert 0.0 <= by_test.loc["lineage-heterogeneity", "p_value"] <= 1.0
    assert (
        by_test.loc["average-and-lineage-joint", "inference_status"]
        == "parametric-bootstrap-required-for-joint-null"
    )
    assert set(sensitivity["group_id"]) == {"lineage1", "lineage2"}
    assert set(sensitivity["inference_status"]) == {"ok"}
    assert set(np.sign(sensitivity["coefficient_change"])) == {-1.0, 1.0}


def test_lineage_joint_parametric_bootstrap_reports_calibrated_p_values():
    rows = []
    for event_index in range(1, 7):
        for gene_index, slope in [(1, 1.5), (2, 2.5)]:
            rows.append(
                _response_row(
                    event_index,
                    slope * event_index + 0.01 * gene_index,
                    gene_index=gene_index,
                )
            )
    result = fit_reconciled_pgls(
        pd.DataFrame(rows),
        _predictor_table(values=tuple(float(value) for value in range(1, 7))),
        ["expression"],
        ["body_size"],
        event_random_effect="no",
        lineage_random_slope="yes",
        lineage_inference="parametric-bootstrap",
        bootstrap_replicates=2,
        seed=4,
    )

    tests = result[result["term_test"].str.contains("lineage")]
    assert set(tests["inference_status"]) == {"ok"}
    assert tests["p_value"].between(0.0, 1.0).all()


def test_hierarchical_model_estimates_shared_species_event_effects():
    offsets = [-1.0, 1.0, -0.7, 0.7, -1.2, 1.2]
    rows = []
    for event_index, offset in enumerate(offsets, start=1):
        for copy, residual in [(1, 0.05), (2, -0.05)]:
            row = _response_row(event_index, 2.0 * event_index + offset + residual)
            row["gene_clade_id"] = "gene{}_{}".format(event_index, copy)
            row["lineage_clade_id"] = "single{}_{}".format(event_index, copy)
            rows.append(row)
    result, random_effects = fit_reconciled_pgls(
        pd.DataFrame(rows),
        _predictor_table(values=tuple(float(value) for value in range(1, 7))),
        ["expression"],
        ["body_size"],
        event_random_effect="yes",
        lineage_random_slope="no",
        return_random_effects=True,
    )

    assert result.iloc[0]["species_event_variance"] > 0.1
    assert result.iloc[0]["event_random_effect"] == "yes"
    assert set(random_effects["effect_type"]) == {"species_event"}
    assert len(random_effects) == 6
    assert set(random_effects["n_observations"]) == {2}


def test_parametric_bootstrap_is_reproducible_and_reports_empirical_inference():
    response = pd.DataFrame(
        [
            _response_row(1, 1.8),
            _response_row(2, 4.3),
            _response_row(3, 5.6),
            _response_row(4, 8.5),
        ]
    )
    predictor = _predictor_table(values=(1.0, 2.0, 3.0, 4.0))
    predictor["contrast_variance"] = [1.0, 1.0, 1.0, 1.0]
    predictor_sampling = _sampling_covariance_table(
        predictor.rename(columns={"branch_clade_id": "gene_clade_id"}),
        np.eye(4) * 0.05,
    )
    arguments = dict(
        response_contrasts=response,
        predictor_contrasts=predictor,
        responses=["expression"],
        predictors=["body_size"],
        predictor_sampling_covariance=predictor_sampling,
        model="replicate-reml",
        inference="parametric-bootstrap",
        bootstrap_replicates=12,
        seed=19,
    )
    first = fit_reconciled_pgls(**arguments).iloc[0]
    second = fit_reconciled_pgls(**arguments).iloc[0]

    for column in [
        "coefficient",
        "standard_error",
        "p_value",
        "confidence_interval_lower",
        "confidence_interval_upper",
    ]:
        assert first[column] == pytest.approx(second[column])
    assert first["inference_method"] == "parametric-bootstrap"
    assert first["measurement_error_model"] == "latent-predictor"
    assert 0.0 < first["p_value"] <= 1.0


def test_model_specific_options_are_validated_instead_of_ignored():
    response = pd.DataFrame(
        [_response_row(1, 2.0), _response_row(2, 4.0), _response_row(3, 6.0)]
    )
    with pytest.raises(ValueError, match="unavailable for legacy"):
        fit_reconciled_pgls(
            response,
            _predictor_table(),
            ["expression"],
            ["body_size"],
            model="legacy",
            inference="parametric-bootstrap",
            bootstrap_replicates=2,
        )
    with pytest.raises(ValueError, match="require the hierarchical"):
        fit_reconciled_pgls(
            response,
            _predictor_table(),
            ["expression"],
            ["body_size"],
            model="replicate-reml",
            event_random_effect="yes",
        )
    with pytest.raises(ValueError, match="requires lineage random slopes"):
        fit_reconciled_pgls(
            response,
            _predictor_table(),
            ["expression"],
            ["body_size"],
            lineage_random_slope="no",
            lineage_inference="likelihood-ratio",
        )
    with pytest.raises(ValueError, match="sequences of names"):
        fit_reconciled_pgls(
            response,
            _predictor_table(),
            "expression",
            ["body_size"],
        )


def test_pgls_rejects_empty_tree_ids_to_prevent_family_pooling():
    response = pd.DataFrame(
        [_response_row(1, 2.0), _response_row(2, 4.0), _response_row(3, 6.0)]
    )
    response["tree_id"] = ""
    with pytest.raises(ValueError, match="non-empty tree_id"):
        fit_reconciled_pgls(
            response,
            _predictor_table(),
            ["expression"],
            ["body_size"],
        )


@pytest.mark.integration
def test_pgls_raw_mode_writes_complete_replicate_aware_bundle_and_audit(tmp_path):
    gene_tree, species_tree, expression, species_traits = _write_raw_regression_inputs(
        tmp_path,
        biological_replicates=True,
    )
    prefix = tmp_path / "analysis"
    audit = tmp_path / "analysis.audit.jsonl"

    assert (
        main(
            [
                "regress",
                "--gene-tree",
                str(gene_tree),
                "--species-tree",
                str(species_tree),
                "--expression",
                str(expression),
                "--species-traits",
                str(species_traits),
                "--responses",
                "expression",
                "--predictors",
                "body_size",
                "--tree-id",
                "OG000001",
                "--biological-id",
                "sample_id",
                "--out-prefix",
                str(prefix),
                "--audit",
                str(audit),
            ]
        )
        is None
    )

    expected_paths = {
        tmp_path / "analysis.reconciliation.tsv",
        tmp_path / "analysis.gene-contrasts.tsv",
        tmp_path / "analysis.species-contrasts.tsv",
        tmp_path / "analysis.response-sampling-covariance.tsv",
        tmp_path / "analysis.response-tip-summary.tsv",
        tmp_path / "analysis.random-effects.tsv",
        tmp_path / "analysis.regression.tsv",
    }
    assert all(path.is_file() for path in expected_paths)
    reconciliation = pd.read_csv(tmp_path / "analysis.reconciliation.tsv", sep="\t")
    gene_contrasts = pd.read_csv(tmp_path / "analysis.gene-contrasts.tsv", sep="\t")
    covariance = pd.read_csv(
        tmp_path / "analysis.response-sampling-covariance.tsv", sep="\t"
    )
    tip_summary = pd.read_csv(tmp_path / "analysis.response-tip-summary.tsv", sep="\t")
    result = pd.read_csv(tmp_path / "analysis.regression.tsv", sep="\t")
    assert set(reconciliation["tree_id"]) == {"OG000001"}
    assert set(gene_contrasts["tree_id"]) == {"OG000001"}
    assert set(gene_contrasts["replicate_model"]) == {"pooled"}
    assert len(covariance) == 10
    assert set(tip_summary["n_biological"]) == {2}
    assert result.iloc[0]["term"] == "body_size"
    assert result.iloc[0]["response"] == "expression"
    assert result.iloc[0]["model"] == "hierarchical"
    assert result.iloc[0]["response_evolution_parameter_status"] == "not-applicable"
    assert result.iloc[0]["response_evolution_optimizer_converged"] == "not-applicable"
    record = json.loads(audit.read_text())
    assert {item["path"] for item in record["inputs"]} == {
        str(path.resolve())
        for path in [gene_tree, species_tree, expression, species_traits]
    }
    assert {item["path"] for item in record["outputs"]} == {
        str(path.resolve()) for path in expected_paths
    }


@pytest.mark.integration
def test_pgls_raw_mode_supports_mixed_response_replication_depth(tmp_path):
    gene_tree, species_tree, expression, species_traits = _write_raw_regression_inputs(
        tmp_path,
        biological_replicates=True,
    )
    frame = pd.read_csv(expression, sep="\t")
    frame["single"] = np.where(
        frame["sample_id"].str.endswith("_1"),
        frame["expression"] + 1.0,
        np.nan,
    )
    frame.to_csv(expression, sep="\t", index=False, na_rep="NA")
    prefix = tmp_path / "mixed-depth"

    main(
        [
            "regress",
            "--gene-tree",
            str(gene_tree),
            "--species-tree",
            str(species_tree),
            "--expression",
            str(expression),
            "--species-traits",
            str(species_traits),
            "--responses",
            "expression,single",
            "--predictors",
            "body_size",
            "--tree-id",
            "OGMIXED",
            "--biological-id",
            "sample_id",
            "--out-prefix",
            str(prefix),
        ]
    )

    result = pd.read_csv(tmp_path / "mixed-depth.regression.tsv", sep="\t")
    assert set(result["response"]) == {"expression", "single"}
    summary = pd.read_csv(tmp_path / "mixed-depth.response-tip-summary.tsv", sep="\t")
    assert set(summary.query("trait == 'expression'")["variance_method"]) == {"pooled"}
    assert set(summary.query("trait == 'single'")["variance_method"]) == {
        "single-observation"
    }
    assert set(summary.query("trait == 'single'")["standard_error"]) == {0.0}


@pytest.mark.integration
def test_reconciled_pgls_accepts_categorical_species_predictor(tmp_path):
    gene_tree, species_tree, expression, species_traits = _write_raw_regression_inputs(
        tmp_path
    )
    traits = pd.read_csv(species_traits, sep="\t")
    traits["habitat"] = [
        "aquatic",
        "terrestrial",
        "arboreal",
        "terrestrial",
        "arboreal",
    ]
    traits.to_csv(species_traits, sep="\t", index=False)
    output = tmp_path / "categorical-pgls.tsv"

    main(
        [
            "regress",
            "--gene-tree",
            str(gene_tree),
            "--species-tree",
            str(species_tree),
            "--expression",
            str(expression),
            "--species-traits",
            str(species_traits),
            "--responses",
            "expression",
            "--predictors",
            "habitat",
            "--categorical-predictors",
            "habitat",
            "--factor-reference",
            "habitat=aquatic",
            "--tree-id",
            "OGCAT",
            "--outfile",
            str(output),
        ]
    )

    result = pd.read_csv(output, sep="\t")
    coefficients = result[result["term_test"] == "coefficient"]
    assert set(coefficients["term"]) == {
        "habitat[arboreal]",
        "habitat[terrestrial]",
    }
    assert set(coefficients["source_term"]) == {"habitat"}
    assert set(coefficients["predictor_type"]) == {"categorical"}
    omnibus = result[result["term_test"] == "omnibus"].iloc[0]
    assert omnibus["term"] == "habitat"
    assert omnibus["degrees_of_freedom"] == 2


@pytest.mark.integration
def test_categorical_origin_mapping_and_origin_leave_one_out_are_auditable(tmp_path):
    gene_tree, species_tree, expression, species_traits = _write_raw_regression_inputs(
        tmp_path
    )
    traits = pd.read_csv(species_traits, sep="\t")
    traits["habitat"] = ["land", "land", "water", "water", "water"]
    traits.to_csv(species_traits, sep="\t", index=False)
    result_path = tmp_path / "origin-pgls.tsv"
    origin_path = tmp_path / "origins.tsv"
    sensitivity_path = tmp_path / "origin-sensitivity.tsv"

    main(
        [
            "regress",
            "--gene-tree",
            str(gene_tree),
            "--species-tree",
            str(species_tree),
            "--expression",
            str(expression),
            "--species-traits",
            str(species_traits),
            "--responses",
            "expression",
            "--predictors",
            "habitat",
            "--categorical-predictors",
            "habitat",
            "--tree-id",
            "OGORIGIN",
            "--categorical-origin-diagnostics",
            "stochastic-map",
            "--origin-map-replicates",
            "10",
            "--origin-min-posterior",
            "0.1",
            "--origin-leave-one-out",
            "yes",
            "--trait-origins-out",
            str(origin_path),
            "--sensitivity-out",
            str(sensitivity_path),
            "--outfile",
            str(result_path),
        ]
    )

    origins = pd.read_csv(origin_path, sep="\t")
    assert set(origins["trait"]) == {"habitat"}
    assert set(origins["mk_model"]) == {"ER"}
    assert set(origins["from_state"]) == {"land", "water"}
    assert origins["posterior_frequency"].between(0.0, 1.0).all()
    sensitivity = pd.read_csv(sensitivity_path, sep="\t")
    assert set(sensitivity["analysis_type"]) == {"trait-origin-leave-one-out"}
    assert (sensitivity["n_omitted_gene_contrasts"] > 0).all()
    assert set(sensitivity["term"]) == {"habitat[water]"}


def test_origin_specific_options_require_stochastic_mapping_mode(tmp_path):
    gene_tree, species_tree, expression, species_traits = _write_raw_regression_inputs(
        tmp_path
    )
    with pytest.raises(ValueError, match="Categorical origin option.*require"):
        main(
            [
                "regress",
                "--gene-tree",
                str(gene_tree),
                "--species-tree",
                str(species_tree),
                "--expression",
                str(expression),
                "--species-traits",
                str(species_traits),
                "--responses",
                "expression",
                "--predictors",
                "body_size",
                "--tree-id",
                "OGINVALID",
                "--origin-map-replicates",
                "10",
            ]
        )


@pytest.mark.integration
def test_reconciled_pgls_propagates_latent_categorical_predictor_replicates(tmp_path):
    gene_tree, species_tree, expression, species_traits = _write_raw_regression_inputs(
        tmp_path
    )
    traits = pd.read_csv(species_traits, sep="\t")
    rows = []
    states = ["water", "land", "water", "land", "water"]
    for leaf_name, state in zip(traits["leaf_name"], states, strict=True):
        rows.append(
            {
                "leaf_name": leaf_name,
                "sample": "{}_1".format(leaf_name),
                "habitat": state,
            }
        )
        rows.append(
            {
                "leaf_name": leaf_name,
                "sample": "{}_2".format(leaf_name),
                "habitat": "land"
                if leaf_name == traits.iloc[0]["leaf_name"]
                else state,
            }
        )
    pd.DataFrame(rows).to_csv(species_traits, sep="\t", index=False)
    prefix = tmp_path / "latent-category"

    main(
        [
            "regress",
            "--gene-tree",
            str(gene_tree),
            "--species-tree",
            str(species_tree),
            "--expression",
            str(expression),
            "--species-traits",
            str(species_traits),
            "--responses",
            "expression",
            "--predictors",
            "habitat",
            "--categorical-predictors",
            "habitat",
            "--predictor-biological-id",
            "sample",
            "--categorical-replicate-policy",
            "latent",
            "--categorical-origin-diagnostics",
            "stochastic-map",
            "--origin-map-replicates",
            "10",
            "--tree-id",
            "OGLATENT",
            "--out-prefix",
            str(prefix),
        ]
    )

    result = pd.read_csv(tmp_path / "latent-category.regression.tsv", sep="\t")
    coefficient = result[result["term_test"] == "coefficient"].iloc[0]
    assert coefficient["measurement_error_model"] == "latent-predictor"
    summary = pd.read_csv(
        tmp_path / "latent-category.predictor-tip-summary.tsv", sep="\t"
    )
    discordant = summary[summary["state"].isna()].iloc[0]
    assert "land:1" in discordant["state_counts"]
    assert "water:1" in discordant["state_counts"]
    origins = pd.read_csv(tmp_path / "latent-category.trait-origins.tsv", sep="\t")
    assert set(origins["trait"]) == {"habitat"}
    assert origins["posterior_frequency"].between(0.0, 1.0).all()


@pytest.mark.integration
def test_reconciled_multilevel_factor_preserves_cross_column_uncertainty(tmp_path):
    gene_tree, species_tree, expression, species_traits = _write_raw_regression_inputs(
        tmp_path
    )
    species_names = pd.read_csv(species_traits, sep="\t")["leaf_name"].tolist()
    replicate_states = [
        ("water", "land"),
        ("land", "air"),
        ("water", "water"),
        ("air", "air"),
        ("land", "land"),
    ]
    rows = []
    for leaf, states in zip(species_names, replicate_states, strict=True):
        for replicate, state in enumerate(states, start=1):
            rows.append(
                {
                    "leaf_name": leaf,
                    "sample": "{}_{}".format(leaf, replicate),
                    "habitat": state,
                }
            )
    pd.DataFrame(rows).to_csv(species_traits, sep="\t", index=False)
    prefix = tmp_path / "multilevel-latent"

    main(
        [
            "regress",
            "--gene-tree",
            str(gene_tree),
            "--species-tree",
            str(species_tree),
            "--expression",
            str(expression),
            "--species-traits",
            str(species_traits),
            "--responses",
            "expression",
            "--predictors",
            "habitat",
            "--categorical-predictors",
            "habitat",
            "--predictor-biological-id",
            "sample",
            "--categorical-replicate-policy",
            "latent",
            "--species-evolution-model",
            "lambda",
            "--tree-id",
            "OGMULTI",
            "--out-prefix",
            str(prefix),
        ]
    )

    result = pd.read_csv(tmp_path / "multilevel-latent.regression.tsv", sep="\t")
    coefficients = result[result["term_test"] == "coefficient"]
    assert len(coefficients) == 2
    assert set(coefficients["measurement_error_model"]) == {"latent-predictor"}
    assert coefficients["predictor_evolution_parameter"].nunique() == 1
    assert set(coefficients["covariance_estimator"]) == {"gaussian-eiv-ML"}
    covariance = pd.read_csv(
        tmp_path / "multilevel-latent.predictor-sampling-covariance.tsv", sep="\t"
    )
    assert "trait_2" in covariance
    cross_terms = covariance[covariance["trait"] != covariance["trait_2"]]
    assert not cross_terms.empty
    assert (cross_terms["sampling_covariance"].abs() > 0.0).any()


@pytest.mark.integration
def test_reconciled_categorical_response_uses_tip_pglmm(tmp_path):
    gene_tree, species_tree, expression, species_traits = _write_raw_regression_inputs(
        tmp_path
    )
    expression_table = pd.read_csv(expression, sep="\t")
    expression_table["state"] = ["low", "middle", "high", "low", "high"]
    expression_table.to_csv(expression, sep="\t", index=False)
    traits = pd.read_csv(species_traits, sep="\t")
    traits["habitat"] = ["water", "land", "water", "land", "water"]
    traits.to_csv(species_traits, sep="\t", index=False)
    prefix = tmp_path / "categorical-response"

    main(
        [
            "regress",
            "--gene-tree",
            str(gene_tree),
            "--species-tree",
            str(species_tree),
            "--expression",
            str(expression),
            "--species-traits",
            str(species_traits),
            "--responses",
            "state",
            "--ordered-responses",
            "state=low|middle|high",
            "--predictors",
            "habitat",
            "--categorical-predictors",
            "habitat",
            "--gene-evolution-model",
            "lambda",
            "--tree-id",
            "OGSTATE",
            "--out-prefix",
            str(prefix),
        ]
    )

    result = pd.read_csv(tmp_path / "categorical-response.regression.tsv", sep="\t")
    coefficients = result[result["term_test"] == "coefficient"]
    assert set(coefficients["response_family"]) == {"ordinal"}
    assert set(coefficients["model"]) == {"reconciled-tip-pglmm"}
    assert set(coefficients["contrast_transform"]) == {"not-applicable-tip-pglmm"}
    assert set(coefficients["response_evolution_parameter_status"]) == {"estimated"}
    assert coefficients["response_evolution_parameter"].between(0.0, 1.0).all()
    assert len(result[result["term_test"] == "threshold"]) == 2
    gene_contrasts = pd.read_csv(
        tmp_path / "categorical-response.gene-contrasts.tsv", sep="\t"
    )
    assert gene_contrasts.empty
    random_effects = pd.read_csv(
        tmp_path / "categorical-response.random-effects.tsv", sep="\t"
    )
    assert set(random_effects["effect_type"]) == {"phylogenetic_tip_logit"}


@pytest.mark.integration
def test_reconciled_categorical_response_and_predictor_replicates_together(tmp_path):
    gene_tree, species_tree, expression, species_traits = _write_raw_regression_inputs(
        tmp_path
    )
    gene_names = pd.read_csv(expression, sep="\t")["leaf_name"].tolist()
    response_states = [
        ("absent", "present"),
        ("absent", "absent"),
        ("present", "absent"),
        ("present", "present"),
        ("absent", "present"),
    ]
    response_rows = [
        {
            "leaf_name": leaf,
            "sample": "{}_{}".format(leaf, replicate),
            "state": state,
        }
        for leaf, states in zip(gene_names, response_states, strict=True)
        for replicate, state in enumerate(states, start=1)
    ]
    pd.DataFrame(response_rows).to_csv(expression, sep="\t", index=False)

    species_names = pd.read_csv(species_traits, sep="\t")["leaf_name"].tolist()
    predictor_states = [
        ("water", "land"),
        ("land", "land"),
        ("water", "water"),
        ("land", "water"),
        ("water", "water"),
    ]
    predictor_rows = [
        {
            "leaf_name": leaf,
            "sample": "{}_{}".format(leaf, replicate),
            "habitat": state,
        }
        for leaf, states in zip(species_names, predictor_states, strict=True)
        for replicate, state in enumerate(states, start=1)
    ]
    pd.DataFrame(predictor_rows).to_csv(species_traits, sep="\t", index=False)
    prefix = tmp_path / "categorical-replicates"

    main(
        [
            "regress",
            "--gene-tree",
            str(gene_tree),
            "--species-tree",
            str(species_tree),
            "--expression",
            str(expression),
            "--species-traits",
            str(species_traits),
            "--responses",
            "state",
            "--predictors",
            "habitat",
            "--biological-id",
            "sample",
            "--predictor-biological-id",
            "sample",
            "--categorical-replicate-policy",
            "latent",
            "--tree-id",
            "OGCATREP",
            "--out-prefix",
            str(prefix),
        ]
    )

    result = pd.read_csv(tmp_path / "categorical-replicates.regression.tsv", sep="\t")
    coefficients = result[result["term_test"] == "coefficient"]
    assert set(coefficients["response_family"]) == {"binomial"}
    assert set(coefficients["predictor_type"]) == {"categorical", "intercept"}
    assert set(coefficients["measurement_error_model"]) == {"latent-predictor"}
    response_summary = pd.read_csv(
        tmp_path / "categorical-replicates.response-tip-summary.tsv", sep="\t"
    )
    predictor_summary = pd.read_csv(
        tmp_path / "categorical-replicates.predictor-tip-summary.tsv", sep="\t"
    )
    assert set(response_summary["variance_method"]) == {"categorical-counts"}
    assert "categorical-latent" in set(predictor_summary["variance_method"])


@pytest.mark.integration
def test_reconciled_negative_binomial_keeps_biological_count_replicates(tmp_path):
    gene_tree, species_tree, expression, species_traits = _write_raw_regression_inputs(
        tmp_path
    )
    gene_names = pd.read_csv(expression, sep="\t")["leaf_name"].tolist()
    rows = []
    for index, leaf_name in enumerate(gene_names, start=1):
        for replicate, count in enumerate((index, index + 2), start=1):
            rows.append(
                {
                    "leaf_name": leaf_name,
                    "sample": "{}_{}".format(leaf_name, replicate),
                    "count": count,
                    "log_library": np.log(1000.0 + 100.0 * replicate),
                }
            )
    pd.DataFrame(rows).to_csv(expression, sep="\t", index=False)
    prefix = tmp_path / "negative-binomial"

    main(
        [
            "regress",
            "--gene-tree",
            str(gene_tree),
            "--species-tree",
            str(species_tree),
            "--expression",
            str(expression),
            "--species-traits",
            str(species_traits),
            "--responses",
            "count",
            "--response-family",
            "count=negative-binomial",
            "--response-offset",
            "count=log_library",
            "--predictors",
            "body_size",
            "--biological-id",
            "sample",
            "--inference",
            "parametric-bootstrap",
            "--bootstrap-replicates",
            "2",
            "--seed",
            "9",
            "--tree-id",
            "OGNB",
            "--out-prefix",
            str(prefix),
        ]
    )

    result = pd.read_csv(tmp_path / "negative-binomial.regression.tsv", sep="\t")
    assert set(result["response_family"]) == {"negative-binomial"}
    assert set(result["link_function"]) == {"log"}
    assert set(result["inference_method"]) == {"parametric-bootstrap"}
    assert np.isfinite(result["response_dispersion"]).all()
    replicate_summary = pd.read_csv(
        tmp_path / "negative-binomial.response-tip-summary.tsv", sep="\t"
    )
    assert set(replicate_summary["variance_method"]) == {"likelihood-replicates"}
    assert set(replicate_summary["n_biological"]) == {2}


@pytest.mark.integration
def test_reconciled_gene_tree_ensemble_combines_tree_uncertainty(tmp_path):
    gene_tree, species_tree, expression, species_traits = _write_raw_regression_inputs(
        tmp_path
    )
    original = gene_tree.read_text()
    alternate = original.replace(
        "Genus_a_g1:1,Genus_b_g1:1", "Genus_b_g1:1,Genus_a_g1:1"
    )
    ensemble = tmp_path / "gene-ensemble.nwk"
    ensemble.write_text(original + "\n" + alternate + "\n")
    output = tmp_path / "ensemble.tsv"
    audit = tmp_path / "ensemble.audit.jsonl"

    main(
        [
            "regress",
            "--gene-tree-ensemble",
            str(ensemble),
            "--species-tree",
            str(species_tree),
            "--expression",
            str(expression),
            "--species-traits",
            str(species_traits),
            "--responses",
            "expression",
            "--predictors",
            "body_size",
            "--tree-id",
            "OGENS",
            "--outfile",
            str(output),
            "--audit",
            str(audit),
        ]
    )

    result = pd.read_csv(output, sep="\t")
    coefficients = result[result["term_test"] == "coefficient"]
    assert set(coefficients["ensemble_size"]) == {2}
    assert set(coefficients["tree_support_fraction"]) == {1.0}
    assert np.isfinite(coefficients["between_tree_variance"]).all()
    assert set(coefficients["inference_method"]) == {"tree-ensemble-rubin"}
    audit_record = json.loads(audit.read_text())
    assert audit_record["primary_input"]["tree_count"] == 2
    assert any(
        record["argument"] == "gene_tree_ensemble" for record in audit_record["inputs"]
    )


@pytest.mark.integration
def test_reconciled_multivariate_pgls_retains_missing_paralog_responses(tmp_path):
    gene_tree, species_tree, expression, species_traits = _write_raw_regression_inputs(
        tmp_path
    )
    gene_names = [
        "Genus_a_g1",
        "Genus_a_g2",
        "Genus_b_g1",
        "Genus_c_g1",
        "Genus_d_g1",
        "Genus_e_g1",
    ]
    gene_tree.write_text(
        "((((Genus_a_g1:1,Genus_a_g2:1):1,Genus_b_g1:2):1,"
        "Genus_c_g1:3):1,(Genus_d_g1:1,Genus_e_g1:1):3);"
    )
    pd.DataFrame(
        {
            "leaf_name": gene_names,
            "expression": [2.0, 2.5, 5.0, 7.0, 8.0, 12.0],
            "expression_alt": [3.0, 3.4, 4.5, np.nan, 9.0, 11.0],
        }
    ).to_csv(expression, sep="\t", index=False)
    output = tmp_path / "multivariate.tsv"

    main(
        [
            "regress",
            "--gene-tree",
            str(gene_tree),
            "--species-tree",
            str(species_tree),
            "--expression",
            str(expression),
            "--species-traits",
            str(species_traits),
            "--responses",
            "expression,expression_alt",
            "--predictors",
            "body_size",
            "--multivariate-responses",
            "yes",
            "--allow-missing-responses",
            "yes",
            "--tree-id",
            "OGMULTIVAR",
            "--outfile",
            str(output),
        ]
    )

    result = pd.read_csv(output, sep="\t")
    coefficients = result[result["term_test"] == "coefficient"]
    assert set(coefficients["response"]) == {"expression", "expression_alt"}
    assert set(coefficients["model"]) == {"reconciled-tip-multivariate-pgls"}
    assert set(coefficients["event_random_effect"]) == {"species-tip"}
    covariance_rows = result[result["term_test"] == "response-covariance"]
    assert set(covariance_rows["source_term"]) == {"(response-covariance)"}
    assert len(covariance_rows) == 6


@pytest.mark.integration
def test_reconciled_categorical_response_shares_species_effect_across_paralogs(
    tmp_path,
):
    gene_tree, species_tree, expression, species_traits = _write_raw_regression_inputs(
        tmp_path
    )
    gene_names = [
        "Genus_a_g1",
        "Genus_a_g2",
        "Genus_b_g1",
        "Genus_c_g1",
        "Genus_d_g1",
        "Genus_e_g1",
    ]
    gene_tree.write_text(
        "((((Genus_a_g1:1,Genus_a_g2:1):1,Genus_b_g1:2):1,"
        "Genus_c_g1:3):1,(Genus_d_g1:1,Genus_e_g1:1):3);"
    )
    pd.DataFrame(
        {
            "leaf_name": gene_names,
            "state": ["absent", "present", "absent", "present", "present", "absent"],
        }
    ).to_csv(expression, sep="\t", index=False)
    output = tmp_path / "paralog-categorical.tsv"
    random_output = tmp_path / "paralog-random.tsv"

    main(
        [
            "regress",
            "--gene-tree",
            str(gene_tree),
            "--species-tree",
            str(species_tree),
            "--expression",
            str(expression),
            "--species-traits",
            str(species_traits),
            "--responses",
            "state",
            "--predictors",
            "body_size",
            "--tree-id",
            "OGPARA",
            "--outfile",
            str(output),
            "--random-effects-out",
            str(random_output),
        ]
    )

    result = pd.read_csv(output, sep="\t")
    assert set(result["response_family"]) == {"binomial"}
    assert set(result["event_random_effect"]) == {"species-tip"}
    random_effects = pd.read_csv(random_output, sep="\t")
    species_effects = random_effects[random_effects["effect_type"] == "species_logit"]
    assert set(species_effects["group_id"]) == {
        "Genus_a",
        "Genus_b",
        "Genus_c",
        "Genus_d",
        "Genus_e",
    }
    assert (
        species_effects.loc[
            species_effects["group_id"] == "Genus_a", "n_observations"
        ].item()
        == 2
    )


@pytest.mark.integration
def test_pgls_raw_mode_propagates_response_and_predictor_replicates_together(
    tmp_path,
):
    gene_tree, species_tree, expression, species_traits = _write_raw_regression_inputs(
        tmp_path,
        biological_replicates=True,
    )
    original = pd.read_csv(species_traits, sep="\t")
    rows = []
    for record in original.to_dict("records"):
        for replicate, offset in [("r1", -0.2), ("r2", 0.2)]:
            rows.append(
                {
                    "leaf_name": record["leaf_name"],
                    "predictor_sample": "{}_{}".format(record["leaf_name"], replicate),
                    "body_size": float(record["body_size"]) + offset,
                }
            )
    pd.DataFrame(rows).to_csv(species_traits, sep="\t", index=False)
    prefix = tmp_path / "predictor-replicates"

    main(
        [
            "regress",
            "--gene-tree",
            str(gene_tree),
            "--species-tree",
            str(species_tree),
            "--expression",
            str(expression),
            "--species-traits",
            str(species_traits),
            "--responses",
            "expression",
            "--predictors",
            "body_size",
            "--tree-id",
            "OG1",
            "--biological-id",
            "sample_id",
            "--predictor-biological-id",
            "predictor_sample",
            "--gene-evolution-model",
            "lambda",
            "--species-evolution-model",
            "lambda",
            "--out-prefix",
            str(prefix),
        ]
    )

    covariance = pd.read_csv(
        tmp_path / "predictor-replicates.predictor-sampling-covariance.tsv",
        sep="\t",
    )
    response_covariance = pd.read_csv(
        tmp_path / "predictor-replicates.response-sampling-covariance.tsv",
        sep="\t",
    )
    summary = pd.read_csv(
        tmp_path / "predictor-replicates.predictor-tip-summary.tsv",
        sep="\t",
    )
    response_summary = pd.read_csv(
        tmp_path / "predictor-replicates.response-tip-summary.tsv",
        sep="\t",
    )
    result = pd.read_csv(tmp_path / "predictor-replicates.regression.tsv", sep="\t")
    assert len(covariance) == 10
    assert len(response_covariance) == 10
    assert set(summary["n_biological"]) == {2}
    assert set(response_summary["n_biological"]) == {2}
    assert set(result["measurement_error_model"]) == {"latent-predictor"}
    assert set(result["reml"]) == {"no"}
    assert set(result["covariance_estimator"]) == {"gaussian-eiv-ML"}
    assert result.iloc[0]["mean_sampling_variance"] > 0.0
    assert result.iloc[0]["mean_predictor_sampling_variance"] > 0.0
    assert set(result["response_evolution_parameter_status"]) == {"estimated"}
    assert set(result["predictor_evolution_parameter_status"]) == {"estimated"}


@pytest.mark.integration
def test_reconciled_poisson_conditions_replicated_predictor_on_species_model(
    tmp_path,
):
    gene_tree, species_tree, expression, species_traits = _write_raw_regression_inputs(
        tmp_path
    )
    expression_frame = pd.read_csv(expression, sep="\t").rename(
        columns={"expression": "count"}
    )
    expression_frame["count"] = expression_frame["count"].astype(int)
    expression_frame.to_csv(expression, sep="\t", index=False)
    original = pd.read_csv(species_traits, sep="\t")
    rows = []
    for record in original.to_dict("records"):
        for replicate, offset in (("r1", -0.5), ("r2", 0.5)):
            rows.append(
                {
                    "leaf_name": record["leaf_name"],
                    "sample": "{}_{}".format(record["leaf_name"], replicate),
                    "body_size": float(record["body_size"]) + offset,
                }
            )
    pd.DataFrame(rows).to_csv(species_traits, sep="\t", index=False)

    coefficients = {}
    for model in ("brownian", "independent"):
        prefix = tmp_path / "tip-latent-{}".format(model)
        main(
            [
                "regress",
                "--gene-tree",
                str(gene_tree),
                "--species-tree",
                str(species_tree),
                "--expression",
                str(expression),
                "--species-traits",
                str(species_traits),
                "--responses",
                "count",
                "--response-family",
                "count=poisson",
                "--predictors",
                "body_size",
                "--predictor-biological-id",
                "sample",
                "--species-evolution-model",
                model,
                "--tree-id",
                "OGLATENT",
                "--out-prefix",
                str(prefix),
            ]
        )
        result = pd.read_csv(str(prefix) + ".regression.tsv", sep="\t")
        coefficient = result[
            (result["term"] == "body_size") & (result["term_test"] == "coefficient")
        ].iloc[0]
        coefficients[model] = float(coefficient["coefficient"])
        assert coefficient["measurement_error_model"] == "latent-predictor"
        assert coefficient["mean_predictor_sampling_variance"] > 0.0
        assert coefficient["mean_latent_predictor_variance"] > 0.0
        assert coefficient["predictor_evolutionary_rate"] > 0.0

    assert coefficients["brownian"] != pytest.approx(coefficients["independent"])


def test_precomputed_reconciled_pgls_accepts_predictor_sampling_covariance():
    response = pd.DataFrame(
        [_response_row(1, 2.0), _response_row(2, 4.0), _response_row(3, 6.0)]
    )
    predictor = _predictor_table()
    predictor["contrast_variance"] = [1.0, 1.5, 2.0]
    covariance = _sampling_covariance_table(
        predictor.rename(columns={"branch_clade_id": "gene_clade_id"}),
        np.diag([0.05, 0.05, 0.05]),
    )

    result = fit_reconciled_pgls(
        response,
        predictor,
        ["expression"],
        ["body_size"],
        predictor_sampling_covariance=covariance,
        event_random_effect="no",
        lineage_random_slope="no",
    )

    assert result.iloc[0]["measurement_error_model"] == "latent-predictor"
    assert result.iloc[0]["mean_predictor_sampling_variance"] > 0.0
    assert result.iloc[0]["predictor_evolutionary_rate"] > 0.0
    assert result.iloc[0]["standard_error"] > 0.0


def test_predictor_sampling_sidecar_must_cover_every_selected_trait():
    response = pd.DataFrame(
        [_response_row(1, 2.0), _response_row(2, 4.0), _response_row(3, 6.0)]
    )
    first = _predictor_table()
    first["contrast_variance"] = 1.0
    second = _predictor_table(values=(0.5, 1.5, 2.5))
    second["trait"] = "temperature"
    second["contrast_variance"] = 1.0
    covariance = _sampling_covariance_table(
        first.rename(columns={"branch_clade_id": "gene_clade_id"}),
        np.eye(3) * 0.05,
    )

    with pytest.raises(ValueError, match="missing selected trait.*temperature"):
        fit_reconciled_pgls(
            response,
            pd.concat([first, second], ignore_index=True),
            ["expression"],
            ["body_size", "temperature"],
            predictor_sampling_covariance=covariance,
            event_random_effect="no",
            lineage_random_slope="no",
        )


def test_grouped_predictor_uncertainty_validation_is_explicit_and_lossless():
    factor_a = sparse.eye(2, format="csr")
    factor_b = sparse.csr_matrix(np.asarray([[0.2, 0.0], [0.0, 0.3]]))
    state = {
        "term_names": ("state[a]", "state[b]"),
        "event_ids": ("event1", "event2"),
        "uncertainty": JointPredictorUncertainty(factors=(factor_a, factor_b)),
    }
    groups = {"state": ("state[a]", "state[b]")}
    normalized = regression_mod._normalize_grouped_uncertainties(
        [state], ["state[a]", "state[b]"], "hierarchical", groups
    )

    with pytest.raises(ValueError, match="missing species event.*event3"):
        regression_mod._grouped_predictor_uncertainties_for_rows(
            pd.DataFrame({"species_event_id": ["event1", "event3"]}), normalized
        )
    with pytest.raises(ValueError, match="assigns term.*more than once"):
        regression_mod._normalize_grouped_uncertainties(
            [state, state], ["state[a]", "state[b]"], "hierarchical", groups
        )
    with pytest.raises(ValueError, match="must match one predictor group"):
        regression_mod._normalize_grouped_uncertainties(
            [
                {
                    **state,
                    "term_names": ("state[a]",),
                    "uncertainty": JointPredictorUncertainty(factors=(factor_a,)),
                }
            ],
            ["state[a]", "state[b]"],
            "hierarchical",
            groups,
        )

    row_scale = np.asarray([2.0, 3.0])
    uncertainty = JointPredictorUncertainty(
        factors=(factor_a, factor_b), row_scale=row_scale
    )
    selected, columns = regression_mod._subset_predictor_uncertainties(
        [uncertainty], [(0, 1)], [0]
    )
    assert columns == [0]
    np.testing.assert_array_equal(selected[0].row_scale, row_scale)


def test_predictor_factor_loading_sidecar_matches_explicit_covariance(monkeypatch):
    response = pd.DataFrame(
        [_response_row(1, 2.0), _response_row(2, 4.0), _response_row(3, 6.0)]
    )
    predictor = _predictor_table()
    predictor["contrast_variance"] = [1.0, 1.5, 2.0]
    diagonal = np.asarray([0.05, 0.10, 0.15])
    explicit = _sampling_covariance_table(
        predictor.rename(columns={"branch_clade_id": "gene_clade_id"}),
        np.diag(diagonal),
    )
    factor_rows = []
    for index, row in predictor.iterrows():
        factor_rows.append(
            {
                "tree_id": row["tree_id"],
                "trait": row["trait"],
                "contrast_id_1": row["branch_clade_id"],
                "contrast_id_2": "latent:{}".format(index),
                "sampling_covariance": np.sqrt(diagonal[index]),
                "covariance_representation": "factor-loading",
            }
        )
    common = dict(
        response_contrasts=response,
        predictor_contrasts=predictor,
        responses=["expression"],
        predictors=["body_size"],
        event_random_effect="no",
        lineage_random_slope="no",
    )
    explicit_result = fit_reconciled_pgls(
        **common, predictor_sampling_covariance=explicit
    )

    def reject_materialization(_covariance):
        raise AssertionError("factor-loading covariance was materialized")

    monkeypatch.setattr(
        regression_mod, "materialize_covariance", reject_materialization
    )
    factor_result = fit_reconciled_pgls(
        **common, predictor_sampling_covariance=pd.DataFrame(factor_rows)
    )

    np.testing.assert_allclose(
        factor_result[["coefficient", "standard_error", "log_likelihood"]],
        explicit_result[["coefficient", "standard_error", "log_likelihood"]],
        rtol=5e-5,
        atol=1e-7,
    )
    np.testing.assert_allclose(
        factor_result[
            ["mean_latent_predictor_variance", "predictor_uncertainty_fraction"]
        ],
        explicit_result[
            ["mean_latent_predictor_variance", "predictor_uncertainty_fraction"]
        ],
        rtol=5e-5,
        atol=1e-7,
    )


@pytest.mark.slow
def test_large_predictor_factor_loading_remains_structured(monkeypatch):
    size = regression_mod.MAX_DENSE_GAUSSIAN_OBSERVATIONS + 1
    predictor_values = np.linspace(-2.0, 2.0, size)
    predictor = _predictor_table(predictor_values)
    predictor["contrast_variance"] = 1.0
    response = pd.DataFrame(
        [
            _response_row(index, 1.5 * value + 0.1 * np.sin(index))
            for index, value in enumerate(predictor_values, start=1)
        ]
    )
    factor = pd.DataFrame(
        {
            "tree_id": "species",
            "trait": "body_size",
            "contrast_id_1": predictor["branch_clade_id"],
            "contrast_id_2": ["latent:{}".format(index) for index in range(size)],
            "sampling_covariance": np.sqrt(0.05),
            "covariance_representation": "factor-loading",
        }
    )

    def reject_materialization(_covariance):
        raise AssertionError("large factor-loading covariance was materialized")

    monkeypatch.setattr(
        regression_mod, "materialize_covariance", reject_materialization
    )
    result = fit_reconciled_pgls(
        response,
        predictor,
        ["expression"],
        ["body_size"],
        predictor_sampling_covariance=factor,
        event_random_effect="no",
        lineage_random_slope="no",
    )

    assert result.iloc[0]["coefficient"] == pytest.approx(1.5, rel=0.02)
    assert result.iloc[0]["measurement_error_model"] == "latent-predictor"


def test_repeated_paralogs_share_one_latent_species_event_uncertainty():
    predictor = _predictor_table()
    predictor["contrast_variance"] = [1.0, 1.5, 2.0]
    covariance = _sampling_covariance_table(
        predictor.rename(columns={"branch_clade_id": "gene_clade_id"}),
        np.diag([0.05, 0.10, 0.15]),
    )
    prepared_predictor = regression_mod._prepare_predictors(predictor, ["body_size"])
    prepared_covariance = regression_mod._prepare_sampling_covariance(
        covariance,
        option_name="--predictor-sampling-covariance",
    )
    posteriors = regression_mod._prepare_predictor_posteriors(
        prepared_predictor,
        ["body_size"],
        prepared_covariance,
    )
    repeated = pd.DataFrame(
        [
            _response_row(1, 2.0, gene_index=1),
            _response_row(1, 2.2, gene_index=2),
            _response_row(2, 4.0),
            _response_row(3, 6.0),
        ]
    )

    uncertainty = regression_mod._predictor_uncertainties_for_rows(
        repeated,
        ["body_size"],
        posteriors,
    )["body_size"]
    loading = continuous_predictor_loading(uncertainty, np.asarray([1.0]))
    covariance = (loading @ loading.T).toarray()

    np.testing.assert_allclose(covariance[0], covariance[1])
    np.testing.assert_allclose(covariance[:, 0], covariance[:, 1])


def test_compatible_grouped_covariance_updates_are_combined_exactly():
    diagonal = np.asarray([0.5, 0.7, 0.9, 1.1])
    first = np.asarray([[1.0, 0.0], [2.0, 0.0], [0.0, 1.5], [0.0, 3.0]])
    second = np.asarray([[0.5, 0.0], [1.0, 0.0], [0.0, -3.0], [0.0, -6.0]])
    factor = factor_diagonal_low_rank_updates(diagonal, [first, second])
    expected = np.diag(diagonal) + first @ first.T + second @ second.T
    values = np.asarray([1.0, -2.0, 3.0, 0.5])

    np.testing.assert_allclose(
        solve_factor(factor, values), np.linalg.solve(expected, values)
    )
    assert factor_logdet(factor) == pytest.approx(np.linalg.slogdet(expected)[1])


def test_zero_predictor_sampling_covariance_recovers_exact_reconciled_pgls():
    response = pd.DataFrame(
        [
            _response_row(1, 2.0),
            _response_row(2, 4.1),
            _response_row(3, 5.8),
            _response_row(4, 8.2),
        ]
    )
    predictor = _predictor_table(values=(1.0, 2.0, 3.0, 4.0))
    predictor["contrast_variance"] = [1.0, 1.0, 1.0, 1.0]
    common = dict(
        response_contrasts=response,
        predictor_contrasts=predictor,
        responses=["expression"],
        predictors=["body_size"],
        model="replicate-reml",
        reml=False,
    )
    exact = fit_reconciled_pgls(**common)
    zero_sampling = _sampling_covariance_table(
        predictor.rename(columns={"branch_clade_id": "gene_clade_id"}),
        np.zeros((4, 4)),
    )
    latent = fit_reconciled_pgls(
        **common,
        predictor_sampling_covariance=zero_sampling,
    )

    assert latent.iloc[0]["coefficient"] == pytest.approx(
        exact.iloc[0]["coefficient"], rel=1e-6
    )
    assert latent.iloc[0]["standard_error"] == pytest.approx(
        exact.iloc[0]["standard_error"], rel=5e-4
    )


def test_legacy_reconciled_pgls_rejects_predictor_sampling_covariance():
    with pytest.raises(ValueError, match="likelihood-based"):
        fit_reconciled_pgls(
            pd.DataFrame([_response_row(1, 2.0)]),
            _predictor_table(values=(1.0,)),
            ["expression"],
            ["body_size"],
            model="legacy",
            predictor_sampling_covariance=pd.DataFrame({"unused": [1]}),
        )


def test_pgls_raw_and_precomputed_modes_reject_mixed_or_incomplete_inputs():
    common = [
        "regress",
        "--responses",
        "expression",
        "--predictors",
        "body_size",
    ]
    with pytest.raises(ValueError, match="require '--gene-tree'"):
        main(common + ["--expression", "expression.tsv"])
    with pytest.raises(ValueError, match="Precomputed-contrast PGLS requires"):
        main(common)
    with pytest.raises(ValueError, match="cannot be combined with precomputed"):
        main(
            common
            + [
                "--gene-tree",
                "gene.nwk",
                "--species-tree",
                "species.nwk",
                "--expression",
                "expression.tsv",
                "--species-traits",
                "traits.tsv",
                "--tree-id",
                "OG1",
                "--infile",
                "gene-contrasts.tsv",
            ]
        )
    with pytest.raises(ValueError, match="Typed response/predictor options"):
        main(
            common
            + [
                "--infile",
                "gene-contrasts.tsv",
                "--predictor-contrasts",
                "species-contrasts.tsv",
                "--response-family",
                "expression=poisson",
            ]
        )
    raw_required = [
        "--species-tree",
        "species.nwk",
        "--expression",
        "expression.tsv",
        "--species-traits",
        "traits.tsv",
        "--tree-id",
        "OG1",
    ]
    with pytest.raises(ValueError, match="exactly one"):
        main(
            common
            + raw_required
            + [
                "--gene-tree",
                "gene.nwk",
                "--gene-tree-ensemble",
                "trees.nwk",
            ]
        )
    with pytest.raises(ValueError, match="one fixed --reconciliation-tree"):
        main(
            common
            + raw_required
            + [
                "--gene-tree-ensemble",
                "trees.nwk",
                "--reconciliation-tree",
                "reconciled.nwk",
            ]
        )


def test_pgls_raw_mode_rejects_mismatched_reconciliation_topology(tmp_path):
    gene_tree, species_tree, expression, species_traits = _write_raw_regression_inputs(
        tmp_path
    )
    reconciliation_tree = tmp_path / "reconciliation.nwk"
    reconciliation_tree.write_text(
        "(((Genus_a_g1:1,Genus_c_g1:1):1,Genus_b_g1:2):1,"
        "(Genus_d_g1:1,Genus_e_g1:1):2);"
    )
    with pytest.raises(ValueError, match="identical rooted topologies"):
        main(
            [
                "regress",
                "--gene-tree",
                str(gene_tree),
                "--reconciliation-tree",
                str(reconciliation_tree),
                "--species-tree",
                str(species_tree),
                "--expression",
                str(expression),
                "--species-traits",
                str(species_traits),
                "--responses",
                "expression",
                "--predictors",
                "body_size",
                "--tree-id",
                "OG1",
            ]
        )


@pytest.mark.integration
def test_pgls_raw_mode_without_prefix_writes_only_requested_primary_output(tmp_path):
    gene_tree, species_tree, expression, species_traits = _write_raw_regression_inputs(
        tmp_path
    )
    outfile = tmp_path / "result.tsv"
    assert (
        main(
            [
                "regress",
                "--gene-tree",
                str(gene_tree),
                "--species-tree",
                str(species_tree),
                "--expression",
                str(expression),
                "--species-traits",
                str(species_traits),
                "--responses",
                "expression",
                "--predictors",
                "body_size",
                "--tree-id",
                "OG1",
                "--outfile",
                str(outfile),
            ]
        )
        is None
    )
    result = pd.read_csv(outfile, sep="\t")
    assert result.iloc[0]["term"] == "body_size"
    assert not list(tmp_path.glob("result.*.tsv"))


@pytest.mark.integration
def test_pgls_raw_mode_applies_gene_and_species_evolution_models(tmp_path):
    gene_tree, species_tree, expression, species_traits = _write_raw_regression_inputs(
        tmp_path
    )
    prefix = tmp_path / "transformed"

    main(
        [
            "regress",
            "--gene-tree",
            str(gene_tree),
            "--species-tree",
            str(species_tree),
            "--expression",
            str(expression),
            "--species-traits",
            str(species_traits),
            "--responses",
            "expression",
            "--predictors",
            "body_size",
            "--tree-id",
            "OG1",
            "--gene-evolution-model",
            "kappa",
            "--gene-evolution-parameter",
            "0.7",
            "--species-evolution-model",
            "delta",
            "--species-evolution-parameter",
            "1.4",
            "--out-prefix",
            str(prefix),
        ]
    )

    gene_contrasts = pd.read_csv(tmp_path / "transformed.gene-contrasts.tsv", sep="\t")
    species_contrasts = pd.read_csv(
        tmp_path / "transformed.species-contrasts.tsv", sep="\t"
    )
    result = pd.read_csv(tmp_path / "transformed.regression.tsv", sep="\t")
    assert set(gene_contrasts["evolution_model"]) == {"kappa"}
    assert set(gene_contrasts["evolution_parameter"]) == {0.7}
    assert set(species_contrasts["evolution_model"]) == {"delta"}
    assert set(species_contrasts["evolution_parameter"]) == {1.4}
    assert set(result["response_evolution_model"]) == {"kappa"}
    assert set(result["response_evolution_parameter_name"]) == {"kappa"}
    assert set(result["response_evolution_parameter_status"]) == {"fixed"}
    assert set(result["response_branch_length_mode"]) == {"original"}
    assert set(result["predictor_evolution_model"]) == {"delta"}
    assert set(result["predictor_evolution_parameter_name"]) == {"delta"}
    assert set(result["predictor_evolution_parameter_status"]) == {"fixed"}
    assert set(result["predictor_branch_length_mode"]) == {"original"}


@pytest.mark.integration
def test_pgls_raw_parameterized_models_are_automatically_estimated(tmp_path):
    gene_tree, species_tree, expression, species_traits = _write_raw_regression_inputs(
        tmp_path
    )
    prefix = tmp_path / "estimated"

    main(
        [
            "regress",
            "--gene-tree",
            str(gene_tree),
            "--species-tree",
            str(species_tree),
            "--expression",
            str(expression),
            "--species-traits",
            str(species_traits),
            "--responses",
            "expression",
            "--predictors",
            "body_size",
            "--tree-id",
            "OG1",
            "--gene-evolution-model",
            "kappa",
            "--gene-evolution-parameter",
            "auto",
            "--species-evolution-model",
            "lambda",
            "--out-prefix",
            str(prefix),
        ]
    )

    gene_contrasts = pd.read_csv(tmp_path / "estimated.gene-contrasts.tsv", sep="\t")
    species_contrasts = pd.read_csv(
        tmp_path / "estimated.species-contrasts.tsv", sep="\t"
    )
    result = pd.read_csv(tmp_path / "estimated.regression.tsv", sep="\t")
    gene_parameter = float(gene_contrasts["evolution_parameter"].iloc[0])
    species_parameter = float(species_contrasts["evolution_parameter"].iloc[0])
    assert 0.0 <= gene_parameter <= 3.0
    assert 0.0 <= species_parameter <= 1.0
    assert gene_contrasts["evolution_parameter"].eq(gene_parameter).all()
    assert species_contrasts["evolution_parameter"].eq(species_parameter).all()
    assert set(result["response_evolution_parameter_status"]) == {"estimated"}
    assert set(result["predictor_evolution_parameter_status"]) == {"estimated"}
    assert np.isfinite(result["predictor_evolution_log_likelihood"]).all()

    identity_out = tmp_path / "identity.tsv"
    main(
        [
            "regress",
            "--gene-tree",
            str(gene_tree),
            "--species-tree",
            str(species_tree),
            "--expression",
            str(expression),
            "--species-traits",
            str(species_traits),
            "--responses",
            "expression",
            "--predictors",
            "body_size",
            "--tree-id",
            "OG1",
            "--gene-evolution-model",
            "kappa",
            "--gene-evolution-parameter",
            "1",
            "--species-evolution-model",
            "lambda",
            "--species-evolution-parameter",
            str(species_parameter),
            "--outfile",
            str(identity_out),
        ]
    )
    identity = pd.read_csv(identity_out, sep="\t")
    assert result.iloc[0]["log_likelihood"] >= identity.iloc[0]["log_likelihood"] - 1e-8


def test_pgls_raw_auto_parameter_is_rejected_for_parameterless_model():
    with pytest.raises(ValueError, match="no shape parameter"):
        main(
            [
                "regress",
                "--gene-tree",
                "gene.nwk",
                "--species-tree",
                "species.nwk",
                "--expression",
                "expression.tsv",
                "--species-traits",
                "traits.tsv",
                "--responses",
                "expression",
                "--predictors",
                "body_size",
                "--tree-id",
                "OG1",
                "--gene-evolution-parameter",
                "auto",
            ]
        )


def test_pgls_raw_auto_gene_parameter_rejects_legacy_model_before_io():
    with pytest.raises(ValueError, match="likelihood-based"):
        main(
            [
                "regress",
                "--gene-tree",
                "gene.nwk",
                "--species-tree",
                "species.nwk",
                "--expression",
                "expression.tsv",
                "--species-traits",
                "traits.tsv",
                "--responses",
                "expression",
                "--predictors",
                "body_size",
                "--tree-id",
                "OG1",
                "--gene-evolution-model",
                "lambda",
                "--model",
                "legacy",
            ]
        )


@pytest.mark.integration
def test_pgls_raw_bootstrap_refits_automatic_gene_parameter(tmp_path):
    gene_tree, species_tree, expression, species_traits = _write_raw_regression_inputs(
        tmp_path, biological_replicates=True
    )
    outfile = tmp_path / "bootstrap.tsv"

    main(
        [
            "regress",
            "--gene-tree",
            str(gene_tree),
            "--species-tree",
            str(species_tree),
            "--expression",
            str(expression),
            "--species-traits",
            str(species_traits),
            "--responses",
            "expression",
            "--predictors",
            "body_size",
            "--tree-id",
            "OG1",
            "--biological-id",
            "sample_id",
            "--gene-evolution-model",
            "lambda",
            "--inference",
            "parametric-bootstrap",
            "--bootstrap-replicates",
            "2",
            "--seed",
            "7",
            "--outfile",
            str(outfile),
        ]
    )

    result = pd.read_csv(outfile, sep="\t")
    assert result.iloc[0]["response_evolution_parameter_status"] == "estimated"
    assert result.iloc[0]["response_evolution_parameter_bootstrap_refit"] == "yes"
    assert result.iloc[0]["inference_method"] == "parametric-bootstrap"
    assert np.isfinite(result.iloc[0]["standard_error"])


@pytest.mark.parametrize("suffix", [".regression.tsv", ".regression-bundle.lock"])
def test_regression_bundle_audit_path_cannot_collide_with_generated_output(
    tmp_path, suffix
):
    prefix = tmp_path / "analysis"
    with pytest.raises(ValueError, match="Output paths must be distinct"):
        main(
            [
                "regress",
                "--gene-tree",
                "gene.nwk",
                "--species-tree",
                "species.nwk",
                "--expression",
                "expression.tsv",
                "--species-traits",
                "traits.tsv",
                "--responses",
                "expression",
                "--predictors",
                "body_size",
                "--tree-id",
                "OG1",
                "--out-prefix",
                str(prefix),
                "--audit",
                str(tmp_path / "analysis{}".format(suffix)),
            ]
        )


def _minimal_pipeline_artifacts(*, replicate_aware=False):
    frame = pd.DataFrame({"value": [1.0]})
    return RegressionPipelineArtifacts(
        reconciliation=frame.copy(),
        gene_contrasts=frame.copy(),
        species_contrasts=frame.copy(),
        response_sampling_covariance=frame.copy() if replicate_aware else None,
        response_tip_summary=frame.copy() if replicate_aware else None,
        results=frame.copy(),
        random_effects=frame.copy(),
    )


def test_regression_bundle_rejects_nonregular_target_before_writing_any_file(tmp_path):
    prefix = tmp_path / "analysis"
    blocked = tmp_path / "analysis.response-sampling-covariance.tsv"
    blocked.mkdir()

    with pytest.raises(ValueError, match="must be a regular file"):
        write_regression_bundle(
            str(prefix),
            _minimal_pipeline_artifacts(replicate_aware=True),
        )

    generated = [
        path
        for path in tmp_path.glob("analysis.*")
        if path != blocked and path.is_file()
    ]
    assert generated == []


def test_regression_bundle_commit_failure_restores_every_existing_output(
    monkeypatch, tmp_path
):
    prefix = tmp_path / "analysis"
    paths = {
        name: path
        for name, path in regression_pipeline_mod.regression_bundle_paths(
            str(prefix)
        ).items()
        if name not in {"response_sampling_covariance_out", "response_tip_summary_out"}
    }
    for path in paths.values():
        with open(path, "w") as handle:
            handle.write("original\n")
    real_replace = regression_pipeline_mod._replace_output
    replace_calls = 0

    def fail_second_replace(source, target):
        nonlocal replace_calls
        replace_calls += 1
        if replace_calls == 2:
            raise OSError("simulated bundle commit failure")
        real_replace(source, target)

    monkeypatch.setattr(regression_pipeline_mod, "_replace_output", fail_second_replace)
    with pytest.raises(OSError, match="bundle commit failure"):
        write_regression_bundle(str(prefix), _minimal_pipeline_artifacts())

    assert all(Path(path).read_text() == "original\n" for path in paths.values())
    assert not [
        path for path in tmp_path.iterdir() if path.name.startswith(".analysis.")
    ]


def test_explicit_output_transactions_are_isolated_across_concurrent_writers(
    monkeypatch, tmp_path
):
    result_path = tmp_path / "result.tsv"
    sidecar_path = tmp_path / "sidecar.tsv"
    first_entered = threading.Event()
    release_first = threading.Event()
    second_entered = threading.Event()
    real_commit = regression_pipeline_mod._commit_regression_outputs
    commit_count = 0
    commit_guard = threading.Lock()

    def controlled_commit(staged_outputs):
        nonlocal commit_count
        with commit_guard:
            commit_count += 1
            position = commit_count
        if position == 1:
            first_entered.set()
            assert release_first.wait(5)
        else:
            second_entered.set()
        real_commit(staged_outputs)

    monkeypatch.setattr(
        regression_pipeline_mod, "_commit_regression_outputs", controlled_commit
    )
    errors = []

    def write(label):
        try:
            regression_pipeline_mod._write_dataframes_transactionally(
                [
                    (str(result_path), pd.DataFrame({"run": [label]})),
                    (str(sidecar_path), pd.DataFrame({"run": [label]})),
                ]
            )
        except BaseException as exc:  # pragma: no cover - asserted below
            errors.append(exc)

    first = threading.Thread(target=write, args=("A",))
    second = threading.Thread(target=write, args=("B",))
    first.start()
    assert first_entered.wait(5)
    second.start()
    assert not second_entered.wait(0.2)
    release_first.set()
    first.join(5)
    second.join(5)

    assert not first.is_alive()
    assert not second.is_alive()
    assert errors == []
    assert pd.read_csv(result_path, sep="\t").iloc[0]["run"] == "B"
    assert pd.read_csv(sidecar_path, sep="\t").iloc[0]["run"] == "B"
    assert not list(tmp_path.glob(".nwkit-output-*.lock"))


@pytest.mark.integration
def test_pgls_failed_bundle_commit_is_rolled_back_and_audits_planned_outputs(
    monkeypatch, tmp_path
):
    gene_tree, species_tree, expression, species_traits = _write_raw_regression_inputs(
        tmp_path,
        biological_replicates=True,
    )
    prefix = tmp_path / "analysis"
    audit = tmp_path / "analysis.audit.jsonl"
    real_replace = regression_pipeline_mod._replace_output
    replace_calls = 0

    def fail_second_replace(source, target):
        nonlocal replace_calls
        replace_calls += 1
        if replace_calls == 2:
            raise OSError("simulated bundle commit failure")
        real_replace(source, target)

    monkeypatch.setattr(regression_pipeline_mod, "_replace_output", fail_second_replace)
    with pytest.raises(OSError, match="bundle commit failure"):
        main(
            [
                "regress",
                "--gene-tree",
                str(gene_tree),
                "--species-tree",
                str(species_tree),
                "--expression",
                str(expression),
                "--species-traits",
                str(species_traits),
                "--responses",
                "expression",
                "--predictors",
                "body_size",
                "--tree-id",
                "OG1",
                "--biological-id",
                "sample_id",
                "--out-prefix",
                str(prefix),
                "--audit",
                str(audit),
            ]
        )

    expected = {str(path.resolve()) for path in tmp_path.glob("analysis.*.tsv")}
    assert expected == set()
    record = json.loads(audit.read_text())
    assert record["status"] == "error"
    assert record["outputs"] == []
    assert {item["path"] for item in record["planned_outputs"]} == {
        str(path.resolve())
        for path in [
            tmp_path / "analysis.reconciliation.tsv",
            tmp_path / "analysis.gene-contrasts.tsv",
            tmp_path / "analysis.species-contrasts.tsv",
            tmp_path / "analysis.response-sampling-covariance.tsv",
            tmp_path / "analysis.response-tip-summary.tsv",
            tmp_path / "analysis.random-effects.tsv",
            tmp_path / "analysis.regression.tsv",
        ]
    }


@pytest.mark.integration
def test_pgls_raw_mode_accepts_separate_nhx_reconciliation_tree(tmp_path):
    gene_tree, species_tree, expression, species_traits = _write_raw_regression_inputs(
        tmp_path
    )
    reconciliation_tree = tmp_path / "gene-reconciled.nhx"
    reconciliation_tree.write_text(
        "(((Genus_a_g1:9,Genus_b_g1:9):9[&&NHX:D=N],Genus_c_g1:9)"
        ":9[&&NHX:D=N],(Genus_d_g1:9,Genus_e_g1:9):9[&&NHX:D=N])"
        ":0[&&NHX:D=N];"
    )
    prefix = tmp_path / "analysis"

    main(
        [
            "regress",
            "--gene-tree",
            str(gene_tree),
            "--reconciliation-tree",
            str(reconciliation_tree),
            "--species-tree",
            str(species_tree),
            "--expression",
            str(expression),
            "--species-traits",
            str(species_traits),
            "--responses",
            "expression",
            "--predictors",
            "body_size",
            "--tree-id",
            "OG1",
            "--event-source",
            "nhx",
            "--out-prefix",
            str(prefix),
        ]
    )

    reconciliation = pd.read_csv(tmp_path / "analysis.reconciliation.tsv", sep="\t")
    internal = reconciliation[reconciliation["node_class"].isin(["root", "intnode"])]
    assert set(internal["event_source"]) == {"nhx"}
    assert set(internal["event_type"]) == {"speciation"}
    contrasts = pd.read_csv(tmp_path / "analysis.gene-contrasts.tsv", sep="\t")
    assert len(contrasts) == 4


@pytest.mark.integration
def test_pgls_raw_known_se_supports_multiple_responses_and_predictors(tmp_path):
    gene_tree, species_tree, expression, species_traits = _write_raw_regression_inputs(
        tmp_path
    )
    expression_table = pd.read_csv(expression, sep="\t")
    expression_table["expression_se"] = [0.2, 0.3, 0.2, 0.4, 0.3]
    expression_table["expression_alt"] = [11.0, 8.0, 10.0, 4.0, 2.0]
    expression_table["expression_alt_se"] = [0.4, 0.2, 0.3, 0.2, 0.4]
    expression_table["expression_n"] = [4, 5, 6, 4, 5]
    expression_table["expression_alt_n"] = [5, 4, 5, 6, 4]
    expression_table.to_csv(expression, sep="\t", index=False)
    species_table = pd.read_csv(species_traits, sep="\t")
    species_table["temperature"] = [5.0, 3.0, 6.0, 2.0, 8.0]
    species_table["body_size_se"] = [0.2, 0.3, 0.2, 0.4, 0.3]
    species_table["temperature_se"] = [0.3, 0.2, 0.4, 0.2, 0.3]
    species_table.to_csv(species_traits, sep="\t", index=False)
    prefix = tmp_path / "analysis"

    main(
        [
            "regress",
            "--gene-tree",
            str(gene_tree),
            "--species-tree",
            str(species_tree),
            "--expression",
            str(expression),
            "--species-traits",
            str(species_traits),
            "--responses",
            "expression,expression_alt",
            "--predictors",
            "body_size,temperature",
            "--tree-id",
            "OG1",
            "--within-variance",
            "known-se",
            "--standard-error-columns",
            "expression_se,expression_alt_se",
            "--sample-size-columns",
            "expression_n,expression_alt_n",
            "--predictor-within-variance",
            "known-se",
            "--predictor-standard-error-columns",
            "body_size_se,temperature_se",
            "--gene-evolution-model",
            "lambda",
            "--species-evolution-model",
            "kappa",
            "--out-prefix",
            str(prefix),
        ]
    )

    result = pd.read_csv(tmp_path / "analysis.regression.tsv", sep="\t")
    assert set(result["response"]) == {"expression", "expression_alt"}
    assert set(result["term"]) == {"body_size", "temperature"}
    assert len(result) == 4
    assert set(result["response_evolution_parameter_status"]) == {"estimated"}
    assert set(result["predictor_evolution_parameter_status"]) == {"estimated"}
    assert set(result["measurement_error_model"]) == {"latent-predictor"}
    assert (
        result.groupby("response")["response_evolution_parameter"].nunique().eq(1).all()
    )
    assert result.groupby("term")["predictor_evolution_parameter"].nunique().eq(1).all()
    covariance = pd.read_csv(
        tmp_path / "analysis.response-sampling-covariance.tsv", sep="\t"
    )
    assert set(covariance["trait"]) == {"expression", "expression_alt"}
    assert len(covariance) == 20
    tip_summary = pd.read_csv(tmp_path / "analysis.response-tip-summary.tsv", sep="\t")
    assert set(tip_summary["variance_method"]) == {"known-se"}
    assert len(tip_summary) == 10
    predictor_covariance = pd.read_csv(
        tmp_path / "analysis.predictor-sampling-covariance.tsv", sep="\t"
    )
    assert set(predictor_covariance["trait"]) == {"body_size", "temperature"}
    assert len(predictor_covariance) == 20
    predictor_summary = pd.read_csv(
        tmp_path / "analysis.predictor-tip-summary.tsv", sep="\t"
    )
    assert set(predictor_summary["variance_method"]) == {"known-se"}
    assert len(predictor_summary) == 10


@pytest.mark.integration
def test_pgls_raw_gene_tree_stdin_has_primary_input_summary_and_hash(
    monkeypatch, tmp_path
):
    gene_tree, species_tree, expression, species_traits = _write_raw_regression_inputs(
        tmp_path
    )
    gene_tree_text = gene_tree.read_text()
    monkeypatch.setattr(sys, "stdin", io.StringIO(gene_tree_text))
    audit = tmp_path / "analysis.audit.jsonl"
    outfile = tmp_path / "analysis.regression.tsv"

    main(
        [
            "regress",
            "--gene-tree",
            "-",
            "--species-tree",
            str(species_tree),
            "--expression",
            str(expression),
            "--species-traits",
            str(species_traits),
            "--responses",
            "expression",
            "--predictors",
            "body_size",
            "--tree-id",
            "OG1",
            "--outfile",
            str(outfile),
            "--audit",
            str(audit),
        ]
    )

    record = json.loads(audit.read_text())
    assert record["stdin"] == {
        "argument": "gene_tree",
        "sha256": hashlib.sha256(gene_tree_text.encode()).hexdigest(),
        "bytes": len(gene_tree_text.encode()),
    }
    assert record["primary_input"]["kind"] == "newick"
    assert record["primary_input"]["first_tree_tip_count"] == 5


def test_expected_scipy_numdiff_warning_is_suppressed_without_hiding_others(
    monkeypatch,
):
    sentinel = object()

    def noisy_minimize(*args, **kwargs):
        warnings.warn_explicit(
            "invalid value encountered in subtract",
            RuntimeWarning,
            filename="scipy/optimize/_numdiff.py",
            lineno=615,
            module="scipy.optimize._numdiff",
        )
        warnings.warn("independent numerical warning", RuntimeWarning, stacklevel=2)
        return sentinel

    monkeypatch.setattr(regression_mod, "minimize", noisy_minimize)
    with pytest.warns(RuntimeWarning, match="independent numerical warning") as caught:
        result = regression_mod._minimize_variance_components(object())
    assert result is sentinel
    assert len(caught) == 1
