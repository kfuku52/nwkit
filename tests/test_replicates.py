import numpy as np
import pandas as pd
import pytest
from ete4 import Tree

from nwkit import contrast as contrast_mod
from nwkit.cli import main
from nwkit.contrast import build_contrast_table
from nwkit.gaussian import DiagonalLowRankCovariance, materialize_covariance
from nwkit.replicates import estimate_likelihood_replicates, estimate_replicate_traits


def _raw_replicates():
    return pd.DataFrame(
        [
            {"leaf_name": "A", "sample": "a1", "expression": 1.0},
            {"leaf_name": "A", "sample": "a2", "expression": 3.0},
            {"leaf_name": "B", "sample": "b1", "expression": 4.0},
            {"leaf_name": "B", "sample": "b2", "expression": 6.0},
            {"leaf_name": "C", "sample": "c1", "expression": 7.0},
            {"leaf_name": "C", "sample": "c2", "expression": 9.0},
        ]
    )


def test_pooled_biological_replicates_estimate_tip_sampling_variance():
    estimates = estimate_replicate_traits(
        _raw_replicates(),
        ["A", "B", "C"],
        ["expression"],
        biological_id="sample",
        tree_id="OG1",
    )

    assert estimates.values_by_trait["expression"] == {
        "A": pytest.approx(2.0),
        "B": pytest.approx(5.0),
        "C": pytest.approx(8.0),
    }
    np.testing.assert_allclose(
        estimates.sampling_covariance_by_trait["expression"], np.ones(3)
    )
    summary = estimates.tip_summary.set_index("leaf_name")
    assert set(summary["n_biological"]) == {2}
    assert set(summary["n_technical"]) == {2}
    np.testing.assert_allclose(summary["within_sd"], np.sqrt(2.0))
    assert set(summary["variance_method"]) == {"pooled"}


def test_technical_replicates_require_explicit_aggregation_and_do_not_add_n():
    dataframe = pd.DataFrame(
        [
            {"leaf_name": "A", "sample": "a1", "tech": "1", "expression": 1.0},
            {"leaf_name": "A", "sample": "a1", "tech": "2", "expression": 3.0},
            {"leaf_name": "A", "sample": "a2", "tech": "1", "expression": 4.0},
            {"leaf_name": "B", "sample": "b1", "tech": "1", "expression": 5.0},
            {"leaf_name": "B", "sample": "b1", "tech": "2", "expression": 7.0},
            {"leaf_name": "B", "sample": "b2", "tech": "1", "expression": 8.0},
        ]
    )
    with pytest.raises(ValueError, match="Technical replicates are present"):
        estimate_replicate_traits(
            dataframe,
            ["A", "B"],
            ["expression"],
            biological_id="sample",
            technical_id="tech",
        )

    estimates = estimate_replicate_traits(
        dataframe,
        ["A", "B"],
        ["expression"],
        biological_id="sample",
        technical_id="tech",
        technical_aggregation="mean",
    )
    summary = estimates.tip_summary.set_index("leaf_name")
    assert summary.loc["A", "n_biological"] == 2
    assert summary.loc["B", "n_biological"] == 2
    assert summary.loc["A", "n_technical"] == 3
    assert summary.loc["B", "n_technical"] == 3
    assert estimates.values_by_trait["expression"] == {
        "A": pytest.approx(3.0),
        "B": pytest.approx(7.0),
    }

    with_missing = dataframe.copy()
    with_missing.loc[
        (with_missing["leaf_name"] == "A") & (with_missing["tech"] == "2"),
        "expression",
    ] = np.nan
    missing_estimates = estimate_replicate_traits(
        with_missing,
        ["A", "B"],
        ["expression"],
        biological_id="sample",
        technical_id="tech",
        technical_aggregation="mean",
    )
    missing_summary = missing_estimates.tip_summary.set_index("leaf_name")
    assert missing_summary.loc["A", "n_biological"] == 2
    assert missing_summary.loc["A", "n_technical"] == 2
    assert missing_summary.loc["B", "n_technical"] == 3


def test_replicate_identifier_options_are_validated_even_without_duplicates():
    dataframe = _raw_replicates().drop_duplicates("leaf_name").copy()
    dataframe["tech"] = ["t1", "", "t3"]
    with pytest.raises(ValueError, match="must not contain missing or empty"):
        estimate_replicate_traits(
            dataframe,
            ["A", "B", "C"],
            ["expression"],
            biological_id="sample",
            technical_id="tech",
        )
    with pytest.raises(ValueError, match="requires '--technical-id'"):
        estimate_replicate_traits(
            dataframe.drop(columns="tech"),
            ["A", "B", "C"],
            ["expression"],
            biological_id="sample",
            technical_aggregation="mean",
        )
    with pytest.raises(ValueError, match="must be distinct"):
        estimate_replicate_traits(
            dataframe,
            ["A", "B", "C"],
            ["expression"],
            biological_id="sample",
            technical_id="sample",
        )


def test_numeric_leaf_identifiers_match_string_tree_tip_names():
    dataframe = pd.DataFrame(
        {
            "leaf_name": [1, 1, 2, 2],
            "sample": ["a", "b", "c", "d"],
            "expression": [1.0, 3.0, 4.0, 6.0],
        }
    )
    estimates = estimate_replicate_traits(
        dataframe,
        ["1", "2"],
        ["expression"],
        biological_id="sample",
    )
    assert estimates.values_by_trait["expression"] == {
        "1": pytest.approx(2.0),
        "2": pytest.approx(5.0),
    }


def test_leaf_specific_variance_rejects_a_leaf_with_one_replicate():
    dataframe = _raw_replicates().query("not (leaf_name == 'C' and sample == 'c2')")
    with pytest.raises(ValueError, match="at least two biological replicates"):
        estimate_replicate_traits(
            dataframe,
            ["A", "B", "C"],
            ["expression"],
            biological_id="sample",
            within_variance="leaf",
        )


def test_batch_adjustment_is_fitted_and_confounding_is_rejected():
    rows = []
    for leaf, baseline in [("A", 1.0), ("B", 4.0)]:
        for batch, shift in [("x", 0.0), ("y", 2.0)]:
            for replicate, residual in [(1, -0.2), (2, 0.2)]:
                rows.append(
                    {
                        "leaf_name": leaf,
                        "sample": "{}{}{}".format(leaf, batch, replicate),
                        "batch": batch,
                        "expression": baseline + shift + residual,
                    }
                )
    estimates = estimate_replicate_traits(
        pd.DataFrame(rows),
        ["A", "B"],
        ["expression"],
        biological_id="sample",
        batch="batch",
    )
    assert estimates.values_by_trait["expression"] == {
        "A": pytest.approx(2.0),
        "B": pytest.approx(5.0),
    }
    assert set(estimates.tip_summary["batch_adjusted"]) == {"yes"}
    assert estimates.model_by_trait["expression"] == "pooled-batch-adjusted"
    covariance = estimates.sampling_covariance_by_trait["expression"]
    assert isinstance(covariance, DiagonalLowRankCovariance)

    selected = pd.DataFrame(rows)
    leaf_design = np.column_stack(
        [(selected["leaf_name"] == leaf).to_numpy(float) for leaf in ["A", "B"]]
    )
    batch_design = (selected["batch"] == "y").to_numpy(float)[:, None]
    design = np.column_stack([leaf_design, batch_design])
    response = selected["expression"].to_numpy(float)
    gram = design.T @ design
    residual = response - design @ np.linalg.solve(gram, design.T @ response)
    variance = float(residual @ residual / (len(response) - design.shape[1]))
    transform = np.column_stack([np.eye(2), np.repeat([[0.5]], 2, axis=0)])
    expected = variance * transform @ np.linalg.inv(gram) @ transform.T
    np.testing.assert_allclose(materialize_covariance(covariance), expected)

    confounded = pd.DataFrame(
        [
            {"leaf_name": "A", "sample": "a1", "batch": "x", "expression": 1},
            {"leaf_name": "A", "sample": "a2", "batch": "x", "expression": 2},
            {"leaf_name": "B", "sample": "b1", "batch": "y", "expression": 3},
            {"leaf_name": "B", "sample": "b2", "batch": "y", "expression": 4},
        ]
    )
    with pytest.raises(ValueError, match="confounded"):
        estimate_replicate_traits(
            confounded,
            ["A", "B"],
            ["expression"],
            biological_id="sample",
            batch="batch",
        )


def test_known_standard_errors_propagate_through_the_pic_transform():
    tree = Tree("((A:1,B:1):1,C:2);", parser=1)
    dataframe = pd.DataFrame(
        {
            "leaf_name": ["A", "B", "C"],
            "expression": [1.0, 3.0, 5.0],
            "expression_se": [0.1, 0.2, 0.3],
        }
    )
    estimates = estimate_replicate_traits(
        dataframe,
        ["A", "B", "C"],
        ["expression"],
        within_variance="known-se",
        se_columns=["expression_se"],
        tree_id="OG1",
    )
    assert estimates.tip_summary["n_biological"].eq("").all()
    assert estimates.tip_summary["within_sd"].eq("").all()

    with_sample_size = dataframe.assign(expression_n=[4, 9, 16])
    sized = estimate_replicate_traits(
        with_sample_size,
        ["A", "B", "C"],
        ["expression"],
        within_variance="known-se",
        se_columns=["expression_se"],
        n_columns=["expression_n"],
    ).tip_summary.set_index("leaf_name")
    assert sized.loc["A", "n_biological"] == 4
    assert sized.loc["A", "within_sd"] == pytest.approx(0.2)
    assert sized["n_technical"].eq("").all()

    contrasts, covariance = build_contrast_table(
        tree,
        estimates.values_by_trait,
        tree_id="OG1",
        sampling_covariance_by_trait=estimates.sampling_covariance_by_trait,
        replicate_model_by_trait=estimates.model_by_trait,
        tip_summary=estimates.tip_summary,
        return_sampling_covariance=True,
    )
    by_taxa = contrasts.set_index("descendant_taxa")
    ab_id = by_taxa.loc["A,B", "branch_clade_id"]
    root_id = by_taxa.loc["A,B,C", "branch_clade_id"]
    assert by_taxa.loc["A,B", "sampling_variance"] == pytest.approx(0.05)
    assert by_taxa.loc["A,B,C", "sampling_variance"] == pytest.approx(0.1025)
    pair = covariance[
        (
            covariance[["contrast_id_1", "contrast_id_2"]]
            .isin([ab_id, root_id])
            .all(axis=1)
        )
        & (covariance["contrast_id_1"] != covariance["contrast_id_2"])
    ]
    assert pair.iloc[0]["sampling_covariance"] == pytest.approx(-0.015)


def test_large_contrast_sampling_covariance_uses_numerically_lossless_sparse_factor(
    monkeypatch,
):
    tree = Tree("((A:1,B:1):1,C:2);", parser=1)
    estimates = estimate_replicate_traits(
        pd.DataFrame(
            {
                "leaf_name": ["A", "B", "C"],
                "expression": [1.0, 2.0, 4.0],
                "expression_se": [0.2, 0.3, 0.4],
            }
        ),
        ["A", "B", "C"],
        ["expression"],
        within_variance="known-se",
        se_columns=["expression_se"],
        tree_id="OG1",
    )
    _, full = build_contrast_table(
        tree,
        estimates.values_by_trait,
        tree_id="OG1",
        sampling_covariance_by_trait=estimates.sampling_covariance_by_trait,
        replicate_model_by_trait=estimates.model_by_trait,
        tip_summary=estimates.tip_summary,
        return_sampling_covariance=True,
        compact_sampling_covariance=False,
    )
    monkeypatch.setattr(contrast_mod, "MAX_FULL_SAMPLING_COVARIANCE_CONTRASTS", 1)
    contrasts, compact = build_contrast_table(
        tree,
        estimates.values_by_trait,
        tree_id="OG1",
        sampling_covariance_by_trait=estimates.sampling_covariance_by_trait,
        replicate_model_by_trait=estimates.model_by_trait,
        tip_summary=estimates.tip_summary,
        return_sampling_covariance=True,
    )

    contrast_ids = contrasts["branch_clade_id"].astype(str).tolist()
    index = {contrast_id: position for position, contrast_id in enumerate(contrast_ids)}
    latent_ids = list(dict.fromkeys(compact["contrast_id_2"].astype(str)))
    latent_index = {
        latent_id: position for position, latent_id in enumerate(latent_ids)
    }
    loading = np.zeros((len(contrast_ids), len(latent_ids)))
    for row in compact.to_dict("records"):
        loading[
            index[str(row["contrast_id_1"])],
            latent_index[str(row["contrast_id_2"])],
        ] = row["sampling_covariance"]
    expected = np.zeros((len(contrast_ids), len(contrast_ids)))
    for row in full.to_dict("records"):
        first = index[str(row["contrast_id_1"])]
        second = index[str(row["contrast_id_2"])]
        expected[first, second] = expected[second, first] = row["sampling_covariance"]

    assert set(compact["covariance_representation"]) == {"factor-loading"}
    np.testing.assert_allclose(loading @ loading.T, expected)


def test_replicate_aware_contrast_cli_writes_all_three_outputs(tmp_path):
    tree_path = tmp_path / "gene.nwk"
    trait_path = tmp_path / "expression.tsv"
    contrast_path = tmp_path / "contrasts.tsv"
    covariance_path = tmp_path / "sampling-covariance.tsv"
    summary_path = tmp_path / "tip-summary.tsv"
    tree_path.write_text("((A:1,B:1):1,C:2);")
    _raw_replicates().to_csv(trait_path, sep="\t", index=False)

    assert (
        main(
            [
                "contrast",
                "--infile",
                str(tree_path),
                "--trait",
                str(trait_path),
                "--columns",
                "expression",
                "--biological-id",
                "sample",
                "--tree-id",
                "OG1",
                "--sampling-covariance-out",
                str(covariance_path),
                "--tip-summary-out",
                str(summary_path),
                "--outfile",
                str(contrast_path),
            ]
        )
        is None
    )
    contrasts = pd.read_csv(contrast_path, sep="\t")
    covariance = pd.read_csv(covariance_path, sep="\t")
    summary = pd.read_csv(summary_path, sep="\t")
    assert set(contrasts["replicate_model"]) == {"pooled"}
    assert sorted(contrasts["sampling_variance"]) == pytest.approx([1.5, 2.0])
    assert len(covariance) == 3
    assert len(summary) == 3


@pytest.mark.parametrize(
    "extra_options",
    [
        ["--technical-id", "tech"],
        ["--batch", "batch"],
        ["--within-variance", "leaf"],
        ["--tip-summary-out", "summary.tsv"],
    ],
)
def test_replicate_only_options_are_never_silently_ignored(tmp_path, extra_options):
    tree_path = tmp_path / "gene.nwk"
    trait_path = tmp_path / "expression.tsv"
    tree_path.write_text("(A:1,B:1);")
    trait_path.write_text("leaf_name\texpression\nA\t1\nB\t2\n")

    with pytest.raises(ValueError, match="Replicate option"):
        main(
            [
                "contrast",
                "--infile",
                str(tree_path),
                "--trait",
                str(trait_path),
                "--columns",
                "expression",
                *extra_options,
            ]
        )


def test_known_se_and_raw_replicate_options_are_mutually_exclusive(tmp_path):
    tree_path = tmp_path / "gene.nwk"
    trait_path = tmp_path / "expression.tsv"
    tree_path.write_text("(A:1,B:1);")
    trait_path.write_text(
        "leaf_name\texpression\texpression_se\nA\t1\t0.1\nB\t2\t0.2\n"
    )
    with pytest.raises(ValueError, match="Known-SE input cannot use"):
        main(
            [
                "contrast",
                "--infile",
                str(tree_path),
                "--trait",
                str(trait_path),
                "--columns",
                "expression",
                "--within-variance",
                "known-se",
                "--batch",
                "batch",
                "--sampling-covariance-out",
                str(tmp_path / "covariance.tsv"),
            ]
        )


def test_multivariate_replicates_allow_trait_specific_missing_tips():
    dataframe = _raw_replicates().rename(columns={"expression": "first"})
    dataframe["second"] = dataframe["first"] * 2.0
    dataframe.loc[dataframe["leaf_name"] == "B", "second"] = np.nan

    estimates = estimate_replicate_traits(
        dataframe,
        ["A", "B", "C"],
        ["first", "second"],
        biological_id="sample",
        allow_missing_traits={"first", "second"},
    )

    assert np.isnan(estimates.values_by_trait["second"]["B"])
    covariance = estimates.sampling_covariance_by_trait["second"]
    assert covariance[1] == 0.0
    summary = estimates.tip_summary.set_index(["trait", "leaf_name"])
    assert summary.loc[("second", "B"), "n_biological"] == 0


def test_known_se_multivariate_missingness_requires_paired_mean_and_se():
    dataframe = pd.DataFrame(
        {
            "leaf_name": ["A", "B", "C"],
            "first": [1.0, 2.0, 3.0],
            "first_se": [0.1, 0.1, 0.1],
            "second": [2.0, np.nan, 6.0],
            "second_se": [0.2, np.nan, 0.2],
        }
    )
    estimates = estimate_replicate_traits(
        dataframe,
        ["A", "B", "C"],
        ["first", "second"],
        within_variance="known-se",
        se_columns=["first_se", "second_se"],
        allow_missing_traits={"first", "second"},
    )

    assert np.isnan(estimates.values_by_trait["second"]["B"])
    assert estimates.sampling_covariance_by_trait["second"][1] == 0.0

    dataframe.loc[dataframe["leaf_name"] == "B", "second_se"] = 0.2
    with pytest.raises(ValueError, match="paired finite means"):
        estimate_replicate_traits(
            dataframe,
            ["A", "B", "C"],
            ["first", "second"],
            within_variance="known-se",
            se_columns=["first_se", "second_se"],
            allow_missing_traits={"first", "second"},
        )

    dataframe.loc[dataframe["leaf_name"] == "B", ["first", "first_se"]] = np.nan
    with pytest.raises(ValueError, match="paired finite means"):
        estimate_replicate_traits(
            dataframe,
            ["A", "B", "C"],
            ["first", "second"],
            within_variance="known-se",
            se_columns=["first_se", "second_se"],
            allow_missing_traits={"first", "second"},
        )


def test_censored_likelihood_replicates_preserve_missing_observations():
    dataframe = pd.DataFrame(
        {
            "leaf_name": ["A", "A", "B", "B"],
            "sample": ["a1", "a2", "b1", "b2"],
            "expression": [1.0, np.nan, 2.0, np.nan],
        }
    )
    estimates = estimate_likelihood_replicates(
        dataframe,
        ["A", "B"],
        ["expression"],
        biological_id="sample",
        allow_missing_traits={"expression"},
    )

    a_values = estimates.values_by_trait["expression"]["A"].values
    b_values = estimates.values_by_trait["expression"]["B"].values
    assert a_values[0] == 1.0 and np.isnan(a_values[1])
    assert b_values[0] == 2.0 and np.isnan(b_values[1])
