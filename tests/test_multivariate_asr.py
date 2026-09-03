import math

import numpy as np
import pandas as pd
import pytest

from nwkit.cli import main
from nwkit.multivariate_asr import compute_mvbm_marginals
from nwkit.util import read_tree


def tree_from(source):
    return read_tree(source, "1", True, quiet=True, rooted="yes")


def test_star_mvbm_matches_sample_covariance_and_root_distribution():
    tree = tree_from("(A:1,B:1,C:1,D:1,E:1,F:1)R;")
    matrix = np.asarray(
        [[0.0, 0.0], [1.0, 2.0], [2.0, 1.0], [3.0, 5.0], [4.0, 3.0], [6.0, 7.0]]
    )
    values = {name: row for name, row in zip("ABCDEF", matrix, strict=True)}
    posterior, fit = compute_mvbm_marginals(tree, values, ("x", "y"))
    expected = np.cov(matrix, rowvar=False, ddof=1)
    np.testing.assert_allclose(fit.sigma, expected, atol=1e-14)
    np.testing.assert_allclose(posterior[tree].mean, matrix.mean(axis=0))
    np.testing.assert_allclose(posterior[tree].covariance, expected / len(matrix))
    assert fit.sigma_rank == 2
    assert fit.fit_status == "ok"
    assert math.isfinite(fit.restricted_log_likelihood)


def test_missing_tip_has_joint_predictive_covariance():
    tree = tree_from("(A:1,B:1,C:2)R;")
    values = {"A": (0.0, 0.0), "B": (2.0, 1.0), "C": None}
    posterior, fit = compute_mvbm_marginals(tree, values, ("x", "y"))
    # Two observed vectors imply a rank-one covariance, which is still usable
    # for conditional moments even though it has no ordinary 2-D density.
    assert fit.fit_status == "singular_covariance"
    assert fit.restricted_log_likelihood is None
    np.testing.assert_allclose(
        posterior[tree["C"]].covariance, 2.5 * fit.sigma, atol=1e-14
    )


def test_mvbm_covariance_and_flat_root_likelihood_match_dense_gls():
    tree = tree_from("((A:0.4,B:1.2):0.7,(C:0.6,D:1.1):0.5,E:1.8)R;")
    values = {
        "A": (-1.0, 2.0),
        "B": (0.5, 4.0),
        "C": (3.0, -2.0),
        "D": (2.0, 5.0),
        "E": (-3.0, 1.0),
    }
    _, fit = compute_mvbm_marginals(tree, values, ("x", "y"))
    tips = [tree[name] for name in values]
    depths = {tree: 0.0}
    for node in tree.traverse(strategy="preorder"):
        if not node.is_root:
            depths[node] = depths[node.up] + float(node.dist)
    phylogenetic = np.asarray(
        [
            [depths[tree.common_ancestor([left, right])] for right in tips]
            for left in tips
        ]
    )
    response = np.asarray(list(values.values()))
    precision = np.linalg.inv(phylogenetic)
    ones = np.ones(len(tips))
    root_precision = float(ones @ precision @ ones)
    root_mean = (ones @ precision @ response) / root_precision
    residual = response - root_mean
    scatter = residual.T @ precision @ residual
    residual_df = len(tips) - 1
    expected_sigma = scatter / residual_df
    sign_tree, logdet_tree = np.linalg.slogdet(phylogenetic)
    sign_sigma, logdet_sigma = np.linalg.slogdet(expected_sigma)
    assert sign_tree == sign_sigma == 1
    quadratic = float(np.trace(np.linalg.solve(expected_sigma, scatter)))
    expected_log_likelihood = -0.5 * (
        residual_df * 2 * math.log(2 * math.pi)
        + 2 * logdet_tree
        + 2 * math.log(root_precision)
        + residual_df * logdet_sigma
        + quadratic
    )
    np.testing.assert_allclose(fit.sigma, expected_sigma, rtol=2e-13, atol=2e-14)
    assert fit.restricted_log_likelihood == pytest.approx(
        expected_log_likelihood, abs=2e-12
    )


def test_mvbm_supports_partial_tip_vectors():
    posterior, fit = compute_mvbm_marginals(
        tree_from("(A:1,B:1,C:1)R;"),
        {"A": (0.0, 1.0), "B": (2.0, None), "C": (1.0, 3.0)},
        ("x", "y"),
    )
    assert fit.num_effective_observations == 3
    imputed = posterior[next(node for node in posterior if node.name == "B")]
    assert np.isfinite(imputed.mean[1])
    assert imputed.covariance[1, 1] > 0.0


def test_mvbm_all_zero_measurement_errors_keep_the_linear_complete_case_path():
    tree = tree_from("((A:1,B:1):1,(C:1,D:1):1,E:2)R;")
    values = {
        "A": (0.0, 1.0),
        "B": (0.5, 1.8),
        "C": (1.2, 0.4),
        "D": (2.1, -0.2),
        "E": (2.8, 0.9),
    }
    expected, expected_fit = compute_mvbm_marginals(tree, values, ("x", "y"))
    zero_errors = {name: (0.0, 0.0) for name in values}
    actual, actual_fit = compute_mvbm_marginals(
        tree,
        values,
        ("x", "y"),
        standard_errors=zero_errors,
    )
    assert type(actual_fit) is type(expected_fit)
    assert actual_fit.num_effective_observations == 5
    assert actual_fit.sigma == pytest.approx(expected_fit.sigma)
    assert actual_fit.restricted_log_likelihood == pytest.approx(
        expected_fit.restricted_log_likelihood
    )
    for node in tree.traverse():
        assert actual[node].mean == pytest.approx(expected[node].mean)
        assert actual[node].covariance == pytest.approx(expected[node].covariance)


@pytest.mark.parametrize(
    "errors, message",
    [
        ({}, "required for observed tip"),
        ({"A": (0.0,)}, "dimension mismatch"),
        (False, "must be supplied as a tip mapping"),
    ],
)
def test_mvbm_validates_all_zero_measurement_error_mappings(errors, message):
    tree = tree_from("(A:1,B:1,C:1,D:1,E:1)R;")
    values = {
        "A": (0.0, 1.0),
        "B": (0.5, 1.8),
        "C": (1.2, 0.4),
        "D": (2.1, -0.2),
        "E": (2.8, 0.9),
    }
    with pytest.raises(ValueError, match=message):
        compute_mvbm_marginals(
            tree,
            values,
            ("x", "y"),
            standard_errors=errors,
        )


def test_mvbm_rejects_incomplete_or_invalid_zero_error_entries():
    tree = tree_from("(A:1,B:1,C:1,D:1,E:1)R;")
    values = {
        name: (float(index), float(index + 1)) for index, name in enumerate("ABCDE")
    }
    incomplete = {name: (0.0, 0.0) for name in "ABCD"}
    with pytest.raises(ValueError, match="required for observed tip 'E'"):
        compute_mvbm_marginals(tree, values, ("x", "y"), standard_errors=incomplete)
    invalid = {name: (0.0, 0.0) for name in "ABCDE"}
    invalid["A"] = (None, 0.0)
    with pytest.raises(ValueError, match="must be numeric and finite"):
        compute_mvbm_marginals(tree, values, ("x", "y"), standard_errors=invalid)


def test_trait_rescaling_and_offsets_transform_mvbm_results():
    tree = tree_from("((A:1,B:2):1,C:2,D:1,E:3)R;")
    values = {
        "A": (-1.0, 2.0),
        "B": (0.5, 4.0),
        "C": (3.0, -2.0),
        "D": (2.0, 5.0),
        "E": (-3.0, 1.0),
    }
    first, first_fit = compute_mvbm_marginals(tree, values, ("x", "y"))
    transform = np.diag([1e20, 0.25])
    offset = np.asarray([1e25, -10.0])
    transformed = {
        name: offset + transform @ np.asarray(value) for name, value in values.items()
    }
    second, second_fit = compute_mvbm_marginals(tree, transformed, ("x", "y"))
    np.testing.assert_allclose(
        second_fit.sigma, transform @ first_fit.sigma @ transform, rtol=1e-10
    )
    for node in tree.traverse():
        np.testing.assert_allclose(
            second[node].mean, offset + transform @ first[node].mean, rtol=2e-11
        )
        np.testing.assert_allclose(
            second[node].covariance,
            transform @ first[node].covariance @ transform,
            rtol=1e-10,
        )


@pytest.mark.integration
def test_mvbm_cli_writes_long_marginals_covariances_and_model(tmp_path):
    trait = tmp_path / "traits.tsv"
    output = tmp_path / "asr.tsv"
    covariance = tmp_path / "covariance.tsv"
    model = tmp_path / "model.tsv"
    tree_output = tmp_path / "asr.nwk"
    trait.write_text(
        "leaf_name\tx\ty\nA\t0\t0\nB\t1\t2\nC\t2\t1\nD\t3\t5\nE\tNA\tNA\nF\t6\t7\n"
    )
    main(
        [
            "asr",
            "-i",
            "[&R]((A:1,B:1):1,(C:1,D:1):1,E:2,F:2)R;",
            "--trait",
            str(trait),
            "--state-column",
            "x,y",
            "--model",
            "MV-BM",
            "--covariance-out",
            str(covariance),
            "--model-out",
            str(model),
            "--tree-out",
            str(tree_output),
            "--tree-annotation",
            "all",
            "-o",
            str(output),
        ]
    )
    result = pd.read_csv(output, sep="\t")
    covariances = pd.read_csv(covariance, sep="\t")
    metadata = pd.read_csv(model, sep="\t").iloc[0]
    assert set(result["trait"]) == {"x", "y"}
    assert len(result) == 18
    assert set(
        covariances[["trait_1", "trait_2"]].itertuples(index=False, name=None)
    ) == {
        ("x", "x"),
        ("x", "y"),
        ("y", "y"),
    }
    assert covariances["correlation"].dropna().between(-1.0, 1.0).all()
    assert metadata["model"] == "MV-BM"
    assert metadata["sigma_rank"] == 2
    assert result.loc[result["name"] == "E", "is_imputed"].all()
    restored = read_tree(str(tree_output), "auto", True, quiet=True)
    assert restored["E"].props["asr_is_imputed"] == "yes"
    assert restored["E"].props["asr_interval_kind"] == "conditional_on_covariance"
    assert float(restored["E"].props["asr_ci_level"]) == pytest.approx(0.95)


def test_mvbm_cli_supports_measurement_error_and_partial_rows(tmp_path):
    trait = tmp_path / "traits.tsv"
    trait.write_text(
        "leaf_name\tx\ty\tse_x\tse_y\n"
        "A\t0\t0\t0.1\t0.2\n"
        "B\t1\tNA\t0.1\tNA\n"
        "C\t2\t1\t0.2\t0.1\n"
    )
    base = [
        "asr",
        "-i",
        "[&R](A:1,B:1,C:1)R;",
        "--trait",
        str(trait),
        "--state-column",
        "x,y",
        "--model",
        "MV-BM",
        "-o",
        str(tmp_path / "out.tsv"),
    ]
    main([*base, "--standard-error-column", "se_x,se_y"])
    result = pd.read_csv(tmp_path / "out.tsv", sep="\t")
    assert set(result["trait"]) == {"x", "y"}
    assert result.loc[
        (result["name"] == "B") & (result["trait"] == "y"), "is_imputed"
    ].all()
    with pytest.raises(ValueError, match="one column per trait"):
        main([*base, "--standard-error-column", "se_x"])


@pytest.mark.integration
def test_mvbm_cli_warns_when_joint_covariance_is_singular(tmp_path, capsys):
    trait = tmp_path / "traits.tsv"
    output = tmp_path / "asr.tsv"
    trait.write_text("leaf_name\tx\ty\nA\t0\t0\nB\t1\t2\nC\tNA\tNA\n")
    main(
        [
            "asr",
            "-i",
            "[&R](A:1,B:1,C:1)R;",
            "--trait",
            str(trait),
            "--state-column",
            "x,y",
            "--model",
            "MV-BM",
            "-o",
            str(output),
        ]
    )
    assert "ordinary multivariate likelihood is not" in capsys.readouterr().err
