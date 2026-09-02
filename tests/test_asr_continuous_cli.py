"""Trait detection and real parser-to-output continuous ASR contracts."""

import math
from io import StringIO

import numpy as np
import pandas as pd
import pytest

from nwkit.asr import asr_main
from nwkit.cli import main, parser
from nwkit.continuous_asr import GaussianMarginal
from nwkit.continuous_asr_io import _summary
from nwkit.rooting_state import get_rooting_info
from nwkit.util import read_tree
from tests.helpers import make_args


def run_asr(tmp_path, text, *options, source="(A:1,B:1,C:1)R;", input_rooted="yes"):
    trait, output = tmp_path / "traits.tsv", tmp_path / "asr.tsv"
    trait.write_text(text)
    main(
        [
            "asr",
            "-i",
            source,
            "--input-rooted",
            input_rooted,
            "--trait",
            str(trait),
            "--state-column",
            "value",
            "-o",
            str(output),
            *options,
        ]
    )
    return pd.read_csv(output, sep="\t")


@pytest.mark.parametrize(
    "values, expected",
    [
        (["1", "2", "3"], "continuous"),
        (["-0.3", "1.4e-1", ".2"], "continuous"),
        (["red", "blue", "red"], "discrete"),
        (["0", "0|1", "1"], "discrete"),
        (["1", "category", "2"], "discrete"),
    ],
)
def test_auto_is_default_and_reports_the_selected_type(
    tmp_path, capsys, values, expected
):
    model = tmp_path / "model.tsv"
    text = "leaf_name\tvalue\n" + "".join(
        f"{name}\t{value}\n" for name, value in zip("ABC", values, strict=True)
    )
    table = run_asr(tmp_path, text, "--model-out", str(model))
    metadata = pd.read_csv(model, sep="\t").iloc[0]
    assert metadata["trait_type"] == expected
    assert metadata["trait_type_requested"] == "auto"
    stderr = capsys.readouterr().err
    assert f"ASR trait type: {expected} (auto)" in stderr
    assert ("numeric category codes" in stderr) == (expected == "continuous")
    assert ("mean" if expected == "continuous" else "map_state") in table.columns


def test_manual_discrete_preserves_numeric_category_spelling(tmp_path):
    table = run_asr(
        tmp_path,
        "leaf_name\tvalue\nA\t001\nB\t1\nC\t01\n",
        "--trait-type",
        "discrete",
        "--rate",
        "0.2",
    )
    assert {"p_001", "p_1", "p_01"}.issubset(table.columns)
    leaves = table[table.node_class == "leaf"].set_index("name")
    assert (
        leaves.loc["A", "p_001"]
        == leaves.loc["B", "p_1"]
        == leaves.loc["C", "p_01"]
        == 1
    )


def test_manual_continuous_fails_on_non_numeric_values(tmp_path):
    with pytest.raises(ValueError, match="Trait value.*numeric and finite"):
        run_asr(
            tmp_path,
            "leaf_name\tvalue\nA\t1\nB\toops\nC\t3\n",
            "--trait-type",
            "continuous",
        )
    assert not (tmp_path / "asr.tsv").exists()


@pytest.mark.parametrize("marker", ["999", "00999", "9.99e2", "0"])
@pytest.mark.parametrize("from_stdin", [False, True])
def test_numeric_se_missing_markers_are_decoded_before_conversion(
    tmp_path, monkeypatch, capsys, marker, from_stdin
):
    text = f"leaf_name\tvalue\tse\nA\t1\t0.1\nB\t2\t{marker}\nC\t3\t0.1\n"
    options = ["--standard-error-column", "se", "--missing-values", marker]
    with pytest.raises(ValueError, match="standard error is required.*'B'"):
        if from_stdin:
            monkeypatch.setattr("sys.stdin", StringIO(text))
            main(
                [
                    "asr",
                    "-i",
                    "[&R](A:1,B:1,C:1);",
                    "--trait",
                    "-",
                    "--state-column",
                    "value",
                    *options,
                ]
            )
        else:
            run_asr(tmp_path, text, *options)
    assert not (tmp_path / "asr.tsv").exists()
    assert capsys.readouterr().out == ""


def test_numeric_se_missing_markers_match_spelling_not_numeric_equality(tmp_path):
    table = run_asr(
        tmp_path,
        "leaf_name\tvalue\tse\nA\t1\t0.1\nB\t2\t999.0\nC\t3\t0.1\n",
        "--standard-error-column",
        "se",
        "--missing-values",
        "999",
        "--sigma2",
        "1",
    )
    assert table.loc[table.name == "B", "observed_se"].item() == 999.0


@pytest.mark.parametrize("bad", ["inf", "-inf", "1e999"])
def test_numeric_overflows_do_not_fall_back_to_categories(tmp_path, bad):
    with pytest.raises(ValueError, match="Trait value.*finite"):
        run_asr(tmp_path, f"leaf_name\tvalue\nA\t1\nB\t{bad}\nC\t3\n")


@pytest.mark.parametrize("bad", ["inf", "-inf", "1e999"])
def test_nonfinite_number_is_not_hidden_by_mixed_text_categories(tmp_path, bad):
    with pytest.raises(ValueError, match="Trait value.*finite"):
        run_asr(tmp_path, f"leaf_name\tvalue\nA\tred\nB\t{bad}\nC\tblue\n")


def test_explicit_discrete_mode_can_use_a_nonfinite_looking_category(tmp_path):
    table = run_asr(
        tmp_path,
        "leaf_name\tvalue\nA\tred\nB\tinf\nC\tblue\n",
        "--trait-type",
        "discrete",
        "--rate",
        "0.2",
    )
    assert "p_inf" in table.columns
    assert table.loc[table.name == "B", "p_inf"].item() == 1.0


def test_custom_missing_policy_exposes_nan_to_auto_finite_validation(tmp_path):
    text = "leaf_name\tvalue\nA\tred\nB\tnan\nC\tblue\n"
    with pytest.raises(ValueError, match="Trait value.*finite"):
        run_asr(tmp_path, text, "--missing-values", "NONE")
    table = run_asr(
        tmp_path,
        text,
        "--missing-values",
        "NONE",
        "--trait-type",
        "discrete",
        "--rate",
        "0.2",
    )
    assert table.loc[table.name == "B", "p_nan"].item() == 1.0


def test_only_matched_nonmissing_tip_values_drive_detection(tmp_path, capsys):
    table = run_asr(
        tmp_path, "leaf_name\tvalue\nA\t0\nB\t2\nC\tNA\nOTHER\ttext\n", "--sigma2", "1"
    )
    leaves = table[table.node_class == "leaf"].set_index("name")
    assert leaves.loc["C", "is_imputed"]
    assert leaves.loc["C", "mean"] == pytest.approx(1)
    assert leaves.loc["C", "variance"] == pytest.approx(1.5)
    stderr = capsys.readouterr().err
    assert stderr.count("Rows in --trait not found in tree") == 1
    assert "continuous (auto)" in stderr


def test_absent_tip_and_custom_missing_markers_are_imputed(tmp_path):
    table = run_asr(
        tmp_path,
        "leaf_name\tvalue\nA\t5\nB\tNONE\n",
        "--missing-values",
        "NONE",
        "--sigma2",
        "2",
        "--target",
        "missing-leaf",
    )
    assert set(table["name"]) == {"B", "C"}
    assert np.all(table["mean"] == 5)
    assert np.all(table["variance"] == 4)
    assert table["is_imputed"].all()


def test_all_missing_auto_needs_an_explicit_type(tmp_path):
    with pytest.raises(ValueError, match="Cannot infer --trait-type"):
        run_asr(tmp_path, "leaf_name\tvalue\nA\tNA\nB\t\nC\t?\n")


def test_all_missing_explicit_discrete_can_use_a_declared_state_space(tmp_path):
    table = run_asr(
        tmp_path,
        "leaf_name\tvalue\nA\tNA\nB\t\nC\t?\n",
        "--trait-type",
        "discrete",
        "--states",
        "0,1",
        "--rate",
        "0.2",
    )
    np.testing.assert_allclose(table[["p_0", "p_1"]], 0.5)


def test_all_missing_explicit_discrete_rejects_a_fitted_process(tmp_path):
    with pytest.raises(ValueError, match="fully fixed transition process"):
        run_asr(
            tmp_path,
            "leaf_name\tvalue\nA\tNA\nB\t\nC\t?\n",
            "--trait-type",
            "discrete",
            "--states",
            "0,1",
        )


def test_all_missing_explicit_continuous_fails_even_with_fixed_rate(tmp_path):
    with pytest.raises(ValueError, match="at least one observed"):
        run_asr(
            tmp_path,
            "leaf_name\tvalue\nA\tNA\nB\t\nC\t?\n",
            "--trait-type",
            "continuous",
            "--sigma2",
            "1",
        )


@pytest.mark.parametrize(
    "options",
    [
        ["--states", "0,1,2"],
        ["--model", "ER"],
        ["--root-prior", "equal"],
        ["--rate", "0.2"],
        ["--rate-bounds", "0.01,1"],
        ["--transition-graph", "ordered"],
        ["--rate-matrix", "q.tsv"],
        ["--output", "map"],
        ["--ambiguous-separator", "|"],
        ["--stochastic-map-out", "maps.tsv"],
        ["--n-sim", "100"],
        ["--threads", "1"],
        ["--hidden-categories", "2"],
        ["--tree-annotation", "map"],
    ],
)
def test_discrete_options_are_not_silently_ignored_in_continuous_mode(
    tmp_path, options
):
    with pytest.raises(ValueError, match="not supported for continuous"):
        run_asr(tmp_path, "leaf_name\tvalue\nA\t0\nB\t1\nC\t2\n", *options)
    assert not (tmp_path / "asr.tsv").exists()


@pytest.mark.parametrize(
    "options",
    [
        ["--sigma2", "1"],
        ["--standard-error-column", "se"],
        ["--ci-level", "0.95"],
        ["--alpha", "0.5"],
        ["--alpha-bounds", "0.01,1"],
        ["--theta", "0"],
        ["--eb-rate", "-0.2"],
        ["--eb-rate-bounds=-1,1"],
        ["--drift", "0.1"],
        ["--covariance-out", "covariance.tsv"],
        ["--model", "BM"],
        ["--root-prior", "flat"],
        ["--output", "summary"],
        ["--tree-annotation", "mean"],
    ],
)
def test_continuous_options_are_not_silently_ignored_in_discrete_mode(
    tmp_path, options
):
    with pytest.raises(ValueError, match="not supported for discrete"):
        run_asr(tmp_path, "leaf_name\tvalue\nA\tx\nB\ty\nC\tx\n", *options)
    assert not (tmp_path / "asr.tsv").exists()


def test_mode_mismatch_guidance_points_to_the_other_trait_type(tmp_path):
    with pytest.raises(ValueError, match="--trait-type discrete"):
        run_asr(
            tmp_path,
            "leaf_name\tvalue\nA\t0\nB\t1\nC\t2\n",
            "--states",
            "0,1,2",
        )
    with pytest.raises(ValueError, match="--trait-type continuous"):
        run_asr(
            tmp_path,
            "leaf_name\tvalue\nA\tx\nB\ty\nC\tx\n",
            "--sigma2",
            "1",
        )


@pytest.mark.parametrize(
    "model,text",
    [
        ("MK-REGIME", "leaf_name\tvalue\nA\tx\nB\ty\nC\tx\n"),
        ("BMS", "leaf_name\tvalue\nA\t0\nB\t1\nC\t2\n"),
        ("OUM", "leaf_name\tvalue\nA\t0\nB\t1\nC\t2\n"),
    ],
)
def test_regime_models_require_an_explicit_branch_map(tmp_path, model, text):
    with pytest.raises(ValueError, match="requires --regime-map"):
        run_asr(tmp_path, text, "--model", model)


def test_regime_map_is_not_silently_accepted_by_nonregime_model(tmp_path):
    with pytest.raises(ValueError, match="requires a regime model"):
        run_asr(
            tmp_path,
            "leaf_name\tvalue\nA\t0\nB\t1\nC\t2\n",
            "--regime-map",
            str(tmp_path / "regimes.tsv"),
        )


def test_known_se_estimates_latent_values_and_records_uncertainty_contract(tmp_path):
    model = tmp_path / "model.tsv"
    table = run_asr(
        tmp_path,
        "leaf_name\tvalue\tse\nA\t0\t1\nB\t2\t1\nC\t4\t1\n",
        "--standard-error-column",
        "se",
        "--model-out",
        str(model),
        "--ci-level",
        "0.9",
    )
    metadata = pd.read_csv(model, sep="\t").iloc[0]
    assert metadata["sigma2"] == pytest.approx(3, abs=1e-7)
    assert metadata["estimation_method"] == "REML"
    assert metadata["root_prior"] == "flat"
    assert metadata["interval_kind"] == "conditional_on_sigma2"
    assert not metadata["parameter_uncertainty_included"]
    assert not metadata["tree_uncertainty_included"]
    a = table[table.name == "A"].iloc[0]
    assert a.observed_value == 0
    assert a.observed_se == 1
    assert 0 < a["mean"] < 2
    assert a.variance > 0
    assert not a.is_imputed
    np.testing.assert_allclose(table.sd**2, table.variance)
    np.testing.assert_allclose((table.ci_upper + table.ci_lower) / 2, table["mean"])
    assert (table.ci_level == 0.9).all()


@pytest.mark.parametrize(
    "text, column, message",
    [
        (
            "leaf_name\tvalue\nA\t1\nB\t2\nC\t3\n",
            "se",
            "Missing required standard-error column",
        ),
        (
            "leaf_name\tvalue\tse\nA\t1\t1\nB\t2\tNA\nC\t3\t0\n",
            "se",
            "standard error is required",
        ),
        (
            "leaf_name\tvalue\tse\nA\t1\t-1\nB\t2\t0\nC\t3\t0\n",
            "se",
            "Standard error.*non-negative",
        ),
        ("leaf_name\tvalue\nA\t1\nB\t2\nC\t3\n", "value", "must differ"),
    ],
)
def test_measurement_error_column_validation(tmp_path, text, column, message):
    with pytest.raises(ValueError, match=message):
        run_asr(tmp_path, text, "--standard-error-column", column)


@pytest.mark.parametrize("level", ["0", "1", "-0.1", "1.1"])
def test_interval_level_must_be_strictly_between_zero_and_one(tmp_path, level):
    with pytest.raises(ValueError, match="--ci-level"):
        run_asr(tmp_path, "leaf_name\tvalue\nA\t1\nB\t2\nC\t3\n", "--ci-level", level)


def test_tiny_positive_interval_level_retains_a_positive_width():
    level = 1e-20
    summary = _summary(GaussianMarginal(0.0, 1.0), level)
    expected_width = level * math.sqrt(math.pi / 2.0)
    assert summary["ci_lower"] == pytest.approx(-expected_width, rel=1e-15)
    assert summary["ci_upper"] == pytest.approx(expected_width, rel=1e-15)


@pytest.mark.parametrize(
    "level", [math.ulp(0.0), 1e-300, 1e-16, math.nextafter(1.0, 0.0)]
)
def test_interval_width_stays_positive_and_finite_across_valid_float_levels(level):
    summary = _summary(GaussianMarginal(0.0, 1.0), level)
    assert 0.0 < summary["ci_upper"] < math.inf
    assert summary["ci_lower"] == -summary["ci_upper"]


@pytest.mark.parametrize("annotation", [None, "mean", "summary", "all"])
def test_nhx_roundtrip_preserves_rooting_numeric_names_and_uncertainties(
    tmp_path, annotation
):
    tree_output = tmp_path / "ancestral.nwk"
    options = ["--sigma2", "1", "--tree-out", str(tree_output)]
    if annotation is not None:
        options += ["--tree-annotation", annotation]
    table = run_asr(
        tmp_path,
        "leaf_name\tvalue\nA\t0\nB\t2\nC\tNA\n",
        *options,
        source="((A:1,B:1)'95':1,C:2)'20';",
    )
    restored = read_tree(str(tree_output), "auto", True, quiet=True)
    assert get_rooting_info(restored).state == "rooted"
    assert {node.name for node in restored.traverse() if not node.is_leaf} == {
        "20",
        "95",
    }
    for node in restored.traverse():
        row = table[table.name.astype(str) == str(node.name)].iloc[0]
        assert float(node.props["asr_mean"]) == pytest.approx(row["mean"], rel=1e-5)
        assert node.props["asr_trait_type"] == "continuous"
        assert "asr_model" not in node.props
        if annotation == "mean":
            assert "asr_sd" not in node.props
        else:
            assert float(node.props["asr_variance"]) == pytest.approx(
                row.variance, rel=1e-5
            )
            assert node.props["asr_interval_kind"] == "conditional_on_sigma2"
    if annotation == "all":
        assert restored["C"].props["asr_is_imputed"] == "yes"


def test_zero_boundary_has_explicit_metadata_and_warning(tmp_path, capsys):
    model = tmp_path / "model.tsv"
    table = run_asr(
        tmp_path, "leaf_name\tvalue\nA\t3\nB\t3\nC\t3\n", "--model-out", str(model)
    )
    metadata = pd.read_csv(model, sep="\t").iloc[0]
    assert metadata.sigma2 == 0
    assert metadata.fit_status == "singular_zero_boundary"
    assert pd.isna(metadata.restricted_log_likelihood)
    assert (table.variance == 0).all()
    assert "exclude rate-estimation uncertainty" in capsys.readouterr().err


def test_fixed_zero_rate_does_not_warn_about_rate_estimation(tmp_path, capsys):
    model = tmp_path / "model.tsv"
    run_asr(
        tmp_path,
        "leaf_name\tvalue\nA\t3\nB\t3\nC\t3\n",
        "--sigma2",
        "0",
        "--model-out",
        str(model),
    )
    metadata = pd.read_csv(model, sep="\t").iloc[0]
    assert not metadata.sigma2_estimated
    assert metadata.fit_status == "singular_zero_boundary"
    assert "rate-estimation uncertainty" not in capsys.readouterr().err


def test_empty_target_selection_still_has_a_header(tmp_path):
    table = run_asr(
        tmp_path, "leaf_name\tvalue\nA\t1\nB\t2\nC\t3\n", "--target", "missing-leaf"
    )
    assert table.empty
    assert {"branch_id", "mean", "variance", "ci_level"}.issubset(table.columns)


def test_primary_stdout_contains_only_the_selected_table(tmp_path, capsys):
    trait = tmp_path / "traits.tsv"
    trait.write_text("leaf_name\tvalue\nA\t0\nB\t2\n")
    main(["asr", "-i", "(A:1,B:1);", "--trait", str(trait), "--state-column", "value"])
    output = capsys.readouterr()
    table = pd.read_csv(StringIO(output.out), sep="\t")
    assert len(table) == 3
    assert "continuous (auto)" in output.err


def test_trait_stdin_is_read_once_for_detection_and_inference(monkeypatch, tmp_path):
    monkeypatch.setattr("sys.stdin", StringIO("leaf_name\tvalue\nA\t0\nB\t2\n"))
    output = tmp_path / "out.tsv"
    main(
        [
            "asr",
            "-i",
            "(A:1,B:1);",
            "--trait",
            "-",
            "--state-column",
            "value",
            "-o",
            str(output),
        ]
    )
    table = pd.read_csv(output, sep="\t")
    assert table[table.node_class == "root"]["mean"].iloc[0] == 1


def test_programmatic_arguments_default_to_auto(tmp_path):
    trait = tmp_path / "traits.tsv"
    trait.write_text("leaf_name\tvalue\nA\t0\nB\t2\n")
    output = tmp_path / "out.tsv"
    asr_main(
        make_args(
            infile="(A:1,B:1);",
            trait=str(trait),
            state_column="value",
            outfile=str(output),
        )
    )
    assert "mean" in pd.read_csv(output, sep="\t")


@pytest.mark.parametrize("input_rooted", ["auto", "no"])
def test_rooting_validation_is_not_skipped_for_continuous_traits(
    tmp_path, input_rooted
):
    with pytest.raises(ValueError, match="ASR requires a rooted tree"):
        run_asr(
            tmp_path, "leaf_name\tvalue\nA\t0\nB\t2\nC\t3\n", input_rooted=input_rooted
        )


def test_cli_parses_auto_and_type_specific_defaults_without_selecting_a_model():
    args = parser.parse_args(
        ["asr", "--trait", "traits.tsv", "--state-column", "value"]
    )
    assert args.trait_type == "auto"
    assert args.model is None
    assert args.root_prior is None
    assert args.output is None
