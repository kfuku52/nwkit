import re
from types import SimpleNamespace

import pandas as pd
import pytest

from nwkit.asr_compare import (
    ComparisonCandidate,
    ComparisonContext,
    _candidate_args,
    _classify_fit,
    _comparison_group,
    _preflight_comparison_figure,
    _successful_row,
    _validate_comparison_options,
    comparison_table,
    evaluate_comparison_candidates,
    resolve_comparison_candidates,
)
from nwkit.asr_compare_figure import (
    _draw_candidate_overview,
    _font_family_for_text,
    _format_number,
    _row_note,
    draw_comparison_figure,
)
from nwkit.asr_comparison import (
    grouped_model_comparison_table,
    model_comparison_table,
    summarize_fit,
)
from nwkit.cli import main


def test_continuous_flat_root_comparison_calculates_weights():
    first = SimpleNamespace(
        sigma2_estimated=True,
        restricted_log_likelihood=-10.0,
        num_effective_observations=20,
        fit_status="ok",
    )
    second = SimpleNamespace(
        sigma2_estimated=True,
        evolution_parameter_estimated=True,
        restricted_log_likelihood=-8.0,
        num_effective_observations=20,
        fit_status="ok",
    )
    table = model_comparison_table(
        [
            summarize_fit("BM", first, trait_type="continuous"),
            summarize_fit("LAMBDA", second, trait_type="continuous"),
        ]
    )
    assert table.iloc[0]["model"] == "LAMBDA"
    assert table["aic_weight"].sum() == pytest.approx(1.0)
    assert table["bic_weight"].sum() == pytest.approx(1.0)


def test_singular_flat_root_fit_keeps_flat_likelihood_contract():
    fit = SimpleNamespace(
        sigma2_estimated=True,
        restricted_log_likelihood=None,
        num_effective_observations=5,
        fit_status="singular_zero_boundary",
    )
    summary = summarize_fit("BM", fit, trait_type="continuous")
    assert summary["likelihood_kind"] == "flat_root_integrated"
    assert summary["log_likelihood"] is None


def test_fitted_model_without_a_likelihood_is_unassigned_from_comparison_sets():
    context = ComparisonContext(
        tree=None,
        trait_df=pd.DataFrame(),
        trait_type="continuous",
        trait_columns=("x", "y"),
        error_columns=None,
        args=SimpleNamespace(trait_type="continuous"),
    )
    fit = SimpleNamespace(
        model="MV-BM",
        trait_names=("x", "y"),
        sigma_estimated=True,
        restricted_log_likelihood=None,
        num_effective_observations=5,
        fit_status="singular_covariance",
        optimizer_success=True,
        optimizer_message="",
    )
    row = _successful_row(
        context,
        ComparisonCandidate("MV-BM", "flat", "model-default"),
        fit,
        0.1,
    )
    assert row["status"] == "no_likelihood"
    assert row["rankable"] == "no"
    assert row["comparison_group"] == ""


def test_comparison_rejects_mixed_likelihood_conventions():
    with pytest.raises(ValueError, match="different likelihood conventions"):
        model_comparison_table(
            [
                {
                    "model": "BM",
                    "likelihood_kind": "flat_root_integrated",
                    "log_likelihood": -1.0,
                    "num_parameters": 1,
                    "sample_size": 10,
                    "fit_status": "ok",
                },
                {
                    "model": "OU",
                    "likelihood_kind": "stationary_root_ml",
                    "log_likelihood": -1.0,
                    "num_parameters": 3,
                    "sample_size": 10,
                    "fit_status": "ok",
                },
            ]
        )


@pytest.mark.parametrize(
    "fit, expected",
    [
        (
            SimpleNamespace(
                sigma2_by_regime={"a": 1.0, "b": 2.0},
                sigma2_estimated=True,
                drift_by_regime={"a": 0.0, "b": 0.1},
                drift_estimated=True,
                restricted_log_likelihood=-3.0,
                num_effective_observations=20,
            ),
            4,
        ),
        (
            SimpleNamespace(
                model="OUMVA",
                alpha=None,
                alpha_by_regime={"a": 0.5, "b": 0.8},
                alpha_estimated=True,
                theta_by_regime={"a": 0.0, "b": 1.0},
                theta_estimated=True,
                sigma2=None,
                sigma2_by_regime={"a": 1.0, "b": 2.0},
                sigma2_estimated=True,
                log_likelihood=-3.0,
                num_effective_observations=20,
                root_prior="stationary",
            ),
            6,
        ),
        (
            SimpleNamespace(
                model="OUM",
                alpha=0.5,
                alpha_by_regime={"a": 0.5, "b": 0.5},
                alpha_estimated=True,
                theta_by_regime={"a": 0.0, "b": 1.0},
                theta_estimated=True,
                sigma2=1.0,
                sigma2_by_regime={"a": 1.0, "b": 1.0},
                sigma2_estimated=True,
                log_likelihood=-3.0,
                num_effective_observations=20,
                root_prior="stationary",
            ),
            4,
        ),
    ],
)
def test_regime_parameter_counts_include_each_free_map_value(fit, expected):
    summary = summarize_fit("candidate", fit, trait_type="continuous")
    assert summary["num_parameters"] == expected


@pytest.mark.parametrize(
    "model,alpha_estimated,theta_estimated,expected",
    [
        ("MV-BM", False, False, 3),
        ("MV-OU", True, True, 6),
        ("MV-OU", False, True, 5),
    ],
)
def test_dense_multivariate_parameter_counts(
    model, alpha_estimated, theta_estimated, expected
):
    fit = SimpleNamespace(
        model=model,
        trait_names=("x", "y"),
        alpha_estimated=alpha_estimated,
        theta_estimated=theta_estimated,
        restricted_log_likelihood=-3.0 if model == "MV-BM" else None,
        log_likelihood=-3.0 if model == "MV-OU" else None,
        num_effective_observations=20,
        root_prior="stationary",
    )
    summary = summarize_fit(model, fit, trait_type="continuous")
    assert summary["num_parameters"] == expected


@pytest.mark.parametrize(
    "field, value, message",
    [
        ("sample_size", 0, "non-positive sample size"),
        ("num_parameters", -1, "negative parameter count"),
    ],
)
def test_comparison_rejects_invalid_dimensions(field, value, message):
    rows = [
        {
            "model": model,
            "likelihood_kind": "discrete_ml",
            "log_likelihood": -1.0,
            "num_parameters": 1,
            "sample_size": 10,
            "fit_status": "ok",
        }
        for model in ("ER", "ARD")
    ]
    rows[0][field] = value
    with pytest.raises(ValueError, match=message):
        model_comparison_table(rows)


@pytest.mark.parametrize(
    "status",
    [
        "singular_support",
        "hidden_rate_homogeneous_boundary",
        "hidden_rate_spread_upper_boundary",
        "base_rate_upper_boundary",
        "switching_rate_lower_boundary",
        "effective_rate_lower_boundary",
        "effective_rate_upper_boundary",
    ],
)
def test_comparison_rejects_singular_or_covarion_boundary_status(status):
    rows = [
        {
            "model": model,
            "likelihood_kind": "discrete_ml",
            "log_likelihood": -1.0,
            "num_parameters": 1,
            "sample_size": 10,
            "fit_status": status if model == "COVARION" else "ok",
        }
        for model in ("COVARION", "second")
    ]
    with pytest.raises(ValueError, match="non-regular fit status"):
        model_comparison_table(rows)


def test_comparison_retains_boundary_status_for_other_regular_model_families():
    rows = [
        {
            "model": model,
            "likelihood_kind": "discrete_ml",
            "log_likelihood": -1.0,
            "num_parameters": 1,
            "sample_size": 10,
            "fit_status": "rate_upper_boundary" if model == "ER" else "ok",
        }
        for model in ("ER", "ARD")
    ]
    table = model_comparison_table(rows)
    assert set(table["fit_status"]) == {"rate_upper_boundary", "ok"}


@pytest.mark.parametrize(
    "status", ["sigma2_lower_boundary", "root_variance_lower_boundary"]
)
def test_zero_variance_boundaries_are_not_rankable(status):
    classification, rankable, message = _classify_fit(
        ComparisonCandidate("BMS", "flat", "model-default"),
        SimpleNamespace(optimizer_success=True),
        {"fit_status": status, "log_likelihood": -1.0},
    )
    assert classification == "nonregular"
    assert rankable == "no"
    assert "variance component reached zero" in message

    rows = [
        {
            "model": model,
            "likelihood_kind": "flat_root_integrated",
            "log_likelihood": -1.0,
            "num_parameters": 1,
            "sample_size": 10,
            "fit_status": status if model == "BMS" else "ok",
        }
        for model in ("BM", "BMS")
    ]
    with pytest.raises(ValueError, match="non-regular fit status"):
        model_comparison_table(rows)


def test_grouped_comparison_never_ranks_different_root_likelihoods_together():
    rows = []
    for model, group, likelihood, count in (
        ("BM", "continuous:scalar:flat_root_integrated", -10.0, 1),
        ("LAMBDA", "continuous:scalar:flat_root_integrated", -9.0, 2),
        ("OU", "continuous:scalar:stationary_root_ml", -8.0, 3),
    ):
        rows.append(
            {
                "model": model,
                "comparison_group": group,
                "rankable": "yes",
                "log_likelihood": likelihood,
                "num_parameters": count,
                "sample_size": 20,
            }
        )
    table = grouped_model_comparison_table(rows, criterion="aic")
    flat = table[table["comparison_group"].str.contains("flat_root")]
    stationary = table[table["comparison_group"].str.contains("stationary_root")]
    assert flat["aic_weight"].sum() == pytest.approx(1.0)
    assert set(flat["criterion_rank"]) == {1}
    assert set(flat["is_best"]) == {"yes"}
    assert stationary.iloc[0]["num_comparable_models"] == 1
    assert pd.isna(stationary.iloc[0]["criterion_rank"])
    assert pd.isna(stationary.iloc[0]["aic_weight"])


def test_grouped_comparison_retains_failed_rows_without_metrics():
    table = grouped_model_comparison_table(
        [
            {
                "model": "BM",
                "comparison_group": "flat",
                "rankable": "yes",
                "log_likelihood": -2.0,
                "num_parameters": 1,
                "sample_size": 10,
            },
            {
                "model": "KAPPA",
                "comparison_group": "",
                "rankable": "no",
                "status": "failed",
                "log_likelihood": None,
                "num_parameters": None,
                "sample_size": None,
            },
        ]
    )
    failed = table[table["model"] == "KAPPA"].iloc[0]
    assert failed["status"] == "failed"
    assert pd.isna(failed["aic"])


def test_grouped_comparison_assigns_dense_ranks_and_joint_winners_to_ties():
    rows = [
        {
            "model": model,
            "comparison_group": "same",
            "rankable": "yes",
            "log_likelihood": likelihood,
            "num_parameters": count,
            "sample_size": 20,
        }
        for model, likelihood, count in (
            ("first", -3.0, 1),
            ("second", -3.0, 1),
            ("third", -3.0, 2),
        )
    ]
    table = grouped_model_comparison_table(rows, criterion="aic")
    tied = table.iloc[:2]
    assert set(tied["criterion_rank"]) == {1}
    assert set(tied["is_best"]) == {"yes"}
    assert set(tied["delta_aic"]) == {0.0}
    assert table.iloc[2]["criterion_rank"] == 2


def test_dense_rank_ties_do_not_chain_across_the_absolute_tolerance():
    rows = [
        {
            "model": model,
            "comparison_group": "same",
            "rankable": "yes",
            "log_likelihood": -value / 2.0,
            "num_parameters": 0,
            "sample_size": 20,
        }
        for model, value in (
            ("first", 0.0),
            ("second", 0.75e-9),
            ("third", 1.5e-9),
        )
    ]
    table = grouped_model_comparison_table(rows, criterion="aic")
    assert list(table["criterion_rank"]) == [1, 1, 2]
    assert list(table["is_best"]) == ["yes", "yes", "no"]
    assert table.iloc[2]["delta_aic"] == pytest.approx(1.5e-9)


def test_information_criterion_ties_are_translation_invariant():
    rows = [
        {
            "model": model,
            "comparison_group": "same",
            "rankable": "yes",
            "log_likelihood": -value / 2.0,
            "num_parameters": 0,
            "sample_size": 20,
        }
        for model, value in (("first", 1e12), ("second", 1e12 + 0.5))
    ]
    table = grouped_model_comparison_table(rows, criterion="aic")
    assert list(table["criterion_rank"]) == [1, 2]
    assert list(table["is_best"]) == ["yes", "no"]
    assert table.iloc[1]["delta_aic"] == pytest.approx(0.5)


def test_selected_criterion_count_excludes_nonfinite_aicc_rows():
    rows = [
        {
            "model": model,
            "comparison_group": "same",
            "rankable": "yes",
            "log_likelihood": -3.0,
            "num_parameters": count,
            "sample_size": sample_size,
        }
        for model, count, sample_size in (
            ("first", 1, 10),
            ("second", 2, 10),
            ("undefined", 2, 3),
        )
    ]
    table = grouped_model_comparison_table(rows, criterion="aicc")
    assert set(table["num_comparable_models"]) == {2}
    assert pd.isna(table.iloc[2]["criterion_rank"])


def test_comparison_table_serializes_nullable_counts_as_integers(monkeypatch):
    rows = [
        {
            "model_id": "BM",
            "model": "BM",
            "comparison_group": "same",
            "rankable": "yes",
            "sample_size": 20,
            "num_parameters": 1,
            "log_likelihood": -3.0,
        },
        {
            "model_id": "failed",
            "model": "failed",
            "comparison_group": "",
            "rankable": "no",
            "sample_size": float("nan"),
            "num_parameters": float("nan"),
            "log_likelihood": float("nan"),
        },
    ]
    monkeypatch.setattr(
        "nwkit.asr_compare.evaluate_comparison_candidates",
        lambda *_args, **_kwargs: rows,
    )
    context = ComparisonContext(
        tree=None,
        trait_df=pd.DataFrame(),
        trait_type="continuous",
        trait_columns=("x",),
        error_columns=None,
        args=SimpleNamespace(),
    )
    table = comparison_table(context, [], automatic=True, criterion="aic")
    assert str(table["sample_size"].dtype) == "Int64"
    assert str(table["num_parameters"].dtype) == "Int64"
    serialized = table.to_csv(sep="\t", index=False, na_rep="")
    assert "\t20\t1\t" in serialized
    assert "\t20.0\t1.0\t" not in serialized


def test_pdf_note_preserves_diagnostic_and_explains_only_the_current_row():
    record = SimpleNamespace(
        message="Boundary fit retained: alpha_upper_boundary.",
        criterion_rank=float("nan"),
        rankable="yes",
        num_comparable_models=2,
        aicc=float("nan"),
    )
    note = _row_note(record, "aicc", "AICc")
    assert "Boundary fit retained" in note
    assert "AICc is unavailable" in note
    assert "Fewer than two" not in note


def test_candidate_overview_draws_every_ranked_and_unranked_model():
    table = pd.DataFrame(
        [
            {
                "model_id": "BM",
                "trait_type": "continuous",
                "trait_columns": "x",
                "status": "ok",
                "rankable": "yes",
                "message": "",
                "comparison_group": "continuous:scalar:flat_root_integrated",
                "criterion_rank": 1,
                "num_comparable_models": 2,
                "is_best": "yes",
                "sample_size": 20,
                "num_parameters": 1,
                "log_likelihood": -10.0,
                "aicc": 22.2,
                "aicc_weight": 0.75,
                "delta_aicc": 0.0,
            },
            {
                "model_id": "LAMBDA",
                "trait_type": "continuous",
                "trait_columns": "x",
                "status": "boundary",
                "rankable": "yes",
                "message": "Boundary fit retained.",
                "comparison_group": "continuous:scalar:flat_root_integrated",
                "criterion_rank": 2,
                "num_comparable_models": 2,
                "is_best": "no",
                "sample_size": 20,
                "num_parameters": 2,
                "log_likelihood": -10.1,
                "aicc": 24.4,
                "aicc_weight": 0.25,
                "delta_aicc": 2.2,
            },
            {
                "model_id": "OU[stationary]",
                "trait_type": "continuous",
                "trait_columns": "x",
                "status": "ok",
                "rankable": "yes",
                "message": "",
                "comparison_group": "continuous:scalar:stationary_root_ml",
                "criterion_rank": float("nan"),
                "num_comparable_models": 1,
                "is_best": "",
                "sample_size": 20,
                "num_parameters": 3,
                "log_likelihood": -9.0,
                "aicc": 25.5,
                "aicc_weight": float("nan"),
                "delta_aicc": float("nan"),
            },
            {
                "model_id": "BMS",
                "trait_type": "continuous",
                "trait_columns": "x",
                "status": "not_applicable",
                "rankable": "no",
                "message": "the model requires --regime-map",
                "comparison_group": "",
                "criterion_rank": float("nan"),
                "num_comparable_models": 0,
                "is_best": "",
                "sample_size": float("nan"),
                "num_parameters": float("nan"),
                "log_likelihood": float("nan"),
                "aicc": float("nan"),
                "aicc_weight": float("nan"),
                "delta_aicc": float("nan"),
            },
        ]
    )

    class CapturePdf:
        def __init__(self):
            self.figures = []

        def savefig(self, figure):
            self.figures.append(figure)

    table = table.iloc[[1, 0, 2, 3]].reset_index(drop=True)
    pdf = CapturePdf()
    _draw_candidate_overview(pdf, table, "aicc")
    assert len(pdf.figures) == 1
    texts = [text.get_text() for text in pdf.figures[0].axes[0].texts]
    assert {"BM", "LAMBDA", "OU[stationary]", "BMS"} <= set(texts)
    assert "#1" in texts
    assert "0.750" in texts
    assert "0.000" in texts
    assert "the model requires --regime-map" in texts
    positions = {
        text.get_text(): text.get_position()[1]
        for text in pdf.figures[0].axes[0].texts
        if text.get_text() in {"BM", "LAMBDA"}
    }
    assert positions["BM"] > positions["LAMBDA"]
    assert any(
        text.startswith("Comparison set 1: Continuous | single trait") for text in texts
    )


def test_candidate_overview_bounds_long_text_and_aligns_numeric_headers():
    long_name = ",".join(f"trait_{index}_{'x' * 35}" for index in range(20))
    table = pd.DataFrame(
        [
            {
                "model_id": "MODEL-WITH-A-VERY-LONG-DISPLAY-NAME",
                "trait_type": "continuous",
                "trait_columns": long_name,
                "status": "failed",
                "rankable": "no",
                "message": "diagnostic " * 120,
                "comparison_group": "",
                "criterion_rank": float("nan"),
                "num_comparable_models": 0,
                "is_best": "",
                "sample_size": 20,
                "num_parameters": 1,
                "log_likelihood": 1e100,
                "aic": 1e100,
                "aic_weight": float("nan"),
                "delta_aic": float("nan"),
            }
        ]
    )

    class CapturePdf:
        def __init__(self):
            self.figures = []

        def savefig(self, figure):
            self.figures.append(figure)

    pdf = CapturePdf()
    _draw_candidate_overview(pdf, table, "aic")
    figure = pdf.figures[0]
    assert len(figure._suptitle.get_text()) <= 110
    notes = [
        text
        for text in figure.axes[0].texts
        if text.get_position()[0] == pytest.approx(0.72) and text.get_text() != "Notes"
    ]
    assert notes
    assert all(len(text.get_text()) <= 76 for text in notes)
    assert _format_number(1e100) == "1.000e+100"

    from matplotlib.backends.backend_agg import FigureCanvasAgg

    canvas = FigureCanvasAgg(figure)
    canvas.draw()
    renderer = canvas.get_renderer()
    width = figure.bbox.width
    for text in [*figure.texts, *figure.axes[0].texts]:
        bounds = text.get_window_extent(renderer=renderer)
        assert bounds.x0 >= -1.0
        assert bounds.x1 <= width + 1.0


def test_comparison_pdf_selects_a_font_covering_japanese_text(tmp_path):
    try:
        family = _font_family_for_text("日本語の連続形質")
    except ValueError:
        pytest.skip("No installed Japanese font is available on this test host.")
    assert family
    table = pd.DataFrame(
        [
            {
                "model_id": "BM",
                "trait_type": "continuous",
                "trait_columns": "日本語の連続形質",
                "status": "ok",
                "rankable": "yes",
                "message": "",
                "comparison_group": "continuous:scalar:flat_root_integrated:flat",
                "criterion_rank": float("nan"),
                "num_comparable_models": 1,
                "is_best": "",
                "sample_size": 20,
                "num_parameters": 1,
                "log_likelihood": -10.0,
                "aic": 22.0,
                "aic_weight": float("nan"),
                "delta_aic": float("nan"),
            }
        ]
    )
    path = tmp_path / "japanese.pdf"
    draw_comparison_figure(table, path, "aic")
    assert path.read_bytes().startswith(b"%PDF")


def test_font_error_reports_the_actually_missing_codepoint(monkeypatch):
    from matplotlib import font_manager, ft2font

    entry = SimpleNamespace(name="ASCII Test", fname="ascii-test.ttf")
    monkeypatch.setattr(font_manager.fontManager, "ttflist", [entry])

    class AsciiFont:
        def get_charmap(self):
            return {codepoint: codepoint for codepoint in range(128)}

    monkeypatch.setattr(ft2font, "FT2Font", lambda _path: AsciiFont())
    with pytest.raises(ValueError) as error:
        _font_family_for_text("A🐍")
    assert "U+1F40D" in str(error.value)
    assert "U+0041" not in str(error.value)


def test_figure_font_preflight_runs_before_model_fitting(monkeypatch):
    context = ComparisonContext(
        tree=None,
        trait_df=pd.DataFrame(),
        trait_type="continuous",
        trait_columns=("trait🐍",),
        error_columns=None,
        args=SimpleNamespace(figure_out="comparison.pdf"),
    )

    def reject_font(text):
        assert "trait🐍" in text
        raise ValueError("missing U+1F40D")

    monkeypatch.setattr("nwkit.asr_compare._font_family_for_text", reject_font)
    with pytest.raises(ValueError, match=r"U\+1F40D"):
        _preflight_comparison_figure(
            context,
            [ComparisonCandidate("BM", "flat", "model-default")],
            "aic",
        )


def test_candidate_overview_numeric_headers_share_value_alignment():
    table = pd.DataFrame(
        [
            {
                "model_id": "BM",
                "trait_type": "continuous",
                "trait_columns": "x",
                "status": "ok",
                "rankable": "yes",
                "message": "",
                "comparison_group": "continuous:scalar:flat_root_integrated:flat",
                "criterion_rank": 1,
                "num_comparable_models": 2,
                "is_best": "yes",
                "sample_size": 20,
                "num_parameters": 1,
                "log_likelihood": -10.0,
                "aic": 22.0,
                "aic_weight": 0.75,
                "delta_aic": 0.0,
            }
        ]
    )

    class CapturePdf:
        def __init__(self):
            self.figure = None

        def savefig(self, figure):
            self.figure = figure

    pdf = CapturePdf()
    _draw_candidate_overview(pdf, table, "aic")
    texts = pdf.figure.axes[0].texts
    expected = {
        "n": "20",
        "k": "1",
        "logL": "-10.000",
        "AIC": "22.000",
        "Delta AIC": "0.000",
        "Weight": "0.750",
    }
    for header, value in expected.items():
        header_text = next(text for text in texts if text.get_text() == header)
        value_text = next(text for text in texts if text.get_text() == value)
        assert header_text.get_position()[0] == value_text.get_position()[0]
        assert header_text.get_ha() == value_text.get_ha() == "right"


def test_candidate_resolution_supports_ou_root_variants_and_exclusions():
    args = SimpleNamespace(
        models="BM,OU[stationary],OU[gaussian]",
        exclude_models="OU[stationary]",
        root_prior=None,
    )
    candidates, automatic = resolve_comparison_candidates("continuous", args)
    assert not automatic
    assert [candidate.model_id for candidate in candidates] == [
        "BM",
        "OU[gaussian]",
    ]
    assert candidates[-1].root_prior_source == "model-token"


def test_discrete_comparison_groups_include_the_root_contract():
    context = ComparisonContext(
        tree=None,
        trait_df=pd.DataFrame(),
        trait_type="discrete",
        trait_columns=("x",),
        error_columns=None,
        args=SimpleNamespace(),
    )
    equal = _comparison_group(
        context,
        ComparisonCandidate("ER", "equal", "model-default"),
        "discrete_ml",
    )
    stationary = _comparison_group(
        context,
        ComparisonCandidate("F81", "stationary", "model-default"),
        "discrete_ml",
    )
    assert equal.endswith(":equal")
    assert stationary.endswith(":stationary")
    assert equal != stationary


@pytest.mark.parametrize(
    "candidate, option, value",
    [
        (ComparisonCandidate("BM", "flat", "model-default"), "alpha", 1.0),
        (
            ComparisonCandidate("MV-OU", "stationary", "model-default"),
            "theta",
            2.0,
        ),
    ],
)
def test_comparison_rejects_options_unused_by_every_selected_model(
    candidate, option, value
):
    args = SimpleNamespace(root_prior=None, **{option: value})
    context = ComparisonContext(
        tree=None,
        trait_df=pd.DataFrame(),
        trait_type="continuous",
        trait_columns=("x", "y") if candidate.model == "MV-OU" else ("x",),
        error_columns=None,
        args=args,
    )
    with pytest.raises(ValueError, match="not used by any selected"):
        _validate_comparison_options(context, [candidate])


def test_comparison_routes_model_specific_option_only_to_consuming_candidates():
    args = SimpleNamespace(root_prior=None, alpha=1.0)
    context = ComparisonContext(
        tree=None,
        trait_df=pd.DataFrame(),
        trait_type="continuous",
        trait_columns=("x",),
        error_columns=None,
        args=args,
    )
    _validate_comparison_options(
        context,
        [
            ComparisonCandidate("BM", "flat", "model-default"),
            ComparisonCandidate("OU", "stationary", "model-default"),
        ],
    )


def test_discrete_candidate_arguments_do_not_leak_fixed_matrix_options():
    args = SimpleNamespace(
        rate=0.5,
        rate_bounds="0.1,2",
        rate_matrix="q.tsv",
        transition_graph="ordered",
    )
    er_args = _candidate_args(args, ComparisonCandidate("ER", "equal", "model-default"))
    custom_args = _candidate_args(
        args, ComparisonCandidate("CUSTOM", "equal", "model-default")
    )
    assert er_args.rate == 0.5
    assert er_args.rate_bounds == "0.1,2"
    assert er_args.transition_graph == "ordered"
    assert er_args.rate_matrix is None
    assert custom_args.rate is None
    assert custom_args.rate_bounds is None
    assert custom_args.transition_graph is None
    assert custom_args.rate_matrix == "q.tsv"


def test_automatic_skipped_hrm_does_not_consume_hidden_category_option():
    context = ComparisonContext(
        tree=None,
        trait_df=pd.DataFrame(),
        trait_type="discrete",
        trait_columns=("x",),
        error_columns=None,
        args=SimpleNamespace(root_prior=None, hidden_categories=3),
    )
    with pytest.raises(ValueError, match="not used by any selected applicable"):
        _validate_comparison_options(
            context,
            [ComparisonCandidate("HRM", "equal", "model-default")],
            automatic=True,
        )


def test_comparison_rejects_fixed_transform_value_invalid_for_any_consumer():
    args = SimpleNamespace(root_prior=None, evolution_parameter=0.5)
    context = ComparisonContext(
        tree=None,
        trait_df=pd.DataFrame(),
        trait_type="continuous",
        trait_columns=("x",),
        error_columns=None,
        args=args,
    )
    with pytest.raises(ValueError, match="invalid for selected model EB"):
        _validate_comparison_options(
            context,
            [
                ComparisonCandidate("LAMBDA", "flat", "model-default"),
                ComparisonCandidate("EB", "flat", "model-default"),
            ],
        )


def test_automatic_candidate_resolution_uses_each_model_default_root():
    args = SimpleNamespace(
        models="all",
        exclude_models="MV-BM",
        root_prior=None,
    )
    candidates, automatic = resolve_comparison_candidates("continuous", args)
    by_model = {candidate.model: candidate for candidate in candidates}
    assert automatic
    assert "MV-BM" not in by_model
    assert by_model["BM"].root_prior == "flat"
    assert by_model["OU"].root_prior == "stationary"


def test_automatic_candidate_failures_are_isolated(monkeypatch):
    args = SimpleNamespace(
        trait="traits.tsv",
        states=None,
        regime_map=None,
        rate_matrix=None,
    )
    context = ComparisonContext(
        tree=None,
        trait_df=pd.DataFrame(),
        trait_type="discrete",
        trait_columns=("x",),
        error_columns=None,
        args=args,
    )
    candidates = [
        ComparisonCandidate("ER", "equal", "model-default"),
        ComparisonCandidate("ARD", "equal", "model-default"),
    ]

    def fake_fit(_context, candidate):
        if candidate.model == "ER":
            raise ValueError("synthetic failure")
        return {
            "model": "ARD",
            "rate_estimated": True,
            "rates": [0.2, 0.3],
            "log_likelihood": -3.0,
            "sample_size": 5,
            "fit_status": "ok",
            "optimizer_success": True,
        }

    monkeypatch.setattr("nwkit.asr_compare._fit_candidate", fake_fit)
    rows = evaluate_comparison_candidates(context, candidates, automatic=True)
    assert [row["status"] for row in rows] == ["failed", "ok"]
    with pytest.raises(ValueError, match="synthetic failure"):
        evaluate_comparison_candidates(context, candidates, automatic=False)


def test_automatic_comparison_skips_nonrankable_hrm_fit(monkeypatch):
    args = SimpleNamespace(
        trait="traits.tsv",
        states=None,
        regime_map=None,
        rate_matrix=None,
        transition_graph=None,
        rate=None,
    )
    context = ComparisonContext(
        tree=None,
        trait_df=pd.DataFrame(),
        trait_type="discrete",
        trait_columns=("x",),
        error_columns=None,
        args=args,
    )

    def unexpected_fit(*_args, **_kwargs):
        raise AssertionError("automatic HRM comparison must not start the fit")

    monkeypatch.setattr("nwkit.asr_compare._fit_candidate", unexpected_fit)
    rows = evaluate_comparison_candidates(
        context,
        [ComparisonCandidate("HRM", "equal", "model-default")],
        automatic=True,
    )
    assert rows[0]["status"] == "not_fitted"
    assert rows[0]["rankable"] == "no"
    assert "request HRM explicitly" in rows[0]["message"]


def test_fixed_eb_and_acdc_are_fitted_once_and_retained_as_alias(monkeypatch):
    args = SimpleNamespace(
        trait_type="continuous",
        evolution_parameter=-0.25,
        evolution_parameter_bounds=None,
        eb_rate=None,
        eb_rate_bounds=None,
        regime_map=None,
    )
    context = ComparisonContext(
        tree=None,
        trait_df=pd.DataFrame(),
        trait_type="continuous",
        trait_columns=("x",),
        error_columns=None,
        args=args,
        cache={"continuous_data": ({}, None)},
    )
    candidates = [
        ComparisonCandidate("EB", "flat", "model-default"),
        ComparisonCandidate("ACDC", "flat", "model-default"),
    ]
    fitted = []

    def fake_fit(_context, candidate):
        fitted.append(candidate.model)
        return SimpleNamespace(
            evolution_model=candidate.model.lower(),
            evolution_parameter_name="rate_change",
            evolution_parameter=-0.25,
            evolution_parameter_estimated=False,
            sigma2_estimated=True,
            restricted_log_likelihood=-3.0,
            num_effective_observations=8,
            fit_status="ok",
            optimizer_success=True,
            optimizer_message="",
        )

    monkeypatch.setattr("nwkit.asr_compare._fit_candidate", fake_fit)
    rows = evaluate_comparison_candidates(context, candidates, automatic=False)
    assert fitted == ["EB"]
    assert rows[1]["status"] == "equivalent"
    assert rows[1]["equivalent_to"] == "EB"
    assert rows[1]["rankable"] == "no"
    assert "fixed exponential-rate parameter" in rows[1]["message"]


def test_fixed_neutral_continuous_models_are_deduplicated(monkeypatch):
    args = SimpleNamespace(
        trait_type="continuous",
        evolution_parameter=1.0,
        evolution_parameter_bounds=None,
        eb_rate=None,
        eb_rate_bounds=None,
        drift=0.0,
        sigma2=None,
        regime_map=None,
        regime_parameters=None,
    )
    context = ComparisonContext(
        tree=None,
        trait_df=pd.DataFrame(),
        trait_type="continuous",
        trait_columns=("x",),
        error_columns=None,
        args=args,
        cache={"continuous_data": ({}, None)},
    )
    candidates = [
        ComparisonCandidate(model, "flat", "model-default")
        for model in ("BM", "LAMBDA", "KAPPA", "DELTA", "BM-DRIFT")
    ]
    fitted = []

    def fake_fit(_context, candidate):
        fitted.append(candidate.model)
        return SimpleNamespace(
            sigma2_estimated=True,
            restricted_log_likelihood=-3.0,
            num_effective_observations=8,
            fit_status="ok",
            optimizer_success=True,
            optimizer_message="",
        )

    monkeypatch.setattr("nwkit.asr_compare._fit_candidate", fake_fit)
    rows = evaluate_comparison_candidates(context, candidates, automatic=False)
    assert fitted == ["BM"]
    assert [row["status"] for row in rows] == ["ok"] + ["equivalent"] * 4
    assert {row["equivalent_to"] for row in rows[1:]} == {"BM"}


def test_one_regime_brownian_and_ou_models_are_deduplicated(monkeypatch):
    assignment = SimpleNamespace(regimes=("only",))
    args = SimpleNamespace(
        trait_type="continuous",
        evolution_parameter=None,
        evolution_parameter_bounds=None,
        eb_rate=None,
        eb_rate_bounds=None,
        drift=None,
        sigma2=None,
        alpha=None,
        theta=None,
        regime_map="regimes.tsv",
        regime_parameters=None,
    )
    context = ComparisonContext(
        tree=None,
        trait_df=pd.DataFrame(),
        trait_type="continuous",
        trait_columns=("x",),
        error_columns=None,
        args=args,
        cache={
            "continuous_data": ({}, None),
            "regime_assignment": assignment,
        },
    )
    candidates = [
        ComparisonCandidate("BM", "flat", "model-default"),
        ComparisonCandidate("BMS", "flat", "model-default"),
        ComparisonCandidate("BM-DRIFT", "flat", "model-default"),
        ComparisonCandidate("BMS-DRIFT", "flat", "model-default"),
        ComparisonCandidate("OU", "stationary", "model-default"),
        ComparisonCandidate("OUM", "stationary", "model-default"),
        ComparisonCandidate("OUMA", "stationary", "model-default"),
        ComparisonCandidate("OUMV", "stationary", "model-default"),
        ComparisonCandidate("OUMVA", "stationary", "model-default"),
    ]
    fitted = []

    def fake_fit(_context, candidate):
        fitted.append(candidate.model)
        common = {
            "num_effective_observations": 8,
            "fit_status": "ok",
            "optimizer_success": True,
            "optimizer_message": "",
            "sigma2_estimated": True,
        }
        if candidate.model == "BM":
            return SimpleNamespace(restricted_log_likelihood=-3.0, **common)
        if candidate.model == "BM-DRIFT":
            return SimpleNamespace(
                restricted_log_likelihood=-2.9,
                drift_estimated=True,
                **common,
            )
        return SimpleNamespace(
            model="OU",
            log_likelihood=-2.5,
            alpha_estimated=True,
            theta_estimated=True,
            root_prior="stationary",
            **common,
        )

    monkeypatch.setattr("nwkit.asr_compare._fit_candidate", fake_fit)
    rows = evaluate_comparison_candidates(context, candidates, automatic=False)
    assert fitted == ["BM", "BM-DRIFT", "OU"]
    aliases = {
        row["model"]: row["equivalent_to"]
        for row in rows
        if row["status"] == "equivalent"
    }
    assert aliases == {
        "BMS": "BM",
        "BMS-DRIFT": "BM-DRIFT",
        "OUM": "OU[stationary]",
        "OUMA": "OU[stationary]",
        "OUMV": "OU[stationary]",
        "OUMVA": "OU[stationary]",
    }


def test_shared_preparation_timing_includes_work_done_before_candidate_fits(
    monkeypatch,
):
    args = SimpleNamespace(
        trait_type="continuous",
        regime_map=None,
        evolution_parameter=None,
        evolution_parameter_bounds=None,
        eb_rate=None,
        eb_rate_bounds=None,
        sigma2=None,
    )
    context = ComparisonContext(
        tree=None,
        trait_df=pd.DataFrame(),
        trait_type="continuous",
        trait_columns=("x",),
        error_columns=None,
        args=args,
        cache={
            "_shared_preparation_started": 10.0,
            "continuous_data": ({}, None),
        },
    )
    monkeypatch.setattr("nwkit.asr_compare.time.perf_counter", lambda: 15.0)
    monkeypatch.setattr(
        "nwkit.asr_compare._fit_candidate",
        lambda *_args: SimpleNamespace(
            sigma2_estimated=True,
            restricted_log_likelihood=-3.0,
            num_effective_observations=8,
            fit_status="ok",
        ),
    )
    rows = evaluate_comparison_candidates(
        context,
        [ComparisonCandidate("BM", "flat", "model-default")],
        automatic=False,
    )
    assert rows[0]["shared_preparation_seconds"] == pytest.approx(5.0)


def test_regime_parameter_table_rejects_incompatible_model_contracts():
    args = SimpleNamespace(
        regime_map="regimes.tsv",
        regime_parameters="parameters.tsv",
        regime_model=None,
        root_prior=None,
        root_mean=None,
        root_variance=None,
        alpha=None,
        alpha_bounds=None,
        evolution_parameter=None,
        evolution_parameter_bounds=None,
        eb_rate=None,
        eb_rate_bounds=None,
    )
    context = ComparisonContext(
        tree=None,
        trait_df=pd.DataFrame(),
        trait_type="continuous",
        trait_columns=("x",),
        error_columns=None,
        args=args,
    )
    candidates = [
        ComparisonCandidate("BMS", "flat", "model-default"),
        ComparisonCandidate("OUM", "stationary", "model-default"),
    ]
    with pytest.raises(ValueError, match=r"BMS=sigma2; OUM=theta"):
        _validate_comparison_options(context, candidates)


def test_continuous_asr_cli_writes_model_comparison(tmp_path):
    trait = tmp_path / "traits.tsv"
    output = tmp_path / "asr.tsv"
    comparison = tmp_path / "comparison.tsv"
    trait.write_text("leaf_name\tx\nA\t0\nB\t1\nC\t2\nD\t3\nE\t5\n")
    main(
        [
            "asr",
            "-i",
            "[&R](((A:1,B:1):1,C:2):1,(D:1,E:1):2)R;",
            "--trait",
            str(trait),
            "--state-column",
            "x",
            "--compare-models",
            "BM,LAMBDA,KAPPA",
            "--model-comparison-out",
            str(comparison),
            "-o",
            str(output),
        ]
    )
    table = pd.read_csv(comparison, sep="\t")
    assert set(table["model"]) == {"BM", "LAMBDA", "KAPPA"}
    assert table["aic_weight"].sum() == pytest.approx(1.0)


def test_discrete_asr_cli_writes_model_comparison(tmp_path):
    trait = tmp_path / "traits.tsv"
    output = tmp_path / "asr.tsv"
    comparison = tmp_path / "comparison.tsv"
    trait.write_text("leaf_name\tx\nA\ta\nB\ta\nC\tb\nD\tb\nE\ta\n")
    main(
        [
            "asr",
            "-i",
            "[&R](((A:1,B:1):1,C:2):1,(D:1,E:1):2)R;",
            "--trait",
            str(trait),
            "--state-column",
            "x",
            "--compare-models",
            "ER,ARD",
            "--model-comparison-out",
            str(comparison),
            "-o",
            str(output),
        ]
    )
    table = pd.read_csv(comparison, sep="\t")
    assert set(table["model"]) == {"ER", "ARD"}
    assert set(table["sample_size"]) == {5}
    assert table["bic_weight"].sum() == pytest.approx(1.0)


def test_legacy_discrete_comparison_excludes_equivalent_alias_weight(tmp_path):
    trait = tmp_path / "traits.tsv"
    output = tmp_path / "asr.tsv"
    comparison = tmp_path / "comparison.tsv"
    trait.write_text("leaf_name\tx\nA\ta\nB\ta\nC\tb\nD\tb\nE\ta\n")
    main(
        [
            "asr",
            "-i",
            "[&R](((A:1,B:1):1,C:2):1,(D:1,E:1):2)R;",
            "--trait",
            str(trait),
            "--state-column",
            "x",
            "--compare-models",
            "ER,SYM,ARD",
            "--model-comparison-out",
            str(comparison),
            "-o",
            str(output),
        ]
    )
    table = pd.read_csv(comparison, sep="\t").set_index("model")
    assert table.loc["SYM", "status"] == "equivalent"
    assert table.loc["SYM", "equivalent_to"] == "ER"
    assert table.loc["SYM", "rankable"] == "no"
    assert pd.isna(table.loc["SYM", "aic_weight"])
    assert table.loc[["ER", "ARD"], "aic_weight"].sum() == pytest.approx(1.0)


def test_asrcompare_cli_groups_flat_and_stationary_root_models(tmp_path):
    trait = tmp_path / "traits.tsv"
    comparison = tmp_path / "comparison.tsv"
    figure = tmp_path / "comparison.pdf"
    trait.write_text("leaf_name\tx\nA\t0\nB\t1\nC\t2\nD\t3\nE\t5\n")
    main(
        [
            "asrcompare",
            "-i",
            "[&R](((A:1,B:1):1,C:2):1,(D:1,E:1):2)R;",
            "--trait",
            str(trait),
            "--state-column",
            "x",
            "--models",
            "BM,LAMBDA,OU",
            "--figure-out",
            str(figure),
            "-o",
            str(comparison),
        ]
    )
    table = pd.read_csv(comparison, sep="\t")
    assert set(table["model"]) == {"BM", "LAMBDA", "OU"}
    assert table["comparison_group"].nunique() == 2
    flat = table[table["likelihood_kind"] == "flat_root_integrated"]
    stationary = table[table["likelihood_kind"] == "stationary_root_ml"]
    assert flat["aic_weight"].sum() == pytest.approx(1.0)
    assert stationary["criterion_rank"].isna().all()
    assert set(table["root_prior"]) == {"flat", "stationary"}
    figure_bytes = figure.read_bytes()
    assert figure_bytes.startswith(b"%PDF")
    assert len(re.findall(rb"/Type /Page\b", figure_bytes)) == 1


def test_asrcompare_matches_legacy_compatible_model_scores(tmp_path):
    trait = tmp_path / "traits.tsv"
    legacy_asr = tmp_path / "asr.tsv"
    legacy_comparison = tmp_path / "legacy.tsv"
    batch_comparison = tmp_path / "batch.tsv"
    trait.write_text("leaf_name\tx\nA\t0\nB\t1\nC\t2\nD\t3\nE\t5\n")
    tree = "[&R](((A:1,B:1):1,C:2):1,(D:1,E:1):2)R;"
    main(
        [
            "asr",
            "-i",
            tree,
            "--trait",
            str(trait),
            "--state-column",
            "x",
            "--compare-models",
            "BM,LAMBDA,KAPPA",
            "--model-comparison-out",
            str(legacy_comparison),
            "-o",
            str(legacy_asr),
        ]
    )
    main(
        [
            "asrcompare",
            "-i",
            tree,
            "--trait",
            str(trait),
            "--state-column",
            "x",
            "--models",
            "BM,LAMBDA,KAPPA",
            "-o",
            str(batch_comparison),
        ]
    )
    legacy = pd.read_csv(legacy_comparison, sep="\t").sort_values("model")
    batch = pd.read_csv(batch_comparison, sep="\t").sort_values("model")
    assert list(legacy["model"]) == list(batch["model"])
    for column in (
        "log_likelihood",
        "num_parameters",
        "sample_size",
        "aic",
        "aicc",
        "bic",
        "aic_weight",
        "aicc_weight",
        "bic_weight",
    ):
        pd.testing.assert_series_equal(
            legacy[column].reset_index(drop=True),
            batch[column].reset_index(drop=True),
            check_names=False,
        )


def test_asrcompare_cli_compares_discrete_models(tmp_path):
    trait = tmp_path / "traits.tsv"
    comparison = tmp_path / "comparison.tsv"
    trait.write_text("leaf_name\tx\nA\ta\nB\ta\nC\tb\nD\tb\nE\ta\n")
    main(
        [
            "asrcompare",
            "-i",
            "[&R](A:1,B:1,C:1,D:1,E:1)R;",
            "--trait",
            str(trait),
            "--state-column",
            "x",
            "--models",
            "ER,ARD",
            "-o",
            str(comparison),
        ]
    )
    table = pd.read_csv(comparison, sep="\t")
    assert list(table["model"]) == ["ER", "ARD"]
    assert table["aic_weight"].sum() == pytest.approx(1.0)
    assert set(table["status"]) <= {"ok", "boundary"}


def test_binary_equivalent_models_are_retained_without_duplicate_ic_weight(tmp_path):
    trait = tmp_path / "traits.tsv"
    comparison = tmp_path / "comparison.tsv"
    trait.write_text("leaf_name\tx\nA\ta\nB\ta\nC\tb\nD\tb\nE\ta\n")
    main(
        [
            "asrcompare",
            "-i",
            "[&R](A:1,B:1,C:1,D:1,E:1)R;",
            "--trait",
            str(trait),
            "--state-column",
            "x",
            "--models",
            "ER,SYM,F81,GTR",
            "-o",
            str(comparison),
        ]
    )
    table = pd.read_csv(comparison, sep="\t")
    by_model = table.set_index("model")
    assert by_model.loc["SYM", "status"] == "equivalent"
    assert by_model.loc["SYM", "equivalent_to"] == "ER"
    assert by_model.loc["GTR", "status"] != "equivalent"
    assert pd.isna(by_model.loc["GTR", "equivalent_to"])
    assert table["rankable"].eq("yes").sum() == 3
    assert table["comparison_group"].nunique() == 2
    stationary = table[table["root_prior"] == "stationary"]
    assert stationary["criterion_rank"].notna().all()
    assert stationary["aic_weight"].sum() == pytest.approx(1.0)


def test_binary_gtr_is_fitted_independently_under_tight_rate_bounds(tmp_path):
    trait = tmp_path / "traits.tsv"
    comparison = tmp_path / "comparison.tsv"
    standalone = tmp_path / "gtr-model.tsv"
    trait.write_text("leaf_name\tx\nA\ta\nB\ta\nC\ta\nD\ta\nE\ta\nF\tb\nG\tb\nH\tb\n")
    tree = (
        "[&R]((((A:0.2,B:0.3):0.4,(C:0.2,D:0.5):0.3):0.6,E:0.8):0.4,"
        "(F:0.2,(G:0.3,H:0.7):0.2):0.9)R;"
    )
    common = [
        "-i",
        tree,
        "--trait",
        str(trait),
        "--state-column",
        "x",
        "--trait-type",
        "discrete",
        "--root-prior",
        "stationary",
        "--rate-bounds",
        "0.2,0.21",
    ]
    main(
        [
            "asrcompare",
            *common,
            "--models",
            "ARD,F81,GTR",
            "-o",
            str(comparison),
        ]
    )
    main(
        [
            "asr",
            *common,
            "--model",
            "GTR",
            "-o",
            str(tmp_path / "gtr-asr.tsv"),
            "--model-out",
            str(standalone),
        ]
    )
    table = pd.read_csv(comparison, sep="\t").set_index("model")
    assert table.loc["F81", "status"] == "equivalent"
    assert table.loc["F81", "equivalent_to"] == "ARD"
    assert table.loc["GTR", "status"] != "equivalent"
    standalone_log_likelihood = pd.read_csv(standalone, sep="\t").loc[
        0, "log_likelihood"
    ]
    assert table.loc["GTR", "log_likelihood"] == pytest.approx(
        standalone_log_likelihood
    )


def test_one_regime_mk_model_reuses_the_matching_base_fit(monkeypatch):
    states = ["a", "b", "c"]
    graph = pd.DataFrame(
        [[False, True, True], [True, False, True], [True, True, False]]
    ).to_numpy(dtype=bool)
    args = SimpleNamespace(
        trait_type="discrete",
        regime_map="regimes.tsv",
        regime_model="GTR",
        rate=0.4,
    )
    context = ComparisonContext(
        tree=None,
        trait_df=pd.DataFrame(),
        trait_type="discrete",
        trait_columns=("x",),
        error_columns=None,
        args=args,
        cache={
            "single_discrete_data": (states, {}, {}),
            "discrete_transition_graph": graph,
            "regime_assignment": SimpleNamespace(regimes=("only",)),
        },
    )
    candidates = [
        ComparisonCandidate("GTR", "stationary", "model-default"),
        ComparisonCandidate("MK-REGIME", "stationary", "--root-prior"),
    ]
    fitted = []

    def fake_fit(_context, candidate):
        fitted.append(candidate.model)
        return {
            "model": candidate.model,
            "rates": (0.2, 0.3, 0.4, 1.0, 1.1),
            "rate_estimated": True,
            "log_likelihood": -3.0,
            "sample_size": 8,
            "fit_status": "ok",
            "optimizer_success": True,
            "optimizer_message": "",
        }

    monkeypatch.setattr("nwkit.asr_compare._fit_candidate", fake_fit)
    rows = evaluate_comparison_candidates(context, candidates, automatic=False)
    assert fitted == ["GTR"]
    assert rows[1]["status"] == "equivalent"
    assert rows[1]["equivalent_to"] == "GTR"


def test_fixed_er_is_not_aliased_to_one_regime_er(monkeypatch):
    states = ["a", "b"]
    graph = pd.DataFrame([[False, True], [True, False]]).to_numpy(dtype=bool)
    args = SimpleNamespace(
        trait_type="discrete",
        regime_map="regimes.tsv",
        regime_model="ER",
        rate=0.4,
    )
    context = ComparisonContext(
        tree=None,
        trait_df=pd.DataFrame(),
        trait_type="discrete",
        trait_columns=("x",),
        error_columns=None,
        args=args,
        cache={
            "single_discrete_data": (states, {}, {}),
            "discrete_transition_graph": graph,
            "regime_assignment": SimpleNamespace(regimes=("only",)),
        },
    )
    candidates = [
        ComparisonCandidate("ER", "equal", "model-default"),
        ComparisonCandidate("MK-REGIME", "equal", "--root-prior"),
    ]
    fitted = []

    def fake_fit(_context, candidate):
        fitted.append(candidate.model)
        return {
            "model": candidate.model,
            "rates": (0.4,),
            "rate_estimated": candidate.model != "ER",
            "log_likelihood": -3.0,
            "sample_size": 8,
            "fit_status": "ok",
            "optimizer_success": True,
            "optimizer_message": "",
        }

    monkeypatch.setattr("nwkit.asr_compare._fit_candidate", fake_fit)
    rows = evaluate_comparison_candidates(context, candidates, automatic=False)
    assert fitted == ["ER", "MK-REGIME"]
    assert all(row["status"] == "ok" for row in rows)


def test_asrcompare_cli_keeps_ou_root_conventions_separate(tmp_path):
    trait = tmp_path / "traits.tsv"
    comparison = tmp_path / "comparison.tsv"
    trait.write_text("leaf_name\tx\nA\t-0.3\nB\t0.1\nC\t0.7\nD\t1.0\nE\t1.4\nF\t1.9\n")
    main(
        [
            "asrcompare",
            "-i",
            "[&R](((A:0.3,B:0.4):0.2,C:0.8):0.5,(D:0.4,(E:0.3,F:0.5):0.2):0.3)R;",
            "--trait",
            str(trait),
            "--state-column",
            "x",
            "--models",
            "OU[stationary],OU[fixed],OU[gaussian]",
            "--root-mean",
            "0",
            "--root-variance",
            "1",
            "-o",
            str(comparison),
        ]
    )
    table = pd.read_csv(comparison, sep="\t")
    assert set(table["model_id"]) == {
        "OU[stationary]",
        "OU[fixed]",
        "OU[gaussian]",
    }
    assert set(table["likelihood_kind"]) == {
        "stationary_root_ml",
        "proper_root_ml",
    }
    assert table["comparison_group"].nunique() == 3
    assert table["criterion_rank"].isna().all()


def test_asrcompare_custom_q_uses_matrix_state_order(tmp_path):
    trait = tmp_path / "traits.tsv"
    matrix = tmp_path / "q.tsv"
    comparison = tmp_path / "comparison.tsv"
    trait.write_text("leaf_name\tx\nA\ta\nB\ta\nC\tb\nD\tb\nE\ta\n")
    matrix.write_text("state\tb\ta\nb\t0\t0.3\na\t0.2\t0\n")
    main(
        [
            "asrcompare",
            "-i",
            "[&R](((A:1,B:1):1,C:2):1,(D:1,E:1):2)R;",
            "--trait",
            str(trait),
            "--state-column",
            "x",
            "--models",
            "ER,CUSTOM",
            "--rate-matrix",
            str(matrix),
            "-o",
            str(comparison),
        ]
    )
    table = pd.read_csv(comparison, sep="\t")
    custom = table[table["model"] == "CUSTOM"].iloc[0]
    assert custom["num_parameters"] == 0
    assert custom["fixed_parameters"] == "Q"
    assert table["aic_weight"].sum() == pytest.approx(1.0)


def test_asrcompare_regime_and_bm_models_share_flat_root_group(tmp_path):
    trait = tmp_path / "traits.tsv"
    regimes = tmp_path / "regimes.tsv"
    comparison = tmp_path / "comparison.tsv"
    trait.write_text("leaf_name\tx\nA\t0\nB\t1\nC\t2\nD\t4\nE\t5\nF\t7\n")
    regimes.write_text(
        "branch_id\tregime\n"
        "0\ta\n1\ta\n2\ta\n3\ta\n4\tb\n5\tb\n"
        "6\tb\n7\tb\n8\tb\n9\tb\n10\tb\n"
    )
    main(
        [
            "asrcompare",
            "-i",
            "[&R](((A:1,B:1):1,C:2):1,(D:1,(E:1,F:1):1):2)R;",
            "--trait",
            str(trait),
            "--state-column",
            "x",
            "--models",
            "BM,BMS",
            "--regime-map",
            str(regimes),
            "-o",
            str(comparison),
        ]
    )
    table = pd.read_csv(comparison, sep="\t")
    assert set(table["model"]) == {"BM", "BMS"}
    assert table["comparison_group"].nunique() == 1
    assert table["aic_weight"].sum() == pytest.approx(1.0)


def test_asrcompare_multivariate_root_conventions_are_not_cross_ranked(tmp_path):
    trait = tmp_path / "traits.tsv"
    comparison = tmp_path / "comparison.tsv"
    trait.write_text(
        "leaf_name\tx\ty\n"
        "A\t0.0\t0.0\nB\t0.9\t1.0\nC\t1.4\t2.3\nD\t2.3\t3.2\n"
        "E\t3.0\t4.5\nF\t3.8\t5.2\nG\t4.1\t6.7\nH\t5.3\t7.5\n"
    )
    main(
        [
            "asrcompare",
            "-i",
            "[&R]((((A:0.5,B:0.7):0.4,(C:0.6,D:0.8):0.3):0.5,"
            "((E:0.4,F:0.9):0.4,G:0.7):0.5):0.2,H:1.2)R;",
            "--trait",
            str(trait),
            "--state-column",
            "x,y",
            "--trait-type",
            "continuous",
            "--models",
            "MV-BM,MV-OU",
            "-o",
            str(comparison),
        ]
    )
    table = pd.read_csv(comparison, sep="\t")
    assert set(table["status"]) == {"ok"}
    assert set(table["sample_size"]) == {8}
    assert table["comparison_group"].nunique() == 2
    assert table["criterion_rank"].isna().all()


def test_asrcompare_reports_single_mk_mixture_without_fake_winner(tmp_path):
    trait = tmp_path / "traits.tsv"
    comparison = tmp_path / "comparison.tsv"
    trait.write_text(
        "leaf_name\tc1\tc2\tc3\n"
        "A\ta\ta\tb\nB\ta\tb\tb\nC\tb\tb\ta\n"
        "D\tb\ta\ta\nE\ta\tb\ta\n"
    )
    main(
        [
            "asrcompare",
            "-i",
            "[&R](((A:1,B:1):1,C:2):1,(D:1,E:1):2)R;",
            "--trait",
            str(trait),
            "--state-column",
            "c1,c2,c3",
            "--trait-type",
            "discrete",
            "--models",
            "MK-MIXTURE",
            "--rate-categories",
            "2",
            "-o",
            str(comparison),
        ]
    )
    row = pd.read_csv(comparison, sep="\t").iloc[0]
    assert row["model"] == "MK-MIXTURE"
    assert row["sample_size"] == 15
    assert pd.isna(row["criterion_rank"])
    assert pd.isna(row["aic_weight"])


def test_asrcompare_all_failures_write_tsv_but_preserve_pdf(tmp_path, monkeypatch):
    trait = tmp_path / "traits.tsv"
    comparison = tmp_path / "comparison.tsv"
    figure = tmp_path / "comparison.pdf"
    trait.write_text("leaf_name\tx\nA\t0\nB\t1\nC\t2\nD\t3\n")
    figure.write_bytes(b"existing-pdf")

    def fail_fit(*_args, **_kwargs):
        raise ValueError("synthetic fit failure")

    monkeypatch.setattr("nwkit.asr_compare._fit_candidate", fail_fit)
    with pytest.raises(ValueError, match="Every selected ASR model failed"):
        main(
            [
                "asrcompare",
                "-i",
                "[&R]((A:1,B:1):1,(C:1,D:1):1)R;",
                "--trait",
                str(trait),
                "--state-column",
                "x",
                "-o",
                str(comparison),
                "--figure-out",
                str(figure),
            ]
        )
    table = pd.read_csv(comparison, sep="\t")
    assert "failed" in set(table["status"])
    assert figure.read_bytes() == b"existing-pdf"


def test_asrcompare_pdf_failure_rolls_back_the_tsv(tmp_path, monkeypatch):
    trait = tmp_path / "traits.tsv"
    comparison = tmp_path / "comparison.tsv"
    figure = tmp_path / "comparison.pdf"
    trait.write_text("leaf_name\tx\nA\t0\nB\t1\nC\t2\nD\t3\nE\t5\n")
    comparison.write_text("existing-table\n")
    figure.write_bytes(b"existing-pdf")

    def fail_figure(*_args, **_kwargs):
        raise RuntimeError("synthetic figure failure")

    monkeypatch.setattr("nwkit.asr_compare.draw_comparison_figure", fail_figure)
    with pytest.raises(RuntimeError, match="synthetic figure failure"):
        main(
            [
                "asrcompare",
                "-i",
                "[&R](((A:1,B:1):1,C:2):1,(D:1,E:1):2)R;",
                "--trait",
                str(trait),
                "--state-column",
                "x",
                "--models",
                "BM,LAMBDA",
                "-o",
                str(comparison),
                "--figure-out",
                str(figure),
            ]
        )
    assert comparison.read_text() == "existing-table\n"
    assert figure.read_bytes() == b"existing-pdf"


def test_model_comparison_rejects_stdout_as_auxiliary_output(tmp_path):
    trait = tmp_path / "traits.tsv"
    trait.write_text("leaf_name\tx\nA\t0\nB\t1\nC\t2\n")
    with pytest.raises(ValueError, match="Auxiliary outputs require file paths"):
        main(
            [
                "asr",
                "-i",
                "[&R](A:1,B:1,C:1)R;",
                "--trait",
                str(trait),
                "--state-column",
                "x",
                "--compare-models",
                "BM,LAMBDA",
                "--model-comparison-out",
                "-",
            ]
        )


@pytest.mark.parametrize(
    "model,state_column,trait_type",
    [
        ("MK-MIXTURE", "x,y", "discrete"),
        ("THRESHOLD", "x", "discrete"),
        ("MV-BM", "x,y", "continuous"),
    ],
)
def test_models_without_compatible_comparison_contract_reject_option(
    tmp_path, model, state_column, trait_type
):
    trait = tmp_path / "traits.tsv"
    trait.write_text("leaf_name\tx\ty\nA\t0\t1\nB\t1\t0\nC\t0\t1\n")
    arguments = [
        "asr",
        "-i",
        "[&R](A:1,B:1,C:1)R;",
        "--trait",
        str(trait),
        "--state-column",
        state_column,
        "--trait-type",
        trait_type,
        "--model",
        model,
        "--compare-models",
        "ER,ARD" if trait_type == "discrete" else "BM,LAMBDA",
        "--model-comparison-out",
        str(tmp_path / "comparison.tsv"),
    ]
    if model == "THRESHOLD":
        arguments.extend(("--states", "0,1"))
    with pytest.raises(ValueError, match="--compare-models"):
        main(arguments)
