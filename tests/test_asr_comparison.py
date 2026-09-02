from types import SimpleNamespace

import pandas as pd
import pytest

from nwkit.asr_comparison import model_comparison_table, summarize_fit
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
