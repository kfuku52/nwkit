"""Validation of simulation truth, failure accounting and estimator routing."""

import sys
from pathlib import Path

import numpy as np
import pytest
from ete4 import Tree

TOOLS = Path(__file__).resolve().parents[1] / "tools"
sys.path.insert(0, str(TOOLS))
from regression_calibration_design import (  # noqa: E402
    Case,
    generate,
    seed_for,
    tree_covariance,
)
from regression_calibration_engine import (  # noqa: E402
    evaluate,
    null_tail,
    record_refits,
)
from validate_regression_calibration import clean, summarize_records  # noqa: E402


def test_shared_edge_covariance_matches_direct_tree_distances():
    for shape in ("balanced", "comb"):
        text, covariance = tree_covariance(9, shape)
        tree = Tree(text, parser=1)
        leaves = {leaf.name: leaf for leaf in tree.leaves()}
        expected = np.empty_like(covariance)
        for i in range(9):
            for j in range(9):
                a, b = leaves[f"S{i}"], leaves[f"S{j}"]
                expected[i, j] = (
                    tree.get_distance(tree, a)
                    + tree.get_distance(tree, b)
                    - tree.get_distance(a, b)
                ) / 2
        expected /= np.mean(np.diag(expected))
        np.testing.assert_allclose(covariance, expected)
        assert np.linalg.eigvalsh(covariance).min() > 0


def test_oracle_calibrates_independent_null_generator():
    case = Case("oracle-test", size=12, copies=2, event_variance=0.6)
    rejected, covered = [], []
    for replicate in range(400):
        data = generate(case, seed_for(9281, case.name, replicate, "data"))
        result = evaluate(case, data, "oracle", 2, 1, 10)
        rejected.append(result["p_value"] < 0.05)
        covered.append(
            result["confidence_interval_lower"]
            <= 0
            <= result["confidence_interval_upper"]
        )
    # A broad deterministic generator smoke check, not a production calibration claim.
    assert 0.01 < np.mean(rejected) < 0.10
    assert 0.90 < np.mean(covered) < 0.99


def test_missingness_and_repeated_events_are_retained_in_truth():
    case = Case(
        "missing",
        size=10,
        copies=5,
        concentrated=True,
        missing=0.3,
        missing_mechanism="clade",
    )
    data = generate(case, 91)
    assert data["generated_rows"] == 18
    assert data["generated_events"] == 10
    assert min(data["events"]) == 3
    assert data["removed_rows"] == 11
    assert len(data["y"]) == 7


def test_conditional_predictor_error_covariance_shared_by_event():
    case = Case("error", size=4, copies=2, beta=0.5, predictor_variance=0.4)
    data = generate(case, 71)
    assert data["true_covariance"][0, 1] == pytest.approx(0.5**2 * 0.4)
    assert data["true_covariance"][0, 2] == 0


def test_production_rsc_bootstrap_observer_counts_reestimation():
    case = Case("small", size=6)
    data = generate(case, 35)
    first = evaluate(case, data, "parametric-bootstrap", 4, 29, 30)
    second = evaluate(case, data, "parametric-bootstrap", 4, 29, 30)
    assert first["status"] == "completed"
    assert first["bootstrap_successes"] == 4
    assert first["bootstrap_attempts"] >= 4
    assert first["p_value"] == second["p_value"]
    assert (
        len(
            {
                tuple(item["variance_components"].values())
                for item in first["bootstrap_attempt_records"]
            }
        )
        > 1
    )


def test_null_refit_failures_give_bounds_not_success_only_pvalue():
    result = null_tail(2.0, [3.0, None, 0.5, None])
    assert result["p_value"] is None
    assert result["p_value_failure_lower"] == 2 / 5
    assert result["p_value_failure_upper"] == 4 / 5
    assert result["bootstrap_attempts"] == 4
    assert result["bootstrap_successes"] == 2


def test_null_bootstrap_refits_no_intercept_null_without_predictors():
    case = Case("null", size=5)
    result = evaluate(case, generate(case, 53), "null-bootstrap", 4, 29, 30)
    assert result["status"] == "completed", result
    assert result["reference_objective"] == "normalized-common-gaussian-ML"
    assert result["bootstrap_attempts"] == 4
    assert result["bootstrap_successes"] == 4
    assert 0 < result["p_value"] <= 1


def test_summary_keeps_failed_and_ineligible_samples_in_denominator():
    prototype = {"method": "wald", "target_beta": 0, "seconds": 1}
    records = []
    for replicate, row in enumerate(
        [
            {
                "status": "completed",
                "inference_status": "ok",
                "p_value": 0.01,
                "coefficient": 2,
                "confidence_interval_lower": 1,
                "confidence_interval_upper": 3,
            },
            {"status": "fit_failed", "error": "optimizer"},
            {"status": "ineligible", "reason": "constant"},
            {"status": "completed", "inference_status": "singular", "p_value": None},
        ]
    ):
        records.append(
            {
                "case": {"name": "test"},
                "replicate": replicate,
                "results": [{**prototype, **row}],
            }
        )
    summary = summarize_records(records)[0]
    assert summary["rejection_all_generated"]["estimate"] == 0.25
    assert summary["rejection_given_available"]["estimate"] == 1
    assert summary["p_available"]["denominator"] == 4
    assert summary["fit_failure_given_eligible"]["denominator"] == 3
    assert summary["interval_available"]["estimate"] == 0.25
    assert summary["coverage_given_interval"]["estimate"] == 0
    with pytest.raises(ValueError, match="Duplicate"):
        summarize_records(records + records)


def test_seed_and_clean_serialization():
    assert seed_for(1, "case", 4, "data") == seed_for(1, "case", 4, "data")
    assert seed_for(1, "case", 4, "data") != seed_for(1, "case", 4, "fit")
    assert clean({"a": np.array([1, np.nan]), "b": np.bool_(True)}) == {
        "a": [1, None],
        "b": True,
    }


def test_observer_restores_estimator_on_failure():
    from nwkit import regress

    original = regress._profile_covariance_fit
    diagnostics = {}
    with (
        pytest.raises(RuntimeError, match="test"),
        record_refits("rsc", True, diagnostics),
    ):
        raise RuntimeError("test")
    assert regress._profile_covariance_fit is original
    assert diagnostics["bootstrap_attempts"] == 0


@pytest.mark.parametrize("family", ["binomial", "poisson", "negative-binomial"])
def test_non_gaussian_fit_uses_no_penalty_and_keeps_diagnostics(family):
    case = Case(
        "glmm",
        engine="glmm",
        size=20,
        family=family,
        baseline=0.5 if family == "binomial" else 2.0,
    )
    result = evaluate(case, generate(case, 18), "wald", 2, 29, 30)
    assert result["status"] == "completed", result
    assert "boundary_warning" in result
    assert "coefficient_covariance_status" in result


def test_nb2_reference_uses_inverse_dispersion_as_size():
    from regression_calibration_laplace import conditional_log_likelihood
    from scipy.stats import nbinom

    y = np.array([0, 1, 2, 8])
    eta = np.array([-1.0, 0.0, 1.0, 2.0])
    alpha = 0.3
    expected = nbinom.logpmf(y, 1 / alpha, 1 / (1 + alpha * np.exp(eta))).sum()
    assert conditional_log_likelihood(
        y, eta, "negative-binomial", alpha
    ) == pytest.approx(expected)


def test_laplace_reconstruction_independent_of_nwkit_likelihood():
    from regression_calibration_laplace import reference

    case = Case("qmc", engine="glmm", size=8)
    result = reference(case, generate(case, 84), powers=(8, 10), scrambles=4, seed=1)
    assert abs(result["reconstruction_error"]) < 1e-5
    assert np.isfinite(result["laplace_minus_reference"])


def test_raw_tip_pipeline_preserves_technical_replicate_invariance():
    from dataclasses import replace

    case = Case(
        "raw", engine="rsc-tips", size=8, biological_replicates=5, sampling_variance=1
    )
    first = evaluate(case, generate(case, 52), "wald", 2, 29, 30)
    technical = replace(case, technical_replicates=2)
    second = evaluate(technical, generate(technical, 52), "wald", 2, 29, 30)
    assert first["status"] == "completed", first
    assert second["status"] == "completed", second
    for field in ("coefficient", "standard_error", "p_value", "n_species_events"):
        assert first[field] == pytest.approx(second[field])
    assert first["n_species_events"] == 7


def test_raw_reconciled_paralogs_do_not_inflate_event_count():
    case = Case("raw", engine="rsc-tips", size=8, copies=3, event_variance=0.5)
    result = evaluate(case, generate(case, 52), "wald", 2, 29, 30)
    assert result["status"] == "completed", result
    assert result["n_species_events"] == 7
    assert result["n_gene_contrasts"] == 21


@pytest.mark.parametrize(
    "name", ["binomial-n60-p0.05", "negative-binomial-n30-boundary"]
)
def test_unpenalized_glmm_checks_beyond_successful_inferior_local_optimum(name):
    from regression_calibration_design import cases
    from regression_calibration_engine import call_glmm

    case = next(case for case in cases() if case.name == name)
    data = generate(case, seed_for(20260910, name, 1, "data"))
    design = np.column_stack([np.ones(len(data["y"])), data["x"]])
    null = call_glmm(case, data, design=design[:, :1])
    full = call_glmm(case, data, design=design)
    assert full.log_likelihood >= null.log_likelihood - 1e-5
    assert full.optimizer_converged


def test_comparison_rejects_changed_input_and_reports_availability_changes():
    from compare_regression_calibration import compare

    before = [
        {
            "case": {"name": "x"},
            "replicate": 0,
            "input_sha256": "same",
            "results": [{"method": "wald", "status": "fit_failed"}],
        }
    ]
    after = [
        {
            "case": {"name": "x"},
            "replicate": 0,
            "input_sha256": "same",
            "results": [
                {
                    "method": "wald",
                    "status": "completed",
                    "inference_status": "ok",
                    "p_value": 0.1,
                }
            ],
        }
    ]
    comparison = compare(before, after)[0]
    assert comparison["became_available"] == 1
    assert comparison["both_available"] == 0
    after[0]["input_sha256"] = "different"
    with pytest.raises(ValueError, match="input mismatch"):
        compare(before, after)
