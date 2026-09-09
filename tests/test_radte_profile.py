"""Profile safety checks, continuation failures, and sequence likelihood endpoints."""

from dataclasses import replace

import numpy as np
import pytest
from scipy.stats import chi2

from nwkit import radte_uncertainty as uncertainty
from nwkit.radte_marginal import refine_quadrature
from nwkit.radte_model import DatingOptimizationError, fit_dates, solve_problem
from nwkit.radte_sequence import SequenceLikelihood, build_quadratic
from tests.test_radte import small_chronology
from tests.test_radte_sequence import simulated_alignment, write_alignment


def test_profile_rejects_inferior_reference_without_publishing_intervals():
    fit, problem = fit_dates(small_chronology(), rate_sd=0.3)
    fit.parameters[0] = 2.5
    fit.ages = problem.unpack_ages(fit.parameters)
    with pytest.raises(ValueError, match="better optimum.*Refit"):
        uncertainty.profile_intervals(fit, problem)
    assert fit.interval_lower is None and fit.interval_upper is None


def test_profile_reduces_failed_steps_and_keeps_diagnostics(monkeypatch):
    fit, problem = fit_dates(small_chronology(), rate_sd=0.3)
    expected, reference = fit_dates(small_chronology(), rate_sd=0.3)
    uncertainty.profile_intervals(expected, reference)
    calls = []

    def intermittent(profile, **kwargs):
        calls.append(profile.chronology.initial.copy())
        if len(calls) <= 3:
            raise DatingOptimizationError(
                "injected numerical failure",
                [dict(success=False, objective=99.0, iterations=1, message="injected")],
            )
        return solve_problem(profile, **kwargs)

    monkeypatch.setattr(uncertainty, "solve_problem", intermittent)
    uncertainty.profile_intervals(fit, problem)
    distances = [
        abs(ages[problem.free[0]] - fit.ages[problem.free[0]]) for ages in calls[:4]
    ]
    assert all(b < a for a, b in zip(distances[:-1], distances[1:], strict=True))
    np.testing.assert_allclose(fit.interval_lower, expected.interval_lower, atol=1e-6)
    np.testing.assert_allclose(fit.interval_upper, expected.interval_upper, atol=1e-6)
    recorded = [a for a in fit.attempts if a.get("phase") == "profile"]
    assert sum(not a["success"] for a in recorded) == 3
    assert all("profile_group" in a and "profile_age" in a for a in recorded)
    assert all("successful_starts" in a and "objective_spread" in a for a in recorded)


def test_permanent_solver_failure_has_bounded_retries(monkeypatch):
    fit, problem = fit_dates(small_chronology(), rate_sd=0.3)
    calls = []

    def fail(profile, **kwargs):
        calls.append(1)
        raise DatingOptimizationError("injected", [])

    monkeypatch.setattr(uncertainty, "solve_problem", fail)
    with pytest.raises(DatingOptimizationError, match="after 32 step reductions"):
        uncertainty.profile_intervals(fit, problem)
    assert len(calls) == 32
    assert fit.interval_lower is None


def test_approximation_errors_are_not_retried(monkeypatch):
    fit, problem = fit_dates(small_chronology(), rate_sd=0.3)
    calls = []

    def invalid(*args):
        calls.append(1)
        raise ValueError("invalid approximation")

    monkeypatch.setattr(uncertainty, "validate_profile_approximation", invalid)
    with pytest.raises(ValueError, match="invalid approximation"):
        uncertainty.profile_intervals(fit, problem)
    assert len(calls) == 1


def test_profile_local_optima_are_reported(monkeypatch):
    fit, problem = fit_dates(small_chronology(), rate_sd=0.3)

    def differing(profile, **kwargs):
        x, attempts = solve_problem(profile, **kwargs)
        attempts.append(
            dict(
                success=True,
                objective=attempts[0]["objective"] + 1,
                iterations=3,
                message="injected second local optimum",
            )
        )
        return x, attempts

    monkeypatch.setattr(uncertainty, "solve_problem", differing)
    uncertainty.profile_intervals(fit, problem)
    assert any(d.startswith("profile_multiple_local_optima:") for d in fit.diagnostics)


@pytest.mark.parametrize("method", ["exact", "marginal"])
def test_sequence_profile_endpoints_have_independently_refitted_likelihood_ratio(
    tmp_path, method
):
    c = small_chronology(max_age=30)
    alignment = write_alignment(tmp_path, simulated_alignment(c, sites=40000))
    likelihood = SequenceLikelihood(c, alignment, model="jc69", gamma_categories=1)
    if method == "marginal":
        likelihood = build_quadratic(likelihood, np.array([n.dist for n in c.edges]))
    fit, problem = fit_dates(c, likelihood=likelihood, rate_sd=0.3, starts=4)
    original = fit.parameters.copy()
    uncertainty.profile_intervals(fit, problem, starts=4)
    np.testing.assert_array_equal(fit.parameters, original)
    assert fit.interval_status == "conditional-profile"
    group = problem.free[0]
    baseline = problem.value_gradient(original)[0]
    for endpoint in (fit.interval_lower[group], fit.interval_upper[group]):
        lo, hi, initial = c.lower.copy(), c.upper.copy(), fit.ages.copy()
        lo[group] = hi[group] = initial[group] = endpoint
        constrained = replace(c, lower=lo, upper=hi, initial=initial)
        extra = (
            {"quadrature_points": problem.quadrature_points}
            if method == "marginal"
            else {}
        )
        profile = type(problem)(
            constrained,
            rho=problem.rho,
            likelihood=likelihood,
            rate_sd=problem.rate_sd,
            **extra,
        )
        # Fresh nuisance initialization and different random starts from the search.
        x, _ = solve_problem(profile, starts=5, seed=91)
        if method == "marginal":
            profile, x, _ = refine_quadrature(
                profile, x, starts=5, maxiter=2000, seed=91
            )
        ratio = 2 * (profile.value_gradient(x)[0] - baseline)
        assert ratio == pytest.approx(chi2.ppf(0.95, 1), abs=1e-4)


def test_profile_attempts_are_serialized_in_manifest(tmp_path):
    import json

    from nwkit.cli import main
    from tests.test_radte import cli_inputs

    prefix = tmp_path / "profile"
    main(
        [
            "radte",
            *cli_inputs(tmp_path),
            "--reconcile",
            "lca",
            "--rate-sd",
            "0.3",
            "--uncertainty",
            "profile",
            "--out-prefix",
            str(prefix),
        ]
    )
    manifest = json.loads(prefix.with_suffix(".manifest.json").read_text())
    attempts = [
        a for a in manifest["optimizer_attempts"] if a.get("phase") == "profile"
    ]
    assert attempts
    assert all(10 < a["profile_age"] < 30 for a in attempts)
    assert all(a["successful_starts"] >= 1 for a in attempts)
    assert all(a["objective_spread"] >= 0 for a in attempts)


def test_sequence_profile_rejects_inaccurate_quadratic_region(tmp_path):
    c = small_chronology()
    alignment = write_alignment(tmp_path, simulated_alignment(c, sites=4000))
    exact = SequenceLikelihood(c, alignment, model="jc69", gamma_categories=1)
    quadratic = build_quadratic(exact, np.array([n.dist for n in c.edges]))
    fit, problem = fit_dates(c, likelihood=quadratic, rate_sd=0.3, starts=4)
    with pytest.raises(ValueError, match="validated quadratic likelihood region"):
        uncertainty.profile_intervals(fit, problem, starts=4)
    assert fit.interval_lower is None


def test_profile_rejects_inferior_fit_even_at_small_rate_variance():
    chronology = small_chronology(
        gene_text="((A_1:0.1,B_1:0.1)S1:0.1,(A_2:0.1,B_2:0.1)S2:0.1)D;"
    )
    fit, problem = fit_dates(chronology, rate_sd=1e-6)
    fit.parameters[0] += 1e-5
    fit.ages = problem.unpack_ages(fit.parameters)
    with pytest.raises(ValueError, match="better optimum"):
        uncertainty.profile_intervals(fit, problem)
    assert fit.interval_lower is None
