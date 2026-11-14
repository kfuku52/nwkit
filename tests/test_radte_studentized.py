import copy
import json
from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.optimize import linprog
from scipy.special import expit, logit
from scipy.stats import norm, t

from nwkit.cli import main
from nwkit.radte_model import DatingProblem, fit_dates, laplace_intervals
from nwkit.radte_sequence import QuadraticLikelihood
from nwkit.radte_studentized import (
    chronology_domain,
    rate_variance_information,
    studentized_intervals,
)
from tests.test_radte import cli_inputs, small_chronology


def quadratic_likelihood(c):
    lengths = np.array([n.dist for n in c.edges])
    mapping = np.array([0, 0, 1, 2, 3, 4])
    return QuadraticLikelihood(
        np.log(np.bincount(mapping, weights=lengths)),
        np.zeros(5),
        np.eye(5) * 100,
        0.0,
        mapping,
    )


def test_studentized_branch_interval_matches_local_gaussian_regression():
    c = small_chronology(max_age=100)
    fit, problem = fit_dates(c)
    original = copy.deepcopy(fit)
    studentized_intervals(fit, problem)
    assert rate_variance_information(fit, problem) == (6, 4)
    group = c.group_by_node[0]
    # Two root and four terminal log rates: Var(log root duration) is
    # sigma^2 * (1/2 + 1/4). Estimate sigma^2 with four residual df.
    sigma_squared = 6 * (np.log(2) / 2) ** 2 / 4
    age_se = (fit.ages[group] - 1) * np.sqrt(sigma_squared * (1 / 2 + 1 / 4))
    low = 1 + c.min_duration
    high = 10
    center = logit((fit.ages[group] - low) / (high - low))
    half_width = (
        t.ppf(0.975, 4)
        * age_se
        * (1 / (fit.ages[group] - low) + 1 / (high - fit.ages[group]))
    )
    expected = low + (high - low) * expit(
        np.array([center - half_width, center + half_width])
    )
    np.testing.assert_allclose(
        [fit.interval_lower[group], fit.interval_upper[group]], expected, rtol=1e-7
    )
    np.testing.assert_array_equal(fit.ages, original.ages)
    np.testing.assert_array_equal(fit.rates, original.rates)
    np.testing.assert_array_equal(fit.parameters, original.parameters)
    assert fit.objective == original.objective
    assert fit.interval_status == "conditional-studentized-curvature"
    assert "studentized_residual_df=4" in fit.diagnostics


@pytest.mark.parametrize(
    "inference,expected", [("marginal", (5, 3)), ("joint-map", (4, 3))]
)
def test_sequence_degrees_of_freedom_count_identifiable_rate_observations(
    inference, expected
):
    c = small_chronology(max_age=100)
    fit, problem = fit_dates(c, likelihood=quadratic_likelihood(c), inference=inference)
    assert problem.rate_variance_estimated
    assert rate_variance_information(fit, problem) == expected
    studentized_intervals(fit, problem)
    assert fit.interval_lower is not None
    fixed, fixed_problem = fit_dates(
        c, likelihood=quadratic_likelihood(c), inference=inference, rate_sd=0.3
    )
    assert not fixed_problem.rate_variance_estimated
    studentized_intervals(fixed, fixed_problem)
    assert fixed.interval_status == "conditional-bounded-normal-curvature"
    assert "studentized_variance_factor=1" in fixed.diagnostics


def test_fixed_rate_sd_uses_normal_critical_value_without_variance_inflation():
    c = small_chronology(max_age=100)
    fit, problem = fit_dates(c, rate_sd=0.3)
    studentized_intervals(fit, problem)
    group = c.group_by_node[0]
    low, high = 1 + c.min_duration, 10
    age = fit.ages[group]
    se = (age - 1) * 0.3 * np.sqrt(1 / 2 + 1 / 4)
    half_width = norm.ppf(0.975) * se * (1 / (age - low) + 1 / (high - age))
    center = logit((age - low) / (high - low))
    expected = low + (high - low) * expit(
        np.array([center - half_width, center + half_width])
    )
    np.testing.assert_allclose(
        [fit.interval_lower[group], fit.interval_upper[group]], expected, rtol=1e-7
    )


def test_feasible_domain_propagation_matches_independent_linear_programs(tmp_path):
    bounds = tmp_path / "bounds.tsv"
    bounds.write_text("node\tage_min\tage_max\nAB\t8\t12\n")
    c = small_chronology(bounds, max_age=30)
    lower = c.lower.copy()
    lower[c.group_by_node[0]] = 0
    c = replace(c, lower=lower)
    lo, hi = chronology_domain(c)
    problem = DatingProblem(c)
    for index in range(len(c.groups)):
        coordinate = np.eye(len(c.groups))[index]
        for sign, endpoint in [(1, lo[index]), (-1, hi[index])]:
            result = linprog(
                sign * coordinate,
                A_ub=-problem.constraints,
                b_ub=np.full(problem.constraints.shape[0], -c.min_duration),
                bounds=list(zip(c.lower, c.upper, strict=True)),
                method="highs",
            )
            assert result.success
            assert endpoint == pytest.approx(result.x[index], abs=1e-8)


def test_bound_transform_returns_feasible_endpoints_without_clipping():
    c = small_chronology(max_age=21)
    fit, problem = fit_dates(c, rate_sd=0.4)
    laplace_intervals(fit, problem)
    assert fit.interval_status == "unavailable-gaussian-interval-crosses-bound"
    studentized_intervals(fit, problem)
    group = c.group_by_node[0]
    assert c.lower[group] < fit.interval_lower[group] < fit.ages[group]
    assert fit.ages[group] < fit.interval_upper[group] < c.upper[group]


def test_studentized_does_not_publish_stale_intervals_after_bound_failure():
    c = small_chronology(max_age=100)
    fit, problem = fit_dates(c)
    studentized_intervals(fit, problem)
    assert fit.interval_lower is not None
    upper = c.upper.copy()
    upper[problem.free] = fit.ages[problem.free] + 1e-7
    constrained = DatingProblem(replace(c, upper=upper))
    studentized_intervals(fit, constrained)
    assert fit.interval_status == "unavailable-active-bound-use-profile-or-bootstrap"
    assert fit.interval_lower is None and fit.interval_upper is None


def test_studentized_refuses_unidentified_age_scale(tmp_path):
    bounds = tmp_path / "bounds.tsv"
    bounds.write_text("node\tage_min\tage_max\nAB\t8\t12\n")
    fit, problem = fit_dates(small_chronology(bounds))
    studentized_intervals(fit, problem)
    assert fit.interval_lower is None and fit.interval_upper is None
    assert fit.interval_status.startswith("unavailable-")


def test_studentized_cli_preserves_point_ages_and_shared_species_intervals(tmp_path):
    inputs = cli_inputs(tmp_path)
    plain, adjusted = tmp_path / "plain", tmp_path / "adjusted"
    for prefix, method in [(plain, "none"), (adjusted, "studentized")]:
        main(
            [
                "radte",
                *inputs,
                "--reconcile",
                "lca",
                "--out-prefix",
                str(prefix),
                "--uncertainty",
                method,
            ]
        )
    first = pd.read_csv(str(plain) + ".nodes.tsv", sep="\t")
    second = pd.read_csv(str(adjusted) + ".nodes.tsv", sep="\t")
    np.testing.assert_array_equal(first.estimated_age, second.estimated_age)
    shared = second[second.event_type == "speciation"]
    assert shared.interval_lower.nunique() == shared.interval_upper.nunique() == 1
    manifest = json.loads(Path(str(adjusted) + ".manifest.json").read_text())
    assert manifest["uncertainty"] == "conditional-studentized-curvature"
    assert manifest["calibration_policy"] == "hard-all-events"
