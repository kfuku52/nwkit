import numpy as np
import pytest
from ete4 import Tree
from numpy.polynomial.hermite import hermgauss
from scipy.optimize import approx_fprime
from scipy.special import logsumexp

from nwkit.radte_inputs import build_chronology
from nwkit.radte_marginal import MarginalDatingProblem
from nwkit.radte_model import fit_dates
from nwkit.radte_sequence import (
    QuadraticLikelihood,
    SequenceLikelihood,
    build_quadratic,
)
from nwkit.reconcile import build_reconciliation_table
from tests.test_radte import small_chronology
from tests.test_radte_sequence import simulated_alignment, write_alignment


def test_root_rate_integration_matches_independent_two_dimensional_quadrature():
    gene = Tree("(A:0.1,B:0.2)Root;", parser=1)
    species = Tree("(A:10,B:10)Root;", parser=1)
    c = build_chronology(
        gene, species, build_reconciliation_table(gene, species, {"A": "A", "B": "B"})
    )
    likelihood = QuadraticLikelihood(
        np.log([0.3]), np.zeros(1), np.array([[80.0]]), 3.0, np.array([0, 0])
    )
    problem = MarginalDatingProblem(
        c, likelihood=likelihood, rate_sd=0.5, quadrature_points=64
    )
    x = np.array([np.log(0.15)])
    value, gradient = problem.value_gradient(x)
    grid, weights = hermgauss(160)
    r = np.exp(x[0] + np.sqrt(2) * 0.5 * grid)
    lengths = r[:, None] + r[None, :]
    log_kernel = -3.0 - 40 * (np.log(lengths) - np.log(0.3)) ** 2
    log_kernel += np.log(weights)[:, None] + np.log(weights)[None, :] - np.log(np.pi)
    reference = -logsumexp(log_kernel)
    assert value == pytest.approx(reference, abs=1e-8)
    posterior_weights = np.exp(log_kernel + reference)
    expected_rate = np.sum(posterior_weights * r[:, None])
    np.testing.assert_allclose(problem.posterior_rates(x), expected_rate, rtol=1e-7)
    np.testing.assert_allclose(
        gradient,
        approx_fprime(x, lambda z: problem.value_gradient(z)[0], 1e-7),
        atol=1e-5,
    )


@pytest.mark.parametrize("rho", [0.0, 0.5])
@pytest.mark.parametrize("fixed_sd", [None, 0.0, 0.4])
def test_marginal_gradient_including_variance_and_correlations(rho, fixed_sd):
    c = small_chronology()
    lengths = np.array([n.dist for n in c.edges])
    # Tree traversal puts the two root edges first.
    mapping = np.array([0, 0, 1, 2, 3, 4])
    center = np.log(np.bincount(mapping, weights=lengths))
    hessian = np.eye(5) * 40 + np.full((5, 5), 0.5)
    likelihood = QuadraticLikelihood(center, np.zeros(5), hessian, 20, mapping)
    problem = MarginalDatingProblem(c, rho=rho, likelihood=likelihood, rate_sd=fixed_sd)
    x = problem.initial_parameters()
    x[0] *= 1.1
    _, gradient = problem.value_gradient(x)
    numeric = approx_fprime(x, lambda z: problem.value_gradient(z)[0], 1e-7)
    np.testing.assert_allclose(gradient, numeric, atol=1e-5, rtol=1e-4)


def test_marginal_sequence_estimates_are_invariant_to_root_branch_split(tmp_path):
    c = small_chronology()
    alignment = write_alignment(tmp_path, simulated_alignment(c, sites=4000))
    exact = SequenceLikelihood(c, alignment, model="jc69", gamma_categories=1)
    quadratic = build_quadratic(exact, np.array([n.dist for n in c.edges]))
    first, problem = fit_dates(c, likelihood=quadratic, starts=3)
    assert problem.marginal
    root_edges = [n for n in c.edges if n.up is c.gene]
    root_edges[0].dist, root_edges[1].dist = 0.29, 0.01
    second, _ = fit_dates(c, likelihood=quadratic, starts=3)
    assert first.objective == pytest.approx(second.objective, abs=1e-7)
    np.testing.assert_allclose(first.ages, second.ages, atol=1e-5)
    assert first.log_rate_sd == pytest.approx(second.log_rate_sd, abs=1e-5)


@pytest.mark.parametrize("rho", [0.0, 0.5])
def test_profile_structure_reuse_preserves_values_and_derivatives(rho, monkeypatch):
    from dataclasses import replace

    c = small_chronology()
    mapping = np.array([0, 0, 1, 2, 3, 4])
    q = QuadraticLikelihood(
        np.log(np.bincount(mapping, weights=[n.dist for n in c.edges])),
        np.arange(5) * 0.01,
        np.eye(5) * 40 + 0.5,
        20,
        mapping,
    )
    parent = MarginalDatingProblem(c, rho=rho, likelihood=q, rate_sd=0.4)
    group = parent.free[0]
    lo, hi = c.lower.copy(), c.upper.copy()
    lo[group] = hi[group] = c.initial[group]
    constrained = replace(c, lower=lo, upper=hi)
    cold = MarginalDatingProblem(constrained, rho=rho, likelihood=q, rate_sd=0.4)

    def forbidden_inverse(*args, **kwargs):
        raise AssertionError("Fixed profile matrices were recomputed")

    with monkeypatch.context() as m:
        m.setattr(np.linalg, "inv", forbidden_inverse)
        shared = MarginalDatingProblem(
            constrained,
            rho=rho,
            likelihood=q,
            rate_sd=0.4,
            shared_structure=parent.shared_structure,
        )
    assert shared.shared_structure is parent.shared_structure
    assert not shared.measurement_covariance.flags.writeable
    x = shared.initial_parameters()
    for candidate in (x, x + 0.001):
        a, ga = cold.value_gradient(candidate)
        b, gb = shared.value_gradient(candidate)
        assert a == pytest.approx(b, abs=1e-12)
        np.testing.assert_allclose(ga, gb, atol=1e-12)
        np.testing.assert_allclose(
            cold.posterior_rates(candidate), shared.posterior_rates(candidate)
        )
    # The input likelihood remains mutable: sharing must detect such changes.
    q.hessian[0, 0] += 1
    changed = MarginalDatingProblem(
        c, rho=rho, likelihood=q, rate_sd=0.4, shared_structure=parent.shared_structure
    )
    assert changed.shared_structure is not parent.shared_structure
    changed_rho = MarginalDatingProblem(
        c, rho=0.2, likelihood=q, rate_sd=0.4, shared_structure=changed.shared_structure
    )
    assert changed_rho.shared_structure is not changed.shared_structure
