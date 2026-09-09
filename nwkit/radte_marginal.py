"""Gaussian branch-likelihood integration with a one-dimensional root integral.

For a reversible model the two root edges have only one observable length.
Condition on their log-rate difference z. All identifiable log lengths are then
affine in the remaining Gaussian log rates, so those rates integrate exactly.
Only z needs Gauss-Hermite quadrature. This is an approximation to the sequence
likelihood, not an approximation that treats a sum of lognormals as lognormal.
"""

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy import sparse
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import Bounds, LinearConstraint
from scipy.special import logsumexp, roots_hermitenorm

from nwkit.radte_model import DatingProblem


@dataclass(frozen=True)
class _MarginalStructure:
    """Immutable model/edge quantities shared across constrained age fits."""

    key: tuple
    root_edges: np.ndarray
    root_row: int
    rows_to_edges: np.ndarray
    contrast_variance: float
    conditional_slope: np.ndarray
    conditional_covariance: np.ndarray
    projected_covariance: np.ndarray
    measurement_covariance: np.ndarray
    observation: np.ndarray
    kernel_minimum: float
    measurement_logdet: float


@lru_cache(maxsize=8)
def _normal_quadrature(points):
    nodes, weights = roots_hermitenorm(points)
    nodes = nodes[weights > 0]
    log_weights = np.log(weights[weights > 0]) - 0.5 * np.log(2 * np.pi)
    nodes.setflags(write=False)
    log_weights.setflags(write=False)
    return nodes, log_weights


class MarginalDatingProblem(DatingProblem):
    marginal = True

    def __init__(
        self,
        chronology,
        *,
        rho=0.0,
        likelihood=None,
        rate_sd=None,
        quadrature_points=32,
        shared_structure=None,
    ):
        super().__init__(chronology, rho=rho, likelihood=likelihood, rate_sd=rate_sd)
        self.quadrature_points = quadrature_points
        self.fixed_sd = rate_sd
        self.normal_nodes, self.log_weights = _normal_quadrature(quadrature_points)
        # Include array contents: callers can mutate a QuadraticLikelihood in
        # place. Bounds/initial ages are deliberately absent from this key.
        key = (
            rho,
            tuple(chronology.edges),
            tuple(chronology.parent),
            tuple(chronology.child),
            id(likelihood),
            likelihood.center.tobytes(),
            likelihood.gradient.tobytes(),
            likelihood.hessian.tobytes(),
            likelihood.mapping.tobytes(),
            likelihood.nll,
        )
        if shared_structure is None or shared_structure.key != key:
            shared_structure = self._build_structure(key)
        self.shared_structure = shared_structure
        self.root_edges = shared_structure.root_edges
        self.root_row = shared_structure.root_row
        self.rows_to_edges = shared_structure.rows_to_edges
        self.contrast_variance = shared_structure.contrast_variance
        self.conditional_slope = shared_structure.conditional_slope
        self.conditional_covariance = shared_structure.conditional_covariance
        self.projected_covariance = shared_structure.projected_covariance
        self.measurement_covariance = shared_structure.measurement_covariance
        self.observation = shared_structure.observation
        self.kernel_minimum = shared_structure.kernel_minimum
        self.measurement_logdet = shared_structure.measurement_logdet
        self._cached_sd = None
        self._cached_covariance = None

    def _build_structure(self, key):
        likelihood = self.likelihood
        root_edges = np.array(
            [
                i
                for i, n in enumerate(self.chronology.edges)
                if n.up is self.chronology.gene
            ]
        )
        root_row = int(likelihood.mapping[root_edges[0]])
        rows_to_edges = np.array(
            [
                np.flatnonzero(likelihood.mapping == row)[-1]
                for row in range(len(likelihood.center))
            ]
        )
        # Choose the second root edge as the Gaussian base; the first enters
        # through log(dt_first * exp(z) + dt_second).
        rows_to_edges[root_row] = root_edges[1]
        contrast = np.zeros(len(self.chronology.edges))
        contrast[root_edges] = [1, -1]
        covariance = np.linalg.inv(self.precision.toarray())
        covariance_contrast = covariance @ contrast
        contrast_variance = float(contrast @ covariance_contrast)
        conditional_slope = covariance_contrast / contrast_variance
        conditional_covariance = (
            covariance
            - np.outer(covariance_contrast, covariance_contrast) / contrast_variance
        )
        projected_covariance = conditional_covariance[
            np.ix_(rows_to_edges, rows_to_edges)
        ]
        measurement_covariance = np.linalg.inv(likelihood.hessian)
        observation = likelihood.center - measurement_covariance @ likelihood.gradient
        kernel_minimum = (
            likelihood.nll
            - 0.5 * likelihood.gradient @ measurement_covariance @ likelihood.gradient
        )
        measurement_logdet = np.linalg.slogdet(measurement_covariance)[1]
        for array in (
            root_edges,
            rows_to_edges,
            conditional_slope,
            conditional_covariance,
            projected_covariance,
            measurement_covariance,
            observation,
        ):
            array.setflags(write=False)
        return _MarginalStructure(
            key,
            root_edges,
            root_row,
            rows_to_edges,
            contrast_variance,
            conditional_slope,
            conditional_covariance,
            projected_covariance,
            measurement_covariance,
            observation,
            kernel_minimum,
            measurement_logdet,
        )

    def _covariance(self, sd):
        if sd != self._cached_sd:
            covariance = self.measurement_covariance + sd**2 * self.projected_covariance
            factor = cho_factor(covariance, lower=True)
            inverse = cho_solve(factor, np.eye(len(covariance)))
            logdet = float(2 * np.log(np.diag(factor[0])).sum())
            self._cached_sd, self._cached_covariance = sd, (factor, inverse, logdet)
        return self._cached_covariance

    def components(self, x):
        ages = self.unpack_ages(x)
        duration = self.chronology.durations(ages)
        mu = x[len(self.free)]
        sd = np.exp(x[-1]) if self.fixed_sd is None else self.fixed_sd
        z = self.normal_nodes * sd * np.sqrt(self.contrast_variance)
        log_duration = np.log(duration)
        means = (
            mu
            + log_duration[self.rows_to_edges][None, :]
            + z[:, None] * self.conditional_slope[self.rows_to_edges][None, :]
        )
        first, second = self.root_edges
        root_terms = np.stack(
            [log_duration[first] + z, np.full_like(z, log_duration[second])], axis=1
        )
        log_total = logsumexp(root_terms, axis=1)
        means[:, self.root_row] += log_total - log_duration[second]
        root_weight = np.exp(root_terms[:, 0] - log_total)
        factor, inverse, logdet = self._covariance(sd)
        residual = means - self.observation[None, :]
        precision_residual = cho_solve(factor, residual.T).T
        log_kernels = self.log_weights - 0.5 * np.sum(
            residual * precision_residual, axis=1
        )
        log_normalizer = logsumexp(log_kernels)
        weights = np.exp(log_kernels - log_normalizer)
        value = (
            self.kernel_minimum
            + 0.5 * (logdet - self.measurement_logdet)
            - log_normalizer
        )
        return value, duration, sd, z, root_weight, precision_residual, inverse, weights

    def value_gradient(self, x):
        duration = self.chronology.durations(self.unpack_ages(x))
        if np.any(duration <= 0):
            violation = np.minimum(duration - self.chronology.min_duration, 0)
            gradient = np.zeros_like(x)
            gradient[: len(self.free)] = self.free_design.T @ (2e12 * violation)
            return float(1e12 * (1 + violation @ violation)), gradient
        value, duration, sd, z, root_weight, v, inverse, weights = self.components(x)
        root = self.root_row
        edge_gradient = np.zeros(len(duration))
        edge_gradient[self.rows_to_edges] = weights @ v / duration[self.rows_to_edges]
        first, second = self.root_edges
        edge_gradient[first] = weights @ (v[:, root] * root_weight) / duration[first]
        edge_gradient[second] = (
            weights @ (v[:, root] * (1 - root_weight)) / duration[second]
        )
        gradient = np.concatenate(
            [self.free_design.T @ edge_gradient, [weights @ v.sum(axis=1)]]
        )
        if self.fixed_sd is None:
            covariance_derivative = 2 * sd**2 * self.projected_covariance
            mean_derivative = (
                z[:, None] * self.conditional_slope[self.rows_to_edges][None, :]
            )
            mean_derivative[:, root] += z * root_weight
            covariance_score = 0.5 * np.sum(
                inverse * covariance_derivative.T
            ) - 0.5 * np.einsum("ki,ij,kj->k", v, covariance_derivative, v)
            sigma_score = covariance_score + np.sum(v * mean_derivative, axis=1)
            gradient = np.append(gradient, weights @ sigma_score)
        return float(value), np.asarray(gradient)

    def initial_parameters(self, initial_ages=None):
        ages = self.chronology.initial if initial_ages is None else initial_ages
        _, mu, _, _, sse = self.branch_statistics(ages)
        initial = np.concatenate([ages[self.free], [mu]])
        if self.fixed_sd is None:
            sd = max(0.1, np.sqrt(sse / len(self.observed)))
            initial = np.append(initial, np.log(sd))
        return initial

    def bounds_and_constraint(self):
        c = self.chronology
        lower = c.lower[self.free].tolist() + [-30.0]
        upper = c.upper[self.free].tolist() + [30.0]
        if self.fixed_sd is None:
            lower.append(-9.0)
            upper.append(2.0)
        extras = len(lower) - len(self.free)
        matrix = sparse.hstack(
            [
                self.free_constraints,
                sparse.csr_matrix((len(self.constraint_lower), extras)),
            ]
        ).toarray()
        return Bounds(lower, upper), LinearConstraint(
            matrix, self.constraint_lower, np.inf
        )

    def rate_posterior(self, x):
        _, _, sd, z, _, precision_residual, inverse, weights = self.components(x)
        cross = sd**2 * self.conditional_covariance[:, self.rows_to_edges]
        means = (
            z[:, None] * self.conditional_slope[None, :] - precision_residual @ cross.T
        )
        covariance = sd**2 * self.conditional_covariance - cross @ inverse @ cross.T
        covariance = (covariance + covariance.T) / 2
        return means, covariance, weights

    def posterior_rates(self, x):
        means, covariance, weights = self.rate_posterior(x)
        mu = x[len(self.free)]
        return weights @ np.exp(mu + means + 0.5 * np.diag(covariance)[None, :])

    def approximation_checks(self, x):
        """Check the approximation at high-posterior-weight log-rate states."""
        means, covariance, weights = self.rate_posterior(x)
        selected = np.argsort(weights)[-min(3, len(weights)) :]
        states = [means[i] for i in selected]
        eigen, vectors = np.linalg.eigh(covariance)
        if eigen[-1] > 0:
            displacement = np.sqrt(eigen[-1]) * vectors[:, -1]
            mode = means[np.argmax(weights)]
            states.extend([mode + displacement, mode - displacement])
        duration = self.chronology.durations(self.unpack_ages(x))
        mu = x[len(self.free)]
        return [
            self.likelihood.check(duration * np.exp(mu + state)) for state in states
        ]


def refine_quadrature(problem, x, *, starts, maxiter, seed):
    from nwkit.radte_model import DatingOptimizationError, solve_problem

    attempts: list[dict] = []
    for points in [64, 128, 256, 512]:
        refined = MarginalDatingProblem(
            problem.chronology,
            rho=problem.rho,
            likelihood=problem.likelihood,
            rate_sd=problem.fixed_sd,
            quadrature_points=points,
            shared_structure=problem.shared_structure,
        )
        value, gradient = problem.value_gradient(x)
        check_value, check_gradient = refined.value_gradient(x)
        if (
            abs(value - check_value) < 1e-5
            and np.max(abs(gradient - check_gradient)) < 1e-4
        ):
            return problem, x, attempts
        try:
            x, new_attempts = solve_problem(
                refined, starts=starts, maxiter=maxiter, seed=seed, initial_parameters=x
            )
        except DatingOptimizationError as exc:
            exc.attempts = attempts + exc.attempts
            raise
        attempts.extend(new_attempts)
        problem = refined
    raise ValueError(
        "Root-rate quadrature did not converge through 512 points; use joint-map/exact or MCMCTree for this family."
    )
