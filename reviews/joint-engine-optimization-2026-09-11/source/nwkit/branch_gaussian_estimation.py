"""Bounded likelihood estimation with fixed BM/OU assignments and fixed jumps."""

import math
from dataclasses import asdict, replace

import numpy as np
from scipy.linalg import null_space, solve_triangular
from scipy.optimize import minimize

from nwkit.branch_gaussian import OUBranch, build_branch_gaussian_process
from nwkit.branch_gaussian_fit_spec import fitted_assignment
from nwkit.gaussian_inference import gaussian_tree_likelihood
from nwkit.gaussian_tree import GaussianRootPrior
from nwkit.optimization import deterministic_multistart


class _BranchFitProblem:
    def __init__(self, tree, observed, errors, assignment, root, parameters):
        self.tree, self.observed, self.errors = tree, observed, errors
        self.parameters = parameters
        # Center biological parameters before constructing OU intercepts: forming
        # (1-exp(-alpha*t))*theta first loses small variations on large offsets.
        center = (
            next((value for value in observed.values() if value is not None), 0.0)
            if root.mode == "flat"
            else root.mean
        )
        self.observed = {
            name: None if value is None else value - center
            for name, value in observed.items()
        }
        self.root = replace(root, mean=0.0)
        self.assignment = replace(
            assignment,
            models_by_branch_id={
                key: replace(
                    model,
                    diffusion=replace(
                        model.diffusion, optimum=model.diffusion.optimum - center
                    ),
                )
                if isinstance(model.diffusion, OUBranch)
                else model
                for key, model in assignment.models_by_branch_id.items()
            },
        )
        self.fit_parameters = tuple(
            replace(
                p,
                initial=p.initial - center,
                lower=p.lower - center,
                upper=p.upper - center,
            )
            if p.parameter == "theta"
            else p
            for p in parameters
        )
        self.initial = np.array([p.coordinate(p.initial) for p in parameters])
        self.nodes = tuple(
            node for node in tree.leaves() if observed.get(str(node.name)) is not None
        )
        dimension = len(self.nodes) - int(root.mode == "flat")
        if not len(parameters) < dimension or len(self.nodes) > 512:
            raise ValueError(
                "Branch fitting requires more independent observations than free "
                "groups (plus the flat root), and at most 512 observed tips."
            )
        if root.mode == "flat" and any(p.parameter == "alpha" for p in parameters):
            raise ValueError(
                "Estimating alpha requires a proper root (fixed/gaussian/stationary); "
                "a flat-root integral changes its root-loading measure with alpha. "
                "Fix alpha or specify a proper root."
            )
        process = self.process(self.initial)
        initial = gaussian_tree_likelihood(
            process, self.observed, standard_errors=errors
        )
        if initial.likelihood_rank != dimension or not math.isfinite(
            initial.log_likelihood
        ):
            raise ValueError(
                "Branch fitting requires a finite, full-rank observation likelihood."
            )
        self.initial_likelihood = initial.log_likelihood
        self.compiled = initial.compiled_tree
        self.dimension = dimension
        self.projection = self._projection(process)
        self.error_variances = np.array(
            [
                0.0 if errors is None else float(errors.get(str(node.name), 0.0)) ** 2
                for node in self.nodes
            ]
        )

    def process(self, coordinates):
        assignment = fitted_assignment(
            self.assignment, self.fit_parameters, coordinates
        )
        return build_branch_gaussian_process(
            self.tree, assignment.models_by_branch_id, root=self.root
        )

    def objective(self, coordinates):
        try:
            result = gaussian_tree_likelihood(
                self.process(coordinates),
                self.observed,
                standard_errors=self.errors,
                compiled_tree=self.compiled,
            )
        except (ValueError, OverflowError, FloatingPointError):
            return math.inf
        if result.likelihood_rank != self.dimension or not math.isfinite(
            result.log_likelihood
        ):
            return math.inf
        # Subtracting a very poor initial likelihood can erase all variation near
        # the optimum through cancellation (e.g. an initial sigma2 of 1e-200).
        return -result.log_likelihood / self.dimension

    def _projection(self, process):
        if self.root.mode != "flat":
            return np.eye(len(self.nodes))
        loadings = {self.tree: 1.0}
        for node in self.tree.traverse("preorder"):
            if not node.is_root:
                loadings[node] = process.transitions[node].slope * loadings[node.up]
        vector = np.array([loadings[node] for node in self.nodes])
        maximum = float(np.max(np.abs(vector)))
        if maximum == 0:
            raise ValueError(
                "The flat root has no representable influence on the observations."
            )
        return null_space((vector / maximum)[None, :]).T

    def moments(self, coordinates):
        process = self.process(coordinates)
        if self.root.mode == "flat":
            process = replace(process, root=GaussianRootPrior("fixed", 0.0, 0.0))
        means, _ = process.marginal_moments()
        covariance = process.covariance(self.nodes)
        covariance.flat[:: len(self.nodes) + 1] += self.error_variances
        mean = self.projection @ np.array([means[node] for node in self.nodes])
        covariance = self.projection @ covariance @ self.projection.T
        if not np.isfinite(mean).all() or not np.isfinite(covariance).all():
            raise ValueError(
                "Branch fitting moments exceed floating-point range; rescale units."
            )
        return mean, 0.5 * covariance + 0.5 * covariance.T


def _difference_points(coordinates, index, step=1e-5):
    left, right = np.array(coordinates, copy=True), np.array(coordinates, copy=True)
    left[index] = max(0.0, left[index] - step)
    right[index] = min(1.0, right[index] + step)
    return left, right, right[index] - left[index]


def _parameter_rank(problem, coordinates):
    """Rank mean/covariance derivatives after whitening the observed distribution."""
    mean, covariance = problem.moments(coordinates)
    try:
        factor = np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as exc:
        raise ValueError(
            "Branch fitting requires positive-definite observed contrast covariance."
        ) from exc
    triangle = np.triu_indices(len(mean))
    columns = []
    for index in range(len(coordinates)):
        left, right, width = _difference_points(coordinates, index)
        left_mean, left_cov = problem.moments(left)
        right_mean, right_cov = problem.moments(right)
        derivative_mean = solve_triangular(
            factor, (right_mean - left_mean) / width, lower=True
        )
        derivative_cov = solve_triangular(
            factor, (right_cov - left_cov) / width, lower=True
        )
        derivative_cov = solve_triangular(factor, derivative_cov.T, lower=True).T
        columns.append(np.concatenate((derivative_mean, derivative_cov[triangle])))
    jacobian = np.column_stack(columns)
    if not np.isfinite(jacobian).all():
        raise ValueError(
            "Branch parameter identifiability exceeds floating-point range."
        )
    scales = np.max(np.abs(jacobian), axis=0)
    if np.any(scales == 0):
        return int(np.count_nonzero(scales))
    normalized = jacobian / scales
    norms = np.linalg.norm(normalized, axis=0)
    informative = scales > 1e-8 / norms
    if not np.all(informative):
        return int(np.count_nonzero(informative))
    singular = np.linalg.svd(normalized / norms, compute_uv=False)
    return int(np.count_nonzero(singular > singular[0] * len(coordinates) * 1e-6))


def _gradient(objective, coordinates):
    result = []
    for index in range(len(coordinates)):
        left, right, width = _difference_points(coordinates, index, step=1e-7)
        result.append((objective(right) - objective(left)) / width)
    return np.array(result)


def _projected_gradient(gradient, coordinates):
    if not np.isfinite(gradient).all():
        return math.inf
    result = np.array(gradient, copy=True)
    result[(coordinates <= 1e-7) & (gradient > 0)] = 0.0
    result[(coordinates >= 1 - 1e-7) & (gradient < 0)] = 0.0
    return float(np.max(np.abs(result)))


def _verified_minimize(objective, initial, *, method, bounds, options):
    result = minimize(objective, initial, method=method, bounds=bounds, options=options)
    if np.isfinite(result.fun) and np.isfinite(result.x).all():
        gradient = _gradient(objective, result.x)
        if _projected_gradient(gradient, result.x) > 1e-4:
            # Polish with central differences when forward-difference stopping is premature.
            result = minimize(
                objective,
                result.x,
                jac=lambda x: _gradient(objective, x),
                method="L-BFGS-B",
                bounds=bounds,
                options={**options, "gtol": 1e-7, "ftol": 1e-14},
            )
        checked = _projected_gradient(_gradient(objective, result.x), result.x)
        result.success = bool(
            result.success and math.isfinite(checked) and checked <= 1e-4
        )
        result.message = f"{result.message}; checked projected gradient={checked:.6g}"
    return result


def _fit_metadata(problem, optimized, rank):
    groups, boundaries = [], []
    for parameter, coordinate in zip(problem.parameters, optimized.x, strict=True):
        boundary = (
            "lower"
            if coordinate <= 1e-6
            else "upper"
            if coordinate >= 1 - 1e-6
            else None
        )
        if boundary:
            boundaries.append(parameter.group)
        groups.append(
            {
                **asdict(parameter),
                "estimate": parameter.value(coordinate),
                "boundary": boundary,
            }
        )
    return {
        "method": "flat_root_integrated" if problem.root.mode == "flat" else "ML",
        "assignment_search": False,
        "root_parameters_estimated": False,
        "jump_parameters_estimated": False,
        "num_parameters_estimated": len(groups),
        "parameter_rank": rank,
        "identifiability": "numerical_local_full_rank",
        "initial_log_likelihood": problem.initial_likelihood,
        "groups": groups,
        "fit_status": "boundary" if boundaries else "ok",
        "optimizer_success": optimized.success,
        "optimizer_message": optimized.message,
        "optimizer_starts": optimized.starts,
        "optimizer_converged_starts": optimized.converged_starts,
        "optimizer_failed_starts": optimized.failed_starts,
        "optimizer_projected_gradient": _projected_gradient(
            _gradient(problem.objective, optimized.x), optimized.x
        ),
        "optimizer_global_optimum_guaranteed": False,
        "boundary_groups": boundaries,
    }


def estimate_branch_parameters(tree, observed, errors, assignment, root, parameters):
    problem = _BranchFitProblem(tree, observed, errors, assignment, root, parameters)
    count = len(parameters)
    # Avoid mistaking a degenerate initial alpha=0 for structural non-identifiability.
    ranks = [_parameter_rank(problem, problem.initial)] + [
        _parameter_rank(problem, np.full(count, fraction))
        for fraction in (0.25, 0.5, 0.75)
    ]
    if max(ranks) < count:
        raise ValueError(
            "Branch-fit parameters are not identifiable from the observed tips; fix or share parameters."
        )
    optimized = deterministic_multistart(
        problem.objective,
        problem.initial,
        [(0.0, 1.0)] * count,
        fallback=False,
        minimizer=_verified_minimize,
        ftol=1e-13,
        fractions=(0.15, 0.5, 0.85),
    )
    rank = _parameter_rank(problem, optimized.x)
    if rank < count:
        raise ValueError(
            "Branch-fit parameters are not identifiable at the fitted boundary/optimum; fix or share parameters."
        )
    if optimized.fun > -problem.initial_likelihood / problem.dimension + 1e-8:
        raise ValueError(
            "Branch parameter optimization is worse than the supplied initial model."
        )
    assignment = fitted_assignment(assignment, parameters, optimized.x)
    return assignment, _fit_metadata(problem, optimized, rank)
