"""OU ancestral reconstruction with fixed or Gaussian nonstationary roots."""

import math
from dataclasses import dataclass

import numpy as np

from nwkit.compiled_tree import CompiledTree
from nwkit.continuous_asr import _ldexp
from nwkit.gaussian_inference import (
    GaussianConditioningResult,
    GaussianTreeLikelihood,
    condition_gaussian_tree,
    gaussian_process_parameter_rank,
    gaussian_tree_likelihood,
)
from nwkit.gaussian_tree import (
    GaussianRootPrior,
    GaussianTreeProcess,
    ou_branch_transition,
)
from nwkit.optimization import MultistartResult, deterministic_multistart
from nwkit.ou_asr import default_alpha_bounds


@dataclass(frozen=True)
class GeneralOuFit:
    alpha: float
    alpha_estimated: bool
    theta: float
    theta_estimated: bool
    sigma2: float
    sigma2_estimated: bool
    root_prior: str
    root_mean: float
    root_variance: float
    log_likelihood: float | None
    num_observed: int
    num_effective_observations: int
    num_observed_positions: int
    likelihood_rank: int
    fit_status: str
    optimizer_success: bool
    optimizer_message: str
    optimizer_starts: int
    optimizer_converged_starts: int
    optimizer_failed_starts: int
    optimizer_grid_evaluations: int
    alpha_bounds: tuple[float, float]
    root_variance_bounds: None = None


def _finite(value, label, *, positive=False):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{label} must be numeric, not boolean.")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{label} must be numeric and finite.") from exc
    if not math.isfinite(result) or (positive and result <= 0.0):
        qualifier = "positive and finite" if positive else "finite"
        raise ValueError(f"{label} must be {qualifier}.")
    return result


def _ou_process(tree, alpha, sigma2, theta, root_prior, root_mean, root_variance):
    transitions = {
        node: ou_branch_transition(float(node.dist), alpha, sigma2, theta)
        for node in tree.traverse()
        if not node.is_root
    }
    root = GaussianRootPrior(
        root_prior,
        root_mean,
        0.0 if root_prior == "fixed" else root_variance,
    )
    return GaussianTreeProcess(tree, transitions, root, "ou", alpha)


def compute_general_ou_marginals(
    tree,
    values_by_leaf,
    *,
    root_prior,
    root_mean,
    root_variance=None,
    alpha=None,
    sigma2=None,
    theta=None,
    alpha_bounds=None,
    standard_errors=None,
    _compute_posterior=True,
):
    """Condition OU1 under a fixed or user-specified Gaussian root prior."""

    if root_prior not in {"fixed", "gaussian"}:
        raise ValueError("General OU roots must be fixed or gaussian.")
    root_mean = _finite(root_mean, "--root-mean")
    if root_prior == "gaussian":
        if root_variance is None:
            raise ValueError("--root-prior gaussian requires --root-variance.")
        root_variance = _finite(root_variance, "--root-variance", positive=True)
    elif root_variance not in {None, 0}:
        raise ValueError("--root-prior fixed cannot take --root-variance.")
    root_variance = 0.0 if root_prior == "fixed" else float(root_variance)

    bounds = default_alpha_bounds(tree) if alpha_bounds is None else alpha_bounds
    bounds = float(bounds[0]), float(bounds[1])
    if (
        not all(math.isfinite(value) and value > 0.0 for value in bounds)
        or bounds[0] >= bounds[1]
    ):
        raise ValueError("OU alpha bounds must be increasing, positive, and finite.")
    fixed_alpha = None if alpha is None else _finite(alpha, "--alpha", positive=True)
    fixed_sigma2 = (
        None if sigma2 is None else _finite(sigma2, "--sigma2", positive=True)
    )
    fixed_theta = None if theta is None else _finite(theta, "--theta")

    observed = np.asarray(
        [float(value) for value in values_by_leaf.values() if value is not None],
        dtype=float,
    )
    if not len(observed) or np.any(~np.isfinite(observed)):
        raise ValueError("OU ancestral reconstruction requires finite observations.")
    compiled = CompiledTree.from_tree(tree)
    maximum_depth = 0.0
    depths = {tree: 0.0}
    for node in compiled.nodes[1:]:
        depth = depths[node.up] + float(node.dist)
        depths[node] = depth
        maximum_depth = max(maximum_depth, depth)
    time_scale = maximum_depth if maximum_depth > 0.0 else 1.0
    trait_size = max(
        float(np.max(np.abs(observed))), abs(root_mean), math.sqrt(root_variance)
    )
    _, trait_exponent = math.frexp(trait_size)
    scaled_observed = np.ldexp(observed, -trait_exponent)
    initial_theta = float(np.mean(scaled_observed))
    scaled_variance = max(
        float(np.var(scaled_observed)),
        float(np.ptp(scaled_observed)) ** 2 / len(observed),
        np.finfo(float).eps,
    )

    names: list[str] = []
    initial: list[float] = []
    parameter_bounds: list[tuple[float | None, float | None]] = []
    if fixed_alpha is None:
        names.append("log_alpha")
        log_alpha_bounds = (math.log(bounds[0]), math.log(bounds[1]))
        initial.append(0.5 * sum(log_alpha_bounds))
        parameter_bounds.append(log_alpha_bounds)
    if fixed_sigma2 is None:
        names.append("log_sigma2")
        min_log_rate = math.log(math.ulp(0.0))
        max_log_rate = math.log(np.finfo(float).max)
        log_rate = (
            math.log(scaled_variance)
            + 2 * trait_exponent * math.log(2.0)
            - math.log(time_scale)
        )
        log_rate = min(max_log_rate, max(min_log_rate, log_rate))
        initial.append(log_rate)
        parameter_bounds.append(
            (
                max(min_log_rate, log_rate - math.log(1e10)),
                min(max_log_rate, log_rate + math.log(1e10)),
            )
        )
    if fixed_theta is None:
        names.append("theta")
        initial.append(initial_theta)
        parameter_bounds.append((None, None))
    if len(observed) < len(names):
        raise ValueError(
            "OU free parameters are not identifiable from the number of observations; "
            "fix alpha, sigma2, or theta."
        )

    def unpack(parameters):
        supplied = dict(zip(names, parameters, strict=True))
        current_alpha = (
            fixed_alpha if fixed_alpha is not None else math.exp(supplied["log_alpha"])
        )
        current_sigma2 = (
            fixed_sigma2
            if fixed_sigma2 is not None
            else math.exp(supplied["log_sigma2"])
        )
        current_theta = (
            fixed_theta
            if fixed_theta is not None
            else _ldexp(float(supplied["theta"]), trait_exponent, "OU theta")
        )
        return current_alpha, current_sigma2, current_theta

    def process_for(parameters):
        current_alpha, current_sigma2, current_theta = unpack(parameters)
        return _ou_process(
            tree,
            current_alpha,
            current_sigma2,
            current_theta,
            root_prior,
            root_mean,
            root_variance,
        )

    observed_nodes = [
        leaf for leaf in tree.leaves() if values_by_leaf.get(str(leaf.name)) is not None
    ]
    if gaussian_process_parameter_rank(
        process_for,
        initial,
        observed_nodes,
        parameter_bounds,
        trait_exponent=trait_exponent,
    ) < len(names):
        raise ValueError(
            "OU free parameters are not jointly identifiable from the observed "
            "mean/covariance design; fix alpha, sigma2, or theta."
        )

    if names:
        prepared_likelihood = GaussianTreeLikelihood.prepare(
            compiled, values_by_leaf, standard_errors
        )

        def objective(parameters):
            try:
                result = prepared_likelihood.evaluate(process_for(parameters))
            except (ValueError, ArithmeticError, OverflowError):
                return 1e100
            return (
                -result.log_likelihood
                if math.isfinite(result.log_likelihood)
                else 1e100
            )

        optimized = deterministic_multistart(
            objective, initial, parameter_bounds, maxiter=1000
        )
    else:
        fixed_parameters = np.asarray(initial, dtype=float)
        optimized = MultistartResult(
            fixed_parameters, 0.0, True, "all parameters fixed", 0, 0, 0
        )
    fitted_alpha, fitted_sigma2, fitted_theta = unpack(optimized.x)
    evaluate = (
        condition_gaussian_tree if _compute_posterior else gaussian_tree_likelihood
    )
    result = evaluate(
        process_for(optimized.x),
        values_by_leaf,
        standard_errors=standard_errors,
        compiled_tree=compiled,
    )
    statuses = []
    tolerance = 1e-5
    if fixed_alpha is None:
        if fitted_alpha <= bounds[0] * (1.0 + tolerance):
            statuses.append("alpha_lower_boundary")
        elif fitted_alpha >= bounds[1] * (1.0 - tolerance):
            statuses.append("alpha_upper_boundary")
    if result.likelihood_rank < result.num_observed:
        statuses.append("singular_support")
        log_likelihood = None
    else:
        log_likelihood = result.log_likelihood
    fit = GeneralOuFit(
        alpha=fitted_alpha,
        alpha_estimated=fixed_alpha is None,
        theta=fitted_theta,
        theta_estimated=fixed_theta is None,
        sigma2=fitted_sigma2,
        sigma2_estimated=fixed_sigma2 is None,
        root_prior=root_prior,
        root_mean=root_mean,
        root_variance=root_variance,
        log_likelihood=log_likelihood,
        num_observed=result.num_observed,
        num_effective_observations=result.num_observed,
        num_observed_positions=result.num_observed_positions,
        likelihood_rank=result.likelihood_rank,
        fit_status="+".join(statuses) if statuses else "ok",
        optimizer_success=optimized.success,
        optimizer_message=optimized.message,
        optimizer_starts=optimized.starts,
        optimizer_converged_starts=optimized.converged_starts,
        optimizer_failed_starts=optimized.failed_starts,
        optimizer_grid_evaluations=0,
        alpha_bounds=bounds,
    )
    posterior = (
        result.marginals if isinstance(result, GaussianConditioningResult) else {}
    )
    return posterior, fit
