"""Memory-linear stationary OU fitting for incomplete multivariate traits."""

import math
from dataclasses import replace

import numpy as np

from nwkit.multivariate_gaussian_asr import (
    DenseMultivariateFit,
    _contracted_positions,
    _decode_cholesky,
    _initial_cholesky,
    _likelihood_scale_adjustment,
    _prepare_observations,
    _restore_sigma,
    _validate_posterior_size,
    parse_alpha_by_trait,
)
from nwkit.optimization import deterministic_multistart
from nwkit.vector_fit import (
    normalized_vector_observations,
    restored_vector_posterior,
)
from nwkit.vector_gaussian import VectorProcess, VectorTransition, condition_vector_tree


def _ou_geometry(data):
    """Observed diameter, without cancellation from a long unobserved stem."""
    observed = set(data.node_indices)
    longest = np.full(len(data.compiled.nodes), -np.inf)
    diameter = 0.0
    for index in data.compiled.postorder:
        if index in observed:
            longest[index] = max(longest[index], 0.0)
        if index == 0:
            continue
        parent = data.compiled.parents[index]
        distance = longest[index] + float(data.compiled.nodes[index].dist)
        diameter = max(diameter, float(distance + longest[parent]))
        longest[parent] = max(longest[parent], distance)
    if not math.isfinite(diameter):
        raise ValueError("Observed OU distances exceed floating-point range.")
    contracted = _contracted_positions(data.compiled)
    distinct = []
    for trait in range(len(data.trait_names)):
        positions = data.node_indices[data.trait_indices == trait]
        distinct.append(len(set(contracted[positions])) > 1)
    return diameter / 2 or 1.0, distinct


def _alpha_spec(data, alpha, alpha_by_trait, alpha_bounds, diagonal):
    scale, distinct = _ou_geometry(data)
    bounds = (
        (1e-6 / scale, 50 / scale)
        if alpha_bounds is None
        else tuple(map(float, alpha_bounds))
    )
    if (
        len(bounds) != 2
        or not all(math.isfinite(v) and v > 0 for v in bounds)
        or bounds[0] >= bounds[1]
    ):
        raise ValueError("OU alpha bounds must be increasing and positive.")
    fixed = None
    if alpha_by_trait is not None:
        if not diagonal or alpha is not None:
            raise ValueError("--alpha-by-trait requires diagonal OU without --alpha.")
        fixed = parse_alpha_by_trait(alpha_by_trait, data.trait_names)
    elif alpha is not None:
        if not math.isfinite(float(alpha)) or float(alpha) <= 0:
            raise ValueError("OU alpha must be finite and positive.")
        fixed = np.full(len(data.trait_names), float(alpha))
    if fixed is None and not (all(distinct) if diagonal else any(distinct)):
        raise ValueError(
            "OU alpha is not identifiable without observations at distinct phylogenetic positions; fix alpha."
        )
    return fixed, bounds


def _process(data, rates, covariance, theta):
    transitions = {}
    for node in data.compiled.nodes[1:]:
        attenuation = np.exp(-rates * float(node.dist))
        innovation = covariance * -np.expm1(
            -(rates[:, None] + rates[None, :]) * float(node.dist)
        )
        transitions[node] = VectorTransition(
            np.diag(attenuation), theta * (1 - attenuation), innovation
        )
    return VectorProcess(data.compiled.tree, len(rates), transitions, theta, covariance)


def _fit_result(
    data,
    result,
    optimized,
    rates,
    sigma,
    diffusion,
    theta,
    fixed,
    bounds,
    diagonal,
    alpha,
):
    dimension = len(data.trait_names)
    eigenvalues = np.linalg.eigvalsh(sigma)
    rank = int(
        np.sum(
            eigenvalues
            > np.finfo(float).eps
            * max(1.0, float(np.max(eigenvalues)))
            * max(100, dimension)
        )
    )
    statuses = []
    if rank < dimension:
        statuses.append("singular_covariance")
    if fixed is None:
        if np.any(rates <= bounds[0] * (1 + 1e-5)):
            statuses.append("alpha_lower_boundary")
        if np.any(rates >= bounds[1] * (1 - 1e-5)):
            statuses.append("alpha_upper_boundary")
    return DenseMultivariateFit(
        trait_names=data.trait_names,
        sigma=_restore_sigma(data, sigma),
        sigma_rank=rank,
        sigma_estimated=True,
        restricted_log_likelihood=None,
        log_likelihood=result.log_likelihood
        - _likelihood_scale_adjustment(data, reml=False),
        num_observed=data.num_observed_tips,
        num_effective_observations=data.num_effective_positions,
        residual_df=len(data.values),
        fit_status="+".join(statuses) if statuses else "ok",
        optimizer_success=optimized.success,
        optimizer_message=optimized.message,
        optimizer_starts=optimized.starts,
        optimizer_converged_starts=optimized.converged_starts,
        optimizer_failed_starts=optimized.failed_starts,
        model="MV-OU-DIAG" if diagonal else "MV-OU",
        alpha=alpha if diagonal else float(rates[0]),
        alpha_by_trait=rates if diagonal else None,
        alpha_estimated=fixed is None,
        diffusion_sigma=_restore_sigma(data, diffusion) if diagonal else None,
        theta=data.centers + data.scales * theta,
        theta_estimated=True,
    )


def fit_pruning_mvou(
    tree,
    values_by_leaf,
    trait_names,
    *,
    alpha=None,
    alpha_by_trait=None,
    alpha_bounds=None,
    standard_errors=None,
    measurement_covariances=None,
    compute_posterior=True,
    diagonal=False,
    _geometry_cache=None,
):
    """Maximize the ordinary stationary likelihood, jointly fitting the optima.

    Storage is linear in tree size. Unlike flat-root BM, optima are fitted
    parameters here, not integrated random variables.
    """
    if measurement_covariances is None:
        data = _prepare_observations(
            tree, values_by_leaf, trait_names, standard_errors, dense_limit=False
        )
        observed, errors = normalized_vector_observations(data)
    else:
        from nwkit.continuous_observation import prepare_correlated_observations

        data, observed, errors, measurement_covariances = (
            prepare_correlated_observations(
                tree,
                values_by_leaf,
                trait_names,
                measurement_covariances,
                standard_errors,
            )
        )
    _validate_posterior_size(data, compute_posterior)
    dimension = len(data.trait_names)
    fixed, physical_bounds = _alpha_spec(
        data, alpha, alpha_by_trait, alpha_bounds, diagonal
    )
    rate_count = (dimension if diagonal else 1) if fixed is None else 0
    if len(data.values) < dimension + dimension * (dimension + 1) // 2 + rate_count:
        raise ValueError(
            "OU has too few effective observed coordinates to estimate all parameters."
        )
    initial, bounds = _initial_cholesky(dimension)
    initial_rate = math.sqrt(physical_bounds[0] * physical_bounds[1])
    initial = [math.log(initial_rate)] * rate_count + initial + [0.0] * dimension
    bounds = (
        [(math.log(physical_bounds[0]), math.log(physical_bounds[1]))] * rate_count
        + bounds
        + [(None, None)] * dimension
    )

    def evaluate(parameters):
        rates = fixed
        if rates is None:
            rates = np.exp(parameters[:rate_count])
            if not diagonal:
                rates = np.full(dimension, rates[0])
        covariance, _ = _decode_cholesky(parameters, dimension, rate_count)
        denominator = rates[:, None] + rates[None, :]
        sigma = covariance / denominator if diagonal else covariance
        diffusion = covariance if diagonal else covariance * denominator
        theta = np.asarray(parameters[-dimension:])
        result = condition_vector_tree(
            _process(data, rates, sigma, theta), observed, error_covariances=errors
        )
        return result, rates, sigma, diffusion, theta

    def objective(parameters):
        try:
            return -evaluate(parameters)[0].log_likelihood
        except (ValueError, ArithmeticError, np.linalg.LinAlgError):
            return 1e100

    optimized = deterministic_multistart(objective, initial, bounds, maxiter=1800)
    result, rates, sigma, diffusion, theta = evaluate(optimized.x)
    fit = _fit_result(
        data,
        result,
        optimized,
        rates,
        sigma,
        diffusion,
        theta,
        fixed,
        physical_bounds,
        diagonal,
        alpha,
    )
    return restored_vector_posterior(
        result, data
    ) if compute_posterior else {}, replace(
        fit, measurement_covariances=measurement_covariances
    )
