"""Stationary full-attraction multivariate OU likelihood fitting."""

import warnings

import numpy as np

from nwkit.full_ou import (
    decode_full_ou,
    full_ou_identifiability,
    full_ou_initial,
    observation_moment_design,
)
from nwkit.multivariate_gaussian_asr import (
    DenseMultivariateFit,
    _likelihood_scale_adjustment,
    _prepare_observations,
    _restore_sigma,
    _validate_posterior_size,
)
from nwkit.optimization import deterministic_multistart
from nwkit.vector_fit import normalized_vector_observations, restored_vector_posterior
from nwkit.vector_gaussian import condition_vector_tree
from nwkit.vector_ou_fit import _ou_geometry
from nwkit.vector_processes import vector_ou_process


def _observations(tree, values, traits, errors, covariances):
    if covariances is not None:
        from nwkit.continuous_observation import prepare_correlated_observations

        return prepare_correlated_observations(
            tree, values, traits, covariances, errors
        )
    data = _prepare_observations(tree, values, traits, errors, dense_limit=False)
    observed, errors = normalized_vector_observations(data)
    return data, observed, errors, None


def _fixed_parameters(data, attraction, diffusion):
    if (attraction is None) != (diffusion is None):
        raise ValueError(
            "Full OU requires both fixed attraction and diffusion, or neither."
        )
    if attraction is None:
        return None
    dimension = len(data.trait_names)
    # Validate physical matrices before broadcasting the trait normalization.
    process = vector_ou_process(
        data.compiled.tree, attraction, diffusion, np.zeros(dimension)
    )
    attraction = np.asarray(attraction) * data.scales[None, :] / data.scales[:, None]
    diffusion = np.asarray(diffusion) / np.outer(data.scales, data.scales)
    covariance = process.root_covariance / np.outer(data.scales, data.scales)
    return attraction, diffusion, covariance


def _result(data, optimized, result, matrices, theta, fixed, diagnostic, boundary):
    attraction, diffusion, covariance = matrices
    status, rank, ratio = diagnostic
    statuses = []
    if status not in {"local_full_rank", "fixed_covariance_parameters"}:
        statuses.append(status)
    if boundary:
        statuses.append("covariance_parameter_boundary")
    return DenseMultivariateFit(
        trait_names=data.trait_names,
        sigma=_restore_sigma(data, covariance),
        sigma_rank=len(data.trait_names),
        sigma_estimated=fixed is None,
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
        model="MV-OU-FULL",
        diffusion_sigma=_restore_sigma(data, diffusion),
        theta=data.centers + data.scales * theta,
        theta_estimated=True,
        attraction_matrix=attraction * data.scales[:, None] / data.scales[None, :],
        attraction_estimated=fixed is None,
        identifiability_status=status,
        identifiability_rank=rank,
        identifiability_ratio=ratio,
    )


def fit_full_mvou(
    tree,
    values_by_leaf,
    trait_names,
    *,
    attraction=None,
    diffusion=None,
    standard_errors=None,
    measurement_covariances=None,
    compute_posterior=True,
):
    """Fit a stationary general OU; fixed matrices must be supplied together.

    C,D positive definite and K skew parameterize all stable A with positive
    definite diffusion. Free covariance parameters are optimized together with
    the trait optima. Numerical local identifiability is reported, not assumed.
    """
    from dataclasses import replace

    data, observed, errors, measurement_covariances = _observations(
        tree, values_by_leaf, trait_names, standard_errors, measurement_covariances
    )
    _validate_posterior_size(data, compute_posterior)
    dimension = len(data.trait_names)
    time_scale, _ = _ou_geometry(data)
    fixed = _fixed_parameters(data, attraction, diffusion)
    initial, bounds = ([], []) if fixed is not None else full_ou_initial(dimension)
    covariance_count = len(initial)
    if len(data.values) < covariance_count + dimension:
        raise ValueError(
            "Full OU has too few effective observed coordinates to estimate all parameters."
        )
    initial += [0.0] * dimension
    bounds += [(None, None)] * dimension

    def evaluate(parameters):
        matrices = (
            fixed
            if fixed is not None
            else decode_full_ou(parameters, dimension, time_scale)
        )
        theta = np.asarray(parameters[-dimension:])
        process = vector_ou_process(tree, matrices[0], matrices[1], theta)
        result = condition_vector_tree(process, observed, error_covariances=errors)
        return result, matrices, theta

    def objective(parameters):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("error", RuntimeWarning)
                return -evaluate(parameters)[0].log_likelihood
        except (ValueError, ArithmeticError, np.linalg.LinAlgError, RuntimeWarning):
            return 1e100

    optimized = deterministic_multistart(objective, initial, bounds, maxiter=1800)
    result, matrices, theta = evaluate(optimized.x)
    diagnostic = ("fixed_covariance_parameters", 0, 1.0)
    if fixed is None:
        design, complete = observation_moment_design(data)
        diagnostic = full_ou_identifiability(
            optimized.x[:covariance_count], dimension, time_scale, design, complete
        )
    boundary = any(
        (low is not None and value <= low + 1e-5)
        or (high is not None and value >= high - 1e-5)
        for value, (low, high) in zip(optimized.x, bounds, strict=True)
    )
    fit = _result(data, optimized, result, matrices, theta, fixed, diagnostic, boundary)
    return (
        restored_vector_posterior(result, data) if compute_posterior else {}
    ), replace(fit, measurement_covariances=measurement_covariances)
