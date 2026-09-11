"""Memory-linear fitting of noisy/incomplete correlated Brownian traits."""

import math

import numpy as np

from nwkit.clade_index import LcaIndex
from nwkit.multivariate_asr import MultivariateGaussianMarginal
from nwkit.multivariate_gaussian_asr import (
    DenseMultivariateFit,
    _decode_cholesky,
    _initial_cholesky,
    _likelihood_scale_adjustment,
    _mvbm_analysis_data,
    _prepare_observations,
    _restore_sigma,
    _trait_covariance_rows,
    _validate_posterior_size,
)
from nwkit.optimization import deterministic_multistart
from nwkit.vector_gaussian import VectorProcess, VectorTransition, condition_vector_tree


def _depths(compiled):
    depths = np.zeros(len(compiled.nodes))
    for index, node in enumerate(compiled.nodes[1:], 1):
        depths[index] = depths[compiled.parents[index]] + float(node.dist)
    if not np.isfinite(depths).all():
        raise ValueError("Multivariate tree depths exceed floating-point range.")
    return depths


def _validate_sparse_bm_design(data, depths):
    """Check covariance identifiability without allocating tip-pair matrices."""
    analysis, fixed = _mvbm_analysis_data(data)
    lca = LcaIndex(data.compiled.tree)
    tolerance = np.finfo(float).eps * float(np.max(depths)) * 100

    def shared(left, right):
        return depths[
            lca.common_ancestor_indices(
                int(analysis.node_indices[left]), int(analysis.node_indices[right])
            )
        ]

    def signal(left, left_reference, right, right_reference):
        terms = [shared(left, right)]
        if left_reference is not None:
            terms.append(-shared(left_reference, right))
        if right_reference is not None:
            terms.append(-shared(left, right_reference))
        if left_reference is not None and right_reference is not None:
            terms.append(shared(left_reference, right_reference))
        return abs(math.fsum(terms)) > tolerance

    rows = [
        _trait_covariance_rows(analysis, k, fixed) for k in range(len(data.trait_names))
    ]
    for first, (left, left_reference) in enumerate(rows):
        for second in range(first, len(rows)):
            right, right_reference = rows[second]
            if not any(
                signal(i, left_reference, j, right_reference)
                for i in left
                for j in right
            ):
                raise ValueError(
                    "MV-BM covariance components are not identifiable from the observed "
                    f"trait/branch overlap: {data.trait_names[first]}, {data.trait_names[second]}."
                )


def normalized_vector_observations(data):
    dimension = len(data.trait_names)
    observed: dict[str, list[float | None]] = {}
    errors: dict[str, np.ndarray] = {}
    for index, trait, value, error in zip(
        data.node_indices, data.trait_indices, data.values, data.errors, strict=True
    ):
        name = str(data.compiled.nodes[index].name)
        observed.setdefault(name, [None] * dimension)[trait] = float(value)
        errors.setdefault(name, np.zeros((dimension, dimension)))[trait, trait] = (
            error**2
        )
    return observed, errors


def _brownian_process(data, sigma, time_scale):
    dimension = len(data.trait_names)
    return VectorProcess(
        data.compiled.tree,
        dimension,
        {
            node: VectorTransition(
                np.eye(dimension),
                np.zeros(dimension),
                sigma * (float(node.dist) / time_scale),
            )
            for node in data.compiled.nodes[1:]
        },
    )


def restored_vector_posterior(result, data):
    return {
        node: MultivariateGaussianMarginal(
            data.centers + data.scales * result.means[index],
            _restore_sigma(data, result.covariances[index]),
        )
        for index, node in enumerate(result.nodes)
    }


def fit_pruning_mvbm(
    tree,
    values_by_leaf,
    trait_names,
    *,
    standard_errors=None,
    measurement_covariances=None,
    compute_posterior=True,
    _geometry_cache=None,
):
    """Fit the same flat-root integrated MV-BM likelihood using vector pruning."""
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
    dimension = len(data.trait_names)
    _validate_posterior_size(data, compute_posterior)
    if np.any(data.count_by_trait < 2):
        raise ValueError("MV-BM covariance needs at least two observations per trait.")
    if len(data.values) - dimension < dimension * (dimension + 1) // 2:
        raise ValueError(
            "MV-BM has too few effective observed coordinates to estimate all trait means and covariance parameters."
        )
    depths = _depths(data.compiled)
    _validate_sparse_bm_design(data, depths)
    time_scale = float(np.max(depths[data.node_indices])) or 1.0
    initial, bounds = _initial_cholesky(dimension)

    def evaluate(parameters):
        sigma, _ = _decode_cholesky(parameters, dimension)
        return sigma, condition_vector_tree(
            _brownian_process(data, sigma, time_scale),
            observed,
            error_covariances=errors,
        )

    def objective(parameters):
        try:
            return -evaluate(parameters)[1].log_likelihood
        except (ValueError, ArithmeticError, OverflowError, np.linalg.LinAlgError):
            return 1e100

    optimized = deterministic_multistart(objective, initial, bounds, maxiter=1200)
    sigma_scaled, result = evaluate(optimized.x)
    eigenvalues = np.linalg.eigvalsh(sigma_scaled)
    tolerance = (
        np.finfo(float).eps * max(1.0, float(np.max(eigenvalues))) * max(100, dimension)
    )
    rank = int(np.sum(eigenvalues > tolerance))
    fit = DenseMultivariateFit(
        trait_names=data.trait_names,
        sigma=_restore_sigma(data, sigma_scaled / time_scale),
        sigma_rank=rank,
        sigma_estimated=True,
        restricted_log_likelihood=result.log_likelihood
        - _likelihood_scale_adjustment(data, reml=True),
        log_likelihood=None,
        num_observed=data.num_observed_tips,
        num_effective_observations=data.num_effective_positions,
        residual_df=len(data.values) - dimension,
        fit_status="ok" if rank == dimension else "singular_covariance",
        optimizer_success=optimized.success,
        optimizer_message=optimized.message,
        optimizer_starts=optimized.starts,
        optimizer_converged_starts=optimized.converged_starts,
        optimizer_failed_starts=optimized.failed_starts,
        model="MV-BM",
        measurement_covariances=measurement_covariances,
    )
    return restored_vector_posterior(result, data) if compute_posterior else {}, fit
