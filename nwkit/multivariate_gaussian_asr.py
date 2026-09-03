"""Dense multivariate Gaussian conditioning for incomplete/noisy trait data."""

import math
from dataclasses import dataclass, replace

import numpy as np
from scipy.linalg import cho_factor, cho_solve

from nwkit.clade_index import LcaIndex
from nwkit.compiled_tree import CompiledTree
from nwkit.multivariate_asr import MultivariateGaussianMarginal
from nwkit.optimization import deterministic_multistart

_LOG_2PI = math.log(2.0 * math.pi)
_MAX_DENSE_OBSERVATIONS = 1000


@dataclass(frozen=True)
class DenseMultivariateFit:
    trait_names: tuple[str, ...]
    sigma: np.ndarray
    sigma_rank: int
    sigma_estimated: bool
    restricted_log_likelihood: float | None
    log_likelihood: float | None
    num_observed: int
    num_effective_observations: int
    residual_df: int
    fit_status: str
    optimizer_success: bool
    optimizer_message: str
    optimizer_starts: int
    optimizer_converged_starts: int
    optimizer_failed_starts: int
    model: str
    alpha: float | None = None
    alpha_by_trait: np.ndarray | None = None
    alpha_estimated: bool = False
    diffusion_sigma: np.ndarray | None = None
    theta: np.ndarray | None = None
    theta_estimated: bool = False


@dataclass(frozen=True)
class _ObservationData:
    compiled: CompiledTree
    trait_names: tuple[str, ...]
    node_indices: np.ndarray
    trait_indices: np.ndarray
    values: np.ndarray
    errors: np.ndarray
    centers: np.ndarray
    scales: np.ndarray
    count_by_trait: np.ndarray
    num_observed_tips: int
    num_effective_positions: int


@dataclass(frozen=True)
class _Geometry:
    compiled: CompiledTree
    depths: np.ndarray
    observed_shared_depth: np.ndarray
    observed_distance: np.ndarray
    node_shared_depth: np.ndarray
    node_distance: np.ndarray


def _finite(value, label, *, nonnegative=False):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{label} must be numeric, not boolean.")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{label} must be numeric and finite.") from exc
    if not math.isfinite(result) or (nonnegative and result < 0.0):
        qualifier = "non-negative and finite" if nonnegative else "finite"
        raise ValueError(f"{label} must be {qualifier}.")
    return result


def _contracted_positions(compiled):
    positions = np.arange(len(compiled.nodes), dtype=int)
    for node_index in range(1, len(compiled.nodes)):
        length = _finite(
            compiled.nodes[node_index].dist,
            "Multivariate Gaussian branch lengths",
            nonnegative=True,
        )
        if length == 0.0:
            positions[node_index] = positions[compiled.parents[node_index]]
    return positions


def _prepare_observations(
    tree, values_by_leaf, trait_names, standard_errors=None
) -> _ObservationData:
    trait_names = tuple(str(value) for value in trait_names)
    dimension = len(trait_names)
    if dimension < 2:
        raise ValueError("A multivariate model requires at least two traits.")
    compiled = CompiledTree.from_tree(tree)
    contracted_positions = _contracted_positions(compiled)
    raw_records = []
    exact_by_position: dict[tuple[int, int], float] = {}
    observed_tips = set()
    observed_positions = set()
    for name, node_index in compiled.leaf_index_by_name.items():
        vector = values_by_leaf.get(name)
        if vector is None:
            continue
        if len(vector) != dimension:
            raise ValueError(f"Trait dimension mismatch for tip '{name}'.")
        error_vector = None if standard_errors is None else standard_errors.get(name)
        if standard_errors is not None and error_vector is None:
            raise ValueError(
                f"Measurement errors are required for observed tip '{name}'."
            )
        if error_vector is not None and len(error_vector) != dimension:
            raise ValueError(f"Measurement-error dimension mismatch for tip '{name}'.")
        for trait_index, raw in enumerate(vector):
            if raw is None:
                continue
            value = _finite(raw, f"Trait '{trait_names[trait_index]}' for '{name}'")
            error = _finite(
                0.0 if error_vector is None else error_vector[trait_index],
                f"Standard error for trait '{trait_names[trait_index]}' at '{name}'",
                nonnegative=True,
            )
            if error == 0.0:
                key = (int(contracted_positions[node_index]), trait_index)
                previous = exact_by_position.get(key)
                if previous is not None:
                    if previous != value:
                        raise ValueError(
                            "Conflicting exact multivariate observations are "
                            "connected by zero-length branches."
                        )
                    observed_tips.add(name)
                    continue
                exact_by_position[key] = value
            raw_records.append((node_index, trait_index, value, error))
            observed_tips.add(name)
            observed_positions.add(int(contracted_positions[node_index]))
    if not raw_records:
        raise ValueError("A multivariate model requires at least one observed value.")
    if len(raw_records) > _MAX_DENSE_OBSERVATIONS:
        raise ValueError(
            "Incomplete/noisy multivariate ASR would require a dense covariance "
            f"larger than {_MAX_DENSE_OBSERVATIONS} observed coordinates; this "
            "limit bounds cubic factorization and all-node solve time."
        )

    centers = np.zeros(dimension, dtype=float)
    scales = np.ones(dimension, dtype=float)
    count_by_trait = np.zeros(dimension, dtype=int)
    for trait_index in range(dimension):
        records = [record for record in raw_records if record[1] == trait_index]
        if not records:
            raise ValueError(
                f"Trait '{trait_names[trait_index]}' has no observed values."
            )
        center = min(records, key=lambda item: (item[3], abs(item[2])))[2]
        differences = np.asarray([item[2] - center for item in records], dtype=float)
        if np.any(~np.isfinite(differences)):
            raise ValueError(
                f"Trait '{trait_names[trait_index]}' range exceeds floating-point range."
            )
        size = max(
            float(np.max(np.abs(differences))),
            max(item[3] for item in records),
        )
        exponent = math.frexp(size)[1] - 1 if size > 0.0 else 0
        centers[trait_index] = center
        scales[trait_index] = math.ldexp(1.0, exponent)
        count_by_trait[trait_index] = len(records)

    node_indices = np.asarray([item[0] for item in raw_records], dtype=int)
    trait_indices = np.asarray([item[1] for item in raw_records], dtype=int)
    values = np.asarray(
        [(item[2] - centers[item[1]]) / scales[item[1]] for item in raw_records],
        dtype=float,
    )
    errors = np.asarray(
        [item[3] / scales[item[1]] for item in raw_records], dtype=float
    )
    return _ObservationData(
        compiled,
        trait_names,
        node_indices,
        trait_indices,
        values,
        errors,
        centers,
        scales,
        count_by_trait,
        len(observed_tips),
        len(observed_positions),
    )


def _geometry(data: _ObservationData) -> _Geometry:
    compiled = data.compiled
    lca = LcaIndex(compiled.tree)
    depths = np.zeros(len(compiled.nodes), dtype=float)
    for index in range(1, len(compiled.nodes)):
        length = _finite(
            compiled.nodes[index].dist,
            "Multivariate Gaussian branch lengths",
            nonnegative=True,
        )
        depths[index] = depths[compiled.parents[index]] + length

    observed_nodes = data.node_indices
    count = len(observed_nodes)
    shared = np.empty((count, count), dtype=float)
    distance = np.empty((count, count), dtype=float)
    for first in range(count):
        first_index = int(observed_nodes[first])
        for second in range(first + 1):
            second_index = int(observed_nodes[second])
            ancestor = lca.common_ancestor_indices(first_index, second_index)
            shared_value = depths[ancestor]
            distance_value = (
                depths[first_index] + depths[second_index] - 2.0 * shared_value
            )
            shared[first, second] = shared[second, first] = shared_value
            distance[first, second] = distance[second, first] = max(0.0, distance_value)

    node_shared = np.empty((len(compiled.nodes), count), dtype=float)
    node_distance = np.empty_like(node_shared)
    for node_index in range(len(compiled.nodes)):
        for observed_index, observed_node in enumerate(observed_nodes):
            observed_node = int(observed_node)
            ancestor = lca.common_ancestor_indices(node_index, observed_node)
            shared_value = depths[ancestor]
            node_shared[node_index, observed_index] = shared_value
            node_distance[node_index, observed_index] = max(
                0.0,
                depths[node_index] + depths[observed_node] - 2.0 * shared_value,
            )
    return _Geometry(compiled, depths, shared, distance, node_shared, node_distance)


def _design(data):
    design = np.zeros((len(data.values), len(data.trait_names)), dtype=float)
    design[np.arange(len(data.values)), data.trait_indices] = 1.0
    return design


def _trait_covariance_rows(data, trait_index, fixed_mean):
    indices = np.flatnonzero(data.trait_indices == trait_index)
    if not np.isnan(fixed_mean[trait_index]):
        return indices, None
    if len(indices) < 2:
        return np.empty(0, dtype=int), None
    return indices[1:], int(indices[0])


def _contrast_covariance_block(
    scalar_covariance, left, left_reference, right, right_reference
):
    block = np.asarray(scalar_covariance[np.ix_(left, right)], dtype=float).copy()
    if left_reference is not None:
        block -= scalar_covariance[left_reference, right][None, :]
    if right_reference is not None:
        block -= scalar_covariance[left, right_reference][:, None]
    if left_reference is not None and right_reference is not None:
        block += scalar_covariance[left_reference, right_reference]
    return block


def _validate_mvbm_covariance_design(data, scalar_covariance, fixed_mean):
    """Require every trait-covariance component to affect an error contrast."""

    scale = float(np.max(np.abs(scalar_covariance), initial=0.0))
    tolerance = np.finfo(float).eps * max(1.0, scale) * max(100, len(data.values))
    rows = [
        _trait_covariance_rows(data, trait_index, fixed_mean)
        for trait_index in range(len(data.trait_names))
    ]
    unidentified = []
    for first in range(len(data.trait_names)):
        left, left_reference = rows[first]
        for second in range(first, len(data.trait_names)):
            right, right_reference = rows[second]
            identifiable = False
            if len(left) and len(right):
                block = _contrast_covariance_block(
                    scalar_covariance,
                    left,
                    left_reference,
                    right,
                    right_reference,
                )
                identifiable = bool(np.any(np.abs(block) > tolerance))
            if not identifiable:
                if first == second:
                    unidentified.append(f"variance({data.trait_names[first]})")
                else:
                    unidentified.append(
                        "covariance({},{})".format(
                            data.trait_names[first], data.trait_names[second]
                        )
                    )
    if unidentified:
        raise ValueError(
            "MV-BM covariance components are not identifiable from the observed "
            "trait/branch overlap after profiling trait means: "
            + ", ".join(unidentified)
            + "."
        )


def _validate_mvou_alpha_design(data, observed_distance):
    """Detect alpha/Sigma confounding when trait-pair distances are constant."""

    scale = float(np.max(np.abs(observed_distance), initial=0.0))
    tolerance = np.finfo(float).eps * max(1.0, scale) * max(100, len(data.values))
    for first in range(len(data.trait_names)):
        left = np.flatnonzero(data.trait_indices == first)
        for second in range(first, len(data.trait_names)):
            right = np.flatnonzero(data.trait_indices == second)
            block = observed_distance[np.ix_(left, right)]
            if block.size and float(np.ptp(block)) > tolerance:
                return
    raise ValueError(
        "MV-OU alpha is confounded with the trait covariance because every "
        "observed pair for each trait combination has one constant phylogenetic "
        "distance; fix --alpha or add observations at distinct positions."
    )


def _validate_mvou_diag_alpha_design(data, geometry):
    """Require each trait-specific decay to span distinct tree positions."""

    positions = _contracted_positions(data.compiled)
    unidentified = []
    for trait_index, trait_name in enumerate(data.trait_names):
        observed = np.flatnonzero(data.trait_indices == trait_index)
        trait_positions = {
            int(positions[int(data.node_indices[index])]) for index in observed
        }
        if len(trait_positions) < 2:
            unidentified.append(trait_name)
            continue
        block = geometry.observed_distance[np.ix_(observed, observed)]
        scale = float(np.max(np.abs(block), initial=0.0))
        tolerance = np.finfo(float).eps * max(1.0, scale) * max(100, len(observed))
        if not np.any(block > tolerance):
            unidentified.append(trait_name)
    if unidentified:
        raise ValueError(
            "MV-OU-DIAG cannot estimate trait-specific alpha without at least "
            "two observations at distinct phylogenetic positions for each trait: "
            + ", ".join(unidentified)
            + ". Fix --alpha/--alpha-by-trait or add observations."
        )


def _observation_covariance(data, scalar_covariance, sigma):
    covariance = (
        scalar_covariance
        * sigma[data.trait_indices[:, None], data.trait_indices[None, :]]
    )
    covariance = np.asarray(covariance, dtype=float)
    covariance[np.diag_indices_from(covariance)] += data.errors * data.errors
    return (covariance + covariance.T) / 2.0


def _factor_covariance(covariance):
    try:
        factor = cho_factor(covariance, lower=True, check_finite=False)
    except np.linalg.LinAlgError as exc:
        raise ValueError(
            "The observed multivariate covariance is singular; positive measurement "
            "errors can separate exact duplicate positions."
        ) from exc
    diagonal = np.diag(factor[0])
    if np.any(diagonal <= 0.0) or np.any(~np.isfinite(diagonal)):
        raise ValueError(
            "The observed multivariate covariance is not positive definite."
        )
    log_determinant = 2.0 * float(np.log(diagonal).sum())
    return factor, log_determinant


def _profile_mean(data, covariance, *, reml, fixed_mean=None):
    full_design = _design(data)
    dimension = len(data.trait_names)
    mean = np.zeros(dimension, dtype=float)
    if fixed_mean is None:
        free_mean_indices = np.arange(dimension, dtype=int)
    else:
        fixed_mean = np.asarray(fixed_mean, dtype=float)
        if fixed_mean.shape != (dimension,):
            raise ValueError("Fixed multivariate means have the wrong dimension.")
        fixed = ~np.isnan(fixed_mean)
        mean[fixed] = fixed_mean[fixed]
        free_mean_indices = np.flatnonzero(~fixed)
    design = full_design[:, free_mean_indices]
    factor, log_determinant = _factor_covariance(covariance)
    adjusted_values = data.values - full_design @ mean
    mean_factor = None
    mean_log_determinant = 0.0
    if len(free_mean_indices):
        solved_design = cho_solve(factor, design, check_finite=False)
        information = design.T @ solved_design
        try:
            mean_factor = cho_factor(information, lower=True, check_finite=False)
        except np.linalg.LinAlgError as exc:
            raise ValueError(
                "Multivariate trait means are not separately identifiable."
            ) from exc
        solved_values = cho_solve(factor, adjusted_values, check_finite=False)
        fitted = cho_solve(mean_factor, design.T @ solved_values, check_finite=False)
        mean[free_mean_indices] = fitted
        mean_log_determinant = 2.0 * float(np.log(np.diag(mean_factor[0])).sum())
    residual = data.values - full_design @ mean
    solved_residual = cho_solve(factor, residual, check_finite=False)
    quadratic = float(residual @ solved_residual)
    if reml:
        rank = len(data.values) - len(free_mean_indices)
        likelihood = -0.5 * (
            rank * _LOG_2PI + log_determinant + mean_log_determinant + quadratic
        )
    else:
        rank = len(data.values)
        likelihood = -0.5 * (rank * _LOG_2PI + log_determinant + quadratic)
    return likelihood, mean, residual, factor, mean_factor, free_mean_indices


def _decode_cholesky(parameters, dimension, offset=0):
    lower = np.zeros((dimension, dimension), dtype=float)
    index = offset
    for row in range(dimension):
        for column in range(row + 1):
            value = float(parameters[index])
            lower[row, column] = math.exp(value) if row == column else value
            index += 1
    sigma = lower @ lower.T
    return sigma, index


def _initial_cholesky(dimension):
    values: list[float] = []
    bounds: list[tuple[float | None, float | None]] = []
    for row in range(dimension):
        for column in range(row + 1):
            if row == column:
                values.append(0.0)
                bounds.append((-20.0, 20.0))
            else:
                values.append(0.0)
                bounds.append((None, None))
    return values, bounds


def _restore_sigma(data, sigma):
    return sigma * data.scales[:, None] * data.scales[None, :]


def _likelihood_scale_adjustment(data, *, reml):
    counts = data.count_by_trait - (1 if reml else 0)
    return float(np.dot(counts, np.log(data.scales)))


def _mvbm_analysis_data(data):
    positions = _contracted_positions(data.compiled)
    root_exact = (positions[data.node_indices] == 0) & (data.errors == 0.0)
    fixed_mean = np.full(len(data.trait_names), np.nan, dtype=float)
    fixed_mean[data.trait_indices[root_exact]] = data.values[root_exact]
    if not np.any(root_exact):
        return data, fixed_mean
    retained = ~root_exact
    if not np.any(retained):
        raise ValueError(
            "MV-BM covariance is not identifiable when every observation is an "
            "exact zero-distance root constraint."
        )
    counts = np.bincount(data.trait_indices[retained], minlength=len(data.trait_names))
    return (
        replace(
            data,
            node_indices=data.node_indices[retained],
            trait_indices=data.trait_indices[retained],
            values=data.values[retained],
            errors=data.errors[retained],
            count_by_trait=counts,
        ),
        fixed_mean,
    )


def _posterior(
    data,
    geometry,
    scalar_observed,
    scalar_cross,
    scalar_variance,
    sigma,
    mean,
    residual,
    factor,
    mean_factor,
    free_mean_indices,
    *,
    flat_root,
):
    dimension = len(data.trait_names)
    design = _design(data)[:, free_mean_indices]
    solved_residual = cho_solve(factor, residual, check_finite=False)
    solved_design = (
        cho_solve(factor, design, check_finite=False)
        if len(free_mean_indices)
        else np.empty((len(data.values), 0), dtype=float)
    )
    mean_covariance = (
        cho_solve(mean_factor, np.eye(len(free_mean_indices)), check_finite=False)
        if flat_root and mean_factor is not None
        else None
    )
    posterior = {}
    for node_index, node in enumerate(geometry.compiled.nodes):
        cross = scalar_cross[node_index][None, :] * sigma[:, data.trait_indices]
        scaled_mean = mean + cross @ solved_residual
        solved_cross = cho_solve(factor, cross.T, check_finite=False)
        covariance = scalar_variance[node_index] * sigma - cross @ solved_cross
        if flat_root and len(free_mean_indices):
            adjustment = np.zeros((dimension, len(free_mean_indices)), dtype=float)
            adjustment[free_mean_indices, np.arange(len(free_mean_indices))] = 1.0
            adjustment -= cross @ solved_design
            assert mean_covariance is not None
            covariance += adjustment @ mean_covariance @ adjustment.T
        covariance = (covariance + covariance.T) / 2.0
        eigenvalues, eigenvectors = np.linalg.eigh(covariance)
        tolerance = (
            np.finfo(float).eps
            * max(1.0, float(np.max(np.abs(eigenvalues))))
            * max(100, dimension)
        )
        if float(np.min(eigenvalues)) < -tolerance:
            raise ValueError(
                "A multivariate posterior covariance is not positive semidefinite."
            )
        covariance = (eigenvectors * np.maximum(eigenvalues, 0.0)) @ eigenvectors.T
        restored_mean = data.centers + data.scales * scaled_mean
        restored_covariance = _restore_sigma(data, covariance)
        posterior[node] = MultivariateGaussianMarginal(
            restored_mean, restored_covariance
        )
    return posterior


def _diagonal_ou_observation_covariance(data, geometry, beta, stationary_sigma):
    observed_depths = geometry.depths[data.node_indices]
    left = np.maximum(0.0, observed_depths[:, None] - geometry.observed_shared_depth)
    right = np.maximum(0.0, observed_depths[None, :] - geometry.observed_shared_depth)
    observed_beta = beta[data.trait_indices]
    decay = np.exp(-observed_beta[:, None] * left - observed_beta[None, :] * right)
    covariance = (
        decay
        * stationary_sigma[data.trait_indices[:, None], data.trait_indices[None, :]]
    )
    covariance[np.diag_indices_from(covariance)] += data.errors * data.errors
    return (covariance + covariance.T) / 2.0


def _diagonal_ou_cross_covariance(data, geometry, node_index, beta, stationary_sigma):
    observed_depths = geometry.depths[data.node_indices]
    shared = geometry.node_shared_depth[node_index]
    left = np.maximum(0.0, geometry.depths[node_index] - shared)
    right = np.maximum(0.0, observed_depths - shared)
    observed_beta = beta[data.trait_indices]
    decay = np.exp(-beta[:, None] * left - observed_beta[None, :] * right)
    return decay * stationary_sigma[:, data.trait_indices]


def _diagonal_ou_posterior(
    data,
    geometry,
    beta,
    stationary_sigma,
    mean,
    residual,
    factor,
):
    solved_residual = cho_solve(factor, residual, check_finite=False)
    posterior = {}
    dimension = len(data.trait_names)
    for node_index, node in enumerate(geometry.compiled.nodes):
        cross = _diagonal_ou_cross_covariance(
            data, geometry, node_index, beta, stationary_sigma
        )
        scaled_mean = mean + cross @ solved_residual
        solved_cross = cho_solve(factor, cross.T, check_finite=False)
        covariance = stationary_sigma - cross @ solved_cross
        covariance = (covariance + covariance.T) / 2.0
        eigenvalues, eigenvectors = np.linalg.eigh(covariance)
        tolerance = (
            np.finfo(float).eps
            * max(1.0, float(np.max(np.abs(eigenvalues))))
            * max(100, dimension)
        )
        if float(np.min(eigenvalues)) < -tolerance:
            raise ValueError(
                "A multivariate posterior covariance is not positive semidefinite."
            )
        covariance = (eigenvectors * np.maximum(eigenvalues, 0.0)) @ eigenvectors.T
        posterior[node] = MultivariateGaussianMarginal(
            data.centers + data.scales * scaled_mean,
            _restore_sigma(data, covariance),
        )
    return posterior


def parse_alpha_by_trait(value, trait_names):
    """Parse one fixed positive attraction rate per named trait."""

    if value is None or (isinstance(value, str) and value == ""):
        return None
    raw = (
        value if isinstance(value, (tuple, list, np.ndarray)) else str(value).split(",")
    )
    if len(raw) != len(trait_names):
        raise ValueError(
            "--alpha-by-trait must contain exactly one value per --state-column "
            f"trait ({len(trait_names)} expected)."
        )
    result = np.asarray(
        [
            _finite(item, f"--alpha-by-trait for '{trait}'")
            for item, trait in zip(raw, trait_names, strict=True)
        ],
        dtype=float,
    )
    if np.any(result <= 0.0):
        raise ValueError("--alpha-by-trait values must be strictly positive.")
    return result


def fit_dense_mvbm(
    tree,
    values_by_leaf,
    trait_names,
    *,
    standard_errors=None,
):
    """Fit flat-root MV-BM with arbitrary per-trait missingness and errors."""

    data = _prepare_observations(
        tree, values_by_leaf, trait_names, standard_errors=standard_errors
    )
    dimension = len(data.trait_names)
    if np.any(data.count_by_trait < 2):
        missing = [
            data.trait_names[index]
            for index, count in enumerate(data.count_by_trait)
            if count < 2
        ]
        raise ValueError(
            "MV-BM covariance needs at least two observations per trait: "
            + ", ".join(missing)
        )
    covariance_parameters = dimension * (dimension + 1) // 2
    if len(data.values) - dimension < covariance_parameters:
        raise ValueError(
            "MV-BM has too few effective observed coordinates to estimate all "
            "trait means and covariance parameters."
        )
    analysis_data, fixed_mean = _mvbm_analysis_data(data)
    geometry = _geometry(analysis_data)
    _validate_mvbm_covariance_design(
        analysis_data, geometry.observed_shared_depth, fixed_mean
    )
    initial, bounds = _initial_cholesky(dimension)

    def evaluate(parameters):
        sigma, _ = _decode_cholesky(parameters, dimension)
        covariance = _observation_covariance(
            analysis_data, geometry.observed_shared_depth, sigma
        )
        return (
            sigma,
            covariance,
            _profile_mean(analysis_data, covariance, reml=True, fixed_mean=fixed_mean),
        )

    def objective(parameters):
        try:
            return -evaluate(parameters)[2][0]
        except (ValueError, ArithmeticError, OverflowError):
            return 1e100

    optimized = deterministic_multistart(objective, initial, bounds, maxiter=1200)
    sigma_scaled, covariance, profile = evaluate(optimized.x)
    likelihood, mean, residual, factor, mean_factor, free_mean_indices = profile
    posterior = _posterior(
        analysis_data,
        geometry,
        geometry.observed_shared_depth,
        geometry.node_shared_depth,
        geometry.depths,
        sigma_scaled,
        mean,
        residual,
        factor,
        mean_factor,
        free_mean_indices,
        flat_root=True,
    )
    sigma = _restore_sigma(data, sigma_scaled)
    likelihood -= _likelihood_scale_adjustment(data, reml=True)
    # Rank is invariant to trait units.  ``sigma_scaled`` is the covariance in
    # independently normalized trait coordinates; using the restored covariance
    # would let one large-unit trait hide valid directions in smaller-unit traits.
    eigenvalues = np.linalg.eigvalsh(sigma_scaled)
    tolerance = (
        np.finfo(float).eps * max(1.0, float(np.max(eigenvalues))) * max(100, dimension)
    )
    rank = int(np.sum(eigenvalues > tolerance))
    fit = DenseMultivariateFit(
        trait_names=data.trait_names,
        sigma=sigma,
        sigma_rank=rank,
        sigma_estimated=True,
        restricted_log_likelihood=likelihood,
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
    )
    return posterior, fit


def fit_dense_mvou(
    tree,
    values_by_leaf,
    trait_names,
    *,
    alpha=None,
    alpha_bounds=None,
    standard_errors=None,
):
    """Fit stationary correlated-trait MV-OU with one shared attraction rate."""

    data = _prepare_observations(
        tree, values_by_leaf, trait_names, standard_errors=standard_errors
    )
    dimension = len(data.trait_names)
    geometry = _geometry(data)
    maximum_depth = float(np.max(geometry.depths))
    time_scale = maximum_depth if maximum_depth > 0.0 else 1.0
    bounds_alpha = (
        (1e-6 / time_scale, 50.0 / time_scale)
        if alpha_bounds is None
        else (float(alpha_bounds[0]), float(alpha_bounds[1]))
    )
    if (
        not all(math.isfinite(value) and value > 0.0 for value in bounds_alpha)
        or bounds_alpha[0] >= bounds_alpha[1]
    ):
        raise ValueError("MV-OU alpha bounds must be increasing and positive.")
    fixed_alpha = None if alpha is None else _finite(alpha, "--alpha")
    if fixed_alpha is not None and fixed_alpha <= 0.0:
        raise ValueError("MV-OU alpha must be positive.")
    if fixed_alpha is None:
        _validate_mvou_alpha_design(data, geometry.observed_distance)
    covariance_parameters = dimension * (dimension + 1) // 2
    free_parameters = dimension + covariance_parameters + int(fixed_alpha is None)
    if len(data.values) < free_parameters:
        raise ValueError(
            "MV-OU has too few effective observed coordinates to estimate all "
            "optima, covariance parameters, and free alpha."
        )
    initial, bounds = _initial_cholesky(dimension)
    if fixed_alpha is None:
        initial.insert(0, math.log(math.sqrt(bounds_alpha[0] * bounds_alpha[1])))
        bounds.insert(0, (math.log(bounds_alpha[0]), math.log(bounds_alpha[1])))

    def evaluate(parameters):
        offset = 0
        current_alpha = fixed_alpha
        if current_alpha is None:
            current_alpha = math.exp(float(parameters[0]))
            offset = 1
        sigma, _ = _decode_cholesky(parameters, dimension, offset)
        scalar = np.exp(-current_alpha * geometry.observed_distance)
        covariance = _observation_covariance(data, scalar, sigma)
        return (
            current_alpha,
            sigma,
            covariance,
            _profile_mean(data, covariance, reml=False),
        )

    def objective(parameters):
        try:
            return -evaluate(parameters)[3][0]
        except (ValueError, ArithmeticError, OverflowError):
            return 1e100

    optimized = deterministic_multistart(objective, initial, bounds, maxiter=1500)
    fitted_alpha, sigma_scaled, covariance, profile = evaluate(optimized.x)
    (
        likelihood,
        theta_scaled,
        residual,
        factor,
        mean_factor,
        free_mean_indices,
    ) = profile
    scalar_observed = np.exp(-fitted_alpha * geometry.observed_distance)
    scalar_cross = np.exp(-fitted_alpha * geometry.node_distance)
    scalar_variance = np.ones(len(geometry.compiled.nodes), dtype=float)
    posterior = _posterior(
        data,
        geometry,
        scalar_observed,
        scalar_cross,
        scalar_variance,
        sigma_scaled,
        theta_scaled,
        residual,
        factor,
        mean_factor,
        free_mean_indices,
        flat_root=False,
    )
    sigma = _restore_sigma(data, sigma_scaled)
    theta = data.centers + data.scales * theta_scaled
    likelihood -= _likelihood_scale_adjustment(data, reml=False)
    # Determine rank before restoring heterogeneous trait units (see MV-BM above).
    eigenvalues = np.linalg.eigvalsh(sigma_scaled)
    tolerance = (
        np.finfo(float).eps * max(1.0, float(np.max(eigenvalues))) * max(100, dimension)
    )
    rank = int(np.sum(eigenvalues > tolerance))
    status = "ok"
    if fixed_alpha is None:
        boundary_tolerance = 1e-5
        if fitted_alpha <= bounds_alpha[0] * (1.0 + boundary_tolerance):
            status = "alpha_lower_boundary"
        elif fitted_alpha >= bounds_alpha[1] * (1.0 - boundary_tolerance):
            status = "alpha_upper_boundary"
    fit = DenseMultivariateFit(
        trait_names=data.trait_names,
        sigma=sigma,
        sigma_rank=rank,
        sigma_estimated=True,
        restricted_log_likelihood=None,
        log_likelihood=likelihood,
        num_observed=data.num_observed_tips,
        num_effective_observations=data.num_effective_positions,
        residual_df=len(data.values),
        fit_status=status if rank == dimension else "singular_covariance",
        optimizer_success=optimized.success,
        optimizer_message=optimized.message,
        optimizer_starts=optimized.starts,
        optimizer_converged_starts=optimized.converged_starts,
        optimizer_failed_starts=optimized.failed_starts,
        model="MV-OU",
        alpha=fitted_alpha,
        alpha_estimated=fixed_alpha is None,
        theta=theta,
        theta_estimated=True,
    )
    return posterior, fit


def fit_dense_mvou_diag(
    tree,
    values_by_leaf,
    trait_names,
    *,
    alpha=None,
    alpha_by_trait=None,
    alpha_bounds=None,
    standard_errors=None,
):
    """Fit stationary MV-OU with diagonal trait-specific attraction rates."""

    data = _prepare_observations(
        tree, values_by_leaf, trait_names, standard_errors=standard_errors
    )
    dimension = len(data.trait_names)
    geometry = _geometry(data)
    maximum_depth = float(np.max(geometry.depths))
    time_scale = maximum_depth if maximum_depth > 0.0 else 1.0
    scaled_geometry = replace(
        geometry,
        depths=geometry.depths / time_scale,
        observed_shared_depth=geometry.observed_shared_depth / time_scale,
        observed_distance=geometry.observed_distance / time_scale,
        node_shared_depth=geometry.node_shared_depth / time_scale,
        node_distance=geometry.node_distance / time_scale,
    )

    has_alpha_by_trait = alpha_by_trait is not None and not (
        isinstance(alpha_by_trait, str) and alpha_by_trait == ""
    )
    if alpha is not None and has_alpha_by_trait:
        raise ValueError("--alpha cannot be combined with --alpha-by-trait.")
    shared_alpha = None
    fixed_alpha = None
    if alpha is not None:
        shared_alpha = _finite(alpha, "--alpha")
        if shared_alpha <= 0.0:
            raise ValueError("MV-OU-DIAG alpha must be positive.")
        fixed_alpha = np.full(dimension, shared_alpha, dtype=float)
    elif has_alpha_by_trait:
        fixed_alpha = parse_alpha_by_trait(alpha_by_trait, data.trait_names)

    physical_bounds = (
        (1e-6 / time_scale, 50.0 / time_scale)
        if alpha_bounds is None
        else (float(alpha_bounds[0]), float(alpha_bounds[1]))
    )
    if (
        not all(math.isfinite(value) and value > 0.0 for value in physical_bounds)
        or physical_bounds[0] >= physical_bounds[1]
    ):
        raise ValueError("MV-OU-DIAG alpha bounds must be increasing and positive.")
    beta_bounds = (
        physical_bounds[0] * time_scale,
        physical_bounds[1] * time_scale,
    )
    if not all(math.isfinite(value) and value > 0.0 for value in beta_bounds):
        raise ValueError("MV-OU-DIAG scaled alpha bounds exceed floating-point range.")

    free_alpha = fixed_alpha is None
    if free_alpha:
        _validate_mvou_diag_alpha_design(data, geometry)
    covariance_parameters = dimension * (dimension + 1) // 2
    free_parameters = (
        dimension + covariance_parameters + (dimension if free_alpha else 0)
    )
    if len(data.values) < free_parameters:
        raise ValueError(
            "MV-OU-DIAG has too few effective observed coordinates to estimate "
            "all optima, diffusion-covariance parameters, and trait-specific "
            "attraction rates."
        )

    if free_alpha:
        initial_beta = np.full(
            dimension, math.sqrt(beta_bounds[0] * beta_bounds[1]), dtype=float
        )
        initial = [math.log(value) for value in initial_beta]
        bounds = [
            (math.log(beta_bounds[0]), math.log(beta_bounds[1]))
            for _ in range(dimension)
        ]
    else:
        assert fixed_alpha is not None
        initial_beta = fixed_alpha * time_scale
        if np.any(~np.isfinite(initial_beta)):
            raise ValueError(
                "MV-OU-DIAG fixed alpha values exceed floating-point range after "
                "tree-time scaling."
            )
        initial = []
        bounds = []
    diffusion_initial, diffusion_bounds = _initial_cholesky(dimension)
    triangular_index = 0
    for row in range(dimension):
        for column in range(row + 1):
            if row == column:
                log_diagonal = 0.5 * math.log(2.0 * initial_beta[row])
                diffusion_initial[triangular_index] = log_diagonal
                lower, upper = diffusion_bounds[triangular_index]
                assert lower is not None and upper is not None
                diffusion_bounds[triangular_index] = (
                    min(lower, log_diagonal - 10.0),
                    max(upper, log_diagonal + 10.0),
                )
            triangular_index += 1
    initial.extend(diffusion_initial)
    bounds.extend(diffusion_bounds)

    def evaluate(parameters):
        offset = 0
        if free_alpha:
            beta = np.exp(np.asarray(parameters[:dimension], dtype=float))
            offset = dimension
        else:
            beta = initial_beta
        diffusion_scaled_time, _ = _decode_cholesky(parameters, dimension, offset)
        stationary_sigma = diffusion_scaled_time / (beta[:, None] + beta[None, :])
        stationary_sigma = (stationary_sigma + stationary_sigma.T) / 2.0
        covariance = _diagonal_ou_observation_covariance(
            data, scaled_geometry, beta, stationary_sigma
        )
        return (
            beta,
            diffusion_scaled_time,
            stationary_sigma,
            covariance,
            _profile_mean(data, covariance, reml=False),
        )

    def objective(parameters):
        try:
            return -evaluate(parameters)[4][0]
        except (ValueError, ArithmeticError, OverflowError):
            return 1e100

    optimized = deterministic_multistart(objective, initial, bounds, maxiter=1800)
    beta, diffusion_scaled_time, sigma_scaled, _covariance, profile = evaluate(
        optimized.x
    )
    (
        likelihood,
        theta_scaled,
        residual,
        factor,
        _mean_factor,
        _free_mean_indices,
    ) = profile
    posterior = _diagonal_ou_posterior(
        data,
        scaled_geometry,
        beta,
        sigma_scaled,
        theta_scaled,
        residual,
        factor,
    )
    sigma = _restore_sigma(data, sigma_scaled)
    diffusion_sigma = _restore_sigma(data, diffusion_scaled_time / time_scale)
    fitted_alpha = beta / time_scale
    theta = data.centers + data.scales * theta_scaled
    likelihood -= _likelihood_scale_adjustment(data, reml=False)
    eigenvalues = np.linalg.eigvalsh(sigma_scaled)
    tolerance = (
        np.finfo(float).eps * max(1.0, float(np.max(eigenvalues))) * max(100, dimension)
    )
    rank = int(np.sum(eigenvalues > tolerance))
    statuses = []
    if rank < dimension:
        statuses.append("singular_covariance")
    if free_alpha:
        boundary_tolerance = 1e-5
        if np.any(fitted_alpha <= physical_bounds[0] * (1.0 + boundary_tolerance)):
            statuses.append("alpha_lower_boundary")
        if np.any(fitted_alpha >= physical_bounds[1] * (1.0 - boundary_tolerance)):
            statuses.append("alpha_upper_boundary")
    fit = DenseMultivariateFit(
        trait_names=data.trait_names,
        sigma=sigma,
        sigma_rank=rank,
        sigma_estimated=True,
        restricted_log_likelihood=None,
        log_likelihood=likelihood,
        num_observed=data.num_observed_tips,
        num_effective_observations=data.num_effective_positions,
        residual_df=len(data.values),
        fit_status="+".join(statuses) if statuses else "ok",
        optimizer_success=optimized.success,
        optimizer_message=optimized.message,
        optimizer_starts=optimized.starts,
        optimizer_converged_starts=optimized.converged_starts,
        optimizer_failed_starts=optimized.failed_starts,
        model="MV-OU-DIAG",
        alpha=shared_alpha,
        alpha_by_trait=fitted_alpha,
        alpha_estimated=free_alpha,
        diffusion_sigma=diffusion_sigma,
        theta=theta,
        theta_estimated=True,
    )
    return posterior, fit
