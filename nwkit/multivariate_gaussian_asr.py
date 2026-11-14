"""Dense multivariate Gaussian conditioning for incomplete/noisy trait data."""

import math
from dataclasses import dataclass, replace

import numpy as np
from scipy.linalg import cho_factor, cho_solve, solve_triangular

from nwkit.clade_index import LcaIndex
from nwkit.compiled_tree import CompiledTree
from nwkit.multivariate_asr import MultivariateGaussianMarginal
from nwkit.optimization import FitResourceError, deterministic_multistart

_LOG_2PI = math.log(2.0 * math.pi)
_MAX_DENSE_OBSERVATIONS = 1000
_MAX_POSTERIOR_BYTES = 256 * 1024**2


def _stationary_diagonal_ou_covariance(covariance, log_alpha):
    # Diffusion Cholesky rows are expressed relative to sqrt(2 * alpha_i).
    # The stationary covariance then needs only the bounded ratio
    # 2 * sqrt(alpha_i * alpha_j) / (alpha_i + alpha_j).
    half_gap = 0.5 * np.abs(log_alpha[:, None] - log_alpha[None, :])
    return covariance * (2.0 * np.exp(-half_gap) / (1.0 + np.exp(-2.0 * half_gap)))


def _diagonal_ou_decay(first_rate, first_path, second_rate, second_path):
    with np.errstate(over="ignore"):
        return np.exp(-first_rate * first_path - second_rate * second_path)


def _restore_diagonal_ou_diffusion(data, covariance, log_alpha):
    log_scale = np.log(data.scales) + 0.5 * (math.log(2.0) + log_alpha)
    with np.errstate(divide="ignore", over="ignore"):
        log_values = (
            np.log(np.abs(covariance)) + log_scale[:, None] + log_scale[None, :]
        )
        restored = np.copysign(np.exp(log_values), covariance)
    if np.any(~np.isfinite(restored)):
        raise ValueError(
            "MV-OU-DIAG diffusion covariance exceeds floating-point range."
        )
    return restored


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
    observed_left: np.ndarray
    observed_right: np.ndarray
    unique_observed_nodes: np.ndarray
    observed_node_inverse: np.ndarray
    lca: LcaIndex
    jump_lengths: tuple[np.ndarray, ...]
    time_scale: float = 1.0

    def path_length(self, node, ancestor):
        steps = self.lca.depth[node] - self.lca.depth[ancestor]
        terms = []
        for level in range(steps.bit_length()):
            if steps & (1 << level):
                terms.append(float(self.jump_lengths[level][node]))
                node = self.lca.ancestors[level][node]
        return math.fsum(terms) / self.time_scale

    def cross(self, node):
        """One node's geometry, using O(observations) temporary storage."""
        shared = np.empty(len(self.unique_observed_nodes))
        left = np.empty_like(shared)
        right = np.empty_like(shared)
        for index, observed in enumerate(self.unique_observed_nodes):
            ancestor = self.lca.common_ancestor_indices(node, int(observed))
            shared[index] = self.depths[ancestor]
            left[index] = self.path_length(node, ancestor)
            right[index] = self.path_length(int(observed), ancestor)
        inverse = self.observed_node_inverse
        return shared[inverse], left[inverse], right[inverse]

    def rescaled(self, scale):
        return replace(
            self,
            depths=self.depths / scale,
            observed_shared_depth=self.observed_shared_depth / scale,
            observed_distance=self.observed_distance / scale,
            observed_left=self.observed_left / scale,
            observed_right=self.observed_right / scale,
            time_scale=self.time_scale * scale,
        )


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
    lengths = np.zeros(len(compiled.nodes), dtype=float)
    for index in range(1, len(compiled.nodes)):
        length = _finite(
            compiled.nodes[index].dist,
            "Multivariate Gaussian branch lengths",
            nonnegative=True,
        )
        depths[index] = depths[compiled.parents[index]] + length
        lengths[index] = length
    if np.any(~np.isfinite(depths)):
        raise ValueError("Multivariate tree depths exceed floating-point range.")
    jumps = [lengths]
    for ancestors in lca.ancestors[:-1]:
        previous = jumps[-1]
        jumps.append(previous + previous[np.asarray(ancestors)])

    unique_nodes, inverse_nodes = np.unique(data.node_indices, return_inverse=True)
    count = len(unique_nodes)
    shared = np.empty((count, count), dtype=float)
    distance = np.empty((count, count), dtype=float)
    left = np.empty_like(distance)
    right = np.empty_like(distance)
    geometry = _Geometry(
        compiled,
        depths,
        shared,
        distance,
        left,
        right,
        unique_nodes,
        inverse_nodes,
        lca,
        tuple(jumps),
    )
    lookup_ancestors = lca.indexed_ancestors()
    path_lengths = _prepare_indexed_path_lengths(lca, jumps)
    batch_size = max(1, 32_768 // max(1, count))
    for start in range(0, count, batch_size):
        stop = min(start + batch_size, count)
        first = unique_nodes[start:stop, None]
        second = unique_nodes[None, :]
        ancestors = lookup_ancestors(first, second)
        shared[start:stop] = depths[ancestors]
        left[start:stop] = path_lengths(first, ancestors)
        right[start:stop] = path_lengths(second, ancestors)
        distance[start:stop] = left[start:stop] + right[start:stop]
    if np.any(~np.isfinite(distance)):
        raise ValueError("Multivariate tree distances exceed floating-point range.")
    return _expand_multivariate_geometry(geometry)


def _prepare_indexed_path_lengths(lca, jumps):
    depths = np.asarray(lca.depth, dtype=np.intp)
    ancestors = np.asarray(lca.ancestors, dtype=np.intp)

    def path_lengths(nodes, ancestor):
        steps = depths[nodes] - depths[ancestor]
        result = np.zeros(steps.shape, dtype=float)
        correction = np.zeros_like(result)
        current = np.broadcast_to(nodes, steps.shape)
        for level, lengths in enumerate(jumps):
            lift = steps & (1 << level)
            length = np.where(lift, lengths[current], 0.0)
            # Sum the selected jumps directly. Taking a difference of node
            # depths can erase a short branch near the end of a long path.
            total = result + length
            correction += np.where(
                result >= length, (result - total) + length, (length - total) + result
            )
            result = total
            current = np.where(lift, ancestors[level, current], current)
        return result + correction

    return path_lengths


def _expand_multivariate_geometry(geometry):
    inverse = geometry.observed_node_inverse
    if np.array_equal(inverse, np.arange(len(geometry.unique_observed_nodes))):
        return geometry
    indices = np.ix_(inverse, inverse)
    return replace(
        geometry,
        observed_shared_depth=geometry.observed_shared_depth[indices],
        observed_distance=geometry.observed_distance[indices],
        observed_left=geometry.observed_left[indices],
        observed_right=geometry.observed_right[indices],
    )


def _cached_geometry(data, cache):
    if cache is None:
        return _geometry(data)
    key = (data.compiled.tree, tuple(data.node_indices))
    if key not in cache:
        cache[key] = _geometry(data)
    return cache[key]


def _validate_posterior_size(data, compute_posterior):
    dimension = len(data.trait_names)
    estimated_bytes = len(data.compiled.nodes) * (256 + 8 * (dimension + dimension**2))
    if compute_posterior and estimated_bytes > _MAX_POSTERIOR_BYTES:
        raise FitResourceError(
            "Multivariate all-node posterior exceeds the 256 MiB estimated "
            "result-storage limit; use asrcompare for likelihood-only fits "
            "or reduce the reconstruction tree."
        )


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
    tolerance = np.finfo(float).eps * scale * 100
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
    tolerance = np.finfo(float).eps * scale * 100
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
        tolerance = np.finfo(float).eps * scale * 100
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


def _prepare_mean_profile(data, *, fixed_mean=None):
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
    adjusted_values = data.values - full_design @ mean
    rhs = np.column_stack([design, adjusted_values])
    return full_design, mean, free_mean_indices, rhs


def _profile_mean(
    data, covariance, *, reml, fixed_mean=None, prepared=None, return_residual=True
):
    if prepared is None:
        prepared = _prepare_mean_profile(data, fixed_mean=fixed_mean)
    full_design, fixed, free_mean_indices, rhs = prepared
    mean = fixed.copy()
    factor, log_determinant = _factor_covariance(covariance)
    whitened = solve_triangular(
        factor[0],
        rhs,
        lower=True,
        check_finite=False,
    )
    whitened_values = whitened[:, -1]
    mean_factor = None
    mean_log_determinant = 0.0
    if len(free_mean_indices):
        whitened_design = whitened[:, :-1]
        information = whitened_design.T @ whitened_design
        try:
            mean_factor = cho_factor(information, lower=True, check_finite=False)
        except np.linalg.LinAlgError as exc:
            raise ValueError(
                "Multivariate trait means are not separately identifiable."
            ) from exc
        fitted = cho_solve(
            mean_factor, whitened_design.T @ whitened_values, check_finite=False
        )
        mean[free_mean_indices] = fitted
        whitened_values = whitened_values - whitened_design @ fitted
        mean_log_determinant = 2.0 * float(np.log(np.diag(mean_factor[0])).sum())
    residual = data.values - full_design @ mean if return_residual else None
    quadratic = float(whitened_values @ whitened_values)
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
    sigma,
    mean,
    residual,
    factor,
    mean_factor,
    free_mean_indices,
    *,
    flat_root,
    alpha=None,
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
    trait_covariance = sigma[:, data.trait_indices]

    def cross_covariance(shared, left, right):
        cross = shared if flat_root else np.exp(-alpha * (left + right))
        return cross[:, None, :] * trait_covariance

    adjustment_base = np.zeros((dimension, len(free_mean_indices)), dtype=float)
    adjustment_base[free_mean_indices, np.arange(len(free_mean_indices))] = 1.0
    prior_scales = (
        geometry.depths if flat_root else np.ones(len(geometry.compiled.nodes))
    )
    posterior = {}
    for start, nodes, cross, conditional_covariance in _posterior_covariance_batches(
        data, geometry, factor, cross_covariance, need_paths=not flat_root
    ):
        scalars = prior_scales[start : start + len(nodes), None, None]
        covariance = scalars * sigma - conditional_covariance
        if mean_covariance is not None:
            adjustment = adjustment_base - cross @ solved_design
            covariance += adjustment @ mean_covariance @ adjustment.swapaxes(1, 2)
        means = mean + cross @ solved_residual
        posterior.update(_restored_posterior_marginals(data, nodes, means, covariance))
    return posterior


def _posterior_covariance_batches(
    data, geometry, factor, cross_covariance, *, need_paths=True
):
    dimension = len(data.trait_names)
    observations = len(data.values)
    batch_size = max(1, 32_768 // (dimension * observations))
    all_nodes = geometry.compiled.nodes
    lookup = geometry.lca.indexed_ancestors()
    path_lengths = _prepare_indexed_path_lengths(geometry.lca, geometry.jump_lengths)
    observed = geometry.unique_observed_nodes[None, :]
    inverse = geometry.observed_node_inverse
    for start in range(0, len(all_nodes), batch_size):
        nodes = all_nodes[start : start + batch_size]
        indices = np.arange(start, start + len(nodes))[:, None]
        ancestor = lookup(indices, observed)
        shared = geometry.depths[ancestor][:, inverse]
        left, right = None, None
        if need_paths:
            left = path_lengths(indices, ancestor)[:, inverse] / geometry.time_scale
            right = path_lengths(observed, ancestor)[:, inverse] / geometry.time_scale
        cross = cross_covariance(shared, left, right)
        solved_cross = cho_solve(
            factor, cross.reshape(-1, observations).T, check_finite=False
        )
        yield (
            start,
            nodes,
            cross,
            cross @ solved_cross.T.reshape(cross.shape).swapaxes(1, 2),
        )


def _restored_posterior_marginals(data, nodes, means, covariance):
    covariance = (covariance + covariance.swapaxes(1, 2)) / 2.0
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    tolerance = (
        np.finfo(float).eps
        * np.maximum(1.0, np.max(np.abs(eigenvalues), axis=1))
        * max(100, len(data.trait_names))
    )
    if np.any(eigenvalues[:, 0] < -tolerance):
        raise ValueError(
            "A multivariate posterior covariance is not positive semidefinite."
        )
    covariance = (eigenvectors * np.maximum(eigenvalues, 0.0)[:, None, :]) @ (
        eigenvectors.swapaxes(1, 2)
    )
    return {
        node: MultivariateGaussianMarginal(
            data.centers + data.scales * node_mean,
            _restore_sigma(data, node_covariance),
        )
        for node, node_mean, node_covariance in zip(
            nodes, means, covariance, strict=True
        )
    }


def _diagonal_ou_observation_covariance(data, geometry, beta, stationary_sigma):
    left = geometry.observed_left
    right = geometry.observed_right
    observed_beta = beta[data.trait_indices]
    decay = _diagonal_ou_decay(
        observed_beta[:, None], left, observed_beta[None, :], right
    )
    covariance = (
        decay
        * stationary_sigma[data.trait_indices[:, None], data.trait_indices[None, :]]
    )
    covariance[np.diag_indices_from(covariance)] += data.errors * data.errors
    return (covariance + covariance.T) / 2.0


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
    observed_beta = beta[data.trait_indices]
    trait_covariance = stationary_sigma[:, data.trait_indices]
    posterior = {}

    def cross_covariance(_, left, right):
        decay = _diagonal_ou_decay(
            beta[:, None], left[:, None, :], observed_beta[None, :], right[:, None, :]
        )
        return decay * trait_covariance

    for _, nodes, cross, conditional_covariance in _posterior_covariance_batches(
        data, geometry, factor, cross_covariance
    ):
        posterior.update(
            _restored_posterior_marginals(
                data,
                nodes,
                mean + cross @ solved_residual,
                stationary_sigma - conditional_covariance,
            )
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
    compute_posterior=True,
    _geometry_cache=None,
):
    """Fit flat-root MV-BM with arbitrary per-trait missingness and errors."""

    data = _prepare_observations(
        tree, values_by_leaf, trait_names, standard_errors=standard_errors
    )
    dimension = len(data.trait_names)
    _validate_posterior_size(data, compute_posterior)
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
    geometry = _cached_geometry(analysis_data, _geometry_cache)
    time_scale = float(np.max(geometry.observed_shared_depth)) or 1.0
    geometry = geometry.rescaled(time_scale)
    _validate_mvbm_covariance_design(
        analysis_data, geometry.observed_shared_depth, fixed_mean
    )
    initial, bounds = _initial_cholesky(dimension)
    mean_profile = _prepare_mean_profile(analysis_data, fixed_mean=fixed_mean)

    def evaluate(parameters, *, return_residual=True):
        sigma, _ = _decode_cholesky(parameters, dimension)
        covariance = _observation_covariance(
            analysis_data, geometry.observed_shared_depth, sigma
        )
        return (
            sigma,
            covariance,
            _profile_mean(
                analysis_data,
                covariance,
                reml=True,
                prepared=mean_profile,
                return_residual=return_residual,
            ),
        )

    def objective(parameters):
        try:
            return -evaluate(parameters, return_residual=False)[2][0]
        except (ValueError, ArithmeticError, OverflowError):
            return 1e100

    optimized = deterministic_multistart(objective, initial, bounds, maxiter=1200)
    sigma_scaled, covariance, profile = evaluate(optimized.x)
    likelihood, mean, residual, factor, mean_factor, free_mean_indices = profile
    posterior = (
        _posterior(
            analysis_data,
            geometry,
            sigma_scaled,
            mean,
            residual,
            factor,
            mean_factor,
            free_mean_indices,
            flat_root=True,
        )
        if compute_posterior
        else {}
    )
    sigma = _restore_sigma(data, sigma_scaled / time_scale)
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
    compute_posterior=True,
    _geometry_cache=None,
):
    """Fit stationary correlated-trait MV-OU with one shared attraction rate."""

    data = _prepare_observations(
        tree, values_by_leaf, trait_names, standard_errors=standard_errors
    )
    dimension = len(data.trait_names)
    _validate_posterior_size(data, compute_posterior)
    geometry = _cached_geometry(data, _geometry_cache)
    # Stationary tip likelihoods do not depend on a shared stem above their
    # MRCA. Use observable distances for units and default attraction bounds.
    time_scale = float(np.max(geometry.observed_distance)) / 2.0 or 1.0
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
    observed_distance = geometry.observed_distance
    if fixed_alpha is None:
        observed_distance = observed_distance / time_scale
        _validate_mvou_alpha_design(data, observed_distance)
    covariance_parameters = dimension * (dimension + 1) // 2
    free_parameters = dimension + covariance_parameters + int(fixed_alpha is None)
    if len(data.values) < free_parameters:
        raise ValueError(
            "MV-OU has too few effective observed coordinates to estimate all "
            "optima, covariance parameters, and free alpha."
        )
    initial, bounds = _initial_cholesky(dimension)
    if fixed_alpha is None:
        log_beta_bounds = (
            math.log(bounds_alpha[0]) + math.log(time_scale),
            math.log(bounds_alpha[1]) + math.log(time_scale),
        )
        initial.insert(0, 0.5 * sum(log_beta_bounds))
        bounds.insert(0, log_beta_bounds)
    mean_profile = _prepare_mean_profile(data)

    def evaluate(parameters, *, return_residual=True):
        offset = 0
        decay_rate = fixed_alpha
        if decay_rate is None:
            decay_rate = math.exp(float(parameters[0]))
            offset = 1
        sigma, _ = _decode_cholesky(parameters, dimension, offset)
        scalar = np.exp(-decay_rate * observed_distance)
        covariance = _observation_covariance(data, scalar, sigma)
        return (
            decay_rate,
            sigma,
            covariance,
            _profile_mean(
                data,
                covariance,
                reml=False,
                prepared=mean_profile,
                return_residual=return_residual,
            ),
        )

    def objective(parameters):
        try:
            return -evaluate(parameters, return_residual=False)[3][0]
        except (ValueError, ArithmeticError, OverflowError):
            return 1e100

    optimized = deterministic_multistart(objective, initial, bounds, maxiter=1500)
    fitted_decay, sigma_scaled, covariance, profile = evaluate(optimized.x)
    fitted_alpha = fixed_alpha if fixed_alpha is not None else fitted_decay / time_scale
    (
        likelihood,
        theta_scaled,
        residual,
        factor,
        mean_factor,
        free_mean_indices,
    ) = profile
    posterior = (
        _posterior(
            data,
            geometry,
            sigma_scaled,
            theta_scaled,
            residual,
            factor,
            mean_factor,
            free_mean_indices,
            flat_root=False,
            alpha=fitted_alpha,
        )
        if compute_posterior
        else {}
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
    compute_posterior=True,
    _geometry_cache=None,
):
    """Fit stationary MV-OU with diagonal trait-specific attraction rates."""

    data = _prepare_observations(
        tree, values_by_leaf, trait_names, standard_errors=standard_errors
    )
    dimension = len(data.trait_names)
    _validate_posterior_size(data, compute_posterior)
    geometry = _cached_geometry(data, _geometry_cache)
    time_scale = float(np.max(geometry.observed_distance)) / 2.0 or 1.0

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
        log_beta_bounds = (
            math.log(physical_bounds[0]) + math.log(time_scale),
            math.log(physical_bounds[1]) + math.log(time_scale),
        )
        log_initial_beta = 0.5 * sum(log_beta_bounds)
        initial = [log_initial_beta] * dimension
        bounds = [log_beta_bounds] * dimension
        analysis_geometry = geometry.rescaled(time_scale)
    else:
        assert fixed_alpha is not None
        initial = []
        bounds = []
        analysis_geometry = geometry
    diffusion_initial, diffusion_bounds = _initial_cholesky(dimension)
    initial.extend(diffusion_initial)
    bounds.extend(diffusion_bounds)
    mean_profile = _prepare_mean_profile(data)

    def evaluate(parameters, *, return_residual=True):
        offset = 0
        if free_alpha:
            log_decay = np.asarray(parameters[:dimension], dtype=float)
            beta = np.exp(log_decay)
            offset = dimension
        else:
            assert fixed_alpha is not None
            beta = fixed_alpha
            log_decay = np.log(beta)
        correlation_scaled_diffusion, _ = _decode_cholesky(
            parameters, dimension, offset
        )
        stationary_sigma = _stationary_diagonal_ou_covariance(
            correlation_scaled_diffusion, log_decay
        )
        covariance = _diagonal_ou_observation_covariance(
            data, analysis_geometry, beta, stationary_sigma
        )
        return (
            beta,
            correlation_scaled_diffusion,
            stationary_sigma,
            covariance,
            _profile_mean(
                data,
                covariance,
                reml=False,
                prepared=mean_profile,
                return_residual=return_residual,
            ),
        )

    def objective(parameters):
        try:
            return -evaluate(parameters, return_residual=False)[4][0]
        except (ValueError, ArithmeticError, OverflowError):
            return 1e100

    optimized = deterministic_multistart(objective, initial, bounds, maxiter=1800)
    beta, normalized_diffusion, sigma_scaled, _covariance, profile = evaluate(
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
    posterior = (
        _diagonal_ou_posterior(
            data,
            analysis_geometry,
            beta,
            sigma_scaled,
            theta_scaled,
            residual,
            factor,
        )
        if compute_posterior
        else {}
    )
    sigma = _restore_sigma(data, sigma_scaled)
    fitted_alpha = beta / time_scale if free_alpha else beta
    diffusion_sigma = _restore_diagonal_ou_diffusion(
        data, normalized_diffusion, np.log(fitted_alpha)
    )
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
