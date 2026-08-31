"""Small linear-algebra helpers for Gaussian covariance calculations."""

import math
from dataclasses import dataclass

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import SuperLU, splu


def center_response(response, design, *, excluded_columns=()):
    """Remove an explicit intercept offset without changing a no-intercept fit.

    Center before solving or taking numerical derivatives: subtracting a large
    fitted intercept afterwards loses useful precision. An uncertain predictor
    cannot be used as the intercept in an errors-in-variables model.
    """
    response = np.asarray(response, dtype=float)
    offset = np.zeros(design.shape[1], dtype=float)
    for column in range(design.shape[1]):
        values = design[:, column]
        if (
            column not in excluded_columns
            and len(values)
            and values[0] != 0.0
            and np.all(values == values[0])
        ):
            mean = float(np.mean(response))
            offset[column] = mean / float(values[0])
            return response - mean, offset
    return response, offset


def residual_variance_scale(response, residual, sampling_diagonal):
    """Choose variance bounds in residual units, independent of a fitted mean."""
    sampling_scale = float(np.mean(sampling_diagonal))
    residual_scale = float(np.mean(np.square(residual)))
    # The floor follows response units and only affects a numerically exact
    # fit. The response passed here has already had its intercept removed.
    return max(
        residual_scale,
        sampling_scale,
        np.finfo(float).eps * float(np.mean(np.square(response))),
        # For an exactly constant response with no sampling error, the ML
        # optimum is zero. Keep the limiting covariance invertible; callers
        # explicitly report the resulting variance-boundary warning.
        math.sqrt(np.finfo(float).tiny),
    )


@dataclass(frozen=True)
class DiagonalLowRankFactor:
    """Factorization of ``diag(diagonal) + low_rank @ low_rank.T``."""

    diagonal: np.ndarray
    low_rank: np.ndarray
    woodbury_cholesky: np.ndarray


@dataclass(frozen=True)
class GroupedDiagonalLowRankFactor:
    """Diagonal covariance with disjoint group and general low-rank updates."""

    diagonal: np.ndarray
    group_index: np.ndarray
    group_loading: np.ndarray
    group_denominator: np.ndarray
    low_rank: np.ndarray
    base_inverse_low_rank: np.ndarray
    woodbury_cholesky: np.ndarray


@dataclass(frozen=True)
class DiagonalLowRankCovariance:
    """Unmaterialized diagonal-plus-low-rank covariance."""

    diagonal: np.ndarray
    low_rank: np.ndarray | sparse.spmatrix


@dataclass(frozen=True)
class DiagonalSparsePrecisionCovariance:
    """Unmaterialized ``diag(d) + L @ Q^-1 @ L.T`` covariance."""

    diagonal: np.ndarray
    loading: sparse.csr_matrix
    precision: sparse.csc_matrix
    marginal_diagonal: np.ndarray | None = None
    grouped_marginal_diagonal: np.ndarray | None = None


@dataclass(frozen=True)
class SparseDiagonalLowRankFactor:
    """Sparse Woodbury factor for ``diag(d) + U @ U.T``."""

    diagonal: np.ndarray
    low_rank: sparse.csr_matrix
    woodbury_factor: SuperLU
    logdet: float


@dataclass(frozen=True)
class SparseCovarianceFactor:
    """Direct sparse factor for ``diag(d) + U @ U.T``.

    This orientation avoids forming the latent-space Woodbury matrix when a
    sparse loading has at least as many columns as observation rows.
    """

    diagonal: np.ndarray
    low_rank: sparse.csr_matrix
    covariance_factor: SuperLU
    logdet: float


@dataclass(frozen=True)
class DiagonalSparsePrecisionFactor:
    """Factor for ``diag(d) + L @ Q^-1 @ L.T`` with sparse ``L`` and ``Q``."""

    diagonal: np.ndarray
    loading: sparse.csr_matrix
    prior_precision: sparse.csc_matrix
    prior_precision_factor: sparse.csr_matrix
    posterior_factor: SuperLU
    prior_factor: SuperLU
    logdet: float


def effective_likelihood_settings(
    n_observations,
    num_parameters,
    reml,
    likelihood_observations=None,
    likelihood_logdet_offset=0.0,
):
    """Validate and normalize event-level Gaussian pseudo-likelihood settings."""
    if likelihood_observations is None:
        likelihood_observations = n_observations
    likelihood_observations = float(likelihood_observations)
    if (
        not math.isfinite(likelihood_observations)
        or likelihood_observations <= 0.0
        or likelihood_observations > n_observations
    ):
        raise ValueError(
            "likelihood_observations must be positive and no larger than the "
            "number of observations."
        )
    if reml and likelihood_observations <= num_parameters:
        raise ValueError(
            "REML requires more effective likelihood observations than fixed-effect "
            "parameters."
        )
    likelihood_logdet_offset = float(likelihood_logdet_offset)
    if not math.isfinite(likelihood_logdet_offset):
        raise ValueError("likelihood_logdet_offset must be finite.")
    effective_count = (
        likelihood_observations - num_parameters if reml else likelihood_observations
    )
    return (
        effective_count,
        likelihood_observations / n_observations,
        likelihood_logdet_offset,
    )


@dataclass(frozen=True)
class NestedLowRankFactor:
    """A factored base covariance followed by one dense low-rank update."""

    base_factor: object
    low_rank: np.ndarray
    base_inverse_low_rank: np.ndarray
    woodbury_cholesky: np.ndarray


def _positive_diagonal(values: np.ndarray) -> np.ndarray:
    diagonal = np.asarray(values, dtype=float)
    if diagonal.ndim != 1 or np.any(diagonal <= 0.0):
        raise np.linalg.LinAlgError("Covariance diagonal is not positive.")
    return diagonal


def _sparse_lu_logdet(solver: SuperLU, label: str) -> float:
    pivots = np.asarray(solver.U.diagonal(), dtype=float)
    if np.any(pivots == 0.0) or not np.isfinite(pivots).all():
        raise np.linalg.LinAlgError("{} is singular.".format(label))
    return float(np.log(np.abs(pivots)).sum())


def is_diagonal(matrix: np.ndarray) -> bool:
    """Return whether a square matrix has exactly zero off-diagonal entries."""
    matrix = np.asarray(matrix)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        return False
    for index, row in enumerate(matrix):
        if np.any(row[:index]) or np.any(row[index + 1 :]):
            return False
    return True


def factor_covariance(covariance: np.ndarray) -> np.ndarray:
    """Factor a positive-definite covariance, preserving diagonal structure."""
    covariance = np.asarray(covariance, dtype=float)
    if is_diagonal(covariance):
        diagonal = np.diag(covariance)
        if np.any(diagonal <= 0.0):
            raise np.linalg.LinAlgError("Covariance is not positive definite.")
        return np.sqrt(diagonal)
    return np.linalg.cholesky(covariance)


def factor_diagonal_low_rank(
    diagonal: np.ndarray, low_rank: np.ndarray
) -> DiagonalLowRankFactor:
    """Factor a positive diagonal covariance plus a low-rank update."""
    diagonal = _positive_diagonal(diagonal)
    low_rank = np.asarray(low_rank, dtype=float)
    if low_rank.ndim != 2 or low_rank.shape[0] != len(diagonal):
        raise ValueError("Low-rank covariance factor has the wrong dimensions.")
    weighted = low_rank / diagonal[:, None]
    woodbury = np.eye(low_rank.shape[1]) + low_rank.T @ weighted
    return DiagonalLowRankFactor(
        diagonal=diagonal,
        low_rank=low_rank,
        woodbury_cholesky=np.linalg.cholesky(woodbury),
    )


def factor_diagonal_sparse_low_rank(
    diagonal: np.ndarray,
    low_rank: sparse.spmatrix,
    *,
    force_observation_space: bool = False,
) -> SparseDiagonalLowRankFactor | SparseCovarianceFactor:
    """Factor diagonal covariance plus a sparse latent loading."""
    diagonal = _positive_diagonal(diagonal)
    loading = sparse.csr_matrix(low_rank, dtype=float)
    if loading.shape[0] != len(diagonal) or not np.isfinite(loading.data).all():
        raise ValueError("Sparse low-rank covariance factor is malformed.")
    if force_observation_space or loading.shape[1] >= loading.shape[0]:
        covariance = (
            sparse.diags(diagonal, format="csc") + loading @ loading.T
        ).tocsc()
        solver = splu(covariance, permc_spec="COLAMD")
        logdet = _sparse_lu_logdet(solver, "Sparse covariance matrix")
        return SparseCovarianceFactor(diagonal, loading, solver, logdet)
    weighted = loading.multiply((1.0 / diagonal)[:, None])
    woodbury = (
        sparse.eye(loading.shape[1], format="csc") + loading.T @ weighted
    ).tocsc()
    solver = splu(woodbury, permc_spec="COLAMD")
    logdet = float(np.log(diagonal).sum()) + _sparse_lu_logdet(
        solver, "Sparse Woodbury matrix"
    )
    return SparseDiagonalLowRankFactor(diagonal, loading, solver, logdet)


def factor_diagonal_sparse_precision_updates(
    diagonal: np.ndarray,
    updates: list[tuple[sparse.spmatrix, sparse.spmatrix, sparse.spmatrix]],
) -> DiagonalSparsePrecisionFactor:
    """Factor a positive diagonal plus sparse latent precision components.

    Each update is ``(loading, precision, precision_factor)`` and represents
    ``loading @ inv(precision) @ loading.T``.  The final element must satisfy
    ``precision_factor.T @ precision_factor == precision``; retaining it also
    permits exact structured Gaussian draws without a dense covariance.
    """
    diagonal = _positive_diagonal(diagonal)
    if not updates:
        raise ValueError("At least one sparse-precision update is required.")
    loadings = []
    precisions = []
    precision_factors = []
    for loading_value, precision_value, precision_factor_value in updates:
        loading, precision, precision_factor = _validated_sparse_precision_update(
            len(diagonal), loading_value, precision_value, precision_factor_value
        )
        loadings.append(loading)
        precisions.append(precision)
        precision_factors.append(precision_factor)
    loading = sparse.hstack(loadings, format="csr")
    prior_precision = sparse.block_diag(precisions, format="csc")
    prior_precision_factor = sparse.block_diag(precision_factors, format="csr")
    prior_factor = splu(prior_precision, permc_spec="COLAMD")
    inverse_diagonal = 1.0 / diagonal
    posterior_precision = (
        prior_precision
        + loading.T @ sparse.diags(inverse_diagonal, format="csc") @ loading
    ).tocsc()
    posterior_factor = splu(posterior_precision, permc_spec="COLAMD")
    logdet = (
        float(np.log(diagonal).sum())
        + _sparse_lu_logdet(posterior_factor, "Sparse posterior precision")
        - _sparse_lu_logdet(prior_factor, "Sparse prior precision")
    )
    return DiagonalSparsePrecisionFactor(
        diagonal=diagonal,
        loading=loading,
        prior_precision=prior_precision,
        prior_precision_factor=prior_precision_factor,
        posterior_factor=posterior_factor,
        prior_factor=prior_factor,
        logdet=logdet,
    )


def _validated_sparse_precision_update(
    observation_count,
    loading_value,
    precision_value,
    precision_factor_value,
):
    loading = sparse.csr_matrix(loading_value, dtype=float)
    precision = sparse.csc_matrix(precision_value, dtype=float)
    precision_factor = sparse.csr_matrix(precision_factor_value, dtype=float)
    if (
        loading.shape[0] != observation_count
        or precision.shape[0] != precision.shape[1]
        or loading.shape[1] != precision.shape[0]
        or precision_factor.shape[1] != precision.shape[0]
        or not np.isfinite(loading.data).all()
        or not np.isfinite(precision.data).all()
        or not np.isfinite(precision_factor.data).all()
    ):
        raise ValueError("Sparse-precision covariance update is malformed.")
    return loading, precision, precision_factor


def _factor_nested_low_rank(base_factor, low_rank: sparse.spmatrix):
    low_rank_array = np.asarray(low_rank.toarray(), dtype=float)
    base_inverse_low_rank = solve_factor(base_factor, low_rank_array)
    woodbury = (
        np.eye(low_rank_array.shape[1]) + low_rank_array.T @ base_inverse_low_rank
    )
    return NestedLowRankFactor(
        base_factor=base_factor,
        low_rank=low_rank_array,
        base_inverse_low_rank=base_inverse_low_rank,
        woodbury_cholesky=np.linalg.cholesky((woodbury + woodbury.T) / 2.0),
    )


def _group_base_solve(
    factor: GroupedDiagonalLowRankFactor, values: np.ndarray
) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        inverse_values = values / factor.diagonal
        group_rhs = np.bincount(
            factor.group_index,
            weights=factor.group_loading * inverse_values,
            minlength=len(factor.group_denominator),
        )
        correction = (
            factor.group_loading
            / factor.diagonal
            * (group_rhs / factor.group_denominator)[factor.group_index]
        )
    else:
        inverse_values = values / factor.diagonal[:, None]
        group_rhs = np.zeros(
            (len(factor.group_denominator), values.shape[1]), dtype=float
        )
        np.add.at(
            group_rhs,
            factor.group_index,
            factor.group_loading[:, None] * inverse_values,
        )
        correction = (
            factor.group_loading[:, None]
            / factor.diagonal[:, None]
            * (group_rhs / factor.group_denominator[:, None])[factor.group_index]
        )
    return inverse_values - correction


def _disjoint_column_groups(factor: np.ndarray) -> tuple[np.ndarray, np.ndarray] | None:
    nonzero = factor != 0.0
    if factor.shape[1] == 0 or np.any(nonzero.sum(axis=1) != 1):
        return None
    if np.any(nonzero.sum(axis=0) == 0):
        return None
    group_index = np.argmax(nonzero, axis=1)
    return group_index, factor[np.arange(len(factor)), group_index]


def _select_grouped_update(
    updates: list[np.ndarray],
) -> tuple[int, np.ndarray, np.ndarray] | None:
    """Select the widest update whose columns partition the observations."""
    candidates = [
        (update.shape[1], index, groups)
        for index, update in enumerate(updates)
        if (groups := _disjoint_column_groups(update)) is not None
    ]
    if not candidates:
        return None
    _, index, (group_index, group_loading) = max(candidates)
    return index, group_index, group_loading


def _merge_compatible_grouped_updates(updates: list[np.ndarray]) -> list[np.ndarray]:
    """Merge collinear updates that encode the same disjoint groups."""
    merged = list(updates)
    reference_index = 0
    while reference_index < len(merged):
        reference_groups = _disjoint_column_groups(merged[reference_index])
        if reference_groups is None:
            reference_index += 1
            continue
        group_index, group_loading = reference_groups
        candidate_index = reference_index + 1
        while candidate_index < len(merged):
            candidate_groups = _disjoint_column_groups(merged[candidate_index])
            if candidate_groups is None or not np.array_equal(
                candidate_groups[0], group_index
            ):
                candidate_index += 1
                continue
            candidate_loading = candidate_groups[1]
            if np.any(group_loading == 0.0):
                candidate_index += 1
                continue
            ratios = candidate_loading / group_loading
            ratio_by_group = np.zeros(merged[reference_index].shape[1], dtype=float)
            compatible = True
            for group in range(len(ratio_by_group)):
                selected = ratios[group_index == group]
                ratio_by_group[group] = float(selected[0])
                if not np.allclose(
                    selected,
                    ratio_by_group[group],
                    rtol=1e-12,
                    atol=1e-14,
                ):
                    compatible = False
                    break
            if not compatible:
                candidate_index += 1
                continue
            group_loading = group_loading * np.sqrt(
                1.0 + np.square(ratio_by_group[group_index])
            )
            combined = merged[reference_index].copy()
            combined[np.arange(len(combined)), group_index] = group_loading
            merged[reference_index] = combined
            del merged[candidate_index]
        reference_index += 1
    return merged


def _factor_sparse_updates(diagonal, updates):
    sparse_updates = [sparse.csr_matrix(update, dtype=float) for update in updates]
    if any(
        update.shape[0] != len(diagonal) or not np.isfinite(update.data).all()
        for update in sparse_updates
    ):
        raise ValueError("Low-rank covariance update has invalid dimensions or values.")
    sparse_base_updates = [
        update
        for update in sparse_updates
        if update.shape[1] == 0
        or update.nnz <= 0.25 * update.shape[0] * update.shape[1]
    ]
    dense_outer_updates = [
        update
        for update in sparse_updates
        if update.shape[1] > 0 and update.nnz > 0.25 * update.shape[0] * update.shape[1]
    ]
    if sparse_base_updates and dense_outer_updates:
        base_factor = factor_diagonal_sparse_low_rank(
            diagonal,
            sparse.hstack(sparse_base_updates, format="csr"),
            force_observation_space=True,
        )
        return _factor_nested_low_rank(
            base_factor, sparse.hstack(dense_outer_updates, format="csr")
        )
    return factor_diagonal_sparse_low_rank(
        diagonal, sparse.hstack(sparse_updates, format="csr")
    )


def _factor_dense_updates(diagonal, updates):
    normalized = [np.asarray(update, dtype=float) for update in updates]
    if any(
        update.ndim != 2
        or update.shape[0] != len(diagonal)
        or not np.isfinite(update).all()
        for update in normalized
    ):
        raise ValueError("Low-rank covariance update has invalid dimensions or values.")
    if not normalized:
        return factor_diagonal_low_rank(
            diagonal, np.zeros((len(diagonal), 0), dtype=float)
        )
    normalized = _merge_compatible_grouped_updates(normalized)
    grouped_update = _select_grouped_update(normalized)
    if grouped_update is None:
        return factor_diagonal_low_rank(diagonal, np.column_stack(normalized))
    grouped_index, group_index, group_loading = grouped_update
    low_rank_updates = [
        update for index, update in enumerate(normalized) if index != grouped_index
    ]
    low_rank = (
        np.column_stack(low_rank_updates)
        if low_rank_updates
        else np.zeros((len(diagonal), 0), dtype=float)
    )
    inverse_group_weight = np.square(group_loading) / diagonal
    group_denominator = 1.0 + np.bincount(
        group_index,
        weights=inverse_group_weight,
        minlength=normalized[grouped_index].shape[1],
    )
    provisional = GroupedDiagonalLowRankFactor(
        diagonal=diagonal,
        group_index=group_index,
        group_loading=group_loading,
        group_denominator=group_denominator,
        low_rank=low_rank,
        base_inverse_low_rank=np.empty_like(low_rank),
        woodbury_cholesky=np.empty((0, 0), dtype=float),
    )
    base_inverse_low_rank = _group_base_solve(provisional, low_rank)
    woodbury = np.eye(low_rank.shape[1]) + low_rank.T @ base_inverse_low_rank
    return GroupedDiagonalLowRankFactor(
        diagonal=diagonal,
        group_index=group_index,
        group_loading=group_loading,
        group_denominator=group_denominator,
        low_rank=low_rank,
        base_inverse_low_rank=base_inverse_low_rank,
        woodbury_cholesky=np.linalg.cholesky(woodbury),
    )


def factor_diagonal_low_rank_updates(
    diagonal: np.ndarray, updates: list[np.ndarray]
) -> (
    DiagonalLowRankFactor
    | GroupedDiagonalLowRankFactor
    | SparseDiagonalLowRankFactor
    | SparseCovarianceFactor
    | NestedLowRankFactor
):
    """Factor diagonal covariance plus low-rank updates, exploiting groups."""
    diagonal = _positive_diagonal(diagonal)
    if (
        any(sparse.issparse(update) for update in updates)
        or sum(update.shape[1] for update in updates) > 512
    ):
        return _factor_sparse_updates(diagonal, updates)
    return _factor_dense_updates(diagonal, updates)


def solve_factor(
    factor: np.ndarray
    | DiagonalLowRankFactor
    | GroupedDiagonalLowRankFactor
    | SparseDiagonalLowRankFactor
    | SparseCovarianceFactor
    | DiagonalSparsePrecisionFactor
    | NestedLowRankFactor,
    values: np.ndarray,
) -> np.ndarray:
    """Apply the inverse covariance represented by a Cholesky-style factor."""
    if isinstance(factor, DiagonalLowRankFactor):
        values = np.asarray(values, dtype=float)
        if values.ndim == 1:
            inverse_values = values / factor.diagonal
        else:
            inverse_values = values / factor.diagonal[:, None]
        correction_rhs = factor.low_rank.T @ inverse_values
        correction = np.linalg.solve(
            factor.woodbury_cholesky.T,
            np.linalg.solve(factor.woodbury_cholesky, correction_rhs),
        )
        weighted_low_rank = factor.low_rank / factor.diagonal[:, None]
        return inverse_values - weighted_low_rank @ correction
    if isinstance(factor, GroupedDiagonalLowRankFactor):
        inverse_values = _group_base_solve(factor, values)
        if factor.low_rank.shape[1] == 0:
            return inverse_values
        correction_rhs = factor.low_rank.T @ inverse_values
        correction = np.linalg.solve(
            factor.woodbury_cholesky.T,
            np.linalg.solve(factor.woodbury_cholesky, correction_rhs),
        )
        return inverse_values - factor.base_inverse_low_rank @ correction
    if isinstance(factor, SparseDiagonalLowRankFactor):
        values = np.asarray(values, dtype=float)
        inverse_values = (
            values / factor.diagonal
            if values.ndim == 1
            else values / factor.diagonal[:, None]
        )
        correction = factor.woodbury_factor.solve(factor.low_rank.T @ inverse_values)
        weighted_loading = factor.low_rank.multiply((1.0 / factor.diagonal)[:, None])
        return inverse_values - np.asarray(weighted_loading @ correction)
    if isinstance(factor, SparseCovarianceFactor):
        return np.asarray(
            factor.covariance_factor.solve(np.asarray(values, dtype=float))
        )
    if isinstance(factor, DiagonalSparsePrecisionFactor):
        return _solve_sparse_precision_factor(factor, values)
    if isinstance(factor, NestedLowRankFactor):
        inverse_values = solve_factor(factor.base_factor, values)
        correction_rhs = factor.low_rank.T @ inverse_values
        correction = np.linalg.solve(
            factor.woodbury_cholesky.T,
            np.linalg.solve(factor.woodbury_cholesky, correction_rhs),
        )
        return inverse_values - factor.base_inverse_low_rank @ correction
    factor = np.asarray(factor, dtype=float)
    values = np.asarray(values, dtype=float)
    if factor.ndim == 1:
        inverse_diagonal = 1.0 / np.square(factor)
        if values.ndim == 1:
            return inverse_diagonal * values
        return inverse_diagonal[:, None] * values
    return np.linalg.solve(factor.T, np.linalg.solve(factor, values))


def _solve_sparse_precision_factor(factor, values):
    values = np.asarray(values, dtype=float)
    inverse_values = (
        values / factor.diagonal
        if values.ndim == 1
        else values / factor.diagonal[:, None]
    )
    correction = factor.posterior_factor.solve(factor.loading.T @ inverse_values)
    projected = np.asarray(factor.loading @ correction)
    return inverse_values - (
        projected / factor.diagonal
        if projected.ndim == 1
        else projected / factor.diagonal[:, None]
    )


def factor_logdet(
    factor: np.ndarray
    | DiagonalLowRankFactor
    | GroupedDiagonalLowRankFactor
    | SparseDiagonalLowRankFactor
    | SparseCovarianceFactor
    | DiagonalSparsePrecisionFactor
    | NestedLowRankFactor,
) -> float:
    """Return the covariance log determinant from its factor."""
    if isinstance(factor, DiagonalLowRankFactor):
        return float(np.log(factor.diagonal).sum()) + 2.0 * float(
            np.log(np.diag(factor.woodbury_cholesky)).sum()
        )
    if isinstance(factor, GroupedDiagonalLowRankFactor):
        return (
            float(np.log(factor.diagonal).sum())
            + float(np.log(factor.group_denominator).sum())
            + 2.0 * float(np.log(np.diag(factor.woodbury_cholesky)).sum())
        )
    if isinstance(
        factor,
        (
            SparseDiagonalLowRankFactor,
            SparseCovarianceFactor,
            DiagonalSparsePrecisionFactor,
        ),
    ):
        return factor.logdet
    if isinstance(factor, NestedLowRankFactor):
        return factor_logdet(factor.base_factor) + 2.0 * float(
            np.log(np.diag(factor.woodbury_cholesky)).sum()
        )
    factor = np.asarray(factor, dtype=float)
    diagonal = factor if factor.ndim == 1 else np.diag(factor)
    return 2.0 * float(np.log(diagonal).sum())


def draw_from_factor(
    factor: np.ndarray
    | DiagonalLowRankFactor
    | GroupedDiagonalLowRankFactor
    | SparseDiagonalLowRankFactor
    | SparseCovarianceFactor
    | DiagonalSparsePrecisionFactor
    | NestedLowRankFactor,
    standard_normal: np.ndarray,
    *,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Transform standard-normal draws by a covariance factor."""
    if isinstance(factor, DiagonalLowRankFactor):
        standard_normal = np.asarray(standard_normal, dtype=float)
        if standard_normal.ndim != 1:
            raise ValueError(
                "Structured covariance draws require one standard-normal vector."
            )
        if rng is None:
            raise ValueError("Structured covariance draws require a random generator.")
        diagonal_draw = np.sqrt(factor.diagonal) * standard_normal
        low_rank_draw = factor.low_rank @ rng.standard_normal(factor.low_rank.shape[1])
        return diagonal_draw + low_rank_draw
    if isinstance(factor, GroupedDiagonalLowRankFactor):
        standard_normal = np.asarray(standard_normal, dtype=float)
        if standard_normal.ndim != 1 or rng is None:
            raise ValueError(
                "Structured covariance draws require a vector and random generator."
            )
        diagonal_draw = np.sqrt(factor.diagonal) * standard_normal
        group_draws = rng.standard_normal(len(factor.group_denominator))
        group_draw = factor.group_loading * group_draws[factor.group_index]
        low_rank_draw = factor.low_rank @ rng.standard_normal(factor.low_rank.shape[1])
        return diagonal_draw + group_draw + low_rank_draw
    if isinstance(factor, (SparseDiagonalLowRankFactor, SparseCovarianceFactor)):
        standard_normal = np.asarray(standard_normal, dtype=float)
        if standard_normal.ndim != 1 or rng is None:
            raise ValueError(
                "Structured covariance draws require a vector and random generator."
            )
        return np.sqrt(factor.diagonal) * standard_normal + np.asarray(
            factor.low_rank @ rng.standard_normal(factor.low_rank.shape[1])
        ).reshape(-1)
    if isinstance(factor, DiagonalSparsePrecisionFactor):
        return _draw_sparse_precision_factor(factor, standard_normal, rng)
    if isinstance(factor, NestedLowRankFactor):
        if rng is None:
            raise ValueError("Structured covariance draws require a random generator.")
        return draw_from_factor(
            factor.base_factor, standard_normal, rng=rng
        ) + factor.low_rank @ rng.standard_normal(factor.low_rank.shape[1])
    factor = np.asarray(factor, dtype=float)
    standard_normal = np.asarray(standard_normal, dtype=float)
    if factor.ndim == 1:
        if standard_normal.ndim == 1:
            return factor * standard_normal
        return factor[:, None] * standard_normal
    return factor @ standard_normal


def _draw_sparse_precision_factor(factor, standard_normal, rng):
    standard_normal = np.asarray(standard_normal, dtype=float)
    if standard_normal.ndim != 1 or rng is None:
        raise ValueError(
            "Structured covariance draws require a vector and random generator."
        )
    perturbation = factor.prior_precision_factor.T @ rng.standard_normal(
        factor.prior_precision_factor.shape[0]
    )
    latent = factor.prior_factor.solve(np.asarray(perturbation).reshape(-1))
    return np.sqrt(factor.diagonal) * standard_normal + np.asarray(
        factor.loading @ latent
    ).reshape(-1)


def materialize_covariance(
    covariance: (
        np.ndarray | DiagonalLowRankCovariance | DiagonalSparsePrecisionCovariance
    ),
) -> np.ndarray:
    """Materialize a diagonal-vector or dense covariance representation."""
    if isinstance(covariance, DiagonalLowRankCovariance):
        update = covariance.low_rank @ covariance.low_rank.T
        if sparse.issparse(update):
            update = update.toarray()
        return np.diag(covariance.diagonal) + np.asarray(update)
    if isinstance(covariance, DiagonalSparsePrecisionCovariance):
        return _materialize_sparse_precision_covariance(covariance)
    covariance = np.asarray(covariance, dtype=float)
    return np.diag(covariance) if covariance.ndim == 1 else covariance


def sparse_precision_update_diagonal(loading, precision, *, solver=None) -> np.ndarray:
    """Return ``diag(L @ inv(Q) @ L.T)`` exactly with bounded workspace.

    A stochastic diagonal estimator is not suitable here: these values enter
    event-balanced Gaussian objectives through a sum of logarithms, so even an
    unbiased estimate of each diagonal element would produce a biased
    likelihood.  Blocked solves keep peak memory bounded without changing the
    objective when an analysis crosses an arbitrary observation threshold.
    """
    loading = sparse.csr_matrix(loading, dtype=float)
    solver = (
        splu(sparse.csc_matrix(precision, dtype=float), permc_spec="COLAMD")
        if solver is None
        else solver
    )
    n_rows = loading.shape[0]
    if n_rows == 0:
        return np.empty(0, dtype=float)
    # A 100k-element dense RHS/output block is under 1 MiB in float64.  Sparse
    # direct solvers need additional work arrays, so this conservative block
    # size keeps 5,000-tip exact diagonals well below dense-covariance memory.
    target_entries = 100_000
    block_size = max(
        1,
        min(n_rows, target_entries // max(1, loading.shape[1])),
    )
    diagonal = np.empty(n_rows, dtype=float)
    for start in range(0, n_rows, block_size):
        stop = min(n_rows, start + block_size)
        block = loading[start:stop]
        solved = solver.solve(block.T.toarray())
        diagonal[start:stop] = np.asarray(
            block.multiply(np.asarray(solved).T).sum(axis=1)
        ).reshape(-1)
    scale = max(1.0, float(np.max(np.abs(diagonal), initial=0.0)))
    tolerance = np.finfo(float).eps * scale * max(1, loading.shape[1]) * 100.0
    if np.any(diagonal < -tolerance):
        raise np.linalg.LinAlgError(
            "Sparse precision covariance has a negative marginal variance."
        )
    diagonal[diagonal < 0.0] = 0.0
    return diagonal


def covariance_marginal_diagonal(covariance, *, precision_factor=None) -> np.ndarray:
    """Return marginal variances without materializing structured covariance."""
    if isinstance(covariance, DiagonalLowRankCovariance):
        loading = covariance.low_rank
        row_squares = (
            np.asarray(loading.multiply(loading).sum(axis=1)).reshape(-1)
            if sparse.issparse(loading)
            else np.einsum("ij,ij->i", loading, loading, optimize=True)
        )
        return np.asarray(covariance.diagonal, dtype=float) + row_squares
    if isinstance(covariance, DiagonalSparsePrecisionCovariance):
        if covariance.marginal_diagonal is not None:
            return np.asarray(covariance.marginal_diagonal, dtype=float)
        return np.asarray(
            covariance.diagonal, dtype=float
        ) + sparse_precision_update_diagonal(
            covariance.loading,
            covariance.precision,
            solver=(
                precision_factor.prior_factor
                if isinstance(precision_factor, DiagonalSparsePrecisionFactor)
                else None
            ),
        )
    covariance = np.asarray(covariance, dtype=float)
    return covariance if covariance.ndim == 1 else np.diag(covariance)


def _covariance_observation_count(covariance) -> int:
    if isinstance(
        covariance, (DiagonalLowRankCovariance, DiagonalSparsePrecisionCovariance)
    ):
        return len(covariance.diagonal)
    return len(np.asarray(covariance))


def _normalized_covariance_groups(groups, n_observations):
    groups = np.asarray(groups)
    if (
        groups.ndim != 1
        or len(groups) != n_observations
        or not np.issubdtype(groups.dtype, np.integer)
        or np.any(groups < 0)
    ):
        raise ValueError("Likelihood groups are malformed.")
    groups = groups.astype(int, copy=False)
    n_groups = int(groups.max()) + 1 if len(groups) else 0
    counts = np.bincount(groups, minlength=n_groups)
    if n_groups == 0 or np.any(counts == 0):
        raise ValueError("Likelihood groups must be non-empty and contiguous.")
    return groups, counts


def grouped_mean_covariance_diagonal(
    covariance, groups, *, precision_factor=None
) -> np.ndarray:
    """Return marginal variances of equally weighted within-group means.

    The calculation preserves diagonal-plus-low-rank and sparse-precision
    representations, so event-level pseudo-likelihoods do not need to
    materialize an observation-by-observation covariance matrix.
    """
    n_observations = _covariance_observation_count(covariance)
    groups, counts = _normalized_covariance_groups(groups, n_observations)
    n_groups = len(counts)
    row_weights = 1.0 / counts[groups].astype(float)
    aggregation = sparse.csr_matrix(
        (row_weights, (groups, np.arange(n_observations))),
        shape=(n_groups, n_observations),
    )

    if isinstance(covariance, DiagonalLowRankCovariance):
        diagonal = np.asarray(covariance.diagonal, dtype=float)
        grouped_diagonal = np.bincount(
            groups,
            weights=np.square(row_weights) * diagonal,
            minlength=n_groups,
        )
        grouped_loading = aggregation @ covariance.low_rank
        return grouped_diagonal + np.asarray(
            grouped_loading.multiply(grouped_loading).sum(axis=1)
            if sparse.issparse(grouped_loading)
            else np.einsum("ij,ij->i", grouped_loading, grouped_loading, optimize=True)
        ).reshape(-1)

    if isinstance(covariance, DiagonalSparsePrecisionCovariance):
        if covariance.grouped_marginal_diagonal is not None:
            grouped = np.asarray(covariance.grouped_marginal_diagonal, dtype=float)
            if grouped.shape != (n_groups,):
                raise ValueError("Cached grouped covariance diagonal is malformed.")
            return grouped.copy()
        grouped_diagonal = np.bincount(
            groups,
            weights=np.square(row_weights) * covariance.diagonal,
            minlength=n_groups,
        )
        grouped_loading = sparse.csr_matrix(aggregation @ covariance.loading)
        return grouped_diagonal + sparse_precision_update_diagonal(
            grouped_loading,
            covariance.precision,
            solver=(
                precision_factor.prior_factor
                if isinstance(precision_factor, DiagonalSparsePrecisionFactor)
                else None
            ),
        )

    covariance = np.asarray(covariance, dtype=float)
    if covariance.ndim == 1:
        return np.bincount(
            groups,
            weights=np.square(row_weights) * covariance,
            minlength=n_groups,
        )
    if covariance.shape != (n_observations, n_observations):
        raise ValueError("Covariance has the wrong dimensions.")
    result = np.empty(n_groups, dtype=float)
    for group in range(n_groups):
        selected = np.flatnonzero(groups == group)
        result[group] = float(
            covariance[np.ix_(selected, selected)].sum() / np.square(len(selected))
        )
    return result


def grouped_average_marginal_logdet(
    covariance, groups, *, precision_factor=None
) -> float:
    """Average marginal log variance within each group, removing row scaling."""
    marginal_variances = covariance_marginal_diagonal(
        covariance, precision_factor=precision_factor
    )
    groups, counts = _normalized_covariance_groups(groups, len(marginal_variances))
    if np.any(marginal_variances <= 0.0) or not np.isfinite(marginal_variances).all():
        raise np.linalg.LinAlgError("Covariance marginal variance is not positive.")
    row_weights = 1.0 / counts[groups].astype(float)
    return float(row_weights @ (np.log(marginal_variances) - np.log(counts[groups])))


def _materialize_sparse_precision_covariance(covariance):
    solver = splu(covariance.precision, permc_spec="COLAMD")
    solved = solver.solve(covariance.loading.T.toarray())
    update = np.asarray(covariance.loading @ solved)
    return np.diag(covariance.diagonal) + (update + update.T) / 2.0
