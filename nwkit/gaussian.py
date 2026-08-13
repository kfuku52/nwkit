"""Small linear-algebra helpers for Gaussian covariance calculations."""

import math
from dataclasses import dataclass

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import SuperLU, splu


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


def factor_logdet(
    factor: np.ndarray
    | DiagonalLowRankFactor
    | GroupedDiagonalLowRankFactor
    | SparseDiagonalLowRankFactor
    | SparseCovarianceFactor
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
    if isinstance(factor, (SparseDiagonalLowRankFactor, SparseCovarianceFactor)):
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


def materialize_covariance(
    covariance: np.ndarray | DiagonalLowRankCovariance,
) -> np.ndarray:
    """Materialize a diagonal-vector or dense covariance representation."""
    if isinstance(covariance, DiagonalLowRankCovariance):
        update = covariance.low_rank @ covariance.low_rank.T
        if sparse.issparse(update):
            update = update.toarray()
        return np.diag(covariance.diagonal) + np.asarray(update)
    covariance = np.asarray(covariance, dtype=float)
    return np.diag(covariance) if covariance.ndim == 1 else covariance
