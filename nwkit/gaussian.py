"""Small linear-algebra helpers for Gaussian covariance calculations."""

from dataclasses import dataclass

import numpy as np


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
    low_rank: np.ndarray


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
    diagonal = np.asarray(diagonal, dtype=float)
    low_rank = np.asarray(low_rank, dtype=float)
    if diagonal.ndim != 1 or np.any(diagonal <= 0.0):
        raise np.linalg.LinAlgError("Covariance diagonal is not positive.")
    if low_rank.ndim != 2 or low_rank.shape[0] != len(diagonal):
        raise ValueError("Low-rank covariance factor has the wrong dimensions.")
    weighted = low_rank / diagonal[:, None]
    woodbury = np.eye(low_rank.shape[1]) + low_rank.T @ weighted
    return DiagonalLowRankFactor(
        diagonal=diagonal,
        low_rank=low_rank,
        woodbury_cholesky=np.linalg.cholesky(woodbury),
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


def factor_diagonal_low_rank_updates(
    diagonal: np.ndarray, updates: list[np.ndarray]
) -> DiagonalLowRankFactor | GroupedDiagonalLowRankFactor:
    """Factor diagonal covariance plus low-rank updates, exploiting groups."""
    diagonal = np.asarray(diagonal, dtype=float)
    if diagonal.ndim != 1 or np.any(diagonal <= 0.0):
        raise np.linalg.LinAlgError("Covariance diagonal is not positive.")
    normalized_updates = [np.asarray(update, dtype=float) for update in updates]
    if any(
        update.ndim != 2
        or update.shape[0] != len(diagonal)
        or not np.isfinite(update).all()
        for update in normalized_updates
    ):
        raise ValueError("Low-rank covariance update has invalid dimensions or values.")
    grouped_update = _select_grouped_update(normalized_updates)
    if grouped_update is None:
        return factor_diagonal_low_rank(diagonal, np.column_stack(normalized_updates))
    grouped_index, group_index, group_loading = grouped_update
    low_rank_updates = [
        update
        for index, update in enumerate(normalized_updates)
        if index != grouped_index
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
        minlength=normalized_updates[grouped_index].shape[1],
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


def solve_factor(
    factor: np.ndarray | DiagonalLowRankFactor | GroupedDiagonalLowRankFactor,
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
    factor = np.asarray(factor, dtype=float)
    values = np.asarray(values, dtype=float)
    if factor.ndim == 1:
        inverse_diagonal = 1.0 / np.square(factor)
        if values.ndim == 1:
            return inverse_diagonal * values
        return inverse_diagonal[:, None] * values
    return np.linalg.solve(factor.T, np.linalg.solve(factor, values))


def factor_logdet(
    factor: np.ndarray | DiagonalLowRankFactor | GroupedDiagonalLowRankFactor,
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
    factor = np.asarray(factor, dtype=float)
    diagonal = factor if factor.ndim == 1 else np.diag(factor)
    return 2.0 * float(np.log(diagonal).sum())


def draw_from_factor(
    factor: np.ndarray | DiagonalLowRankFactor | GroupedDiagonalLowRankFactor,
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
        return (
            np.diag(covariance.diagonal) + covariance.low_rank @ covariance.low_rank.T
        )
    covariance = np.asarray(covariance, dtype=float)
    return np.diag(covariance) if covariance.ndim == 1 else covariance
