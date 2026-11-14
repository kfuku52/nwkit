"""Numerical error contracts shared by dense and sparse Gaussian factors."""

import math

import numpy as np
import pytest
from scipy import sparse

from nwkit import gaussian


@pytest.mark.parametrize("groups", [[2, 0, 1, 4, 3], [1, 0, 1, 2, 0]])
def test_grouped_dense_covariance_uses_the_correct_rows_in_group_order(groups):
    loading = np.array([[1, 1], [2, -1], [3, 0], [-1, 2], [0, -1]], dtype=float)
    covariance = loading @ loading.T
    membership = np.asarray(groups) == np.arange(max(groups) + 1)[:, None]
    weights = membership / membership.sum(axis=1, keepdims=True)
    expected = np.diag(weights @ covariance @ weights.T)
    np.testing.assert_allclose(
        gaussian.grouped_mean_covariance_diagonal(covariance, groups), expected
    )


def test_grouped_dense_covariance_keeps_large_finite_mean_variance():
    covariance = np.array([[1.6e308, 8e307], [8e307, 1.2e308]])
    with np.errstate(over="raise", invalid="raise"):
        observed = gaussian.grouped_mean_covariance_diagonal(covariance, [0, 0])
    np.testing.assert_allclose(observed, [1.1e308], rtol=1e-15)


@pytest.mark.parametrize("columns", [None, 2])
def test_diagonal_cholesky_solve_keeps_finite_results_at_subnormal_variance(columns):
    diagonal = np.ldexp(np.array([1.0, 2.0, 4.0]), -1060)
    values = np.array([-2.0, 0.0, 3.0])
    rhs = diagonal * values
    if columns:
        rhs = np.column_stack([rhs, -rhs])
        values = np.column_stack([values, -values])
    factor = gaussian.factor_covariance(np.diag(diagonal))
    with np.errstate(over="raise", invalid="raise"):
        solution = gaussian.solve_factor(factor, rhs)
    np.testing.assert_allclose(solution, values, rtol=1e-15, atol=0.0)


@pytest.mark.parametrize("kind", ["dense", "low-rank", "grouped", "nested"])
@pytest.mark.parametrize("vector", [True, False])
def test_cholesky_solves_match_dense_oracle_without_refactorization(
    kind, vector, monkeypatch
):
    rng = np.random.default_rng(241)
    diagonal = np.linspace(0.5, 2.0, 8)
    loading = rng.normal(size=(8, 3))
    covariance = np.diag(diagonal) + loading @ loading.T
    if kind == "dense":
        factor = gaussian.factor_covariance(covariance)
    elif kind == "low-rank":
        factor = gaussian.factor_diagonal_low_rank(diagonal, loading)
    elif kind == "grouped":
        groups = np.eye(2)[np.arange(8) % 2]
        covariance += groups @ groups.T
        factor = gaussian.factor_diagonal_low_rank_updates(diagonal, [groups, loading])
        assert isinstance(factor, gaussian.GroupedDiagonalLowRankFactor)
    else:
        base = gaussian.factor_diagonal_low_rank(diagonal, loading[:, :1])
        factor = gaussian._factor_nested_low_rank(
            base, sparse.csr_matrix(loading[:, 1:])
        )
    rhs = rng.normal(size=8 if vector else (8, 3))
    expected = np.linalg.solve(covariance, rhs)

    def forbidden(*args, **kwargs):
        raise AssertionError("refactorized a Cholesky factor")

    monkeypatch.setattr(np.linalg, "solve", forbidden)
    np.testing.assert_allclose(
        gaussian.solve_factor(factor, rhs), expected, rtol=1e-12, atol=1e-12
    )


@pytest.mark.parametrize(
    "columns, force_observation_space, label",
    [(1, True, "covariance"), (3, False, "covariance"), (2, False, "Woodbury")],
)
def test_sparse_singular_candidates_raise_linear_algebra_errors(
    columns, force_observation_space, label
):
    # Although the abstract covariance is positive definite, adding this tiny
    # diagonal to a rank-one update loses it in floating-point arithmetic.
    with pytest.raises(
        np.linalg.LinAlgError, match=f"Sparse {label} matrix is singular"
    ):
        gaussian.factor_diagonal_sparse_low_rank(
            np.full(3, 1e-20),
            sparse.csr_matrix(np.ones((3, columns))),
            force_observation_space=force_observation_space,
        )


def test_sparse_singular_prior_precision_uses_the_same_error_contract():
    with pytest.raises(
        np.linalg.LinAlgError,
        match="Sparse prior precision is (singular|not positive definite)",
    ):
        gaussian.factor_diagonal_sparse_precision_updates(
            np.ones(3),
            [
                (
                    sparse.csr_matrix(np.ones((3, 1))),
                    sparse.csc_matrix((1, 1)),
                    sparse.csr_matrix((1, 1)),
                )
            ],
        )


@pytest.mark.parametrize("precision_update", [False, True])
def test_sparse_covariance_solve_and_determinant_at_subnormal_variances(
    precision_update,
):
    exponent = -1060
    relative_diagonal = np.array([1.0, 2.0, 4.0, 8.0])
    relative_loading = np.array([0.25, -0.5, 0.0, 0.5])
    diagonal = np.ldexp(relative_diagonal, exponent)
    loading = sparse.csr_matrix(np.ldexp(relative_loading[:, None], exponent // 2))
    precision = sparse.eye(1, format="csc")
    base_covariance = np.diag(relative_diagonal) + np.outer(
        relative_loading, relative_loading
    )
    expected = np.linspace(-1.0, 2.0, len(diagonal))
    base_rhs = np.stack([base_covariance @ expected, relative_diagonal], axis=1)
    rhs = np.ldexp(base_rhs, exponent)
    independent_solution = np.linalg.solve(base_covariance, np.ldexp(rhs, -exponent))
    expected_logdet = np.linalg.slogdet(base_covariance)[1] + len(
        diagonal
    ) * exponent * math.log(2.0)
    with np.errstate(over="raise", invalid="raise", divide="raise"):
        if precision_update:
            factor = gaussian.factor_diagonal_sparse_precision_updates(
                diagonal, [(loading, precision, precision.tocsr())]
            )
        else:
            factor = gaussian.factor_diagonal_sparse_low_rank(diagonal, loading)
        observed = gaussian.solve_factor(factor, rhs)
    np.testing.assert_allclose(observed, independent_solution, rtol=2e-14, atol=2e-14)
    assert gaussian.factor_logdet(factor) == pytest.approx(
        expected_logdet, rel=0.0, abs=1e-12
    )


def test_sparse_roundoff_cannot_turn_rank_deficiency_into_an_indefinite_factor():
    loading = sparse.csr_matrix(
        [[-20.749325931793987], [167.02490282686168], [165.85872694090048]]
    )
    with pytest.raises(np.linalg.LinAlgError, match="singular|positive definite"):
        gaussian.factor_diagonal_sparse_low_rank(
            np.full(3, 1e-16), loading, force_observation_space=True
        )


def test_sparse_precision_helpers_reject_indefinite_matrices_before_solving():
    precision = sparse.csc_matrix([[1.0, 2.0], [2.0, 1.0]])
    with pytest.raises(np.linalg.LinAlgError, match="singular|positive definite"):
        gaussian.sparse_precision_update_diagonal(sparse.eye(2), precision)


@pytest.mark.parametrize("scale", [2.0**-900, 2.0**1023])
def test_sparse_precision_validation_and_solve_use_relative_units(scale):
    reference = np.array([[1.0, 0.25], [0.25, 1.0]])
    matrix = sparse.csc_matrix(scale * reference)
    with np.errstate(over="raise", invalid="raise"):
        factor = gaussian._factor_sparse_lu(matrix, "Sparse precision")
        np.testing.assert_allclose(factor.solve(scale * reference[:, 1]), [0.0, 1.0])
        asymmetric = sparse.csc_matrix(scale * np.array([[1.0, 0.5], [0.25, 1.0]]))
        with pytest.raises(np.linalg.LinAlgError, match="not symmetric"):
            gaussian._factor_sparse_lu(asymmetric, "Sparse precision")


def test_sparse_precision_diagonal_reuses_states_with_distinct_observation_loadings():
    precision = np.array([[2.0, -0.3], [-0.3, 1.5]])
    loading = np.array([[0.0, 2.0], [0.5, 0.0], [0.0, 0.0], [0.0, -1.5]])
    expected = np.diag(loading @ np.linalg.solve(precision, loading.T))
    observed = gaussian.sparse_precision_update_diagonal(
        sparse.csr_matrix(loading), sparse.csc_matrix(precision)
    )
    np.testing.assert_allclose(observed, expected)


@pytest.mark.parametrize("blocks", [1, 2])
def test_sparse_inertia_check_rejects_zero_diagonal_pivoting(blocks):
    # This positive-diagonal matrix develops a zero Schur pivot. SuperLU can
    # row-swap it into an all-positive U diagonal unless we certify that the
    # symmetric row and column permutations remained identical. Two blocks
    # also cover an even number of negative eigenvalues.
    zero_pivot = sparse.csc_matrix(
        [
            [1.0, 1.0, 0.0, 0.0],
            [1.0, 1.0, 1.0, 0.0],
            [0.0, 1.0, 1.0, 1.0],
            [0.0, 0.0, 1.0, 1.0],
        ]
    )
    indefinite = sparse.kron(sparse.eye(blocks), zero_pivot, format="csc")
    with pytest.raises(np.linalg.LinAlgError, match="positive definite"):
        gaussian._factor_sparse_lu(indefinite, "Sparse probe matrix")


def test_sparse_spd_factor_retains_the_symmetric_ldlt_structure():
    loading = sparse.csr_matrix([[1.0, 0.2], [0.0, -0.7], [0.5, 0.0], [0.3, 1.1]])
    matrix = sparse.eye(4, format="csc") + loading @ loading.T
    solver = gaussian._factor_sparse_lu(matrix, "Sparse probe matrix")
    assert np.array_equal(solver.perm_r, solver.perm_c)
    diagonal = sparse.diags(solver.U.diagonal(), format="csc")
    difference = (solver.U.T - solver.L @ diagonal).tocoo()
    assert np.max(np.abs(difference.data), initial=0.0) < 1e-12


def test_unexpected_sparse_backend_errors_are_not_reclassified(monkeypatch):
    error = RuntimeError("unexpected backend failure")

    def fail(*args, **kwargs):
        raise error

    monkeypatch.setattr(gaussian, "splu", fail)
    with pytest.raises(RuntimeError) as caught:
        gaussian.factor_diagonal_sparse_low_rank(
            np.ones(3), sparse.eye(3, format="csr")
        )
    assert caught.value is error
