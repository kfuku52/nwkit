"""Numerical error contracts shared by dense and sparse Gaussian factors."""

import numpy as np
import pytest
from scipy import sparse

from nwkit import gaussian


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
