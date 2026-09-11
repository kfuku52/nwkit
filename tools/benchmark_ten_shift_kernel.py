"""Experimental dense likelihood kernel for the noisy/missing 100-tip workload.

Benchmarks a proposal; does not modify NWKIT or benchmark detection outcomes.
"""

import json
import time
from functools import partial
from pathlib import Path

import numpy as np
from benchmark_shift_alpha_models import clean
from benchmark_ten_shifts_missing import fixture
from scipy.linalg import cholesky, solve_triangular

from nwkit.shift_joint_model import evaluate_joint, integral_decay, joint_design
from nwkit.shift_native_model import ShiftLayout


def shared_times(tree):
    n = len(tree.leaf_names)
    shared = np.zeros((n, n))
    for i, (first, last) in enumerate(tree.tip_intervals):
        shared[first:last, first:last] += tree.times[i]
    return shared


def dense_evaluate(data, layout, alpha, covariance, noise, shared):
    # Independent closed-form fixed-root covariance on an ultrametric tree.
    p = len(alpha)
    if np.isinf(alpha).all():
        full = np.einsum("ij,ab->iajb", np.eye(len(shared)), covariance)
    else:
        marginal = integral_decay(2 * alpha, 1.0)
        diffusion = covariance / np.sqrt(marginal[:, None] * marginal[None, :])
        sums = alpha[:, None] + alpha[None, :]
        cov = diffusion[None, None] * integral_decay(
            sums[None, None], shared[:, :, None, None]
        )
        cov *= np.exp(-(1 - shared[:, :, None, None]) * sums[None, None])
        full = cov.transpose(0, 2, 1, 3)
    matrix = full.reshape(len(shared) * p, len(shared) * p).copy()
    matrix.flat[:: len(matrix) + 1] += (data.variances + noise).ravel()
    mask = np.isfinite(data.values)
    matrix = matrix[np.ix_(mask.ravel(), mask.ravel())]
    factor = cholesky(matrix, lower=True, check_finite=False)
    full_design = joint_design(data, layout, alpha)
    design = full_design[mask]
    white = solve_triangular(
        factor,
        np.column_stack((data.values[mask], design)),
        lower=True,
        check_finite=False,
    )
    y, x = white[:, 0], white[:, 1:]
    norms = np.linalg.norm(x, axis=0)
    q, r = np.linalg.qr(x / norms, mode="reduced")
    beta = solve_triangular(r, q.T @ y, check_finite=False) / norms
    residual = y - q @ (q.T @ y)
    inverse = (
        solve_triangular(r, np.eye(len(norms)), check_finite=False) / norms[:, None]
    )
    ll = -0.5 * (
        len(y) * np.log(2 * np.pi)
        + 2 * np.log(np.diag(factor)).sum()
        + residual @ residual
    )
    predicted = np.einsum("ijk,k->ij", full_design, beta)
    return ll, beta, inverse @ inverse.T, predicted


def timed(function):
    for _ in range(3):
        function()
    values = []
    for _ in range(7):
        started = time.perf_counter()
        for _ in range(10):
            function()
        values.append((time.perf_counter() - started) / 10)
    return values


def main():
    rows = []
    for truth in ("shared", "different"):
        data, generating = fixture(
            dict(traits=2, truth=truth, replicate=999, missing_rate=0.2)
        )
        started = time.perf_counter()
        shared = shared_times(data.tree)
        setup = time.perf_counter() - started
        for count in (0, 10):
            layout = ShiftLayout.build(
                data.tree, generating["true_branches"] if count else []
            )
            for alpha in (
                np.ones(2),
                np.array([0.25, 4]),
                np.zeros(2),
                np.full(2, np.inf),
            ):
                covariance = (
                    np.array([[1.0, 0.8], [0.8, 1.0]])
                    / data.scales[:, None]
                    / data.scales[None, :]
                )
                noise = np.full(2, 0.04) / data.scales**2
                prune = partial(evaluate_joint, data, layout, alpha, covariance, noise)
                dense = partial(
                    dense_evaluate, data, layout, alpha, covariance, noise, shared
                )
                baseline = prune()
                experimental = dense()
                errors = dict(
                    predicted=float(
                        np.max(np.abs(baseline.predicted - experimental[3]))
                    ),
                    ll=abs(baseline.log_likelihood - experimental[0]),
                    beta=float(
                        np.max(
                            np.abs(baseline.coefficients.T.ravel() - experimental[1])
                        )
                    ),
                    covariance=float(
                        np.max(
                            np.abs(baseline.coefficient_covariance - experimental[2])
                        )
                    ),
                )
                assert max(errors.values()) < 1e-8, errors
                old, new = timed(prune), timed(dense)
                rows.append(
                    dict(
                        truth=truth,
                        shifts=count,
                        data_sha256=generating["data_sha256"],
                        observed_coordinates=int(np.isfinite(data.values).sum()),
                        covariance_coordinate=covariance,
                        measurement_variance=noise,
                        alpha=alpha.tolist(),
                        errors=errors,
                        pruning_seconds=old,
                        dense_seconds=new,
                        median_ratio=float(np.median(old) / np.median(new)),
                        dense_setup_seconds=setup,
                    )
                )
    Path("/bench/results/kernel.json").write_text(
        json.dumps(clean(rows), indent=2, allow_nan=False) + "\n"
    )
    print(json.dumps(clean(rows), allow_nan=False), flush=True)


if __name__ == "__main__":
    main()
