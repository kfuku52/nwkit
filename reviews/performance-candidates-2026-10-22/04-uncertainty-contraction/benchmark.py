"""Prototype grouped predictor covariance projection without production edits."""

import argparse
import fcntl
import inspect
import json
import os
import platform
import statistics
import sys
import time
import tracemalloc
from contextlib import contextmanager
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
LOCK = Path("/tmp/nwkit-perf-20261022-9fYwTD/measurement.lock")
sys.path.insert(0, str(ROOT))

LOCAL_LOOP = """            for observation in range(n_observations):
                block = slopes.T @ covariance[observation] @ slopes
                start = observation * dimension
                projected[start : start + dimension, start : start + dimension] += block
"""
LOCAL_BATCH = """            blocks = slopes.T @ covariance @ slopes
            rows = np.arange(n_observations)
            view = projected.reshape(n_observations, dimension, n_observations, dimension)
            view[rows, :, rows, :] += blocks
"""
JOINT_LOOP = """        for first in range(n_observations):
            first_start = first * dimension
            for second in range(n_observations):
                second_start = second * dimension
                block = slopes.T @ covariance[:, :, first, second] @ slopes
                projected[
                    first_start : first_start + dimension,
                    second_start : second_start + dimension,
                ] += block
"""
JOINT_BATCH = """        blocks = slopes.T @ covariance.transpose(2, 3, 0, 1) @ slopes
        projected.reshape(n_observations, dimension, n_observations, dimension)[:] += (
            blocks.transpose(0, 2, 1, 3)
        )
"""
JOINT_CHUNKS = """        view = projected.reshape(n_observations, dimension, n_observations, dimension)
        for start in range(0, n_observations, 64):
            end = min(start + 64, n_observations)
            blocks = slopes.T @ covariance[:, :, start:end, :].transpose(2, 3, 0, 1) @ slopes
            view[start:end] += blocks.transpose(0, 2, 1, 3)
"""


def build_candidate(batch_all=False):
    import nwkit.phylogenetic_glmm as glmm

    source = inspect.getsource(glmm._project_predictor_uncertainty)
    assert source.count(LOCAL_LOOP) == 1
    assert source.count(JOINT_LOOP) == 1
    source = source.replace(LOCAL_LOOP, LOCAL_BATCH)
    source = source.replace(JOINT_LOOP, JOINT_BATCH if batch_all else JOINT_CHUNKS)
    namespace = {}
    exec(compile(source, str(HERE / "benchmark.py"), "exec"), glmm.__dict__, namespace)
    return namespace["_project_predictor_uncertainty"]


@contextmanager
def prototype():
    import nwkit.phylogenetic_glmm as glmm

    original = glmm._project_predictor_uncertainty
    glmm._project_predictor_uncertainty = build_candidate()
    try:
        yield
    finally:
        glmm._project_predictor_uncertainty = original


def stats(samples):
    return {
        "samples_seconds_per_call": samples,
        "median_seconds": statistics.median(samples),
        "min_seconds": min(samples),
        "max_seconds": max(samples),
        "relative_range": (max(samples) - min(samples)) / statistics.median(samples),
    }


def timed(functions, count=None):
    # Batch very short helpers to amortize timer/scheduler noise. Count includes
    # every call, and all variants use the same inputs and number of invocations.
    warm = {}
    for name, function in functions.items():
        start = time.perf_counter()
        function()
        warm[name] = time.perf_counter() - start
    if count is None:
        count = max(1, min(1000, int(0.035 / min(warm.values())) + 1))
    samples = {name: [] for name in functions}
    for repetition in range(7):
        names = list(functions)
        if repetition % 2:
            names.reverse()
        for name in names:
            start = time.perf_counter()
            for _ in range(count):
                functions[name]()
            samples[name].append((time.perf_counter() - start) / count)
    output = {name: stats(values) for name, values in samples.items()}
    for record in output.values():
        record["speedup"] = (
            output["baseline"]["median_seconds"] / record["median_seconds"]
        )
    return {"batch_calls": count, "timings": output}


def peak(function):
    tracemalloc.start()
    try:
        function()
        return tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()


def fixture(n, k, dimension, kind, scale=1.0):
    import numpy as np

    rng = np.random.default_rng(481 + n + k + dimension)
    coefficients = rng.normal(size=(k + 1, dimension))
    if kind == "local":
        loading = rng.normal(size=(n, k, k)) * 0.03
        covariance = loading @ loading.transpose(0, 2, 1)
    elif kind == "joint":
        loading = rng.normal(size=(k, n, min(8, n * k))) * 0.03
        covariance = np.einsum("anr,bmr->abnm", loading, loading)
    else:
        raise ValueError(kind)
    return coefficients, [covariance * scale], [tuple(range(1, k + 1))], n


def equivalence(first, second):
    import numpy as np

    np.testing.assert_allclose(first, second, rtol=1e-12, atol=1e-13, equal_nan=True)
    finite = np.isfinite(first) & np.isfinite(second)
    delta = np.abs(first[finite] - second[finite])
    scale = np.maximum(np.maximum(np.abs(first[finite]), np.abs(second[finite])), 1.0)
    return {
        "maximum_absolute_difference": float(delta.max(initial=0.0)),
        "maximum_scaled_difference": float((delta / scale).max(initial=0.0)),
        "bitwise_equal": bool(np.array_equal(first, second, equal_nan=True)),
    }


def helper_case(n, k, dimension, kind):
    import nwkit.phylogenetic_glmm as glmm

    arguments = fixture(n, k, dimension, kind)
    baseline = glmm._project_predictor_uncertainty
    functions = {
        "baseline": lambda: baseline(*arguments),
        "batched": lambda: build_all(*arguments),
        "chunked": lambda: build_chunks(*arguments),
    }
    build_all, build_chunks = build_candidate(batch_all=True), build_candidate()
    expected = functions["baseline"]()
    comparisons = {
        name: equivalence(expected, function()) for name, function in functions.items()
    }
    result = timed(functions)
    result.update(
        observations=n,
        predictor_columns=k,
        response_linear_predictors=dimension,
        kind=kind,
        output_shape=list(expected.shape),
        equivalence=comparisons,
        peak_python_bytes={
            name: peak(function) for name, function in functions.items()
        },
    )
    return result


def edge_checks():
    import numpy as np
    from scipy import sparse

    import nwkit.phylogenetic_glmm as glmm
    from nwkit.sparse_laplace import (
        ContinuousPredictorUncertainty,
        GroupedPredictorUncertainty,
        JointPredictorUncertainty,
    )

    baseline, candidate = glmm._project_predictor_uncertainty, build_candidate()
    result = []
    for kind in ("local", "joint"):
        for n, k, dimension in ((1, 2, 1), (7, 3, 2), (65, 2, 3)):
            for scale in (0.0, 1e-200, 1e200):
                args = fixture(n, k, dimension, kind, scale)
                result.append(
                    {
                        "case": [kind, n, k, dimension, scale],
                        **equivalence(baseline(*args), candidate(*args)),
                    }
                )
    # Non-contiguous grouped covariance, multiple groups, repeated selected
    # columns, mixed continuous and structured cases retain original dispatch.
    coefficients, terms, columns, n = fixture(9, 3, 2, "joint")
    cases = {
        "reversed_noncontiguous": (
            coefficients,
            [terms[0][:, :, ::-1, ::-1]],
            columns,
            n,
        ),
        "duplicate_selected_columns": (coefficients, terms, [(2, 2, 1)], n),
        "multiple_terms": (coefficients, [terms[0], terms[0] * 0.5], columns * 2, n),
        "continuous": (coefficients, [np.eye(n) * 0.01], [1], n),
        "continuous_structured": (
            coefficients,
            [ContinuousPredictorUncertainty(np.eye(n) * 0.1, np.arange(n))],
            [1],
            n,
        ),
        "empty": (coefficients, [], [], n),
        "joint_structured": (
            coefficients[:, :1],
            [
                JointPredictorUncertainty(
                    factors=tuple(sparse.eye(n, format="csr") for _ in range(3))
                )
            ],
            [(1, 2, 3)],
            n,
        ),
        "local_structured": (
            coefficients,
            [
                GroupedPredictorUncertainty(
                    tuple(np.eye(3) * 0.1 for _ in range(n)), np.arange(n)
                )
            ],
            [(1, 2, 3)],
            n,
        ),
    }
    for label, args in cases.items():
        result.append({"case": label, **equivalence(baseline(*args), candidate(*args))})
    # An independently assembled linear operator verifies the observation-major
    # block order, instead of comparing only two implementations of the loop.
    for kind in ("local", "joint"):
        from scipy.linalg import block_diag

        coefficients, terms, columns, n = fixture(7, 3, 2, kind)
        slopes = coefficients[np.asarray(columns[0])]
        operator = np.kron(np.eye(n), slopes.T)
        if kind == "local":
            covariance = block_diag(*terms[0])
        else:
            covariance = terms[0].transpose(2, 0, 3, 1).reshape(n * 3, n * 3)
        oracle = operator @ covariance @ operator.T
        oracle = (oracle + oracle.T) / 2.0
        result.append(
            {
                "case": f"independent_dense_operator_{kind}",
                **equivalence(oracle, candidate(coefficients, terms, columns, n)),
            }
        )
    invalid = {
        "term_count": (coefficients, terms, [], n),
        "wrong_group_shape": (coefficients, [np.zeros((3, 3, n))], [(1, 2, 3)], n),
        "wrong_continuous_shape": (coefficients, [np.zeros((3, 3))], [1], n),
    }
    for label, args in invalid.items():
        errors = []
        for function in (baseline, candidate):
            try:
                function(*args)
            except Exception as error:
                errors.append([type(error).__name__, str(error)])
        assert len(errors) == 2 and errors[0] == errors[1]
        result.append({"case": label, "matching_errors": errors[0]})
    # The function's baseline does not validate PSD/finite covariance itself;
    # nonfinite blocks must propagate identically rather than become zeros.
    args = fixture(3, 2, 2, "joint")
    args[1][0][0, 0, 1, 2] = np.nan
    result.append(
        {"case": "nan_propagation", **equivalence(baseline(*args), candidate(*args))}
    )
    return result


def fit_case():
    import numpy as np

    import nwkit.phylogenetic_glmm as glmm

    n = 12
    rng = np.random.default_rng(749)
    design = np.column_stack([np.ones(n), rng.normal(size=(n, 2))])
    response = rng.poisson(np.exp(0.6 + design[:, 1] * 0.3)).astype(float)
    loading = rng.normal(size=(2, n, 4)) * 0.03
    uncertainty = np.einsum("anr,bmr->abnm", loading, loading)
    covariance = 0.3 * np.ones((n, n)) + 0.7 * np.eye(n)

    def fit():
        return glmm.fit_phylogenetic_glmm(
            response,
            design,
            lambda _: covariance,
            family="poisson",
            predictor_uncertainties=[uncertainty],
            predictor_columns=[(1, 2)],
        )

    def fit_candidate():
        with prototype():
            return fit()

    baseline_fit, candidate_fit = fit(), fit_candidate()
    comparisons = {}
    for name in baseline_fit.__dataclass_fields__:
        first, second = getattr(baseline_fit, name), getattr(candidate_fit, name)
        if isinstance(first, np.ndarray):
            np.testing.assert_allclose(
                first, second, rtol=2e-4, atol=2e-6, equal_nan=True
            )
            comparisons[name] = float(np.max(np.abs(first - second), initial=0.0))
        elif isinstance(first, dict):
            assert first.keys() == second.keys()
            for key in first:
                np.testing.assert_allclose(
                    first[key], second[key], rtol=2e-4, atol=2e-6
                )
        elif isinstance(first, float):
            if np.isnan(first):
                assert np.isnan(second)
            else:
                np.testing.assert_allclose(first, second, rtol=2e-4, atol=2e-6)
                comparisons[name] = abs(first - second)
        else:
            assert first == second, (name, first, second)
    return {
        "observations": n,
        "family": "poisson",
        "predictor_columns": 2,
        "equivalence_max_abs": comparisons,
        "baseline_log_likelihood": baseline_fit.log_likelihood,
        **timed({"baseline": fit, "chunked": fit_candidate}, count=1),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fit-only", action="store_true")
    parser.add_argument("--checks-only", action="store_true")
    args = parser.parse_args()
    LOCK.parent.mkdir(exist_ok=True, parents=True)
    with LOCK.open("w") as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        import numpy as np
        import scipy

        from nwkit import __version__

        record = {
            "baseline_commit": "d7116f80b9979739c69a5431a5e34b822795408b",
            "environment": {
                "platform": platform.platform(),
                "python": platform.python_version(),
                "numpy": np.__version__,
                "scipy": scipy.__version__,
                "nwkit": __version__,
                "threads": {
                    name: os.getenv(name)
                    for name in (
                        "OPENBLAS_NUM_THREADS",
                        "OMP_NUM_THREADS",
                        "VECLIB_MAXIMUM_THREADS",
                        "MKL_NUM_THREADS",
                    )
                },
                "host_contention": "Unrelated user CPU jobs were active; raw spreads retained.",
            },
        }
        if args.checks_only:
            record["edge_checks"] = edge_checks()
            output = HERE / "additional-checks.json"
            output.write_text(json.dumps(record, indent=2) + "\n")
            print(json.dumps({"output": str(output)}))
            return
        if not args.fit_only:
            record["edge_checks"] = edge_checks()
            record["helper_cases"] = [
                helper_case(*case)
                for case in (
                    (16, 2, 1, "joint"),
                    (64, 2, 2, "joint"),
                    (128, 4, 3, "joint"),
                    (256, 3, 2, "joint"),
                    (64, 2, 2, "local"),
                    (256, 4, 3, "local"),
                )
            ]
        record["fit_case"] = fit_case()
        output = HERE / ("fit-only-results.json" if args.fit_only else "results.json")
        output.write_text(json.dumps(record, indent=2) + "\n")
        print(
            json.dumps(
                {"output": str(output), "cases": len(record.get("helper_cases", []))}
            )
        )


if __name__ == "__main__":
    main()
