"""Measure local factor reuse in the complete K permutation operation.

No production module is modified. Set all BLAS thread variables to one before
running. The main function serializes work with the shared measurement lock.
"""

import fcntl
import json
import os
import platform
import statistics
import sys
import time
import tracemalloc
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import cho_solve

from nwkit import __version__
from nwkit.signal_stats import k_statistic
from nwkit.signal_stats import permutation_test as baseline_permutation_test

DIRECTORY = Path(__file__).resolve().parent
LOCK = Path("/tmp/nwkit-perf-20261022-9fYwTD/measurement.lock")


def prepared_k(covariance, count):
    factor = np.linalg.cholesky(covariance)
    weights = cho_solve((factor, True), np.ones(count))
    precision = float(weights.sum())
    expected = (np.trace(covariance) - count / precision) / (count - 1)

    def statistic(values):
        mean = float(weights @ values / precision)
        residual = values - mean
        quadratic = float(residual @ cho_solve((factor, True), residual))
        if quadratic <= 0 or expected <= 0:
            raise ValueError(
                "K is undefined for constant traits or degenerate covariance."
            )
        return float(np.sum((values - mean) ** 2) / quadratic / expected)

    return statistic


def candidate_permutation_test(covariance, values, errors, observed, simulations, rng):
    # The zero-draw API returns without examining covariance in the baseline.
    # No cached state survives a call. Error-aware permutations must refit.
    if simulations <= 0 or np.any(errors):
        return baseline_permutation_test(
            covariance, values, errors, observed, simulations, rng
        )
    statistic = prepared_k(covariance, len(values))
    exceed = 0
    for _ in range(simulations):
        order = rng.permutation(len(values))
        value = statistic(values[order])
        exceed += value >= observed - 1e-12 * max(1.0, abs(observed))
    return (1 + exceed) / (1 + simulations)


def covariance_fixture(count, shape):
    if shape == "star":
        return np.eye(count)
    if shape == "ladder":
        positions = np.arange(count)
        return np.minimum(positions[:, None], positions[None, :]) + np.eye(count) + 1
    covariance = np.eye(count)
    block = 2
    while block < count:
        for start in range(0, count, block):
            covariance[start : start + block, start : start + block] += 1
        block *= 2
    return covariance


def summarize(times):
    return {
        "samples_seconds": times,
        "median_seconds": statistics.median(times),
        "min_seconds": min(times),
        "max_seconds": max(times),
        "relative_range": (max(times) - min(times)) / statistics.median(times),
    }


def peak_allocation(function):
    tracemalloc.start()
    try:
        function()
        return tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()


def benchmark_case(count, simulations, shape):
    covariance = covariance_fixture(count, shape)
    values = np.random.default_rng(10 + count).normal(size=count)
    errors = np.zeros(count)
    observed = k_statistic(covariance, values)
    arguments = covariance, values, errors, observed, simulations
    functions = {
        "baseline": baseline_permutation_test,
        "candidate": candidate_permutation_test,
    }

    def call(name):
        rng = np.random.default_rng(974)
        result = functions[name](*arguments, rng)
        return result, rng.bit_generator.state

    expected, expected_state = call("baseline")
    actual, actual_state = call("candidate")
    assert actual == expected
    assert actual_state == expected_state
    # Warm up each complete operation, including local factor preparation.
    for name in functions:
        call(name)
    samples = {name: [] for name in functions}
    for repetition in range(7):
        order = (
            ["baseline", "candidate"]
            if repetition % 2 == 0
            else ["candidate", "baseline"]
        )
        for name in order:
            start = time.perf_counter()
            value, state = call(name)
            samples[name].append(time.perf_counter() - start)
            assert value == expected and state == expected_state
    timings = {name: summarize(times) for name, times in samples.items()}
    prepared = prepared_k(covariance, count)
    rng = np.random.default_rng(105)
    maximum_difference = 0.0
    for _ in range(100):
        permuted = values[rng.permutation(count)]
        maximum_difference = max(
            maximum_difference,
            abs(k_statistic(covariance, permuted) - prepared(permuted)),
        )
    return {
        "tips": count,
        "simulations": simulations,
        "shape": shape,
        "p_value": expected,
        "p_value_exact_match": True,
        "rng_state_exact_match": True,
        "maximum_k_absolute_difference_100_permutations": maximum_difference,
        "timings": timings,
        "median_speedup": timings["baseline"]["median_seconds"]
        / timings["candidate"]["median_seconds"],
        "peak_python_bytes": {
            name: peak_allocation(lambda name=name: call(name)) for name in functions
        },
    }


def edge_checks():
    checks = {}
    covariance = covariance_fixture(8, "balanced")
    values = np.array([1, 1.5, 2, 1.8, 4, 4.2, 3.8, 4.5])
    for label, errors in (
        ("positive_errors", np.array([0.1, 0.2] * 4)),
        ("mixed_errors", np.array([0.0, 0.2] * 4)),
        ("zero_errors", np.zeros(8)),
    ):
        from nwkit.signal_stats import k_fit

        observed = k_fit(covariance, values, errors)[0]
        results, states = [], []
        for function in (baseline_permutation_test, candidate_permutation_test):
            rng = np.random.default_rng(47)
            results.append(function(covariance, values, errors, observed, 9, rng))
            states.append(rng.bit_generator.state)
        assert results[0] == results[1] and states[0] == states[1]
        checks[label] = {
            "p_value": results[0],
            "exact_match": True,
            "rng_exact_match": True,
        }
    # Keep the baseline's zero-draw behavior even for an invalid covariance.
    for function in (baseline_permutation_test, candidate_permutation_test):
        assert function(None, values, np.zeros(8), 1, 0, np.random.default_rng(1)) == 1
    checks["zero_draw_invalid_covariance"] = True
    failure_checks = {}
    for label, c, y in (
        ("singular_covariance", np.zeros((8, 8)), values),
        ("constant_trait", covariance, np.ones(8)),
    ):
        failures = []
        for function in (baseline_permutation_test, candidate_permutation_test):
            try:
                function(c, y, np.zeros(8), 1, 2, np.random.default_rng(1))
            except (ValueError, np.linalg.LinAlgError) as exc:
                failures.append((type(exc).__name__, str(exc)))
        assert len(failures) == 2 and failures[0] == failures[1]
        failure_checks[label] = {"exception": failures[0], "exact_match": True}
    checks["failures"] = failure_checks
    # Deliberately put the observed threshold on each of several permuted K's.
    # Reusing the very same factor and operations preserves all finite ties here.
    prepared = prepared_k(covariance, 8)
    rng = np.random.default_rng(3)
    differences = []
    for _ in range(100):
        permuted = values[rng.permutation(8)]
        differences.append(abs(k_statistic(covariance, permuted) - prepared(permuted)))
    assert max(differences) == 0
    checks["ties_max_k_difference"] = max(differences)
    # Different covariance on every invocation proves there is no stale global cache.
    checks["local_cache_only"] = "factor is held in a closure owned by one invocation"
    return checks


def main():
    with LOCK.open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        result = {
            "candidate": "01-signal-factor-reuse",
            "baseline_commit": "d7116f80b9979739c69a5431a5e34b822795408b",
            "environment": {
                "python": sys.version,
                "platform": platform.platform(),
                "numpy": np.__version__,
                "scipy": scipy.__version__,
                "nwkit": __version__,
                "thread_environment": {
                    name: os.environ.get(name)
                    for name in (
                        "OPENBLAS_NUM_THREADS",
                        "OMP_NUM_THREADS",
                        "VECLIB_MAXIMUM_THREADS",
                        "MKL_NUM_THREADS",
                    )
                },
            },
            "scope": "Complete permutation_test operation; includes candidate factor setup, excludes CLI/tree I/O and observed K calculation common to both",
            "edge_checks": edge_checks(),
            "cases": [],
        }
        for arguments in (
            (64, 999, "balanced"),
            (256, 999, "balanced"),
            (512, 199, "balanced"),
            (128, 199, "ladder"),
            (128, 199, "star"),
        ):
            case = benchmark_case(*arguments)
            result["cases"].append(case)
            print(
                json.dumps(
                    {
                        key: case[key]
                        for key in ("tips", "shape", "median_speedup", "p_value")
                    }
                ),
                flush=True,
            )
        (DIRECTORY / "results.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
