"""Compare exact diagonal detectors, preserving the production implementation."""

import fcntl
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
sys.path.insert(0, str(ROOT))
LOCK = Path("/tmp/nwkit-perf-20261022-9fYwTD/measurement.lock")

import numpy as np  # noqa: E402

from nwkit import gaussian  # noqa: E402

BASELINE = gaussian.is_diagonal


def bulk_is_diagonal(matrix):
    matrix = np.asarray(matrix)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        return False
    return bool(np.count_nonzero(matrix) == np.count_nonzero(np.diag(matrix)))


def guarded_is_diagonal(matrix):
    matrix = np.asarray(matrix)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        return False
    if matrix.dtype.kind not in "biufc":
        return BASELINE(matrix)
    if len(matrix) == 0:
        return True
    if np.any(matrix[0, 1:]):
        return False
    return bool(np.count_nonzero(matrix[1:]) == np.count_nonzero(np.diag(matrix)[1:]))


@contextmanager
def prototype():
    gaussian.is_diagonal = guarded_is_diagonal
    try:
        yield
    finally:
        gaussian.is_diagonal = BASELINE


def whole_factor(matrix, detector):
    previous = gaussian.is_diagonal
    gaussian.is_diagonal = detector
    try:
        return gaussian.factor_covariance(matrix)
    finally:
        gaussian.is_diagonal = previous


def summary(values):
    median = statistics.median(values)
    return {
        "samples_seconds_per_call": values,
        "median_seconds": median,
        "min_seconds": min(values),
        "max_seconds": max(values),
        "relative_range": (max(values) - min(values)) / median,
    }


def allocation(function):
    tracemalloc.start()
    try:
        function()
        return tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()


def timed_case(name, matrix, factor=False):
    detectors = {
        "baseline": BASELINE,
        "bulk": bulk_is_diagonal,
        "guarded": guarded_is_diagonal,
    }
    functions = {
        key: (
            (lambda detector=detector: whole_factor(matrix, detector))
            if factor
            else (lambda detector=detector: detector(matrix))
        )
        for key, detector in detectors.items()
    }
    expected = functions["baseline"]()
    for function in functions.values():
        actual = function()
        np.testing.assert_array_equal(actual, expected)
    loops = {}
    # Calibrate each version separately so even microsecond early exits have
    # substantial sample duration. Calibrations are outside recorded timings.
    for key, function in functions.items():
        count = 1
        while True:
            start = time.perf_counter()
            for _ in range(count):
                function()
            elapsed = time.perf_counter() - start
            if elapsed >= 0.03:
                break
            count *= 2
        loops[key] = count
    samples = {key: [] for key in functions}
    for repeat in range(7):
        order = list(functions)
        if repeat % 2:
            order.reverse()
        for key in order:
            start = time.perf_counter()
            for _ in range(loops[key]):
                actual = functions[key]()
            samples[key].append((time.perf_counter() - start) / loops[key])
            np.testing.assert_array_equal(actual, expected)
    records = {key: summary(values) for key, values in samples.items()}
    for key in records:
        records[key]["calls_per_sample"] = loops[key]
        records[key]["peak_python_bytes"] = allocation(functions[key])
        records[key]["speedup_over_baseline"] = (
            records["baseline"]["median_seconds"] / records[key]["median_seconds"]
        )
    return {
        "name": name,
        "shape": matrix.shape,
        "dtype": str(matrix.dtype),
        "strides": matrix.strides,
        "scope": "factor_covariance" if factor else "is_diagonal",
        "equivalent_exactly": True,
        "measurements": records,
    }


def edge_checks():
    cases = {
        "empty": np.empty((0, 0)),
        "scalar": np.asarray(1),
        "vector": np.ones(3),
        "nonsquare": np.ones((2, 3)),
        "zeros": np.zeros((3, 3)),
        "signed_zero": np.diag([0.0, -0.0, 2.0]),
        "nan_diagonal": np.diag([np.nan, 2.0]),
        "inf_diagonal": np.diag([np.inf, -np.inf]),
        "nan_offdiagonal": np.array([[1, np.nan], [0, 2]]),
        "inf_offdiagonal": np.array([[1, 0], [np.inf, 2]]),
        "tiny_offdiagonal": np.array([[1, 0], [np.nextafter(0.0, 1.0), 2]]),
        "complex_diagonal": np.diag([1j, 2 + 3j]),
        "complex_offdiagonal": np.array([[1, 1j], [0, 2]]),
        "integer": np.diag(np.arange(4)),
        "bool": np.eye(3, dtype=bool),
        "object": np.diag(np.array([None, 1, 2], dtype=object)),
        "strings": np.array([["a", ""], ["", "b"]]),
        "fortran": np.asfortranarray(np.diag(np.arange(4))),
        "strided": np.diag(np.arange(8))[::2, ::2],
        "negative_stride": np.diag(np.arange(8))[::-1, ::-1],
    }
    results = {}
    for name, matrix in cases.items():
        expected = BASELINE(matrix)
        assert guarded_is_diagonal(matrix) == expected, name
        assert bulk_is_diagonal(matrix) == expected, name
        results[name] = {"baseline": bool(expected), "guarded_exact": True}

    class NeverBoolean:
        def __bool__(self):
            raise RuntimeError("diagonal truth value must not be inspected")

    objects = np.zeros((2, 2), dtype=object)
    objects[0, 0] = NeverBoolean()
    assert BASELINE(objects) is True
    assert guarded_is_diagonal(objects) is True
    try:
        bulk_is_diagonal(objects)
    except RuntimeError as error:
        results["object_diagonal_truth_contract"] = {
            "baseline": True,
            "guarded_exact": True,
            "bulk_error": str(error),
        }
    else:
        raise AssertionError("bulk should inspect the forbidden diagonal object")
    return results


def main():
    import scipy

    from nwkit import __version__

    with LOCK.open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        output = {
            "baseline_commit": "d7116f80b9979739c69a5431a5e34b822795408b",
            "environment": {
                "platform": platform.platform(),
                "python": platform.python_version(),
                "numpy": np.__version__,
                "scipy": scipy.__version__,
                "nwkit": __version__,
                "threads": {
                    name: os.environ.get(name)
                    for name in (
                        "OPENBLAS_NUM_THREADS",
                        "OMP_NUM_THREADS",
                        "VECLIB_MAXIMUM_THREADS",
                        "MKL_NUM_THREADS",
                    )
                },
                "host_contention": "Other user-owned CPU workloads were active.",
            },
            "edge_checks": edge_checks(),
            "cases": [],
        }
        for size in (64, 512, 2048):
            diagonal = np.diag(np.linspace(0.5, 2.0, size))
            dense = np.full((size, size), 0.1) + diagonal
            late = diagonal.copy()
            late[-2, -1] = late[-1, -2] = 0.1
            for label, matrix in (
                ("diagonal", diagonal),
                ("dense_early_exit", dense),
                ("late_offdiagonal", late),
            ):
                result = timed_case(f"{label}_{size}", matrix)
                output["cases"].append(result)
                print(result["name"], flush=True)
            if size <= 512:
                for label, matrix in (
                    ("diagonal", diagonal),
                    ("dense", dense),
                    ("late_offdiagonal", late),
                ):
                    result = timed_case(f"factor_{label}_{size}", matrix, factor=True)
                    output["cases"].append(result)
                    print(result["name"], flush=True)
        view = np.diag(np.linspace(0.5, 2.0, 1024))[::2, ::2]
        output["cases"].append(timed_case("strided_diagonal_512", view))
        output["cases"].append(
            timed_case("fortran_diagonal_512", np.asfortranarray(view))
        )
        (HERE / "results.json").write_text(json.dumps(output, indent=2) + "\n")


if __name__ == "__main__":
    main()
