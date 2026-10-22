"""Benchmark grouped scores for the existing cluster-HC1 estimator."""

import argparse
import fcntl
import inspect
import json
import math
import os
import platform
import statistics
import subprocess
import sys
import time
import tracemalloc
from functools import partial
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
LOCK = Path("/tmp/nwkit-perf-20261022-9fYwTD/measurement.lock")
for key in (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "MKL_NUM_THREADS",
):
    os.environ[key] = "1"
sys.path.insert(0, str(ROOT))

OLD = """    meat = np.zeros((num_parameters, num_parameters), dtype=float)
    for event_index in range(n_events):
        selected = inverse == event_index
        score = design[selected].T @ (weights[selected] * residual[selected])
        meat += np.outer(score, score)
"""
NEW = """    scores = np.column_stack([
        np.bincount(inverse, weights=design[:, column] * (weights * residual), minlength=n_events)
        for column in range(num_parameters)
    ])
    meat = scores.T @ scores
"""


def candidate_function():
    from nwkit import regress

    source = inspect.getsource(regress._fit_model)
    if source.count(OLD) != 1:
        raise RuntimeError("Baseline score loop changed; review prototype first.")
    namespace = dict(regress.__dict__)
    exec(compile(source.replace(OLD, NEW), str(HERE / "prototype"), "exec"), namespace)
    return namespace["_fit_model"]


def kernel_baseline(design, weights, residual, inverse, n_events):
    import numpy as np

    meat = np.zeros((design.shape[1],) * 2)
    for event in range(n_events):
        selected = inverse == event
        score = design[selected].T @ (weights[selected] * residual[selected])
        meat += np.outer(score, score)
    return meat


def kernel_candidate(design, weights, residual, inverse, n_events):
    import numpy as np

    weighted_residual = weights * residual
    scores = np.column_stack(
        [
            np.bincount(
                inverse,
                weights=design[:, column] * weighted_residual,
                minlength=n_events,
            )
            for column in range(design.shape[1])
        ]
    )
    return scores.T @ scores


def kernel_add_at(design, weights, residual, inverse, n_events):
    import numpy as np

    scores = np.zeros((n_events, design.shape[1]))
    np.add.at(scores, inverse, design * (weights * residual)[:, None])
    return scores.T @ scores


def measure(functions, repeat=7):
    records = {name: [] for name in functions}
    batches = {}
    for name, function in functions.items():
        start = time.perf_counter()
        function()
        elapsed = time.perf_counter() - start
        batches[name] = min(10000, max(1, math.ceil(0.04 / elapsed)))
    for repetition in range(repeat):
        names = list(functions)
        if repetition % 2:
            names.reverse()
        for name in names:
            start = time.perf_counter()
            for _ in range(batches[name]):
                functions[name]()
            records[name].append((time.perf_counter() - start) / batches[name])
    result = {}
    for name, function in functions.items():
        tracemalloc.start()
        function()
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        samples = records[name]
        result[name] = {
            "seconds": samples,
            "operations_per_sample": batches[name],
            "median_seconds": statistics.median(samples),
            "min_seconds": min(samples),
            "max_seconds": max(samples),
            "peak_python_bytes": peak,
        }
    for name in functions:
        result[name]["speedup_over_baseline"] = (
            result["baseline"]["median_seconds"] / result[name]["median_seconds"]
        )
    return result


def compare_rows(first, second):
    max_absolute = 0.0
    max_scaled = 0.0
    checked = 0
    assert len(first) == len(second)
    for left, right in zip(first, second, strict=True):
        assert left.keys() == right.keys()
        for key in left:
            a, b = left[key], right[key]
            if isinstance(a, (float, int)):
                checked += 1
                if math.isnan(a):
                    assert math.isnan(b)
                    continue
                difference = abs(a - b)
                max_absolute = max(max_absolute, difference)
                max_scaled = max(max_scaled, difference / max(1.0, abs(a)))
                assert math.isclose(a, b, rel_tol=2e-11, abs_tol=2e-12), (key, a, b)
            else:
                assert a == b, (key, a, b)
    return {
        "numeric_fields": checked,
        "max_absolute": max_absolute,
        "max_scaled": max_scaled,
    }


def fixture(n, g, p, weighting, seed=442):
    import numpy as np
    import pandas as pd

    from nwkit.model_matrix import PredictorTerm

    rng = np.random.default_rng(seed)
    inverse = np.arange(n) % g
    rng.shuffle(inverse)
    design = rng.normal(size=(n, p))
    beta = rng.normal(size=p)
    residual = rng.normal(size=n) + rng.normal(size=g)[inverse]
    variance = np.exp(rng.normal(scale=0.5, size=n))
    response = design @ beta + residual
    counts = np.bincount(inverse, minlength=g)
    weights = 1.0 / counts[inverse] if weighting == "event" else np.ones(n)
    terms = [f"x{i}" for i in range(p)]
    metadata = {name: PredictorTerm(name, name, "continuous") for name in terms}
    dataframe = pd.DataFrame(
        {
            "tree_id": "synthetic-family",
            "trait": "response",
            "raw_contrast": response * np.sqrt(variance),
            "contrast_variance": variance,
            "_predictor_values": list(design * np.sqrt(variance)[:, None]),
            "species_event_id": [f"event-{index}" for index in inverse],
            "lineage_clade_id": [f"lineage-{index % 17}" for index in inverse],
        }
    )
    arguments = dict(
        confidence_level=0.95,
        event_weighting=weighting,
        coverage_policy="complete",
        excluded_ineligible=0,
        excluded_coverage=0,
        predictor_metadata=metadata,
        predictor_groups={name: (name,) for name in terms},
    )
    return dataframe, terms, arguments, (design, weights, residual, inverse, g)


def check_failures(baseline, candidate):
    import numpy as np

    rows = []
    for name in ("few-events", "rank-deficient", "zero-variance"):
        frame, terms, kwargs, _ = fixture(40, 10, 3, "event")
        if name == "few-events":
            frame["species_event_id"] = "single"
        elif name == "rank-deficient":
            frame["_predictor_values"] = [np.ones(3) for _ in range(len(frame))]
        else:
            frame["contrast_variance"] = 0.0
        errors = []
        for function in (baseline, candidate):
            try:
                with np.errstate(all="ignore"):
                    function(frame, terms, **kwargs)
            except Exception as error:
                errors.append((type(error).__name__, str(error)))
            else:
                errors.append(None)
        assert errors[0] is not None and errors[0] == errors[1], (name, errors)
        rows.append({"case": name, "exception": errors[0]})
    return rows


def run_benchmarks(short_only=False):
    import numpy as np
    import pandas as pd
    import scipy

    from nwkit import __version__, regress

    candidate = candidate_function()
    baseline = regress._fit_model
    results = {
        "baseline_commit": "d7116f80b9979739c69a5431a5e34b822795408b",
        "environment": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "pandas": pd.__version__,
            "nwkit": __version__,
            "threads": 1,
            "exclusive_lock": str(LOCK),
        },
        "whole_helper": [],
        "kernel": [],
    }
    workloads = ((200, 50, 3), (3000, 600, 8), (12000, 3000, 4))
    if short_only:
        previous = json.loads((HERE / "results-initial.json").read_text())
        for field in ("whole_helper", "kernel"):
            results[field] = [
                row for row in previous[field] if row["observations"] != 200
            ]
        results["measurement_notes"] = (
            "200-row cases remeasured with 40ms target batches; larger cases retained from initial run. External CPU contention was observed by the parent outside our advisory lock; effects remain exploratory."
        )
        workloads = workloads[:1]
    for n, g, p in workloads:
        for weighting in ("event", "contrast"):
            frame, terms, kwargs, kernel_args = fixture(n, g, p, weighting)
            first = baseline(frame, terms, **kwargs)
            second = candidate(frame, terms, **kwargs)
            result = {
                "observations": n,
                "events": g,
                "predictors": p,
                "event_weighting": weighting,
                "equivalence": compare_rows(first, second),
                "measurements": measure(
                    {
                        "baseline": partial(baseline, frame, terms, **kwargs),
                        "candidate": partial(candidate, frame, terms, **kwargs),
                    }
                ),
            }
            results["whole_helper"].append(result)
            expected = kernel_baseline(*kernel_args)
            residuals = {}
            for name, function in (
                ("candidate", kernel_candidate),
                ("add_at", kernel_add_at),
            ):
                actual = function(*kernel_args)
                np.testing.assert_allclose(actual, expected, rtol=2e-12, atol=2e-12)
                residuals[name] = float(np.max(np.abs(actual - expected)))
            results["kernel"].append(
                {
                    "observations": n,
                    "events": g,
                    "predictors": p,
                    "event_weighting": weighting,
                    "max_absolute_differences": residuals,
                    "measurements": measure(
                        {
                            "baseline": partial(kernel_baseline, *kernel_args),
                            "candidate": partial(kernel_candidate, *kernel_args),
                            "add_at": partial(kernel_add_at, *kernel_args),
                        }
                    ),
                }
            )
    results["failure_equivalence"] = check_failures(baseline, candidate)
    # Categorical omnibus consumer and exact-fit degeneracy also preserve rows.
    frame, terms, kwargs, _ = fixture(200, 50, 3, "event")
    from nwkit.model_matrix import PredictorTerm

    kwargs["predictor_metadata"] = {
        term: PredictorTerm(term, "factor", "categorical", level=term) for term in terms
    }
    kwargs["predictor_groups"] = {"factor": tuple(terms)}
    results["omnibus_equivalence"] = compare_rows(
        baseline(frame, terms, **kwargs), candidate(frame, terms, **kwargs)
    )
    original_response = frame["raw_contrast"].copy()
    for scale in (1e-12, 1e12):
        frame["raw_contrast"] = original_response * scale
        results[f"response_rescale_{scale}"] = compare_rows(
            baseline(frame, terms, **kwargs), candidate(frame, terms, **kwargs)
        )
    (HERE / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(
        json.dumps(
            {
                "whole_helper": [
                    {
                        k: row[k]
                        for k in (
                            "observations",
                            "events",
                            "predictors",
                            "event_weighting",
                        )
                    }
                    | {
                        "speedup": row["measurements"]["candidate"][
                            "speedup_over_baseline"
                        ]
                    }
                    for row in results["whole_helper"]
                ]
            },
            indent=2,
        )
    )


def run_tests():
    environment = dict(os.environ)
    environment["PYTHONPATH"] = os.pathsep.join((str(HERE), str(ROOT)))
    command = [
        sys.executable,
        "tools/check.py",
        "test",
        "--",
        "-p",
        "candidate03_plugin",
        "tests/test_regress.py",
        "tests/test_numerical_invariance.py",
        "tests/test_regression_inference_contract.py",
        "-k",
        "cluster or event_weighting or gaussian_fit_is_invariant",
        "-rs",
    ]
    completed = subprocess.run(
        command, cwd=ROOT, env=environment, text=True, capture_output=True
    )
    (HERE / "tests.log").write_text(
        "Command: " + " ".join(command) + "\n" + completed.stdout + completed.stderr
    )
    print(completed.stdout)
    print(completed.stderr)
    if completed.returncode:
        raise SystemExit(completed.returncode)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tests", action="store_true")
    parser.add_argument("--short-only", action="store_true")
    args = parser.parse_args()
    with LOCK.open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if args.tests:
            run_tests()
        else:
            run_benchmarks(short_only=args.short_only)


if __name__ == "__main__":
    main()
