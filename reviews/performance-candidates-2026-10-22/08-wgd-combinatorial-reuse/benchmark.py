"""Evaluate local log-gamma reuse in the unbounded count transition kernel."""

import argparse
import fcntl
import inspect
import json
import os
import platform
import statistics
import subprocess
import sys
import time
import tracemalloc
from contextlib import contextmanager
from functools import partial
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
LOCK = Path("/tmp/nwkit-perf-20261022-9fYwTD/measurement.lock")
sys.path.insert(0, str(ROOT))
for key in (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "MKL_NUM_THREADS",
):
    os.environ[key] = "1"


def build_candidate():
    from nwkit import wgd_count_model as model

    source = inspect.getsource(model._birth_death_log_transition)
    marker = "    descendants = np.arange(1, max_count + 1)[None, :]\n"
    insertion = (
        marker
        + """    gamma_table = gammaln(np.arange(max_count + 2))
    survivor_axis = np.arange(1, max_count + 1)[:, None]
    gamma_difference = gamma_table[np.maximum(descendants - survivor_axis + 1, 0)]
"""
    )
    assert source.count(marker) == 1
    source = source.replace(marker, insertion)
    replacements = {
        "gammaln(count + 1)": "gamma_table[count + 1]",
        "gammaln(surviving + 1)": "gamma_table[2 : count + 2, None]",
        "gammaln(count - surviving + 1)": "gamma_table[1 : count + 1, None][::-1]",
        "gammaln(descendants)": "gamma_table[None, 1 : max_count + 1]",
        "gammaln(surviving)": "gamma_table[1 : count + 1, None]",
        "gammaln(descendants - surviving + 1)": "gamma_difference[:count]",
    }
    for old, new in replacements.items():
        assert source.count(old) == 1, old
        source = source.replace(old, new)
    namespace = {}
    exec(compile(source, str(HERE / "benchmark.py"), "exec"), model.__dict__, namespace)
    return namespace["_birth_death_log_transition"]


@contextmanager
def prototype():
    from nwkit import wgd_count_model as model

    baseline = model._birth_death_log_transition
    model._birth_death_log_transition = build_candidate()
    try:
        yield
    finally:
        model._birth_death_log_transition = baseline


def equivalent(first, second):
    import numpy as np

    first, second = np.asarray(first), np.asarray(second)
    assert first.shape == second.shape
    np.testing.assert_array_equal(np.isneginf(first), np.isneginf(second))
    np.testing.assert_array_equal(np.isposinf(first), np.isposinf(second))
    np.testing.assert_array_equal(np.isnan(first), np.isnan(second))
    finite = np.isfinite(first) & np.isfinite(second)
    maximum = float(np.max(np.abs(first[finite] - second[finite]), initial=0))
    np.testing.assert_allclose(first, second, rtol=1e-13, atol=1e-12)
    return {"array_equal": bool(np.array_equal(first, second)), "max_abs": maximum}


def summarize(values):
    median = statistics.median(values)
    return {
        "samples_seconds_per_call": values,
        "median_seconds": median,
        "min_seconds": min(values),
        "max_seconds": max(values),
        "relative_range": (max(values) - min(values)) / median,
    }


def peak(function):
    tracemalloc.start()
    try:
        function()
        return tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()


def measure(functions):
    warmup = {}
    outputs = {}
    for name, function in functions.items():
        start = time.perf_counter()
        outputs[name] = function()
        warmup[name] = time.perf_counter() - start
    comparisons = {
        name: equivalent(outputs["baseline"], value) for name, value in outputs.items()
    }
    count = max(1, min(1000, int(0.03 / min(warmup.values())) + 1))
    samples = {name: [] for name in functions}
    for repetition in range(7):
        names = list(functions)
        if repetition % 2:
            names.reverse()
        for name in names:
            start = time.perf_counter()
            for _ in range(count):
                value = functions[name]()
            samples[name].append((time.perf_counter() - start) / count)
            equivalent(outputs["baseline"], value)
    timings = {name: summarize(values) for name, values in samples.items()}
    for record in timings.values():
        record["speedup"] = (
            timings["baseline"]["median_seconds"] / record["median_seconds"]
        )
    return {
        "batch_calls": count,
        "equivalence": comparisons,
        "timings": timings,
        "peak_python_bytes": {
            name: peak(function) for name, function in functions.items()
        },
    }


def tree_fixture(counts, detection=None):
    import numpy as np

    from nwkit.wgd_count_model import CountLikelihood, CountTree

    tree = CountTree(
        (-1, 0, 0, 1, 1),
        (0.0, 0.4, 0.7, 0.25, 0.45),
        (2, 3, 4),
        ("A", "B", "C"),
        (0, 1, 2, 3, 4),
        ("ABC", "BC", "A", "B", "C"),
    )
    return CountLikelihood(tree, np.asarray(counts), detection=detection)


def tree_call(model, rates, bound, event, kernel):
    from nwkit import wgd_count_model as module

    previous = module._birth_death_log_transition
    module._birth_death_log_transition = kernel
    try:
        return model.family_log_likelihoods(rates, 1.0, bound, event)
    finally:
        module._birth_death_log_transition = previous


def edge_checks(baseline, candidate):
    import numpy as np

    from nwkit.wgd_count_model import birth_death_transition

    cases = []
    for duplication, loss, time_value in (
        (0, 0, 1),
        (0.2, 0.3, 0),
        (0, 0.3, 0.9),
        (0.2, 0, 0.9),
        (0.2, 0.2, 0.9),
        (0.2, 0.200000000001, 0.9),
        (1e-8, 0, 0.4),
        (0.2, 0.4, 100),
        (1e-200, 2e-200, 1e200),
        (2e100, 1e100, 1e-100),
    ):
        for bound in (1, 20, 86):
            args = duplication, loss, time_value, bound
            first, second = baseline(*args), candidate(*args)
            result = equivalent(first, second)
            probability = birth_death_transition(*args)
            np.testing.assert_allclose(
                np.exp(second), probability, rtol=2e-11, atol=2e-13
            )
            cases.append({"args": args, **result})
    failures = []
    for args in ((-0.2, 0.3, 1, 4), (0.2, np.nan, 1, 4), (0.2, 0.3, np.inf, 4)):
        errors = []
        for function in (baseline, candidate):
            try:
                function(*args)
            except Exception as exc:
                errors.append((type(exc).__name__, str(exc)))
            else:
                raise AssertionError("Expected invalid rate/time to fail")
        assert errors[0] == errors[1]
        failures.append({"args": [str(value) for value in args], "error": errors[0]})
    return {"kernels": cases, "invalid_inputs": failures}


def run_benchmarks():
    import numpy as np
    import scipy

    from nwkit import __version__
    from nwkit import wgd_count_model as module
    from nwkit.wgd_count_model import MultiplicationEvent

    baseline = module._birth_death_log_transition
    candidate = build_candidate()
    result = {
        "baseline_commit": "d7116f80b9979739c69a5431a5e34b822795408b",
        "environment": {
            "platform": platform.platform(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "nwkit": __version__,
            "threads": 1,
            "exclusive_lock": str(LOCK),
            "contention": "Other user jobs active outside this advisory lock",
        },
        "edges": edge_checks(baseline, candidate),
        "kernels": [],
        "tree_likelihoods": [],
    }
    for bound in (20, 60, 120):
        args = (0.2, 0.3, 0.9, bound)
        result["kernels"].append(
            {
                "args": args,
                **measure(
                    {
                        "baseline": lambda args=args: baseline(*args),
                        "candidate": lambda args=args: candidate(*args),
                    }
                ),
            }
        )
    for name, counts, rates, bound, event, detection in (
        (
            "ordinary_probability_control",
            [[1, 2, 3], [4, 2, 1], [0, 1, np.nan]],
            [[0.2, 0.3]],
            40,
            None,
            None,
        ),
        (
            "rare_yule_missing",
            [[84, 1, 1], [84, 1, np.nan], [84, np.nan, 1]],
            [[1e-8, 0]],
            86,
            None,
            None,
        ),
        (
            "rare_yule_pulse_detection",
            [[84, 1, 1], [84, 1, np.nan]],
            [[1e-8, 0]],
            86,
            MultiplicationEvent(2, 0.3),
            np.array([0.8, 0.9, 1.0]),
        ),
    ):
        model = tree_fixture(counts, detection)
        counters = {}
        for is_candidate in (False, True):
            original = module._birth_death_log_transition
            calls = [0]
            chosen = candidate if is_candidate else baseline

            def counted(*args, chosen=chosen, calls=calls):
                calls[0] += 1
                return chosen(*args)

            module._birth_death_log_transition = counted
            try:
                model.family_log_likelihoods(rates, 1.0, bound, event)
            finally:
                module._birth_death_log_transition = original
            counters["candidate" if is_candidate else "baseline"] = calls[0]
        result["tree_likelihoods"].append(
            {
                "name": name,
                "families": len(counts),
                "tips": 3,
                "max_count": bound,
                "rates": rates,
                "event": None if event is None else vars(event),
                "log_kernel_calls": counters,
                **measure(
                    {
                        "baseline": partial(
                            tree_call, model, rates, bound, event, baseline
                        ),
                        "candidate": partial(
                            tree_call, model, rates, bound, event, candidate
                        ),
                    }
                ),
            }
        )
    path = HERE / "results.json"
    path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    for record in result["kernels"]:
        print("kernel", record["args"], record["timings"]["candidate"]["speedup"])
    for record in result["tree_likelihoods"]:
        print(record["name"], record["timings"]["candidate"]["speedup"])


def run_tests(consumer=False):
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(HERE)
    environment["NWKIT_CI_DIAGNOSTICS"] = "0"
    command = [
        sys.executable,
        "tools/check.py",
        "test",
        "--",
        "-p",
        "candidate08_plugin",
        (
            "tests/test_wgd_count.py::test_real_cli_and_metadata_no_uncalibrated_probability"
            if consumer
            else "tests/test_wgd_count_model.py"
        ),
        "-rs",
    ]
    completed = subprocess.run(
        command, cwd=ROOT, env=environment, capture_output=True, text=True
    )
    (HERE / ("consumer-tests.log" if consumer else "tests.log")).write_text(
        "+ " + " ".join(command) + "\n" + completed.stdout + completed.stderr
    )
    print(completed.stdout[-3000:])
    if completed.returncode:
        raise SystemExit(completed.returncode)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tests", action="store_true")
    parser.add_argument("--consumer-tests", action="store_true")
    args = parser.parse_args()
    LOCK.parent.mkdir(parents=True, exist_ok=True)
    with LOCK.open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if args.tests or args.consumer_tests:
            run_tests(consumer=args.consumer_tests)
        else:
            run_benchmarks()


if __name__ == "__main__":
    main()
