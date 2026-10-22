"""Compare complete covariance-profile fits with combined design/response RHS."""

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
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
LOCK = Path("/tmp/nwkit-perf-20261022-9fYwTD/measurement.lock")
for key in (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "MKL_NUM_THREADS",
):
    os.environ[key] = "1"
sys.path.insert(0, str(ROOT))


@contextmanager
def prototype():
    import nwkit.regress as regress

    original = regress._profile_covariance_fit
    source = inspect.getsource(original)
    before = "inverse_design = _solve_positive_definite(cholesky, design)"
    after = """inverse_inputs = _solve_positive_definite(
                    cholesky, np.column_stack((design, response))
                )
                inverse_design = inverse_inputs[:, :-1]
                inverse_response = inverse_inputs[:, -1]"""
    assert source.count(before) == 2
    source = source.replace(before, after)
    before = "design.T @ _solve_positive_definite(cholesky, response)"
    assert source.count(before) == 2
    source = source.replace(before, "design.T @ inverse_response")
    namespace = {}
    exec(
        compile(source, str(HERE / "benchmark.py"), "exec"), regress.__dict__, namespace
    )
    candidate = namespace["_profile_covariance_fit"]
    regress._profile_covariance_fit = candidate
    # Preserve imported aliases when an existing consumer was imported earlier.
    aliases = []
    for module in list(sys.modules.values()):
        if module is None or not getattr(module, "__name__", "").startswith("nwkit."):
            continue
        if getattr(module, "_profile_covariance_fit", None) is original:
            aliases.append(module)
            module._profile_covariance_fit = candidate
    try:
        yield
    finally:
        regress._profile_covariance_fit = original
        for module in aliases:
            module._profile_covariance_fit = original


def fixture(kind, n, p, reml, seed=606, near=False, scale=1.0):
    import numpy as np

    rng = np.random.default_rng(seed)
    design = np.column_stack((np.ones(n), rng.normal(size=(n, p - 1))))
    response = design @ np.linspace(0.2, 1.0, p) + rng.normal(size=n) * (
        1e-8 if near else 1.0
    )
    options = {"reml": reml}
    if kind.startswith("dense"):
        positions = np.arange(n)
        component = 0.7 ** np.abs(
            positions[:, None] - positions[None, :]
        ) + 0.2 * np.eye(n)
        fixed = np.repeat(0.2 if kind == "dense-fixed" else 0.0, n)
        components = [("tree", component)]
    elif kind == "diagonal":
        fixed = np.zeros(n)
        components = [("tree", rng.uniform(0.5, 1.5, n))]
    elif kind == "lowrank":
        loading = rng.normal(size=(n, 16)) / 4
        fixed = rng.uniform(0.1, 0.3, n)
        components = [("tree", None)]
        options["component_factors"] = {"tree": loading}
    else:
        raise ValueError(kind)
    # Scaling both response and fixed observation variance preserves units.
    return response * scale, design, fixed * scale**2, components, options


def call_fit(arguments, candidate=False, count=False):
    import nwkit.regress as regress

    response, design, fixed, components, options = arguments
    calls = []
    original_solve = regress._solve_positive_definite
    if count:

        def counting(factor, values):
            calls.append(list(values.shape))
            return original_solve(factor, values)

        regress._solve_positive_definite = counting
    try:
        if candidate:
            with prototype():
                result = regress._profile_covariance_fit(
                    response, design, fixed, components, **options
                )
        else:
            result = regress._profile_covariance_fit(
                response, design, fixed, components, **options
            )
        return result, calls
    finally:
        regress._solve_positive_definite = original_solve


def compare(first, second):
    import numpy as np

    from nwkit.gaussian import materialize_covariance

    assert first.keys() == second.keys()
    differences = {}
    for name, value in first.items():
        other = second[name]
        if name == "cholesky":
            # Factor dataclasses have multiple internal representations; fitted
            # covariance and inverse residual are independently compared below.
            assert type(value) is type(other)
            continue
        if name == "covariance":
            value, other = materialize_covariance(value), materialize_covariance(other)
        if isinstance(value, dict):
            assert value.keys() == other.keys()
            value, other = list(value.values()), list(other.values())
        if name == "optimizer_message":
            differences[name] = {
                "baseline": value,
                "candidate": other,
                "exact_match": value == other,
            }
            continue
        if isinstance(value, (bool, str)):
            assert value == other, (name, value, other)
            continue
        a, b = np.asarray(value), np.asarray(other)
        absolute = float(np.max(np.abs(a - b), initial=0.0))
        scaled = float(np.max(np.abs(a - b) / np.maximum(1, np.abs(a)), initial=0.0))
        np.testing.assert_allclose(a, b, rtol=1e-4, atol=1e-6)
        differences[name] = {"max_absolute": absolute, "max_scaled": scaled}
    return differences


def peak(function):
    tracemalloc.start()
    try:
        function()
        return tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()


def benchmark_case(kind, n, p, reml):
    arguments = fixture(kind, n, p, reml)
    first, baseline_calls = call_fit(arguments, False, True)
    second, candidate_calls = call_fit(arguments, True, True)
    differences = compare(first, second)
    # Prototype compilation and installation are not part of the numerical
    # algorithm; keep them outside timed batches.
    import nwkit.regress as regress

    baseline = regress._profile_covariance_fit
    with prototype():
        candidate = regress._profile_covariance_fit
    response, design, fixed, components, options = arguments
    functions = {
        "baseline": lambda: baseline(response, design, fixed, components, **options),
        "candidate": lambda: candidate(response, design, fixed, components, **options),
    }
    elapsed = []
    for function in functions.values():
        start = time.perf_counter()
        function()
        elapsed.append(time.perf_counter() - start)
    batch = max(1, min(100, int(0.04 / min(elapsed)) + 1))
    samples = {name: [] for name in functions}
    for repetition in range(7):
        order = list(functions) if repetition % 2 == 0 else list(reversed(functions))
        for name in order:
            start = time.perf_counter()
            for _ in range(batch):
                result = functions[name]()
            samples[name].append((time.perf_counter() - start) / batch)
            compare(first, result)
    timings = {
        name: {
            "samples_seconds": times,
            "median_seconds": statistics.median(times),
            "min_seconds": min(times),
            "max_seconds": max(times),
            "peak_python_bytes": peak(functions[name]),
        }
        for name, times in samples.items()
    }
    return {
        "kind": kind,
        "observations": n,
        "predictors": p,
        "reml": reml,
        "batch": batch,
        "timings": timings,
        "speedup": timings["baseline"]["median_seconds"]
        / timings["candidate"]["median_seconds"],
        "equivalence": differences,
        "solve_counts": {
            "baseline": len(baseline_calls),
            "candidate": len(candidate_calls),
        },
        "solve_shapes": {
            "baseline_first_six": baseline_calls[:6],
            "candidate_first_six": candidate_calls[:6],
        },
    }


def edge_checks():
    import numpy as np

    records = []
    cases = [
        ("dense", 32, 3, True, False, 1e-6),
        ("dense", 32, 3, False, False, 1e6),
        ("dense", 32, 3, True, True, 1.0),
        ("diagonal", 32, 3, False, True, 1.0),
        ("lowrank", 32, 3, True, False, 1.0),
    ]
    for kind, n, p, reml, near, scale in cases:
        args = fixture(kind, n, p, reml, near=near, scale=scale)
        first, _ = call_fit(args)
        second, _ = call_fit(args, True)
        records.append(
            {
                "kind": kind,
                "reml": reml,
                "near_exact": near,
                "response_scale": scale,
                "equivalence": compare(first, second),
            }
        )
    for label in ("singular_component", "rank_deficient_design", "invalid_sampling"):
        args = list(fixture("dense", 16, 3, True))
        if label == "singular_component":
            args[3] = [("tree", np.ones((16, 16)))]
        elif label == "rank_deficient_design":
            args[1][:, 1] = args[1][:, 0]
        else:
            args[2][0] = -1
        outcomes = []
        for candidate in (False, True):
            try:
                call_fit(args, candidate)
                outcomes.append({"exception": None})
            except Exception as error:
                outcomes.append(
                    {"exception": type(error).__name__, "message": str(error)}
                )
        assert outcomes[0] == outcomes[1], (label, outcomes)
        records.append({"case": label, "outcomes": outcomes})
    args = list(fixture("dense", 32, 3, True))
    component = np.ones((32, 32)) + 1e-8 * np.eye(32)
    args[3] = [("tree", component)]
    first, _ = call_fit(args)
    second, _ = call_fit(args, True)
    records.append(
        {
            "case": "ill_conditioned_covariance",
            "condition_number": float(np.linalg.cond(component)),
            "equivalence": compare(first, second),
        }
    )
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tests", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    LOCK.parent.mkdir(parents=True, exist_ok=True)
    with LOCK.open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if args.tests:
            environment = os.environ.copy()
            environment["PYTHONPATH"] = str(HERE) + os.pathsep + str(ROOT)
            environment["PYTEST_PLUGINS"] = "candidate06_plugin"
            environment["NWKIT_CI_DIAGNOSTICS"] = "0"
            command = [
                sys.executable,
                "tools/check.py",
                "test",
                "--",
                "tests/test_numerical_invariance.py",
                "tests/test_regress.py::test_large_dense_gaussian_fit_is_rejected_before_optimization",
                "tests/test_regress.py::test_large_structured_gaussian_fit_remains_available",
                "tests/test_regress.py::test_diagonal_profile_fit_matches_dense_reference",
                "tests/test_regress.py::test_grouped_low_rank_factor_matches_dense_linear_algebra",
                "tests/test_regress.py::test_sparse_base_and_dense_update_use_nested_factorization",
                "tests/test_regress.py::test_replicate_reml_matches_gls_with_evolutionary_and_sampling_covariance",
                "-rs",
            ]
            result = subprocess.run(
                command, cwd=ROOT, env=environment, capture_output=True, text=True
            )
            if (HERE / "tests.log").exists():
                (HERE / "tests-initial.log").write_text(
                    (HERE / "tests.log").read_text()
                )
            (HERE / "tests.log").write_text(
                " ".join(command) + "\n" + result.stdout + result.stderr
            )
            print(result.stdout, result.stderr)
            raise SystemExit(result.returncode)
        import numpy as np
        import scipy

        from nwkit import __version__

        workloads = [
            ("dense", 128, 3, False),
            ("dense", 512, 8, True),
            ("dense-fixed", 128, 3, True),
            ("diagonal", 2000, 5, True),
            ("lowrank", 512, 4, True),
        ]
        results = {
            "baseline_commit": "d7116f80b9979739c69a5431a5e34b822795408b",
            "environment": {
                "platform": platform.platform(),
                "python": platform.python_version(),
                "numpy": np.__version__,
                "scipy": scipy.__version__,
                "nwkit": __version__,
                "threads": 1,
            },
            "scope": "complete _profile_covariance_fit excluding fixture/import/prototype installation",
            "workloads": [],
        }
        if args.resume and (HERE / "results.json").exists():
            results = json.loads((HERE / "results.json").read_text())
        completed = {
            (row["kind"], row["observations"], row["predictors"], row["reml"])
            for row in results["workloads"]
        }
        for case in workloads:
            if case in completed:
                continue
            result = benchmark_case(*case)
            results["workloads"].append(result)
            (HERE / "results.json").write_text(json.dumps(results, indent=2) + "\n")
            print(case, result["speedup"], flush=True)
        results["edge_checks"] = edge_checks()
        (HERE / "results.json").write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
