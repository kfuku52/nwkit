"""Observed-only intermediate PGLS covariance; retain full final output."""

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
sys.path.insert(0, str(ROOT))
for key in (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "MKL_NUM_THREADS",
):
    os.environ[key] = "1"

OLD_BLOCK = """        full_covariance = np.zeros((n_tips * n_traits,) * 2, dtype=float)
        for name, component in covariance_components.items():
            tip_covariance = dense_component(
                name, decoded if callable(component) else None
            )
            full_covariance += np.kron(trait_covariances[name], tip_covariance)
        if fixed_dense is not None:
            full_covariance += fixed_dense
        observed_covariance = full_covariance[
            np.ix_(observed_indices, observed_indices)
        ]
"""
NEW_BLOCK = """        if full_output or len(observed_indices) == n_tips * n_traits:
            full_covariance = np.zeros((n_tips * n_traits,) * 2, dtype=float)
            for name, component in covariance_components.items():
                tip_covariance = dense_component(
                    name, decoded if callable(component) else None
                )
                full_covariance += np.kron(trait_covariances[name], tip_covariance)
            if fixed_dense is not None:
                full_covariance += fixed_dense
            observed_covariance = full_covariance[
                np.ix_(observed_indices, observed_indices)
            ]
        else:
            full_covariance = None
            observed_covariance = np.zeros((len(observed_indices),) * 2, dtype=float)
            for name, component in covariance_components.items():
                tip_covariance = dense_component(
                    name, decoded if callable(component) else None
                )
                observed_covariance += (
                    trait_covariances[name][trait_pairs] * tip_covariance[tip_pairs]
                )
            if fixed_dense is not None:
                observed_covariance += fixed_dense[
                    np.ix_(observed_indices, observed_indices)
                ]
"""


def build_candidate():
    import nwkit.multivariate_pgls as module

    source = inspect.getsource(module.fit_multivariate_pgls)
    assert source.count(OLD_BLOCK) == 1
    source = source.replace(OLD_BLOCK, NEW_BLOCK)
    old = "    def state(parameters: np.ndarray):\n"
    new = """    observed_traits = observed_indices // n_tips
    observed_tips = observed_indices % n_tips
    trait_pairs = np.ix_(observed_traits, observed_traits)
    tip_pairs = np.ix_(observed_tips, observed_tips)

    def state(parameters: np.ndarray, *, full_output=False):
"""
    assert source.count(old) == 1
    source = source.replace(old, new)
    old = """    fitted_objective, coefficients, coefficient_covariance, fitted_covariance = state(
        result.x
    )
"""
    new = """    fitted_objective, coefficients, coefficient_covariance, fitted_covariance = state(
        result.x, full_output=True
    )
"""
    assert source.count(old) == 1
    source = source.replace(old, new)
    namespace = {}
    exec(
        compile(source, str(HERE / "benchmark.py"), "exec"), module.__dict__, namespace
    )
    return namespace["fit_multivariate_pgls"]


@contextmanager
def prototype():
    import nwkit.multivariate_pgls as module

    original = module.fit_multivariate_pgls
    candidate = build_candidate()
    aliases = []
    for name, loaded in tuple(sys.modules.items()):
        if name.startswith("nwkit.") and loaded is not None:
            for key, value in tuple(vars(loaded).items()):
                if value is original:
                    aliases.append((loaded, key, value))
                    setattr(loaded, key, candidate)
    try:
        yield
    finally:
        for loaded, key, value in aliases:
            setattr(loaded, key, value)


def compare(first, second):
    import numpy as np

    record = {}
    for key in first.__dataclass_fields__:
        a, b = getattr(first, key), getattr(second, key)
        if isinstance(a, dict):
            assert a.keys() == b.keys()
            record[key] = {name: compare_array(a[name], b[name]) for name in a}
        elif isinstance(a, np.ndarray):
            record[key] = compare_array(a, b)
        elif isinstance(a, float):
            np.testing.assert_allclose(a, b, rtol=1e-6, atol=1e-8)
            record[key] = {"absolute_difference": abs(a - b), "exact": a == b}
        else:
            assert a == b, (key, a, b)
            record[key] = {"exact": True, "value": a}
    return record


def compare_array(first, second):
    import numpy as np

    np.testing.assert_allclose(first, second, rtol=1e-6, atol=1e-8)
    return {
        "max_absolute": float(np.max(np.abs(first - second), initial=0)),
        "exact": bool(np.array_equal(first, second)),
        "shape": list(first.shape),
    }


def stats(samples):
    return {
        "samples_seconds_per_call": samples,
        "median_seconds": statistics.median(samples),
        "min_seconds": min(samples),
        "max_seconds": max(samples),
        "relative_range": (max(samples) - min(samples)) / statistics.median(samples),
    }


def peak(function):
    tracemalloc.start()
    try:
        function()
        return tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()


def timed(functions, repetitions=7):
    warm = {}
    for name, function in functions.items():
        start = time.perf_counter()
        function()
        warm[name] = time.perf_counter() - start
    batch = max(1, min(500, int(0.025 / min(warm.values())) + 1))
    samples = {name: [] for name in functions}
    for repetition in range(repetitions):
        names = list(functions)
        if repetition % 2:
            names.reverse()
        for name in names:
            start = time.perf_counter()
            for _ in range(batch):
                functions[name]()
            samples[name].append((time.perf_counter() - start) / batch)
    output = {name: stats(sample) for name, sample in samples.items()}
    return {
        "batch_calls": batch,
        "timings": output,
        "speedup": output["baseline"]["median_seconds"]
        / output["candidate"]["median_seconds"],
        "peak_python_bytes": {name: peak(fn) for name, fn in functions.items()},
    }


def fixture(n, d, fraction, *, fixed=False, multiple=False, reml=True):
    import numpy as np

    rng = np.random.default_rng(1010 + n + d)
    x = rng.normal(size=n)
    design = np.column_stack([np.ones(n), x])
    covariance = 0.65 ** np.abs(np.arange(n)[:, None] - np.arange(n)[None, :])
    trait_factor = np.eye(d) + 0.15 * np.tril(np.ones((d, d)), -1)
    noise = np.linalg.cholesky(covariance) @ rng.normal(size=(n, d)) @ trait_factor.T
    responses = 0.3 + x[:, None] * np.linspace(0.5, 1.2, d) + noise
    mask = rng.uniform(size=(n, d)) < fraction
    # Every retained tip has at least one observed response.
    mask[np.arange(n), np.arange(n) % d] = False
    responses[mask] = np.nan
    components = {"tree": covariance}
    if multiple:
        components["independent"] = np.eye(n)
    options = {"reml": reml}
    if fixed:
        options["fixed_covariance"] = np.full(n * d, 0.05)
    return responses, design, components, options


def full_case(n, d, fraction, *, fixed=False, multiple=False, reml=True):
    import numpy as np

    import nwkit.multivariate_pgls as module

    responses, design, components, options = fixture(
        n, d, fraction, fixed=fixed, multiple=multiple, reml=reml
    )
    candidate = build_candidate()
    functions = {
        "baseline": lambda: module.fit_multivariate_pgls(
            responses, design, components, **options
        ),
        "candidate": lambda: candidate(responses, design, components, **options),
    }
    first, second = functions["baseline"](), functions["candidate"]()
    differences = compare(first, second)
    measured = timed(functions)
    return {
        "tips": n,
        "traits": d,
        "missing_fraction_requested": fraction,
        "missing_fraction_actual": float(np.isnan(responses).mean()),
        "reml": reml,
        "fixed_noise": fixed,
        "multiple_components": multiple,
        "equivalence": differences,
        **measured,
    }


def kernel_case(n, d, fraction):
    import numpy as np

    responses, _, components, _ = fixture(n, d, fraction)
    covariance = components["tree"]
    observed = np.flatnonzero(np.isfinite(responses).ravel(order="F"))
    trait_indices, tip_indices = observed // n, observed % n
    trait_pairs = np.ix_(trait_indices, trait_indices)
    tip_pairs = np.ix_(tip_indices, tip_indices)
    trait_covariance = np.eye(d) + 0.2 * np.ones((d, d))

    def baseline():
        full = np.zeros((n * d, n * d))
        full += np.kron(trait_covariance, covariance)
        return full[np.ix_(observed, observed)]

    def candidate():
        selected = np.zeros((len(observed), len(observed)))
        selected += trait_covariance[trait_pairs] * covariance[tip_pairs]
        return selected

    return {
        "tips": n,
        "traits": d,
        "observed": len(observed),
        "missing_fraction_actual": float(np.isnan(responses).mean()),
        "equivalence": compare_array(baseline(), candidate()),
        **timed({"baseline": baseline, "candidate": candidate}),
    }


def edge_checks():
    import numpy as np

    import nwkit.multivariate_pgls as module

    candidate = build_candidate()
    records = []
    responses, design, components, options = fixture(20, 2, 0.5, fixed=True, reml=False)
    base = components["tree"]
    components = {"tree": lambda parameter: base + float(parameter) * np.eye(20)}
    options.update(
        evolution_parameter_bounds=(0.01, 0.5), evolution_parameter_initial=0.2
    )
    first = module.fit_multivariate_pgls(responses, design, components, **options)
    second = candidate(responses, design, components, **options)
    records.append(
        {"case": "callable_estimated_shape", "equivalence": compare(first, second)}
    )
    observed = np.flatnonzero(np.isfinite(responses).ravel(order="F"))
    full_design = np.kron(np.eye(responses.shape[1]), design)
    residual = (
        responses.ravel(order="F")[observed]
        - full_design[observed] @ second.coefficients.ravel()
    )
    selected = second.fitted_covariance[np.ix_(observed, observed)]
    sign, logdet = np.linalg.slogdet(selected)
    assert sign == 1
    expected = -0.5 * (
        len(observed) * np.log(2 * np.pi)
        + logdet
        + residual @ np.linalg.solve(selected, residual)
    )
    np.testing.assert_allclose(second.log_likelihood, expected, rtol=0, atol=1e-10)
    records.append(
        {
            "case": "independent_observed_likelihood",
            "absolute_difference": abs(second.log_likelihood - float(expected)),
        }
    )
    for label in ("all_missing_tip", "unidentified_response", "singular_covariance"):
        responses, design, components, options = fixture(16, 2, 0.5)
        if label == "all_missing_tip":
            responses[0] = np.nan
        elif label == "unidentified_response":
            responses[:, 1] = np.nan
            responses[0, 1] = 1.0
        else:
            components = {"tree": np.zeros((16, 16))}
        outcomes = []
        for function in (module.fit_multivariate_pgls, candidate):
            try:
                function(responses, design, components, **options)
                outcomes.append({"exception": None})
            except Exception as exc:
                outcomes.append({"exception": type(exc).__name__, "message": str(exc)})
        assert outcomes[0] == outcomes[1], (label, outcomes)
        records.append({"case": label, "outcomes": outcomes})
    return records


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tests", action="store_true")
    parser.add_argument("--additional", action="store_true")
    args = parser.parse_args()
    LOCK.parent.mkdir(parents=True, exist_ok=True)
    with LOCK.open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if args.tests:
            command = [
                sys.executable,
                "tools/check.py",
                "test",
                "--",
                "tests/test_numerical_invariance.py",
                "tests/test_multivariate_numerics.py",
                "tests/test_ordinary_regression.py::test_sparse_multivariate_pgls_matches_dense_fit",
                "tests/test_ordinary_regression.py::test_sparse_multivariate_pgls_matches_dense_fixed_factor",
                "tests/test_ordinary_regression.py::test_sparse_multivariate_pgls_matches_dense_with_missing_response",
                "tests/test_ordinary_regression.py::test_multivariate_gaussian_pgls_retains_partially_observed_tips",
                "tests/test_ordinary_regression.py::test_multivariate_cli_combines_biological_replicates_with_partial_missingness",
                "-rs",
            ]
            environment = os.environ.copy()
            environment.update(
                PYTHONPATH=str(HERE) + os.pathsep + str(ROOT),
                PYTEST_PLUGINS="candidate10_plugin",
                NWKIT_CI_DIAGNOSTICS="0",
            )
            result = subprocess.run(
                command, cwd=ROOT, env=environment, capture_output=True, text=True
            )
            (HERE / "tests.log").write_text(
                " ".join(command) + "\n" + result.stdout + result.stderr
            )
            print(result.stdout, end="", flush=True)
            print(result.stderr, end="", flush=True)
            raise SystemExit(result.returncode)

        import numpy as np
        import scipy

        from nwkit import __version__

        if args.additional:
            output = json.loads((HERE / "results.json").read_text())
            print("high missing full fit", flush=True)
            output["additional_high_missing_fit"] = full_case(128, 4, 0.95)
            output["edge_checks"] = edge_checks()
            (HERE / "results.json").write_text(json.dumps(output, indent=2) + "\n")
            print("additional completed", flush=True)
            return

        output = {
            "baseline_commit": "d7116f80b9979739c69a5431a5e34b822795408b",
            "environment": {
                "platform": platform.platform(),
                "python": platform.python_version(),
                "numpy": np.__version__,
                "scipy": scipy.__version__,
                "nwkit": __version__,
                "threads": 1,
            },
            "full_fits": [],
            "covariance_kernels": [],
        }
        for case in ((32, 2, 0.0), (64, 2, 0.5), (64, 3, 0.75), (192, 2, 0.75)):
            print("full", case, flush=True)
            output["full_fits"].append(full_case(*case))
            (HERE / "results.json").write_text(json.dumps(output, indent=2) + "\n")
        print("multiple/fixed", flush=True)
        output["full_fits"].append(
            full_case(24, 2, 0.5, fixed=True, multiple=True, reml=False)
        )
        for case in ((256, 4, 0.5), (256, 4, 0.95), (512, 4, 0.95)):
            print("kernel", case, flush=True)
            output["covariance_kernels"].append(kernel_case(*case))
        output["edge_checks"] = edge_checks()
        (HERE / "results.json").write_text(json.dumps(output, indent=2) + "\n")
        print("completed", flush=True)


if __name__ == "__main__":
    main()
