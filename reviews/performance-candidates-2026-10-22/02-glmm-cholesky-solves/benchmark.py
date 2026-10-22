"""Isolated dense GLMM Cholesky-solve experiment; never edits nwkit sources."""

import argparse
import ast
import fcntl
import inspect
import json
import math
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
NAMES = (
    "_random_mode",
    "_scalar_random_mode",
    "_polish_scalar_random_mode",
    "_marginal_coefficient_covariance",
    "_marginal_scalar_coefficient_covariance",
)


class ReplaceCholeskySolve(ast.NodeTransformer):
    def visit_Call(self, node):
        if (
            isinstance(node.func, ast.Attribute)
            and ast.unparse(node.func) == "np.linalg.solve"
            and len(node.args) == 2
            and isinstance(node.args[0], ast.Attribute)
            and node.args[0].attr == "T"
            and isinstance(node.args[1], ast.Call)
            and ast.unparse(node.args[1].func) == "np.linalg.solve"
            and len(node.args[1].args) == 2
            and ast.dump(node.args[0].value) == ast.dump(node.args[1].args[0])
        ):
            return ast.copy_location(
                ast.Call(
                    func=ast.Name(id="_candidate02_cho_solve", ctx=ast.Load()),
                    args=[
                        ast.Tuple(
                            elts=[node.args[0].value, ast.Constant(value=True)],
                            ctx=ast.Load(),
                        ),
                        node.args[1].args[1],
                    ],
                    keywords=[],
                ),
                node,
            )
        return self.generic_visit(node)


@contextmanager
def prototype():
    from scipy.linalg import cho_solve

    import nwkit.phylogenetic_glmm as glmm

    originals = {name: getattr(glmm, name) for name in NAMES}
    old_alias = glmm.__dict__.get("_candidate02_cho_solve")
    glmm._candidate02_cho_solve = cho_solve
    for name, function in originals.items():
        source = inspect.getsource(function)
        tree = ast.fix_missing_locations(
            ReplaceCholeskySolve().visit(ast.parse(source))
        )
        assert ast.unparse(tree) != ast.unparse(ast.parse(source)), name
        namespace = {}
        exec(
            compile(tree, str(HERE / "benchmark.py"), "exec"), glmm.__dict__, namespace
        )
        setattr(glmm, name, namespace[name])
    try:
        yield
    finally:
        for name, function in originals.items():
            setattr(glmm, name, function)
        if old_alias is None:
            del glmm._candidate02_cho_solve
        else:
            glmm._candidate02_cho_solve = old_alias


def compare_outputs(first, second):
    import numpy as np

    if hasattr(first, "__dataclass_fields__"):
        numeric_fields = [
            "coefficients",
            "coefficient_covariance",
            "random_modes",
            "component_variances",
            "log_likelihood",
            "dispersion",
            "thresholds",
            "threshold_covariance",
            "objective_value",
        ]
        differences = {}
        for name in numeric_fields:
            if hasattr(first, name):
                differences[name] = compare_outputs(
                    getattr(first, name), getattr(second, name)
                )
        assert first.optimizer_converged == second.optimizer_converged
        assert first.boundary_warning == second.boundary_warning
        return differences
    if isinstance(first, dict):
        assert first.keys() == second.keys()
        return {key: compare_outputs(first[key], second[key]) for key in first}
    if isinstance(first, tuple):
        return [compare_outputs(a, b) for a, b in zip(first, second, strict=True)]
    if first is None:
        assert second is None
        return None
    first, second = np.asarray(first), np.asarray(second)
    difference = float(np.max(np.abs(first - second), initial=0))
    np.testing.assert_allclose(first, second, rtol=2e-4, atol=2e-6, equal_nan=True)
    return {"max_abs_difference": difference}


def measure(function, repetitions=7):
    outputs = []
    warmup_times = []
    for candidate in (False, True):
        if candidate:
            with prototype():
                start = time.perf_counter()
                outputs.append(function())
                warmup_times.append(time.perf_counter() - start)
        else:
            start = time.perf_counter()
            outputs.append(function())
            warmup_times.append(time.perf_counter() - start)
    equivalence = compare_outputs(*outputs)
    batch = min(500, max(1, math.ceil(0.05 / min(warmup_times))))
    timings = {"baseline": [], "candidate": []}
    for repetition in range(repetitions):
        for candidate in (False, True) if repetition % 2 == 0 else (True, False):
            label = "candidate" if candidate else "baseline"
            if candidate:
                with prototype():
                    start = time.perf_counter()
                    for _ in range(batch):
                        result = function()
                    elapsed = (time.perf_counter() - start) / batch
            else:
                start = time.perf_counter()
                for _ in range(batch):
                    result = function()
                elapsed = (time.perf_counter() - start) / batch
            compare_outputs(outputs[int(candidate)], result)
            timings[label].append(elapsed)
    peak = {}
    for candidate in (False, True):
        context = prototype() if candidate else contextmanager(lambda: (yield))()
        with context:
            tracemalloc.start()
            function()
            peak["candidate" if candidate else "baseline"] = (
                tracemalloc.get_traced_memory()[1]
            )
            tracemalloc.stop()
    medians = {key: statistics.median(values) for key, values in timings.items()}
    return {
        "seconds": timings,
        "calls_per_timed_trial": batch,
        "median_seconds": medians,
        "range_seconds": {
            key: [min(values), max(values)] for key, values in timings.items()
        },
        "speedup": medians["baseline"] / medians["candidate"],
        "peak_python_bytes": peak,
        "equivalence": equivalence,
    }


def workloads():
    import numpy as np
    from scipy.linalg import toeplitz

    import nwkit.phylogenetic_glmm as glmm

    results = {}
    for n in (64, 256, 512):
        rng = np.random.default_rng(2202 + n)
        covariance = toeplitz(0.8 ** np.arange(n)) + 0.2 * np.eye(n)
        precision = np.linalg.inv(covariance)
        design = np.column_stack([np.ones(n), rng.normal(size=n)])
        values = rng.poisson(np.exp(0.2 + 0.4 * design[:, 1]))
        kwargs = dict(
            family="poisson",
            dispersion=None,
            zero_probability=None,
            offset=np.zeros(n),
            trials=None,
            censor_lower=None,
            censor_upper=None,
        )
        results[f"poisson_mode_{n}"] = measure(
            lambda values=values, design=design, precision=precision, n=n, kwargs=kwargs: (
                glmm._scalar_random_mode(
                    values,
                    design @ np.array([0.2, 0.4]),
                    precision,
                    np.arange(n),
                    **kwargs,
                )
            )
        )
        weights = np.exp(0.2 + design[:, 1] * 0.4)
        results[f"scalar_marginal_covariance_{n}"] = measure(
            lambda design=design, n=n, weights=weights, precision=precision: (
                glmm._marginal_scalar_coefficient_covariance(
                    design, np.arange(n), weights, precision
                )
            )
        )
        counts = rng.multinomial(1, [0.3, 0.4, 0.3], size=n)
        linear = np.zeros((n, 2))
        categorical_precision = np.kron(precision, np.eye(2))
        results[f"multinomial_mode_{n}"] = measure(
            lambda counts=counts, linear=linear, categorical_precision=categorical_precision: (
                glmm._random_mode(
                    counts,
                    linear,
                    categorical_precision,
                    family="multinomial",
                    thresholds=np.empty(0),
                )
            )
        )
    for n in (24, 64):
        rng = np.random.default_rng(2203 + n)
        covariance = toeplitz(0.8 ** np.arange(n)) + 0.2 * np.eye(n)
        design = np.column_stack([np.ones(n), rng.normal(size=n)])
        values = rng.poisson(np.exp(0.2 + 0.4 * design[:, 1]))
        results[f"poisson_complete_fit_{n}"] = measure(
            lambda values=values, design=design, covariance=covariance: (
                glmm.fit_phylogenetic_glmm(
                    values,
                    design,
                    lambda parameter, covariance=covariance: covariance,
                    family="poisson",
                )
            )
        )
    return results


def edge_cases():
    import numpy as np
    from scipy.linalg import cho_solve

    import nwkit.phylogenetic_glmm as glmm

    rng = np.random.default_rng(2204)
    records = []
    for exponent in (0, 4, 8, 12):
        orthogonal, _ = np.linalg.qr(rng.normal(size=(32, 32)))
        matrix = (orthogonal * np.geomspace(1.0, 10.0**exponent, 32)) @ orthogonal.T
        factor = np.linalg.cholesky(matrix)
        rhs = rng.normal(size=(32, 3))
        first = np.linalg.solve(factor.T, np.linalg.solve(factor, rhs))
        second = cho_solve((factor, True), rhs)
        records.append(
            {
                "condition_number": float(np.linalg.cond(matrix)),
                "relative_solution_difference": float(
                    np.linalg.norm(first - second) / np.linalg.norm(first)
                ),
                "candidate_relative_residual": float(
                    np.linalg.norm(matrix @ second - rhs) / np.linalg.norm(rhs)
                ),
            }
        )
    invalid = []
    for matrix in (np.array([[1.0, 2.0], [2.0, 1.0]]), -np.eye(2)):
        errors = []
        for candidate in (False, True):
            context = prototype() if candidate else contextmanager(lambda: (yield))()
            with context:
                try:
                    glmm._marginal_scalar_coefficient_covariance(
                        np.eye(2), np.arange(2), np.zeros(2), matrix
                    )
                except Exception as exc:
                    errors.append(type(exc).__name__)
        assert errors == ["LinAlgError", "LinAlgError"]
        invalid.append(errors)
    nonfinite = []
    for value in (np.nan, np.inf):
        errors = []
        for candidate in (False, True):
            context = prototype() if candidate else contextmanager(lambda: (yield))()
            with context:
                try:
                    glmm._marginal_scalar_coefficient_covariance(
                        np.eye(2), np.arange(2), np.ones(2), np.eye(2) * value
                    )
                except Exception as exc:
                    errors.append(type(exc).__name__)
                else:
                    errors.append("returned")
        nonfinite.append({"precision_value": str(value), "baseline_candidate": errors})
    return {
        "condition_cases": records,
        "invalid_precision_exceptions": invalid,
        "nonfinite_internal_helper_behavior": nonfinite,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=HERE / "results.json")
    args = parser.parse_args()
    for name in (
        "OPENBLAS_NUM_THREADS",
        "OMP_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "MKL_NUM_THREADS",
    ):
        assert os.environ.get(name) == "1", name
    with LOCK.open("a") as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        import ete4
        import numpy as np
        import pandas as pd
        import scipy

        results = {
            "environment": {
                "python": sys.version,
                "platform": platform.platform(),
                "numpy": np.__version__,
                "scipy": scipy.__version__,
                "pandas": pd.__version__,
                "ete4": ete4.__version__,
                "threads": {
                    name: os.environ[name]
                    for name in (
                        "OPENBLAS_NUM_THREADS",
                        "OMP_NUM_THREADS",
                        "VECLIB_MAXIMUM_THREADS",
                        "MKL_NUM_THREADS",
                    )
                },
            },
            "workloads": workloads(),
            "edge_cases": edge_cases(),
        }
        args.output.write_text(json.dumps(results, indent=2) + "\n")
        print(
            json.dumps(
                {
                    name: result["speedup"]
                    for name, result in results["workloads"].items()
                },
                indent=2,
            )
        )


if __name__ == "__main__":
    main()
