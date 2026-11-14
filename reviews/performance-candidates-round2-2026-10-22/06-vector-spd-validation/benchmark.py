"""Discard-free branch covariance validation for vector Gaussian pruning."""

import argparse
import hashlib
import inspect
import json
import os
import subprocess
import sys
import textwrap
from contextlib import contextmanager
from dataclasses import fields, replace
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = Path("/tmp/nwkit-perf-round2-20261022-rkzoGe/baseline")
for variable in (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "MKL_NUM_THREADS",
):
    os.environ[variable] = "1"
os.environ["NWKIT_CI_DIAGNOSTICS"] = "0"
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(SOURCE))

import numpy as np  # noqa: E402
from support import environment, measurement_lock, paired, peak_bytes  # noqa: E402

import nwkit  # noqa: E402
from nwkit import vector_gaussian as vg  # noqa: E402
from nwkit.compiled_tree import CompiledTree  # noqa: E402
from nwkit.util import read_tree  # noqa: E402

ORIGINAL = vg._propagate


def validate_branch_covariance(matrix):
    """Retain the original checks without constructing an unused inverse."""
    label = "Branch covariance"
    if not np.allclose(matrix, matrix.T, rtol=1e-12, atol=1e-14):
        raise ValueError(f"{label} must be symmetric.")
    try:
        factor = vg.cho_factor(matrix, lower=True, check_finite=False)
    except np.linalg.LinAlgError as exc:
        raise ValueError(f"{label} must be positive definite.") from exc
    # cho_solve originally performs this finite-factor check. The earlier
    # _array check still validates the complete input before factorization.
    np.asarray_chkfinite(factor[0])


def build_candidate():
    source = textwrap.dedent(inspect.getsource(ORIGINAL))
    old = '    _inverse(covariance, "Branch covariance")\n'
    assert source.count(old) == 1
    source = source.replace(old, "    _validate_branch_covariance(covariance)\n")
    namespace = dict(
        vg.__dict__, _validate_branch_covariance=validate_branch_covariance
    )
    exec(compile(source, "<candidate06_propagate>", "exec"), namespace)
    return namespace["_propagate"]


CANDIDATE = build_candidate()


@contextmanager
def patched(function):
    previous = vg._propagate
    vg._propagate = function
    try:
        yield
    finally:
        vg._propagate = previous


def fixture(tips, dimension, *, mode="partial", shared=False, flat=False, zero=False):
    rng = np.random.default_rng(600611 + tips + dimension)
    parts = [f"T{i}:{0.3 + (i % 7) * 0.1}" for i in range(tips)]
    while len(parts) > 1:
        parts = [
            "(" + ",".join(parts[i : i + 2]) + "):0.4"
            if i + 1 < len(parts)
            else parts[i]
            for i in range(0, len(parts), 2)
        ]
    tree = read_tree(parts[0] + ";", "1", True, quiet=True, rooted="yes")
    compiled = CompiledTree.from_tree(tree)
    raw = rng.normal(size=(dimension, dimension)) / np.sqrt(dimension)
    sigma = raw @ raw.T + np.eye(dimension) * 0.4
    slope = np.eye(dimension) * 0.85
    if dimension > 1:
        slope[0, 1] = 0.05
    intercept = rng.normal(scale=0.1, size=dimension)
    common = vg.VectorTransition(slope, intercept, sigma)
    transitions = {
        node: (
            vg.VectorTransition(
                np.eye(dimension), np.zeros(dimension), np.zeros_like(sigma)
            )
            if zero
            else common
            if shared
            else vg.VectorTransition(slope, intercept, sigma * float(node.dist))
        )
        for node in compiled.nodes[1:]
    }
    process = vg.VectorProcess(
        tree,
        dimension,
        transitions,
        None if flat else np.linspace(-0.2, 0.2, dimension),
        None if flat else sigma * 1.2,
    )
    observed, errors = {}, {}
    for i in range(tips):
        values = rng.normal(size=dimension)
        if zero:
            values[:] = np.linspace(-0.5, 0.5, dimension)
        observed[f"T{i}"] = [
            None if mode == "partial" and (i + j) % 5 == 0 else float(value)
            for j, value in enumerate(values)
        ]
        if mode != "exact":
            if mode == "partial":
                scale = np.array(
                    [0.0 if (i + j) % 4 == 0 else 0.1 for j in range(dimension)]
                )
            else:
                scale = np.full(dimension, 0.1)
            errors[f"T{i}"] = np.diag(scale)
    return process, observed, errors if mode != "exact" else None


def call(function, inputs):
    process, observed, errors = inputs
    with patched(function):
        return vg.condition_vector_tree(process, observed, error_covariances=errors)


def compare_results(before, after, *, draws=False):
    result = {}
    for field in fields(before):
        name = field.name
        left, right = getattr(before, name), getattr(after, name)
        if name in {"nodes", "parents"}:
            assert left == right, name
        else:
            np.testing.assert_array_equal(left, right, err_msg=name)
        result[name] = "exact"
    if draws:
        np.testing.assert_array_equal(
            before.sample(19, seed=617), after.sample(19, seed=617)
        )
        result["seeded_joint_draws"] = "exact"
    return result


def measure(name, inputs):
    def baseline():
        return call(ORIGINAL, inputs)

    def candidate():
        return call(CANDIDATE, inputs)

    process = inputs[0]
    return {
        "case": name,
        "tips": len(inputs[1]),
        "traits": process.dimension,
        "scope": "complete condition_vector_tree including compiled-tree, observation, pruning and smoothing setup",
        "comparison": compare_results(baseline(), candidate(), draws=True),
        "timing": paired(baseline, candidate),
        "peak_bytes": {
            "baseline": peak_bytes(baseline),
            "candidate": peak_bytes(candidate),
        },
    }


def exception_result(function, inputs):
    try:
        result = call(function, inputs)
    except Exception as exc:
        return {"type": type(exc).__name__, "message": str(exc)}
    return {"type": None, "finite": bool(np.isfinite(result.log_likelihood))}


def edge_checks():
    checks = []
    for dimension in (1, 2, 8):
        for mode in ("exact", "noisy", "partial"):
            for flat in (False, True):
                inputs = fixture(8, dimension, mode=mode, flat=flat)
                compare_results(
                    call(ORIGINAL, inputs), call(CANDIDATE, inputs), draws=True
                )
                checks.append(
                    {"case": f"d{dimension}-{mode}-flat{flat}", "exact": True}
                )
    inputs = fixture(4, 2, mode="noisy")
    process, observed, errors = inputs
    node = next(iter(process.transitions))
    original_transition = process.transitions[node]
    invalid = {
        "asymmetric": np.array([[1.0, 0.2], [0.1, 1.0]]),
        "indefinite": np.array([[1.0, 2.0], [2.0, 1.0]]),
        "negative": -np.eye(2),
        "singular_nonzero": np.ones((2, 2)),
        "nan_upper": np.array([[1.0, np.nan], [0.0, 1.0]]),
        "inf_diagonal": np.diag([np.inf, 1.0]),
        "wrong_shape": np.ones((1, 1)),
    }
    for name, covariance in invalid.items():
        transitions = dict(process.transitions)
        transitions[node] = replace(original_transition, covariance=covariance)
        broken = replace(process, transitions=transitions), observed, errors
        before, after = (
            exception_result(ORIGINAL, broken),
            exception_result(CANDIDATE, broken),
        )
        assert before == after, (name, before, after)
        assert before["type"] is not None
        checks.append({"case": name, "exception": before})
    # Validate tiny/huge units, nearly singular SPD, and accepted symmetry noise.
    for name, covariance in {
        "tiny_spd": np.diag([1e-120, 2e-120]),
        "huge_spd": np.diag([1e120, 2e120]),
        "nearly_singular_spd": np.array([[1.0, 1.0 - 1e-12], [1.0 - 1e-12, 1.0]]),
        "within_symmetry_tolerance": np.array([[1.0, 0.2 + 1e-14], [0.2, 1.0]]),
    }.items():
        transitions = dict(process.transitions)
        transitions[node] = replace(original_transition, covariance=covariance)
        altered = replace(process, transitions=transitions), observed, errors
        before, after = (
            exception_result(ORIGINAL, altered),
            exception_result(CANDIDATE, altered),
        )
        assert before == after, (name, before, after)
        if before["type"] is None:
            compare_results(call(ORIGINAL, altered), call(CANDIDATE, altered))
        checks.append({"case": name, "result": before, "exact": True})
    return checks


def run_tests():
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join((str(SOURCE), str(HERE)))
    command = [
        sys.executable,
        "tools/check.py",
        "test",
        "--",
        "-p",
        "candidate06_plugin",
        "tests/test_vector_gaussian.py",
        "tests/test_vector_fit.py",
        "tests/test_asr_multivariate_diagnostics.py",
        "-rs",
    ]
    with measurement_lock(), (HERE / "tests.log").open("w") as log:
        subprocess.run(
            command,
            cwd=SOURCE,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=True,
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tests", action="store_true")
    args = parser.parse_args()
    if args.tests:
        run_tests()
        return
    assert Path(nwkit.__file__).resolve().parent.parent == SOURCE.resolve()
    with measurement_lock():
        results = {
            "environment": environment(),
            "source_sha256": hashlib.sha256(
                inspect.getsource(ORIGINAL).encode()
            ).hexdigest(),
            "source_scope": "immutable baseline vector_gaussian._propagate",
            "cases": [
                measure(
                    "8tips-2traits-exact-shared",
                    fixture(8, 2, mode="exact", shared=True),
                ),
                measure("64tips-4traits-partial-unique", fixture(64, 4)),
                measure(
                    "64tips-16traits-noisy-shared",
                    fixture(64, 16, mode="noisy", shared=True),
                ),
                measure("128tips-32traits-partial-unique", fixture(128, 32)),
                measure(
                    "32tips-64traits-noisy-flat",
                    fixture(32, 64, mode="noisy", flat=True),
                ),
                measure(
                    "64tips-4traits-zero-branch-control",
                    fixture(64, 4, mode="noisy", zero=True),
                ),
            ],
            "edge_checks": edge_checks(),
        }
        (HERE / "results.json").write_text(
            json.dumps(results, indent=2, allow_nan=False) + "\n"
        )
        for case in results["cases"]:
            print(
                case["case"], case["timing"]["speedup"], case["peak_bytes"], flush=True
            )
        print("edge_checks", len(results["edge_checks"]), flush=True)


if __name__ == "__main__":
    main()
