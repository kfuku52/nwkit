"""Reuse final PCA factor/whitening while retaining the likelihood SVD."""

import argparse
import hashlib
import inspect
import json
import os
import subprocess
import sys
import textwrap
from contextlib import contextmanager
from dataclasses import fields
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
os.environ["DYLD_FALLBACK_LIBRARY_PATH"] = "/Users/kf/miniforge3-arm64/lib"
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(SOURCE))

import numpy as np  # noqa: E402
from support import environment, measurement_lock, paired, peak_bytes  # noqa: E402

import nwkit  # noqa: E402
from nwkit import phylogenetic_pca as module  # noqa: E402

ORIGINAL = module.fit_pca
ORIGINAL_PROFILE = module._profile_likelihood


def build_candidate():
    profile = textwrap.dedent(inspect.getsource(ORIGINAL_PROFILE))
    old_head = """def _profile_likelihood(covariance, values):
    n, p = values.shape
    factor, _, residual = _centered(covariance, values)
    whitened = solve_triangular(factor, residual, lower=True)
"""
    new_head = "def _likelihood_from_white(factor, whitened, n, p):\n"
    assert profile.count(old_head) == 1
    profile = profile.replace(old_head, new_head)
    namespace = dict(module.__dict__)
    exec(compile(profile, "<candidate07_profile>", "exec"), namespace)
    source = textwrap.dedent(inspect.getsource(ORIGINAL))
    old_white = "    white = white / scale\n"
    new_white = "    unscaled_white = white\n" + old_white
    old_call = "float(_profile_likelihood(covariance, values))"
    new_call = "float(_likelihood_from_white(factor, unscaled_white, n, p))"
    assert source.count(old_white) == source.count(old_call) == 1
    source = source.replace(old_white, new_white).replace(old_call, new_call)
    exec(compile(source, "<candidate07_fit>", "exec"), namespace)
    return namespace["fit_pca"]


CANDIDATE = build_candidate()


@contextmanager
def prototype():
    from nwkit import pca

    previous = module.fit_pca
    consumer = pca.fit_pca
    module.fit_pca = CANDIDATE
    pca.fit_pca = CANDIDATE
    try:
        yield
    finally:
        module.fit_pca = previous
        pca.fit_pca = consumer


def fixture(n, p, *, seed=7125, nonultrametric=False):
    covariance = np.eye(n) * 0.75
    width = 2
    while width < n:
        for start in range(0, n, width):
            covariance[start : start + width, start : start + width] += 0.25
        width *= 2
    if nonultrametric:
        covariance[np.diag_indices(n)] += np.linspace(0.01, 0.9, n)
    rng = np.random.default_rng(seed + n + p)
    values = np.linalg.cholesky(covariance) @ rng.normal(size=(n, p))
    values += np.linspace(-1, 2, p)
    return covariance, values


def compare(first, second):
    result = {}
    for field in fields(first):
        left, right = getattr(first, field.name), getattr(second, field.name)
        if isinstance(left, np.ndarray):
            np.testing.assert_array_equal(left, right, err_msg=field.name)
            result[field.name] = {
                "exact": True,
                "shape": list(left.shape),
                "max_abs_difference": float(np.max(np.abs(left - right), initial=0)),
            }
        else:
            assert left == right, (field.name, left, right)
            result[field.name] = {"exact": True, "value": left}
    return result


def outcome(function, covariance, values, **kwargs):
    try:
        return function(covariance, values, **kwargs)
    except Exception as error:
        return type(error).__name__, str(error)


def edge_checks():
    checks = []

    def check(name, covariance, values, **kwargs):
        before = outcome(ORIGINAL, covariance, values, **kwargs)
        after = outcome(CANDIDATE, covariance, values, **kwargs)
        if isinstance(before, tuple):
            assert before == after, (name, before, after)
            checks.append({"case": name, "error": list(before), "exact": True})
        else:
            compare(before, after)
            checks.append(
                {
                    "case": name,
                    "status": before.status,
                    "likelihood_is_none": before.log_likelihood is None,
                    "exact": True,
                }
            )

    covariance, values = fixture(16, 4)
    for mode in ("cov", "corr"):
        for model, fixed in (
            ("BM", None),
            ("LAMBDA", 0),
            ("LAMBDA", 0.4),
            ("LAMBDA", 1),
        ):
            check(
                f"{mode}-{model}-{fixed}",
                covariance,
                values,
                model=model,
                mode=mode,
                lambda_value=fixed,
            )
        check(
            f"rank-deficient-{mode}",
            covariance,
            np.column_stack((values[:, 0], 2 * values[:, 0])),
            mode=mode,
        )
        check(
            f"almost-dependent-{mode}",
            covariance,
            np.column_stack((values[:, 0], values[:, 0] + 1e-15 * values[:, 1])),
            mode=mode,
        )
        for exponent in (-200, 200):
            check(
                f"tree-units-{exponent}-{mode}",
                covariance * 2.0**exponent,
                values,
                mode=mode,
            )
        check(
            f"mixed-trait-units-{mode}",
            covariance,
            values * 2.0 ** np.array([100, -100, 40, -40]),
            mode=mode,
        )
        check(f"large-offset-{mode}", covariance, values + 1e10, mode=mode)
    check(
        "free-lambda-nonultrametric",
        *fixture(32, 3, nonultrametric=True),
        model="LAMBDA",
        mode="corr",
    )
    wide_c, wide_y = fixture(8, 24)
    check("wide-fixed", wide_c, wide_y, mode="corr")
    check("wide-free-rejected", wide_c, wide_y, model="LAMBDA")
    repeated = np.array([[1, 0], [-1, 0], [0, 1], [0, -1], [0, 0]], dtype=float)
    check("repeated-eigenvalues", np.eye(5), repeated)
    check("constant", covariance, np.ones_like(values))
    bad_values = values.copy()
    bad_values[0, 0] = np.nan
    check("nonfinite", covariance, bad_values)
    bad_covariance = covariance.copy()
    bad_covariance[0, 1] += 0.3
    check("asymmetric", bad_covariance, values)
    check("singular", np.ones_like(covariance), values)
    check("bad-mode", covariance, values, mode="bad")
    check("bad-model", covariance, values, model="bad")
    check("bad-fixed", covariance, values, model="BM", lambda_value=0.4)
    # Every fixed-rank likelihood is still checked by its original SVD.
    inverse = np.linalg.inv(covariance)
    fit = CANDIDATE(covariance, values)
    center = (np.ones(len(values)) @ inverse @ values) / inverse.sum()
    residual = values - center
    trait = residual.T @ inverse @ residual / len(values)
    direct = -0.5 * (
        values.size * (np.log(2 * np.pi) + 1)
        + values.shape[1] * np.linalg.slogdet(covariance)[1]
        + len(values) * np.linalg.slogdet(trait)[1]
    )
    difference = abs(fit.log_likelihood - direct)
    assert difference < 1e-7  # Existing independent matrix-normal test tolerance.
    return checks, {
        "absolute_log_likelihood_difference": difference,
        "existing_test_atol": 1e-7,
    }


def measure_case(name, n, p, **kwargs):
    covariance, values = fixture(n, p, nonultrametric=name.startswith("nonultra"))
    baseline = lambda: ORIGINAL(covariance, values, **kwargs)  # noqa: E731
    candidate = lambda: CANDIDATE(covariance, values, **kwargs)  # noqa: E731
    comparison = compare(baseline(), candidate())
    print(name, flush=True)
    timing = paired(baseline, candidate)
    return {
        "case": name,
        "tips": n,
        "traits": p,
        "kwargs": kwargs,
        "scope": "complete fit_pca including lambda search when requested; fixture/import excluded",
        "comparison": comparison,
        "timing": timing,
        "baseline_peak_bytes": peak_bytes(baseline),
        "candidate_peak_bytes": peak_bytes(candidate),
    }


def run_tests():
    environment_variables = os.environ.copy()
    environment_variables["PYTHONPATH"] = os.pathsep.join(
        [str(SOURCE), str(HERE), str(HERE.parent)]
    )
    command = [
        sys.executable,
        "tools/check.py",
        "test",
        "--",
        "-p",
        "candidate07_plugin",
        "tests/test_pca.py",
        "tests/test_numerical_invariance.py",
        "-rs",
    ]
    print(" ".join(command), flush=True)
    with (HERE / "tests.log").open("w") as handle:
        completed = subprocess.run(
            command,
            cwd=SOURCE,
            env=environment_variables,
            stdout=handle,
            stderr=subprocess.STDOUT,
        )
    return completed.returncode


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tests", action="store_true")
    args = parser.parse_args()
    assert Path(nwkit.__file__).resolve().parent.parent == SOURCE.resolve()
    with measurement_lock():
        if args.tests:
            return run_tests()
        cases = [
            measure_case("small-bm-cov", 8, 2),
            measure_case("bm-cov", 128, 4),
            measure_case("large-bm-cov", 512, 4),
            measure_case("large-bm-corr", 512, 16, mode="corr"),
            measure_case("many-trait-bm-corr", 128, 32, mode="corr"),
            measure_case("wide-bm-corr", 32, 64, mode="corr"),
            measure_case(
                "fixed-lambda-corr",
                128,
                4,
                model="LAMBDA",
                lambda_value=0.4,
                mode="corr",
            ),
            measure_case(
                "fixed-zero-lambda-cov", 512, 4, model="LAMBDA", lambda_value=0.0
            ),
            measure_case("estimated-lambda-corr", 64, 4, model="LAMBDA", mode="corr"),
            measure_case("nonultra-estimated-lambda-cov", 128, 4, model="LAMBDA"),
        ]
        edges, direct = edge_checks()
        result = {
            "environment": environment(),
            "target_source_sha256": hashlib.sha256(
                inspect.getsource(ORIGINAL).encode()
            ).hexdigest(),
            "profile_source_sha256": hashlib.sha256(
                inspect.getsource(ORIGINAL_PROFILE).encode()
            ).hexdigest(),
            "cases": cases,
            "edge_checks": edges,
            "independent_matrix_normal": direct,
        }
        (HERE / "results.json").write_text(
            json.dumps(result, indent=2, allow_nan=False) + "\n"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
