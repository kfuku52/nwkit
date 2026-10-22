"""Reuse one consecutive full-OU distance-pair product without reassociation."""

import argparse
import hashlib
import inspect
import json
import os
import subprocess
import sys
from contextlib import contextmanager
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
from nwkit import full_ou as module  # noqa: E402
from nwkit.multivariate_gaussian_asr import _prepare_observations  # noqa: E402
from nwkit.util import read_tree  # noqa: E402

ORIGINAL = module._moment_vector


def consecutive_moments(parameters, dimension, time_scale, design):
    attraction, _, covariance = module.decode_full_ou(parameters, dimension, time_scale)
    distances = {distance for row in design for distance in row[:2]}
    slopes = {distance: module.expm(-attraction * distance) for distance in distances}
    output = []
    last_pair = None
    product = None
    for left, right, first, second in design:
        pair = (left, right)
        if pair != last_pair:
            product = slopes[left] @ covariance @ slopes[right].T
            last_pair = pair
        output.append(product[first, second])
    return np.array(output)


@contextmanager
def patched(function):
    previous = module._moment_vector
    module._moment_vector = function
    try:
        yield
    finally:
        module._moment_vector = previous


def identify(parameters, dimension, time_scale, design, complete, function):
    with patched(function):
        return module.full_ou_identifiability(
            parameters, dimension, time_scale, design, complete
        )


def jacobian(parameters, dimension, time_scale, design, function):
    columns = []
    for index, value in enumerate(parameters):
        step = 1e-4 * max(1.0, abs(float(value)))
        delta = np.zeros(len(parameters))
        delta[index] = step
        columns.append(
            (
                function(parameters + delta, dimension, time_scale, design)
                - function(parameters - delta, dimension, time_scale, design)
            )
            / (2 * step)
        )
    return np.column_stack(columns)


def compare(parameters, dimension, time_scale, design, complete):
    baseline = ORIGINAL(parameters, dimension, time_scale, design)
    candidate = consecutive_moments(parameters, dimension, time_scale, design)
    np.testing.assert_array_equal(baseline, candidate)
    before = jacobian(parameters, dimension, time_scale, design, ORIGINAL)
    after = jacobian(parameters, dimension, time_scale, design, consecutive_moments)
    np.testing.assert_array_equal(before, after)
    first = identify(parameters, dimension, time_scale, design, complete, ORIGINAL)
    second = identify(
        parameters, dimension, time_scale, design, complete, consecutive_moments
    )
    assert first == second, (first, second)
    return {
        "all_moments_exact": True,
        "all_jacobian_entries_exact": True,
        "all_diagnostic_fields_exact": True,
        "maximum_moment_difference": float(
            np.max(np.abs(baseline - candidate), initial=0)
        ),
        "status": first[0],
        "rank": first[1],
        "singular_ratio": first[2],
    }


def parameters_for(dimension, seed=1010, *, initial=False):
    values, _ = module.full_ou_initial(dimension)
    values = np.asarray(values)
    if not initial:
        values = values + np.random.default_rng(seed).normal(0, 0.12, len(values))
    return values


def native_design(dimension, tips, *, equal=False, gap=None, limit=4096):
    if gap is not None:
        lengths = [1 + gap * i for i in range(tips)]
    elif equal:
        lengths = [1.0] * tips
    else:
        lengths = [0.25 + (i % 8) * 0.125 for i in range(tips)]
    leaves = [f"T{i}:{length:.16g}" for i, length in enumerate(lengths)]
    tree = read_tree(
        "(" + ",".join(leaves) + ")R;", "1", True, quiet=True, rooted="yes"
    )
    values = {
        f"T{i}": list(np.random.default_rng(710 + i).normal(size=dimension))
        for i in range(tips)
    }
    data = _prepare_observations(
        tree, values, tuple(f"trait{i}" for i in range(dimension)), None
    )
    return module.observation_moment_design(data, limit=limit)


def pair_runs(design):
    previous = None
    count = 0
    for left, right, _, _ in design:
        current = (left, right)
        count += current != previous
        previous = current
    return count


def measure_case(name, parameters, dimension, design, complete, *, time_scale=1.0):
    def first():
        return identify(parameters, dimension, time_scale, design, complete, ORIGINAL)

    def second():
        return identify(
            parameters, dimension, time_scale, design, complete, consecutive_moments
        )

    comparison = compare(parameters, dimension, time_scale, design, complete)
    timing = paired(first, second)
    before_peak, after_peak = peak_bytes(first), peak_bytes(second)
    print(name, f"{timing['speedup']:.3f}x", comparison, flush=True)
    return {
        "case": name,
        "scope": "complete full_ou_identifiability including all finite differences, factors, expm and SVD",
        "dimension": dimension,
        "parameters": len(parameters),
        "design_entries": len(design),
        "design_complete": complete,
        "consecutive_distance_pairs": pair_runs(design),
        "unique_distances": len({value for row in design for value in row[:2]}),
        "time_scale": time_scale,
        "comparison": comparison,
        "timing": timing,
        "traced_peak_bytes": {"baseline": before_peak, "candidate": after_peak},
    }


def edge_checks():
    results = []
    for gap in (0.0, 1e-9, 1e-7, 1e-5, 1e-3):
        design, complete = native_design(2, 8, gap=gap)
        comparison = compare(parameters_for(2, initial=True), 2, 1, design, complete)
        results.append({"case": "near_ultrametric", "gap": gap, **comparison})
    for dimension in (2, 4, 8):
        design, complete = native_design(dimension, 12)
        parameters = parameters_for(dimension)
        for scale in (1e-5, 1e5):
            scaled = tuple(
                (left * scale, right * scale, first, second)
                for left, right, first, second in design
            )
            comparison = compare(parameters, dimension, scale, scaled, complete)
            results.append(
                {
                    "case": "time_units",
                    "dimension": dimension,
                    "scale": scale,
                    **comparison,
                }
            )
        shuffled = list(design)
        np.random.default_rng(8010).shuffle(shuffled)
        shuffled = tuple(shuffled + shuffled[:5])
        comparison = compare(parameters, dimension, 1, shuffled, complete)
        results.append(
            {"case": "shuffled_duplicate_design", "dimension": dimension, **comparison}
        )
    design, complete = native_design(2, 8, limit=2)
    results.append(
        {
            "case": "incomplete_subset",
            **compare(parameters_for(2), 2, 1, design, complete),
        }
    )
    # A nonnormal attraction with complex eigenvalues and a negative diagonal.
    from scipy.linalg import solve_continuous_lyapunov

    attraction = np.array([[-1.0, 10.0], [-1.0, 3.0]])
    diffusion = np.eye(2)
    covariance = solve_continuous_lyapunov(attraction, diffusion)
    parameters = []
    for matrix in (covariance, diffusion):
        lower = np.linalg.cholesky(matrix)
        parameters += [np.log(lower[0, 0]), lower[1, 0], np.log(lower[1, 1])]
    parameters += [(attraction @ covariance - diffusion / 2)[1, 0]]
    design, complete = native_design(2, 8)
    results.append(
        {
            "case": "nonnormal_complex",
            **compare(np.asarray(parameters), 2, 1, design, complete),
        }
    )
    empty_before = ORIGINAL(parameters_for(2), 2, 1, ())
    empty_after = consecutive_moments(parameters_for(2), 2, 1, ())
    np.testing.assert_array_equal(empty_before, empty_after)
    results.append({"case": "empty_moment_design", "exact": True})
    return results


def run_measurements():
    cases = []
    for dimension, tips in ((2, 12), (4, 16), (8, 32)):
        design, complete = native_design(dimension, tips)
        cases.append(
            measure_case(
                f"native_star_d{dimension}_n{tips}",
                parameters_for(dimension),
                dimension,
                design,
                complete,
            )
        )
    design, complete = native_design(2, 8, equal=True)
    cases.append(
        measure_case(
            "ultrametric_deficient_d2",
            parameters_for(2, initial=True),
            2,
            design,
            complete,
        )
    )
    design, complete = native_design(4, 16)
    shuffled = list(design)
    np.random.default_rng(1037).shuffle(shuffled)
    cases.append(
        measure_case(
            "shuffled_external_d4", parameters_for(4), 4, tuple(shuffled), complete
        )
    )
    unique = tuple((0.5, 0.2 + i / 127, i % 4, (i // 4) % 4) for i in range(64))
    cases.append(
        measure_case("unique_distance_pairs_d4", parameters_for(4), 4, unique, False)
    )
    design, complete = native_design(2, 8, limit=2)
    cases.append(
        measure_case("incomplete_subset_d2", parameters_for(2), 2, design, complete)
    )
    output = {
        "environment": environment(),
        "source_sha256": hashlib.sha256(
            inspect.getsource(ORIGINAL).encode()
        ).hexdigest(),
        "source_path": module.__file__,
        "strategy": "Reuse only the immediately preceding identical distance pair's complete two-matmul product; no reassociation or long-lived cache",
        "cases": cases,
        "edge_checks": edge_checks(),
    }
    (HERE / "results.json").write_text(
        json.dumps(output, indent=2, allow_nan=False) + "\n"
    )


def run_tests():
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join((str(SOURCE), str(HERE), str(HERE.parent)))
    command = [
        sys.executable,
        "tools/check.py",
        "test",
        "--",
        "-p",
        "candidate10_plugin",
        "tests/test_full_ou.py",
        "-rs",
    ]
    with (HERE / "tests.log").open("w") as handle:
        print("command:", " ".join(command), file=handle, flush=True)
        subprocess.run(
            command,
            cwd=SOURCE,
            env=env,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=True,
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tests", action="store_true")
    args = parser.parse_args()
    assert Path(nwkit.__file__).resolve().parent.parent == SOURCE.resolve()
    with measurement_lock():
        if args.tests:
            run_tests()
        else:
            run_measurements()


if __name__ == "__main__":
    main()
