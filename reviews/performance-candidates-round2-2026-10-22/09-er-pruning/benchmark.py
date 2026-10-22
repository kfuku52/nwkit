"""Operation-local ER diagonal-plus-constant log pruning prototype."""

import argparse
import hashlib
import inspect
import itertools
import json
import math
import os
import subprocess
import sys
import textwrap
from contextlib import contextmanager
from functools import partial
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
from scipy.special import logsumexp  # noqa: E402
from support import environment, measurement_lock, paired, peak_bytes  # noqa: E402

from nwkit import asr  # noqa: E402
from nwkit.util import read_tree  # noqa: E402

ORIGINAL = asr._compute_log_inside_likelihoods
MEASURED_MINIMUM_STATES = 16
PROTOTYPE_CONFIGURATION = {
    "minimum_states": MEASURED_MINIMUM_STATES,
    "rate_dispatch": "Fixed and free rates, intentionally exploratory",
    "pytest_minimum_states": 2,
    "proposed_integration_scope": {
        "minimum_states": 32,
        "fixed_rate_only": True,
        "implemented_in_this_prototype": False,
    },
}


def er_descriptor(matrix, cache):
    key = id(matrix)
    if key in cache:
        return cache[key][1]
    description = None
    if (
        matrix.ndim == 2
        and matrix.shape[0] == matrix.shape[1]
        and len(matrix) >= 2
        and matrix.dtype.kind in "biuf"
    ):
        diagonal = np.diag(matrix)
        a, b = float(diagonal[0]), float(matrix[0, 1])
        if (
            math.isfinite(a)
            and math.isfinite(b)
            and 0.0 <= b <= a
            and np.all(diagonal == a)
            and np.count_nonzero(matrix != b) == (0 if a == b else len(matrix))
        ):
            delta = a - b
            description = (
                -math.inf if delta == 0 else math.log(delta),
                -math.inf if b == 0 else math.log(b),
            )
    # Keep the array alive throughout this call; identity cannot be recycled.
    cache[key] = matrix, description
    return description


def generic_log_probabilities(matrix, cache):
    return (
        None
        if er_descriptor(matrix, cache) is not None
        else asr._log_probabilities(matrix)
    )


def candidate_message(matrix, log_matrix, inside, cache):
    descriptor = er_descriptor(matrix, cache)
    if descriptor is None:
        if log_matrix is None:
            log_matrix = asr._log_probabilities(matrix)
        return logsumexp(log_matrix + inside[None, :], axis=1)
    log_delta, log_off_diagonal = descriptor
    return np.logaddexp(log_delta + inside, log_off_diagonal + logsumexp(inside))


def build_candidate(minimum_states=16):
    source = textwrap.dedent(inspect.getsource(ORIGINAL))
    before = "    transition_matrix_cache = dict()\n"
    source = source.replace(before, before + "    er_descriptor_cache = {}\n", 1)
    source = source.replace(
        "                log_matrix = _log_probabilities(matrix)\n",
        "                log_matrix = None\n",
        1,
    )
    source = source.replace(
        "                        _log_probabilities(matrix),\n",
        "                        _generic_log_probabilities(matrix, er_descriptor_cache),\n",
        1,
    )
    before = "            term = logsumexp(log_matrix + inside[child][None, :], axis=1)"
    after = "            term = _candidate_message(matrix, log_matrix, inside[child], er_descriptor_cache)"
    assert source.count(before) == 1
    source = source.replace(before, after)
    namespace = dict(asr.__dict__)
    namespace.update(
        _candidate_message=candidate_message,
        _generic_log_probabilities=generic_log_probabilities,
    )
    exec(compile(source, "<candidate09_inside>", "exec"), namespace)
    compiled = namespace["_compute_log_inside_likelihoods"]

    def wrapped(
        tree, likelihood_by_leaf, rate_matrix, *, fixed_transition_matrices=None
    ):
        sample = (
            next(iter(rate_matrix.values()))
            if isinstance(rate_matrix, dict)
            else rate_matrix
        )
        function = ORIGINAL if len(sample) < minimum_states else compiled
        return function(
            tree,
            likelihood_by_leaf,
            rate_matrix,
            fixed_transition_matrices=fixed_transition_matrices,
        )

    wrapped.compiled = compiled
    return wrapped


CANDIDATE = build_candidate(MEASURED_MINIMUM_STATES)


@contextmanager
def patched(function):
    previous = asr._compute_log_inside_likelihoods
    asr._compute_log_inside_likelihoods = function
    try:
        yield
    finally:
        asr._compute_log_inside_likelihoods = previous


def fixture(tips, states, *, unique=False):
    rng = np.random.default_rng(9090 + tips + states)
    count = 0

    def length():
        return float(rng.uniform(0.2, 0.8)) if unique else 0.4

    parts = [f"T{i}:{length():.17g}" for i in range(tips)]
    while len(parts) > 1:
        following = []
        for start in range(0, len(parts), 2):
            if start + 1 == len(parts):
                following.append(parts[start])
            else:
                following.append(
                    f"({parts[start]},{parts[start + 1]})I{count}:{length():.17g}"
                )
                count += 1
        parts = following
    tree = read_tree(parts[0] + ";", "1", True, quiet=True, rooted="yes")
    labels = tuple(f"s{i}" for i in range(states))
    observed, likelihood = {}, {}
    for i, leaf in enumerate(tree.leaves()):
        row = np.zeros(states)
        row[i % states] = 1.0
        value = labels[i % states]
        if i % 11 == 0:
            row[:] = 1.0
            value = None
        elif i % 7 == 0:
            row[(i + 1) % states] = 1.0
            value = (labels[i % states], labels[(i + 1) % states])
        observed[leaf.name] = value
        likelihood[leaf.name] = row
    return tree, labels, observed, likelihood


def full_call(arguments, function, *, fitted=False):
    options = {} if fitted else {"rate": 0.15 / len(arguments[1])}
    with patched(function):
        return asr.compute_mk_marginals(*arguments, model="ER", **options)


def compare(first, second):
    differences = []
    exact = True

    def visit(left, right, path):
        nonlocal exact
        if isinstance(left, dict):
            assert set(left) == set(right), path
            for key in left:
                visit(left[key], right[key], path + "/" + str(key))
        elif isinstance(left, (tuple, list)):
            assert len(left) == len(right), path
            for index, (a, b) in enumerate(zip(left, right, strict=True)):
                visit(a, b, path + "/" + str(index))
        elif isinstance(left, (float, int, np.ndarray, np.number)) and not isinstance(
            left, bool
        ):
            a, b = np.asarray(left), np.asarray(right)
            assert a.shape == b.shape, path
            for predicate in (np.isnan, np.isneginf, np.isposinf):
                np.testing.assert_array_equal(predicate(a), predicate(b), err_msg=path)
            finite = np.isfinite(a) & np.isfinite(b)
            delta = np.abs(a[finite] - b[finite])
            absolute = float(np.max(delta, initial=0.0))
            scaled = float(
                np.max(delta / np.maximum(1, np.abs(a[finite])), initial=0.0)
            )
            exact &= bool(np.array_equal(a, b, equal_nan=True))
            if absolute:
                differences.append(
                    {"path": path, "max_absolute": absolute, "max_scaled": scaled}
                )
        else:
            assert left == right, path

    visit(first, second, "result")
    return {
        "all_fields_exact": exact,
        "max_absolute": max((x["max_absolute"] for x in differences), default=0.0),
        "max_scaled": max((x["max_scaled"] for x in differences), default=0.0),
        "different_fields": differences,
    }


def edge_checks():
    output = []
    rng = np.random.default_rng(911)
    for states in (2, 8, 16, 64):
        for rate in (0.0, 1e-300, 0.01, 10.0):
            for kind in ("one", "ambiguous", "missing", "impossible", "dynamic"):
                matrix = asr._er_transition_matrix(0.7, rate, states)
                inside = rng.normal(-20, 4, states)
                if kind == "one":
                    inside[:] = -np.inf
                    inside[0] = 0.0
                elif kind == "ambiguous":
                    inside[:] = -np.inf
                    inside[:2] = 0.0
                elif kind == "missing":
                    inside[:] = 0.0
                elif kind == "impossible":
                    inside[:] = -np.inf
                expected = logsumexp(
                    asr._log_probabilities(matrix) + inside[None, :], axis=1
                )
                actual = candidate_message(matrix, None, inside, {})
                np.testing.assert_allclose(expected, actual, rtol=1e-13, atol=1e-13)
                output.append(
                    {
                        "states": states,
                        "rate": rate,
                        "kind": kind,
                        **compare(expected, actual),
                    }
                )
    # Exact structure checking must reject even a one-ULP nonER perturbation.
    matrix = asr._er_transition_matrix(1, 0.1, 16)
    for kind in ("one_ulp", "unequal_diagonal", "negative", "diagonal_smaller", "nan"):
        changed = matrix.copy()
        if kind == "one_ulp":
            changed[0, 1] = np.nextafter(changed[0, 1], np.inf)
        elif kind == "unequal_diagonal":
            changed[0, 0] = np.nextafter(changed[0, 0], np.inf)
        elif kind == "negative":
            changed[0, 1] = -0.1
        elif kind == "diagonal_smaller":
            np.fill_diagonal(changed, 0.0)
        else:
            changed[0, 1] = np.nan
        assert er_descriptor(changed, {}) is None
        output.append({"fallback": kind})
    return output


def independent_enumeration():
    tree, states, observed, likelihood = fixture(4, 4, unique=True)
    rate = 0.07
    matrix = asr._build_rate_matrix("ER", states, [rate])
    transitions = asr._transition_matrices_for_tree(tree, matrix)
    nodes = tuple(tree.traverse("preorder"))
    indices = {node: index for index, node in enumerate(nodes)}
    marginal = np.zeros((len(nodes), len(states)))
    total = 0.0
    for assignment in itertools.product(range(len(states)), repeat=len(nodes)):
        weight = 1.0 / len(states)
        for index, node in enumerate(nodes):
            if not node.is_root:
                weight *= transitions[node][
                    assignment[indices[node.up]], assignment[index]
                ]
            if node.is_leaf:
                weight *= likelihood[node.name][assignment[index]]
        total += weight
        for index, value in enumerate(assignment):
            marginal[index, value] += weight
    forced = build_candidate(2)
    with patched(forced):
        posterior, fit = asr.compute_mk_marginals(
            tree, states, observed, likelihood, rate=rate
        )
    actual = np.asarray([posterior[node] for node in nodes])
    np.testing.assert_allclose(actual, marginal / total, rtol=1e-12, atol=1e-13)
    assert abs(fit["log_likelihood"] - math.log(total)) < 1e-12
    return {
        "assignments": len(states) ** len(nodes),
        "posterior_max_absolute": float(np.max(np.abs(actual - marginal / total))),
        "log_likelihood_absolute": abs(fit["log_likelihood"] - math.log(total)),
    }


def measure():
    with measurement_lock():
        result = {
            "environment": environment(),
            "prototype_configuration": PROTOTYPE_CONFIGURATION,
            "source_sha256": hashlib.sha256(
                inspect.getsource(ORIGINAL).encode()
            ).hexdigest(),
            "edges": edge_checks(),
            "independent": independent_enumeration(),
            "cases": [],
        }
        for tips, states, unique, fitted in (
            (32, 2, False, False),
            (64, 8, False, False),
            (128, 16, False, False),
            (128, 32, False, False),
            (128, 64, False, False),
            (128, 64, True, False),
            (32, 16, False, True),
        ):
            arguments = fixture(tips, states, unique=unique)
            before = partial(full_call, arguments, ORIGINAL, fitted=fitted)
            after = partial(full_call, arguments, CANDIDATE, fitted=fitted)
            entry = {
                "tips": tips,
                "states": states,
                "unique_lengths": unique,
                "fitted_rate": fitted,
                "comparison": compare(before(), after()),
                "timing": paired(before, after),
                "baseline_peak_bytes": peak_bytes(before),
                "candidate_peak_bytes": peak_bytes(after),
            }
            print(
                {
                    key: entry[key]
                    for key in ("tips", "states", "unique_lengths", "fitted_rate")
                },
                entry["timing"]["speedup"],
                entry["comparison"]["max_scaled"],
                flush=True,
            )
            result["cases"].append(entry)
        q = asr._build_rate_matrix("ER", tuple(str(i) for i in range(64)), [0.002])
        matrix = asr._transition_matrix(q, 0.4)
        inside = np.linspace(-20, 0, 64)
        cache = {}

        def before():
            return logsumexp(asr._log_probabilities(matrix) + inside[None, :], axis=1)

        after = partial(candidate_message, matrix, None, inside, cache)
        result["kernel64"] = paired(before, after)
        (HERE / "results.json").write_text(
            json.dumps(result, indent=2, allow_nan=False) + "\n"
        )


def tests():
    with measurement_lock():
        env = os.environ.copy()
        env["PYTHONPATH"] = os.pathsep.join((str(SOURCE), str(HERE), str(HERE.parent)))
        env["DYLD_FALLBACK_LIBRARY_PATH"] = "/Users/kf/miniforge3-arm64/lib"
        command = [
            sys.executable,
            "tools/check.py",
            "test",
            "--",
            "-p",
            "candidate09_plugin",
            "tests/test_asr.py",
            "tests/test_discrete_asr_models.py",
            "-k",
            "not stochastic and not map and not rsc and not origin",
            "-rs",
        ]
        print(" ".join(command), flush=True)
        subprocess.run(command, cwd=SOURCE, env=env, check=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tests", action="store_true")
    options = parser.parse_args()
    tests() if options.tests else measure()
