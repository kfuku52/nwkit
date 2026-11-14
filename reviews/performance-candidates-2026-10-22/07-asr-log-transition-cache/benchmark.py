"""Evaluate call-local caching of repeated fixed transition logarithms."""

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
from collections import Counter
from contextlib import contextmanager
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
LOCK = Path("/tmp/nwkit-perf-20261022-9fYwTD/measurement.lock")
sys.path.insert(0, str(ROOT))


def build_candidate(policy="shared_only"):
    import nwkit.asr as asr

    source = inspect.getsource(asr._compute_log_inside_likelihoods)
    initialization = "    transition_matrix_cache = dict()\n"
    assert source.count(initialization) == 1
    addition = "    fixed_log_cache = {}\n"
    if policy == "shared_only":
        addition += """    fixed_shared_ids = set()
    if fixed_transition_matrices is not None:
        fixed_shared_ids = {
            key for key, count in _candidate07_counter(
                id(matrix) for matrix in fixed_transition_matrices.values()
            ).items() if count > 1
        }
"""
    source = source.replace(initialization, initialization + addition)
    previous = """                matrix = fixed_transition_matrices[child]
                log_matrix = _log_probabilities(matrix)
"""
    common = """                matrix = fixed_transition_matrices[child]
                matrix_id = id(matrix)
"""
    if policy == "shared_only":
        replacement = (
            common
            + """                if matrix_id in fixed_shared_ids:
                    if matrix_id not in fixed_log_cache:
                        fixed_log_cache[matrix_id] = _log_probabilities(matrix)
                    log_matrix = fixed_log_cache[matrix_id]
                else:
                    log_matrix = _log_probabilities(matrix)
"""
        )
    else:
        replacement = (
            common
            + """                if matrix_id not in fixed_log_cache:
                    fixed_log_cache[matrix_id] = _log_probabilities(matrix)
                log_matrix = fixed_log_cache[matrix_id]
"""
        )
    assert source.count(previous) == 1
    source = source.replace(previous, replacement)
    namespace = {}
    globals_ = asr.__dict__
    globals_["_candidate07_counter"] = Counter
    exec(compile(source, str(HERE / "benchmark.py"), "exec"), globals_, namespace)
    return namespace["_compute_log_inside_likelihoods"]


@contextmanager
def prototype(policy="shared_only"):
    import nwkit.asr as asr

    original = asr._compute_log_inside_likelihoods
    previous_counter = asr.__dict__.get("_candidate07_counter")
    asr._compute_log_inside_likelihoods = build_candidate(policy)
    try:
        yield
    finally:
        asr._compute_log_inside_likelihoods = original
        if previous_counter is None:
            asr.__dict__.pop("_candidate07_counter", None)
        else:
            asr._candidate07_counter = previous_counter


def equivalent(first, second):
    import numpy as np

    if isinstance(first, dict):
        assert first.keys() == second.keys()
        for key in first:
            equivalent(first[key], second[key])
    elif isinstance(first, (tuple, list)):
        assert len(first) == len(second)
        for a, b in zip(first, second, strict=True):
            equivalent(a, b)
    elif isinstance(first, (np.ndarray, float, int, np.number)):
        np.testing.assert_array_equal(first, second, strict=False)
    else:
        assert first == second


def fixture(tips, states, unique=False, generator="symmetric", shape="balanced"):
    import numpy as np
    from ete4 import Tree

    if shape == "star":
        tree = Tree()
        for index in range(tips):
            tree.add_child(name=f"T{index}", dist=0.2)
    else:
        nodes = []
        for index in range(tips):
            leaf = Tree()
            leaf.name, leaf.dist = f"T{index}", 0.2
            nodes.append(leaf)
        while len(nodes) > 1:
            parents = []
            for index in range(0, len(nodes), 2):
                if index + 1 == len(nodes):
                    parents.append(nodes[index])
                    continue
                parent = Tree()
                parent.dist = 0.2
                parent.add_child(nodes[index])
                parent.add_child(nodes[index + 1])
                parents.append(parent)
            nodes = parents
        tree = nodes[0]
    tree.dist = 0
    for index, node in enumerate(tree.traverse()):
        if not node.is_root:
            node.dist = 0.02 + index * 0.001 if unique else 0.2
    rng = np.random.default_rng(6701 + tips + states)
    matrix = rng.uniform(0.05, 0.2, size=(states, states))
    if generator in ("symmetric", "disconnected"):
        matrix = (matrix + matrix.T) / 2
    if generator == "disconnected":
        half = states // 2
        matrix[:half, half:] = matrix[half:, :half] = 0
    np.fill_diagonal(matrix, 0)
    np.fill_diagonal(matrix, -matrix.sum(axis=1))
    observations = {
        name: index % states for index, name in enumerate(tree.leaf_names())
    }
    likelihood = {name: np.eye(states)[state] for name, state in observations.items()}
    return tree, matrix, observations, likelihood


def measure(functions):
    warm = {}
    for name, function in functions.items():
        start = time.perf_counter()
        function()
        warm[name] = time.perf_counter() - start
    calls = max(1, min(500, int(0.04 / min(warm.values())) + 1))
    samples = {name: [] for name in functions}
    for repetition in range(7):
        names = list(functions)
        if repetition % 2:
            names.reverse()
        for name in names:
            start = time.perf_counter()
            for _ in range(calls):
                functions[name]()
            samples[name].append((time.perf_counter() - start) / calls)
    result = {}
    for name, values in samples.items():
        result[name] = {
            "samples_seconds_per_call": values,
            "median_seconds": statistics.median(values),
            "min_seconds": min(values),
            "max_seconds": max(values),
            "relative_range": (max(values) - min(values)) / statistics.median(values),
        }
        result[name]["speedup"] = statistics.median(
            samples["baseline"]
        ) / statistics.median(values)
    return {"batch_calls": calls, "timings": result}


def peak(function):
    tracemalloc.start()
    try:
        function()
        return tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()


def call_with(policy, function):
    if policy == "baseline":
        return function()
    with prototype(policy):
        return function()


def count_logs(policy, function):
    import nwkit.asr as asr

    original = asr._log_probabilities
    calls = []

    def counting(value):
        calls.append(value.ndim)
        return original(value)

    asr._log_probabilities = counting
    try:
        call_with(policy, function)
    finally:
        asr._log_probabilities = original
    return {
        "total": len(calls),
        "matrix_calls": calls.count(2),
        "vector_calls": calls.count(1),
    }


def inside_case(
    tips, states, unique=False, generator="symmetric", shape="balanced", regime=False
):
    import nwkit.asr as asr

    tree, matrix, _, likelihood = fixture(tips, states, unique, generator, shape)
    matrices = asr._transition_matrices_for_tree(tree, matrix)
    rate = matrix
    if regime:
        fast = matrix * 2
        rate = {
            node: matrix if index % 2 else fast for index, node in enumerate(matrices)
        }
        matrices = None
    baseline = asr._compute_log_inside_likelihoods
    candidates = {
        policy: build_candidate(policy) for policy in ("cached", "shared_only")
    }
    functions = {
        "baseline": lambda: baseline(
            tree, likelihood, rate, fixed_transition_matrices=matrices
        ),
        **{
            policy: lambda candidate=candidate: candidate(
                tree, likelihood, rate, fixed_transition_matrices=matrices
            )
            for policy, candidate in candidates.items()
        },
    }
    output = functions["baseline"]()
    for function in functions.values():
        equivalent(output, function())
    result = {
        "tips": tips,
        "states": states,
        "unique_lengths": unique,
        "generator": generator,
        "shape": shape,
        "regime_path": regime,
        "equivalence": "all nested numeric arrays and scalars exactly equal",
        **measure(functions),
        "peak_python_bytes": {
            name: peak(function) for name, function in functions.items()
        },
    }
    # Count through the same original call in a scoped prototype.
    for policy in ("baseline", "cached", "shared_only"):

        def dynamic_call():
            return asr._compute_log_inside_likelihoods(
                tree, likelihood, rate, fixed_transition_matrices=matrices
            )

        result.setdefault("log_calls", {})[policy] = count_logs(policy, dynamic_call)
    print(
        json.dumps(
            {
                "inside": {
                    key: value
                    for key, value in result.items()
                    if key not in ("peak_python_bytes", "log_calls")
                }
            }
        ),
        flush=True,
    )
    return result


def marginal_case(tips, states, unique=False):
    import nwkit.asr as asr

    tree, matrix, observations, likelihood = fixture(tips, states, unique)
    state_names = [str(index) for index in range(states)]
    observed_names = {name: str(state) for name, state in observations.items()}

    def fit():
        return asr.compute_mk_marginals(
            tree,
            state_names,
            observed_names,
            likelihood,
            model="CUSTOM",
            fixed_rate_matrix=matrix,
        )

    expected = fit()
    equivalent(expected, call_with("shared_only", fit))
    # Install once outside timing: no source parsing/context overhead in fit time.
    baseline = asr._compute_log_inside_likelihoods
    candidate = build_candidate()

    def run(function):
        old = asr._compute_log_inside_likelihoods
        asr._compute_log_inside_likelihoods = function
        try:
            return fit()
        finally:
            asr._compute_log_inside_likelihoods = old

    functions = {
        "baseline": lambda: run(baseline),
        "shared_only": lambda: run(candidate),
    }
    result = {
        "tips": tips,
        "states": states,
        "unique_lengths": unique,
        "equivalence": "posterior and every fit field exactly equal",
        **measure(functions),
        "peak_python_bytes": {
            name: peak(function) for name, function in functions.items()
        },
    }
    print(json.dumps({"marginal": result}), flush=True)
    return result


def edges():
    import numpy as np

    import nwkit.asr as asr

    records = []
    for generator in ("symmetric", "asymmetric", "disconnected"):
        tree, matrix, _, likelihood = fixture(12, 4, generator=generator)
        nodes = [node for node in tree.traverse() if not node.is_root]
        nodes[0].dist = 0
        for node in tree.leaves():
            likelihood[node.name] = np.ones(4)
        likelihood[next(tree.leaf_names())] = np.array([0.0, 1.0, 0.0, 0.0])
        matrices = asr._transition_matrices_for_tree(tree, matrix)
        # A deliberately impossible child state tests -inf propagation.
        for impossible in (False, True):
            if impossible:
                likelihood[next(tree.leaf_names())] = np.zeros(4)
            baseline = asr._compute_log_inside_likelihoods(
                tree, likelihood, matrix, fixed_transition_matrices=matrices
            )
            for policy in ("cached", "shared_only"):
                with prototype(policy):
                    equivalent(
                        baseline,
                        asr._compute_log_inside_likelihoods(
                            tree, likelihood, matrix, fixed_transition_matrices=matrices
                        ),
                    )
            records.append(
                {"generator": generator, "impossible": impossible, "exact": True}
            )
        # Mutate existing arrays between calls; local caches must not survive.
        likelihood = {name: np.ones(4) for name in tree.leaf_names()}
        with prototype():
            asr._compute_log_inside_likelihoods(
                tree, likelihood, matrix, fixed_transition_matrices=matrices
            )
            for transition in {
                id(value): value for value in matrices.values()
            }.values():
                transition[:] = np.eye(4)
            actual = asr._compute_log_inside_likelihoods(
                tree, likelihood, matrix, fixed_transition_matrices=matrices
            )
        expected = asr._compute_log_inside_likelihoods(
            tree, likelihood, matrix, fixed_transition_matrices=matrices
        )
        equivalent(expected, actual)
        records.append(
            {"generator": generator, "mutated_between_calls": True, "exact": True}
        )
    return records


def tests():
    env = os.environ.copy()
    env["PYTHONPATH"] = str(HERE) + os.pathsep + str(ROOT)
    env["PYTEST_PLUGINS"] = "candidate07_plugin"
    env["NWKIT_CI_DIAGNOSTICS"] = "0"
    command = [
        sys.executable,
        "tools/check.py",
        "test",
        "--",
        "tests/test_asr.py",
        "tests/test_asr_polytomy.py",
        "tests/test_asr_regimes.py",
        "tests/test_discrete_asr_models.py",
        "-rs",
    ]
    with (HERE / "prototype-tests.log").open("w") as stream:
        stream.write(" ".join(command) + "\n")
        stream.flush()
        completed = subprocess.run(
            command, cwd=ROOT, env=env, stdout=stream, stderr=subprocess.STDOUT
        )
    if completed.returncode:
        raise SystemExit(completed.returncode)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tests", action="store_true")
    args = parser.parse_args()
    for key in (
        "OPENBLAS_NUM_THREADS",
        "OMP_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "MKL_NUM_THREADS",
    ):
        os.environ[key] = "1"
    LOCK.parent.mkdir(parents=True, exist_ok=True)
    with LOCK.open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if args.tests:
            tests()
            return
        import numpy as np
        import scipy

        from nwkit import __version__

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
            },
            "inside": [
                inside_case(*case)
                for case in (
                    (128, 4),
                    (512, 4),
                    (256, 16),
                    (128, 32),
                    (256, 8, True),
                    (128, 8, False, "asymmetric"),
                    (128, 8, False, "disconnected"),
                    (128, 4, False, "symmetric", "star"),
                    (128, 4, False, "symmetric", "balanced", True),
                )
            ],
            "marginals": [
                marginal_case(*case) for case in ((128, 4), (128, 16), (128, 8, True))
            ],
            "edge_checks": edges(),
        }
        (HERE / "results.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
