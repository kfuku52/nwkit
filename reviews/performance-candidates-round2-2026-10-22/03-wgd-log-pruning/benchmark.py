"""Bounded batching of independent rare-family log-pruning reductions.

Run this file with the Python environment and immutable source in PROTOCOL.md.
No production file is modified. The complete consuming operation is reported
separately from the transition-prepared pruning kernel.
"""

import argparse
import hashlib
import inspect
import json
import os
import subprocess
import sys
import textwrap
from contextlib import contextmanager
from functools import partial
from pathlib import Path

HERE = Path(__file__).resolve().parent
COMMON = HERE.parent
SOURCE = Path(
    os.environ.get(
        "NWKIT_BENCH_SOURCE", "/tmp/nwkit-perf-round2-20261022-rkzoGe/baseline"
    )
)
for key in (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "MKL_NUM_THREADS",
):
    os.environ[key] = "1"
os.environ["NWKIT_CI_DIAGNOSTICS"] = "0"
sys.path.insert(0, str(COMMON))
sys.path.insert(0, str(SOURCE))

from support import environment, measurement_lock, paired, peak_bytes  # noqa: E402

DEFAULT_WORKSPACE_BYTES = 1024**2
WORKSPACE_BYTES = DEFAULT_WORKSPACE_BYTES
MAX_BATCH = 64


def batched_contribution(matrix, values):
    import numpy as np
    from scipy.special import logsumexp

    # Four float-sized temporaries conservatively estimate the broadcast and
    # stabilized logsumexp workspace. Large state spaces retain the old loop.
    batch = max(1, min(MAX_BATCH, WORKSPACE_BYTES // (32 * matrix.size)))
    if len(values) <= 1 or batch == 1:
        return np.array([logsumexp(matrix + row[None, :], axis=1) for row in values])
    contribution = np.empty((len(values), matrix.shape[0]))
    for start in range(0, len(values), batch):
        end = min(start + batch, len(values))
        contribution[start:end] = logsumexp(
            matrix[None, :, :] + values[start:end, None, :], axis=2
        )
    return contribution


def build_candidate():
    from nwkit import wgd_count_model as module

    source = textwrap.dedent(inspect.getsource(module.CountLikelihood._prune_log))
    old = """            contribution = np.array(
                [
                    logsumexp(matrices[child] + row[None, :], axis=1)
                    for row in likelihoods.pop(child)
                ]
            )"""
    new = """            contribution = _candidate_contribution(
                matrices[child], likelihoods.pop(child)
            )"""
    assert source.count(old) == 1, "Prototype requires the pinned baseline body"
    source = source.replace(old, new)
    namespace = dict(module.__dict__, _candidate_contribution=batched_contribution)
    exec(compile(source, str(HERE / "benchmark.py"), "exec"), namespace)
    return namespace["_prune_log"]


@contextmanager
def prototype(function=None):
    from nwkit.wgd_count_model import CountLikelihood

    baseline = CountLikelihood._prune_log
    CountLikelihood._prune_log = build_candidate() if function is None else function
    try:
        yield
    finally:
        CountLikelihood._prune_log = baseline


def equivalent(first, second):
    import numpy as np

    first, second = np.asarray(first), np.asarray(second)
    assert first.shape == second.shape
    for predicate in (np.isnan, np.isneginf, np.isposinf):
        np.testing.assert_array_equal(predicate(first), predicate(second))
    finite = np.isfinite(first) & np.isfinite(second)
    difference = float(np.max(np.abs(first[finite] - second[finite]), initial=0.0))
    np.testing.assert_allclose(first, second, rtol=1e-13, atol=1e-12)
    return {
        "array_equal": bool(np.array_equal(first, second, equal_nan=True)),
        "max_abs_difference": difference,
        "negative_infinity_count": int(np.sum(np.isneginf(second))),
        "nan_count": int(np.sum(np.isnan(second))),
    }


def tree_fixture(counts, *, detection=None, scales=None, ascertainment="observed"):
    import numpy as np

    from nwkit.wgd_count_model import CountLikelihood, CountTree

    tree = CountTree(
        (-1, 0, 1, 1, 0),
        (0.0, 0.4, 0.25, 0.45, 0.7),
        (2, 3, 4),
        ("A", "B", "C"),
        (0, 1, 2, 3, 4),
        ("ABC", "AB", "A", "B", "C"),
    )
    return CountLikelihood(
        tree,
        np.asarray(counts, dtype=float),
        detection=detection,
        rate_scales=scales,
        ascertainment=ascertainment,
    )


def rare_counts(families, maximum=84, *, missing=True):
    import numpy as np

    counts = np.ones((families, 3))
    counts[:, 0] = maximum - np.arange(families) % min(16, maximum - 1)
    if missing:
        counts[1::4, 1] = np.nan
        counts[2::4, 2] = np.nan
    return counts


def full_call(model, rates, root_mean, bound, event, implementation):
    with prototype(implementation):
        return model.family_log_likelihoods(rates, root_mean, bound, event)


def rare_switches(model, rates, root_mean, bound, event, implementation):
    from nwkit.wgd_count_model import CountLikelihood

    rows = []

    def counted(self, matrices, mean, max_count, selected):
        rows.append(len(selected))
        return implementation(self, matrices, mean, max_count, selected)

    original = CountLikelihood._prune_log
    CountLikelihood._prune_log = counted
    try:
        model.family_log_likelihoods(rates, root_mean, bound, event)
    finally:
        CountLikelihood._prune_log = original
    return rows


def measure_case(name, baseline, candidate, metadata):
    result = {
        "name": name,
        **metadata,
        "equivalence": equivalent(baseline(), candidate()),
        "timing": paired(baseline, candidate),
        "peak_python_bytes": {
            "baseline": peak_bytes(baseline),
            "candidate": peak_bytes(candidate),
        },
    }
    print(name, result["timing"]["speedup"], flush=True)
    return result


def check_edges(baseline, candidate):
    import numpy as np

    from nwkit.wgd_count_model import MultiplicationEvent

    checks = []
    rng = np.random.default_rng(73019)
    counts = rng.integers(0, 7, (19, 3)).astype(float)
    counts[:, 0] += 1
    counts[1::3, 1:] = np.nan
    counts[2::4, 0] = np.nan
    # The missing patterns leave an observed copy in every row.
    counts[2::4, 2] = 1
    for bound in (8, 20, 86):
        for detection in (np.ones(3), np.array([0.65, 0.9, 0.8])):
            model = tree_fixture(counts, detection=detection)
            for rates in ([[0.2, 0.3]], [[1e-8, 0]], [[0, 0.3]], [[0, 0]]):
                for retention in (None, 0.0, 0.3, 1.0):
                    event = (
                        None if retention is None else MultiplicationEvent(1, retention)
                    )
                    matrices = model._log_transitions(rates, 1, bound, event)
                    rows = np.array([18, 1, 2, 9, 1, 0, 17])
                    for mean in (1.0, 1.7):
                        checks.append(
                            {
                                "bound": bound,
                                "detection": detection.tolist(),
                                "rates": rates,
                                "retention": retention,
                                "root_mean": mean,
                                **equivalent(
                                    baseline(model, matrices, mean, bound, rows),
                                    candidate(model, matrices, mean, bound, rows),
                                ),
                            }
                        )
    # Low budget and state-size fallback exercise both paths and last short tiles.
    global WORKSPACE_BYTES
    for state_count in (1, 2, 17, 87, 513):
        matrix = -rng.uniform(0, 1e6, (state_count, state_count))
        matrix[rng.random(matrix.shape) < 0.25] = -np.inf
        values = -rng.uniform(0, 1e5, (67, state_count))
        values[1, :] = -np.inf
        values[2, ::2] = -np.inf
        from scipy.special import logsumexp

        expected = np.array(
            [logsumexp(matrix + row[None, :], axis=1) for row in values]
        )
        for budget in (1, 16384, 1024**2, 4 * 1024**2):
            WORKSPACE_BYTES = budget
            checks.append(
                {
                    "states": state_count,
                    "workspace_bytes": budget,
                    **equivalent(expected, batched_contribution(matrix, values)),
                }
            )
    WORKSPACE_BYTES = DEFAULT_WORKSPACE_BYTES
    convergence = []
    for retention in (None, 0.0, 0.3):
        model = tree_fixture([[40, 1, 1], [39, np.nan, 1]])
        event = None if retention is None else MultiplicationEvent(1, retention)
        with prototype(baseline):
            first = model.convergence([[1e-8, 0]], 1, 42, event)
        with prototype(candidate):
            second = model.convergence([[1e-8, 0]], 1, 42, event)
        convergence.append({"retention": retention, **equivalent(first, second)})
    return {"pruning_and_helper_cases": checks, "state_doubling": convergence}


def run_benchmarks():
    import numpy as np

    import nwkit
    from nwkit import wgd_count_model as module
    from nwkit.wgd_count_model import MultiplicationEvent

    baseline = module.CountLikelihood._prune_log
    candidate = build_candidate()
    source_file = Path(module.__file__)
    results = {
        "environment": environment(),
        "source": {
            "nwkit_file": nwkit.__file__,
            "target_file": str(source_file),
            "target_sha256": hashlib.sha256(source_file.read_bytes()).hexdigest(),
        },
        "workspace_estimate_bytes": WORKSPACE_BYTES,
        "maximum_family_batch": MAX_BATCH,
        "numerical_checks": check_edges(baseline, candidate),
        "pruning_kernel": [],
        "complete_family_likelihood": [],
    }
    for families, bound in (
        (1, 20),
        (16, 20),
        (128, 20),
        (128, 86),
        (512, 86),
        (128, 160),
    ):
        counts = rare_counts(families, maximum=min(84, bound))
        model = tree_fixture(counts, detection=np.array([0.8, 0.9, 1.0]))
        matrices = model._log_transitions([[1e-8, 0]], 1, bound, None)
        rows = np.arange(families)
        results["pruning_kernel"].append(
            measure_case(
                f"prune_{families}x{bound}",
                partial(baseline, model, matrices, 1, bound, rows),
                partial(candidate, model, matrices, 1, bound, rows),
                {
                    "families": families,
                    "max_count": bound,
                    "setup": "precomputed transitions; includes emissions and complete log pruning",
                },
            )
        )
    cases = [
        (
            "ordinary_control",
            np.tile([[1, 2, 3], [4, 2, 1], [0, 1, np.nan]], (43, 1)),
            [[0.2, 0.3]],
            20,
            None,
            None,
        ),
        ("rare_1", rare_counts(1), [[1e-8, 0]], 86, None, None),
        ("rare_16", rare_counts(16), [[1e-8, 0]], 86, None, None),
        ("rare_128", rare_counts(128), [[1e-8, 0]], 86, None, None),
        ("rare_512", rare_counts(512), [[1e-8, 0]], 86, None, None),
        (
            "rare_128_larger_states",
            rare_counts(128, maximum=150),
            [[1e-5, 0]],
            160,
            None,
            None,
        ),
        (
            "rare_pulse_detection",
            rare_counts(64),
            [[1e-8, 0]],
            86,
            MultiplicationEvent(1, 0.3),
            np.array([0.8, 0.9, 1.0]),
        ),
    ]
    for name, counts, rates, bound, event, detection in cases:
        model = tree_fixture(counts, detection=detection)
        metadata = {
            "families": len(counts),
            "tips": 3,
            "max_count": bound,
            "rates": rates,
            "event": None if event is None else vars(event),
            "setup": "complete family_log_likelihoods; includes emissions, probability pruning, transitions, log recovery and ascertainment",
            "rare_rows_per_category": rare_switches(
                model, rates, 1, bound, event, baseline
            ),
        }
        results["complete_family_likelihood"].append(
            measure_case(
                name,
                partial(full_call, model, rates, 1, bound, event, baseline),
                partial(full_call, model, rates, 1, bound, event, candidate),
                metadata,
            )
        )
    (HERE / "results.json").write_text(
        json.dumps(results, indent=2, allow_nan=False) + "\n"
    )


def run_tests():
    env = os.environ.copy()
    env["PYTHONPATH"] = str(HERE) + os.pathsep + str(SOURCE)
    command = [
        sys.executable,
        "tools/check.py",
        "test",
        "--",
        "-p",
        "candidate03_plugin",
        "tests/test_wgd_count_model.py",
        "tests/test_wgd_count.py",
        "-rs",
    ]
    completed = subprocess.run(
        command, cwd=SOURCE, env=env, text=True, capture_output=True
    )
    (HERE / "tests.log").write_text(
        "COMMAND: "
        + " ".join(command)
        + "\nCWD: "
        + str(SOURCE)
        + "\n"
        + completed.stdout
        + completed.stderr
    )
    print(completed.stdout[-2000:])
    completed.check_returncode()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tests", action="store_true")
    args = parser.parse_args()
    with measurement_lock():
        if args.tests:
            run_tests()
        else:
            run_benchmarks()


if __name__ == "__main__":
    main()
