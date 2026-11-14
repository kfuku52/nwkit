"""Memory-bounded row batching of the WGD log-semiring matrix product."""

import argparse
import hashlib
import inspect
import json
import os
import sys
from contextlib import contextmanager
from functools import partial
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = Path("/tmp/nwkit-perf-round2-20261022-rkzoGe/baseline")
for key in (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "MKL_NUM_THREADS",
):
    os.environ[key] = "1"
os.environ["NWKIT_CI_DIAGNOSTICS"] = "0"
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(SOURCE))

import numpy as np  # noqa: E402
from scipy.special import logsumexp  # noqa: E402
from scipy.stats import binom  # noqa: E402
from support import environment, measurement_lock, paired, peak_bytes  # noqa: E402

from nwkit import wgd_count_model as model  # noqa: E402

ORIGINAL = model._log_matrix_product
WORKSPACE_BYTES = 1024**2
MAX_BATCH = 32


def batched_product(left, right):
    """Retain the k-axis order while batching independent output rows."""
    batch = max(1, min(MAX_BATCH, WORKSPACE_BYTES // max(1, 32 * right.size)))
    if len(left) <= 1 or batch == 1:
        return ORIGINAL(left, right)
    stop = min(len(left), batch)
    first = logsumexp(left[:stop, :, None] + right[None, :, :], axis=1)
    output = np.empty((len(left), right.shape[1]), dtype=first.dtype)
    output[:stop] = first
    del first
    for start in range(stop, len(left), batch):
        end = min(len(left), start + batch)
        output[start:end] = logsumexp(
            left[start:end, :, None] + right[None, :, :], axis=1
        )
    return output


@contextmanager
def prototype(implementation):
    previous = model._log_matrix_product
    model._log_matrix_product = implementation
    try:
        yield
    finally:
        model._log_matrix_product = previous


def equivalent(before, after):
    before, after = np.asarray(before), np.asarray(after)
    np.testing.assert_array_equal(before, after)
    assert before.dtype == after.dtype
    finite = np.isfinite(before) & np.isfinite(after)
    return {
        "exact": True,
        "shape": list(before.shape),
        "dtype": str(before.dtype),
        "max_absolute_difference": float(
            np.max(np.abs(before[finite] - after[finite]), initial=0)
        ),
        "negative_infinity_count": int(np.isneginf(before).sum()),
        "nan_count": int(np.isnan(before).sum()),
    }


def measured(name, baseline, candidate, scope, **metadata):
    return {
        "case": name,
        "scope": scope,
        **metadata,
        "comparison": equivalent(baseline(), candidate()),
        "timing": paired(baseline, candidate),
        "peak_bytes": {
            "baseline": peak_bytes(baseline),
            "candidate": peak_bytes(candidate),
        },
    }


def pulse(bound, retention=0.6, multiplicity=2):
    value = np.full((bound + 1, bound + 1), -np.inf)
    for count in range(bound + 1):
        value[count, count:] = binom.logpmf(
            np.arange(bound - count + 1), count * (multiplicity - 1), retention
        )
    return value


def count_fixture(families=1, maximum=84, *, rare=True, missing=False):
    tree = model.CountTree(
        (-1, 0, 1, 1, 0),
        (0.0, 0.4, 0.25, 0.45, 0.7),
        (2, 3, 4),
        ("A", "B", "C"),
        (0, 1, 2, 3, 4),
        ("ABC", "AB", "A", "B", "C"),
    )
    counts = np.ones((families, 3))
    counts[:, 0] = maximum if rare else 2
    if missing:
        counts[::2, 1] = np.nan
    return model.CountLikelihood(tree, counts, ascertainment="root-clades")


def call_with(implementation, function):
    with prototype(implementation):
        return function()


def count_calls(function):
    calls = []

    def counted(left, right):
        calls.append([list(left.shape), list(right.shape)])
        return ORIGINAL(left, right)

    call_with(counted, function)
    return calls


def edge_checks():
    global WORKSPACE_BYTES
    rng = np.random.default_rng(1948)
    checks = []
    for size in (1, 2, 3, 7, 17, 41):
        for dtype in (np.float32, np.float64):
            left = (-rng.exponential(12, (size, size))).astype(dtype)
            right = (-rng.exponential(12, (size, size))).astype(dtype)
            left[rng.random(left.shape) < 0.25] = -np.inf
            right[rng.random(right.shape) < 0.25] = -np.inf
            checks.append(
                equivalent(ORIGINAL(left, right), batched_product(left, right))
            )
    for shape in ((0, 5, 7), (1, 5, 7), (9, 20, 11), (30, 5, 9)):
        rows, middle, columns = shape
        left = -rng.random((rows, middle))
        right = -rng.random((middle, columns))
        checks.append(equivalent(ORIGINAL(left, right), batched_product(left, right)))
    left, right = -rng.random((21, 21)), -rng.random((21, 21))
    left[0] = -np.inf
    right[:, -1] = -np.inf
    left[1] = -1e300
    right[2] = -1e300
    for budget in (1, 2 * 32 * right.size, 1024**2):
        WORKSPACE_BYTES = budget
        checks.append(equivalent(ORIGINAL(left, right), batched_product(left, right)))
    WORKSPACE_BYTES = 1024**2
    backing = -rng.random((42, 42))
    checks.append(
        equivalent(
            ORIGINAL(backing[::2, ::2], backing[1::2, 1::2]),
            batched_product(backing[::2, ::2], backing[1::2, 1::2]),
        )
    )
    for retention in (0.0, 0.4, 1.0):
        for multiplicity in (2, 3):
            before = model._birth_death_log_transition(0.2, 0.3, 0.2, 20)
            checks.append(
                equivalent(
                    ORIGINAL(before, pulse(20, retention, multiplicity)),
                    batched_product(before, pulse(20, retention, multiplicity)),
                )
            )
    return checks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "results.json")
    args = parser.parse_args()
    with measurement_lock():
        import nwkit

        assert Path(nwkit.__file__).resolve().parent.parent == SOURCE.resolve()
        results = {
            "environment": environment(),
            "prototype": {
                "workspace_bytes_estimate": WORKSPACE_BYTES,
                "max_row_batch": MAX_BATCH,
                "baseline_function_sha256": hashlib.sha256(
                    inspect.getsource(ORIGINAL).encode()
                ).hexdigest(),
            },
            "edge_checks": edge_checks(),
            "cases": [],
        }
        for bound in (4, 20, 40, 86, 128):
            left = model._birth_death_log_transition(0.2, 0.3, 0.2, bound)
            right = pulse(bound)
            case = measured(
                f"kernel_bound_{bound}",
                partial(ORIGINAL, left, right),
                partial(batched_product, left, right),
                "log-semiring matrix-product kernel; inputs prepared outside timing",
                bound=bound,
            )
            results["cases"].append(case)
            print(case["case"], case["timing"]["speedup"], flush=True)
        for bound in (20, 40, 86, 128):
            instance = count_fixture(maximum=min(18, bound))
            event = model.MultiplicationEvent(2, 0.6, fraction=0.4)
            function = partial(
                instance._log_transitions, [[0.2, 0.3]], 1.0, bound, event
            )

            def transitions_with(implementation, function=function):
                with prototype(implementation):
                    values = function()
                    return np.stack([values[key] for key in sorted(values)])

            case = measured(
                f"complete_log_transitions_bound_{bound}",
                partial(transitions_with, ORIGINAL),
                partial(transitions_with, batched_product),
                "complete _log_transitions including all birth/death kernels and pulse construction",
                bound=bound,
                product_calls=count_calls(function),
            )
            results["cases"].append(case)
            print(case["case"], case["timing"]["speedup"], flush=True)
        for name, families, bound, rate, rare, retention in (
            ("rare_small", 1, 20, 1e-50, True, 0.6),
            ("rare_1", 1, 86, 1e-12, True, 0.6),
            ("rare_16_missing", 16, 86, 1e-12, True, 0.6),
            ("ordinary_control", 16, 40, 0.2, False, 0.6),
            ("zero_retention_control", 1, 86, 1e-12, True, 0.0),
        ):
            instance = count_fixture(
                families, min(bound - 2, 84), rare=rare, missing=True
            )
            event = model.MultiplicationEvent(2, retention, fraction=0.4)
            function = partial(
                instance.family_log_likelihoods, [[rate, 0]], 1.0, bound, event
            )
            case = measured(
                "complete_family_likelihood_" + name,
                partial(call_with, ORIGINAL, function),
                partial(call_with, batched_product, function),
                "complete family_log_likelihoods including emissions, transitions, pruning and ascertainment",
                families=families,
                bound=bound,
                product_calls=count_calls(function),
            )
            results["cases"].append(case)
            print(case["case"], case["timing"]["speedup"], flush=True)
        args.output.write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
