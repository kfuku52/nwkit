"""Compare complete operations in fresh baseline/current checkout processes."""

import argparse
import dataclasses
import fcntl
import importlib.util
import json
import os
import platform
import statistics
import subprocess
import sys
import tempfile
import time
import tracemalloc
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CASES = (
    "signal",
    "cluster",
    "glmm",
    "projection",
    "diagonal",
    "wgd",
    "dtt",
    "multivariate",
)


def fixture_module(directory):
    path = (
        ROOT / "reviews/performance-candidates-2026-10-22" / directory / "benchmark.py"
    )
    spec = importlib.util.spec_from_file_location("fixture_" + directory, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def operation(case):
    import numpy as np

    from nwkit import disparity, gaussian, phylogenetic_glmm, regress, signal_stats
    from nwkit import wgd_count_model as wgd

    if case == "signal":
        fixture = fixture_module("01-signal-factor-reuse")
        covariance = fixture.covariance_fixture(256, "balanced")
        values = np.random.default_rng(123).normal(size=256)
        observed = signal_stats.k_statistic(covariance, values)

        def run():
            rng = np.random.default_rng(974)
            result = signal_stats.permutation_test(
                covariance, values, np.zeros(256), observed, 199, rng
            )
            return result, rng.bit_generator.state

        return run
    if case == "cluster":
        fixture = fixture_module("03-cluster-score-aggregation")
        frame, terms, kwargs, _ = fixture.fixture(6000, 1500, 5, "event")
        return lambda: regress._fit_model(frame, terms, **kwargs)
    if case == "glmm":
        n = 64
        rng = np.random.default_rng(2203 + n)
        covariance = 0.8 ** np.abs(np.arange(n)[:, None] - np.arange(n))
        covariance += 0.2 * np.eye(n)
        design = np.column_stack([np.ones(n), rng.normal(size=n)])
        values = rng.poisson(np.exp(0.2 + 0.4 * design[:, 1]))
        return lambda: phylogenetic_glmm.fit_phylogenetic_glmm(
            values, design, lambda parameter: covariance, family="poisson"
        )
    if case == "projection":
        fixture = fixture_module("04-uncertainty-contraction")
        args = fixture.fixture(128, 4, 3, "joint")
        return lambda: phylogenetic_glmm._project_predictor_uncertainty(*args)
    if case == "diagonal":
        covariance = np.diag(np.linspace(0.5, 2, 512))
        return lambda: gaussian.factor_covariance(covariance)
    if case == "wgd":
        fixture = fixture_module("08-wgd-combinatorial-reuse")
        model = fixture.tree_fixture(
            [[84, 1, 1], [84, 1, np.nan]], np.array([0.8, 0.9, 1])
        )
        event = wgd.MultiplicationEvent(2, 0.3)
        return lambda: model.family_log_likelihoods(np.array([[1e-8, 0]]), 1, 86, event)
    if case == "multivariate":
        from nwkit.multivariate_pgls import fit_multivariate_pgls

        fixture = fixture_module("10-observed-covariance")
        responses, design, components, options = fixture.fixture(128, 4, 0.9)
        return lambda: fit_multivariate_pgls(responses, design, components, **options)
    fixture = fixture_module("09-dtt-batched-simulation")
    design, factor = fixture.fixture("balanced", 256, 3)
    return lambda: disparity.simulate_curves(design, factor, 128, seed=42, threads=1)


def flatten(value, key="result", output=None):
    import numpy as np

    if output is None:
        output = {}
    if dataclasses.is_dataclass(value):
        value = {
            field.name: getattr(value, field.name)
            for field in dataclasses.fields(value)
        }
    if isinstance(value, dict):
        for name, item in value.items():
            flatten(item, key + "/" + str(name), output)
    elif isinstance(value, (tuple, list)):
        for index, item in enumerate(value):
            flatten(item, key + "/" + str(index), output)
    elif value is None or isinstance(value, int) and abs(value) >= 2**63:
        output[key] = np.asarray(str(value))
    else:
        array = np.asarray(value)
        if array.dtype.kind == "O":
            raise TypeError((key, type(value)))
        output[key] = array
    return output


def worker(args):
    sys.path.insert(0, str(args.checkout))
    import numpy as np
    import scipy

    from nwkit import __version__

    function = operation(args.case)
    start = time.perf_counter()
    function()
    warmup = time.perf_counter() - start
    loops = min(1000, max(1, int(0.05 / max(warmup, 1e-9))))
    start = time.perf_counter()
    for _ in range(loops):
        result = function()
    seconds = (time.perf_counter() - start) / loops
    np.savez_compressed(args.signature, **flatten(result))
    tracemalloc.start()
    function()
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()
    print(
        json.dumps(
            {
                "seconds": seconds,
                "loops": loops,
                "peak_python_bytes": peak,
                "python": platform.python_version(),
                "numpy": np.__version__,
                "scipy": scipy.__version__,
                "nwkit": __version__,
            }
        )
    )


def compare(first_path, second_path, case):
    import numpy as np

    exact, maximum_scaled = True, 0.0
    with np.load(first_path) as first, np.load(second_path) as second:
        assert set(first.files) == set(second.files)
        for key in first.files:
            left, right = first[key], second[key]
            assert left.shape == right.shape, key
            same = (
                np.array_equal(left, right, equal_nan=True)
                if left.dtype.kind in "fciub"
                else np.array_equal(left, right)
            )
            exact &= same
            if left.dtype.kind in "fciu" and case in ("glmm", "cluster"):
                np.testing.assert_allclose(
                    left,
                    right,
                    rtol=2e-4 if case == "glmm" else 1e-11,
                    atol=2e-6 if case == "glmm" else 1e-12,
                    equal_nan=True,
                    err_msg=key,
                )
                finite = np.isfinite(left) & np.isfinite(right)
                if np.any(finite):
                    delta = np.abs(left[finite] - right[finite]) / np.maximum(
                        1, np.abs(left[finite])
                    )
                    maximum_scaled = max(maximum_scaled, float(np.max(delta)))
            elif case == "glmm" and key.endswith("optimizer_message"):
                # Different successful finite-difference optimizer stopping reasons
                # remain visible in the signatures rather than hidden.
                continue
            else:
                assert same, key
    return {"exact": exact, "maximum_scaled_numeric_difference": maximum_scaled}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--checkout", type=Path, default=ROOT)
    parser.add_argument("--case", choices=CASES, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--repeat", type=int, default=7)
    parser.add_argument("--signature", type=Path)
    args = parser.parse_args()
    if args.signature:
        worker(args)
        return
    if args.baseline is None or args.output is None or args.repeat < 3:
        parser.error("provide --baseline, --output and at least three repetitions")
    environment = os.environ.copy()
    for name in (
        "OPENBLAS_NUM_THREADS",
        "OMP_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "MKL_NUM_THREADS",
    ):
        environment[name] = "1"
    environment.pop("PYTEST_PLUGINS", None)
    records, comparisons = {"baseline": [], "current": []}, []
    with Path(tempfile.gettempdir(), "nwkit-performance-adoption.lock").open(
        "a"
    ) as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        with tempfile.TemporaryDirectory(prefix="nwkit-adoption-") as directory:
            paths = {name: Path(directory, name + ".npz") for name in records}
            for repetition in range(args.repeat):
                names = (
                    list(records) if repetition % 2 == 0 else list(reversed(records))
                )
                for name in names:
                    checkout = args.baseline if name == "baseline" else args.checkout
                    command = [
                        sys.executable,
                        str(Path(__file__).resolve()),
                        "--case",
                        args.case,
                        "--checkout",
                        str(checkout.resolve()),
                        "--signature",
                        str(paths[name]),
                    ]
                    completed = subprocess.run(
                        command,
                        capture_output=True,
                        text=True,
                        env=environment,
                        check=True,
                        timeout=300,
                    )
                    records[name].append(json.loads(completed.stdout))
                comparisons.append(
                    compare(paths["baseline"], paths["current"], args.case)
                )
    summary = {
        name: {
            "median_seconds": statistics.median(row["seconds"] for row in rows),
            "min_seconds": min(row["seconds"] for row in rows),
            "max_seconds": max(row["seconds"] for row in rows),
        }
        for name, rows in records.items()
    }
    result = {
        "case": args.case,
        "baseline": str(args.baseline.resolve()),
        "checkout": str(args.checkout.resolve()),
        "environment": platform.platform(),
        "summary": summary,
        "records": records,
        "comparisons": comparisons,
        "median_speedup": summary["baseline"]["median_seconds"]
        / summary["current"]["median_seconds"],
        "contention": "Other user-owned CPU jobs may remain outside the measurement lock",
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {"case": args.case, "speedup": result["median_speedup"], "summary": summary}
        )
    )


if __name__ == "__main__":
    main()
