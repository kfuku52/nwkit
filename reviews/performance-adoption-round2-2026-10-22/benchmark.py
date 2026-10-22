"""Compare complete operations in fresh baseline/current checkout processes."""

import argparse
import dataclasses
import fcntl
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
    "lambda",
    "pca",
    "pgls",
    "covariance",
    "wgd",
    "fullou",
    "maps",
    "histories",
    "er",
    "er-free-control",
    "pca-search-control",
    "maps-default-control",
    "covariance-zero-control",
)


def operation(case):
    import fixtures
    import numpy as np

    from nwkit import asr, full_ou, phylogenetic_pca, signal_stats, stochastic_history
    from nwkit import wgd_count_model as wgd
    from nwkit.ordinary_regression import _fit_ordinary_gaussian

    if case == "lambda":
        arguments = fixtures.signal_fixture(256)
        return lambda: signal_stats.lambda_fit(*arguments)
    if case in ("pca", "pca-search-control"):
        arguments = fixtures.pca_fixture(
            512 if case == "pca" else 64, 16 if case == "pca" else 4
        )
        options = {"mode": "corr"} if case == "pca" else {"model": "LAMBDA"}
        return lambda: phylogenetic_pca.fit_pca(*arguments, **options)
    if case in ("covariance", "covariance-zero-control"):
        tree, process, leaves = fixtures.gaussian_fixture(
            512, model="ou", parameter=0.4
        )
        if case == "covariance-zero-control":
            from nwkit.gaussian_tree import GaussianTransition, GaussianTreeProcess

            process = GaussianTreeProcess(
                tree,
                {
                    node: GaussianTransition(0 if i % 7 == 0 else 0.8, 0, 0.4)
                    for i, node in enumerate(tuple(tree.traverse("preorder"))[1:])
                },
                process.root,
                "custom",
            )
        return lambda: process.covariance(leaves)
    if case == "pgls":
        tree, _, leaves = fixtures.gaussian_fixture(256)
        names = [str(leaf.name) for leaf in leaves]
        rng = np.random.default_rng(379)
        design = np.column_stack([np.ones(256), rng.normal(size=256)])
        values = design @ np.array([0.5, -0.8]) + rng.normal(size=256)
        return lambda: _fit_ordinary_gaussian(
            values,
            design,
            np.zeros(256),
            tree,
            names,
            evolution_model="brownian",
            evolution_parameter=None,
            branch_length="original",
            custom_covariance=None,
            reml=True,
        )
    if case == "wgd":
        model = fixtures.count_fixture(fixtures.rare_counts(128))
        event = wgd.MultiplicationEvent(2, 0.3)
        return lambda: model.family_log_likelihoods(np.array([[1e-8, 0]]), 1, 86, event)
    if case == "fullou":
        parameters = fixtures.ou_parameters(8)
        design, complete = fixtures.ou_design(8, 32)
        return lambda: full_ou.full_ou_identifiability(
            parameters, 8, 1, design, complete
        )
    if case in ("maps", "histories", "maps-default-control"):
        tree, states, fit = fixtures.map_fixture(32, unique=True)
        if case != "histories":
            count = 100 if case == "maps-default-control" else 1
            return lambda: asr._simulate_stochastic_maps(
                tree, states, fit, count, seed=8081
            ).to_dict("list")

        def history():
            sample = stochastic_history.sample_histories(
                tree, states, fit, 1, seed=8081
            )
            nodes = list(tree.traverse("preorder"))
            return (
                sample.states,
                sample.draws,
                [sample.depths[n] for n in nodes],
                [sample.lengths[sample.branch_ids[n]] for n in nodes[1:]],
            )

        return history
    arguments = fixtures.er_fixture(
        128 if case == "er" else 32, 64 if case == "er" else 16
    )
    tree = arguments[0]
    nodes = list(tree.traverse("preorder"))
    options = {"rate": 0.15 / 64} if case == "er" else {}

    def er():
        posterior, fit = asr.compute_mk_marginals(*arguments, model="ER", **options)
        matrices = {"inside", "log_inside", "transition_matrices", "posterior_by_node"}
        result = {key: value for key, value in fit.items() if key not in matrices}
        result["posterior"] = np.array([posterior[n] for n in nodes])
        result["inside"] = np.array([fit["inside"][n] for n in nodes])
        result["log_inside"] = np.array([fit["log_inside"][n] for n in nodes])
        result["transitions"] = np.array(
            [fit["transition_matrices"][n] for n in nodes[1:]]
        )
        return result

    return er


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
                "source": __import__("nwkit").__file__,
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
            if left.dtype.kind in "fciu" and case in ("lambda", "er"):
                tolerance = (
                    1e-5
                    if case == "lambda" and key.endswith("estimate")
                    else 1e-7
                    if case == "lambda"
                    else 1e-13
                )
                np.testing.assert_allclose(
                    left, right, rtol=0, atol=tolerance, equal_nan=True, err_msg=key
                )
                finite = np.isfinite(left) & np.isfinite(right)
                if np.any(finite):
                    delta = np.abs(left[finite] - right[finite]) / np.maximum(
                        1, np.abs(left[finite])
                    )
                    maximum_scaled = max(maximum_scaled, float(np.max(delta)))
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
                        check=False,
                        timeout=300,
                    )
                    if completed.returncode:
                        raise RuntimeError(completed.stderr)
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
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(
        json.dumps(
            {"case": args.case, "speedup": result["median_speedup"], "summary": summary}
        )
    )


if __name__ == "__main__":
    main()
