"""Bounded replica batches for DTT; no production sources are edited."""

import argparse
import fcntl
import json
import os
import platform
import statistics
import subprocess
import sys
import time
import tracemalloc
from contextlib import contextmanager
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
LOCK = Path("/tmp/nwkit-perf-20261022-9fYwTD/measurement.lock")
sys.path.insert(0, str(ROOT))
for variable in (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "MKL_NUM_THREADS",
):
    os.environ[variable] = "1"

BATCH_SIZE = 32


def row_squared(delta):
    import numpy as np

    output = []
    for row in delta:
        copied = row.copy()
        output.append(np.dot(copied, copied))
    return np.asarray(output)


def candidate_chunk(payload):
    import numpy as np

    from nwkit.disparity import _simulate_one

    design, factor, seeds = payload
    if len(seeds) <= 1:
        return [(index, _simulate_one(design, factor, seed)) for index, seed in seeds]
    nodes, traits = len(design.parents), factor.shape[0]
    counts = np.zeros(nodes, dtype=int)
    for node in range(nodes - 1, -1, -1):
        counts[node] = (
            1
            if design.tip_rows[node] >= 0
            else sum(counts[child] for child in design.children[node])
        )
    square_roots = np.sqrt(design.lengths)
    indexed = []
    for start in range(0, len(seeds), BATCH_SIZE):
        selected = seeds[start : start + BATCH_SIZE]
        replicas = len(selected)
        innovations = np.empty((nodes - 1, replicas, traits))
        # Same RNG object, draw dimensions and GEMM call as each baseline draw.
        for column, (_, seed) in enumerate(selected):
            rng = np.random.default_rng(seed)
            innovations[:, column, :] = (
                rng.standard_normal((nodes - 1, traits)) @ factor.T
            )
        values = np.zeros((nodes, replicas, traits))
        for node in range(1, nodes):
            values[node] = (
                values[design.parents[node]]
                + square_roots[node] * innovations[node - 1]
            )
        means = np.zeros_like(values)
        sums = np.zeros((nodes, replicas))
        for node in range(nodes - 1, -1, -1):
            if design.tip_rows[node] >= 0:
                means[node] = values[node]
                continue
            merged = 0
            for child in design.children[node]:
                total = merged + counts[child]
                delta = means[child] - means[node]
                squared = row_squared(delta)
                sums[node] += sums[child] + squared * merged * counts[child] / np.full(
                    replicas, total
                )
                means[node] += delta * (counts[child] / total)
                merged = total
        disparities = np.divide(
            2 * sums,
            np.repeat((counts - 1)[:, None], replicas, axis=1),
            out=np.zeros_like(sums),
            where=counts[:, None] > 1,
        )
        if not np.isfinite(disparities).all() or np.any(disparities[0] <= 0):
            raise ValueError(
                "DTT requires finite, positive total disparity; rescale nonconstant traits."
            )
        # Retain the baseline sparse matvec kernel for each replica while
        # concentrating the change in propagation and moment accumulation.
        curves = np.asarray(
            [
                np.asarray(
                    design.weights @ (disparities[:, column] / disparities[0, column])
                )
                for column in range(replicas)
            ]
        )
        indexed.extend(
            (index, curve.copy())
            for (index, _), curve in zip(selected, curves, strict=True)
        )
    return indexed


@contextmanager
def prototype(enabled=True):
    import nwkit.disparity as disparity

    original = disparity._simulate_chunk
    if enabled:
        disparity._simulate_chunk = candidate_chunk
    try:
        yield
    finally:
        disparity._simulate_chunk = original


def fixture(shape, tips, traits):
    import numpy as np
    from ete4 import Tree

    from nwkit.disparity import make_design

    if shape == "star":
        tree = Tree()
        for index in range(tips):
            tree.add_child(name=f"t{index}", dist=1)
    elif shape == "ladder":
        tree = Tree()
        current = tree
        for index in range(tips - 2):
            current.add_child(name=f"t{index}", dist=1)
            current = current.add_child(dist=1)
        current.add_child(name=f"t{tips - 2}", dist=1)
        current.add_child(name=f"t{tips - 1}", dist=1)
    else:
        degree = 3 if shape == "polytomy" else 2
        nodes = []
        for index in range(tips):
            node = Tree()
            node.name = f"t{index}"
            node.dist = 1
            nodes.append(node)
        while len(nodes) > 1:
            parents = []
            for offset in range(0, len(nodes), degree):
                group = nodes[offset : offset + degree]
                if len(group) == 1:
                    parents.extend(group)
                    continue
                parent = Tree()
                parent.dist = 1
                for child in group:
                    parent.add_child(child)
                parents.append(parent)
            nodes = parents
        tree = nodes[0]
    depths = {tree: 0.0}
    for node in tree.traverse("preorder"):
        if not node.is_root:
            depths[node] = depths[node.up] + float(node.dist)
    height = max(depths[leaf] for leaf in tree.leaves())
    for leaf in tree.leaves():
        leaf.dist = height - depths[leaf.up]
    for node in tree.traverse("preorder"):
        if not node.is_root:
            node.dist /= height
            depths[node] = depths[node.up] + float(node.dist)
    design = make_design(tree, sorted(tree.leaf_names()), depths)
    rng = np.random.default_rng(985 + tips + traits)
    factor = np.tril(rng.normal(scale=0.2, size=(traits, traits)))
    factor[np.diag_indices(traits)] = np.linspace(0.7, 1.3, traits)
    return design, factor


def compare(first, second, times):
    import numpy as np

    from nwkit.disparity import curve_area

    difference = float(np.max(np.abs(first - second), initial=0))
    summaries = {}
    for label, function in (
        ("mean", lambda values: np.mean(values, axis=0)),
        ("median", lambda values: np.median(values, axis=0)),
        ("envelope", lambda values: np.quantile(values, [0.025, 0.975], axis=0)),
        ("area", lambda values: curve_area(times, values)),
        ("partial_area", lambda values: curve_area(times, values, (0.13, 0.87))),
    ):
        left, right = function(first), function(second)
        summaries[label] = {
            "exact": bool(np.array_equal(left, right)),
            "max_absolute_difference": float(np.max(np.abs(left - right), initial=0)),
        }
    np.testing.assert_array_equal(first, second)
    return {
        "exact_curves": True,
        "max_absolute_difference": difference,
        "summaries": summaries,
    }


def call(design, factor, count, candidate=False, threads=1):
    from nwkit.disparity import simulate_curves

    with prototype(candidate):
        return simulate_curves(design, factor, count, seed=42, threads=threads)


def peak(function):
    tracemalloc.start()
    try:
        function()
        return tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()


def benchmark(shape, tips, traits, simulations):
    design, factor = fixture(shape, tips, traits)
    baseline = call(design, factor, simulations)
    candidate = call(design, factor, simulations, True)
    equivalence = compare(baseline, candidate, design.times)
    warm = []
    for enabled in (False, True):
        start = time.perf_counter()
        call(design, factor, simulations, enabled)
        warm.append(time.perf_counter() - start)
    batch = max(1, min(1000, int(0.025 / min(warm)) + 1))
    timings = {"baseline": [], "candidate": []}
    for repetition in range(7):
        choices = (False, True) if repetition % 2 == 0 else (True, False)
        for enabled in choices:
            with prototype(enabled):
                from nwkit.disparity import simulate_curves

                start = time.perf_counter()
                for _ in range(batch):
                    result = simulate_curves(design, factor, simulations, seed=42)
                elapsed = (time.perf_counter() - start) / batch
            name = "candidate" if enabled else "baseline"
            timings[name].append(elapsed)
            compare(baseline, result, design.times)
    records = {}
    for name, enabled in (("baseline", False), ("candidate", True)):
        samples = timings[name]
        median = statistics.median(samples)
        records[name] = {
            "samples_seconds": samples,
            "median_seconds": median,
            "min_seconds": min(samples),
            "max_seconds": max(samples),
            "relative_range": (max(samples) - min(samples)) / median,
            "peak_python_bytes": peak(
                lambda enabled=enabled: call(design, factor, simulations, enabled)
            ),
        }
    return {
        "shape": shape,
        "tips": tips,
        "traits": traits,
        "simulations": simulations,
        "nodes": len(design.parents),
        "times": len(design.times),
        "trial_batch_calls": batch,
        "replica_batch_size": BATCH_SIZE,
        "timings": records,
        "speedup": records["baseline"]["median_seconds"]
        / records["candidate"]["median_seconds"],
        "equivalence": equivalence,
    }


def direct_checks(parallel=True):
    import numpy as np

    from nwkit.disparity import _simulate_one, simulate_curves

    checks = []
    for shape, tips, traits, count in (
        ("balanced", 8, 2, 4),
        ("balanced", 16, 1, 33),
        ("balanced", 16, 3, 65),
        ("ladder", 20, 4, 35),
        ("star", 12, 2, 12),
        ("polytomy", 27, 5, 33),
    ):
        design, factor = fixture(shape, tips, traits)
        first = call(design, factor, count)
        second = call(design, factor, count, True)
        record = {
            "shape": shape,
            "tips": tips,
            "traits": traits,
            "count": count,
            "equivalence": compare(first, second, design.times),
        }
        for prefix in (1, 4):
            np.testing.assert_array_equal(
                second[:prefix], call(design, factor, prefix, True)
            )
        if parallel:
            np.testing.assert_array_equal(second, call(design, factor, count, True, 2))
            record["serial_parallel_exact"] = True
        # Check the unbatched per-draw oracle and the final random state.
        seeds = np.random.SeedSequence(42).spawn(count)
        for index in (0, count - 1):
            np.testing.assert_array_equal(
                second[index], _simulate_one(design, factor, seeds[index])
            )
        record["individual_draw_exact"] = True
        checks.append(record)
    design, factor = fixture("balanced", 8, 2)
    for count in (0, 1):
        first = call(design, factor, count)
        second = call(design, factor, count, True)
        np.testing.assert_array_equal(first, second)
        checks.append({"simulation_count": count, "exact": True})
    for label, bad_factor, seed in (
        ("constant", np.zeros((2, 2)), 42),
        ("nonfinite", np.full((2, 2), np.nan), 42),
        ("invalid_seed", factor, -1),
    ):
        outcomes = []
        for enabled in (False, True):
            try:
                with prototype(enabled):
                    simulate_curves(design, bad_factor, 4, seed=seed)
            except Exception as error:
                outcomes.append([type(error).__name__, str(error)])
            else:
                outcomes.append([None])
        assert outcomes[0] == outcomes[1], (label, outcomes)
        checks.append({"case": label, "outcomes": outcomes})
    return checks


def additional_checks():
    import numpy as np
    from ete4 import Tree

    from nwkit.disparity import make_design

    checks = []
    design, factor = fixture("polytomy", 27, 5)
    for label, selected in (
        ("tiny_factor", factor * 1e-140),
        ("large_factor", factor * 1e140),
        ("strided_factor", factor[::-1, ::-1]),
    ):
        first = call(design, selected, 35)
        second = call(design, selected, 35, True)
        checks.append(
            {"case": label, "equivalence": compare(first, second, design.times)}
        )
    tree = Tree(
        "((((A:.25,B:.25):0):.25,C:.5):.5,((D:.4,E:.4):0,F:.4):.6,G:1);",
        parser=1,
    )
    depths = {tree: 0.0}
    for node in tree.traverse("preorder"):
        if not node.is_root:
            depths[node] = depths[node.up] + node.dist
    unusual = make_design(tree, sorted(tree.leaf_names(), reverse=True), depths)
    _, factor = fixture("balanced", 8, 3)
    first = call(unusual, factor, 35)
    second = call(unusual, factor, 35, True)
    checks.append(
        {
            "case": "zero_edges_unary_polytomy",
            "equivalence": compare(first, second, unusual.times),
        }
    )

    # Compare the actual generator calls and post-draw states in serial mode.
    original_rng = np.random.default_rng
    histories = []
    for enabled in (False, True):
        history = []

        class GeneratorProbe:
            def __init__(self, seed, history=history):
                self.generator = original_rng(seed)
                self.history = history

            def standard_normal(self, size):
                values = self.generator.standard_normal(size)
                self.history.append(
                    {
                        "shape": list(size),
                        "post_draw_state": json.loads(
                            json.dumps(self.generator.bit_generator.state)
                        ),
                    }
                )
                return values

        np.random.default_rng = GeneratorProbe
        try:
            call(unusual, factor, 35, enabled)
        finally:
            np.random.default_rng = original_rng
        histories.append(history)
    assert histories[0] == histories[1]
    checks.append(
        {
            "case": "rng_calls_and_final_states",
            "replicas": len(histories[0]),
            "exact": True,
            "draw_shape": histories[0][0]["shape"],
        }
    )
    return checks


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checks-only", action="store_true")
    parser.add_argument("--no-parallel", action="store_true")
    parser.add_argument("--tests", action="store_true")
    args = parser.parse_args()
    LOCK.parent.mkdir(parents=True, exist_ok=True)
    with LOCK.open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if args.tests:
            environment = os.environ.copy()
            environment["PYTHONPATH"] = str(HERE) + os.pathsep + str(ROOT)
            environment["PYTEST_PLUGINS"] = "candidate09_plugin"
            environment["NWKIT_CI_DIAGNOSTICS"] = "0"
            command = [
                sys.executable,
                "tools/check.py",
                "test",
                "--",
                "tests/test_dtt.py",
                "tests/test_dtt_audit.py",
                "-rs",
            ]
            result = subprocess.run(
                command, cwd=ROOT, env=environment, capture_output=True, text=True
            )
            (HERE / "tests.log").write_text(
                " ".join(command) + "\n" + result.stdout + result.stderr
            )
            print(result.stdout + result.stderr)
            raise SystemExit(result.returncode)
        checks = direct_checks(not args.no_parallel)
        checks.extend(additional_checks())
        if args.checks_only:
            (HERE / "direct-checks.json").write_text(
                json.dumps(checks, indent=2) + "\n"
            )
            print("All direct checks passed.")
            return
        import numpy as np
        import scipy

        from nwkit import __version__

        output = {
            "baseline_commit": "d7116f80b9979739c69a5431a5e34b822795408b",
            "scope": "complete serial simulate_curves; setup included, imports/fixtures excluded",
            "environment": {
                "platform": platform.platform(),
                "python": platform.python_version(),
                "numpy": np.__version__,
                "scipy": scipy.__version__,
                "nwkit": __version__,
                "threads": 1,
            },
            "direct_checks": checks,
            "workloads": [],
        }
        for parameters in (
            ("balanced", 8, 2, 4),
            ("balanced", 64, 2, 128),
            ("balanced", 256, 3, 128),
            ("ladder", 128, 1, 96),
            ("star", 128, 2, 96),
            ("polytomy", 81, 4, 96),
            ("balanced", 64, 2, 1),
        ):
            record = benchmark(*parameters)
            output["workloads"].append(record)
            (HERE / "results.json").write_text(json.dumps(output, indent=2) + "\n")
            print(parameters, "speedup", record["speedup"], flush=True)


if __name__ == "__main__":
    main()
