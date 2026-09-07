"""Reproduce core tree-operation workloads, with warmup + three-run medians.

Set OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1.
Use --checkout to run this same workload against a different source checkout.
Memory is peak tracemalloc allocation, not process RSS.
"""

import argparse
import json
import os
import platform
import statistics
import sys
import time
import tracemalloc
from pathlib import Path


def measure(function):
    function()
    times = []
    result = None
    for _ in range(3):
        start = time.perf_counter()
        result = function()
        times.append(time.perf_counter() - start)
    return statistics.median(times), result


def peak_allocation(function):
    tracemalloc.start()
    try:
        function()
        return tracemalloc.get_traced_memory()[1] / 1024**2
    finally:
        tracemalloc.stop()


def ladder(tips):
    from ete4 import Tree

    root = Tree()
    current = root
    for index in range(tips - 2):
        current.add_child(name=f"T{index}", dist=1)
        current = current.add_child(dist=1)
    current.add_child(name=f"T{tips - 2}", dist=1)
    current.add_child(name=f"T{tips - 1}", dist=1)
    return root


def tree_workloads():
    from ete4 import Tree

    from nwkit.consensus import _build_consensus_tree
    from nwkit.nwk2table import _sister_branch_id
    from nwkit.rf_distance import robinson_foulds
    from nwkit.util import assign_branch_ids, get_subtree_leaf_name_sets

    first, second = ladder(1600), ladder(1600)
    rf_seconds, rf_value = measure(lambda: robinson_foulds(first, second))
    rf_memory = peak_allocation(lambda: robinson_foulds(first, second))
    star = Tree()
    for index in range(8000):
        star.add_child(name=f"t{index}", dist=1)
    ids = assign_branch_ids(star)
    nodes = list(star.traverse())
    sister_seconds, sisters = measure(
        lambda: [_sister_branch_id(node, ids) for node in nodes]
    )
    expected = [-1] + [
        ids[star.children[1] if node is star.children[0] else star.children[0]]
        for node in star.children
    ]
    assert sisters == expected
    intervals = list(get_subtree_leaf_name_sets(star).values())
    interval_seconds, values = measure(lambda: [tuple(value) for value in intervals])
    assert sum(map(len, values)) == 16000
    n = 400
    all_mask = (1 << n) - 1
    selected = [all_mask ^ ((1 << index) - 1) for index in range(1, n - 1)]
    consensus_seconds, consensus = measure(
        lambda: _build_consensus_tree(
            [f"T{i}" for i in range(n)], all_mask, selected, {}, {}
        )
    )
    assert robinson_foulds(consensus, ladder(n))[0] == 0
    return {
        "rf_1600": {
            "seconds": rf_seconds,
            "peak_python_MiB": rf_memory,
            "result": rf_value,
        },
        "sister_8000": {"seconds": sister_seconds},
        "intervals_8000": {"seconds": interval_seconds},
        "consensus_400": {
            "seconds": consensus_seconds,
            "tips": len(list(consensus.leaves())),
        },
    }


def gaussian_workload():
    import numpy as np

    from nwkit.regress import _profile_covariance_fit

    rng = np.random.default_rng(43)
    n = 200
    response = rng.normal(size=n)
    design = np.column_stack([np.ones(n), rng.normal(size=n)])
    covariance = 0.2 * np.ones((n, n)) + 0.8 * np.eye(n)
    seconds, fit = measure(
        lambda: _profile_covariance_fit(
            response, design, np.zeros(n), [("tree", covariance)], reml=True
        )
    )
    return {
        "seconds": seconds,
        "objective": fit["objective"],
        "beta": fit["beta"].tolist(),
        "beta_covariance": fit["beta_covariance"].tolist(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--checkout", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    sys.path.insert(0, str(args.checkout.resolve()))
    import numpy as np
    import scipy

    from nwkit import __version__

    result = {
        "environment": {
            "platform": platform.platform(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "nwkit": __version__,
            "threads": {
                key: os.environ.get(key)
                for key in (
                    "OPENBLAS_NUM_THREADS",
                    "OMP_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS",
                )
            },
        },
        "tree": tree_workloads(),
        "gaussian_200": gaussian_workload(),
    }
    text = json.dumps(result, indent=2) + "\n"
    print(text, end="")
    if args.output:
        args.output.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    main()
