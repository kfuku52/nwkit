"""Dispatch, workspace and repeated-node rounding controls for candidate 02."""

import json
from pathlib import Path
from types import FunctionType

import benchmark as bench
import numpy as np


def main():
    with bench.measurement_lock():
        result = {
            "environment": bench.environment(),
            "controls": [],
            "production_proposal_minimum_requested": 64,
        }
        small_dispatch = bench._candidate_function(32)
        for count in (32,):
            _, process, leaves = bench.fixture(count)

            def baseline(process=process, leaves=leaves):
                return bench.ORIGINAL(process, leaves)

            def candidate(process=process, leaves=leaves):
                return small_dispatch(process, leaves)

            result["controls"].append(
                {
                    "case": "small_dispatch",
                    "minimum_requested_for_this_control": 32,
                    "size": count,
                    "equivalence": bench.compare_arrays(baseline(), candidate()),
                    "timing": bench.paired(baseline, candidate),
                }
            )
        tree, _, leaves = bench.fixture(64)
        process = bench.gt.GaussianTreeProcess(
            tree,
            {
                node: bench.gt.GaussianTransition(0 if index % 7 == 0 else 0.8, 0, 0.4)
                for index, node in enumerate(tuple(tree.traverse("preorder"))[1:])
            },
            bench.gt.GaussianRootPrior("gaussian", 0, 0.6),
            "custom",
        )
        result["controls"].append(
            {
                "case": "zero_slope_fallback",
                "size": 64,
                "equivalence": bench.compare_arrays(
                    bench.ORIGINAL(process, leaves), bench.CANDIDATE(process, leaves)
                ),
                "timing": bench.paired(
                    lambda: bench.ORIGINAL(process, leaves),
                    lambda: bench.CANDIDATE(process, leaves),
                ),
            }
        )
        _, process, leaves = bench.fixture(48, model="ou", parameter=0.43)
        baseline = bench.ORIGINAL(process, leaves)
        result["workspace_checks"] = []
        for budget in (24, 192, bench.WORKSPACE_BYTES):
            previous = bench.WORKSPACE_BYTES
            bench.WORKSPACE_BYTES = budget
            try:
                actual = bench._candidate_function(0)(process, leaves)
                result["workspace_checks"].append(
                    {"budget": budget, **bench.compare_arrays(baseline, actual)}
                )
            finally:
                bench.WORKSPACE_BYTES = previous

        def incorrect_own_blocks(covariance, lca, requested, variances, *args):
            bench._fill_blocks(covariance, lca, requested, variances, *args)
            own = {}
            for output, index in enumerate(requested):
                own.setdefault(index, []).append(output)
            for index, indices in own.items():
                covariance[np.ix_(indices, indices)] = variances[index]

        scope = dict(bench.CANDIDATE.__globals__)
        scope["_fill_blocks"] = incorrect_own_blocks
        initial = FunctionType(bench.CANDIDATE.__code__, scope)
        rng = np.random.default_rng(7128)
        nodes = tuple(tree.traverse("preorder"))
        process = bench.gt.GaussianTreeProcess(
            tree,
            {
                node: bench.gt.GaussianTransition(
                    float(rng.uniform(0.7, 1.2)), 0, float(rng.uniform(0.1, 2))
                )
                for node in nodes[1:]
            },
            bench.gt.GaussianRootPrior("gaussian", 0, 0.6),
            "custom",
        )
        requested = list(nodes) * 2
        rng.shuffle(requested)
        reference = bench.ORIGINAL(process, requested)
        corrected = bench.CANDIDATE(process, requested)
        wrong = initial(process, requested)
        differences = np.abs(reference - wrong)
        result["duplicate_arithmetic_fix"] = {
            "initial_shortcut_differing_entries": int(np.count_nonzero(differences)),
            "initial_shortcut_max_absolute_difference": float(differences.max()),
            "corrected": bench.compare_arrays(reference, corrected),
        }
        output = Path(__file__).with_name("additional-results.json")
        output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
        print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
