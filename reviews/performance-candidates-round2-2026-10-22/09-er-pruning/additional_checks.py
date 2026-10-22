"""Focused fixed-P, regime, lifetime and downstream seed checks."""

import json
from pathlib import Path

import benchmark as study

np, asr = study.np, study.asr


def run():
    with study.measurement_lock():
        rows = []
        arguments = study.fixture(8, 16, unique=True)
        tree, states, observed, likelihood = arguments
        q = asr._build_rate_matrix("ER", states, [0.02])
        original_p = asr._transition_matrices_for_tree(tree, q)
        for mode in ("original", "one-ulp-P", "ordered-P", "mixed-regime"):
            fixed = {node: matrix.copy() for node, matrix in original_p.items()}
            supplied_q = q
            if mode == "one-ulp-P":
                for matrix in fixed.values():
                    matrix[0, 1] = np.nextafter(matrix[0, 1], np.inf)
            elif mode in {"ordered-P", "mixed-regime"}:
                ordered = np.zeros_like(q)
                for i in range(len(states) - 1):
                    ordered[i, i + 1] = 0.03
                    ordered[i + 1, i] = 0.07
                np.fill_diagonal(ordered, -ordered.sum(axis=1))
                fixed = {
                    node: asr._transition_matrix(ordered, float(node.dist))
                    for node in fixed
                }
                if mode == "mixed-regime":
                    supplied_q = {
                        node: q if i % 2 else ordered
                        for i, node in enumerate(original_p)
                    }
                    fixed = None
            before = study.ORIGINAL(
                tree, likelihood, supplied_q, fixed_transition_matrices=fixed
            )
            after = study.CANDIDATE(
                tree, likelihood, supplied_q, fixed_transition_matrices=fixed
            )
            comparison = study.compare(before, after)
            if mode in {"one-ulp-P", "ordered-P"}:
                assert comparison["all_fields_exact"]
            rows.append({"case": mode, **comparison})
        # A new call must not retain a description of a mutated matrix.
        fixed = original_p
        study.CANDIDATE(tree, likelihood, q, fixed_transition_matrices=fixed)
        for matrix in {id(x): x for x in fixed.values()}.values():
            matrix[0, 1] = np.nextafter(matrix[0, 1], np.inf)
        before = study.ORIGINAL(tree, likelihood, q, fixed_transition_matrices=fixed)
        after = study.CANDIDATE(tree, likelihood, q, fixed_transition_matrices=fixed)
        comparison = study.compare(before, after)
        assert comparison["all_fields_exact"]
        rows.append({"case": "mutated-input-new-call", **comparison})

        failures = []
        for rate in (0.0, -1.0, np.nan, np.inf):
            outcomes = []
            for implementation in (study.ORIGINAL, study.CANDIDATE):
                try:
                    with study.patched(implementation):
                        asr.compute_mk_marginals(*arguments, rate=rate)
                except Exception as error:
                    outcomes.append((type(error).__name__, str(error)))
                else:
                    outcomes.append(("success", ""))
            assert outcomes[0] == outcomes[1]
            failures.append({"rate": str(rate), "outcome": outcomes[0]})

        # Shared fixed rate and branch transitions; compare downstream draw counts.
        arguments = study.fixture(8, 32)
        baseline_posterior, baseline_fit = study.full_call(arguments, study.ORIGINAL)
        candidate_posterior, candidate_fit = study.full_call(arguments, study.CANDIDATE)
        baseline_maps = asr._simulate_stochastic_maps(
            arguments[0], arguments[1], baseline_fit, 20, seed=971, threads=1
        )
        candidate_maps = asr._simulate_stochastic_maps(
            arguments[0], arguments[1], candidate_fit, 20, seed=971, threads=1
        )
        import pandas as pd

        pd.testing.assert_frame_equal(baseline_maps, candidate_maps)
        result = {
            "environment": study.environment(),
            "fixed_P_regime_lifetime": rows,
            "failures": failures,
            "seeded_map_counts": {
                "exact": True,
                "seed": 971,
                "draws": 20,
                "states": 32,
            },
        }
        (Path(__file__).parent / "additional-results.json").write_text(
            json.dumps(result, indent=2, allow_nan=False) + "\n"
        )
        print("5 fixed-P/regime/lifetime controls; 4 failures; seeded maps: passed")


if __name__ == "__main__":
    run()
