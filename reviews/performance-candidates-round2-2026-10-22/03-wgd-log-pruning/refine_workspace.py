"""Compare smaller family tiles after observing the initial memory penalty."""

import json
from functools import partial
from pathlib import Path

import benchmark as candidate
from support import environment, measurement_lock


def run():
    import numpy as np

    from nwkit.wgd_count_model import CountLikelihood

    baseline = CountLikelihood._prune_log
    prototype = candidate.build_candidate()
    results = {"environment": environment(), "cases": []}
    for families, bound in ((128, 86), (512, 86), (128, 160)):
        counts = candidate.rare_counts(families, maximum=84 if bound == 86 else 150)
        model = candidate.tree_fixture(counts)
        rates = [[1e-8, 0]] if bound == 86 else [[1e-5, 0]]
        for budget in (256 * 1024, 1024**2, 4 * 1024**2):
            candidate.WORKSPACE_BYTES = budget
            record = candidate.measure_case(
                f"families{families}_bound{bound}_budget{budget}",
                partial(candidate.full_call, model, rates, 1, bound, None, baseline),
                partial(candidate.full_call, model, rates, 1, bound, None, prototype),
                {
                    "families": families,
                    "max_count": bound,
                    "workspace_estimate_bytes": budget,
                },
            )
            assert record["equivalence"]["array_equal"]
            assert np.isfinite(record["timing"]["speedup"])
            results["cases"].append(record)
    Path(__file__).with_name("workspace-results.json").write_text(
        json.dumps(results, indent=2, allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    with measurement_lock():
        run()
