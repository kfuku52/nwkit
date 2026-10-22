"""One longer-sample check after the initial full-operation noise."""

import json
from functools import partial

import benchmark as candidate


def main():
    with candidate.measurement_lock():
        rows = []
        for name, families, bound, rare in (
            ("rare_1", 1, 86, True),
            ("rare_16", 16, 86, True),
            ("ordinary_control", 16, 40, False),
        ):
            instance = candidate.count_fixture(
                families, min(bound - 2, 84), rare=rare, missing=True
            )
            event = candidate.model.MultiplicationEvent(2, 0.6, fraction=0.4)
            function = partial(
                instance.family_log_likelihoods,
                [[1e-12 if rare else 0.2, 0]],
                1.0,
                bound,
                event,
            )
            baseline = partial(candidate.call_with, candidate.ORIGINAL, function)
            proposed = partial(candidate.call_with, candidate.batched_product, function)
            row = {
                "case": name,
                "scope": "complete family_log_likelihoods, longer samples",
                "comparison": candidate.equivalent(baseline(), proposed()),
                "timing": candidate.paired(baseline, proposed, minimum_seconds=0.3),
                "product_calls": candidate.count_calls(function),
            }
            rows.append(row)
            print(name, row["timing"]["speedup"], flush=True)
        instance = candidate.count_fixture(maximum=18)
        event = candidate.model.MultiplicationEvent(2, 0.6, fraction=0.4)

        def transition_call(implementation):
            with candidate.prototype(implementation):
                values = instance._log_transitions([[0.2, 0.3]], 1.0, 40, event)
                return candidate.np.stack([values[key] for key in sorted(values)])

        baseline = partial(transition_call, candidate.ORIGINAL)
        proposed = partial(transition_call, candidate.batched_product)
        rows.append(
            {
                "case": "complete_log_transitions_40",
                "comparison": candidate.equivalent(baseline(), proposed()),
                "timing": candidate.paired(baseline, proposed, minimum_seconds=0.3),
            }
        )
        print(rows[-1]["case"], rows[-1]["timing"]["speedup"], flush=True)
        result = {"environment": candidate.environment(), "cases": rows}
        (candidate.HERE / "followup-results.json").write_text(
            json.dumps(result, indent=2, allow_nan=False) + "\n"
        )


if __name__ == "__main__":
    main()
