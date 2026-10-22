"""Failure/lifetime controls and one non-timed explanatory profile."""

import cProfile
import io
import json
import pstats

import benchmark as study


def failure_result(function, arguments, kwargs):
    try:
        function(*arguments, **kwargs)
    except Exception as error:
        return {"type": type(error).__name__, "message": str(error)}
    raise AssertionError("Failure control unexpectedly succeeded")


def failures():
    cases = []
    tree, states, observed, weights = study.fixture(4)
    tree.children[0].dist = 0.0
    cases.append(("zero_branch", (tree, states, observed, weights), {}))
    tree, states, observed, weights = study.fixture(4)
    tree.children[0].dist = -1.0
    cases.append(("negative_branch", (tree, states, observed, weights), {}))
    tree, states, observed, weights = study.fixture(4)
    tree.children[0].dist = float("nan")
    cases.append(("nonfinite_branch", (tree, states, observed, weights), {}))
    cases.append(
        ("threshold_count", study.fixture(4, ordinal=True), {"thresholds": "0"})
    )
    cases.append(
        ("unordered_thresholds", study.fixture(4, ordinal=True), {"thresholds": "1,0"})
    )
    tree, states, observed, weights = study.fixture(4)
    observed = {
        name: None if index else value
        for index, (name, value) in enumerate(observed.items())
    }
    cases.append(("one_observation", (tree, states, observed, weights), {}))
    tree, states, observed, weights = study.fixture(6, ordinal=True)
    observed = {name: "low" for name in observed}
    weights = {name: study.np.array([1.0, 0.0, 0.0]) for name in weights}
    cases.append(("missing_ordinal_category", (tree, states, observed, weights), {}))
    tree, states, observed, weights = study.fixture(4)
    weights["T0"] = study.np.zeros(len(states))
    cases.append(("no_allowed_state", (tree, states, observed, weights), {}))
    results = []
    for name, arguments, kwargs in cases:
        before = failure_result(study.ORIGINAL, arguments, kwargs)
        after = failure_result(study.CANDIDATE, arguments, kwargs)
        assert before == after, (name, before, after)
        results.append({"case": name, "same_failure": before})
    return results


def lifetime():
    arguments = study.fixture(4)
    settings = dict(num_samples=12, burnin=8, thin=1, chains=2, seed=419)
    old = study.ORIGINAL(*arguments, **settings)
    first = study.CANDIDATE(*arguments, **settings)
    study.compare_results(old, first)
    arguments[0].children[0].dist *= 1.7
    new = study.ORIGINAL(*arguments, **settings)
    second = study.CANDIDATE(*arguments, **settings)
    study.compare_results(new, second)
    study.compare_results(second, study.CANDIDATE(*arguments, **settings))
    assert first[1].liability_marginals != second[1].liability_marginals
    return {"new_call_reprepares_changed_branch": True, "repeat_call_exact": True}


def profile():
    arguments = study.fixture(32)
    settings = dict(num_samples=32, burnin=16, thin=1, chains=2, seed=522)
    profiler = cProfile.Profile()
    profiler.runcall(study.ORIGINAL, *arguments, **settings)
    output = io.StringIO()
    stats = pstats.Stats(profiler, stream=output).sort_stats("cumulative")
    stats.print_stats(35)
    study.HERE.joinpath("profile.txt").write_text(output.getvalue())
    selected = []
    for (filename, line, name), (
        _primitive,
        total,
        own,
        cumulative,
        _,
    ) in stats.stats.items():
        if name in {
            "_conditional_parameters",
            "_truncated_draw",
            "diagnose_threshold_draws",
            "rvs",
        }:
            selected.append(
                {
                    "file": filename,
                    "line": line,
                    "function": name,
                    "calls": total,
                    "own_seconds": own,
                    "cumulative_seconds": cumulative,
                }
            )
    return {
        "scope": "one explanatory cProfile call, not a timing comparison",
        "selected_functions": selected,
    }


def main():
    with study.measurement_lock():
        data = {"failures": failures(), "lifetime": lifetime(), "profile": profile()}
        study.HERE.joinpath("additional-results.json").write_text(
            json.dumps(data, indent=2, allow_nan=False) + "\n"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
