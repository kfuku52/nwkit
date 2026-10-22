"""Default-draw, random-stream, guard-order and cache-lifetime controls."""

import copy
import json
from functools import partial
from types import SimpleNamespace

import benchmark as bench

asr, history, np = bench.asr, bench.history, bench.np


def recursive_equal(left, right):
    if isinstance(left, np.ndarray):
        np.testing.assert_array_equal(left, right)
    elif isinstance(left, dict):
        assert set(left) == set(right)
        for key in left:
            recursive_equal(left[key], right[key])
    elif isinstance(left, (tuple, list)):
        assert len(left) == len(right)
        for a, b in zip(left, right, strict=True):
            recursive_equal(a, b)
    else:
        assert left == right, (left, right)


def final_states(function):
    original = np.random.default_rng
    generators = []

    def tracked(*args, **kwargs):
        generator = original(*args, **kwargs)
        generators.append(generator)
        return generator

    np.random.default_rng = tracked
    try:
        output = function()
        return output, [copy.deepcopy(x.bit_generator.state) for x in generators]
    finally:
        np.random.default_rng = original


def guard_case(problem, *, histories=False, limit=None, count=1, fit_override=None):
    tree, states, fit = problem
    if fit_override is not None:
        fit = fit_override
    results = []
    for enabled in (False, True):
        with bench.patched(enabled):
            original_context = asr._build_uniformization_context
            context_calls = []

            def forbidden(*args, _calls=context_calls, **kwargs):
                _calls.append(1)
                raise AssertionError(
                    "Context allocation began before the guard rejected work"
                )

            asr._build_uniformization_context = forbidden
            changed = []
            if limit:
                module, name, value = limit
                changed.append((module, name, getattr(module, name)))
                setattr(module, name, value)
            try:
                entry = (
                    history.sample_histories
                    if histories
                    else asr._simulate_stochastic_maps
                )
                try:
                    entry(tree, states, fit, count, seed=8081)
                except Exception as exc:
                    outcome = (type(exc).__name__, str(exc))
                else:
                    raise AssertionError("Expected guard failure")
                assert not context_calls
                results.append(outcome)
            finally:
                asr._build_uniformization_context = original_context
                for module, name, value in changed:
                    setattr(module, name, value)
    assert results[0] == results[1], results
    return {
        "same_exception": True,
        "exception": results[0],
        "contexts_before_failure": 0,
    }


def custom_problem(q, *, length=1.0, target=0):
    tree = bench.read_tree(f"(A:{length})R;", "1", True, quiet=True, rooted="yes")
    labels = tuple(str(i) for i in range(len(q)))
    likelihood = {"A": np.eye(len(q))[target]}
    inside, _, _, transitions = asr._compute_log_inside_likelihoods(tree, likelihood, q)
    root = np.zeros(len(q))
    root[0] = 1
    fit = {
        "rate_matrix": q,
        "log_inside": inside,
        "transition_matrices": transitions,
        "posterior_by_node": {tree: root},
    }
    return tree, labels, fit


def main():
    with bench.measurement_lock():
        unique = bench.fixture(32, unique=True)
        regime = bench.fixture(16, states=4, unique=True, regime=True, projection=True)
        full = [
            bench.measure("unique-count-default100", unique, 100),
            bench.measure("unique-history-default100", unique, 100, histories=True),
        ]
        streams = []
        for histories in (False, True):
            for problem in (unique, regime):
                before, before_states = final_states(
                    partial(bench.call, problem, 4, enabled=False, histories=histories)
                )
                after, after_states = final_states(
                    partial(bench.call, problem, 4, enabled=True, histories=histories)
                )
                same = bench.equivalent(before, after, histories)
                recursive_equal(before_states, after_states)
                streams.append(
                    {
                        "histories": histories,
                        "regime_projection": problem is regime,
                        "comparison": same,
                        "generator_count": len(before_states),
                        "all_final_rng_states_exact": True,
                    }
                )
        guards = {
            "count_work": guard_case(
                unique, limit=(asr, "_MAX_STOCHASTIC_MAP_WORK", 1)
            ),
            "history_work": guard_case(
                unique, histories=True, limit=(asr, "_MAX_STOCHASTIC_MAP_WORK", 1)
            ),
            "history_segments": guard_case(
                unique, histories=True, limit=(history, "_MAX_HISTORY_SEGMENTS", 1)
            ),
            "history_records": guard_case(unique, histories=True, count=100000),
            "history_dense_bytes": guard_case(
                unique,
                histories=True,
                fit_override={"rate_matrix": SimpleNamespace(shape=(1024, 1024))},
            ),
        }
        metadata = []
        for step in range(3):
            if step:
                unique[2]["rate_matrix"] *= 1.1
                for node in unique[0].traverse():
                    if not node.is_root:
                        node.dist *= 1.01
            for histories in (False, True):
                before = bench.setup(unique, 4, enabled=False, histories=histories)
                after = bench.setup(unique, 4, enabled=True, histories=histories)
                recursive_equal(before, after)
                metadata.append(
                    {
                        "mutation_step": step,
                        "histories": histories,
                        "all_prepared_spec_fields_exact": True,
                    }
                )
        zeros = custom_problem(np.zeros((2, 2)), length=0.0)
        q = np.zeros((65, 65))
        for state in range(63):
            q[state, state], q[state, state + 1] = -0.5, 0.5
        q[64, 64], q[64, 63] = -20, 20
        rare = custom_problem(q, target=63)
        edge = []
        for label, problem in (
            ("zero_rate_zero_length", zeros),
            ("rare_endpoint_chain", rare),
        ):
            for histories in (False, True):
                left, left_states = final_states(
                    partial(bench.call, problem, 2, enabled=False, histories=histories)
                )
                right, right_states = final_states(
                    partial(bench.call, problem, 2, enabled=True, histories=histories)
                )
                comparison = bench.equivalent(left, right, histories)
                recursive_equal(left_states, right_states)
                edge.append(
                    {
                        "case": label,
                        "histories": histories,
                        "comparison": comparison,
                        "rng_exact": True,
                    }
                )
        branches = list(unique[0].leaves())
        with bench.patched(True):
            assert asr._validate_stochastic_map_work(unique[2], branches, 1) is None
            assert history._history_work(unique[2], branches, 1) is None
        output = {
            "environment": bench.environment(),
            "default_draw_operations": full,
            "streams": streams,
            "guard_order": guards,
            "metadata_mutation_lifetime": metadata,
            "edges": edge,
            "default_guard_returns_none": True,
        }
        (bench.HERE / "additional-results.json").write_text(
            json.dumps(output, indent=2, allow_nan=False) + "\n"
        )
        for item in full:
            print(item["case"], item["timing"]["speedup"], flush=True)


if __name__ == "__main__":
    main()
