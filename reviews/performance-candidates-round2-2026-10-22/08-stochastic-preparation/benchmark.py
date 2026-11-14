"""Operation-local uniformization-parameter handoff; no production edits."""

import argparse
import hashlib
import inspect
import json
import os
import subprocess
import sys
import textwrap
from contextlib import contextmanager
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = Path("/tmp/nwkit-perf-round2-20261022-rkzoGe/baseline")
for variable in (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "MKL_NUM_THREADS",
):
    os.environ[variable] = "1"
os.environ["NWKIT_CI_DIAGNOSTICS"] = "0"
os.environ["DYLD_FALLBACK_LIBRARY_PATH"] = "/Users/kf/miniforge3-arm64/lib"
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(SOURCE))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from support import environment, measurement_lock, paired, peak_bytes  # noqa: E402

import nwkit  # noqa: E402
from nwkit import asr  # noqa: E402
from nwkit import stochastic_history as history  # noqa: E402
from nwkit.stochastic_map_tables import history_table  # noqa: E402
from nwkit.util import assign_branch_ids, read_tree  # noqa: E402

NAMES = (
    "_validate_stochastic_map_work",
    "_build_uniformization_context",
    "_build_uniformization_contexts",
    "_stochastic_uniformization_contexts",
    "_simulate_stochastic_maps",
)
ORIGINAL = {(asr, name): getattr(asr, name) for name in NAMES}
ORIGINAL.update(
    {
        (history, name): getattr(history, name)
        for name in ("_history_work", "sample_histories")
    }
)


def clone(module, name, changes, append=""):
    source = textwrap.dedent(inspect.getsource(ORIGINAL[module, name]))
    for old, new in changes:
        assert source.count(old) == 1, (name, old)
        source = source.replace(old, new, 1)
    new_name = "_candidate08" + name
    source = source.replace("def " + name + "(", "def " + new_name + "(", 1)
    source += append
    # Actual module globals retain runtime failure-injection/limit monkeypatches.
    exec(compile(source, "<candidate08_" + name + ">", "exec"), module.__dict__)
    return module.__dict__.pop(new_name)


def build_candidates():
    functions = {}
    functions[asr, "_validate_stochastic_map_work"] = clone(
        asr,
        "_validate_stochastic_map_work",
        [
            (
                "def _validate_stochastic_map_work(fit, branches, num_simulations):",
                "def _validate_stochastic_map_work(fit, branches, num_simulations, *, _return_parameters=False):",
            ),
            ("    cached_terms = {}\n", "    cached_terms = {}\n    parameters = {}\n"),
            (
                "            max_n = _uniformization_parameters(matrix, node.dist)[2]",
                "            parameters[key] = _uniformization_parameters(matrix, node.dist)\n            max_n = parameters[key][2]",
            ),
        ],
        "    if _return_parameters:\n        return parameters\n",
    )
    functions[asr, "_build_uniformization_context"] = clone(
        asr,
        "_build_uniformization_context",
        [
            (
                "def _build_uniformization_context(rate_matrix, branch_length):",
                "def _build_uniformization_context(rate_matrix, branch_length, *, _parameters=None):",
            ),
            (
                "    omega, lam, max_n = _uniformization_parameters(rate_matrix, branch_length)",
                "    omega, lam, max_n = (\n        _uniformization_parameters(rate_matrix, branch_length)\n        if _parameters is None else _parameters\n    )",
            ),
        ],
    )
    functions[asr, "_build_uniformization_contexts"] = clone(
        asr,
        "_build_uniformization_contexts",
        [
            (
                "def _build_uniformization_contexts(rate_matrix, branch_lengths):",
                "def _build_uniformization_contexts(rate_matrix, branch_lengths, *, _parameters=None):",
            ),
            (
                "            rate_matrix, branch_length\n",
                "            rate_matrix, branch_length,\n            _parameters=None if _parameters is None else _parameters[(id(rate_matrix), branch_length)],\n",
            ),
        ],
    )
    functions[asr, "_stochastic_uniformization_contexts"] = clone(
        asr,
        "_stochastic_uniformization_contexts",
        [
            (
                "def _stochastic_uniformization_contexts(fit, branches):",
                "def _stochastic_uniformization_contexts(fit, branches, *, _parameters=None):",
            ),
            (
                "                cache[key] = _build_uniformization_context(matrix, node.dist)",
                "                cache[key] = _build_uniformization_context(\n                    matrix, node.dist,\n                    _parameters=None if _parameters is None else _parameters[key],\n                )",
            ),
            (
                '        fit["rate_matrix"], [float(node.dist) for node in branches]\n',
                '        fit["rate_matrix"], [float(node.dist) for node in branches],\n        _parameters=_parameters,\n',
            ),
        ],
    )
    functions[asr, "_simulate_stochastic_maps"] = clone(
        asr,
        "_simulate_stochastic_maps",
        [
            (
                "    _validate_stochastic_map_work(fit, branches, num_simulations)\n    contexts = _stochastic_uniformization_contexts(fit, branches)",
                "    parameters = _validate_stochastic_map_work(\n        fit, branches, num_simulations, _return_parameters=True\n    )\n    contexts = _stochastic_uniformization_contexts(fit, branches, _parameters=parameters)\n    del parameters",
            )
        ],
    )
    functions[history, "_history_work"] = clone(
        history,
        "_history_work",
        [
            (
                "def _history_work(fit, branches, count):",
                "def _history_work(fit, branches, count, *, _return_parameters=False):",
            ),
            (
                "    _validate_stochastic_map_work(fit, branches, count)",
                "    parameters = _validate_stochastic_map_work(fit, branches, count, _return_parameters=True)",
            ),
            (
                "            terms[key] = _uniformization_parameters(matrix, node.dist)[2]",
                "            terms[key] = parameters[key][2]",
            ),
        ],
        "    if _return_parameters:\n        return parameters\n",
    )
    functions[history, "sample_histories"] = clone(
        history,
        "sample_histories",
        [
            (
                "    _history_work(fit, branches, count)\n    contexts = _stochastic_uniformization_contexts(fit, branches)",
                "    parameters = _history_work(fit, branches, count, _return_parameters=True)\n    contexts = _stochastic_uniformization_contexts(fit, branches, _parameters=parameters)\n    del parameters",
            )
        ],
    )
    return functions


CANDIDATE = build_candidates()


@contextmanager
def patched(enabled):
    for (module, name), function in (CANDIDATE if enabled else ORIGINAL).items():
        setattr(module, name, function)
    try:
        yield
    finally:
        for (module, name), function in ORIGINAL.items():
            setattr(module, name, function)


def fixture(tips, states=3, *, unique=False, regime=False, projection=False):
    lengths = iter(np.linspace(0.15, 0.8, 4 * tips))
    leaves = [f"T{i}:{next(lengths) if unique else 0.4:.16g}" for i in range(tips)]
    while len(leaves) > 1:
        leaves = [
            "("
            + ",".join(leaves[i : i + 2])
            + f"):{next(lengths) if unique else 0.4:.16g}"
            if i + 1 < len(leaves)
            else leaves[i]
            for i in range(0, len(leaves), 2)
        ]
    tree = read_tree(leaves[0] + ";", "1", True, quiet=True, rooted="yes")
    labels = tuple(str(i) for i in range(states))
    observed = {
        f"T{i}": None if i % 7 == 0 else labels[i % states] for i in range(tips)
    }
    evidence = {
        name: np.ones(states)
        if value is None
        else np.array([float(value == x) for x in labels])
        for name, value in observed.items()
    }
    _, fit = asr.compute_mk_marginals(tree, labels, observed, evidence, rate=0.45)
    if regime:
        matrices = [fit["rate_matrix"].copy(), fit["rate_matrix"].copy() * 1.7]
        by_node = {
            node: matrices[i % 2]
            for i, node in enumerate(tree.traverse("preorder"))
            if not node.is_root
        }
        inside, _, terms, transitions = asr._compute_log_inside_likelihoods(
            tree, evidence, by_node
        )
        outside = asr._compute_log_outside_likelihoods(
            tree, terms, transitions, fit["root_prior"]
        )
        fit.update(
            rate_matrix_by_node=by_node,
            log_inside=inside,
            transition_matrices=transitions,
            posterior_by_node={
                node: asr._probabilities_from_logs(inside[node] + outside[node])
                for node in tree.traverse()
            },
        )
    if projection:
        fit["state_projection"] = np.arange(states) % 2
        labels = ("0", "1")
    return tree, labels, fit


def call(problem, count, *, enabled, histories=False, threads=1):
    tree, states, fit = problem
    with patched(enabled):
        if histories:
            return history.sample_histories(
                tree, states, fit, count, seed=8081, threads=threads
            )
        return asr._simulate_stochastic_maps(
            tree, states, fit, count, seed=8081, threads=threads
        )


def equivalent(before, after, histories):
    if not histories:
        pd.testing.assert_frame_equal(before, after, check_exact=True)
        return {"all_fields_exact": True, "rows": len(before)}
    assert before.states == after.states
    assert before.branch_ids == after.branch_ids
    assert before.depths == after.depths
    assert before.lengths == after.lengths
    assert len(before.draws) == len(after.draws)
    for left, right in zip(before.draws, after.draws, strict=True):
        assert left.simulation == right.simulation and left.branches == right.branches
        assert left.counts == right.counts
        np.testing.assert_array_equal(left.node_states, right.node_states)
    pd.testing.assert_frame_equal(
        history_table(before), history_table(after), check_exact=True
    )
    return {
        "all_fields_exact": True,
        "draws": len(before.draws),
        "history_rows": len(history_table(before)),
    }


def parameter_calls(function):
    original = asr._uniformization_parameters
    calls = 0

    def counted(*args, **kwargs):
        nonlocal calls
        calls += 1
        return original(*args, **kwargs)

    asr._uniformization_parameters = counted
    try:
        function()
    finally:
        asr._uniformization_parameters = original
    return calls


def setup(problem, count, *, enabled, histories):
    tree, _, fit = problem
    branches = [node for node in tree.traverse("preorder") if not node.is_root]
    with patched(enabled):
        if histories:
            values = history._history_work(
                fit,
                branches,
                count,
                **({"_return_parameters": True} if enabled else {}),
            )
        else:
            values = asr._validate_stochastic_map_work(
                fit,
                branches,
                count,
                **({"_return_parameters": True} if enabled else {}),
            )
        contexts = asr._stochastic_uniformization_contexts(
            fit, branches, **({"_parameters": values} if enabled else {})
        )
        return asr._build_stochastic_map_spec(
            tree, fit, assign_branch_ids(tree), contexts
        )


def measure(name, problem, count, *, histories=False):
    def baseline():
        return call(problem, count, enabled=False, histories=histories)

    def candidate():
        return call(problem, count, enabled=True, histories=histories)

    print(name, flush=True)
    return {
        "case": name,
        "scope": "complete count-table operation"
        if not histories
        else "complete conditional-history sample",
        "draws": count,
        "branches": len(problem[2]["transition_matrices"]),
        "states": problem[2]["rate_matrix"].shape[0],
        "comparison": equivalent(baseline(), candidate(), histories),
        "uniformization_calls": {
            "baseline": parameter_calls(baseline),
            "candidate": parameter_calls(candidate),
        },
        "timing": paired(baseline, candidate, minimum_seconds=0.08),
        "peak_bytes": {
            "baseline": peak_bytes(baseline),
            "candidate": peak_bytes(candidate),
        },
    }


def tests():
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join([str(SOURCE), str(HERE)])
    command = [
        sys.executable,
        "tools/check.py",
        "test",
        "--",
        "-p",
        "candidate08_plugin",
        "tests/test_asr.py",
        "tests/test_stochastic_history.py",
        "tests/test_stochastic_history_audit.py",
        "-k",
        "stochastic or uniformization or bridge or history or segment",
        "-rs",
    ]
    with measurement_lock(), (HERE / "tests.log").open("w") as output:
        output.write(" ".join(command) + "\n")
        output.flush()
        subprocess.run(
            command,
            cwd=SOURCE,
            env=env,
            stdout=output,
            stderr=subprocess.STDOUT,
            check=True,
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tests", action="store_true")
    args = parser.parse_args()
    if args.tests:
        tests()
        return
    with measurement_lock():
        assert Path(nwkit.__file__).resolve().parent.parent == SOURCE.resolve()
        repeated = fixture(32)
        unique = fixture(32, unique=True)
        small = fixture(4)
        regime = fixture(16, states=4, unique=True, regime=True, projection=True)
        cases = [
            measure("small-repeated-count-1", small, 1),
            measure("unique-count-1", unique, 1),
            measure("unique-count-20", unique, 20),
            measure("repeated-count-20", repeated, 20),
            measure("unique-history-1", unique, 1, histories=True),
            measure("unique-history-20", unique, 20, histories=True),
            measure("repeated-history-20", repeated, 20, histories=True),
            measure("regime-hidden-history-4", regime, 4, histories=True),
        ]
        setup_cases = []
        for histories in (False, True):

            def baseline(histories=histories):
                return setup(unique, 20, enabled=False, histories=histories)

            def candidate(histories=histories):
                return setup(unique, 20, enabled=True, histories=histories)

            setup_cases.append(
                {
                    "histories": histories,
                    "timing": paired(baseline, candidate, minimum_seconds=0.08),
                    "uniformization_calls": {
                        "baseline": parameter_calls(baseline),
                        "candidate": parameter_calls(candidate),
                    },
                }
            )
        result = {
            "environment": environment(),
            "target_sources": {
                name: hashlib.sha256(inspect.getsource(value).encode()).hexdigest()
                for (module, name), value in ORIGINAL.items()
            },
            "complete_operations": cases,
            "setup_only": setup_cases,
        }
        (HERE / "results.json").write_text(
            json.dumps(result, indent=2, allow_nan=False) + "\n"
        )
        for case in cases:
            print(
                case["case"],
                case["timing"]["speedup"],
                case["uniformization_calls"],
                flush=True,
            )


if __name__ == "__main__":
    main()
