"""Prepare invariant threshold-Gibbs conditionals without changing arithmetic."""

import argparse
import dataclasses
import hashlib
import inspect
import json
import math
import os
import subprocess
import sys
import textwrap
from contextlib import contextmanager
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCE = Path(
    os.environ.get(
        "NWKIT_BENCH_SOURCE", "/tmp/nwkit-perf-round2-20261022-rkzoGe/baseline"
    )
)
for variable in (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "MKL_NUM_THREADS",
):
    os.environ[variable] = "1"
os.environ["NWKIT_CI_DIAGNOSTICS"] = "0"
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(SOURCE))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from support import environment, measurement_lock, paired, peak_bytes  # noqa: E402

import nwkit  # noqa: E402
from nwkit import threshold_asr as module  # noqa: E402
from nwkit.compiled_tree import CompiledTree  # noqa: E402
from nwkit.gaussian_tree import (  # noqa: E402
    GaussianRootPrior,
    GaussianTransition,
    GaussianTreeProcess,
)
from nwkit.util import read_tree  # noqa: E402

ORIGINAL = module.compute_threshold_marginals
ORIGINAL_CONDITIONAL = module._conditional_parameters


def prepare_conditionals(compiled, process):
    """Per-fit metadata; no mutable tree/process is retained after the call."""
    prepared = []
    for index, node in enumerate(compiled.nodes):
        if index == 0:
            assert process.root.variance is not None
            precision = 1.0 / process.root.variance
            parent = (-1, 0.0, 0.0, 1.0)
            root_numerator = process.root.mean / process.root.variance
        else:
            transition = process.transitions[node]
            precision = 1.0 / transition.variance
            parent = (
                compiled.parents[index],
                transition.slope,
                transition.intercept,
                transition.variance,
            )
            root_numerator = 0.0
        children = []
        for child in compiled.children[index]:
            transition = process.transitions[compiled.nodes[child]]
            precision += transition.slope * transition.slope / transition.variance
            children.append(
                (child, transition.slope, transition.intercept, transition.variance)
            )
        variance = 1.0 / precision
        prepared.append(
            (parent, root_numerator, tuple(children), variance, math.sqrt(variance))
        )
    return tuple(prepared)


def conditional_prepared(index, values, prepared):
    parent, root_numerator, children, variance, sd = prepared[index]
    parent_index, slope, intercept, branch_variance = parent
    if parent_index < 0:
        numerator = root_numerator
    else:
        numerator = (slope * values[parent_index] + intercept) / branch_variance
    for child, slope, intercept, branch_variance in children:
        numerator += slope * (values[child] - intercept) / branch_variance
    return numerator * variance, sd


def build_candidate():
    source = textwrap.dedent(inspect.getsource(ORIGINAL))
    old_prepare = "    process = build_threshold_process(tree)\n"
    new_prepare = (
        old_prepare + "    prepared = _prepare_conditionals(compiled, process)\n"
    )
    old_call = """                    mean, sd = _conditional_parameters(
                        int(node_index), values, compiled, process
                    )"""
    new_call = """                    mean, sd = _conditional_prepared(
                        int(node_index), values, prepared
                    )"""
    assert source.count(old_prepare) == source.count(old_call) == 1
    source = source.replace(old_prepare, new_prepare).replace(old_call, new_call)
    namespace = dict(module.__dict__)
    namespace.update(
        _prepare_conditionals=prepare_conditionals,
        _conditional_prepared=conditional_prepared,
    )
    exec(compile(source, "<candidate05_threshold>", "exec"), namespace)
    return namespace["compute_threshold_marginals"], namespace


COMPILED_CANDIDATE, CANDIDATE_GLOBALS = build_candidate()


def candidate(*args, **kwargs):
    # Honour runtime monkeypatches used by existing memory/disk and failure
    # tests. The production integration would not need this prototype wrapper.
    CANDIDATE_GLOBALS.update(module.__dict__)
    return COMPILED_CANDIDATE(*args, **kwargs)


CANDIDATE = candidate


def fixture(size, *, ordinal=False, sparse=False, ambiguous=False, topology="balanced"):
    leaves = [f"T{i}:{0.15 + (i % 9) * 0.13:.12g}" for i in range(size)]
    if topology == "star":
        text = "(" + ",".join(leaves) + ")R;"
    elif topology == "pectinate":
        text = f"({leaves[0]},{leaves[1]})I0:0.27"
        for index, leaf in enumerate(leaves[2:], 1):
            text = f"({text},{leaf})I{index}:0.27"
        text = text.rsplit(":", 1)[0] + ";"
    else:
        current = leaves
        index = 0
        while len(current) > 1:
            following = []
            for offset in range(0, len(current), 2):
                if offset + 1 == len(current):
                    following.append(current[offset])
                else:
                    following.append(
                        f"({current[offset]},{current[offset + 1]})I{index}:0.31"
                    )
                    index += 1
            current = following
        text = current[0].rsplit(":", 1)[0] + ";"
    tree = read_tree(text, "1", True, quiet=True, rooted="yes")
    states = ("low", "mid", "high") if ordinal else ("low", "high")
    observed, likelihoods = {}, {}
    for index, leaf in enumerate(tree.leaves()):
        name = str(leaf.name)
        category = index % len(states)
        observed[name] = states[category]
        likelihoods[name] = np.asarray(
            [float(i == category) for i in range(len(states))]
        )
        if sparse and index >= 2:
            observed[name] = None
            likelihoods[name] = np.ones(len(states))
        if ambiguous and index % 7 == 4:
            observed[name] = "low|high"
            likelihoods[name] = np.asarray([1.0, 0.0, 1.0])
        elif ambiguous and index % 7 == 5:
            observed[name] = None
            likelihoods[name] = np.ones(len(states))
    return tree, states, observed, likelihoods


@contextmanager
def capture(function):
    """Record traces and final chain RNG states outside timed calls."""
    traces, generators = [], []
    original_storage = module.trace_storage
    original_rng = np.random.default_rng

    @contextmanager
    def recorded_storage(*args, **kwargs):
        with original_storage(*args, **kwargs) as array:
            yield array
            traces.append(np.array(array, copy=True))

    def recorded_rng(*args, **kwargs):
        generator = original_rng(*args, **kwargs)
        generators.append(generator)
        return generator

    module.trace_storage = recorded_storage
    np.random.default_rng = recorded_rng
    try:
        yield traces, generators
    finally:
        module.trace_storage = original_storage
        np.random.default_rng = original_rng


def compare_results(first, second):
    before, old_fit = first
    after, new_fit = second
    assert before.keys() == after.keys()
    for node in before:
        np.testing.assert_array_equal(before[node], after[node])
    checked = []
    for field in dataclasses.fields(old_fit):
        left, right = getattr(old_fit, field.name), getattr(new_fit, field.name)
        if isinstance(left, pd.DataFrame):
            pd.testing.assert_frame_equal(left, right, check_exact=True)
        elif isinstance(left, dict):
            assert left == right, field.name
        elif isinstance(left, (float, np.number)):
            np.testing.assert_array_equal(left, right)
        else:
            assert left == right, field.name
        checked.append(field.name)
    return {
        "posterior_exact": True,
        "all_fit_fields_exact": True,
        "fit_fields": checked,
    }


def whole_case(name, arguments, kwargs):
    def baseline():
        return ORIGINAL(*arguments, **kwargs)

    def candidate():
        return CANDIDATE(*arguments, **kwargs)

    with capture(ORIGINAL) as (old_traces, old_rng):
        first = baseline()
    with capture(CANDIDATE) as (new_traces, new_rng):
        second = candidate()
    comparison = compare_results(first, second)
    assert len(old_traces) == len(new_traces) == 1
    np.testing.assert_array_equal(old_traces[0], new_traces[0])
    old_states = [generator.bit_generator.state for generator in old_rng]
    new_states = [generator.bit_generator.state for generator in new_rng]
    assert old_states == new_states
    comparison.update(
        trace_exact=True,
        trace_sha256=hashlib.sha256(old_traces[0].tobytes()).hexdigest(),
        chain_rng_states_exact=True,
        chain_count=len(old_states),
    )
    timing = paired(baseline, candidate)
    left = timing["baseline"]["raw_seconds"]
    right = timing["candidate"]["raw_seconds"]
    return {
        "case": name,
        "tips": len(arguments[2]),
        "observed_tips": sum(item is not None for item in arguments[2].values()),
        "settings": kwargs,
        "scope": "complete compute_threshold_marginals including preparation, sampler and diagnostics",
        "comparison": comparison,
        "timing": timing,
        "candidate_faster_pairs": sum(
            one > two for one, two in zip(left, right, strict=True)
        ),
        "raw_ranges_separate": max(right) < min(left),
        "peak_traced_bytes": {
            "baseline": peak_bytes(baseline),
            "candidate": peak_bytes(candidate),
        },
    }


def helper_checks():
    results = []
    for topology in ("balanced", "pectinate", "star"):
        tree, *_ = fixture(24, topology=topology)
        compiled = CompiledTree.from_tree(tree)
        for general in (False, True):
            if general:
                rng = np.random.default_rng(120)
                process = GaussianTreeProcess(
                    tree,
                    {
                        node: GaussianTransition(
                            float(rng.choice([-1.0, 0.0, 0.5, 1.25])),
                            float(rng.normal()),
                            float(rng.lognormal(0, 2)),
                        )
                        for node in compiled.nodes[1:]
                    },
                    GaussianRootPrior("gaussian", 1.7, 0.3),
                    "test",
                )
            else:
                process = module.build_threshold_process(tree)
            prepared = prepare_conditionals(compiled, process)
            count = 0
            for seed in range(4):
                values = np.random.default_rng(seed).normal(size=len(compiled.nodes))
                for index in range(len(compiled.nodes)):
                    expected = ORIGINAL_CONDITIONAL(index, values, compiled, process)
                    actual = conditional_prepared(index, values, prepared)
                    np.testing.assert_array_equal(expected, actual)
                    count += 1
            results.append(
                {
                    "topology": topology,
                    "general_affine": general,
                    "exact_comparisons": count,
                }
            )
    return results


def conditional_case(size):
    tree, *_ = fixture(size)
    compiled = CompiledTree.from_tree(tree)
    process = module.build_threshold_process(tree)
    prepared = prepare_conditionals(compiled, process)
    values = np.random.default_rng(11).normal(size=len(compiled.nodes))

    def baseline():
        return [
            ORIGINAL_CONDITIONAL(i, values, compiled, process)
            for i in range(len(values))
        ]

    def candidate():
        return [conditional_prepared(i, values, prepared) for i in range(len(values))]

    np.testing.assert_array_equal(baseline(), candidate())
    return {
        "case": f"conditional_sweep_{size}_tips",
        "scope": "all node conditional parameters, excludes preparation and all RNG/draw/diagnostic work",
        "exact": True,
        "timing": paired(baseline, candidate),
        "preparation_timing": paired(
            lambda: prepare_conditionals(compiled, process),
            lambda: prepare_conditionals(compiled, process),
        ),
        "prepared_peak_bytes": peak_bytes(
            lambda: prepare_conditionals(compiled, process)
        ),
    }


def run_tests():
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join((str(SOURCE), str(HERE), str(HERE.parent)))
    command = [
        sys.executable,
        "tools/check.py",
        "test",
        "--",
        "-p",
        "candidate05_plugin",
        "tests/test_threshold_asr.py",
        "tests/test_threshold_reference.py",
        "-rs",
    ]
    print("COMMAND:", " ".join(command), flush=True)
    return subprocess.run(command, cwd=SOURCE, env=env, check=False).returncode


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tests", action="store_true")
    parser.add_argument("--output", type=Path, default=HERE / "results.json")
    arguments = parser.parse_args()
    with measurement_lock():
        assert Path(nwkit.__file__).resolve().parent.parent == SOURCE.resolve()
        if arguments.tests:
            return run_tests()
        results = {
            "environment": environment(),
            "source_sha256": hashlib.sha256(
                Path(module.__file__).read_bytes()
            ).hexdigest(),
            "helper_checks": helper_checks(),
            "whole_operations": [],
            "conditional_kernels": [],
        }
        cases = [
            (
                "binary4_observed",
                fixture(4),
                dict(num_samples=32, burnin=16, thin=1, chains=2, seed=521),
            ),
            (
                "binary32_observed",
                fixture(32),
                dict(num_samples=32, burnin=16, thin=1, chains=2, seed=522),
            ),
            (
                "binary128_sparse",
                fixture(128, sparse=True),
                dict(num_samples=32, burnin=64, thin=8, chains=2, seed=523),
            ),
            (
                "ordinal24_estimated_ambiguous",
                fixture(24, ordinal=True, ambiguous=True),
                dict(num_samples=32, burnin=24, thin=1, chains=2, seed=527),
            ),
            (
                "ordinal24_fixed_ambiguous",
                fixture(24, ordinal=True, ambiguous=True),
                dict(
                    num_samples=32,
                    burnin=24,
                    thin=1,
                    chains=2,
                    seed=527,
                    thresholds="0,1",
                ),
            ),
        ]
        for name, inputs, settings in cases:
            result = whole_case(name, inputs, settings)
            results["whole_operations"].append(result)
            arguments.output.write_text(
                json.dumps(results, indent=2, allow_nan=False) + "\n"
            )
            print(
                name,
                result["timing"]["speedup"],
                result["comparison"]["trace_exact"],
                flush=True,
            )
        for size in (4, 32, 128):
            result = conditional_case(size)
            results["conditional_kernels"].append(result)
        arguments.output.write_text(
            json.dumps(results, indent=2, allow_nan=False) + "\n"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
