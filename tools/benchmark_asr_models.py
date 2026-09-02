#!/usr/bin/env python3
"""Reproducible ASR timing and equivalent-output benchmark.

Run from a checkout root.  The JSON includes robust medians, peak traced Python
memory, and numerical digests so two revisions can be compared without retaining
large posterior objects.
"""

import argparse
import json
import math
import statistics
import time
import tracemalloc

import numpy as np

from nwkit.asr_regimes import RegimeAssignment
from nwkit.continuous_asr import compute_bm_marginals
from nwkit.ou_asr import compute_ou_marginals
from nwkit.regime_continuous_asr import compute_bms_marginals
from nwkit.util import read_tree


def _tree_text(num_tips):
    rng = np.random.default_rng(8128)
    nodes = [f"T{index}:{rng.uniform(0.2, 1.8):.12g}" for index in range(num_tips)]
    internal = 0
    while len(nodes) > 1:
        next_nodes = []
        for index in range(0, len(nodes), 2):
            if index + 1 == len(nodes):
                next_nodes.append(nodes[index])
                continue
            length = rng.uniform(0.2, 1.8)
            next_nodes.append(
                f"({nodes[index]},{nodes[index + 1]})N{internal}:{length:.12g}"
            )
            internal += 1
        nodes = next_nodes
    return nodes[0].rsplit(":", 1)[0] + "R;"


def _inputs(num_tips):
    tree = read_tree(_tree_text(num_tips), "1", True, quiet=True, rooted="yes")
    rng = np.random.default_rng(271828)
    values = {str(leaf.name): float(rng.normal()) for leaf in tree.leaves()}
    nodes = list(tree.traverse())
    labels = ["slow" if index % 3 else "fast" for index in range(len(nodes))]
    assignment = RegimeAssignment(
        tuple(dict.fromkeys(labels)), dict(zip(nodes, labels, strict=True)), "memory"
    )
    return tree, values, assignment


def _digest(result):
    posterior, fit = result
    likelihood = getattr(fit, "restricted_log_likelihood", None)
    if likelihood is None:
        likelihood = getattr(fit, "log_likelihood", math.nan)
    return [
        float(math.fsum(item.mean for item in posterior.values())),
        float(math.fsum(item.variance for item in posterior.values())),
        float(likelihood),
    ]


def _measure(function, repeats):
    function()
    timings = []
    peak = 0
    digest = None
    for _ in range(repeats):
        tracemalloc.start()
        start = time.perf_counter()
        result = function()
        timings.append(time.perf_counter() - start)
        _, current_peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        peak = max(peak, current_peak)
        digest = _digest(result)
    return {
        "median_seconds": statistics.median(timings),
        "min_seconds": min(timings),
        "peak_mib": peak / 1024.0 / 1024.0,
        "digest": digest,
    }


def _fixed_oum(tree, values, assignment):
    theta = {"slow": -0.3, "fast": 0.6}
    try:
        from nwkit.regime_gaussian_asr import compute_regime_ou_marginals
    except ImportError:
        from nwkit.regime_continuous_asr import compute_oum_marginals

        return compute_oum_marginals(
            tree,
            values,
            assignment,
            alpha=0.7,
            sigma2=1.1,
            theta_by_regime=theta,
            _tree_validated=True,
        )
    return compute_regime_ou_marginals(
        tree,
        values,
        assignment,
        model="OUM",
        alpha=0.7,
        sigma2=1.1,
        regime_parameters={regime: {"theta": value} for regime, value in theta.items()},
        _tree_validated=True,
    )


def _mixture_cache_benchmark(num_tips, repeats):
    try:
        from nwkit.asr import (
            _log_likelihood,
            _transition_matrices_for_tree,
        )
    except ImportError:
        return None
    tree, _, _ = _inputs(num_tips)
    rng = np.random.default_rng(314159)
    characters = []
    for _ in range(16):
        characters.append(
            {str(leaf.name): np.eye(3)[rng.integers(0, 3)] for leaf in tree.leaves()}
        )
    matrix = np.full((3, 3), 0.4, dtype=float)
    np.fill_diagonal(matrix, -0.8)
    rates = np.asarray([0.2, 0.65, 1.25, 2.1])
    prior = np.full(3, 1.0 / 3.0)
    cached = [_transition_matrices_for_tree(tree, matrix * rate) for rate in rates]

    def score(use_cache):
        total = 0.0
        for likelihoods in characters:
            for rate, transitions in zip(rates, cached, strict=True):
                total += _log_likelihood(
                    tree,
                    likelihoods,
                    prior,
                    matrix * rate,
                    **({"fixed_transition_matrices": transitions} if use_cache else {}),
                )
        return total

    expected = score(False)
    observed = score(True)
    if not math.isclose(expected, observed, rel_tol=0.0, abs_tol=1e-12):
        raise RuntimeError("Cached and uncached Mk likelihoods differ.")

    def timing(use_cache):
        values = []
        for _ in range(repeats):
            start = time.perf_counter()
            score(use_cache)
            values.append(time.perf_counter() - start)
        return statistics.median(values)

    uncached = timing(False)
    cached_time = timing(True)
    return {
        "uncached_median_seconds": uncached,
        "cached_median_seconds": cached_time,
        "speedup": uncached / cached_time,
        "equivalent_score": observed,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tips", type=int, default=1024)
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--mixture-tips", type=int, default=128)
    args = parser.parse_args()
    tree, values, assignment = _inputs(args.tips)
    models = {
        "BM": lambda: compute_bm_marginals(
            tree, values, sigma2=0.9, _tree_validated=True
        ),
        "OU": lambda: compute_ou_marginals(
            tree,
            values,
            alpha=0.7,
            sigma2=1.1,
            theta=0.2,
            _tree_validated=True,
        ),
        "BMS": lambda: compute_bms_marginals(
            tree,
            values,
            assignment,
            sigma2_by_regime={"slow": 0.5, "fast": 1.4},
            _tree_validated=True,
        ),
        "OUM": lambda: _fixed_oum(tree, values, assignment),
    }
    output = {
        "tips": args.tips,
        "repeats": args.repeats,
        "models": {
            name: _measure(function, args.repeats) for name, function in models.items()
        },
    }
    mixture = _mixture_cache_benchmark(args.mixture_tips, args.repeats)
    if mixture is not None:
        output["mk_mixture_transition_cache"] = mixture
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
