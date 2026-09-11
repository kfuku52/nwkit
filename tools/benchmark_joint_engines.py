"""Matched fixed-parameter joint GLS workloads for engine optimization review.

Run baseline with frozen pre-change NWKIT, optimized with the candidate source.
All inputs are deterministic. Numerical outputs accompany timings for comparison.
"""

import argparse
import json
import resource
import time
from pathlib import Path

import numpy as np

from nwkit.shift_joint_model import evaluate_joint
from nwkit.shift_native_model import ShiftData, ShiftLayout
from nwkit.shift_native_provenance import native_implementation


def workload(tips, traits):
    # Native imports above must precede legacy helpers that prepend their
    # checkout path. Preserve the explicitly selected installed/frozen package.
    from benchmark_shift_covariance import balanced_tree
    from benchmark_ten_shifts_missing import fixture

    if tips == 100 and traits == 2:
        data, truth = fixture(
            dict(truth="different", replicate=999, traits=2, missing_rate=0.2)
        )
        layout = ShiftLayout.build(data.tree, truth["true_branches"])
    else:
        tree = balanced_tree(tips)
        rng = np.random.default_rng(62100 + tips + traits)
        values = rng.normal(size=(tips, traits))
        values[rng.random(values.shape) < 0.2] = np.nan
        data = ShiftData.build(
            tree,
            values,
            tuple(f"x{i}" for i in range(traits)),
            np.full(values.shape, 0.01),
        )
        layout = ShiftLayout.build(tree)
    S = (
        (0.2 * np.eye(traits) + 0.8 * np.ones((traits, traits)))
        / data.scales[:, None]
        / data.scales[None]
    )
    return (
        data,
        layout,
        np.geomspace(0.25, 4, traits),
        S,
        np.full(traits, 0.04) / data.scales**2,
    )


def measure(function):
    for _ in range(3):
        function()
    samples = []
    for _ in range(7):
        start = time.perf_counter()
        for _ in range(5):
            function()
        samples.append((time.perf_counter() - start) / 5)
    return samples


def main():
    from benchmark_shift_alpha_models import clean

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--engine", choices=["baseline", "optimized"], required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    records = []
    for tips, traits in [(100, 2), (100, 5), (1000, 2), (1000, 5)]:
        data, layout, alpha, S, noise = workload(tips, traits)
        for root in ["OUfixedRoot", "OUrandomRoot"]:
            setup = time.perf_counter()
            if args.engine == "optimized":
                from nwkit.shift_joint_dense import DenseJointContext, dense_eligible
                from nwkit.vector_whitening import VectorObservationPlan

                context = (
                    DenseJointContext(data, root) if dense_eligible(data) else None
                )
                plan = VectorObservationPlan.build(
                    data.tree.compiled, np.isfinite(data.values)
                )

                def evaluate(
                    data=data,
                    layout=layout,
                    alpha=alpha,
                    S=S,
                    noise=noise,
                    root=root,
                    context=context,
                    plan=plan,
                ):
                    if context is not None:
                        return context.evaluate(layout, alpha, S, noise)
                    return evaluate_joint(
                        data, layout, alpha, S, noise, root_model=root, _plan=plan
                    )
            else:

                def evaluate(
                    data=data, layout=layout, alpha=alpha, S=S, noise=noise, root=root
                ):
                    return evaluate_joint(
                        data, layout, alpha, S, noise, root_model=root
                    )

            setup_seconds = time.perf_counter() - setup
            fit = evaluate()
            records.append(
                dict(
                    tips=tips,
                    traits=traits,
                    root=root,
                    engine=fit.engine,
                    observed=fit.num_observations,
                    setup_seconds=setup_seconds,
                    seconds=measure(evaluate),
                    log_likelihood=fit.log_likelihood,
                    coefficients=fit.coefficients,
                    coefficient_covariance=fit.coefficient_covariance,
                    predicted=fit.predicted,
                )
            )
            print(tips, traits, root, records[-1]["seconds"], flush=True)
    args.output.write_text(
        json.dumps(
            clean(
                dict(
                    implementation=native_implementation(),
                    records=records,
                    peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                )
            ),
            indent=2,
            allow_nan=False,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
