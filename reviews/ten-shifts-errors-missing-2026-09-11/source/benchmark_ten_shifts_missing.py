"""Paired noisy incomplete 100-tip, ten-shift native alpha benchmark."""

import argparse
import cProfile
import hashlib
import json
import resource
import signal
import time
from pathlib import Path

import numpy as np
from benchmark_shift_alpha_models import alarm_handler, clean, describe_fit
from benchmark_shift_covariance import balanced_tree

import nwkit.shift_native_search as search_module
from nwkit.shift_native_fit import NativeFitOptions
from nwkit.shift_native_heuristic import NativeSearchOptions, heuristic_native_search
from nwkit.shift_native_model import ShiftData, ShiftLayout
from nwkit.shift_native_provenance import native_implementation_sha256
from nwkit.shift_simulation import simulate_shift
from nwkit.shift_simulation_cli import explicit_simulation


def fixture(job):
    tree = balanced_tree(100)
    p = job.get("traits", 2)
    seed = 921000 + job["replicate"]
    rng = np.random.default_rng(seed)
    # Sixteen disjoint 6/7-tip clades; leave six as observed background.
    eligible = [
        i
        for i in range(1, len(tree.branch_ids))
        if tree.tip_intervals[i][1] - tree.tip_intervals[i][0] in (6, 7)
    ]
    assert len(eligible) == 16
    chosen = sorted(rng.choice(eligible, size=10, replace=False).tolist())
    branches = sorted(tree.branch_ids[i] for i in chosen)
    layout = ShiftLayout.build(tree, branches)
    alpha = np.ones(p) if job["truth"] == "shared" else np.geomspace(0.25, 4, p)
    corr = 0.2 * np.eye(p) + 0.8 * np.ones((p, p))
    marginal = -np.expm1(-2 * alpha) / (2 * alpha)
    diffusion = corr / np.sqrt(marginal[:, None] * marginal[None, :])
    optima = np.zeros((11, p))
    for k, group in enumerate(layout.groups[1:], 1):
        branch = group[0]
        index = tree.branch_ids.index(branch)
        age = tree.remaining_times[tree.compiled.parents[index]]
        direction = rng.normal(size=p)
        direction /= np.linalg.norm(direction)
        optima[k] = 3 * direction / -np.expm1(-alpha * age)
    missing = rng.random((100, p)) < job.get("missing_rate", 0.2)
    parameters = dict(
        trait_names=[f"x{i}" for i in range(p)],
        alpha=alpha.tolist(),
        diffusion_covariance=diffusion.tolist(),
        regime_optima=optima.tolist(),
        shift_branch_ids=branches,
        sampling_standard_errors=[0.1] * p,
        measurement_covariance=(0.04 * np.eye(p)).tolist(),
        missing=missing.tolist(),
    )
    spec = explicit_simulation(tree, parameters)
    values, _ = simulate_shift(spec, seed=seed + 100000)
    data = ShiftData.build(tree, values[0], spec.trait_names, spec.sampling_variances)
    truth = dict(
        parameters=parameters,
        seed=seed,
        data_sha256=hashlib.sha256(values.tobytes()).hexdigest(),
        missing_count=int(missing.sum()),
        tip_names=list(tree.leaf_names),
        true_branches=branches,
    )
    return data, truth


def metrics(branches, truth):
    found, actual = set(branches), set(truth)
    tp, fp, fn = len(found & actual), len(found - actual), len(actual - found)
    return dict(
        tp=tp,
        fp=fp,
        fn=fn,
        precision=tp / len(found) if found else 0.0,
        recall=tp / len(actual),
        f1=2 * tp / (len(found) + len(actual)),
        estimated_shifts=len(found),
        exact_set=found == actual,
    )


def worker(job, path):
    signal.signal(signal.SIGALRM, alarm_handler)
    result = dict(job=job, implementation_sha256=native_implementation_sha256())
    profile = cProfile.Profile() if job.get("profile") else None
    started = time.perf_counter()
    fit_times = []
    original_fit = search_module.fit_native_layout

    def measured_fit(data, layout, **kwargs):
        start = time.perf_counter()
        record = dict(branches=list(layout.shifts))
        try:
            fitted = original_fit(data, layout, **kwargs)
            record["status"] = "complete"
            return fitted
        except Exception:
            record["status"] = "failed_or_interrupted"
            raise
        finally:
            record["seconds"] = time.perf_counter() - start
            fit_times.append(record)
            print(json.dumps(record), flush=True)

    search_module.fit_native_layout = measured_fit
    try:
        data, truth = fixture(job)
        result["generating"] = truth
        options = NativeFitOptions(
            trait_covariance="full",
            alpha_model=job["mode"],
            estimate_measurement_error=True,
        )
        budget = NativeSearchOptions(
            max_shifts=12,
            candidate_pool=32,
            refit_budget=26,
            beam_width=2,
            screening_budget=10000,
            lasso_iterations=200,
        )
        signal.alarm(job.get("timeout", 1800))
        if profile:
            profile.enable()
        t = time.perf_counter()
        search = heuristic_native_search(
            data, options=budget, fit_arguments={"options": options}, criterion="AIC"
        )
        result["search_seconds"] = time.perf_counter() - t
        if profile:
            profile.disable()
        signal.alarm(0)
        retained = [
            describe_fit(data, fit) for fit in search.best_by_complexity.values()
        ]
        selected = {}
        for criterion in ("AIC", "BIC"):
            best = min(retained, key=lambda f: f["criteria"][criterion]["score"])
            selected[criterion] = dict(
                fit=best, **metrics(best["shift_branch_ids"], truth["true_branches"])
            )
        result.update(
            status="complete",
            selected=selected,
            retained=retained,
            records=search.records,
            metadata=search.metadata,
        )
    except TimeoutError as exc:
        result.update(status="timeout", error=str(exc))
    except Exception as exc:
        result.update(status="failed", error=f"{type(exc).__name__}: {exc}")
    finally:
        signal.alarm(0)
        if profile:
            profile.disable()
            profile.dump_stats(str(path.with_suffix(".prof")))
    search_module.fit_native_layout = original_fit
    result["fit_times"] = fit_times
    result["elapsed_seconds"] = time.perf_counter() - started
    result["peak_rss_kib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    path.write_text(json.dumps(clean(result), indent=2, allow_nan=False) + "\n")
    print(
        json.dumps(
            dict(
                job=job,
                status=result["status"],
                elapsed=result["elapsed_seconds"],
                error=result.get("error"),
            )
        ),
        flush=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    worker(json.loads(args.job), args.output)
