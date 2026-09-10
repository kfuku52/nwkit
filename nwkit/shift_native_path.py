"""AIC candidates from covariance-updated paths of OU optimum increments.

The path defines the candidate set, not the final estimates or AIC penalty.
Every retained layout is refitted without shrinkage. Unlike branch-pool beam
search, the path retains joint configurations and does not standardize each
branch column. Consequently the two searches need not select the same model,
and neither certifies a global discrete AIC minimum.
"""

import numpy as np

from nwkit.shift_native_fit import NativeFitOptions
from nwkit.shift_native_heuristic import NativeSearchOptions
from nwkit.shift_native_model import ShiftLayout
from nwkit.shift_native_screen import (
    _kkt_residual,
    _proximal_step,
    _whitened_matrices,
)
from nwkit.shift_native_search import NativeLayoutEvaluator, observable_layout


def _path_layout(data, coefficients, branches, max_shifts):
    magnitudes = np.linalg.norm(coefficients, axis=1)
    active = np.flatnonzero(magnitudes > 1e-8 * max(1.0, float(max(magnitudes))))
    active_count = len(active)
    if active_count > max_shifts:
        # A finite lambda grid can jump across the cap. Retain its strongest
        # effects as an explicitly approximate proposal, then refit without lasso.
        active = sorted(active, key=lambda j: (-magnitudes[j], branches[j]))[
            :max_shifts
        ]
    if not len(active):
        return None, active_count
    try:
        layout = ShiftLayout.build(data.tree, [branches[j] for j in active])
    except ValueError:
        return None, active_count
    return (layout if observable_layout(data, layout) else None), active_count


def _candidate_path(data, fit, options, point_budget):
    matrices, responses, branches = _whitened_matrices(
        data, fit, options.memory_limit, optimum_increments=True
    )
    gradient = np.column_stack(
        [x.T @ y for x, y in zip(matrices, responses, strict=True)]
    )
    maximum = float(np.max(np.linalg.norm(gradient, axis=1)))
    coefficients = np.zeros((len(branches), len(responses)))
    step = 1 / max(1.0, max(float(np.sum(x * x)) for x in matrices))
    layouts, records = [], []
    for fraction in np.geomspace(0.999, 0.005, 80)[:point_budget]:
        strength = maximum * fraction
        for _iteration in range(options.lasso_iterations):
            coefficients, gradient, step = _proximal_step(
                matrices, responses, coefficients, strength, step
            )
            residual = _kkt_residual(coefficients, gradient, strength)
            if residual <= 1e-5 * max(1.0, maximum):
                break
            step *= 1.05
        layout, active = _path_layout(data, coefficients, branches, options.max_shifts)
        records.append(
            {
                "relative_strength": float(fraction),
                "active_branches": active,
                "truncated_at_shift_cap": active > options.max_shifts,
                "iterations": _iteration + 1,
                "kkt_residual": residual,
                "converged": residual <= 1e-5 * max(1.0, maximum),
            }
        )
        if layout is not None and layout not in layouts:
            layouts.append(layout)
        if active > 2 * options.max_shifts:
            break
    return layouts, records


def _initial_covariance(data, evaluator, fit_arguments):
    null_fit = evaluator.best[0]
    options = fit_arguments.get("options", NativeFitOptions())
    if (
        fit_arguments.get("alpha_height") is not None
        or options.root_model != "OUfixedRoot"
    ):
        return null_fit, 0
    # The first path uses the exact fixed-root Brownian limit. This fit supplies
    # covariance only; it is not scored as a fixed-alpha candidate in selection.
    seed = NativeLayoutEvaluator(data, {**fit_arguments, "alpha_height": 0.0})
    seed.evaluate(null_fit["layout"])
    return seed.best[0], 1


def sparse_native_search(data, *, options=None, fit_arguments=None, criterion="AIC"):
    options = NativeSearchOptions() if options is None else options
    options.validate(data, uses_candidate_pool=False)
    if criterion not in {"AIC", "AICc"} or options.convergence:
        raise ValueError("Native path search requires AIC or AICc without convergence.")
    fit_arguments = fit_arguments or {}
    evaluator = NativeLayoutEvaluator(data, fit_arguments, criterion)
    evaluator.evaluate(ShiftLayout.build(data.tree))
    if evaluator.best_information is None:
        raise ValueError("No candidate has finite native information criterion.")
    paths: list[dict] = []
    covariance_fits = 0
    if options.max_shifts:
        fit, covariance_fits = _initial_covariance(data, evaluator, fit_arguments)
        for cycle in range(2):
            remaining_points = options.screening_budget - sum(
                len(record["path"]) for record in paths
            )
            remaining_refits = (
                options.refit_budget - len(evaluator.records) - covariance_fits
            )
            if remaining_points <= 0 or remaining_refits <= 0:
                break
            layouts, records = _candidate_path(data, fit, options, remaining_points)
            fresh = [layout for layout in layouts if layout not in evaluator.scores]
            limit = max(1, remaining_refits // (2 - cycle))
            if len(fresh) > limit:
                indices = np.linspace(0, len(fresh) - 1, limit).round().astype(int)
                fresh = [fresh[i] for i in indices]
            for layout in fresh:
                evaluator.evaluate(layout)
            paths.append(
                {
                    "covariance_source_shifts": list(fit["layout"].shifts),
                    "candidate_layouts": len(layouts),
                    "refitted_candidates": len(fresh),
                    "path": records,
                }
            )
            fit = evaluator.best_information
    points = sum(len(record["path"]) for record in paths)
    return evaluator.finish(
        {
            "strategy": "covariance_updated_optimum_path",
            "complete_discrete_enumeration": False,
            "continuous_global_optimum_certified": False,
            "candidate_parameterization": "unstandardized_OU_optimum_increments",
            "candidate_generation_only": True,
            "covariance_seed_fits": covariance_fits,
            "refitted_candidates": len(evaluator.records),
            "refit_budget": options.refit_budget,
            "screening_evaluations": points,
            "screening_budget": options.screening_budget,
            "budget_exhausted": len(evaluator.records) + covariance_fits
            >= options.refit_budget
            or points >= options.screening_budget,
            "all_paths_converged": all(
                row["converged"] for path in paths for row in path["path"]
            ),
            "paths": paths,
        }
    )
