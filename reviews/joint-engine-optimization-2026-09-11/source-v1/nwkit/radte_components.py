"""Separate input-refit variation from conditional interval widths and variances."""

from copy import deepcopy

import numpy as np
import pandas as pd

CONDITIONAL_COLUMNS = [
    "sample",
    "shared_age_id",
    "estimated_age",
    "interval_lower",
    "interval_upper",
    "conditional_variance",
    "interval_status",
    "method",
    "error",
]
COMPONENT_COLUMNS = [
    "shared_age_id",
    "input_samples",
    "conditional_samples",
    "input_refit_sd",
    "input_refit_variance",
    "mean_conditional_interval_width",
    "mean_conditional_bootstrap_variance",
    "conditional_method",
    "status",
]


def conditional_intervals(
    chronology, fit, problem, args, sample_id, method, *, audit=None
):
    from nwkit.radte_uncertainty import bootstrap_intervals, profile_intervals

    requested = getattr(args, "ensemble_within_uncertainty", "none")
    if requested == "none":
        return []
    result = deepcopy(fit)
    error = ""
    try:
        if requested == "profile":
            profile_intervals(
                result,
                problem,
                level=args.interval_level,
                starts=args.starts,
                maxiter=args.maxiter,
                seed=args.seed,
            )
        else:
            bootstrap_intervals(
                result,
                problem,
                replicates=args.bootstrap_replicates,
                level=args.interval_level,
                rho=args.rate_correlation,
                starts=args.starts,
                maxiter=args.maxiter,
                seed=args.seed,
                rate_sd=args.rate_sd,
            )
    except (ValueError, FloatingPointError) as exc:
        result.interval_lower = result.interval_upper = None
        result.interval_status = "unavailable-conditional-fit-failed"
        error = str(exc)
    if audit is not None:
        attempts = [
            {
                key: None
                if isinstance(value, (float, np.floating)) and not np.isfinite(value)
                else value
                for key, value in attempt.items()
            }
            for attempt in result.attempts[len(fit.attempts) :]
        ]
        audit.append(
            dict(
                sample=sample_id,
                interval_status=result.interval_status,
                error=error,
                diagnostics=result.diagnostics,
                optimizer_attempts=attempts,
            )
        )
    variance = (
        np.var(result.samples * chronology.scale, axis=0, ddof=1)
        if requested == "bootstrap"
        and result.samples is not None
        and len(result.samples) >= 20
        and result.interval_lower is not None
        else None
    )
    return [
        dict(
            sample=sample_id,
            shared_age_id=key,
            estimated_age=float(result.ages[i] * chronology.scale),
            interval_lower=np.nan
            if result.interval_lower is None
            else float(result.interval_lower[i] * chronology.scale),
            interval_upper=np.nan
            if result.interval_upper is None
            else float(result.interval_upper[i] * chronology.scale),
            conditional_variance=np.nan if variance is None else float(variance[i]),
            interval_status=result.interval_status,
            method=method,
            error=error,
        )
        for i, key in enumerate(chronology.groups)
    ]


def component_tables(chronology, fit):
    conditional = pd.DataFrame(fit.conditional_intervals, columns=CONDITIONAL_COLUMNS)
    rows = []
    metadata = fit.ensemble_metadata
    if metadata and fit.samples is not None:
        enough_outer = len(fit.samples) >= max(
            20, int(np.ceil(0.9 * metadata["sample_count"]))
        )
        requested = metadata.get("within_fit_method", "none")
        same_method = len(set(metadata.get("sample_methods", []))) == 1
        for j, key in enumerate(chronology.groups):
            values = fit.samples[:, j] * chronology.scale
            values = values[np.isfinite(values)]
            part = conditional.loc[conditional.shared_age_id == key]
            good = part.loc[
                np.isfinite(part.interval_lower.to_numpy(dtype=float))
                & np.isfinite(part.interval_upper.to_numpy(dtype=float))
            ]
            enough_input = (
                same_method
                and enough_outer
                and len(values) >= max(20, int(np.ceil(0.9 * len(fit.samples))))
            )
            enough_within = enough_input and len(good) >= max(
                20, int(np.ceil(0.9 * len(values)))
            )
            variance = float(np.var(values, ddof=1)) if enough_input else np.nan
            mixed = len(set(good.method)) > 1
            enough_within = enough_within and not mixed
            bootstrap = good.conditional_variance.dropna()
            rows.append(
                dict(
                    shared_age_id=key,
                    input_samples=len(values),
                    conditional_samples=len(good),
                    input_refit_sd=np.sqrt(variance),
                    input_refit_variance=variance,
                    mean_conditional_interval_width=float(
                        (good.interval_upper - good.interval_lower).mean()
                    )
                    if enough_within
                    else np.nan,
                    mean_conditional_bootstrap_variance=float(bootstrap.mean())
                    if enough_within and len(bootstrap) == len(good)
                    else np.nan,
                    conditional_method=requested,
                    status="unavailable-mixed-input-methods"
                    if not same_method
                    else "unavailable-input-coverage"
                    if not enough_input
                    else "input-only"
                    if requested == "none"
                    else "unavailable-mixed-conditional-methods"
                    if mixed
                    else "separate-components"
                    if enough_within
                    else "unavailable-conditional-coverage",
                )
            )
    return {
        "conditional_intervals": conditional,
        "uncertainty_components": pd.DataFrame(rows, columns=COMPONENT_COLUMNS),
    }
