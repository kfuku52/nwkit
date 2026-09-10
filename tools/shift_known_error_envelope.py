"""Research-only finite nuisance-grid envelope for the known-error no-shift test.

This is not connected to the CLI. It covers only the supplied fitted grid,
not continuous alpha/variance values, and does not address later-stage means.
"""

import numpy as np


def known_error_envelope(
    search, values, seed=1, replicates=199, level=0.05, max_evaluations=None
):
    from nwkit.shift_calibration import validate_calibration_options

    validate_calibration_options(replicates, level)
    if not search.known_error:
        raise ValueError(
            "The research envelope requires known positive observation errors"
        )
    if max_evaluations is not None and (
        isinstance(max_evaluations, bool)
        or not isinstance(max_evaluations, int)
        or max_evaluations < 1
    ):
        raise ValueError("max_evaluations must be a positive integer or None")
    values = np.asarray(values, dtype=float)
    if values.shape != (search.n,) or not np.all(np.isfinite(values)):
        raise ValueError("Finite tip-aligned observations are required")
    family = search.families[0]
    z = search.q @ (values - values.mean())
    best, at = search.profile(z)
    winner = int(family[np.argmax(best[family, 0])])
    observed = float(2 * (best[:, 0].max() - best[winner, 0]))
    fitted_index = int(at[winner, 0])
    # Ordering may depend on the observation, but rejection requires the entire
    # fixed grid. Every remaining point stays in the set after deduplication.
    fitted_variance = search.cache[fitted_index][1]
    order = [
        fitted_index,
        *(i for i, item in enumerate(search.cache) if item[1] == fitted_variance),
        *range(len(search.cache)),
    ]
    unique = []
    seen = set()
    for index in order:
        item = search.cache[index]
        key = (None if item[1] == 0 else item[0], item[1])
        if key not in seen:
            seen.add(key)
            unique.append(index)
    noise = np.random.default_rng(seed).normal(size=(search.d, replicates))
    lower = 0.0
    evaluated = []
    failure = None
    for index in unique:
        if max_evaluations is not None and len(evaluated) == max_evaluations:
            break
        item = search.cache[index]
        try:
            p = search._probability_from_noise(
                np.zeros(search.d), item[2], 1, family, observed, noise
            )
        except (ValueError, np.linalg.LinAlgError) as exc:
            alpha = search.grid[item[0]]
            failure = {
                "alpha_height": None if np.isinf(alpha) else float(alpha),
                "process_variance": float(item[1]),
                "error": str(exc),
            }
            break
        lower = max(lower, p)
        alpha = search.grid[item[0]]
        evaluated.append(
            dict(
                alpha_height=None if np.isinf(alpha) else float(alpha),
                process_variance=float(item[1]),
                probability=p,
            )
        )
        if lower > level:
            break
    complete = len(evaluated) == len(unique)
    result = {
        "method": "research_known_error_finite_grid_envelope",
        "search_backend": type(search).__name__,
        "statistic": observed,
        "level": level,
        "replicates": replicates,
        "seed": seed,
        "p_value_lower_bound": lower,
        "p_value_upper_bound": lower if complete else 1.0,
        "reject": False if lower > level else (True if complete else None),
        "grid_complete": complete,
        "grid_point_count": len(unique),
        "evaluations": evaluated,
        "scope": "No-shift family and declared finite alpha/process-variance grid only; no continuous-parameter guarantee.",
    }

    if failure is not None:
        result["failed_evaluation"] = failure
    return result
