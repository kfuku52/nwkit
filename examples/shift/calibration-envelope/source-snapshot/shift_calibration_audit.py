"""Reconstruct calibration selection metrics instead of trusting saved flags."""

from collections import defaultdict

import numpy as np
from shift_simulation_cases import rate

from nwkit.shift_calibration import CalibratedSearch
from nwkit.shift_candidates import tip_groups
from nwkit.util import read_tree


def check_probability_metadata(test, index, search, B, level):
    p = test["p_value"]
    if index or search.known_error:
        if test["p_value_kind"] != "plugin" or test["p_value_lower_bound"] != p:
            raise ValueError("Invalid plug-in probability metadata")
        return
    evaluated = test["null_alpha_evaluations"]
    grid = [None if np.isinf(a) else float(a) for a in search.grid]
    if (
        not evaluated
        or [r["alpha_height"] for r in evaluated] != grid[: len(evaluated)]
    ):
        raise ValueError("Invalid null nuisance evaluation grid")
    values = [r["p_value"] for r in evaluated]
    if any(
        not 0 < value <= 1 or abs(value * (B + 1) - round(value * (B + 1))) > 1e-10
        for value in values
    ):
        raise ValueError("Invalid nuisance Monte Carlo probability")
    lower = max(values)
    if test["p_value_lower_bound"] != lower:
        raise ValueError("Incorrect probability lower bound")
    if len(evaluated) == len(grid):
        if test["p_value_kind"] != "grid_supremum" or p != lower:
            raise ValueError("Incorrect full-grid supremum probability")
    elif test["p_value_kind"] != "conservative_upper_bound" or p != 1 or lower <= level:
        raise ValueError("Incorrect early-acceptance probability bound")


def audit_record(row, design, search, *, replay_bootstrap=False):
    fit, truth = row["fit"], row["truth"]
    B = design["calibration_replicates"]
    level = design["calibration_level"]
    if (
        fit["calibration_replicates"] != B
        or fit["calibration_level"] != level
        or fit["seed"] != row["case"]["seed"] + 9000000000
    ):
        raise ValueError("Fit calibration options disagree with protocol")
    winner = fit["winner"]
    if (
        not 0 <= winner < len(search.models)
        or fit["model"] != search.models[winner]
        or fit["candidate_count"] != len(search.models)
    ):
        raise ValueError("Invalid selected candidate")
    values = np.asarray(truth["observations"])
    best, at = search.profile(search.q @ (values - values.mean()))
    best = best[:, 0]
    tests = fit["tests"]
    for index, test in enumerate(tests):
        if index >= len(search.families) - 1 or test["stage"] != index:
            raise ValueError("Invalid calibration stage sequence")
        family = search.families[index]
        statistic = 2 * (best.max() - best[family].max())
        if test["candidate_count"] != len(family) or not np.isclose(
            test["statistic"], statistic, rtol=1e-9, atol=1e-8
        ):
            raise ValueError(
                "Calibration statistic does not match refitted candidate profile"
            )
        p = test["p_value"]
        if not 0 < p <= 1 or abs(p * (B + 1) - round(p * (B + 1))) > 1e-10:
            raise ValueError("Invalid plus-one Monte Carlo probability")
        check_probability_metadata(test, index, search, B, level)
        if index < len(tests) - 1 and p > level:
            raise ValueError("Calibration continued after an accepted family")
    if tests and tests[-1]["p_value"] > level:
        selected_family = search.families[len(tests) - 1]
    elif len(tests) == len(search.families) - 1:
        selected_family = search.families[-1]
    else:
        raise ValueError("Calibration stopped before testing all required families")
    if winner not in selected_family or not np.isclose(
        best[winner], best[selected_family].max(), rtol=1e-9, atol=1e-8
    ):
        raise ValueError("Winner does not maximize the accepted family")
    if not np.isclose(
        fit["contrast_log_likelihood"], best[winner], rtol=1e-9, atol=1e-8
    ):
        raise ValueError("Saved likelihood disagrees with direct refit")
    alpha = search.grid[search.cache[at[winner, 0]][0]]
    expected_height = None if np.isinf(alpha) else float(alpha)
    expected_alpha = None if np.isinf(alpha) else float(alpha / search.height)
    expected_status = (
        "brownian_limit"
        if alpha == 0
        else ("independent_limit" if np.isinf(alpha) else "finite")
    )
    if (
        fit["alpha_height"] != expected_height
        or fit["alpha"] != expected_alpha
        or fit["alpha_status"] != expected_status
    ):
        raise ValueError("Alpha estimate/status disagrees with refitted grid maximum")
    if replay_bootstrap:
        replay = search.fit(values, seed=fit["seed"], replicates=B, level=level)
        if [test["p_value"] for test in replay["tests"]] != [
            test["p_value"] for test in tests
        ]:
            raise ValueError(
                "Saved bootstrap probabilities disagree with seeded replay"
            )
    profile = fit["alpha_profile"]
    expected_grid = [None if np.isinf(a) else float(a) for a in search.grid]
    if [p["alpha_height"] for p in profile] != expected_grid:
        raise ValueError("Alpha grid disagrees with fit settings")
    ll = np.array([p["log_likelihood"] for p in profile])
    refitted = search.alpha_profile(search.q @ (values - values.mean()), winner)
    if not np.allclose(ll, refitted, rtol=1e-9, atol=1e-8):
        raise ValueError("Alpha profile disagrees with direct refit")
    supported = ll >= ll.max() - 1.920729410347062
    limit_supported = bool(
        np.any(supported & ((search.grid == 0) | np.isinf(search.grid)))
    )
    if (
        not np.all(np.isfinite(ll))
        or not np.isclose(ll.max(), best[winner], rtol=1e-9, atol=1e-8)
        or [p["supported"] for p in profile] != supported.tolist()
        or fit["alpha_limit_supported"] != limit_supported
    ):
        raise ValueError("Inconsistent alpha support diagnostics")
    groups = [
        list(g)
        for g in tip_groups(
            search.tree, fit["model"]["shift_branch_ids"], fit["model"]["groups"]
        )
    ]
    metrics = {
        "any_shift": bool(fit["model"]["shift_branch_ids"]),
        "partition_recovered": sorted(groups) == truth["shared_partition"],
        "mean_rmse": float(
            np.sqrt(
                np.mean(
                    (
                        np.asarray(fit["predicted"])
                        - np.asarray(list(truth["tip_mean"].values()))
                    )
                    ** 2
                )
            )
        ),
    }
    for key, value in metrics.items():
        if not np.isclose(row[key], value, rtol=1e-12, atol=1e-12):
            raise ValueError(f"Saved {key} disagrees with reconstructed fit metrics")
    return {**metrics, "alpha_limit_supported": limit_supported}


def audit_records(rows, design, summary, *, replay_bootstrap=False):
    engines, grouped = {}, defaultdict(list)
    for row in rows:
        c = row["case"]
        key = row["tree"], c["standard_error"]
        if key not in engines:
            engines[key] = CalibratedSearch(
                read_tree(row["tree"], "auto", True, quiet=True),
                variances=np.full(c["tips"], c["standard_error"] ** 2),
            )
        result = audit_record(
            row, design, engines[key], replay_bootstrap=replay_bootstrap
        )
        grouped[c["family"], c["scenario"], c["root_model"]].append(result)
    recomputed = []
    for (family, scenario, root), members in grouped.items():
        recomputed.append(
            {
                "family": family,
                "scenario": scenario,
                "root_model": root,
                "attempted": len(members),
                "completed": len(members),
                **{
                    key: rate(sum(r[key] for r in members), len(members))
                    for key in (
                        "any_shift",
                        "partition_recovered",
                        "alpha_limit_supported",
                    )
                },
            }
        )
    if summary != recomputed:
        raise ValueError("Summary disagrees with reconstructed record metrics")
    return recomputed
