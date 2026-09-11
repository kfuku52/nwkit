"""Full-search plug-in bootstrap; no uniform composite-null guarantee."""

import itertools
import math

import numpy as np

from nwkit.shift_native_model import ShiftData, covariance_geometry


def simulate_native_data(data, result, rng):
    """Simulate in original units, preserving errors and observation masks."""
    if "joint_fit" in result:
        from nwkit.shift_simulation import simulate_joint_data

        return simulate_joint_data(data, result, rng)
    values = np.empty_like(data.values)
    tree = data.tree
    for trait, fit in enumerate(result["fits"]):
        slopes, innovations, root_variance = covariance_geometry(
            tree, fit.alpha_height, fit.process_variance, fit.root_model
        )
        state = np.zeros(len(tree.branch_ids))
        state[0] = rng.normal() * math.sqrt(root_variance)
        noise = rng.normal(size=len(state))
        parents = np.asarray(tree.compiled.parents)
        for indices in tree.levels:
            state[indices] = (
                slopes[indices] * state[parents[indices]]
                + np.sqrt(innovations[indices]) * noise[indices]
            )
        observed = (
            fit.predicted
            + state[list(tree.compiled.leaf_indices)]
            + rng.normal(size=len(tree.leaf_names))
            * np.sqrt(data.variances[:, trait] + fit.measurement_variance)
        )
        values[:, trait] = data.centers[trait] + data.scales[trait] * observed
    values[~np.isfinite(data.values)] = np.nan
    errors = data.variances * data.scales[None, :] ** 2
    return ShiftData.build(tree, values, data.trait_names, errors)


def _family_statistic(search, complexity):
    restricted = [result for cost, result in search.families() if cost <= complexity]
    if not restricted:
        raise ValueError("Bootstrap search did not retain its null family.")
    return max(
        0.0,
        2
        * (
            search.families()[-1][1]["log_likelihood"]
            - restricted[-1]["log_likelihood"]
        ),
    )


def calibrate_native_search(data, search, run_search, *, replicates, seed, level=0.05):
    """Repeat the complete search for each draw; stop at first accepted family.

    run_search must capture exactly the same configuration used for search.
    Numerical failures propagate, including their family/replicate coordinates.
    """
    if (
        isinstance(replicates, bool)
        or not isinstance(replicates, int)
        or replicates < 1
    ):
        raise ValueError("Native calibration requires a positive integer draw count.")
    if not 0 < level < 1 or 1 / (replicates + 1) > level:
        raise ValueError("Calibration draws cannot resolve the requested test level.")
    families = search.families()
    streams = np.random.SeedSequence(seed).spawn(max(0, len(families) - 1))
    tests = []
    selected = families[-1][1]
    for family_index, (complexity, restricted) in enumerate(families[:-1]):
        statistic = _family_statistic(search, complexity)
        samples = []
        for draw, stream in enumerate(streams[family_index].spawn(replicates)):
            try:
                simulated = simulate_native_data(
                    data, restricted, np.random.default_rng(stream)
                )
                bootstrap = run_search(simulated)
                samples.append(_family_statistic(bootstrap, complexity))
            except (ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
                raise ValueError(
                    f"Native calibration failed at complexity {complexity}, draw {draw}; no draws discarded: {exc}"
                ) from exc
        exceedances = sum(
            value >= statistic - 1e-10 * max(1, statistic) for value in samples
        )
        p_value = (1 + exceedances) / (replicates + 1)
        tests.append(
            {
                "complexity": complexity,
                "statistic": statistic,
                "exceedances": exceedances,
                "replicates": replicates,
                "p_value": p_value,
                "rejected": p_value <= level,
                "bootstrap_statistics": samples,
            }
        )
        if p_value > level:
            selected = restricted
            break
    return selected, {
        "method": "full_search_plugin_parametric_bootstrap_fixed_sequence",
        "research_only": True,
        "uniform_composite_null_control_proven": False,
        "full_search_repeated": True,
        "failed_draw_policy": "abort_without_discarding",
        "seed": seed,
        "replicates_per_test": replicates,
        "level": level,
        "tests": tests,
    }


def _aic_gain(search):
    null = search.best_by_complexity.get(0)
    selected = search.best_information
    if null is None or selected is None:
        raise ValueError(
            "Global-null calibration requires the null fit and AIC winner."
        )
    scores = []
    for result in (null, selected):
        record = result.get("information_criterion", {})
        score = record.get("score")
        if (
            record.get("criterion") != "AIC"
            or score is None
            or not math.isfinite(score)
        ):
            raise ValueError("Global-null calibration requires finite AIC scores.")
        scores.append(score)
    return max(0.0, scores[0] - scores[1])


def gate_native_aic(data, search, run_search, *, replicates, seed, level=0.05):
    """Test the global no-shift null by replaying the entire AIC search.

    The statistic is max(0, AIC(null) - min AIC). Passing the gate retains
    the original AIC winner; this does not test individual selected branches.
    """
    if (
        isinstance(replicates, bool)
        or not isinstance(replicates, int)
        or replicates < 1
    ):
        raise ValueError("Native calibration requires a positive integer draw count.")
    if not 0 < level < 1 or 1 / (replicates + 1) > level:
        raise ValueError("Calibration draws cannot resolve the requested test level.")
    statistic = _aic_gain(search)
    null = search.best_by_complexity[0]
    samples = []
    for draw, stream in enumerate(np.random.SeedSequence(seed).spawn(replicates)):
        try:
            simulated = simulate_native_data(data, null, np.random.default_rng(stream))
            samples.append(_aic_gain(run_search(simulated)))
        except (ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
            raise ValueError(
                f"Native global-null calibration failed at draw {draw}; no draws discarded: {exc}"
            ) from exc
    exceedances = sum(
        value >= statistic - 1e-10 * max(1, statistic) for value in samples
    )
    p_value = (1 + exceedances) / (replicates + 1)
    rejected = p_value <= level
    return (search.best_information if rejected else null), {
        "method": "full_search_plugin_parametric_bootstrap_global_null_aic_gate",
        "research_only": True,
        "uniform_composite_null_control_proven": False,
        "controls_false_branches_under_nonnull": False,
        "full_search_repeated": True,
        "failed_draw_policy": "abort_without_discarding",
        "statistic_name": "max_zero_null_aic_minus_minimum_search_aic",
        "seed": seed,
        "replicates": replicates,
        "level": level,
        "statistic": statistic,
        "exceedances": exceedances,
        "p_value": p_value,
        "rejected": rejected,
        "bootstrap_statistics": samples,
        "ungated_shift_branch_ids": list(search.best_information["layout"].shifts),
    }


def native_selection_support(data, selected, select_data, *, replicates, seed):
    """Stability frequencies; select_data replays the entire selection procedure."""
    if (
        isinstance(replicates, bool)
        or not isinstance(replicates, int)
        or replicates < 1
    ):
        raise ValueError("Native support requires a positive integer draw count.")
    branches = {branch: 0 for branch in data.tree.branch_ids if branch}
    grouping = 0
    pairs = {
        pair: [0, 0]
        for pair in itertools.combinations((0, *selected["layout"].shifts), 2)
    }
    streams = np.random.SeedSequence(seed).spawn(replicates)
    for draw, stream in enumerate(streams):
        simulation_stream, selection_stream = stream.spawn(2)
        try:
            simulated = simulate_native_data(
                data, selected, np.random.default_rng(simulation_stream)
            )
            result = select_data(simulated, int(selection_stream.generate_state(1)[0]))
        except (ValueError, ArithmeticError, np.linalg.LinAlgError) as exc:
            raise ValueError(
                f"Native support failed at draw {draw}; no draws discarded: {exc}"
            ) from exc
        for branch in result["layout"].shifts:
            branches[branch] += 1
        grouping += result["layout"] == selected["layout"]
        labels = {
            branch: group
            for group, members in enumerate(result["layout"].groups)
            for branch in members
        }
        for (first, second), counts in pairs.items():
            if first in labels and second in labels:
                counts[0] += 1
                counts[1] += labels[first] == labels[second]
    return {
        "method": "parametric_full_selection_stability",
        "interpretation": "selection_frequency_under_fitted_model_not_posterior_or_p_value",
        "seed": seed,
        "replicates": replicates,
        "branch_frequencies": [
            {"branch_id": b, "count": count, "frequency": count / replicates}
            for b, count in branches.items()
        ],
        "exact_layout_frequency": grouping / replicates,
        "regime_pair_frequencies": [
            {
                "first_branch_id": first,
                "second_branch_id": second,
                "cooccurrences": cooccurrences,
                "shared_regime_count": shared,
                "shared_regime_frequency": shared / replicates,
                "shared_given_both_selected": shared / cooccurrences
                if cooccurrences
                else None,
            }
            for (first, second), (cooccurrences, shared) in pairs.items()
        ],
    }
