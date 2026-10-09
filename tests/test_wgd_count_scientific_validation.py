"""Independent count-model validation, including a bounded replicated study.

Use the existing checker, for example:
  python tools/check.py test -- tests/test_wgd_count_scientific_validation.py --run-studies -s

The simulator below draws individual birth/death waiting times and Bernoulli
retentions/detections; it does not use production transitions or simulation.
All families originate at the root. Missing masks and detection are known and
independent of latent counts. Family rate categories are an explicit discrete
mixture, not continuous gamma. Local SSD deliberately misspecifies the fitted
homogeneous null. All four branches are searched on a fixed midpoint grid.

The opt-in study defaults to eight datasets each for null/local SSD and four for
the other regimes, with 19 refitted null draws per dataset (minimum p=.05).
NWKIT_WGD_COUNT_STUDY_REPLICATES changes the null/SSD replication count;
NWKIT_WGD_COUNT_STUDY_BOOTSTRAP changes draws for all regimes. Settings, seeds,
all bootstrap statistics, failures, and exact 95% binomial intervals are saved
in pytest's temporary directory and printed. No nominal power/FPR target is
asserted: intervals describe this simulation experiment, not empirical genomes.
Failed fits are recorded and fail validation, never silently discarded.
"""

import json
import os
import platform
import sys
import time
from dataclasses import dataclass

import numpy as np
import pytest
import scipy
from scipy.stats import beta, binomtest

from nwkit.wgd_count_fit import _burst_model, calibrate_scan, scan_counts
from nwkit.wgd_count_model import CountLikelihood, CountTree, MultiplicationEvent


def _tree():
    return CountTree(
        (-1, 0, 1, 1, 0),
        (0.0, 0.6, 0.4, 0.4, 1.0),
        (2, 3, 4),
        ("A", "B", "C"),
        (0, 1, 2, 3, 4),
        ("ABC", "AB", "A", "B", "C"),
    )


def _gillespie_edge(count, duplication, loss, duration, rng):
    elapsed = 0.0
    for _ in range(100000):
        intensity = count * (duplication + loss)
        if intensity == 0:
            return count
        elapsed += rng.exponential(1 / intensity)
        if elapsed >= duration:
            return count
        count += 1 if rng.random() < duplication / (duplication + loss) else -1
        if count > 10000:
            raise ValueError(
                "Independent Gillespie simulation exceeded copy safety bound."
            )
    raise ValueError("Independent Gillespie simulation exceeded jump safety bound.")


def _independent_counts(
    families,
    rates,
    root_mean,
    seed,
    *,
    event=None,
    groups=(0, 0, 0, 0, 0),
    scales=(1.0,),
    detection=(1.0, 1.0, 1.0),
    missing=False,
    ascertainment="observed",
):
    design = _tree()
    rng = np.random.default_rng(seed)
    result = np.full((families, 3), np.nan)
    masks = np.zeros_like(result, dtype=bool)
    if missing:
        masks[::4, 0] = True
        masks[1::4, 1] = True
    for row in range(families):
        for _ in range(10000):
            root = 1
            while rng.random() > 1 / root_mean:
                root += 1
            latent = [root]
            scale = scales[int(rng.integers(len(scales)))]
            for node in range(1, len(design.parents)):
                count = latent[design.parents[node]]
                duplication, loss = np.asarray(rates[groups[node]]) * scale
                if event is not None and node == event.node:
                    count = _gillespie_edge(
                        count,
                        duplication,
                        loss,
                        design.lengths[node] * event.fraction,
                        rng,
                    )
                    extra = np.count_nonzero(
                        rng.random(count * (event.multiplicity - 1)) < event.retention
                    )
                    count += int(extra)
                    duration = design.lengths[node] * (1 - event.fraction)
                else:
                    duration = design.lengths[node]
                latent.append(_gillespie_edge(count, duplication, loss, duration, rng))
            observed = np.array(
                [
                    np.nan
                    if masks[row, column]
                    else np.count_nonzero(rng.random(latent[node]) < detection[column])
                    for column, node in enumerate(design.tip_nodes)
                ]
            )
            selected = np.nansum(observed) > 0
            if ascertainment == "root-clades":
                selected = np.nansum(observed[:2]) > 0 and observed[2] > 0
            if selected:
                result[row] = observed
                break
        else:
            raise ValueError(
                "Independent simulation cannot sample an ascertained family."
            )
    np.testing.assert_array_equal(np.isnan(result), masks)
    return result


@pytest.mark.parametrize("ascertainment", ["observed", "root-clades"])
def test_independent_history_pattern_probabilities_match_pruning(ascertainment):
    rates = np.array([[0.2, 0.45]])
    event = MultiplicationEvent(1, 0.65)
    scales = (0.15, 0.55, 1.05, 2.25)
    detection = (0.7, 0.8, 0.9)
    counts = _independent_counts(
        16000,
        rates,
        1.4,
        68113,
        event=event,
        scales=scales,
        detection=detection,
        missing=True,
        ascertainment=ascertainment,
    )
    # Compare complete patterns within each fixed missingness stratum, not
    # frequencies pooled across different ascertainment denominators.
    for residue in range(4):
        rows = counts[residue::4]
        present = ~np.isnan(rows[0])
        for target in (np.ones(3), np.array([2.0, 2.0, 1.0])):
            target = target.copy()
            target[~present] = np.nan
            model = CountLikelihood(
                _tree(),
                target[None, :],
                detection=np.array(detection),
                rate_scales=np.array(scales),
                ascertainment=ascertainment,
            )
            probability = np.exp(model.log_likelihood(rates, 1.4, 64, event))
            matches = np.count_nonzero(
                np.all(rows[:, present] == target[present], axis=1)
            )
            assert binomtest(matches, len(rows), probability).pvalue > 1e-5
            assert model.convergence(rates, 1.4, 64, event) < 1e-7


@dataclass(frozen=True)
class _Scenario:
    name: str
    rates: tuple[tuple[float, float], ...] = ((0.12, 0.22),)
    root_mean: float = 1.2
    event: MultiplicationEvent | None = None
    groups: tuple[int, ...] = (0, 0, 0, 0, 0)
    scales: tuple[float, ...] = (1.0,)
    detection: tuple[float, ...] = (1.0, 1.0, 1.0)
    missing: bool = False
    misspecified_groups: bool = False


_SCENARIOS = (
    _Scenario("null"),
    _Scenario("pure-wgd", ((0.0, 0.0),), 1.0, MultiplicationEvent(1, 1.0)),
    _Scenario("strong-wgd", event=MultiplicationEvent(1, 0.85)),
    _Scenario(
        "local-ssd",
        ((0.08, 0.15), (0.9, 0.15)),
        groups=(0, 1, 0, 0, 0),
        misspecified_groups=True,
    ),
    _Scenario(
        "heterogeneous",
        ((0.08, 0.15), (0.25, 0.35)),
        groups=(0, 1, 0, 0, 0),
        scales=(0.15, 0.55, 1.05, 2.25),
    ),
    _Scenario("missing-annotation", detection=(0.65, 0.8, 0.9), missing=True),
)


def _scenario_model(scenario, seed, families=64):
    counts = _independent_counts(
        families,
        scenario.rates,
        scenario.root_mean,
        seed,
        event=scenario.event,
        groups=scenario.groups,
        scales=scenario.scales,
        detection=scenario.detection,
        missing=scenario.missing,
    )
    return CountLikelihood(
        _tree(),
        counts,
        detection=np.array(scenario.detection),
        rate_scales=np.array(scenario.scales),
        branch_groups=None if scenario.misspecified_groups else scenario.groups,
    )


def _verified_scan(model):
    result = scan_counts(model, fractions=(0.5,), max_states=128)
    for fit in (
        result.background,
        *(c.event_fit for c in result.candidates),
        *(c.burst_fit for c in result.candidates),
    ):
        assert fit.converged and fit.state_error <= 1e-7
        assert np.isfinite(fit.log_likelihood)
    assert {c.node for c in result.candidates} == {1, 2, 3, 4}
    return result


@pytest.mark.parametrize(
    "index,scenario", list(enumerate(_SCENARIOS)), ids=[s.name for s in _SCENARIOS]
)
def test_independent_scenarios_do_not_lose_known_feasible_likelihood(index, scenario):
    model = _scenario_model(scenario, 73131 + index)
    scan = _verified_scan(model)
    if scenario.misspecified_groups:
        candidate = next(c for c in scan.candidates if c.node == 1)
        fit = candidate.burst_fit
        reference_model = _burst_model(model, 1)
    elif scenario.event is not None:
        fit = next(
            c.event_fit for c in scan.candidates if c.node == scenario.event.node
        )
        reference_model = model
    else:
        fit, reference_model = scan.background, model
    # Exact zero rates are outside the optimizer's documented numerical box.
    # The pure-WGD identity is checked below, but its feasible fit comparator
    # uses the same lower bound rather than claiming zero is an allowed MLE.
    time_unit = np.median([t for t in model.tree.lengths[1:] if t > 0])
    feasible_rates = np.maximum(scenario.rates, np.exp(-16.0) / time_unit)
    feasible = reference_model.log_likelihood(
        feasible_rates, scenario.root_mean, fit.max_count, scenario.event
    )
    if scenario.name == "pure-wgd":
        np.testing.assert_array_equal(model.counts, np.tile([2, 2, 1], (64, 1)))
        assert model.log_likelihood([[0.0, 0.0]], 1.0, 8, scenario.event) == 0.0
    assert fit.log_likelihood >= feasible - 1e-6
    print(
        json.dumps(
            {
                "pilot": scenario.name,
                "seed": 73131 + index,
                "families": len(model.counts),
                "top_node": scan.candidates[0].node,
                "statistic": scan.candidates[0].improvement,
                "burst_minus_event_aic": scan.candidates[0].burst_aic_difference,
            }
        )
    )


def _binomial_interval(successes, total):
    return [
        0.0
        if successes == 0
        else float(beta.ppf(0.025, successes, total - successes + 1)),
        1.0
        if successes == total
        else float(beta.ppf(0.975, successes + 1, total - successes)),
    ]


@pytest.mark.slow
@pytest.mark.study
@pytest.mark.parametrize(
    "index,scenario", list(enumerate(_SCENARIOS)), ids=[s.name for s in _SCENARIOS]
)
def test_replicated_independent_bootstrap_study(index, scenario, tmp_path):
    replicates = int(os.environ.get("NWKIT_WGD_COUNT_STUDY_REPLICATES", "8"))
    if scenario.name not in {"null", "local-ssd"}:
        replicates = min(replicates, 4)
    draws = int(os.environ.get("NWKIT_WGD_COUNT_STUDY_BOOTSTRAP", "19"))
    assert replicates >= 1 and draws >= 19
    records, failures = [], []
    start = time.perf_counter()
    path = tmp_path / f"{scenario.name}.json"
    protocol = {
        "regime": scenario.name,
        "replicates": replicates,
        "bootstrap_draws": draws,
        "families": 64,
        "alpha": 0.05,
        "fractions": [0.5],
        "candidate_nodes": [1, 2, 3, 4],
        "max_states": 128,
        "state_tolerance": 1e-7,
        "max_iterations": 200,
        "root_mean": scenario.root_mean,
        "rates": scenario.rates,
        "groups": scenario.groups,
        "scales": scenario.scales,
        "detection": scenario.detection,
        "missing": scenario.missing,
        "fitted_null_groups_misspecified": scenario.misspecified_groups,
        "retention": None if scenario.event is None else scenario.event.retention,
        "python": sys.version,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
    }
    for replicate in range(replicates):
        data_seed = 170000 + 1000 * index + replicate
        bootstrap_seed = 910000 + 1000 * index + replicate
        model = _scenario_model(scenario, data_seed)
        serialized_counts = [
            [None if np.isnan(value) else int(value) for value in row]
            for row in model.counts
        ]
        try:
            observed = _verified_scan(model)
            calibrated = calibrate_scan(
                model, observed, draws, bootstrap_seed, max_states=128
            )
            candidates = [
                {
                    "node": c.node,
                    "statistic": c.improvement,
                    "p_value": c.p_value,
                    "p_value_mc_se": c.p_value_mc_se,
                    "burst_minus_event_aic": c.burst_aic_difference,
                    "retention": c.event_fit.event.retention,
                    "state_error": c.event_fit.state_error,
                    "nuisance_bound_reached": c.event_fit.boundary,
                }
                for c in calibrated.candidates
            ]
            record = {
                "replicate": replicate,
                "data_seed": data_seed,
                "bootstrap_seed": bootstrap_seed,
                "counts": serialized_counts,
                "background": {
                    "rates": observed.background.rates.tolist(),
                    "root_mean": observed.background.root_mean,
                    "log_likelihood": observed.background.log_likelihood,
                    "state_error": observed.background.state_error,
                    "nuisance_bound_reached": observed.background.boundary,
                },
                "candidates": candidates,
                "bootstrap_max_statistics": calibrated.bootstrap_statistics,
                "raw_rejection": any(c.p_value <= 0.05 for c in calibrated.candidates),
                "conditional_support": any(
                    c.p_value <= 0.05 and c.burst_aic_difference > 0
                    for c in calibrated.candidates
                ),
                "true_branch_top_rank": calibrated.candidates[0].node == 1,
                "true_branch_supported": any(
                    c.node == 1 and c.p_value <= 0.05 and c.burst_aic_difference > 0
                    for c in calibrated.candidates
                ),
            }
            records.append(record)
            print(
                json.dumps(
                    {
                        "regime": scenario.name,
                        "replicate": replicate,
                        "data_seed": data_seed,
                        "raw_rejection": record["raw_rejection"],
                        "conditional_support": record["conditional_support"],
                        "top_p": calibrated.candidates[0].p_value,
                        "elapsed_seconds": time.perf_counter() - start,
                    }
                ),
                flush=True,
            )
        except (ValueError, AssertionError) as exc:
            failures.append(
                {
                    "replicate": replicate,
                    "data_seed": data_seed,
                    "bootstrap_seed": bootstrap_seed,
                    "counts": serialized_counts,
                    "error": str(exc),
                }
            )
        path.write_text(
            json.dumps(
                {"protocol": protocol, "records": records, "failures": failures},
                indent=2,
                allow_nan=False,
            ),
            encoding="utf-8",
        )
    summary = {
        "regime": scenario.name,
        "completed": len(records),
        "failed": len(failures),
        "output": str(path),
        "elapsed_seconds": time.perf_counter() - start,
    }
    for field in (
        "raw_rejection",
        "conditional_support",
        "true_branch_top_rank",
        "true_branch_supported",
    ):
        successes = sum(record[field] for record in records)
        summary[field] = {
            "count": successes,
            "denominator": len(records),
            "exact_95_percent_interval": _binomial_interval(successes, len(records))
            if records
            else None,
        }
    path.write_text(
        json.dumps(
            {
                "protocol": protocol,
                "records": records,
                "failures": failures,
                "summary": summary,
            },
            indent=2,
            allow_nan=False,
        ),
        encoding="utf-8",
    )
    print(json.dumps(summary), flush=True)
    assert not failures, (
        f"{scenario.name}: failed refits retained in {path}: {failures}"
    )
