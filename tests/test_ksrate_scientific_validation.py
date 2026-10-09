"""Additive-tree references and repeated-dataset bootstrap coverage."""

import itertools
import json

import numpy as np
import pytest
from ete4 import Tree
from scipy.stats import binom, binomtest

from nwkit.ksrate import make_trios
from nwkit.ksrate_model import (
    KsObservation,
    PairKs,
    bootstrap_corrections,
    grouped_corrections,
    median_order_statistic_interval,
    simultaneous_median_corrections,
    trio_values,
)


def balanced_newick(names):
    if len(names) == 1:
        return names[0] + ":1"
    midpoint = len(names) // 2
    return (
        "("
        + balanced_newick(names[:midpoint])
        + ","
        + balanced_newick(names[midpoint:])
        + "):1"
    )


@pytest.mark.parametrize("seed", range(20))
def test_all_trios_recover_focal_path_on_asymmetric_twenty_species_trees(seed):
    tree = Tree(balanced_newick([f"S{i:02d}" for i in range(20)]) + ";", parser=1)
    rng = np.random.default_rng(seed)
    for node in tree.traverse():
        node.dist = 0.0 if node.is_root else float(rng.uniform(0.01, 2))
    leaves = tuple(tree.leaves())
    records = [
        KsObservation(a.name, b.name, "reference", tree.get_distance(a, b))
        for a, b in itertools.combinations(leaves, 2)
    ]
    trios, events = make_trios(tree, outgroup_policy="all")
    result = trio_values(
        trios, PairKs(records, [leaf.name for leaf in leaves]).estimates()
    )
    for trio, correction in zip(trios, result, strict=True):
        ancestor = tree.common_ancestor([trio.focal, trio.sister])
        expected = 2 * tree.get_distance(ancestor, trio.focal)
        assert correction == pytest.approx(expected, rel=1e-12, abs=1e-12)
    assert all(
        event["trio_indices"] == []
        for event in events
        if event["descendant_taxa"].count(",") == 19
    )


def _bootstrap_reference(scales, trios, draws, seed):
    key = (trios[0].focal, trios[0].species_event_id)
    observations = [
        KsObservation(a, b, f"f{i:03d}", distance * scale)
        for a, b, distance in (("A", "B", 0.4), ("A", "C", 1.2), ("B", "C", 1.4))
        for i, scale in enumerate(scales)
    ]
    pairs = PairKs(observations, "ABC")
    result = bootstrap_corrections(pairs, trios, draws, seed=seed)[key]
    # Direct iid family-index resampling is an independent bootstrap reference.
    direct_rng = np.random.default_rng(seed)
    selected = direct_rng.integers(0, len(scales), size=(draws, len(scales)))
    direct = 0.2 * np.median(scales[selected], axis=1)
    np.testing.assert_allclose(result, direct, rtol=2e-14, atol=2e-14)
    point = grouped_corrections(trios, trio_values(trios, pairs.estimates()))[key]
    assert point == pytest.approx(0.2 * np.median(scales), rel=2e-14)
    return pairs, result


@pytest.mark.parametrize("seed", [616, 617, 618])
def test_family_bootstrap_matches_independent_resampling(seed):
    tree = Tree("((A:0.1,B:0.3):0.4,C:0.7);", parser=1)
    trios, _ = make_trios(tree, ["A"])
    scales = np.random.default_rng(seed).lognormal(mean=0, sigma=0.6, size=80)
    _bootstrap_reference(scales, trios, 199, seed=10000 + seed)


@pytest.mark.slow
@pytest.mark.study
def test_repeated_dataset_coverage_and_independent_bootstrap_reference():
    tree = Tree("((A:0.1,B:0.3):0.4,C:0.7);", parser=1)
    trios, _ = make_trios(tree, ["A"])
    key = (trios[0].focal, trios[0].species_event_id)
    true_focal_ks = 0.2
    rng = np.random.default_rng(616)
    experiments, families, draws = 5000, 80, 199
    covered = conservative_covered = 0
    widths, conservative_widths = [], []
    for experiment in range(experiments):
        scales = rng.lognormal(mean=0, sigma=0.6, size=families)
        pairs, result = _bootstrap_reference(scales, trios, draws, 10000 + experiment)
        low, high = np.quantile(result, [0.025, 0.975])
        covered += int(low <= true_focal_ks <= high)
        widths.append(float(high - low))
        conservative, _ = simultaneous_median_corrections(pairs, trios, 0.95)
        low, high = conservative[key]
        conservative_covered += int(low <= true_focal_ks <= high)
        conservative_widths.append(high - low)
    # Legacy percentile coverage is measured, not assumed adequate. The
    # alternative's finite-sample claim is separately checked by exact ranks.
    interval = binomtest(covered, experiments).proportion_ci(confidence_level=0.95)
    conservative_interval = binomtest(conservative_covered, experiments).proportion_ci(
        confidence_level=0.95
    )
    assert (
        binomtest(conservative_covered, experiments, 0.95, alternative="less").pvalue
        > 1e-6
    )
    print(
        json.dumps(
            {
                "study": "ks_family_bootstrap_coverage",
                "seed": 616,
                "datasets": experiments,
                "families_per_dataset": families,
                "bootstrap_draws": draws,
                "nominal_coverage": 0.95,
                "covered": covered,
                "observed_coverage": covered / experiments,
                "coverage_ci_lower": interval.low,
                "coverage_ci_upper": interval.high,
                "median_interval_width": float(np.median(widths)),
                "conservative_covered": conservative_covered,
                "conservative_observed_coverage": conservative_covered / experiments,
                "conservative_coverage_ci_lower": conservative_interval.low,
                "conservative_coverage_ci_upper": conservative_interval.high,
                "conservative_median_interval_width": float(
                    np.median(conservative_widths)
                ),
                "assumptions": "iid lognormal shared family multipliers; additive true distances; all pairs observed; population family median multiplier one",
            },
            sort_keys=True,
        )
    )


@pytest.mark.parametrize("n", [1, 5, 6, 10, 30, 80, 199])
@pytest.mark.parametrize("error", [0.05, 0.05 / 3, 0.01])
def test_finite_sample_median_ranks_meet_exact_binomial_coverage(n, error):
    lower, upper = median_order_statistic_interval(np.arange(n), error)
    if np.isnan(lower):
        assert 2 * 0.5**n > error
        return
    k = int(lower) + 1
    assert int(upper) == n - k
    assert 2 * binom.cdf(k - 1, n, 0.5) <= error
    assert 2 * binom.cdf(k, n, 0.5) > error
