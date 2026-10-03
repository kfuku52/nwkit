"""Independent full-forest/mask checks for the research integration prototype."""

import itertools
import math
from collections import Counter

import numpy as np
import pytest

from nwkit.mul_locus import Locus, LocusParameters, sample_locus_tree
from nwkit.mul_locus_integral import (
    OVERFLOW,
    IntegratedBank,
    Moments,
    WorkBudget,
    accumulate,
    bounded_interval,
    build_paired_banks,
    chernoff_interval,
    conditional_observations,
    detection_distribution,
    genealogy_distribution,
    integrated_probability,
    score_integrated_bank,
)
from nwkit.mul_locus_mc import LocusBank, calibrate, make_tasks, score_bank
from nwkit.mul_msc_fit import species_topology_signature
from tests.test_mul_coalescent import oracle_ancestral
from tests.test_mul_locus import (
    SPECIES,
    color,
    joint_locus_forest,
    model,
    nested_locus_with_loss,
    strip_ids,
)
from tests.test_mul_reconcile import parser, tree


def joint_oracle(root, ne):
    pending, tips = [root], []
    while pending:
        node = pending.pop()
        if node.kind == "tip":
            tips.append(node)
        pending.extend(node.children)
    labels = {tip: f"g{i}_{tip.species}" for i, tip in enumerate(tips)}
    forest = joint_locus_forest(root, None, ne, labels)
    normalizer = math.fsum(forest.values())
    result = Counter()
    for state, p in forest.items():
        if not state:
            result[None] += p / normalizer
        for topology, q in oracle_ancestral(state).items() if state else ():
            result[topology] += p * q / normalizer
    return result


def mask_oracle(topologies, detection):
    result = Counter()
    for text, mass in topologies.items():
        if text is None:
            result[None] += mass
            continue
        gene = tree(text + ";")
        names = list(gene.leaf_names())
        for mask in itertools.product((False, True), repeat=len(names)):
            kept = [name for name, present in zip(names, mask, strict=True) if present]
            p = math.prod(
                detection[name.rsplit("_", 1)[1]]
                if present
                else 1 - detection[name.rsplit("_", 1)[1]]
                for name, present in zip(names, mask, strict=True)
            )
            if not p:
                continue
            if not kept:
                signature = None
            elif len(kept) == 1:
                signature = ("tip", kept[0].rsplit("_", 1)[1])
            else:
                pruned = gene.copy()
                pruned.prune(kept)
                signature = species_topology_signature(pruned, parser())
            result[signature] += mass * p
    return result


@pytest.mark.parametrize("candidate", [0, 1, 2])
@pytest.mark.parametrize("ne", [0.1, 0.7, 3.0])
def test_all_genealogies_and_detection_against_independent_forest(candidate, ne):
    point = LocusParameters(0, 0, ne, 0.5)
    population = next(
        t[3]
        for t in make_tasks(tree(SPECIES), "X", "A B", [point], 100)[0]
        if t[0] == candidate
    )
    locus = sample_locus_tree(population, point, 0.5, np.random.default_rng(1))
    expected = joint_oracle(locus, ne)
    colored = Counter()
    for text, p in expected.items():
        colored[color(text)] += p
    actual = Counter()
    for gene, p in genealogy_distribution(locus, ne, max_states=1000000).items():
        actual[strip_ids(gene)] += p
    assert dict(actual) == pytest.approx(dict(colored), abs=3e-13, rel=2e-12)
    detection = {"A": 0.8, "X": 0.6, "B": 1.0}
    assert conditional_observations(
        locus, ne, detection, max_states=1000000
    ) == pytest.approx(dict(mask_oracle(expected, detection)), abs=3e-13, rel=2e-12)


def test_nested_daughter_bounds_loss_and_undetected_copies_are_retained():
    locus, _ = nested_locus_with_loss()
    detection = {"A": 0.85, "B": 0.8, "X": 0.65, "C": 0.0}
    expected = mask_oracle(joint_oracle(locus, 0.7), detection)
    actual = conditional_observations(locus, 0.7, detection, max_states=1000000)
    assert actual == pytest.approx(dict(expected), abs=3e-13, rel=2e-12)
    assert math.fsum(actual.values()) == pytest.approx(1, abs=3e-13)
    assert any(key is not None and str(key).count("'X'") == 2 for key in actual)


@pytest.mark.parametrize("duration", [1e-100, 0.001, 2.0])
def test_tiny_daughter_normalizer_is_not_rounded_to_zero(duration):
    daughter = Locus(
        0,
        "speciation",
        children=[Locus(0, "tip", "A"), Locus(0, "tip", "B")],
        daughter=True,
    )
    locus = Locus(duration, "origin", children=[daughter])
    assert conditional_observations(locus, 1, {"A": 1, "B": 1}) == pytest.approx(
        {("node", ("tip", "A"), ("tip", "B")): 1.0}
    )


def test_extinct_and_zero_duration_histories():
    assert conditional_observations(
        Locus(1, "origin", children=[Locus(0, "loss")]), 1, {}
    ) == {None: 1.0}
    locus = Locus(
        0,
        "origin",
        children=[
            Locus(
                0, "speciation", children=[Locus(0, "tip", "A"), Locus(0, "tip", "B")]
            )
        ],
    )
    assert len(genealogy_distribution(locus, 1)) == 1


def test_detection_distribution_keeps_all_sizes_and_coalesces_identical_colors():
    gene = ("node", ("tip", "X", 0), ("node", ("tip", "X", 1), ("tip", "X", 2)))
    actual = detection_distribution(gene, {"X": 0.5})
    assert actual[None] == 0.125
    assert actual[("tip", "X")] == 0.375
    assert actual[("node", ("tip", "X"), ("tip", "X"))] == 0.375
    assert actual[strip_ids(gene)] == 0.125


def test_detection_overflow_is_exact_probability_aggregation_not_hidden_tip_truncation():
    gene = (
        "node",
        ("tip", "A", 0),
        ("node", ("tip", "X", 1), ("node", ("tip", "B", 2), ("tip", "X", 3))),
    )
    full = detection_distribution(gene, {"A": 0.8, "X": 0.6, "B": 0.7})
    bounded = detection_distribution(
        gene, {"A": 0.8, "X": 0.6, "B": 0.7}, max_observed_tips=2
    )
    from nwkit.mul_locus import signature_size

    expected = {key: p for key, p in full.items() if signature_size(key) <= 2}
    expected[OVERFLOW] = math.fsum(
        p for key, p in full.items() if signature_size(key) > 2
    )
    assert bounded == pytest.approx(expected, abs=3e-14)
    stratum = {"patterns": {}, "selected": Moments()}
    accumulate(stratum, bounded, 2)
    assert stratum["selected"].mean == pytest.approx(
        math.fsum(p for key, p in full.items() if signature_size(key) == 2)
    )
    assert OVERFLOW not in stratum["patterns"]


def test_large_hidden_genealogy_detection_retains_mass_within_observation_budget():
    gene = ("tip", "X", 0)
    for i in range(1, 20):
        gene = ("node", gene, ("tip", "X", i))
    result = detection_distribution(
        gene, {"X": 0.6}, max_observed_tips=4, max_states=1000
    )
    assert result[OVERFLOW] == pytest.approx(
        1 - sum(math.comb(20, k) * 0.6**k * 0.4 ** (20 - k) for k in range(5))
    )
    assert math.fsum(result.values()) == pytest.approx(1)


@pytest.mark.parametrize("value", [-1, float("nan"), float("inf"), True])
def test_invalid_detection_probabilities(value):
    with pytest.raises(ValueError, match="Detection"):
        detection_distribution(("tip", "A", 0), {"A": value})


def test_caps_abort_instead_of_selecting_or_retrying_histories():
    locus, _ = nested_locus_with_loss()
    with pytest.raises(ValueError, match="cap"):
        conditional_observations(
            locus, 0.7, {"A": 1, "B": 1, "X": 1, "C": 1}, max_states=5
        )
    with pytest.raises(ValueError, match="cap"):
        detection_distribution(
            ("node", ("tip", "X", 0), ("tip", "X", 1)), {"X": 0.5}, max_states=2
        )
    with pytest.raises(ValueError):
        WorkBudget(0)


def test_sparse_moments_include_unsampled_patterns_as_zeros():
    values = [0.4, 0.2, 0.7, 0, 0, 0, 0]
    moments = Moments()
    for value in values[:3]:
        moments.add(value)
    final = moments.with_zeros(len(values))
    assert final.mean == pytest.approx(np.mean(values))
    assert final.m2 / (final.n - 1) == pytest.approx(np.var(values, ddof=1))
    assert moments.n == 3


def test_empirical_bernstein_formula_and_missing_support():
    moments = Moments(100, 0.2, 0.5)
    alpha = 0.01
    radius = math.sqrt(2 * 0.5 / 99 * math.log(4 / alpha) / 100) + 7 * math.log(
        4 / alpha
    ) / (3 * 99)
    assert bounded_interval(moments, alpha) == pytest.approx(
        (max(0, 0.2 - radius), 0.2 + radius)
    )
    lo, hi = bounded_interval(Moments(100000, 0, 0), alpha)
    assert lo == 0 < hi < 0.001
    with pytest.raises(ValueError, match="two"):
        bounded_interval(Moments(1, 0, 0), alpha)


@pytest.mark.parametrize("mean", [0, 0.001, 0.2, 0.5, 0.999, 1])
def test_chernoff_inversion_endpoints_and_complement_symmetry(mean):
    lo, hi = chernoff_interval(Moments(500, mean, 0), 1e-7)
    inverse = chernoff_interval(Moments(500, 1 - mean, 0), 1e-7)
    assert (lo, hi) == pytest.approx((1 - inverse[1], 1 - inverse[0]))
    assert 0 <= lo <= mean <= hi <= 1
    if mean == 0:
        assert hi == pytest.approx(-math.expm1(-math.log(2e7) / 500))
        assert hi < bounded_interval(Moments(500, mean, 0), 1e-7)[1]
    if 0 < mean < 1:
        for endpoint in (lo, hi):
            if endpoint in (0, 1) or 1 - endpoint < 1e-12:
                continue  # The conservative outer bracket can be a float endpoint.
            kl = mean * math.log(mean / endpoint) + (1 - mean) * math.log(
                (1 - mean) / (1 - endpoint)
            )
            assert kl == pytest.approx(math.log(2e7) / 500)


def test_chernoff_exact_nonbernoulli_failure_mass():
    from scipy.stats import multinomial

    n, alpha = 32, 0.25
    probabilities, support = [0.8, 0.15, 0.05], [0, 0.3, 1]
    truth = sum(p * x for p, x in zip(probabilities, support, strict=True))
    failures = 0.0
    for a in range(n + 1):
        for b in range(n - a + 1):
            counts = [a, b, n - a - b]
            mean = sum(k * x for k, x in zip(counts, support, strict=True)) / n
            lo, hi = chernoff_interval(Moments(n, mean, 0), alpha)
            if not lo <= truth <= hi:
                failures += multinomial.pmf(counts, n, probabilities)
    assert 0 < failures <= alpha


@pytest.mark.parametrize("n", [True, 2.5, float("nan"), float("inf"), None])
@pytest.mark.parametrize("interval", [bounded_interval, chernoff_interval])
def test_invalid_moment_count_cannot_return_a_spurious_interval(n, interval):
    with pytest.raises(ValueError, match="two IID"):
        interval(Moments(n, 0.5, 0.1), 0.01)


@pytest.mark.parametrize("method", ["detection-rb", "hybrid-rb"])
def test_integrated_cli_bank_dispatch_and_simultaneous_budget(method):
    from nwkit.mul_locus_cli import bank_record
    from nwkit.mul_locus_mc import build_bank, integration_alpha, validate_model

    config = {**model(), "samples": 16, "integration": method}
    species = tree(SPECIES)
    parameters = validate_model(config, species)
    task = make_tasks(species, "X", "A B", parameters, 100)[0][0]
    bank = build_bank(task, config)
    assert isinstance(bank, IntegratedBank)
    assert bank.samples == bank.attempts == 16
    assert bank.interval_method == "chernoff-kl"
    assert bank_record(bank)["strata"][0]["selected"]["n"] > 0
    assert integration_alpha(config, 20) == integration_alpha(
        {**config, "integration": "ancestral-stratified"}, 20
    )


@pytest.mark.parametrize("value", [-0.1, 1.1, float("nan"), float("inf"), True])
def test_invalid_probability_moments(value):
    with pytest.raises(ValueError):
        Moments().add(value)


def test_selection_and_pattern_use_same_prior_weighted_denominator():
    signature = ("node", ("tip", "A"), ("tip", "B"))
    strata = []
    for weight, probability, selected in ((0.8, 0.2, 0.4), (0.1, 0.3, 0.9)):
        strata.append(
            {
                "weight": weight,
                "samples": 1000,
                "patterns": {signature: Moments(1000, probability, 1)},
                "selected": Moments(1000, selected, 1),
            }
        )
    bank = IntegratedBank(
        0, "NA", 0, None, LocusParameters(0.1, 0.1, 1, 0.5), strata, "detection"
    )
    estimate, lo, hi = integrated_probability(bank, signature, 0.01)
    assert estimate == pytest.approx((0.8 * 0.2 + 0.1 * 0.3) / (0.8 * 0.4 + 0.1 * 0.9))
    assert lo < estimate < hi
    absent = ("node", ("tip", "X"), ("tip", "X"))
    assert integrated_probability(bank, absent, 0.01)[0:2] == (0, 0)
    assert score_integrated_bank(bank, [absent], 0.01)["log_likelihood"] == -math.inf


def test_accumulation_excludes_large_observations_not_large_hidden_genealogies():
    signature = ("node", ("tip", "A"), ("tip", "B"))
    large = ("node", signature, ("tip", "X"))
    stratum = {"patterns": {}, "selected": Moments()}
    accumulate(stratum, {None: 0.1, signature: 0.4, large: 0.5}, 2)
    assert stratum["selected"].mean == 0.4
    assert set(stratum["patterns"]) == {signature}


@pytest.mark.parametrize("mode", ["plug-in", "grid-supremum"])
def test_custom_scorer_used_in_every_observed_and_null_search(mode):
    signature = ("node", ("tip", "A"), ("tip", "B"))
    point = LocusParameters(0, 0, 1, 0.5)
    banks = [
        LocusBank(c, "NA", i, None, point, Counter({signature: 100}), 100, 100)
        for i, c in enumerate((0, 0, 1))
    ]
    calls = []

    def scorer(bank, observations, alpha):
        calls.append(bank.grid)
        return score_bank(bank, observations, alpha)

    _, result = calibrate(
        banks,
        [signature],
        model(),
        0.01,
        2,
        null_calibration=mode,
        scorer=scorer,
        sampler=lambda *args, **kwargs: ([signature], 1),
    )
    assert len(calls) == (1 + result.get("total_replicates", result["replicates"])) * 3
    assert result["p_value"] == 1


def test_bounded_intervals_cover_exact_discrete_iid_law():
    # Enumerate all IID draws, including non-Bernoulli values, not a coverage simulation.
    n, alpha = 64, 0.25
    support, probabilities = (0.0, 0.2, 1.0), (0.1, 0.7, 0.2)
    expectation = math.fsum(x * p for x, p in zip(support, probabilities, strict=True))
    failure = 0.0
    for a in range(n + 1):
        for b in range(n - a + 1):
            c = n - a - b
            m = Moments()
            for x, k in zip(support, (a, b, c), strict=True):
                for _ in range(k):
                    m.add(x)
            lo, hi = bounded_interval(m, alpha)
            if not lo <= expectation <= hi:
                coefficient = math.factorial(n) / math.prod(
                    math.factorial(k) for k in (a, b, c)
                )
                failure += coefficient * math.prod(
                    p**k for p, k in zip(probabilities, (a, b, c), strict=True)
                )
    assert 0 < failure <= alpha


def test_paired_banks_retain_large_hidden_histories_and_repeat_exactly():
    point = LocusParameters(0, 0, 1, 0.5)
    task = make_tasks(tree(SPECIES), "X", "A B", [point], 100)[0][0]
    config = {
        **model(20),
        "max_observed_tips": 2,
        "detection": {"A": 0.7, "X": 0.7, "B": 0.7},
    }
    first, audit = build_paired_banks(task, config, exact_tip_limit=2)
    second, other = build_paired_banks(task, config, exact_tip_limit=2)
    assert audit == other == {"exact_histories": 0, "sampled_histories": 20}
    assert first["histogram"].strata == second["histogram"].strata
    assert (
        first["detection"].strata
        == first["hybrid"].strata
        == second["detection"].strata
    )
    assert first["hybrid"].strata[0]["selected"].mean > 0


def test_zero_dl_integrated_means_equal_exact_law_without_topology_mc_noise():
    point = LocusParameters(0, 0, 1, 0.5)
    task = make_tasks(tree(SPECIES), "X", "A B", [point], 100)[0][0]
    config = model(20)
    banks, audit = build_paired_banks(task, config)
    assert audit == {"exact_histories": 20, "sampled_histories": 0}
    locus = sample_locus_tree(task[3], point, 0.5, np.random.default_rng(7))
    expected = conditional_observations(locus, point.ne, config["detection"])
    stratum = banks["hybrid"].strata[0]
    for key, moment in stratum["patterns"].items():
        assert moment.mean == pytest.approx(expected[key], abs=3e-14)
        assert moment.m2 == 0
    assert stratum["selected"].n == 20


@pytest.mark.parametrize(
    "update", [{"samples": 1}, {"samples": 3}, {"samples": 20, "max_attempts": 1}]
)
def test_paired_bank_budget_rejected_before_sampling(update, monkeypatch):
    import nwkit.mul_locus_integral as module

    point = LocusParameters(0.1, 0.1, 1, 0.5)
    task = make_tasks(tree(SPECIES), "X", "A B", [point], 100)[0][0]
    monkeypatch.setattr(
        module,
        "sample_locus_tree",
        lambda *args, **kwargs: pytest.fail("invalid budget sampled a history"),
    )
    with pytest.raises(ValueError, match="budget"):
        build_paired_banks(task, {**model(), **update})


def test_accumulate_rejects_invalid_probability_mass():
    with pytest.raises(ArithmeticError, match="mass"):
        accumulate({"selected": Moments(), "patterns": {}}, {None: 1.2}, 4)
