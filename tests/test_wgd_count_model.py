import itertools
import math

import numpy as np
import pytest
from scipy.linalg import expm

from nwkit.wgd_count_model import (
    CountLikelihood,
    CountTree,
    MultiplicationEvent,
    _birth_death_log_transition,
    birth_death_parameters,
    birth_death_transition,
    multiplication_transition,
    rate_categories,
    root_probabilities,
)


def tree():
    return CountTree(
        (-1, 0, 0), (0.0, 0.4, 0.7), (1, 2), ("A", "B"), (0, 1, 2), ("AB", "A", "B")
    )


@pytest.mark.parametrize(
    "duplication,loss",
    [
        (0, 0),
        (0, 0.3),
        (0.2, 0),
        (0.2, 0.2),
        (0.3, 0.2),
        (0.2, 0.3),
        (0.2, 0.200000000001),
    ],
)
def test_branching_probabilities_against_independent_matrix_exponential(
    duplication, loss
):
    # A much larger killed CTMC is an independent numerical reference.
    bound = 120
    generator = np.zeros((bound + 1, bound + 1))
    for n in range(1, bound + 1):
        generator[n, n] = -n * (duplication + loss)
        generator[n, n - 1] = n * loss
        if n < bound:
            generator[n, n + 1] = n * duplication
    reference = expm(generator * 0.9)
    actual = birth_death_transition(duplication, loss, 0.9, 20)
    np.testing.assert_allclose(actual[:6], reference[:6, :21], atol=2e-14, rtol=2e-10)
    assert np.all(actual >= 0)
    assert np.all(actual.sum(axis=1) <= 1 + 1e-13)


@pytest.mark.parametrize(
    "duplication,loss", [(0, 0), (0, 0.3), (0.2, 0), (0.2, 0.2), (0.3, 0.2), (0.2, 0.3)]
)
def test_log_branching_kernel_matches_probability_kernel(duplication, loss):
    actual = _birth_death_log_transition(duplication, loss, 0.9, 20)
    expected = birth_death_transition(duplication, loss, 0.9, 20)
    np.testing.assert_allclose(np.exp(actual), expected, atol=2e-14, rtol=2e-12)
    assert not np.any(np.isnan(actual))


def test_high_copy_low_birth_likelihood_remains_finite_against_yule_law():
    duplication = 1e-8
    model = CountLikelihood(tree(), np.array([[84, 1]]), ascertainment="root-clades")
    assert birth_death_transition(duplication, 0, 0.4, 86)[1, 84] == 0
    expected = -duplication * 1.1 + 83 * np.log(-np.expm1(-duplication * 0.4))
    assert model.log_likelihood([[duplication, 0]], 1, 86) == pytest.approx(
        expected, abs=1e-11
    )
    assert model.convergence([[duplication, 0]], 1, 86) < 1e-11


def test_high_copy_rare_pulse_likelihood_matches_independent_yule_sum():
    duplication = 1e-8
    model = CountLikelihood(tree(), np.array([[84, 1]]), ascertainment="root-clades")
    event = MultiplicationEvent(1, 1)
    half_time = duplication * 0.2
    log_extra = math.log(-math.expm1(-half_time))
    terms = []
    for before in range(1, 43):
        after = 2 * before
        terms.append(
            -half_time
            + (before - 1) * log_extra
            + math.log(math.comb(83, after - 1))
            - after * half_time
            + (84 - after) * log_extra
        )
    peak = max(terms)
    expected = peak + math.log(sum(math.exp(term - peak) for term in terms))
    expected -= duplication * 0.7
    assert model.log_likelihood([[duplication, 0]], 1, 86, event) == pytest.approx(
        expected, abs=1e-11
    )


@pytest.mark.parametrize("duplication", [1e-8, 2.5652894719336476, 6.155496485662881])
@pytest.mark.parametrize("mean", [1.0, 2.5])
@pytest.mark.parametrize("multiplicity", [2, 3])
def test_loss_free_fully_detected_pulse_has_certain_selection(
    duplication, mean, multiplicity
):
    # Births and a pulse retain original copies. With no loss or detection error,
    # every positive root state observes both root clades with probability one.
    model = CountLikelihood(tree(), np.array([[1, 1]]), ascertainment="root-clades")
    event = MultiplicationEvent(1, 1.0, multiplicity=multiplicity)
    with np.errstate(invalid="raise"):
        actual = model._selection_logs([[duplication, 0.0]], 1.0, mean, event)
    np.testing.assert_array_equal(actual, [0.0])


@pytest.mark.parametrize("retention", [0.0, 0.3, 1.0])
def test_log_pruning_matches_scaled_pruning_with_pulse_and_missingness(retention):
    model = CountLikelihood(
        tree(), np.array([[1, 2], [0, 1], [1, np.nan]]), detection=np.array([0.8, 0.9])
    )
    rates, root_mean, bound = [[0.2, 0.3]], 1.5, 40
    event = MultiplicationEvent(1, retention)
    matrices = model._log_transitions(rates, 1, bound, event)
    unconditional = model._prune_log(matrices, root_mean, bound, np.arange(3))
    expected = model.family_log_likelihoods(rates, root_mean, bound, event)
    selection = model._selection_logs(rates, 1, root_mean, event)
    np.testing.assert_allclose(
        unconditional - selection[model.pattern_index], expected, atol=2e-12
    )


def test_time_units_and_semigroup():
    p = birth_death_transition(0.2, 0.4, 1.3, 100)
    np.testing.assert_allclose(
        p, birth_death_transition(0.0002, 0.0004, 1300, 100), atol=2e-14
    )
    composed = birth_death_transition(0.2, 0.4, 0.5, 100) @ birth_death_transition(
        0.2, 0.4, 0.8, 100
    )
    np.testing.assert_allclose(composed[:5], p[:5], atol=2e-14)


def test_event_limits_and_tail_not_renormalized():
    np.testing.assert_array_equal(multiplication_transition(0, 2, 8), np.eye(9))
    full = multiplication_transition(1, 3, 8)
    assert full[2, 6] == 1
    assert full[3].sum() == 0
    partial = multiplication_transition(0.5, 2, 8)
    np.testing.assert_allclose(partial[2, 2:5], [0.25, 0.5, 0.25])
    assert partial[8].sum() == pytest.approx(0.5**8)
    assert root_probabilities(4, 3).sum() == pytest.approx(1 - 0.75**3)


def test_pruning_matches_brute_force_and_ascertainment():
    model = CountLikelihood(
        tree(), np.array([[1, 2], [0, 1], [1, np.nan]]), detection=np.array([0.8, 0.9])
    )
    rates = [[0.2, 0.3]]
    bound = 40
    prior = root_probabilities(1.5, bound)
    matrices = model.transitions(rates, 1, bound)
    from scipy.stats import binom

    expected = []
    for observations in model.counts:
        probability = zero = 0.0
        for root in range(bound + 1):
            numerator = denominator = prior[root]
            for tip, value, detection in zip(
                (1, 2), observations, model.detection, strict=True
            ):
                if np.isnan(value):
                    numerator *= matrices[tip][root].sum()
                    denominator *= matrices[tip][root].sum()
                else:
                    numerator *= sum(
                        matrices[tip][root, n] * binom.pmf(value, n, detection)
                        for n in range(bound + 1)
                    )
                    denominator *= sum(
                        matrices[tip][root, n] * (1 - detection) ** n
                        for n in range(bound + 1)
                    )
            probability += numerator
            zero += denominator
        expected.append(np.log(probability / (1 - zero)))
    np.testing.assert_allclose(
        model.family_log_likelihoods(rates, 1.5, bound), expected, atol=1e-12
    )


def test_missing_is_not_zero_and_observation_conditioning_normalizes():
    bound = 40
    counts = np.array(list(itertools.product(range(10), repeat=2))[1:], dtype=float)
    model = CountLikelihood(tree(), counts)
    probabilities = np.exp(model.family_log_likelihoods([[0.1, 0.3]], 1, bound))
    assert probabilities.sum() == pytest.approx(1, abs=1e-9)
    missing = CountLikelihood(tree(), np.array([[1, np.nan], [1, 0]]))
    assert (
        missing.family_log_likelihoods([[0.1, 0.3]], 1, bound)[0]
        > missing.family_log_likelihoods([[0.1, 0.3]], 1, bound)[1]
    )


def test_rate_mixture_conditioned_after_mixing():
    model = CountLikelihood(tree(), np.array([[1, 0]]), rate_scales=np.array([0.1, 3]))
    components = [
        CountLikelihood(tree(), model.counts, rate_scales=np.array([scale]))
        for scale in model.rate_scales
    ]
    naive = np.log(
        np.mean([np.exp(m.log_likelihood([[0.2, 0.8]], 1, 40)) for m in components])
    )
    actual = model.log_likelihood([[0.2, 0.8]], 1, 40)
    assert abs(naive - actual) > 0.01


def test_zero_retention_matches_null_and_state_convergence():
    model = CountLikelihood(tree(), np.array([[1, 2], [2, 3]]))
    null = model.log_likelihood([[0.2, 0.3]], 1.5, 40)
    assert (
        model.log_likelihood([[0.2, 0.3]], 1.5, 40, MultiplicationEvent(1, 0)) == null
    )
    event = MultiplicationEvent(1, 0.7)
    assert model.convergence([[0.2, 0.3]], 1.5, 3, event) > 1e-7
    assert model.convergence([[0.2, 0.3]], 1.5, 40, event) < 1e-10


def test_simulation_seed_missingness_and_known_zero_rate_event():
    model = CountLikelihood(tree(), np.array([[1, np.nan], [1, 1]]))
    event = MultiplicationEvent(1, 1)
    one = model.simulate([[0, 0]], 1, np.random.default_rng(19), event)
    two = model.simulate([[0, 0]], 1, np.random.default_rng(19), event)
    np.testing.assert_array_equal(one, two)
    np.testing.assert_array_equal(one, [[2, np.nan], [2, 1]])
    simulated_model = CountLikelihood(tree(), one)
    assert simulated_model.log_likelihood([[0, 0]], 1, 8, event) == pytest.approx(0)
    assert model.log_likelihood([[0, 0]], 1, 8, event) == -np.inf


def test_simulated_frequency_matches_conditional_probability():
    model = CountLikelihood(tree(), np.ones((20000, 2)), detection=np.array([0.7, 0.8]))
    samples = model.simulate(
        [[0.2, 0.3]], 1, np.random.default_rng(531), MultiplicationEvent(1, 0.6)
    )
    target = CountLikelihood(tree(), np.array([[1, 1]]), detection=model.detection)
    expected = np.exp(
        target.log_likelihood([[0.2, 0.3]], 1, 40, MultiplicationEvent(1, 0.6))
    )
    assert np.mean(np.all(samples == [1, 1], axis=1)) == pytest.approx(
        expected, abs=0.015
    )
    assert np.all(samples.sum(axis=1) > 0)


@pytest.mark.parametrize(
    "counts",
    [[[0, 0]], [[np.nan, np.nan]], [[1, -1]], [[1, 0.5]], [[1, np.inf]], [[1]]],
)
def test_bad_counts(counts):
    with pytest.raises(ValueError):
        CountLikelihood(tree(), np.array(counts))


def test_invalid_model_arguments():
    with pytest.raises(ValueError):
        birth_death_parameters(-1, 0, 1)
    with pytest.raises(ValueError):
        birth_death_transition(1, 1, 1, 0)
    with pytest.raises(ValueError):
        multiplication_transition(1.01, 2, 4)
    with pytest.raises(ValueError):
        multiplication_transition(0.5, 2.5, 4)
    with pytest.raises(ValueError):
        root_probabilities(0.5, 4)
    with pytest.raises(ValueError):
        MultiplicationEvent(1, 0.5, 1.1)
    with pytest.raises(ValueError):
        CountLikelihood(tree(), np.array([[1, 2]]), detection=np.array([0, 1]))
    with pytest.raises(ValueError):
        CountLikelihood(tree(), np.array([[1, 2]]), branch_groups=(0, 0, 2))
    model = CountLikelihood(tree(), np.array([[1, 2]]))
    with pytest.raises(ValueError):
        model.log_likelihood([[0.2, 0.3]], 1, 1)
    with pytest.raises(ValueError):
        model.log_likelihood([[0.2, 0.3]], 1, 4, MultiplicationEvent(0, 0.5))


def test_gamma_categories():
    np.testing.assert_array_equal(rate_categories(None), [1])
    scales = rate_categories(1, 4)
    assert scales.mean() == pytest.approx(1)
    assert np.all(np.diff(scales) > 0)
    with pytest.raises(ValueError):
        rate_categories(0)


def test_root_clade_ascertainment_normalizes_and_matches_independent_selection():
    counts = np.array(list(itertools.product(range(1, 10), repeat=2)), dtype=float)
    model = CountLikelihood(tree(), counts, ascertainment="root-clades")
    rates = [[0.2, 0.4]]
    probabilities = np.exp(model.family_log_likelihoods(rates, 1, 40))
    assert probabilities.sum() == pytest.approx(1, abs=2e-8)
    branches = model.transitions(rates, 1, 40)
    selection = (1 - branches[1][1, 0]) * (1 - branches[2][1, 0])
    expected = branches[1][1, 1] * branches[2][1, 1] / selection
    assert probabilities[0] == pytest.approx(expected, abs=1e-13)
    with pytest.raises(ValueError, match="every root-child"):
        CountLikelihood(tree(), np.array([[1, 0]]), ascertainment="root-clades")


def test_root_clade_simulation_observation_and_missingness_selection():
    design = CountTree(
        (-1, 0, 1, 1, 0),
        (0.0, 0.5, 0.5, 0.5, 1.0),
        (2, 3, 4),
        ("A", "B", "C"),
        (0, 1, 2, 3, 4),
        ("ABC", "AB", "A", "B", "C"),
    )
    model = CountLikelihood(
        design,
        np.tile([1, np.nan, 1], (1000, 1)),
        ascertainment="root-clades",
        detection=np.array([0.7, 1, 0.7]),
    )
    result = model.simulate([[0.2, 0.5]], 1.2, np.random.default_rng(901))
    assert np.all(result[:, 0] > 0)
    assert np.all(result[:, 2] > 0)
    assert np.all(np.isnan(result[:, 1]))


@pytest.mark.parametrize("tips", [2, 3, 4])
@pytest.mark.parametrize("mean", [1.0, 2.0, 4.0])
@pytest.mark.parametrize("loss", [0.4, 40.0])
def test_root_clade_selection_matches_positive_geometric_sum(tips, mean, loss):
    from scipy.special import logsumexp

    design = CountTree(
        (-1, *(0 for _ in range(tips))),
        (0.0, *(1.0 for _ in range(tips))),
        tuple(range(1, tips + 1)),
        tuple(f"S{i}" for i in range(tips)),
        tuple(range(tips + 1)),
        tuple(f"clade{i}" for i in range(tips + 1)),
    )
    detection = np.linspace(0.6, 0.9, tips)
    model = CountLikelihood(
        design,
        np.ones((1, tips)),
        detection=detection,
        ascertainment="root-clades",
    )
    survival = detection * np.exp(-loss)
    copies = np.arange(1, 2001)
    if mean == 1:
        expected = np.log(survival).sum()
    else:
        log_prior = -np.log(mean) + (copies - 1) * np.log1p(-1 / mean)
        log_selected = np.log(-np.expm1(copies[:, None] * np.log1p(-survival))).sum(
            axis=1
        )
        expected = logsumexp(log_prior + log_selected)
    actual = model._selection_logs([[0.0, loss]], 1.0, mean, None)[0]
    assert actual == pytest.approx(expected, abs=1e-12)
    if loss == 40:
        # Conditional observation is almost surely one copy in every clade;
        # neither the root prior nor the selection denominator may be clipped.
        assert model.log_likelihood([[0.0, loss]], mean, 128) == pytest.approx(
            0.0, abs=1e-10
        )


@pytest.mark.parametrize(
    "nodes,names",
    [((1, 2, 2), ("A", "B", "B")), ((1, 1, 2), ("A", "A", "B")), ((1, 2), ("A", "A"))],
)
def test_duplicate_tip_mappings_are_rejected(nodes, names):
    with pytest.raises(ValueError, match="uniquely|unique"):
        CountTree(
            (-1, 0, 0), (0.0, 1.0, 1.0), nodes, names, (0, 1, 2), ("AB", "A", "B")
        )


@pytest.mark.parametrize("bound", [8, 16])
def test_fully_missing_subtree_integrates_exactly_without_truncated_tail(bound):
    design = CountTree(
        (-1, 0, 0),
        (0.0, 1.0, 100.0),
        (1, 2),
        ("A", "B"),
        (0, 1, 2),
        ("AB", "A", "B"),
    )
    model = CountLikelihood(design, np.array([[1.0, np.nan]]))
    # Independent Yule probability: one root copy is still one with exp(-lambda*t).
    for event in (None, MultiplicationEvent(2, 1.0)):
        assert model.log_likelihood([[0.2, 0.0]], 1.0, bound, event) == pytest.approx(
            -0.2, abs=1e-14
        )


def test_missing_subtree_is_marginalized_per_family_not_per_tip():
    design = CountTree(
        (-1, 0, 1, 1, 0),
        (0.0, 100.0, 1.0, 1.0, 1.0),
        (2, 3, 4),
        ("A", "B", "C"),
        (0, 1, 2, 3, 4),
        ("ABC", "AB", "A", "B", "C"),
    )
    counts = np.array([[np.nan, np.nan, 1.0], [1.0, np.nan, 1.0]])
    model = CountLikelihood(design, counts)
    logs = model.family_log_likelihoods([[0.2, 0.0]], 1.0, 8)
    assert logs[0] == pytest.approx(-0.2, abs=1e-14)
    assert logs[1] < -20


@pytest.mark.parametrize("node", [0, -1, True, 1.5])
def test_invalid_event_node_is_not_silently_ignored(node):
    with pytest.raises(ValueError, match="positive integer"):
        MultiplicationEvent(node, 0.5)


def test_zero_retention_state_check_cannot_hide_split_edge_truncation():
    design = CountTree(
        (-1, 0, 0), (0.0, 1.0, 1.0), (1, 2), ("A", "B"), (0, 1, 2), ("AB", "A", "B")
    )
    model = CountLikelihood(design, np.array([[1.0, 1.0]]))
    rates = [[20.0, 20.0]]
    zero = MultiplicationEvent(1, 0.0)
    tiny = MultiplicationEvent(1, 1e-12)
    # With one root copy and exact detection, endpoint counts are represented
    # even at bound eight, but intermediate pulse counts can exceed that bound.
    assert model.convergence(rates, 1.0, 8) < 1e-12
    assert model.convergence(rates, 1.0, 8, zero) > 0.5
    assert model.convergence(rates, 1.0, 128, zero) < 1e-7
    assert model.log_likelihood(rates, 1.0, 256, zero) == pytest.approx(
        model.log_likelihood(rates, 1.0, 256, tiny), abs=1e-10
    )
