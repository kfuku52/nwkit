import itertools
import math

import numpy as np
import pytest
from scipy.linalg import expm
from scipy.special import logsumexp

from nwkit import wgd_count_model
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
    "retention,split_zero,expected_calls",
    [(0.0, False, 0), (0.0, True, 1), (0.4, False, 1), (1.0, False, 1)],
)
def test_likelihood_reuses_the_same_event_matrix_across_rate_scales(
    monkeypatch, retention, split_zero, expected_calls
):
    model = CountLikelihood(tree(), [[1, 2], [1, np.nan]], rate_scales=[0.25, 1, 2])
    event = MultiplicationEvent(1, retention, 0.5, 3)
    rates, root, bound = [[0.12, 0.2]], 1.2, 8
    prior = root_probabilities(root, bound)
    emissions = model._emissions(bound)
    observed = [
        model._prune(
            emissions,
            model.transitions(rates, scale, bound, event, _split_zero_event=split_zero),
            prior,
        )
        for scale in model.rate_scales
    ]
    selected = [
        model._selection_logs(rates, scale, root, event) for scale in model.rate_scales
    ]
    expected = (logsumexp(observed, axis=0) - np.log(3)) - (
        logsumexp(selected, axis=0) - np.log(3)
    )[model.pattern_index]
    calls = []

    def counted(*args):
        calls.append(args)
        return multiplication_transition(*args)

    monkeypatch.setattr(wgd_count_model, "multiplication_transition", counted)
    actual = model.family_log_likelihoods(
        rates, root, bound, event, _split_zero_event=split_zero
    )
    np.testing.assert_array_equal(actual, expected)
    assert calls == [(retention, 3, bound)] * expected_calls


def test_likelihood_transition_cache_is_bounded_and_preserves_results(monkeypatch):
    model = CountLikelihood(tree(), np.array([[2.0, 3.0], [1.0, np.nan]]))
    rates = np.array([[0.12, 0.17]])
    event = MultiplicationEvent(1, 0.3, 0.5, 2)
    bound = 8
    matrix_bytes = (bound + 1) ** 2 * np.dtype(float).itemsize
    monkeypatch.setattr(wgd_count_model, "_TRANSITION_CACHE_BYTES", 0)
    expected = model.family_log_likelihoods(rates, 1.4, bound, event)
    assert not model._transition_cache
    monkeypatch.setattr(wgd_count_model, "_TRANSITION_CACHE_BYTES", 2 * matrix_bytes)
    np.testing.assert_array_equal(
        model.family_log_likelihoods(rates, 1.4, bound, event), expected
    )
    assert len(model._transition_cache) == 2
    assert model._transition_cache_bytes == 2 * matrix_bytes
    # Public kernels remain independent mutable arrays, even on a cache hit.
    matrix = model.transitions(rates, 1.0, bound)[2]
    matrix[:] = np.nan
    np.testing.assert_array_equal(
        model.family_log_likelihoods(rates, 1.4, bound, event), expected
    )
    # Oversized entries do not displace the cached working set.
    model._cached_transition(*rates[0], 0.7, 20)
    assert model._transition_cache_bytes == 2 * matrix_bytes


def test_transition_cache_separates_log_kernels_and_count_bounds(monkeypatch):
    model = CountLikelihood(tree(), np.array([[1.0, 1.0]]))
    monkeypatch.setattr(wgd_count_model, "_TRANSITION_CACHE_ENTRIES", 2)
    probability = model._cached_transition(0.2, 0.3, 0.5, 8)
    log = model._cached_transition(0.2, 0.3, 0.5, 8, log=True)
    np.testing.assert_allclose(np.exp(log), probability, rtol=1e-12, atol=1e-14)
    assert model._cached_transition(0.2, 0.3, 0.5, 8) is probability
    larger = model._cached_transition(0.2, 0.3, 0.5, 16)
    assert larger.shape == (17, 17)
    assert list(model._transition_cache.values())[0] is probability
    assert model._transition_cache_bytes == probability.nbytes + larger.nbytes


def test_emissions_reuse_one_bound_and_refresh_on_observation_changes(monkeypatch):
    counts = np.array([[1.0, 2.0], [np.nan, 1.0]])
    detection = np.array([0.9, 1.0])
    model = CountLikelihood(tree(), counts, detection=detection)
    first = model._emissions(8)
    assert model._emissions(8) is first
    assert model._emissions(16) is not first
    second = model._emissions(8)
    assert second is not first
    np.testing.assert_array_equal(second[1], first[1])
    counts[0, 0] = 2.0
    detection[1] = 0.8
    updated = model._emissions(8)
    assert updated is not second
    assert np.any(updated[1] != second[1])
    np.testing.assert_array_equal(second[2], [2, 1])
    assert updated[2].shape == (2, 9)
    assert np.all(updated[2][np.arange(2), second[2]] < 1)
    monkeypatch.setattr(wgd_count_model, "_EMISSION_CACHE_BYTES", 0)
    uncached = model._emissions(8)
    assert model._emission_cache is None
    for node in model.tree.tip_nodes:
        np.testing.assert_array_equal(uncached[node], updated[node])


@pytest.mark.parametrize("detection", [[1, 1], [0.82, 1], [1, 0.9]])
@pytest.mark.parametrize("retention", [0.0, 0.3])
def test_deterministic_tip_indices_match_dense_emissions(detection, retention):
    counts = np.array([[0, 1], [3, 0], [np.nan, 2], [2, np.nan], [1, 3]])
    model = CountLikelihood(tree(), counts, detection=np.array(detection))
    bound = 8
    rates = np.array([[0.12, 0.2]])
    event = MultiplicationEvent(2, retention, 0.5, 2)
    transitions = model.transitions(rates, 1, bound, event)
    emissions = model._emissions(bound)
    assert [emissions[i].ndim for i in tree().tip_nodes] == [
        1 if p == 1 else 2 for p in detection
    ]
    dense = {}
    for column, node in enumerate(tree().tip_nodes):
        observed_counts = counts[:, column]
        emission = (observed_counts[:, None] == np.arange(bound + 1)).astype(float)
        emission[np.isnan(observed_counts)] = 1
        dense[node] = emission if detection[column] == 1 else emissions[node]
    prior = root_probabilities(1.3, bound)
    np.testing.assert_array_equal(
        model._prune(emissions, transitions, prior),
        model._prune(dense, transitions, prior),
    )


@pytest.mark.parametrize("ascertainment", ["observed", "root-clades"])
@pytest.mark.parametrize("detection", [[1, 0.93, 0.8], [1e-180, 1e-190, 1e-170]])
@pytest.mark.parametrize("loss", [0.15, 1700])
def test_shared_missingness_selection_matches_batched_path(
    ascertainment, detection, loss
):
    design = CountTree(
        (-1, 0, 1, 1, 0),
        (0.0, 0.1, 0.4, 0.5, 0.6),
        (2, 3, 4),
        ("A", "B", "C"),
        (0, 1, 2, 3, 4),
        ("ABC", "AB", "A", "B", "C"),
    )
    observations = [[np.nan, 1.0, 1.0], [np.nan, 2.0, 1.0]]
    scalar = CountLikelihood(
        design, observations, detection=detection, ascertainment=ascertainment
    )
    batched = CountLikelihood(
        design,
        [*observations, [1.0, 2.0, 1.0]],
        detection=detection,
        ascertainment=ascertainment,
    )
    which = batched.pattern_index[0]
    for event in (
        None,
        MultiplicationEvent(1, 0.4, 0.25, 3),
        MultiplicationEvent(3, 1.0, 0, 2),
    ):
        args = (np.array([[0.05, loss]]), 1.8, 1.2, event)
        np.testing.assert_array_equal(
            scalar._selection_logs(*args), batched._selection_logs(*args)[[which]]
        )


@pytest.mark.parametrize("budget", [1, 1024**2])
@pytest.mark.parametrize("families", [1, 67])
def test_log_pruning_tiles_preserve_impossible_states_and_reduction_order(
    budget, families, monkeypatch
):
    rng = np.random.default_rng(303)
    matrix = rng.normal(-100, 20, (23, 23))
    values = rng.normal(-500, 200, (families, 23))
    matrix[0] = -np.inf
    values[0] = -np.inf
    expected = np.array([logsumexp(matrix + row[None, :], axis=1) for row in values])
    monkeypatch.setattr(wgd_count_model, "_LOG_PRUNING_WORKSPACE_BYTES", budget)
    np.testing.assert_array_equal(
        wgd_count_model._log_pruning_contribution(matrix, values), expected
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


@pytest.mark.parametrize(
    "rates_time",
    [
        (0, 0, 10),
        (0, 0.02, 10),
        (0.03, 0, 10),
        (0.02, 0.02, 40),
        (0.01, 0.03, 100),
        (0.03, 0.01, 100),
        (1e-8, 0.4, 500),
        (0.4, 1e-8, 500),
    ],
)
def test_large_count_geometric_recurrence_matches_convolution_and_log_kernel(
    rates_time, monkeypatch
):
    bound = 275
    actual = birth_death_transition(*rates_time, bound)
    monkeypatch.setattr(wgd_count_model, "_GEOMETRIC_RECURRENCE_MIN_COUNT", bound + 1)
    expected = birth_death_transition(*rates_time, bound)
    assert np.all(actual >= 0)
    assert np.all(actual.sum(axis=1) <= 1 + 1e-12)
    np.testing.assert_allclose(actual, expected, rtol=4e-12, atol=5e-323)
    representable = actual > 1e-280
    log_probs = _birth_death_log_transition(*rates_time, bound)
    np.testing.assert_allclose(
        np.log(actual[representable]), log_probs[representable], rtol=0, atol=2e-11
    )


@pytest.mark.parametrize(
    "duplication,loss", [(0.2, 0.3), (0.2, 0.2), (0.2, 0), (0, 0.3), (1e-8, 0)]
)
def test_log_branching_prefix_does_not_depend_on_global_state_bound(duplication, loss):
    larger = _birth_death_log_transition(duplication, loss, 0.9, 90)
    for bound in (1, 4, 9):
        smaller = _birth_death_log_transition(duplication, loss, 0.9, bound)
        np.testing.assert_allclose(
            smaller, larger[: bound + 1, : bound + 1], rtol=1e-13, atol=1e-12
        )


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
@pytest.mark.parametrize("detection", [(0.8, 0.9), (1.0, 0.9), (1.0, 1.0)])
def test_log_pruning_matches_scaled_pruning_with_pulse_and_missingness(
    retention, detection
):
    model = CountLikelihood(
        tree(),
        np.array([[1, 2], [0, 1], [1, np.nan], [np.nan, 1]]),
        detection=np.array(detection),
    )
    rates, root_mean, bound = [[0.2, 0.3]], 1.5, 40
    event = MultiplicationEvent(1, retention)
    matrices = model._log_transitions(rates, 1, bound, event)
    unconditional = model._prune_log(matrices, root_mean, bound, np.arange(4))
    expected = model.family_log_likelihoods(rates, root_mean, bound, event)
    selection = model._selection_logs(rates, 1, root_mean, event)
    np.testing.assert_allclose(
        unconditional - selection[model.pattern_index], expected, atol=2e-12
    )


@pytest.mark.parametrize("fraction", [0.0, 0.5, 1.0])
@pytest.mark.parametrize("retention,multiplicity", [(0, 2), (0.4, 3), (1, 2), (1, 12)])
def test_log_pulse_transition_matches_explicit_truncated_matrices(
    fraction, retention, multiplicity
):
    model = CountLikelihood(tree(), np.array([[1, 1]]))
    rates, bound = [[0.6, 0.2]], 7
    event = MultiplicationEvent(1, retention, fraction, multiplicity)
    actual = model._log_transitions(rates, 1, bound, event, _split_zero_event=True)
    time = model.tree.lengths[1]
    before = birth_death_transition(*rates[0], time * fraction, bound)
    after = birth_death_transition(*rates[0], time * (1 - fraction), bound)
    pulse = np.zeros((bound + 1, bound + 1))
    for count in range(bound + 1):
        trials = count * (multiplicity - 1)
        for extra in range(min(trials, bound - count) + 1):
            pulse[count, count + extra] = (
                math.comb(trials, extra)
                * retention**extra
                * (1 - retention) ** (trials - extra)
            )
    expected = before @ pulse @ after
    np.testing.assert_array_equal(np.isfinite(actual[1]), expected > 0)
    np.testing.assert_allclose(np.exp(actual[1]), expected, rtol=2e-13, atol=1e-15)


def test_time_units_and_semigroup():
    p = birth_death_transition(0.2, 0.4, 1.3, 100)
    np.testing.assert_allclose(
        p, birth_death_transition(0.0002, 0.0004, 1300, 100), atol=2e-14
    )
    composed = birth_death_transition(0.2, 0.4, 0.5, 100) @ birth_death_transition(
        0.2, 0.4, 0.8, 100
    )
    np.testing.assert_allclose(composed[:5], p[:5], atol=2e-14)


@pytest.mark.parametrize("bound", [8, 70])
def test_event_limits_and_tail_not_renormalized(bound):
    np.testing.assert_array_equal(
        multiplication_transition(0, 2, bound), np.eye(bound + 1)
    )
    full = multiplication_transition(1, 3, bound)
    assert full[2, 6] == 1
    assert full[bound // 3 + 1].sum() == 0
    partial = multiplication_transition(0.5, 2, bound)
    np.testing.assert_allclose(partial[2, 2:5], [0.25, 0.5, 0.25])
    assert partial[bound].sum() == pytest.approx(0.5**bound, rel=1e-12, abs=0)
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


def test_missing_patterns_keep_tip_order_across_byte_boundary():
    design = CountTree(
        (-1, 0, *([1] * 8), 0),
        (0.0, *([1.0] * 10)),
        tuple(range(10, 1, -1)),
        tuple(f"T{i}" for i in range(9)),
        tuple(range(11)),
        tuple(f"C{i}" for i in range(11)),
    )
    rows = np.ones((5, 9))
    rows[[0, 2], -1] = np.nan
    rows[1, 0] = np.nan
    rows[4, [3, 8]] = np.nan
    missing = np.isnan(rows)
    expected_patterns, expected_index = np.unique(missing, axis=0, return_inverse=True)
    model = CountLikelihood(design, rows)

    assert model.root_clade_columns == [list(range(1, 9)), [0]]
    assert model.patterns.dtype == bool
    np.testing.assert_array_equal(model.patterns, expected_patterns)
    np.testing.assert_array_equal(model.pattern_index, expected_index)


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


@pytest.mark.parametrize("loss", [740.0, 1000.0])
def test_rare_observed_selection_preserves_pulse_mixture_and_missingness(loss):
    design = CountTree(
        (-1, 0, 0, 0),
        (0.0, 0.0, 1.0, 1.0),
        (1, 2, 3),
        ("A", "B", "C"),
        (0, 1, 2, 3),
        ("ABC", "A", "B", "C"),
    )
    model = CountLikelihood(
        design,
        np.array([[np.nan, 1, 0], [1, np.nan, np.nan]]),
        rate_scales=np.array([1.0, 1.5]),
    )
    event = MultiplicationEvent(3, 0.5)
    actual = model.family_log_likelihoods([[0.0, loss]], 2.0, 96, event)
    # In this rare-survival limit a single copy survives. At either rate,
    # its relative chances in B and C are 1 and 1.5 because C has a pulse.
    # With only the zero-length A observed, the likelihood is the root prior.
    np.testing.assert_allclose(actual, np.log([1 / 2.5, 1 / 2]), atol=1e-12)


@pytest.mark.parametrize("rare_clades,root_moment", [(1, 2), (2, 6)])
def test_rare_root_clade_selection_preserves_size_biased_root(rare_clades, root_moment):
    tips = 1 + rare_clades
    design = CountTree(
        (-1,) + (0,) * tips,
        (0.0, 0.0) + (1.0,) * rare_clades,
        tuple(range(1, tips + 1)),
        tuple(f"S{i}" for i in range(tips)),
        tuple(range(tips + 1)),
        tuple(f"clade{i}" for i in range(tips + 1)),
    )
    counts = np.ones((2, tips))
    counts[1, 0] = 2
    model = CountLikelihood(design, counts, ascertainment="root-clades")
    # Observing all rare clades weights each possible root count n by
    # n**rare_clades. The first two moments of Geometric(1/2) are 2 and 6.
    expected = np.log([0.5 / root_moment, 0.25 * 2**rare_clades / root_moment])
    actual = model.family_log_likelihoods([[0.0, 1000.0]], 2.0, 8)
    np.testing.assert_allclose(actual, expected, atol=1e-12)


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
