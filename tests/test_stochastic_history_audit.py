"""Numerical-scale, coordinate, resource and independent CTMC map regressions."""

from types import SimpleNamespace

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.linalg import expm

from nwkit import asr
from nwkit.cli import main
from nwkit.stochastic_history import MapDraw, MapSample, _segments, sample_histories
from nwkit.stochastic_map_io import write_extended_maps
from nwkit.stochastic_map_tables import (
    duration_table,
    history_table,
    probability_table,
    time_table,
)
from nwkit.util import assign_branch_ids, read_tree


def single_branch_sample(length, fractions=(0.25, 0.75)):
    tree = read_tree("(A:1)R;", 1, False)
    node = tree.children[0]
    ids = assign_branch_ids(tree)
    draws = tuple(
        MapDraw(
            index,
            np.array([0, 1]),
            {ids[node]: ((0.0, fraction * length, 0), (fraction * length, length, 1))},
            {(ids[node], 0, 1): 1},
        )
        for index, fraction in enumerate(fractions, 1)
    )
    return tree, MapSample(
        ("a", "b"), draws, ids, {tree: 0.0, node: length}, {ids[node]: length}
    )


@pytest.mark.parametrize("length", [1e200, 1e-200])
def test_duration_statistics_are_scale_equivariant(length):
    _, sample = single_branch_sample(length)
    table = duration_table(sample)
    assert table.duration_sd.to_numpy() / length == pytest.approx([np.sqrt(0.125)] * 2)
    assert table.mean_duration.to_numpy() / length == pytest.approx([0.5, 0.5])
    assert table.total_duration.to_numpy() / length == pytest.approx([1, 1])
    assert table.duration_q025.to_numpy() / length == pytest.approx([0.2625, 0.2625])
    assert table.duration_q975.to_numpy() / length == pytest.approx([0.7375, 0.7375])


def test_unrepresentable_duration_totals_are_rejected():
    _, sample = single_branch_sample(1e308, (0.75,) * 4)
    with pytest.raises(ValueError, match="total duration overflows"):
        duration_table(sample)


def deep_sample():
    tree = read_tree("((A:4):10000000000000000)R;", 1, False)
    middle, node = tree.children[0], tree.children[0].children[0]
    ids = assign_branch_ids(tree)
    draw = MapDraw(
        1,
        np.array([0, 0, 1]),
        {ids[middle]: ((0.0, 1e16, 0),), ids[node]: ((0.0, 0.5, 0), (0.5, 4.0, 1))},
        {(ids[node], 0, 1): 1},
    )
    return MapSample(
        ("a", "b"),
        (draw,),
        ids,
        {tree: 0.0, middle: 1e16, node: 1e16 + 4},
        {ids[middle]: 1e16, ids[node]: 4.0},
    )


def test_time_bin_durations_use_local_coordinates_to_avoid_cancellation():
    sample = deep_sample()
    durations = time_table(sample, 1).query("quantity == 'duration' and state == 'b'")
    assert durations.total.iloc[0] == 3.5
    assert duration_table(sample).query("state == 'b'").total_duration.sum() == 3.5


@pytest.mark.parametrize(
    "builder", [history_table, lambda sample: probability_table(sample, 9)]
)
def test_root_coordinate_exports_reject_collapsed_times(builder):
    with pytest.raises(ValueError, match="not distinguishable"):
        builder(deep_sample())


def test_branch_depth_precision_is_checked_before_sampling():
    tree = read_tree("((A:3):10000000000000000)R;", 1, False)
    with pytest.raises(ValueError, match="loses precision"):
        sample_histories(tree, ["a", "b"], {}, 1)


def test_probability_grid_work_is_bounded_before_sampling(monkeypatch):
    tree = read_tree("(A:1)R;", 1, False)
    args = SimpleNamespace(
        n_sim=40000,
        threads=1,
        map_probabilities_out="unused.tsv",
        map_grid_points=10001,
    )

    def fail(*args, **kwargs):
        pytest.fail("Sampling started before the probability work budget was checked")

    monkeypatch.setattr("nwkit.stochastic_map_io.sample_histories", fail)
    with pytest.raises(ValueError, match="draw/grid evaluations"):
        write_extended_maps(tree, ["a", "b"], {}, args)


def test_figure_uses_literal_state_and_tip_labels(tmp_path):
    traits = tmp_path / "traits.tsv"
    traits.write_text("leaf_name\tstate\n$\\notacommand$\t$\\unknown$\nB\tblue\n")
    figure = tmp_path / "map.svg"
    main(
        [
            "asr",
            "-i",
            r"($\notacommand$:1,B:1);",
            "--trait",
            str(traits),
            "--state-column",
            "state",
            "--model",
            "ER",
            "--rate",
            ".5",
            "--n-sim",
            "2",
            "--seed",
            "4",
            "-o",
            str(tmp_path / "nodes.tsv"),
            "--map-figure-out",
            str(figure),
        ]
    )
    text = figure.read_text()
    assert r"$\notacommand$" in text and r"$\unknown$" in text


def test_asymmetric_three_state_bridge_matches_independent_integrals():
    q = np.array([[-1.1, 0.9, 0.2], [0.1, -0.7, 0.6], [0.4, 0.2, -0.6]])
    length, time = 1.7, 0.65
    start, end = 0, 2
    denominator = expm(q * length)[start, end]

    def conditional_probability(t, state):
        return (
            expm(q * t)[start, state] * expm(q * (length - t))[state, end] / denominator
        )

    expected_probability = [conditional_probability(time, state) for state in range(3)]
    expected_duration = [
        quad(conditional_probability, 0, length, args=(state,))[0] for state in range(3)
    ]
    context = asr._build_uniformization_context(q, length)
    rng, times_rng = np.random.default_rng(451), np.random.default_rng(671)
    probabilities, durations = np.zeros(3), np.zeros(3)
    count = 5000
    for _ in range(count):
        path = [start]
        asr._sample_bridge_transition_counts(
            start, end, length, q, rng, uniformization_context=context, history=path
        )
        segments = _segments(path, length, np.arange(3), times_rng)
        for left, right, state in segments:
            durations[state] += right - left
            if left <= time < right:
                probabilities[state] += 1
    assert probabilities / count == pytest.approx(expected_probability, abs=0.025)
    assert durations / count == pytest.approx(expected_duration, abs=0.025)


def test_expanded_state_simulation_arrays_are_bounded():
    from nwkit.stochastic_history import _history_work

    # No events: a loop/segment budget alone misses the dense CTMC arrays.
    fit = {"rate_matrix": np.zeros((64, 64))}
    branches = [SimpleNamespace(dist=0.0)] * 3000
    with pytest.raises(ValueError, match="simulation arrays exceed 256 MiB"):
        _history_work(fit, branches, 1)


def test_unrepresentable_time_bin_rates_are_rejected():
    _, sample = single_branch_sample(1e-310, (0.5,))
    with pytest.raises(ValueError, match="time-bin rate overflows"):
        time_table(sample, 1)


def test_numerical_failure_preserves_existing_asr_outputs(tmp_path, capsys):
    traits = tmp_path / "traits.tsv"
    traits.write_text("leaf_name\tstate\nA\tred\nB\tred\n")
    output, history, model = (
        tmp_path / name for name in ("nodes.tsv", "history.tsv", "model.tsv")
    )
    for path in (output, history, model):
        path.write_text("original")
    with pytest.raises(ValueError, match="loses precision"):
        main(
            [
                "asr",
                "-i",
                "((A:3):10000000000000000,B:10000000000000000);",
                "--trait",
                str(traits),
                "--state-column",
                "state",
                "--states",
                "red,blue",
                "--model",
                "ER",
                "--rate",
                "0",
                "--n-sim",
                "1",
                "--seed",
                "5",
                "-o",
                str(output),
                "--model-out",
                str(model),
                "--map-history-out",
                str(history),
            ]
        )
    assert all(path.read_text() == "original" for path in (output, history, model))
    assert capsys.readouterr().out == ""


def test_uniformization_cutoff_controls_rare_endpoint_error():
    from scipy.special import gammainc, gammaln, logsumexp
    from scipy.stats import poisson

    # Terminal hitting time is Erlang(63, .5). An unreachable fast state sets
    # omega=20, introducing many virtual self events along the slow chain.
    # The former unconditional cutoff (64) missed most conditional N mass.
    steps, slow, omega = 63, 0.5, 20.0
    q = np.zeros((65, 65))
    for state in range(steps):
        q[state, state] = -slow
        q[state, state + 1] = slow
    q[64, 64], q[64, 63] = -omega, omega
    _, lam, limit = asr._uniformization_parameters(q, 1.0)
    endpoint_probability = gammainc(steps, slow)
    assert poisson.sf(limit, lam) / endpoint_probability <= 1e-12
    expected = 1 - (steps / slow) * gammainc(steps + 1, slow) / endpoint_probability
    logs, occupations = [], []
    for n in range(steps, limit + 1):
        k = np.arange(steps, n + 1)
        log_weights = (
            poisson.logpmf(n, lam)
            + gammaln(k)
            - gammaln(steps)
            - gammaln(k - steps + 1)
            + steps * np.log(slow / omega)
            + (k - steps) * np.log1p(-slow / omega)
        )
        logs.extend(log_weights)
        occupations.extend(1 - k / (n + 1))
    actual = np.exp(np.asarray(logs) - logsumexp(logs)) @ occupations
    assert actual == pytest.approx(expected, rel=1e-10)
