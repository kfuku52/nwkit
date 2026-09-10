"""Independent CTMC bridge expectations, conservation and CLI map exports."""

from io import StringIO

import numpy as np
import pandas as pd
import pytest
from scipy.linalg import expm

from nwkit import asr
from nwkit.cli import main
from nwkit.stochastic_history import _segments, sample_histories
from nwkit.stochastic_map_tables import (
    duration_table,
    history_table,
    probability_table,
    time_table,
)
from nwkit.util import read_tree


@pytest.fixture(scope="module")
def problem():
    tree = read_tree("((A:1,B:0.7):0.6,C:1.8);", 0, False)
    states = ["red", "blue"]
    observed = {"A": "red", "B": "blue", "C": None}
    likelihoods = {
        name: np.ones(2)
        if value is None
        else np.array([float(value == state) for state in states])
        for name, value in observed.items()
    }
    _, fit = asr.compute_mk_marginals(
        tree, states, observed, likelihoods, model="ER", rate=0.8
    )
    return tree, states, fit


@pytest.fixture(scope="module")
def sample(problem):
    tree, states, fit = problem
    return sample_histories(tree, states, fit, 80, seed=192)


def test_segment_coverage_endpoints_and_shared_node_states(problem, sample):
    nodes = list(problem[0].traverse("preorder"))
    index = {node: i for i, node in enumerate(nodes)}
    for draw in sample.draws:
        for node in nodes[1:]:
            segments = draw.branches[sample.branch_ids[node]]
            assert segments[0][0] == 0 and segments[-1][1] == node.dist
            assert segments[0][2] == draw.node_states[index[node.up]]
            assert segments[-1][2] == draw.node_states[index[node]]
            assert sum(end - start for start, end, _ in segments) == pytest.approx(
                node.dist
            )
            for before, after in zip(segments, segments[1:], strict=False):
                assert before[1] == after[0] and before[2] != after[2]


def test_counts_are_identical_to_legacy_random_stream(problem, sample):
    tree, states, fit = problem
    original = asr._simulate_stochastic_maps(tree, states, fit, 80, seed=192)
    total, any_count = asr._merge_simulation_counts(
        draw.counts for draw in sample.draws
    )
    current = pd.DataFrame(
        asr._stochastic_map_rows(tree, states, sample.branch_ids, total, any_count, 80)
    )
    pd.testing.assert_frame_equal(original, current)


def test_thread_workers_are_reproducible(problem, sample):
    tree, states, fit = problem
    parallel = sample_histories(tree, states, fit, 80, seed=192, threads=2)
    pd.testing.assert_frame_equal(history_table(sample), history_table(parallel))


def test_duration_and_time_conservation(problem, sample):
    duration = duration_table(sample)
    expected = sum(sample.lengths.values())
    assert duration.mean_duration.sum() == pytest.approx(expected)
    assert duration.groupby(
        "branch_id"
    ).mean_fraction.sum().to_numpy() == pytest.approx(np.ones(4))
    bins = time_table(sample, 7)
    occupancy = bins[bins.quantity == "duration"]
    assert occupancy["mean"].sum() == pytest.approx(expected)
    assert occupancy.groupby("bin").per_lineage_time.sum().to_numpy() == pytest.approx(
        np.ones(7)
    )
    transitions = bins[bins.quantity == "transitions"]
    assert transitions.total.sum() == sum(
        sum(draw.counts.values()) for draw in sample.draws
    )


def test_probability_normalization_and_tip_conditioning(sample):
    probabilities = probability_table(sample, 11)
    assert probabilities.groupby(
        ["branch_id", "position"]
    ).probability.sum().to_numpy() == pytest.approx(np.ones(44))
    a = probabilities[(probabilities.name == "A") & (probabilities.position == 1)]
    assert a.loc[a.state == "red", "probability"].iloc[0] == 1
    assert a.loc[a.state == "blue", "probability"].iloc[0] == 0
    assert np.isfinite(probabilities.mc_se).all()


def test_conditional_bridge_against_matrix_exponential():
    # Condition on 0 at the parent and 1 at the child in a symmetric CTMC.
    # Independently compute the state probability at t using matrix exponentials.
    q = np.array([[-0.8, 0.8], [0.8, -0.8]])
    length, time = 1.3, 0.4
    denominator = expm(q * length)[0, 1]
    expected = expm(q * time)[0, 0] * expm(q * (length - time))[0, 1] / denominator
    context = asr._build_uniformization_context(q, length)
    rng, times = np.random.default_rng(837), np.random.default_rng(995)
    in_zero, durations = [], []
    for _ in range(4000):
        path = [0]
        counts = asr._sample_bridge_transition_counts(
            0, 1, length, q, rng, uniformization_context=context, history=path
        )
        segments = _segments(path, length, np.arange(2), times)
        assert sum(counts.values()) % 2 == 1
        in_zero.append(
            next(state for start, end, state in segments if start <= time < end) == 0
        )
        durations.append(
            sum(end - start for start, end, state in segments if state == 0)
        )
    assert np.mean(in_zero) == pytest.approx(expected, abs=0.025)
    assert np.mean(durations) == pytest.approx(length / 2, abs=0.025)


def test_hidden_projection_merges_unobservable_transitions():
    segments = _segments(
        [0, 1, 1, 2, 3, 0], 2.0, np.array([0, 0, 1, 1]), np.random.default_rng(4)
    )
    assert [s[2] for s in segments] == [0, 1, 0]
    assert sum(s[1] - s[0] for s in segments) == pytest.approx(2)
    assert _segments([1], 0.0, np.array([0, 0]), np.random.default_rng(3)) == (
        (0.0, 0.0, 0),
    )


def cli_args(tmp_path, **options):
    tree, traits = tmp_path / "tree.nwk", tmp_path / "traits.tsv"
    tree.write_text("((A:1,B:0.7):0.6,C:1.8);\n")
    traits.write_text("leaf_name\tstate\nA\tred\nB\tblue\nC\tNA\n")
    result = [
        "asr",
        "--infile",
        str(tree),
        "--trait",
        str(traits),
        "--state-column",
        "state",
        "--model",
        "ER",
        "--rate",
        ".8",
        "--states",
        "red,blue",
        "--n-sim",
        "20",
        "--seed",
        "13",
    ]
    for name, value in options.items():
        result += ["--" + name.replace("_", "-"), str(value)]
    return result


@pytest.mark.parametrize("extension", ["png", "svg", "pdf"])
def test_cli_all_outputs_and_figure(tmp_path, extension, capsys):
    options = {
        key: tmp_path / (key + ".tsv")
        for key in (
            "outfile",
            "model_out",
            "stochastic_map_out",
            "map_history_out",
            "map_summary_out",
            "map_time_out",
            "map_probabilities_out",
        )
    }
    options["map_figure_out"] = tmp_path / ("probabilities." + extension)
    main(cli_args(tmp_path, **options))
    for key, path in options.items():
        assert path.stat().st_size > 0
        if key != "map_figure_out":
            assert len(pd.read_csv(path, sep="\t")) > 0
    assert capsys.readouterr().out == ""
    if extension == "png":
        assert options["map_figure_out"].read_bytes().startswith(b"\x89PNG")


def test_map_export_failure_preserves_all_outputs(tmp_path, monkeypatch, capsys):
    import nwkit.stochastic_map_io as module

    output, history, figure = (
        tmp_path / name for name in ("nodes.tsv", "history.tsv", "figure.png")
    )
    for path in (output, history, figure):
        path.write_text("original")

    def fail(*args, **kwargs):
        raise OSError("figure export failed")

    monkeypatch.setattr(module, "draw_map_figure", fail)
    with pytest.raises(OSError, match="figure export failed"):
        main(
            cli_args(
                tmp_path, outfile=output, map_history_out=history, map_figure_out=figure
            )
        )
    assert all(path.read_text() == "original" for path in (output, history, figure))
    assert capsys.readouterr().out == ""


def test_history_only_allows_simulation_controls_and_stdout(tmp_path, capsys):
    history = tmp_path / "history.tsv"
    main(cli_args(tmp_path, map_history_out=history))
    output = capsys.readouterr().out
    assert output.count("branch_id\t") == 1
    assert len(pd.read_csv(StringIO(output), sep="\t")) > 0


@pytest.mark.parametrize(
    "key,value,match",
    [
        ("n_sim", 0, "positive"),
        ("threads", 0, "positive"),
        ("map_time_bins", 0, "between"),
        ("map_grid_points", 1, "between"),
        ("map_figure_out", "bad.txt", "PNG"),
    ],
)
def test_invalid_controls_rejected_before_output(tmp_path, key, value, match):
    options = {
        "map_history_out": tmp_path / "history.tsv",
        "map_time_out": tmp_path / "time.tsv",
        "map_probabilities_out": tmp_path / "prob.tsv",
        key: value,
    }
    with pytest.raises(ValueError, match=match):
        main(cli_args(tmp_path, **options))
    assert not (tmp_path / "history.tsv").exists()


def test_auxiliary_path_cannot_replace_input(tmp_path):
    args = cli_args(tmp_path, map_history_out=tmp_path / "traits.tsv")
    before = (tmp_path / "traits.tsv").read_bytes()
    with pytest.raises(ValueError, match="overwrite input"):
        main(args)
    assert (tmp_path / "traits.tsv").read_bytes() == before


@pytest.mark.parametrize("model", ["HRM", "COVARION", "MK-REGIME", "PAGEL-INDEPENDENT"])
def test_extended_models_preserve_legacy_counts(tmp_path, model):
    from nwkit.util import assign_branch_ids

    tree_text = "[&R]((A:1,B:1):1,(C:1,D:1):1,E:2,F:2)R;"
    trait = tmp_path / "traits.tsv"
    trait.write_text(
        "leaf_name\tstate\tsecond\nA\t0\t0\nB\t0\t1\nC\t1\t0\nD\t1\t1\nE\t0\t0\nF\t1\t1\n"
    )
    base = [
        "asr",
        "-i",
        tree_text,
        "--trait-type",
        "discrete",
        "--trait",
        str(trait),
        "--state-column",
        "state,second" if model.startswith("PAGEL") else "state",
        "--model",
        model,
        "--rate-bounds",
        "0.01,3",
        "--n-sim",
        "12",
        "--seed",
        "37",
        "-o",
        str(tmp_path / "nodes.tsv"),
    ]
    if model == "MK-REGIME":
        tree = read_tree(tree_text, 1, False)
        regimes = tmp_path / "regimes.tsv"
        regimes.write_text(
            "branch_id\tregime\n"
            + "".join(
                f"{identifier}\t{'foreground' if node.name in {'C', 'D'} else 'background'}\n"
                for node, identifier in assign_branch_ids(tree).items()
            )
        )
        base += ["--regime-map", str(regimes)]
    legacy, current, history = (
        tmp_path / name for name in ("legacy.tsv", "current.tsv", "history.tsv")
    )
    main(base + ["--stochastic-map-out", str(legacy)])
    main(
        base + ["--stochastic-map-out", str(current), "--map-history-out", str(history)]
    )
    assert legacy.read_bytes() == current.read_bytes()
    records = pd.read_csv(history, sep="\t", dtype={"state": str})
    if model in {"HRM", "COVARION"}:
        assert set(records.state) == {"0", "1"}
    changes = records.groupby(["simulation", "branch_id"]).size() - 1
    assert changes.sum() == pd.read_csv(current, sep="\t").total_count.sum()


def test_zero_rate_zero_length_and_root_only():
    for source in ("(A:0,B:1)R;", "A;"):
        tree = read_tree(source, 1, False)
        observed = {tip.name: "a" for tip in tree.leaves()}
        likelihoods = {name: np.array([1.0, 0.0]) for name in observed}
        _, fit = asr.compute_mk_marginals(
            tree, ["a", "b"], observed, likelihoods, rate=0
        )
        sample = sample_histories(tree, ["a", "b"], fit, 1, seed=7)
        history, duration, probabilities = (
            history_table(sample),
            duration_table(sample),
            probability_table(sample, 3),
        )
        assert len(history) == len(sample.lengths)
        assert sum(sample.draws[0].counts.values()) == 0
        if sample.lengths:
            assert set(history.state) == {"a"}
            assert duration[duration.branch_length == 0].mean_fraction.tolist() == [
                "",
                "",
            ]
            assert set(duration.duration_sd) == {""}
            assert time_table(sample, 3).query("quantity == 'duration'")[
                "mean"
            ].sum() == pytest.approx(1)
        else:
            assert history.empty and duration.empty and probabilities.empty
            with pytest.raises(ValueError, match="positive tree time span"):
                time_table(sample, 3)


def test_time_boundary_is_assigned_to_following_bin():
    from nwkit.stochastic_history import MapDraw, MapSample
    from nwkit.util import assign_branch_ids

    tree = read_tree("(A:2)R;", 1, False)
    child = tree.children[0]
    ids = assign_branch_ids(tree)
    draw = MapDraw(
        1,
        np.array([0, 1]),
        {ids[child]: ((0.0, 1.0, 0), (1.0, 2.0, 1))},
        {(ids[child], 0, 1): 1},
    )
    sample = MapSample(
        ("a", "b"), (draw,), ids, {tree: 0.0, child: 2.0}, {ids[child]: 2.0}
    )
    transitions = time_table(sample, 2).query(
        "quantity == 'transitions' and state == 'a'"
    )
    assert transitions.total.tolist() == [0, 1]
    assert (
        probability_table(sample, 3)
        .query("position == 0.5 and state == 'b'")
        .probability.iloc[0]
        == 1
    )


@pytest.mark.parametrize(
    "model,trait_type",
    [("BM", "continuous"), ("MK-MIXTURE", "discrete"), ("THRESHOLD", "discrete")],
)
def test_unsupported_models_reject_map_outputs(tmp_path, model, trait_type):
    history = tmp_path / "history.tsv"
    with pytest.raises(ValueError):
        main(
            cli_args(
                tmp_path, model=model, trait_type=trait_type, map_history_out=history
            )
        )
    assert not history.exists()
