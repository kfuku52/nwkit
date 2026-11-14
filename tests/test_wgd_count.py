import json
from dataclasses import replace
from itertools import product

import numpy as np
import pandas as pd
import pytest

from nwkit.cli import main
from nwkit.wgd_count import read_counts


def inputs(tmp_path):
    tree = tmp_path / "species.nwk"
    tree.write_text("((A:0.5,B:0.5):0.5,C:1);")
    counts = tmp_path / "counts.tsv"
    counts.write_text(
        "family_id\tA\tB\tC\n" + "".join(f"f{i}\t2\t2\t1\n" for i in range(10))
    )
    return tree, counts


def command(tree, counts):
    return [
        "wgd-count",
        "-i",
        str(tree),
        "--counts",
        str(counts),
        "--candidate-branches",
        "1",
        "--event-fractions",
        "0.5",
        "--family-gamma-shape",
        "none",
        "--rate-model",
        "homogeneous",
        "--max-states",
        "64",
    ]


def test_real_cli_and_metadata_no_uncalibrated_probability(tmp_path):
    tree, counts = inputs(tmp_path)
    out = tmp_path / "events.tsv"
    model = tmp_path / "model.json"
    main(command(tree, counts) + ["-o", str(out), "--model-out", str(model)])
    table = pd.read_csv(out, sep="\t", keep_default_na=False)
    assert list(table.branch_id) == [1]
    assert table.retention.iloc[0] > 0.99
    assert table.count_support.iloc[0] == "not_calibrated"
    assert table.p_value.iloc[0] == "NA"
    assert table.species_event_id.iloc[0].startswith("clade-sha256:")
    metadata = json.loads(model.read_text())
    assert metadata["schema_version"] == 1
    assert metadata["calibration"] == "not-run"
    assert metadata["background"]["state_converged"]
    assert "after mixing" in metadata["ascertainment"]
    assert "not posterior" in metadata["uncertainty"]
    assert len(metadata["candidates"][0]["branch_burst_rate_groups"]) == len(
        metadata["species_event_ids"]
    )
    assert (
        bool(table.nuisance_bound_reached.iloc[0])
        == metadata["candidates"][0]["event"]["nuisance_bound_reached"]
    )
    assert (
        bool(table.background_nuisance_bound_reached.iloc[0])
        == metadata["background"]["nuisance_bound_reached"]
    )
    assert (
        bool(table.branch_burst_nuisance_bound_reached.iloc[0])
        == metadata["candidates"][0]["branch_burst"]["nuisance_bound_reached"]
    )


@pytest.mark.parametrize(
    "event_bound,background_bound,burst_bound", list(product((False, True), repeat=3))
)
def test_cli_exposes_each_fit_boundary_without_changing_support_status(
    tmp_path, monkeypatch, event_bound, background_bound, burst_bound
):
    from nwkit.wgd_count_fit import CountFit, CountScan, ScanCandidate
    from nwkit.wgd_count_model import MultiplicationEvent

    background = CountFit(
        np.array([[0.1, 0.2]]),
        1.0,
        -0.1,
        3,
        8,
        0.0,
        True,
        background_bound,
        "Verified test fit",
    )
    event = replace(
        background,
        event=MultiplicationEvent(1, 0.8),
        log_likelihood=0.0,
        num_parameters=4,
        boundary=event_bound,
    )
    burst = replace(
        background,
        rates=np.tile(background.rates, (2, 1)),
        log_likelihood=-0.3,
        num_parameters=5,
        boundary=burst_bound,
    )
    candidate = ScanCandidate(1, event, burst, 0.2, burst.aic - event.aic)
    scan = CountScan(background, (candidate,), fractions=(0.5,))
    monkeypatch.setattr("nwkit.wgd_count.scan_counts", lambda *args, **kwargs: scan)
    calibrated = replace(
        scan,
        candidates=(replace(candidate, p_value=0.01, p_value_mc_se=0.01),),
        bootstrap_statistics=(0.0,) * 99,
        calibration="plugin-parametric-bootstrap-search-maximum",
    )
    monkeypatch.setattr(
        "nwkit.wgd_count.calibrate_scan", lambda *args, **kwargs: calibrated
    )
    tree, counts = inputs(tmp_path)
    out, model = tmp_path / "events.tsv", tmp_path / "model.json"
    main(
        command(tree, counts)
        + ["--bootstrap", "99", "-o", str(out), "--model-out", str(model)]
    )
    row = pd.read_csv(out, sep="\t").iloc[0]
    assert bool(row.nuisance_bound_reached) is event_bound
    assert bool(row.background_nuisance_bound_reached) is background_bound
    assert bool(row.branch_burst_nuisance_bound_reached) is burst_bound
    assert row.count_support == "count_supported_conditional"
    metadata = json.loads(model.read_text())
    assert metadata["background"]["nuisance_bound_reached"] is background_bound
    assert metadata["candidates"][0]["event"]["nuisance_bound_reached"] is event_bound
    assert (
        metadata["candidates"][0]["branch_burst"]["nuisance_bound_reached"]
        is burst_bound
    )


def test_real_bootstrap_cli_respects_seed_and_search(tmp_path):
    tree, counts = inputs(tmp_path)
    first = tmp_path / "first.tsv"
    second = tmp_path / "second.tsv"
    metadata = tmp_path / "model.json"
    args = command(tree, counts) + ["--bootstrap", "2", "--seed", "815"]
    main(args + ["-o", str(first), "--model-out", str(metadata)])
    main(args + ["-o", str(second)])
    assert first.read_text() == second.read_text()
    table = pd.read_csv(first, sep="\t")
    assert table.p_value.iloc[0] == pytest.approx(1 / 3)
    assert table.p_value_method.iloc[0] == "plugin-parametric-bootstrap-search-maximum"
    assert table.count_support.iloc[0] == "background_compatible"
    audit = json.loads(metadata.read_text())["bootstrap_tail_audit"]
    assert audit["method"] == "pointwise-clopper-pearson"
    assert audit["minimum_plus_one_p"] == pytest.approx(1 / 3)
    (candidate,) = audit["candidates"]
    assert (
        candidate["num_exceedances"] == 0 and candidate["tail_probability_lower"] == 0
    )
    assert candidate["alpha_resolution"] == "overlaps"


def test_bootstrap_tail_audit_distinguishes_both_limits():
    from types import SimpleNamespace

    from nwkit.wgd_count import bootstrap_tail_audit

    result = SimpleNamespace(
        bootstrap_statistics=(0.0,) * 100,
        candidates=[
            SimpleNamespace(node=1, improvement=0.0),
            SimpleNamespace(node=2, improvement=1.0),
        ],
    )
    audit = bootstrap_tail_audit(result, SimpleNamespace(branch_ids=(0, 10, 20)), 0.05)
    first, last = audit["candidates"]
    assert first == {
        "branch_id": 10,
        "num_exceedances": 100,
        "tail_probability_lower": pytest.approx(0.025**0.01),
        "tail_probability_upper": 1.0,
        "alpha_resolution": "above",
    }
    assert last == {
        "branch_id": 20,
        "num_exceedances": 0,
        "tail_probability_lower": 0.0,
        "tail_probability_upper": pytest.approx(1 - 0.025**0.01),
        "alpha_resolution": "below",
    }


def test_real_bootstrap_checkpoint_resumes_holes_and_matches_outputs(
    tmp_path, monkeypatch
):
    import nwkit.wgd_count_fit as fitting
    from nwkit.wgd_count_checkpoint import _statistics_digest

    tree, counts = inputs(tmp_path)
    checkpoint = tmp_path / "bootstrap.json"
    first, second = tmp_path / "first.tsv", tmp_path / "second.tsv"
    first_model, second_model = tmp_path / "first.json", tmp_path / "second.json"
    args = command(tree, counts) + ["--bootstrap", "3", "--seed", "815"]
    options = ["--bootstrap-checkpoint", str(checkpoint)]
    main(args + options + ["-o", str(first), "--model-out", str(first_model)])
    previous = json.loads(checkpoint.read_text())
    expected = list(previous["statistics"])
    assert len(expected) == 3 and all(value is not None for value in expected)
    previous["statistics"][1] = None
    previous["statistics_sha256"] = _statistics_digest(previous["statistics"])
    checkpoint.write_text(json.dumps(previous))
    calls = []
    refit = fitting._bootstrap_scan_statistic

    def record(*args):
        calls.append(args[1])
        return refit(*args)

    monkeypatch.setattr(fitting, "_bootstrap_scan_statistic", record)
    main(args + options + ["-o", str(second), "--model-out", str(second_model)])
    assert len(calls) == 1
    assert json.loads(checkpoint.read_text())["statistics"] == expected
    assert first.read_bytes() == second.read_bytes()
    assert first_model.read_bytes() == second_model.read_bytes()
    corrupted = json.loads(checkpoint.read_text())
    corrupted["statistics"][0] += 1
    checkpoint.write_text(json.dumps(corrupted))
    with pytest.raises(ValueError, match="checksum does not match"):
        main(args + options + ["-o", str(second)])
    previous["statistics"][0] = -1
    checkpoint.write_text(json.dumps(previous))
    with pytest.raises(ValueError, match="invalid statistics"):
        main(args + options + ["-o", str(second)])
    with pytest.raises(ValueError, match="does not match"):
        main(command(tree, counts) + ["--bootstrap", "3", "--seed", "816", *options])


def test_bootstrap_checkpoint_has_distinct_file_inputs_and_outputs(tmp_path):
    tree, counts = inputs(tmp_path)
    with pytest.raises(ValueError, match="positive --bootstrap"):
        main(
            command(tree, counts) + ["--bootstrap-checkpoint", str(tmp_path / "b.json")]
        )
    with pytest.raises(ValueError, match="not STDIN"):
        main(
            command(tree, "-")
            + ["--bootstrap", "1", "--bootstrap-checkpoint", str(tmp_path / "b.json")]
        )
    for path in (tree, counts):
        with pytest.raises(ValueError):
            main(
                command(tree, counts)
                + ["--bootstrap", "1", "--bootstrap-checkpoint", str(path)]
            )


def test_tip_and_family_labels_remain_literal_and_missing_is_explicit(tmp_path):
    path = tmp_path / "counts.tsv"
    path.write_text(
        "\ufefffamily_id\tNA\t001\nf001\t1\tNA\nNA\t2\t0\n", encoding="utf-8"
    )
    families, counts = read_counts(str(path), ("001", "NA"), ",NA")
    assert families == ["f001", "NA"]
    assert pd.isna(counts[0, 0])
    assert counts[0, 1] == 1
    assert counts[1, 0] == 0


@pytest.mark.parametrize(
    "content",
    [
        "family_id\tA\tA\nf\t1\t1\n",
        "family_id\tA\tB\nf\t1\n",
        "family_id\tA\tB\nf\t1\t1\nf\t1\t1\n",
        "family_id\tA\tB\nf\t1\t-1\n",
        "family_id\tA\tB\nf\t1\t0.5\n",
        "family_id\tA\tB\nf\t1\tinf\n",
    ],
)
def test_invalid_counts_fail(content, tmp_path):
    path = tmp_path / "counts.tsv"
    path.write_text(content)
    with pytest.raises(ValueError):
        read_counts(str(path), ("A", "B"), ",NA")


def test_output_cannot_replace_declared_input_or_other_output(tmp_path):
    tree, counts = inputs(tmp_path)
    original = counts.read_text()
    with pytest.raises(ValueError):
        main(command(tree, counts) + ["-o", str(counts)])
    assert counts.read_text() == original
    out = tmp_path / "result.tsv"
    out.write_text("keep")
    with pytest.raises(ValueError):
        main(command(tree, counts) + ["-o", str(out), "--model-out", str(out)])
    assert out.read_text() == "keep"


def test_single_stdin_owner_and_stdout_primary(tmp_path, monkeypatch, capsys):
    import io

    tree, counts = inputs(tmp_path)
    monkeypatch.setattr("sys.stdin", io.StringIO(counts.read_text()))
    main(command(tree, "-") + ["-o", "-"])
    output = capsys.readouterr().out
    assert output.startswith("rank\tevent_id\tbranch_id")
    with pytest.raises(ValueError, match="STDIN"):
        main(["wgd-count", "-i", "-", "--counts", "-"])


def test_unknown_root_and_missing_lengths_are_not_invented(tmp_path):
    tree, counts = inputs(tmp_path)
    tree.write_text("(A:1,B:1,C:1);")
    with pytest.raises(ValueError, match="rooted"):
        main(command(tree, counts))
    tree.write_text("((A,B:0.5):0.5,C:1);")
    with pytest.raises(ValueError, match="branch length"):
        main(command(tree, counts))
