import json

import numpy as np
import pandas as pd
import pytest

from nwkit.cli import main
from nwkit.ksrate import make_trios
from nwkit.ksrate_model import (
    KsObservation,
    KsTrio,
    PairKs,
    bootstrap_corrections,
    corrected_ks,
    grouped_corrections,
    median_order_statistic_interval,
    simultaneous_median_corrections,
    trio_values,
)
from nwkit.util import read_tree


def observations():
    values = {
        ("A", "B"): 0.3,
        ("A", "C"): 1.5,
        ("A", "D"): 1.6,
        ("B", "C"): 1.6,
        ("B", "D"): 1.7,
        ("C", "D"): 0.7,
    }
    return [
        KsObservation(a, b, f"f{i}", value * multiplier)
        for (a, b), value in values.items()
        for i, multiplier in enumerate([0.8, 0.9, 1, 1.1, 1.2])
    ]


def inputs(tmp_path):
    tree = tmp_path / "species.nwk"
    tree.write_text("((A:1,B:1):1,(C:1,D:1):1);")
    records = observations()
    ks = tmp_path / "ks.tsv"
    pd.DataFrame([record.__dict__ for record in records]).to_csv(
        ks, sep="\t", index=False
    )
    return tree, ks


def test_formula_on_asymmetric_rates_and_negative_not_clipped():
    assert corrected_ks(0.3, 1.5, 1.6) == pytest.approx(0.2)
    assert corrected_ks(0.3, 1.6, 1.5) == pytest.approx(0.4)
    assert corrected_ks(0.1, 0.1, 1) == pytest.approx(-0.8)


def test_finite_ks_cancellation_does_not_overflow_intermediate_sum():
    assert corrected_ks(1e308, 1e308, 1e308) == 1e308
    assert corrected_ks(5e-324, 0, 0) == 5e-324
    with pytest.raises(ValueError, match="finite"):
        corrected_ks(1e308, 1e308, 0)


@pytest.mark.parametrize("value", [1e308, 5e-324])
def test_pair_and_group_medians_are_stable_at_float_extremes(value):
    pairs = PairKs([KsObservation("A", "B", f"f{i}", value) for i in range(2)], "AB")
    assert pairs.estimates()[("A", "B")] == value
    trios = [KsTrio("A", "B", "C", "AB", 1)] * 2
    assert grouped_corrections(trios, np.array([value, value]))[("A", "AB")] == value
    assert grouped_corrections(trios, np.array([-1e308, 1e308]))[("A", "AB")] == 0


def test_simultaneous_endpoints_avoid_cancellable_intermediate_overflow():
    pairs = PairKs(
        [
            KsObservation(a, b, f"f{i}", 1e308)
            for a, b in [("A", "B"), ("A", "C"), ("B", "C")]
            for i in range(8)
        ],
        "ABC",
    )
    trios = [KsTrio("A", "B", "C", "AB", 1)]
    intervals, _ = simultaneous_median_corrections(pairs, trios, 0.95)
    assert intervals[("A", "AB")] == (1e308, 1e308)


def test_weighted_median_and_joint_family_resampling():
    pairs = PairKs(observations(), "ABCD")
    estimate = pairs.estimates(np.array([2, 1, 0, 1, 0]))
    assert estimate[("A", "B")] == pytest.approx(0.3 * 0.85)
    trios = [KsTrio("A", "B", "C", "AB", 1), KsTrio("B", "A", "C", "AB", 1)]
    one = bootstrap_corrections(pairs, trios, 99, 198)
    two = bootstrap_corrections(pairs, trios, 99, 198)
    for key in one:
        np.testing.assert_array_equal(one[key], two[key])
    np.testing.assert_allclose(one[("B", "AB")], one[("A", "AB")] * 2, atol=1e-14)


def test_cli_known_corrections_multiple_outgroups_and_unresolved_root(tmp_path):
    tree, ks = inputs(tmp_path)
    out, trios, model = (
        tmp_path / name for name in ("nodes.tsv", "trios.tsv", "model.json")
    )
    main(
        [
            "ksrate",
            "-i",
            str(tree),
            "--ks-tsv",
            str(ks),
            "-o",
            str(out),
            "--trios-out",
            str(trios),
            "--model-out",
            str(model),
            "--bootstrap",
            "99",
            "--seed",
            "198",
        ]
    )
    nodes = pd.read_csv(out, sep="\t")
    for focal, expected in zip("ABCD", [0.2, 0.4, 0.6, 0.8], strict=True):
        rows = nodes[nodes.focal == focal]
        assert rows.iloc[0].corrected_ks == pytest.approx(expected)
        assert rows.iloc[0].num_complete_trios == 2
        assert rows.iloc[0].interval_status == "ok"
        assert rows.iloc[0].ci_lower <= expected <= rows.iloc[0].ci_upper
        assert rows.iloc[1].status == "no_external_outgroup"
        assert pd.isna(rows.iloc[1].corrected_ks)
    metadata = json.loads(model.read_text())
    assert metadata["formula"] == "Ks(F,S)+Ks(F,O)-Ks(S,O)"
    assert "not posterior" in metadata["uncertainty"]
    assert len(pd.read_csv(trios, sep="\t")) == 8


def test_nonmonotone_corrections_are_reported_not_projected(tmp_path):
    tree = tmp_path / "tree.nwk"
    tree.write_text("(((A:1,B:1):1,C:2):1,D:3);")
    values = {
        ("A", "B"): 2,
        ("A", "C"): 2,
        ("B", "C"): 1,
        ("A", "D"): 1,
        ("B", "D"): 1,
        ("C", "D"): 2,
    }
    ks = tmp_path / "ks.tsv"
    pd.DataFrame(
        [
            {"species_a": a, "species_b": b, "family_id": "f", "ks": value}
            for (a, b), value in values.items()
        ]
    ).to_csv(ks, sep="\t", index=False)
    out = tmp_path / "out.tsv"
    main(
        [
            "ksrate",
            "-i",
            str(tree),
            "--ks-tsv",
            str(ks),
            "--focals",
            "A",
            "--bootstrap",
            "0",
            "-o",
            str(out),
        ]
    )
    nodes = pd.read_csv(out, sep="\t")
    assert list(nodes.raw_corrected_ks[:2]) == [3, 1]
    assert nodes.monotone_from_younger_node.iloc[1] == "no"


def test_missing_bootstrap_comparisons_make_interval_unavailable(tmp_path):
    tree, ks = inputs(tmp_path)
    table = pd.read_csv(ks, sep="\t")
    table["family_id"] = [f"unique{i}" for i in range(len(table))]
    table.iloc[::5].to_csv(ks, sep="\t", index=False)
    out = tmp_path / "out.tsv"
    main(
        [
            "ksrate",
            "-i",
            str(tree),
            "--ks-tsv",
            str(ks),
            "--bootstrap",
            "99",
            "-o",
            str(out),
        ]
    )
    nodes = pd.read_csv(out, sep="\t")
    assert nodes.interval_status.iloc[0] == "unavailable_missing_bootstrap_comparisons"
    assert nodes.num_bootstrap_estimable.iloc[0] < 99
    assert pd.isna(nodes.ci_lower.iloc[0])


def test_partial_pairs_and_exact_keys(tmp_path):
    pairs = PairKs([KsObservation("001", "NA", "NA", 0.3)], ["001", "NA", "C"])
    trios = [KsTrio("001", "NA", "C", "id", 1)]
    assert np.isnan(trio_values(trios, pairs.estimates())[0])
    tree = read_tree("((A:1,B:1):1,(C:1,D:1):1);", 1, False)
    trios, events = make_trios(tree, ["A"])
    assert {(trio.sister, trio.outgroup) for trio in trios} == {("B", "C"), ("B", "D")}
    assert events[-1]["trio_indices"] == []


def test_duplicate_family_unknown_species_and_bad_distance_rejected():
    record = KsObservation("A", "B", "f", 1)
    for records in (
        [record, record],
        [KsObservation("A", "C", "f", 1)],
        [KsObservation("A", "B", "f", -1)],
        [KsObservation("A", "A", "f", 1)],
    ):
        with pytest.raises(ValueError):
            PairKs(records, "AB")


@pytest.mark.parametrize("weight", [np.inf, np.nan, -np.inf, 1e20])
def test_nonfinite_family_weights_are_rejected(weight):
    pairs = PairKs([KsObservation("A", "B", "f", 1.0)], "AB")
    with pytest.raises(ValueError, match="weights"):
        pairs.estimates([weight])


def test_family_weight_sum_cannot_overflow_integer_accumulation():
    pairs = PairKs(
        [KsObservation("A", "B", "f", 1.0), KsObservation("A", "B", "g", 2.0)],
        "AB",
    )
    with pytest.raises(ValueError, match="accumulation"):
        pairs.estimates(np.array([2**62, 2**62], dtype=np.int64))


def test_median_order_statistic_interval_has_exact_ranks_and_no_short_sample_fallback():
    assert median_order_statistic_interval(np.arange(80), 0.05) == (30, 49)
    low, high = median_order_statistic_interval(np.arange(5), 0.05)
    assert np.isnan(low) and np.isnan(high)
    for error in (0, 1, np.nan):
        with pytest.raises(ValueError):
            median_order_statistic_interval([1, 2], error)


def test_simultaneous_bounds_do_not_drop_unbounded_complete_trios():
    pairs = PairKs(observations(), "ABCD")
    trios = [KsTrio("A", "B", "C", "AB", 1), KsTrio("A", "B", "D", "AB", 1)]
    intervals, comparisons = simultaneous_median_corrections(pairs, trios, 0.95)
    assert comparisons == 5
    assert all(np.isnan(value) for value in intervals[("A", "AB")])


def test_conservative_cli_interval_does_not_require_bootstrap_draws(tmp_path):
    tree, ks = inputs(tmp_path)
    records = [
        KsObservation(record.species_a, record.species_b, f"f{i:03d}", record.ks)
        for i, record in enumerate(observations() * 4)
    ]
    pd.DataFrame([record.__dict__ for record in records]).to_csv(
        ks, sep="\t", index=False
    )
    out, report = tmp_path / "out.tsv", tmp_path / "model.json"
    main(
        [
            "ksrate",
            "-i",
            str(tree),
            "--ks-tsv",
            str(ks),
            "--focals",
            "A",
            "--ci-method",
            "pair-median-bonferroni",
            "--bootstrap",
            "0",
            "-o",
            str(out),
            "--model-out",
            str(report),
        ]
    )
    nodes = pd.read_csv(out, sep="\t")
    first = nodes.iloc[0]
    assert first.ci_method == "pair-median-bonferroni"
    assert first.interval_status == "ok"
    assert first.bootstrap_interval_status == "not_run"
    assert pd.isna(first.bootstrap_ci_lower)
    assert first.ci_lower <= 0.2 <= first.ci_upper
    assert first.num_simultaneous_pairs == 5
    assert "iid independent families" in json.loads(report.read_text())["uncertainty"]


def test_inputs_and_related_outputs_are_protected(tmp_path):
    tree, ks = inputs(tmp_path)
    original = ks.read_text()
    with pytest.raises(ValueError):
        main(["ksrate", "-i", str(tree), "--ks-tsv", str(ks), "-o", str(ks)])
    assert ks.read_text() == original
    with pytest.raises(ValueError, match="STDIN"):
        main(["ksrate", "-i", "-", "--ks-tsv", "-"])
