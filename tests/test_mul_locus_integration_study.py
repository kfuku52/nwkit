import importlib.util
import json
import sys
from argparse import Namespace
from collections import Counter
from pathlib import Path

import pytest

from nwkit.mul_locus import LocusParameters
from nwkit.mul_locus_mc import (
    LocusBank,
    build_bank,
    make_tasks,
    pattern_probability,
    validate_model,
)


@pytest.fixture
def probe(monkeypatch):
    pytest.importorskip("msprime", reason="optional independent research simulator")
    pytest.importorskip("Bio", reason="optional independent research tree estimator")
    directory = Path(__file__).resolve().parents[1] / "examples" / "mul-locus"
    monkeypatch.syspath_prepend(str(directory))
    spec = importlib.util.spec_from_file_location(
        "nwkit_locus_integration_probe", directory / "integration.py"
    )
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, spec.name, module)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def args(tmp_path):
    return Namespace(
        output=tmp_path,
        samples=20,
        bank_seed=20261102,
        seed=20261101,
        calibration_seed=20261103,
        exact_tip_limit=4,
        max_states=100000,
        families=2,
        bootstrap=3,
        replicates=2,
    )


def test_copy_vector_distinguishes_topology_and_post_detection_copy_support(probe):
    a, b = ("tip", "A"), ("tip", "B")
    first = ("node", ("node", a, a), b)
    second = ("node", ("node", a, b), a)
    absent = ("node", a, a)
    bank = LocusBank(
        0, "NA", 0, None, LocusParameters(0, 0, 1, 0.5), Counter({first: 5}), 5, 5
    )
    rows = probe.support_rows([bank], [first, second, absent])
    assert [r["absence"] for r in rows] == [
        "none",
        "topology",
        "copy-count-or-detection",
    ]
    assert [r["copy_vector_hits"] for r in rows] == [5, 5, 0]


@pytest.mark.parametrize("integration", ["selected-histogram", "ancestral-stratified"])
def test_saved_histogram_banks_restore_strata_and_exact_probabilities(
    probe, tmp_path, integration
):
    config = probe.scenario_model("missing-ils", 40, 20261112) | {
        "integration": integration
    }
    species = probe.tree(probe.SPECIES)
    task = make_tasks(species, "X", "A B", validate_model(config, species), 100)[0][0]
    original = build_bank(task, config)
    path = tmp_path / "banks.json"
    path.write_text(json.dumps([probe.bank_record(original)]))
    restored = probe.read_banks(path, config)[0]
    assert restored.counts == original.counts
    assert restored.strata == original.strata
    assert (
        restored.samples == original.samples and restored.attempts == original.attempts
    )
    for signature in original.counts:
        assert pattern_probability(restored, signature, 0.001) == pattern_probability(
            original, signature, 0.001
        )


def test_seed_namespaces_canonical_subset_independent_and_shared_across_methods(
    probe, args
):
    seen = set()
    for name in reversed(probe.SCENARIOS):
        for i, case in enumerate(probe.scenario_cases(name)):
            for r in range(args.replicates):
                rows = probe.case_records(args, name, i, r, case)
                assert len({row["calibration_seed"] for row in rows}) == 1
                seed = rows[0]["calibration_seed"]
                assert seed not in seen
                seen.add(seed)
                assert rows[0]["data_seed_namespace"][2] == tuple(
                    probe.SCENARIOS
                ).index(name)
                assert int(seed) >= 2**32


def test_bank_failure_keeps_every_planned_trial_with_seed_and_truth(
    probe, args, monkeypatch
):
    def fail(*a, **kw):
        raise ValueError("work cap; no history discarded")

    monkeypatch.setattr(probe, "build_paired_banks", fail)
    rows = probe.run_probe(args, "baseline")
    assert len(rows) == 7 * 2 * 3
    assert all(row["status"] == "bank-failed" for row in rows)
    assert all("calibration_seed" in row and "true_parameters" in row for row in rows)
    assert json.loads((args.output / "baseline/failed-trials.json").read_text()) == rows


@pytest.mark.parametrize("mode", ["exception", "unexpected", "wrong-count"])
def test_independent_generation_failure_keeps_all_methods(
    probe, args, monkeypatch, mode
):
    (args.output / "baseline").mkdir()

    def fail(*a, **kw):
        if mode == "exception":
            raise ArithmeticError("reference failure")
        if mode == "unexpected":
            raise UnboundLocalError("reference estimator failed unexpectedly")
        return [], 0

    monkeypatch.setattr(probe, "reference_selected", fail)
    case = probe.scenario_cases("baseline")[0]
    rows = probe.evaluate_case(
        args,
        "baseline",
        0,
        0,
        case,
        {},
        probe.scenario_model("baseline", 20, args.bank_seed),
    )
    assert len(rows) == 3
    assert all(row["status"] == "generation-failed" for row in rows)
    assert all("calibration_seed" in row for row in rows)
    for row in rows:
        saved = json.loads(
            (
                args.output
                / "baseline/null-d0-ne0-r1"
                / f"{row['method']}-summary.json"
            ).read_text()
        )
        assert saved == row


def test_calibration_failure_preserves_methods_and_shared_generation(
    probe, args, monkeypatch
):
    (args.output / "baseline").mkdir()
    observations = [("node", ("tip", "A"), ("tip", "B"))] * args.families
    calls = []
    monkeypatch.setattr(probe, "reference_selected", lambda *a, **kw: (observations, 3))
    monkeypatch.setattr(probe, "integration_tables", lambda *a: [])
    monkeypatch.setattr(probe, "integration_alpha", lambda *a: 0.001)

    def fail(*a, **kw):
        calls.append((a[2]["seed"], kw["scorer"], kw["sampler"]))
        raise ValueError("support failure")

    monkeypatch.setattr(probe, "calibrate", fail)
    case = probe.scenario_cases("baseline")[0]
    rows = probe.evaluate_case(
        args,
        "baseline",
        0,
        0,
        case,
        {method: [] for method in probe.METHODS},
        probe.scenario_model("baseline", 20, args.bank_seed),
    )
    assert len(rows) == 3
    assert all(row["status"] == "calibration-failed" for row in rows)
    assert len({c[0] for c in calls}) == 1
    assert calls[1][1] == calls[2][1] == probe.score_integrated_bank
    assert all(c[2] is probe.reference_selected for c in calls)


def test_true_only_reference_avoids_sequence_and_estimated_tree_pipeline(
    probe, monkeypatch
):
    import numpy as np
    import pilot

    def forbidden(*a, **kw):
        pytest.fail("true-only sampling invoked the sequence/NJ pipeline")

    monkeypatch.setattr(pilot, "estimate_tree", forbidden)
    monkeypatch.setattr(pilot.msprime, "sim_mutations", forbidden)
    config = probe.scenario_model("baseline", 20, 123)
    point = probe.scenario_cases("baseline")[0][2]
    population = probe.tree(probe.SPECIES)
    values, attempts = pilot.reference_true_selected(
        population, point, config, np.random.default_rng(20261110), count=3
    )
    assert len(values) == 3
    assert attempts >= 3
    with pytest.raises(ValueError, match="True-only"):
        pilot.reference_family(
            population,
            point,
            config,
            np.random.default_rng(1),
            estimated=True,
            true_only=True,
        )
