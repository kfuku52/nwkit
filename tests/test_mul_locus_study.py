import importlib.util
import json
import os
import sys
from argparse import Namespace
from concurrent.futures.process import BrokenProcessPool
from pathlib import Path

import pandas as pd
import pytest

from nwkit.mul_locus_mc import make_tasks, validate_model
from tests.test_mul_reconcile import tree


def abrupt_study_worker(*args):
    os._exit(17)


def empty_worker_init(*args):
    pass


@pytest.fixture
def study(monkeypatch):
    pytest.importorskip("msprime", reason="optional independent research simulator")
    pytest.importorskip("Bio", reason="optional independent research tree estimator")
    directory = Path(__file__).resolve().parents[1] / "examples" / "mul-locus"
    monkeypatch.syspath_prepend(str(directory))
    spec = importlib.util.spec_from_file_location(
        "nwkit_locus_grid_study", directory / "calibration.py"
    )
    module = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, spec.name, module)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def study_args(tmp_path):
    return Namespace(
        output=tmp_path,
        cpus=1,
        replicates=2,
        families=30,
        samples=10,
        bootstrap=3,
        seed=20261033,
        bank_seed=20261030,
        calibration_seed=20261034,
        event_alpha=0.05,
        scenarios=["baseline", "turnover", "missing-ils"],
    )


def test_protocol_contains_all_on_grid_nulls_and_off_grid_diagnostics(study):
    for i, name in enumerate(study.SCENARIOS):
        model = study.scenario_model(name, 100000, 20261030 + i)
        tasks, _ = make_tasks(
            tree(study.SPECIES),
            "X",
            "A B",
            validate_model(model, tree(study.SPECIES)),
            100,
        )
        assert len(tasks) == 20
        null = {(t[4].duplication, t[4].loss, t[4].ne) for t in tasks if t[0] == 0}
        cases = study.scenario_cases(name)
        assert len(cases) == 7
        assert {
            (p.duplication, p.loss, p.ne)
            for _, truth, p, g in cases
            if truth == 0 and g
        } == null
        assert all(
            (p.duplication, p.loss, p.ne) not in null
            for _, truth, p, g in cases
            if truth == 0 and not g
        )
        assert {truth for _, truth, _, _ in cases} == {0, 1, 2}
        assert model["integration"] == "ancestral-stratified"


def test_failure_aware_rates_keep_every_planned_analysis(study):
    rows = [
        {
            "scenario": "baseline",
            "case": "null",
            "true_candidate": 0,
            "null_on_grid": True,
            "status": "completed",
            "reported_event": True,
            "point_event": True,
            "fitted_null_event_same_draws": True,
            "parent_correct": False,
        },
        {
            "scenario": "baseline",
            "case": "null",
            "true_candidate": 0,
            "null_on_grid": True,
            "status": "calibration-failed",
        },
    ]
    result = study.summarize(rows)[0]
    assert result["planned"] == 2 and result["completed"] == 1 and result["failed"] == 1
    assert result["failure_aware_rate_lower"] == 0.5
    assert result["failure_aware_rate_upper"] == 1
    assert result["failure_aware_95_lower"] < 0.5
    assert result["failure_aware_95_upper"] == 1


def test_bank_failure_is_retained_for_all_planned_datasets(
    study, study_args, tmp_path, monkeypatch
):
    monkeypatch.setattr(
        study,
        "build_bank",
        lambda *a: (_ for _ in ()).throw(ValueError("injected bank cap")),
    )
    rows = study.run_scenario(
        study_args, "baseline", 0, study.scenario_model("baseline", 10, 3)
    )
    assert len(rows) == 14
    assert all(r["status"] == "bank-failed" for r in rows)
    saved = json.loads((tmp_path / "baseline" / "bank-failures.json").read_text())
    assert len(saved) == 20
    assert all("injected bank cap" in r["error"] for r in saved)
    assert sum(r["failed"] for r in study.summarize(rows)) == 14
    assert len(pd.read_csv(tmp_path / "baseline" / "summary.tsv", sep="\t")) == 14
    assert len({r["calibration_seed"] for r in rows}) == 14
    assert all(r["families"] == 30 and r["true_parameters"] for r in rows)


def test_scenario_subset_retains_full_protocol_streams(study, study_args, tmp_path):
    full = study.write_protocol(study_args)
    study_args.output = tmp_path / "subset"
    study_args.output.mkdir()
    study_args.scenarios = ["missing-ils"]
    subset = study.write_protocol(study_args)
    assert subset["models"]["missing-ils"] == full["models"]["missing-ils"]


def test_main_uses_fixed_scenario_indices(study, monkeypatch, tmp_path):
    calls = []

    def run(args, name, index, config):
        calls.append((name, index))
        return [
            {
                "scenario": name,
                "case": "null",
                "true_candidate": 0,
                "null_on_grid": True,
                "status": "bank-failed",
            }
        ]

    monkeypatch.setattr(study, "run_scenario", run)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "calibration.py",
            "--output",
            str(tmp_path / "run"),
            "--scenarios",
            "missing-ils",
            "baseline",
        ],
    )
    assert study.main() == 1
    assert calls == [("missing-ils", 2), ("baseline", 0)]


def test_calibration_seed_collision_fails_before_protocol(
    study, study_args, monkeypatch
):
    class RepeatedSeed:
        def __init__(self, *a):
            pass

        def generate_state(self, *a):
            return [7]

    monkeypatch.setattr(study.np.random, "SeedSequence", RepeatedSeed)
    with pytest.raises(ValueError, match="Calibration seed collision.*baseline"):
        study.write_protocol(study_args)
    assert not (study_args.output / "protocol.json").exists()


def test_real_cross_case_seed_collision_is_rejected(study, study_args, monkeypatch):
    cases = study.scenario_cases("baseline")
    jobs = [
        ("baseline", *cases[i][:3], cases[i][3], 0, i, r, study_args)
        for i, r in ((3, 11575), (5, 13904))
    ]
    assert [study.trial_record(j)["calibration_seed"] for j in jobs] == [3590509910] * 2
    study_args.scenarios = ["baseline"]
    monkeypatch.setattr(study, "trial_jobs", lambda *a: jobs)
    with pytest.raises(ValueError, match="Calibration seed collision"):
        study.validate_calibration_seeds(study_args)


@pytest.mark.parametrize("phase", ["bank", "evaluation", "saved-evaluation"])
def test_worker_crash_preserves_planned_trials_and_partial_results(
    study, study_args, monkeypatch, phase
):
    class FailedPool:
        def __init__(self, **kw):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return None

        def map(self, function, jobs, *a, **kw):
            if function is study.integration_bank:
                if phase != "bank":
                    return iter([object() for _ in jobs])
                return self.crash(object())
            job = jobs[0]
            row = {
                **study.trial_record(job),
                "status": "completed",
                "reported_event": False,
                "point_event": False,
                "fitted_null_event_same_draws": False,
                "parent_correct": False,
            }
            if phase == "saved-evaluation":
                saved = {
                    **study.trial_record(jobs[1]),
                    **{
                        k: v for k, v in row.items() if k not in study.trial_record(job)
                    },
                }
                directory = study_args.output / "baseline" / "null-d0-ne0-r2"
                directory.mkdir()
                (directory / "summary.json").write_text(json.dumps(saved))
            return self.crash(row)

        def crash(self, first):
            yield first
            raise BrokenProcessPool("injected native worker crash")

    study_args.cpus = 2
    monkeypatch.setattr(study, "ProcessPoolExecutor", FailedPool)
    monkeypatch.setattr(study, "bank_record", lambda b: {"samples": 1})
    rows = study.run_scenario(
        study_args, "baseline", 0, study.scenario_model("baseline", 10, 3)
    )
    assert len(rows) == 14
    directory = study_args.output / "baseline"
    if phase == "bank":
        assert all(r["status"] == "bank-failed" for r in rows)
        assert len(json.loads((directory / "partial-banks.json").read_text())) == 1
        assert len(json.loads((directory / "bank-failures.json").read_text())) == 19
    else:
        assert rows[0]["status"] == "completed"
        completed = 2 if phase == "saved-evaluation" else 1
        assert all(r["status"] == "completed" for r in rows[:completed])
        assert all(r["status"] == "worker-failed" for r in rows[completed:])
        assert len(list(directory.glob("*/failure.json"))) == 14 - completed
    assert all(
        "injected native worker crash" in r["error"]
        for r in rows
        if r["status"] != "completed"
    )
    assert sum(r["planned"] for r in study.summarize(rows)) == 14


@pytest.mark.parametrize(
    "saved",
    [
        "{",
        {"status": "completed"},
        {"status": "generation-failed", "error": "saved generation error"},
    ],
)
def test_worker_recovery_preserves_unreadable_and_failed_records(
    study, study_args, saved
):
    job = study.trial_jobs(study_args, "baseline", 0)[0]
    directory = study_args.output / "baseline" / "null-d0-ne0-r1"
    directory.mkdir(parents=True)
    record = (
        saved
        if isinstance(saved, str)
        else json.dumps({**study.trial_record(job), **saved})
    )
    (directory / "failure.json").write_text(record)
    row = study.failed_trial(job, "worker-failed", "native worker error")
    assert (directory / "failure.json").read_text() == record
    if isinstance(saved, dict) and saved["status"] == "generation-failed":
        assert (
            row["status"] == "generation-failed"
            and row["error"] == "saved generation error"
        )
    else:
        assert row["status"] == "worker-failed"
        assert (directory / "worker-failure.json").exists()


@pytest.mark.slow
@pytest.mark.parametrize("phase", ["bank", "evaluation"])
def test_real_process_death_retains_all_planned_trials(
    study, study_args, monkeypatch, phase
):
    study_args.cpus = 2
    target = "integration_bank" if phase == "bank" else "evaluate"
    monkeypatch.setattr(study, target, abrupt_study_worker)
    monkeypatch.setattr(study, "worker_init", empty_worker_init)
    if phase == "evaluation":
        monkeypatch.setattr(study, "integration_bank", native_integration_bank)
    rows = study.run_scenario(
        study_args, "baseline", 0, study.scenario_model("baseline", 2, 3)
    )
    status = "bank-failed" if phase == "bank" else "worker-failed"
    assert len(rows) == 14 and all(row["status"] == status for row in rows)
    assert all("BrokenProcessPool" in row["error"] for row in rows)
    assert sum(row["planned"] for row in study.summarize(rows)) == 14
    assert len(list((study_args.output / "baseline").glob("*/failure.json"))) == 14


def native_integration_bank(task, config):
    from nwkit.mul_locus_mc import build_bank

    return build_bank(task, config)
