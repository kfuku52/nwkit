"""Independent evidence checks reject numerically plausible corruptions."""

import copy
import gzip
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from compare_shift_null_calibration import summarize  # noqa: E402
from shift_calibration_audit import check_probability_metadata  # noqa: E402
from validate_shift_null_contract import tree_text  # noqa: E402
from verify_shift_null_contract import (  # noqa: E402
    audit,
    check_stages,
    check_winner,
    same_generated_array,
    same_replayed_model,
    same_replayed_summary,
    same_replayed_tests,
)

from nwkit.shift_calibration import CalibratedSearch  # noqa: E402
from nwkit.shift_candidates import enumerate_candidates  # noqa: E402
from nwkit.util import read_tree  # noqa: E402


@pytest.fixture(scope="module")
def evidence():
    tree = read_tree(tree_text(4, "pectinate"), "auto", True, quiet=True)
    search = CalibratedSearch(tree)
    y = np.array([0.1, -0.3, 1.2, 0.7])
    return tree, search.models, y, search.fit(y, seed=183, replicates=19)


def test_valid_winner_and_stage_contract(evidence):
    tree, models, y, fit = evidence
    assert check_winner(tree, y, np.zeros(4), fit, models) < 1e-10
    assert check_stages(fit, False, 19, 0.05) == bool(fit["model"]["shift_branch_ids"])


def test_nuisance_grid_accepts_only_machine_roundoff():
    search = SimpleNamespace(known_error=False, grid=[0.1, np.inf])
    test = dict(
        p_value=1.0,
        p_value_lower_bound=1.0,
        p_value_kind="grid_supremum",
        null_alpha_evaluations=[
            dict(alpha_height=np.nextafter(0.1, np.inf), p_value=1.0),
            dict(alpha_height=None, p_value=1.0),
        ],
    )
    check_probability_metadata(test, 0, search, 19, 0.05)
    for value in [0.100001, None, np.nan]:
        changed = copy.deepcopy(test)
        changed["null_alpha_evaluations"][0]["alpha_height"] = value
        with pytest.raises(ValueError, match="evaluation grid"):
            check_probability_metadata(changed, 0, search, 19, 0.05)


def test_winner_grid_membership_uses_the_metadata_roundoff_contract():
    with gzip.open(
        "examples/shift/null-contract-timing-pilot/records.jsonl.gz", "rt"
    ) as stream:
        row = next(
            row
            for row in map(json.loads, stream)
            if row["lanes"][0]["fit"]["alpha_height"] not in (None, 0)
        )
    specification = json.loads(
        Path("examples/shift/null-contract-timing-pilot/protocol.json").read_text()
    )
    cell = next(
        cell for cell in specification["cells"] if cell["cell_id"] == row["cell_id"]
    )
    tree = read_tree(cell["tree"], "auto", True, quiet=True)
    lane = row["lanes"][0]
    models, _ = enumerate_candidates(tree, convergence=lane["convergence"])
    fit = copy.deepcopy(lane["fit"])
    fit["alpha_height"] = np.nextafter(fit["alpha_height"], np.inf)
    assert (
        check_winner(
            tree,
            np.asarray(row["observations"]),
            np.asarray(row["known_variances"]),
            fit,
            models,
        )
        < 1e-10
    )
    fit["alpha_height"] += 1e-6
    with pytest.raises(ValueError, match="outside the declared grid"):
        check_winner(
            tree,
            np.asarray(row["observations"]),
            np.asarray(row["known_variances"]),
            fit,
            models,
        )


@pytest.mark.parametrize(
    "field", ["predicted", "contrast_log_likelihood", "candidate_count"]
)
def test_independent_winner_detects_changed_evidence(evidence, field):
    tree, models, y, fit = evidence
    changed = copy.deepcopy(fit)
    if field == "predicted":
        changed[field][0] += 0.01
    else:
        changed[field] += 1
    with pytest.raises(ValueError):
        check_winner(tree, y, np.zeros(4), changed, models)


def test_stage_rejects_invalid_monte_carlo_resolution(evidence):
    fit = copy.deepcopy(evidence[-1])
    fit["tests"][0]["p_value"] = 0.123
    with pytest.raises(ValueError, match="Monte Carlo"):
        check_stages(fit, False, 19, 0.05)


def test_paired_power_failure_remains_a_loss():
    case = dict(family="primary", scenario="single", root_model="OUfixedRoot")
    fit = dict(
        status="completed", any_shift=True, partition_recovered=True, mean_rmse=0.2
    )
    rows = [dict(case=case, plugin=fit, envelope=fit)] * 99
    rows.append(dict(case=case, plugin=dict(status="failed"), envelope=fit))
    result = summarize(rows)[0]
    assert result["attempted"] == 100
    assert result["failed"] == 1
    assert result["power_loss_upper_one_sided_95"] > 0.04


def test_replay_roundoff_does_not_relax_discrete_decisions(evidence):
    saved = evidence[-1]["tests"]
    replayed = copy.deepcopy(saved)
    replayed[0]["statistic"] = np.nextafter(saved[0]["statistic"], np.inf)
    assert same_replayed_tests(replayed, saved)
    for field in ("stage", "candidate_count", "p_value", "statistic"):
        changed = copy.deepcopy(replayed)
        changed[0][field] += 1e-6
        assert not same_replayed_tests(changed, saved)
    changed = copy.deepcopy(replayed)
    changed[0]["null_alpha_evaluations"][0]["p_value"] += 1e-6
    assert not same_replayed_tests(changed, saved)


def test_replay_grid_roundoff_keeps_probabilities_and_grid_identity_exact():
    # This archived dataset exposed a one-ulp geomspace difference in
    # Linux CI. All its statistics and probabilities otherwise matched exactly.
    with gzip.open(
        "examples/shift/null-contract-timing-pilot/records.jsonl.gz", "rt"
    ) as stream:
        row = next(
            row
            for row in map(json.loads, stream)
            if row["cell_id"] == 0 and row["replicate"] == 4
        )
    saved = row["lanes"][0]["fit"]["tests"]
    replayed = copy.deepcopy(saved)
    grid = replayed[0]["null_alpha_evaluations"]
    for entry in grid:
        alpha = entry["alpha_height"]
        if alpha is not None and alpha > 0:
            entry["alpha_height"] = np.nextafter(alpha, np.inf)
    assert same_replayed_tests(replayed, saved)
    assert replayed != saved

    for field in ("p_value", "p_value_lower_bound"):
        changed = copy.deepcopy(replayed)
        changed[0][field] = np.nextafter(changed[0][field], np.inf)
        assert not same_replayed_tests(changed, saved)
    changed = copy.deepcopy(replayed)
    changed[0]["null_alpha_evaluations"][0]["p_value"] = np.nextafter(
        grid[0]["p_value"], np.inf
    )
    assert not same_replayed_tests(changed, saved)
    for value in (0.100001, None, np.nan):
        changed = copy.deepcopy(replayed)
        changed[0]["null_alpha_evaluations"][10]["alpha_height"] = value
        assert not same_replayed_tests(changed, saved)
    for alter in (lambda entries: entries.pop(), lambda entries: entries.reverse()):
        changed = copy.deepcopy(replayed)
        alter(changed[0]["null_alpha_evaluations"])
        assert not same_replayed_tests(changed, saved)


def test_replay_accepts_only_same_tip_partition_representation(evidence):
    tree, _, _, _ = evidence
    models, _ = enumerate_candidates(tree, convergence=False)
    saved = models[8]
    equivalent = models[7]
    changed_partition = models[9]
    assert same_replayed_model(tree, equivalent, saved)
    assert not same_replayed_model(tree, changed_partition, saved)


def test_summary_roundoff_cannot_hide_changed_counts_or_acceptance():
    saved = json.loads(
        Path("examples/shift/null-contract-timing-pilot/summary.json").read_text()
    )
    replayed = copy.deepcopy(saved)
    field = "simultaneous_one_sided_95_upper"
    replayed["cell_results"][0][field] = np.nextafter(
        saved["cell_results"][0][field], np.inf
    )
    assert same_replayed_summary(replayed, saved)
    for key in (field, "attempted", "any_shift", "acceptance_upper"):
        changed = copy.deepcopy(replayed)
        changed["cell_results"][0][key] += 1e-6
        assert not same_replayed_summary(changed, saved)


def test_generated_array_roundoff_does_not_allow_broadcasting_or_changed_data():
    observed = np.array([0.1, 0.2])
    assert same_generated_array(observed, np.nextafter(observed, np.inf))
    assert not same_generated_array(observed, observed + 1e-6)
    assert not same_generated_array(observed, observed[None, :])
    assert not same_generated_array(observed, [np.nan, 0.2])


def test_all_failed_pairs_have_missing_rmse_and_worst_case_bound():
    row = dict(
        case=dict(family="primary", scenario="single", root_model="OUfixedRoot"),
        plugin=dict(status="failed"),
        envelope=dict(status="completed"),
    )
    result = summarize([row])[0]
    assert result["failed"] == 1
    assert result["power_loss_upper_one_sided_95"] == 1
    assert result["plugin_mean_rmse"] is None
    assert result["envelope_mean_rmse"] is None
    json.dumps(result, allow_nan=False)


def test_archived_engine_requires_explicit_scope_and_intact_snapshot(tmp_path):
    source = Path("examples/shift/null-contract-timing-pilot")
    bundle = tmp_path / "archived"
    shutil.copytree(source, bundle)
    with pytest.raises(ValueError, match="Active generator/fitting hash mismatch"):
        audit(bundle)
    before = {
        str(p.relative_to(bundle)): p.read_bytes()
        for p in bundle.rglob("*")
        if p.is_file()
    }
    # The saved experiment declares one BLAS/OpenMP thread. Apply that protocol
    # before importing NumPy in a child; changing os.environ after import would
    # leave the parent's already initialized BLAS pools unchanged.
    threads = json.loads((bundle / "protocol.json").read_text())["thread_environment"]
    replay = subprocess.run(
        [
            sys.executable,
            str(
                Path(__file__).resolve().parents[1]
                / "tools/verify_shift_null_contract.py"
            ),
            str(bundle),
            "--replay-stride",
            "1",
            "--frozen-engine",
        ],
        env={**os.environ, **threads},
        capture_output=True,
        text=True,
    )
    assert replay.returncode == 0, replay.stdout + replay.stderr
    result = json.loads(replay.stdout)
    assert {
        str(p.relative_to(bundle)): p.read_bytes()
        for p in bundle.rglob("*")
        if p.is_file()
    } == before
    assert result["status"] == "passed"
    assert result["complete_search_replays"] == 480
    assert result["engine_scope"] == "archived engine; not current CLI validation"
    snapshot = bundle / "source-snapshot" / "shift_calibration.py"
    snapshot.write_text(snapshot.read_text() + "\n# altered snapshot\n")
    with pytest.raises(ValueError, match="Snapshot hash mismatch"):
        audit(bundle, frozen_engine=True)


def test_replay_mismatch_reports_identity_without_changing_rejection(monkeypatch):
    import verify_shift_null_contract as verifier

    monkeypatch.setattr(verifier, "same_replayed_model", lambda *_: False)
    with pytest.raises(
        ValueError, match="Seeded complete-search replay disagrees"
    ) as exc:
        audit(
            Path("examples/shift/null-contract-timing-pilot"),
            replay_stride=1,
            frozen_engine=True,
        )
    detail = json.loads(str(exc.value).split(": ", 1)[1])
    assert set(detail) == {
        "cell_id",
        "replicate",
        "convergence",
        "bootstrap_seed",
        "model_matches",
        "tests_match",
        "replayed_tests",
        "saved_tests",
    }
    assert detail["cell_id"] == 0 and detail["replicate"] == 0
    assert detail["model_matches"] is False
    assert type(detail["bootstrap_seed"]) is int


def test_replay_uses_exact_saved_inputs_after_tolerant_regeneration(
    monkeypatch, tmp_path
):
    import verify_shift_null_contract as verifier

    bundle = tmp_path / "one-dataset"
    shutil.copytree(Path("examples/shift/null-contract-timing-pilot"), bundle)
    with gzip.open(bundle / "records.jsonl.gz", "rt") as stream:
        record = json.loads(next(stream))
    # A single four-tip dataset is enough to check replay input fidelity. The
    # separate archived-engine test retains the full 480 real replays.
    specification = json.loads((bundle / "protocol.json").read_text())
    specification["cells"] = specification["cells"][:1]
    specification["replicates"] = 1
    (bundle / "protocol.json").write_text(json.dumps(specification))
    with gzip.open(bundle / "records.jsonl.gz", "wt") as stream:
        stream.write(json.dumps(record) + "\n")
    (bundle / "summary.json").write_text(
        json.dumps(verifier.summarize([record], specification))
    )
    compare = verifier.same_generated_array

    def regenerated_with_roundoff(actual, expected):
        if np.ndim(actual) == 1:
            # Simulate another BLAS's valid, near-machine-precision input.
            actual[:] = np.nextafter(np.asarray(expected), np.inf)
        return compare(actual, expected)

    class Replay:
        def __init__(self, tree, convergence, variances):
            self.convergence = convergence

        def fit(self, values, seed, replicates):
            assert seed == record["bootstrap_seed"]
            np.testing.assert_array_equal(values, record["observations"])
            lane = next(
                lane
                for lane in record["lanes"]
                if lane["convergence"] == self.convergence
            )
            assert replicates == lane["fit"]["calibration_replicates"]
            return lane["fit"]

    class ReplayModule(ModuleType):
        @property
        def CalibratedSearch(self):
            # The trusted snapshot is still compiled and verified; replace
            # only bootstrap fitting with an oracle for the supplied inputs.
            return Replay

    monkeypatch.setattr(verifier, "same_generated_array", regenerated_with_roundoff)
    monkeypatch.setattr(verifier, "ModuleType", ReplayModule)
    result = verifier.audit(bundle, replay_stride=1, frozen_engine=True)
    assert result["complete_search_replays"] == 2


def test_archived_mode_still_rejects_revised_generator(tmp_path):
    bundle = tmp_path / "archived"
    shutil.copytree(Path("examples/shift/null-contract-timing-pilot"), bundle)
    snapshot = bundle / "source-snapshot" / "shift_continuous_reference.py"
    snapshot.write_text(snapshot.read_text() + "\n# different generator\n")
    spec = json.loads((bundle / "protocol.json").read_text())
    spec["source_sha256"]["tools/shift_continuous_reference.py"] = hashlib.sha256(
        snapshot.read_bytes()
    ).hexdigest()
    (bundle / "protocol.json").write_text(json.dumps(spec))
    with pytest.raises(ValueError, match="Active generator/fitting hash mismatch"):
        audit(bundle, frozen_engine=True)


def test_archived_mode_never_executes_bundle_supplied_engine(tmp_path):
    bundle = tmp_path / "archived"
    shutil.copytree(Path("examples/shift/null-contract-timing-pilot"), bundle)
    marker = tmp_path / "executed"
    snapshot = bundle / "source-snapshot" / "shift_calibration.py"
    snapshot.write_text(f"from pathlib import Path\nPath({str(marker)!r}).touch()\n")
    specification = json.loads((bundle / "protocol.json").read_text())
    specification["source_sha256"]["nwkit/shift_calibration.py"] = hashlib.sha256(
        snapshot.read_bytes()
    ).hexdigest()
    (bundle / "protocol.json").write_text(json.dumps(specification))
    with pytest.raises(ValueError, match="not a trusted checked-in snapshot"):
        audit(bundle, frozen_engine=True)
    assert not marker.exists()
