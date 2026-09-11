"""Independent evidence checks reject numerically plausible corruptions."""

import copy
import hashlib
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from compare_shift_null_calibration import summarize  # noqa: E402
from validate_shift_null_contract import tree_text  # noqa: E402
from verify_shift_null_contract import (  # noqa: E402
    audit,
    check_stages,
    check_winner,
    same_generated_array,
    same_replayed_summary,
    same_replayed_tests,
)

from nwkit.shift_calibration import CalibratedSearch  # noqa: E402
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
    result = audit(bundle, replay_stride=1, frozen_engine=True)
    assert {
        str(p.relative_to(bundle)): p.read_bytes()
        for p in bundle.rglob("*")
        if p.is_file()
    } == before
    assert result["status"] == "passed"
    assert result["engine_scope"] == "archived engine; not current CLI validation"
    snapshot = bundle / "source-snapshot" / "shift_calibration.py"
    snapshot.write_text(snapshot.read_text() + "\n# altered snapshot\n")
    with pytest.raises(ValueError, match="Snapshot hash mismatch"):
        audit(bundle, frozen_engine=True)


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
