"""Independent-generator and replay checks for the signal calibration study."""

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

from nwkit.signal_stats import lambda_fit

PATH = Path(__file__).resolve().parents[1] / "tools/validate_signal_calibration.py"
SPEC = importlib.util.spec_from_file_location("validate_signal_calibration", PATH)
assert SPEC is not None and SPEC.loader is not None
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)
VERIFY_PATH = PATH.with_name("verify_signal_calibration.py")
VERIFY_SPEC = importlib.util.spec_from_file_location(
    "verify_signal_calibration", VERIFY_PATH
)
assert VERIFY_SPEC is not None and VERIFY_SPEC.loader is not None
verifier = importlib.util.module_from_spec(VERIFY_SPEC)
VERIFY_SPEC.loader.exec_module(verifier)


def test_study_covariances_represent_distinct_tree_shapes():
    balanced = module._covariance(8, "balanced")
    pectinate = module._covariance(8, "pectinate")
    assert balanced[0, 1] == 0.5
    assert balanced[0, 2] == 0.25
    assert balanced[0, 4] == 0
    assert pectinate[6, 7] > pectinate[1, 7] > 0
    assert np.min(np.linalg.eigvalsh(balanced)) > 0
    assert np.min(np.linalg.eigvalsh(pectinate)) > 0
    covariance, errors, truth = module._scenario("balanced-8-missing-null")
    assert covariance.shape == (6, 6)
    assert errors.shape == (6,)
    assert truth == 0
    _, known_errors, _ = module._scenario("balanced-8-known-se-null")
    assert np.all(known_errors > 0)


def test_no_se_lambda_ratio_is_location_scale_invariant():
    for name in ("balanced-8-null", "pectinate-8-null", "balanced-32-null"):
        covariance, errors, _ = module._scenario(name)
        values = np.random.default_rng(module._seed(173, name, 0, 7)).normal(
            size=len(errors)
        )
        original = lambda_fit(covariance, values, errors, ci_level=None)
        transformed = lambda_fit(covariance, 7 - 3 * values, errors, ci_level=None)
        np.testing.assert_allclose(
            transformed["likelihood_ratio"],
            original["likelihood_ratio"],
            rtol=0,
            atol=1e-8,
        )


def test_signal_study_replays_independently_of_scenario_order(tmp_path):
    first, second = tmp_path / "first", tmp_path / "second"
    options = ["--outer", "2", "--inner", "9", "--seed", "173"]
    module.main(
        [
            "--output",
            str(first),
            "--scenarios",
            "balanced-8-null,pectinate-8-null",
            *options,
        ]
    )
    module.main(
        [
            "--output",
            str(second),
            "--scenarios",
            "pectinate-8-null,balanced-8-null",
            *options,
        ]
    )
    rows = []
    for output in (first, second):
        protocol = json.loads((output / "protocol.json").read_text())
        assert protocol["completed"]
        records = [
            json.loads(line)
            for line in (output / "records.jsonl").read_text().splitlines()
        ]
        assert len(records) == 4
        assert all(0 <= row["bootstrap_p_value"] <= 1 for row in records)
        rows.append(
            sorted(records, key=lambda row: (row["scenario"], row["replicate"]))
        )
    assert rows[0] == rows[1]


def test_bootstrap_failure_keeps_chi2_and_profile_result(monkeypatch):
    def fail(*args):
        raise ValueError("injected bootstrap failure")

    monkeypatch.setattr(module, "lambda_parametric_bootstrap", fail)
    row = module._one(("balanced-8-null", 0, 173, 9, True))
    assert row["status"] in {"ok", "boundary"}
    assert "chi2_p_value" in row and "bootstrap_p_value" not in row
    assert row["bootstrap_error"] == "injected bootstrap failure"
    summary = module._summarize([row], ["balanced-8-null"], 0.05)[0]
    assert summary["chi2_p_available"]["count"] == 1
    assert summary["bootstrap_p_available"]["count"] == 0
    assert summary["profile_interval_available"]["count"] == 1


def test_interrupted_study_keeps_incomplete_protocol_and_finished_rows(
    monkeypatch, tmp_path
):
    original = module._one
    calls = 0

    def stop_after_first(task):
        nonlocal calls
        calls += 1
        if calls > 1:
            raise RuntimeError("injected interruption")
        return original(task)

    monkeypatch.setattr(module, "_one", stop_after_first)
    output = tmp_path / "interrupted"
    with pytest.raises(RuntimeError, match="injected interruption"):
        module.main(
            [
                "--output",
                str(output),
                "--scenarios",
                "balanced-8-null",
                "--outer",
                "2",
                "--inner",
                "1",
            ]
        )
    assert not json.loads((output / "protocol.json").read_text())["completed"]
    assert len((output / "records.jsonl").read_text().splitlines()) == 1


def test_signal_evidence_audit_rejects_tampered_summary(tmp_path):
    output = tmp_path / "study"
    module.main(
        [
            "--output",
            str(output),
            "--scenarios",
            "balanced-8-null",
            "--outer",
            "2",
            "--inner",
            "9",
        ]
    )
    assert verifier.verify(output) == 2
    summary = json.loads((output / "summary.json").read_text())
    summary[0]["bootstrap_rejection_all"]["count"] += 1
    (output / "summary.json").write_text(json.dumps(summary))
    with pytest.raises(ValueError, match="Rate numerator or denominator mismatch"):
        verifier.verify(output)


def test_profile_only_study_skips_bootstrap_and_keeps_chi2(tmp_path):
    output = tmp_path / "profile"
    module.main(
        [
            "--output",
            str(output),
            "--scenarios",
            "balanced-8-alternative",
            "--outer",
            "2",
            "--inner",
            "0",
            "--profile-ci",
        ]
    )
    assert verifier.verify(output) == 2
    rows = [
        json.loads(line) for line in (output / "records.jsonl").read_text().splitlines()
    ]
    assert all("ci_lower" in row and "chi2_p_value" in row for row in rows)
    assert all("bootstrap_p_value" not in row for row in rows)
    summary = json.loads((output / "summary.json").read_text())[0]
    assert summary["bootstrap_p_available"]["count"] == 0
    assert summary["profile_interval_available"]["count"] == 2
