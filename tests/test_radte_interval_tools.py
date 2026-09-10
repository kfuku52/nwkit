import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

from nwkit.cli import main
from tests.test_radte import cli_inputs


@pytest.fixture
def runner(monkeypatch):
    tools = Path(__file__).resolve().parents[1] / "tools"
    monkeypatch.syspath_prepend(str(tools))
    spec = importlib.util.spec_from_file_location(
        "interval_runner", tools / "validate_radte_intervals.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_summary_does_not_hide_unavailable_intervals(runner):
    rows = [
        dict(status="completed", interval_available=True, covered=True, width=2),
        dict(status="completed", interval_available=False, covered=False),
        dict(status="failed", interval_available=False, covered=False),
    ]
    result = runner.interval_summary(rows)
    assert result["families"] == 3
    assert result["coverage_among_available"] == 1
    assert result["interval_return_fraction"] == pytest.approx(1 / 3)
    assert result["correct_interval_returned_fraction"] == pytest.approx(1 / 3)
    assert "availability_wilson_95" in result and "correct_return_wilson_95" in result


def test_interval_exception_does_not_erase_other_methods(runner, monkeypatch, tmp_path):
    from radte_interval_simulation import simulate

    directory = tmp_path / "input"
    simulate(directory, 2, 50, 10, 0.3)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "validate",
            "--output",
            str(tmp_path / "out"),
            "--branch-only",
            "--targets",
            "D",
            "absent",
        ],
    )
    options = runner.options()

    def fail(*args, **kwargs):
        raise ValueError("interval failure")

    monkeypatch.setattr(runner, "laplace_intervals", fail)
    rows = runner.evaluate_family(directory, options, 0)
    indexed = {(row["method"], row["target"]): row for row in rows}
    assert (
        indexed["laplace", "D"]["interval_status"] == "unavailable-interval-exception"
    )
    assert indexed["studentized", "D"]["status"] == "completed"
    assert indexed["studentized", "D"]["error"] is None
    assert indexed["studentized", "absent"]["status"] == "target-unmatched"


def test_validation_role_requires_frozen_protocol(runner, monkeypatch, tmp_path):
    monkeypatch.setattr(
        sys,
        "argv",
        ["validate", "--output", str(tmp_path), "--study-role", "validation"],
    )
    with pytest.raises(SystemExit):
        runner.options()


def test_exact_interval_cli_preserves_points_and_records_method(tmp_path):
    inputs = cli_inputs(tmp_path)
    for method in ["none", "exact-log-duration"]:
        main(
            [
                "radte",
                *inputs,
                "--reconcile",
                "lca",
                "--uncertainty",
                method,
                "--out-prefix",
                str(tmp_path / method),
            ]
        )
    original = pd.read_csv(tmp_path / "none.nodes.tsv", sep="\t")
    exact = pd.read_csv(tmp_path / "exact-log-duration.nodes.tsv", sep="\t")
    pd.testing.assert_series_equal(original.estimated_age, exact.estimated_age)
    manifest = json.loads((tmp_path / "exact-log-duration.manifest.json").read_text())
    assert manifest["uncertainty"] == "conditional-exact-log-duration-t"
    assert (exact.interval_lower <= exact.estimated_age).all()
    assert (exact.interval_upper >= exact.estimated_age).all()
