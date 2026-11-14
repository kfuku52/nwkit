"""Failure checks for the native optimization measurement protocol."""

import copy
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

RUNNER_PATH = (
    Path(__file__).resolve().parents[1]
    / "reviews/native-optimization-2026-09-11/run_benchmarks.py"
)
SPEC = importlib.util.spec_from_file_location(
    "native_optimization_benchmark", RUNNER_PATH
)
RUNNER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RUNNER)


def manifest():
    return {
        "python": "3.10",
        "numpy": "1.26",
        "scipy": "1.15",
        "platform": "test",
        "thread_limits": {"OMP_NUM_THREADS": "1"},
        "harness_sha256": {"runner": "abc"},
        "source_sha256": {"before": {"nwkit/__init__.py": "abc"}},
    }


@pytest.mark.parametrize("field", list(manifest()))
def test_reuse_rejects_mismatched_or_missing_provenance(field):
    current = manifest()
    RUNNER.validate_baseline(copy.deepcopy(current), current)
    prior = copy.deepcopy(current)
    prior[field] = {} if field == "source_sha256" else "changed"
    with pytest.raises(ValueError, match="mismatch"):
        RUNNER.validate_baseline(prior, current)
    del prior[field]
    with pytest.raises(ValueError, match="mismatch"):
        RUNNER.validate_baseline(prior, current)


def test_checks_are_not_disabled_by_python_optimization():
    command = (
        "import runpy; m=runpy.run_path(" + repr(str(RUNNER_PATH)) + "); "
        "m['validate_baseline']({}, " + repr(manifest()) + ")"
    )
    result = subprocess.run(
        [sys.executable, "-O", "-c", command], capture_output=True, text=True
    )
    assert result.returncode != 0
    assert "Baseline measurement mismatch" in result.stderr


@pytest.mark.parametrize("value", [np.nan, np.inf, -np.inf])
def test_matching_nonfinite_arrays_are_not_equivalent_evidence(value):
    arrays = {"x0": np.array([[value]])}
    with pytest.raises(ValueError, match="Non-finite"):
        RUNNER.compare_screen_arrays(arrays, arrays)


def test_screening_comparison_checks_keys_and_values():
    arrays = {"x0": np.array([[1.0, 2.0]])}
    assert RUNNER.compare_screen_arrays(arrays, arrays) == 0
    with pytest.raises(ValueError, match="keys"):
        RUNNER.compare_screen_arrays(arrays, {})
    with pytest.raises(AssertionError):
        RUNNER.compare_screen_arrays(arrays, {"x0": np.array([[1.0, 3.0]])})


def test_search_comparison_checks_budgets_even_when_results_match():
    path = RUNNER_PATH.parent / "final-measurements/balanced128-2-est-before-0.json"
    reference = json.loads(path.read_text())
    result = copy.deepcopy(reference)
    result["configuration"]["output"] = "another-output.json"
    RUNNER.compare_search_results(reference, result, "case")
    result["configuration"]["refit_budget"] += 1
    with pytest.raises(ValueError, match="configuration"):
        RUNNER.compare_search_results(reference, result, "case")
    result = copy.deepcopy(reference)
    result["best_log_likelihood"] += 1
    with pytest.raises(ValueError, match="best_log_likelihood"):
        RUNNER.compare_search_results(reference, result, "case")


def test_missing_source_cannot_silently_use_installed_package(tmp_path):
    output = tmp_path / "output"
    result = subprocess.run(
        [
            sys.executable,
            str(RUNNER_PATH),
            "--before",
            str(tmp_path / "missing"),
            "--after",
            str(tmp_path / "missing"),
            "--output",
            str(output),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "must contain a nwkit source package" in result.stderr
    assert not output.exists()


@pytest.mark.parametrize(
    "extra,suffix",
    [([], ".npz"), (["--tips", "0"], ".json"), (["--traits", "0"], ".json")],
)
@pytest.mark.skipif(
    sys.platform == "win32", reason="The resource-based benchmark is POSIX-only."
)
def test_screen_harness_rejects_invalid_sizes_and_colliding_outputs(
    tmp_path, extra, suffix
):
    tool = RUNNER_PATH.parents[2] / "tools/benchmark_native_screen.py"
    output = tmp_path / ("result" + suffix)
    result = subprocess.run(
        [sys.executable, str(tool), "--output", str(output), *extra],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 2
    assert not output.exists()
