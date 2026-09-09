import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest


def load_tool(name):
    spec = importlib.util.spec_from_file_location(
        f"nwkit_check_{name}",
        Path(__file__).resolve().parents[1] / "tools" / f"{name}.py",
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


check = load_tool("check")
maintainability = load_tool("check_maintainability")
ci_matrix = load_tool("ci_matrix")
ete_build = load_tool("build_ete_windows")


def test_check_runner_bootstraps_from_source_without_installed_dependencies(tmp_path):
    result = subprocess.run(
        [sys.executable, "-S", str(Path(check.__file__)), "--help"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "release" in result.stdout


def test_quick_runner_passes_targets_and_uses_incremental_types(monkeypatch):
    commands = []
    monkeypatch.setattr(
        check, "run", lambda *command, **kwargs: commands.append(command)
    )
    check.main(["quick", "--", "tests/test_asr.py", "-k", "time_units"])
    assert (check.PYTHON, "-m", "mypy") in commands
    assert commands[-1][-5:] == (
        "-m",
        "not slow",
        "tests/test_asr.py",
        "-k",
        "time_units",
    )


def test_quick_runner_allows_explicit_slow_selection(monkeypatch):
    commands = []
    monkeypatch.setattr(
        check, "run", lambda *command, **kwargs: commands.append(command)
    )
    check.main(["quick", "--", "-m", "slow", "tests/test_asr.py"])
    assert "not slow" not in commands[-1]
    assert commands[-1][-3:] == ("-m", "slow", "tests/test_asr.py")


@pytest.mark.parametrize("mode", ["full", "release", "dist"])
def test_full_validation_cannot_be_mistaken_for_a_selected_suite(monkeypatch, mode):
    monkeypatch.setattr(
        check,
        "run",
        lambda *args, **kwargs: pytest.fail(
            "validation started before rejecting selection"
        ),
    )
    with pytest.raises(SystemExit) as exc:
        check.main([mode, "--", "-k", "one_test"])
    assert exc.value.code == 2


def test_full_checks_keep_the_complete_suite_and_uncached_type_checks(monkeypatch):
    commands = []
    monkeypatch.setattr(
        check, "run", lambda *command, **kwargs: commands.append(command)
    )
    check.main(["full"])
    assert (check.PYTHON, "-m", "mypy", "--no-incremental") in commands
    assert (
        check.PYTHON,
        "-m",
        "coverage",
        "run",
        "-m",
        "pytest",
        "tests/",
        "-q",
    ) in commands


def test_complexity_cleanup_is_not_penalized_for_raising_the_average():
    baseline = {"module:large": 20, "module:small": 1}
    assert maintainability.complexity_increases({"module:large": 19}, baseline) == []
    assert (
        maintainability.complexity_increases(
            {"module:large": 20, "module:new": 2}, baseline
        )
        == []
    )


def test_complexity_growth_warns_but_common_limit_still_rejects():
    assert maintainability.complexity_increases(
        {"module:existing": 4}, {"module:existing": 3}
    )
    assert maintainability.complexity_violations({"module:existing": 4}) == []
    assert maintainability.complexity_violations({"module:existing": 41})
    assert maintainability.complexity_violations({"module:new": 40}) == []
    assert maintainability.complexity_violations({"module:new": 41}, {})


def test_documented_exception_has_a_hard_ceiling():
    exception = {
        "module:legacy": {
            "limit": 50,
            "reason": "Legacy orchestration awaiting responsibility separation.",
            "tests": ["tests/test_draw.py"],
        }
    }
    maintainability.validate_exceptions(exception, {"module:legacy": 50})
    assert maintainability.complexity_violations({"module:legacy": 50}, exception) == []
    assert maintainability.complexity_violations({"module:legacy": 51}, exception)


@pytest.mark.parametrize(
    "record",
    [
        {"limit": 50, "reason": "", "tests": ["tests/test_draw.py"]},
        {"limit": 40, "reason": "Reason", "tests": ["tests/test_draw.py"]},
        {"limit": 50, "reason": "Reason", "tests": []},
        {"limit": 50, "reason": "Reason", "tests": ["tests/missing_test.py"]},
        {"limit": 50, "reason": "Reason", "tests": ["tests/../setup.py"]},
    ],
)
def test_complexity_exception_requires_a_rationale_and_real_tests(record):
    with pytest.raises(ValueError):
        maintainability.validate_exceptions(
            {"module:legacy": record}, {"module:legacy": 50}
        )


def test_complexity_baseline_update_cannot_bypass_hard_limits(
    tmp_path, monkeypatch, capsys
):
    baseline = tmp_path / "baseline.json"
    exceptions = tmp_path / "exceptions.json"
    baseline.write_text('{"module:existing": 16}')
    exceptions.write_text("{}")
    monkeypatch.setattr(maintainability, "BASELINE_PATH", baseline)
    monkeypatch.setattr(maintainability, "EXCEPTIONS_PATH", exceptions)
    current = {"module:existing": 20}
    monkeypatch.setattr(maintainability, "collect_complexities", lambda: current)
    assert maintainability.main([]) == 0
    assert "increased from 16 to 20" in capsys.readouterr().err
    assert json.loads(baseline.read_text())["module:existing"] == 16
    assert maintainability.main(["--update-baseline"]) == 0
    assert json.loads(baseline.read_text())["module:existing"] == 20
    current["module:existing"] = 41
    with pytest.raises(RuntimeError, match="exceeds limit 40"):
        maintainability.main(["--update-baseline"])
    assert json.loads(baseline.read_text())["module:existing"] == 20
    assert exceptions.read_text() == "{}"


def test_complexity_keys_distinguish_methods_and_nested_functions():
    results = {
        "nwkit/example.py": [
            {
                "type": "class",
                "name": "Model",
                "complexity": 3,
                "methods": [
                    {
                        "type": "method",
                        "name": "fit",
                        "complexity": 2,
                        "closures": [
                            {"type": "function", "name": "objective", "complexity": 1}
                        ],
                    }
                ],
            },
            {"type": "function", "name": "fit", "complexity": 3},
        ]
    }
    assert maintainability.function_complexities(results) == {
        "nwkit/example.py:Model.fit": 2,
        "nwkit/example.py:Model.fit.<locals>.objective": 1,
        "nwkit/example.py:fit": 3,
    }


def test_ci_keeps_minimum_and_latest_python_for_numerical_changes():
    selected = ci_matrix.select_coverage(
        ["nwkit/gaussian.py", "tests/test_numerical_invariance.py"]
    )
    assert selected["source_checks"]  # quality job runs all tests on Python 3.14
    assert selected["matrix"]["include"] == [
        {"os": "ubuntu-latest", "python-version": "3.10", "extras": "test,image"}
    ]
    assert not selected["macos_clean"]


def test_ci_docs_only_skips_numerical_suites_but_dependencies_restore_all_versions():
    docs = ci_matrix.select_coverage(["README.md", "PHYLOGENETIC_REGRESSION.md"])
    assert not docs["source_checks"] and docs["matrix"]["include"] == []
    dependencies = ci_matrix.select_coverage(["pyproject.toml"])
    assert len(dependencies["matrix"]["include"]) == 6
    assert dependencies["macos_clean"]
    assert ci_matrix.select_coverage([], full=True) == dependencies


def test_ci_filesystem_changes_keep_native_os_coverage():
    selected = ci_matrix.select_coverage(["nwkit/output_transaction.py"])
    assert {row["os"] for row in selected["matrix"]["include"]} == {
        "ubuntu-latest",
        "macos-latest",
        "windows-latest",
    }
    assert selected["macos_clean"]


def test_ci_ignores_only_a_plain_version_change():
    assert ci_matrix.version_only_change(
        '__version__ = "0.40.5"', '__version__ = "0.40.6"'
    ) == (True, "0.40.6")
    assert not ci_matrix.version_only_change(
        '__version__ = "0.40.5"', '__version__ = "0.40.6"\nimport numpy'
    )[0]
    assert not ci_matrix.version_only_change('__version__ = "0.40.5"', "")[0]


@pytest.mark.parametrize("path", ["ete4/core/tree.pyx", r"ete4\core\tree.pyx"])
def test_verified_ete_patch_produces_importable_module_names(path):
    namespace = {"path": path}
    exec(
        ete_build.patch_setup("name = path.replace('/', '.')[:-len('.pyx')]"), namespace
    )
    assert namespace["name"] == "ete4.core.tree"


def test_ete_patch_refuses_an_unreviewed_source_change():
    with pytest.raises(ValueError, match="revalidate or remove"):
        ete_build.patch_setup("new upstream implementation")
