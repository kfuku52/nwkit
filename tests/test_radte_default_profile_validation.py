"""Keep validation denominators honest when estimation/intervals fail."""

import importlib.util
from pathlib import Path

import pytest


@pytest.fixture
def runner(monkeypatch):
    root = Path(__file__).resolve().parents[1]
    monkeypatch.syspath_prepend(str(root / "tools"))
    spec = importlib.util.spec_from_file_location(
        "default_profile_validation", root / "tools/validate_radte_default_profile.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_failure_denominators_are_distinct(runner):
    rows = [
        dict(
            point_success=True,
            interval_available=True,
            covered=True,
            bias=1,
            width=4,
            status="complete",
            profile_estimator="map",
        ),
        dict(
            point_success=True,
            interval_available=True,
            covered=False,
            bias=-1,
            width=2,
            status="complete",
            profile_estimator="map",
        ),
        dict(
            point_success=True,
            interval_available=False,
            covered=False,
            bias=3,
            status="profile-timeout",
            none_estimator="map",
        ),
        dict(
            point_success=False,
            interval_available=False,
            covered=False,
            status="none-exit-1",
        ),
    ]
    result = runner.summarize(rows)
    assert result["point_successes"] == 3
    assert result["availability"] == 0.5
    assert result["coverage_among_returned"] == 0.5
    assert result["correct_return_fraction"] == 0.25
    assert result["bias"] == 1
    assert result["rmse"] == pytest.approx((11 / 3) ** 0.5)
    assert result["median_width"] == 3
    assert result["statuses"]["profile-timeout"] == 1


def test_no_intervals_is_missing_coverage_not_zero_width(runner):
    result = runner.summarize(
        [
            dict(
                point_success=False,
                interval_available=False,
                covered=False,
                status="none-timeout",
            )
        ]
    )
    assert result["coverage_among_returned"] is None
    assert result["median_width"] is None
    assert result["correct_return_fraction"] == 0
    assert result["coverage_wilson95"] is None


def test_generated_chronology_retains_independent_truth(runner, tmp_path):
    from ete4 import Tree

    truth = runner.generate_trees(
        tmp_path, dict(species=4, rate_sd=0.3, scenario="nested"), 171
    )
    assert truth == 7.5
    genealogy = Tree(str(tmp_path / "truth.nwk"), parser=1)
    duplication = next(n for n in genealogy.traverse() if n.name == "D")
    assert {genealogy.get_distance(duplication, t) for t in duplication.leaves()} == {
        truth
    }


@pytest.mark.parametrize("failed_method", ["none", "profile"])
def test_nonfinite_estimate_is_retained_as_failure(
    runner, monkeypatch, tmp_path, failed_method
):
    import json
    from types import SimpleNamespace

    def generate(directory, case, seed):
        for name in ["gene.nwk", "species.nwk", "mapping.tsv", "truth.json"]:
            (directory / name).write_text("fixture")
        (directory / "alignment.fa").write_text(">A\nAAA\n")
        return 20.0

    def read_target(prefix):
        value = "nan" if prefix.name == failed_method else "21.0"
        return {"method": "map", "diagnostics": []}, {"estimated_age": value}

    monkeypatch.setattr(runner, "generate_trees", generate)
    monkeypatch.setattr(runner, "run", lambda *args: ("complete", 0.0))
    monkeypatch.setattr(runner, "read_target", read_target)
    row = runner.evaluate(
        {"name": "root", "seed": 12, "codons": 1},
        0,
        SimpleNamespace(output=tmp_path, iqtree="unused", timeout=1),
    )
    assert row["status"] == failed_method + "-nonfinite-estimate"
    assert row["point_success"] == (failed_method == "profile")
    assert not row["interval_available"]
    assert not row["covered"]
    if failed_method == "profile":
        assert row["estimate"] == 21.0
        assert row["bias"] == 1.0
    assert json.loads((tmp_path / "root/f0000/result.json").read_text()) == row
    assert runner.summarize([row])["correct_return_fraction"] == 0


def test_oracle_nonfinite_estimate_is_a_failure(runner, monkeypatch, tmp_path):
    import json
    from types import SimpleNamespace

    import diagnose_radte_default_profile as oracle

    source = tmp_path / "source"
    source.mkdir()
    (source / "result.json").write_text(json.dumps({"case": "root", "family": 0}))
    (source / "truth.json").write_text(
        json.dumps({"age": 20, "case": {"rate_sd": 0.3}})
    )
    (source / "profile.command.json").write_text(
        json.dumps(["python", "--uncertainty", "profile", "--out-prefix", "old"])
    )
    monkeypatch.setattr(
        oracle.subprocess, "run", lambda *a, **k: SimpleNamespace(returncode=0)
    )
    monkeypatch.setattr(
        oracle,
        "read_target",
        lambda p: ({"method": "map", "diagnostics": []}, {"estimated_age": "nan"}),
    )
    result = oracle.replay(source, tmp_path / "output", 1)
    assert result["status"] == "nonfinite-estimate"
    assert not result["point_success"]
    assert not result["interval_available"]
    assert json.loads((tmp_path / "output/result.json").read_text()) == result
