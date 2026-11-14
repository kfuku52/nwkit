"""Check independent simulation truth, failure denominators and label invariance."""

import importlib.util
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from nwkit.shift_reference import evaluate_shift_model
from nwkit.util import assign_branch_ids, read_tree

TOOLS = Path(__file__).resolve().parents[1] / "tools"
spec = importlib.util.spec_from_file_location(
    "shift_simulation_cases", TOOLS / "shift_simulation_cases.py"
)
cases = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = cases
spec.loader.exec_module(cases)
spec_runner = importlib.util.spec_from_file_location(
    "validate_shift_simulations", TOOLS / "validate_shift_simulations.py"
)
runner = importlib.util.module_from_spec(spec_runner)
spec_runner.loader.exec_module(runner)


@pytest.mark.parametrize("root_model", ["OUfixedRoot", "OUrandomRoot"])
@pytest.mark.parametrize("scenario", cases.SCENARIOS)
def test_independent_truth_matches_gaussian_reference(root_model, scenario):
    newick, truth = cases.simulate(
        tips=16, scenario=scenario, root_model=root_model, se=0.2, effect=2, seed=1
    )
    tree = read_tree(newick, "auto", True, quiet=True)
    ids = assign_branch_ids(tree)
    alpha = truth["alpha"]
    effects = {}
    for node, branch in ids.items():
        if str(branch) in truth["shift_optima"]:
            # These disconnected shifts begin at the same depth. The reference
            # accepts induced tip-mean changes, whereas the generator uses optima.
            duration = node.dist + node.get_farthest_leaf()[1]
            effects[node] = truth["shift_optima"][str(branch)] * (
                -np.expm1(-alpha * duration)
            )
    mean, cov, _ = evaluate_shift_model(
        tree,
        observations=truth["observations"],
        standard_errors=[0.2] * 16,
        alpha=alpha,
        sigma2=truth["sigma2"],
        intercept=0,
        mean_effects=effects,
        root_model=root_model,
    )
    np.testing.assert_allclose(mean, list(truth["tip_mean"].values()), atol=1e-14)
    np.testing.assert_allclose(cov, truth["covariance"], atol=1e-14)
    assert (
        len(truth["shift_branch_ids"])
        == {"null": 0, "single": 1, "distinct": 2, "convergent": 2}[scenario]
    )
    assert (
        len(truth["shared_partition"])
        == {"null": 1, "single": 2, "distinct": 3, "convergent": 2}[scenario]
    )


@pytest.mark.parametrize("root_model", ["OUfixedRoot", "OUrandomRoot"])
def test_innovation_draws_have_expected_moments(root_model):
    draws = []
    for seed in range(2000):
        _, truth = cases.simulate(
            tips=8,
            scenario="convergent",
            root_model=root_model,
            se=0.2,
            effect=2,
            seed=seed,
        )
        draws.append(truth["observations"])
    np.testing.assert_allclose(
        np.mean(draws, axis=0), list(truth["tip_mean"].values()), atol=0.04
    )
    np.testing.assert_allclose(
        np.cov(draws, rowvar=False), truth["covariance"], atol=0.035
    )


def test_failures_are_not_removed_from_returned_recovery_rate():
    good = {
        "status": "completed",
        "any_shift": True,
        "any_merge": False,
        "exact_edges": True,
        "ancestry_recovered": True,
        "shared_recovered": True,
        "tip_mean_rmse": 0.1,
    }
    result = cases.summarize([good, {"status": "failed"}])
    assert result["shared_recovered"]["rate"] == 1
    assert result["shared_recovered_returned"]["rate"] == 0.5
    assert result["failure"]["count"] == 1
    assert result["shared_recovered"]["wilson_95"][0] < 0.3
    empty = cases.summarize([{"status": "failed"}])
    assert empty["shared_recovered"]["rate"] is None
    assert empty["shared_recovered_returned"]["rate"] == 0
    assert empty["bootstrap_truth_shared_partition"]["mean_frequency"] is None


def test_partition_ignores_regime_labels_and_input_order():
    assert cases.partition({"b": 1, "a": 1, "c": 0}) == cases.partition(
        {"c": "x", "a": "y", "b": "y"}
    )


def test_fit_preserves_failure_record(tmp_path, monkeypatch):
    def fail(args):
        raise ValueError("all bootstrap replicates failed.")

    monkeypatch.setattr(
        runner.parser, "parse_args", lambda args: SimpleNamespace(handler=fail)
    )
    _, truth = cases.simulate(
        tips=8, scenario="null", root_model="OUfixedRoot", se=0, effect=2, seed=4
    )
    row = runner.fit(
        tmp_path,
        truth,
        SimpleNamespace(rscript="Rscript", criterion="BIC", bootstrap=3),
        "bootstrap",
    )
    assert row["status"] == "failed"
    assert "all bootstrap" in row["error"]
    assert (tmp_path / "bootstrap" / "result.json").exists()


@pytest.mark.parametrize("tips,se", [(7, 0), (12, 0), (8, -1), (8, float("nan"))])
def test_invalid_generating_parameters(tips, se):
    with pytest.raises(ValueError):
        cases.simulate(
            tips=tips,
            scenario="null",
            root_model="OUfixedRoot",
            se=se,
            effect=2,
            seed=0,
        )


@pytest.mark.parametrize("converged", [False, True])
def test_metrics_distinguish_ancestry_from_shared_optima(converged):
    newick, truth = cases.simulate(
        tips=8, scenario="convergent", root_model="OUfixedRoot", se=0, effect=2, seed=2
    )
    tree = read_tree(newick, "auto", True, quiet=True)
    ids = assign_branch_ids(tree)
    labels = {}
    rows = []
    for node in tree.traverse("preorder"):
        labels[node] = (
            ids[node]
            if node.is_root or ids[node] in truth["shift_branch_ids"]
            else labels[node.up]
        )
        if node.is_leaf:
            rows.append(
                {
                    "leaf_name": node.name,
                    "branch_id": ids[node],
                    "regime": truth["tip_optimum"][node.name]
                    if converged
                    else labels[node],
                    "predicted": truth["tip_mean"][node.name],
                }
            )
    model = {
        "tip_predictions": rows,
        "branches": [
            {"branch_id": ids[n], "parent": -1 if n.is_root else ids[n.up]} for n in ids
        ],
        "shift_branch_ids": truth["shift_branch_ids"],
        "parameters": {"alpha": 0.7},
        "convergence": {"merges": 1} if converged else None,
        "backend_version": "test",
    }
    result = cases.selection_metrics(model, truth)
    assert result["exact_edges"] and result["ancestry_recovered"]
    assert result["shared_recovered"] == converged
    assert result["tip_mean_rmse"] == 0
    # Matching regime labels at the BM boundary cannot recover OU optima.
    model["parameters"]["alpha"] = 0
    assert not cases.selection_metrics(model, truth)["shared_recovered"]


@pytest.mark.parametrize("effect", [0, -1, float("inf")])
def test_effect_is_a_positive_magnitude(effect):
    with pytest.raises(ValueError):
        cases.simulate(
            tips=8,
            scenario="single",
            root_model="OUfixedRoot",
            se=0,
            effect=effect,
            seed=0,
        )


def test_export_rejects_missing_or_duplicate_trials(monkeypatch):
    monkeypatch.syspath_prepend(str(TOOLS))
    spec_export = importlib.util.spec_from_file_location(
        "summarize_shift_validation", TOOLS / "summarize_shift_validation.py"
    )
    exporter = importlib.util.module_from_spec(spec_export)
    spec_export.loader.exec_module(exporter)
    manifest = {
        "options": {
            "tips": [8],
            "standard_errors": [0],
            "effects": [2],
            "replicates": 1,
            "bootstrap": 0,
        }
    }
    rows = [
        {"case": f"c{i:03d}-r0000", "mode": mode, "status": "failed"}
        for i in range(8)
        for mode in ("shift", "convergence")
    ]
    exporter.validate_records(manifest, rows)
    with pytest.raises(ValueError, match="incomplete or duplicated"):
        exporter.validate_records(manifest, rows[:-1])
    with pytest.raises(ValueError, match="incomplete or duplicated"):
        exporter.validate_records(manifest, rows[:-1] + rows[:1])


@pytest.fixture
def audited_run(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(TOOLS))
    import summarize_shift_validation as exporter

    def fail(args):
        raise ValueError("intentional refit failure")

    monkeypatch.setattr(
        runner.parser, "parse_args", lambda args: SimpleNamespace(handler=fail)
    )
    options = SimpleNamespace(
        output=tmp_path / "run",
        tips=[8],
        standard_errors=[0.0],
        effects=[2.0],
        replicates=1,
        bootstrap=0,
        seed=1,
        rscript="unused",
        criterion="BIC",
    )
    runner.run(options)
    return exporter, options.output


def rewrite_json(path, mutate):
    import json

    value = json.loads(path.read_text())
    mutate(value)
    path.write_text(json.dumps(value))


@pytest.mark.parametrize(
    "artifact", ["summary", "truth", "traits", "record", "result", "error"]
)
def test_export_rejects_inconsistent_evidence(audited_run, tmp_path, artifact):
    import json

    exporter, run = audited_run
    case = run / "c000-r0000"
    if artifact == "summary":
        rewrite_json(run / "summary.json", lambda rows: rows[0].update(attempted=999))
    elif artifact == "truth":
        rewrite_json(case / "truth.json", lambda truth: truth.update(seed=99))
    elif artifact == "traits":
        path = case / "traits.tsv"
        path.write_text(path.read_text().replace("t0\t", "wrong\t"))
    elif artifact == "result":
        rewrite_json(case / "shift" / "result.json", lambda row: row.update(seed=99))
    else:
        path = run / "records.jsonl"
        rows = [json.loads(line) for line in path.read_text().splitlines()]
        rows[0].update(seed=99) if artifact == "record" else rows[0].update(error="")
        path.write_text("\n".join(json.dumps(row) for row in rows) + "\n")
        if artifact == "error":
            (case / "shift" / "result.json").write_text(json.dumps(rows[0]))
    with pytest.raises(ValueError):
        exporter.export([run], tmp_path / "evidence")
    assert not (tmp_path / "evidence").exists()


def test_export_keeps_failed_trials_and_publishes_atomically(
    audited_run, tmp_path, monkeypatch
):
    import json

    exporter, run = audited_run
    output = tmp_path / "evidence"
    original = exporter.write_json

    def fail_write(path, value):
        if path.name == "records.json":
            raise OSError("simulated disk failure")
        original(path, value)

    monkeypatch.setattr(exporter, "write_json", fail_write)
    with pytest.raises(OSError, match="disk failure"):
        exporter.export([run], output)
    assert not output.exists()
    monkeypatch.setattr(exporter, "write_json", original)
    exporter.export([run], output)
    rows = json.loads((output / "records.json").read_text())
    assert len(rows) == 16 and all(row["status"] == "failed" for row in rows)
    with pytest.raises(FileExistsError):
        exporter.export([run], output)
    assert len(json.loads((output / "records.json").read_text())) == 16


@pytest.mark.parametrize(
    "corruption", ["cycle", "tips", "observed", "root", "search", "bootstrap", "metric"]
)
def test_export_checks_completed_models(audited_run, tmp_path, corruption):
    import json

    exporter, run = audited_run
    directory = run / "c000-r0000"
    truth = json.loads((directory / "truth.json").read_text())
    tree = read_tree((directory / "tree.nwk").read_text(), "auto", True, quiet=True)
    ids = assign_branch_ids(tree)
    observations = dict(zip(truth["tip_names"], truth["observations"], strict=True))
    model = {
        "criterion": "BIC",
        "root_model": "OUfixedRoot",
        "convergence_searched": False,
        "max_shifts": 2,
        "search_strategy": "exhaustive",
        "trait": "value",
        "standard_error_column": "se",
        "branches": [
            {"branch_id": ids[n], "parent": -1 if n.is_root else ids[n.up]} for n in ids
        ],
        "tip_predictions": [
            {
                "leaf_name": n.name,
                "branch_id": ids[n],
                "predicted": 0,
                "observed": observations[n.name],
                "standard_error": 0,
                "regime": "baseline",
            }
            for n in tree.leaves()
        ],
        "shift_branch_ids": [],
        "parameters": {"alpha": 0.7},
        "convergence": None,
        "backend_version": "test",
        "bootstrap": None,
    }
    path = run / "records.jsonl"
    rows = [json.loads(line) for line in path.read_text().splitlines()]
    row = rows[0]
    row.pop("error")
    row.update(status="completed", **cases.selection_metrics(model, truth))
    if corruption == "metric":
        row["tip_mean_rmse"] = 99
    elif corruption == "cycle":
        model["branches"][1]["parent"] = model["branches"][1]["branch_id"]
    elif corruption == "tips":
        model["tip_predictions"].pop()
    elif corruption == "observed":
        model["tip_predictions"][0]["observed"] += 1
    elif corruption == "root":
        model["root_model"] = "OUrandomRoot"
    elif corruption == "search":
        model["max_shifts"] = 3
    else:
        model["bootstrap"] = {"attempted": 3}
    path.write_text("\n".join(json.dumps(r) for r in rows) + "\n")
    (directory / "shift" / "result.json").write_text(json.dumps(row))
    (directory / "shift" / "model.json").write_text(json.dumps(model))
    with pytest.raises(ValueError):
        exporter.export([run], tmp_path / "evidence")
    assert not (tmp_path / "evidence").exists()


@pytest.mark.parametrize(
    "changes",
    [
        {"tips": [8, 8]},
        {"effects": [2, 2]},
        {"standard_errors": [0, 0]},
        {"bootstrap": 2147483648},
        {"replicates": 0},
    ],
)
def test_invalid_grid_is_rejected_before_creating_output(tmp_path, changes):
    options = SimpleNamespace(
        output=tmp_path / "run",
        tips=[8],
        effects=[2],
        standard_errors=[0],
        bootstrap=0,
        replicates=1,
        seed=1,
    )
    vars(options).update(changes)
    with pytest.raises(ValueError):
        runner.run(options)
    assert not options.output.exists()
