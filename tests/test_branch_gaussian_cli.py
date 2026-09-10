"""Model-table validation and real CLI round trips through Gaussian consumers."""

import json
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from nwkit.branch_gaussian import build_branch_gaussian_process
from nwkit.branch_gaussian_input import load_branch_gaussian_assignment
from nwkit.cli import main
from nwkit.gaussian_inference import condition_gaussian_tree, simulate_gaussian_process
from nwkit.gaussian_tree import GaussianRootPrior
from nwkit.util import assign_branch_ids, read_tree

TREE = "((A:0.3,B:0.7)I:0.4,C:1.1)R;"
HEADER = "branch_id\tmodel\tsigma2\talpha\ttheta\tjump_mean\tjump_variance\n"
MODELS = (
    HEADER
    + "1\tBM\t0.8\t\t\t\t\n2\tBM\t0.7\t\t\t\t\n3\tOU\t1.3\t0.6\t0.9\t0.2\t0.15\n4\tJUMP\t\t\t\t-0.4\t0.5\n"
)
ROOT = ["--root-prior", "gaussian", "--root-mean", "0.3", "--root-variance", "0.9"]


@pytest.fixture
def inputs(tmp_path):
    paths = {name: tmp_path / name for name in ("tree.nwk", "models.tsv", "traits.tsv")}
    paths["tree.nwk"].write_text(TREE)
    paths["models.tsv"].write_text(MODELS)
    paths["traits.tsv"].write_text(
        "leaf_name\tx\tse\nA\t1.2\t0.2\nB\t-0.4\t0.1\nC\t2.1\t0.4\n"
    )
    return paths


def command(inputs, *, observe=True, root=ROOT):
    arguments = [
        "asr",
        "--model",
        "BRANCH-GAUSSIAN",
        "-i",
        str(inputs["tree.nwk"]),
        "--input-rooted",
        "yes",
        "--branch-models",
        str(inputs["models.tsv"]),
        *root,
    ]
    if observe:
        arguments += [
            "--trait",
            str(inputs["traits.tsv"]),
            "--state-column",
            "x",
            "--standard-error-column",
            "se",
        ]
    return arguments


def process_from(inputs):
    tree = read_tree(str(inputs["tree.nwk"]), "1", True, quiet=True, rooted="yes")
    assignment = load_branch_gaussian_assignment(
        tree, branch_models=inputs["models.tsv"]
    )
    return build_branch_gaussian_process(
        tree,
        assignment.models_by_branch_id,
        root=GaussianRootPrior("gaussian", 0.3, 0.9),
    ), assignment


def test_cli_asr_likelihood_and_model_roundtrip(inputs, tmp_path):
    out, model_json, model_tsv = [
        tmp_path / name for name in ("asr.tsv", "model.json", "normalized.tsv")
    ]
    assert (
        main(
            [
                *command(inputs),
                "-o",
                str(out),
                "--process-out",
                str(model_json),
                "--branch-models-out",
                str(model_tsv),
            ]
        )
        is None
    )
    process, assignment = process_from(inputs)
    result = condition_gaussian_tree(
        process,
        {"A": 1.2, "B": -0.4, "C": 2.1},
        standard_errors={"A": 0.2, "B": 0.1, "C": 0.4},
    )
    expected = {
        identifier: result.marginals[node]
        for node, identifier in assign_branch_ids(process.tree).items()
    }
    frame = pd.read_csv(out, sep="\t")
    for row in frame.itertuples():
        assert row.mean == pytest.approx(expected[row.branch_id].mean)
        assert row.variance == pytest.approx(expected[row.branch_id].variance)
    saved = json.loads(model_json.read_text())
    assert saved["log_likelihood"] == pytest.approx(result.log_likelihood)
    assert saved["root"] == {"mode": "gaussian", "mean": 0.3, "variance": 0.9}
    assert saved["jump_position"] == "branch_end_after_diffusion"
    assert len(saved["nodes"]) == 5
    assert len(saved["transitions"]) == 4
    restored = load_branch_gaussian_assignment(process.tree, branch_models=model_tsv)
    assert restored.models_by_branch_id == assignment.models_by_branch_id
    likelihood = tmp_path / "likelihood.tsv"
    assert (
        main([*command(inputs), "--output", "likelihood", "-o", str(likelihood)])
        is None
    )
    row = pd.read_csv(likelihood, sep="\t").iloc[0]
    assert row.log_likelihood == pytest.approx(result.log_likelihood)
    assert row.num_observed == 3
    args = command(inputs)
    args[args.index(str(inputs["models.tsv"]))] = str(model_tsv)
    repeated = tmp_path / "repeat.tsv"
    main([*args, "-o", str(repeated)])
    assert repeated.read_bytes() == out.read_bytes()


def test_regime_tables_are_equivalent_and_preserve_names(inputs, tmp_path):
    regimes, definitions = tmp_path / "regimes.tsv", tmp_path / "definitions.tsv"
    regimes.write_text(
        "branch_id\tregime\n1\t背景\n2\tother\n3\tselected\n4\tjump\n",
        encoding="utf-8-sig",
    )
    definitions.write_text(
        "regime\tmodel\tsigma2\talpha\ttheta\tjump_mean\tjump_variance\n背景\tBM\t0.8\t\t\t\t\nother\tBM\t0.7\t\t\t\t\nselected\tOU\t1.3\t0.6\t0.9\t0.2\t0.15\njump\tJUMP\t\t\t\t-0.4\t0.5\n"
    )
    args = command(inputs)
    position = args.index("--branch-models")
    args[position : position + 2] = [
        "--branch-regimes",
        str(regimes),
        "--regime-models",
        str(definitions),
    ]
    out, saved, audit = (
        tmp_path / "out.tsv",
        tmp_path / "saved.json",
        tmp_path / "audit.jsonl",
    )
    exported = tmp_path / "exported.tsv"
    main(
        [
            *args,
            "-o",
            str(out),
            "--process-out",
            str(saved),
            "--branch-models-out",
            str(exported),
            "--audit",
            str(audit),
        ]
    )
    metadata = json.loads(saved.read_text())
    assert metadata["branch_regimes"][0] == {"branch_id": 1, "regime": "背景"}
    process, direct = process_from(inputs)
    parsed = load_branch_gaussian_assignment(
        process.tree, branch_regimes=regimes, regime_models=definitions
    )
    assert parsed.models_by_branch_id == direct.models_by_branch_id
    record = json.loads(audit.read_text())
    assert record["status"] == "ok"
    recorded = {item["path"] for item in record["inputs"]}
    assert {
        str(path.resolve())
        for path in (regimes, definitions, inputs["tree.nwk"], inputs["traits.tsv"])
    } == recorded
    assert str(exported.resolve()) in {item["path"] for item in record["outputs"]}


def test_prior_simulation_matches_api_and_seed(inputs, tmp_path):
    outfile = tmp_path / "draws.tsv"
    main(
        [
            *command(inputs, observe=False),
            "--output",
            "prior-samples",
            "--prior-samples",
            "7",
            "--seed",
            "42",
            "-o",
            str(outfile),
        ]
    )
    process, _ = process_from(inputs)
    samples = simulate_gaussian_process(process, num_samples=7, seed=42)
    table = pd.read_csv(outfile, sep="\t")
    np.testing.assert_allclose(table.value.to_numpy().reshape(7, 5), samples.values)
    assert table.simulation.tolist() == [
        identifier for identifier in range(1, 8) for _ in range(5)
    ]
    ids = assign_branch_ids(process.tree)
    assert table.branch_id.tolist() == [ids[node] for node in samples.nodes] * 7


def test_flat_root_simulation_and_missing_observations(inputs, tmp_path):
    out = tmp_path / "out.tsv"
    args = command(inputs, observe=False, root=["--root-prior", "flat"])
    with pytest.raises(ValueError, match="root-value"):
        main([*args, "--output", "prior-samples", "-o", str(out)])
    main(
        [
            *args,
            "--output",
            "prior-samples",
            "--prior-root-value",
            "2",
            "--prior-samples",
            "2",
            "-o",
            str(out),
        ]
    )
    frame = pd.read_csv(out, sep="\t")
    assert frame.loc[frame.branch_id == 0, "value"].tolist() == [2, 2]
    inputs["traits.tsv"].write_text("leaf_name\tx\tse\nA\t1.2\t0.2\nC\t2.1\t0.4\n")
    main([*command(inputs), "-o", str(out)])
    frame = pd.read_csv(out, sep="\t").set_index("name")
    assert frame.loc["B", "is_imputed"]
    assert np.isfinite(frame.loc["B", "mean"])


@pytest.mark.parametrize(
    "source,match",
    [
        ("", "Empty"),
        ("\n", "Empty"),
        ("branch_id\tmodel\tmodel\n1\tBM\tBM\n", "Duplicate TSV"),
        ("branch_id\tmodel\tdrift\n1\tBM\t0\n", "columns"),
        ("branch_id\tmodel\tsigma2\n1\tBM\n", "width"),
        ("branch_id\tmodel\tsigma2\n\n", "Empty TSV row"),
        ("branch_id\tmodel\tsigma2\n1\tBM\t1\n1\tBM\t1\n", "Duplicate branch_id"),
        ("branch_id\tmodel\tsigma2\n0\tBM\t1\n", "root 0"),
        ("branch_id\tmodel\tsigma2\n1.0\tBM\t1\n", "branch_id"),
        ("branch_id\tmodel\tsigma2\n01\tBM\t1\n", "branch_id"),
        ("branch_id\tmodel\tsigma2\nTrue\tBM\t1\n", "branch_id"),
        ("branch_id\tmodel\tsigma2\n1\tBM\t1\n", "cover every"),
        ("branch_id\tmodel\tsigma2\n1\tBM\t-1\n", "non-negative"),
        ("branch_id\tmodel\tsigma2\n1\tBM\tNaN\n", "Non-finite"),
        ("branch_id\tmodel\tsigma2\n1\tBM\tx\n", "numeric"),
        ("branch_id\tmodel\tsigma2\n1\tBM\t\n", "Missing sigma2"),
        ("branch_id\tmodel\tsigma2\n1\tbad\t1\n", "model must"),
        (HEADER + "1\tBM\t1\t1\t\t\t\n", "not used"),
        (HEADER + "1\tJUMP\t1\t\t\t0\t1\n", "not used"),
        (HEADER + "1\tOU\t1\t1\t\t\t\n", "Missing theta"),
        (HEADER + "1\tBM\t1\t\t\t0\t\n", "Missing jump_variance"),
        (HEADER + "1\tBM\t1\t\t\t\t1\n", "Missing jump_mean"),
        (HEADER + "1\tJUMP\t\t\t\t0\t-1\n", "non-negative"),
    ],
)
def test_invalid_model_tsv_leaves_outputs_untouched(inputs, tmp_path, source, match):
    inputs["models.tsv"].write_text(source)
    out, model = tmp_path / "out.tsv", tmp_path / "model.json"
    out.write_text("original table")
    model.write_text("original model")
    with pytest.raises(ValueError, match=match):
        main([*command(inputs), "-o", str(out), "--process-out", str(model)])
    assert out.read_text() == "original table"
    assert model.read_text() == "original model"


@pytest.mark.parametrize(
    "extra",
    [
        ["--ci-level", "1"],
        ["--ci-level", "nan"],
        ["--prior-samples", "2"],
        ["--prior-root-value", "1"],
        ["--root-prior", "fixed"],
        ["--root-prior", "flat"],
        ["--root-variance", "0"],
        ["--root-mean", "nan"],
        ["--state-column", "x,se"],
        ["--standard-error-column", "x"],
        ["--process-out", "-"],
        ["--branch-models-out", "-"],
        ["--output", "prior-samples"],
        ["--regime-models", "unused.tsv"],
    ],
)
def test_incompatible_options_are_rejected(inputs, tmp_path, extra):
    with pytest.raises((ValueError, SystemExit)):
        main([*command(inputs), "-o", str(tmp_path / "out.tsv"), *extra])


@pytest.mark.parametrize(
    "extra", [["--prior-samples", "0"], ["--prior-samples", "10001"], ["--seed", "-1"]]
)
def test_invalid_simulation_options(inputs, extra):
    with pytest.raises(ValueError):
        main([*command(inputs, observe=False), "--output", "prior-samples", *extra])


def test_stdin_and_stdout_roundtrip(inputs, monkeypatch, capsys):
    args = command(inputs)
    args[args.index(str(inputs["models.tsv"]))] = "-"
    monkeypatch.setattr("sys.stdin", StringIO("\ufeff" + MODELS))
    main(args)
    output = capsys.readouterr().out
    assert len(pd.read_csv(StringIO(output), sep="\t")) == 5
    args[args.index(str(inputs["tree.nwk"]))] = "-"
    with pytest.raises(ValueError, match="STDIN"):
        main(args)


@pytest.mark.parametrize("input_name", ["tree.nwk", "models.tsv", "traits.tsv"])
def test_outputs_cannot_replace_inputs(inputs, input_name):
    original = inputs[input_name].read_bytes()
    with pytest.raises(ValueError):
        main([*command(inputs), "-o", str(inputs[input_name])])
    assert inputs[input_name].read_bytes() == original


def test_output_aliases_and_failed_staging_preserve_existing_files(
    inputs, tmp_path, monkeypatch
):
    out, saved = tmp_path / "out.tsv", tmp_path / "model.json"
    out.write_text("old table")
    saved.write_text("old model")
    alias = tmp_path / "alias.tsv"
    alias.symlink_to(out)
    with pytest.raises(ValueError):
        main([*command(inputs), "-o", str(out), "--branch-models-out", str(alias)])
    original = Path.write_text

    def fail_model(self, *args, **kwargs):
        if "model.json.stage" in self.name:
            raise OSError("injected staging failure")
        return original(self, *args, **kwargs)

    monkeypatch.setattr(Path, "write_text", fail_model)
    with pytest.raises(OSError, match="staging failure"):
        main([*command(inputs), "-o", str(out), "--process-out", str(saved)])
    assert out.read_text() == "old table"
    assert saved.read_text() == "old model"


def test_unknown_duplicate_and_unused_regimes(inputs, tmp_path):
    process, _ = process_from(inputs)
    mapping, definitions = tmp_path / "mapping.tsv", tmp_path / "definitions.tsv"
    definitions.write_text("regime\tmodel\tsigma2\nx\tBM\t1\n")
    mapping.write_text("branch_id\tregime\n1\tx\n2\tx\n3\tx\n4\tunknown\n")
    with pytest.raises(ValueError, match="exactly match"):
        load_branch_gaussian_assignment(
            process.tree, branch_regimes=mapping, regime_models=definitions
        )
    mapping.write_text("branch_id\tregime\n1\tx\n1\tx\n")
    with pytest.raises(ValueError, match="Duplicate branch_id"):
        load_branch_gaussian_assignment(
            process.tree, branch_regimes=mapping, regime_models=definitions
        )
    mapping.write_text("branch_id\tregime\n1\tx\n2\tx\n3\tx\n4\tx\n")
    definitions.write_text("regime\tmodel\tsigma2\nx\tBM\t1\ny\tBM\t1\n")
    with pytest.raises(ValueError, match="exactly match"):
        load_branch_gaussian_assignment(
            process.tree, branch_regimes=mapping, regime_models=definitions
        )


@pytest.mark.parametrize(
    "root",
    [
        ["--root-prior", "fixed"],
        ["--root-prior", "gaussian", "--root-mean", "0"],
        ["--root-prior", "stationary", "--root-mean", "0"],
    ],
)
def test_proper_root_parameters_are_required(inputs, root):
    with pytest.raises(ValueError):
        main(command(inputs, root=root))


@pytest.mark.parametrize(
    "data",
    [
        "leaf_name\tx\tse\nA\t1\t\n",
        "leaf_name\tx\tse\nA\t1\t-1\n",
        "leaf_name\tx\tse\nA\tNA\t0\nB\tNA\t0\nC\tNA\t0\n",
    ],
)
def test_missing_or_invalid_observation_information(inputs, data):
    inputs["traits.tsv"].write_text(data)
    with pytest.raises(ValueError):
        main(command(inputs))


def test_simulation_budget_is_checked_before_drawing(inputs, monkeypatch):
    monkeypatch.setattr("nwkit.branch_gaussian_output._MAX_SIMULATION_VALUES", 4)
    with pytest.raises(ValueError, match="node-value limit"):
        main(
            [
                *command(inputs, observe=False),
                "--output",
                "prior-samples",
                "--prior-samples",
                "1",
            ]
        )


def test_loader_does_not_consume_stdin_twice(inputs, monkeypatch):
    process, _ = process_from(inputs)
    stream = StringIO("not read")
    monkeypatch.setattr("sys.stdin", stream)
    with pytest.raises(ValueError, match="STDIN"):
        load_branch_gaussian_assignment(
            process.tree, branch_regimes="-", regime_models="-"
        )
    assert stream.tell() == 0


def test_checked_in_example_matches_documented_python_result(tmp_path):
    example = Path(__file__).resolve().parents[1] / "examples" / "branch_gaussian"
    out = tmp_path / "likelihood.tsv"
    main(
        [
            "asr",
            "--model",
            "BRANCH-GAUSSIAN",
            "-i",
            str(example / "tree.nwk"),
            "--input-rooted",
            "yes",
            "--branch-regimes",
            str(example / "branch_regimes.tsv"),
            "--regime-models",
            str(example / "regime_models.tsv"),
            *ROOT,
            "--output",
            "likelihood",
            "--trait",
            str(example / "traits.tsv"),
            "--state-column",
            "x",
            "--standard-error-column",
            "se",
            "-o",
            str(out),
        ]
    )
    assert pd.read_csv(out, sep="\t").iloc[0].log_likelihood == pytest.approx(
        -5.487990618, abs=5e-10
    )
