"""Input preservation and complete publication for tree-ensemble runs."""

import os
from contextlib import redirect_stdout

import pandas as pd
import pytest

from nwkit.cli import main

TREE = "((A:1,B:1):1,(C:1,D:1):1);"


def arguments(tmp_path):
    traits = tmp_path / "traits.tsv"
    traits.write_text("leaf_name\tstate\nA\t1\nB\t2\nC\t4\nD\t5\n")
    return [
        "asr",
        "-i",
        TREE,
        "--input-rooted",
        "yes",
        "--trait",
        str(traits),
        "--state-column",
        "state",
        "--trait-type",
        "continuous",
    ]


@pytest.mark.parametrize("option", ["-o", "--model-out", "--tree-out"])
@pytest.mark.parametrize(
    "alias",
    [
        "same",
        pytest.param(
            "symlink",
            marks=pytest.mark.skipif(
                os.name == "nt", reason="symlink creation may need privileges"
            ),
        ),
        "hardlink",
    ],
)
def test_regime_input_cannot_be_replaced(tmp_path, option, alias):
    args = arguments(tmp_path)
    regimes = tmp_path / "regimes.tsv"
    regimes.write_text("branch_id\tregime\n" + "".join(f"{i}\tr\n" for i in range(7)))
    args += ["--regime-map", str(regimes)]
    parameters = tmp_path / "parameters.tsv"
    original = "regime\tsigma2\nr\t1\n"
    parameters.write_text(original)
    destination = parameters
    if alias != "same":
        destination = tmp_path / "alias.tsv"
        if alias == "symlink":
            destination.symlink_to(parameters)
        else:
            os.link(parameters, destination)
    with pytest.raises(ValueError, match="regime-parameters"):
        main(
            args
            + [
                "--model",
                "BMS",
                "--regime-parameters",
                str(parameters),
                option,
                str(destination),
            ]
        )
    assert parameters.read_text() == original


@pytest.mark.parametrize("stdout", [False, True])
@pytest.mark.parametrize("failure", ["tips", "late"])
@pytest.mark.parametrize("existing", [False, True])
def test_ensemble_failure_preserves_every_output(
    tmp_path, monkeypatch, capsys, stdout, failure, existing
):
    args = arguments(tmp_path)
    trees = tmp_path / "trees.nwk"
    trees.write_text(TREE if failure == "late" else TREE.replace("D:", "E:"))
    paths = {
        name: tmp_path / name for name in ["reference.tsv", "ensemble.tsv", "model.tsv"]
    }
    if existing:
        for path in paths.values():
            path.write_text("previous " + path.name)
    if failure == "late":
        import nwkit.asr_tree_ensemble as module

        original = module.write_tree_ensemble

        def fail_after_fit(*args):
            original(*args)
            raise OSError("late ensemble failure")

        monkeypatch.setattr(module, "write_tree_ensemble", fail_after_fit)
    with pytest.raises((ValueError, OSError), match="tip set|late ensemble failure"):
        main(
            args
            + [
                "--model",
                "BM",
                "--sigma2",
                "1",
                "-o",
                "-" if stdout else str(paths["reference.tsv"]),
                "--model-out",
                str(paths["model.tsv"]),
                "--tree-ensemble",
                str(trees),
                "--tree-ensemble-out",
                str(paths["ensemble.tsv"]),
            ]
        )
    assert capsys.readouterr().out == ""
    for path in paths.values():
        if existing:
            assert path.read_text() == "previous " + path.name
        else:
            assert not path.exists()
    assert not list(tmp_path.glob(".*.stage.*"))


def test_ensemble_publishes_stdout_and_figure(tmp_path, capsys):
    args = arguments(tmp_path)
    trees = tmp_path / "trees.nwk"
    trees.write_text(TREE)
    ensemble = tmp_path / "ensemble.tsv"
    figure = tmp_path / "figure.svg"
    main(
        args
        + [
            "--model",
            "BM",
            "--sigma2",
            "1",
            "-o",
            "-",
            "--tree-ensemble",
            str(trees),
            "--tree-ensemble-out",
            str(ensemble),
            "--figure-out",
            str(figure),
        ]
    )
    assert "branch_id\t" in capsys.readouterr().out
    assert len(pd.read_csv(ensemble, sep="\t")) == 7
    assert "<svg" in figure.read_text()


@pytest.mark.parametrize("failure", ["write", "flush"])
def test_ensemble_stdout_failure_restores_previous_file(tmp_path, failure):
    args = arguments(tmp_path)
    trees = tmp_path / "trees.nwk"
    trees.write_text(TREE)
    ensemble = tmp_path / "ensemble.tsv"
    ensemble.write_text("previous ensemble")

    class BrokenOutput:
        def write(self, value):
            if failure == "write":
                raise BrokenPipeError("closed consumer")
            return len(value)

        def flush(self):
            raise BrokenPipeError("closed consumer")

    with (
        redirect_stdout(BrokenOutput()),
        pytest.raises(BrokenPipeError, match="closed consumer"),
    ):
        main(
            args
            + [
                "--model",
                "BM",
                "--sigma2",
                "1",
                "-o",
                "-",
                "--tree-ensemble",
                str(trees),
                "--tree-ensemble-out",
                str(ensemble),
            ]
        )
    assert ensemble.read_text() == "previous ensemble"
