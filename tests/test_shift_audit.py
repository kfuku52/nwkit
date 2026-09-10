"""Regression checks for scale, output integrity and known observation errors."""

import json
import os

import numpy as np
import pandas as pd
import pytest

from nwkit.cli import main
from nwkit.shift_math import close_in_units, observation_variances, remaining_heights
from nwkit.util import read_tree
from tests import test_shift as shift_support
from tests.test_shift import fake_backend

shift_inputs = shift_support.shift_inputs


def test_checks_have_no_absolute_trait_unit_floor():
    assert not close_in_units(1e-30, 2e-30)
    assert close_in_units(1e-30, 1e-30 * (1 + 1e-9))
    assert close_in_units(1e16 - (1e16 + 2), -3, operands=(1e16, 1e16 + 2))


def test_short_terminal_branch_geometry_avoids_root_depth_subtraction():
    tree = read_tree(
        "((A:0.0000000000000001,B:0.0000000000000001):1,C:1);", "auto", True, quiet=True
    )
    clade = next(n for n in tree.traverse() if set(n.leaf_names()) == {"A", "B"})
    assert remaining_heights(tree)[clade] == 1e-16


@pytest.mark.parametrize("error", [-1, float("inf"), float("nan"), 1e200, 1e-200])
def test_unrepresentable_or_invalid_error_variance_is_rejected(error):
    with pytest.raises(ValueError):
        observation_variances([error])


def test_se_alignment(shift_inputs, tmp_path, monkeypatch):
    frame = pd.read_csv(tmp_path / "trait.tsv", sep="\t")
    frame["se"] = [0.4, 0.2, 0.1, 0]
    frame.to_csv(tmp_path / "trait.tsv", sep="\t", index=False)

    def backend(directory, args):
        data = pd.read_csv(directory / "trait.tsv", sep="\t")
        assert data.standard_error.tolist() == pytest.approx([0.1, 0.2, 0, 0.4])
        assert data.observation_variance.tolist() == pytest.approx(
            [0.01, 0.04, 0, 0.16]
        )
        return fake_backend(directory, args)

    monkeypatch.setattr("nwkit.shift.run_backend", backend)
    main([*shift_inputs, "--standard-error-column", "se"])
    saved = json.loads((tmp_path / "model.json").read_text())
    assert saved["observation_error"] == "known_independent_variances"
    assert [r["standard_error"] for r in saved["tip_predictions"]] == [0.1, 0.2, 0, 0.4]


@pytest.mark.parametrize("error", ["NA", "-1", "inf", "1e200", "1e-200"])
def test_invalid_se_fails_before_backend(shift_inputs, tmp_path, monkeypatch, error):
    table = pd.read_csv(tmp_path / "trait.tsv", sep="\t")
    table["se"] = error
    table.to_csv(tmp_path / "trait.tsv", sep="\t", index=False)

    def fail(*args):
        pytest.fail("backend must not run")

    monkeypatch.setattr("nwkit.shift.run_backend", fail)
    with pytest.raises(ValueError):
        main([*shift_inputs, "--standard-error-column", "se"])


@pytest.mark.parametrize(
    "field,value",
    [
        ("coverage", 1.5),
        ("ensemble_failed", -1),
        ("ensemble_attempted", 2),
        ("evaluated_configurations", 1.5),
        ("alpha_upper", -1),
    ],
)
def test_invalid_search_diagnostics_fail_closed(
    shift_inputs, tmp_path, monkeypatch, field, value
):
    def corrupt(directory, args):
        result = fake_backend(directory, args)
        table = pd.read_csv(directory / "search.tsv", sep="\t")
        table[field] = value
        table.to_csv(directory / "search.tsv", sep="\t", index=False)
        return result

    monkeypatch.setattr("nwkit.shift.run_backend", corrupt)
    with pytest.raises(ValueError):
        main(shift_inputs)
    assert not (tmp_path / "model.json").exists()


def test_failed_output_install_restores_all_outputs(
    shift_inputs, tmp_path, monkeypatch
):
    import nwkit.output_transaction as transaction

    monkeypatch.setattr("nwkit.shift.run_backend", fake_backend)
    model, effects = tmp_path / "model.json", tmp_path / "effects.tsv"
    model.write_text("old model")
    effects.write_text("old effects")
    original = transaction.replace_output

    def replace(source, target):
        if str(target) == str(effects):
            raise OSError("injected install failure")
        return original(source, target)

    monkeypatch.setattr(transaction, "replace_output", replace)
    with pytest.raises(OSError, match="injected"):
        main([*shift_inputs, "--effects-out", str(effects)])
    assert model.read_text() == "old model"
    assert effects.read_text() == "old effects"


@pytest.mark.integration
@pytest.mark.skipif(
    not os.environ.get("NWKIT_TEST_RSCRIPT"), reason="requires R and kfl1ou"
)
@pytest.mark.parametrize("mode", ["none", "zero", "known"])
def test_real_se_search(tmp_path, mode):
    from tests.test_shift_reference import TREE, VALUES

    (tmp_path / "tree.nwk").write_text(TREE)
    frame = pd.DataFrame(
        {
            "leaf_name": [f"t{i}" for i in range(8)],
            "x": VALUES,
            "se": np.linspace(0.01, 0.15, 8) if mode == "known" else np.zeros(8),
        }
    )
    # Reverse order to exercise the named R input-error matrix.
    frame.iloc[::-1].to_csv(tmp_path / "traits.tsv", sep="\t", index=False)
    args = [
        "shift",
        "--selection",
        "ic",
        "-i",
        str(tmp_path / "tree.nwk"),
        "--trait",
        str(tmp_path / "traits.tsv"),
        "--state-column",
        "x",
        "--max-shifts",
        "1",
        "--criterion",
        "BIC",
        "--model-out",
        str(tmp_path / "model.json"),
        "--rscript",
        os.environ["NWKIT_TEST_RSCRIPT"],
    ]
    main(args + (["--standard-error-column", "se"] if mode != "none" else []))
    model = json.loads((tmp_path / "model.json").read_text())
    assert np.isfinite(model["parameters"]["log_likelihood"])
    assert [r["standard_error"] for r in model["tip_predictions"]] == pytest.approx(
        frame.se
    )
    if mode == "zero":
        main(args)
        no_error = json.loads((tmp_path / "model.json").read_text())
        assert model["parameters"] == pytest.approx(
            no_error["parameters"], rel=1e-7, abs=1e-8
        )
        assert model["shift_branch_ids"] == no_error["shift_branch_ids"]


def test_empty_se_column_is_not_silently_ignored(shift_inputs):
    with pytest.raises(ValueError, match="nonempty"):
        main([*shift_inputs, "--standard-error-column", ""])
