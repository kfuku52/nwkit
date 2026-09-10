import json

import numpy as np
import pandas as pd
import pytest

from nwkit.cli import main
from nwkit.shift_native_model import ShiftLayout, ShiftTree
from nwkit.util import read_tree
from tests.test_shift_native_model import TREE


@pytest.fixture
def native_inputs(tmp_path):
    tree = tmp_path / "tree.nwk"
    tree.write_text(TREE)
    data = pd.DataFrame(
        np.random.default_rng(1852).normal(size=(8, 2)), columns=["x", "y"]
    )
    data.insert(0, "leaf_name", list("abcdefgh"))
    data.loc[2, "y"] = np.nan
    data["x_se"], data["y_se"] = 0.1, 0.2
    data.to_csv(tmp_path / "traits.tsv", sep="\t", index=False, na_rep="NA")
    parsed = ShiftTree.build(read_tree(TREE, "auto", True, quiet=True))
    layout = ShiftLayout.build(parsed, [1, 7], [[0, 7], [1]])
    groups = layout.node_groups(parsed)
    pd.DataFrame(
        {
            "branch_id": parsed.branch_ids,
            "regime": ["background" if x == 0 else "new regime" for x in groups],
        }
    ).to_csv(tmp_path / "map.tsv", sep="\t", index=False)
    return [
        "shift",
        "--selection",
        "native",
        "-i",
        str(tree),
        "--trait",
        str(tmp_path / "traits.tsv"),
        "--state-column",
        "x,y",
        "--standard-error-column",
        "x_se,y_se",
        "--regime-map",
        str(tmp_path / "map.tsv"),
        "--model-out",
        str(tmp_path / "model.json"),
        "--alpha",
        "0.2,0.4",
        "--process-tip-variance",
        "1,2",
        "--effects-out",
        str(tmp_path / "effects.tsv"),
        "--regime-parameters-out",
        str(tmp_path / "regimes.tsv"),
        "--tip-summary-out",
        str(tmp_path / "tips.tsv"),
        "-o",
        str(tmp_path / "output-map.tsv"),
    ]


def test_native_cli_never_calls_r_and_exports_original_regime_names(
    native_inputs, tmp_path, monkeypatch
):
    def forbidden(*args, **kwargs):
        raise AssertionError("Native inference attempted to call kfl1ou")

    monkeypatch.setattr("nwkit.shift.run_backend", forbidden)
    main(native_inputs)
    model = json.loads((tmp_path / "model.json").read_text())
    assert model["backend"] == "nwkit"
    assert model["inference_role"] == "fixed_layout_parameter_estimation"
    assert model["trait_names"] == ["x", "y"]
    assert model["selection_calibration"] is None
    assert model["shift_branch_ids"] == [1, 7]
    assert [r["num_observations"] for r in model["traits"]] == [8, 7]
    assert {row["regime"] for row in model["branches"]} == {"background", "new regime"}
    assert model["traits"][0]["alpha"] == pytest.approx(0.2)
    tips = pd.read_csv(tmp_path / "tips.tsv", sep="\t")
    assert len(tips) == 16 and tips.observed.isna().sum() == 1
    assert tips.predicted.notna().all()
    assert list(pd.read_csv(tmp_path / "output-map.tsv", sep="\t").columns) == [
        "branch_id",
        "regime",
    ]


def test_native_cli_missing_observation_error_and_output_collision_preserve_files(
    native_inputs, tmp_path
):
    target = tmp_path / "model.json"
    target.write_text("previous output")
    table = pd.read_csv(tmp_path / "traits.tsv", sep="\t")
    table.loc[0, "x_se"] = np.nan
    table.to_csv(tmp_path / "traits.tsv", sep="\t", index=False, na_rep="NA")
    with pytest.raises(ValueError, match="standard error"):
        main(native_inputs)
    assert target.read_text() == "previous output"
    native_inputs[native_inputs.index("--model-out") + 1] = str(tmp_path / "traits.tsv")
    original = (tmp_path / "traits.tsv").read_bytes()
    with pytest.raises(ValueError, match="input|same|replace"):
        main(native_inputs)
    assert (tmp_path / "traits.tsv").read_bytes() == original


def test_native_options_cannot_be_silently_ignored_by_other_backends(native_inputs):
    native_inputs[native_inputs.index("--selection") + 1] = "ic"
    with pytest.raises(ValueError, match="require --selection native"):
        main(native_inputs)


def test_native_alpha_limit_is_not_serialized_as_infinity(native_inputs, tmp_path):
    native_inputs[native_inputs.index("--alpha") + 1] = "inf"
    main(native_inputs)
    model = json.loads((tmp_path / "model.json").read_text())
    assert all(row["alpha"] is None for row in model["traits"])
    assert all(row["optimum"] is None for row in model["regime_parameters"])
    assert "Infinity" not in (tmp_path / "model.json").read_text()


def test_native_search_cli_and_resume_reject_changed_input(
    native_inputs, tmp_path, monkeypatch
):
    index = native_inputs.index("--regime-map")
    del native_inputs[index : index + 2]
    native_inputs.extend(["--max-shifts", "1", "--calibration-replicates", "19"])
    main(native_inputs)
    model = json.loads((tmp_path / "model.json").read_text())
    assert model["inference_role"] == "research_joint_shift_selection"
    assert model["selection_calibration"]["full_search_repeated"]
    # A terminal shift affecting only the missing trait coordinate is unobservable.
    assert len(model["candidates"]) == 14
    assert model["completion_status"] == "complete"
    assert model["search"]["complete_discrete_enumeration"]
    native_inputs[native_inputs.index("--model-out") + 1] = str(
        tmp_path / "resumed.json"
    )
    native_inputs.extend(["--resume-model", str(tmp_path / "model.json")])

    def forbidden(*args, **kwargs):
        raise AssertionError("Resume unexpectedly reran inference")

    monkeypatch.setattr("nwkit.shift_native_selection.select_native", forbidden)
    main(native_inputs)
    assert json.loads((tmp_path / "resumed.json").read_text()) == model
    table = pd.read_csv(tmp_path / "traits.tsv", sep="\t")
    table.loc[0, "x"] += 0.01
    table.to_csv(tmp_path / "traits.tsv", sep="\t", index=False)
    with pytest.raises(ValueError, match="analysis_input_sha256 differs"):
        main(native_inputs)
    assert json.loads((tmp_path / "resumed.json").read_text()) == model


def test_native_resume_rejects_changed_configuration(native_inputs, tmp_path):
    main(native_inputs)
    native_inputs[native_inputs.index("--model-out") + 1] = str(
        tmp_path / "resumed.json"
    )
    native_inputs.extend(["--resume-model", str(tmp_path / "model.json")])
    native_inputs[native_inputs.index("--alpha") + 1] = "0.8"
    with pytest.raises(ValueError, match="configuration_sha256 differs"):
        main(native_inputs)
    assert not (tmp_path / "resumed.json").exists()
