import json

import numpy as np
import pandas as pd
import pytest

from nwkit.cli import main, parser
from nwkit.shift_native_model import ShiftLayout, ShiftTree
from nwkit.util import read_tree
from tests.test_shift_native_model import TREE


@pytest.mark.parametrize("selection", ["native", "calibrated", "ic"])
def test_search_budget_defaults_agree_across_cli_api_and_legacy_dispatch(
    selection, monkeypatch
):
    from nwkit.shift_cli import _command_shift
    from nwkit.shift_native_heuristic import NativeSearchOptions
    from nwkit.shift_native_limits import NATIVE_SEARCH_DEFAULTS

    args = parser.parse_args(
        [
            "shift",
            "--selection",
            selection,
            "--trait",
            "traits.tsv",
            "--state-column",
            "x",
            "--model-out",
            "model.json",
        ]
    )
    options = NativeSearchOptions()
    for name, expected in NATIVE_SEARCH_DEFAULTS.items():
        assert getattr(args, name) == expected
        actual = (
            options.memory_limit // 1024**2
            if name == "search_memory_mb"
            else getattr(options, name)
        )
        assert actual == expected
    if selection != "native":
        sentinel = object()
        monkeypatch.setattr("nwkit.shift.shift_main", lambda _: sentinel)
        assert _command_shift(args) is sentinel


def test_auto_shift_cap_and_convergent_search_match_explicit_cap(
    native_inputs, tmp_path
):
    args = native_inputs.copy()
    index = args.index("--regime-map")
    del args[index : index + 2]
    args.extend(
        [
            "--max-shifts",
            "auto",
            "--criterion",
            "AICc",
            "--convergence",
            "--search-strategy",
            "lasso",
            "--candidate-pool",
            "3",
            "--refit-budget",
            "8",
            "--bootstrap",
            "1",
        ]
    )
    main(args)
    auto = json.loads((tmp_path / "model.json").read_text())
    assert auto["configuration"]["max_shifts"] == "auto"
    assert auto["configuration"]["convergence"]
    assert auto["search"]["shift_limit"] == {
        "requested": "auto",
        "resolved": 3,
        "constraints": {"tree": 6, "refit_budget": 7, "candidate_pool": 3},
        "budget_limited": True,
    }
    assert any(
        len(row["groups"]) < len(row["shift_branch_ids"]) + 1
        for row in auto["candidates"]
    )
    args[args.index("--max-shifts") + 1] = "3"
    main(args)
    explicit = json.loads((tmp_path / "model.json").read_text())
    assert auto["candidates"] == explicit["candidates"]
    assert auto["shift_branch_ids"] == explicit["shift_branch_ids"]
    assert auto["information_criterion"] == explicit["information_criterion"]
    assert auto["selection_support"] == explicit["selection_support"]


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


@pytest.mark.parametrize("criterion", ["AIC", "AICc", "BIC", "pBIC"])
def test_native_fixed_layout_exports_information_criterion(
    native_inputs, tmp_path, criterion
):
    main(native_inputs + ["--criterion", criterion])
    model = json.loads((tmp_path / "model.json").read_text())
    assert model["information_criterion"]["criterion"] == criterion
    assert model["information_criterion"]["status"] == "ok"


@pytest.mark.parametrize("criterion", ["AIC", "AICc", "BIC", "pBIC"])
def test_native_information_selection_replays_support_without_calibration(
    native_inputs, tmp_path, criterion, monkeypatch
):
    def unexpected_calibration(*args, **kwargs):
        raise AssertionError("Information criteria must not run calibration")

    monkeypatch.setattr(
        "nwkit.shift_native_selection.calibrate_native_search", unexpected_calibration
    )
    arguments = list(native_inputs)
    index = arguments.index("--regime-map")
    del arguments[index : index + 2]
    main(
        arguments
        + [
            "--criterion",
            criterion,
            "--max-shifts",
            "1",
            "--bootstrap",
            "2",
            "--calibration-replicates",
            "0",
        ]
    )
    model = json.loads((tmp_path / "model.json").read_text())
    assert model["information_criterion"]["criterion"] == criterion
    assert model["selection_calibration"] is None
    assert model["selection_support"] is not None
    arguments[arguments.index("--model-out") + 1] = str(tmp_path / "resumed.json")
    with pytest.raises(ValueError, match="configuration"):
        main(
            arguments
            + [
                "--criterion",
                "BIC" if criterion == "AIC" else "AIC",
                "--resume-model",
                str(tmp_path / "model.json"),
                "--max-shifts",
                "1",
                "--bootstrap",
                "2",
                "--calibration-replicates",
                "0",
            ]
        )


@pytest.mark.parametrize(
    "strategy, expected",
    [
        ("auto", "group_lasso_beam"),
        ("native-path", "covariance_updated_optimum_path"),
        ("lasso", "group_lasso_beam"),
    ],
)
@pytest.mark.parametrize("criterion", ["AIC", "AICc"])
def test_native_aic_strategy_dispatch_and_support_replay(
    native_inputs, tmp_path, strategy, expected, criterion
):
    arguments = list(native_inputs)
    index = arguments.index("--regime-map")
    del arguments[index : index + 2]
    main(
        arguments
        + [
            "--criterion",
            criterion,
            "--max-shifts",
            "1",
            "--search-strategy",
            strategy,
            "--exhaustive-max-configurations",
            "1",
            "--bootstrap",
            "2",
        ]
    )
    model = json.loads((tmp_path / "model.json").read_text())
    assert model["search"]["strategy"] == expected
    assert model["selection_support"] is not None
    assert model["selection_calibration"] is None


def test_native_aic_global_null_gate_cli_and_resume(native_inputs, tmp_path):
    index = native_inputs.index("--regime-map")
    del native_inputs[index : index + 2]
    native_inputs.extend(
        [
            "--criterion",
            "AIC",
            "--search-strategy",
            "native-path",
            "--global-null-gate",
            "--max-shifts",
            "1",
            "--calibration-replicates",
            "19",
        ]
    )
    main(native_inputs)
    model = json.loads((tmp_path / "model.json").read_text())
    calibration = model["selection_calibration"]
    assert calibration["full_search_repeated"]
    assert calibration["replicates"] == 19
    assert model["configuration"]["global_null_gate"] is True
    if not calibration["rejected"]:
        assert model["shift_branch_ids"] == []
    native_inputs[native_inputs.index("--model-out") + 1] = str(
        tmp_path / "resumed.json"
    )
    native_inputs.extend(["--resume-model", str(tmp_path / "model.json")])
    main(native_inputs)
    native_inputs.remove("--global-null-gate")
    with pytest.raises(ValueError, match="configuration_sha256 differs"):
        main(native_inputs)


@pytest.mark.parametrize(
    "selection,criterion", [("native", "BIC"), ("ic", "AIC"), ("calibrated", "AIC")]
)
def test_global_null_gate_rejects_other_selection_modes(
    native_inputs, selection, criterion
):
    index = native_inputs.index("--regime-map")
    del native_inputs[index : index + 2]
    native_inputs[native_inputs.index("--selection") + 1] = selection
    native_inputs.extend(["--criterion", criterion, "--global-null-gate"])
    with pytest.raises(ValueError, match="requires native AIC search"):
        main(native_inputs)


def test_global_null_gate_rejects_fixed_layout(native_inputs):
    native_inputs.extend(["--criterion", "AIC", "--global-null-gate"])
    with pytest.raises(ValueError, match="without --regime-map"):
        main(native_inputs)


def test_native_support_repeats_global_gate_with_distinct_seeds(
    native_inputs, tmp_path, monkeypatch
):
    import nwkit.shift_native_selection as selection

    index = native_inputs.index("--regime-map")
    del native_inputs[index : index + 2]
    native_inputs.extend(
        [
            "--criterion",
            "AIC",
            "--global-null-gate",
            "--max-shifts",
            "1",
            "--calibration-replicates",
            "19",
            "--bootstrap",
            "2",
        ]
    )
    original = selection.gate_native_aic
    seeds = []

    def gate(*args, **kwargs):
        seeds.append(kwargs["seed"])
        return original(*args, **kwargs)

    monkeypatch.setattr(selection, "gate_native_aic", gate)
    main(native_inputs)
    assert len(seeds) == len(set(seeds)) == 3
    model = json.loads((tmp_path / "model.json").read_text())
    assert model["selection_support"]["replicates"] == 2
