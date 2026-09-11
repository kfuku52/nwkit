import json
from dataclasses import replace

import numpy as np
import pandas as pd
import pytest

from nwkit.cli import main
from nwkit.shift_joint_model import evaluate_joint
from nwkit.shift_joint_screen import JointQuickProfile, joint_group_lasso_screen
from nwkit.shift_native_bootstrap import simulate_native_data
from nwkit.shift_native_fit import NativeFitOptions, fit_native_layout
from nwkit.shift_native_heuristic import NativeSearchOptions, heuristic_native_search
from nwkit.shift_native_model import ShiftLayout
from nwkit.shift_simulation import simulate_shift, simulation_from_fit
from nwkit.shift_simulation_cli import explicit_simulation
from tests.test_shift_joint_covariance import dense_ou, fixture_data


@pytest.mark.parametrize("root_model", ["OUfixedRoot", "OUrandomRoot"])
def test_correlated_simulation_moments_and_observation_mask(root_model):
    data = fixture_data()
    parameters = {
        "trait_names": ["x", "y"],
        "root_model": root_model,
        "alpha": [0.5, 2.0],
        "process_tip_covariance": [[1.0, 0.5], [0.5, 1.0]],
        "shift_branch_ids": [2],
        "regime_optima": [[0.2, -0.3], [1.0, -1.0]],
        "sampling_standard_errors": [0.1, 0.2],
        "measurement_covariance": [[0.1, 0.03], [0.03, 0.2]],
    }
    spec = explicit_simulation(data.tree, parameters)
    values, latent = simulate_shift(spec, 30000, seed=7)
    # Independent physical-OU mean propagation from the explicit optima.
    expected = np.zeros((len(data.tree.branch_ids), 2))
    expected[0] = parameters["regime_optima"][0]
    labels = spec.layout.node_groups(data.tree)
    for i in range(1, len(expected)):
        attenuation = np.exp(-spec.alpha_height * data.tree.times[i])
        expected[i] = attenuation * expected[data.tree.compiled.parents[i]] + (
            1 - attenuation
        ) * np.asarray(parameters["regime_optima"][labels[i]])
    np.testing.assert_allclose(latent.mean(axis=0), expected, atol=0.025)
    dense = dense_ou(
        data, spec.alpha_height, spec.covariance_coordinate, spec.root_model
    )
    for i in range(8):
        dense[i * 2 : (i + 1) * 2, i * 2 : (i + 1) * 2] += (
            spec.measurement_covariance + np.diag(spec.sampling_variances[i])
        )
    np.testing.assert_allclose(
        np.cov(values.reshape(30000, -1), rowvar=False), dense, atol=0.04
    )
    mask = np.zeros((8, 2), bool)
    mask[2, 1] = True
    missing = replace(spec, missing=mask)
    a, _ = simulate_shift(missing, 3, seed=8)
    b, _ = simulate_shift(missing, 3, seed=8)
    np.testing.assert_array_equal(a, b)
    assert np.isnan(a[:, 2, 1]).all()
    assert np.isfinite(a[:, 0]).all()


def test_fitted_simulation_and_joint_quick_profile():
    data = fixture_data(0.02, True)
    options = NativeFitOptions(
        trait_covariance="full", alpha_model="shared", optimizer_starts=2
    )
    fitted = fit_native_layout(
        data, ShiftLayout.build(data.tree), options=options, alpha_height=0.8
    )
    branches = list(range(1, len(data.tree.branch_ids)))
    profile = JointQuickProfile(data, fitted, branches, 0.8)
    parent = data.tree.branch_ids[1]
    child = data.tree.branch_ids[data.tree.compiled.children[1][0]]
    for layout in [
        ShiftLayout.build(data.tree, [2]),
        ShiftLayout.build(data.tree, [parent, child], [(0, child), (parent,)]),
    ]:
        joint = fitted["joint_fit"]
        direct = evaluate_joint(
            data,
            layout,
            joint.alpha_height,
            joint.covariance_coordinate,
            joint.measurement_variance,
        )
        adjustment = sum(
            np.isfinite(data.values[:, j]).sum() * np.log(data.scales[j])
            for j in range(2)
        )
        assert profile.score(layout) == pytest.approx(
            direct.log_likelihood - adjustment, abs=1e-9
        )
    simulated = simulate_native_data(data, fitted, np.random.default_rng(21))
    np.testing.assert_array_equal(np.isnan(simulated.values), np.isnan(data.values))
    np.testing.assert_allclose(
        simulated.variances * simulated.scales**2, data.variances * data.scales**2
    )
    spec = simulation_from_fit(data, fitted)
    assert spec.coefficients.shape == (1, 2)
    with pytest.raises(ValueError, match="memory"):
        joint_group_lasso_screen(data, fitted, memory_limit=1)


def test_joint_heuristic_search_uses_correlated_screen():
    data = fixture_data()
    options = NativeFitOptions(trait_covariance="full", alpha_model="shared")
    search = heuristic_native_search(
        data,
        options=NativeSearchOptions(
            max_shifts=1, candidate_pool=4, refit_budget=7, lasso_iterations=20
        ),
        fit_arguments={"options": options, "alpha_height": 0.7},
        criterion="AIC",
    )
    assert (
        search.metadata["screening"]["method"]
        == "joint_covariance_whitened_group_lasso"
    )
    assert search.best_information["trait_covariance"] == "full"


def test_simulation_cli_and_fit_model_roundtrip(tmp_path):
    data = fixture_data()
    tree = tmp_path / "tree.nwk"
    tree.write_text(data.tree.compiled.tree.write(format_root_node=True))
    parameters = tmp_path / "parameters.json"
    parameters.write_text(
        json.dumps(
            {
                "trait_names": ["x", "y"],
                "alpha": 0.7,
                "process_tip_covariance": [[1.0, 0.8], [0.8, 1.0]],
                "regime_optima": [[0.0, 0.0]],
            }
        )
    )
    output = tmp_path / "traits.tsv"
    truth = tmp_path / "truth.json"
    common = ["--infile", str(tree), "--input-rooted", "yes"]
    main(
        [
            "shift-simulate",
            *common,
            "--parameters",
            str(parameters),
            "--outfile",
            str(output),
            "--truth-out",
            str(truth),
            "--latent-out",
            str(tmp_path / "latent.tsv"),
            "--seed",
            "2",
        ]
    )
    before = output.read_bytes()
    main(
        [
            "shift-simulate",
            *common,
            "--parameters",
            str(parameters),
            "--outfile",
            str(output),
            "--truth-out",
            str(truth),
            "--seed",
            "2",
        ]
    )
    assert output.read_bytes() == before
    np.testing.assert_allclose(
        json.loads(truth.read_text())["process_tip_covariance"],
        [[1.0, 0.8], [0.8, 1.0]],
        atol=1e-14,
    )
    model = tmp_path / "fit.json"
    main(
        [
            "shift",
            *common,
            "--trait",
            str(output),
            "--state-column",
            "x,y",
            "--selection",
            "native",
            "--trait-covariance",
            "full",
            "--alpha-model",
            "shared",
            "--alpha",
            ".7",
            "--max-shifts",
            "0",
            "--criterion",
            "AIC",
            "--model-out",
            str(model),
            "--outfile",
            str(tmp_path / "map.tsv"),
        ]
    )
    fitted = json.loads(model.read_text())
    assert fitted["trait_covariance"] == "full"
    assert fitted["configuration"]["alpha_model"] == "shared"
    main(
        [
            "shift-simulate",
            *common,
            "--model-in",
            str(model),
            "--outfile",
            str(tmp_path / "again.tsv"),
            "--truth-out",
            str(tmp_path / "again.json"),
            "--replicates",
            "2",
        ]
    )
    again = pd.read_csv(tmp_path / "again.tsv", sep="\t")
    assert len(again) == 16 and again["replicate"].nunique() == 2
    with pytest.raises(ValueError, match="input"):
        main(
            [
                "shift-simulate",
                *common,
                "--parameters",
                str(parameters),
                "--outfile",
                str(parameters),
                "--truth-out",
                str(truth),
            ]
        )


def test_full_bootstrap_cli_uses_complete_search(tmp_path):
    data = fixture_data()
    tree = tmp_path / "tree.nwk"
    tree.write_text(data.tree.compiled.tree.write(format_root_node=True))
    table = tmp_path / "traits.tsv"
    pd.DataFrame(
        {
            "leaf_name": data.tree.leaf_names,
            "x": data.values[:, 0],
            "y": data.values[:, 1],
        }
    ).to_csv(table, sep="\t", index=False)
    model = tmp_path / "fit.json"
    main(
        [
            "shift",
            "--infile",
            str(tree),
            "--input-rooted",
            "yes",
            "--trait",
            str(table),
            "--state-column",
            "x,y",
            "--selection",
            "native",
            "--trait-covariance",
            "full",
            "--alpha-model",
            "shared",
            "--alpha",
            ".7",
            "--max-shifts",
            "1",
            "--calibration-replicates",
            "19",
            "--model-out",
            str(model),
            "--outfile",
            str(tmp_path / "map.tsv"),
        ]
    )
    result = json.loads(model.read_text())
    assert result["selection_calibration"]["full_search_repeated"]
    assert result["search"]["strategy"] == "exhaustive"
