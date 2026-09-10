"""Shared ASR consumers and exact end-jump path refinement."""

import json

import numpy as np
import pandas as pd
import pytest

from nwkit.asr_paths import _grid_counts, simulate_fitted_paths
from nwkit.branch_gaussian_asr import compute_branch_marginals
from nwkit.branch_gaussian_paths import refine_branch_process
from nwkit.cli import main
from nwkit.gaussian_inference import gaussian_tree_likelihood
from tests.test_branch_gaussian_cli import command, process_from
from tests.test_branch_gaussian_cli import inputs as branch_inputs


@pytest.fixture(name="inputs")
def _inputs(tmp_path):
    return branch_inputs.__wrapped__(tmp_path)


def fitted(inputs):
    process, assignment = process_from(inputs)
    observations = {"A": 1.2, "B": -0.4, "C": 2.1}
    errors = dict.fromkeys(observations, 0.0)
    posterior, fit = compute_branch_marginals(
        process.tree, observations, errors, assignment=assignment, root=process.root
    )
    return process.tree, observations, errors, posterior, fit


@pytest.mark.parametrize("steps", [1, 7, 41])
def test_refining_diffusion_keeps_endpoint_joint_law_and_one_jump(inputs, steps):
    tree, observed, errors, _, fit = fitted(inputs)
    counts = _grid_counts(tree, steps, 1, 1)
    refined, chains = refine_branch_process(tree, fit, counts)
    original_nodes = list(tree.traverse())
    refined_nodes = [
        refined.tree if node.is_root else chains[node][0][-1] for node in original_nodes
    ]
    np.testing.assert_allclose(
        refined.covariance(refined_nodes),
        fit.process.covariance(original_nodes),
        atol=2e-14,
    )
    expected = gaussian_tree_likelihood(fit.process, observed, standard_errors=errors)
    actual = gaussian_tree_likelihood(refined, observed, standard_errors=errors)
    assert actual.log_likelihood == pytest.approx(expected.log_likelihood, abs=2e-13)
    assert len(refined.transitions) == sum(counts.values()) + len(fit.jump_nodes)
    for node in fit.jump_nodes:
        chain, times = chains[node]
        assert times[-1] == times[-2] == node.dist
        assert sum(point.dist == 0 for point in chain[1:]) == 1


@pytest.mark.parametrize("mode", ["conditional", "unconditional"])
def test_branch_paths_are_seeded_and_share_endpoints(inputs, mode):
    tree, observed, errors, posterior, fit = fitted(inputs)
    options = dict(
        fit=fit, model="BRANCH-GAUSSIAN", count=4, steps=9, mode=mode, seed=5
    )
    paths = simulate_fitted_paths(tree, observed, errors, posterior, **options)
    repeated = simulate_fitted_paths(tree, observed, errors, posterior, **options)
    for node, path in paths.branches.items():
        np.testing.assert_array_equal(path.values, repeated.branches[node].values)
        parent = (
            paths.root_values
            if node.up.is_root
            else paths.branches[node.up].values[:, -1, :]
        )
        np.testing.assert_array_equal(path.values[:, 0, :], parent)
        if mode == "conditional" and node.is_leaf:
            np.testing.assert_allclose(
                path.values[:, -1, 0], observed[node.name], atol=1e-14
            )


def test_zero_length_pure_jump_survives_grid_and_conditioning(inputs):
    inputs["tree.nwk"].write_text("((A:0.3,B:0)I:0.4,C:1.1)R;")
    tree, observed, errors, posterior, fit = fitted(inputs)
    paths = simulate_fitted_paths(
        tree,
        observed,
        errors,
        posterior,
        fit=fit,
        model="BRANCH-GAUSSIAN",
        count=3,
        steps=7,
        mode="conditional",
        seed=1,
    )
    node = next(node for node in tree.leaves() if node.name == "B")
    assert paths.branches[node].elapsed.tolist() == [0, 0, 0]
    np.testing.assert_allclose(paths.branches[node].values[:, -1, 0], -0.4, atol=1e-14)
    assert np.any(paths.branches[node].values[:, 0, 0] != -0.4)


def test_grid_budget_includes_added_jumps(inputs, monkeypatch):
    tree, observed, errors, posterior, fit = fitted(inputs)
    monkeypatch.setattr("nwkit.branch_gaussian_paths._MAX_GRID_NODES", 5)
    with pytest.raises(ValueError, match="including end jumps"):
        simulate_fitted_paths(
            tree,
            observed,
            errors,
            posterior,
            fit=fit,
            model="BRANCH-GAUSSIAN",
            steps=1,
            seed=1,
        )


def test_summary_tree_model_samples_and_diagnostics_share_process(inputs, tmp_path):
    paths = {
        key: tmp_path / filename
        for key, filename in {
            "out": "summary.tsv",
            "model": "model.tsv",
            "tree": "asr.nwk",
            "samples": "samples.tsv",
            "ppc": "ppc.tsv",
            "cv": "cv.tsv",
            "process": "process.json",
        }.items()
    }
    args = command(inputs)
    for option, key in [
        ("-o", "out"),
        ("--model-out", "model"),
        ("--tree-out", "tree"),
        ("--posterior-samples-out", "samples"),
        ("--posterior-predictive-out", "ppc"),
        ("--cross-validation-out", "cv"),
        ("--process-out", "process"),
    ]:
        args.extend([option, str(paths[key])])
    main(
        [
            *args,
            "--posterior-samples",
            "4",
            "--posterior-predictive-simulations",
            "3",
            "--seed",
            "17",
        ]
    )
    assert len(pd.read_csv(paths["samples"], sep="\t")) == 20
    assert "asr_model=BRANCH-GAUSSIAN" in paths["tree"].read_text()
    model = pd.read_csv(paths["model"], sep="\t").iloc[0]
    assert model.num_parameters_estimated == 0
    assert model.log_likelihood == pytest.approx(
        json.loads(paths["process"].read_text())["log_likelihood"]
    )
    assert len(pd.read_csv(paths["cv"], sep="\t")) == 3
    assert not pd.read_csv(paths["ppc"], sep="\t").empty


@pytest.mark.parametrize("prior", [False, True])
@pytest.mark.parametrize("extension", ["png", "svg", "pdf"])
def test_cli_figures_and_prior_labels(inputs, tmp_path, extension, prior):
    figure = tmp_path / f"figure.{extension}"
    extra = (
        ["--output", "prior-samples", "--prior-samples", "3"]
        if prior
        else ["--figure-simulation-mode", "conditional"]
    )
    main(
        [
            *command(inputs, observe=not prior),
            *extra,
            "--figure-out",
            str(figure),
            "--figure-simulations",
            "2",
            "--figure-simulation-steps",
            "5",
            "--species-overlap-node-plot",
            "no",
            "--seed",
            "8",
            "-o",
            str(tmp_path / "out.tsv"),
        ]
    )
    assert figure.stat().st_size > 2000
    if extension == "svg":
        svg = figure.read_text()
        assert "Prescribed Gaussian end jump" in svg
        assert (
            "Prior trait distribution" if prior else "Ancestral trait reconstruction"
        ) in svg
        assert ("Observed tip" not in svg) == prior
        assert "Imputed tip" not in svg


def test_flat_prior_figure_uses_supplied_root_and_no_inferred_label(inputs, tmp_path):
    figure = tmp_path / "prior.svg"
    main(
        [
            *command(inputs, observe=False, root=["--root-prior", "flat"]),
            "--output",
            "prior-samples",
            "--prior-samples",
            "2",
            "--prior-root-value",
            "2.5",
            "--figure-out",
            str(figure),
            "--figure-simulations",
            "1",
            "--figure-simulation-steps",
            "3",
            "-o",
            str(tmp_path / "out.tsv"),
        ]
    )
    svg = figure.read_text()
    assert "Root fixed at the supplied starting value" in svg
    assert "inferred" not in svg


def test_figure_failure_preserves_tables_and_suppresses_stdout(
    inputs, tmp_path, monkeypatch, capsys
):
    figure, table, model = [
        tmp_path / name for name in ("figure.png", "model.tsv", "model.json")
    ]
    for path in (figure, table, model):
        path.write_text("original")

    def fail(*args, **kwargs):
        raise OSError("injected figure failure")

    monkeypatch.setattr("matplotlib.figure.Figure.savefig", fail)
    with pytest.raises(OSError, match="injected figure"):
        main(
            [
                *command(inputs),
                "--figure-out",
                str(figure),
                "--model-out",
                str(table),
                "--process-out",
                str(model),
            ]
        )
    assert capsys.readouterr().out == ""
    assert all(path.read_text() == "original" for path in (figure, table, model))


@pytest.mark.parametrize(
    "extra",
    [
        ["--output", "likelihood", "--figure-out", "unused.png"],
        ["--output", "prior-samples"],
        ["--sigma2", "1"],
        ["--bootstrap-out", "unused.tsv"],
        ["--model", "BM"],
        ["--model", "OU"],
    ],
)
def test_undefined_combinations_fail_before_outputs(inputs, tmp_path, extra):
    out = tmp_path / "out.tsv"
    out.write_text("original")
    with pytest.raises(ValueError):
        main([*command(inputs), "-o", str(out), *extra])
    assert out.read_text() == "original"


def test_root_is_required_explicitly(inputs):
    with pytest.raises(ValueError, match="explicit --root-prior"):
        main(command(inputs, root=[]))


def test_comparison_marks_fixed_branch_model_inapplicable(inputs, tmp_path):
    from nwkit.asr_models import model_names

    output = tmp_path / "comparison.tsv"
    args = [
        "asrcompare",
        "-i",
        str(inputs["tree.nwk"]),
        "--input-rooted",
        "yes",
        "--trait",
        str(inputs["traits.tsv"]),
        "--state-column",
        "x",
        "-o",
        str(output),
    ]
    excluded = ",".join(
        model
        for model in model_names("continuous")
        if model not in {"BRANCH-GAUSSIAN", "BM"}
    )
    main([*args, "--models", "all", "--exclude-models", excluded])
    row = pd.read_csv(output, sep="\t").set_index("model").loc["BRANCH-GAUSSIAN"]
    assert row.status == "not_applicable"
    previous = output.read_bytes()
    with pytest.raises(ValueError, match="not applicable"):
        main([*args, "--models", "BRANCH-GAUSSIAN"])
    assert output.read_bytes() == previous
