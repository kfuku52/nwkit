"""Independent likelihood oracles and ASR contracts for fixed-layout estimation."""

import json
import math
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy.optimize import differential_evolution

from nwkit.branch_gaussian import BranchGaussianModel, BrownianBranch, OUBranch
from nwkit.branch_gaussian_asr import compute_branch_marginals
from nwkit.branch_gaussian_fit_spec import BranchFitParameter, load_branch_fit_spec
from nwkit.branch_gaussian_input import (
    BranchGaussianAssignment,
    load_branch_gaussian_assignment,
)
from nwkit.cli import main
from nwkit.gaussian_tree import GaussianRootPrior
from nwkit.util import assign_branch_ids, read_tree

FIT_HEADER = "regime\tparameter\tgroup\tlower\tupper\n"
EXAMPLE = Path(__file__).resolve().parents[1] / "examples/branch_gaussian/plot"


def brownian_problem(tmp_path):
    path = tmp_path / "tree.nwk"
    path.write_text("((A:0.8,B:1.2)I:0.5,(C:0.7,D:1.3)J:0.6,(E:0.9,F:1.1)K:0.4)R;")
    tree = read_tree(str(path), "1", True, quiet=True, rooted="yes")
    ids = assign_branch_ids(tree)
    models = {
        i: BranchGaussianModel(BrownianBranch(1.0))
        for node, i in ids.items()
        if not node.is_root
    }
    regimes = {i: "baseline" for i in models}
    values = dict(zip("ABCDEF", (1.0, -0.5, 2.3, 0.4, -1.2, 0.7), strict=True))
    return tree, BranchGaussianAssignment(models, regimes), values


def dense_moments(tree, assignment, root):
    """Assemble independent noise loadings from biological BM/OU parameters."""
    nodes = list(tree.traverse("preorder"))
    ids = assign_branch_ids(tree)
    means = {tree: root.mean}
    loadings = {tree: np.zeros(len(nodes))}
    loadings[tree][0] = math.sqrt(root.variance or 0.0)
    for index, node in enumerate(nodes[1:], start=1):
        model = assignment.models_by_branch_id[ids[node]]
        diffusion = model.diffusion
        if isinstance(diffusion, OUBranch) and diffusion.alpha > 0:
            a = math.exp(-diffusion.alpha * node.dist)
            b = -math.expm1(-diffusion.alpha * node.dist) * diffusion.optimum
            q = (
                diffusion.variance_rate
                * -math.expm1(-2 * diffusion.alpha * node.dist)
                / (2 * diffusion.alpha)
            )
        else:
            a, b = 1.0, 0.0
            q = 0.0 if diffusion is None else diffusion.variance_rate * node.dist
        if model.jump is not None:
            b += model.jump.mean
            q += model.jump.variance
        means[node] = a * means[node.up] + b
        loadings[node] = a * loadings[node.up]
        loadings[node][index] = math.sqrt(q)
    leaves = list(tree.leaves())
    matrix = np.array([loadings[node] for node in leaves])
    return np.array([means[node] for node in leaves]), matrix @ matrix.T


@pytest.mark.parametrize("mode", ["fixed", "flat"])
def test_bm_rate_matches_closed_form_ml_or_flat_root_integral(tmp_path, mode):
    tree, assignment, observed = brownian_problem(tmp_path)
    root = GaussianRootPrior(mode, 0.0, None if mode == "flat" else 0.0)
    spec = (BranchFitParameter("rate", "sigma2", ("baseline",), 1.0, 0.01, 20.0),)
    posterior, fit = compute_branch_marginals(
        tree, observed, None, assignment=assignment, root=root, fit_spec=spec
    )
    _, covariance = dense_moments(tree, assignment, GaussianRootPrior("fixed"))
    y = np.array([observed[str(n.name)] for n in tree.leaves()])
    inverse = np.linalg.inv(covariance)
    mu = (np.ones(len(y)) @ inverse @ y) / inverse.sum() if mode == "flat" else 0.0
    residual = y - mu
    expected = residual @ inverse @ residual / (len(y) - int(mode == "flat"))
    assert fit.estimation["groups"][0]["estimate"] == pytest.approx(expected, rel=2e-6)
    assert fit.estimation["parameter_rank"] == 1
    assert fit.estimation["optimizer_projected_gradient"] <= 1e-4
    assert fit.estimation["initial_log_likelihood"] <= fit.log_likelihood
    assert posterior[tree].mean == pytest.approx(mu, abs=1e-8)


@pytest.fixture
def mixed_inputs(tmp_path):
    for name in ("tree.nwk", "traits.tsv", "branch_regimes.tsv"):
        (tmp_path / name).write_bytes((EXAMPLE / name).read_bytes())
    definitions = pd.read_csv(EXAMPLE / "regime_models.tsv", sep="\t")
    definitions["sigma2"] = 0.3
    definitions.loc[definitions.model == "OU", "alpha"] = 1.8
    definitions.to_csv(tmp_path / "regime_models.tsv", sep="\t", index=False)
    (tmp_path / "fit.tsv").write_text(
        FIT_HEADER + "Background\tsigma2\trate\t0.01\t3\n"
        "Low optimum\tsigma2\trate\t0.01\t3\n"
        "High optimum\tsigma2\trate\t0.01\t3\n"
        "High optimum + event\tsigma2\trate\t0.01\t3\n"
        "Low optimum\talpha\tpull\t0\t8\n"
        "High optimum\talpha\tpull\t0\t8\n"
        "High optimum + event\talpha\tpull\t0\t8\n"
        "High optimum\ttheta\thigh\t0\t5\n"
        "High optimum + event\ttheta\thigh\t0\t5\n"
    )
    return tmp_path


def mixed_problem(directory):
    tree = read_tree(str(directory / "tree.nwk"), "1", True, quiet=True, rooted="yes")
    assignment = load_branch_gaussian_assignment(
        tree,
        branch_regimes=directory / "branch_regimes.tsv",
        regime_models=directory / "regime_models.tsv",
    )
    spec = load_branch_fit_spec(directory / "fit.tsv", assignment, tree)
    data = pd.read_csv(directory / "traits.tsv", sep="\t").set_index("leaf_name")
    observed = {
        name: None if pd.isna(row.x) else row.x for name, row in data.iterrows()
    }
    errors = {name: row.se for name, row in data.iterrows()}
    return tree, assignment, spec, observed, errors


def mixed_command(directory):
    return [
        "asr",
        "--model",
        "BRANCH-GAUSSIAN",
        "-i",
        str(directory / "tree.nwk"),
        "--input-rooted",
        "yes",
        "--trait",
        str(directory / "traits.tsv"),
        "--state-column",
        "x",
        "--standard-error-column",
        "se",
        "--branch-regimes",
        str(directory / "branch_regimes.tsv"),
        "--regime-models",
        str(directory / "regime_models.tsv"),
        "--branch-fit",
        str(directory / "fit.tsv"),
        "--root-prior",
        "gaussian",
        "--root-mean",
        "0",
        "--root-variance",
        "0.2",
    ]


def test_mixed_shared_fit_matches_independent_dense_optimization(mixed_inputs):
    tree, assignment, spec, observed, errors = mixed_problem(mixed_inputs)
    root = GaussianRootPrior("gaussian", 0.0, 0.2)
    _, fit = compute_branch_marginals(
        tree, observed, errors, assignment=assignment, root=root, fit_spec=spec
    )
    leaf_names = list(tree.leaf_names())
    keep = [i for i, name in enumerate(leaf_names) if observed[name] is not None]
    y = np.array([observed[leaf_names[i]] for i in keep])

    def objective(values):
        rate, alpha, theta = math.exp(values[0]), values[1], values[2]
        models = {}
        for identifier, model in assignment.models_by_branch_id.items():
            diffusion = replace(model.diffusion, variance_rate=rate)
            if isinstance(diffusion, OUBranch):
                diffusion = replace(diffusion, alpha=alpha)
                if assignment.regime_by_branch_id[identifier].startswith("High"):
                    diffusion = replace(diffusion, optimum=theta)
            models[identifier] = replace(model, diffusion=diffusion)
        means, covariance = dense_moments(tree, BranchGaussianAssignment(models), root)
        covariance = covariance[np.ix_(keep, keep)] + np.diag(
            [errors[leaf_names[i]] ** 2 for i in keep]
        )
        residual = y - means[keep]
        return 0.5 * (
            len(y) * math.log(2 * math.pi)
            + np.linalg.slogdet(covariance)[1]
            + residual @ np.linalg.solve(covariance, residual)
        )

    oracle = differential_evolution(
        objective,
        [(math.log(0.01), math.log(3)), (0, 8), (0, 5)],
        seed=17,
        tol=1e-9,
        polish=True,
    )
    assert oracle.success
    assert fit.log_likelihood == pytest.approx(-oracle.fun, abs=2e-6)
    estimates = {row["group"]: row["estimate"] for row in fit.estimation["groups"]}
    np.testing.assert_allclose(
        [math.log(estimates["rate"]), estimates["pull"], estimates["high"]],
        oracle.x,
        rtol=2e-4,
        atol=2e-4,
    )
    assert fit.branch_assignment.regime_by_branch_id == assignment.regime_by_branch_id
    for identifier, model in fit.branch_assignment.models_by_branch_id.items():
        assert model.jump == assignment.models_by_branch_id[identifier].jump
    assert fit.estimation["num_parameters_estimated"] == 3
    assert fit.estimation["root_parameters_estimated"] is False


def test_cli_fit_exports_replayable_models_and_plot(mixed_inputs):
    output = mixed_inputs / "asr.tsv"
    normalized = mixed_inputs / "fitted.tsv"
    main(
        [
            *mixed_command(mixed_inputs),
            "-o",
            str(output),
            "--model-out",
            str(mixed_inputs / "model.tsv"),
            "--branch-models-out",
            str(normalized),
            "--process-out",
            str(mixed_inputs / "process.json"),
            "--figure-out",
            str(mixed_inputs / "fit.png"),
            "--figure-simulations",
            "2",
            "--figure-simulation-mode",
            "conditional",
            "--seed",
            "12",
        ]
    )
    metadata = json.loads((mixed_inputs / "process.json").read_text())
    assert metadata["schema_version"] == 2
    assert metadata["estimation"]["num_parameters_estimated"] == 3
    assert (mixed_inputs / "fit.png").read_bytes().startswith(b"\x89PNG")
    table = pd.read_csv(mixed_inputs / "model.tsv", sep="\t").iloc[0]
    assert table.estimation_method == "ML"
    assert json.loads(table.parameter_estimation)["assignment_search"] is False
    command = mixed_command(mixed_inputs)
    for flag in ("--branch-fit", "--branch-regimes", "--regime-models"):
        index = command.index(flag)
        del command[index : index + 2]
    main(
        [
            *command,
            "--branch-models",
            str(normalized),
            "-o",
            str(mixed_inputs / "replay.tsv"),
        ]
    )
    pd.testing.assert_frame_equal(
        pd.read_csv(output, sep="\t"),
        pd.read_csv(mixed_inputs / "replay.tsv", sep="\t"),
        rtol=1e-12,
        atol=1e-12,
    )


@pytest.mark.parametrize(
    "row,match",
    [
        ("unknown\tsigma2\tr\t0.01\t3", "Invalid fit regime"),
        ("Background\talpha\tr\t0\t3", "not used"),
        ("Background\tjump_mean\tr\t0\t3", "only sigma2"),
        ("Background\tsigma2\t\t0.01\t3", "empty group"),
        ("Background\tsigma2\tr\t0\t3", "positive lower"),
        ("High optimum\talpha\tr\t-1\t3", "nonnegative lower"),
        ("Background\tsigma2\tr\t1\t3", "outside"),
        ("Background\tsigma2\tr\t0.01\tinf", "finite"),
        ("Background\tsigma2\tr\tnan\t3", "finite"),
        ("Background\tsigma2\tr\ta\t3", "Non-numeric"),
        ("Background\tsigma2\tr\t3\t0.01", "increasing"),
        ("Background\tsigma2\tr\t0.01\t3\nBackground\tsigma2\tr\t0.01\t3", "Duplicate"),
        (
            "Background\tsigma2\tr\t0.01\t3\nHigh optimum\talpha\tr\t0\t3",
            "same parameter",
        ),
        (
            "Background\tsigma2\tr\t0.01\t3\nHigh optimum\tsigma2\tr\t0.02\t3",
            "same parameter",
        ),
        ("", "between 1 and 20"),
    ],
)
def test_invalid_fit_spec_preserves_all_outputs(mixed_inputs, row, match, capsys):
    (mixed_inputs / "fit.tsv").write_text(FIT_HEADER + row + ("\n" if row else ""))
    out, metadata, figure = [
        mixed_inputs / name for name in ("output.tsv", "model.json", "old.png")
    ]
    for path in (out, metadata, figure):
        path.write_bytes(b"preserve old output")
    with pytest.raises(ValueError, match=match):
        main(
            [
                *mixed_command(mixed_inputs),
                "-o",
                str(out),
                "--process-out",
                str(metadata),
                "--figure-out",
                str(figure),
            ]
        )
    assert capsys.readouterr().out == ""
    for path in (out, metadata, figure):
        assert path.read_bytes() == b"preserve old output"


def test_fit_input_cannot_be_an_output(mixed_inputs):
    before = (mixed_inputs / "fit.tsv").read_bytes()
    with pytest.raises(ValueError, match="input|overwrit|replace"):
        main([*mixed_command(mixed_inputs), "-o", str(mixed_inputs / "fit.tsv")])
    assert (mixed_inputs / "fit.tsv").read_bytes() == before


def test_sharing_requires_matching_initial_values(mixed_inputs):
    text = (
        (mixed_inputs / "regime_models.tsv")
        .read_text()
        .replace("Background\tBM\t0.3", "Background\tBM\t0.4")
    )
    (mixed_inputs / "regime_models.tsv").write_text(text)
    with pytest.raises(ValueError, match="initial value"):
        mixed_problem(mixed_inputs)


def test_alpha_estimation_rejects_flat_root(mixed_inputs):
    tree, assignment, spec, observed, errors = mixed_problem(mixed_inputs)
    with pytest.raises(ValueError, match="alpha requires a proper root"):
        compute_branch_marginals(
            tree,
            observed,
            errors,
            assignment=assignment,
            root=GaussianRootPrior("flat", variance=None),
            fit_spec=spec,
        )


def test_indistinguishable_star_alpha_and_rate_are_rejected(tmp_path):
    tree_path = tmp_path / "star.nwk"
    tree_path.write_text("(A:1,B:1,C:1,D:1,E:1,F:1)R;")
    tree = read_tree(str(tree_path), "1", True, quiet=True, rooted="yes")
    models = {
        i: BranchGaussianModel(OUBranch(1.0, 1.0, 0.0))
        for node, i in assign_branch_ids(tree).items()
        if not node.is_root
    }
    assignment = BranchGaussianAssignment(models, dict.fromkeys(models, "all"))
    spec = (
        BranchFitParameter("rate", "sigma2", ("all",), 1.0, 0.01, 10),
        BranchFitParameter("alpha", "alpha", ("all",), 1.0, 0.1, 3),
    )
    with pytest.raises(ValueError, match="not identifiable"):
        compute_branch_marginals(
            tree,
            dict(zip("ABCDEF", range(6), strict=True)),
            None,
            assignment=assignment,
            root=GaussianRootPrior("fixed"),
            fit_spec=spec,
        )


def test_rate_on_unobserved_branch_is_rejected(tmp_path):
    tree, assignment, observed = brownian_problem(tmp_path)
    identifier = next(
        i for node, i in assign_branch_ids(tree).items() if node.name == "F"
    )
    regimes = {**assignment.regime_by_branch_id, identifier: "missing"}
    assignment = replace(assignment, regime_by_branch_id=regimes)
    observed["F"] = None
    spec = (BranchFitParameter("rate", "sigma2", ("missing",), 1.0, 0.01, 20),)
    with pytest.raises(ValueError, match="not identifiable"):
        compute_branch_marginals(
            tree,
            observed,
            None,
            assignment=assignment,
            root=GaussianRootPrior("fixed"),
            fit_spec=spec,
        )


def test_zero_alpha_makes_free_optimum_unidentifiable(tmp_path):
    tree, assignment, observed = brownian_problem(tmp_path)
    assignment = replace(
        assignment,
        models_by_branch_id={
            i: BranchGaussianModel(OUBranch(0, 1, 0))
            for i in assignment.models_by_branch_id
        },
    )
    spec = (BranchFitParameter("theta", "theta", ("baseline",), 0.0, -3, 3),)
    with pytest.raises(ValueError, match="not identifiable"):
        compute_branch_marginals(
            tree,
            observed,
            None,
            assignment=assignment,
            root=GaussianRootPrior("fixed"),
            fit_spec=spec,
        )


def test_bound_is_reported_as_constrained_fit(tmp_path):
    tree, assignment, observed = brownian_problem(tmp_path)
    observed = dict.fromkeys(observed, 0.0)
    spec = (BranchFitParameter("rate", "sigma2", ("baseline",), 1.0, 0.01, 20),)
    _, fit = compute_branch_marginals(
        tree,
        observed,
        None,
        assignment=assignment,
        root=GaussianRootPrior("fixed"),
        fit_spec=spec,
    )
    assert fit.fit_status == "boundary"
    assert fit.estimation["groups"][0]["boundary"] == "lower"
    assert fit.estimation["groups"][0]["estimate"] == 0.01


def test_claimed_optimizer_success_needs_independent_gradient_check(
    tmp_path, monkeypatch
):
    from scipy.optimize import OptimizeResult

    from nwkit import branch_gaussian_estimation

    def stopped(objective, initial, **kwargs):
        point = np.full_like(initial, 0.9)
        return OptimizeResult(
            x=point, fun=objective(point), success=True, message="premature stopping"
        )

    monkeypatch.setattr(branch_gaussian_estimation, "minimize", stopped)
    tree, assignment, observed = brownian_problem(tmp_path)
    spec = (BranchFitParameter("rate", "sigma2", ("baseline",), 1.0, 0.01, 20),)
    with pytest.raises(ValueError, match="No optimizer start"):
        compute_branch_marginals(
            tree,
            observed,
            None,
            assignment=assignment,
            root=GaussianRootPrior("fixed"),
            fit_spec=spec,
        )


def test_diagnostics_refit_estimated_parameters_with_held_out_data(tmp_path):
    from nwkit.asr_continuous_diagnostics import _refit

    tree, assignment, observed = brownian_problem(tmp_path)
    spec = (BranchFitParameter("rate", "sigma2", ("baseline",), 1.0, 0.01, 20),)
    _, fit = compute_branch_marginals(
        tree,
        observed,
        None,
        assignment=assignment,
        root=GaussianRootPrior("fixed"),
        fit_spec=spec,
    )
    training = {**observed, "A": None}
    refitted = _refit(
        tree, training, None, None, SimpleNamespace(model="BRANCH-GAUSSIAN"), fit, None
    )
    assert refitted.num_observed == 5
    assert refitted.estimation["num_parameters_estimated"] == 1
    assert refitted.estimation["groups"][0]["estimate"] != pytest.approx(
        fit.estimation["groups"][0]["estimate"]
    )


@pytest.mark.parametrize(
    "trait_scale,time_scale", [(1e-4, 1.0), (1e4, 1e-4), (1.0, 1e4)]
)
def test_estimation_respects_trait_and_time_units(tmp_path, trait_scale, time_scale):
    tree, assignment, observed = brownian_problem(tmp_path)
    spec = (BranchFitParameter("rate", "sigma2", ("baseline",), 1.0, 0.01, 20),)
    _, original = compute_branch_marginals(
        tree,
        observed,
        None,
        assignment=assignment,
        root=GaussianRootPrior("fixed"),
        fit_spec=spec,
    )
    rate_scale = trait_scale**2 / time_scale
    for node in tree.traverse():
        if not node.is_root:
            node.dist *= time_scale
    assignment = replace(
        assignment,
        models_by_branch_id={
            i: BranchGaussianModel(BrownianBranch(rate_scale))
            for i in assignment.models_by_branch_id
        },
    )
    spec = (
        replace(
            spec[0], initial=rate_scale, lower=0.01 * rate_scale, upper=20 * rate_scale
        ),
    )
    _, scaled = compute_branch_marginals(
        tree,
        {name: value * trait_scale for name, value in observed.items()},
        None,
        assignment=assignment,
        root=GaussianRootPrior("fixed"),
        fit_spec=spec,
    )
    assert scaled.estimation["groups"][0]["estimate"] / rate_scale == pytest.approx(
        original.estimation["groups"][0]["estimate"], rel=3e-6
    )
    assert scaled.log_likelihood + len(observed) * math.log(
        trait_scale
    ) == pytest.approx(original.log_likelihood, abs=2e-8)


@pytest.mark.parametrize(
    "change,match",
    [
        ("prior", "prior-samples"),
        ("direct", "requires --branch-regimes"),
        ("model", "require --model BRANCH-GAUSSIAN"),
    ],
)
def test_fit_options_reject_incompatible_modes(mixed_inputs, change, match):
    command = mixed_command(mixed_inputs)
    if change == "prior":
        command += ["--output", "prior-samples"]
    elif change == "direct":
        for flag in ("--branch-regimes", "--regime-models"):
            index = command.index(flag)
            del command[index : index + 2]
        command += ["--branch-models", str(mixed_inputs / "unused.tsv")]
    else:
        command[command.index("BRANCH-GAUSSIAN")] = "OU"
    with pytest.raises(ValueError, match=match):
        main(command)


def test_flat_root_optimum_is_rejected_when_confounded_with_root(mixed_inputs):
    tree, assignment, _, observed, errors = mixed_problem(mixed_inputs)
    assignment = BranchGaussianAssignment(
        {
            i: BranchGaussianModel(OUBranch(1.0, 0.3, 0.0))
            for i in assignment.models_by_branch_id
        },
        dict.fromkeys(assignment.models_by_branch_id, "all"),
    )
    spec = (BranchFitParameter("theta", "theta", ("all",), 0.0, -3, 3),)
    with pytest.raises(ValueError, match="not identifiable"):
        compute_branch_marginals(
            tree,
            observed,
            errors,
            assignment=assignment,
            root=GaussianRootPrior("flat", variance=None),
            fit_spec=spec,
        )


def test_mixed_fixed_alpha_flat_root_can_estimate_rates_and_relative_optimum(
    mixed_inputs,
):
    tree, assignment, spec, observed, errors = mixed_problem(mixed_inputs)
    spec = tuple(p for p in spec if p.parameter != "alpha")
    _, fit = compute_branch_marginals(
        tree,
        observed,
        errors,
        assignment=assignment,
        root=GaussianRootPrior("flat", variance=None),
        fit_spec=spec,
    )
    assert fit.estimation["parameter_rank"] == 2
    assert fit.estimation["method"] == "flat_root_integrated"
    assert fit.estimation["initial_log_likelihood"] <= fit.log_likelihood


def test_fit_table_can_be_read_from_stdin(mixed_inputs, monkeypatch):
    import io

    command = mixed_command(mixed_inputs)
    command[command.index(str(mixed_inputs / "fit.tsv"))] = "-"
    monkeypatch.setattr(
        "sys.stdin", io.StringIO((mixed_inputs / "fit.tsv").read_text())
    )
    main(
        [*command, "--output", "likelihood", "-o", str(mixed_inputs / "likelihood.tsv")]
    )
    assert (
        pd.read_csv(mixed_inputs / "likelihood.tsv", sep="\t")
        .iloc[0]
        .num_parameters_estimated
        == 3
    )


@pytest.mark.parametrize("time_scale", [1e-6, 1.0, 1e6])
def test_alpha_coordinates_include_zero_and_respect_time_units(time_scale):
    parameter = BranchFitParameter(
        "alpha",
        "alpha",
        ("OU",),
        1.8 / time_scale,
        0.0,
        8.0 / time_scale,
        0.9 * time_scale,
    )
    for value in (0.0, 1e-7, 0.4, 1.8, 7.9, 8.0):
        coordinate = parameter.coordinate(value / time_scale)
        assert parameter.value(coordinate) * time_scale == pytest.approx(
            value, rel=1e-12, abs=1e-16
        )
        reference = BranchFitParameter("alpha", "alpha", ("OU",), 1.8, 0.0, 8.0, 0.9)
        assert coordinate == pytest.approx(reference.coordinate(value), abs=2e-15)


def test_mixed_fit_and_sharing_are_invariant_to_units_and_table_row_order(mixed_inputs):
    tree, assignment, spec, observed, errors = mixed_problem(mixed_inputs)
    root = GaussianRootPrior("gaussian", 0.0, 0.2)
    _, original = compute_branch_marginals(
        tree, observed, errors, assignment=assignment, root=root, fit_spec=spec
    )
    trait_scale, time_scale = 1e3, 1e-4
    for node in tree.traverse():
        if not node.is_root:
            node.dist *= time_scale
    models = {}
    for i, model in assignment.models_by_branch_id.items():
        diffusion = replace(
            model.diffusion,
            variance_rate=model.diffusion.variance_rate * trait_scale**2 / time_scale,
        )
        if isinstance(diffusion, OUBranch):
            diffusion = replace(
                diffusion,
                alpha=diffusion.alpha / time_scale,
                optimum=diffusion.optimum * trait_scale,
            )
        jump = (
            None
            if model.jump is None
            else replace(
                model.jump,
                mean=model.jump.mean * trait_scale,
                variance=model.jump.variance * trait_scale**2,
            )
        )
        models[i] = replace(model, diffusion=diffusion, jump=jump)
    scales = {
        "sigma2": trait_scale**2 / time_scale,
        "alpha": 1 / time_scale,
        "theta": trait_scale,
    }
    scaled_spec = tuple(
        replace(
            p,
            initial=p.initial * scales[p.parameter],
            lower=p.lower * scales[p.parameter],
            upper=p.upper * scales[p.parameter],
            time_scale=p.time_scale * time_scale
            if p.parameter == "alpha"
            else p.time_scale,
        )
        for p in spec
    )
    _, scaled = compute_branch_marginals(
        tree,
        {k: None if v is None else v * trait_scale for k, v in observed.items()},
        {k: v * trait_scale for k, v in errors.items()},
        assignment=replace(assignment, models_by_branch_id=models),
        root=replace(root, variance=root.variance * trait_scale**2),
        fit_spec=scaled_spec,
    )
    for first, second in zip(
        original.estimation["groups"], scaled.estimation["groups"], strict=True
    ):
        assert second["estimate"] / scales[first["parameter"]] == pytest.approx(
            first["estimate"], rel=2e-5
        )
    assert scaled.log_likelihood + original.num_observed * math.log(
        trait_scale
    ) == pytest.approx(original.log_likelihood, abs=1e-7)
    # Parsing order must not change the specified groups or their starting points.
    lines = (mixed_inputs / "fit.tsv").read_text().splitlines()
    (mixed_inputs / "fit.tsv").write_text(
        "\n".join([lines[0], *reversed(lines[1:])]) + "\n"
    )
    _, _, reordered, _, _ = mixed_problem(mixed_inputs)
    assert reordered == spec


def test_replicate_constant_is_restored_in_initial_and_fitted_likelihood(mixed_inputs):
    from scipy.stats import norm

    main(
        [
            *mixed_command(mixed_inputs),
            "-o",
            str(mixed_inputs / "base.tsv"),
            "--process-out",
            str(mixed_inputs / "base.json"),
        ]
    )
    original = json.loads((mixed_inputs / "base.json").read_text())
    replicates = mixed_inputs / "replicates.tsv"
    se = math.sqrt(2) * 0.1
    replicates.write_text(
        f"leaf_name\ttrait\tvalue\tstandard_error\nA\tx\t-0.45\t{se:.17g}\nA\tx\t-0.35\t{se:.17g}\n"
    )
    main(
        [
            *mixed_command(mixed_inputs),
            "--replicate-observations",
            str(replicates),
            "-o",
            str(mixed_inputs / "replicated.tsv"),
            "--process-out",
            str(mixed_inputs / "replicated.json"),
        ]
    )
    replicated = json.loads((mixed_inputs / "replicated.json").read_text())
    constant = norm.logpdf(0.1, scale=0.2)
    assert replicated["replicate_log_constant"] == pytest.approx(constant)
    assert replicated["log_likelihood"] == pytest.approx(
        original["log_likelihood"] + constant, abs=1e-8
    )
    assert replicated["estimation"]["initial_log_likelihood"] == pytest.approx(
        original["estimation"]["initial_log_likelihood"] + constant, abs=1e-8
    )
    for first, second in zip(
        original["estimation"]["groups"],
        replicated["estimation"]["groups"],
        strict=True,
    ):
        assert first["estimate"] == pytest.approx(second["estimate"], rel=2e-5)


def test_fitted_asr_exports_tree_posterior_predictive_and_refitted_cv(mixed_inputs):
    lines = (mixed_inputs / "fit.tsv").read_text().splitlines()
    (mixed_inputs / "fit.tsv").write_text(
        "\n".join([lines[0], *[line for line in lines[1:] if "\tsigma2\t" in line]])
        + "\n"
    )
    main(
        [
            *mixed_command(mixed_inputs),
            "-o",
            str(mixed_inputs / "summary.tsv"),
            "--tree-out",
            str(mixed_inputs / "asr.nwk"),
            "--posterior-samples-out",
            str(mixed_inputs / "samples.tsv"),
            "--posterior-samples",
            "2",
            "--posterior-predictive-out",
            str(mixed_inputs / "ppc.tsv"),
            "--posterior-predictive-simulations",
            "3",
            "--cross-validation-out",
            str(mixed_inputs / "cv.tsv"),
            "--seed",
            "15",
        ]
    )
    assert "asr_model=BRANCH-GAUSSIAN" in (mixed_inputs / "asr.nwk").read_text()
    assert len(pd.read_csv(mixed_inputs / "samples.tsv", sep="\t")) == 30
    assert not pd.read_csv(mixed_inputs / "ppc.tsv", sep="\t").empty
    validation = pd.read_csv(mixed_inputs / "cv.tsv", sep="\t")
    assert len(validation) == 7
    assert validation.log_score.notna().all()


def test_extremely_poor_initial_rate_preserves_likelihood_resolution(tmp_path):
    tree, assignment, observed = brownian_problem(tmp_path)
    _, covariance = dense_moments(tree, assignment, GaussianRootPrior("fixed"))
    values = np.array([observed[str(node.name)] for node in tree.leaves()])
    expected = values @ np.linalg.solve(covariance, values) / len(values)
    assignment = replace(
        assignment,
        models_by_branch_id={
            key: replace(model, diffusion=BrownianBranch(1e-200))
            for key, model in assignment.models_by_branch_id.items()
        },
    )
    _, fit = compute_branch_marginals(
        tree,
        observed,
        None,
        assignment=assignment,
        root=GaussianRootPrior("fixed"),
        fit_spec=(
            BranchFitParameter("rate", "sigma2", ("baseline",), 1e-200, 1e-250, 1e250),
        ),
    )
    assert fit.estimation["groups"][0]["estimate"] == pytest.approx(expected, rel=2e-6)


@pytest.mark.parametrize(
    "gradient,coordinate", [(np.inf, 0.0), (-np.inf, 1.0), (np.nan, 0.5)]
)
def test_nonfinite_gradient_cannot_certify_boundary_convergence(gradient, coordinate):
    from nwkit.branch_gaussian_estimation import _projected_gradient

    assert math.isinf(_projected_gradient(np.array([gradient]), np.array([coordinate])))


@pytest.mark.parametrize("mode", ["gaussian", "flat"])
def test_fitting_is_invariant_to_large_trait_origin(mixed_inputs, mode):
    tree, assignment, spec, observed, errors = mixed_problem(mixed_inputs)
    root = GaussianRootPrior(mode, 0.0, None if mode == "flat" else 0.2)
    if mode == "flat":
        spec = tuple(p for p in spec if p.parameter != "alpha")
    _, original = compute_branch_marginals(
        tree, observed, errors, assignment=assignment, root=root, fit_spec=spec
    )
    offset = 1e10
    translated = replace(
        assignment,
        models_by_branch_id={
            key: replace(
                model,
                diffusion=replace(
                    model.diffusion, optimum=model.diffusion.optimum + offset
                ),
            )
            if isinstance(model.diffusion, OUBranch)
            else model
            for key, model in assignment.models_by_branch_id.items()
        },
    )
    translated_spec = tuple(
        replace(
            p,
            initial=p.initial + offset,
            lower=p.lower + offset,
            upper=p.upper + offset,
        )
        if p.parameter == "theta"
        else p
        for p in spec
    )
    _, fitted = compute_branch_marginals(
        tree,
        {
            name: None if value is None else value + offset
            for name, value in observed.items()
        },
        errors,
        assignment=translated,
        root=replace(root, mean=offset),
        fit_spec=translated_spec,
    )
    assert fitted.log_likelihood == pytest.approx(original.log_likelihood, abs=2e-5)
    for actual, expected in zip(
        fitted.estimation["groups"], original.estimation["groups"], strict=True
    ):
        estimate = actual["estimate"] - (
            offset if actual["parameter"] == "theta" else 0
        )
        assert estimate == pytest.approx(expected["estimate"], rel=2e-5, abs=2e-5)


def test_informative_initial_model_is_used_for_rank_check(mixed_inputs):
    tree, assignment, spec, observed, errors = mixed_problem(mixed_inputs)
    spec = tuple(replace(p, upper=1e100) if p.parameter == "alpha" else p for p in spec)
    _, fit = compute_branch_marginals(
        tree,
        observed,
        errors,
        assignment=assignment,
        root=GaussianRootPrior("gaussian", 0.0, 0.2),
        fit_spec=spec,
    )
    assert fit.log_likelihood == pytest.approx(-2.77486466116334, abs=2e-7)
    assert fit.estimation["parameter_rank"] == 3
