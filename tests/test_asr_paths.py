"""Check simulation laws, continuity, finite-grid geometry and CLI rendering."""

import math
from types import SimpleNamespace

import numpy as np
import pytest
from ete4 import Tree

from nwkit.asr_paths import _grid_counts, _refined_tree, simulate_fitted_paths
from nwkit.asr_regimes import RegimeAssignment
from nwkit.cli import main
from nwkit.continuous_asr import GaussianMarginal, compute_bm_marginals
from nwkit.continuous_asr_process import fitted_scalar_process
from nwkit.rooting_state import set_rooting_info
from tests.test_asr_figure import _command


def rooted(source):
    tree = Tree(source)
    set_rooting_info(tree, True)
    return tree


def ou_fit():
    return SimpleNamespace(
        alpha=0.8,
        sigma2=1.2,
        theta=2.0,
        root_prior="fixed",
        root_mean=0.0,
        root_variance=0.0,
    )


def test_unconditional_ou_endpoint_moments():
    tree = rooted("(A:2);")
    fit = ou_fit()
    result = simulate_fitted_paths(
        tree,
        {"A": 99.0},
        None,
        {tree: GaussianMarginal(0.0, 0.0)},
        fit=fit,
        model="OU",
        count=16000,
        steps=4,
        seed=321,
    )
    path = result.branches[next(tree.leaves())]
    values = path.values[:, -1, 0]
    decay = math.exp(-fit.alpha * 2)
    expected_mean = fit.theta * (1 - decay)
    expected_variance = fit.sigma2 / (2 * fit.alpha) * (1 - decay**2)
    assert abs(values.mean() - expected_mean) < 6 * math.sqrt(
        expected_variance / len(values)
    )
    assert abs(values.var() - expected_variance) < 0.04
    assert np.all(result.root_values == 0)
    assert not np.any(values == 99.0)


@pytest.mark.parametrize("error", [0.0, 0.3])
def test_conditional_ou_midpoint_matches_gaussian_bridge(error):
    tree = rooted("(A:2);")
    fit = ou_fit()
    observation = 3.0
    paths = simulate_fitted_paths(
        tree,
        {"A": observation},
        {"A": error},
        {tree: GaussianMarginal(0.0, 0.0)},
        fit=fit,
        model="OU",
        count=10000,
        steps=2,
        mode="conditional",
        seed=71,
    )
    path = paths.branches[next(tree.leaves())]
    decay = math.exp(-fit.alpha)
    variance = fit.sigma2 / (2 * fit.alpha)
    m1, m2 = fit.theta * (1 - decay), fit.theta * (1 - decay**2)
    v1, v2 = variance * (1 - decay**2), variance * (1 - decay**4)
    cross = decay * v1
    mean = m1 + cross / (v2 + error**2) * (observation - m2)
    conditional_variance = v1 - cross**2 / (v2 + error**2)
    middle = path.values[:, 1, 0]
    assert abs(middle.mean() - mean) < 6 * math.sqrt(conditional_variance / len(middle))
    assert abs(middle.var() - conditional_variance) < 0.035
    if not error:
        np.testing.assert_allclose(path.values[:, -1, 0], observation, atol=1e-13)
    else:
        assert path.values[:, -1, 0].var() > 0.03


@pytest.mark.parametrize("mode", ["unconditional", "conditional"])
def test_paths_are_joint_reproducible_and_do_not_modify_the_input(mode):
    tree = rooted("((A:1,B:0.5):1,C:0,D:0.75):80;")
    observed = {"A": 1.0, "B": 2.0, "C": 0.0, "D": None}
    posterior, fit = compute_bm_marginals(tree, observed, sigma2=0.5)
    before = [(node, node.up, dict(node.props)) for node in tree.traverse()]
    arguments = dict(fit=fit, model="BM", count=3, steps=12, mode=mode, seed=84)
    first = simulate_fitted_paths(tree, observed, None, posterior, **arguments)
    second = simulate_fitted_paths(tree, observed, None, posterior, **arguments)
    assert [(node, node.up, dict(node.props)) for node in tree.traverse()] == before
    for node, path in first.branches.items():
        assert path.elapsed[0] == 0
        assert path.elapsed[-1] == node.dist
        parent = (
            first.root_values
            if node.up is tree
            else first.branches[node.up].values[:, -1, :]
        )
        np.testing.assert_array_equal(path.values[:, 0, :], parent)
        np.testing.assert_array_equal(path.values, second.branches[node].values)
        if node.dist == 0:
            np.testing.assert_array_equal(path.values[:, -1, :], parent)
    if mode == "conditional":
        for node in tree.leaves():
            if observed[node.name] is not None:
                np.testing.assert_allclose(
                    first.branches[node].values[:, -1, 0],
                    observed[node.name],
                    atol=1e-12,
                )
    else:
        assert "flat root" in first.root_description


@pytest.mark.parametrize(
    "model",
    ["BM", "BM-DRIFT", "BMS", "BMS-DRIFT", "OU", "OUM", "OUMA", "OUMV", "OUMVA"],
)
def test_refinement_preserves_original_scalar_transitions(model):
    tree = rooted("((A:1,B:0.5):0.8,C:0);")
    by_node = {node: "high" if node.name == "B" else "low" for node in tree.traverse()}
    assignment = RegimeAssignment(("low", "high"), by_node, "test")
    fit = ou_fit()
    fit.drift = -0.3
    fit.sigma2_by_regime = {"low": 0.7, "high": 1.4}
    fit.drift_by_regime = {"low": -0.3, "high": 0.6}
    fit.alpha_by_regime = {"low": 0.8, "high": 1.7}
    fit.theta_by_regime = {"low": 2.0, "high": -1.0}
    process = fitted_scalar_process(tree, model, fit, regime_assignment=assignment)
    refined, chains, mapped = _refined_tree(
        tree, _grid_counts(tree, 13, 1, 1), assignment
    )
    expanded = fitted_scalar_process(refined, model, fit, regime_assignment=mapped)
    for node, (chain, _) in chains.items():
        slope, intercept, variance = 1.0, 0.0, 0.0
        for child in chain[1:]:
            edge = expanded.transitions[child]
            slope, intercept, variance = (
                edge.slope * slope,
                edge.slope * intercept + edge.intercept,
                edge.slope**2 * variance + edge.variance,
            )
        expected = process.transitions[node]
        np.testing.assert_allclose(
            [slope, intercept, variance],
            [expected.slope, expected.intercept, expected.variance],
            rtol=1e-12,
            atol=1e-14,
        )


def test_vector_paths_retain_trait_covariance():
    tree = rooted("(A:2);")
    sigma = np.array([[1.0, 0.65], [0.65, 2.0]])
    fit = SimpleNamespace(trait_names=("x", "y"), sigma=sigma)
    result = simulate_fitted_paths(
        tree,
        {"A": [1.0, 2.0]},
        None,
        {tree: SimpleNamespace(mean=np.zeros(2))},
        fit=fit,
        model="MV-BM",
        count=12000,
        steps=4,
        seed=121,
    )
    values = result.branches[next(tree.leaves())].values[:, -1, :]
    np.testing.assert_allclose(np.cov(values.T), 2 * sigma, atol=0.09)


@pytest.mark.parametrize("mode", ["unconditional", "conditional"])
def test_simulation_cli_adds_panel_with_shared_scales_and_keeps_tsv(
    tmp_path, monkeypatch, mode
):
    import nwkit.asr_figure as plotting

    command = _command(tmp_path, "--sigma2", "0.5")
    main(command)
    before = (tmp_path / "nodes.tsv").read_bytes()
    captured = []
    original = plotting.build_continuous_asr_figure

    def capture(tree, table, **kwargs):
        figure = original(tree, table, **kwargs)
        assert len(figure.axes) == 3
        assert len({ax.get_ylim() for ax in figure.axes}) == 1
        assert figure.axes[1].get_xlim() == figure.axes[2].get_xlim()
        captured.append(kwargs["simulation"])
        return figure

    monkeypatch.setattr(plotting, "build_continuous_asr_figure", capture)
    main(
        [
            *command,
            "--figure-out",
            str(tmp_path / "paths.svg"),
            "--figure-simulations",
            "2",
            "--figure-simulation-mode",
            mode,
            "--figure-simulation-steps",
            "8",
            "--seed",
            "11",
        ]
    )
    assert captured[0].count == 2
    assert captured[0].mode == mode
    assert (tmp_path / "nodes.tsv").read_bytes() == before
    svg = (tmp_path / "paths.svg").read_text()
    assert "simulation" in svg if mode == "unconditional" else "posterior paths" in svg


@pytest.mark.parametrize(
    "options,match",
    [
        (["--figure-simulations", "-1"], "integer"),
        (["--figure-simulations", "1"], "requires --figure-out"),
        (["--figure-simulation-steps", "0"], "integer"),
        (["--figure-simulation-steps", "20"], "require --figure-simulations"),
        (["--figure-simulation-mode", "conditional"], "require --figure-simulations"),
        (
            [
                "--figure-out",
                "plot.pdf",
                "--figure-simulations",
                "1",
                "--model",
                "KAPPA",
            ],
            "does not support",
        ),
    ],
)
def test_simulation_options_fail_before_output(tmp_path, options, match):
    with pytest.raises(ValueError, match=match):
        main(_command(tmp_path, *options))
    assert not (tmp_path / "nodes.tsv").exists()


def test_resource_guard_precedes_grid_allocation():
    tree = rooted("(A:1,B:1);")
    with pytest.raises(ValueError, match="grid is too large"):
        _grid_counts(tree, 100, 10000, 2)
    with pytest.raises(ValueError, match="grid is too large"):
        _grid_counts(tree, 10**500, 1, 1)
    with pytest.raises(ValueError, match="branch traces"):
        _grid_counts(tree, 1, 30000, 1)


@pytest.mark.parametrize("model", ["MV-BM", "MV-OU", "MV-OU-DIAG", "MV-OU-FULL"])
def test_refinement_preserves_vector_transitions(model):
    from nwkit.asr_multivariate_diagnostics import fitted_vector_process

    tree = rooted("(A:1.7,B:0);")
    fit = SimpleNamespace(
        trait_names=("x", "y"),
        sigma=np.array([[1.0, 0.4], [0.4, 1.5]]),
        alpha=0.7,
        alpha_by_trait=np.array([0.5, 1.2]),
        attraction_matrix=np.array([[0.7, -0.2], [0.3, 1.1]]),
        diffusion_sigma=np.array([[1.0, 0.4], [0.4, 1.5]]),
        theta=np.array([1.0, -2.0]),
    )
    process = fitted_vector_process(tree, model, fit)
    refined, chains, _ = _refined_tree(tree, _grid_counts(tree, 12, 1, 2), None)
    expanded = fitted_vector_process(refined, model, fit)
    for node, (chain, _) in chains.items():
        slope, intercept, covariance = np.eye(2), np.zeros(2), np.zeros((2, 2))
        for point in chain[1:]:
            edge = expanded.transitions[point]
            slope = edge.slope @ slope
            intercept = edge.slope @ intercept + edge.intercept
            covariance = edge.slope @ covariance @ edge.slope.T + edge.covariance
        expected = process.transitions[node]
        np.testing.assert_allclose(slope, expected.slope, atol=1e-12)
        np.testing.assert_allclose(intercept, expected.intercept, atol=1e-12)
        np.testing.assert_allclose(covariance, expected.covariance, atol=1e-12)


@pytest.mark.parametrize("model", ["MV-BM", "MV-OU-FULL"])
@pytest.mark.parametrize("mode", ["unconditional", "conditional"])
def test_multivariate_simulation_cli(tmp_path, monkeypatch, mode, model):
    import nwkit.asr_figure as plotting

    traits = tmp_path / "multi.tsv"
    traits.write_text(
        "leaf_name\tx\ty\nA\t0.1\t2.4\nB\t0.8\t1.7\nC\t1.2\t3.2\nD\t2.1\t2.1\nE\t-0.7\t0.4\nF\t-0.5\t2.0\n"
    )
    original = plotting.build_continuous_asr_figure
    captured = []

    def capture(tree, table, **kwargs):
        figure = original(tree, table, **kwargs)
        assert len(figure.axes) == 5
        simulation = kwargs["simulation"]
        assert simulation.root_values.shape == (2, 2)
        if mode == "conditional":
            for node in tree.leaves():
                values = table[table.name == node.name]["observed_value"].to_numpy(
                    dtype=float
                )
                np.testing.assert_allclose(
                    simulation.branches[node].values[0, -1, :], values, atol=1e-10
                )
        captured.append(simulation)
        return figure

    monkeypatch.setattr(plotting, "build_continuous_asr_figure", capture)
    main(
        [
            "asr",
            "-i",
            "[&R]((A:1,B:1):1,(C:1,D:1):1,(E:1,F:1):1);",
            "--trait",
            str(traits),
            "--state-column",
            "x,y",
            "--model",
            model,
            *(
                [
                    "--attraction-matrix",
                    "0.7,-0.2;0.3,1.1",
                    "--diffusion-matrix",
                    "1,0.3;0.3,1.2",
                ]
                if model == "MV-OU-FULL"
                else []
            ),
            "-o",
            str(tmp_path / "multi.tsv.out"),
            "--figure-out",
            str(tmp_path / "multi.svg"),
            "--figure-simulations",
            "2",
            "--figure-simulation-steps",
            "4",
            "--figure-simulation-mode",
            mode,
            "--seed",
            "34",
        ]
    )
    assert len(captured) == 1


def test_correlated_measurement_error_paths_match_node_moments():
    from nwkit.full_ou_fit import fit_full_mvou

    tree = rooted("(A:1,B:2,C:3,D:4);")
    observed = {"A": [0.1, 2.4], "B": [0.8, 1.7], "C": [1.2, 3.2], "D": [2.1, 2.1]}
    covariance = np.array([[0.5, 0.35], [0.35, 0.7]])
    posterior, fit = fit_full_mvou(
        tree,
        observed,
        ("x", "y"),
        attraction=np.array([[0.7, -0.2], [0.3, 1.1]]),
        diffusion=np.array([[1.0, 0.3], [0.3, 1.2]]),
        measurement_covariances=dict.fromkeys(observed, covariance),
    )
    paths = simulate_fitted_paths(
        tree,
        observed,
        None,
        posterior,
        fit=fit,
        model="MV-OU-FULL",
        count=8000,
        steps=1,
        mode="conditional",
        seed=67,
    )
    tip = next(tree.leaves())
    samples = paths.branches[tip].values[:, -1, :]
    np.testing.assert_allclose(samples.mean(axis=0), posterior[tip].mean, atol=0.025)
    np.testing.assert_allclose(np.cov(samples.T), posterior[tip].covariance, atol=0.025)
