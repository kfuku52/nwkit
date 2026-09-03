import math

import numpy as np
import pandas as pd
import pytest

from nwkit.asr_regimes import RegimeAssignment
from nwkit.cli import main
from nwkit.continuous_asr import compute_bm_marginals
from nwkit.ou_asr import compute_ou_marginals
from nwkit.regime_continuous_asr import (
    compute_bms_marginals,
    compute_oum_marginals,
)
from nwkit.util import read_tree


def tree_from(source):
    return read_tree(source, "1", True, quiet=True, rooted="yes")


def assignment(tree, values):
    nodes = list(tree.traverse())
    by_node = dict(zip(nodes, values, strict=True))
    regimes = tuple(dict.fromkeys(values))
    return RegimeAssignment(regimes, by_node, "memory")


def path_to_root(node):
    result = {node: 0.0}
    distance = 0.0
    while not node.is_root:
        distance += float(node.dist)
        node = node.up
        result[node] = distance
    return result


def patristic_distance(first, second):
    first_path = path_to_root(first)
    second_path = path_to_root(second)
    common = next(node for node in first_path if node in second_path)
    return first_path[common] + second_path[common]


def dense_oum(tree, observed, alpha, sigma2, theta_by_regime, regimes):
    nodes = list(tree.traverse())
    index = {node: position for position, node in enumerate(nodes)}
    root_variance = sigma2 / (2 * alpha)
    means = np.zeros(len(nodes))
    means[0] = theta_by_regime[regimes.by_node[tree]]
    for node in nodes[1:]:
        decay = math.exp(-alpha * float(node.dist))
        theta = theta_by_regime[regimes.by_node[node]]
        means[index[node]] = decay * means[index[node.up]] + (1 - decay) * theta
    covariance = np.asarray(
        [
            [
                root_variance * math.exp(-alpha * patristic_distance(left, right))
                for right in nodes
            ]
            for left in nodes
        ]
    )
    observed_nodes = [
        tree[name] for name, value in observed.items() if value is not None
    ]
    observed_indices = [index[node] for node in observed_nodes]
    response = np.asarray([observed[node.name] for node in observed_nodes])
    observed_covariance = covariance[np.ix_(observed_indices, observed_indices)]
    cross = covariance[:, observed_indices]
    solved = np.linalg.solve(observed_covariance, response - means[observed_indices])
    conditional_mean = means + cross @ solved
    conditional_covariance = covariance - cross @ np.linalg.solve(
        observed_covariance, cross.T
    )
    residual = response - means[observed_indices]
    sign, log_determinant = np.linalg.slogdet(observed_covariance)
    assert sign == 1
    log_likelihood = -0.5 * (
        len(response) * math.log(2 * math.pi)
        + log_determinant
        + residual @ np.linalg.solve(observed_covariance, residual)
    )
    return {
        node: (
            conditional_mean[index[node]],
            conditional_covariance[index[node], index[node]],
        )
        for node in nodes
    }, log_likelihood


def test_fixed_equal_bms_is_exactly_the_existing_bm_model():
    tree = tree_from("((A:0.4,B:1.2)I:0.7,C:1.5,D:0.8)R;")
    values = {"A": -1.0, "B": 2.0, "C": 4.0, "D": None}
    errors = {"A": 0.2, "B": 0.3, "C": 0.4}
    regimes = assignment(tree, ["r0", "r0", "r1", "r1", "r0", "r1"])
    expected, expected_fit = compute_bm_marginals(
        tree, values, sigma2=0.8, standard_errors=errors
    )
    observed, fit = compute_bms_marginals(
        tree,
        values,
        regimes,
        sigma2_by_regime={"r0": 0.8, "r1": 0.8},
        standard_errors=errors,
    )
    assert fit.restricted_log_likelihood == pytest.approx(
        expected_fit.restricted_log_likelihood, abs=1e-12
    )
    for node in tree.traverse():
        assert observed[node] == pytest.approx(expected[node])


def test_fixed_bms_star_uses_each_incoming_branch_rate():
    tree = tree_from("(A:1,B:1,C:2)R;")
    regimes = assignment(tree, ["slow", "slow", "slow", "fast"])
    posterior, fit = compute_bms_marginals(
        tree,
        {"A": 0.0, "B": 2.0, "C": None},
        regimes,
        sigma2_by_regime={"slow": 1.0, "fast": 4.0},
    )
    assert posterior[tree].mean == pytest.approx(1.0)
    assert posterior[tree].variance == pytest.approx(0.5)
    assert posterior[tree["C"]].mean == pytest.approx(1.0)
    assert posterior[tree["C"]].variance == pytest.approx(8.5)
    assert fit.sigma2_by_regime == {"slow": 1.0, "fast": 4.0}


def test_fixed_zero_bms_preserves_singular_boundary_semantics():
    tree = tree_from("((A:1,B:2):1,C:2,D:1)R;")
    values = {"A": 3.0, "B": 3.0, "C": 3.0, "D": None}
    regimes = assignment(tree, ["r0", "r0", "r1", "r1", "r0", "r1"])
    expected, expected_fit = compute_bm_marginals(tree, values, sigma2=0.0)
    observed, fit = compute_bms_marginals(
        tree,
        values,
        regimes,
        sigma2_by_regime={"r0": 0.0, "r1": 0.0},
    )
    assert fit.fit_status == expected_fit.fit_status == "singular_zero_boundary"
    assert (
        fit.restricted_log_likelihood is expected_fit.restricted_log_likelihood is None
    )
    assert fit.residual_df == expected_fit.residual_df
    for node in tree.traverse():
        assert observed[node] == expected[node]


def test_estimated_single_regime_bms_matches_existing_bm_fit():
    tree = tree_from("((A:0.4,B:1.2)I:0.7,C:1.5,D:0.8,E:2)R;")
    values = {"A": -1.0, "B": 2.0, "C": 4.0, "D": 1.0, "E": 3.0}
    regimes = assignment(tree, ["shared"] * len(list(tree.traverse())))
    expected, expected_fit = compute_bm_marginals(tree, values)
    observed, fit = compute_bms_marginals(tree, values, regimes)
    assert fit.sigma2_by_regime["shared"] == pytest.approx(
        expected_fit.sigma2, rel=1e-6
    )
    assert fit.restricted_log_likelihood == pytest.approx(
        expected_fit.restricted_log_likelihood, abs=1e-10
    )
    for node in tree.traverse():
        assert observed[node].mean == pytest.approx(expected[node].mean, abs=1e-10)
        assert observed[node].variance == pytest.approx(
            expected[node].variance, rel=1e-6, abs=1e-10
        )


def test_bms_rejects_rate_on_stem_common_to_every_observation():
    tree = tree_from("((A:1,B:1,C:1,D:1)I:2,E:1)R;")
    regimes = assignment(
        tree,
        [
            "terminal",
            "stem",
            "terminal",
            "terminal",
            "terminal",
            "terminal",
            "terminal",
        ],
    )
    with pytest.raises(ValueError, match="not separately identifiable"):
        compute_bms_marginals(
            tree,
            {"A": 0.0, "B": 1.0, "C": 2.0, "D": 3.0, "E": None},
            regimes,
        )


def test_equal_optimum_oum_matches_single_optimum_ou():
    tree = tree_from("((A:0.4,B:1.2)I:0.7,C:1.5,D:0.8)R;")
    values = {"A": -1.0, "B": 2.0, "C": 4.0, "D": None}
    errors = {"A": 0.2, "B": 0.3, "C": 0.4}
    regimes = assignment(tree, ["r0", "r0", "r1", "r1", "r0", "r1"])
    expected, expected_fit = compute_ou_marginals(
        tree,
        values,
        alpha=0.7,
        sigma2=1.4,
        theta=0.3,
        standard_errors=errors,
    )
    observed, fit = compute_oum_marginals(
        tree,
        values,
        regimes,
        alpha=0.7,
        sigma2=1.4,
        theta_by_regime={"r0": 0.3, "r1": 0.3},
        standard_errors=errors,
    )
    assert fit.log_likelihood == pytest.approx(expected_fit.log_likelihood, abs=1e-12)
    for node in tree.traverse():
        assert observed[node].mean == pytest.approx(expected[node].mean, abs=1e-12)
        assert observed[node].variance == pytest.approx(
            expected[node].variance, abs=1e-12
        )


def test_estimated_single_regime_theta_matches_existing_ou_fit():
    tree = tree_from("((A:0.4,B:1.2)I:0.7,C:1.5,D:0.8,E:2)R;")
    values = {"A": -1.0, "B": 2.0, "C": 4.0, "D": 1.0, "E": 3.0}
    regimes = assignment(tree, ["shared"] * len(list(tree.traverse())))
    expected, expected_fit = compute_ou_marginals(tree, values, alpha=0.7, sigma2=1.4)
    observed, fit = compute_oum_marginals(tree, values, regimes, alpha=0.7, sigma2=1.4)
    assert fit.theta_by_regime["shared"] == pytest.approx(expected_fit.theta, abs=2e-7)
    assert fit.optimizer_starts == expected_fit.optimizer_starts == 0
    assert fit.log_likelihood == pytest.approx(expected_fit.log_likelihood, abs=1e-10)
    for node in tree.traverse():
        assert observed[node].mean == pytest.approx(expected[node].mean, abs=2e-7)
        assert observed[node].variance == pytest.approx(
            expected[node].variance, abs=2e-12
        )


def test_estimated_single_regime_oum_covariance_matches_existing_ou_fit():
    tree = tree_from("((A:0.4,B:1.2)I:0.7,C:1.5,D:0.8,E:2)R;")
    values = {"A": -1.0, "B": 2.0, "C": 4.0, "D": 1.0, "E": 3.0}
    regimes = assignment(tree, ["shared"] * len(list(tree.traverse())))
    _, expected = compute_ou_marginals(tree, values, theta=0.3)
    _, observed = compute_oum_marginals(
        tree, values, regimes, theta_by_regime={"shared": 0.3}
    )
    assert observed.alpha == pytest.approx(expected.alpha, rel=2e-6)
    assert observed.sigma2 == pytest.approx(expected.sigma2, rel=2e-6)
    assert observed.log_likelihood == pytest.approx(expected.log_likelihood, abs=1e-9)
    assert observed.fit_status == expected.fit_status


def test_oum_rejects_collinear_regime_optimum_design():
    tree = tree_from("(A:1,B:1,C:1,D:1)R;")
    regimes = assignment(
        tree, ["root_only", "terminal", "terminal", "terminal", "terminal"]
    )
    with pytest.raises(ValueError, match="not separately identifiable"):
        compute_oum_marginals(
            tree,
            {"A": 0.0, "B": 1.0, "C": 2.0, "D": 3.0},
            regimes,
            alpha=0.7,
            sigma2=1.2,
        )


def test_fixed_multi_optimum_oum_matches_independent_dense_conditioning():
    tree = tree_from("((A:0.4,B:1.2)I:0.7,C:1.5,D:0.8)R;")
    values = {"A": -1.0, "B": 2.0, "C": 4.0, "D": None}
    regimes = assignment(tree, ["cold", "cold", "warm", "cold", "warm", "warm"])
    theta = {"cold": -0.5, "warm": 2.0}
    posterior, fit = compute_oum_marginals(
        tree,
        values,
        regimes,
        alpha=0.8,
        sigma2=1.6,
        theta_by_regime=theta,
    )
    expected, expected_log_likelihood = dense_oum(
        tree, values, 0.8, 1.6, theta, regimes
    )
    assert fit.log_likelihood == pytest.approx(expected_log_likelihood, abs=2e-12)
    for node, marginal in posterior.items():
        assert (marginal.mean, marginal.variance) == pytest.approx(
            expected[node], abs=2e-12
        )


@pytest.mark.integration
@pytest.mark.parametrize("model", ["BMS", "OUM"])
def test_regime_continuous_cli_fixed_parameters(model, tmp_path):
    traits = tmp_path / "traits.tsv"
    regime_map = tmp_path / "regimes.tsv"
    parameters = tmp_path / "parameters.tsv"
    output = tmp_path / "asr.tsv"
    metadata = tmp_path / "model.tsv"
    traits.write_text("leaf_name\tvalue\nA\t0\nB\t1\nC\t3\nD\t4\n")
    regime_map.write_text(
        "branch_id\tregime\n0\tlow\n1\tlow\n2\thigh\n3\tlow\n4\tlow\n5\thigh\n6\thigh\n"
    )
    column = "sigma2" if model == "BMS" else "theta"
    first, second = (0.5, 2.0) if model == "BMS" else (0.0, 3.0)
    parameters.write_text(f"regime\t{column}\nlow\t{first}\nhigh\t{second}\n")
    options = [] if model == "BMS" else ["--alpha", "0.7", "--sigma2", "1.2"]
    main(
        [
            "asr",
            "-i",
            "[&R]((A:1,B:1):1,(C:1,D:1):1)R;",
            "--trait",
            str(traits),
            "--state-column",
            "value",
            "--model",
            model,
            "--regime-map",
            str(regime_map),
            "--regime-parameters",
            str(parameters),
            "--model-out",
            str(metadata),
            "-o",
            str(output),
            *options,
        ]
    )
    result = pd.read_csv(output, sep="\t")
    fit = pd.read_csv(metadata, sep="\t").iloc[0]
    assert np.all(np.isfinite(result[["mean", "variance"]]))
    assert fit["model"] == model
    assert fit["regimes"] == "low,high"
    assert fit[f"{column}_low"] == first
    assert fit[f"{column}_high"] == second


@pytest.mark.integration
def test_regime_name_cannot_overwrite_model_metadata_columns(tmp_path):
    traits = tmp_path / "traits.tsv"
    regime_map = tmp_path / "regimes.tsv"
    parameters = tmp_path / "parameters.tsv"
    output = tmp_path / "asr.tsv"
    metadata = tmp_path / "model.tsv"
    traits.write_text("leaf_name\tvalue\nA\t0\nB\t1\nC\t3\n")
    regime_map.write_text(
        "branch_id\tregime\n0\testimated\n1\testimated\n2\testimated\n3\testimated\n"
    )
    parameters.write_text("regime\tsigma2\nestimated\t0.5\n")
    main(
        [
            "asr",
            "-i",
            "[&R](A:1,B:1,C:1)R;",
            "--trait",
            str(traits),
            "--state-column",
            "value",
            "--model",
            "BMS",
            "--regime-map",
            str(regime_map),
            "--regime-parameters",
            str(parameters),
            "--model-out",
            str(metadata),
            "-o",
            str(output),
        ]
    )
    fit = pd.read_csv(metadata, sep="\t").iloc[0]
    assert not fit["sigma2_estimated"]
    escaped = "regime_" + "estimated".encode().hex()
    assert fit[f"sigma2_{escaped}"] == 0.5
