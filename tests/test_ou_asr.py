import math
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from scipy.optimize import minimize

from nwkit import ou_asr
from nwkit.cli import main
from nwkit.ou_asr import compute_ou_marginals
from nwkit.util import read_tree


def tree_from(source):
    return read_tree(source, "1", True, quiet=True, rooted="yes")


def ancestor_distances(node):
    result, distance = {node: 0.0}, 0.0
    while not node.is_root:
        distance += float(node.dist)
        node = node.up
        result[node] = distance
    return result


def patristic_distance(left, right):
    left_distances = ancestor_distances(left)
    right_distances = ancestor_distances(right)
    return min(
        left_distances[node] + right_distances[node]
        for node in left_distances.keys() & right_distances.keys()
    )


def dense_ou(tree, values, errors, alpha, sigma2, theta):
    nodes = list(tree.traverse())
    observed = [node for node in tree.leaves() if values.get(node.name) is not None]
    root_variance = sigma2 / (2.0 * alpha)
    covariance = np.array(
        [
            [
                root_variance * math.exp(-alpha * patristic_distance(left, right))
                for right in nodes
            ]
            for left in nodes
        ]
    )
    indices = [nodes.index(node) for node in observed]
    observed_covariance = covariance[np.ix_(indices, indices)] + np.diag(
        [errors.get(node.name, 0.0) ** 2 for node in observed]
    )
    response = np.array([values[node.name] for node in observed], dtype=float)
    centered = response - theta
    solved = np.linalg.solve(observed_covariance, centered)
    cross = covariance[:, indices]
    means = theta + cross @ solved
    conditional = covariance - cross @ np.linalg.solve(observed_covariance, cross.T)
    likelihood = -0.5 * (
        len(observed) * math.log(2.0 * math.pi)
        + np.linalg.slogdet(observed_covariance)[1]
        + centered @ solved
    )
    return nodes, means, np.diag(conditional), likelihood


def test_fixed_ou_polytomy_matches_independent_dense_covariance_oracle():
    tree = tree_from("((A:0.3,B:0.7)I:0.4,C:1.1,D:0.2)R;")
    values = {"A": 1.2, "B": -0.4, "C": 2.1, "D": None}
    errors = {"A": 0.2, "B": 0.0, "C": 0.4}
    alpha, sigma2, theta = 0.8, 1.7, 0.3
    posterior, fit = compute_ou_marginals(
        tree,
        values,
        alpha=alpha,
        sigma2=sigma2,
        theta=theta,
        standard_errors=errors,
    )
    nodes, means, variances, likelihood = dense_ou(
        tree, values, errors, alpha, sigma2, theta
    )
    assert fit.log_likelihood == pytest.approx(likelihood, abs=1e-12)
    for index, node in enumerate(nodes):
        assert posterior[node].mean == pytest.approx(means[index], abs=1e-12)
        assert posterior[node].variance == pytest.approx(
            max(0.0, variances[index]), abs=1e-12
        )


def test_fixed_alpha_sigma_estimates_dense_gls_theta():
    tree = tree_from("((A:0.2,B:0.6):0.3,C:0.9,D:1.2)R;")
    values = {"A": 0.1, "B": 1.7, "C": -0.3, "D": 0.8}
    errors = {name: 0.1 for name in values}
    alpha, sigma2 = 0.9, 1.3
    _, fit = compute_ou_marginals(
        tree,
        values,
        alpha=alpha,
        sigma2=sigma2,
        standard_errors=errors,
    )
    leaves = list(tree.leaves())
    variance = sigma2 / (2.0 * alpha)
    covariance = (
        np.array(
            [
                [variance * math.exp(-alpha * patristic_distance(a, b)) for b in leaves]
                for a in leaves
            ]
        )
        + np.eye(len(leaves)) * 0.01
    )
    inverse_one = np.linalg.solve(covariance, np.ones(len(leaves)))
    expected = (
        inverse_one
        @ np.array([values[node.name] for node in leaves])
        / sum(inverse_one)
    )
    assert fit.theta == pytest.approx(expected, rel=1e-7, abs=1e-8)
    assert fit.theta_estimated
    assert not fit.alpha_estimated
    assert not fit.sigma2_estimated


def test_profiled_theta_and_global_variance_search_find_boundary_optimum():
    tree = tree_from("[&R](A:0.2,B:0.7,C:1.3)R;")
    _, fit = compute_ou_marginals(
        tree,
        {"A": -3.0, "B": -4.8, "C": 0.2},
        alpha=1.0,
        standard_errors={"A": 0.6, "B": 0.0005, "C": 1.5},
    )
    assert fit.theta == pytest.approx(-4.7999981943, abs=2e-9)
    assert fit.root_variance == pytest.approx(fit.root_variance_bounds[0])
    assert fit.log_likelihood == pytest.approx(-5.1061516568, abs=2e-9)
    assert fit.fit_status == "root_variance_lower_boundary"
    assert fit.optimizer_success
    assert fit.optimizer_grid_evaluations == 49


def test_one_dimensional_ou_search_polishes_just_inside_a_bound():
    source = "[&R]((A:0.2,B:0.6):0.3,C:0.9,D:1.2)R;"
    values = {"A": 0.1, "B": 1.7, "C": -0.3, "D": 0.8}
    _, broad = compute_ou_marginals(
        tree_from(source),
        values,
        sigma2=1.3,
        alpha_bounds=(0.001, 100.0),
    )
    _, narrowed = compute_ou_marginals(
        tree_from(source),
        values,
        sigma2=1.3,
        alpha_bounds=(broad.alpha / 1.01, 100.0),
    )
    assert narrowed.alpha == pytest.approx(broad.alpha, rel=3e-7)
    assert narrowed.log_likelihood == pytest.approx(broad.log_likelihood, abs=1e-12)
    assert narrowed.fit_status == "ok"
    assert narrowed.optimizer_starts >= 1


def test_constant_exact_trait_returns_explicit_ou_boundary_fit():
    tree = tree_from("[&R]((A:1,B:1):1,C:2,D:1.5)R;")
    _, fit = compute_ou_marginals(tree, {"A": 2.0, "B": 2.0, "C": 2.0, "D": 2.0})
    assert fit.theta == 2.0
    assert fit.root_variance == pytest.approx(fit.root_variance_bounds[0])
    assert "root_variance_lower_boundary" in fit.fit_status
    assert fit.optimizer_success


@pytest.mark.parametrize("alpha", [350.0, 360.0, 400.0, 700.0, 745.0, 1000.0])
def test_strong_ou_attraction_uses_the_independent_branch_limit(alpha):
    tree = tree_from("[&R](A:1,B:1)R;")
    posterior, fit = compute_ou_marginals(
        tree,
        {"A": -1.0, "B": 1.0},
        alpha=alpha,
        sigma2=2.0 * alpha,
        theta=0.0,
    )
    expected = -math.log(2.0 * math.pi) - 1.0
    assert fit.log_likelihood == pytest.approx(expected, abs=2e-13)
    assert posterior[tree].mean == pytest.approx(0.0, abs=1e-300)
    assert posterior[tree].variance == pytest.approx(1.0)


def test_ou_observation_and_position_counts_have_distinct_meanings():
    tree = tree_from("[&R]((A:0,B:0)I:0,C:1)R;")
    _, fit = compute_ou_marginals(
        tree,
        {"A": 1.0, "B": 1.5, "C": 2.0},
        alpha=0.7,
        sigma2=1.2,
        theta=0.5,
        standard_errors={"A": 0.2, "B": 0.3, "C": 0.1},
    )
    assert fit.num_observed == 3
    assert fit.num_effective_observations == 3
    assert fit.num_observed_positions == 2


def test_ou_rejects_more_free_parameters_than_observed_positions():
    tree = tree_from("[&R](A:1,B:1)R;")
    with pytest.raises(ValueError, match="at least 3 distinct observed tree positions"):
        compute_ou_marginals(tree, {"A": 0.0, "B": 1.0})


def test_ou_rejects_standard_error_lost_during_power_of_two_scaling():
    tree = tree_from("[&R](A:1,B:1)R;")
    with pytest.raises(ValueError, match="OU standard error exceeds"):
        compute_ou_marginals(
            tree,
            {"A": 0.0, "B": 1e308},
            alpha=1.0,
            sigma2=1.0,
            theta=0.0,
            standard_errors={"A": math.ulp(0.0), "B": 1.0},
        )


def test_fitted_ou_matches_independent_dense_multistart_likelihood():
    tree = tree_from(
        "(((A:0.2,B:0.5):0.3,(C:0.4,D:0.7):0.2):0.4,"
        "((E:0.3,F:0.6):0.5,(G:0.2,H:0.8):0.3):0.2)R;"
    )
    values = dict(
        zip(
            "ABCDEFGH",
            [
                1.2043053758055426,
                1.0660051721718675,
                -0.7080006226350897,
                -0.016925310991197984,
                -0.6361758193419043,
                0.2614600133921077,
                -0.5070819256214758,
                -0.27089908894694165,
            ],
            strict=True,
        )
    )
    errors = {name: 0.1 + (index % 3) * 0.03 for index, name in enumerate("ABCDEFGH")}
    _, fit = compute_ou_marginals(
        tree,
        values,
        standard_errors=errors,
        alpha_bounds=(0.05, 5.0),
    )

    def objective(parameters):
        alpha = math.exp(parameters[0])
        root_variance = math.exp(parameters[1])
        sigma2 = 2.0 * alpha * root_variance
        return -dense_ou(tree, values, errors, alpha, sigma2, parameters[2])[3]

    starts = []
    initial_variance = np.var(list(values.values()))
    for alpha in np.geomspace(0.05, 5.0, 7):
        result = minimize(
            objective,
            [
                math.log(alpha),
                math.log(initial_variance),
                np.mean(list(values.values())),
            ],
            method="L-BFGS-B",
            bounds=[(math.log(0.05), math.log(5.0)), (-30, 30), (None, None)],
        )
        if result.success:
            starts.append(result)
    oracle = min(starts, key=lambda result: result.fun)
    assert fit.log_likelihood == pytest.approx(-oracle.fun, abs=2e-7)
    assert fit.alpha == pytest.approx(math.exp(oracle.x[0]), rel=2e-5)
    assert fit.root_variance == pytest.approx(math.exp(oracle.x[1]), rel=2e-5)
    assert fit.theta == pytest.approx(oracle.x[2], abs=2e-6)


def test_ou_multistart_skips_uncompetitive_powell_fallback(monkeypatch):
    results = iter(
        [
            SimpleNamespace(
                success=True, fun=1.0, x=np.array([0.0, 0.0]), message="ok"
            ),
            SimpleNamespace(
                success=False, fun=1.1, x=np.array([0.5, 0.5]), message="failed"
            ),
        ]
    )
    monkeypatch.setattr(
        ou_asr,
        "_two_dimensional_starts",
        lambda axes, objectives, initial: [(0, 0), (1, 1)],
    )
    monkeypatch.setattr(
        ou_asr,
        "_optimize_ou_start",
        lambda objective, start, bounds: (
            (result := next(results)),
            bool(result.success),
            [result.message],
        ),
    )

    def unexpected_powell(*args):
        raise AssertionError("uncompetitive failed endpoints must not run Powell")

    monkeypatch.setattr(ou_asr, "_powell_ou_start", unexpected_powell)
    candidates = []
    stats = ou_asr._search_ou_two_dimensions(
        [(-1.0, 1.0), (-1.0, 1.0)],
        [0.0, 0.0],
        lambda values: float(np.sum(np.square(values))),
        lambda *candidate: candidates.append(candidate),
    )
    assert stats == (2, 1, 1, 49)
    assert any(candidate[1] == "l-bfgs-b-endpoint" for candidate in candidates)


def test_ou_multistart_uses_powell_when_all_primary_starts_fail(monkeypatch):
    results = iter(
        [
            SimpleNamespace(
                success=False, fun=2.0, x=np.array([0.5, 0.5]), message="failed"
            ),
            SimpleNamespace(
                success=False, fun=3.0, x=np.array([0.8, 0.8]), message="failed"
            ),
        ]
    )
    monkeypatch.setattr(
        ou_asr,
        "_two_dimensional_starts",
        lambda axes, objectives, initial: [(0, 0), (1, 1)],
    )
    monkeypatch.setattr(
        ou_asr,
        "_optimize_ou_start",
        lambda objective, start, bounds: (
            (result := next(results)),
            False,
            [result.message],
        ),
    )
    fallback_calls = []

    def successful_powell(objective, start, bounds, messages):
        fallback_calls.append(tuple(start))
        return (
            SimpleNamespace(
                success=True, fun=1.0, x=np.array([0.0, 0.0]), message="ok"
            ),
            messages + ["Powell fallback: ok"],
        )

    monkeypatch.setattr(ou_asr, "_powell_ou_start", successful_powell)
    stats = ou_asr._search_ou_two_dimensions(
        [(-1.0, 1.0), (-1.0, 1.0)],
        [0.0, 0.0],
        lambda values: float(np.sum(np.square(values))),
        lambda *candidate: None,
    )
    assert stats == (2, 1, 1, 49)
    assert len(fallback_calls) == 1


def test_zero_length_nodes_share_the_same_ou_marginal():
    tree = tree_from("((A:0,B:0)I:0,C:1)R;")
    posterior, _ = compute_ou_marginals(
        tree,
        {"A": 2.0, "B": 2.0, "C": -1.0},
        alpha=0.7,
        sigma2=1.2,
        theta=0.4,
    )
    assert posterior[tree] == posterior[tree["I"]]
    assert posterior[tree] == posterior[tree["A"]]
    assert posterior[tree] == posterior[tree["B"]]
    assert posterior[tree].mean == 2.0
    assert posterior[tree].variance == 0.0


def test_ou_time_unit_change_preserves_fixed_model():
    first = tree_from("(A:0.3,B:0.8,C:1.1)R;")
    second = tree_from("(A:3,B:8,C:11)R;")
    values = {"A": -1.0, "B": 0.5, "C": 2.0}
    original, fit1 = compute_ou_marginals(
        first, values, alpha=0.6, sigma2=1.4, theta=0.2
    )
    rescaled, fit2 = compute_ou_marginals(
        second, values, alpha=0.06, sigma2=0.14, theta=0.2
    )
    assert fit1.log_likelihood == pytest.approx(fit2.log_likelihood, abs=1e-12)
    for name in ("A", "B", "C"):
        assert original[first[name]].mean == pytest.approx(
            rescaled[second[name]].mean, abs=1e-12
        )
        assert original[first[name]].variance == pytest.approx(
            rescaled[second[name]].variance, abs=1e-12
        )
    assert original[first].mean == pytest.approx(rescaled[second].mean, abs=1e-12)
    assert original[first].variance == pytest.approx(
        rescaled[second].variance, abs=1e-12
    )


def test_ou_trait_affine_change_has_expected_moments_and_density_jacobian():
    first = tree_from("((A:0.3,B:0.8):0.2,C:1.1)R;")
    second = tree_from("((A:0.3,B:0.8):0.2,C:1.1)R;")
    values = {"A": -1.0, "B": 0.5, "C": 2.0}
    errors = {"A": 0.1, "B": 0.2, "C": 0.3}
    scale, offset = 8.0, 1e12
    transformed_values = {
        name: offset + scale * value for name, value in values.items()
    }
    transformed_errors = {name: scale * error for name, error in errors.items()}
    original, fit1 = compute_ou_marginals(
        first,
        values,
        alpha=0.6,
        sigma2=1.4,
        theta=0.2,
        standard_errors=errors,
    )
    transformed, fit2 = compute_ou_marginals(
        second,
        transformed_values,
        alpha=0.6,
        sigma2=1.4 * scale**2,
        theta=offset + scale * 0.2,
        standard_errors=transformed_errors,
    )
    assert fit2.log_likelihood == pytest.approx(
        fit1.log_likelihood - len(values) * math.log(scale), abs=2e-6
    )
    for left, right in zip(first.traverse(), second.traverse(), strict=True):
        assert transformed[right].mean == pytest.approx(
            offset + scale * original[left].mean, abs=2e-4
        )
        assert transformed[right].variance == pytest.approx(
            scale**2 * original[left].variance, rel=2e-12
        )


def test_ou_power_of_two_scaling_handles_large_finite_trait_units():
    first = tree_from("(A:0.3,B:0.8,C:1.1)R;")
    second = tree_from("(A:0.3,B:0.8,C:1.1)R;")
    base, unit = math.ldexp(1.0, 500), math.ldexp(1.0, 452)
    requested_values = {"A": -1.0, "B": 0.5, "C": 2.0}
    large_values = {
        name: base + unit * value for name, value in requested_values.items()
    }
    small_values = {name: (value - base) / unit for name, value in large_values.items()}
    large_theta = base + unit * 0.2
    small_theta = (large_theta - base) / unit
    small, fit1 = compute_ou_marginals(
        first, small_values, alpha=0.6, sigma2=1.4, theta=small_theta
    )
    large, fit2 = compute_ou_marginals(
        second,
        large_values,
        alpha=0.6,
        sigma2=math.ldexp(1.4, 904),
        theta=large_theta,
    )
    assert math.isfinite(fit2.log_likelihood)
    for left, right in zip(first.traverse(), second.traverse(), strict=True):
        normalized_mean = (large[right].mean - base) / unit
        normalized_variance = math.ldexp(large[right].variance, -904)
        assert normalized_mean == pytest.approx(small[left].mean, abs=0.04)
        assert normalized_variance == pytest.approx(small[left].variance, rel=1e-12)
    assert fit2.log_likelihood == pytest.approx(
        fit1.log_likelihood - len(small_values) * math.log(unit), abs=1e-10
    )


def test_ou_requires_positive_alpha_and_sigma2():
    tree = tree_from("(A:1,B:1)R;")
    values = {"A": 0.0, "B": 1.0}
    with pytest.raises(ValueError, match="--alpha must be strictly positive"):
        compute_ou_marginals(tree, values, alpha=0, sigma2=1, theta=0)
    with pytest.raises(ValueError, match="--sigma2 must be strictly positive"):
        compute_ou_marginals(tree, values, alpha=1, sigma2=0, theta=0)


def test_ou_rejects_joint_parameter_fit_at_one_observed_position():
    tree = tree_from("(A:1,B:1)R;")
    with pytest.raises(ValueError, match="not separately identifiable"):
        compute_ou_marginals(tree, {"A": 1.0, "B": None}, theta=0.0)
    posterior, fit = compute_ou_marginals(
        tree,
        {"A": 1.0, "B": None},
        alpha=0.5,
        sigma2=1.0,
        theta=0.0,
    )
    assert fit.log_likelihood == pytest.approx(-0.5 * (math.log(2 * math.pi) + 1))
    assert math.isfinite(posterior[tree["B"]].variance)


def test_ou_cli_fixed_and_fitted_metadata(tmp_path):
    trait = tmp_path / "trait.tsv"
    trait.write_text("leaf_name\tvalue\nA\t1.2\nB\t-0.4\nC\t2.1\nD\t0.8\nE\t1.4\n")
    source = "[&R]((A:0.3,B:0.7):0.4,(C:0.4,D:0.8):0.3,E:1.1)R;"
    for label, options, estimated in (
        ("fixed", ["--alpha", "0.8", "--sigma2", "1.7", "--theta", "0.3"], False),
        ("fitted", ["--alpha-bounds", "0.01,10"], True),
    ):
        output = tmp_path / f"{label}.tsv"
        model = tmp_path / f"{label}_model.tsv"
        main(
            [
                "asr",
                "-i",
                source,
                "--trait",
                str(trait),
                "--state-column",
                "value",
                "--model",
                "OU",
                "--model-out",
                str(model),
                "-o",
                str(output),
                *options,
            ]
        )
        metadata = pd.read_csv(model, sep="\t").iloc[0]
        assert metadata.model == "OU"
        assert metadata.root_prior == "stationary"
        assert bool(metadata.alpha_estimated) is estimated
        assert bool(metadata.sigma2_estimated) is estimated
        assert bool(metadata.theta_estimated) is estimated
        assert metadata.likelihood_kind == "stationary_root_ml"
        assert metadata.interval_kind == "conditional_on_parameters"
        assert metadata.num_effective_observations == 5
        assert metadata.num_observed_positions == 5
        assert metadata.optimizer_starts >= 0
        assert metadata.optimizer_converged_starts <= metadata.optimizer_starts
        assert metadata.optimizer_failed_starts <= metadata.optimizer_starts
        assert np.isfinite(metadata.log_likelihood)
        assert {"mean", "variance", "ci_lower", "ci_upper"}.issubset(
            pd.read_csv(output, sep="\t").columns
        )


def test_ou_cli_rejects_alpha_bounds_when_alpha_is_fixed(tmp_path):
    trait = tmp_path / "trait.tsv"
    trait.write_text("leaf_name\tvalue\nA\t0\nB\t1\n")
    with pytest.raises(ValueError, match="--alpha-bounds cannot be combined"):
        main(
            [
                "asr",
                "-i",
                "[&R](A:1,B:1)R;",
                "--trait",
                str(trait),
                "--state-column",
                "value",
                "--model",
                "OU",
                "--alpha",
                "0.5",
                "--alpha-bounds",
                "0.01,1",
            ]
        )
