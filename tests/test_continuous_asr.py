"""Independent Gaussian oracles for continuous ancestral reconstruction."""

import math

import numpy as np
import pytest
from scipy.linalg import eigh
from scipy.optimize import brentq, minimize_scalar

import nwkit.continuous_asr as continuous
from nwkit.continuous_asr import compute_bm_marginals
from nwkit.util import read_tree
from tests.helpers import make_deep_ladder_tree


def tree_from(text):
    return read_tree(text, "1", True, rooted="yes", quiet=True)


def dense_marginals(tree, values, sigma2, errors=None):
    """Invert the full graph Laplacian after conditioning on exact tips."""
    nodes = list(tree.traverse())
    index = {node: i for i, node in enumerate(nodes)}
    precision = np.zeros((len(nodes), len(nodes)))
    linear = np.zeros(len(nodes))
    locked = {}
    for node in nodes:
        i = index[node]
        if not node.is_root:
            j = index[node.up]
            edge_precision = 1.0 / (sigma2 * node.dist)
            precision[i, i] += edge_precision
            precision[j, j] += edge_precision
            precision[i, j] -= edge_precision
            precision[j, i] -= edge_precision
        if node.is_leaf and values.get(node.name) is not None:
            error = 0.0 if errors is None else errors[node.name]
            if error == 0.0:
                locked[i] = values[node.name]
            else:
                precision[i, i] += 1.0 / error**2
                linear[i] += values[node.name] / error**2
    free = [i for i in range(len(nodes)) if i not in locked]
    fixed = list(locked)
    fixed_values = np.asarray(list(locked.values()))
    mean, variance = np.zeros(len(nodes)), np.zeros(len(nodes))
    covariance = np.linalg.inv(precision[np.ix_(free, free)])
    mean[free] = covariance @ (
        linear[free] - precision[np.ix_(free, fixed)] @ fixed_values
    )
    variance[free] = covariance.diagonal()
    mean[fixed] = fixed_values
    return {node: (mean[i], variance[i]) for node, i in index.items()}


def dense_tip_covariance(tree, values):
    depths = {tree: 0.0}
    for node in tree.traverse(strategy="preorder"):
        if not node.is_root:
            depths[node] = depths[node.up] + node.dist
    leaves = [leaf for leaf in tree.leaves() if values.get(leaf.name) is not None]
    covariance = np.array(
        [[depths[tree.common_ancestor([a, b])] for b in leaves] for a in leaves]
    )
    return leaves, covariance


def dense_reml(tree, values, sigma2, errors=None):
    leaves, covariance = dense_tip_covariance(tree, values)
    y = np.array([values[leaf.name] for leaf in leaves])
    noise = np.array(
        [0.0 if errors is None else errors[leaf.name] ** 2 for leaf in leaves]
    )
    covariance = sigma2 * covariance + np.diag(noise)
    inverse = np.linalg.inv(covariance)
    root_precision = inverse.sum()
    mean = (inverse @ y).sum() / root_precision
    residual = y - mean
    return -0.5 * (
        (len(leaves) - 1) * math.log(2 * math.pi)
        + np.linalg.slogdet(covariance)[1]
        + math.log(root_precision)
        + residual @ inverse @ residual
    )


def _group_root_covariance(data):
    """Unit-rate covariance with the contracted root state fixed at zero."""
    count = len(data.parents)
    depths = np.zeros(count)
    ancestors = []
    for group in range(count):
        lineage = {group}
        current = group
        while data.parents[current] >= 0:
            current = data.parents[current]
            lineage.add(current)
        ancestors.append(lineage)
        if group:
            depths[group] = depths[data.parents[group]] + data.lengths[group]
    covariance = np.empty((count, count))
    for first in range(count):
        for second in range(count):
            common = ancestors[first] & ancestors[second]
            covariance[first, second] = max(depths[group] for group in common)
    return covariance


def test_production_noise_shift_is_below_every_generalized_noise_eigenvalue():
    """Check the curvature-bound premise using independent dense matrices."""
    rng = np.random.default_rng(872613)
    sources = [
        "(A:1,B:1,C:1,D:1,E:1,F:1)R;",
        "((A:1,B:1)I:1,(C:1,(D:1,E:1)J:1)K:1,F:1)R;",
        "(A:1,(B:1,C:1)I:1,D:1,(E:1,F:1)J:1)R;",
    ]
    for case in range(240):
        tree = tree_from(sources[case % len(sources)])
        for node in tree.traverse():
            if not node.is_root:
                node.dist = (
                    0.0 if rng.random() < 0.08 else float(10 ** rng.uniform(-2.0, 2.0))
                )
        values = {name: float(rng.normal()) for name in "ABCDEF"}
        errors = {name: float(10 ** rng.uniform(-2.0, 2.0)) for name in values}
        if case % 5 == 0:
            errors["A"] = 0.0
        data = continuous._prepare_data(tree, values, errors, None)
        observed = [
            group
            for group, message in enumerate(data.observations)
            if message is not None
        ]
        if len(observed) < 2:
            continue
        upper = continuous._rate_upper_bound(data)
        if upper == 0.0:
            continue
        shift = continuous._noise_rate_shift(data, upper)
        anchor = min(observed, key=lambda group: data.observations[group].variance)
        ordered = [anchor, *(group for group in observed if group != anchor)]
        contrasts = np.column_stack(
            (-np.ones(len(ordered) - 1), np.eye(len(ordered) - 1))
        )
        root_covariance = _group_root_covariance(data)[np.ix_(ordered, ordered)]
        process_covariance = contrasts @ root_covariance @ contrasts.T
        observation_variances = np.diag(
            [data.observations[group].variance for group in ordered]
        )
        noise_covariance = contrasts @ observation_variances @ contrasts.T
        minimum_eigenvalue = float(
            eigh(noise_covariance, process_covariance, eigvals_only=True)[0]
        )
        assert 0.0 < shift <= minimum_eigenvalue * (1.0 + 1e-12)


def test_two_tip_root_mean_variance_and_residual_likelihood():
    tree = tree_from("(A:1,B:3)R;")
    posterior, fit = compute_bm_marginals(tree, {"A": 0.0, "B": 4.0}, sigma2=2.0)
    assert posterior[tree].mean == pytest.approx(1.0)
    assert posterior[tree].variance == pytest.approx(1.5)
    assert fit.restricted_log_likelihood == pytest.approx(
        -0.5 * (math.log(2 * math.pi * 8) + 16 / 8)
    )
    assert posterior[tree["A"]].variance == posterior[tree["B"]].variance == 0.0


@pytest.mark.parametrize("errors", [None, {"A": 0.2, "B": 0.0, "C": 0.7, "D": 1.2}])
@pytest.mark.parametrize(
    "source",
    [
        "((A:1,B:2)I:0.5,C:3,D:0.25)R;",
        "(((A:1)J:2,B:3)I:1,(C:2,D:4)K:0.1)R;",
        "(A:0.3,B:0.7,C:1.1,D:0.5)R;",
    ],
)
def test_fixed_rate_all_node_marginals_match_full_precision_inverse(source, errors):
    tree = tree_from(source)
    values = {"A": -1.0, "B": 2.0, "C": 4.0, "D": None}
    posterior, fit = compute_bm_marginals(
        tree, values, sigma2=0.4, standard_errors=errors
    )
    expected = dense_marginals(tree, values, 0.4, errors)
    for node, marginal in posterior.items():
        np.testing.assert_allclose(
            (marginal.mean, marginal.variance), expected[node], atol=1e-11
        )
    assert fit.restricted_log_likelihood == pytest.approx(
        dense_reml(tree, values, 0.4, errors)
    )


def test_outside_observations_affect_an_internal_node_and_missing_tip():
    tree = tree_from("((A:1,B:1)I:1,C:2,D:1)R;")
    values = {"A": 0.0, "B": 0.0, "C": 8.0}
    posterior, _ = compute_bm_marginals(tree, values, sigma2=1.0)
    assert posterior[tree["I"]].mean == pytest.approx(8 / 7)
    assert posterior[tree["D"]].mean == pytest.approx(posterior[tree].mean)
    assert posterior[tree["D"]].variance == pytest.approx(posterior[tree].variance + 1)


def test_matches_phytools_fastanc_means_and_variances_fixture():
    # Independently generated with phytools 2.3.0 / ape 5.8:
    # fastAnc(tree, c(A=-2, B=.5, C=1.1, D=7.3), vars=TRUE).
    # R is not needed to run this regression test.
    tree = tree_from("((A:0.2,B:1.3)I:0.6,(C:1.1,D:0.8)J:0.4)R;")
    expected = {
        "R": (1.33696397941681, 6.67383557818386),
        "I": (-0.99343910806175, 2.53566461097103),
        "J": (2.89056603773585, 5.43336655986709),
    }
    posterior, _ = compute_bm_marginals(tree, {"A": -2.0, "B": 0.5, "C": 1.1, "D": 7.3})
    for name, moments in expected.items():
        marginal = posterior[tree[name]]
        np.testing.assert_allclose(
            (marginal.mean, marginal.variance), moments, atol=1e-13
        )


def test_exact_bm_reml_matches_independent_tip_covariance():
    tree = tree_from("((A:0.2,B:1.3):0.6,C:1.1,D:0.8);")
    values = {"A": -2.0, "B": 0.5, "C": 1.1, "D": 7.3}
    leaves, covariance = dense_tip_covariance(tree, values)
    inverse = np.linalg.inv(covariance)
    y = np.array([values[leaf.name] for leaf in leaves])
    residual = y - (inverse @ y).sum() / inverse.sum()
    expected_rate = (residual @ inverse @ residual) / (len(y) - 1)
    _, fit = compute_bm_marginals(tree, values)
    assert fit.sigma2 == pytest.approx(expected_rate)
    assert fit.restricted_log_likelihood == pytest.approx(
        dense_reml(tree, values, expected_rate)
    )


@pytest.mark.parametrize("error, expected_rate", [(1.0, 3.0), (2.0, 0.0), (3.0, 0.0)])
def test_known_se_star_reml_has_analytic_interior_and_zero_boundaries(
    error, expected_rate
):
    tree = tree_from("(A:1,B:1,C:1)R;")
    values = {"A": 0.0, "B": 2.0, "C": 4.0}
    errors = dict.fromkeys(values, error)
    posterior, fit = compute_bm_marginals(tree, values, standard_errors=errors)
    assert fit.sigma2 == pytest.approx(expected_rate, abs=1e-7)
    assert posterior[tree].variance == pytest.approx((expected_rate + error**2) / 3)
    assert fit.restricted_log_likelihood == pytest.approx(
        dense_reml(tree, values, expected_rate, errors)
    )
    if expected_rate == 0:
        assert fit.fit_status == "zero_boundary"
        assert all(marginal.mean == pytest.approx(2) for marginal in posterior.values())
        assert all(marginal.variance > 0 for marginal in posterior.values())


def test_heteroscedastic_reml_matches_dense_optimization():
    tree = tree_from("((A:0.2,B:2.3)I:0.6,C:1.1,D:0.8)R;")
    values = {"A": -3.0, "B": 0.5, "C": 1.1, "D": 7.3}
    errors = {"A": 0.02, "B": 0.3, "C": 2.0, "D": 5.0}
    expected = minimize_scalar(
        lambda log_rate: -dense_reml(tree, values, math.exp(log_rate), errors),
        bounds=(-15, 10),
        method="bounded",
        options={"xatol": 1e-10},
    )
    posterior, fit = compute_bm_marginals(tree, values, standard_errors=errors)
    assert fit.sigma2 == pytest.approx(math.exp(expected.x), rel=1e-6)
    for node, expected_values in dense_marginals(
        tree, values, fit.sigma2, errors
    ).items():
        np.testing.assert_allclose(
            (posterior[node].mean, posterior[node].variance), expected_values, rtol=1e-9
        )


@pytest.mark.parametrize("source", ["(O:0.1,A:1,B:1,C:1)R;", "(C:1,B:1,A:1,O:0.1)R;"])
def test_heteroscedastic_reml_finds_the_best_of_multiple_likelihood_modes(source):
    tree = tree_from(source)
    values = {"O": 0.0, "A": 2.2, "B": 7.7, "C": 74.2}
    errors = {"O": 0.0, "A": 0.55, "B": 2.6, "C": 23.3}
    # Independent contrasts to O have covariance r*C+S. Diagonalization gives
    # an independent score with three stationary points, not a unimodal fit.
    covariance = np.eye(3) + 0.1 * np.ones((3, 3))
    noise = np.diag([errors[name] ** 2 for name in "ABC"])
    eigenvalues, eigenvectors = eigh(noise, covariance)
    squares = (eigenvectors.T @ [values[name] for name in "ABC"]) ** 2

    def derivative(rate):
        variances = rate + eigenvalues
        return np.sum((variances - squares) / variances**2)

    def log_likelihood(rate):
        variances = rate + eigenvalues
        return -0.5 * (
            3 * math.log(2 * math.pi)
            + np.linalg.slogdet(covariance)[1]
            + np.sum(np.log(variances) + squares / variances)
        )

    candidates = [0.0] + [
        brentq(derivative, left, right, xtol=1e-12)
        for left, right in [(20.0, 50.0), (150.0, 250.0), (500.0, 800.0)]
    ]
    expected_rate = max(candidates, key=log_likelihood)
    assert expected_rate == pytest.approx(30.120948020760228)
    posterior, fit = compute_bm_marginals(tree, values, standard_errors=errors)
    assert fit.sigma2 == pytest.approx(expected_rate, rel=2e-6)
    assert fit.restricted_log_likelihood == pytest.approx(
        log_likelihood(expected_rate), abs=1e-10
    )
    expected = dense_marginals(tree, values, expected_rate, errors)
    for node, marginal in posterior.items():
        assert (marginal.mean, marginal.variance) == pytest.approx(
            expected[node], rel=3e-6
        )


@pytest.mark.parametrize("error", [1e-6, 1e-8, 1e-100])
def test_rate_independent_measurement_residuals_do_not_change_the_optimum(error):
    tree = tree_from("((A:0,B:0)I:1,C:1)R;")
    posterior, fit = compute_bm_marginals(
        tree,
        {"A": 0.0, "B": 2.0, "C": 4.0},
        standard_errors={"A": error, "B": error, "C": 1.0},
    )
    # A/B reduce to mean=1, variance=v. Their within-position residual is
    # independent of r and may dwarf every between-position likelihood change.
    variance = error**2 / 2
    expected_rate = (3.0**2 - 1.0 - variance) / 2
    assert fit.sigma2 == pytest.approx(expected_rate, rel=1e-6)
    assert posterior[tree].mean == pytest.approx(
        (5 * expected_rate + 1 + 4 * variance) / 9
    )
    assert posterior[tree].variance == pytest.approx(
        (expected_rate + variance) * (expected_rate + 1) / 9
    )
    # The discarded optimization constant still belongs in the reported density.
    assert fit.residual_df == 2
    assert fit.num_effective_observations == 3
    assert fit.restricted_log_likelihood == pytest.approx(-1 / error**2, rel=1e-9)


@pytest.mark.parametrize("size", [1e160, 1e200, 1e300])
def test_unresolvable_rate_is_rejected_instead_of_returning_a_numerical_boundary(size):
    # A/B determine r=0.5, but that rate cannot be resolved after scaling to the
    # extremely uncertain outlier. Returning the first evaluable rate is wrong.
    with pytest.raises(ValueError, match="range|underflow"):
        compute_bm_marginals(
            tree_from("(A:1,B:1,C:1)R;"),
            {"A": 1.0, "B": 2.0, "C": size},
            standard_errors={"A": 0.0, "B": 0.0, "C": size},
        )


def test_underflowing_probe_edge_does_not_abort_a_resolvable_rate_fit():
    tree = tree_from("(A:1,B:1e-320,C:1)R;")
    values = {"A": 0.0, "B": 1.0, "C": 2.0}
    errors = dict.fromkeys(values, 0.01)
    _, fit = compute_bm_marginals(tree, values, standard_errors=errors)
    oracle = minimize_scalar(
        lambda log_rate: -dense_reml(tree, values, math.exp(log_rate), errors=errors),
        bounds=(-20.0, 20.0),
        method="bounded",
        options={"xatol": 1e-12},
    )
    assert fit.sigma2 == pytest.approx(math.exp(oracle.x), rel=1e-6)
    assert fit.restricted_log_likelihood == pytest.approx(-oracle.fun, abs=1e-10)


def test_time_scaling_does_not_erase_a_representable_subnormal_branch():
    tree = tree_from("(A:4e-323,B:64,C:1)R;")
    posterior, fit = compute_bm_marginals(
        tree, {"A": 0.0, "B": 1.0, "C": 2.0}, sigma2=2.0
    )
    assert fit.sigma2 == 2.0
    assert fit.fit_status == "ok"
    assert posterior[tree["A"]].variance == 0.0
    assert math.isfinite(fit.restricted_log_likelihood)


def test_unrepresentable_variance_in_the_final_fixed_model_is_still_rejected():
    with pytest.raises(ValueError, match="variance underflowed"):
        compute_bm_marginals(
            tree_from("(A:1,B:1e-310,C:1)R;"),
            {"A": 0.0, "B": 1.0, "C": 2.0},
            sigma2=1e-20,
        )


@pytest.mark.parametrize("size", [1e20, 1e100, 1e150])
def test_uncertain_outlier_does_not_erase_differences_between_exact_observations(size):
    tree = tree_from("(A:1,B:1,C:1)R;")
    values = {"A": 1.0, "B": 2.0, "C": size}
    errors = {"A": 0.0, "B": 0.0, "C": size}
    posterior, fit = compute_bm_marginals(tree, values, standard_errors=errors)
    assert fit.sigma2 == pytest.approx(0.5, rel=1e-5)
    assert posterior[tree].mean == pytest.approx(1.5)
    assert posterior[tree].variance == pytest.approx(0.25, rel=1e-5)
    with pytest.raises(ValueError, match="Conflicting exact observations"):
        compute_bm_marginals(tree, values, sigma2=0.0, standard_errors=errors)


def test_zero_rate_retains_an_exact_anchor_despite_an_uncertain_outlier():
    tree = tree_from("(A:1,B:1,C:1)R;")
    posterior, fit = compute_bm_marginals(
        tree,
        {"A": 2.0, "B": 2.0, "C": 1e100},
        standard_errors={"A": 0.0, "B": 0.0, "C": 1e100},
    )
    assert fit.sigma2 == 0.0
    assert all((item.mean, item.variance) == (2.0, 0.0) for item in posterior.values())


@pytest.mark.parametrize("fixed_rate", [None, 0.3])
def test_polytomy_is_equivalent_to_exact_zero_edge_resolution(fixed_rate):
    values = {"A": -1.0, "B": 3.0, "C": 5.0}
    star = tree_from("(A:0.3,B:0.7,C:1.1,D:0.5)R;")
    resolved = tree_from("(((A:0.3,B:0.7)I:0,C:1.1)J:0,D:0.5)R;")
    first, first_fit = compute_bm_marginals(star, values, sigma2=fixed_rate)
    second, second_fit = compute_bm_marginals(resolved, values, sigma2=fixed_rate)
    assert first_fit == second_fit
    for node, marginal in first.items():
        assert marginal == second[resolved[node.name]]
    assert second[resolved["I"]] == second[resolved["J"]] == second[resolved]


@pytest.mark.parametrize("fixed_rate", [None, 0.8])
def test_flat_root_brownian_inference_is_invariant_to_root_placement(fixed_rate):
    sources = [
        "[&R]((A:.3,B:.7)I:.4,(C:1.1,D:.8)J:.6)R;",
        "[&R](A:.3,B:.7,(C:1.1,D:.8)J:1.0)I;",
    ]
    values = {"A": -3.0, "B": 0.5, "C": 1.1, "D": 7.3}
    errors = {"A": 0.02, "B": 0.3, "C": 2.0, "D": 5.0}
    trees = [tree_from(source) for source in sources]
    results = [
        compute_bm_marginals(tree, values, sigma2=fixed_rate, standard_errors=errors)
        for tree in trees
    ]
    (first, first_fit), (second, second_fit) = results
    assert second_fit.sigma2 == pytest.approx(first_fit.sigma2, rel=1e-6)
    assert second_fit.restricted_log_likelihood == pytest.approx(
        first_fit.restricted_log_likelihood, abs=1e-12
    )
    for name in "ABCDJ":
        observed, expected = second[trees[1][name]], first[trees[0][name]]
        assert observed.mean == pytest.approx(expected.mean, rel=1e-6, abs=1e-12)
        assert observed.variance == pytest.approx(
            expected.variance, rel=1e-6, abs=1e-12
        )
    observed, expected = second[trees[1]], first[trees[0]["I"]]
    assert observed.mean == pytest.approx(expected.mean, rel=1e-6, abs=1e-12)
    assert observed.variance == pytest.approx(expected.variance, rel=1e-6, abs=1e-12)


def test_zero_edge_duplicate_exact_observations_count_once():
    tree = tree_from("(A:0,B:0,C:1)R;")
    posterior, fit = compute_bm_marginals(tree, {"A": 3.0, "B": 3.0, "C": 5.0})
    assert fit.sigma2 == pytest.approx(4)
    assert fit.num_observed == 3
    assert fit.num_effective_observations == 2
    assert fit.residual_df == 1
    assert posterior[tree] == posterior[tree["A"]] == posterior[tree["B"]]


def test_zero_edge_duplicate_likelihood_matches_independent_reduced_space():
    contracted = tree_from("((A:0,B:0)I:1,C:1,D:2)R;")
    reduced = tree_from("(I:1,C:1,D:2)R;")
    contracted_values = {"A": 1.0, "B": 1.0, "C": 3.0, "D": -2.0}
    reduced_values = {"I": 1.0, "C": 3.0, "D": -2.0}
    _, fixed = compute_bm_marginals(contracted, contracted_values, sigma2=0.7)
    assert fixed.restricted_log_likelihood == pytest.approx(
        dense_reml(reduced, reduced_values, 0.7), abs=1e-12
    )

    _, fitted = compute_bm_marginals(contracted, contracted_values)
    oracle = minimize_scalar(
        lambda log_rate: -dense_reml(reduced, reduced_values, math.exp(log_rate)),
        bounds=(-20.0, 20.0),
        method="bounded",
    )
    assert fitted.sigma2 == pytest.approx(math.exp(oracle.x), rel=1e-5)
    assert fitted.restricted_log_likelihood == pytest.approx(-oracle.fun, abs=1e-10)


def test_zero_edge_noisy_observations_are_independent_measurements():
    tree = tree_from("(A:0,B:0,C:1)R;")
    posterior, fit = compute_bm_marginals(
        tree, {"A": 1.0, "B": 3.0}, sigma2=0.5, standard_errors={"A": 1.0, "B": 1.0}
    )
    assert posterior[tree].mean == 2
    assert posterior[tree].variance == pytest.approx(0.5)
    assert posterior[tree["C"]].variance == pytest.approx(1.0)
    assert fit.num_effective_observations == 2


@pytest.mark.parametrize(
    "source, rate", [("(A:0,B:0,C:1);", 1.0), ("(A:1,B:1,C:1);", 0.0)]
)
def test_conflicting_exact_observations_fail(source, rate):
    with pytest.raises(ValueError, match="Conflicting exact observations"):
        compute_bm_marginals(tree_from(source), {"A": 1.0, "B": 2.0}, sigma2=rate)


def test_identical_exact_data_report_singular_boundary_and_zero_uncertainty():
    tree = tree_from("(A:1,B:2,C:3)R;")
    posterior, fit = compute_bm_marginals(tree, {"A": 12.0, "B": 12.0, "C": 12.0})
    assert fit.sigma2 == 0
    assert fit.restricted_log_likelihood is None
    assert fit.fit_status == "singular_zero_boundary"
    assert all((item.mean, item.variance) == (12.0, 0.0) for item in posterior.values())


@pytest.mark.parametrize(
    "source, values",
    [
        ("(A:1,B:1);", {"A": 1.0}),
        ("(A:0,B:0);", {"A": 1.0, "B": 1.0}),
    ],
)
def test_unidentifiable_rate_requires_fixed_sigma2(source, values):
    tree = tree_from(source)
    with pytest.raises(ValueError, match="not identifiable.*--sigma2"):
        compute_bm_marginals(tree, values)
    posterior, _ = compute_bm_marginals(tree, values, sigma2=2.0)
    assert all(item.mean == 1.0 for item in posterior.values())


def test_no_observations_cannot_anchor_a_flat_root_prior():
    with pytest.raises(ValueError, match="at least one observed"):
        compute_bm_marginals(tree_from("(A:1,B:1);"), {"A": None}, sigma2=1.0)


@pytest.mark.parametrize("error", [-1.0, math.inf, math.nan])
def test_invalid_standard_errors_are_rejected(error):
    with pytest.raises(ValueError, match="Standard error"):
        compute_bm_marginals(
            tree_from("(A:1,B:1);"), {"A": 1.0}, sigma2=1, standard_errors={"A": error}
        )


@pytest.mark.parametrize("rate", [-1.0, math.inf, math.nan, True])
def test_invalid_fixed_rates_are_rejected(rate):
    with pytest.raises(ValueError, match="--sigma2"):
        compute_bm_marginals(tree_from("(A:1,B:1);"), {"A": 1.0}, sigma2=rate)


@pytest.mark.parametrize("time_scale", [1e-100, 1e-12, 1e12, 1e100])
@pytest.mark.parametrize("fixed", [False, True])
def test_time_scaling_preserves_marginals_and_likelihood(time_scale, fixed):
    tree = tree_from("((A:0.5,B:2)I:1,C:3,D:2)R;")
    values = {"A": -2.0, "B": 0.5, "C": 3.0}
    errors = {"A": 0.1, "B": 0.2, "C": 0.4}
    first, first_fit = compute_bm_marginals(
        tree, values, sigma2=0.7 if fixed else None, standard_errors=errors
    )
    for node in tree.traverse():
        if not node.is_root:
            node.dist *= time_scale
    second, second_fit = compute_bm_marginals(
        tree, values, sigma2=0.7 / time_scale if fixed else None, standard_errors=errors
    )
    assert second_fit.sigma2 * time_scale == pytest.approx(first_fit.sigma2, rel=2e-6)
    assert second_fit.restricted_log_likelihood == pytest.approx(
        first_fit.restricted_log_likelihood, abs=1e-10
    )
    for node in first:
        np.testing.assert_allclose(
            (second[node].mean, second[node].variance),
            (first[node].mean, first[node].variance),
            rtol=3e-6,
        )


@pytest.mark.parametrize("scale, offset", [(1e-50, 0), (1e50, 0), (2.0, 1e10)])
def test_trait_units_and_offsets_transform_marginals_and_rate(scale, offset):
    tree = tree_from("((A:1,B:2)I:1,C:3)R;")
    values = {"A": -2.0, "B": 1.0, "C": 4.0}
    first, fit1 = compute_bm_marginals(tree, values)
    second, fit2 = compute_bm_marginals(
        tree, {name: offset + scale * value for name, value in values.items()}
    )
    assert fit2.sigma2 == pytest.approx(fit1.sigma2 * scale**2)
    assert fit2.restricted_log_likelihood == pytest.approx(
        fit1.restricted_log_likelihood - fit1.residual_df * math.log(scale)
    )
    for node in first:
        assert second[node].mean == pytest.approx(offset + scale * first[node].mean)
        assert second[node].variance == pytest.approx(scale**2 * first[node].variance)


def test_large_star_matches_sample_variance_and_does_not_resolve_tree():
    count = 4096
    tree = tree_from("(" + ",".join(f"T{i}:1" for i in range(count)) + ")R;")
    values = {f"T{i}": float(i % 7) for i in range(count)}
    expected_variance = np.var(list(values.values()), ddof=1)
    posterior, fit = compute_bm_marginals(tree, values)
    assert fit.sigma2 == pytest.approx(expected_variance)
    assert posterior[tree].mean == pytest.approx(np.mean(list(values.values())))
    assert posterior[tree].variance == pytest.approx(expected_variance / count)
    assert len(tree.children) == count


def test_deep_tree_does_not_use_python_recursion():
    tree = make_deep_ladder_tree(1500)
    values = {leaf.name: float(index % 3) for index, leaf in enumerate(tree.leaves())}
    posterior, fit = compute_bm_marginals(tree, values, sigma2=1.0)
    assert len(posterior) == 2999
    assert fit.fit_status == "ok"


def test_optimizer_failure_is_not_silently_accepted(monkeypatch):
    class Failed:
        success = False
        fun = 1.0

    monkeypatch.setattr(continuous, "minimize_scalar", lambda *args, **kwargs: Failed())
    with pytest.raises(ValueError, match="failed to converge"):
        compute_bm_marginals(
            tree_from("(A:1,B:1,C:1);"),
            {"A": 0.0, "B": 2.0, "C": 4.0},
            standard_errors={"A": 1.0, "B": 1.0, "C": 1.0},
        )
