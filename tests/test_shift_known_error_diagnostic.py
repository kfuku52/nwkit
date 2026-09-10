"""Validate paired generator comparisons against the public fit and dense GLS."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from diagnose_shift_known_error import diagnose, summarize  # noqa: E402

from nwkit.shift_calibration import CalibratedSearch  # noqa: E402
from nwkit.util import read_tree  # noqa: E402


def test_plugin_matches_fit_and_oracle_matches_independent_dense_gls():
    tree = read_tree("((A:1,B:1):1,(C:1,D:1):1);", "auto", True, quiet=True)
    search = CalibratedSearch(
        tree,
        convergence=False,
        variances=[0.01, 0.03, 0.1, 0.25],
        alpha_grid=[0, np.inf],
    )
    y = np.array([0.2, -0.4, 1.1, 0.9])
    covariance = search.geometry(0.3)[0] + np.diag(search.variances)
    row = diagnose(search, y, np.zeros(4), covariance, search.families[0], 41, 19)
    fit = search.fit(y, seed=41, replicates=19)
    assert row["plugin_p"] == fit["tests"][0]["p_value"]
    noise = np.random.default_rng(41).normal(size=(3, 19))
    simulations = np.linalg.cholesky(search.q @ covariance @ search.q.T) @ noise
    z = np.column_stack((search.q @ (y - y.mean()), simulations))
    scores = np.full((len(search.models), 20), -np.inf)
    for ai, variance, _, _, _, _ in search.cache:
        K, weights = search.geometry(search.grid[ai])
        C = search.q @ (variance * K + np.diag(search.variances)) @ search.q.T
        inverse = np.linalg.inv(C)
        _, logdet = np.linalg.slogdet(C)
        for m in range(len(search.models)):
            X = search.q @ ((search.loads[m] * weights[m]) @ search.transform[m])
            X = X[:, : search.dim[m]]
            beta = np.linalg.solve(X.T @ inverse @ X, X.T @ inverse @ z)
            residual = z - X @ beta
            score = -0.5 * (
                3 * np.log(2 * np.pi)
                + logdet
                + np.sum(residual * (inverse @ residual), axis=0)
            )
            scores[m] = np.maximum(scores[m], score)
    statistics = 2 * (scores.max(axis=0) - scores[search.families[0]].max(axis=0))
    assert row["statistic"] == pytest.approx(statistics[0], abs=1e-9)
    assert (
        row["oracle_p"]
        == (1 + np.count_nonzero(statistics[1:] >= statistics[0] - 1e-10)) / 20
    )


def test_summary_retains_failed_attempts_and_paired_disagreement():
    record = dict(plugin_reject=True, oracle_reject=False)
    rows = [
        dict(
            case={"alpha_height": 0}, status="completed", null=record, one_shift=record
        ),
        dict(case={"alpha_height": 0}, status="failed"),
    ]
    summary = summarize(rows)
    assert summary["attempted_pairs"] == 2
    assert summary["completed_pairs"] == 1
    cell = summary["cells"][0]
    assert cell["attempted"] == 2
    assert cell["completed"] == 1
    assert cell["plugin_only"] == 1
    assert cell["oracle_only"] == 0
    assert cell["plugin"]["denominator"] == 1


@pytest.fixture(scope="module")
def audited_pair():
    from diagnose_shift_known_error import execute

    return execute(
        dict(
            case_id=0,
            tips=4,
            alpha_height=0.01,
            known_error=True,
            seed=101,
            replicates=19,
        )
    )


def test_auditor_refits_and_rejects_changed_observed_statistic(audited_pair):
    import copy

    from verify_shift_known_error import check_row

    assert audited_pair["status"] == "completed"
    assert check_row(audited_pair) == 0
    altered = copy.deepcopy(audited_pair)
    altered["one_shift"]["statistic"] += 1
    with pytest.raises(AssertionError):
        check_row(altered)


def test_auditor_regenerates_true_covariance(audited_pair):
    import copy

    from verify_shift_known_error import check_row

    altered = copy.deepcopy(audited_pair)
    altered["generating_covariance"][0][0] += 0.1
    with pytest.raises(AssertionError):
        check_row(altered)


@pytest.mark.parametrize(
    "probabilities,budget,reject,complete,lower",
    [
        ([0.05], 1, None, False, 0.05),
        ([0.05, 0.2], 2, False, False, 0.2),
        ([0.05] * 50, None, True, True, 0.05),
    ],
)
def test_nuisance_envelope_never_rejects_an_incomplete_grid(
    monkeypatch, probabilities, budget, reject, complete, lower
):
    from shift_known_error_envelope import known_error_envelope

    tree = read_tree("((A:1,B:1):1,(C:1,D:1):1);", "auto", True, quiet=True)
    search = CalibratedSearch(
        tree, convergence=False, variances=np.ones(4) * 0.1, alpha_grid=[0]
    )
    values = iter(probabilities)
    monkeypatch.setattr(search, "_probability_from_noise", lambda *args: next(values))
    result = known_error_envelope(
        search, [0.2, -0.4, 1.1, 0.9], replicates=19, max_evaluations=budget
    )
    assert result["reject"] is reject
    assert result["grid_complete"] is complete
    assert result["p_value_lower_bound"] == lower
    assert result["p_value_upper_bound"] == (lower if complete else 1)
    assert result["grid_point_count"] == 50


def test_envelope_first_point_reproduces_plugin_probability(audited_pair):
    from shift_known_error_envelope import known_error_envelope

    search = CalibratedSearch(
        read_tree(audited_pair["tree"], "auto", True, quiet=True),
        convergence=False,
        variances=audited_pair["variances"],
    )
    result = known_error_envelope(
        search,
        audited_pair["values"],
        seed=audited_pair["null"]["bootstrap_seed"],
        replicates=19,
        max_evaluations=1,
    )
    assert result["p_value_lower_bound"] == audited_pair["null"]["plugin_p"]
    assert result["grid_complete"] is False
    assert result["reject"] is (
        False if audited_pair["null"]["plugin_p"] > 0.05 else None
    )


def test_envelope_reports_bootstrap_boundary_failure_as_unresolved(monkeypatch):
    from shift_known_error_envelope import known_error_envelope

    tree = read_tree("((A:1,B:1):1,(C:1,D:1):1);", "auto", True, quiet=True)
    search = CalibratedSearch(
        tree, convergence=False, variances=np.ones(4) * 0.1, alpha_grid=[0]
    )
    original = search._probability_from_noise
    calls = []

    def reaches_boundary(*args):
        calls.append(1)
        if len(calls) == 1:
            return 0.05
        # A genuine numerical guard from the production profile, not a fake p.
        return original(
            np.zeros(3),
            search.cache[-1][2],
            1,
            search.families[0],
            0,
            np.random.default_rng(1).normal(size=(3, 19)),
        )

    monkeypatch.setattr(search, "_probability_from_noise", reaches_boundary)
    result = known_error_envelope(search, [0.2, -0.4, 1.1, 0.9], replicates=19)
    assert result["reject"] is None
    assert result["p_value_lower_bound"] == 0.05
    assert result["p_value_upper_bound"] == 1
    assert "upper grid boundary" in result["failed_evaluation"]["error"]
    assert len(result["evaluations"]) == 1


def test_variance_tail_certification_matches_a_wider_independent_grid():
    from shift_variance_tail import TailCertifiedSearch

    tree = read_tree("((A:1,B:1):1,(C:1,D:1):1);", "auto", True, quiet=True)
    search = TailCertifiedSearch(
        tree, convergence=False, variances=np.ones(4) * 0.1, alpha_grid=[0, np.inf]
    )
    Z = search.cache[-1][2] @ np.random.default_rng(1).normal(size=(3, 19))
    with pytest.raises(ValueError, match="upper grid boundary"):
        CalibratedSearch.profile(search, Z)
    scores, selected = search.profile(Z)
    for m in range(len(search.models)):
        for column in range(Z.shape[1]):
            item = search.cache[selected[m, column]]
            K, weights = search.geometry(search.grid[item[0]])
            C = search.q @ (item[1] * K + np.diag(search.variances)) @ search.q.T
            inverse = np.linalg.inv(C)
            X = search.q @ ((search.loads[m] * weights[m]) @ search.transform[m])
            X = X[:, : search.dim[m]]
            z = Z[:, column]
            beta = np.linalg.solve(X.T @ inverse @ X, X.T @ inverse @ z)
            residual = z - X @ beta
            expected = -0.5 * (
                3 * np.log(2 * np.pi)
                + np.linalg.slogdet(C)[1]
                + residual @ inverse @ residual
            )
            assert scores[m, column] == pytest.approx(expected, abs=1e-9)
    assert search.extension_count > 0
    assert search.tail_log_likelihood_upper_bound < scores.min()
    # Evaluate additional, otherwise unnecessary variances with a separate
    # dense inverse GLS formula, checking the claimed upper-tail exclusion.
    ceiling = search.variance_grid[-1]
    for variance in (ceiling, ceiling * 2, ceiling * 100):
        for alpha in search.grid:
            K, weights = search.geometry(alpha)
            C = search.q @ (variance * K + np.diag(search.variances)) @ search.q.T
            inverse = np.linalg.inv(C)
            _, logdet = np.linalg.slogdet(C)
            for m in range(len(search.models)):
                X = search.q @ ((search.loads[m] * weights[m]) @ search.transform[m])
                X = X[:, : search.dim[m]]
                beta = np.linalg.solve(X.T @ inverse @ X, X.T @ inverse @ Z)
                residual = Z - X @ beta
                ll = -0.5 * (
                    3 * np.log(2 * np.pi)
                    + logdet
                    + np.sum(residual * (inverse @ residual), axis=0)
                )
                assert np.all(ll < scores[m])
    # Previously certified fits must remain unchanged after extending cache
    # for another data batch; this is essential for bootstrap batching.
    search.profile(Z * 20)
    replay, _ = search.profile(Z)
    np.testing.assert_allclose(replay, scores, rtol=0, atol=1e-10)


def test_tail_backend_completes_a_full_nuisance_grid():
    from shift_known_error_envelope import known_error_envelope
    from shift_variance_tail import TailCertifiedSearch

    tree = read_tree("((A:1,B:1):1,(C:1,D:1):1);", "auto", True, quiet=True)
    search = TailCertifiedSearch(
        tree, convergence=False, variances=np.ones(4) * 0.1, alpha_grid=[0]
    )
    result = known_error_envelope(search, [1e4, -1e4, 1e4, -1e4], replicates=19)
    assert result["grid_complete"] is True
    assert result["reject"] is True
    assert result["p_value_lower_bound"] == result["p_value_upper_bound"] == 0.05
    assert len(result["evaluations"]) == result["grid_point_count"]
    assert "failed_evaluation" not in result
    assert search.extension_count > 0
