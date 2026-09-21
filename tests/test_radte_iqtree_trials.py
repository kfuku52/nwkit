"""Reject poor IQ2MC trial scores without consuming their undefined derivatives."""

import gzip
import shutil
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pytest

from nwkit.radte_iqtree import IQTreeLikelihood, read_export
from nwkit.radte_model import DatingProblem, fit_dates, solve_problem
from tests.test_radte import small_chronology


def _export(tmp_path, score="-123.5", gradient="1 2 inf 3 4"):
    prefix = tmp_path / "trial"
    prefix.with_suffix(".mcmctree.hessian").write_text(
        "4\n(A:1,B:2,(C:3,D:4):5);\n1 2 5 3 4\n"
        + gradient
        + "\nHessian\n"
        + "1 0 0 0 1\n" * 5
    )
    with gzip.open(prefix.with_suffix(".ckp.gz"), "wt") as handle:
        handle.write(f"CandidateSet:\n 0: {score} (A,B,C,D);\n")
    return prefix


def test_score_only_export_does_not_accept_undefined_derivatives(tmp_path):
    prefix = _export(tmp_path)
    result = read_export(prefix, derivatives=False)
    assert result[3] == 123.5
    assert not np.isfinite(result[4]).all()
    with pytest.raises(ValueError, match="Nonfinite IQ-TREE gradient"):
        read_export(prefix)


@pytest.mark.parametrize("derivatives", [False, True])
@pytest.mark.parametrize("score", ["nan", "inf", "-inf"])
def test_nonfinite_score_is_never_accepted(tmp_path, derivatives, score):
    prefix = _export(tmp_path, score=score, gradient="1 2 3 4 5")
    with pytest.raises(ValueError, match="Nonfinite IQ-TREE likelihood"):
        read_export(prefix, derivatives=derivatives)


def test_cached_score_cannot_bypass_derivative_validation():
    exact = IQTreeLikelihood.__new__(IQTreeLikelihood)
    lengths = np.ones(5)
    exact.cache = OrderedDict(
        {
            lengths.tobytes(): (
                lengths,
                123.5,
                np.full(5, np.inf),
                np.eye(5),
                np.arange(5),
                {},
            )
        }
    )
    assert exact.value(lengths) == 123.5
    with pytest.raises(ValueError, match="Nonfinite IQ-TREE gradient"):
        exact.value_gradient(lengths)
    exact.cache[lengths.tobytes()] = (
        lengths,
        123.5,
        np.ones(5),
        np.full((5, 5), np.nan),
        np.arange(5),
        {},
    )
    with pytest.raises(ValueError, match="Nonfinite IQ-TREE Hessian"):
        exact.value_gradient(lengths)


class _TrialLikelihood:
    def __init__(self):
        self.rejected_scores = 0
        self.gradient_calls = 0

    def value(self, lengths):
        if np.min(lengths) < 1e-9:
            self.rejected_scores += 1
        return 1000 * np.sum((lengths - 0.05) ** 2)

    def value_gradient(self, lengths):
        self.gradient_calls += 1
        if np.min(lengths) < 1e-9:
            raise ValueError("Derivative undefined at trial point")
        return self.value(lengths), 2000 * (lengths - 0.05)


def test_solver_rejects_trial_without_requesting_its_derivative():
    likelihood = _TrialLikelihood()
    problem = DatingProblem(small_chronology(), likelihood=likelihood, rate_sd=0.3)
    x, attempts = solve_problem(problem, starts=1, maxiter=500)
    assert likelihood.rejected_scores > 0
    assert likelihood.gradient_calls > 0
    assert attempts[0]["success"]
    assert problem.feasible(x)
    # The exact analytical optimum has every branch length .05 and no rate
    # deviations; both the likelihood and rate penalty are zero there.
    assert problem.value_gradient(x)[0] == pytest.approx(0, abs=1e-9)
    assert problem.value(x) == pytest.approx(problem.value_gradient(x)[0], abs=1e-12)


def test_solver_still_rejects_undefined_derivative_at_accepted_point():
    class InvalidLikelihood(_TrialLikelihood):
        def value_gradient(self, lengths):
            raise ValueError("Derivative undefined at accepted point")

    problem = DatingProblem(
        small_chronology(), likelihood=InvalidLikelihood(), rate_sd=0.3
    )
    with pytest.raises(ValueError, match="Derivative undefined at accepted point"):
        solve_problem(problem, starts=1, maxiter=500)


@pytest.mark.integration
@pytest.mark.skipif(shutil.which("iqtree3") is None, reason="IQ-TREE runtime required")
def test_freerate_clock_fit_with_tiny_rejected_iq2mc_trials():
    c = small_chronology(
        gene_text="((A_1:0.12,A_2:0.15)nA:0.17,(B_1:0.13,B_2:0.11)nB:0.18)Root;"
    )
    alignment = Path(__file__).parent / "data/radte-iqtree-trials/alignment.fa"
    exact = IQTreeLikelihood(
        c, alignment, "GY+F3X4+R4", executable="iqtree3", interface="cli"
    )
    try:
        tiny = {
            "A_1": 5.812731384787386e-13,
            "A_2": 1.906806182439666e-14,
            "B_1": 5.604098119492752e-13,
            "B_2": 1.814881853565924e-12,
        }
        lengths = np.array(
            [tiny.get(n.name, 0.004035275727651381 / 2) for n in c.edges]
        )
        assert np.isfinite(exact.value(lengths))
        assert exact.value(lengths) > exact.prefit_nll
        # IQ-TREE must still evaluate the supplied lengths, without clipping.
        np.testing.assert_allclose(
            exact.evaluate(lengths, derivatives=False)[0], lengths, rtol=1e-5, atol=0
        )
        fit, problem = fit_dates(c, likelihood=exact, rate_sd=0.3, maxiter=500)
        assert problem.feasible(fit.parameters)
        assert np.isfinite(problem.value_gradient(fit.parameters)[1]).all()
        assert fit.objective == pytest.approx(problem.value(fit.parameters), abs=1e-8)
    finally:
        exact.close()
