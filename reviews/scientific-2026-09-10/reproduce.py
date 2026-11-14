"""Small independent diagnostics for the scientific review; no production edits.

Run from the repository root with:
    PYTHONPATH=. python reviews/scientific-2026-09-10/reproduce.py

These are counterexamples and population calculations, not an error-rate study
of the complete application. Existing archived shift/RADTE simulations are not
rerun by this script.
"""

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.stats import chi2

from nwkit.gaussian import grouped_average_marginal_logdet
from nwkit.phylogenetic_glmm import (
    _censored_gaussian_log_likelihood,
    _coefficient_likelihood_inference,
    _draw_censored_gaussian,
)
from nwkit.regress import _profile_covariance_fit
from nwkit.threshold_asr import _rhat_and_ess


def event_objective():
    # Independent events, two paralogs per event. The stated generative model
    # has residual variance s=1 and a shared event effect with variance t=1.
    # The production event balancing doubles the row-specific covariance.
    m = 20
    groups = np.repeat(np.arange(m), 2)
    z = np.eye(m)[groups]
    identity = np.eye(2 * m)
    shared = z @ z.T
    truth = identity + shared

    def expected_objective(st):
        s, t = st
        working = 2 * s * identity + t * shared
        return 0.5 * (
            grouped_average_marginal_logdet(working, groups)
            + np.trace(np.linalg.solve(working, truth))
        )

    step = 1e-5
    gradient = [
        (
            expected_objective(np.ones(2) + step * direction)
            - expected_objective(np.ones(2) - step * direction)
        )
        / (2 * step)
        for direction in np.eye(2)
    ]
    fit = minimize(lambda log_st: expected_objective(np.exp(log_st)), [0.0, 0.0])
    return {
        "interpretation": "Expected ML objective, known zero mean; no simulation error",
        "true_residual_and_event_variance": [1.0, 1.0],
        "expected_gradient_at_truth": gradient,
        "population_optimum": np.exp(fit.x).tolist(),
        "optimizer_success": bool(fit.success),
    }


def event_bootstrap_scale():
    m = 10
    n = 2 * m
    groups = np.repeat(np.arange(m), 2)
    x = np.ones((n, 1))
    y = np.random.default_rng(19).normal(size=n)
    fit = _profile_covariance_fit(
        y,
        x,
        np.zeros(n),
        [("evolutionary_rate", np.full(n, 2.0))],
        reml=False,
        likelihood_observations=m,
        likelihood_groups=groups,
    )
    rate = fit["component_variances"]["evolutionary_rate"]
    # Closed form, not Monte Carlo: bootstrap draws use the fitted working
    # covariance 2*rate*I. Refitting gives sum(centered_y**2)/(2*m).
    return {
        "fitted_rate": rate,
        "working_row_variance": float(np.asarray(fit["covariance"])[0]),
        "exact_expected_bootstrap_refit_rate": (n - 1) / m * rate,
        "expected_bootstrap_rate_ratio": (n - 1) / m,
        "scope": "replicate-reml covariance without random effects, ML fit, intercept only",
    }


def threshold_diagnostics():
    rng = np.random.default_rng(72)
    traces = rng.normal(size=(4, 1000, 1))
    traces += np.linspace(-2, 2, 1000)[None, :, None]
    rhat, ess = _rhat_and_ess(traces)
    # A plain split R-hat suffices for this counterexample; it is not a full
    # implementation of rank-normalized/folded R-hat or multi-lag ESS.
    split = np.concatenate([traces[:, :500, 0], traces[:, 500:, 0]], axis=0)
    within = np.var(split, axis=1, ddof=1).mean()
    between = 500 * np.var(split.mean(axis=1), ddof=1)
    split_rhat = np.sqrt(((499 / 500) * within + between / 500) / within)
    constant_rhat, constant_ess = _rhat_and_ess(np.zeros((4, 1000, 1)))
    return {
        "scope": "diagnostic helper counterexample, not an actual fitted threshold chain",
        "drifting_chain_rhat": float(rhat[0]),
        "drifting_chain_ess": float(ess[0]),
        "passes_current_status_cutoffs": bool(rhat[0] <= 1.05 and ess[0] >= 400),
        "plain_split_rhat": float(split_rhat),
        "constant_chain_rhat": float(constant_rhat[0]),
        "constant_chain_ess": float(constant_ess[0]),
    }


def censored_bootstrap():
    # Example observation mechanism: left censoring at 0, 50 censored and
    # 50 exact observations. A conditional bootstrap must respect Y>0 for the
    # exact rows; an unconditional bootstrap must regenerate censoring labels.
    lower = np.full(100, np.nan)
    upper = np.r_[np.zeros(50), np.full(50, np.nan)]
    generated = _draw_censored_gaussian(
        np.random.default_rng(93), np.zeros(100), 1.0, lower, upper
    )
    # E[uncensored synthetic Y]=0. The likelihood score at mu=0 therefore
    # contains 50*(-phi(0)/Phi(0)) from permanently censored rows.
    expected_values = np.r_[np.full(50, np.nan), np.zeros(50)]
    step = 1e-5

    def ll(mu):
        return float(
            np.sum(
                _censored_gaussian_log_likelihood(
                    expected_values, np.full(100, mu), 1.0, lower, upper
                )
            )
        )

    return {
        "scope": "known SD, zero latent variance limit, common detection threshold 0",
        "censored_rows_remaining_censored": int(np.isnan(generated[:50]).sum()),
        "exact_rows_generated_below_detection_threshold": int(
            (generated[50:] < 0).sum()
        ),
        "expected_refit_loglikelihood_score_at_generating_mean": (ll(step) - ll(-step))
        / (2 * step),
    }


def penalized_likelihood_ratio():
    # Exact one-dimensional Gaussian likelihood plus Gaussian penalty. The
    # generic inference helper still assigns an ordinary chi-square tail.
    y = 2.0
    penalty_precision = 3.0
    optimum = np.array([y / (1 + penalty_precision), 0.0])

    def objective(beta):
        return (
            0.5 * (beta[0] - y) ** 2
            + 0.5 * penalty_precision * beta[0] ** 2
            + 0.5 * beta[1] ** 2
        )

    stats, pvalues, _, _ = _coefficient_likelihood_inference(
        objective,
        optimum,
        [(-20.0, 20.0), (-20.0, 20.0)],
        1,
        np.array([[1 / (1 + penalty_precision)]]),
        inference="likelihood-ratio",
        confidence_level=0.95,
    )
    return {
        "scope": "generic helper counterexample; not an end-to-end GLMM calibration",
        "reported_statistic": float(stats[0]),
        "reported_p_value": float(pvalues[0]),
        "exact_null_tail_for_this_statistic": float(chi2.sf(y**2, 1)),
        "null_distribution": "chi-square(1) / 4, not chi-square(1)",
    }


def main():
    root = Path(__file__).resolve().parents[2]
    sources = [
        "nwkit/regress.py",
        "nwkit/gaussian.py",
        "nwkit/threshold_asr.py",
        "nwkit/phylogenetic_glmm.py",
        "nwkit/shift_backend.py",
        "nwkit/shift_cli.py",
    ]
    result = {
        "event_objective": event_objective(),
        "event_bootstrap_scale": event_bootstrap_scale(),
        "threshold_diagnostics": threshold_diagnostics(),
        "censored_bootstrap": censored_bootstrap(),
        "penalized_likelihood_ratio": penalized_likelihood_ratio(),
        "source_sha256": {
            source: hashlib.sha256((root / source).read_bytes()).hexdigest()
            for source in sources
        },
    }
    output = Path(__file__).with_name("results.json")
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
