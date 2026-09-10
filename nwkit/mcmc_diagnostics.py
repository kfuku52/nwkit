"""Rank/split diagnostics following Vehtari et al. (2021).

Arrays have shape (independent chains, retained draws). ESS uses FFT
autocovariances and Geyer's initial positive/monotone sequence, including the
antithetic correction used by Stan posterior. No sampler dependencies.
"""

import math

import numpy as np
from scipy.fft import next_fast_len
from scipy.special import ndtri
from scipy.stats import rankdata

DIAGNOSTIC_VERSION = "rank_split_v1"


def _split(x):
    half = x.shape[1] // 2
    return np.concatenate((x[:, :half], x[:, -half:]), axis=0)


def _rank(x):
    ranks = rankdata(x, method="average").reshape(x.shape)
    return ndtri((ranks - 0.375) / (x.size + 0.25))


def _constant(x):
    return bool(np.any(np.ptp(x, axis=1) == 0))


def _rhat(x):
    if np.ptp(x) == 0:
        return math.nan
    n = x.shape[1]
    within = np.var(x, axis=1, ddof=1).mean()
    if within == 0:
        return math.inf
    return float(np.sqrt((n - 1) / n + np.var(x.mean(axis=1), ddof=1) / within))


def _ess(x):
    """Multi-chain ESS on already transformed/split draws."""
    m, n = x.shape
    if n < 3 or np.ptp(x) == 0:
        return math.nan
    centered = x - x.mean(axis=1, keepdims=True)
    length = next_fast_len(2 * n)
    spectrum = np.fft.rfft(centered, n=length, axis=1)
    acov = np.fft.irfft(spectrum * spectrum.conjugate(), n=length, axis=1)
    acov = acov[:, :n].mean(axis=0) / n
    within = acov[0] * n / (n - 1)
    pooled = acov[0] + np.var(x.mean(axis=1), ddof=1)
    rho = np.zeros(n)
    rho[0] = even = 1.0
    rho[1] = odd = 1 - (within - acov[1]) / pooled
    t = 0
    while t < n - 5 and even + odd > 0:
        t += 2
        even = 1 - (within - acov[t]) / pooled
        odd = 1 - (within - acov[t + 1]) / pooled
        if even + odd >= 0:
            rho[t : t + 2] = (even, odd)
    if even > 0:
        rho[t] = even
    for index in range(2, t - 1, 2):
        if rho[index] + rho[index + 1] > rho[index - 2] + rho[index - 1]:
            rho[index : index + 2] = (rho[index - 2] + rho[index - 1]) / 2
    tau = -1 + 2 * rho[:t].sum() + rho[t]
    return float(m * n / max(tau, 1 / math.log10(m * n)))


def diagnose(values, *, indicator=False, structural=False):
    """Return diagnostics and explicit reasons; unavailable values never pass.

    Tail ESS is inapplicable to Bernoulli indicators; their mean ESS and MCSE
    directly assess the reported probability. Fixed quantities are excluded.
    """
    x = np.asarray(values, dtype=float)
    if x.ndim != 2 or not all(x.shape):
        raise ValueError("MCMC diagnostics require a nonempty chains-by-draws array.")
    result = dict.fromkeys(
        ("rhat", "ess_bulk", "ess_tail", "ess_mean", "mcse_mean"), math.nan
    )
    if structural:
        return {**result, "status": "structural_constant"}
    if not np.all(np.isfinite(x)):
        return {**result, "status": "nonfinite_trace"}
    if x.shape[0] < 2:
        return {**result, "status": "mcmc_rhat_unavailable"}
    if x.shape[1] < 8:
        return {**result, "status": "insufficient_draws"}
    if _constant(x) or _constant(_split(x)):
        reason = "unresolved_rare_category" if indicator else "constant_trace"
        return {**result, "status": reason}
    split = _split(x)
    ranked = _rank(split)
    rank_rhat = _rhat(ranked)
    if indicator:
        # Folding a balanced binary variable can be identically constant.
        result["rhat"] = rank_rhat
    else:
        folded_rhat = _rhat(_rank(_split(np.abs(x - np.median(x)))))
        result["rhat"] = float(np.maximum(rank_rhat, folded_rhat))
        tails = [
            _ess(_split((x <= np.quantile(x, q)).astype(float))) for q in (0.05, 0.95)
        ]
        result["ess_tail"] = float(np.min(tails))
    result["ess_bulk"] = _ess(ranked)
    result["ess_mean"] = _ess(split)
    result["mcse_mean"] = float(np.sqrt(np.var(x, ddof=1) / result["ess_mean"]))
    required = [result[k] for k in ("rhat", "ess_bulk", "ess_mean", "mcse_mean")]
    if not indicator:
        required.append(result["ess_tail"])
    reasons = []
    if not all(math.isfinite(v) for v in required):
        reasons.append("diagnostic_unavailable")
    if result["rhat"] > 1.01:
        reasons.append("mcmc_rhat")
    if result["ess_bulk"] < 400 or (not indicator and result["ess_tail"] < 400):
        reasons.append("mcmc_low_ess")
    if indicator and (result["ess_mean"] < 400 or result["mcse_mean"] > 0.01):
        reasons.append("mcmc_probability_precision")
    return {**result, "status": "+".join(reasons) if reasons else "ok"}
