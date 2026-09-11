"""Phylogenetic signal from Gaussian likelihoods and label permutations.

Independent implementation of Blomberg's K and Pagel's lambda. Covariances
are relative to the retained tips' MRCA; the root mean is profiled by GLS.
"""

import numpy as np
from scipy.linalg import cho_solve, solve_triangular
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import chi2


def gls(covariance, values):
    factor = np.linalg.cholesky(covariance)
    weights = cho_solve((factor, True), np.ones(len(values)))
    mean = float(weights @ values / weights.sum())
    residual = values - mean
    quadratic = float(residual @ cho_solve((factor, True), residual))
    logdet = float(2 * np.log(np.diag(factor)).sum())
    return mean, quadratic, logdet, float(weights.sum())


def k_statistic(covariance, values):
    mean, quadratic, _, precision = gls(covariance, values)
    expected = (np.trace(covariance) - len(values) / precision) / (len(values) - 1)
    if quadratic <= 0 or expected <= 0:
        raise ValueError("K is undefined for constant traits or degenerate covariance.")
    return float(np.sum((values - mean) ** 2) / quadratic / expected)


def _best_scalar(function, grid):
    """Search each grid interval and explicitly retain boundary optima."""
    candidates = [(float(x), function(float(x))) for x in grid]
    for left, right in zip(grid[:-1], grid[1:], strict=True):
        result = minimize_scalar(
            lambda x: -function(x)[0],
            bounds=(left, right),
            method="bounded",
            options={"xatol": 1e-9},
        )
        if not result.success:
            raise ValueError("Signal likelihood optimization did not converge.")
        candidates.append((float(result.x), function(float(result.x))))
    return max(candidates, key=lambda item: item[1][0])


def rate_fit(covariance, values, errors):
    """Profile the root mean and ML diffusion rate, retaining rate=0."""
    n = len(values)
    if not np.any(errors):
        mean, quadratic, logdet, _ = gls(covariance, values)
        rate = quadratic / n
        if rate <= 0:
            raise ValueError("Zero residual variance has no finite Gaussian density.")
        likelihood = -0.5 * (n * (1 + np.log(2 * np.pi * rate)) + logdet)
        return float(likelihood), float(rate), mean
    evaluate, lower, upper = _error_rate_profile(covariance, values, errors)
    if upper <= 0:
        return evaluate(0.0)
    # Search in log(rate): exact observations can identify rates many orders of
    # magnitude smaller than the marginal variation of noisy observations.
    logs = np.linspace(
        np.log(lower),
        np.log(upper),
        max(2, int(np.ceil(np.log10(upper) - np.log10(lower))) + 1),
    )
    _, fit = _best_scalar(lambda z: evaluate(np.exp(z)), logs)
    return max([fit, evaluate(0.0)], key=lambda item: item[0])


def _error_rate_profile(covariance, values, errors):
    exact = np.flatnonzero(errors == 0)
    noisy = np.flatnonzero(errors > 0)
    if len(exact) and np.ptp(values[exact]) == 0:
        raise ValueError("Unbounded likelihood at singular zero diffusion rate.")
    order = np.r_[exact, noisy]
    factor = np.linalg.cholesky(covariance[np.ix_(order, order)])
    k = len(exact)
    centered = values - values[0]
    full_y = solve_triangular(factor, centered[order], lower=True)
    upper = float(full_y @ full_y)
    exact_y = full_y[:k]
    exact_one = solve_triangular(factor[:k, :k], np.ones(k), lower=True)
    noise_y = centered[noisy] - factor[k:, :k] @ exact_y
    noise_one = np.ones(len(noisy)) - factor[k:, :k] @ exact_one
    # Whiten by positive SEs, then use an SVD of the conditional BM factor.
    # Unlike diagonalizing C^-1/2 E C^-1/2, this never invents variance for
    # exact observations or loses the rank of E through eigensolver roundoff.
    left, singular, _ = np.linalg.svd(
        factor[k:, k:] / errors[noisy, None], full_matrices=False
    )
    noise_y = left.T @ (noise_y / errors[noisy])
    noise_one = left.T @ (noise_one / errors[noisy])
    log_eigen = 2 * np.log(singular)
    logdet = 2 * np.log(np.diag(factor)[:k]).sum() + 2 * np.log(errors[noisy]).sum()
    if k:
        mean = exact_one @ exact_y / (exact_one @ exact_one)
        lower = float(np.sum((exact_y - mean * exact_one) ** 2) / len(values))
    else:
        lower = (
            min(upper, float(np.min(errors) ** 2 / np.max(np.diag(covariance)))) * 1e-10
        )
    if (
        lower <= 0
        or not np.isfinite(lower + upper)
        or not np.all(np.isfinite(log_eigen))
    ):
        raise ValueError(
            "Sampling-error dynamic range exceeds floating-point resolution; rescale units."
        )

    def evaluate(rate):
        if rate == 0 and k:
            return (-np.inf, 0.0, np.nan)
        log_rate = np.log(rate) if rate > 0 else -np.inf
        log_diagonal = np.logaddexp(0, log_rate + log_eigen)
        weights = np.exp(-log_diagonal / 2)
        y = np.r_[exact_y / np.sqrt(rate) if k else [], noise_y * weights]
        ones = np.r_[exact_one / np.sqrt(rate) if k else [], noise_one * weights]
        mean = float(ones @ y / (ones @ ones))
        q = float(np.sum((y - mean * ones) ** 2))
        determinant = logdet + (k * log_rate if k else 0) + log_diagonal.sum()
        likelihood = -0.5 * (len(values) * np.log(2 * np.pi) + determinant + q)
        return float(likelihood), float(rate), float(mean + values[0])

    return evaluate, lower, upper


def k_fit(covariance, values, errors):
    if not np.any(errors):
        return k_statistic(covariance, values), None
    _, rate, _ = rate_fit(covariance, values, errors)
    # Finite-sample correction used by the Ives/phytools error-aware K.
    corrected = rate * len(values) / (len(values) - 1)
    effective = corrected * covariance + np.diag(errors**2)
    return k_statistic(effective, values), corrected


def permutation_test(covariance, values, errors, observed, simulations, rng):
    exceed = 0
    for _ in range(simulations):
        order = rng.permutation(len(values))
        statistic, _ = k_fit(covariance, values[order], errors[order])
        exceed += statistic >= observed - 1e-12 * max(1.0, abs(observed))
    return (1 + exceed) / (1 + simulations)


def lambda_fit(covariance, values, errors, ci_level=0.95):
    diagonal = np.diag(np.diag(covariance))
    shared = covariance - diagonal
    if not np.any(shared):
        return {"status": "unidentifiable_lambda"}
    cache = {}

    def profile(value):
        value = float(value)
        if value not in cache:
            cache[value] = rate_fit(diagonal + value * shared, values, errors)
        return cache[value]

    estimate, fit = _best_scalar(profile, np.linspace(0, 1, 11))
    likelihood, rate, mean = fit
    null = profile(0)[0]
    if (
        rate == 0
        or max(v[0] for v in cache.values()) - min(v[0] for v in cache.values()) < 1e-9
    ):
        return {
            "status": "unidentifiable_lambda",
            "sigma2": rate,
            "log_likelihood": likelihood,
        }
    lr = max(0.0, 2 * (likelihood - null))
    lower, upper = profile_interval(profile, estimate, likelihood, ci_level)
    return {
        "status": "boundary" if estimate in (0.0, 1.0) else "ok",
        "estimate": estimate,
        "sigma2": rate,
        "root_mean": mean,
        "log_likelihood": likelihood,
        "null_log_likelihood": null,
        "likelihood_ratio": lr,
        "p_value": float(chi2.sf(lr, 1)),
        "ci_lower": lower,
        "ci_upper": upper,
    }


def profile_interval(profile, estimate, likelihood, level):
    """Connected profile-likelihood interval around the selected maximum."""
    cutoff = likelihood - chi2.ppf(level, 1) / 2

    def difference(value):
        return profile(value)[0] - cutoff

    limits = []
    for endpoint in (0.0, 1.0):
        previous = estimate
        limit = endpoint
        for value in np.linspace(estimate, endpoint, 65)[1:]:
            if difference(value) < 0:
                limit = brentq(difference, min(previous, value), max(previous, value))
                break
            previous = value
        limits.append(float(limit))
    return limits


def bh_adjust(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    adjusted = values[order] * len(values) / np.arange(1, len(values) + 1)
    result = np.empty(len(values))
    result[order] = np.minimum(1, np.minimum.accumulate(adjusted[::-1])[::-1])
    return result
