"""Phylogenetic signal from Gaussian likelihoods and label permutations.

Independent implementation of Blomberg's K and Pagel's lambda. Covariances
are relative to the retained tips' MRCA; the root mean is profiled by GLS.
"""

from functools import cache

import numpy as np
from scipy.linalg import cho_solve, eigh, solve_triangular
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import chi2


def _prepare_gls(covariance, count):
    factor = np.linalg.cholesky(covariance)
    weights = cho_solve((factor, True), np.ones(count))
    precision = float(weights.sum())
    logdet = float(2 * np.log(np.diag(factor)).sum())
    return factor, weights, precision, logdet


def _gls_from_preparation(prepared, values):
    factor, weights, precision, logdet = prepared
    mean = float(weights @ values / precision)
    residual = values - mean
    quadratic = float(residual @ cho_solve((factor, True), residual))
    return mean, quadratic, logdet, precision


def gls(covariance, values):
    return _gls_from_preparation(_prepare_gls(covariance, len(values)), values)


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


def _zero_error_rate(gls_result, count):
    mean, quadratic, logdet, _ = gls_result
    rate = quadratic / count
    if rate <= 0:
        raise ValueError("Zero residual variance has no finite Gaussian density.")
    likelihood = -0.5 * (count * (1 + np.log(2 * np.pi * rate)) + logdet)
    return float(likelihood), float(rate), mean


def rate_fit(covariance, values, errors):
    """Profile the root mean and ML diffusion rate, retaining rate=0."""
    if not np.any(errors):
        return _zero_error_rate(gls(covariance, values), len(values))
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
    whitened_y = np.empty(len(values), dtype=float)
    whitened_one = np.empty(len(values), dtype=float)

    def evaluate(rate):
        if rate == 0 and k:
            return (-np.inf, 0.0, np.nan)
        log_rate = np.log(rate) if rate > 0 else -np.inf
        log_diagonal = np.logaddexp(0, log_rate + log_eigen)
        weights = np.exp(-log_diagonal / 2)
        if k:
            whitened_y[:k] = exact_y / np.sqrt(rate)
            whitened_one[:k] = exact_one / np.sqrt(rate)
        whitened_y[k:] = noise_y * weights
        whitened_one[k:] = noise_one * weights
        mean = float(whitened_one @ whitened_y / (whitened_one @ whitened_one))
        q = float(np.sum((whitened_y - mean * whitened_one) ** 2))
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


def _prepared_k_statistic(covariance, count):
    """Reuse one covariance factor while only the tip labels change."""
    factor = np.linalg.cholesky(covariance)
    weights = cho_solve((factor, True), np.ones(count))
    precision = float(weights.sum())
    expected = (np.trace(covariance) - count / precision) / (count - 1)

    def statistic(values):
        mean = float(weights @ values / precision)
        residual = values - mean
        quadratic = float(residual @ cho_solve((factor, True), residual))
        if quadratic <= 0 or expected <= 0:
            raise ValueError(
                "K is undefined for constant traits or degenerate covariance."
            )
        return float(np.sum((values - mean) ** 2) / quadratic / expected)

    return statistic


def permutation_test(covariance, values, errors, observed, simulations, rng):
    # Sampling errors move with the labels and require a rate refit per draw.
    # Prepare only after the first permutation, preserving the random stream
    # even when an invalid covariance/constant trait fails on that draw.
    statistic = None
    exceed = 0
    for _ in range(simulations):
        order = rng.permutation(len(values))
        if statistic is None and not np.any(errors):
            statistic = _prepared_k_statistic(covariance, len(values))
        if statistic is None:
            value, _ = k_fit(covariance, values[order], errors[order])
        else:
            value = statistic(values[order])
        exceed += value >= observed - 1e-12 * max(1.0, abs(observed))
    return (1 + exceed) / (1 + simulations)


def _zero_error_lambda_spectrum(covariance):
    sd = np.sqrt(np.diag(covariance))
    correlation = covariance / sd[:, None] / sd[None, :]
    eigenvalues, eigenvectors = eigh(correlation, driver="evd")
    # A resolved eigenvalue alone does not ensure accurate log likelihoods:
    # near lambda=1, spectral errors grow with conditioning and sample size.
    # Keep Cholesky profiles in that regime without changing scientific rank.
    accuracy_floor = (
        np.sqrt(np.finfo(float).eps) * len(sd) * np.max(np.abs(eigenvalues))
    )
    if eigenvalues[0] <= accuracy_floor:
        return None
    projected_one = eigenvectors.T @ (np.ones(len(sd)) / sd)
    diagonal_logdet = float(2 * np.log(sd).sum())
    return sd, eigenvalues, eigenvectors, projected_one, diagonal_logdet


def _zero_error_lambda_profile(covariance, values, errors, prepare_profile=None):
    """Prepare an exact-observation profile from the covariance spectrum."""
    if np.any(errors) or np.ptp(values) == 0:
        return None
    # Check the endpoints first, including rejection of a singular covariance.
    if prepare_profile is None:
        boundary = {
            0.0: rate_fit(np.diag(np.diag(covariance)), values, errors),
            1.0: rate_fit(covariance, values, errors),
        }
        spectrum = _zero_error_lambda_spectrum(covariance)
    else:
        boundary, spectrum = prepare_profile(values)
    if spectrum is None:
        return None
    sd, eigenvalues, eigenvectors, projected_one, diagonal_logdet = spectrum
    reference = values[0]
    projected_y = eigenvectors.T @ ((values - reference) / sd)

    def profile(value):
        if value in boundary:
            return boundary[value]
        variance = (1 - value) + value * eigenvalues
        weight = projected_one / np.sqrt(variance)
        response = projected_y / np.sqrt(variance)
        mean = float(weight @ response / (weight @ weight))
        residual = response - mean * weight
        rate = float(residual @ residual / len(values))
        if rate <= 0:
            raise ValueError("Zero residual variance has no finite Gaussian density.")
        logdet = diagonal_logdet + float(np.log(variance).sum())
        likelihood = -0.5 * (len(values) * (1 + np.log(2 * np.pi * rate)) + logdet)
        return float(likelihood), rate, float(mean + reference)

    return profile


def lambda_fit(covariance, values, errors, ci_level=0.95):
    """Fit lambda; a None interval level skips profile limits for resampling."""
    return _fit_lambda(covariance, values, errors, ci_level)


def _fit_lambda(covariance, values, errors, ci_level, prepare_profile=None):
    diagonal = np.diag(np.diag(covariance))
    shared = covariance - diagonal
    if not np.any(shared):
        return {"status": "unidentifiable_lambda"}
    cache = {}
    spectral_profile = _zero_error_lambda_profile(
        covariance, values, errors, prepare_profile
    )

    def profile(value):
        value = float(value)
        if value not in cache:
            cache[value] = (
                rate_fit(diagonal + value * shared, values, errors)
                if spectral_profile is None
                else spectral_profile(value)
            )
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
    result = {
        "status": "boundary" if estimate in (0.0, 1.0) else "ok",
        "estimate": estimate,
        "sigma2": rate,
        "root_mean": mean,
        "log_likelihood": likelihood,
        "null_log_likelihood": null,
        "likelihood_ratio": lr,
        "p_value": float(chi2.sf(lr, 1)),
    }
    if ci_level is not None:
        result["ci_lower"], result["ci_upper"] = profile_interval(
            profile, estimate, likelihood, ci_level
        )
    return result


def lambda_parametric_bootstrap(
    covariance, values, errors, observed_lr, simulations, rng
):
    """Simulate at the fitted lambda=0 null and refit both models each time."""
    diagonal = np.diag(covariance)
    _, null_rate, null_mean = rate_fit(np.diag(diagonal), values, errors)
    standard_deviation = np.sqrt(null_rate * diagonal + errors**2)
    if not np.all(np.isfinite(standard_deviation)):
        raise ValueError("Non-finite lambda null simulation variance.")
    # Every refitted likelihood ratio is nonnegative. At observed LR=0 all
    # simulations would count as exceedances, so the Monte Carlo P is 1.
    if observed_lr == 0:
        return 1.0

    @cache
    def prepare_spectrum():
        return _zero_error_lambda_spectrum(covariance)

    @cache
    def prepare_null():
        return _prepare_gls(np.diag(diagonal), len(diagonal))

    @cache
    def prepare_bm():
        return _prepare_gls(covariance, len(diagonal))

    cache_endpoints = False

    def prepare_profile(simulated):
        nonlocal cache_endpoints
        # Finish the first eigensolve before retaining the endpoint factors.
        if not cache_endpoints:
            boundary = {
                0.0: rate_fit(np.diag(diagonal), simulated, errors),
                1.0: rate_fit(covariance, simulated, errors),
            }
            spectrum = prepare_spectrum()
            cache_endpoints = spectrum is not None
            return boundary, spectrum
        count = len(diagonal)
        boundary = {
            0.0: _zero_error_rate(
                _gls_from_preparation(prepare_null(), simulated), count
            ),
            1.0: _zero_error_rate(
                _gls_from_preparation(prepare_bm(), simulated), count
            ),
        }
        return boundary, prepare_spectrum()

    exceed = 0
    tolerance = 1e-12 * max(1.0, abs(observed_lr))
    for _ in range(simulations):
        simulated = null_mean + standard_deviation * rng.standard_normal(len(values))
        fit = _fit_lambda(covariance, simulated, errors, None, prepare_profile)
        # A zero-rate or flat replicate has no detectable lambda effect (LR=0).
        lr = fit.get("likelihood_ratio", 0.0)
        exceed += lr >= observed_lr - tolerance
    return (1 + exceed) / (1 + simulations)


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
