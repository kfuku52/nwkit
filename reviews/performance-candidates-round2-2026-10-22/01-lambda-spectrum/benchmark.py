"""Scoped zero-error lambda spectral prototype and complete-fit comparisons."""

import argparse
import json
import os
import subprocess
import sys
from contextlib import contextmanager
from pathlib import Path

SOURCE = Path("/tmp/nwkit-perf-round2-20261022-rkzoGe/baseline")
HERE = Path(__file__).resolve().parent
for variable in (
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "MKL_NUM_THREADS",
):
    os.environ[variable] = "1"
os.environ["NWKIT_CI_DIAGNOSTICS"] = "0"
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(SOURCE))

import numpy as np  # noqa: E402
from scipy.linalg import eigh  # noqa: E402
from scipy.stats import chi2  # noqa: E402
from support import environment, measurement_lock, paired, peak_bytes  # noqa: E402

from nwkit import signal_stats as stats  # noqa: E402

ORIGINAL = stats.lambda_fit


def candidate_environment():
    return {
        **environment(),
        "nwkit_source": str(Path(stats.__file__).resolve()),
        "nwkit_package": str(Path(sys.modules["nwkit"].__file__).resolve()),
    }


def spectral_lambda_fit(covariance, values, errors, ci_level=0.95):
    """Evaluate only the exact-observation, spectrally resolved regime."""
    if np.any(errors) or np.ptp(values) == 0:
        return ORIGINAL(covariance, values, errors, ci_level)
    diagonal = np.diag(np.diag(covariance))
    shared = covariance - diagonal
    if not np.any(shared):
        return {"status": "unidentifiable_lambda"}

    # Evaluate both boundaries through the original rate calculation. In
    # particular, a singular supplied C must still fail without invented rank.
    boundary = {
        0.0: stats.rate_fit(diagonal, values, errors),
        1.0: stats.rate_fit(covariance, values, errors),
    }
    standard_deviation = np.sqrt(np.diag(covariance))
    correlation = covariance / standard_deviation[:, None] / standard_deviation[None, :]
    eigenvalues, eigenvectors = eigh(correlation, driver="evd")
    resolution = 64 * np.finfo(float).eps * len(values) * np.max(np.abs(eigenvalues))
    if eigenvalues[0] <= resolution:
        # Eigendecomposition cannot resolve the endpoint spectrum reliably.
        # Retain original Cholesky profiles rather than altering a rank/tolerance.
        return ORIGINAL(covariance, values, errors, ci_level)
    reference = values[0]
    projected_one = eigenvectors.T @ (np.ones(len(values)) / standard_deviation)
    projected_y = eigenvectors.T @ ((values - reference) / standard_deviation)
    diagonal_logdet = float(2 * np.log(standard_deviation).sum())
    cache = dict(boundary)

    def profile(value):
        value = float(value)
        if value not in cache:
            variance = (1 - value) + value * eigenvalues
            weight = projected_one / np.sqrt(variance)
            response = projected_y / np.sqrt(variance)
            mean = float(weight @ response / (weight @ weight))
            residual = response - mean * weight
            rate = float(residual @ residual / len(values))
            if rate <= 0:
                raise ValueError(
                    "Zero residual variance has no finite Gaussian density."
                )
            logdet = diagonal_logdet + float(np.log(variance).sum())
            likelihood = -0.5 * (len(values) * (1 + np.log(2 * np.pi * rate)) + logdet)
            cache[value] = float(likelihood), rate, float(mean + reference)
        return cache[value]

    # Identical optimization, statuses, inference and interval logic.
    estimate, fit = stats._best_scalar(profile, np.linspace(0, 1, 11))
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
        result["ci_lower"], result["ci_upper"] = stats.profile_interval(
            profile, estimate, likelihood, ci_level
        )
    return result


@contextmanager
def patch_fit(function):
    previous = stats.lambda_fit
    stats.lambda_fit = function
    try:
        yield
    finally:
        stats.lambda_fit = previous


def fixture(n, *, shape="balanced", seed=71, truth=0.65, nonultrametric=False):
    rng = np.random.default_rng(seed)
    covariance = np.zeros((n, n))
    if shape == "balanced":
        size = 2
        while size < n:
            for start in range(0, n, size):
                covariance[start : start + size, start : start + size] += 0.25
            size *= 2
    elif shape == "pectinate":
        for start in range(1, n - 1):
            covariance[start:, start:] += 1 / n
    covariance[np.diag_indices(n)] += rng.lognormal(0, 0.7, n) if nonultrametric else 1
    diagonal = np.diag(np.diag(covariance))
    generating = diagonal + truth * (covariance - diagonal)
    values = 1.7 + np.linalg.cholesky(generating) @ rng.standard_normal(n)
    return covariance, values, np.zeros(n)


def compare_fits(first, second):
    assert set(first) == set(second), (first, second)
    assert first["status"] == second["status"], (first, second)
    differences = {}
    for name in first:
        if name == "status":
            continue
        left, right = first[name], second[name]
        assert np.isfinite(left) and np.isfinite(right), (name, left, right)
        differences[name] = {
            "baseline": left,
            "candidate": right,
            "absolute_difference": abs(left - right),
            "scaled_difference": abs(left - right) / max(1, abs(left)),
        }
    return {"status": first["status"], "fields": differences, "exact": first == second}


def measure_case(name, args, *, ci_level=0.95):
    def baseline():
        return ORIGINAL(*args, ci_level=ci_level)

    def candidate():
        return spectral_lambda_fit(*args, ci_level=ci_level)

    return {
        "case": name,
        "n_tips": len(args[1]),
        "ci_level": ci_level,
        "scope": "complete lambda_fit, including preparation, search and optional intervals",
        "comparison": compare_fits(baseline(), candidate()),
        "timing": paired(baseline, candidate),
        "peak_traced_bytes": {
            "baseline": peak_bytes(baseline),
            "candidate": peak_bytes(candidate),
        },
    }


def bootstrap_comparison(args, *, simulations=19, seed=443):
    observed = ORIGINAL(*args)["likelihood_ratio"]

    def run(function):
        rng = np.random.default_rng(seed)
        with patch_fit(function):
            p_value = stats.lambda_parametric_bootstrap(
                *args, observed, simulations, rng
            )
        return {"p_value": p_value, "rng_state": rng.bit_generator.state}

    def baseline():
        return run(ORIGINAL)

    def candidate():
        return run(spectral_lambda_fit)

    a, b = baseline(), candidate()
    return {
        "case": f"seeded_{simulations}_draw_lambda_bootstrap_{len(args[1])}tips",
        "observed_lr": observed,
        "baseline": a,
        "candidate": b,
        "exact": a == b,
        "scope": (
            "zero-LR early-exit control; no simulation is performed"
            if observed == 0
            else "complete bootstrap; each fit includes its own spectral preparation"
        ),
        "performed_draws": 0 if observed == 0 else simulations,
        "timing": paired(baseline, candidate),
    }


def profile_audit(args):
    """Capture identical grid points and check the Gaussian formula directly."""
    from scipy.stats import multivariate_normal

    c, values, errors = args
    records = []
    original_search = stats._best_scalar
    for fit_function in (ORIGINAL, spectral_lambda_fit):
        captured = {}

        def grid_only(profile, grid, _captured=captured):
            for value in grid:
                _captured[float(value)] = profile(float(value))
            return max(_captured.items(), key=lambda item: item[1][0])

        stats._best_scalar = grid_only
        try:
            fit_function(*args, ci_level=None)
        finally:
            stats._best_scalar = original_search
        records.append(captured)
    output = []
    diagonal = np.diag(np.diag(c))
    for value in records[0]:
        left, right = records[0][value], records[1][value]
        scaled = diagonal + value * (c - diagonal)
        independent_one = np.linalg.solve(scaled, np.ones(len(values)))
        independent_mean = float(independent_one @ values / independent_one.sum())
        residual = values - independent_mean
        independent_rate = float(
            residual @ np.linalg.solve(scaled, residual) / len(values)
        )
        direct_logpdf = float(
            multivariate_normal.logpdf(
                values,
                mean=np.full(len(values), independent_mean),
                cov=independent_rate * scaled,
            )
        )
        output.append(
            {
                "lambda": value,
                "baseline": list(left),
                "candidate": list(right),
                "candidate_likelihood_minus_direct_logpdf": right[0] - direct_logpdf,
                "independent_direct_mle": [
                    direct_logpdf,
                    independent_rate,
                    independent_mean,
                ],
                "candidate_minus_independent_mle": [
                    right[0] - direct_logpdf,
                    right[1] - independent_rate,
                    right[2] - independent_mean,
                ],
                "absolute_differences": [
                    abs(a - b) for a, b in zip(left, right, strict=True)
                ],
            }
        )
    return output


def additional():
    c, y, e = fixture(64)
    guarded = np.ones((8, 8)) + 1e-14 * np.eye(8)
    guard_args = (guarded, y[:8], e[:8])
    comparisons = []
    for name, args in (
        ("interior_large_offset", (c, y + 1e9, e)),
        ("interior_small_time_large_traits", (c * 2.0**-300, y * 2.0**100, e)),
        ("unresolved_spectrum_guard", guard_args),
    ):
        comparisons.append(
            {
                "case": name,
                "comparison": compare_fits(ORIGINAL(*args), spectral_lambda_fit(*args)),
            }
        )
    direct_profiles = profile_audit(fixture(64, nonultrametric=True))
    bootstrap = bootstrap_comparison(fixture(64))
    assert bootstrap["observed_lr"] > 0
    assert bootstrap["performed_draws"] == 19
    null_rng = np.random.default_rng(912)
    null_rate_fit = stats.rate_fit(np.diag(np.diag(c)), y, e)
    means, rates = null_rate_fit[2], null_rate_fit[1]
    replicate_lrs = []
    for _ in range(19):
        draw = means + np.sqrt(rates * np.diag(c)) * null_rng.standard_normal(len(y))
        fits = [
            function(c, draw, e, ci_level=None)
            for function in (ORIGINAL, spectral_lambda_fit)
        ]
        replicate_lrs.append([fit.get("likelihood_ratio", 0.0) for fit in fits])
    cutoffs = [0.0, 0.1, 0.5, 2.0, 6.0, *[pair[0] for pair in replicate_lrs]]
    tails = []
    for cutoff in cutoffs:
        tolerance = 1e-12 * max(1.0, abs(cutoff))
        counts = [
            sum(pair[column] >= cutoff - tolerance for pair in replicate_lrs)
            for column in (0, 1)
        ]
        tails.append(
            {"cutoff": cutoff, "exceedances": counts, "exact": counts[0] == counts[1]}
        )
    return {
        "environment": candidate_environment(),
        "edge_comparisons": comparisons,
        "direct_profiles": direct_profiles,
        "positive_lr_bootstrap": bootstrap,
        "bootstrap_tail_audit": {"replicate_lrs": replicate_lrs, "cutoffs": tails},
    }


def edge_cases():
    cases = {}
    c, y, e = fixture(8)
    cases["positive_errors_unchanged"] = (c, y, np.full(8, 0.25))
    cases["mixed_errors_unchanged"] = (c, y, np.array([0.0, 0.1] * 4))
    cases["constant_unchanged"] = (c, np.ones(8), e)
    cases["star"] = (np.eye(8), y, e)
    cases["singular"] = (np.ones((8, 8)), y, e)
    cases["nonpositive_diagonal"] = (c - 10 * np.eye(8), y, e)
    near = np.ones((8, 8)) + 1e-12 * np.eye(8)
    cases["near_singular_1e-12_control"] = (near, y, e)
    cases["large_trait_offset"] = (c, y + 1e9, e)
    cases["small_time_units"] = (c * 2.0**-300, y * 2.0**100, e)
    results = []
    for name, args in cases.items():
        outcomes = []
        for function in (ORIGINAL, spectral_lambda_fit):
            try:
                outcomes.append({"fit": function(*args)})
            except (ValueError, np.linalg.LinAlgError) as exc:
                outcomes.append({"exception": type(exc).__name__, "message": str(exc)})
        if "fit" in outcomes[0] and "fit" in outcomes[1]:
            comparison = compare_fits(outcomes[0]["fit"], outcomes[1]["fit"])
        else:
            comparison = {
                "same_exception": outcomes[0] == outcomes[1],
                "outcomes": outcomes,
            }
        results.append({"case": name, "comparison": comparison})
    return results


def run_tests():
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join((str(HERE), str(SOURCE)))
    command = [
        sys.executable,
        "tools/check.py",
        "test",
        "--",
        "-p",
        "candidate01_plugin",
        "tests/test_signal.py",
        "tests/test_signal_calibration.py",
        "-rs",
    ]
    log = HERE / "prototype-tests.log"
    with measurement_lock(), log.open("w") as handle:
        handle.write("Prototype scoped monkeypatch; cwd=" + str(SOURCE) + "\n")
        handle.write(" ".join(command) + "\n")
        handle.flush()
        outcome = subprocess.run(
            command, cwd=SOURCE, env=env, stdout=handle, stderr=subprocess.STDOUT
        )
    return outcome.returncode


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tests", action="store_true")
    parser.add_argument("--additional", action="store_true")
    parser.add_argument("--output", type=Path, default=HERE / "results.json")
    options = parser.parse_args()
    if options.tests:
        return run_tests()
    if options.additional:
        with measurement_lock():
            result = additional()
            (HERE / "additional-results.json").write_text(
                json.dumps(result, indent=2, allow_nan=False) + "\n"
            )
        return 0
    with measurement_lock():
        results = {
            "candidate": "01 zero-error lambda spectral profile",
            "environment": candidate_environment(),
            "workloads": [
                measure_case("balanced_8", fixture(8)),
                measure_case("balanced_64", fixture(64)),
                measure_case("balanced_256", fixture(256)),
                measure_case("balanced_512_no_interval", fixture(512), ci_level=None),
                measure_case(
                    "pectinate_128_nonultrametric",
                    fixture(128, shape="pectinate", nonultrametric=True),
                ),
            ],
            "edge_cases": edge_cases(),
            "bootstrap": bootstrap_comparison(fixture(32, truth=1)),
        }
        options.output.write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
