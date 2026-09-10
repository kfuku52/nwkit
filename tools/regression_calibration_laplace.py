#!/usr/bin/env python3
"""Independent importance-QMC check of Laplace likelihoods at fitted parameters.

This measures integration error, not bias in the optimum or CI calibration.
Agreement across randomized Sobol scrambles is required before interpreting a
reference value. A posterior Gaussian proposal keeps this practical for small
correlated phylogenetic latent vectors without changing the target integral.
"""

import argparse
import hashlib
import sys
from pathlib import Path

import numpy as np
from scipy.special import expit, gammaln, logsumexp
from scipy.stats import norm, qmc

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from regression_calibration_design import Case, seed_for  # noqa: E402
from regression_calibration_engine import call_glmm  # noqa: E402
from validate_regression_calibration import (  # noqa: E402
    encode,
    read_records,
    source_files,
)


def conditional_log_likelihood(y, eta, family, alpha=None):
    if family == "binomial":
        values = y * eta - np.logaddexp(0, eta)
    elif family == "poisson":
        values = y * eta - np.exp(eta) - gammaln(y + 1)
    elif family == "negative-binomial":
        size = 1 / alpha
        values = (
            gammaln(y + size)
            - gammaln(size)
            - gammaln(y + 1)
            + size * np.log(size)
            + y * eta
            - (y + size) * np.logaddexp(np.log(size), eta)
        )
    else:
        raise ValueError(family)
    return np.sum(values, axis=-1)


def reference(case, data, *, powers=(10, 12, 14), scrambles=8, seed=1):
    fit = call_glmm(case, data)
    y = np.asarray(data["y"])
    n = len(y)
    design = np.column_stack([np.ones(n), data["x"]])
    fixed = design @ fit.coefficients.reshape(-1)
    covariance = np.asarray(data["covariance"])
    covariance = (
        covariance
        / np.mean(np.diag(covariance))
        * fit.component_variances["phylogenetic"]
    )
    precision = np.linalg.inv(covariance)
    mode = np.asarray(fit.random_modes).reshape(-1)
    eta = fixed + mode
    if case.family == "binomial":
        weights = expit(eta) * expit(-eta)
    elif case.family == "poisson":
        weights = np.exp(eta)
    else:
        weights = (
            np.exp(eta)
            * (1 + fit.dispersion * y)
            / (1 + fit.dispersion * np.exp(eta)) ** 2
        )
    posterior_precision = precision + np.diag(weights)
    proposal_covariance = np.linalg.inv(posterior_precision)
    proposal_factor = np.linalg.cholesky(proposal_covariance)
    prior_logdet = np.linalg.slogdet(covariance)[1]
    proposal_logdet = np.linalg.slogdet(proposal_covariance)[1]
    conditional_mode = conditional_log_likelihood(y, eta, case.family, fit.dispersion)
    reconstructed_laplace = (
        conditional_mode
        - 0.5 * mode @ precision @ mode
        - 0.5 * prior_logdet
        - 0.5 * np.linalg.slogdet(posterior_precision)[1]
    )
    convergence = []
    for power in powers:
        log_integrals, effective_samples = [], []
        for scramble in range(scrambles):
            uniforms = qmc.Sobol(
                n,
                scramble=True,
                seed=np.random.default_rng(seed_for(seed, case.name, scramble, "qmc")),
            ).random_base2(power)
            z = norm.ppf(
                np.clip(uniforms, np.nextafter(0.0, 1.0), np.nextafter(1.0, 0.0))
            )
            latent = mode + z @ proposal_factor.T
            log_weights = (
                conditional_log_likelihood(
                    y, fixed + latent, case.family, fit.dispersion
                )
                - 0.5 * np.einsum("ij,jk,ik->i", latent, precision, latent)
                + 0.5 * np.sum(z * z, axis=1)
                + 0.5 * (proposal_logdet - prior_logdet)
            )
            log_integrals.append(
                float(logsumexp(log_weights) - np.log(len(log_weights)))
            )
            normalized = np.exp(log_weights - logsumexp(log_weights))
            effective_samples.append(float(1 / np.sum(normalized**2)))
        values = np.exp(np.array(log_integrals) - max(log_integrals))
        estimate = float(max(log_integrals) + np.log(np.mean(values)))
        log_mcse = float(np.std(values, ddof=1) / np.sqrt(scrambles) / np.mean(values))
        convergence.append(
            {
                "power": power,
                "draws_per_scramble": 2**power,
                "scrambles": scrambles,
                "log_integral": estimate,
                "log_mcse": log_mcse,
                "min_importance_ess": min(effective_samples),
            }
        )
    stable = (
        len(convergence) >= 2
        and convergence[-1]["log_mcse"] < 0.005
        and abs(convergence[-1]["log_integral"] - convergence[-2]["log_integral"])
        < 0.01
    )
    return {
        "fit_log_likelihood": fit.log_likelihood,
        "reconstructed_laplace": reconstructed_laplace,
        "reconstruction_error": reconstructed_laplace - fit.log_likelihood,
        "reference_stable_at_0.01_log_units": stable,
        "laplace_minus_reference": fit.log_likelihood - convergence[-1]["log_integral"],
        "convergence": convergence,
        "boundary_warning": fit.boundary_warning,
        "separation_warning": fit.separation_warning,
        "coefficient": fit.coefficients.reshape(-1),
        "component_variances": fit.component_variances,
        "dispersion_alpha": fit.dispersion,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--cases",
        default="binomial-n8-p0.5,binomial-n30-boundary,poisson-n30,negative-binomial-n30",
    )
    parser.add_argument("--per-case", type=int, default=3)
    args = parser.parse_args()
    selected = set(args.cases.split(","))
    counts = {name: 0 for name in selected}
    results = []
    for row in read_records(args.input):
        name = row["case"]["name"]
        if name not in selected or counts[name] >= args.per_case:
            continue
        case = Case(**row["case"])
        data = {
            key: np.asarray(value) if isinstance(value, list) else value
            for key, value in row["data"].items()
        }
        if len(np.unique(data["y"])) < 2:
            continue
        result = reference(case, data, seed=row["fit_seed"])
        results.append(
            {
                "case": name,
                "replicate": row["replicate"],
                "input_sha256": row["input_sha256"],
                **result,
            }
        )
        counts[name] += 1
        print(
            f"{name} replicate {row['replicate']}: Laplace-reference={result['laplace_minus_reference']:.5f}; stable={result['reference_stable_at_0.01_log_units']}",
            flush=True,
        )
        if all(value >= args.per_case for value in counts.values()):
            break
    source_root, paths = source_files()
    report = {
        "source_sha256": {
            str(path.relative_to(source_root)): hashlib.sha256(
                path.read_bytes()
            ).hexdigest()
            for path in paths
        },
        "scope": "likelihood at fitted parameters only, not coefficient/SE reference optimization",
        "input_protocol_sha256": hashlib.sha256(
            (args.input / "protocol.json").read_bytes()
        ).hexdigest(),
        "reference_source_sha256": hashlib.sha256(
            Path(__file__).read_bytes()
        ).hexdigest(),
        "results": results,
    }
    with args.output.open("x") as handle:
        handle.write(encode(report) + "\n")


if __name__ == "__main__":
    main()
