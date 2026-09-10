"""Information criteria at native maximum-likelihood fits.

pBIC uses conditional information in OU-optimum coordinates. Boundary fits
retain the declared parameter count; these scores are not calibrated tests.
"""

import math

import numpy as np

CRITERIA = {"AIC", "AICc", "BIC", "pBIC"}


def native_information_criterion(data, result, criterion):
    if criterion not in CRITERIA:
        raise ValueError(
            "Native information criterion must be AIC, AICc, BIC, or pBIC."
        )
    shifts = len(result["layout"].shifts)
    n = len(data.tree.leaf_names)
    penalty = 2 * shifts if criterion in {"AIC", "AICc"} else shifts * math.log(n)
    if criterion == "pBIC":
        penalty = 2 * shifts * math.log(len(data.tree.branch_ids) - 2)
    parameters = shifts
    for j, (fit, record) in enumerate(
        zip(result["fits"], result["traits"], strict=True)
    ):
        optimizer = record["optimizer"]
        covariance_parameters = sum(
            int(optimizer[key])
            for key in (
                "alpha_estimated",
                "process_variance_estimated",
                "measurement_variance_estimated",
            )
        )
        dimension = fit.num_mean_parameters
        parameters += dimension + covariance_parameters
        if criterion in {"AIC", "AICc"}:
            penalty += 2 * (dimension + covariance_parameters)
        elif criterion == "BIC":
            penalty += (dimension + covariance_parameters) * math.log(
                fit.num_observations
            )
        else:
            sign, logdet = np.linalg.slogdet(fit.coefficient_covariance)
            variance = float(np.nanvar(data.values[:, j], ddof=1))
            if sign <= 0 or variance <= 0 or (dimension > 1 and fit.alpha_height == 0):
                return {
                    "criterion": criterion,
                    "score": None,
                    "status": "singular_optimum_information",
                }
            # Native offset coefficients equal optimum offsets times (1-exp(-alpha*H)).
            jacobian = (
                0.0
                if dimension == 1
                else 2 * (dimension - 1) * math.log(-math.expm1(-fit.alpha_height))
            )
            penalty += (
                covariance_parameters * math.log(fit.num_observations)
                - logdet
                + jacobian
                + dimension * math.log(variance)
            )
    correction = {}
    if criterion == "AICc":
        # Match kfl1ou's independent-trait convention: shared locations count
        # once, and observed scalar coordinates sum across traits.
        observations = sum(fit.num_observations for fit in result["fits"])
        denominator = observations - parameters - 1
        correction = {"sample_size": observations, "small_sample_correction": None}
        if denominator <= 0:
            return {
                "criterion": criterion,
                "score": None,
                "penalty": None,
                "parameter_count": parameters,
                "status": "insufficient_aicc_sample_size",
                **correction,
            }
        extra = 2 * parameters * (parameters + 1) / denominator
        penalty += extra
        correction["small_sample_correction"] = extra
    score = -2 * result["log_likelihood"] + penalty
    return {
        "criterion": criterion,
        "score": float(score) if math.isfinite(score) else None,
        "penalty": float(penalty) if math.isfinite(penalty) else None,
        "parameter_count": parameters,
        "status": "ok" if math.isfinite(score) else "nonfinite_information",
        **correction,
    }
