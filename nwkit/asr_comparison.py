"""Information-criterion summaries for compatible ASR likelihood fits."""

import math

import numpy as np
import pandas as pd


def _estimated_mapping_count(fit, value_name, flag_name, *, shared_name=None):
    if not getattr(fit, flag_name, False):
        return 0
    values = getattr(fit, value_name, None)
    if values is None:
        return 1
    if shared_name is not None and getattr(fit, shared_name, None) is not None:
        return 1
    return len(values)


def _continuous_parameter_count(fit):
    model = getattr(fit, "model", "")
    if model in {"MV-BM", "MV-OU"}:
        dimension = len(fit.trait_names)
        covariance_count = dimension * (dimension + 1) // 2
        if model == "MV-BM":
            return covariance_count
        theta_count = dimension if getattr(fit, "theta_estimated", False) else 0
        return (
            covariance_count + theta_count + int(getattr(fit, "alpha_estimated", False))
        )
    if hasattr(fit, "evolution_parameter_estimated"):
        return int(fit.sigma2_estimated) + int(fit.evolution_parameter_estimated)
    if hasattr(fit, "eb_rate_estimated"):
        return int(fit.sigma2_estimated) + int(fit.eb_rate_estimated)
    if hasattr(fit, "drift_by_regime"):
        return _estimated_mapping_count(
            fit, "sigma2_by_regime", "sigma2_estimated"
        ) + _estimated_mapping_count(fit, "drift_by_regime", "drift_estimated")
    if hasattr(fit, "drift_estimated"):
        return int(fit.sigma2_estimated) + int(fit.drift_estimated)
    if hasattr(fit, "theta_by_regime"):
        return sum(
            (
                _estimated_mapping_count(
                    fit,
                    "alpha_by_regime",
                    "alpha_estimated",
                    shared_name="alpha",
                ),
                _estimated_mapping_count(fit, "theta_by_regime", "theta_estimated"),
                _estimated_mapping_count(
                    fit,
                    "sigma2_by_regime",
                    "sigma2_estimated",
                    shared_name="sigma2",
                ),
            )
        )
    if hasattr(fit, "alpha_estimated") and hasattr(fit, "theta_estimated"):
        sigma_count = int(getattr(fit, "sigma2_estimated", False))
        return int(fit.alpha_estimated) + sigma_count + int(fit.theta_estimated)
    if hasattr(fit, "sigma2_by_regime"):
        return len(fit.sigma2_by_regime) if fit.sigma2_estimated else 0
    if hasattr(fit, "sigma"):
        dimension = len(fit.trait_names)
        return dimension * (dimension + 1) // 2
    return int(getattr(fit, "sigma2_estimated", False))


def _discrete_parameter_count(fit):
    model = fit["model"]
    if model == "COVARION":
        return 2 + int(fit.get("base_rate_estimated", True))
    if model == "MK-MIXTURE":
        base = len(fit["rates"]) if fit.get("base_rate_estimated", True) else 0
        categories = int(fit["rate_categories"])
        mixture = 1 if fit["rate_mixture"] == "gamma" else 2 * (categories - 1)
        return base + mixture
    if not fit.get("rate_estimated", False):
        return 0
    return len(fit.get("rates", ()))


def summarize_fit(model, fit, *, trait_type):
    if trait_type == "discrete":
        likelihood = fit.get("log_likelihood")
        return {
            "model": model,
            "likelihood_kind": "discrete_ml",
            "log_likelihood": likelihood,
            "num_parameters": _discrete_parameter_count(fit),
            "sample_size": int(fit.get("sample_size", fit.get("num_characters", 1))),
            "fit_status": fit.get("fit_status", "ok"),
        }
    restricted = getattr(fit, "restricted_log_likelihood", None)
    if restricted is not None:
        likelihood = restricted
        likelihood_kind = "flat_root_integrated"
        sample_size = int(fit.num_effective_observations)
    else:
        likelihood = getattr(fit, "log_likelihood", None)
        likelihood_kind = (
            "stationary_root_ml"
            if getattr(fit, "root_prior", "stationary") == "stationary"
            else "proper_root_ml"
        )
        sample_size = int(fit.num_effective_observations)
    return {
        "model": model,
        "likelihood_kind": likelihood_kind,
        "log_likelihood": likelihood,
        "num_parameters": _continuous_parameter_count(fit),
        "sample_size": sample_size,
        "fit_status": getattr(fit, "fit_status", "ok"),
    }


def model_comparison_table(summaries):
    """Calculate AIC/AICc/BIC and weights for one likelihood convention."""

    rows = [dict(summary) for summary in summaries]
    if len(rows) < 2:
        raise ValueError("ASR model comparison requires at least two fitted models.")
    kinds = {row["likelihood_kind"] for row in rows}
    if len(kinds) != 1:
        raise ValueError(
            "ASR models with different likelihood conventions cannot be compared: "
            + ", ".join(sorted(kinds))
        )
    for row in rows:
        likelihood = row["log_likelihood"]
        fit_status = str(row.get("fit_status", "ok"))
        singular = "singular" in fit_status
        covarion_boundary = row["model"] == "COVARION" and "boundary" in fit_status
        if singular or covarion_boundary:
            raise ValueError(
                f"Model {row['model']} has non-regular fit status {fit_status!r} "
                "and cannot be ranked by regular-model information criteria."
            )
        if likelihood is None or not math.isfinite(float(likelihood)):
            raise ValueError(
                f"Model {row['model']} has no finite ordinary likelihood for comparison."
            )
        count = int(row["num_parameters"])
        sample_size = int(row["sample_size"])
        if count < 0:
            raise ValueError(f"Model {row['model']} has a negative parameter count.")
        if sample_size <= 0:
            raise ValueError(f"Model {row['model']} has a non-positive sample size.")
        row["aic"] = 2.0 * count - 2.0 * float(likelihood)
        denominator = sample_size - count - 1
        row["aicc"] = (
            row["aic"] + 2.0 * count * (count + 1) / denominator
            if denominator > 0
            else np.nan
        )
        row["bic"] = math.log(sample_size) * count - 2.0 * float(likelihood)
    for criterion in ("aic", "aicc", "bic"):
        finite = [row[criterion] for row in rows if math.isfinite(row[criterion])]
        if not finite:
            for row in rows:
                row[f"delta_{criterion}"] = np.nan
                row[f"{criterion}_weight"] = np.nan
            continue
        minimum = min(finite)
        raw = [
            math.exp(-0.5 * (row[criterion] - minimum))
            if math.isfinite(row[criterion])
            else 0.0
            for row in rows
        ]
        total = math.fsum(raw)
        for row, value in zip(rows, raw, strict=True):
            row[f"delta_{criterion}"] = (
                row[criterion] - minimum if math.isfinite(row[criterion]) else np.nan
            )
            row[f"{criterion}_weight"] = value / total if value else 0.0
    return pd.DataFrame(rows).sort_values("aic", kind="stable").reset_index(drop=True)
