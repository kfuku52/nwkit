"""Information-criterion summaries for compatible ASR likelihood fits."""

import math

import numpy as np
import pandas as pd

from nwkit.asr_models import model_definition

INFORMATION_CRITERIA = ("aic", "aicc", "bic")
IC_TIE_ABS_TOLERANCE = 1e-9
IC_TIE_REL_TOLERANCE = 0.0


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


def _continuous_root_prior(model, fit, explicit_root_prior):
    if explicit_root_prior is not None:
        return explicit_root_prior
    fit_root_prior = getattr(fit, "root_prior", None)
    if fit_root_prior is not None:
        return str(fit_root_prior)
    try:
        return model_definition(model).default_root_prior
    except ValueError:
        # Synthetic/test fit objects may not carry a registered model name.  In
        # that case the presence of an ordinary likelihood distinguishes the
        # proper-root contract from a flat-root integrated fit.
        return (
            "stationary" if getattr(fit, "log_likelihood", None) is not None else "flat"
        )


def summarize_fit(model, fit, *, trait_type, root_prior=None):
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
    resolved_root_prior = _continuous_root_prior(model, fit, root_prior)
    if resolved_root_prior == "flat":
        likelihood = getattr(fit, "restricted_log_likelihood", None)
        likelihood_kind = "flat_root_integrated"
        sample_size = int(fit.num_effective_observations)
    else:
        likelihood = getattr(fit, "log_likelihood", None)
        likelihood_kind = (
            "stationary_root_ml"
            if resolved_root_prior == "stationary"
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


def _information_criteria(row):
    likelihood = row["log_likelihood"]
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
    aic = 2.0 * count - 2.0 * float(likelihood)
    denominator = sample_size - count - 1
    return {
        "aic": aic,
        "aicc": (
            aic + 2.0 * count * (count + 1) / denominator if denominator > 0 else np.nan
        ),
        "bic": math.log(sample_size) * count - 2.0 * float(likelihood),
    }


def _criterion_weights(rows, criterion, *, minimum_finite=1):
    finite = [
        (index, float(row[criterion]))
        for index, row in enumerate(rows)
        if math.isfinite(float(row[criterion]))
    ]
    if len(finite) < minimum_finite:
        return None
    minimum = min(value for _, value in finite)
    deltas = {
        index: (
            0.0
            if math.isclose(
                value,
                minimum,
                rel_tol=IC_TIE_REL_TOLERANCE,
                abs_tol=IC_TIE_ABS_TOLERANCE,
            )
            else value - minimum
        )
        for index, value in finite
    }
    raw = {index: math.exp(-0.5 * deltas[index]) for index, _value in finite}
    total = math.fsum(raw.values())
    return {index: (deltas[index], raw[index] / total) for index, _value in finite}


def _assign_dense_ranks(scored, finite_selected):
    rank = 0
    cluster_value = None
    for value, _order, index in finite_selected:
        if cluster_value is None or not math.isclose(
            value,
            cluster_value,
            rel_tol=IC_TIE_REL_TOLERANCE,
            abs_tol=IC_TIE_ABS_TOLERANCE,
        ):
            rank += 1
            cluster_value = value
        scored[index]["criterion_value"] = value
        scored[index]["criterion_rank"] = rank
        scored[index]["is_best"] = "yes" if rank == 1 else "no"


def _score_group(scored, group, indices, criterion):
    group_rows = [scored[index] for index in indices]
    for name in INFORMATION_CRITERIA:
        weights = _criterion_weights(group_rows, name, minimum_finite=2)
        if weights is None:
            continue
        for relative_index, (delta, weight) in weights.items():
            row = scored[indices[relative_index]]
            row[f"delta_{name}"] = delta
            row[f"{name}_weight"] = weight

    finite_selected = sorted(
        (
            (float(scored[index][criterion]), order, index)
            for order, index in enumerate(indices)
            if math.isfinite(float(scored[index][criterion]))
        ),
        key=lambda item: (item[0], item[1]),
    )
    for row in scored:
        if str(row.get("comparison_group", "")) == group:
            row["num_comparable_models"] = len(finite_selected)
    if len(finite_selected) >= 2:
        _assign_dense_ranks(scored, finite_selected)


def grouped_model_comparison_table(rows, *, criterion="aic"):
    """Score rankable rows within explicitly compatible comparison groups.

    Unlike :func:`model_comparison_table`, this batch-oriented interface retains
    failed and non-rankable candidates.  Deltas, weights, and a selected winner
    are emitted only when at least two finite candidates share a group.
    """

    if criterion not in INFORMATION_CRITERIA:
        raise ValueError(
            "Unsupported model-comparison criterion: {}. Choose {}.".format(
                criterion, ", ".join(INFORMATION_CRITERIA)
            )
        )
    scored = [dict(row) for row in rows]
    metric_columns = [
        *INFORMATION_CRITERIA,
        *(f"delta_{name}" for name in INFORMATION_CRITERIA),
        *(f"{name}_weight" for name in INFORMATION_CRITERIA),
    ]
    for row in scored:
        for name in metric_columns:
            row[name] = np.nan
        row["num_comparable_models"] = 0
        row["criterion"] = criterion
        row["criterion_value"] = np.nan
        row["criterion_rank"] = np.nan
        row["is_best"] = ""

    grouped: dict[str, list[int]] = {}
    for index, row in enumerate(scored):
        if str(row.get("rankable", "no")).lower() not in {"yes", "true", "1"}:
            continue
        row.update(_information_criteria(row))
        group = str(row.get("comparison_group", ""))
        if group == "":
            raise ValueError(f"Rankable model {row['model']} has no comparison_group.")
        grouped.setdefault(group, []).append(index)

    for group, indices in grouped.items():
        _score_group(scored, group, indices, criterion)

    return pd.DataFrame(scored)


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
        fit_status = str(row.get("fit_status", "ok"))
        singular = "singular" in fit_status
        covarion_boundary = row["model"] == "COVARION" and "boundary" in fit_status
        if singular or covarion_boundary:
            raise ValueError(
                f"Model {row['model']} has non-regular fit status {fit_status!r} "
                "and cannot be ranked by regular-model information criteria."
            )
        row.update(_information_criteria(row))
    for criterion in INFORMATION_CRITERIA:
        weights = _criterion_weights(rows, criterion)
        if weights is None:
            for row in rows:
                row[f"delta_{criterion}"] = np.nan
                row[f"{criterion}_weight"] = np.nan
            continue
        for index, row in enumerate(rows):
            delta, weight = weights.get(index, (np.nan, 0.0))
            row[f"delta_{criterion}"] = delta
            row[f"{criterion}_weight"] = weight
    return pd.DataFrame(rows).sort_values("aic", kind="stable").reset_index(drop=True)
