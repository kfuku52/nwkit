import math
from dataclasses import dataclass

import numpy as np
import pandas as pd

TIP_SUMMARY_COLUMNS = [
    "tree_id",
    "leaf_name",
    "trait",
    "n_biological",
    "n_technical",
    "mean",
    "within_sd",
    "standard_error",
    "variance_method",
    "batch_adjusted",
]


@dataclass
class ReplicateEstimates:
    values_by_trait: dict[str, dict[str, float]]
    sampling_covariance_by_trait: dict[str, pd.DataFrame]
    tip_summary: pd.DataFrame
    model_by_trait: dict[str, str]


def _require_nonempty_ids(dataframe, column, option_name):
    if column not in dataframe.columns:
        raise ValueError("Column '{}' is required by '{}'.".format(column, option_name))
    values = dataframe[column]
    invalid = values.isna() | values.astype(str).str.strip().eq("")
    if invalid.any():
        raise ValueError(
            "Column '{}' must not contain missing or empty values.".format(column)
        )


def _numeric_trait(dataframe, trait):
    original = dataframe[trait]
    converted = pd.to_numeric(original, errors="coerce")
    invalid = original.notna() & converted.isna()
    if invalid.any():
        leaves = sorted(set(dataframe.loc[invalid, "leaf_name"].astype(str)))
        raise ValueError(
            "Trait column '{}' has non-numeric values for tree tips: {}.".format(
                trait, ", ".join(leaves)
            )
        )
    finite = converted.notna() & ~np.isfinite(
        converted.to_numpy(dtype=float, na_value=np.nan)
    )
    if finite.any():
        leaves = sorted(set(dataframe.loc[finite, "leaf_name"].astype(str)))
        raise ValueError(
            "Trait column '{}' has non-finite values for tree tips: {}.".format(
                trait, ", ".join(leaves)
            )
        )
    return converted.astype(float)


def _validate_leaf_coverage(dataframe, leaf_names, trait):
    observed = set(dataframe.loc[dataframe[trait].notna(), "leaf_name"].astype(str))
    missing = sorted(set(leaf_names) - observed)
    if missing:
        raise ValueError(
            "Trait column '{}' has no biological observations for tree tips: {}.".format(
                trait, ", ".join(missing)
            )
        )


def _aggregate_technical_replicates(
    dataframe,
    traits,
    biological_id,
    technical_id,
    batch,
    technical_aggregation,
):
    key = ["leaf_name", biological_id]
    duplicated = dataframe.duplicated(subset=key, keep=False)
    if not duplicated.any():
        result = dataframe.copy()
        for trait in traits:
            result["_n_technical_{}".format(trait)] = result[trait].notna().astype(int)
        return result
    if technical_id is None:
        examples = dataframe.loc[duplicated, key].drop_duplicates().head(5)
        labels = ["{}:{}".format(*row) for row in examples.to_numpy()]
        raise ValueError(
            "Repeated leaf_name/biological ID rows require '--technical-id' "
            "or prior aggregation (examples: {}).".format(", ".join(labels))
        )
    _require_nonempty_ids(dataframe, technical_id, "--technical-id")
    duplicated_technical = dataframe.duplicated(subset=key + [technical_id], keep=False)
    if duplicated_technical.any():
        raise ValueError(
            "Technical IDs must be unique within each leaf_name/biological ID."
        )
    if technical_aggregation == "error":
        raise ValueError(
            "Technical replicates are present; select '--technical-aggregation mean' "
            "to aggregate transformed continuous values explicitly."
        )
    if technical_aggregation != "mean":
        raise ValueError(
            "Unsupported technical replicate aggregation: {}.".format(
                technical_aggregation
            )
        )
    if batch is not None:
        batch_counts = dataframe.groupby(key, sort=False, dropna=False)[batch].nunique()
        if (batch_counts > 1).any():
            raise ValueError(
                "Technical replicates for one biological observation cannot span batches."
            )
    aggregations = {trait: "mean" for trait in traits}
    if batch is not None:
        aggregations[batch] = "first"
    result = dataframe.groupby(key, sort=False, as_index=False).agg(aggregations)
    counts = dataframe.groupby(key, sort=False)[traits].count().reset_index()
    counts = counts.rename(
        columns={trait: "_n_technical_{}".format(trait) for trait in traits}
    )
    result = result.merge(counts, on=key, how="left", validate="one_to_one")
    return result


def _known_se_estimates(
    dataframe,
    leaf_names,
    traits,
    se_columns,
    n_columns,
    tree_id,
):
    dataframe = dataframe[dataframe[traits].notna().any(axis=1)].copy()
    if dataframe["leaf_name"].duplicated().any():
        raise ValueError("Known-SE input requires exactly one row per leaf_name.")
    by_leaf = dataframe.set_index("leaf_name").reindex(leaf_names)
    values_by_trait = {}
    covariance_by_trait = {}
    models = {}
    rows = []
    for index, trait in enumerate(traits):
        se_column = se_columns[index]
        if se_column not in dataframe.columns:
            raise ValueError("Known-SE column '{}' is absent.".format(se_column))
        values = pd.to_numeric(by_leaf[trait], errors="coerce")
        ses = pd.to_numeric(by_leaf[se_column], errors="coerce")
        if values.isna().any() or ses.isna().any():
            raise ValueError(
                "Known-SE input requires finite means and standard errors for trait '{}'.".format(
                    trait
                )
            )
        values_array = values.to_numpy(dtype=float)
        se_array = ses.to_numpy(dtype=float)
        if not np.isfinite(values_array).all() or not np.isfinite(se_array).all():
            raise ValueError(
                "Known-SE input requires finite means and standard errors for trait '{}'.".format(
                    trait
                )
            )
        if np.any(se_array < 0.0):
            raise ValueError("Known standard errors must be non-negative.")
        if n_columns:
            n_column = n_columns[index]
            if n_column not in dataframe.columns:
                raise ValueError("Sample-size column '{}' is absent.".format(n_column))
            n_values = pd.to_numeric(by_leaf[n_column], errors="coerce")
            if n_values.isna().any() or not np.isfinite(n_values.to_numpy(float)).all():
                raise ValueError("Known sample sizes must be positive integers.")
            n_array = n_values.to_numpy(dtype=float)
            if np.any(n_array <= 0.0) or np.any(n_array != np.floor(n_array)):
                raise ValueError("Known sample sizes must be positive integers.")
            n_array = n_array.astype(int)
        else:
            n_array = None
        covariance = pd.DataFrame(
            np.diag(se_array**2), index=leaf_names, columns=leaf_names
        )
        values_by_trait[trait] = dict(zip(leaf_names, values_array, strict=True))
        covariance_by_trait[trait] = covariance
        models[trait] = "known-se"
        for leaf_index, leaf_name in enumerate(leaf_names):
            sample_size = "" if n_array is None else int(n_array[leaf_index])
            within_sd = (
                ""
                if n_array is None
                else se_array[leaf_index] * math.sqrt(n_array[leaf_index])
            )
            rows.append(
                {
                    "tree_id": tree_id,
                    "leaf_name": leaf_name,
                    "trait": trait,
                    "n_biological": sample_size,
                    "n_technical": "",
                    "mean": values_array[leaf_index],
                    "within_sd": within_sd,
                    "standard_error": se_array[leaf_index],
                    "variance_method": "known-se",
                    "batch_adjusted": "no",
                }
            )
    return ReplicateEstimates(
        values_by_trait=values_by_trait,
        sampling_covariance_by_trait=covariance_by_trait,
        tip_summary=pd.DataFrame(rows, columns=TIP_SUMMARY_COLUMNS),
        model_by_trait=models,
    )


def _pooled_no_batch(dataframe, leaf_names, trait):
    selected = dataframe[dataframe[trait].notna()].copy()
    grouped = selected.groupby("leaf_name", sort=False)[trait]
    means = grouped.mean().reindex(leaf_names)
    counts = grouped.size().reindex(leaf_names).astype(int)
    residuals = selected[trait] - selected["leaf_name"].map(means)
    degrees_of_freedom = len(selected) - len(leaf_names)
    if degrees_of_freedom <= 0:
        raise ValueError(
            "Pooled within-leaf variance for trait '{}' needs at least one "
            "residual biological-replicate degree of freedom.".format(trait)
        )
    variance = float(np.dot(residuals, residuals) / degrees_of_freedom)
    covariance = np.diag(variance / counts.to_numpy(dtype=float))
    within_sd = np.repeat(math.sqrt(max(variance, 0.0)), len(leaf_names))
    return means.to_numpy(float), covariance, counts.to_numpy(int), within_sd


def _leaf_specific_no_batch(dataframe, leaf_names, trait):
    selected = dataframe[dataframe[trait].notna()].copy()
    grouped = selected.groupby("leaf_name", sort=False)[trait]
    means = grouped.mean().reindex(leaf_names)
    counts = grouped.size().reindex(leaf_names).astype(int)
    if (counts < 2).any():
        missing = counts.index[counts < 2].tolist()
        raise ValueError(
            "Leaf-specific variance for trait '{}' needs at least two biological "
            "replicates per leaf (insufficient: {}).".format(trait, ", ".join(missing))
        )
    variances = grouped.var(ddof=1).reindex(leaf_names).to_numpy(float)
    covariance = np.diag(variances / counts.to_numpy(dtype=float))
    return (
        means.to_numpy(float),
        covariance,
        counts.to_numpy(int),
        np.sqrt(np.maximum(variances, 0.0)),
    )


def _pooled_with_batch(dataframe, leaf_names, trait, batch):
    selected = dataframe[dataframe[trait].notna()].copy()
    if (
        selected[batch].isna().any()
        or selected[batch].astype(str).str.strip().eq("").any()
    ):
        raise ValueError("Batch values must be present for every observed trait value.")
    leaf_levels = list(leaf_names)
    batch_levels = sorted(set(selected[batch].astype(str)))
    leaf_design = np.column_stack(
        [
            (selected["leaf_name"].astype(str) == leaf).to_numpy(float)
            for leaf in leaf_levels
        ]
    )
    if len(batch_levels) > 1:
        batch_design = np.column_stack(
            [
                (selected[batch].astype(str) == level).to_numpy(float)
                for level in batch_levels[1:]
            ]
        )
        design = np.column_stack([leaf_design, batch_design])
    else:
        batch_design = np.empty((len(selected), 0), dtype=float)
        design = leaf_design
    rank = int(np.linalg.matrix_rank(design))
    if rank != design.shape[1]:
        raise ValueError(
            "Batch and leaf effects are confounded for trait '{}'; the observation "
            "design matrix is rank deficient.".format(trait)
        )
    degrees_of_freedom = len(selected) - rank
    if degrees_of_freedom <= 0:
        raise ValueError(
            "Batch-adjusted variance for trait '{}' needs positive residual degrees "
            "of freedom.".format(trait)
        )
    response = selected[trait].to_numpy(float)
    gram = design.T @ design
    coefficients = np.linalg.solve(gram, design.T @ response)
    residuals = response - design @ coefficients
    variance = float(np.dot(residuals, residuals) / degrees_of_freedom)
    coefficient_covariance = variance * np.linalg.inv(gram)
    transform = np.zeros((len(leaf_levels), design.shape[1]), dtype=float)
    transform[:, : len(leaf_levels)] = np.eye(len(leaf_levels))
    if batch_design.shape[1]:
        average_batch = batch_design.mean(axis=0)
        transform[:, len(leaf_levels) :] = average_batch
    adjusted_means = transform @ coefficients
    mean_covariance = transform @ coefficient_covariance @ transform.T
    counts = (
        selected.groupby("leaf_name", sort=False)
        .size()
        .reindex(leaf_levels)
        .to_numpy(int)
    )
    within_sd = np.repeat(math.sqrt(max(variance, 0.0)), len(leaf_levels))
    return adjusted_means, mean_covariance, counts, within_sd


def _validate_replicate_request(
    dataframe,
    tree_leaf_names,
    traits,
    biological_id,
    technical_id,
    batch,
    within_variance,
):
    leaf_names = [str(name) for name in tree_leaf_names]
    traits = list(traits)
    if "leaf_name" not in dataframe.columns:
        raise ValueError("Replicate input requires a 'leaf_name' column.")
    if not leaf_names or any(name.strip() == "" for name in leaf_names):
        raise ValueError("Tree leaf names must be non-empty.")
    if len(leaf_names) != len(set(leaf_names)):
        raise ValueError("Tree leaf names must be unique.")
    if not traits or any(
        not isinstance(trait, str) or trait.strip() == "" for trait in traits
    ):
        raise ValueError("Traits must be a non-empty sequence of non-empty names.")
    if len(traits) != len(set(traits)) or "leaf_name" in traits:
        raise ValueError("Trait names must be unique and cannot be 'leaf_name'.")
    if within_variance not in {"known-se", "leaf", "pooled"}:
        raise ValueError(
            "Unsupported within-leaf variance method: {}.".format(within_variance)
        )
    for trait in traits:
        if trait not in dataframe.columns:
            raise ValueError("Trait column '{}' is absent.".format(trait))
    identifiers = [
        value for value in [biological_id, technical_id, batch] if value is not None
    ]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError(
            "Biological-ID, technical-ID, and batch columns must be distinct."
        )
    conflicting_roles = sorted(set(identifiers).intersection(traits + ["leaf_name"]))
    if conflicting_roles:
        raise ValueError(
            "Identifier/batch columns cannot also be trait or leaf-name columns: {}.".format(
                ", ".join(conflicting_roles)
            )
        )
    return leaf_names, traits


def _validate_known_se_request(
    traits,
    biological_id,
    technical_id,
    batch,
    technical_aggregation,
    se_columns,
    n_columns,
):
    incompatible = [
        label
        for label, value in [
            ("--biological-id", biological_id),
            ("--technical-id", technical_id),
            ("--batch", batch),
        ]
        if value is not None
    ]
    if technical_aggregation != "error":
        incompatible.append("--technical-aggregation")
    if incompatible:
        raise ValueError(
            "Known-SE input cannot use raw-replicate option(s): {}.".format(
                ", ".join(incompatible)
            )
        )
    if len(se_columns) != len(traits):
        raise ValueError("Known-SE input requires one standard-error column per trait.")
    if n_columns and len(n_columns) != len(traits):
        raise ValueError("Provide one sample-size column per trait or none.")


def _prepare_raw_replicates(
    dataframe,
    leaf_names,
    traits,
    biological_id,
    technical_id,
    batch,
    within_variance,
    technical_aggregation,
    se_columns,
    n_columns,
):
    if se_columns or n_columns:
        raise ValueError(
            "Standard-error and sample-size columns require known-SE input."
        )
    if biological_id is None:
        raise ValueError(
            "Raw replicate input requires '--biological-id', or select "
            "'--within-variance known-se'."
        )
    for trait in traits:
        dataframe[trait] = _numeric_trait(dataframe, trait)
    dataframe = dataframe[dataframe[traits].notna().any(axis=1)].copy()
    if dataframe.empty:
        raise ValueError("Replicate input contains no observed trait values.")
    _require_nonempty_ids(dataframe, biological_id, "--biological-id")
    if technical_aggregation not in {"error", "mean"}:
        raise ValueError(
            "Unsupported technical replicate aggregation: {}.".format(
                technical_aggregation
            )
        )
    if technical_id is None and technical_aggregation != "error":
        raise ValueError("'--technical-aggregation' requires '--technical-id'.")
    if technical_id is not None:
        _require_nonempty_ids(dataframe, technical_id, "--technical-id")
    if batch is not None:
        _require_nonempty_ids(dataframe, batch, "--batch")
        if within_variance == "leaf":
            raise ValueError(
                "Leaf-specific variance with batch adjustment is not supported; "
                "use pooled variance or known standard errors."
            )
    dataframe = _aggregate_technical_replicates(
        dataframe,
        traits,
        biological_id,
        technical_id,
        batch,
        technical_aggregation,
    )
    for trait in traits:
        _validate_leaf_coverage(dataframe, leaf_names, trait)
    return dataframe


def _estimate_one_trait(dataframe, leaf_names, trait, batch, within_variance):
    if batch is not None:
        estimates = _pooled_with_batch(dataframe, leaf_names, trait, batch)
        method = "pooled-batch-adjusted"
    elif within_variance == "leaf":
        estimates = _leaf_specific_no_batch(dataframe, leaf_names, trait)
        method = "leaf"
    else:
        estimates = _pooled_no_batch(dataframe, leaf_names, trait)
        method = "pooled"
    means, covariance, counts, within_sd = estimates
    covariance = (np.asarray(covariance) + np.asarray(covariance).T) / 2.0
    if not np.isfinite(means).all() or not np.isfinite(covariance).all():
        raise ValueError("Replicate estimation produced non-finite values.")
    minimum_eigenvalue = float(np.linalg.eigvalsh(covariance).min())
    tolerance = np.finfo(float).eps * max(1.0, float(np.max(np.abs(covariance))))
    if minimum_eigenvalue < -tolerance:
        raise ValueError("Replicate estimation produced an invalid covariance matrix.")
    covariance[np.abs(covariance) < tolerance] = 0.0
    return means, covariance, counts, within_sd, method


def _replicate_tip_rows(
    dataframe,
    leaf_names,
    trait,
    means,
    covariance,
    counts,
    within_sd,
    method,
    batch,
    tree_id,
):
    selected = dataframe[dataframe[trait].notna()]
    technical_counts = (
        selected.groupby("leaf_name", sort=False)["_n_technical_{}".format(trait)]
        .sum()
        .reindex(leaf_names)
        .to_numpy(int)
    )
    standard_errors = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    rows = []
    for leaf_index, leaf_name in enumerate(leaf_names):
        rows.append(
            {
                "tree_id": tree_id,
                "leaf_name": leaf_name,
                "trait": trait,
                "n_biological": int(counts[leaf_index]),
                "n_technical": int(technical_counts[leaf_index]),
                "mean": float(means[leaf_index]),
                "within_sd": float(within_sd[leaf_index]),
                "standard_error": float(standard_errors[leaf_index]),
                "variance_method": method,
                "batch_adjusted": "yes" if batch is not None else "no",
            }
        )
    return rows


def estimate_replicate_traits(
    dataframe,
    tree_leaf_names,
    traits,
    *,
    biological_id=None,
    technical_id=None,
    batch=None,
    within_variance="pooled",
    technical_aggregation="error",
    se_columns=(),
    n_columns=(),
    tree_id="",
):
    leaf_names, traits = _validate_replicate_request(
        dataframe,
        tree_leaf_names,
        traits,
        biological_id,
        technical_id,
        batch,
        within_variance,
    )
    dataframe = dataframe.copy()
    _require_nonempty_ids(dataframe, "leaf_name", "replicate input")
    dataframe["leaf_name"] = dataframe["leaf_name"].astype(str)
    dataframe = dataframe[dataframe["leaf_name"].isin(set(leaf_names))].copy()
    if within_variance == "known-se":
        _validate_known_se_request(
            traits,
            biological_id,
            technical_id,
            batch,
            technical_aggregation,
            se_columns,
            n_columns,
        )
        return _known_se_estimates(
            dataframe, leaf_names, traits, list(se_columns), list(n_columns), tree_id
        )
    dataframe = _prepare_raw_replicates(
        dataframe,
        leaf_names,
        traits,
        biological_id,
        technical_id,
        batch,
        within_variance,
        technical_aggregation,
        se_columns,
        n_columns,
    )
    values_by_trait = {}
    covariance_by_trait = {}
    models = {}
    rows = []
    for trait in traits:
        means, covariance, counts, within_sd, method = _estimate_one_trait(
            dataframe, leaf_names, trait, batch, within_variance
        )
        covariance_frame = pd.DataFrame(
            covariance, index=leaf_names, columns=leaf_names
        )
        values_by_trait[trait] = dict(zip(leaf_names, means, strict=True))
        covariance_by_trait[trait] = covariance_frame
        models[trait] = method
        rows.extend(
            _replicate_tip_rows(
                dataframe,
                leaf_names,
                trait,
                means,
                covariance,
                counts,
                within_sd,
                method,
                batch,
                tree_id,
            )
        )
    return ReplicateEstimates(
        values_by_trait=values_by_trait,
        sampling_covariance_by_trait=covariance_by_trait,
        tip_summary=pd.DataFrame(rows, columns=TIP_SUMMARY_COLUMNS),
        model_by_trait=models,
    )
