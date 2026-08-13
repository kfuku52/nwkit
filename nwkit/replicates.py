import math
from dataclasses import dataclass

import numpy as np
import pandas as pd

from nwkit.model_matrix import CategoricalObservation, ReplicatedObservation

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
    "state",
    "state_counts",
    "state_probabilities",
    "n_levels",
]


@dataclass
class ReplicateEstimates:
    values_by_trait: dict[str, dict[str, object]]
    sampling_covariance_by_trait: dict[str, pd.DataFrame | np.ndarray]
    tip_summary: pd.DataFrame
    model_by_trait: dict[str, str]


def _allowed_missing_trait_set(traits, allow_missing_traits):
    allowed = set(allow_missing_traits)
    if allowed - set(traits):
        raise ValueError("Allowed-missing traits must be selected response traits.")
    return allowed


def _categorical_technical_observations(dataframe, traits, biological_id, technical_id):
    key = ["leaf_name", biological_id]
    duplicated = dataframe.duplicated(subset=key, keep=False)
    if duplicated.any() and technical_id is None:
        raise ValueError(
            "Repeated categorical leaf/biological-ID rows require '--technical-id'."
        )
    if technical_id is not None:
        _require_nonempty_ids(dataframe, technical_id, "--technical-id")
        if dataframe.duplicated(subset=key + [technical_id], keep=False).any():
            raise ValueError(
                "Technical IDs must be unique within each leaf_name/biological ID."
            )
    rows = []
    for group_key, group in dataframe.groupby(key, sort=False, dropna=False):
        row: dict[str, object] = {
            "leaf_name": str(group_key[0]),
            biological_id: str(group_key[1]),
        }
        for trait in traits:
            states = group[trait].dropna().astype(str)
            unique = states.unique().tolist()
            if len(unique) > 1:
                raise ValueError(
                    "Technical replicates disagree for categorical trait '{}' at "
                    "tree tip '{}'.".format(trait, group_key[0])
                )
            row[trait] = None if not unique else unique[0]
            row["_n_technical_{}".format(trait)] = len(states)
        rows.append(row)
    return pd.DataFrame(rows)


def estimate_categorical_traits(
    dataframe,
    tree_leaf_names,
    traits,
    *,
    biological_id,
    technical_id=None,
    policy="error",
    tree_id="",
):
    """Aggregate categorical biological replicates without numeric averaging."""
    if policy not in {"counts", "error", "latent"}:
        raise ValueError("Unsupported categorical replicate policy: {}.".format(policy))
    if biological_id is None:
        raise ValueError("Categorical replicate input requires a biological ID.")
    leaf_names = [str(name) for name in tree_leaf_names]
    traits = list(traits)
    required = ["leaf_name", biological_id]
    for column in required:
        _require_nonempty_ids(dataframe, column, "categorical replicate input")
    dataframe = dataframe.copy()
    dataframe["leaf_name"] = dataframe["leaf_name"].astype(str)
    dataframe = dataframe[dataframe["leaf_name"].isin(set(leaf_names))]
    observations = _categorical_technical_observations(
        dataframe, traits, biological_id, technical_id
    )
    values_by_trait: dict[str, dict[str, object]] = {}
    models = {}
    summary_rows = []
    for trait in traits:
        levels = sorted(set(observations[trait].dropna().astype(str)))
        if len(levels) < 2:
            raise ValueError(
                "Categorical trait '{}' needs at least two observed levels.".format(
                    trait
                )
            )
        values_by_leaf = {}
        for leaf_name in leaf_names:
            selected = observations.loc[
                observations["leaf_name"] == leaf_name, trait
            ].dropna()
            if selected.empty:
                raise ValueError(
                    "Categorical trait '{}' has no observations for tree tip '{}'.".format(
                        trait, leaf_name
                    )
                )
            counts = selected.astype(str).value_counts().reindex(levels, fill_value=0)
            observed_levels = counts[counts > 0].index.tolist()
            if policy == "error" and len(observed_levels) != 1:
                raise ValueError(
                    "Biological replicates disagree for categorical trait '{}' at "
                    "tree tip '{}'; use '--categorical-replicate-policy latent'.".format(
                        trait, leaf_name
                    )
                )
            probabilities = counts.to_numpy(dtype=float) / float(counts.sum())
            values_by_leaf[leaf_name] = (
                CategoricalObservation(
                    dict(zip(levels, probabilities, strict=True)),
                    int(counts.sum()),
                )
                if policy == "counts" or len(observed_levels) > 1
                else observed_levels[0]
            )
            technical_count = int(
                observations.loc[
                    observations["leaf_name"] == leaf_name,
                    "_n_technical_{}".format(trait),
                ].sum()
            )
            summary_rows.append(
                {
                    "tree_id": tree_id,
                    "leaf_name": leaf_name,
                    "trait": trait,
                    "n_biological": int(counts.sum()),
                    "n_technical": technical_count,
                    "mean": "",
                    "within_sd": "",
                    "standard_error": "",
                    "variance_method": "categorical-{}".format(policy),
                    "batch_adjusted": "no",
                    "state": observed_levels[0] if len(observed_levels) == 1 else "",
                    "state_counts": ";".join(
                        "{}:{}".format(level, int(counts[level])) for level in levels
                    ),
                    "state_probabilities": ";".join(
                        "{}:{:.17g}".format(level, probabilities[index])
                        for index, level in enumerate(levels)
                    ),
                    "n_levels": len(levels),
                }
            )
        values_by_trait[trait] = values_by_leaf
        models[trait] = "categorical-{}".format(policy)
    return ReplicateEstimates(
        values_by_trait=values_by_trait,
        sampling_covariance_by_trait={},
        tip_summary=pd.DataFrame(summary_rows, columns=TIP_SUMMARY_COLUMNS),
        model_by_trait=models,
    )


def estimate_likelihood_replicates(
    dataframe,
    tree_leaf_names,
    traits,
    *,
    biological_id,
    technical_id=None,
    technical_aggregation="error",
    allow_missing_traits=(),
    tree_id="",
):
    """Retain numeric replicates, including permitted censored observations."""
    if biological_id is None:
        raise ValueError("Non-Gaussian replicate input requires a biological ID.")
    if technical_aggregation not in {"error", "mean"}:
        raise ValueError(
            "Unsupported technical replicate aggregation: {}.".format(
                technical_aggregation
            )
        )
    leaf_names = [str(name) for name in tree_leaf_names]
    traits = list(traits)
    allow_missing_traits = _allowed_missing_trait_set(traits, allow_missing_traits)
    dataframe = dataframe.copy()
    _require_nonempty_ids(dataframe, "leaf_name", "non-Gaussian replicate input")
    _require_nonempty_ids(dataframe, biological_id, "--biological-id")
    dataframe["leaf_name"] = dataframe["leaf_name"].astype(str)
    dataframe = dataframe[dataframe["leaf_name"].isin(set(leaf_names))]
    for trait in traits:
        dataframe[trait] = _numeric_trait(dataframe, trait)
        if trait not in allow_missing_traits:
            _validate_leaf_coverage(dataframe, leaf_names, trait)
    observations = _aggregate_technical_replicates(
        dataframe,
        traits,
        biological_id,
        technical_id,
        None,
        technical_aggregation,
    )
    values_by_trait: dict[str, dict[str, object]] = {}
    summary_rows = []
    for trait in traits:
        values_by_leaf: dict[str, object] = {}
        for leaf_name in leaf_names:
            selected = observations.loc[observations["leaf_name"] == leaf_name, trait]
            if trait not in allow_missing_traits:
                selected = selected.dropna()
            values = tuple(float(value) for value in selected)
            if not values:
                raise ValueError(
                    "Trait '{}' has no observations for tree tip '{}'.".format(
                        trait, leaf_name
                    )
                )
            values_by_leaf[leaf_name] = ReplicatedObservation(values)
            technical_count = int(
                observations.loc[
                    observations["leaf_name"] == leaf_name,
                    "_n_technical_{}".format(trait),
                ].sum()
            )
            summary_rows.append(
                _likelihood_summary_row(
                    tree_id, leaf_name, trait, values, technical_count
                )
            )
        values_by_trait[trait] = values_by_leaf
    return ReplicateEstimates(
        values_by_trait=values_by_trait,
        sampling_covariance_by_trait={},
        tip_summary=pd.DataFrame(summary_rows, columns=TIP_SUMMARY_COLUMNS),
        model_by_trait={trait: "likelihood-replicates" for trait in traits},
    )


def _likelihood_summary_row(tree_id, leaf_name, trait, values, technical_count):
    finite_values = np.asarray(values)[np.isfinite(values)]
    return {
        "tree_id": tree_id,
        "leaf_name": leaf_name,
        "trait": trait,
        "n_biological": len(values),
        "n_technical": technical_count,
        "mean": float(np.mean(finite_values)) if finite_values.size else "",
        "within_sd": (
            float(np.std(finite_values, ddof=1)) if finite_values.size > 1 else 0.0
        ),
        "standard_error": "",
        "variance_method": "likelihood-replicates",
        "batch_adjusted": "no",
    }


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


def _known_se_trait_arrays(by_leaf, trait, se_column, allow_missing):
    if se_column not in by_leaf.columns:
        raise ValueError("Known-SE column '{}' is absent.".format(se_column))
    values = pd.to_numeric(by_leaf[trait], errors="coerce")
    ses = pd.to_numeric(by_leaf[se_column], errors="coerce")
    if (values.isna() != ses.isna()).any() or (
        not allow_missing and values.isna().any()
    ):
        raise ValueError(
            "Known-SE input requires paired finite means and standard errors "
            "for trait '{}'.".format(trait)
        )
    values_array = values.to_numpy(dtype=float)
    se_array = ses.to_numpy(dtype=float)
    observed = np.isfinite(values_array)
    if not observed.any():
        raise ValueError("Known-SE trait '{}' has no observations.".format(trait))
    if np.any(np.isinf(values_array)) or np.any(np.isinf(se_array)):
        raise ValueError(
            "Known-SE input requires finite means and standard errors for trait '{}'.".format(
                trait
            )
        )
    if np.any(se_array < 0.0):
        raise ValueError("Known standard errors must be non-negative.")
    return values_array, se_array, observed


def _known_sample_sizes(by_leaf, n_columns, index, observed):
    if not n_columns:
        return None
    n_column = n_columns[index]
    if n_column not in by_leaf.columns:
        raise ValueError("Sample-size column '{}' is absent.".format(n_column))
    n_array = pd.to_numeric(by_leaf[n_column], errors="coerce").to_numpy(float)
    if np.any(np.isinf(n_array)) or np.any(np.isnan(n_array) != ~observed):
        raise ValueError("Known sample sizes must be positive integers.")
    if np.any(n_array[observed] <= 0.0) or np.any(
        n_array[observed] != np.floor(n_array[observed])
    ):
        raise ValueError("Known sample sizes must be positive integers.")
    return n_array


def _known_se_summary_rows(
    leaf_names, trait, values, ses, observed, sample_sizes, tree_id
):
    rows = []
    for leaf_index, leaf_name in enumerate(leaf_names):
        has_sample_size = sample_sizes is not None and observed[leaf_index]
        rows.append(
            {
                "tree_id": tree_id,
                "leaf_name": leaf_name,
                "trait": trait,
                "n_biological": (
                    int(sample_sizes[leaf_index]) if has_sample_size else ""
                ),
                "n_technical": "",
                "mean": values[leaf_index],
                "within_sd": (
                    ses[leaf_index] * math.sqrt(sample_sizes[leaf_index])
                    if has_sample_size
                    else ""
                ),
                "standard_error": ses[leaf_index] if observed[leaf_index] else "",
                "variance_method": "known-se",
                "batch_adjusted": "no",
            }
        )
    return rows


def _known_se_observation_rows(dataframe, traits, se_columns, n_columns):
    for index, trait in enumerate(traits):
        se_column = se_columns[index]
        if se_column not in dataframe.columns:
            raise ValueError("Known-SE column '{}' is absent.".format(se_column))
        values = pd.to_numeric(dataframe[trait], errors="coerce")
        ses = pd.to_numeric(dataframe[se_column], errors="coerce")
        invalid_numeric = (dataframe[trait].notna() & values.isna()) | (
            dataframe[se_column].notna() & ses.isna()
        )
        if invalid_numeric.any() or (values.isna() != ses.isna()).any():
            raise ValueError(
                "Known-SE input requires paired finite means and standard errors "
                "for trait '{}'.".format(trait)
            )
        if n_columns:
            n_column = n_columns[index]
            if n_column not in dataframe.columns:
                raise ValueError("Sample-size column '{}' is absent.".format(n_column))
            sample_sizes = pd.to_numeric(dataframe[n_column], errors="coerce")
            invalid_sample_size = dataframe[n_column].notna() & sample_sizes.isna()
            if (
                invalid_sample_size.any()
                or (sample_sizes.isna() != values.isna()).any()
            ):
                raise ValueError("Known sample sizes must be positive integers.")
    return dataframe[dataframe[traits].notna().any(axis=1)].copy()


def _known_se_estimates(
    dataframe,
    leaf_names,
    traits,
    se_columns,
    n_columns,
    tree_id,
    allow_missing_traits,
):
    dataframe = _known_se_observation_rows(dataframe, traits, se_columns, n_columns)
    if dataframe["leaf_name"].duplicated().any():
        raise ValueError("Known-SE input requires exactly one row per leaf_name.")
    by_leaf = dataframe.set_index("leaf_name").reindex(leaf_names)
    values_by_trait = {}
    covariance_by_trait = {}
    models = {}
    rows = []
    for index, trait in enumerate(traits):
        values, ses, observed = _known_se_trait_arrays(
            by_leaf,
            trait,
            se_columns[index],
            trait in allow_missing_traits,
        )
        sample_sizes = _known_sample_sizes(by_leaf, n_columns, index, observed)
        covariance_diagonal = np.zeros(len(leaf_names), dtype=float)
        covariance_diagonal[observed] = ses[observed] ** 2
        covariance_by_trait[trait] = covariance_diagonal
        values_by_trait[trait] = dict(zip(leaf_names, values, strict=True))
        models[trait] = "known-se"
        rows.extend(
            _known_se_summary_rows(
                leaf_names,
                trait,
                values,
                ses,
                observed,
                sample_sizes,
                tree_id,
            )
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
    covariance = variance / counts.to_numpy(dtype=float)
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
    covariance = variances / counts.to_numpy(dtype=float)
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
    allow_missing_traits,
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
        if trait not in allow_missing_traits:
            _validate_leaf_coverage(dataframe, leaf_names, trait)
    return dataframe


def _estimate_one_trait(
    dataframe, leaf_names, trait, batch, within_variance, *, allow_missing=False
):
    observed_leaf_names = [
        leaf_name
        for leaf_name in leaf_names
        if dataframe.loc[dataframe["leaf_name"] == leaf_name, trait].notna().any()
    ]
    if not observed_leaf_names:
        raise ValueError("Trait '{}' contains no observations.".format(trait))
    fitted_leaf_names = observed_leaf_names if allow_missing else leaf_names
    if batch is not None:
        estimates = _pooled_with_batch(dataframe, fitted_leaf_names, trait, batch)
        method = "pooled-batch-adjusted"
    elif within_variance == "leaf":
        estimates = _leaf_specific_no_batch(dataframe, fitted_leaf_names, trait)
        method = "leaf"
    else:
        estimates = _pooled_no_batch(dataframe, fitted_leaf_names, trait)
        method = "pooled"
    means, covariance, counts, within_sd = estimates
    if fitted_leaf_names != leaf_names:
        means, covariance, counts, within_sd = _expand_trait_estimates(
            leaf_names,
            fitted_leaf_names,
            means,
            covariance,
            counts,
            within_sd,
        )
    covariance = _validated_sampling_covariance(covariance, means, allow_missing)
    return means, covariance, counts, within_sd, method


def _expand_trait_estimates(
    leaf_names, fitted_leaf_names, means, covariance, counts, within_sd
):
    expanded_means = np.full(len(leaf_names), np.nan, dtype=float)
    covariance = np.asarray(covariance, dtype=float)
    expanded_covariance = np.zeros(
        len(leaf_names) if covariance.ndim == 1 else (len(leaf_names), len(leaf_names)),
        dtype=float,
    )
    expanded_counts = np.zeros(len(leaf_names), dtype=int)
    expanded_within_sd = np.full(len(leaf_names), np.nan, dtype=float)
    selected = [leaf_names.index(name) for name in fitted_leaf_names]
    expanded_means[selected] = means
    if covariance.ndim == 1:
        expanded_covariance[selected] = covariance
    else:
        expanded_covariance[np.ix_(selected, selected)] = covariance
    expanded_counts[selected] = counts
    expanded_within_sd[selected] = within_sd
    return expanded_means, expanded_covariance, expanded_counts, expanded_within_sd


def _validated_sampling_covariance(covariance, means, allow_missing):
    covariance = np.asarray(covariance, dtype=float)
    if covariance.ndim == 1:
        if np.any(covariance < 0.0):
            raise ValueError(
                "Replicate estimation produced an invalid covariance diagonal."
            )
        symmetric = covariance.copy()
        minimum_eigenvalue = float(covariance.min())
    else:
        symmetric = (covariance + covariance.T) / 2.0
        minimum_eigenvalue = float(np.linalg.eigvalsh(symmetric).min())
    invalid_means = np.isinf(means).any() or (
        not allow_missing and not np.isfinite(means).all()
    )
    if invalid_means or not np.isfinite(symmetric).all():
        raise ValueError("Replicate estimation produced non-finite values.")
    tolerance = np.finfo(float).eps * max(1.0, float(np.max(np.abs(symmetric))))
    if minimum_eigenvalue < -tolerance:
        raise ValueError("Replicate estimation produced an invalid covariance matrix.")
    symmetric[np.abs(symmetric) < tolerance] = 0.0
    return symmetric


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
        .fillna(0)
        .to_numpy(int)
    )
    covariance = np.asarray(covariance, dtype=float)
    diagonal = covariance if covariance.ndim == 1 else np.diag(covariance)
    standard_errors = np.sqrt(np.maximum(diagonal, 0.0))
    rows = []
    for leaf_index, leaf_name in enumerate(leaf_names):
        observed = counts[leaf_index] > 0
        rows.append(
            {
                "tree_id": tree_id,
                "leaf_name": leaf_name,
                "trait": trait,
                "n_biological": int(counts[leaf_index]),
                "n_technical": int(technical_counts[leaf_index]),
                "mean": float(means[leaf_index]) if observed else "",
                "within_sd": float(within_sd[leaf_index]) if observed else "",
                "standard_error": (
                    float(standard_errors[leaf_index]) if observed else ""
                ),
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
    allow_missing_traits=(),
    tree_id="",
):
    """Estimate tip means and sampling covariance from response replicates."""
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
    allow_missing_traits = _allowed_missing_trait_set(traits, allow_missing_traits)
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
            dataframe,
            leaf_names,
            traits,
            list(se_columns),
            list(n_columns),
            tree_id,
            allow_missing_traits,
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
        allow_missing_traits,
    )
    values_by_trait = {}
    covariance_by_trait = {}
    models = {}
    rows = []
    for trait in traits:
        means, covariance, counts, within_sd, method = _estimate_one_trait(
            dataframe,
            leaf_names,
            trait,
            batch,
            within_variance,
            allow_missing=trait in allow_missing_traits,
        )
        covariance_frame = (
            covariance
            if np.asarray(covariance).ndim == 1
            else pd.DataFrame(covariance, index=leaf_names, columns=leaf_names)
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
