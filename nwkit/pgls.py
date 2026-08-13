import math
import os
import sys
import warnings
from io import StringIO
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import minimize
from scipy.stats import chi2
from scipy.stats import t as student_t

from nwkit.evolution import evolution_model_spec, validate_evolution_parameter
from nwkit.gaussian import (
    DiagonalLowRankCovariance,
    draw_from_factor,
    factor_diagonal_low_rank_updates,
    factor_logdet,
    is_diagonal,
    materialize_covariance,
    solve_factor,
)
from nwkit.measurement_error import (
    _predictor_error_variance_diagonal,
    fit_conditional_eiv_gaussian,
    fit_latent_predictor,
)
from nwkit.model_matrix import PredictorTerm
from nwkit.sparse_laplace import (
    ContinuousPredictorUncertainty,
    JointPredictorUncertainty,
)
from nwkit.util import (
    normalized_missing_path_key,
    read_input_text,
    validate_distinct_output_paths,
)


def _minimize_variance_components(*args, **kwargs):
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="invalid value encountered in subtract",
            category=RuntimeWarning,
            module=r"scipy\.optimize\._numdiff",
        )
        return minimize(*args, **kwargs)


RESPONSE_REQUIRED_COLUMNS = {
    "tree_id",
    "gene_clade_id",
    "lineage_clade_id",
    "event_type",
    "eligible",
    "coverage_status",
    "species_event_id",
    "species_event_taxa",
    "species_numerator_event_id",
    "species_denominator_event_id",
    "trait",
    "evolution_model",
    "evolution_parameter_name",
    "evolution_parameter",
    "branch_length_mode",
    "raw_contrast",
    "contrast_variance",
}

MAX_DENSE_GAUSSIAN_OBSERVATIONS = 2000

PREDICTOR_REQUIRED_COLUMNS = {
    "tree_id",
    "branch_clade_id",
    "descendant_taxa",
    "numerator_clade_id",
    "denominator_clade_id",
    "trait",
    "evolution_model",
    "evolution_parameter_name",
    "evolution_parameter",
    "branch_length_mode",
    "raw_contrast",
}

RESULT_COLUMNS = [
    "model_id",
    "tree_id",
    "predictor_tree_id",
    "response",
    "response_type",
    "response_family",
    "response_level",
    "response_reference",
    "link_function",
    "response_dispersion",
    "zero_probability",
    "coefficient_penalty",
    "coefficient_prior_sd",
    "separation_warning",
    "term",
    "source_term",
    "predictor_type",
    "predictor_level",
    "predictor_reference",
    "factor_coding",
    "term_test",
    "coefficient",
    "standard_error",
    "statistic",
    "degrees_of_freedom",
    "p_value",
    "confidence_level",
    "confidence_interval_lower",
    "confidence_interval_upper",
    "n_gene_contrasts",
    "n_species_events",
    "n_repeated_gene_contrasts",
    "n_lineages",
    "n_excluded_ineligible",
    "n_excluded_coverage",
    "num_parameters",
    "matrix_rank",
    "condition_number",
    "weighted_residual_sum_squares",
    "residual_scale",
    "r_squared_uncentered",
    "intercept",
    "event_weighting",
    "covariance_estimator",
    "contrast_transform",
    "response_evolution_model",
    "response_evolution_parameter_name",
    "response_evolution_parameter",
    "response_evolution_parameter_status",
    "response_evolution_optimizer_converged",
    "response_evolution_optimizer_message",
    "response_evolution_boundary_warning",
    "response_evolution_parameter_bootstrap_refit",
    "response_branch_length_mode",
    "predictor_evolution_model",
    "predictor_evolution_parameter_name",
    "predictor_evolution_parameter",
    "predictor_evolution_parameter_status",
    "predictor_evolution_optimizer_converged",
    "predictor_evolution_optimizer_message",
    "predictor_evolution_boundary_warning",
    "predictor_evolution_log_likelihood",
    "predictor_branch_length_mode",
    "coverage_policy",
    "small_sample_warning",
    "inference_status",
    "model",
    "inference_method",
    "reml",
    "evolutionary_rate",
    "species_event_variance",
    "lineage_slope_variance",
    "mean_sampling_variance",
    "sampling_variance_fraction",
    "mean_predictor_sampling_variance",
    "mean_latent_predictor_variance",
    "predictor_uncertainty_fraction",
    "predictor_evolutionary_rate",
    "predictor_rate_optimizer_converged",
    "predictor_rate_optimizer_message",
    "predictor_rate_boundary_warning",
    "measurement_error_model",
    "log_likelihood",
    "optimizer_converged",
    "boundary_warning",
    "event_random_effect",
    "lineage_random_slope",
    "ensemble_size",
    "between_tree_variance",
    "tree_support_fraction",
]

SAMPLING_COVARIANCE_REQUIRED_COLUMNS = {
    "tree_id",
    "trait",
    "contrast_id_1",
    "contrast_id_2",
    "sampling_covariance",
}

RANDOM_EFFECT_COLUMNS = [
    "model_id",
    "tree_id",
    "response",
    "effect_type",
    "group_id",
    "term",
    "conditional_mode",
    "conditional_standard_error",
    "conditional_interval_lower",
    "conditional_interval_upper",
    "total_coefficient",
    "total_standard_error",
    "total_interval_lower",
    "total_interval_upper",
    "reliability",
    "variance_component",
    "n_observations",
    "inference_status",
]

SENSITIVITY_COLUMNS = [
    "analysis_type",
    "model_id",
    "tree_id",
    "response",
    "group_id",
    "group_label",
    "term",
    "source_term",
    "n_omitted_gene_contrasts",
    "n_omitted_species_events",
    "n_retained_gene_contrasts",
    "n_retained_species_events",
    "full_coefficient",
    "omitted_coefficient",
    "coefficient_change",
    "relative_change",
    "sign_changed",
    "inference_status",
    "message",
]


def _parse_names(value, option_name):
    if value in (None, ""):
        raise ValueError("'{}' is required.".format(option_name))
    names = [item.strip() for item in str(value).split(",")]
    if any(name == "" for name in names):
        raise ValueError("'{}' contains an empty name.".format(option_name))
    if len(names) != len(set(names)):
        raise ValueError("'{}' contains duplicated names.".format(option_name))
    return names


def _read_tsv(path, option_name):
    text = read_input_text(path)
    if text.strip() == "":
        raise ValueError("'{}' is empty.".format(option_name))
    header = text.splitlines()[0].split("\t")
    duplicated = sorted({name for name in header if header.count(name) > 1})
    if duplicated:
        raise ValueError(
            "'{}' contains duplicated column name(s): {}.".format(
                option_name, ", ".join(duplicated)
            )
        )
    try:
        dataframe = pd.read_csv(
            StringIO(text),
            sep="\t",
            dtype=str,
            keep_default_na=False,
            na_filter=False,
        )
    except Exception as exc:
        raise ValueError("'{}' is not a valid TSV table.".format(option_name)) from exc
    if dataframe.empty:
        raise ValueError("'{}' contains no data rows.".format(option_name))
    return dataframe


def _require_columns(dataframe, required, option_name):
    missing = sorted(required - set(dataframe.columns))
    if missing:
        raise ValueError(
            "'{}' is missing required column(s): {}.".format(
                option_name, ", ".join(missing)
            )
        )


def _numeric_column(dataframe, column, option_name, *, positive=False):
    try:
        values = pd.to_numeric(dataframe[column], errors="raise").to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "'{}' column '{}' must contain numeric values.".format(option_name, column)
        ) from exc
    if not np.isfinite(values).all():
        raise ValueError(
            "'{}' column '{}' must contain finite values.".format(option_name, column)
        )
    if positive and np.any(values <= 0.0):
        raise ValueError(
            "'{}' column '{}' must contain positive values.".format(option_name, column)
        )
    return values


def _selected_traits(dataframe, requested, selector_option, table_option):
    available = set(dataframe["trait"])
    missing = [trait for trait in requested if trait not in available]
    if missing:
        raise ValueError(
            "{} trait(s) are absent from '{}': {}.".format(
                selector_option, table_option, ", ".join(missing)
            )
        )


EVOLUTION_METADATA_COLUMNS = [
    "evolution_model",
    "evolution_parameter_name",
    "evolution_parameter",
    "branch_length_mode",
]


def _prepare_evolution_metadata(dataframe, option_name, group_columns):
    """Validate and normalize the transform recorded for every contrast group."""
    normalized = dataframe.copy()
    for column in ["evolution_model", "evolution_parameter_name", "branch_length_mode"]:
        normalized[column] = normalized[column].fillna("").astype(str)
    normalized_parameters = []
    for row_index, row in normalized.iterrows():
        model = row["evolution_model"]
        try:
            spec = evolution_model_spec(model)
        except ValueError as exc:
            raise ValueError(
                "'{}' contains unsupported evolution_model '{}' at row {}.".format(
                    option_name, model, row_index + 2
                )
            ) from exc
        if not spec.contrast_supported:
            raise ValueError(
                "'{}' evolution_model '{}' cannot produce local contrasts.".format(
                    option_name, model
                )
            )
        expected_name = spec.parameter_name or ""
        if row["evolution_parameter_name"] != expected_name:
            raise ValueError(
                "'{}' evolution_parameter_name is inconsistent with evolution_model "
                "'{}' at row {}.".format(option_name, model, row_index + 2)
            )
        raw_parameter = row["evolution_parameter"]
        parameter_missing = (
            raw_parameter is None
            or (isinstance(raw_parameter, str) and raw_parameter == "")
            or bool(pd.isna(raw_parameter))
        )
        if spec.parameter_name is None:
            if not parameter_missing:
                raise ValueError(
                    "'{}' evolution_model '{}' must have an empty "
                    "evolution_parameter at row {}.".format(
                        option_name, model, row_index + 2
                    )
                )
            parameter: float | str = ""
        else:
            if parameter_missing:
                raise ValueError(
                    "'{}' evolution_model '{}' requires evolution_parameter at row "
                    "{}.".format(option_name, model, row_index + 2)
                )
            try:
                parameter = float(raw_parameter)
                validate_evolution_parameter(model, parameter)
            except (TypeError, ValueError, OverflowError) as exc:
                raise ValueError(
                    "'{}' contains an invalid evolution_parameter for model '{}' at "
                    "row {}: {}".format(option_name, model, row_index + 2, exc)
                ) from exc
        expected_modes = (
            {"original", "unit"} if spec.branch_lengths_used else {"not-applicable"}
        )
        if row["branch_length_mode"] not in expected_modes:
            raise ValueError(
                "'{}' branch_length_mode is inconsistent with evolution_model '{}' "
                "at row {}.".format(option_name, model, row_index + 2)
            )
        normalized_parameters.append(parameter)
    normalized["evolution_parameter"] = normalized_parameters
    for key, rows in normalized.groupby(group_columns, sort=False, dropna=False):
        metadata = rows[EVOLUTION_METADATA_COLUMNS].drop_duplicates()
        if len(metadata) != 1:
            key_values = key if isinstance(key, tuple) else (key,)
            description = ", ".join(
                "{}={}".format(column, value)
                for column, value in zip(group_columns, key_values, strict=True)
            )
            raise ValueError(
                "'{}' mixes evolutionary transforms within contrast group {}.".format(
                    option_name, description
                )
            )
    return normalized


def _evolution_metadata_by_group(dataframe, group_columns):
    metadata = {}
    for key, rows in dataframe.groupby(group_columns, sort=False, dropna=False):
        normalized_key = key if isinstance(key, tuple) else (key,)
        metadata[tuple(str(value) for value in normalized_key)] = {
            column: rows.iloc[0][column] for column in EVOLUTION_METADATA_COLUMNS
        }
    return metadata


def _prepare_responses(dataframe, responses, coverage_policy):
    _require_columns(dataframe, RESPONSE_REQUIRED_COLUMNS, "--infile")
    _selected_traits(dataframe, responses, "--responses", "--infile")
    selected = dataframe[dataframe["trait"].isin(responses)].copy()
    for column in [
        "tree_id",
        "gene_clade_id",
        "lineage_clade_id",
        "event_type",
        "eligible",
        "coverage_status",
        "species_event_id",
        "species_event_taxa",
        "species_numerator_event_id",
        "species_denominator_event_id",
        "trait",
    ]:
        selected[column] = selected[column].fillna("").astype(str)
    invalid_event = selected[~selected["event_type"].isin(["speciation"])]
    invalid_eligible = selected[~selected["eligible"].isin(["yes", "no"])]
    invalid_coverage = selected[
        ~selected["coverage_status"].isin(["complete", "partial", "not-applicable"])
    ]
    if not invalid_event.empty:
        raise ValueError(
            "'--infile' selected response rows must contain only speciation events."
        )
    if not invalid_eligible.empty:
        raise ValueError("'--infile' contains invalid eligible values.")
    if not invalid_coverage.empty:
        raise ValueError("'--infile' contains invalid coverage_status values.")
    inconsistent_coverage = selected[
        (
            (selected["eligible"] == "yes")
            & (selected["coverage_status"] == "not-applicable")
        )
        | (
            (selected["eligible"] == "no")
            & (selected["coverage_status"] != "not-applicable")
        )
    ]
    if not inconsistent_coverage.empty:
        raise ValueError(
            "'--infile' eligible and coverage_status values are inconsistent."
        )
    if (selected["gene_clade_id"] == "").any():
        raise ValueError("'--infile' selected rows require gene_clade_id.")
    if (selected["tree_id"] == "").any():
        raise ValueError(
            "'--infile' selected rows require a non-empty tree_id so gene "
            "families cannot be pooled accidentally."
        )
    eligible_rows = selected[selected["eligible"] == "yes"]
    if (eligible_rows["lineage_clade_id"] == "").any():
        raise ValueError("'--infile' eligible rows require lineage_clade_id.")
    if (eligible_rows["species_event_id"] == "").any():
        raise ValueError("'--infile' eligible rows require species_event_id.")
    orientation_columns = [
        "species_numerator_event_id",
        "species_denominator_event_id",
    ]
    if eligible_rows[orientation_columns].eq("").any(axis=None):
        raise ValueError("'--infile' eligible rows require species-event orientation.")
    if (
        eligible_rows["species_numerator_event_id"]
        == eligible_rows["species_denominator_event_id"]
    ).any():
        raise ValueError(
            "'--infile' species numerator and denominator event IDs must differ."
        )
    duplicated = selected.duplicated(
        subset=["tree_id", "trait", "gene_clade_id"], keep=False
    )
    if duplicated.any():
        raise ValueError(
            "'--infile' contains duplicated tree_id/trait/gene_clade_id rows."
        )
    selected["raw_contrast"] = _numeric_column(selected, "raw_contrast", "--infile")
    selected["contrast_variance"] = _numeric_column(
        selected, "contrast_variance", "--infile", positive=True
    )
    selected = _prepare_evolution_metadata(
        selected,
        "--infile",
        ["tree_id", "trait"],
    )
    selected["_eligible_for_model"] = selected["eligible"] == "yes"
    selected["_coverage_for_model"] = (
        True if coverage_policy == "any" else selected["coverage_status"] == "complete"
    )
    return selected


def _prepare_predictors(dataframe, predictors):
    _require_columns(dataframe, PREDICTOR_REQUIRED_COLUMNS, "--predictor-contrasts")
    _selected_traits(dataframe, predictors, "--predictors", "--predictor-contrasts")
    selected = dataframe[dataframe["trait"].isin(predictors)].copy()
    for column in [
        "tree_id",
        "branch_clade_id",
        "descendant_taxa",
        "numerator_clade_id",
        "denominator_clade_id",
        "trait",
    ]:
        selected[column] = selected[column].fillna("").astype(str)
    if (selected["tree_id"] == "").any():
        raise ValueError(
            "'--predictor-contrasts' selected rows require a non-empty tree_id."
        )
    if selected["tree_id"].nunique(dropna=False) != 1:
        raise ValueError(
            "'--predictor-contrasts' selected rows must come from exactly one tree_id."
        )
    if (selected["branch_clade_id"] == "").any():
        raise ValueError(
            "'--predictor-contrasts' selected rows require branch_clade_id."
        )
    if selected[["numerator_clade_id", "denominator_clade_id"]].eq("").any(axis=None):
        raise ValueError(
            "'--predictor-contrasts' selected rows require contrast orientation."
        )
    if (selected["numerator_clade_id"] == selected["denominator_clade_id"]).any():
        raise ValueError(
            "'--predictor-contrasts' numerator and denominator clade IDs must differ."
        )
    duplicated = selected.duplicated(subset=["trait", "branch_clade_id"], keep=False)
    if duplicated.any():
        raise ValueError(
            "'--predictor-contrasts' contains duplicated trait/branch_clade_id rows."
        )
    selected["raw_contrast"] = _numeric_column(
        selected, "raw_contrast", "--predictor-contrasts"
    )
    selected = _prepare_evolution_metadata(
        selected,
        "--predictor-contrasts",
        ["tree_id", "trait"],
    )
    return selected


def _predictor_by_event(dataframe, predictors):
    by_event: dict[str, dict[str, Any]] = dict()
    for record in dataframe.to_dict("records"):
        event_id = record["branch_clade_id"]
        event = by_event.setdefault(
            event_id,
            {
                "descendant_taxa": record["descendant_taxa"],
                "numerator_clade_id": record["numerator_clade_id"],
                "denominator_clade_id": record["denominator_clade_id"],
                "values": {},
            },
        )
        metadata = (
            record["descendant_taxa"],
            record["numerator_clade_id"],
            record["denominator_clade_id"],
        )
        expected = (
            event["descendant_taxa"],
            event["numerator_clade_id"],
            event["denominator_clade_id"],
        )
        if metadata != expected:
            raise ValueError(
                "'--predictor-contrasts' has inconsistent metadata for species event '{}'.".format(
                    event_id
                )
            )
        event["values"][record["trait"]] = float(record["raw_contrast"])
    complete = {
        event_id: event
        for event_id, event in by_event.items()
        if all(predictor in event["values"] for predictor in predictors)
    }
    return complete


def _validate_and_attach_predictors(
    responses,
    predictor_events,
    predictors,
    predictor_posteriors=None,
):
    missing_events = sorted(set(responses["species_event_id"]) - set(predictor_events))
    if missing_events:
        raise ValueError(
            "'--predictor-contrasts' lacks selected predictor rows for species event(s): {}.".format(
                ", ".join(missing_events)
            )
        )
    predictor_values = []
    for record in responses.to_dict("records"):
        event = predictor_events[record["species_event_id"]]
        if record["species_event_taxa"] != event["descendant_taxa"]:
            raise ValueError(
                "Species-event taxa disagree for species_event_id '{}'.".format(
                    record["species_event_id"]
                )
            )
        if (
            record["species_numerator_event_id"] != event["numerator_clade_id"]
            or record["species_denominator_event_id"] != event["denominator_clade_id"]
        ):
            raise ValueError(
                "Contrast orientation disagrees for species_event_id '{}'.".format(
                    record["species_event_id"]
                )
            )
        predictor_values.append(
            [
                (
                    predictor_posteriors[name]["posterior"].mean[
                        predictor_posteriors[name]["event_index"][
                            record["species_event_id"]
                        ]
                    ]
                    if predictor_posteriors and name in predictor_posteriors
                    else event["values"][name]
                )
                for name in predictors
            ]
        )
    attached = responses.copy()
    attached["_predictor_values"] = predictor_values
    return attached


def _prepare_sampling_covariance(
    dataframe,
    *,
    option_name="--response-sampling-covariance",
):
    _require_columns(
        dataframe,
        SAMPLING_COVARIANCE_REQUIRED_COLUMNS,
        option_name,
    )
    prepared = dataframe.copy()
    for column in ["tree_id", "trait", "contrast_id_1", "contrast_id_2"]:
        prepared[column] = prepared[column].fillna("").astype(str)
        if prepared[column].eq("").any():
            raise ValueError(
                "'{}' column '{}' cannot be empty.".format(option_name, column)
            )
    prepared["sampling_covariance"] = _numeric_column(
        prepared,
        "sampling_covariance",
        option_name,
    )
    if "covariance_representation" not in prepared.columns:
        prepared["covariance_representation"] = "covariance"
    prepared["covariance_representation"] = (
        prepared["covariance_representation"].fillna("covariance").astype(str)
    )
    supported_representations = {"covariance", "factor-loading"}
    unsupported = sorted(
        set(prepared["covariance_representation"]) - supported_representations
    )
    if unsupported:
        raise ValueError(
            "'{}' contains unsupported covariance representation(s): {}.".format(
                option_name, ", ".join(unsupported)
            )
        )
    keys = []
    for record in prepared.to_dict("records"):
        first = str(record["contrast_id_1"])
        second = str(record["contrast_id_2"])
        keys.append(
            (
                str(record["tree_id"]),
                str(record["trait"]),
                min(first, second),
                max(first, second),
            )
            if record["covariance_representation"] == "covariance"
            else (
                str(record["tree_id"]),
                str(record["trait"]),
                "factor:{}".format(first),
                second,
            )
        )
    if len(keys) != len(set(keys)):
        raise ValueError(
            "'{}' contains duplicated covariance pairs.".format(option_name)
        )
    prepared["_pair_key"] = keys
    return prepared


def _validate_covariance_matrix(matrix, label):
    matrix = np.asarray(matrix, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("{} must be a square matrix.".format(label))
    if not np.isfinite(matrix).all():
        raise ValueError("{} must contain finite values.".format(label))
    symmetric = (matrix + matrix.T) / 2.0
    asymmetry = float(np.max(np.abs(matrix - matrix.T))) if matrix.size else 0.0
    scale = max(1.0, float(np.max(np.abs(symmetric))) if matrix.size else 1.0)
    tolerance = np.finfo(float).eps * scale * max(1, len(matrix)) * 100.0
    if asymmetry > tolerance:
        raise ValueError("{} must be symmetric.".format(label))
    eigenvalues = (
        np.diag(symmetric) if is_diagonal(symmetric) else np.linalg.eigvalsh(symmetric)
    )
    if eigenvalues.size and float(eigenvalues.min()) < -tolerance:
        raise ValueError("{} must be positive semidefinite.".format(label))
    symmetric[np.abs(symmetric) < tolerance] = 0.0
    return symmetric


def _validate_factor_covariance(covariance, size, label):
    diagonal = np.asarray(covariance.diagonal, dtype=float)
    factor = covariance.low_rank
    factor_finite = (
        np.isfinite(factor.data).all()
        if sparse.issparse(factor)
        else np.isfinite(np.asarray(factor)).all()
    )
    if (
        diagonal.shape != (size,)
        or np.any(diagonal < 0.0)
        or not np.isfinite(diagonal).all()
        or factor.shape[0] != size
        or not factor_finite
    ):
        raise ValueError("{} factor representation is malformed.".format(label))
    return covariance


def _validate_covariance_representation(covariance, size, label):
    """Validate either a diagonal vector or a dense covariance matrix."""
    if isinstance(covariance, DiagonalLowRankCovariance):
        return _validate_factor_covariance(covariance, size, label)
    covariance = np.asarray(covariance, dtype=float)
    if covariance.ndim == 1:
        if covariance.shape != (size,) or not np.isfinite(covariance).all():
            raise ValueError(
                "{} diagonal has invalid dimensions or values.".format(label)
            )
        if np.any(covariance < 0.0):
            raise ValueError("{} diagonal must be non-negative.".format(label))
        return covariance
    covariance = _validate_covariance_matrix(covariance, label)
    if covariance.shape != (size, size):
        raise ValueError("{} has the wrong dimensions.".format(label))
    return covariance


def _covariance_diagonal(covariance, factor=None):
    if covariance is None:
        if factor is None:
            raise ValueError("A factor is required for an implicit covariance.")
        if sparse.issparse(factor):
            return np.asarray(factor.multiply(factor).sum(axis=1)).reshape(-1)
        return np.sum(np.square(np.asarray(factor, dtype=float)), axis=1)
    if isinstance(covariance, DiagonalLowRankCovariance):
        update = covariance.low_rank
        update_diagonal = (
            np.asarray(update.multiply(update).sum(axis=1)).reshape(-1)
            if sparse.issparse(update)
            else np.sum(np.square(update), axis=1)
        )
        return np.asarray(covariance.diagonal, dtype=float) + update_diagonal
    covariance = np.asarray(covariance, dtype=float)
    return covariance if covariance.ndim == 1 else np.diag(covariance)


def _sampling_factor_for_model(selected, contrast_ids, label, model_label):
    row_index = {contrast_id: index for index, contrast_id in enumerate(contrast_ids)}
    missing = sorted(set(contrast_ids) - set(selected["contrast_id_1"].astype(str)))
    if missing:
        raise ValueError(
            "{} factor is incomplete for model '{}' (missing: {}).".format(
                label, model_label, ", ".join(missing[:10])
            )
        )
    latent_ids = list(dict.fromkeys(selected["contrast_id_2"].astype(str)))
    latent_index = {latent_id: index for index, latent_id in enumerate(latent_ids)}
    rows = []
    columns = []
    values = []
    for record in selected.to_dict("records"):
        contrast_id = str(record["contrast_id_1"])
        if contrast_id not in row_index:
            continue
        rows.append(row_index[contrast_id])
        columns.append(latent_index[str(record["contrast_id_2"])])
        values.append(float(record["sampling_covariance"]))
    factor = sparse.csr_matrix(
        (values, (rows, columns)), shape=(len(contrast_ids), len(latent_ids))
    )
    return DiagonalLowRankCovariance(np.zeros(len(contrast_ids)), factor)


def _sampling_explicit_covariance_for_model(
    selected, tree_id, response, contrast_ids, label, model_label
):
    pair_values = {
        record["_pair_key"]: float(record["sampling_covariance"])
        for record in selected.to_dict("records")
    }
    size = len(contrast_ids)
    matrix = np.zeros((size, size), dtype=float)
    missing = []
    for index1, first in enumerate(contrast_ids):
        for index2 in range(index1, size):
            second = contrast_ids[index2]
            key = (tree_id, response, min(first, second), max(first, second))
            if key not in pair_values:
                missing.append("{}|{}".format(first, second))
                continue
            value = pair_values[key]
            matrix[index1, index2] = value
            matrix[index2, index1] = value
    if missing:
        raise ValueError(
            "{} is incomplete for model '{}' (missing: {}).".format(
                label, model_label, ", ".join(missing[:10])
            )
        )
    return _validate_covariance_matrix(matrix, label)


def _sampling_matrix_for_model(
    dataframe,
    tree_id,
    response,
    contrast_ids,
    response_rows,
    *,
    option_name="--response-sampling-covariance",
    label="Response sampling covariance",
    model_label=None,
):
    if dataframe is None:
        replicate_columns = {"sampling_variance", "replicate_model"}.intersection(
            response_rows.columns
        )
        if replicate_columns:
            raise ValueError(
                "Replicate-aware contrasts require the full "
                "'{}' table; diagonal variances alone "
                "cannot represent covariance introduced by the PIC transform.".format(
                    option_name
                )
            )
        return np.zeros(len(response_rows), dtype=float)
    selected = dataframe[
        (dataframe["tree_id"] == tree_id) & (dataframe["trait"] == response)
    ]
    resolved_model_label = model_label or _model_id(tree_id, response)
    if selected.empty:
        raise ValueError(
            "{} has no rows for model '{}'.".format(label, resolved_model_label)
        )
    representations = set(selected["covariance_representation"])
    if representations == {"factor-loading"}:
        return _sampling_factor_for_model(
            selected, contrast_ids, label, resolved_model_label
        )
    if representations != {"covariance"}:
        raise ValueError(
            "{} mixes covariance representations within one model.".format(label)
        )
    return _sampling_explicit_covariance_for_model(
        selected,
        tree_id,
        response,
        contrast_ids,
        label,
        resolved_model_label,
    )


def _prepare_predictor_posteriors(
    predictor_contrasts,
    predictors,
    sampling_covariance,
):
    if sampling_covariance is None:
        return {}
    if "contrast_variance" not in predictor_contrasts.columns:
        raise ValueError(
            "Predictor sampling covariance requires 'contrast_variance' in "
            "'--predictor-contrasts'."
        )
    prepared = predictor_contrasts.copy()
    prepared["contrast_variance"] = _numeric_column(
        prepared,
        "contrast_variance",
        "--predictor-contrasts",
        positive=True,
    )
    posteriors = {}
    available_traits = set(sampling_covariance["trait"])
    for predictor in predictors:
        if predictor not in available_traits:
            continue
        rows = prepared[prepared["trait"] == predictor].copy()
        rows = rows.sort_values("branch_clade_id", kind="stable")
        tree_id = str(rows.iloc[0]["tree_id"])
        event_ids = rows["branch_clade_id"].astype(str).tolist()
        covariance = _sampling_matrix_for_model(
            sampling_covariance,
            tree_id,
            predictor,
            event_ids,
            rows,
            option_name="--predictor-sampling-covariance",
            label="Predictor sampling covariance",
            model_label=predictor,
        )
        conditioning_covariance = (
            materialize_covariance(covariance)
            if isinstance(covariance, DiagonalLowRankCovariance)
            else covariance
        )
        posterior = fit_latent_predictor(
            rows["raw_contrast"].to_numpy(dtype=float),
            np.diag(rows["contrast_variance"].to_numpy(dtype=float)),
            conditioning_covariance,
            include_intercept=False,
        )
        prior_covariance = posterior.evolutionary_rate * np.diag(
            rows["contrast_variance"].to_numpy(dtype=float)
        )
        posteriors[predictor] = {
            "event_ids": event_ids,
            "event_index": {
                event_id: index for index, event_id in enumerate(event_ids)
            },
            "posterior": posterior,
            "posterior_factor": _covariance_factor(posterior.covariance),
            "sampling_covariance": conditioning_covariance,
            "prior_covariance": prior_covariance,
        }
    return posteriors


def _covariance_factor(covariance):
    symmetric = (np.asarray(covariance, dtype=float) + np.asarray(covariance).T) / 2.0
    diagonal = np.diag(symmetric)
    if np.array_equal(symmetric, np.diag(diagonal)):
        if np.any(diagonal < 0.0):
            raise ValueError(
                "Predictor posterior covariance is not positive semidefinite."
            )
        return np.sqrt(diagonal)
    eigenvalues, eigenvectors = np.linalg.eigh(symmetric)
    tolerance = (
        np.finfo(float).eps
        * max(1.0, float(np.max(np.abs(symmetric))))
        * max(1, len(symmetric))
        * 100.0
    )
    if float(eigenvalues.min()) < -tolerance:
        raise ValueError("Predictor posterior covariance is not positive semidefinite.")
    positive = eigenvalues > tolerance
    return eigenvectors[:, positive] * np.sqrt(eigenvalues[positive])


def _predictor_uncertainties_for_rows(rows, predictors, posteriors):
    if not posteriors:
        return {}
    event_ids = rows["species_event_id"].astype(str).tolist()
    uncertainties = {}
    for predictor in predictors:
        if predictor not in posteriors:
            continue
        state = posteriors[predictor]
        indices = [state["event_index"][event_id] for event_id in event_ids]
        uncertainties[predictor] = ContinuousPredictorUncertainty(
            factor=state["posterior_factor"],
            observation_index=np.asarray(indices, dtype=int),
        )
    return uncertainties


def _grouped_predictor_uncertainties_for_rows(rows, states):
    event_ids = rows["species_event_id"].astype(str).tolist()
    selected = []
    for state in states or ():
        indices = [state["event_index"][event_id] for event_id in event_ids]
        uncertainty = state["uncertainty"]
        selected.append(
            {
                **state,
                "uncertainty": JointPredictorUncertainty(
                    factors=tuple(factor[indices] for factor in uncertainty.factors)
                ),
            }
        )
    return selected


def _model_id(tree_id, response):
    return "{}:{}".format(tree_id, response) if tree_id else response


def _predictor_term_fields(metadata):
    return {
        "source_term": metadata.source,
        "predictor_type": metadata.kind,
        "predictor_level": metadata.level,
        "predictor_reference": metadata.reference,
        "factor_coding": metadata.coding,
        "term_test": "coefficient",
    }


def _append_predictor_omnibus_rows(
    rows, beta, covariance, predictors, metadata, groups
):
    index_by_term = {term: index for index, term in enumerate(predictors)}
    for source, terms in groups.items():
        if not terms or metadata[terms[0]].kind not in {"categorical", "ordered"}:
            continue
        indices = [index_by_term[term] for term in terms]
        coefficients = np.asarray(beta, dtype=float)[indices]
        selected_covariance = np.asarray(covariance, dtype=float)[
            np.ix_(indices, indices)
        ]
        try:
            statistic = float(
                coefficients @ np.linalg.solve(selected_covariance, coefficients)
            )
        except np.linalg.LinAlgError:
            statistic = float(
                coefficients @ np.linalg.pinv(selected_covariance) @ coefficients
            )
        template = next(row for row in rows if row["term"] == terms[0]).copy()
        template.update(
            {
                "term": source,
                "source_term": source,
                "predictor_level": "",
                "term_test": "omnibus",
                "coefficient": "",
                "standard_error": "",
                "statistic": statistic,
                "degrees_of_freedom": len(indices),
                "p_value": float(chi2.sf(statistic, len(indices))),
                "confidence_interval_lower": "",
                "confidence_interval_upper": "",
                "inference_status": "ok",
            }
        )
        rows.append(template)


def _fit_model(
    dataframe,
    predictors,
    *,
    confidence_level,
    event_weighting,
    coverage_policy,
    excluded_ineligible,
    excluded_coverage,
    predictor_metadata,
    predictor_groups,
):
    tree_id = str(dataframe.iloc[0]["tree_id"])
    response = str(dataframe.iloc[0]["trait"])
    y_raw = dataframe["raw_contrast"].to_numpy(dtype=float)
    variances = dataframe["contrast_variance"].to_numpy(dtype=float)
    predictor_raw = np.asarray(dataframe["_predictor_values"].tolist(), dtype=float)
    scale = np.sqrt(variances)
    y = y_raw / scale
    design = predictor_raw / scale[:, None]
    event_ids = dataframe["species_event_id"].to_numpy(dtype=str)
    unique_events, inverse, counts = np.unique(
        event_ids, return_inverse=True, return_counts=True
    )
    n_observations = len(dataframe)
    n_events = len(unique_events)
    num_parameters = len(predictors)
    if n_events <= num_parameters:
        raise ValueError(
            "Model '{}' needs more unique species events than predictors (events={}; predictors={}).".format(
                _model_id(tree_id, response), n_events, num_parameters
            )
        )
    weights = np.ones(n_observations, dtype=float)
    if event_weighting == "equal":
        weights = 1.0 / counts[inverse]
    sqrt_weights = np.sqrt(weights)
    weighted_design = design * sqrt_weights[:, None]
    weighted_y = y * sqrt_weights
    matrix_rank = int(np.linalg.matrix_rank(weighted_design))
    if matrix_rank != num_parameters:
        raise ValueError(
            "Model '{}' predictor matrix is rank deficient (rank={}; predictors={}).".format(
                _model_id(tree_id, response), matrix_rank, num_parameters
            )
        )
    gram = weighted_design.T @ weighted_design
    condition_number = float(np.linalg.cond(gram))
    if (
        not math.isfinite(condition_number)
        or condition_number > 1.0 / np.finfo(float).eps
    ):
        raise ValueError(
            "Model '{}' predictor matrix is numerically singular.".format(
                _model_id(tree_id, response)
            )
        )
    beta = np.linalg.solve(gram, weighted_design.T @ weighted_y)
    residual = y - design @ beta
    weighted_residual = sqrt_weights * residual
    weighted_rss = float(weighted_residual @ weighted_residual)
    weighted_total = float(weighted_y @ weighted_y)
    r_squared = (
        float("nan") if weighted_total == 0.0 else 1.0 - weighted_rss / weighted_total
    )
    bread = np.linalg.inv(gram)
    meat = np.zeros((num_parameters, num_parameters), dtype=float)
    for event_index in range(n_events):
        selected = inverse == event_index
        score = design[selected].T @ (weights[selected] * residual[selected])
        meat += np.outer(score, score)
    degrees_of_freedom = n_events - num_parameters
    covariance = (n_events / degrees_of_freedom) * (bread @ meat @ bread)
    covariance = (covariance + covariance.T) / 2.0
    diagonal = np.diag(covariance)
    tolerance = np.finfo(float).eps * max(1.0, float(np.max(np.abs(diagonal))))
    if np.any(diagonal < -tolerance):
        raise ValueError(
            "Model '{}' produced an invalid covariance estimate.".format(
                _model_id(tree_id, response)
            )
        )
    standard_errors = np.sqrt(np.maximum(diagonal, 0.0))
    critical = float(student_t.ppf(0.5 + confidence_level / 2.0, degrees_of_freedom))
    residual_scale = weighted_rss / degrees_of_freedom
    lineages = set(dataframe["lineage_clade_id"])
    rows = []
    for index, predictor in enumerate(predictors):
        coefficient = float(beta[index])
        standard_error = float(standard_errors[index])
        if standard_error == 0.0:
            statistic: float | str = ""
            p_value: float | str = ""
            inference_status = "zero-cluster-variance"
        else:
            statistic_value = coefficient / standard_error
            statistic = statistic_value
            p_value = float(
                2.0 * student_t.sf(abs(statistic_value), degrees_of_freedom)
            )
            inference_status = "ok"
        rows.append(
            {
                "model_id": _model_id(tree_id, response),
                "tree_id": tree_id,
                "response": response,
                "response_type": "continuous",
                "response_family": "gaussian",
                "response_level": "",
                "response_reference": "",
                "link_function": "identity",
                "term": predictor,
                **_predictor_term_fields(predictor_metadata[predictor]),
                "coefficient": coefficient,
                "standard_error": standard_error,
                "statistic": statistic,
                "degrees_of_freedom": degrees_of_freedom,
                "p_value": p_value,
                "confidence_level": confidence_level,
                "confidence_interval_lower": coefficient - critical * standard_error,
                "confidence_interval_upper": coefficient + critical * standard_error,
                "n_gene_contrasts": n_observations,
                "n_species_events": n_events,
                "n_repeated_gene_contrasts": n_observations - n_events,
                "n_lineages": len(lineages),
                "n_excluded_ineligible": excluded_ineligible,
                "n_excluded_coverage": excluded_coverage,
                "num_parameters": num_parameters,
                "matrix_rank": matrix_rank,
                "condition_number": condition_number,
                "weighted_residual_sum_squares": weighted_rss,
                "residual_scale": residual_scale,
                "r_squared_uncentered": r_squared,
                "intercept": "no",
                "event_weighting": event_weighting,
                "covariance_estimator": "species-event-cluster-HC1",
                "contrast_transform": "gene-contrast-variance",
                "coverage_policy": coverage_policy,
                "small_sample_warning": "yes" if n_events < 20 else "no",
                "inference_status": inference_status,
            }
        )
    _append_predictor_omnibus_rows(
        rows, beta, covariance, predictors, predictor_metadata, predictor_groups
    )
    return rows


def _solve_positive_definite(cholesky, values):
    return solve_factor(cholesky, values)


def _prepare_component_matrices(n_observations, components, component_factors):
    component_matrices = []
    component_scales = []
    for name, matrix in components:
        factor = component_factors.get(name)
        if matrix is None:
            if factor is None:
                raise ValueError(
                    "Implicit variance component '{}' requires a factor.".format(name)
                )
            prepared_matrix = None
        else:
            prepared_matrix = np.asarray(matrix, dtype=float)
            valid_shape = prepared_matrix.shape in {
                (n_observations,),
                (n_observations, n_observations),
            }
            if not valid_shape:
                raise ValueError(
                    "Variance component '{}' has the wrong dimensions.".format(name)
                )
            if not np.isfinite(prepared_matrix).all():
                raise ValueError(
                    "Variance component '{}' contains non-finite values.".format(name)
                )
        diagonal = _covariance_diagonal(prepared_matrix, factor)
        positive = diagonal[diagonal > 0.0]
        if not len(positive):
            raise ValueError("Variance component '{}' has zero diagonal.".format(name))
        component_matrices.append((name, prepared_matrix))
        component_scales.append(float(np.mean(positive)))
    return component_matrices, component_scales


def _normalize_component_factors(
    n_observations, component_matrices, component_scales, component_factors
):
    unknown_factor_names = set(component_factors) - {
        name for name, _ in component_matrices
    }
    if unknown_factor_names:
        raise ValueError("Low-rank factors reference unknown variance components.")
    normalized_factors = {}
    for (name, matrix), scale in zip(component_matrices, component_scales, strict=True):
        if name not in component_factors:
            continue
        factor_value = component_factors[name]
        factor = (
            sparse.csr_matrix(factor_value, dtype=float)
            if sparse.issparse(factor_value)
            else np.asarray(factor_value, dtype=float)
        )
        finite = (
            np.isfinite(factor.data).all()
            if sparse.issparse(factor)
            else np.isfinite(factor).all()
        )
        row_squares = (
            np.asarray(factor.multiply(factor).sum(axis=1)).reshape(-1)
            if sparse.issparse(factor)
            else np.sum(np.square(factor), axis=1)
        )
        valid_factor = (
            factor.ndim == 2
            and factor.shape[0] == n_observations
            and finite
            and np.allclose(
                row_squares,
                _covariance_diagonal(matrix, factor),
            )
        )
        if not valid_factor:
            raise ValueError(
                "Low-rank factor for variance component '{}' is inconsistent.".format(
                    name
                )
            )
        normalized_factors[name] = factor / math.sqrt(scale)
    return normalized_factors


def _profile_component_representation(
    fixed_covariance, component_matrices, component_scales, normalized_factors
):
    fixed_is_structured = isinstance(fixed_covariance, DiagonalLowRankCovariance) or (
        fixed_covariance.ndim == 1 or is_diagonal(fixed_covariance)
    )
    structured_model = fixed_is_structured and all(
        matrix is None
        or np.asarray(matrix).ndim == 1
        or is_diagonal(matrix)
        or name in normalized_factors
        for name, matrix in component_matrices
    )
    if structured_model:
        if isinstance(fixed_covariance, DiagonalLowRankCovariance):
            working_fixed_covariance = np.asarray(
                fixed_covariance.diagonal, dtype=float
            ).copy()
            fixed_updates = [fixed_covariance.low_rank]
        else:
            working_fixed_covariance = (
                fixed_covariance.copy()
                if fixed_covariance.ndim == 1
                else np.diag(fixed_covariance).copy()
            )
            fixed_updates = []
        normalized_components = [
            (
                name,
                None
                if name in normalized_factors
                else _covariance_diagonal(matrix) / scale,
            )
            for (name, matrix), scale in zip(
                component_matrices, component_scales, strict=True
            )
        ]
    else:
        working_fixed_covariance = materialize_covariance(fixed_covariance)
        normalized_components = [
            (
                name,
                None if matrix is None else materialize_covariance(matrix) / scale,
            )
            for (name, matrix), scale in zip(
                component_matrices, component_scales, strict=True
            )
        ]
        if any(component is None for _, component in normalized_components):
            raise ValueError(
                "Implicit low-rank components require diagonal fixed covariance."
            )
        fixed_updates = []
    return (
        structured_model,
        working_fixed_covariance,
        normalized_components,
        fixed_updates,
    )


def _covariance_is_zero(covariance):
    if isinstance(covariance, DiagonalLowRankCovariance):
        loading = covariance.low_rank
        loading_is_zero = (
            loading.nnz == 0
            if sparse.issparse(loading)
            else not np.any(np.asarray(loading, dtype=float))
        )
        return not np.any(covariance.diagonal) and loading_is_zero
    return not np.any(covariance)


def _prepare_profile_covariance_components(
    n_observations, fixed_covariance, components, component_factors
):
    component_factors = {} if component_factors is None else component_factors
    component_matrices, component_scales = _prepare_component_matrices(
        n_observations, components, component_factors
    )
    fixed_covariance = _validate_covariance_representation(
        fixed_covariance, n_observations, "Fixed sampling covariance"
    )
    normalized_factors = _normalize_component_factors(
        n_observations,
        component_matrices,
        component_scales,
        component_factors,
    )
    (
        structured_model,
        working_fixed_covariance,
        normalized_components,
        fixed_updates,
    ) = _profile_component_representation(
        fixed_covariance,
        component_matrices,
        component_scales,
        normalized_factors,
    )
    return (
        fixed_covariance,
        component_matrices,
        component_scales,
        normalized_factors,
        structured_model,
        working_fixed_covariance,
        normalized_components,
        fixed_updates,
    )


def _requires_dense_profile_covariance(fixed_covariance, components, component_factors):
    if isinstance(fixed_covariance, DiagonalLowRankCovariance):
        fixed_is_dense = False
    else:
        fixed_covariance = np.asarray(fixed_covariance)
        fixed_is_dense = fixed_covariance.ndim == 2 and not is_diagonal(
            fixed_covariance
        )
    if fixed_is_dense:
        return True
    factor_names = set() if component_factors is None else set(component_factors)
    return any(
        matrix is not None
        and np.asarray(matrix).ndim == 2
        and not is_diagonal(np.asarray(matrix))
        and name not in factor_names
        for name, matrix in components
    )


def _profile_covariance_fit(
    y,
    design,
    fixed_covariance,
    components,
    *,
    reml,
    starting_log_variances=None,
    component_factors=None,
    allow_large_dense=False,
):
    n_observations = len(y)
    num_parameters = design.shape[1]
    if n_observations > MAX_DENSE_GAUSSIAN_OBSERVATIONS and (
        _requires_dense_profile_covariance(
            fixed_covariance, components, component_factors
        )
    ):
        message = (
            "Dense Gaussian covariance fitting is limited to {} observations "
            "by default (received {}); use diagonal sampling covariance and "
            "low-rank variance components for larger analyses."
        ).format(MAX_DENSE_GAUSSIAN_OBSERVATIONS, n_observations)
        if not allow_large_dense:
            raise ValueError(
                message + " Pass allow_large_dense=True to attempt the allocation."
            )
        warnings.warn(
            message + " Large dense allocation was explicitly enabled; attempting it.",
            RuntimeWarning,
            stacklevel=2,
        )
    (
        fixed_covariance,
        component_matrices,
        component_scales,
        normalized_factors,
        structured_model,
        working_fixed_covariance,
        normalized_components,
        fixed_updates,
    ) = _prepare_profile_covariance_components(
        n_observations, fixed_covariance, components, component_factors
    )
    ordinary_beta = np.linalg.lstsq(design, y, rcond=None)[0]
    ordinary_residual = y - design @ ordinary_beta
    response_scale = max(
        float(np.mean(y**2)),
        float(np.mean(ordinary_residual**2)),
        float(np.mean(_covariance_diagonal(fixed_covariance))) if len(y) else 0.0,
        np.finfo(float).tiny,
    )
    lower_variance = max(response_scale * 1e-12, np.finfo(float).tiny)
    upper_variance = max(response_scale * 1e6, lower_variance * 1e6)
    bounds = [(math.log(lower_variance), math.log(upper_variance))] * len(
        normalized_components
    )

    def evaluate(log_variances, response=y, return_details=False):
        variances = np.exp(np.asarray(log_variances, dtype=float))
        covariance = working_fixed_covariance.copy()
        if structured_model:
            low_rank_updates = list(fixed_updates)
            for variance, (name, component) in zip(
                variances, normalized_components, strict=True
            ):
                if name in normalized_factors:
                    low_rank_updates.append(
                        math.sqrt(float(variance)) * normalized_factors[name]
                    )
                else:
                    covariance = covariance + variance * component
            if np.any(covariance <= 0.0) or not np.isfinite(covariance).all():
                return float("inf")
            try:
                if low_rank_updates:
                    cholesky = factor_diagonal_low_rank_updates(
                        covariance, low_rank_updates
                    )
                    low_rank = (
                        sparse.hstack(
                            [sparse.csr_matrix(update) for update in low_rank_updates],
                            format="csr",
                        )
                        if any(sparse.issparse(update) for update in low_rank_updates)
                        else np.column_stack(low_rank_updates)
                    )
                    covariance_representation = DiagonalLowRankCovariance(
                        covariance, low_rank
                    )
                else:
                    cholesky = np.sqrt(covariance)
                    covariance_representation = covariance
                inverse_design = _solve_positive_definite(cholesky, design)
                gram = design.T @ inverse_design
                gram_sign, gram_logdet = np.linalg.slogdet(gram)
                if gram_sign <= 0.0:
                    return float("inf")
                beta = np.linalg.solve(
                    gram,
                    design.T @ _solve_positive_definite(cholesky, response),
                )
                residual = response - design @ beta
                inverse_residual = _solve_positive_definite(cholesky, residual)
                quadratic = float(residual @ inverse_residual)
                covariance_logdet = factor_logdet(cholesky)
            except np.linalg.LinAlgError:
                return float("inf")
        else:
            for variance, (_, component) in zip(
                variances, normalized_components, strict=True
            ):
                assert component is not None
                covariance = covariance + variance * component
            covariance = (covariance + covariance.T) / 2.0
            covariance_representation = covariance
            try:
                cholesky = np.linalg.cholesky(covariance)
                inverse_design = _solve_positive_definite(cholesky, design)
                gram = design.T @ inverse_design
                gram_sign, gram_logdet = np.linalg.slogdet(gram)
                if gram_sign <= 0.0:
                    return float("inf")
                beta = np.linalg.solve(
                    gram,
                    design.T @ _solve_positive_definite(cholesky, response),
                )
                residual = response - design @ beta
                inverse_residual = _solve_positive_definite(cholesky, residual)
                quadratic = float(residual @ inverse_residual)
                covariance_logdet = 2.0 * float(np.log(np.diag(cholesky)).sum())
            except np.linalg.LinAlgError:
                return float("inf")
        effective_n = n_observations - num_parameters if reml else n_observations
        objective = 0.5 * (
            effective_n * math.log(2.0 * math.pi)
            + covariance_logdet
            + quadratic
            + (gram_logdet if reml else 0.0)
        )
        if not math.isfinite(objective):
            return float("inf")
        if not return_details:
            return objective
        beta_covariance = np.linalg.inv(gram)
        component_variances = {
            name: float(variance / scale)
            for variance, (name, _), scale in zip(
                variances,
                normalized_components,
                component_scales,
                strict=True,
            )
        }
        return {
            "objective": objective,
            "beta": beta,
            "beta_covariance": beta_covariance,
            "residual": residual,
            "inverse_residual": inverse_residual,
            "covariance": covariance_representation,
            "cholesky": cholesky,
            "quadratic": quadratic,
            "component_variances": component_variances,
            "log_variances": np.asarray(log_variances, dtype=float),
            "lower_variance": lower_variance,
            "upper_variance": upper_variance,
        }

    if len(normalized_components) == 1 and _covariance_is_zero(fixed_covariance):
        unit_fit = evaluate(np.asarray([0.0]), return_details=True)
        if not isinstance(unit_fit, dict):
            raise ValueError(
                "Variance-component closed-form fit produced an invalid covariance."
            )
        effective_n = n_observations - num_parameters if reml else n_observations
        optimum = float(unit_fit["quadratic"]) / effective_n
        optimum = float(np.clip(optimum, lower_variance, upper_variance))
        details = evaluate(np.log(np.asarray([optimum])), return_details=True)
        if not isinstance(details, dict):
            raise ValueError(
                "Variance-component closed-form fit produced an invalid fit."
            )
        details["optimizer_converged"] = True
        details["optimizer_message"] = "closed-form single-scale covariance fit"
        details["reml"] = bool(reml)
        details["boundary_warning"] = bool(
            optimum <= lower_variance * 10.0 or optimum >= upper_variance / 10.0
        )
        return details

    if starting_log_variances is None:
        starts = [
            np.log(
                np.asarray(
                    [response_scale]
                    + [max(response_scale * 0.1, lower_variance)]
                    * (len(normalized_components) - 1)
                )
            ),
            np.log(np.repeat(max(response_scale, lower_variance), len(components))),
        ]
    else:
        starts = [np.asarray(starting_log_variances, dtype=float)]
    candidates = []
    for start in starts:
        result = _minimize_variance_components(
            evaluate,
            np.clip(start, bounds[0][0], bounds[0][1]),
            method="L-BFGS-B",
            bounds=bounds,
        )
        if math.isfinite(float(result.fun)):
            candidates.append(result)
    if not candidates:
        raise ValueError("Variance-component optimization failed to find a finite fit.")
    result = min(candidates, key=lambda candidate: float(candidate.fun))
    if not result.success:
        fallback = _minimize_variance_components(
            evaluate,
            result.x,
            method="Powell",
            bounds=bounds,
            options={"maxiter": 5000},
        )
        if math.isfinite(float(fallback.fun)) and float(fallback.fun) <= float(
            result.fun
        ):
            result = fallback
    if not result.success:
        raise ValueError(
            "Variance-component optimization did not converge: {}".format(
                result.message
            )
        )
    details = evaluate(result.x, return_details=True)
    if not isinstance(details, dict):
        raise ValueError("Variance-component optimization produced an invalid fit.")
    details["optimizer_converged"] = bool(result.success)
    details["optimizer_message"] = str(result.message)
    details["reml"] = bool(reml)
    details["boundary_warning"] = bool(
        np.any(np.exp(details["log_variances"]) <= lower_variance * 10.0)
        or np.any(np.exp(details["log_variances"]) >= upper_variance / 10.0)
    )
    return details


def _random_effect_policy(policy, identifiable, label):
    if policy not in {"auto", "no", "yes"}:
        raise ValueError("Unsupported {} policy: {}.".format(label, policy))
    if policy == "yes" and not identifiable:
        raise ValueError("{} was requested but is not identifiable.".format(label))
    return identifiable if policy == "auto" else policy == "yes"


def _covariance_component_inner(left, right):
    """Frobenius inner product without materializing low-rank covariances."""
    left_kind, left_values = left
    right_kind, right_values = right
    if left_kind == "diagonal" and right_kind == "diagonal":
        return float(np.dot(left_values, right_values))
    if left_kind == "factor" and right_kind == "diagonal":
        left_kind, right_kind = right_kind, left_kind
        left_values, right_values = right_values, left_values
    if left_kind == "diagonal" and right_kind == "factor":
        row_squares = (
            np.asarray(right_values.multiply(right_values).sum(axis=1)).reshape(-1)
            if sparse.issparse(right_values)
            else np.einsum("ij,ij->i", right_values, right_values)
        )
        return float(np.dot(left_values, row_squares))
    if left_kind == "factor" and right_kind == "factor":
        cross = sparse.csr_matrix(left_values).T @ sparse.csr_matrix(right_values)
        return float(np.dot(cross.data, cross.data))
    raise ValueError("Unknown covariance-component representation.")


def _component_adds_rank(components, candidate):
    """Test covariance-component independence using only a small Gram matrix."""
    representations = [*components, candidate]
    count = len(representations)
    gram = np.empty((count, count), dtype=float)
    for row in range(count):
        for column in range(row + 1):
            value = _covariance_component_inner(
                representations[row], representations[column]
            )
            gram[row, column] = value
            gram[column, row] = value
    norms = np.sqrt(np.maximum(np.diag(gram), 0.0))
    if norms[-1] <= np.finfo(float).tiny:
        return False
    normalized = gram / np.outer(norms, norms)
    normalized = np.nan_to_num(normalized, nan=0.0, posinf=0.0, neginf=0.0)
    tolerance = np.finfo(float).eps * max(1, count) * 1000.0
    existing_rank = int(np.linalg.matrix_rank(normalized[:-1, :-1], tol=tolerance))
    return int(np.linalg.matrix_rank(normalized, tol=tolerance)) > existing_rank


def _parametric_bootstrap_coefficients(
    fit,
    design,
    fixed_covariance,
    components,
    *,
    reml,
    replicates,
    seed,
    component_factors=None,
    allow_large_dense=False,
):
    if replicates < 2:
        raise ValueError("Parametric bootstrap requires at least two replicates.")
    rng = np.random.default_rng(seed)
    coefficients: list[np.ndarray] = []
    mean = design @ fit["beta"]
    maximum_attempts = max(replicates * 3, replicates + 10)
    attempts = 0
    while len(coefficients) < replicates and attempts < maximum_attempts:
        attempts += 1
        response = mean + draw_from_factor(
            fit["cholesky"], rng.standard_normal(len(mean)), rng=rng
        )
        try:
            bootstrap_fit = _profile_covariance_fit(
                response,
                design,
                fixed_covariance,
                components,
                reml=reml,
                starting_log_variances=fit["log_variances"],
                component_factors=component_factors,
                allow_large_dense=allow_large_dense,
            )
        except ValueError:
            continue
        if not bootstrap_fit["optimizer_converged"]:
            continue
        coefficients.append(bootstrap_fit["beta"])
    if len(coefficients) < replicates:
        raise ValueError(
            "Parametric bootstrap produced only {} successful fits in {} attempts.".format(
                len(coefficients), attempts
            )
        )
    return np.asarray(coefficients, dtype=float)


def _parametric_bootstrap_eiv_coefficients(
    fit,
    design,
    predictor_uncertainties,
    predictor_columns,
    fixed_covariance,
    components,
    *,
    reml,
    replicates,
    seed,
    component_factors=None,
    allow_large_dense=False,
):
    if replicates < 2:
        raise ValueError("Parametric bootstrap requires at least two replicates.")
    rng = np.random.default_rng(seed)
    coefficients: list[np.ndarray] = []
    mean = design @ fit["beta"]
    maximum_attempts = max(replicates * 3, replicates + 10)
    attempts = 0
    starting = np.concatenate([fit["beta"], fit["log_variances"]])
    while len(coefficients) < replicates and attempts < maximum_attempts:
        attempts += 1
        response = mean + draw_from_factor(
            fit["cholesky"], rng.standard_normal(len(mean)), rng=rng
        )
        try:
            bootstrap_fit = fit_conditional_eiv_gaussian(
                response,
                design,
                predictor_uncertainties,
                predictor_columns,
                fixed_covariance,
                components,
                reml=False,
                starting_parameters=starting,
                component_factors=component_factors,
                allow_large_dense=allow_large_dense,
            )
        except ValueError:
            continue
        if not bootstrap_fit["optimizer_converged"]:
            continue
        coefficients.append(bootstrap_fit["beta"])
    if len(coefficients) < replicates:
        raise ValueError(
            "Parametric bootstrap produced only {} successful errors-in-variables "
            "fits in {} attempts.".format(len(coefficients), attempts)
        )
    return np.asarray(coefficients, dtype=float)


def _fit_profile_or_eiv(
    response,
    design,
    fixed_covariance,
    components,
    predictor_uncertainties,
    predictor_columns,
    *,
    reml,
    component_factors,
    allow_large_dense,
    starting_fit=None,
):
    if predictor_uncertainties:
        starting_parameters = None
        if starting_fit is not None:
            starting_parameters = np.concatenate(
                [starting_fit["beta"], starting_fit["log_variances"]]
            )
        return fit_conditional_eiv_gaussian(
            response,
            design,
            predictor_uncertainties,
            predictor_columns,
            fixed_covariance,
            components,
            reml=False,
            starting_parameters=starting_parameters,
            component_factors=component_factors,
            allow_large_dense=allow_large_dense,
        )
    return _profile_covariance_fit(
        response,
        design,
        fixed_covariance,
        components,
        reml=reml,
        starting_log_variances=(
            None if starting_fit is None else starting_fit["log_variances"]
        ),
        component_factors=component_factors,
        allow_large_dense=allow_large_dense,
    )


def _without_lineage_component(components, component_factors):
    retained = [
        component
        for component in components
        if component[0] != "lineage_slope_variance"
    ]
    factors = {
        name: factor
        for name, factor in component_factors.items()
        if name != "lineage_slope_variance"
    }
    return retained, factors


def _likelihood_ratio(null_fit, full_fit):
    return max(
        0.0,
        2.0 * (float(null_fit["objective"]) - float(full_fit["objective"])),
    )


def _empirical_tail_probability(statistics, observed_statistic):
    return float(
        (1 + np.sum(np.asarray(statistics) >= observed_statistic))
        / (len(statistics) + 1)
    )


def _lineage_test_significance(term_test, statistic):
    if term_test == "lineage-heterogeneity":
        return (
            1.0 if statistic <= 0.0 else float(0.5 * chi2.sf(statistic, 1)),
            "boundary-mixture-chi-square",
        )
    return "", "parametric-bootstrap-required-for-joint-null"


def _normal_critical_value(confidence_level):
    return float(student_t.ppf(0.5 + confidence_level / 2.0, math.inf))


def _total_lineage_slope(fixed_coefficient, random_mode):
    return float(fixed_coefficient) + float(random_mode)


def _lineage_reliability(conditional_variance, lineage_variance):
    if lineage_variance <= 0.0:
        return 0.0
    return max(0.0, min(1.0, 1.0 - conditional_variance / lineage_variance))


def _subset_predictor_uncertainties(uncertainties, columns, retained_indices):
    retained_set = set(retained_indices)
    new_index = {old: new for new, old in enumerate(retained_indices)}
    selected_uncertainties = []
    selected_columns: list[int | tuple[int, ...]] = []
    for uncertainty, uncertainty_columns in zip(uncertainties, columns, strict=True):
        old_columns = (
            (uncertainty_columns,)
            if isinstance(uncertainty_columns, int)
            else tuple(uncertainty_columns)
        )
        selected = [
            (position, column)
            for position, column in enumerate(old_columns)
            if column in retained_set
        ]
        if not selected:
            continue
        if len(selected) != len(old_columns):
            if not isinstance(uncertainty, JointPredictorUncertainty):
                raise RuntimeError("Cannot partially subset predictor uncertainty.")
            uncertainty = JointPredictorUncertainty(
                factors=tuple(uncertainty.factors[position] for position, _ in selected)
            )
        remapped = tuple(new_index[column] for _, column in selected)
        selected_uncertainties.append(uncertainty)
        selected_columns.append(remapped[0] if len(remapped) == 1 else remapped)
    return selected_uncertainties, selected_columns


def _lineage_test_models(
    design,
    predictors,
    predictor_groups,
    components,
    component_factors,
    random_designs,
    predictor_uncertainties,
    predictor_columns,
):
    null_components, null_factors = _without_lineage_component(
        components, component_factors
    )
    models = [
        {
            "term_test": "lineage-heterogeneity",
            "term": "(lineage-slope heterogeneity)",
            "source_term": "(all predictors)",
            "design": design,
            "components": null_components,
            "component_factors": null_factors,
            "predictor_uncertainties": predictor_uncertainties,
            "predictor_columns": predictor_columns,
            "degrees_of_freedom": 1,
        }
    ]
    index_by_predictor = {
        predictor: index for index, predictor in enumerate(predictors)
    }
    lineage_design_by_predictor = dict(random_designs.get("lineage", ()))
    for source, source_terms in predictor_groups.items():
        omitted = {index_by_predictor[term] for term in source_terms}
        retained = [index for index in range(len(predictors)) if index not in omitted]
        selected_uncertainties, selected_columns = _subset_predictor_uncertainties(
            predictor_uncertainties, predictor_columns, retained
        )
        retained_lineage = [
            lineage_design_by_predictor[predictors[index]]
            for index in retained
            if predictors[index] in lineage_design_by_predictor
        ]
        joint_components = list(null_components)
        joint_factors = dict(null_factors)
        if retained_lineage:
            joint_components.append(("lineage_slope_variance", None))
            joint_factors["lineage_slope_variance"] = sparse.hstack(
                retained_lineage, format="csr"
            )
        models.append(
            {
                "term_test": "average-and-lineage-joint",
                "term": source,
                "source_term": source,
                "design": design[:, retained],
                "components": joint_components,
                "component_factors": joint_factors,
                "predictor_uncertainties": selected_uncertainties,
                "predictor_columns": selected_columns,
                "degrees_of_freedom": len(source_terms) + 1,
            }
        )
    return models


def _parametric_bootstrap_likelihood_ratio(
    null_fit,
    null_model,
    full_model,
    observed_statistic,
    *,
    replicates,
    seed,
    fixed_covariance,
    allow_large_dense,
):
    rng = np.random.default_rng(seed)
    mean = null_model["design"] @ null_fit["beta"]
    statistics: list[float] = []
    attempts = 0
    maximum_attempts = max(replicates * 3, replicates + 10)
    while len(statistics) < replicates and attempts < maximum_attempts:
        attempts += 1
        response = mean + draw_from_factor(
            null_fit["cholesky"], rng.standard_normal(len(mean)), rng=rng
        )
        try:
            bootstrap_null = _fit_profile_or_eiv(
                response,
                null_model["design"],
                fixed_covariance,
                null_model["components"],
                null_model["predictor_uncertainties"],
                null_model["predictor_columns"],
                reml=False,
                component_factors=null_model["component_factors"],
                allow_large_dense=allow_large_dense,
                starting_fit=null_fit,
            )
            bootstrap_full = _fit_profile_or_eiv(
                response,
                full_model["design"],
                fixed_covariance,
                full_model["components"],
                full_model["predictor_uncertainties"],
                full_model["predictor_columns"],
                reml=False,
                component_factors=full_model["component_factors"],
                allow_large_dense=allow_large_dense,
                starting_fit=full_model["fit"],
            )
        except ValueError:
            continue
        statistic = _likelihood_ratio(bootstrap_null, bootstrap_full)
        statistics.append(statistic)
    if len(statistics) < replicates:
        raise ValueError(
            "Lineage parametric bootstrap produced only {} successful fits in "
            "{} attempts.".format(len(statistics), attempts)
        )
    return _empirical_tail_probability(statistics, observed_statistic)


def _build_covariance_components(
    design,
    evolutionary_variances,
    sampling_covariance,
    event_inverse,
    event_counts,
    lineage_inverse,
    lineage_counts,
    unique_events,
    unique_lineages,
    predictors,
    *,
    n_events,
    num_parameters,
    event_weighting,
    model,
    event_random_effect,
    lineage_random_slope,
):
    n_observations = len(design)
    balance = np.ones(n_observations, dtype=float)
    if event_weighting == "equal":
        balance = np.sqrt(event_counts[event_inverse].astype(float))
    fixed_covariance = (
        DiagonalLowRankCovariance(
            np.asarray(sampling_covariance.diagonal, dtype=float).copy(),
            sampling_covariance.low_rank.copy(),
        )
        if isinstance(sampling_covariance, DiagonalLowRankCovariance)
        else sampling_covariance.copy()
    )
    if not np.all(balance == 1.0):
        if isinstance(fixed_covariance, DiagonalLowRankCovariance):
            fixed_covariance = DiagonalLowRankCovariance(
                fixed_covariance.diagonal * np.square(balance),
                fixed_covariance.low_rank.multiply(balance[:, None])
                if sparse.issparse(fixed_covariance.low_rank)
                else fixed_covariance.low_rank * balance[:, None],
            )
        elif fixed_covariance.ndim == 1:
            fixed_covariance *= np.square(balance)
        else:
            fixed_covariance *= np.outer(balance, balance)
    raw_evolutionary_component = evolutionary_variances.copy()
    evolutionary_component = evolutionary_variances * np.square(balance)
    components = [("evolutionary_rate", evolutionary_component)]
    component_factors = {}
    raw_components = {"evolutionary_rate": raw_evolutionary_component}
    random_designs = {}
    use_event = False
    use_lineage = False
    if model == "hierarchical":
        covariance_representations = [("diagonal", evolutionary_component)]
        event_may_be_identifiable = (
            np.any(event_counts > 1)
            and int(np.sum(event_counts > 1)) >= 2
            and n_events > num_parameters
        )
        if event_may_be_identifiable:
            base_event_design = sparse.csr_matrix(
                (
                    np.ones(n_observations),
                    (np.arange(n_observations), event_inverse),
                ),
                shape=(n_observations, len(unique_events)),
            )
            event_design = base_event_design.multiply(balance[:, None]).tocsr()
            event_identifiable = _component_adds_rank(
                covariance_representations, ("factor", event_design)
            )
        else:
            event_identifiable = False
        use_event = _random_effect_policy(
            event_random_effect, event_identifiable, "species-event random effect"
        )
        if use_event:
            covariance_representations.append(("factor", event_design))
        repeated_lineages = lineage_counts >= 2
        lineage_may_be_identifiable = (
            int(np.sum(repeated_lineages)) >= 2
            and n_events > num_parameters
            and np.any(np.abs(design) > 0.0)
        )
        lineage_designs = []
        lineage_factor = sparse.csr_matrix((n_observations, 0), dtype=float)
        if lineage_may_be_identifiable:
            base_lineage_design = sparse.csr_matrix(
                (
                    balance,
                    (np.arange(n_observations), lineage_inverse),
                ),
                shape=(n_observations, len(unique_lineages)),
            )
            for predictor_index, predictor in enumerate(predictors):
                random_design = base_lineage_design.multiply(
                    design[:, predictor_index, None]
                ).tocsr()
                lineage_designs.append((predictor, random_design))
            lineage_factor = sparse.hstack(
                [random_design for _, random_design in lineage_designs], format="csr"
            )
        lineage_identifiable = lineage_may_be_identifiable and _component_adds_rank(
            covariance_representations, ("factor", lineage_factor)
        )
        use_lineage = _random_effect_policy(
            lineage_random_slope,
            lineage_identifiable,
            "lineage random slope",
        )
    if use_event:
        components.append(("species_event_variance", None))
        component_factors["species_event_variance"] = event_design
        raw_components["species_event_variance"] = np.asarray(
            base_event_design.multiply(base_event_design).sum(axis=1)
        ).reshape(-1)
        random_designs["species_event"] = event_design
    if use_lineage:
        components.append(("lineage_slope_variance", None))
        component_factors["lineage_slope_variance"] = lineage_factor
        unbalanced_lineage = lineage_factor.multiply((1.0 / balance)[:, None])
        raw_components["lineage_slope_variance"] = np.asarray(
            unbalanced_lineage.multiply(unbalanced_lineage).sum(axis=1)
        ).reshape(-1)
        random_designs["lineage"] = lineage_designs
    return (
        fixed_covariance,
        components,
        component_factors,
        raw_components,
        random_designs,
        use_event,
        use_lineage,
    )


def _predictor_measurement_result_fields(
    predictor,
    predictor_uncertainties,
    predictor_posteriors,
    grouped_uncertainties_by_term,
    mean_response_sampling_variance,
):
    if predictor in grouped_uncertainties_by_term:
        state, term_index = grouped_uncertainties_by_term[predictor]
        factor = state["uncertainty"].factors[term_index]
        mean_variance = float(
            np.mean(np.asarray(factor.multiply(factor).sum(axis=1)).reshape(-1))
        )
        return {
            "mean_predictor_sampling_variance": mean_variance,
            "mean_latent_predictor_variance": mean_variance,
            "predictor_uncertainty_fraction": "",
            "predictor_evolutionary_rate": "",
            "predictor_rate_optimizer_converged": "not-applicable",
            "predictor_rate_optimizer_message": "categorical-state uncertainty",
            "predictor_rate_boundary_warning": "not-applicable",
            "measurement_error_model": "latent-predictor",
        }
    if predictor not in predictor_posteriors:
        return {
            "mean_predictor_sampling_variance": "",
            "mean_latent_predictor_variance": "",
            "predictor_uncertainty_fraction": "",
            "predictor_evolutionary_rate": "",
            "predictor_rate_optimizer_converged": "",
            "predictor_rate_optimizer_message": "",
            "predictor_rate_boundary_warning": "",
            "measurement_error_model": (
                "response-only" if mean_response_sampling_variance > 0.0 else "none"
            ),
        }
    state = predictor_posteriors[predictor]
    posterior = state["posterior"]
    uncertainty = predictor_uncertainties[predictor]
    if isinstance(uncertainty, ContinuousPredictorUncertainty):
        selected_factor = uncertainty.factor[uncertainty.observation_index]
        if np.ndim(selected_factor) == 1:
            mean_posterior_variance = float(np.mean(np.square(selected_factor)))
        else:
            mean_posterior_variance = float(
                np.mean(np.sum(np.square(selected_factor), axis=1))
            )
    else:
        mean_posterior_variance = float(np.mean(np.diag(uncertainty)))
    mean_prior_variance = float(np.mean(np.diag(state["prior_covariance"])))
    return {
        "mean_predictor_sampling_variance": float(
            np.mean(np.diag(state["sampling_covariance"]))
        ),
        "mean_latent_predictor_variance": mean_posterior_variance,
        "predictor_uncertainty_fraction": (
            mean_posterior_variance / mean_prior_variance
        ),
        "predictor_evolutionary_rate": posterior.evolutionary_rate,
        "predictor_rate_optimizer_converged": (
            "yes" if posterior.optimizer_converged else "no"
        ),
        "predictor_rate_optimizer_message": posterior.optimizer_message,
        "predictor_rate_boundary_warning": (
            "yes" if posterior.boundary_warning else "no"
        ),
        "measurement_error_model": "latent-predictor",
    }


def _weighted_predictor_uncertainties(
    predictor_uncertainties,
    grouped_predictor_uncertainties,
    index_by_predictor,
    balance,
):
    raw: list[Any] = []
    weighted: list[Any] = []
    columns: list[int | tuple[int, ...]] = []
    grouped_by_term = {}
    outer_balance = None
    for predictor, uncertainty in predictor_uncertainties.items():
        if isinstance(uncertainty, ContinuousPredictorUncertainty):
            raw.append(uncertainty)
            weighted.append(
                ContinuousPredictorUncertainty(
                    factor=uncertainty.factor,
                    observation_index=uncertainty.observation_index,
                    row_scale=balance,
                )
            )
        else:
            if outer_balance is None:
                outer_balance = np.outer(balance, balance)
            covariance = np.asarray(uncertainty, dtype=float)
            raw.append(covariance)
            weighted.append(covariance * outer_balance)
        columns.append(index_by_predictor[predictor])
    for state in grouped_predictor_uncertainties:
        uncertainty = state["uncertainty"]
        raw.append(uncertainty)
        weighted.append(
            JointPredictorUncertainty(
                factors=uncertainty.factors,
                row_scale=balance,
            )
        )
        term_names = tuple(state["term_names"])
        columns.append(tuple(index_by_predictor[term] for term in term_names))
        for term_index, term in enumerate(term_names):
            grouped_by_term[term] = (state, term_index)
    return raw, weighted, columns, grouped_by_term


def _append_lineage_inference_rows(
    rows,
    response_values,
    design,
    predictors,
    predictor_groups,
    fixed_covariance,
    components,
    component_factors,
    random_designs,
    predictor_uncertainties,
    predictor_columns,
    *,
    lineage_inference,
    bootstrap_replicates,
    seed,
    allow_large_dense,
):
    if lineage_inference == "none" or "lineage" not in random_designs:
        return
    full_model = {
        "design": design,
        "components": components,
        "component_factors": component_factors,
        "predictor_uncertainties": predictor_uncertainties,
        "predictor_columns": predictor_columns,
    }
    full_fit = _fit_profile_or_eiv(
        response_values,
        design,
        fixed_covariance,
        components,
        predictor_uncertainties,
        predictor_columns,
        reml=False,
        component_factors=component_factors,
        allow_large_dense=allow_large_dense,
    )
    full_model["fit"] = full_fit
    test_models = _lineage_test_models(
        design,
        predictors,
        predictor_groups,
        components,
        component_factors,
        random_designs,
        predictor_uncertainties,
        predictor_columns,
    )
    for test_index, null_model in enumerate(test_models):
        null_fit = _fit_profile_or_eiv(
            response_values,
            null_model["design"],
            fixed_covariance,
            null_model["components"],
            null_model["predictor_uncertainties"],
            null_model["predictor_columns"],
            reml=False,
            component_factors=null_model["component_factors"],
            allow_large_dense=allow_large_dense,
        )
        statistic = _likelihood_ratio(null_fit, full_fit)
        if lineage_inference == "parametric-bootstrap":
            p_value: float | str = _parametric_bootstrap_likelihood_ratio(
                null_fit,
                null_model,
                full_model,
                statistic,
                replicates=bootstrap_replicates,
                seed=seed + test_index,
                fixed_covariance=fixed_covariance,
                allow_large_dense=allow_large_dense,
            )
            status = "ok"
        else:
            p_value, status = _lineage_test_significance(
                null_model["term_test"], statistic
            )
        template = rows[0].copy()
        full_variances = full_fit["component_variances"]
        template.update(
            {
                "term": null_model["term"],
                "source_term": null_model["source_term"],
                "predictor_type": "",
                "predictor_level": "",
                "predictor_reference": "",
                "factor_coding": "",
                "term_test": null_model["term_test"],
                "coefficient": "",
                "standard_error": "",
                "statistic": statistic,
                "degrees_of_freedom": null_model["degrees_of_freedom"],
                "p_value": p_value,
                "confidence_interval_lower": "",
                "confidence_interval_upper": "",
                "inference_status": status,
                "inference_method": (
                    "parametric-bootstrap-likelihood-ratio"
                    if lineage_inference == "parametric-bootstrap"
                    else "likelihood-ratio"
                ),
                "reml": "no",
                "covariance_estimator": "gaussian{}-ML".format(
                    "-eiv" if predictor_uncertainties else ""
                ),
                "weighted_residual_sum_squares": full_fit["quadratic"],
                "residual_scale": full_variances["evolutionary_rate"],
                "evolutionary_rate": full_variances["evolutionary_rate"],
                "species_event_variance": full_variances.get(
                    "species_event_variance", 0.0
                ),
                "lineage_slope_variance": full_variances.get(
                    "lineage_slope_variance", 0.0
                ),
                "log_likelihood": -float(full_fit["objective"]),
                "optimizer_converged": (
                    "yes" if full_fit["optimizer_converged"] else "no"
                ),
                "boundary_warning": ("yes" if full_fit["boundary_warning"] else "no"),
            }
        )
        rows.append(template)


def _lineage_random_effect_rows(
    fit,
    design,
    random_designs,
    predictors,
    unique_lineages,
    lineage_counts,
    *,
    model_id,
    tree_id,
    response,
    variance,
    confidence_level,
):
    rows = []
    if variance <= 0.0:
        critical = float("nan")
    else:
        critical = _normal_critical_value(confidence_level)
    predictor_index = {predictor: index for index, predictor in enumerate(predictors)}
    inverse_design = _solve_positive_definite(fit["cholesky"], design)
    for predictor, random_design in random_designs["lineage"]:
        fixed_index = predictor_index[predictor]
        for start in range(0, random_design.shape[1], 64):
            stop = min(start + 64, random_design.shape[1])
            design_batch = np.asarray(
                random_design[:, start:stop].toarray(), dtype=float
            )
            inverse_batch = _solve_positive_definite(fit["cholesky"], design_batch)
            modes = variance * (design_batch.T @ fit["inverse_residual"])
            conditional_variance = variance - np.square(variance) * np.einsum(
                "ij,ij->j", design_batch, inverse_batch
            )
            fixed_random_information = inverse_design.T @ design_batch
            cross = fit["beta_covariance"] @ fixed_random_information
            integrated_variance = conditional_variance + np.square(
                variance
            ) * np.einsum("ij,ij->j", fixed_random_information, cross)
            covariance_fixed_random = -variance * cross[fixed_index]
            total_variance = (
                float(fit["beta_covariance"][fixed_index, fixed_index])
                + integrated_variance
                + 2.0 * covariance_fixed_random
            )
            for local_index, lineage_index in enumerate(range(start, stop)):
                mode = float(modes[local_index])
                conditional_se = math.sqrt(
                    max(float(integrated_variance[local_index]), 0.0)
                )
                total = _total_lineage_slope(fit["beta"][fixed_index], mode)
                total_se = math.sqrt(max(float(total_variance[local_index]), 0.0))
                reliability = _lineage_reliability(
                    float(conditional_variance[local_index]), variance
                )
                rows.append(
                    {
                        "model_id": model_id,
                        "tree_id": tree_id,
                        "response": response,
                        "effect_type": "lineage_slope",
                        "group_id": unique_lineages[lineage_index],
                        "term": predictor,
                        "conditional_mode": mode,
                        "conditional_standard_error": conditional_se,
                        "conditional_interval_lower": mode - critical * conditional_se,
                        "conditional_interval_upper": mode + critical * conditional_se,
                        "total_coefficient": total,
                        "total_standard_error": total_se,
                        "total_interval_lower": total - critical * total_se,
                        "total_interval_upper": total + critical * total_se,
                        "reliability": reliability,
                        "variance_component": variance,
                        "n_observations": int(lineage_counts[lineage_index]),
                        "inference_status": "empirical-bayes-conditional-on-variance",
                    }
                )
    return rows


def _fit_covariance_model(
    dataframe,
    predictors,
    sampling_covariance,
    predictor_uncertainties,
    predictor_posteriors,
    grouped_predictor_uncertainties,
    *,
    confidence_level,
    event_weighting,
    coverage_policy,
    excluded_ineligible,
    excluded_coverage,
    model,
    inference,
    bootstrap_replicates,
    seed,
    reml,
    event_random_effect,
    lineage_random_slope,
    lineage_inference,
    return_fit_state,
    predictor_metadata,
    predictor_groups,
    allow_large_dense,
):
    tree_id = str(dataframe.iloc[0]["tree_id"])
    response = str(dataframe.iloc[0]["trait"])
    model_id = _model_id(tree_id, response)
    response_values = dataframe["raw_contrast"].to_numpy(dtype=float)
    design = np.asarray(dataframe["_predictor_values"].tolist(), dtype=float)
    evolutionary_variances = dataframe["contrast_variance"].to_numpy(dtype=float)
    event_ids = dataframe["species_event_id"].to_numpy(dtype=str)
    lineage_ids = dataframe["lineage_clade_id"].to_numpy(dtype=str)
    unique_events, event_inverse, event_counts = np.unique(
        event_ids, return_inverse=True, return_counts=True
    )
    unique_lineages, lineage_inverse, lineage_counts = np.unique(
        lineage_ids, return_inverse=True, return_counts=True
    )
    n_observations = len(dataframe)
    n_events = len(unique_events)
    num_parameters = len(predictors)
    if n_events <= num_parameters:
        raise ValueError(
            "Model '{}' needs more unique species events than predictors "
            "(events={}; predictors={}).".format(model_id, n_events, num_parameters)
        )
    if int(np.linalg.matrix_rank(design)) != num_parameters:
        raise ValueError(
            "Model '{}' predictor matrix is rank deficient.".format(model_id)
        )
    sampling_covariance = _validate_covariance_representation(
        sampling_covariance, n_observations, "Response sampling covariance"
    )
    (
        fixed_covariance,
        components,
        component_factors,
        raw_components,
        random_designs,
        use_event,
        use_lineage,
    ) = _build_covariance_components(
        design,
        evolutionary_variances,
        sampling_covariance,
        event_inverse,
        event_counts,
        lineage_inverse,
        lineage_counts,
        unique_events,
        unique_lineages,
        predictors,
        n_events=n_events,
        num_parameters=num_parameters,
        event_weighting=event_weighting,
        model=model,
        event_random_effect=event_random_effect,
        lineage_random_slope=lineage_random_slope,
    )
    if lineage_inference != "none" and not use_lineage:
        raise ValueError(
            "Lineage inference was requested but lineage random slopes are not "
            "identifiable for model '{}'.".format(model_id)
        )
    balance = np.ones(n_observations, dtype=float)
    if event_weighting == "equal":
        balance = np.sqrt(event_counts[event_inverse].astype(float))
    index_by_predictor = {
        predictor: index for index, predictor in enumerate(predictors)
    }
    (
        raw_predictor_uncertainties,
        balanced_predictor_uncertainties,
        predictor_uncertainty_columns,
        grouped_uncertainties_by_term,
    ) = _weighted_predictor_uncertainties(
        predictor_uncertainties,
        grouped_predictor_uncertainties,
        index_by_predictor,
        balance,
    )
    if balanced_predictor_uncertainties:
        eiv_components = components
        eiv_fixed_covariance = fixed_covariance
        fit = fit_conditional_eiv_gaussian(
            response_values,
            design,
            balanced_predictor_uncertainties,
            predictor_uncertainty_columns,
            eiv_fixed_covariance,
            eiv_components,
            reml=False,
            component_factors=component_factors,
            allow_large_dense=allow_large_dense,
        )
    else:
        fit = _profile_covariance_fit(
            response_values,
            design,
            fixed_covariance,
            components,
            reml=reml,
            component_factors=component_factors,
            allow_large_dense=allow_large_dense,
        )
    effective_reml = bool(fit.get("reml", reml))
    beta = fit["beta"]
    beta_covariance = fit["beta_covariance"]
    bootstrap_coefficients = None
    if inference == "parametric-bootstrap":
        if balanced_predictor_uncertainties:
            bootstrap_coefficients = _parametric_bootstrap_eiv_coefficients(
                fit,
                design,
                balanced_predictor_uncertainties,
                predictor_uncertainty_columns,
                eiv_fixed_covariance,
                eiv_components,
                reml=effective_reml,
                replicates=bootstrap_replicates,
                seed=seed,
                component_factors=component_factors,
                allow_large_dense=allow_large_dense,
            )
        else:
            bootstrap_coefficients = _parametric_bootstrap_coefficients(
                fit,
                design,
                fixed_covariance,
                components,
                reml=effective_reml,
                replicates=bootstrap_replicates,
                seed=seed,
                component_factors=component_factors,
                allow_large_dense=allow_large_dense,
            )
        standard_errors = np.std(bootstrap_coefficients, axis=0, ddof=1)
    elif inference == "wald":
        standard_errors = np.sqrt(np.maximum(np.diag(beta_covariance), 0.0))
    else:
        raise ValueError("Unsupported inference method: {}.".format(inference))
    degrees_of_freedom = n_events - num_parameters
    critical = float(student_t.ppf(0.5 + confidence_level / 2.0, degrees_of_freedom))
    inverse_y = _solve_positive_definite(fit["cholesky"], response_values)
    total_quadratic = float(response_values @ inverse_y)
    r_squared = (
        float("nan")
        if total_quadratic == 0.0
        else 1.0 - fit["quadratic"] / total_quadratic
    )
    component_variances = fit["component_variances"]
    evolutionary_rate = component_variances["evolutionary_rate"]
    event_variance = component_variances.get("species_event_variance", 0.0)
    lineage_variance = component_variances.get("lineage_slope_variance", 0.0)
    mean_sampling_variance = float(np.mean(_covariance_diagonal(sampling_covariance)))
    fitted_variance = mean_sampling_variance
    for uncertainty, columns in zip(
        raw_predictor_uncertainties,
        predictor_uncertainty_columns,
        strict=True,
    ):
        fitted_variance += float(
            np.mean(
                _predictor_error_variance_diagonal(
                    beta, uncertainty, columns, n_observations
                )
            )
        )
    for name, variance in component_variances.items():
        fitted_variance += variance * float(
            np.mean(_covariance_diagonal(raw_components[name]))
        )
    sampling_fraction = (
        0.0 if fitted_variance == 0.0 else mean_sampling_variance / fitted_variance
    )
    condition_number = float(np.linalg.cond(beta_covariance))
    rows = []
    for index, predictor in enumerate(predictors):
        coefficient = float(beta[index])
        standard_error = float(standard_errors[index])
        if bootstrap_coefficients is None:
            lower = coefficient - critical * standard_error
            upper = coefficient + critical * standard_error
        else:
            alpha = (1.0 - confidence_level) / 2.0
            lower, upper = np.quantile(
                bootstrap_coefficients[:, index], [alpha, 1.0 - alpha]
            )
            lower = float(lower)
            upper = float(upper)
        if standard_error == 0.0:
            statistic: float | str = ""
            p_value: float | str = ""
            inference_status = "zero-model-variance"
        else:
            statistic_value = coefficient / standard_error
            statistic = statistic_value
            if bootstrap_coefficients is None:
                p_value = float(
                    2.0 * student_t.sf(abs(statistic_value), degrees_of_freedom)
                )
            else:
                centered = bootstrap_coefficients[:, index] - coefficient
                p_value = float(
                    (1 + np.sum(np.abs(centered) >= abs(coefficient)))
                    / (len(centered) + 1)
                )
            inference_status = "ok"
        rows.append(
            {
                "model_id": model_id,
                "tree_id": tree_id,
                "response": response,
                "response_type": "continuous",
                "response_family": "gaussian",
                "response_level": "",
                "response_reference": "",
                "link_function": "identity",
                "term": predictor,
                **_predictor_term_fields(predictor_metadata[predictor]),
                "coefficient": coefficient,
                "standard_error": standard_error,
                "statistic": statistic,
                "degrees_of_freedom": degrees_of_freedom,
                "p_value": p_value,
                "confidence_level": confidence_level,
                "confidence_interval_lower": lower,
                "confidence_interval_upper": upper,
                "n_gene_contrasts": n_observations,
                "n_species_events": n_events,
                "n_repeated_gene_contrasts": n_observations - n_events,
                "n_lineages": len(unique_lineages),
                "n_excluded_ineligible": excluded_ineligible,
                "n_excluded_coverage": excluded_coverage,
                "num_parameters": num_parameters,
                "matrix_rank": int(np.linalg.matrix_rank(design)),
                "condition_number": condition_number,
                "weighted_residual_sum_squares": fit["quadratic"],
                "residual_scale": evolutionary_rate,
                "r_squared_uncentered": r_squared,
                "intercept": "no",
                "event_weighting": event_weighting,
                "covariance_estimator": "gaussian{}-{}".format(
                    "-eiv" if predictor_uncertainties else "",
                    "REML" if effective_reml else "ML",
                ),
                "contrast_transform": "gene-evolutionary-plus-sampling-covariance",
                "coverage_policy": coverage_policy,
                "small_sample_warning": "yes" if n_events < 20 else "no",
                "inference_status": inference_status,
                "model": model,
                "inference_method": inference,
                "reml": "yes" if effective_reml else "no",
                "evolutionary_rate": evolutionary_rate,
                "species_event_variance": event_variance,
                "lineage_slope_variance": lineage_variance,
                "mean_sampling_variance": mean_sampling_variance,
                "sampling_variance_fraction": sampling_fraction,
                **_predictor_measurement_result_fields(
                    predictor,
                    predictor_uncertainties,
                    predictor_posteriors,
                    grouped_uncertainties_by_term,
                    mean_sampling_variance,
                ),
                "log_likelihood": -float(fit["objective"]),
                "optimizer_converged": "yes" if fit["optimizer_converged"] else "no",
                "boundary_warning": "yes" if fit["boundary_warning"] else "no",
                "event_random_effect": "yes" if use_event else "no",
                "lineage_random_slope": "yes" if use_lineage else "no",
            }
        )
    _append_predictor_omnibus_rows(
        rows,
        beta,
        beta_covariance,
        predictors,
        predictor_metadata,
        predictor_groups,
    )
    _append_lineage_inference_rows(
        rows,
        response_values,
        design,
        predictors,
        predictor_groups,
        fixed_covariance,
        components,
        component_factors,
        random_designs,
        balanced_predictor_uncertainties,
        predictor_uncertainty_columns,
        lineage_inference=lineage_inference,
        bootstrap_replicates=bootstrap_replicates,
        seed=seed,
        allow_large_dense=allow_large_dense,
    )
    random_effect_rows = []
    if use_event:
        variance = event_variance
        modes = variance * random_designs["species_event"].T @ fit["inverse_residual"]
        for index, event_id in enumerate(unique_events):
            random_effect_rows.append(
                {
                    "model_id": model_id,
                    "tree_id": tree_id,
                    "response": response,
                    "effect_type": "species_event",
                    "group_id": event_id,
                    "term": "(intercept)",
                    "conditional_mode": float(modes[index]),
                    "variance_component": variance,
                    "n_observations": int(event_counts[index]),
                    "inference_status": "empirical-bayes-conditional-on-variance",
                }
            )
    if use_lineage:
        random_effect_rows.extend(
            _lineage_random_effect_rows(
                fit,
                design,
                random_designs,
                predictors,
                unique_lineages,
                lineage_counts,
                model_id=model_id,
                tree_id=tree_id,
                response=response,
                variance=lineage_variance,
                confidence_level=confidence_level,
            )
        )
    if return_fit_state:
        fit_state = {
            "beta": np.asarray(beta, dtype=float),
            "design": np.asarray(design, dtype=float),
            "contrast_ids": dataframe["gene_clade_id"].astype(str).tolist(),
            "evolutionary_rate": float(evolutionary_rate),
            "fitted_covariance_factor": fit["cholesky"],
        }
    else:
        fit_state = None
    return rows, random_effect_rows, fit_state


def _validate_reconciled_pgls_options(
    *,
    confidence_level=0.95,
    event_weighting="equal",
    coverage_policy="complete",
    model="hierarchical",
    response_sampling_covariance=None,
    predictor_sampling_covariance=None,
    inference="wald",
    bootstrap_replicates=1000,
    seed=1,
    reml=True,
    event_random_effect="auto",
    lineage_random_slope="auto",
    lineage_inference="none",
    lineage_leave_one_out=False,
    return_fit_state=False,
):
    if not 0.0 < confidence_level < 1.0:
        raise ValueError("confidence_level must be between zero and one.")
    if event_weighting not in {"equal", "observation"}:
        raise ValueError("Unsupported event weighting: {}.".format(event_weighting))
    if coverage_policy not in {"complete", "any"}:
        raise ValueError("Unsupported coverage policy: {}.".format(coverage_policy))
    if model not in {"hierarchical", "legacy", "replicate-reml"}:
        raise ValueError("Unsupported PGLS model: {}.".format(model))
    if inference not in {"parametric-bootstrap", "wald"}:
        raise ValueError("Unsupported inference method: {}.".format(inference))
    if (
        not isinstance(bootstrap_replicates, int)
        or isinstance(bootstrap_replicates, bool)
        or bootstrap_replicates < 2
    ):
        raise ValueError("bootstrap_replicates must be an integer of at least two.")
    if not isinstance(seed, int) or isinstance(seed, bool) or seed < 0:
        raise ValueError("Bootstrap seed must be a non-negative integer.")
    if not isinstance(reml, bool):
        raise ValueError("reml must be a boolean.")
    if not isinstance(lineage_leave_one_out, bool):
        raise ValueError("lineage_leave_one_out must be a boolean.")
    for policy, label in [
        (event_random_effect, "event_random_effect"),
        (lineage_random_slope, "lineage_random_slope"),
    ]:
        if policy not in {"auto", "yes", "no"}:
            raise ValueError("{} must be auto, yes, or no.".format(label))
    if lineage_inference not in {
        "none",
        "likelihood-ratio",
        "parametric-bootstrap",
    }:
        raise ValueError(
            "lineage_inference must be none, likelihood-ratio, or parametric-bootstrap."
        )
    if model == "legacy" and (
        response_sampling_covariance is not None
        or predictor_sampling_covariance is not None
    ):
        raise ValueError("Sampling covariance requires a likelihood-based PGLS model.")
    if model == "legacy" and inference != "wald":
        raise ValueError(
            "Parametric-bootstrap inference is unavailable for legacy PGLS."
        )
    if model == "legacy" and return_fit_state:
        raise ValueError("Fit-state output requires a likelihood-based PGLS model.")
    if model != "hierarchical" and (
        event_random_effect == "yes" or lineage_random_slope == "yes"
    ):
        raise ValueError(
            "Requested random effects require the hierarchical PGLS model."
        )
    if model != "hierarchical" and lineage_inference != "none":
        raise ValueError("Lineage inference requires the hierarchical PGLS model.")
    if model == "legacy" and lineage_leave_one_out:
        raise ValueError(
            "Lineage leave-one-out requires a likelihood-based PGLS model."
        )
    if lineage_random_slope == "no" and lineage_inference != "none":
        raise ValueError("Lineage inference requires lineage random slopes.")


def _unique_trait_names(values, label):
    if isinstance(values, (str, bytes)):
        raise ValueError("Responses and predictors must be sequences of names.")
    names = list(values)
    if not names or len(names) != len(set(names)):
        raise ValueError(
            "{} must be a non-empty sequence of unique names.".format(label)
        )
    return names


def _normalize_predictor_metadata(predictors, metadata):
    normalized = (
        {
            predictor: PredictorTerm(predictor, predictor, "continuous")
            for predictor in predictors
        }
        if metadata is None
        else dict(metadata)
    )
    if set(normalized) != set(predictors):
        raise ValueError("Predictor metadata must describe every encoded predictor.")
    return normalized


def _normalize_predictor_groups(predictors, groups):
    normalized = (
        {predictor: (predictor,) for predictor in predictors}
        if groups is None
        else {str(source): tuple(terms) for source, terms in groups.items()}
    )
    grouped_terms = [term for terms in normalized.values() for term in terms]
    if sorted(grouped_terms) != sorted(predictors):
        raise ValueError("Predictor groups must partition the encoded predictors.")
    return normalized


def _normalize_one_grouped_uncertainty(state, predictors):
    normalized = dict(state)
    term_names = tuple(normalized.get("term_names", ()))
    event_ids = [str(value) for value in normalized.get("event_ids", ())]
    uncertainty = normalized.get("uncertainty")
    factors = (
        ()
        if not isinstance(uncertainty, JointPredictorUncertainty)
        else tuple(
            sparse.csr_matrix(factor, dtype=float) for factor in uncertainty.factors
        )
    )
    expected_shape = None if not factors else factors[0].shape
    valid = (
        bool(term_names)
        and set(term_names) <= set(predictors)
        and len(event_ids) == len(set(event_ids))
        and len(factors) == len(term_names)
        and expected_shape is not None
        and expected_shape[0] == len(event_ids)
        and all(
            factor.shape == expected_shape and np.isfinite(factor.data).all()
            for factor in factors
        )
    )
    if not valid:
        raise ValueError("Grouped predictor uncertainty is malformed.")
    normalized.update(
        {
            "term_names": term_names,
            "event_ids": event_ids,
            "event_index": {
                event_id: index for index, event_id in enumerate(event_ids)
            },
            "uncertainty": JointPredictorUncertainty(factors=factors),
        }
    )
    return normalized


def _normalize_grouped_uncertainties(states, predictors, model):
    normalized = [
        _normalize_one_grouped_uncertainty(state, predictors) for state in states or ()
    ]
    if model == "legacy" and normalized:
        raise ValueError(
            "Grouped predictor uncertainty requires a likelihood-based PGLS model."
        )
    return normalized


def _precomputed_transform_fields(response_metadata, predictor_metadata, tree_id):
    return {
        "predictor_tree_id": tree_id,
        "response_evolution_model": response_metadata["evolution_model"],
        "response_evolution_parameter_name": response_metadata[
            "evolution_parameter_name"
        ],
        "response_evolution_parameter": response_metadata["evolution_parameter"],
        "response_evolution_parameter_status": (
            "recorded"
            if response_metadata["evolution_parameter_name"]
            else "not-applicable"
        ),
        "response_evolution_optimizer_converged": "not-applicable",
        "response_evolution_optimizer_message": "precomputed transform",
        "response_evolution_boundary_warning": "not-applicable",
        "response_evolution_parameter_bootstrap_refit": "no",
        "response_branch_length_mode": response_metadata["branch_length_mode"],
        "predictor_evolution_model": predictor_metadata["evolution_model"],
        "predictor_evolution_parameter_name": predictor_metadata[
            "evolution_parameter_name"
        ],
        "predictor_evolution_parameter": predictor_metadata["evolution_parameter"],
        "predictor_evolution_parameter_status": (
            "recorded"
            if predictor_metadata["evolution_parameter_name"]
            else "not-applicable"
        ),
        "predictor_evolution_optimizer_converged": "not-applicable",
        "predictor_evolution_optimizer_message": "precomputed transform",
        "predictor_evolution_boundary_warning": "not-applicable",
        "predictor_evolution_log_likelihood": "",
        "predictor_branch_length_mode": predictor_metadata["branch_length_mode"],
    }


def _attach_precomputed_transform_fields(
    rows,
    response_metadata,
    predictor_metadata,
    predictor_groups,
    predictor_tree_id,
):
    for row in rows:
        response_transform = response_metadata[
            (str(row["tree_id"]), str(row["response"]))
        ]
        evolution_term = (
            str(row["term"])
            if row["term_test"] == "coefficient"
            else (
                next(iter(predictor_groups.values()))[0]
                if row["term_test"] == "lineage-heterogeneity"
                else predictor_groups[str(row["source_term"])][0]
            )
        )
        row.update(
            _precomputed_transform_fields(
                response_transform,
                predictor_metadata[(evolution_term,)],
                predictor_tree_id,
            )
        )


def _reconciled_return_value(
    result,
    random_effects,
    sensitivity,
    fit_states,
    *,
    return_random_effects,
    return_sensitivity,
    return_fit_state,
):
    outputs = [result]
    if return_random_effects:
        outputs.append(random_effects)
    if return_sensitivity:
        outputs.append(sensitivity)
    if return_fit_state:
        outputs.append(fit_states)
    return outputs[0] if len(outputs) == 1 else tuple(outputs)


def _sensitivity_model_groups(prepared_responses, omissions, lineage_leave_one_out):
    groups = []
    for (tree_id, response), model_rows in prepared_responses.groupby(
        ["tree_id", "trait"], sort=True, dropna=False
    ):
        if lineage_leave_one_out:
            for lineage_id in sorted(
                model_rows["lineage_clade_id"].astype(str).unique()
            ):
                groups.append(
                    {
                        "analysis_type": "lineage-leave-one-out",
                        "group_id": lineage_id,
                        "group_label": lineage_id,
                        "tree_id": str(tree_id),
                        "response": str(response),
                        "lineage_ids": {lineage_id},
                        "event_ids": set(),
                    }
                )
        for omission in omissions or ():
            groups.append(
                {
                    **omission,
                    "tree_id": str(tree_id),
                    "response": str(response),
                    "lineage_ids": set(omission.get("lineage_ids", ())),
                    "event_ids": set(omission.get("event_ids", ())),
                }
            )
    return groups


def _refit_random_effect_policy(policy):
    return "auto" if policy == "yes" else policy


def _relative_coefficient_change(full_coefficient, coefficient_change):
    if full_coefficient == 0.0:
        return float("nan")
    return coefficient_change / abs(full_coefficient)


def _coefficient_sign_changed(full_coefficient, omitted_coefficient):
    return bool(
        np.signbit(omitted_coefficient) != np.signbit(full_coefficient)
        and omitted_coefficient != 0.0
        and full_coefficient != 0.0
    )


def _coefficient_result_rows(result):
    return result[result["term_test"] == "coefficient"]


def _compute_reconciled_sensitivity(
    prepared_responses,
    prepared_predictors,
    responses,
    predictors,
    full_result,
    *,
    omissions,
    lineage_leave_one_out,
    fit_options,
):
    sensitivity_rows = []
    model_groups = _sensitivity_model_groups(
        prepared_responses, omissions, lineage_leave_one_out
    )
    for omission_index, omission in enumerate(model_groups):
        selected = prepared_responses[
            (prepared_responses["tree_id"].astype(str) == omission["tree_id"])
            & (prepared_responses["trait"].astype(str) == omission["response"])
        ].copy()
        remove = selected["lineage_clade_id"].astype(str).isin(
            omission["lineage_ids"]
        ) | selected["species_event_id"].astype(str).isin(omission["event_ids"])
        if not remove.any():
            continue
        omitted_rows = selected.loc[~remove].copy()
        coefficient_results = _coefficient_result_rows(full_result)
        full_coefficients = coefficient_results[
            (coefficient_results["tree_id"].astype(str) == omission["tree_id"])
            & (coefficient_results["response"].astype(str) == omission["response"])
        ]
        omitted_event_count = len(
            set(selected["species_event_id"].astype(str))
            - set(omitted_rows["species_event_id"].astype(str))
        )
        try:
            omitted_result = fit_reconciled_pgls(
                omitted_rows,
                prepared_predictors,
                [omission["response"]],
                predictors,
                **{
                    **fit_options,
                    "inference": "wald",
                    "lineage_inference": "none",
                    "lineage_leave_one_out": False,
                    "sensitivity_omissions": None,
                    "return_random_effects": False,
                    "return_sensitivity": False,
                    "return_fit_state": False,
                    "seed": fit_options["seed"] + omission_index + 1,
                    "event_random_effect": _refit_random_effect_policy(
                        fit_options["event_random_effect"]
                    ),
                    "lineage_random_slope": _refit_random_effect_policy(
                        fit_options["lineage_random_slope"]
                    ),
                },
            )
        except ValueError as exc:
            for record in full_coefficients.to_dict("records"):
                sensitivity_rows.append(
                    {
                        "analysis_type": omission["analysis_type"],
                        "model_id": record["model_id"],
                        "tree_id": omission["tree_id"],
                        "response": omission["response"],
                        "group_id": omission["group_id"],
                        "group_label": omission.get(
                            "group_label", omission["group_id"]
                        ),
                        "term": record["term"],
                        "source_term": record["source_term"],
                        "n_omitted_gene_contrasts": int(remove.sum()),
                        "n_omitted_species_events": omitted_event_count,
                        "n_retained_gene_contrasts": len(omitted_rows),
                        "n_retained_species_events": omitted_rows["species_event_id"]
                        .astype(str)
                        .nunique(),
                        "full_coefficient": record["coefficient"],
                        "inference_status": "refit-failed",
                        "message": str(exc),
                    }
                )
            continue
        omitted_coefficients = _coefficient_result_rows(omitted_result).set_index(
            "term"
        )
        for record in full_coefficients.to_dict("records"):
            term = str(record["term"])
            full_coefficient = float(record["coefficient"])
            if term not in omitted_coefficients.index:
                omitted_coefficient = float("nan")
                change = float("nan")
                relative_change = float("nan")
                sign_changed: bool | str = ""
                status = "coefficient-unavailable"
            else:
                omitted_coefficient = float(
                    omitted_coefficients.loc[term, "coefficient"]
                )
                change = omitted_coefficient - full_coefficient
                relative_change = _relative_coefficient_change(full_coefficient, change)
                sign_changed = _coefficient_sign_changed(
                    full_coefficient, omitted_coefficient
                )
                status = "ok"
            sensitivity_rows.append(
                {
                    "analysis_type": omission["analysis_type"],
                    "model_id": record["model_id"],
                    "tree_id": omission["tree_id"],
                    "response": omission["response"],
                    "group_id": omission["group_id"],
                    "group_label": omission.get("group_label", omission["group_id"]),
                    "term": term,
                    "source_term": record["source_term"],
                    "n_omitted_gene_contrasts": int(remove.sum()),
                    "n_omitted_species_events": omitted_event_count,
                    "n_retained_gene_contrasts": len(omitted_rows),
                    "n_retained_species_events": omitted_rows["species_event_id"]
                    .astype(str)
                    .nunique(),
                    "full_coefficient": full_coefficient,
                    "omitted_coefficient": omitted_coefficient,
                    "coefficient_change": change,
                    "relative_change": relative_change,
                    "sign_changed": sign_changed,
                    "inference_status": status,
                    "message": "",
                }
            )
    return pd.DataFrame(sensitivity_rows, columns=SENSITIVITY_COLUMNS)


def fit_reconciled_pgls(
    response_contrasts,
    predictor_contrasts,
    responses,
    predictors,
    *,
    confidence_level=0.95,
    event_weighting="equal",
    coverage_policy="complete",
    model="hierarchical",
    response_sampling_covariance=None,
    predictor_sampling_covariance=None,
    inference="wald",
    bootstrap_replicates=1000,
    seed=1,
    reml=True,
    event_random_effect="auto",
    lineage_random_slope="auto",
    lineage_inference="none",
    lineage_leave_one_out=False,
    sensitivity_omissions=None,
    return_random_effects=False,
    return_sensitivity=False,
    return_fit_state=False,
    predictor_metadata=None,
    predictor_groups=None,
    predictor_group_uncertainties=None,
    allow_large_dense=False,
):
    _validate_reconciled_pgls_options(
        confidence_level=confidence_level,
        event_weighting=event_weighting,
        coverage_policy=coverage_policy,
        model=model,
        response_sampling_covariance=response_sampling_covariance,
        predictor_sampling_covariance=predictor_sampling_covariance,
        inference=inference,
        bootstrap_replicates=bootstrap_replicates,
        seed=seed,
        reml=reml,
        event_random_effect=event_random_effect,
        lineage_random_slope=lineage_random_slope,
        lineage_inference=lineage_inference,
        lineage_leave_one_out=lineage_leave_one_out,
        return_fit_state=return_fit_state,
    )
    responses = _unique_trait_names(responses, "responses")
    predictors = _unique_trait_names(predictors, "predictors")
    predictor_metadata = _normalize_predictor_metadata(predictors, predictor_metadata)
    predictor_groups = _normalize_predictor_groups(predictors, predictor_groups)
    predictor_group_uncertainties = _normalize_grouped_uncertainties(
        predictor_group_uncertainties, predictors, model
    )
    prepared_responses = _prepare_responses(
        response_contrasts, responses, coverage_policy
    )
    prepared_predictors = _prepare_predictors(predictor_contrasts, predictors)
    response_evolution_metadata = _evolution_metadata_by_group(
        prepared_responses,
        ["tree_id", "trait"],
    )
    predictor_evolution_metadata = _evolution_metadata_by_group(
        prepared_predictors,
        ["trait"],
    )
    predictor_tree_id = str(prepared_predictors.iloc[0]["tree_id"])
    prepared_sampling = (
        None
        if response_sampling_covariance is None
        else _prepare_sampling_covariance(response_sampling_covariance)
    )
    prepared_predictor_sampling = (
        None
        if predictor_sampling_covariance is None
        else _prepare_sampling_covariance(
            predictor_sampling_covariance,
            option_name="--predictor-sampling-covariance",
        )
    )
    predictor_posteriors = _prepare_predictor_posteriors(
        prepared_predictors,
        predictors,
        prepared_predictor_sampling,
    )
    predictor_events = _predictor_by_event(prepared_predictors, predictors)
    rows = []
    random_effect_rows = []
    fit_states = {}
    group_columns = ["tree_id", "trait"]
    for model_index, ((tree_id, response), unfiltered) in enumerate(
        prepared_responses.groupby(group_columns, sort=True, dropna=False)
    ):
        ineligible = ~unfiltered["_eligible_for_model"]
        coverage_excluded = (
            unfiltered["_eligible_for_model"] & ~unfiltered["_coverage_for_model"]
        )
        filtered = unfiltered[
            unfiltered["_eligible_for_model"] & unfiltered["_coverage_for_model"]
        ].copy()
        if filtered.empty:
            raise ValueError(
                "Model '{}' has no rows after eligibility and coverage filtering.".format(
                    _model_id(str(tree_id), str(response))
                )
            )
        filtered = _validate_and_attach_predictors(
            filtered,
            predictor_events,
            predictors,
            predictor_posteriors,
        )
        if model == "legacy":
            model_rows = _fit_model(
                filtered,
                predictors,
                confidence_level=confidence_level,
                event_weighting=event_weighting,
                coverage_policy=coverage_policy,
                excluded_ineligible=int(ineligible.sum()),
                excluded_coverage=int(coverage_excluded.sum()),
                predictor_metadata=predictor_metadata,
                predictor_groups=predictor_groups,
            )
            for row in model_rows:
                row.update(
                    {
                        "model": "legacy",
                        "inference_method": "species-event-cluster-HC1",
                        "reml": "no",
                        "optimizer_converged": "not-applicable",
                        "boundary_warning": "not-applicable",
                        "event_random_effect": "no",
                        "lineage_random_slope": "no",
                    }
                )
            rows.extend(model_rows)
        else:
            contrast_ids = filtered["gene_clade_id"].astype(str).tolist()
            sampling_matrix = _sampling_matrix_for_model(
                prepared_sampling,
                str(tree_id),
                str(response),
                contrast_ids,
                filtered,
            )
            model_rows, model_random_effects, model_fit_state = _fit_covariance_model(
                filtered,
                predictors,
                sampling_matrix,
                _predictor_uncertainties_for_rows(
                    filtered,
                    predictors,
                    predictor_posteriors,
                ),
                predictor_posteriors,
                _grouped_predictor_uncertainties_for_rows(
                    filtered, predictor_group_uncertainties
                ),
                confidence_level=confidence_level,
                event_weighting=event_weighting,
                coverage_policy=coverage_policy,
                excluded_ineligible=int(ineligible.sum()),
                excluded_coverage=int(coverage_excluded.sum()),
                model=model,
                inference=inference,
                bootstrap_replicates=bootstrap_replicates,
                seed=seed + model_index,
                reml=reml,
                event_random_effect=event_random_effect,
                lineage_random_slope=lineage_random_slope,
                lineage_inference=lineage_inference,
                return_fit_state=return_fit_state,
                predictor_metadata=predictor_metadata,
                predictor_groups=predictor_groups,
                allow_large_dense=allow_large_dense,
            )
            rows.extend(model_rows)
            random_effect_rows.extend(model_random_effects)
            if return_fit_state:
                fit_states[(str(tree_id), str(response))] = model_fit_state
    if not rows:
        raise ValueError("No PGLS models were fitted.")
    _attach_precomputed_transform_fields(
        rows,
        response_evolution_metadata,
        predictor_evolution_metadata,
        predictor_groups,
        predictor_tree_id,
    )
    result = pd.DataFrame(rows, columns=RESULT_COLUMNS)
    random_effects = pd.DataFrame(random_effect_rows, columns=RANDOM_EFFECT_COLUMNS)
    if lineage_leave_one_out or sensitivity_omissions:
        sensitivity = _compute_reconciled_sensitivity(
            prepared_responses,
            prepared_predictors,
            responses,
            predictors,
            result,
            omissions=sensitivity_omissions,
            lineage_leave_one_out=lineage_leave_one_out,
            fit_options={
                "confidence_level": confidence_level,
                "event_weighting": event_weighting,
                "coverage_policy": coverage_policy,
                "model": model,
                "response_sampling_covariance": prepared_sampling,
                "predictor_sampling_covariance": prepared_predictor_sampling,
                "bootstrap_replicates": bootstrap_replicates,
                "seed": seed,
                "reml": reml,
                "event_random_effect": event_random_effect,
                "lineage_random_slope": lineage_random_slope,
                "predictor_metadata": predictor_metadata,
                "predictor_groups": predictor_groups,
                "predictor_group_uncertainties": predictor_group_uncertainties,
                "allow_large_dense": allow_large_dense,
            },
        )
    else:
        sensitivity = pd.DataFrame(columns=SENSITIVITY_COLUMNS)
    return _reconciled_return_value(
        result,
        random_effects,
        sensitivity,
        fit_states,
        return_random_effects=return_random_effects,
        return_sensitivity=return_sensitivity,
        return_fit_state=return_fit_state,
    )


RAW_PGLS_REQUIRED_ARGUMENTS = {
    "expression": "--expression",
    "species_traits": "--species-traits",
    "species_tree": "--species-tree",
    "tree_id": "--tree-id",
}

ORDINARY_PGLS_REQUIRED_ARGUMENTS = {
    "data": "--data",
    "tree": "--tree",
}

ORDINARY_PGLS_ONLY_ARGUMENTS = {
    "branch_length": "--branch-length",
    "compare_evolution_models": "--compare-evolution-models",
    "data": "--data",
    "evolution_covariance": "--evolution-covariance",
    "evolution_model": "--evolution-model",
    "evolution_parameter": "--evolution-parameter",
    "intercept": "--intercept",
    "model_comparison_out": "--model-comparison-out",
    "predictor_branch_length": "--predictor-branch-length",
    "predictor_evolution_model": "--predictor-evolution-model",
    "predictor_evolution_parameter": "--predictor-evolution-parameter",
    "predictor_sampling_covariance_out": "--predictor-sampling-covariance-out",
    "predictor_tip_summary_out": "--predictor-tip-summary-out",
    "sampling_covariance_out": "--sampling-covariance-out",
    "tip_summary_out": "--tip-summary-out",
    "tree": "--tree",
    "tree_format": "--tree-format",
}

RAW_PGLS_ONLY_ARGUMENTS = {
    "batch": "--batch",
    "biological_id": "--biological-id",
    "event_source": "--event-source",
    "expression": "--expression",
    "gene_branch_length": "--gene-branch-length",
    "gene_evolution_model": "--gene-evolution-model",
    "gene_evolution_parameter": "--gene-evolution-parameter",
    "gene_tree_format": "--gene-tree-format",
    "gene_tree_ensemble": "--gene-tree-ensemble",
    "missing_values": "--missing-values",
    "out_prefix": "--out-prefix",
    "categorical_origin_diagnostics": "--categorical-origin-diagnostics",
    "origin_map_replicates": "--origin-map-replicates",
    "origin_map_threads": "--origin-map-threads",
    "origin_min_posterior": "--origin-min-posterior",
    "origin_leave_one_out": "--origin-leave-one-out",
    "quoted_node_names": "--quoted-node-names",
    "reconciliation_tree": "--reconciliation-tree",
    "reconciliation_tree_format": "--reconciliation-tree-format",
    "sample_size_columns": "--sample-size-columns",
    "species_branch_length": "--species-branch-length",
    "species_evolution_model": "--species-evolution-model",
    "species_evolution_parameter": "--species-evolution-parameter",
    "species_map_tsv": "--species-map-tsv",
    "species_parser": "--species-parser",
    "species_regex": "--species-regex",
    "species_traits": "--species-traits",
    "species_tree": "--species-tree",
    "species_tree_format": "--species-tree-format",
    "predictor_batch": "--predictor-batch",
    "predictor_biological_id": "--predictor-biological-id",
    "predictor_sample_size_columns": "--predictor-sample-size-columns",
    "predictor_standard_error_columns": "--predictor-standard-error-columns",
    "predictor_technical_aggregation": "--predictor-technical-aggregation",
    "predictor_technical_id": "--predictor-technical-id",
    "predictor_within_variance": "--predictor-within-variance",
    "standard_error_columns": "--standard-error-columns",
    "technical_aggregation": "--technical-aggregation",
    "technical_id": "--technical-id",
    "tree_id": "--tree-id",
    "trait_origins_out": "--trait-origins-out",
    "unmatched": "--unmatched",
    "within_variance": "--within-variance",
}


def _nonempty_argument(args, name):
    value = getattr(args, name, None)
    return value is not None and (not isinstance(value, str) or value.strip() != "")


def _validate_evolution_cli_pair(
    args,
    model_name,
    parameter_name,
    model_option,
    parameter_option,
    *,
    parameter_required,
):
    model = getattr(args, model_name, None) or "brownian"
    parameter = getattr(args, parameter_name, None)
    spec = evolution_model_spec(model)
    if parameter == "auto":
        if spec.parameter_name is None:
            raise ValueError(
                "'{} {}' has no shape parameter to estimate.".format(
                    model_option, model
                )
            )
    elif parameter is not None:
        try:
            validate_evolution_parameter(model, parameter)
        except ValueError as exc:
            raise ValueError(
                "Invalid {} for {}: {}".format(parameter_option, model_option, exc)
            ) from exc
    if parameter_required and spec.parameter_name is not None and parameter is None:
        raise ValueError(
            "'{} {}' requires '{}'.".format(
                model_option,
                model,
                parameter_option,
            )
        )
    return model


def _require_pgls_arguments(args, required, description):
    missing = [
        option
        for name, option in required.items()
        if not _nonempty_argument(args, name)
    ]
    if missing:
        raise ValueError(
            "{} requires: {}.".format(description, ", ".join(sorted(missing)))
        )


def _require_biological_id_for_variance(args):
    within_variance = getattr(args, "within_variance", None)
    if within_variance in {"pooled", "leaf"} and not _nonempty_argument(
        args, "biological_id"
    ):
        raise ValueError(
            "'--within-variance {}' requires '--biological-id'.".format(within_variance)
        )
    predictor_within = getattr(args, "predictor_within_variance", None)
    if predictor_within in {"pooled", "leaf"} and not _nonempty_argument(
        args, "predictor_biological_id"
    ):
        raise ValueError(
            "'--predictor-within-variance {}' requires "
            "'--predictor-biological-id'.".format(predictor_within)
        )


def _ordinary_incompatible_options(args):
    incompatible = [
        option
        for name, option in {
            "event_source": "--event-source",
            "expression": "--expression",
            "gene_branch_length": "--gene-branch-length",
            "gene_evolution_model": "--gene-evolution-model",
            "gene_evolution_parameter": "--gene-evolution-parameter",
            "gene_tree": "--gene-tree",
            "gene_tree_ensemble": "--gene-tree-ensemble",
            "gene_tree_format": "--gene-tree-format",
            "infile": "--infile",
            "out_prefix": "--out-prefix",
            "predictor_contrasts": "--predictor-contrasts",
            "random_effects_out": "--random-effects-out",
            "sensitivity_out": "--sensitivity-out",
            "trait_origins_out": "--trait-origins-out",
            "reconciliation_tree": "--reconciliation-tree",
            "reconciliation_tree_format": "--reconciliation-tree-format",
            "response_sampling_covariance": "--response-sampling-covariance",
            "predictor_sampling_covariance": "--predictor-sampling-covariance",
            "species_branch_length": "--species-branch-length",
            "species_evolution_model": "--species-evolution-model",
            "species_evolution_parameter": "--species-evolution-parameter",
            "species_map_tsv": "--species-map-tsv",
            "species_parser": "--species-parser",
            "species_regex": "--species-regex",
            "species_traits": "--species-traits",
            "species_tree": "--species-tree",
            "species_tree_format": "--species-tree-format",
            "tree_id": "--tree-id",
        }.items()
        if _nonempty_argument(args, name)
    ]
    incompatible.extend(
        option
        for name, option in [
            ("event_weighting", "--event-weighting"),
            ("speciation_coverage", "--speciation-coverage"),
            ("model", "--model"),
            ("event_random_effect", "--event-random-effect"),
            ("lineage_random_slope", "--lineage-random-slope"),
            ("lineage_inference", "--lineage-inference"),
            ("lineage_leave_one_out", "--lineage-leave-one-out"),
            (
                "categorical_origin_diagnostics",
                "--categorical-origin-diagnostics",
            ),
            ("origin_map_replicates", "--origin-map-replicates"),
            ("origin_map_threads", "--origin-map-threads"),
            ("origin_min_posterior", "--origin-min-posterior"),
            ("origin_leave_one_out", "--origin-leave-one-out"),
        ]
        if _nonempty_argument(args, name)
    )
    return incompatible


def _validate_ordinary_custom_model_options(args, evolution_model):
    for name, option in [
        ("compare_evolution_models", "--compare-evolution-models"),
        ("model_comparison_out", "--model-comparison-out"),
    ]:
        if getattr(args, name, None) is not None and not _nonempty_argument(args, name):
            raise ValueError("'{}' must not be empty.".format(option))
    comparison_value = getattr(args, "compare_evolution_models", None)
    comparison_out = getattr(args, "model_comparison_out", None)
    if (comparison_value is None) != (comparison_out is None):
        raise ValueError(
            "'--compare-evolution-models' and '--model-comparison-out' must be used together."
        )
    comparison_models = (
        []
        if comparison_value is None
        else [item.strip() for item in str(comparison_value).split(",")]
    )
    uses_custom = (
        evolution_model == "custom"
        or getattr(args, "predictor_evolution_model", None) == "custom"
        or "custom" in comparison_models
    )
    covariance_supplied = _nonempty_argument(args, "evolution_covariance")
    if uses_custom and not covariance_supplied:
        raise ValueError(
            "The custom evolution model requires '--evolution-covariance'."
        )
    if covariance_supplied and not uses_custom:
        raise ValueError(
            "'--evolution-covariance' requires a custom selected or comparison model."
        )


def _validate_ordinary_pgls_mode(args):
    _require_pgls_arguments(args, ORDINARY_PGLS_REQUIRED_ARGUMENTS, "Conventional PGLS")
    incompatible = _ordinary_incompatible_options(args)
    if incompatible:
        raise ValueError(
            "Conventional PGLS cannot use reconciled/precomputed option(s): {}.".format(
                ", ".join(sorted(incompatible))
            )
        )
    evolution_model = _validate_evolution_cli_pair(
        args,
        "evolution_model",
        "evolution_parameter",
        "--evolution-model",
        "--evolution-parameter",
        parameter_required=False,
    )
    _validate_ordinary_custom_model_options(args, evolution_model)
    predictor_model = getattr(args, "predictor_evolution_model", None)
    predictor_parameter = getattr(args, "predictor_evolution_parameter", None)
    if predictor_model is not None:
        _validate_evolution_cli_pair(
            args,
            "predictor_evolution_model",
            "predictor_evolution_parameter",
            "--predictor-evolution-model",
            "--predictor-evolution-parameter",
            parameter_required=False,
        )
    elif predictor_parameter is not None:
        try:
            validate_evolution_parameter(evolution_model, predictor_parameter)
        except ValueError as exc:
            raise ValueError(
                "Invalid --predictor-evolution-parameter for inherited "
                "--evolution-model: {}".format(exc)
            ) from exc
    _require_biological_id_for_variance(args)


def _validate_raw_output_options(args):
    if getattr(args, "out_prefix", None) is not None and not _nonempty_argument(
        args, "out_prefix"
    ):
        raise ValueError("'--out-prefix' must not be empty.")
    if not _nonempty_argument(args, "out_prefix"):
        return
    if getattr(args, "outfile", "-") != "-":
        raise ValueError("'--out-prefix' cannot be combined with '--outfile'.")
    if _nonempty_argument(args, "random_effects_out"):
        raise ValueError(
            "'--out-prefix' cannot be combined with '--random-effects-out'; "
            "the bundle includes random effects."
        )
    if _nonempty_argument(args, "sensitivity_out"):
        raise ValueError(
            "'--out-prefix' cannot be combined with '--sensitivity-out'; "
            "the bundle includes sensitivity diagnostics."
        )
    if _nonempty_argument(args, "trait_origins_out"):
        raise ValueError(
            "'--out-prefix' cannot be combined with '--trait-origins-out'; "
            "the bundle includes trait-origin diagnostics."
        )


def _validate_raw_pgls_mode(args):
    _require_pgls_arguments(args, RAW_PGLS_REQUIRED_ARGUMENTS, "Raw-input PGLS")
    gene_tree = _nonempty_argument(args, "gene_tree")
    ensemble = _nonempty_argument(args, "gene_tree_ensemble")
    if gene_tree == ensemble:
        raise ValueError(
            "Raw-input PGLS requires exactly one of '--gene-tree' or "
            "'--gene-tree-ensemble'."
        )
    if ensemble and _nonempty_argument(args, "reconciliation_tree"):
        raise ValueError(
            "--gene-tree-ensemble cannot use one fixed --reconciliation-tree; "
            "embed reconciliation annotations in each sampled tree."
        )
    incompatible = [
        option
        for name, option in [
            ("infile", "--infile"),
            ("predictor_contrasts", "--predictor-contrasts"),
            ("response_sampling_covariance", "--response-sampling-covariance"),
            ("predictor_sampling_covariance", "--predictor-sampling-covariance"),
        ]
        if _nonempty_argument(args, name)
    ]
    if incompatible:
        raise ValueError(
            "Raw-input PGLS cannot be combined with precomputed input(s): {}.".format(
                ", ".join(incompatible)
            )
        )
    if _nonempty_argument(
        args, "reconciliation_tree_format"
    ) and not _nonempty_argument(args, "reconciliation_tree"):
        raise ValueError(
            "'--reconciliation-tree-format' requires '--reconciliation-tree'."
        )
    gene_model = _validate_evolution_cli_pair(
        args,
        "gene_evolution_model",
        "gene_evolution_parameter",
        "--gene-evolution-model",
        "--gene-evolution-parameter",
        parameter_required=False,
    )
    _validate_evolution_cli_pair(
        args,
        "species_evolution_model",
        "species_evolution_parameter",
        "--species-evolution-model",
        "--species-evolution-parameter",
        parameter_required=False,
    )
    gene_parameter = getattr(args, "gene_evolution_parameter", None)
    if (
        (getattr(args, "model", None) or "hierarchical") == "legacy"
        and evolution_model_spec(gene_model).parameter_name is not None
        and gene_parameter in {None, "auto"}
    ):
        raise ValueError(
            "Automatic gene evolution-parameter estimation requires a likelihood-based "
            "reconciled model, not '--model legacy'."
        )
    _require_biological_id_for_variance(args)
    origin_mode = getattr(args, "categorical_origin_diagnostics", None) or "none"
    origin_options = [
        option
        for name, option in [
            ("origin_map_replicates", "--origin-map-replicates"),
            ("origin_map_threads", "--origin-map-threads"),
            ("origin_min_posterior", "--origin-min-posterior"),
            ("trait_origins_out", "--trait-origins-out"),
        ]
        if _nonempty_argument(args, name)
    ]
    if bool(getattr(args, "origin_leave_one_out", False)):
        origin_options.append("--origin-leave-one-out")
    if origin_mode == "none" and origin_options:
        raise ValueError(
            "Categorical origin option(s) require "
            "'--categorical-origin-diagnostics stochastic-map': {}.".format(
                ", ".join(origin_options)
            )
        )
    for name, option in [
        ("origin_map_replicates", "--origin-map-replicates"),
        ("origin_map_threads", "--origin-map-threads"),
    ]:
        value = getattr(args, name, None)
        if value is not None and (isinstance(value, bool) or value <= 0):
            raise ValueError("'{}' must be a positive integer.".format(option))
    minimum = getattr(args, "origin_min_posterior", None)
    if minimum is not None and not 0.0 <= minimum <= 1.0:
        raise ValueError("'--origin-min-posterior' must lie in [0, 1].")
    _validate_raw_output_options(args)


def _pgls_input_mode(args):
    if any(_nonempty_argument(args, name) for name in ORDINARY_PGLS_ONLY_ARGUMENTS):
        _validate_ordinary_pgls_mode(args)
        return "ordinary"
    if _nonempty_argument(args, "gene_tree") or _nonempty_argument(
        args, "gene_tree_ensemble"
    ):
        _validate_raw_pgls_mode(args)
        return "raw"
    raw_options = [
        option
        for name, option in RAW_PGLS_ONLY_ARGUMENTS.items()
        if _nonempty_argument(args, name)
    ]
    if raw_options:
        raise ValueError(
            "Raw-input option(s) require '--gene-tree' or '--gene-tree-ensemble': {}.".format(
                ", ".join(sorted(raw_options))
            )
        )
    missing = [
        option
        for name, option in [
            ("infile", "--infile"),
            ("predictor_contrasts", "--predictor-contrasts"),
        ]
        if not _nonempty_argument(args, name)
    ]
    if missing:
        raise ValueError(
            "Precomputed-contrast PGLS requires: {}.".format(", ".join(missing))
        )
    return "contrasts"


def _warn_pgls_diagnostics(results):
    if (results["small_sample_warning"] == "yes").any():
        sys.stderr.write(
            "Warning: at least one model has fewer than 20 unique species events; "
            "small-sample inference may be unstable.\n"
        )
    if (results["boundary_warning"] == "yes").any():
        sys.stderr.write(
            "Warning: at least one variance component is near the optimization "
            "boundary; inspect variance-component diagnostics.\n"
        )
    if (results["optimizer_converged"] == "no").any():
        sys.stderr.write(
            "Warning: at least one variance-component optimizer did not converge; "
            "do not interpret that model without further diagnosis.\n"
        )
    if (
        "predictor_rate_boundary_warning" in results
        and (results["predictor_rate_boundary_warning"] == "yes").any()
    ):
        sys.stderr.write(
            "Warning: at least one latent-predictor evolutionary rate is near "
            "its optimization boundary.\n"
        )
    if (
        "predictor_rate_optimizer_converged" in results
        and (results["predictor_rate_optimizer_converged"] == "no").any()
    ):
        sys.stderr.write(
            "Warning: at least one latent-predictor evolutionary-rate optimizer "
            "did not converge.\n"
        )


def _warn_ordinary_pgls_diagnostics(results):
    if (results["small_sample_warning"] == "yes").any():
        sys.stderr.write(
            "Warning: at least one conventional PGLS model has fewer than 20 "
            "species; small-sample inference may be unstable.\n"
        )
    if (results["boundary_warning"] == "yes").any():
        sys.stderr.write(
            "Warning: at least one conventional PGLS variance or evolution-model "
            "parameter is near its optimization boundary.\n"
        )
    if (results["optimizer_converged"] == "no").any():
        sys.stderr.write(
            "Warning: at least one conventional PGLS optimizer did not converge; "
            "do not interpret that model without further diagnosis.\n"
        )
    if (
        "predictor_evolution_boundary_warning" in results
        and (results["predictor_evolution_boundary_warning"] == "yes").any()
    ):
        sys.stderr.write(
            "Warning: at least one conventional-PGLS latent predictor is near "
            "an evolution-model optimization boundary.\n"
        )


def _warn_ordinary_model_comparison(results):
    if results.empty:
        return
    if (results["boundary_warning"] == "yes").any():
        sys.stderr.write(
            "Warning: at least one compared evolutionary model has a parameter "
            "near its optimization boundary.\n"
        )
    if (results["optimizer_converged"] == "no").any():
        sys.stderr.write(
            "Warning: at least one compared evolutionary model did not converge.\n"
        )
    if pd.to_numeric(results["aicc"], errors="coerce").isna().any():
        sys.stderr.write(
            "Warning: AICc and AICc weights are unavailable where the sample "
            "size is not larger than the likelihood parameter count plus one.\n"
        )


def _write_pgls_outputs(
    args, results, random_effects, sensitivity=None, trait_origins=None
):
    _validate_pgls_file_output_paths(args)
    random_effects_path = getattr(args, "random_effects_out", None)
    file_outputs = []
    if random_effects_path is not None:
        file_outputs.append((random_effects_path, random_effects))
    sensitivity_path = getattr(args, "sensitivity_out", None)
    if sensitivity_path is not None:
        file_outputs.append(
            (
                sensitivity_path,
                pd.DataFrame(columns=SENSITIVITY_COLUMNS)
                if sensitivity is None
                else sensitivity,
            )
        )
    trait_origins_path = getattr(args, "trait_origins_out", None)
    if trait_origins_path is not None:
        from nwkit.rsc_diagnostics import ORIGIN_DIAGNOSTIC_COLUMNS

        file_outputs.append(
            (
                trait_origins_path,
                pd.DataFrame(columns=ORIGIN_DIAGNOSTIC_COLUMNS)
                if trait_origins is None
                else trait_origins,
            )
        )
    if args.outfile != "-":
        file_outputs.append((args.outfile, results))
    if file_outputs:
        from nwkit.pgls_pipeline import _write_dataframes_transactionally

        _write_dataframes_transactionally(file_outputs)
    if args.outfile == "-":
        print(results.to_csv(sep="\t", index=False), end="")


def _paths_identify_same_file(input_path, output_path):
    if normalized_missing_path_key(input_path) == normalized_missing_path_key(
        output_path
    ):
        return True
    if not os.path.exists(input_path) or not os.path.exists(output_path):
        return False
    try:
        return os.path.samefile(input_path, output_path)
    except OSError:
        return False


def _validate_pgls_file_output_paths(args):
    random_effects_path = getattr(args, "random_effects_out", None)
    sensitivity_path = getattr(args, "sensitivity_out", None)
    trait_origins_path = getattr(args, "trait_origins_out", None)
    if random_effects_path == "-":
        raise ValueError("'--random-effects-out' cannot be STDOUT.")
    if sensitivity_path == "-":
        raise ValueError("'--sensitivity-out' cannot be STDOUT.")
    if trait_origins_path == "-":
        raise ValueError("'--trait-origins-out' cannot be STDOUT.")
    output_paths = [
        args.outfile,
        random_effects_path,
        sensitivity_path,
        trait_origins_path,
    ]
    validate_distinct_output_paths(
        [
            ("--outfile", args.outfile),
            ("--random-effects-out", random_effects_path),
            ("--sensitivity-out", sensitivity_path),
            ("--trait-origins-out", trait_origins_path),
        ]
    )
    input_paths = [
        getattr(args, name, None)
        for name in [
            "expression",
            "gene_tree",
            "gene_tree_ensemble",
            "infile",
            "predictor_contrasts",
            "predictor_sampling_covariance",
            "reconciliation_tree",
            "response_sampling_covariance",
            "species_map_tsv",
            "species_traits",
            "species_tree",
        ]
    ]
    for input_path in input_paths:
        if input_path in (None, "", "-"):
            continue
        for output_path in output_paths:
            if output_path in (None, "", "-"):
                continue
            if _paths_identify_same_file(input_path, output_path):
                raise ValueError(
                    "PGLS output must not overwrite an input file: '{}'".format(
                        os.path.realpath(output_path)
                    )
                )


def pgls_main(args):
    responses = _parse_names(args.responses, "--responses")
    predictors = _parse_names(args.predictors, "--predictors")
    input_mode = _pgls_input_mode(args)
    if (
        input_mode != "ordinary"
        and (
            bool(getattr(args, "lineage_leave_one_out", False))
            or bool(getattr(args, "origin_leave_one_out", False))
        )
        and not _nonempty_argument(args, "out_prefix")
        and not _nonempty_argument(args, "sensitivity_out")
    ):
        raise ValueError(
            "Leave-one-out diagnostics require '--sensitivity-out' or '--out-prefix'."
        )
    if (
        input_mode == "raw"
        and getattr(args, "categorical_origin_diagnostics", None) == "stochastic-map"
        and not _nonempty_argument(args, "out_prefix")
        and not _nonempty_argument(args, "trait_origins_out")
    ):
        raise ValueError(
            "Categorical origin diagnostics require '--trait-origins-out' or "
            "'--out-prefix'."
        )
    if input_mode == "ordinary":
        from nwkit.ordinary_pgls import (
            build_ordinary_pgls,
            validate_ordinary_pgls_output_paths,
            write_ordinary_pgls_outputs,
        )

        validate_ordinary_pgls_output_paths(args)
        ordinary_artifacts = build_ordinary_pgls(args, responses, predictors)
        _warn_ordinary_pgls_diagnostics(ordinary_artifacts.results)
        _warn_ordinary_model_comparison(ordinary_artifacts.model_comparison)
        write_ordinary_pgls_outputs(args, ordinary_artifacts)
        return
    if input_mode == "raw":
        from nwkit.pgls_pipeline import (
            _active_pgls_bundle_paths,
            build_pgls_ensemble_pipeline,
            build_pgls_pipeline,
            validate_pgls_bundle_target,
            write_pgls_bundle,
        )

        out_prefix = getattr(args, "out_prefix", None)
        if out_prefix is not None:
            validate_pgls_bundle_target(
                out_prefix,
                protected_inputs=[
                    getattr(args, name, None)
                    for name in [
                        "expression",
                        "gene_tree",
                        "gene_tree_ensemble",
                        "reconciliation_tree",
                        "species_map_tsv",
                        "species_traits",
                        "species_tree",
                    ]
                ],
            )
        else:
            _validate_pgls_file_output_paths(args)
        pipeline_artifacts = (
            build_pgls_ensemble_pipeline(args, responses, predictors)
            if _nonempty_argument(args, "gene_tree_ensemble")
            else build_pgls_pipeline(args, responses, predictors)
        )
        _warn_pgls_diagnostics(pipeline_artifacts.results)
        if out_prefix is None:
            _write_pgls_outputs(
                args,
                pipeline_artifacts.results,
                pipeline_artifacts.random_effects,
                pipeline_artifacts.sensitivity,
                pipeline_artifacts.trait_origins,
            )
        else:
            written_paths = _active_pgls_bundle_paths(out_prefix, pipeline_artifacts)
            for argument, path in written_paths.items():
                setattr(args, argument, path)
            write_pgls_bundle(out_prefix, pipeline_artifacts)
        return
    unsupported_typed_options = [
        option
        for name, option in [
            ("categorical_responses", "--categorical-responses"),
            ("ordered_responses", "--ordered-responses"),
            ("response_reference", "--response-reference"),
            ("response_family", "--response-family"),
            ("response_offset", "--response-offset"),
            ("response_trials", "--response-trials"),
            ("response_censor_lower", "--response-censor-lower"),
            ("response_censor_upper", "--response-censor-upper"),
            ("response_dispersion", "--response-dispersion"),
            ("response_zero_probability", "--response-zero-probability"),
            ("categorical_predictors", "--categorical-predictors"),
            ("ordered_predictors", "--ordered-predictors"),
            ("factor_reference", "--factor-reference"),
        ]
        if _nonempty_argument(args, name)
    ]
    unsupported_typed_options.extend(
        option
        for name, option in [
            ("multivariate_responses", "--multivariate-responses"),
            ("allow_missing_responses", "--allow-missing-responses"),
        ]
        if bool(getattr(args, name, False))
    )
    if unsupported_typed_options:
        raise ValueError(
            "Typed response/predictor options require raw-input or conventional PGLS, "
            "not precomputed contrasts: {}.".format(
                ", ".join(unsupported_typed_options)
            )
        )
    _validate_pgls_file_output_paths(args)
    response_table = _read_tsv(args.infile, "--infile")
    predictor_table = _read_tsv(args.predictor_contrasts, "--predictor-contrasts")
    covariance_path = getattr(args, "response_sampling_covariance", None)
    sampling_covariance = (
        None
        if covariance_path is None
        else _read_tsv(covariance_path, "--response-sampling-covariance")
    )
    predictor_covariance_path = getattr(args, "predictor_sampling_covariance", None)
    predictor_sampling_covariance = (
        None
        if predictor_covariance_path is None
        else _read_tsv(
            predictor_covariance_path,
            "--predictor-sampling-covariance",
        )
    )
    results, random_effects, sensitivity = fit_reconciled_pgls(
        response_table,
        predictor_table,
        responses,
        predictors,
        confidence_level=args.confidence_level,
        event_weighting=getattr(args, "event_weighting", None) or "equal",
        coverage_policy=getattr(args, "speciation_coverage", None) or "complete",
        model=getattr(args, "model", None) or "hierarchical",
        response_sampling_covariance=sampling_covariance,
        predictor_sampling_covariance=predictor_sampling_covariance,
        inference=getattr(args, "inference", "wald"),
        bootstrap_replicates=getattr(args, "bootstrap_replicates", 1000),
        seed=getattr(args, "seed", 1),
        reml=getattr(args, "reml", True),
        event_random_effect=getattr(args, "event_random_effect", None) or "auto",
        lineage_random_slope=getattr(args, "lineage_random_slope", None) or "auto",
        lineage_inference=getattr(args, "lineage_inference", None) or "none",
        lineage_leave_one_out=bool(getattr(args, "lineage_leave_one_out", False)),
        allow_large_dense=getattr(args, "allow_large_dense", False),
        return_random_effects=True,
        return_sensitivity=True,
    )
    _warn_pgls_diagnostics(results)
    _write_pgls_outputs(args, results, random_effects, sensitivity)
