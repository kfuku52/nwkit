import math
import os
import sys
import warnings
from io import StringIO
from typing import Any

import numpy as np
import pandas as pd
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
    _predictor_error_covariance,
    fit_conditional_eiv_gaussian,
    fit_latent_predictor,
)
from nwkit.model_matrix import PredictorTerm
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
    "variance_component",
    "n_observations",
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
        return np.zeros((len(response_rows), len(response_rows)), dtype=float)
    selected = dataframe[
        (dataframe["tree_id"] == tree_id) & (dataframe["trait"] == response)
    ]
    if selected.empty:
        raise ValueError(
            "{} has no rows for model '{}'.".format(
                label,
                model_label or _model_id(tree_id, response),
            )
        )
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
                label,
                model_label or _model_id(tree_id, response),
                ", ".join(missing[:10]),
            )
        )
    return _validate_covariance_matrix(matrix, label)


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
        posterior = fit_latent_predictor(
            rows["raw_contrast"].to_numpy(dtype=float),
            np.diag(rows["contrast_variance"].to_numpy(dtype=float)),
            covariance,
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
            "sampling_covariance": covariance,
            "prior_covariance": prior_covariance,
        }
    return posteriors


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
        covariance = state["posterior"].covariance
        uncertainties[predictor] = covariance[np.ix_(indices, indices)]
    return uncertainties


def _grouped_predictor_uncertainties_for_rows(rows, states):
    event_ids = rows["species_event_id"].astype(str).tolist()
    selected = []
    for state in states or ():
        indices = [state["event_index"][event_id] for event_id in event_ids]
        covariance = np.asarray(state["covariance"], dtype=float)
        selected.append(
            {
                **state,
                "covariance": covariance[:, :, indices, :][:, :, :, indices],
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


def _profile_covariance_fit(
    y,
    design,
    fixed_covariance,
    components,
    *,
    reml,
    starting_log_variances=None,
    component_factors=None,
):
    n_observations = len(y)
    num_parameters = design.shape[1]
    component_matrices = []
    component_scales = []
    for name, matrix in components:
        matrix = np.asarray(matrix, dtype=float)
        diagonal = np.diag(matrix)
        positive = diagonal[diagonal > 0.0]
        if not len(positive):
            raise ValueError("Variance component '{}' has zero diagonal.".format(name))
        scale = float(np.mean(positive))
        component_matrices.append((name, matrix))
        component_scales.append(scale)
    fixed_covariance = _validate_covariance_matrix(
        fixed_covariance, "Fixed sampling covariance"
    )
    component_factors = {} if component_factors is None else component_factors
    unknown_factor_names = set(component_factors) - {
        name for name, _ in component_matrices
    }
    if unknown_factor_names:
        raise ValueError("Low-rank factors reference unknown variance components.")
    normalized_factors = {}
    for (name, matrix), scale in zip(component_matrices, component_scales, strict=True):
        if name not in component_factors:
            continue
        factor = np.asarray(component_factors[name], dtype=float)
        if factor.ndim != 2 or factor.shape[0] != n_observations:
            raise ValueError(
                "Low-rank factor for variance component '{}' has the wrong "
                "dimensions.".format(name)
            )
        if not np.isfinite(factor).all() or not np.allclose(
            np.sum(np.square(factor), axis=1), np.diag(matrix)
        ):
            raise ValueError(
                "Low-rank factor for variance component '{}' is inconsistent.".format(
                    name
                )
            )
        normalized_factors[name] = factor / math.sqrt(scale)
    structured_model = is_diagonal(fixed_covariance) and all(
        is_diagonal(matrix) or name in normalized_factors
        for name, matrix in component_matrices
    )
    if structured_model:
        working_fixed_covariance = np.diag(fixed_covariance).copy()
        normalized_components = [
            (
                name,
                None if name in normalized_factors else np.diag(matrix) / scale,
            )
            for (name, matrix), scale in zip(
                component_matrices, component_scales, strict=True
            )
        ]
    else:
        working_fixed_covariance = fixed_covariance
        normalized_components = [
            (name, matrix / scale)
            for (name, matrix), scale in zip(
                component_matrices, component_scales, strict=True
            )
        ]
    ordinary_beta = np.linalg.lstsq(design, y, rcond=None)[0]
    ordinary_residual = y - design @ ordinary_beta
    response_scale = max(
        float(np.mean(y**2)),
        float(np.mean(ordinary_residual**2)),
        float(np.mean(np.diag(fixed_covariance))) if len(y) else 0.0,
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
            low_rank_updates = []
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
                    low_rank = np.column_stack(low_rank_updates)
                    cholesky = factor_diagonal_low_rank_updates(
                        covariance, low_rank_updates
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

    if len(normalized_components) == 1 and not np.any(fixed_covariance):
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
    details = evaluate(result.x, return_details=True)
    if not isinstance(details, dict):
        raise ValueError("Variance-component optimization produced an invalid fit.")
    details["optimizer_converged"] = bool(result.success)
    details["optimizer_message"] = str(result.message)
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


def _component_adds_rank(components, candidate):
    size = candidate.shape[0]
    upper = np.triu_indices(size)
    existing = np.column_stack([matrix[upper] for _, matrix in components])
    candidate_vector = candidate[upper][:, None]
    scale = max(1.0, float(np.max(np.abs(candidate_vector))))
    tolerance = np.finfo(float).eps * scale * max(1, len(candidate_vector)) * 100.0
    augmented_rank = int(
        np.linalg.matrix_rank(
            np.column_stack([existing, candidate_vector]), tol=tolerance
        )
    )
    return augmented_rank > int(np.linalg.matrix_rank(existing, tol=tolerance))


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
            )
        except ValueError:
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
                reml=reml,
                starting_parameters=starting,
            )
        except ValueError:
            continue
        coefficients.append(bootstrap_fit["beta"])
    if len(coefficients) < replicates:
        raise ValueError(
            "Parametric bootstrap produced only {} successful errors-in-variables "
            "fits in {} attempts.".format(len(coefficients), attempts)
        )
    return np.asarray(coefficients, dtype=float)


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
    fixed_covariance = sampling_covariance.copy()
    if not np.all(balance == 1.0):
        fixed_covariance *= np.outer(balance, balance)
    raw_evolutionary_component = np.diag(evolutionary_variances)
    evolutionary_component = np.diag(evolutionary_variances * np.square(balance))
    components = [("evolutionary_rate", evolutionary_component)]
    component_factors = {}
    raw_components = {"evolutionary_rate": raw_evolutionary_component}
    random_designs = {}
    use_event = False
    use_lineage = False
    if model == "hierarchical":
        event_may_be_identifiable = (
            np.any(event_counts > 1)
            and int(np.sum(event_counts > 1)) >= 2
            and n_events > num_parameters
        )
        if event_may_be_identifiable:
            base_event_design = np.zeros(
                (n_observations, len(unique_events)), dtype=float
            )
            base_event_design[np.arange(n_observations), event_inverse] = 1.0
            event_design = balance[:, None] * base_event_design
            event_component = event_design @ event_design.T
            event_identifiable = _component_adds_rank(components, event_component)
        else:
            event_identifiable = False
        use_event = _random_effect_policy(
            event_random_effect, event_identifiable, "species-event random effect"
        )
        repeated_lineages = lineage_counts >= 2
        lineage_identifiable = (
            int(np.sum(repeated_lineages)) >= 2
            and n_events > num_parameters
            and np.any(np.abs(design) > 0.0)
        )
        use_lineage = _random_effect_policy(
            lineage_random_slope,
            lineage_identifiable,
            "lineage random slope",
        )
    if use_event:
        components.append(("species_event_variance", event_component))
        component_factors["species_event_variance"] = event_design
        raw_components["species_event_variance"] = (
            base_event_design @ base_event_design.T
        )
        random_designs["species_event"] = event_design
    if use_lineage:
        base_lineage_design = np.zeros(
            (n_observations, len(unique_lineages)), dtype=float
        )
        base_lineage_design[np.arange(n_observations), lineage_inverse] = balance
        lineage_designs = []
        lineage_component = np.zeros((n_observations, n_observations), dtype=float)
        raw_lineage_component = np.zeros((n_observations, n_observations), dtype=float)
        for predictor_index, predictor in enumerate(predictors):
            random_design = design[:, predictor_index, None] * base_lineage_design
            lineage_designs.append((predictor, random_design))
            lineage_component += random_design @ random_design.T
            raw_random_design = design[:, predictor_index, None] * (
                base_lineage_design / balance[:, None]
            )
            raw_lineage_component += raw_random_design @ raw_random_design.T
        if np.max(np.abs(lineage_component)) == 0.0:
            raise ValueError("Lineage random-slope component is identically zero.")
        if not _component_adds_rank(components, lineage_component):
            if lineage_random_slope == "yes":
                raise ValueError(
                    "Lineage random slope was requested but its covariance "
                    "component is not separately identifiable."
                )
            use_lineage = False
    if use_lineage:
        components.append(("lineage_slope_variance", lineage_component))
        component_factors["lineage_slope_variance"] = np.column_stack(
            [design for _, design in lineage_designs]
        )
        raw_components["lineage_slope_variance"] = raw_lineage_component
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
        covariance = state["covariance"][term_index, term_index]
        mean_variance = float(np.mean(np.diag(covariance)))
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
    mean_posterior_variance = float(
        np.mean(np.diag(predictor_uncertainties[predictor]))
    )
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
    raw: list[np.ndarray] = []
    weighted: list[np.ndarray] = []
    columns: list[int | tuple[int, ...]] = []
    grouped_by_term = {}
    outer_balance = np.outer(balance, balance)
    for predictor, uncertainty in predictor_uncertainties.items():
        covariance = np.asarray(uncertainty, dtype=float)
        raw.append(covariance)
        weighted.append(covariance * outer_balance)
        columns.append(index_by_predictor[predictor])
    for state in grouped_predictor_uncertainties:
        covariance = np.asarray(state["covariance"], dtype=float)
        raw.append(covariance)
        weighted.append(
            covariance * balance[None, None, :, None] * balance[None, None, None, :]
        )
        term_names = tuple(state["term_names"])
        columns.append(tuple(index_by_predictor[term] for term in term_names))
        for term_index, term in enumerate(term_names):
            grouped_by_term[term] = (state, term_index)
    return raw, weighted, columns, grouped_by_term


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
    return_fit_state,
    predictor_metadata,
    predictor_groups,
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
    sampling_covariance = _validate_covariance_matrix(
        sampling_covariance, "Response sampling covariance"
    )
    if sampling_covariance.shape != (n_observations, n_observations):
        raise ValueError("Response sampling covariance has the wrong dimensions.")
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
        fit = fit_conditional_eiv_gaussian(
            response_values,
            design,
            balanced_predictor_uncertainties,
            predictor_uncertainty_columns,
            fixed_covariance,
            components,
            reml=reml,
        )
    else:
        fit = _profile_covariance_fit(
            response_values,
            design,
            fixed_covariance,
            components,
            reml=reml,
            component_factors=component_factors,
        )
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
                fixed_covariance,
                components,
                reml=reml,
                replicates=bootstrap_replicates,
                seed=seed,
            )
        else:
            bootstrap_coefficients = _parametric_bootstrap_coefficients(
                fit,
                design,
                fixed_covariance,
                components,
                reml=reml,
                replicates=bootstrap_replicates,
                seed=seed,
                component_factors=component_factors,
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
    mean_sampling_variance = float(np.mean(np.diag(sampling_covariance)))
    fitted_variance = mean_sampling_variance
    for uncertainty, columns in zip(
        raw_predictor_uncertainties,
        predictor_uncertainty_columns,
        strict=True,
    ):
        fitted_variance += float(
            np.mean(
                np.diag(
                    _predictor_error_covariance(
                        beta, uncertainty, columns, n_observations
                    )
                )
            )
        )
    for name, variance in component_variances.items():
        fitted_variance += variance * float(np.mean(np.diag(raw_components[name])))
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
                    "REML" if reml else "ML",
                ),
                "contrast_transform": "gene-evolutionary-plus-sampling-covariance",
                "coverage_policy": coverage_policy,
                "small_sample_warning": "yes" if n_events < 20 else "no",
                "inference_status": inference_status,
                "model": model,
                "inference_method": inference,
                "reml": "yes" if reml else "no",
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
                }
            )
    if use_lineage:
        variance = lineage_variance
        for predictor, random_design in random_designs["lineage"]:
            modes = variance * random_design.T @ fit["inverse_residual"]
            for index, lineage_id in enumerate(unique_lineages):
                random_effect_rows.append(
                    {
                        "model_id": model_id,
                        "tree_id": tree_id,
                        "response": response,
                        "effect_type": "lineage_slope",
                        "group_id": lineage_id,
                        "term": predictor,
                        "conditional_mode": float(modes[index]),
                        "variance_component": variance,
                        "n_observations": int(lineage_counts[index]),
                    }
                )
    if return_fit_state:
        fit_state = {
            "beta": np.asarray(beta, dtype=float),
            "design": np.asarray(design, dtype=float),
            "contrast_ids": dataframe["gene_clade_id"].astype(str).tolist(),
            "evolutionary_rate": float(evolutionary_rate),
            "fitted_covariance": materialize_covariance(fit["covariance"]),
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
    for policy, label in [
        (event_random_effect, "event_random_effect"),
        (lineage_random_slope, "lineage_random_slope"),
    ]:
        if policy not in {"auto", "yes", "no"}:
            raise ValueError("{} must be auto, yes, or no.".format(label))
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
    covariance = np.asarray(normalized.get("covariance"), dtype=float)
    expected = (len(term_names), len(term_names), len(event_ids), len(event_ids))
    valid = (
        bool(term_names)
        and set(term_names) <= set(predictors)
        and len(event_ids) == len(set(event_ids))
        and covariance.shape == expected
        and np.isfinite(covariance).all()
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
            "covariance": covariance,
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
            else predictor_groups[str(row["source_term"])][0]
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
    fit_states,
    *,
    return_random_effects,
    return_fit_state,
):
    if return_random_effects and return_fit_state:
        return result, random_effects, fit_states
    if return_random_effects:
        return result, random_effects
    if return_fit_state:
        return result, fit_states
    return result


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
    return_random_effects=False,
    return_fit_state=False,
    predictor_metadata=None,
    predictor_groups=None,
    predictor_group_uncertainties=None,
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
                return_fit_state=return_fit_state,
                predictor_metadata=predictor_metadata,
                predictor_groups=predictor_groups,
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
    return _reconciled_return_value(
        result,
        random_effects,
        fit_states,
        return_random_effects=return_random_effects,
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


def _write_pgls_outputs(args, results, random_effects):
    _validate_pgls_file_output_paths(args)
    random_effects_path = getattr(args, "random_effects_out", None)
    file_outputs = []
    if random_effects_path is not None:
        file_outputs.append((random_effects_path, random_effects))
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
    if random_effects_path == "-":
        raise ValueError("'--random-effects-out' cannot be STDOUT.")
    output_paths = [args.outfile, random_effects_path]
    validate_distinct_output_paths(
        [
            ("--outfile", args.outfile),
            ("--random-effects-out", random_effects_path),
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
    results, random_effects = fit_reconciled_pgls(
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
        return_random_effects=True,
    )
    _warn_pgls_diagnostics(results)
    _write_pgls_outputs(args, results, random_effects)
