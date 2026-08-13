"""End-to-end reconciliation, contrast, and PGLS orchestration."""

import hashlib
import math
import os
import secrets
import shutil
import stat
import tempfile
from contextlib import ExitStack
from dataclasses import dataclass
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import norm, t

from nwkit.clade_index import CladeIndex
from nwkit.contrast import (
    _orient_children,
    _read_mixed_replicate_traits,
    _read_typed_traits,
    _validate_replicate_options,
    build_contrast_table,
    calculate_contrasts,
)
from nwkit.conventions import (
    DEFAULT_TABLE_MISSING_VALUES_CSV,
    pgls_bundle_lock_path,
    pgls_bundle_paths,
)
from nwkit.evolution import (
    build_evolutionary_covariance,
    build_sparse_evolutionary_model,
    evolution_model_spec,
    evolutionary_covariance_factory,
    optimization_parameterization,
    parameter_near_boundary,
    transformed_edge_variances,
)
from nwkit.gaussian import DiagonalLowRankCovariance, draw_from_factor
from nwkit.measurement_error import fit_latent_predictor, fit_sparse_latent_predictor
from nwkit.model_matrix import (
    PredictorTerm,
    ReplicatedObservation,
    encode_predictors,
    parse_key_values,
    parse_name_list,
    parse_numeric_key_values,
    parse_ordered_levels,
    resolve_response_specs,
    response_family_link,
    validate_response_auxiliaries,
)
from nwkit.multivariate_pgls import fit_multivariate_pgls
from nwkit.ordinary_pgls import (
    _categorical_shape_settings,
    _global_bounded_scalar_minimize,
    estimate_marginal_evolution_parameter,
)
from nwkit.pgls import (
    RANDOM_EFFECT_COLUMNS,
    RESPONSE_REQUIRED_COLUMNS,
    RESULT_COLUMNS,
    SENSITIVITY_COLUMNS,
    fit_reconciled_pgls,
)
from nwkit.phylogenetic_glmm import (
    SCALAR_RESPONSE_FAMILIES,
    fit_phylogenetic_glmm,
    summarize_glmm_coefficient,
    summarize_glmm_omnibus,
    summarize_glmm_threshold,
)
from nwkit.reconcile import _report_unmatched_species, build_reconciliation_table
from nwkit.rsc_diagnostics import (
    ORIGIN_DIAGNOSTIC_COLUMNS,
    build_categorical_origin_diagnostics,
)
from nwkit.sparse_laplace import (
    ContinuousPredictorUncertainty,
    GmrfPredictorUncertainty,
    GroupedPredictorUncertainty,
    JointPredictorUncertainty,
    SparseCovarianceModel,
    sparse_group_covariance,
)
from nwkit.species_parser import get_species_parser
from nwkit.util import (
    acquire_exclusive_lock,
    normalized_missing_path_key,
    read_tip_table,
    read_tree,
    read_tree_strings,
    validate_distinct_output_paths,
)


@dataclass
class PglsPipelineArtifacts:
    reconciliation: pd.DataFrame
    gene_contrasts: pd.DataFrame
    species_contrasts: pd.DataFrame
    response_sampling_covariance: pd.DataFrame | None
    response_tip_summary: pd.DataFrame | None
    results: pd.DataFrame
    random_effects: pd.DataFrame
    sensitivity: pd.DataFrame | None = None
    trait_origins: pd.DataFrame | None = None
    predictor_sampling_covariance: pd.DataFrame | None = None
    predictor_tip_summary: pd.DataFrame | None = None


@dataclass
class _GeneResponseInputs:
    values_by_trait: dict[str, dict[str, Any]]
    sampling_covariance_by_trait: (
        dict[str, pd.DataFrame | np.ndarray | DiagonalLowRankCovariance] | None
    )
    replicate_model_by_trait: dict[str, str] | None
    tip_summary: pd.DataFrame | None


@dataclass
class _SpeciesPredictorInputs:
    values_by_trait: dict[str, dict[str, Any]]
    sampling_covariance_by_trait: (
        dict[str, pd.DataFrame | np.ndarray | DiagonalLowRankCovariance] | None
    )
    replicate_model_by_trait: dict[str, str] | None
    tip_summary: pd.DataFrame | None
    sparse_posterior_by_trait: dict[str, SparseCovarianceModel] | None = None


def _ensemble_result_key_columns(results: pd.DataFrame) -> list[str]:
    candidates = [
        "response",
        "response_family",
        "response_level",
        "term",
        "source_term",
        "predictor_type",
        "predictor_level",
        "predictor_reference",
        "factor_coding",
        "term_test",
    ]
    return [column for column in candidates if column in results.columns]


def _combine_ensemble_results(
    frames: list[pd.DataFrame], base_tree_id: str, confidence_level: float
) -> pd.DataFrame:
    combined = pd.concat(frames, ignore_index=True)
    ensemble_size = len(frames)
    key_columns = _ensemble_result_key_columns(combined)
    rows = []
    for _key, group in combined.groupby(key_columns, sort=False, dropna=False):
        if group["tree_id"].duplicated().any():
            raise RuntimeError(
                "Tree-ensemble result keys are not unique within each tree."
            )
        row = group.iloc[0].to_dict()
        row["tree_id"] = base_tree_id
        row["model_id"] = "ensemble:{}:{}".format(base_tree_id, row.get("response", ""))
        row["ensemble_size"] = ensemble_size
        row["tree_support_fraction"] = group["tree_id"].nunique() / ensemble_size
        coefficients = pd.to_numeric(group["coefficient"], errors="coerce")
        standard_errors = pd.to_numeric(group["standard_error"], errors="coerce")
        valid = coefficients.notna() & standard_errors.notna()
        if valid.any():
            estimates = coefficients[valid].to_numpy(float)
            within_variance = float(
                np.mean(standard_errors[valid].to_numpy(float) ** 2)
            )
            between_variance = (
                0.0 if len(estimates) == 1 else float(np.var(estimates, ddof=1))
            )
            total_variance = (
                within_variance + (1.0 + 1.0 / len(estimates)) * between_variance
            )
            estimate = float(np.mean(estimates))
            standard_error = math.sqrt(max(total_variance, 0.0))
            statistic = "" if standard_error == 0.0 else estimate / standard_error
            extra_variance = (1.0 + 1.0 / len(estimates)) * between_variance
            degrees_of_freedom = (
                math.inf
                if len(estimates) <= 1 or extra_variance <= 0.0
                else (len(estimates) - 1)
                * (1.0 + within_variance / extra_variance) ** 2
            )
            critical = float(
                norm.ppf(0.5 + confidence_level / 2.0)
                if math.isinf(degrees_of_freedom)
                else t.ppf(0.5 + confidence_level / 2.0, degrees_of_freedom)
            )
            row.update(
                {
                    "coefficient": estimate,
                    "standard_error": standard_error,
                    "statistic": statistic,
                    "p_value": (
                        ""
                        if standard_error == 0.0
                        else float(
                            2.0
                            * (
                                norm.sf(abs(float(statistic)))
                                if math.isinf(degrees_of_freedom)
                                else t.sf(abs(float(statistic)), degrees_of_freedom)
                            )
                        )
                    ),
                    "degrees_of_freedom": degrees_of_freedom,
                    "confidence_interval_lower": estimate - critical * standard_error,
                    "confidence_interval_upper": estimate + critical * standard_error,
                    "between_tree_variance": between_variance,
                    "inference_method": "tree-ensemble-rubin",
                }
            )
        elif coefficients.notna().any():
            estimates = coefficients.dropna().to_numpy(float)
            row["coefficient"] = float(np.mean(estimates))
            row["standard_error"] = ""
            row["statistic"] = ""
            row["p_value"] = ""
            row["confidence_interval_lower"] = ""
            row["confidence_interval_upper"] = ""
            row["between_tree_variance"] = (
                0.0 if len(estimates) == 1 else float(np.var(estimates, ddof=1))
            )
            row["inference_method"] = "tree-ensemble-descriptive"
            row["inference_status"] = "no-within-tree-standard-error"
        else:
            row["statistic"] = ""
            row["p_value"] = ""
            row["between_tree_variance"] = ""
            row["inference_method"] = "tree-ensemble-unpooled"
            row["inference_status"] = "omnibus-requires-coefficient-covariance"
        for column in [
            "response_evolution_parameter",
            "predictor_evolution_parameter",
            "evolutionary_rate",
            "species_event_variance",
            "lineage_slope_variance",
        ]:
            numeric = pd.to_numeric(group[column], errors="coerce").dropna()
            if len(numeric):
                row[column] = float(numeric.mean())
        optimizer_values = set(group["optimizer_converged"].astype(str))
        row["optimizer_converged"] = "yes" if optimizer_values == {"yes"} else "no"
        row["boundary_warning"] = (
            "yes" if (group["boundary_warning"].astype(str) == "yes").any() else "no"
        )
        rows.append(row)
    return pd.DataFrame(rows).reindex(columns=RESULT_COLUMNS)


def _concat_optional_frames(frames):
    available = [frame for frame in frames if frame is not None and not frame.empty]
    return None if not available else pd.concat(available, ignore_index=True)


def _active_pgls_bundle_paths(
    prefix: str,
    artifacts: PglsPipelineArtifacts,
) -> dict[str, str]:
    paths = pgls_bundle_paths(prefix)
    inactive = set()
    if artifacts.response_sampling_covariance is None:
        inactive.add("response_sampling_covariance_out")
    elif artifacts.response_tip_summary is None:
        raise RuntimeError(
            "Replicate-aware covariance requires a response tip summary."
        )
    if artifacts.response_tip_summary is None:
        inactive.add("response_tip_summary_out")
    if artifacts.predictor_sampling_covariance is None:
        inactive.add("predictor_sampling_covariance_out")
    elif artifacts.predictor_tip_summary is None:
        raise RuntimeError(
            "Replicate-aware predictor covariance requires a predictor tip summary."
        )
    if artifacts.predictor_tip_summary is None:
        inactive.add("predictor_tip_summary_out")
    if artifacts.sensitivity is None or artifacts.sensitivity.empty:
        inactive.add("sensitivity_out")
    if artifacts.trait_origins is None or artifacts.trait_origins.empty:
        inactive.add("trait_origins_out")
    return {name: path for name, path in paths.items() if name not in inactive}


def _regular_output_mode(path: str) -> int | None:
    try:
        path_stat = os.lstat(path)
    except FileNotFoundError:
        return None
    if not stat.S_ISREG(path_stat.st_mode):
        raise ValueError(
            "Existing PGLS bundle target must be a regular file: '{}'.".format(path)
        )
    return stat.S_IMODE(path_stat.st_mode)


def _new_output_mode(directory: str) -> int:
    for _ in range(100):
        probe = os.path.join(
            directory,
            ".nwkit-pgls-mode-probe-{}".format(secrets.token_hex(16)),
        )
        try:
            descriptor = os.open(
                probe,
                os.O_CREAT | os.O_EXCL | os.O_WRONLY,
                0o666,
            )
        except FileExistsError:
            continue
        try:
            return stat.S_IMODE(os.fstat(descriptor).st_mode)
        finally:
            os.close(descriptor)
            os.remove(probe)
    raise FileExistsError("Could not allocate a PGLS output-mode probe.")


def _stage_dataframe(path: str, dataframe: pd.DataFrame, output_mode: int) -> str:
    absolute_path = os.path.abspath(path)
    directory = os.path.dirname(absolute_path)
    descriptor, staged_path = tempfile.mkstemp(
        prefix=".{}.stage.".format(os.path.basename(absolute_path)),
        dir=directory,
    )
    descriptor_open = True
    staged_stat = os.fstat(descriptor)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            descriptor_open = False
            dataframe.to_csv(handle, sep="\t", index=False)
            handle.flush()
            os.fsync(handle.fileno())
            if hasattr(os, "fchmod"):
                os.fchmod(handle.fileno(), output_mode)
        if not hasattr(os, "fchmod"):
            os.chmod(staged_path, output_mode)
        path_stat = os.lstat(staged_path)
        if (
            not stat.S_ISREG(path_stat.st_mode)
            or path_stat.st_dev != staged_stat.st_dev
            or path_stat.st_ino != staged_stat.st_ino
        ):
            raise RuntimeError("A PGLS staging file was replaced before commit.")
        return staged_path
    except BaseException:
        if descriptor_open:
            os.close(descriptor)
        if os.path.lexists(staged_path):
            os.remove(staged_path)
        raise


def _backup_regular_output(path: str) -> str:
    absolute_path = os.path.abspath(path)
    directory = os.path.dirname(absolute_path)
    descriptor, backup_path = tempfile.mkstemp(
        prefix=".{}.backup.".format(os.path.basename(absolute_path)),
        dir=directory,
    )
    descriptor_open = True
    flags = os.O_RDONLY
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    if hasattr(os, "O_NONBLOCK"):
        flags |= os.O_NONBLOCK
    source_descriptor = None
    try:
        source_descriptor = os.open(absolute_path, flags)
        source_stat = os.fstat(source_descriptor)
        if not stat.S_ISREG(source_stat.st_mode):
            raise ValueError(
                "Existing PGLS bundle target must be a regular file: '{}'.".format(path)
            )
        with os.fdopen(source_descriptor, "rb") as source_handle:
            source_descriptor = None
            with os.fdopen(descriptor, "wb") as backup_handle:
                descriptor_open = False
                shutil.copyfileobj(source_handle, backup_handle, length=1024 * 1024)
                backup_handle.flush()
                os.fsync(backup_handle.fileno())
                if hasattr(os, "fchmod"):
                    os.fchmod(backup_handle.fileno(), stat.S_IMODE(source_stat.st_mode))
        if not hasattr(os, "fchmod"):
            os.chmod(backup_path, stat.S_IMODE(source_stat.st_mode))
        return backup_path
    except BaseException:
        if source_descriptor is not None:
            os.close(source_descriptor)
        if descriptor_open:
            os.close(descriptor)
        if os.path.lexists(backup_path):
            os.remove(backup_path)
        raise


def _replace_output(source: str, target: str) -> None:
    os.replace(source, target)


def _restore_pgls_outputs(transactions: list[dict[str, Any]]) -> None:
    for transaction in reversed(transactions):
        target = transaction["target"]
        backup = transaction["backup"]
        if transaction["installed"]:
            if backup is None:
                if os.path.lexists(target):
                    os.remove(target)
            elif os.path.lexists(backup):
                _replace_output(backup, target)
        elif backup is not None and os.path.lexists(backup):
            os.remove(backup)


def _commit_pgls_outputs(staged_outputs: list[tuple[str, str]]) -> None:
    transactions: list[dict[str, Any]] = []
    commit_succeeded = False
    restoration_succeeded = False
    try:
        for target, staged_path in staged_outputs:
            transaction: dict[str, Any] = {
                "target": target,
                "staged_path": staged_path,
                "backup": None,
                "installed": False,
            }
            transactions.append(transaction)
            if os.path.lexists(target):
                transaction["backup"] = _backup_regular_output(target)
            _replace_output(staged_path, target)
            transaction["installed"] = True
        commit_succeeded = True
    except BaseException:
        for transaction in transactions:
            if not transaction["installed"] and not os.path.lexists(
                transaction["staged_path"]
            ):
                transaction["installed"] = True
        try:
            _restore_pgls_outputs(transactions)
            restoration_succeeded = True
        except BaseException as restore_exc:
            raise RuntimeError(
                "Failed to restore PGLS bundle outputs after a commit error; "
                "backup files were preserved."
            ) from restore_exc
        raise
    finally:
        if commit_succeeded or restoration_succeeded:
            for transaction in transactions:
                backup = transaction["backup"]
                if backup is not None and os.path.lexists(backup):
                    os.remove(backup)
        for _, staged_path in staged_outputs:
            if os.path.lexists(staged_path):
                os.remove(staged_path)


def _transaction_output_modes(outputs):
    output_modes: dict[str, int] = {}
    for path, _ in outputs:
        absolute_path = os.path.abspath(path)
        directory = os.path.dirname(absolute_path)
        mode = _regular_output_mode(absolute_path)
        output_modes[absolute_path] = (
            _new_output_mode(directory) if mode is None else mode
        )
    return output_modes


def _transaction_output_lock_path(absolute_path):
    directory = os.path.realpath(os.path.dirname(absolute_path))
    identity = hashlib.sha256(
        os.path.join(directory, os.path.basename(absolute_path)).encode("utf-8")
    ).hexdigest()
    return os.path.join(directory, ".nwkit-output-{}.lock".format(identity))


def _transaction_output_lock_paths(output_modes):
    return [
        _transaction_output_lock_path(absolute_path)
        for absolute_path in sorted(output_modes)
    ]


def _write_dataframes_transactionally(
    outputs: list[tuple[str, pd.DataFrame]],
) -> None:
    output_modes = _transaction_output_modes(outputs)
    lock_paths = _transaction_output_lock_paths(output_modes)
    with ExitStack() as locks:
        for lock_path in lock_paths:
            locks.enter_context(
                acquire_exclusive_lock(lock_path, lock_label="NWKIT output")
            )
        staged_outputs: list[tuple[str, str]] = []
        try:
            for path, dataframe in outputs:
                absolute_path = os.path.abspath(path)
                staged_outputs.append(
                    (
                        absolute_path,
                        _stage_dataframe(
                            absolute_path,
                            dataframe,
                            output_modes[absolute_path],
                        ),
                    )
                )
            _commit_pgls_outputs(staged_outputs)
        except BaseException:
            for _, staged_path in staged_outputs:
                if os.path.lexists(staged_path):
                    os.remove(staged_path)
            raise


def validate_pgls_bundle_target(
    prefix: str,
    protected_inputs: list[str | None] | None = None,
) -> dict[str, str]:
    """Validate deterministic bundle targets before an expensive model fit."""
    if not isinstance(prefix, str) or prefix.strip() in {"", "-"}:
        raise ValueError("'--out-prefix' must be a non-empty filesystem prefix.")
    paths = pgls_bundle_paths(prefix)
    lock_path = pgls_bundle_lock_path(prefix)
    protected_paths = {**paths, "transaction_lock": lock_path}
    validate_distinct_output_paths(
        [
            ("--out-prefix {}".format(name), path)
            for name, path in protected_paths.items()
        ]
    )
    for input_path in protected_inputs or []:
        if input_path in (None, "", "-"):
            continue
        for output_path in protected_paths.values():
            same_path = normalized_missing_path_key(
                input_path
            ) == normalized_missing_path_key(output_path)
            same_file = False
            if (
                not same_path
                and os.path.exists(input_path)
                and os.path.exists(output_path)
            ):
                try:
                    same_file = os.path.samefile(input_path, output_path)
                except OSError:
                    same_file = False
            if same_path or same_file:
                raise ValueError(
                    "PGLS bundle path must not overwrite an input file: '{}'.".format(
                        os.path.realpath(output_path)
                    )
                )
    return paths


def _effective_raw_args(args: Any) -> SimpleNamespace:
    values = vars(args).copy()
    defaults = {
        "batch": None,
        "biological_id": None,
        "event_source": "lca",
        "event_random_effect": "auto",
        "event_weighting": "equal",
        "gene_branch_length": "original",
        "gene_evolution_model": "brownian",
        "gene_evolution_parameter": None,
        "gene_tree_format": "auto",
        "missing_values": DEFAULT_TABLE_MISSING_VALUES_CSV,
        "model": "hierarchical",
        "quoted_node_names": True,
        "predictor_batch": None,
        "predictor_biological_id": None,
        "predictor_sample_size_columns": None,
        "predictor_standard_error_columns": None,
        "predictor_technical_aggregation": "error",
        "predictor_technical_id": None,
        "predictor_within_variance": "pooled",
        "categorical_predictors": None,
        "categorical_responses": None,
        "categorical_replicate_policy": "error",
        "factor_coding": "treatment",
        "factor_reference": None,
        "ordered_predictors": None,
        "ordered_responses": None,
        "response_reference": None,
        "response_family": None,
        "response_offset": None,
        "response_trials": None,
        "response_censor_lower": None,
        "response_censor_upper": None,
        "response_dispersion": None,
        "response_zero_probability": None,
        "coefficient_penalty": "student-t",
        "coefficient_prior_sd": 2.5,
        "multivariate_responses": False,
        "allow_missing_responses": False,
        "allow_large_dense": False,
        "sample_size_columns": None,
        "speciation_coverage": "complete",
        "species_branch_length": "original",
        "species_evolution_model": "brownian",
        "species_evolution_parameter": None,
        "species_tree_format": "auto",
        "standard_error_columns": None,
        "technical_aggregation": "error",
        "technical_id": None,
        "lineage_random_slope": "auto",
        "lineage_inference": "none",
        "lineage_leave_one_out": False,
        "categorical_origin_diagnostics": "none",
        "origin_map_replicates": 200,
        "origin_map_threads": 1,
        "origin_min_posterior": 0.5,
        "origin_leave_one_out": False,
        "unmatched": "error",
        "within_variance": "pooled",
    }
    for name, default in defaults.items():
        if values.get(name) is None:
            values[name] = default
    return SimpleNamespace(**values)


def _species_labels(gene_tree, args: Any) -> dict[str, str | None]:
    parser = get_species_parser(args=args)
    return {
        str(leaf.name): parser.parse(leaf.name).species_label
        for leaf in gene_tree.leaves()
    }


def _tree_clade_ids(tree) -> set[str]:
    clades = CladeIndex(tree)
    return {clades.clade_id_for_node(node) for node in tree.traverse()}


def _validate_matching_gene_topologies(gene_tree, reconciliation_tree) -> None:
    if _tree_clade_ids(gene_tree) != _tree_clade_ids(reconciliation_tree):
        raise ValueError(
            "'--gene-tree' and '--reconciliation-tree' must have identical rooted "
            "topologies and tip names; branch lengths and annotations may differ."
        )


def _read_gene_response_inputs(
    args: Any,
    gene_tree,
    responses: list[str],
    auxiliary_columns: list[str] | None = None,
) -> _GeneResponseInputs:
    categorical_responses = parse_name_list(args.categorical_responses)
    ordered_responses = parse_ordered_levels(
        args.ordered_responses, "--ordered-responses"
    )
    response_families = parse_key_values(args.response_family, "--response-family")
    auxiliary_columns = [] if auxiliary_columns is None else auxiliary_columns
    selected_columns = list(dict.fromkeys(responses + auxiliary_columns))
    replicate_requested = _validate_replicate_options(args)
    replicate_estimates = (
        _read_mixed_replicate_traits(
            args,
            gene_tree,
            selected_columns,
            set(categorical_responses) | set(ordered_responses),
            args.tree_id,
            categorical_policy="counts",
            non_gaussian_columns={
                response
                for response, family in response_families.items()
                if family in SCALAR_RESPONSE_FAMILIES
            },
            allow_missing_columns=(
                (
                    {response for response in responses}
                    if args.multivariate_responses and args.allow_missing_responses
                    else set()
                )
                | {
                    response
                    for response, family in response_families.items()
                    if family == "censored-gaussian"
                }
            ),
            option_name="--expression",
        )
        if replicate_requested
        else None
    )
    values_by_trait = (
        replicate_estimates.values_by_trait
        if replicate_estimates is not None
        else _read_typed_traits(
            args,
            gene_tree,
            selected_columns,
            categorical=categorical_responses,
            ordered=ordered_responses,
            allow_missing={
                response
                for response, family in response_families.items()
                if family == "censored-gaussian"
            }
            | (
                set(responses)
                if args.multivariate_responses and args.allow_missing_responses
                else set()
            ),
            option_name="--expression",
        )
    )
    tip_summary = (
        None if replicate_estimates is None else replicate_estimates.tip_summary
    )
    if tip_summary is not None:
        tip_summary = tip_summary[tip_summary["trait"].isin(responses)].copy()
    return _GeneResponseInputs(
        values_by_trait=values_by_trait,
        sampling_covariance_by_trait=(
            None
            if replicate_estimates is None
            else replicate_estimates.sampling_covariance_by_trait
        ),
        replicate_model_by_trait=(
            None if replicate_estimates is None else replicate_estimates.model_by_trait
        ),
        tip_summary=tip_summary,
    )


def _read_gene_auxiliary_values(args, gene_tree, columns, *, allow_missing=False):
    if not columns:
        return {}
    dataframe, _, _ = read_tip_table(
        args.expression,
        option_name="--expression",
        tree_leaf_names=list(gene_tree.leaf_names()),
        required_columns=columns,
        unmatched=args.unmatched,
        missing_values=args.missing_values,
        duplicate_leaf_names="allow",
    )
    values_by_column = {}
    for column in columns:
        values_by_leaf = {}
        for leaf_name in gene_tree.leaf_names():
            selected_rows = dataframe.loc[dataframe["leaf_name"] == str(leaf_name)]
            biological_id = getattr(args, "biological_id", None)
            if biological_id is None:
                numeric = pd.to_numeric(
                    selected_rows[column], errors="coerce"
                ).to_numpy(float)
            else:
                observations = []
                for _identifier, group in selected_rows.groupby(
                    biological_id, sort=False, dropna=False
                ):
                    group_values = pd.to_numeric(
                        group[column], errors="coerce"
                    ).to_numpy(float)
                    finite_group = group_values[np.isfinite(group_values)]
                    if len(finite_group) and np.any(finite_group != finite_group[0]):
                        raise ValueError(
                            "Auxiliary response column '{}' differs among technical "
                            "rows for gene tip '{}'.".format(column, leaf_name)
                        )
                    observations.append(
                        float("nan")
                        if not len(finite_group)
                        else float(finite_group[0])
                    )
                numeric = np.asarray(observations, dtype=float)
            finite = numeric[np.isfinite(numeric)]
            if not allow_missing and (not len(numeric) or len(finite) != len(numeric)):
                raise ValueError(
                    "Auxiliary response column '{}' must be finite for gene tip '{}'.".format(
                        column, leaf_name
                    )
                )
            if biological_id is None and len(finite) and np.any(finite != finite[0]):
                raise ValueError(
                    "Auxiliary response column '{}' differs among rows for gene tip '{}'.".format(
                        column, leaf_name
                    )
                )
            values_by_leaf[str(leaf_name)] = (
                ReplicatedObservation(tuple(numeric))
                if biological_id is not None
                else float("nan")
                if not len(finite)
                else float(finite[0])
            )
        values_by_column[column] = values_by_leaf
    return values_by_column


def _map_response_auxiliary(column_by_response, values_by_column):
    return {
        response: values_by_column[column]
        for response, column in column_by_response.items()
    }


def _read_species_predictor_inputs(
    args: Any,
    species_tree,
    predictors: list[str],
) -> _SpeciesPredictorInputs:
    categorical_predictors = parse_name_list(args.categorical_predictors)
    ordered_predictors = parse_ordered_levels(args.ordered_predictors)
    values = vars(args).copy()
    values.update(
        {
            "batch": args.predictor_batch,
            "biological_id": args.predictor_biological_id,
            "sample_size_columns": args.predictor_sample_size_columns,
            "sampling_covariance_out": None,
            "standard_error_columns": args.predictor_standard_error_columns,
            "technical_aggregation": args.predictor_technical_aggregation,
            "technical_id": args.predictor_technical_id,
            "tip_summary_out": None,
            "within_variance": args.predictor_within_variance,
        }
    )
    predictor_args = SimpleNamespace(**values)
    replicate_requested = _validate_replicate_options(predictor_args)
    discrete_predictors = set(categorical_predictors) | set(ordered_predictors)
    replicate_estimates = (
        _read_mixed_replicate_traits(
            predictor_args,
            species_tree,
            predictors,
            discrete_predictors,
            "species",
            categorical_policy=args.categorical_replicate_policy,
            option_name="--species-traits",
        )
        if replicate_requested
        else None
    )
    values_by_trait = (
        replicate_estimates.values_by_trait
        if replicate_estimates is not None
        else _read_typed_traits(
            predictor_args,
            species_tree,
            predictors,
            categorical=categorical_predictors,
            ordered=ordered_predictors,
            option_name="--species-traits",
        )
    )
    return _SpeciesPredictorInputs(
        values_by_trait=values_by_trait,
        sampling_covariance_by_trait=(
            None
            if replicate_estimates is None
            else replicate_estimates.sampling_covariance_by_trait
        ),
        replicate_model_by_trait=(
            None if replicate_estimates is None else replicate_estimates.model_by_trait
        ),
        tip_summary=(
            None if replicate_estimates is None else replicate_estimates.tip_summary
        ),
    )


def _build_gene_response_contrasts(
    args: Any,
    gene_tree,
    reconciliation_by_id: dict[str, dict[str, Any]],
    response_inputs: _GeneResponseInputs,
    response: str,
    evolution_parameter: float | None,
) -> tuple[pd.DataFrame, pd.DataFrame | None]:
    sampling_by_trait = response_inputs.sampling_covariance_by_trait
    replicate_models = response_inputs.replicate_model_by_trait
    response_tip_summary = (
        None
        if response_inputs.tip_summary is None
        else response_inputs.tip_summary[
            response_inputs.tip_summary["trait"] == response
        ].copy()
    )
    output = build_contrast_table(
        gene_tree,
        {response: response_inputs.values_by_trait[response]},
        branch_length=args.gene_branch_length,
        evolution_model=args.gene_evolution_model,
        evolution_parameter=evolution_parameter,
        reconciliation_by_id=reconciliation_by_id,
        event_type="speciation",
        eligible_only=False,
        speciation_coverage=args.speciation_coverage,
        tree_id=args.tree_id,
        sampling_covariance_by_trait=(
            None
            if sampling_by_trait is None
            else {response: sampling_by_trait[response]}
        ),
        replicate_model_by_trait=(
            None if replicate_models is None else {response: replicate_models[response]}
        ),
        tip_summary=response_tip_summary,
        return_sampling_covariance=sampling_by_trait is not None,
        tree_option_name="--gene-tree",
    )
    if sampling_by_trait is None:
        return output, None
    contrasts, covariance = output
    return contrasts, covariance


def _fixed_evolution_diagnostics(model: str, parameter: float | None) -> dict[str, Any]:
    parameterless = evolution_model_spec(model).parameter_name is None
    status = "not-applicable" if parameterless else "fixed"
    return {
        "parameter": parameter,
        "parameter_status": status,
        "log_likelihood": "",
        "optimizer_converged": None if parameterless else True,
        "optimizer_message": (
            "model has no shape parameter"
            if parameterless
            else "fixed evolutionary transform"
        ),
        "boundary_warning": None if parameterless else False,
    }


def _fit_candidate_reconciled_model(
    args: Any,
    gene_contrasts: pd.DataFrame,
    species_contrasts: pd.DataFrame,
    response: str,
    predictors: list[str],
    sampling_covariance: pd.DataFrame | None,
    predictor_sampling_covariance: pd.DataFrame | None,
    predictor_group_uncertainties,
) -> pd.DataFrame:
    return fit_reconciled_pgls(
        gene_contrasts,
        species_contrasts,
        [response],
        predictors,
        confidence_level=args.confidence_level,
        event_weighting=args.event_weighting,
        coverage_policy=args.speciation_coverage,
        model=args.model,
        response_sampling_covariance=sampling_covariance,
        predictor_sampling_covariance=predictor_sampling_covariance,
        predictor_group_uncertainties=predictor_group_uncertainties,
        inference="wald",
        bootstrap_replicates=2,
        seed=args.seed,
        reml=args.reml,
        event_random_effect=args.event_random_effect,
        lineage_random_slope=args.lineage_random_slope,
        lineage_inference="none",
        allow_large_dense=args.allow_large_dense,
    )


def _estimate_gene_response_parameter(
    args: Any,
    gene_tree,
    reconciliation_by_id: dict[str, dict[str, Any]],
    response_inputs: _GeneResponseInputs,
    response: str,
    species_contrasts: pd.DataFrame,
    predictor_sampling_covariance: pd.DataFrame | None,
    predictors: list[str],
    predictor_group_uncertainties,
) -> tuple[pd.DataFrame, pd.DataFrame | None, dict[str, Any]]:
    bounds, decode = optimization_parameterization(
        gene_tree,
        args.gene_evolution_model,
        branch_length=args.gene_branch_length,
    )
    cache: dict[float, dict[str, Any] | None] = {}

    def cached_candidate(encoded_value: float) -> dict[str, Any] | None:
        parameter = float(decode(float(encoded_value)))
        if parameter not in cache:
            try:
                contrasts, covariance = _build_gene_response_contrasts(
                    args,
                    gene_tree,
                    reconciliation_by_id,
                    response_inputs,
                    response,
                    parameter,
                )
                result = _fit_candidate_reconciled_model(
                    args,
                    contrasts,
                    species_contrasts,
                    response,
                    predictors,
                    covariance,
                    predictor_sampling_covariance,
                    predictor_group_uncertainties,
                )
                likelihoods = pd.to_numeric(
                    result["log_likelihood"], errors="coerce"
                ).to_numpy(dtype=float)
                if not np.isfinite(likelihoods).all() or not np.allclose(
                    likelihoods, likelihoods[0]
                ):
                    raise ValueError(
                        "Candidate reconciled fit returned an inconsistent likelihood."
                    )
                cache[parameter] = {
                    "parameter": parameter,
                    "contrasts": contrasts,
                    "covariance": covariance,
                    "result": result,
                    "objective": -float(likelihoods[0]),
                }
            except ValueError:
                cache[parameter] = None
        return cache[parameter]

    def objective(encoded_value: float) -> float:
        candidate = cached_candidate(encoded_value)
        return float("inf") if candidate is None else float(candidate["objective"])

    optimized = _global_bounded_scalar_minimize(objective, bounds)
    encoded_candidates = [bounds[0], bounds[1]]
    if math.isfinite(float(optimized.fun)):
        encoded_candidates.append(float(optimized.x))
    candidates = [cached_candidate(value) for value in encoded_candidates]
    finite_candidates = [candidate for candidate in candidates if candidate is not None]
    if not finite_candidates:
        raise ValueError(
            "Gene evolution-parameter optimization failed to find a finite "
            "reconciled PGLS fit for response '{}'.".format(response)
        )
    selected = min(
        finite_candidates, key=lambda candidate: float(candidate["objective"])
    )
    parameter = float(selected["parameter"])
    inner_converged = selected["result"]["optimizer_converged"].eq("yes").all()
    diagnostics = {
        "parameter": parameter,
        "parameter_status": "estimated",
        "log_likelihood": -float(selected["objective"]),
        "optimizer_converged": bool(optimized.success and inner_converged),
        "optimizer_message": str(optimized.message),
        "boundary_warning": parameter_near_boundary(
            gene_tree,
            args.gene_evolution_model,
            parameter,
            branch_length=args.gene_branch_length,
        ),
    }
    return selected["contrasts"], selected["covariance"], diagnostics


def _build_gene_contrasts(
    args: Any,
    gene_tree,
    reconciliation: pd.DataFrame,
    response_inputs: _GeneResponseInputs,
    responses: list[str],
    species_contrasts: pd.DataFrame,
    predictor_sampling_covariance: pd.DataFrame | None,
    predictors: list[str],
    predictor_group_uncertainties=(),
) -> tuple[
    pd.DataFrame,
    pd.DataFrame | None,
    dict[str, dict[str, Any]],
]:
    reconciliation_by_id = {
        str(row["gene_clade_id"]): row for row in reconciliation.to_dict("records")
    }
    spec = evolution_model_spec(args.gene_evolution_model)
    auto_parameter = (
        spec.parameter_name is not None
        and args.gene_evolution_parameter
        in {
            None,
            "auto",
        }
    )
    if auto_parameter and args.model == "legacy":
        raise ValueError(
            "Automatic gene evolution-parameter estimation requires a likelihood-based "
            "reconciled model, not '--model legacy'."
        )
    contrast_frames = []
    covariance_frames = []
    diagnostics_by_response = {}
    for response in responses:
        if auto_parameter:
            contrasts, covariance, diagnostics = _estimate_gene_response_parameter(
                args,
                gene_tree,
                reconciliation_by_id,
                response_inputs,
                response,
                species_contrasts,
                predictor_sampling_covariance,
                predictors,
                predictor_group_uncertainties,
            )
        else:
            parameter = (
                None
                if spec.parameter_name is None
                else float(args.gene_evolution_parameter)
            )
            contrasts, covariance = _build_gene_response_contrasts(
                args,
                gene_tree,
                reconciliation_by_id,
                response_inputs,
                response,
                parameter,
            )
            diagnostics = _fixed_evolution_diagnostics(
                args.gene_evolution_model, parameter
            )
        contrast_frames.append(contrasts)
        if covariance is not None:
            covariance_frames.append(covariance)
        diagnostics_by_response[response] = diagnostics
    covariance_table = (
        None
        if response_inputs.sampling_covariance_by_trait is None
        else pd.concat(covariance_frames, ignore_index=True)
    )
    return (
        pd.concat(contrast_frames, ignore_index=True),
        covariance_table,
        diagnostics_by_response,
    )


def _build_species_contrasts(
    args: Any,
    species_tree,
    predictor_inputs: _SpeciesPredictorInputs,
    predictors: list[str],
    predictor_groups: dict[str, tuple[str, ...]] | None = None,
    predictor_group_uncertainties: dict[str, np.ndarray] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame | None, dict[str, dict[str, Any]]]:
    spec = evolution_model_spec(args.species_evolution_model)
    auto_parameter = (
        spec.parameter_name is not None
        and args.species_evolution_parameter in {None, "auto"}
    )
    frames = []
    covariance_frames = []
    diagnostics_by_predictor = {}
    sampling_by_trait = predictor_inputs.sampling_covariance_by_trait
    replicate_models = predictor_inputs.replicate_model_by_trait
    predictor_groups = (
        {predictor: (predictor,) for predictor in predictors}
        if predictor_groups is None
        else predictor_groups
    )
    predictor_group_uncertainties = predictor_group_uncertainties or {}
    shared_diagnostics = {}
    if auto_parameter:
        for source, terms in predictor_groups.items():
            if len(terms) == 1:
                continue
            shared_diagnostics[source] = _estimate_factor_evolution_parameter(
                species_tree,
                predictor_inputs.values_by_trait,
                terms,
                evolution_model=args.species_evolution_model,
                branch_length=args.species_branch_length,
                sampling_covariance_by_observation=predictor_group_uncertainties.get(
                    source
                ),
            )
    source_by_term = {
        term: source for source, terms in predictor_groups.items() for term in terms
    }
    for predictor in predictors:
        predictor_sampling = (
            None
            if not sampling_by_trait or predictor not in sampling_by_trait
            else sampling_by_trait[predictor]
        )
        parameter: float | None
        source = source_by_term[predictor]
        if source in shared_diagnostics:
            diagnostics = shared_diagnostics[source]
            parameter = float(diagnostics["parameter"])
        elif auto_parameter:
            diagnostics = estimate_marginal_evolution_parameter(
                species_tree,
                predictor_inputs.values_by_trait[predictor],
                predictor,
                evolution_model=args.species_evolution_model,
                branch_length=args.species_branch_length,
                sampling_covariance=predictor_sampling,
                allow_large_dense=args.allow_large_dense,
            )
            parameter = float(diagnostics["parameter"])
        else:
            parameter = (
                None
                if spec.parameter_name is None
                else float(args.species_evolution_parameter)
            )
            diagnostics = _fixed_evolution_diagnostics(
                args.species_evolution_model, parameter
            )
        output = build_contrast_table(
            species_tree,
            {predictor: predictor_inputs.values_by_trait[predictor]},
            branch_length=args.species_branch_length,
            evolution_model=args.species_evolution_model,
            evolution_parameter=parameter,
            tree_id="species",
            sampling_covariance_by_trait=(
                None if predictor_sampling is None else {predictor: predictor_sampling}
            ),
            replicate_model_by_trait=(
                None
                if predictor_sampling is None or replicate_models is None
                else {predictor: replicate_models[predictor]}
            ),
            tip_summary=(
                None
                if predictor_inputs.tip_summary is None
                else predictor_inputs.tip_summary[
                    predictor_inputs.tip_summary["trait"] == predictor
                ].copy()
            ),
            return_sampling_covariance=predictor_sampling is not None,
            tree_option_name="--species-tree",
        )
        if predictor_sampling is None:
            frames.append(output)
        else:
            contrasts, covariance = output
            frames.append(contrasts)
            covariance_frames.append(covariance)
        diagnostics_by_predictor[predictor] = diagnostics
    covariance_table = (
        None
        if not covariance_frames
        else pd.concat(covariance_frames, ignore_index=True)
    )
    return (
        pd.concat(frames, ignore_index=True),
        covariance_table,
        diagnostics_by_predictor,
    )


def _prepare_reconciled_tip_predictors(
    args: Any,
    species_tree,
    raw_inputs: _SpeciesPredictorInputs,
    predictors: list[str],
    categorical_predictors: list[str],
    ordered_predictors: dict[str, tuple[str, ...]],
    factor_references: dict[str, str],
    predictor_diagnostics: dict[str, dict[str, Any]],
):
    """Condition continuous replicated predictors before gene-tip PGLMM use."""
    leaf_names = [str(leaf.name) for leaf in species_tree.leaves()]
    values_by_trait = {
        trait: dict(values) for trait, values in raw_inputs.values_by_trait.items()
    }
    posterior_covariances: dict[str, pd.DataFrame] = {}
    sparse_posteriors: dict[str, SparseCovarianceModel] = {}
    diagnostics = {
        trait: dict(values) for trait, values in predictor_diagnostics.items()
    }
    sampling_by_trait = raw_inputs.sampling_covariance_by_trait or {}
    for predictor in predictors:
        if predictor in categorical_predictors or predictor not in sampling_by_trait:
            continue
        observed = np.asarray(
            [float(values_by_trait[predictor][name]) for name in leaf_names],
            dtype=float,
        )
        sampling_value = sampling_by_trait[predictor]
        if isinstance(sampling_value, DiagonalLowRankCovariance):
            sampling = sampling_value
        else:
            sampling = (
                sampling_value.loc[leaf_names, leaf_names].to_numpy(dtype=float)
                if isinstance(sampling_value, pd.DataFrame)
                else np.asarray(sampling_value, dtype=float)
            )
        predictor_diagnostic = diagnostics[predictor]
        parameter = predictor_diagnostic.get("parameter")
        sampling_diagonal = (
            sampling.diagonal
            if isinstance(sampling, DiagonalLowRankCovariance)
            else sampling
            if sampling.ndim == 1
            else np.diag(sampling)
        )
        if np.all(sampling_diagonal > 0.0) and (
            isinstance(sampling, DiagonalLowRankCovariance)
            or sampling.ndim == 1
            or np.array_equal(sampling, np.diag(sampling_diagonal))
        ):
            sparse_prior = build_sparse_evolutionary_model(
                species_tree,
                leaf_names,
                model=args.species_evolution_model,
                parameter=None if parameter in {None, ""} else float(parameter),
                branch_length=args.species_branch_length,
            )
            posterior = fit_sparse_latent_predictor(
                observed,
                sparse_prior,
                sampling,
                include_intercept=True,
            )
            values_by_trait[predictor] = dict(
                zip(leaf_names, posterior.mean, strict=True)
            )
            sparse_posteriors[predictor] = posterior.covariance_model
            mean_posterior_variance = posterior.mean_posterior_variance
            mean_prior_variance = (
                posterior.evolutionary_rate * sparse_prior.covariance_scale
            )
        else:
            prior_covariance = build_evolutionary_covariance(
                species_tree,
                leaf_names,
                model=args.species_evolution_model,
                parameter=None if parameter in {None, ""} else float(parameter),
                branch_length=args.species_branch_length,
            )
            posterior = fit_latent_predictor(
                observed,
                prior_covariance,
                sampling,
                include_intercept=True,
            )
            values_by_trait[predictor] = dict(
                zip(leaf_names, posterior.mean, strict=True)
            )
            posterior_covariances[predictor] = pd.DataFrame(
                posterior.covariance,
                index=leaf_names,
                columns=leaf_names,
            )
            mean_posterior_variance = float(np.mean(np.diag(posterior.covariance)))
            mean_prior_variance = float(
                np.mean(np.diag(posterior.evolutionary_rate * prior_covariance))
            )
        diagnostics[predictor].update(
            {
                "optimizer_converged": bool(
                    predictor_diagnostic.get("optimizer_converged", True)
                    and posterior.optimizer_converged
                ),
                "optimizer_message": "{}; {}".format(
                    predictor_diagnostic.get(
                        "optimizer_message", "fixed evolutionary model"
                    ),
                    posterior.optimizer_message,
                ),
                "boundary_warning": bool(
                    predictor_diagnostic.get("boundary_warning", False)
                    or posterior.boundary_warning
                ),
                "log_likelihood": posterior.log_likelihood,
                "evolutionary_rate": posterior.evolutionary_rate,
                "rate_optimizer_converged": posterior.optimizer_converged,
                "rate_optimizer_message": posterior.optimizer_message,
                "rate_boundary_warning": posterior.boundary_warning,
                "mean_sampling_variance": float(np.mean(sampling_diagonal)),
                "mean_posterior_variance": mean_posterior_variance,
                "uncertainty_fraction": float(
                    mean_posterior_variance
                    / max(mean_prior_variance, np.finfo(float).tiny)
                ),
            }
        )
    encoded = encode_predictors(
        values_by_trait,
        predictors,
        leaf_names,
        categorical=categorical_predictors,
        ordered_levels=ordered_predictors,
        factor_references=factor_references,
        factor_coding=args.factor_coding,
    )
    inputs = _SpeciesPredictorInputs(
        values_by_trait=encoded.values_by_trait,
        sampling_covariance_by_trait=posterior_covariances or None,
        replicate_model_by_trait=raw_inputs.replicate_model_by_trait,
        tip_summary=raw_inputs.tip_summary,
        sparse_posterior_by_trait=sparse_posteriors or None,
    )
    return encoded, inputs, diagnostics


def _positive_semidefinite_loading(matrix, dimension):
    scale = max(1.0, float(np.max(np.abs(matrix))))
    tolerance = np.finfo(float).eps * scale * max(1, dimension) * 100.0
    if float(np.max(np.abs(matrix - matrix.T))) > tolerance:
        raise ValueError("Categorical predictor sampling covariance must be symmetric.")
    eigenvalues, eigenvectors = np.linalg.eigh((matrix + matrix.T) / 2.0)
    if float(np.min(eigenvalues)) < -tolerance:
        raise ValueError(
            "Categorical predictor sampling covariance must be positive semidefinite."
        )
    positive = eigenvalues > tolerance
    return eigenvectors[:, positive] * np.sqrt(eigenvalues[positive])


def _categorical_sampling_factor(sampling_covariance, n_tips, n_terms):
    sampling = np.asarray(sampling_covariance, dtype=float)
    expected = (n_tips, n_terms, n_terms)
    if sampling.shape != expected or not np.isfinite(sampling).all():
        raise ValueError(
            "Categorical predictor sampling covariance has invalid dimensions or values."
        )
    loading_rows: list[int] = []
    loading_columns: list[int] = []
    loading_values: list[float] = []
    latent_offset = 0
    for tip_index, matrix in enumerate(sampling):
        block = _positive_semidefinite_loading(matrix, n_terms)
        block_rows, block_columns = np.nonzero(block)
        loading_rows.extend(int(row) * n_tips + tip_index for row in block_rows)
        loading_columns.extend(latent_offset + int(column) for column in block_columns)
        loading_values.extend(
            float(block[row, column])
            for row, column in zip(block_rows, block_columns, strict=True)
        )
        latent_offset += block.shape[1]
    loading = sparse.csr_matrix(
        (loading_values, (loading_rows, loading_columns)),
        shape=(n_tips * n_terms, latent_offset),
    )
    return DiagonalLowRankCovariance(np.zeros(n_tips * n_terms, dtype=float), loading)


def _factor_evolution_diagnostics(tree, model, parameter, branch_length, fit):
    return {
        "parameter": parameter,
        "parameter_status": "estimated",
        "log_likelihood": fit.log_likelihood,
        "optimizer_converged": fit.optimizer_converged,
        "optimizer_message": fit.optimizer_message,
        "boundary_warning": parameter_near_boundary(
            tree, model, parameter, branch_length=branch_length
        ),
    }


def _estimate_factor_evolution_parameter(
    tree,
    values_by_trait,
    terms,
    *,
    evolution_model,
    branch_length,
    sampling_covariance_by_observation=None,
):
    """Estimate one coding-invariant shape parameter for a factor's columns."""
    leaf_names = [str(leaf.name) for leaf in tree.leaves()]
    response = np.column_stack(
        [[float(values_by_trait[term][name]) for name in leaf_names] for term in terms]
    )
    bounds, decode = optimization_parameterization(
        tree, evolution_model, branch_length=branch_length
    )
    fixed_covariance = (
        None
        if sampling_covariance_by_observation is None
        else _categorical_sampling_factor(
            sampling_covariance_by_observation, len(leaf_names), len(terms)
        )
    )
    factory = evolutionary_covariance_factory(
        tree,
        leaf_names,
        model=evolution_model,
        branch_length=branch_length,
    )
    encoded_bounds = (float(bounds[0]), float(bounds[1]))
    fit = fit_multivariate_pgls(
        response,
        np.ones((len(leaf_names), 1), dtype=float),
        {"phylogenetic": factory},
        fixed_covariance=fixed_covariance,
        evolution_parameter_bounds=encoded_bounds,
        evolution_parameter_decoder=decode,
        evolution_parameter_initial=sum(encoded_bounds) / 2.0,
        reml=False,
    )
    if fit.evolution_parameter is None:
        raise RuntimeError(
            "Factor evolutionary shape estimation returned no parameter."
        )
    parameter = float(fit.evolution_parameter)
    return _factor_evolution_diagnostics(
        tree, evolution_model, parameter, branch_length, fit
    )


def _categorical_contrast_uncertainties(
    species_tree,
    species_contrasts,
    encoded_predictors,
    predictor_diagnostics,
    *,
    evolution_model,
    branch_length,
):
    clades = CladeIndex(species_tree)
    node_by_id = {
        clades.clade_id_for_node(node): node for node in species_tree.traverse()
    }
    states = []
    audit_rows = []
    for uncertainty in encoded_predictors.uncertainties:
        term_names = tuple(uncertainty.term_names)
        event_ids = (
            species_contrasts[species_contrasts["trait"] == term_names[0]][
                "branch_clade_id"
            ]
            .astype(str)
            .tolist()
        )
        coefficient_matrices = []
        for term_name in term_names:
            diagnostics = predictor_diagnostics[term_name]
            parameter = diagnostics["parameter"]
            _, coefficients_by_node, leaf_names = calculate_contrasts(
                species_tree,
                encoded_predictors.values_by_trait[term_name],
                branch_length=branch_length,
                evolution_model=evolution_model,
                evolution_parameter=parameter,
                return_coefficients=True,
                sparse_coefficients=True,
            )
            if leaf_names != [str(leaf.name) for leaf in species_tree.leaves()]:
                raise RuntimeError("Categorical contrast leaf ordering changed.")
            coefficient_matrices.append(
                sparse.vstack(
                    [
                        coefficients_by_node[node_by_id[event_id]]
                        for event_id in event_ids
                    ],
                    format="csr",
                )
            )
        event_count = len(event_ids)
        dimension = len(term_names)
        local_factors = []
        for matrix in uncertainty.covariance_by_observation:
            symmetric = (matrix + matrix.T) / 2.0
            eigenvalues, eigenvectors = np.linalg.eigh(symmetric)
            tolerance = (
                np.finfo(float).eps
                * max(1.0, float(np.max(np.abs(symmetric))))
                * max(1, dimension)
                * 100.0
            )
            if float(eigenvalues.min()) < -tolerance:
                raise ValueError(
                    "Categorical predictor sampling covariance is not positive "
                    "semidefinite."
                )
            positive = eigenvalues > tolerance
            local_factors.append(
                eigenvectors[:, positive] * np.sqrt(eigenvalues[positive])
            )
        latent_factor = sparse.block_diag(local_factors, format="csr")
        factor_by_term = tuple(
            (
                coefficient_matrices[term_index]
                @ latent_factor[np.arange(len(local_factors)) * dimension + term_index]
            ).tocsr()
            for term_index in range(dimension)
        )
        joint_uncertainty = JointPredictorUncertainty(factors=factor_by_term)
        states.append(
            {
                "source": uncertainty.source,
                "term_names": term_names,
                "event_ids": event_ids,
                "event_index": {
                    event_id: index for index, event_id in enumerate(event_ids)
                },
                "uncertainty": joint_uncertainty,
            }
        )
        full_audit = event_count * dimension <= 500
        for first, first_term in enumerate(term_names):
            for second, second_term in enumerate(term_names):
                if full_audit:
                    cross_covariance = (
                        factor_by_term[first] @ factor_by_term[second].T
                    ).tocsr()
                    marginal_covariance = np.empty(0, dtype=float)
                else:
                    cross_covariance = None
                    marginal_covariance = np.asarray(
                        factor_by_term[first]
                        .multiply(factor_by_term[second])
                        .sum(axis=1)
                    ).reshape(-1)
                for first_event, first_id in enumerate(event_ids):
                    second_events = (
                        range(first_event, event_count)
                        if full_audit
                        else range(first_event, first_event + 1)
                    )
                    for second_event in second_events:
                        audit_rows.append(
                            {
                                "tree_id": "species",
                                "trait": first_term,
                                "trait_2": second_term,
                                "contrast_id_1": first_id,
                                "contrast_id_2": event_ids[second_event],
                                "sampling_covariance": float(
                                    cross_covariance[first_event, second_event]
                                    if cross_covariance is not None
                                    else marginal_covariance[first_event]
                                ),
                                "audit_scope": (
                                    "full" if full_audit else "marginal-only"
                                ),
                            }
                        )
    return states, pd.DataFrame(audit_rows)


def _yes_no(value: bool) -> str:
    return "yes" if value else "no"


def _diagnostic_flag(value: bool | None) -> str:
    return "not-applicable" if value is None else _yes_no(value)


def _attach_evolution_diagnostics(
    results: pd.DataFrame,
    response_diagnostics: dict[str, dict[str, Any]],
    predictor_diagnostics: dict[str, dict[str, Any]],
    predictor_groups: dict[str, tuple[str, ...]] | None = None,
) -> pd.DataFrame:
    updated = results.copy()
    diagnostic_columns = [
        "response_evolution_parameter_status",
        "response_evolution_optimizer_converged",
        "response_evolution_optimizer_message",
        "response_evolution_boundary_warning",
        "response_evolution_parameter_bootstrap_refit",
        "predictor_evolution_parameter_status",
        "predictor_evolution_optimizer_converged",
        "predictor_evolution_optimizer_message",
        "predictor_evolution_boundary_warning",
        "predictor_evolution_log_likelihood",
    ]
    for column in diagnostic_columns:
        updated[column] = updated[column].astype(object)
    for row_index, row in updated.iterrows():
        response = response_diagnostics[str(row["response"])]
        predictor_term = str(row["term"])
        if row.get("term_test", "coefficient") != "coefficient":
            if predictor_groups is None:
                raise ValueError("Grouped predictor diagnostics require term groups.")
            predictor_term = predictor_groups[str(row["source_term"])][0]
        predictor = predictor_diagnostics[predictor_term]
        updated.loc[row_index, "response_evolution_parameter_status"] = response[
            "parameter_status"
        ]
        updated.loc[row_index, "response_evolution_optimizer_converged"] = (
            _diagnostic_flag(response["optimizer_converged"])
        )
        updated.loc[row_index, "response_evolution_optimizer_message"] = response[
            "optimizer_message"
        ]
        updated.loc[row_index, "response_evolution_boundary_warning"] = (
            _diagnostic_flag(response["boundary_warning"])
        )
        updated.loc[row_index, "response_evolution_parameter_bootstrap_refit"] = "no"
        updated.loc[row_index, "predictor_evolution_parameter_status"] = predictor[
            "parameter_status"
        ]
        updated.loc[row_index, "predictor_evolution_optimizer_converged"] = (
            _diagnostic_flag(predictor["optimizer_converged"])
        )
        updated.loc[row_index, "predictor_evolution_optimizer_message"] = predictor[
            "optimizer_message"
        ]
        updated.loc[row_index, "predictor_evolution_boundary_warning"] = (
            _diagnostic_flag(predictor["boundary_warning"])
        )
        updated.loc[row_index, "predictor_evolution_log_likelihood"] = predictor[
            "log_likelihood"
        ]
        if (
            response["optimizer_converged"] is False
            or predictor["optimizer_converged"] is False
        ):
            updated.loc[row_index, "optimizer_converged"] = "no"
        if (
            response["boundary_warning"] is True
            or predictor["boundary_warning"] is True
        ):
            updated.loc[row_index, "boundary_warning"] = "yes"
    return updated


def _positive_semidefinite_factor(matrix: np.ndarray, label: str) -> np.ndarray:
    symmetric = (
        np.asarray(matrix, dtype=float) + np.asarray(matrix, dtype=float).T
    ) / 2.0
    eigenvalues, eigenvectors = np.linalg.eigh(symmetric)
    tolerance = (
        np.finfo(float).eps
        * max(1.0, float(np.max(np.abs(symmetric))))
        * max(1, len(symmetric))
        * 100.0
    )
    if float(eigenvalues.min()) < -tolerance:
        raise ValueError("{} is not positive semidefinite.".format(label))
    positive = eigenvalues > tolerance
    if not positive.any():
        return np.zeros((len(symmetric), 0), dtype=float)
    return eigenvectors[:, positive] * np.sqrt(eigenvalues[positive])


def _prepare_response_tip_simulator(
    args: Any,
    gene_tree,
    reconciliation: pd.DataFrame,
    response_inputs: _GeneResponseInputs,
    response: str,
    parameter: float,
    contrast_ids: list[str],
) -> dict[str, Any]:
    """Prepare an O(nodes)-storage PIC forward/inverse bootstrap transform."""
    clades = CladeIndex(gene_tree)
    reconciliation_by_id = {
        str(row["gene_clade_id"]): row for row in reconciliation.to_dict("records")
    }
    orientation = _orient_children(gene_tree, clades, reconciliation_by_id)
    leaf_names = [str(leaf.name) for leaf in gene_tree.leaves()]
    edge_variances = transformed_edge_variances(
        gene_tree,
        model=args.gene_evolution_model,
        parameter=parameter,
        branch_length=args.gene_branch_length,
    )
    variance_by_node = {}
    weights_by_node = {}
    for node in gene_tree.traverse(strategy="postorder"):
        if node.is_leaf:
            variance_by_node[node] = float(edge_variances[node])
            continue
        numerator, denominator = orientation[node]
        numerator_variance = variance_by_node[numerator]
        denominator_variance = variance_by_node[denominator]
        scale = max(numerator_variance, denominator_variance)
        scaled_numerator = numerator_variance / scale
        scaled_denominator = denominator_variance / scale
        total = scaled_numerator + scaled_denominator
        weights_by_node[node] = (
            scaled_denominator / total,
            scaled_numerator / total,
        )
        smaller = min(numerator_variance, denominator_variance)
        larger = max(numerator_variance, denominator_variance)
        variance_by_node[node] = smaller / (1.0 + smaller / larger) + float(
            edge_variances[node]
        )
    node_by_id = {
        clades.clade_id_for_node(node): node
        for node in gene_tree.traverse()
        if not node.is_leaf
    }
    try:
        selected_nodes = [node_by_id[contrast_id] for contrast_id in contrast_ids]
    except KeyError as exc:
        raise ValueError(
            "Bootstrap fit state references an unknown response contrast."
        ) from exc
    if len(set(selected_nodes)) != len(selected_nodes):
        raise ValueError("Bootstrap response contrasts must be unique.")

    evolutionary_model = build_sparse_evolutionary_model(
        gene_tree,
        leaf_names,
        model=args.gene_evolution_model,
        parameter=parameter,
        branch_length=args.gene_branch_length,
    )
    sampling_factor = None
    sampling_by_trait = response_inputs.sampling_covariance_by_trait or {}
    if response in sampling_by_trait:
        sampling_value = sampling_by_trait[response]
        if isinstance(sampling_value, DiagonalLowRankCovariance):
            diagonal_factor = sparse.diags(
                np.sqrt(sampling_value.diagonal), format="csr"
            )
            sampling_factor = sparse.hstack(
                [diagonal_factor, sparse.csr_matrix(sampling_value.low_rank)],
                format="csr",
            )
            sampling_covariance = None
        else:
            sampling_covariance = (
                sampling_value.reindex(index=leaf_names, columns=leaf_names).to_numpy(
                    dtype=float
                )
                if isinstance(sampling_value, pd.DataFrame)
                else np.asarray(sampling_value, dtype=float)
            )
        if sampling_covariance is None:
            pass
        else:
            sampling_diagonal = (
                sampling_covariance
                if sampling_covariance.ndim == 1
                else np.diag(sampling_covariance)
            )
            if sampling_covariance.ndim == 1 or np.array_equal(
                sampling_covariance, np.diag(sampling_diagonal)
            ):
                if np.any(sampling_diagonal < 0.0):
                    raise ValueError("Bootstrap tip sampling covariance is not PSD.")
                sampling_factor = np.sqrt(sampling_diagonal)
            else:
                if len(leaf_names) > 2000 and not args.allow_large_dense:
                    raise ValueError(
                        "Shape-refitted bootstrap has dense tip sampling covariance; "
                        "pass '--allow-large-dense yes' to attempt it."
                    )
                sampling_factor = _positive_semidefinite_factor(
                    sampling_covariance, "Bootstrap tip sampling covariance"
                )
    return {
        "leaf_names": leaf_names,
        "orientation": orientation,
        "weights": weights_by_node,
        "selected_nodes": selected_nodes,
        "evolutionary_model": evolutionary_model,
        "sampling_factor": sampling_factor,
    }


def _simulate_response_tip_values(
    gene_tree,
    fit_state: dict[str, Any],
    simulator: dict[str, Any],
    rng: np.random.Generator,
) -> dict[str, float]:
    leaf_names = simulator["leaf_names"]
    evolutionary_model = simulator["evolutionary_model"]
    base = evolutionary_model.sample(
        rng,
        variance=float(fit_state["evolutionary_rate"])
        * float(evolutionary_model.covariance_scale),
    )
    sampling_factor = simulator["sampling_factor"]
    if sampling_factor is not None:
        if sparse.issparse(sampling_factor):
            base += np.asarray(
                sampling_factor @ rng.standard_normal(sampling_factor.shape[1])
            ).reshape(-1)
        elif np.asarray(sampling_factor).ndim == 1:
            base += sampling_factor * rng.standard_normal(len(leaf_names))
        else:
            base += sampling_factor @ rng.standard_normal(sampling_factor.shape[1])

    leaves = list(gene_tree.leaves())
    estimate_by_node = {leaf: float(base[index]) for index, leaf in enumerate(leaves)}
    contrast_by_node = {}
    orientation = simulator["orientation"]
    weights = simulator["weights"]
    for node in gene_tree.traverse(strategy="postorder"):
        if node.is_leaf:
            continue
        numerator, denominator = orientation[node]
        first = estimate_by_node[numerator]
        second = estimate_by_node[denominator]
        contrast_by_node[node] = first - second
        first_weight, second_weight = weights[node]
        estimate_by_node[node] = first_weight * first + second_weight * second

    contrast_mean = fit_state["design"] @ fit_state["beta"]
    contrast_error = draw_from_factor(
        fit_state["fitted_covariance_factor"],
        rng.standard_normal(len(contrast_mean)),
        rng=rng,
    )
    target = contrast_mean + contrast_error
    for node, value in zip(simulator["selected_nodes"], target, strict=True):
        contrast_by_node[node] = float(value)

    reconstructed = {gene_tree: estimate_by_node[gene_tree]}
    for node in gene_tree.traverse(strategy="preorder"):
        if node.is_leaf:
            continue
        numerator, denominator = orientation[node]
        first_weight, second_weight = weights[node]
        contrast = contrast_by_node[node]
        reconstructed[numerator] = reconstructed[node] + second_weight * contrast
        reconstructed[denominator] = reconstructed[node] - first_weight * contrast
    return {str(leaf.name): float(reconstructed[leaf]) for leaf in leaves}


def _bootstrap_shape_refitted_coefficients(
    args: Any,
    gene_tree,
    reconciliation: pd.DataFrame,
    response_inputs: _GeneResponseInputs,
    response: str,
    species_contrasts: pd.DataFrame,
    predictor_sampling_covariance: pd.DataFrame | None,
    predictors: list[str],
    predictor_group_uncertainties,
    diagnostics: dict[str, Any],
    fit_state: dict[str, Any],
    seed: int,
) -> np.ndarray:
    parameter = float(diagnostics["parameter"])
    simulator = _prepare_response_tip_simulator(
        args,
        gene_tree,
        reconciliation,
        response_inputs,
        response,
        parameter,
        fit_state["contrast_ids"],
    )
    reconciliation_by_id = {
        str(row["gene_clade_id"]): row for row in reconciliation.to_dict("records")
    }
    rng = np.random.default_rng(seed)
    coefficients: list[np.ndarray] = []
    maximum_attempts = max(
        args.bootstrap_replicates * 3, args.bootstrap_replicates + 10
    )
    attempts = 0
    while len(coefficients) < args.bootstrap_replicates and attempts < maximum_attempts:
        attempts += 1
        simulated_values = _simulate_response_tip_values(
            gene_tree,
            fit_state,
            simulator,
            rng,
        )
        simulated_inputs = _GeneResponseInputs(
            values_by_trait={
                **response_inputs.values_by_trait,
                response: simulated_values,
            },
            sampling_covariance_by_trait=response_inputs.sampling_covariance_by_trait,
            replicate_model_by_trait=response_inputs.replicate_model_by_trait,
            tip_summary=response_inputs.tip_summary,
        )
        try:
            contrasts, covariance, _ = _estimate_gene_response_parameter(
                args,
                gene_tree,
                reconciliation_by_id,
                simulated_inputs,
                response,
                species_contrasts,
                predictor_sampling_covariance,
                predictors,
                predictor_group_uncertainties,
            )
            bootstrap_result = _fit_candidate_reconciled_model(
                args,
                contrasts,
                species_contrasts,
                response,
                predictors,
                covariance,
                predictor_sampling_covariance,
                predictor_group_uncertainties,
            ).set_index("term")
        except ValueError:
            continue
        if (bootstrap_result["optimizer_converged"] != "yes").any():
            continue
        coefficients.append(
            bootstrap_result.loc[predictors, "coefficient"].to_numpy(dtype=float)
        )
    if len(coefficients) < args.bootstrap_replicates:
        raise ValueError(
            "Shape-parameter-refitted bootstrap produced only {} successful fits "
            "in {} attempts for response '{}'.".format(
                len(coefficients), attempts, response
            )
        )
    return np.asarray(coefficients, dtype=float)


def _apply_shape_refitted_bootstrap(
    args: Any,
    results: pd.DataFrame,
    gene_tree,
    reconciliation: pd.DataFrame,
    response_inputs: _GeneResponseInputs,
    responses: list[str],
    species_contrasts: pd.DataFrame,
    predictor_sampling_covariance: pd.DataFrame | None,
    predictors: list[str],
    predictor_group_uncertainties,
    response_diagnostics: dict[str, dict[str, Any]],
    fit_states: dict[tuple[str, str], dict[str, Any]],
) -> pd.DataFrame:
    updated = results.copy()
    alpha = (1.0 - args.confidence_level) / 2.0
    for response_index, response in enumerate(responses):
        diagnostics = response_diagnostics[response]
        if diagnostics["parameter_status"] != "estimated":
            continue
        fit_state = fit_states[(str(args.tree_id), response)]
        bootstrap = _bootstrap_shape_refitted_coefficients(
            args,
            gene_tree,
            reconciliation,
            response_inputs,
            response,
            species_contrasts,
            predictor_sampling_covariance,
            predictors,
            predictor_group_uncertainties,
            diagnostics,
            fit_state,
            args.seed + response_index,
        )
        for predictor_index, predictor in enumerate(predictors):
            selected = (updated["response"] == response) & (
                updated["term"] == predictor
            )
            row_index = updated.index[selected].item()
            coefficient = float(updated.loc[row_index, "coefficient"])
            samples = bootstrap[:, predictor_index]
            standard_error = float(np.std(samples, ddof=1))
            lower, upper = np.quantile(samples, [alpha, 1.0 - alpha])
            updated.loc[row_index, "standard_error"] = standard_error
            updated.loc[row_index, "confidence_interval_lower"] = float(lower)
            updated.loc[row_index, "confidence_interval_upper"] = float(upper)
            updated.loc[row_index, "inference_method"] = "parametric-bootstrap"
            updated.loc[row_index, "response_evolution_parameter_bootstrap_refit"] = (
                "yes"
            )
            if standard_error == 0.0:
                updated.loc[row_index, "statistic"] = ""
                updated.loc[row_index, "p_value"] = ""
                updated.loc[row_index, "inference_status"] = "zero-model-variance"
            else:
                updated.loc[row_index, "statistic"] = coefficient / standard_error
                centered = samples - coefficient
                updated.loc[row_index, "p_value"] = float(
                    (1 + np.sum(np.abs(centered) >= abs(coefficient)))
                    / (len(centered) + 1)
                )
                updated.loc[row_index, "inference_status"] = "ok"
    return updated


def _lift_predictor_uncertainties_to_gene_tips(
    encoded_predictors,
    predictor_inputs,
    species_leaf_names,
    gene_species_indices,
    *,
    coefficient_offset,
):
    uncertainties: list[
        ContinuousPredictorUncertainty
        | GmrfPredictorUncertainty
        | GroupedPredictorUncertainty
    ] = []
    columns: list[int | tuple[int, ...]] = []
    term_index = {
        term.name: index + coefficient_offset
        for index, term in enumerate(encoded_predictors.terms)
    }
    sampling = predictor_inputs.sampling_covariance_by_trait or {}
    sparse_posteriors = predictor_inputs.sparse_posterior_by_trait or {}
    for term in encoded_predictors.terms:
        if term.kind != "continuous" or (
            term.name not in sampling and term.name not in sparse_posteriors
        ):
            continue
        observation_index = np.asarray(gene_species_indices, dtype=int)
        if term.name in sparse_posteriors:
            uncertainties.append(
                GmrfPredictorUncertainty(
                    model=sparse_posteriors[term.name],
                    observation_index=observation_index,
                )
            )
        else:
            covariance = sampling[term.name]
            if isinstance(covariance, pd.DataFrame):
                matrix = covariance.loc[
                    species_leaf_names, species_leaf_names
                ].to_numpy(dtype=float)
            else:
                matrix = np.asarray(covariance, dtype=float)
            uncertainties.append(
                ContinuousPredictorUncertainty(
                    factor=_positive_semidefinite_factor(
                        matrix, "Continuous predictor posterior covariance"
                    ),
                    observation_index=observation_index,
                )
            )
        columns.append(term_index[term.name])
    for uncertainty in encoded_predictors.uncertainties:
        selected = np.asarray(gene_species_indices, dtype=int)
        species_covariance = uncertainty.covariance_by_observation
        factors = tuple(
            _positive_semidefinite_factor(
                covariance,
                "Categorical predictor covariance",
            )
            for covariance in species_covariance
        )
        uncertainties.append(
            GroupedPredictorUncertainty(
                factors=factors,
                observation_index=selected,
            )
        )
        columns.append(tuple(term_index[name] for name in uncertainty.term_names))
    return uncertainties, columns


def _reconciled_pglmm_base_row(
    args,
    response,
    response_spec,
    fit,
    *,
    n_gene_tips,
    n_species,
    num_parameters,
    matrix_rank,
    predictor_uncertainties,
):
    evolution_spec = evolution_model_spec(args.gene_evolution_model)
    parameter_status = (
        "not-applicable"
        if evolution_spec.parameter_name is None
        else fit.evolution_parameter_status
    )
    return {
        "model_id": "{}:{}".format(args.tree_id, response),
        "tree_id": args.tree_id,
        "predictor_tree_id": "species",
        "response": response,
        "response_type": response_spec.kind,
        "response_family": response_spec.family,
        "response_level": "",
        "response_reference": fit.reference,
        "link_function": response_family_link(response_spec.family),
        "response_dispersion": "" if fit.dispersion is None else fit.dispersion,
        "zero_probability": (
            "" if fit.zero_probability is None else fit.zero_probability
        ),
        "coefficient_penalty": fit.coefficient_penalty,
        "coefficient_prior_sd": (
            "" if fit.coefficient_prior_sd is None else fit.coefficient_prior_sd
        ),
        "separation_warning": "yes" if fit.separation_warning else "no",
        "degrees_of_freedom": 1,
        "n_gene_contrasts": 0,
        "n_species_events": 0,
        "n_repeated_gene_contrasts": 0,
        "n_lineages": n_gene_tips,
        "n_excluded_ineligible": 0,
        "n_excluded_coverage": 0,
        "num_parameters": num_parameters,
        "matrix_rank": matrix_rank,
        "condition_number": (
            float(np.linalg.cond(fit.coefficient_covariance))
            if np.isfinite(fit.coefficient_covariance).all()
            else ""
        ),
        "weighted_residual_sum_squares": "",
        "residual_scale": "",
        "r_squared_uncentered": "",
        "intercept": "no" if response_spec.family == "ordinal" else "yes",
        "event_weighting": "not-applicable",
        "covariance_estimator": "laplace-ML",
        "contrast_transform": "not-applicable-tip-pglmm",
        "response_evolution_model": args.gene_evolution_model,
        "response_evolution_parameter_name": evolution_spec.parameter_name or "",
        "response_evolution_parameter": (
            "" if fit.evolution_parameter is None else fit.evolution_parameter
        ),
        "response_evolution_parameter_status": parameter_status,
        "response_evolution_optimizer_converged": (
            "yes" if fit.optimizer_converged else "no"
        ),
        "response_evolution_optimizer_message": fit.optimizer_message,
        "response_evolution_boundary_warning": (
            "yes" if fit.boundary_warning else "no"
        ),
        "response_evolution_parameter_bootstrap_refit": "no",
        "response_branch_length_mode": args.gene_branch_length,
        "coverage_policy": args.speciation_coverage,
        "small_sample_warning": "yes" if n_species < 20 else "no",
        "inference_status": "ok",
        "model": "reconciled-tip-pglmm",
        "inference_method": fit.coefficient_inference,
        "reml": "no",
        "evolutionary_rate": fit.component_variances["phylogenetic"],
        "species_event_variance": fit.component_variances.get("group", 0.0),
        "lineage_slope_variance": 0.0,
        "mean_sampling_variance": "",
        "sampling_variance_fraction": "",
        "mean_predictor_sampling_variance": "",
        "mean_latent_predictor_variance": "",
        "predictor_uncertainty_fraction": "",
        "predictor_evolutionary_rate": "",
        "predictor_rate_optimizer_converged": "",
        "predictor_rate_optimizer_message": "",
        "predictor_rate_boundary_warning": "",
        "measurement_error_model": (
            "latent-predictor" if predictor_uncertainties else "none"
        ),
        "log_likelihood": fit.log_likelihood,
        "optimizer_converged": "yes" if fit.optimizer_converged else "no",
        "boundary_warning": "yes" if fit.boundary_warning else "no",
        "event_random_effect": "species-tip"
        if "group" in fit.component_variances
        else "no",
        "lineage_random_slope": "no",
    }


def _reconciled_pglmm_rows(
    args,
    response,
    response_spec,
    fit,
    terms,
    predictor_groups,
    predictor_diagnostics,
    *,
    confidence_level,
    n_gene_tips,
    n_species,
    matrix_rank,
    predictor_uncertainties,
):
    coefficient_count = fit.coefficients.size
    num_parameters = (
        coefficient_count
        + len(fit.thresholds)
        + len(fit.component_variances)
        + (1 if fit.evolution_parameter_status == "estimated" else 0)
        + (1 if fit.dispersion is not None else 0)
        + (1 if fit.zero_probability is not None else 0)
    )
    base = _reconciled_pglmm_base_row(
        args,
        response,
        response_spec,
        fit,
        n_gene_tips=n_gene_tips,
        n_species=n_species,
        num_parameters=num_parameters,
        matrix_rank=matrix_rank,
        predictor_uncertainties=predictor_uncertainties,
    )
    critical = float(norm.ppf(0.5 + confidence_level / 2.0))
    dimension = fit.coefficients.shape[1]
    response_levels = (
        [""]
        if response_spec.family == "ordinal"
        or response_spec.family in SCALAR_RESPONSE_FAMILIES
        else list(fit.levels[:-1])
    )
    rows = []
    row_by_flat_index = {}
    for term_index, term in enumerate(terms):
        for level_index, level in enumerate(response_levels):
            flat_index = term_index * dimension + level_index
            coefficient = float(fit.coefficients[term_index, level_index])
            (
                standard_error,
                statistic,
                p_value,
                lower,
                upper,
                inference_status,
            ) = summarize_glmm_coefficient(fit, flat_index, coefficient, critical)
            row = base.copy()
            row.update(
                {
                    "response_level": level,
                    "term": term.name,
                    "source_term": term.source,
                    "predictor_type": term.kind,
                    "predictor_level": term.level,
                    "predictor_reference": term.reference,
                    "factor_coding": term.coding,
                    "term_test": "coefficient",
                    "coefficient": coefficient,
                    "standard_error": standard_error,
                    "statistic": statistic,
                    "p_value": p_value,
                    "confidence_level": confidence_level,
                    "confidence_interval_lower": lower,
                    "confidence_interval_upper": upper,
                    "inference_status": inference_status,
                }
            )
            if term.kind != "intercept":
                diagnostics = predictor_diagnostics[term.name]
                predictor_spec = evolution_model_spec(args.species_evolution_model)
                row.update(
                    {
                        "predictor_evolution_model": args.species_evolution_model,
                        "predictor_evolution_parameter_name": (
                            predictor_spec.parameter_name or ""
                        ),
                        "predictor_evolution_parameter": diagnostics["parameter"],
                        "predictor_evolution_parameter_status": diagnostics[
                            "parameter_status"
                        ],
                        "predictor_evolution_optimizer_converged": (
                            _diagnostic_flag(diagnostics["optimizer_converged"])
                        ),
                        "predictor_evolution_optimizer_message": diagnostics[
                            "optimizer_message"
                        ],
                        "predictor_evolution_boundary_warning": (
                            _diagnostic_flag(diagnostics["boundary_warning"])
                        ),
                        "predictor_evolution_log_likelihood": diagnostics[
                            "log_likelihood"
                        ],
                        "predictor_branch_length_mode": args.species_branch_length,
                    }
                )
                if "mean_posterior_variance" in diagnostics:
                    row.update(
                        {
                            "mean_predictor_sampling_variance": diagnostics[
                                "mean_sampling_variance"
                            ],
                            "mean_latent_predictor_variance": diagnostics[
                                "mean_posterior_variance"
                            ],
                            "predictor_uncertainty_fraction": diagnostics[
                                "uncertainty_fraction"
                            ],
                            "predictor_evolutionary_rate": diagnostics[
                                "evolutionary_rate"
                            ],
                            "predictor_rate_optimizer_converged": _diagnostic_flag(
                                diagnostics["rate_optimizer_converged"]
                            ),
                            "predictor_rate_optimizer_message": diagnostics[
                                "rate_optimizer_message"
                            ],
                            "predictor_rate_boundary_warning": _diagnostic_flag(
                                diagnostics["rate_boundary_warning"]
                            ),
                        }
                    )
            row_by_flat_index[flat_index] = row
            rows.append(row)
    term_index_by_name = {term.name: index for index, term in enumerate(terms)}
    for source, names in predictor_groups.items():
        if dimension == 1 and len(names) == 1:
            continue
        indices = [
            term_index_by_name[name] * dimension + response_index
            for name in names
            for response_index in range(dimension)
        ]
        omnibus_statistic, omnibus_p_value, omnibus_status = summarize_glmm_omnibus(
            fit, indices
        )
        row = row_by_flat_index[indices[0]].copy()
        row.update(
            {
                "response_level": "",
                "term": source,
                "source_term": source,
                "predictor_level": "",
                "term_test": "omnibus",
                "coefficient": "",
                "standard_error": "",
                "statistic": omnibus_statistic,
                "degrees_of_freedom": len(indices),
                "p_value": omnibus_p_value,
                "confidence_interval_lower": "",
                "confidence_interval_upper": "",
                "inference_method": "wald",
                "inference_status": omnibus_status,
            }
        )
        rows.append(row)
    for threshold_index, threshold in enumerate(fit.thresholds):
        (
            threshold_standard_error,
            lower,
            upper,
            threshold_status,
        ) = summarize_glmm_threshold(fit, threshold_index, critical)
        row = base.copy()
        row.update(
            {
                "response_level": fit.levels[threshold_index],
                "term": "threshold[{}]".format(fit.levels[threshold_index]),
                "source_term": "(threshold)",
                "predictor_type": "threshold",
                "term_test": "threshold",
                "coefficient": float(threshold),
                "standard_error": threshold_standard_error,
                "confidence_level": confidence_level,
                "confidence_interval_lower": lower,
                "confidence_interval_upper": upper,
                "inference_status": threshold_status,
            }
        )
        rows.append(row)
    return rows


def _reconciled_categorical_context(gene_tree, species_tree, species_labels):
    gene_tip_names = [str(leaf.name) for leaf in gene_tree.leaves()]
    species_leaf_names = [str(leaf.name) for leaf in species_tree.leaves()]
    species_index = {name: index for index, name in enumerate(species_leaf_names)}
    unmapped = [name for name in gene_tip_names if species_labels[name] is None]
    if unmapped:
        raise ValueError(
            "Categorical reconciled responses require every gene tip to map to a "
            "species (unmapped: {}).".format(", ".join(unmapped))
        )
    gene_species = [str(species_labels[name]) for name in gene_tip_names]
    gene_species_indices = [species_index[name] for name in gene_species]
    return gene_tip_names, species_leaf_names, gene_species, gene_species_indices


def _reconciled_categorical_design(gene_tip_names, gene_species, encoded_predictors):
    predictor_columns = [
        [
            encoded_predictors.values_by_trait[term.name][species]
            for species in gene_species
        ]
        for term in encoded_predictors.terms
    ]
    design = np.column_stack(
        [np.ones(len(gene_tip_names), dtype=float)] + predictor_columns
    )
    terms = [PredictorTerm("(intercept)", "(intercept)", "intercept")] + list(
        encoded_predictors.terms
    )
    group_covariance = (
        None
        if len(set(gene_species)) == len(gene_species)
        else sparse_group_covariance(gene_species)
    )
    return design, terms, group_covariance


def _gene_covariance_factory(args, gene_tree, gene_tip_names):
    return evolutionary_covariance_factory(
        gene_tree,
        gene_tip_names,
        model=args.gene_evolution_model,
        branch_length=args.gene_branch_length,
    )


def _categorical_random_row(
    args,
    response,
    effect_type,
    group_id,
    level,
    conditional_mode,
    variance_component,
    n_observations,
):
    return {
        "model_id": "{}:{}".format(args.tree_id, response),
        "tree_id": args.tree_id,
        "response": response,
        "effect_type": effect_type,
        "group_id": group_id,
        "term": level,
        "conditional_mode": float(conditional_mode),
        "variance_component": variance_component,
        "n_observations": n_observations,
    }


def _categorical_tip_random_rows(
    args, response, fit, gene_tip_names, level_names, effect_type
):
    modes = fit.component_random_modes["phylogenetic"]
    return [
        _categorical_random_row(
            args,
            response,
            effect_type,
            tip_name,
            level,
            modes[tip_index, level_index],
            fit.component_variances["phylogenetic"],
            1,
        )
        for tip_index, tip_name in enumerate(gene_tip_names)
        for level_index, level in enumerate(level_names)
    ]


def _categorical_group_random_rows(
    args, response, fit, gene_species, level_names, effect_type
):
    if "group" not in fit.component_random_modes:
        return []
    modes = fit.component_random_modes["group"]
    rows = []
    for species in sorted(set(gene_species)):
        members = [
            index for index, value in enumerate(gene_species) if value == species
        ]
        for level_index, level in enumerate(level_names):
            rows.append(
                _categorical_random_row(
                    args,
                    response,
                    effect_type,
                    species,
                    level,
                    np.mean(modes[members, level_index]),
                    fit.component_variances["group"],
                    len(members),
                )
            )
    return rows


def _categorical_random_rows(
    args, response, response_spec, fit, gene_tip_names, gene_species
):
    scalar = response_spec.family in SCALAR_RESPONSE_FAMILIES
    level_names = (
        [""] if response_spec.family == "ordinal" or scalar else list(fit.levels[:-1])
    )
    tip_effect = (
        "phylogenetic_tip_linear_predictor" if scalar else "phylogenetic_tip_logit"
    )
    group_effect = "species_linear_predictor" if scalar else "species_logit"
    return _categorical_tip_random_rows(
        args, response, fit, gene_tip_names, level_names, tip_effect
    ) + _categorical_group_random_rows(
        args, response, fit, gene_species, level_names, group_effect
    )


def _fit_one_reconciled_categorical_response(
    args,
    response,
    response_spec,
    response_inputs,
    gene_tree,
    gene_tip_names,
    species_leaf_names,
    gene_species_indices,
    base_design,
    base_terms,
    group_covariance,
    encoded_predictors,
    predictor_inputs,
    response_offsets,
    response_trials,
    response_censor_lower,
    response_censor_upper,
    response_dispersions,
    response_zero_probabilities,
):
    ordinal = response_spec.family == "ordinal"
    design = base_design[:, 1:] if ordinal else base_design
    terms = base_terms[1:] if ordinal else base_terms
    if (
        len(design) <= design.shape[1]
        or np.linalg.matrix_rank(design) < design.shape[1]
    ):
        raise ValueError(
            "Reconciled categorical model '{}' has a rank-deficient or saturated "
            "design.".format(response)
        )
    predictor_uncertainties, predictor_columns = (
        _lift_predictor_uncertainties_to_gene_tips(
            encoded_predictors,
            predictor_inputs,
            species_leaf_names,
            gene_species_indices,
            coefficient_offset=0 if ordinal else 1,
        )
    )
    fixed_parameter = (
        None
        if args.gene_evolution_parameter in {None, "auto"}
        else float(args.gene_evolution_parameter)
    )
    shape_bounds, shape_decoder, shape_initial = _categorical_shape_settings(
        gene_tree,
        args.gene_evolution_model,
        fixed_parameter,
        args.gene_branch_length,
    )
    fit = fit_phylogenetic_glmm(
        [response_inputs.values_by_trait[response][name] for name in gene_tip_names],
        design,
        _gene_covariance_factory(args, gene_tree, gene_tip_names),
        family=response_spec.family,
        levels=response_spec.levels,
        reference=response_spec.reference,
        group_covariance=group_covariance,
        evolution_parameter=fixed_parameter,
        evolution_parameter_bounds=shape_bounds,
        evolution_parameter_decoder=shape_decoder,
        evolution_parameter_initial=shape_initial,
        predictor_uncertainties=predictor_uncertainties,
        predictor_columns=predictor_columns,
        offset=(
            None
            if response not in response_offsets
            else [response_offsets[response][name] for name in gene_tip_names]
        ),
        trials=(
            None
            if response not in response_trials
            else [response_trials[response][name] for name in gene_tip_names]
        ),
        censor_lower=(
            None
            if response not in response_censor_lower
            else [response_censor_lower[response][name] for name in gene_tip_names]
        ),
        censor_upper=(
            None
            if response not in response_censor_upper
            else [response_censor_upper[response][name] for name in gene_tip_names]
        ),
        dispersion=response_dispersions.get(response),
        zero_probability=response_zero_probabilities.get(response),
        coefficient_penalty=args.coefficient_penalty,
        coefficient_prior_sd=args.coefficient_prior_sd,
        inference=args.inference,
        confidence_level=args.confidence_level,
        bootstrap_replicates=args.bootstrap_replicates,
        seed=args.seed,
        allow_large_dense=args.allow_large_dense,
    )
    return fit, design, terms, predictor_uncertainties


def _fit_reconciled_categorical_responses(
    args,
    gene_tree,
    species_tree,
    species_labels,
    response_inputs,
    response_specs,
    responses,
    encoded_predictors,
    predictor_inputs,
    predictor_diagnostics,
    response_offsets,
    response_trials,
    response_censor_lower,
    response_censor_upper,
    response_dispersions,
    response_zero_probabilities,
):
    gene_tip_names, species_leaf_names, gene_species, gene_species_indices = (
        _reconciled_categorical_context(gene_tree, species_tree, species_labels)
    )
    base_design, base_terms, group_covariance = _reconciled_categorical_design(
        gene_tip_names, gene_species, encoded_predictors
    )
    rows = []
    random_rows = []
    for response in responses:
        response_spec = response_specs[response]
        fit, design, terms, predictor_uncertainties = (
            _fit_one_reconciled_categorical_response(
                args,
                response,
                response_spec,
                response_inputs,
                gene_tree,
                gene_tip_names,
                species_leaf_names,
                gene_species_indices,
                base_design,
                base_terms,
                group_covariance,
                encoded_predictors,
                predictor_inputs,
                response_offsets,
                response_trials,
                response_censor_lower,
                response_censor_upper,
                response_dispersions,
                response_zero_probabilities,
            )
        )
        rows.extend(
            _reconciled_pglmm_rows(
                args,
                response,
                response_spec,
                fit,
                terms,
                encoded_predictors.groups,
                predictor_diagnostics,
                confidence_level=args.confidence_level,
                n_gene_tips=len(gene_tip_names),
                n_species=len(set(gene_species)),
                matrix_rank=int(np.linalg.matrix_rank(design)),
                predictor_uncertainties=predictor_uncertainties,
            )
        )
        random_rows.extend(
            _categorical_random_rows(
                args, response, response_spec, fit, gene_tip_names, gene_species
            )
        )
    return (
        pd.DataFrame(rows, columns=RESULT_COLUMNS),
        pd.DataFrame(random_rows, columns=RANDOM_EFFECT_COLUMNS),
    )


def _reconciled_multivariate_rows(
    args,
    responses,
    fit,
    terms,
    predictor_diagnostics,
    *,
    n_gene_tips,
    n_species,
    matrix_rank,
):
    critical = float(norm.ppf(0.5 + args.confidence_level / 2.0))
    evolution_spec = evolution_model_spec(args.gene_evolution_model)
    phylogenetic_trait_covariance = fit.component_trait_covariances["phylogenetic"]
    species_trait_covariance = fit.component_trait_covariances.get(
        "species", np.zeros_like(phylogenetic_trait_covariance)
    )
    coefficient_count = len(terms)
    num_parameters = (
        fit.coefficients.size
        + sum(
            covariance.shape[0] * (covariance.shape[0] + 1) // 2
            for covariance in fit.component_trait_covariances.values()
        )
        + (1 if fit.evolution_parameter_status == "estimated" else 0)
    )
    base = {
        "model_id": "{}:multivariate".format(args.tree_id),
        "tree_id": args.tree_id,
        "predictor_tree_id": "species",
        "response_type": "multivariate-continuous",
        "response_family": "multivariate-gaussian",
        "response_level": "",
        "response_reference": "",
        "link_function": "identity",
        "response_dispersion": "",
        "zero_probability": "",
        "coefficient_penalty": "none",
        "coefficient_prior_sd": "",
        "separation_warning": "no",
        "degrees_of_freedom": fit.n_observations - fit.coefficients.size,
        "n_gene_contrasts": 0,
        "n_species_events": 0,
        "n_repeated_gene_contrasts": 0,
        "n_lineages": n_gene_tips,
        "n_excluded_ineligible": 0,
        "n_excluded_coverage": 0,
        "num_parameters": num_parameters,
        "matrix_rank": matrix_rank,
        "condition_number": float(np.linalg.cond(fit.coefficient_covariance)),
        "intercept": "yes",
        "event_weighting": "not-applicable",
        "covariance_estimator": (
            "multivariate-gaussian-REML" if fit.reml else "multivariate-gaussian-ML"
        ),
        "contrast_transform": "not-applicable-tip-multivariate-pgls",
        "response_evolution_model": args.gene_evolution_model,
        "response_evolution_parameter_name": evolution_spec.parameter_name or "",
        "response_evolution_parameter": (
            "" if fit.evolution_parameter is None else fit.evolution_parameter
        ),
        "response_evolution_parameter_status": (
            "not-applicable"
            if evolution_spec.parameter_name is None
            else fit.evolution_parameter_status
        ),
        "response_evolution_optimizer_converged": (
            "yes" if fit.optimizer_converged else "no"
        ),
        "response_evolution_optimizer_message": fit.optimizer_message,
        "response_evolution_boundary_warning": (
            "yes" if fit.boundary_warning else "no"
        ),
        "response_evolution_parameter_bootstrap_refit": "no",
        "response_branch_length_mode": args.gene_branch_length,
        "coverage_policy": args.speciation_coverage,
        "small_sample_warning": "yes" if n_species < 20 else "no",
        "inference_status": "ok",
        "model": "reconciled-tip-multivariate-pgls",
        "inference_method": "wald",
        "reml": "yes" if fit.reml else "no",
        "measurement_error_model": "response-covariance",
        "log_likelihood": fit.log_likelihood,
        "optimizer_converged": "yes" if fit.optimizer_converged else "no",
        "boundary_warning": "yes" if fit.boundary_warning else "no",
        "event_random_effect": "species-tip"
        if "species" in fit.component_trait_covariances
        else "no",
        "lineage_random_slope": "no",
    }
    rows = []
    for response_index, response in enumerate(responses):
        for term_index, term in enumerate(terms):
            flat_index = response_index * coefficient_count + term_index
            coefficient = float(fit.coefficients[response_index, term_index])
            standard_error = math.sqrt(
                max(float(fit.coefficient_covariance[flat_index, flat_index]), 0.0)
            )
            statistic = "" if standard_error == 0.0 else coefficient / standard_error
            row = base.copy()
            row.update(
                {
                    "response": response,
                    "term": term.name,
                    "source_term": term.source,
                    "predictor_type": term.kind,
                    "predictor_level": term.level,
                    "predictor_reference": term.reference,
                    "factor_coding": term.coding,
                    "term_test": "coefficient",
                    "coefficient": coefficient,
                    "standard_error": standard_error,
                    "statistic": statistic,
                    "p_value": (
                        ""
                        if standard_error == 0.0
                        else float(2.0 * norm.sf(abs(float(statistic))))
                    ),
                    "confidence_level": args.confidence_level,
                    "confidence_interval_lower": coefficient
                    - critical * standard_error,
                    "confidence_interval_upper": coefficient
                    + critical * standard_error,
                    "evolutionary_rate": phylogenetic_trait_covariance[
                        response_index, response_index
                    ],
                    "species_event_variance": species_trait_covariance[
                        response_index, response_index
                    ],
                }
            )
            if term.kind != "intercept":
                diagnostics = predictor_diagnostics[term.name]
                row["predictor_evolution_model"] = args.species_evolution_model
                row["predictor_evolution_parameter"] = diagnostics["parameter"]
            rows.append(row)
    for component_name, covariance in fit.component_trait_covariances.items():
        for first, first_response in enumerate(responses):
            for second in range(first, len(responses)):
                row = base.copy()
                second_response = responses[second]
                row.update(
                    {
                        "response": "{}|{}".format(first_response, second_response),
                        "term": "{}-covariance[{},{}]".format(
                            component_name, first_response, second_response
                        ),
                        "source_term": "(response-covariance)",
                        "predictor_type": "response-covariance",
                        "term_test": "response-covariance",
                        "coefficient": covariance[first, second],
                        "confidence_level": args.confidence_level,
                    }
                )
                rows.append(row)
    return pd.DataFrame(rows, columns=RESULT_COLUMNS)


def _fit_reconciled_multivariate_responses(
    args,
    gene_tree,
    species_tree,
    species_labels,
    response_inputs,
    responses,
    encoded_predictors,
    predictor_diagnostics,
):
    gene_tip_names, _species_names, gene_species, _indices = (
        _reconciled_categorical_context(gene_tree, species_tree, species_labels)
    )
    design, terms, group_covariance = _reconciled_categorical_design(
        gene_tip_names, gene_species, encoded_predictors
    )
    if encoded_predictors.uncertainties:
        raise ValueError(
            "Reconciled multivariate PGLS does not yet combine categorical "
            "predictor uncertainty."
        )
    response_matrix = np.asarray(
        [
            [response_inputs.values_by_trait[response][name] for response in responses]
            for name in gene_tip_names
        ],
        dtype=float,
    )
    fixed_diagonal = np.zeros(response_matrix.size, dtype=float)
    dense_sampling_blocks: list[tuple[int, np.ndarray]] = []
    fixed_loading_blocks = []
    for response_index, response in enumerate(responses):
        sampling = response_inputs.sampling_covariance_by_trait or {}
        if response not in sampling:
            fixed_loading_blocks.append(
                sparse.csr_matrix((len(gene_tip_names), 0), dtype=float)
            )
            continue
        sampling_value = sampling[response]
        covariance = (
            sampling_value
            if isinstance(sampling_value, DiagonalLowRankCovariance)
            else sampling_value.loc[gene_tip_names, gene_tip_names].to_numpy(float)
            if isinstance(sampling_value, pd.DataFrame)
            else np.asarray(sampling_value, dtype=float)
        )
        start = response_index * len(gene_tip_names)
        if isinstance(covariance, DiagonalLowRankCovariance):
            fixed_diagonal[start : start + len(gene_tip_names)] = covariance.diagonal
            fixed_loading_blocks.append(sparse.csr_matrix(covariance.low_rank))
        elif covariance.ndim == 1 or np.array_equal(
            covariance, np.diag(np.diag(covariance))
        ):
            fixed_diagonal[start : start + len(gene_tip_names)] = (
                covariance if covariance.ndim == 1 else np.diag(covariance)
            )
            fixed_loading_blocks.append(
                sparse.csr_matrix((len(gene_tip_names), 0), dtype=float)
            )
        else:
            dense_sampling_blocks.append((start, covariance))
            fixed_loading_blocks.append(None)
    if dense_sampling_blocks:
        fixed_covariance = np.diag(fixed_diagonal)
        for response_index, loading in enumerate(fixed_loading_blocks):
            if loading is not None and loading.shape[1]:
                update = loading @ loading.T
                start = response_index * len(gene_tip_names)
                fixed_covariance[
                    start : start + len(gene_tip_names),
                    start : start + len(gene_tip_names),
                ] += update.toarray()
        for start, covariance in dense_sampling_blocks:
            fixed_covariance[
                start : start + len(gene_tip_names),
                start : start + len(gene_tip_names),
            ] = covariance
    elif any(loading.shape[1] for loading in fixed_loading_blocks):
        fixed_covariance = DiagonalLowRankCovariance(
            fixed_diagonal,
            sparse.block_diag(fixed_loading_blocks, format="csr"),
        )
    else:
        fixed_covariance = fixed_diagonal
    fixed_parameter = (
        None
        if args.gene_evolution_parameter in {None, "auto"}
        else float(args.gene_evolution_parameter)
    )
    shape_bounds, shape_decoder, shape_initial = _categorical_shape_settings(
        gene_tree,
        args.gene_evolution_model,
        fixed_parameter,
        args.gene_branch_length,
    )
    components = {
        "phylogenetic": _gene_covariance_factory(args, gene_tree, gene_tip_names)
    }
    if group_covariance is not None:
        components["species"] = group_covariance
    fit = fit_multivariate_pgls(
        response_matrix,
        design,
        components,
        fixed_covariance=fixed_covariance,
        evolution_parameter=fixed_parameter,
        evolution_parameter_bounds=shape_bounds,
        evolution_parameter_decoder=shape_decoder,
        evolution_parameter_initial=shape_initial,
        reml=args.reml,
        allow_large_dense=args.allow_large_dense,
    )
    return _reconciled_multivariate_rows(
        args,
        responses,
        fit,
        terms,
        predictor_diagnostics,
        n_gene_tips=len(gene_tip_names),
        n_species=len(set(gene_species)),
        matrix_rank=int(np.linalg.matrix_rank(design)),
    )


def _validate_reconciled_response_modes(raw_args, responses, response_specs):
    continuous = [
        response
        for response in responses
        if response_specs[response].family == "gaussian"
    ]
    non_gaussian = [
        response
        for response in responses
        if response_specs[response].family != "gaussian"
    ]
    if raw_args.multivariate_responses:
        if len(responses) < 2:
            raise ValueError("Multivariate PGLS requires at least two responses.")
        if non_gaussian:
            raise ValueError("Multivariate PGLS currently requires Gaussian responses.")
        if raw_args.inference != "wald":
            raise ValueError("Multivariate PGLS currently supports Wald inference.")
        if raw_args.model != "hierarchical":
            raise ValueError(
                "Reconciled multivariate responses require '--model hierarchical'."
            )
    if non_gaussian and raw_args.model != "hierarchical":
        raise ValueError(
            "Non-Gaussian reconciled responses require '--model hierarchical'."
        )
    if continuous and raw_args.inference in {
        "likelihood-ratio",
        "profile-likelihood",
    }:
        raise ValueError(
            "Gaussian reconciled PGLS supports Wald or parametric-bootstrap inference."
        )
    return continuous, non_gaussian


def _merge_predictor_covariance_audit(
    predictor_sampling_covariance, grouped_predictor_covariance_audit
):
    fit_covariance = predictor_sampling_covariance
    if grouped_predictor_covariance_audit.empty:
        return predictor_sampling_covariance, fit_covariance
    if predictor_sampling_covariance is None:
        return grouped_predictor_covariance_audit, fit_covariance
    scalar_audit = predictor_sampling_covariance.copy()
    scalar_audit["trait_2"] = scalar_audit["trait"]
    audit = pd.concat(
        [scalar_audit, grouped_predictor_covariance_audit],
        ignore_index=True,
        sort=False,
    )
    return audit, fit_covariance


def _reconciled_multivariate_artifacts(
    raw_args,
    gene_tree,
    species_tree,
    species_labels,
    reconciliation,
    response_inputs,
    responses,
    encoded_predictors,
    predictor_inputs,
    predictor_diagnostics,
    species_contrasts,
    predictor_sampling_covariance,
    predictor_sampling_covariance_for_fit,
    predictor_group_uncertainties,
    trait_origins,
):
    if (
        predictor_sampling_covariance_for_fit is not None
        or predictor_group_uncertainties
    ):
        raise ValueError(
            "Reconciled multivariate PGLS does not yet combine predictor "
            "measurement error; fit separate response models for that combination."
        )
    results = _fit_reconciled_multivariate_responses(
        raw_args,
        gene_tree,
        species_tree,
        species_labels,
        response_inputs,
        responses,
        encoded_predictors,
        predictor_diagnostics,
    )
    return PglsPipelineArtifacts(
        reconciliation=reconciliation,
        gene_contrasts=pd.DataFrame(columns=sorted(RESPONSE_REQUIRED_COLUMNS)),
        species_contrasts=species_contrasts,
        response_sampling_covariance=None,
        response_tip_summary=response_inputs.tip_summary,
        predictor_sampling_covariance=predictor_sampling_covariance,
        predictor_tip_summary=predictor_inputs.tip_summary,
        results=results.reindex(columns=RESULT_COLUMNS),
        random_effects=pd.DataFrame(columns=RANDOM_EFFECT_COLUMNS),
        trait_origins=trait_origins,
    )


def _fit_reconciled_gaussian_responses(
    raw_args,
    gene_tree,
    reconciliation,
    response_inputs,
    continuous_responses,
    species_contrasts,
    predictor_sampling_covariance_for_fit,
    encoded_predictor_names,
    predictor_group_uncertainties,
    encoded_predictors,
    predictor_diagnostics,
    sensitivity_omissions,
):
    if not continuous_responses:
        return (
            pd.DataFrame(columns=sorted(RESPONSE_REQUIRED_COLUMNS)),
            None,
            {},
            pd.DataFrame(columns=RESULT_COLUMNS),
            pd.DataFrame(columns=RANDOM_EFFECT_COLUMNS),
            pd.DataFrame(columns=SENSITIVITY_COLUMNS),
            {},
            False,
        )
    gene_contrasts, sampling_covariance, response_diagnostics = _build_gene_contrasts(
        raw_args,
        gene_tree,
        reconciliation,
        response_inputs,
        continuous_responses,
        species_contrasts,
        predictor_sampling_covariance_for_fit,
        encoded_predictor_names,
        predictor_group_uncertainties,
    )
    refit_shape = raw_args.inference == "parametric-bootstrap" and any(
        diagnostics["parameter_status"] == "estimated"
        for diagnostics in response_diagnostics.values()
    )
    fitted = fit_reconciled_pgls(
        gene_contrasts,
        species_contrasts,
        continuous_responses,
        encoded_predictor_names,
        confidence_level=raw_args.confidence_level,
        event_weighting=raw_args.event_weighting,
        coverage_policy=raw_args.speciation_coverage,
        model=raw_args.model,
        response_sampling_covariance=sampling_covariance,
        predictor_sampling_covariance=predictor_sampling_covariance_for_fit,
        predictor_group_uncertainties=predictor_group_uncertainties,
        inference="wald" if refit_shape else raw_args.inference,
        bootstrap_replicates=raw_args.bootstrap_replicates,
        seed=raw_args.seed,
        reml=raw_args.reml,
        event_random_effect=raw_args.event_random_effect,
        lineage_random_slope=raw_args.lineage_random_slope,
        lineage_inference=raw_args.lineage_inference,
        lineage_leave_one_out=raw_args.lineage_leave_one_out,
        sensitivity_omissions=sensitivity_omissions,
        return_random_effects=True,
        return_sensitivity=True,
        return_fit_state=refit_shape,
        predictor_metadata=encoded_predictors.metadata_by_term,
        predictor_groups=encoded_predictors.groups,
        allow_large_dense=raw_args.allow_large_dense,
    )
    if refit_shape:
        results, random_effects, sensitivity, fit_states = fitted
    else:
        results, random_effects, sensitivity = fitted
        fit_states = {}
    results = _attach_evolution_diagnostics(
        results,
        response_diagnostics,
        predictor_diagnostics,
        encoded_predictors.groups,
    )
    return (
        gene_contrasts,
        sampling_covariance,
        response_diagnostics,
        results,
        random_effects,
        sensitivity,
        fit_states,
        refit_shape,
    )


def _merge_reconciled_result_frames(continuous, non_gaussian, columns):
    frames = [frame for frame in [continuous, non_gaussian] if not frame.empty]
    if not frames:
        return pd.DataFrame(columns=columns)
    return pd.concat(frames, ignore_index=True).reindex(columns=columns)


def build_pgls_pipeline(
    args: Any,
    responses: list[str],
    predictors: list[str],
) -> PglsPipelineArtifacts:
    """Run reconciliation, both PIC transforms, and hierarchical PGLS in memory."""
    raw_args = _effective_raw_args(args)
    if raw_args.allow_missing_responses and not raw_args.multivariate_responses:
        raise ValueError(
            "--allow-missing-responses requires '--multivariate-responses yes'."
        )
    gene_tree = read_tree(
        raw_args.gene_tree,
        raw_args.gene_tree_format,
        raw_args.quoted_node_names,
    )
    if raw_args.reconciliation_tree is None:
        reconciliation_tree = gene_tree
    else:
        reconciliation_tree = read_tree(
            raw_args.reconciliation_tree,
            raw_args.reconciliation_tree_format or raw_args.gene_tree_format,
            raw_args.quoted_node_names,
        )
        _validate_matching_gene_topologies(gene_tree, reconciliation_tree)
    species_tree = read_tree(
        raw_args.species_tree,
        raw_args.species_tree_format,
        raw_args.quoted_node_names,
    )
    species_labels = _species_labels(reconciliation_tree, raw_args)
    _report_unmatched_species(
        species_labels,
        species_tree,
        policy=raw_args.unmatched,
    )
    reconciliation = build_reconciliation_table(
        reconciliation_tree,
        species_tree,
        species_labels,
        event_source=raw_args.event_source,
        tree_id=raw_args.tree_id,
    )
    categorical_responses = parse_name_list(raw_args.categorical_responses)
    ordered_responses = parse_ordered_levels(
        raw_args.ordered_responses, "--ordered-responses"
    )
    response_references = parse_key_values(
        raw_args.response_reference, "--response-reference"
    )
    response_families = parse_key_values(raw_args.response_family, "--response-family")
    response_offset_columns = parse_key_values(
        raw_args.response_offset, "--response-offset"
    )
    response_trial_columns = parse_key_values(
        raw_args.response_trials, "--response-trials"
    )
    response_lower_columns = parse_key_values(
        raw_args.response_censor_lower, "--response-censor-lower"
    )
    response_upper_columns = parse_key_values(
        raw_args.response_censor_upper, "--response-censor-upper"
    )
    response_dispersions = parse_numeric_key_values(
        raw_args.response_dispersion, "--response-dispersion", lower=0.0
    )
    response_zero_probabilities = parse_numeric_key_values(
        raw_args.response_zero_probability,
        "--response-zero-probability",
        lower=0.0,
        upper=1.0,
    )
    raw_args.trait = raw_args.expression
    response_inputs = _read_gene_response_inputs(
        raw_args,
        gene_tree,
        responses,
    )
    regular_auxiliary = _read_gene_auxiliary_values(
        raw_args,
        gene_tree,
        sorted(
            set(response_offset_columns.values()) | set(response_trial_columns.values())
        ),
    )
    censor_auxiliary = _read_gene_auxiliary_values(
        raw_args,
        gene_tree,
        sorted(
            set(response_lower_columns.values()) | set(response_upper_columns.values())
        ),
        allow_missing=True,
    )
    response_offsets = _map_response_auxiliary(
        response_offset_columns, regular_auxiliary
    )
    response_trials = _map_response_auxiliary(response_trial_columns, regular_auxiliary)
    response_censor_lower = _map_response_auxiliary(
        response_lower_columns, censor_auxiliary
    )
    response_censor_upper = _map_response_auxiliary(
        response_upper_columns, censor_auxiliary
    )
    response_specs = resolve_response_specs(
        response_inputs.values_by_trait,
        responses,
        [str(leaf.name) for leaf in gene_tree.leaves()],
        categorical=categorical_responses,
        ordered_levels=ordered_responses,
        references=response_references,
        families=response_families,
        allow_missing=(
            raw_args.multivariate_responses and raw_args.allow_missing_responses
        )
        or any(family == "censored-gaussian" for family in response_families.values()),
    )
    validate_response_auxiliaries(
        response_specs,
        offsets=response_offset_columns,
        trials=response_trial_columns,
        censor_lower=response_lower_columns,
        censor_upper=response_upper_columns,
        dispersions=response_dispersions,
        zero_probabilities=response_zero_probabilities,
    )
    if raw_args.coefficient_prior_sd <= 0.0:
        raise ValueError("--coefficient-prior-sd must be positive.")
    continuous_responses, non_gaussian_response_names = (
        _validate_reconciled_response_modes(raw_args, responses, response_specs)
    )
    if not continuous_responses and (
        raw_args.lineage_inference != "none" or raw_args.lineage_leave_one_out
    ):
        raise ValueError(
            "Lineage-slope inference and leave-one-out require a continuous "
            "reconciled-contrast response."
        )
    if raw_args.multivariate_responses and (
        raw_args.lineage_inference != "none" or raw_args.lineage_leave_one_out
    ):
        raise ValueError(
            "Lineage-slope inference and leave-one-out are not available for "
            "the multivariate tip-level response model."
        )
    raw_args.trait = raw_args.species_traits
    raw_predictor_inputs = _read_species_predictor_inputs(
        raw_args,
        species_tree,
        predictors,
    )
    categorical_predictors = parse_name_list(raw_args.categorical_predictors)
    ordered_predictors = parse_ordered_levels(raw_args.ordered_predictors)
    factor_references = parse_key_values(
        raw_args.factor_reference, "--factor-reference"
    )
    encoded_predictors = encode_predictors(
        raw_predictor_inputs.values_by_trait,
        predictors,
        [str(leaf.name) for leaf in species_tree.leaves()],
        categorical=categorical_predictors,
        ordered_levels=ordered_predictors,
        factor_references=factor_references,
        factor_coding=raw_args.factor_coding,
    )
    encoded_sampling_covariance = dict(
        raw_predictor_inputs.sampling_covariance_by_trait or {}
    )
    encoded_replicate_models = dict(raw_predictor_inputs.replicate_model_by_trait or {})
    predictor_inputs = _SpeciesPredictorInputs(
        values_by_trait=encoded_predictors.values_by_trait,
        sampling_covariance_by_trait=(encoded_sampling_covariance or None),
        replicate_model_by_trait=encoded_replicate_models or None,
        tip_summary=raw_predictor_inputs.tip_summary,
    )
    encoded_predictor_names = encoded_predictors.term_names
    (
        species_contrasts,
        predictor_sampling_covariance,
        predictor_diagnostics,
    ) = _build_species_contrasts(
        raw_args,
        species_tree,
        predictor_inputs,
        encoded_predictor_names,
        encoded_predictors.groups,
        {
            uncertainty.source: uncertainty.covariance_by_observation
            for uncertainty in encoded_predictors.uncertainties
        },
    )
    predictor_group_uncertainties, grouped_predictor_covariance_audit = (
        _categorical_contrast_uncertainties(
            species_tree,
            species_contrasts,
            encoded_predictors,
            predictor_diagnostics,
            evolution_model=raw_args.species_evolution_model,
            branch_length=raw_args.species_branch_length,
        )
    )
    predictor_sampling_covariance, predictor_sampling_covariance_for_fit = (
        _merge_predictor_covariance_audit(
            predictor_sampling_covariance, grouped_predictor_covariance_audit
        )
    )
    if raw_args.origin_leave_one_out and not continuous_responses:
        raise ValueError(
            "Origin leave-one-out requires at least one continuous response."
        )
    if (
        raw_args.origin_leave_one_out
        and raw_args.categorical_origin_diagnostics == "none"
    ):
        raise ValueError(
            "'--origin-leave-one-out yes' requires "
            "'--categorical-origin-diagnostics stochastic-map'."
        )
    if raw_args.categorical_origin_diagnostics == "stochastic-map":
        if not categorical_predictors:
            raise ValueError(
                "Categorical origin diagnostics require at least one "
                "--categorical-predictors trait."
            )
        trait_origins, origin_omissions = build_categorical_origin_diagnostics(
            species_tree,
            raw_predictor_inputs.values_by_trait,
            categorical_predictors,
            species_contrasts["branch_clade_id"].astype(str).unique(),
            num_simulations=raw_args.origin_map_replicates,
            minimum_posterior=raw_args.origin_min_posterior,
            seed=raw_args.seed,
            threads=raw_args.origin_map_threads,
        )
    else:
        trait_origins = pd.DataFrame(columns=ORIGIN_DIAGNOSTIC_COLUMNS)
        origin_omissions = []
    (
        tip_encoded_predictors,
        tip_predictor_inputs,
        tip_predictor_diagnostics,
    ) = _prepare_reconciled_tip_predictors(
        raw_args,
        species_tree,
        raw_predictor_inputs,
        predictors,
        categorical_predictors,
        ordered_predictors,
        factor_references,
        predictor_diagnostics,
    )
    if raw_args.multivariate_responses:
        if raw_args.origin_leave_one_out:
            raise ValueError(
                "Origin leave-one-out is not available for multivariate responses."
            )
        if (
            tip_predictor_inputs.sampling_covariance_by_trait
            or tip_predictor_inputs.sparse_posterior_by_trait
        ):
            raise ValueError(
                "Reconciled multivariate PGLS does not yet support predictor "
                "measurement uncertainty; use separate univariate models."
            )
        return _reconciled_multivariate_artifacts(
            raw_args,
            gene_tree,
            species_tree,
            species_labels,
            reconciliation,
            response_inputs,
            responses,
            tip_encoded_predictors,
            tip_predictor_inputs,
            tip_predictor_diagnostics,
            species_contrasts,
            predictor_sampling_covariance,
            predictor_sampling_covariance_for_fit,
            predictor_group_uncertainties,
            trait_origins,
        )
    (
        gene_contrasts,
        sampling_covariance,
        response_diagnostics,
        continuous_results,
        continuous_random_effects,
        sensitivity,
        fit_states,
        refit_shape_in_bootstrap,
    ) = _fit_reconciled_gaussian_responses(
        raw_args,
        gene_tree,
        reconciliation,
        response_inputs,
        continuous_responses,
        species_contrasts,
        predictor_sampling_covariance_for_fit,
        encoded_predictor_names,
        predictor_group_uncertainties,
        encoded_predictors,
        predictor_diagnostics,
        origin_omissions if raw_args.origin_leave_one_out else [],
    )
    if non_gaussian_response_names:
        categorical_results, categorical_random_effects = (
            _fit_reconciled_categorical_responses(
                raw_args,
                gene_tree,
                species_tree,
                species_labels,
                response_inputs,
                response_specs,
                non_gaussian_response_names,
                tip_encoded_predictors,
                tip_predictor_inputs,
                tip_predictor_diagnostics,
                response_offsets,
                response_trials,
                response_censor_lower,
                response_censor_upper,
                response_dispersions,
                response_zero_probabilities,
            )
        )
    else:
        categorical_results = pd.DataFrame(columns=RESULT_COLUMNS)
        categorical_random_effects = pd.DataFrame(columns=RANDOM_EFFECT_COLUMNS)
    results = _merge_reconciled_result_frames(
        continuous_results, categorical_results, RESULT_COLUMNS
    )
    random_effects = _merge_reconciled_result_frames(
        continuous_random_effects,
        categorical_random_effects,
        RANDOM_EFFECT_COLUMNS,
    )
    if refit_shape_in_bootstrap:
        results = _apply_shape_refitted_bootstrap(
            raw_args,
            results,
            gene_tree,
            reconciliation,
            response_inputs,
            continuous_responses,
            species_contrasts,
            predictor_sampling_covariance_for_fit,
            encoded_predictor_names,
            predictor_group_uncertainties,
            response_diagnostics,
            fit_states,
        )
    return PglsPipelineArtifacts(
        reconciliation=reconciliation,
        gene_contrasts=gene_contrasts,
        species_contrasts=species_contrasts,
        response_sampling_covariance=sampling_covariance,
        response_tip_summary=response_inputs.tip_summary,
        predictor_sampling_covariance=predictor_sampling_covariance,
        predictor_tip_summary=predictor_inputs.tip_summary,
        results=results,
        random_effects=random_effects,
        sensitivity=sensitivity,
        trait_origins=trait_origins,
    )


def build_pgls_ensemble_pipeline(
    args: Any,
    responses: list[str],
    predictors: list[str],
) -> PglsPipelineArtifacts:
    """Fit a Newick sample and combine within- and between-tree uncertainty."""
    tree_strings = read_tree_strings(args.gene_tree_ensemble)
    if len(tree_strings) < 2:
        raise ValueError("--gene-tree-ensemble requires at least two Newick trees.")
    artifacts = []
    for tree_index, tree_string in enumerate(tree_strings, start=1):
        values = vars(args).copy()
        values["gene_tree"] = tree_string
        values["gene_tree_ensemble"] = None
        values["tree_id"] = "{}#{}".format(args.tree_id, tree_index)
        artifacts.append(
            build_pgls_pipeline(SimpleNamespace(**values), responses, predictors)
        )
    results = _combine_ensemble_results(
        [artifact.results for artifact in artifacts],
        str(args.tree_id),
        args.confidence_level,
    )
    return PglsPipelineArtifacts(
        reconciliation=pd.concat(
            [artifact.reconciliation for artifact in artifacts], ignore_index=True
        ),
        gene_contrasts=pd.concat(
            [artifact.gene_contrasts for artifact in artifacts], ignore_index=True
        ),
        species_contrasts=artifacts[0].species_contrasts,
        response_sampling_covariance=_concat_optional_frames(
            [artifact.response_sampling_covariance for artifact in artifacts]
        ),
        response_tip_summary=artifacts[0].response_tip_summary,
        predictor_sampling_covariance=artifacts[0].predictor_sampling_covariance,
        predictor_tip_summary=artifacts[0].predictor_tip_summary,
        results=results,
        random_effects=pd.concat(
            [artifact.random_effects for artifact in artifacts], ignore_index=True
        ),
        sensitivity=_concat_optional_frames(
            [artifact.sensitivity for artifact in artifacts]
        ),
        trait_origins=artifacts[0].trait_origins,
    )


def write_pgls_bundle(
    prefix: str,
    artifacts: PglsPipelineArtifacts,
) -> dict[str, str]:
    """Transactionally write an end-to-end bundle and return committed paths."""
    validate_pgls_bundle_target(prefix)
    written = _active_pgls_bundle_paths(prefix, artifacts)
    validate_distinct_output_paths(
        [
            ("--out-prefix {}".format(name.replace("_", " ")), path)
            for name, path in written.items()
        ]
    )
    frames = {
        "reconciliation_out": artifacts.reconciliation,
        "gene_contrasts_out": artifacts.gene_contrasts,
        "species_contrasts_out": artifacts.species_contrasts,
        "response_sampling_covariance_out": artifacts.response_sampling_covariance,
        "response_tip_summary_out": artifacts.response_tip_summary,
        "predictor_sampling_covariance_out": (artifacts.predictor_sampling_covariance),
        "predictor_tip_summary_out": artifacts.predictor_tip_summary,
        "random_effects_out": artifacts.random_effects,
        "sensitivity_out": artifacts.sensitivity,
        "trait_origins_out": artifacts.trait_origins,
        "outfile": artifacts.results,
    }
    with acquire_exclusive_lock(
        pgls_bundle_lock_path(prefix),
        lock_label="PGLS output bundle",
    ):
        _write_dataframes_transactionally(
            [(path, frames[name]) for name, path in written.items()]
        )
    return written
