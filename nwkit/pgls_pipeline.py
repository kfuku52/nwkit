"""End-to-end reconciliation, contrast, and PGLS orchestration."""

import math
import os
import secrets
import shutil
import stat
import tempfile
from dataclasses import dataclass
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd

from nwkit.clade_index import CladeIndex
from nwkit.contrast import (
    _orient_children,
    _read_numeric_traits,
    _read_replicate_traits,
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
    evolution_model_spec,
    optimization_parameterization,
    parameter_near_boundary,
)
from nwkit.ordinary_pgls import (
    _global_bounded_scalar_minimize,
    estimate_marginal_evolution_parameter,
)
from nwkit.pgls import fit_reconciled_pgls
from nwkit.reconcile import _report_unmatched_species, build_reconciliation_table
from nwkit.species_parser import get_species_parser
from nwkit.util import (
    acquire_exclusive_lock,
    normalized_missing_path_key,
    read_tree,
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
    predictor_sampling_covariance: pd.DataFrame | None = None
    predictor_tip_summary: pd.DataFrame | None = None


@dataclass
class _GeneResponseInputs:
    values_by_trait: dict[str, dict[str, float]]
    sampling_covariance_by_trait: dict[str, pd.DataFrame] | None
    replicate_model_by_trait: dict[str, str] | None
    tip_summary: pd.DataFrame | None


@dataclass
class _SpeciesPredictorInputs:
    values_by_trait: dict[str, dict[str, float]]
    sampling_covariance_by_trait: dict[str, pd.DataFrame] | None
    replicate_model_by_trait: dict[str, str] | None
    tip_summary: pd.DataFrame | None


def _active_pgls_bundle_paths(
    prefix: str,
    artifacts: PglsPipelineArtifacts,
) -> dict[str, str]:
    paths = pgls_bundle_paths(prefix)
    inactive = set()
    if artifacts.response_sampling_covariance is None:
        if artifacts.response_tip_summary is not None:
            raise RuntimeError(
                "A response tip summary requires response sampling covariance."
            )
        inactive.update(
            {"response_sampling_covariance_out", "response_tip_summary_out"}
        )
    elif artifacts.response_tip_summary is None:
        raise RuntimeError(
            "Replicate-aware covariance requires a response tip summary."
        )
    if artifacts.predictor_sampling_covariance is None:
        if artifacts.predictor_tip_summary is not None:
            raise RuntimeError(
                "A predictor tip summary requires predictor sampling covariance."
            )
        inactive.update(
            {"predictor_sampling_covariance_out", "predictor_tip_summary_out"}
        )
    elif artifacts.predictor_tip_summary is None:
        raise RuntimeError(
            "Replicate-aware predictor covariance requires a predictor tip summary."
        )
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


def _write_dataframes_transactionally(
    outputs: list[tuple[str, pd.DataFrame]],
) -> None:
    output_modes: dict[str, int] = {}
    for path, _ in outputs:
        absolute_path = os.path.abspath(path)
        directory = os.path.dirname(absolute_path)
        mode = _regular_output_mode(absolute_path)
        output_modes[absolute_path] = (
            _new_output_mode(directory) if mode is None else mode
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
) -> _GeneResponseInputs:
    replicate_requested = _validate_replicate_options(args)
    replicate_estimates = (
        _read_replicate_traits(
            args,
            gene_tree,
            responses,
            args.tree_id,
            option_name="--expression",
        )
        if replicate_requested
        else None
    )
    values_by_trait = (
        replicate_estimates.values_by_trait
        if replicate_estimates is not None
        else _read_numeric_traits(
            args, gene_tree, responses, option_name="--expression"
        )
    )
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
        tip_summary=(
            None if replicate_estimates is None else replicate_estimates.tip_summary
        ),
    )


def _read_species_predictor_inputs(
    args: Any,
    species_tree,
    predictors: list[str],
) -> _SpeciesPredictorInputs:
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
    replicate_estimates = (
        _read_replicate_traits(
            predictor_args,
            species_tree,
            predictors,
            "species",
            option_name="--species-traits",
        )
        if replicate_requested
        else None
    )
    values_by_trait = (
        replicate_estimates.values_by_trait
        if replicate_estimates is not None
        else _read_numeric_traits(
            predictor_args,
            species_tree,
            predictors,
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
        inference="wald",
        bootstrap_replicates=2,
        seed=args.seed,
        reml=args.reml,
        event_random_effect=args.event_random_effect,
        lineage_random_slope=args.lineage_random_slope,
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
    for predictor in predictors:
        parameter: float | None
        if auto_parameter:
            diagnostics = estimate_marginal_evolution_parameter(
                species_tree,
                predictor_inputs.values_by_trait[predictor],
                predictor,
                evolution_model=args.species_evolution_model,
                branch_length=args.species_branch_length,
                sampling_covariance=(
                    None if sampling_by_trait is None else sampling_by_trait[predictor]
                ),
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
                None
                if sampling_by_trait is None
                else {predictor: sampling_by_trait[predictor]}
            ),
            replicate_model_by_trait=(
                None
                if replicate_models is None
                else {predictor: replicate_models[predictor]}
            ),
            tip_summary=(
                None
                if predictor_inputs.tip_summary is None
                else predictor_inputs.tip_summary[
                    predictor_inputs.tip_summary["trait"] == predictor
                ].copy()
            ),
            return_sampling_covariance=sampling_by_trait is not None,
            tree_option_name="--species-tree",
        )
        if sampling_by_trait is None:
            frames.append(output)
        else:
            contrasts, covariance = output
            frames.append(contrasts)
            covariance_frames.append(covariance)
        diagnostics_by_predictor[predictor] = diagnostics
    covariance_table = (
        None
        if sampling_by_trait is None
        else pd.concat(covariance_frames, ignore_index=True)
    )
    return (
        pd.concat(frames, ignore_index=True),
        covariance_table,
        diagnostics_by_predictor,
    )


def _yes_no(value: bool) -> str:
    return "yes" if value else "no"


def _diagnostic_flag(value: bool | None) -> str:
    return "not-applicable" if value is None else _yes_no(value)


def _attach_evolution_diagnostics(
    results: pd.DataFrame,
    response_diagnostics: dict[str, dict[str, Any]],
    predictor_diagnostics: dict[str, dict[str, Any]],
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
        predictor = predictor_diagnostics[str(row["term"])]
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


def _selected_response_contrast_matrix(
    args: Any,
    gene_tree,
    reconciliation: pd.DataFrame,
    response_inputs: _GeneResponseInputs,
    response: str,
    parameter: float,
    contrast_ids: list[str],
) -> tuple[np.ndarray, list[str]]:
    clades = CladeIndex(gene_tree)
    reconciliation_by_id = {
        str(row["gene_clade_id"]): row for row in reconciliation.to_dict("records")
    }
    orientation = _orient_children(gene_tree, clades, reconciliation_by_id)
    _, coefficients_by_node, leaf_names = calculate_contrasts(
        gene_tree,
        response_inputs.values_by_trait[response],
        branch_length=args.gene_branch_length,
        evolution_model=args.gene_evolution_model,
        evolution_parameter=parameter,
        orientation_by_node=orientation,
        return_coefficients=True,
    )
    coefficient_by_id = {
        clades.clade_id_for_node(node): coefficients
        for node, coefficients in coefficients_by_node.items()
    }
    try:
        matrix = np.vstack(
            [coefficient_by_id[contrast_id] for contrast_id in contrast_ids]
        )
    except KeyError as exc:
        raise ValueError(
            "Bootstrap fit state references an unknown response contrast."
        ) from exc
    if int(np.linalg.matrix_rank(matrix)) != len(contrast_ids):
        raise ValueError("Selected response contrasts are not linearly independent.")
    return matrix, leaf_names


def _simulate_response_tip_values(
    args: Any,
    gene_tree,
    response_inputs: _GeneResponseInputs,
    response: str,
    diagnostics: dict[str, Any],
    fit_state: dict[str, Any],
    contrast_matrix: np.ndarray,
    leaf_names: list[str],
    rng: np.random.Generator,
) -> dict[str, float]:
    parameter = float(diagnostics["parameter"])
    phylogenetic_covariance = build_evolutionary_covariance(
        gene_tree,
        leaf_names,
        model=args.gene_evolution_model,
        parameter=parameter,
        branch_length=args.gene_branch_length,
    )
    if response_inputs.sampling_covariance_by_trait is None:
        tip_sampling_covariance = np.zeros_like(phylogenetic_covariance)
    else:
        tip_sampling_covariance = (
            response_inputs.sampling_covariance_by_trait[response]
            .reindex(index=leaf_names, columns=leaf_names)
            .to_numpy(dtype=float)
        )
    tip_covariance = (
        float(fit_state["evolutionary_rate"]) * phylogenetic_covariance
        + tip_sampling_covariance
    )
    tip_factor = _positive_semidefinite_factor(
        tip_covariance, "Bootstrap tip covariance"
    )
    contrast_factor = _positive_semidefinite_factor(
        fit_state["fitted_covariance"],
        "Bootstrap fitted contrast covariance",
    )
    contrast_inverse = np.linalg.pinv(contrast_matrix)
    if not np.allclose(
        contrast_matrix @ contrast_inverse,
        np.eye(len(contrast_matrix)),
        rtol=1e-8,
        atol=1e-10,
    ):
        raise ValueError("Selected response contrasts cannot be mapped back to tips.")
    contrast_mean = fit_state["design"] @ fit_state["beta"]
    contrast_error = contrast_factor @ rng.standard_normal(contrast_factor.shape[1])
    null_projection = np.eye(len(leaf_names)) - contrast_inverse @ contrast_matrix
    null_error = null_projection @ (
        tip_factor @ rng.standard_normal(tip_factor.shape[1])
    )
    simulated = contrast_inverse @ (contrast_mean + contrast_error) + null_error
    if not np.allclose(
        contrast_matrix @ simulated,
        contrast_mean + contrast_error,
        rtol=1e-8,
        atol=1e-10,
    ):
        raise RuntimeError(
            "Bootstrap tip simulation did not preserve fitted contrasts."
        )
    return dict(zip(leaf_names, simulated, strict=True))


def _bootstrap_shape_refitted_coefficients(
    args: Any,
    gene_tree,
    reconciliation: pd.DataFrame,
    response_inputs: _GeneResponseInputs,
    response: str,
    species_contrasts: pd.DataFrame,
    predictor_sampling_covariance: pd.DataFrame | None,
    predictors: list[str],
    diagnostics: dict[str, Any],
    fit_state: dict[str, Any],
    seed: int,
) -> np.ndarray:
    parameter = float(diagnostics["parameter"])
    contrast_matrix, leaf_names = _selected_response_contrast_matrix(
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
            args,
            gene_tree,
            response_inputs,
            response,
            diagnostics,
            fit_state,
            contrast_matrix,
            leaf_names,
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
            )
            bootstrap_result = _fit_candidate_reconciled_model(
                args,
                contrasts,
                species_contrasts,
                response,
                predictors,
                covariance,
                predictor_sampling_covariance,
            ).set_index("term")
        except ValueError:
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


def build_pgls_pipeline(
    args: Any,
    responses: list[str],
    predictors: list[str],
) -> PglsPipelineArtifacts:
    """Run reconciliation, both PIC transforms, and hierarchical PGLS in memory."""
    raw_args = _effective_raw_args(args)
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
    raw_args.trait = raw_args.expression
    response_inputs = _read_gene_response_inputs(
        raw_args,
        gene_tree,
        responses,
    )
    raw_args.trait = raw_args.species_traits
    predictor_inputs = _read_species_predictor_inputs(
        raw_args,
        species_tree,
        predictors,
    )
    (
        species_contrasts,
        predictor_sampling_covariance,
        predictor_diagnostics,
    ) = _build_species_contrasts(
        raw_args,
        species_tree,
        predictor_inputs,
        predictors,
    )
    gene_contrasts, sampling_covariance, response_diagnostics = _build_gene_contrasts(
        raw_args,
        gene_tree,
        reconciliation,
        response_inputs,
        responses,
        species_contrasts,
        predictor_sampling_covariance,
        predictors,
    )
    refit_shape_in_bootstrap = raw_args.inference == "parametric-bootstrap" and any(
        diagnostics["parameter_status"] == "estimated"
        for diagnostics in response_diagnostics.values()
    )
    fitted = fit_reconciled_pgls(
        gene_contrasts,
        species_contrasts,
        responses,
        predictors,
        confidence_level=raw_args.confidence_level,
        event_weighting=raw_args.event_weighting,
        coverage_policy=raw_args.speciation_coverage,
        model=raw_args.model,
        response_sampling_covariance=sampling_covariance,
        predictor_sampling_covariance=predictor_sampling_covariance,
        inference="wald" if refit_shape_in_bootstrap else raw_args.inference,
        bootstrap_replicates=raw_args.bootstrap_replicates,
        seed=raw_args.seed,
        reml=raw_args.reml,
        event_random_effect=raw_args.event_random_effect,
        lineage_random_slope=raw_args.lineage_random_slope,
        return_random_effects=True,
        return_fit_state=refit_shape_in_bootstrap,
    )
    if refit_shape_in_bootstrap:
        results, random_effects, fit_states = fitted
    else:
        results, random_effects = fitted
        fit_states = {}
    results = _attach_evolution_diagnostics(
        results,
        response_diagnostics,
        predictor_diagnostics,
    )
    if refit_shape_in_bootstrap:
        results = _apply_shape_refitted_bootstrap(
            raw_args,
            results,
            gene_tree,
            reconciliation,
            response_inputs,
            responses,
            species_contrasts,
            predictor_sampling_covariance,
            predictors,
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
