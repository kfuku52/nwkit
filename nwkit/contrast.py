import csv
import math
from io import StringIO

import numpy as np
import pandas as pd
from scipy import sparse

from nwkit.clade_index import CladeIndex
from nwkit.evolution import (
    evolution_model_spec,
    transformed_edge_variances,
    validate_evolution_parameter,
)
from nwkit.gaussian import DiagonalLowRankCovariance
from nwkit.util import (
    assign_branch_ids,
    get_node_class,
    is_rooted,
    read_input_text,
    read_tip_table,
    read_tree,
    validate_distinct_output_paths,
    validate_outputs_do_not_replace_inputs,
    validate_unique_named_leaves,
)

BASE_CONTRAST_COLUMNS = [
    "tree_id",
    "branch_id",
    "branch_clade_id",
    "node_class",
    "descendant_taxa",
    "num_taxa",
    "trait",
    "evolution_model",
    "evolution_parameter_name",
    "evolution_parameter",
    "branch_length_mode",
    "numerator_branch_id",
    "numerator_clade_id",
    "denominator_branch_id",
    "denominator_clade_id",
    "numerator_taxa",
    "denominator_taxa",
    "raw_contrast",
    "contrast_variance",
    "standardized_contrast",
    "ancestral_estimate",
]

RECONCILIATION_CONTEXT_COLUMNS = [
    "gene_branch_id",
    "gene_clade_id",
    "species_branch_id",
    "species_event_id",
    "species_numerator_branch_id",
    "species_numerator_event_id",
    "species_denominator_branch_id",
    "species_denominator_event_id",
    "lineage_id",
    "lineage_clade_id",
    "event_type",
    "event_source",
    "event_status",
    "collapsed_event_boundary",
    "transfer_source_species",
    "transfer_destination_species",
    "mapping_status",
    "eligible",
    "coverage_status",
    "reason",
    "descendant_species_taxa",
    "num_descendant_species_taxa",
    "species_event_taxa",
    "num_species_event_taxa",
    "numerator_observed_species_taxa",
    "num_numerator_observed_species_taxa",
    "numerator_species_clade_taxa",
    "num_numerator_species_clade_taxa",
    "numerator_species_coverage",
    "denominator_observed_species_taxa",
    "num_denominator_observed_species_taxa",
    "denominator_species_clade_taxa",
    "num_denominator_species_clade_taxa",
    "denominator_species_coverage",
]

RECONCILIATION_REQUIRED_COLUMNS = set(RECONCILIATION_CONTEXT_COLUMNS).union(
    {
        "tree_id",
        "contrast_numerator_gene_clade_id",
        "contrast_denominator_gene_clade_id",
    }
)

RECONCILIATION_ENUMS = {
    "event_type": {"leaf", "speciation", "duplication", "transfer", "unresolved"},
    "event_source": {"lca", "nhx", "species-overlap"},
    "event_status": {"resolved", "unresolved", "not-applicable"},
    "mapping_status": {"mapped", "unmapped"},
    "eligible": {"yes", "no"},
    "coverage_status": {"complete", "partial", "not-applicable"},
    "collapsed_event_boundary": {"yes", "no"},
}

COVERAGE_SIDES = ("numerator", "denominator")

REPLICATE_CONTRAST_COLUMNS = [
    "sampling_variance",
    "evolutionary_variance",
    "replicate_model",
    "min_n_biological",
]

SAMPLING_COVARIANCE_COLUMNS = [
    "tree_id",
    "trait",
    "contrast_id_1",
    "contrast_id_2",
    "sampling_covariance",
    "covariance_representation",
]

MAX_FULL_SAMPLING_COVARIANCE_CONTRASTS = 500


def _parse_columns(value):
    if value in (None, ""):
        raise ValueError("'--columns' is required.")
    columns = [column.strip() for column in str(value).split(",")]
    if any(column == "" for column in columns):
        raise ValueError("'--columns' contains an empty column name.")
    if len(columns) != len(set(columns)):
        raise ValueError("'--columns' contains duplicated column names.")
    if "leaf_name" in columns:
        raise ValueError("'leaf_name' cannot be used as a numeric trait column.")
    return columns


def _validate_contrast_tree(tree, branch_length, option_name="--infile"):
    validate_unique_named_leaves(
        tree, option_name=option_name, context=" for 'contrast'"
    )
    if len(list(tree.leaves())) < 2:
        raise ValueError(
            "'{}' must contain at least two tips for contrasts.".format(option_name)
        )
    if not is_rooted(tree):
        raise ValueError("Contrasts require a rooted tree.")
    nonbinary = [
        node for node in tree.traverse() if not node.is_leaf and len(node.children) != 2
    ]
    if nonbinary:
        raise ValueError(
            "Contrasts require a strictly bifurcating tree; found {} non-binary internal node(s).".format(
                len(nonbinary)
            )
        )
    if branch_length == "unit":
        return
    for node in tree.traverse():
        if node.is_root:
            continue
        if node.dist is None:
            raise ValueError(
                "Contrasts require branch lengths on every non-root branch; "
                "use '--branch-length unit' to ignore input lengths."
            )
        try:
            distance = float(node.dist)
        except (TypeError, ValueError) as exc:
            raise ValueError("Tree branch lengths must be numeric.") from exc
        if not math.isfinite(distance) or distance <= 0.0:
            raise ValueError(
                "Contrasts require positive finite branch lengths; "
                "use '--branch-length unit' to ignore input lengths."
            )


def _read_numeric_traits(args, tree, columns, option_name="--trait"):
    dataframe, _, _ = read_tip_table(
        args.trait,
        option_name=option_name,
        tree_leaf_names=list(tree.leaf_names()),
        required_columns=columns,
        unmatched=args.unmatched,
        missing_values=args.missing_values,
    )
    dataframe = dataframe[dataframe["leaf_name"].isin(set(tree.leaf_names()))]
    row_by_leaf = dataframe.set_index("leaf_name", drop=False).to_dict("index")
    values_by_trait = dict()
    for column in columns:
        values_by_leaf = dict()
        missing = list()
        nonnumeric = list()
        for leaf_name in tree.leaf_names():
            row = row_by_leaf.get(str(leaf_name))
            value = None if row is None else row[column]
            if value is None or pd.isna(value):
                missing.append(str(leaf_name))
                continue
            try:
                numeric_value = float(value)
            except (TypeError, ValueError):
                nonnumeric.append(str(leaf_name))
                continue
            if not math.isfinite(numeric_value):
                nonnumeric.append(str(leaf_name))
                continue
            values_by_leaf[str(leaf_name)] = numeric_value
        if missing:
            raise ValueError(
                "Trait column '{}' has missing values for tree tips: {}.".format(
                    column, ", ".join(sorted(missing))
                )
            )
        if nonnumeric:
            raise ValueError(
                "Trait column '{}' has non-numeric or non-finite values for tree tips: {}.".format(
                    column, ", ".join(sorted(nonnumeric))
                )
            )
        values_by_trait[column] = values_by_leaf
    return values_by_trait


def _read_typed_traits(
    args,
    tree,
    columns,
    *,
    categorical=(),
    ordered=(),
    allow_missing=(),
    option_name="--trait",
):
    """Read continuous or categorical tip traits with global type inference."""
    dataframe, _, _ = read_tip_table(
        args.trait,
        option_name=option_name,
        tree_leaf_names=list(tree.leaf_names()),
        required_columns=columns,
        unmatched=args.unmatched,
        missing_values=args.missing_values,
    )
    dataframe = dataframe[dataframe["leaf_name"].isin(set(tree.leaf_names()))]
    row_by_leaf = dataframe.set_index("leaf_name", drop=False).to_dict("index")
    discrete = set(categorical) | set(ordered)
    allow_missing = set(allow_missing)
    values_by_trait = {}
    for column in columns:
        raw_values = [
            None
            if row_by_leaf.get(str(name)) is None
            else row_by_leaf[str(name)][column]
            for name in tree.leaf_names()
        ]
        missing = [
            str(name)
            for name, value in zip(tree.leaf_names(), raw_values, strict=True)
            if value is None or pd.isna(value)
        ]
        if missing and column not in allow_missing:
            raise ValueError(
                "Trait column '{}' has missing values for tree tips: {}.".format(
                    column, ", ".join(sorted(missing))
                )
            )
        numeric = pd.to_numeric(pd.Series(raw_values), errors="coerce")
        categorical_column = column in discrete or (
            column not in allow_missing and numeric.isna().any()
        )
        if categorical_column and missing:
            raise ValueError(
                "Categorical trait '{}' cannot contain missing values.".format(column)
            )
        numeric_values = numeric.to_numpy(float)
        if not categorical_column:
            invalid = (
                np.isinf(numeric_values)
                if column in allow_missing
                else ~np.isfinite(numeric_values)
            )
            if invalid.any():
                raise ValueError(
                    "Trait column '{}' has non-finite values.".format(column)
                )
        values_by_trait[column] = {
            str(name): str(value) if categorical_column else float(number)
            for name, value, number in zip(
                tree.leaf_names(), raw_values, numeric, strict=True
            )
        }
    return values_by_trait


def _parse_optional_columns(value, option_name, expected_length):
    if value in (None, ""):
        return []
    columns = [column.strip() for column in str(value).split(",")]
    if any(column == "" for column in columns):
        raise ValueError("'{}' contains an empty column name.".format(option_name))
    if len(columns) != expected_length:
        raise ValueError(
            "'{}' must contain exactly {} column name(s).".format(
                option_name, expected_length
            )
        )
    if len(columns) != len(set(columns)):
        raise ValueError("'{}' contains duplicated column names.".format(option_name))
    return columns


def _validate_known_se_options(args):
    biological_id = getattr(args, "biological_id", None)
    technical_id = getattr(args, "technical_id", None)
    batch = getattr(args, "batch", None)
    technical_aggregation = getattr(args, "technical_aggregation", "error")
    incompatible = [
        option
        for option, value in [
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


def _validate_raw_replicate_options(args):
    technical_id = getattr(args, "technical_id", None)
    technical_aggregation = getattr(args, "technical_aggregation", "error")
    se_columns = getattr(args, "standard_error_columns", None)
    n_columns = getattr(args, "sample_size_columns", None)
    if se_columns is not None or n_columns is not None:
        raise ValueError(
            "Standard-error and sample-size columns require "
            "'--within-variance known-se'."
        )
    if technical_aggregation != "error" and technical_id is None:
        raise ValueError("'--technical-aggregation' requires '--technical-id'.")


def _validate_nonreplicate_options(args):
    technical_id = getattr(args, "technical_id", None)
    batch = getattr(args, "batch", None)
    within_variance = getattr(args, "within_variance", "pooled")
    technical_aggregation = getattr(args, "technical_aggregation", "error")
    se_columns = getattr(args, "standard_error_columns", None)
    n_columns = getattr(args, "sample_size_columns", None)
    covariance_out = getattr(args, "sampling_covariance_out", None)
    tip_summary_out = getattr(args, "tip_summary_out", None)
    incompatible = [
        option
        for option, enabled in [
            ("--technical-id", technical_id is not None),
            ("--batch", batch is not None),
            ("--within-variance leaf", within_variance == "leaf"),
            ("--technical-aggregation", technical_aggregation != "error"),
            ("--standard-error-columns", se_columns is not None),
            ("--sample-size-columns", n_columns is not None),
            ("--sampling-covariance-out", covariance_out is not None),
            ("--tip-summary-out", tip_summary_out is not None),
        ]
        if enabled
    ]
    if incompatible:
        raise ValueError(
            "Replicate option(s) require '--biological-id' or "
            "'--within-variance known-se': {}.".format(", ".join(incompatible))
        )


def _validate_replicate_options(args):
    biological_id = getattr(args, "biological_id", None)
    known_se = getattr(args, "within_variance", "pooled") == "known-se"
    if known_se:
        _validate_known_se_options(args)
        return True
    if biological_id is not None:
        _validate_raw_replicate_options(args)
        return True
    _validate_nonreplicate_options(args)
    return False


def _read_replicate_traits(
    args,
    tree,
    columns,
    tree_id,
    option_name="--trait",
    *,
    allow_missing_columns=(),
):
    from nwkit.replicates import estimate_replicate_traits

    se_columns = _parse_optional_columns(
        getattr(args, "standard_error_columns", None),
        "--standard-error-columns",
        len(columns),
    )
    n_columns = _parse_optional_columns(
        getattr(args, "sample_size_columns", None),
        "--sample-size-columns",
        len(columns),
    )
    within_variance = getattr(args, "within_variance", "pooled")
    if within_variance == "known-se" and not se_columns:
        se_columns = ["{}_se".format(column) for column in columns]
    required_columns = list(columns) + list(se_columns) + list(n_columns)
    for optional_column in [
        getattr(args, "biological_id", None),
        getattr(args, "technical_id", None),
        getattr(args, "batch", None),
    ]:
        if optional_column is not None:
            required_columns.append(optional_column)
    dataframe, _, _ = read_tip_table(
        args.trait,
        option_name=option_name,
        tree_leaf_names=list(tree.leaf_names()),
        required_columns=required_columns,
        unmatched=args.unmatched,
        missing_values=args.missing_values,
        duplicate_leaf_names="allow",
    )
    return estimate_replicate_traits(
        dataframe,
        list(tree.leaf_names()),
        columns,
        biological_id=getattr(args, "biological_id", None),
        technical_id=getattr(args, "technical_id", None),
        batch=getattr(args, "batch", None),
        within_variance=within_variance,
        technical_aggregation=getattr(args, "technical_aggregation", "error"),
        se_columns=se_columns,
        n_columns=n_columns,
        allow_missing_traits=allow_missing_columns,
        tree_id=tree_id,
    )


def _read_mixed_replicate_traits(
    args,
    tree,
    columns,
    categorical_columns,
    tree_id,
    *,
    categorical_policy="error",
    non_gaussian_columns=(),
    allow_missing_columns=(),
    option_name="--trait",
):
    """Read a mixture of continuous and categorical replicate traits."""
    from nwkit.replicates import (
        ReplicateEstimates,
        estimate_categorical_traits,
        estimate_likelihood_replicates,
    )

    categorical_columns = set(categorical_columns)
    non_gaussian_columns = set(non_gaussian_columns)
    allow_missing_columns = set(allow_missing_columns)
    inspection, _, _ = read_tip_table(
        args.trait,
        option_name=option_name,
        tree_leaf_names=list(tree.leaf_names()),
        required_columns=list(columns),
        unmatched=args.unmatched,
        missing_values=args.missing_values,
        duplicate_leaf_names="allow",
    )
    for column in columns:
        observed = inspection[column].dropna()
        numeric = pd.to_numeric(observed, errors="coerce")
        if numeric.isna().any():
            categorical_columns.add(column)
    overlap = categorical_columns & non_gaussian_columns
    if overlap:
        raise ValueError(
            "Traits cannot be both categorical and scalar non-Gaussian: {}.".format(
                ", ".join(sorted(overlap))
            )
        )
    continuous_columns = [
        column
        for column in columns
        if column not in categorical_columns and column not in non_gaussian_columns
    ]
    estimates = []
    if continuous_columns:
        estimates.append(
            _read_replicate_traits(
                args,
                tree,
                continuous_columns,
                tree_id,
                option_name=option_name,
                allow_missing_columns=allow_missing_columns & set(continuous_columns),
            )
        )
    if categorical_columns:
        if getattr(args, "within_variance", "pooled") == "known-se":
            raise ValueError(
                "Known standard errors do not apply to categorical traits."
            )
        if getattr(args, "batch", None) is not None:
            raise ValueError(
                "Batch adjustment for categorical traits is not supported."
            )
        biological_id = getattr(args, "biological_id", None)
        if biological_id is None:
            raise ValueError("Categorical replicate input requires '--biological-id'.")
        required_columns = list(categorical_columns) + [biological_id]
        technical_id = getattr(args, "technical_id", None)
        if technical_id is not None:
            required_columns.append(technical_id)
        dataframe, _, _ = read_tip_table(
            args.trait,
            option_name=option_name,
            tree_leaf_names=list(tree.leaf_names()),
            required_columns=required_columns,
            unmatched=args.unmatched,
            missing_values=args.missing_values,
            duplicate_leaf_names="allow",
        )
        estimates.append(
            estimate_categorical_traits(
                dataframe,
                list(tree.leaf_names()),
                [column for column in columns if column in categorical_columns],
                biological_id=biological_id,
                technical_id=technical_id,
                policy=categorical_policy,
                tree_id=tree_id,
            )
        )
    if non_gaussian_columns:
        if getattr(args, "within_variance", "pooled") == "known-se":
            raise ValueError(
                "Known standard errors do not apply to non-Gaussian likelihood replicates."
            )
        if getattr(args, "batch", None) is not None:
            raise ValueError(
                "Batch adjustment for non-Gaussian likelihood replicates is not yet supported."
            )
        biological_id = getattr(args, "biological_id", None)
        if biological_id is None:
            raise ValueError("Non-Gaussian replicate input requires '--biological-id'.")
        required_columns = list(non_gaussian_columns) + [biological_id]
        technical_id = getattr(args, "technical_id", None)
        if technical_id is not None:
            required_columns.append(technical_id)
        dataframe, _, _ = read_tip_table(
            args.trait,
            option_name=option_name,
            tree_leaf_names=list(tree.leaf_names()),
            required_columns=required_columns,
            unmatched=args.unmatched,
            missing_values=args.missing_values,
            duplicate_leaf_names="allow",
        )
        estimates.append(
            estimate_likelihood_replicates(
                dataframe,
                list(tree.leaf_names()),
                [column for column in columns if column in non_gaussian_columns],
                biological_id=biological_id,
                technical_id=technical_id,
                technical_aggregation=getattr(args, "technical_aggregation", "error"),
                allow_missing_traits=(allow_missing_columns & non_gaussian_columns),
                tree_id=tree_id,
            )
        )
    return ReplicateEstimates(
        values_by_trait={
            trait: values
            for estimate in estimates
            for trait, values in estimate.values_by_trait.items()
        },
        sampling_covariance_by_trait={
            trait: covariance
            for estimate in estimates
            for trait, covariance in estimate.sampling_covariance_by_trait.items()
        },
        tip_summary=pd.concat(
            [estimate.tip_summary for estimate in estimates], ignore_index=True
        ),
        model_by_trait={
            trait: model
            for estimate in estimates
            for trait, model in estimate.model_by_trait.items()
        },
    )


def _validate_reconciliation_domains(dataframe):
    for column, allowed_values in RECONCILIATION_ENUMS.items():
        invalid = sorted(set(dataframe[column]) - allowed_values)
        if invalid:
            raise ValueError(
                "'--reconciliation' contains invalid {} value(s): {}.".format(
                    column, ", ".join(invalid)
                )
            )
    mapped_missing_event = dataframe[
        (dataframe["mapping_status"] == "mapped")
        & (
            (dataframe["species_event_id"] == "")
            | (dataframe["species_branch_id"] == "")
        )
    ]
    if not mapped_missing_event.empty:
        raise ValueError(
            "'--reconciliation' mapped rows must define species_event_id and species_branch_id."
        )
    unmapped_with_event = dataframe[
        (dataframe["mapping_status"] == "unmapped")
        & (
            (dataframe["species_event_id"] != "")
            | (dataframe["species_branch_id"] != "")
        )
    ]
    if not unmapped_with_event.empty:
        raise ValueError(
            "'--reconciliation' unmapped rows cannot define species event identifiers."
        )
    invalid_eligible = dataframe[
        (dataframe["eligible"] == "yes")
        & (
            (dataframe["event_type"] != "speciation")
            | (dataframe["event_status"] != "resolved")
            | (dataframe["mapping_status"] != "mapped")
            | ~dataframe["coverage_status"].isin(["complete", "partial"])
        )
    ]
    if not invalid_eligible.empty:
        raise ValueError(
            "'--reconciliation' eligible=yes is only valid for resolved, mapped speciation rows."
        )
    invalid_status = dataframe[
        (
            (dataframe["event_type"] == "leaf")
            & (dataframe["event_status"] != "not-applicable")
        )
        | (
            (dataframe["event_type"] == "unresolved")
            & (dataframe["event_status"] != "unresolved")
        )
        | (
            dataframe["event_type"].isin(["speciation", "duplication", "transfer"])
            & (dataframe["event_status"] != "resolved")
        )
    ]
    if not invalid_status.empty:
        raise ValueError(
            "'--reconciliation' event_type and event_status values are inconsistent."
        )
    invalid_coverage = dataframe[
        (
            (dataframe["eligible"] == "yes")
            & (dataframe["coverage_status"] == "not-applicable")
        )
        | (
            (dataframe["eligible"] == "no")
            & (dataframe["coverage_status"] != "not-applicable")
        )
    ]
    if not invalid_coverage.empty:
        raise ValueError(
            "'--reconciliation' eligible and coverage_status values are inconsistent."
        )
    species_orientation_columns = [
        "species_numerator_branch_id",
        "species_numerator_event_id",
        "species_denominator_branch_id",
        "species_denominator_event_id",
    ]
    missing_species_orientation = dataframe[
        (dataframe["eligible"] == "yes")
        & dataframe[species_orientation_columns].eq("").any(axis=1)
    ]
    unexpected_species_orientation = dataframe[
        (dataframe["eligible"] == "no")
        & dataframe[species_orientation_columns].ne("").any(axis=1)
    ]
    if (
        not missing_species_orientation.empty
        or not unexpected_species_orientation.empty
    ):
        raise ValueError(
            "'--reconciliation' species contrast orientation and eligible values are inconsistent."
        )
    nontransfer_endpoints = dataframe[
        (dataframe["event_type"] != "transfer")
        & (
            (dataframe["transfer_source_species"] != "")
            | (dataframe["transfer_destination_species"] != "")
        )
    ]
    if not nontransfer_endpoints.empty:
        raise ValueError(
            "'--reconciliation' transfer endpoints are only valid for transfer events."
        )
    transfer_endpoints = dataframe[dataframe["event_type"] == "transfer"]
    incomplete_transfer = transfer_endpoints[
        (transfer_endpoints["transfer_source_species"] == "")
        != (transfer_endpoints["transfer_destination_species"] == "")
    ]
    if not incomplete_transfer.empty:
        raise ValueError(
            "'--reconciliation' transfer endpoints must define both source and destination."
        )


def _csv_items(value):
    if value == "":
        return tuple()
    try:
        return tuple(next(csv.reader([value], strict=True)))
    except csv.Error as exc:
        raise ValueError(
            "'--reconciliation' contains malformed descendant-taxon CSV."
        ) from exc


def _validated_int(value, column, *, positive=False):
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            "'--reconciliation' {} values must be integers.".format(column)
        ) from exc
    minimum = 1 if positive else 0
    if parsed < minimum:
        qualifier = "positive" if positive else "non-negative"
        raise ValueError(
            "'--reconciliation' {} values must be {} integers.".format(
                column, qualifier
            )
        )
    return parsed


def _validate_reconciliation_coverage(dataframe):
    detail_columns = list()
    for side in COVERAGE_SIDES:
        detail_columns.extend(
            [
                "{}_observed_species_taxa".format(side),
                "num_{}_observed_species_taxa".format(side),
                "{}_species_clade_taxa".format(side),
                "num_{}_species_clade_taxa".format(side),
                "{}_species_coverage".format(side),
            ]
        )
    not_applicable = dataframe[dataframe["coverage_status"] == "not-applicable"]
    if any((not_applicable[column] != "").any() for column in detail_columns):
        raise ValueError(
            "'--reconciliation' not-applicable coverage rows must have empty coverage details."
        )
    for _, row in dataframe[
        dataframe["coverage_status"] != "not-applicable"
    ].iterrows():
        is_complete = True
        for side in COVERAGE_SIDES:
            observed_column = "{}_observed_species_taxa".format(side)
            full_column = "{}_species_clade_taxa".format(side)
            observed_count_column = "num_{}_observed_species_taxa".format(side)
            full_count_column = "num_{}_species_clade_taxa".format(side)
            fraction_column = "{}_species_coverage".format(side)
            observed = _csv_items(row[observed_column])
            full = _csv_items(row[full_column])
            observed_count = _validated_int(
                row[observed_count_column], observed_count_column, positive=True
            )
            full_count = _validated_int(
                row[full_count_column], full_count_column, positive=True
            )
            if (
                len(observed) != len(set(observed))
                or len(full) != len(set(full))
                or len(observed) != observed_count
                or len(full) != full_count
                or not set(observed).issubset(full)
            ):
                raise ValueError(
                    "'--reconciliation' {} coverage taxa and counts are inconsistent.".format(
                        side
                    )
                )
            try:
                fraction = float(row[fraction_column])
            except (TypeError, ValueError) as exc:
                raise ValueError(
                    "'--reconciliation' {} values must be finite fractions.".format(
                        fraction_column
                    )
                ) from exc
            expected_fraction = observed_count / full_count
            if (
                not math.isfinite(fraction)
                or not (0.0 < fraction <= 1.0)
                or not math.isclose(
                    fraction, expected_fraction, rel_tol=1e-12, abs_tol=0.0
                )
            ):
                raise ValueError(
                    "'--reconciliation' {} values must match observed/full counts.".format(
                        fraction_column
                    )
                )
            is_complete &= set(observed) == set(full)
        expected_status = "complete" if is_complete else "partial"
        if row["coverage_status"] != expected_status:
            raise ValueError(
                "'--reconciliation' coverage_status does not match its coverage details."
            )


def _validate_reconciliation_structure(dataframe, clades):
    branch_ids = [
        _validated_int(value, "gene_branch_id") for value in dataframe["gene_branch_id"]
    ]
    if len(branch_ids) != len(set(branch_ids)):
        raise ValueError(
            "'--reconciliation' contains duplicated gene_branch_id values."
        )
    dataframe["gene_branch_id"] = branch_ids
    lineage_ids = [
        _validated_int(value, "lineage_id") for value in dataframe["lineage_id"]
    ]
    for value in dataframe.loc[
        dataframe["mapping_status"] == "mapped", "species_branch_id"
    ]:
        _validated_int(value, "species_branch_id")
    for column in ("species_numerator_branch_id", "species_denominator_branch_id"):
        for value in dataframe.loc[dataframe["eligible"] == "yes", column]:
            _validated_int(value, column)
    if len(set(dataframe["event_source"])) != 1:
        raise ValueError(
            "'--reconciliation' must contain exactly one event_source value."
        )
    reconciliation_ids = list(dataframe["gene_clade_id"])
    if any(value == "" for value in reconciliation_ids):
        raise ValueError("'--reconciliation' gene_clade_id values cannot be empty.")
    if len(reconciliation_ids) != len(set(reconciliation_ids)):
        raise ValueError("'--reconciliation' contains duplicated gene_clade_id values.")
    node_by_id = {clades.clade_id_for_node(node): node for node in clades.mask_by_node}
    if set(reconciliation_ids) != set(node_by_id):
        raise ValueError(
            "'--reconciliation' gene_clade_id values do not exactly cover --infile."
        )
    record_by_id = dataframe.set_index("gene_clade_id", drop=False).to_dict("index")
    for row_number, (_, row) in enumerate(dataframe.iterrows()):
        clade_id = row["gene_clade_id"]
        node = node_by_id[clade_id]
        if (row["event_type"] == "leaf") != node.is_leaf:
            raise ValueError(
                "'--reconciliation' event_type=leaf does not match --infile topology."
            )
        numerator_id = row["contrast_numerator_gene_clade_id"]
        denominator_id = row["contrast_denominator_gene_clade_id"]
        orientation_present = numerator_id != "" and denominator_id != ""
        if (numerator_id == "") != (denominator_id == ""):
            raise ValueError(
                "'--reconciliation' contrast orientation must define both child clade IDs."
            )
        if orientation_present != (row["eligible"] == "yes"):
            raise ValueError(
                "'--reconciliation' contrast orientation and eligible values are inconsistent."
            )
        if orientation_present:
            child_ids = {clades.clade_id_for_node(child) for child in node.children}
            if {numerator_id, denominator_id} != child_ids:
                raise ValueError(
                    "'--reconciliation' contrast orientation must reference immediate children."
                )
        lineage_clade_id = row["lineage_clade_id"]
        if lineage_clade_id not in node_by_id:
            raise ValueError(
                "'--reconciliation' lineage_clade_id references an unknown clade."
            )
        lineage_record = record_by_id[lineage_clade_id]
        if lineage_ids[row_number] != lineage_record["gene_branch_id"]:
            raise ValueError(
                "'--reconciliation' lineage_id and lineage_clade_id are inconsistent."
            )
        node_mask = clades.mask_by_node[node]
        lineage_mask = clades.mask_by_node[node_by_id[lineage_clade_id]]
        if node_mask & ~lineage_mask:
            raise ValueError(
                "'--reconciliation' lineage_clade_id must identify an ancestor clade."
            )
    _validate_reconciliation_coverage(dataframe)


def _read_reconciliation(path, clades):
    if path in (None, ""):
        return None, ""
    dataframe = pd.read_csv(
        StringIO(read_input_text(path)),
        sep="\t",
        keep_default_na=False,
        dtype=str,
    )
    missing_columns = sorted(RECONCILIATION_REQUIRED_COLUMNS - set(dataframe.columns))
    if missing_columns:
        raise ValueError(
            "'--reconciliation' is missing required column(s): {}.".format(
                ", ".join(missing_columns)
            )
        )
    _validate_reconciliation_domains(dataframe)
    _validate_reconciliation_structure(dataframe, clades)
    tree_ids = set(dataframe["tree_id"])
    if len(tree_ids) > 1:
        raise ValueError("'--reconciliation' must contain exactly one tree_id value.")
    record_by_id = {
        str(record["gene_clade_id"]): record for record in dataframe.to_dict("records")
    }
    return record_by_id, next(iter(tree_ids), "")


def _orient_children(tree, clades, reconciliation_by_id):
    node_by_clade_id = {
        clades.clade_id_for_node(node): node for node in tree.traverse()
    }
    orientation_by_node = dict()
    for node in tree.traverse():
        if node.is_leaf:
            continue
        children = sorted(node.children, key=clades.names_for_node)
        if reconciliation_by_id is not None:
            record = reconciliation_by_id[clades.clade_id_for_node(node)]
            numerator_id = record["contrast_numerator_gene_clade_id"]
            denominator_id = record["contrast_denominator_gene_clade_id"]
            if (numerator_id == "") != (denominator_id == ""):
                raise ValueError(
                    "'--reconciliation' contrast orientation is incomplete at gene_clade_id {}.".format(
                        clades.clade_id_for_node(node)
                    )
                )
            if numerator_id != "":
                if (
                    numerator_id not in node_by_clade_id
                    or denominator_id not in node_by_clade_id
                ):
                    raise ValueError(
                        "'--reconciliation' contrast orientation references an unknown clade."
                    )
                oriented = [
                    node_by_clade_id[numerator_id],
                    node_by_clade_id[denominator_id],
                ]
                if set(oriented) != set(node.children):
                    raise ValueError(
                        "'--reconciliation' contrast orientation must reference the two immediate children at gene_clade_id {}.".format(
                            clades.clade_id_for_node(node)
                        )
                    )
                children = oriented
        orientation_by_node[node] = tuple(children)
    return orientation_by_node


def _edge_variance(node, edge_variances):
    return float(edge_variances[node])


def _validated_leaf_values(tree, values_by_leaf):
    expected = {str(leaf.name) for leaf in tree.leaves()}
    observed = {str(name) for name in values_by_leaf}
    missing = sorted(expected - observed)
    extra = sorted(observed - expected)
    if missing or extra:
        raise ValueError(
            "Trait values must exactly cover tree tips (missing={}; extra={}).".format(
                ",".join(missing), ",".join(extra)
            )
        )
    validated = dict()
    for name in sorted(expected):
        try:
            value = float(values_by_leaf[name])
        except (TypeError, ValueError) as exc:
            raise ValueError("Trait values must be numeric and finite.") from exc
        if not math.isfinite(value):
            raise ValueError("Trait values must be numeric and finite.")
        validated[name] = value
    return validated


def _require_finite(value, label):
    if not math.isfinite(value):
        raise ValueError("A contrast produced a non-finite {}.".format(label))
    return value


def _stable_ancestral_estimate(value1, variance1, value2, variance2):
    weight1, weight2 = _stable_ancestral_weights(variance1, variance2)
    return value1 * weight1 + value2 * weight2


def _stable_ancestral_weights(variance1, variance2):
    scale = max(variance1, variance2)
    scaled1 = variance1 / scale
    scaled2 = variance2 / scale
    denominator = scaled1 + scaled2
    return scaled2 / denominator, scaled1 / denominator


def _stable_adjusted_variance(variance1, variance2):
    smaller = min(variance1, variance2)
    larger = max(variance1, variance2)
    return smaller / (1.0 + smaller / larger)


def calculate_contrasts(
    tree,
    values_by_leaf,
    branch_length="original",
    evolution_model="brownian",
    evolution_parameter=None,
    orientation_by_node=None,
    return_coefficients=False,
    sparse_coefficients=False,
):
    if branch_length not in {"original", "unit"}:
        raise ValueError("Unsupported contrast branch-length mode.")
    evolution_spec = evolution_model_spec(evolution_model)
    _validate_contrast_tree(
        tree,
        branch_length if evolution_spec.branch_lengths_used else "unit",
    )
    edge_variances = transformed_edge_variances(
        tree,
        model=evolution_model,
        parameter=evolution_parameter,
        branch_length=branch_length,
    )
    values_by_leaf = _validated_leaf_values(tree, values_by_leaf)
    if orientation_by_node is None:
        clades = CladeIndex(tree)
        orientation_by_node = _orient_children(tree, clades, reconciliation_by_id=None)
    estimate_by_node = dict()
    variance_by_node = dict()
    contrast_by_node = dict()
    coefficient_by_node = dict()
    contrast_coefficient_by_node = dict()
    leaf_names = [str(leaf.name) for leaf in tree.leaves()]
    leaf_index = {name: index for index, name in enumerate(leaf_names)}
    for node in tree.traverse(strategy="postorder"):
        if node.is_leaf:
            estimate_by_node[node] = values_by_leaf[str(node.name)]
            variance_by_node[node] = _edge_variance(node, edge_variances)
            if return_coefficients:
                if sparse_coefficients:
                    coefficients = sparse.csr_matrix(
                        (
                            [1.0],
                            ([0], [leaf_index[str(node.name)]]),
                        ),
                        shape=(1, len(leaf_names)),
                    )
                else:
                    coefficients = np.zeros(len(leaf_names), dtype=float)
                    coefficients[leaf_index[str(node.name)]] = 1.0
                coefficient_by_node[node] = coefficients
            continue
        numerator, denominator = orientation_by_node[node]
        numerator_variance = variance_by_node[numerator]
        denominator_variance = variance_by_node[denominator]
        contrast_variance = numerator_variance + denominator_variance
        if contrast_variance <= 0.0 or not math.isfinite(contrast_variance):
            raise ValueError("A contrast has non-positive or non-finite variance.")
        numerator_estimate = estimate_by_node[numerator]
        denominator_estimate = estimate_by_node[denominator]
        raw_contrast = _require_finite(
            numerator_estimate - denominator_estimate, "raw contrast"
        )
        ancestral_estimate = _require_finite(
            _stable_ancestral_estimate(
                numerator_estimate,
                numerator_variance,
                denominator_estimate,
                denominator_variance,
            ),
            "ancestral estimate",
        )
        weight1, weight2 = _stable_ancestral_weights(
            numerator_variance, denominator_variance
        )
        if return_coefficients:
            contrast_coefficients = (
                coefficient_by_node[numerator] - coefficient_by_node[denominator]
            )
            ancestral_coefficients = (
                weight1 * coefficient_by_node[numerator]
                + weight2 * coefficient_by_node[denominator]
            )
        adjusted_variance = _stable_adjusted_variance(
            numerator_variance, denominator_variance
        )
        node_variance = _require_finite(
            adjusted_variance + _edge_variance(node, edge_variances),
            "adjusted variance",
        )
        standardized_contrast = _require_finite(
            raw_contrast / math.sqrt(contrast_variance),
            "standardized contrast",
        )
        estimate_by_node[node] = ancestral_estimate
        variance_by_node[node] = node_variance
        if return_coefficients:
            coefficient_by_node[node] = ancestral_coefficients
            contrast_coefficient_by_node[node] = contrast_coefficients
        contrast_by_node[node] = {
            "raw_contrast": raw_contrast,
            "contrast_variance": contrast_variance,
            "standardized_contrast": standardized_contrast,
            "ancestral_estimate": ancestral_estimate,
        }
    if return_coefficients:
        return contrast_by_node, contrast_coefficient_by_node, leaf_names
    return contrast_by_node


def _structured_tip_sampling_factor(coefficients, leaf_names, covariance):
    tip_diagonal = np.asarray(covariance.diagonal, dtype=float)
    tip_low_rank = covariance.low_rank
    finite_loading = (
        np.isfinite(tip_low_rank.data).all()
        if sparse.issparse(tip_low_rank)
        else np.isfinite(np.asarray(tip_low_rank, dtype=float)).all()
    )
    if (
        tip_diagonal.shape != (len(leaf_names),)
        or np.any(tip_diagonal < 0.0)
        or not np.isfinite(tip_diagonal).all()
        or len(tip_low_rank.shape) != 2
        or tip_low_rank.shape[0] != len(leaf_names)
        or not finite_loading
    ):
        raise ValueError("Tip sampling covariance factor is incomplete.")
    coefficient_matrix = sparse.csr_matrix(coefficients)
    factor = sparse.hstack(
        [
            coefficient_matrix.multiply(np.sqrt(tip_diagonal)[None, :]),
            coefficient_matrix @ sparse.csr_matrix(tip_low_rank),
        ],
        format="csr",
    )
    latent_ids = ["latent:tip:{}".format(name) for name in leaf_names] + [
        "latent:low-rank:{}".format(index) for index in range(tip_low_rank.shape[1])
    ]
    return factor, latent_ids


def _tip_sampling_array(value, leaf_names, trait):
    if isinstance(value, pd.DataFrame):
        covariance = value.reindex(index=leaf_names, columns=leaf_names)
        if covariance.isna().any(axis=None):
            raise ValueError(
                "Tip sampling covariance is incomplete for trait '{}'.".format(trait)
            )
        return covariance.to_numpy(dtype=float)
    covariance = np.asarray(value, dtype=float)
    if covariance.shape != (len(leaf_names),):
        raise ValueError("Tip sampling covariance diagonal is incomplete.")
    if not np.isfinite(covariance).all() or np.any(covariance < 0.0):
        raise ValueError(
            "Tip sampling covariance diagonal for trait '{}' must contain "
            "finite non-negative values.".format(trait)
        )
    return covariance


def _dense_contrast_sampling_covariance(coefficients, tip_covariance):
    dense_coefficients = (
        coefficients.toarray() if sparse.issparse(coefficients) else coefficients
    )
    covariance = dense_coefficients @ tip_covariance @ dense_coefficients.T
    covariance = (covariance + covariance.T) / 2.0
    if not np.isfinite(covariance).all():
        raise ValueError("Contrast sampling covariance contains non-finite values.")
    eigenvalues = np.linalg.eigvalsh(covariance)
    tolerance = (
        np.finfo(float).eps
        * max(1.0, float(np.max(np.abs(covariance))))
        * max(1, len(covariance))
    )
    if eigenvalues.size and float(eigenvalues.min()) < -tolerance:
        raise ValueError("Contrast sampling covariance is not positive semidefinite.")
    covariance[np.abs(covariance) < tolerance] = 0.0
    return covariance


def _contrast_sampling_representation(coefficients, leaf_names, value, trait):
    if isinstance(value, DiagonalLowRankCovariance):
        factor, latent_ids = _structured_tip_sampling_factor(
            coefficients, leaf_names, value
        )
    else:
        tip_covariance = _tip_sampling_array(value, leaf_names, trait)
        if tip_covariance.ndim != 1:
            covariance = _dense_contrast_sampling_covariance(
                coefficients, tip_covariance
            )
            return False, None, None, np.diag(covariance), covariance
        factor = sparse.csr_matrix(coefficients).multiply(
            np.sqrt(tip_covariance)[None, :]
        )
        latent_ids = ["latent:{}".format(name) for name in leaf_names]
    factor = _prune_sparse_factor_loadings(factor)
    diagonal = np.asarray(factor.multiply(factor).sum(axis=1)).reshape(-1)
    return True, factor, latent_ids, diagonal, None


def _sampling_covariance_table(
    table,
    coefficient_by_trait_and_clade,
    leaf_names_by_trait,
    sampling_covariance_by_trait,
    replicate_model_by_trait,
    tip_summary,
    compact_sampling_covariance,
):
    covariance_rows = []
    table = table.copy()
    for column in REPLICATE_CONTRAST_COLUMNS:
        table[column] = pd.Series([""] * len(table), index=table.index, dtype=object)
    n_by_trait_and_leaf = {}
    if tip_summary is not None and not tip_summary.empty:
        n_by_trait_and_leaf = {
            (str(row["trait"]), str(row["leaf_name"])): int(row["n_biological"])
            for row in tip_summary.to_dict("records")
            if row["n_biological"] not in (None, "")
            and not pd.isna(row["n_biological"])
        }
    for trait, trait_rows in table.groupby("trait", sort=False):
        if trait not in sampling_covariance_by_trait:
            raise ValueError(
                "Sampling covariance is absent for trait '{}'.".format(trait)
            )
        row_indices = list(trait_rows.index)
        contrast_ids = [
            str(table.loc[index, "branch_clade_id"]) for index in row_indices
        ]
        leaf_names = leaf_names_by_trait[trait]
        coefficient_values = [
            coefficient_by_trait_and_clade[(trait, contrast_id)]
            for contrast_id in contrast_ids
        ]
        coefficients = (
            sparse.vstack(coefficient_values, format="csr")
            if any(sparse.issparse(value) for value in coefficient_values)
            else np.vstack(coefficient_values)
        )
        (
            diagonal_sampling,
            contrast_factor,
            factor_latent_ids,
            sampling_diagonal,
            contrast_covariance,
        ) = _contrast_sampling_representation(
            coefficients,
            leaf_names,
            sampling_covariance_by_trait[trait],
            trait,
        )
        for position, row_index in enumerate(row_indices):
            coefficient = coefficients[position]
            if sparse.issparse(coefficient):
                coefficient = coefficient.tocsr()
                supported_indices = coefficient.indices[
                    np.abs(coefficient.data) > np.finfo(float).eps
                ]
            else:
                supported_indices = np.flatnonzero(
                    np.abs(coefficient) > np.finfo(float).eps
                )
            supported_leaves = [leaf_names[index] for index in supported_indices]
            sample_sizes = [
                n_by_trait_and_leaf[(trait, leaf)]
                for leaf in supported_leaves
                if (trait, leaf) in n_by_trait_and_leaf
            ]
            table.loc[row_index, "sampling_variance"] = float(
                max(sampling_diagonal[position], 0.0)
            )
            table.loc[row_index, "evolutionary_variance"] = float(
                table.loc[row_index, "contrast_variance"]
            )
            table.loc[row_index, "replicate_model"] = replicate_model_by_trait[trait]
            table.loc[row_index, "min_n_biological"] = (
                min(sample_sizes) if sample_sizes else ""
            )
        if (
            compact_sampling_covariance
            and diagonal_sampling
            and len(contrast_ids) > MAX_FULL_SAMPLING_COVARIANCE_CONTRASTS
        ):
            assert factor_latent_ids is not None
            coo = contrast_factor.tocoo()
            represented_rows = set(coo.row)
            for row, column, value in zip(coo.row, coo.col, coo.data, strict=True):
                covariance_rows.append(
                    {
                        "tree_id": str(table.loc[row_indices[row], "tree_id"]),
                        "trait": trait,
                        "contrast_id_1": contrast_ids[row],
                        "contrast_id_2": factor_latent_ids[column],
                        "sampling_covariance": float(value),
                        "covariance_representation": "factor-loading",
                    }
                )
            for row in sorted(set(range(len(contrast_ids))) - represented_rows):
                covariance_rows.append(
                    {
                        "tree_id": str(table.loc[row_indices[row], "tree_id"]),
                        "trait": trait,
                        "contrast_id_1": contrast_ids[row],
                        "contrast_id_2": "latent:zero",
                        "sampling_covariance": 0.0,
                        "covariance_representation": "factor-loading",
                    }
                )
        else:
            if contrast_covariance is None:
                contrast_covariance = np.asarray(
                    (contrast_factor @ contrast_factor.T).toarray()
                )
            for row1, contrast_id_1 in enumerate(contrast_ids):
                for row2 in range(row1, len(contrast_ids)):
                    covariance_rows.append(
                        {
                            "tree_id": str(table.loc[row_indices[row1], "tree_id"]),
                            "trait": trait,
                            "contrast_id_1": contrast_id_1,
                            "contrast_id_2": contrast_ids[row2],
                            "sampling_covariance": float(
                                contrast_covariance[row1, row2]
                            ),
                            "covariance_representation": "covariance",
                        }
                    )
    return table, pd.DataFrame(covariance_rows, columns=SAMPLING_COVARIANCE_COLUMNS)


def _prune_sparse_factor_loadings(factor):
    """Drop roundoff-scale loadings while retaining row covariance numerically."""
    factor = sparse.csr_matrix(factor, dtype=float).copy()
    relative_tolerance = np.finfo(float).eps * 16.0
    for row in range(factor.shape[0]):
        start, stop = factor.indptr[row : row + 2]
        if start == stop:
            continue
        values = factor.data[start:stop]
        tolerance = relative_tolerance * float(np.max(np.abs(values)))
        values[np.abs(values) <= tolerance] = 0.0
    factor.eliminate_zeros()
    return factor


def build_contrast_table(
    tree,
    values_by_trait,
    branch_length="original",
    evolution_model="brownian",
    evolution_parameter=None,
    reconciliation_by_id=None,
    event_type="all",
    eligible_only=True,
    speciation_coverage="complete",
    tree_id="",
    sampling_covariance_by_trait=None,
    replicate_model_by_trait=None,
    tip_summary=None,
    return_sampling_covariance=False,
    compact_sampling_covariance=True,
    tree_option_name="--infile",
):
    if branch_length not in {"original", "unit"}:
        raise ValueError("Unsupported contrast branch-length mode.")
    evolution_spec = evolution_model_spec(evolution_model)
    if not evolution_spec.contrast_supported:
        raise ValueError(
            "Evolutionary model '{}' cannot be represented by local contrasts.".format(
                evolution_model
            )
        )
    evolution_parameter = validate_evolution_parameter(
        evolution_model, evolution_parameter
    )
    _validate_contrast_tree(
        tree,
        branch_length if evolution_spec.branch_lengths_used else "unit",
        option_name=tree_option_name,
    )
    if event_type not in {"all", "speciation", "duplication", "transfer", "unresolved"}:
        raise ValueError("Unsupported event type: {}".format(event_type))
    if speciation_coverage not in {"complete", "any"}:
        raise ValueError(
            "Unsupported speciation coverage policy: {}".format(speciation_coverage)
        )
    branch_ids = assign_branch_ids(tree)
    clades = CladeIndex(tree)
    orientation_by_node = _orient_children(tree, clades, reconciliation_by_id)
    rows = list()
    coefficient_by_trait_and_clade = {}
    leaf_names_by_trait = {}
    for trait, values_by_leaf in values_by_trait.items():
        contrast_result = calculate_contrasts(
            tree,
            values_by_leaf,
            branch_length=branch_length,
            evolution_model=evolution_model,
            evolution_parameter=evolution_parameter,
            orientation_by_node=orientation_by_node,
            return_coefficients=sampling_covariance_by_trait is not None,
            sparse_coefficients=(
                sampling_covariance_by_trait is not None
                and len(list(tree.leaves())) > MAX_FULL_SAMPLING_COVARIANCE_CONTRASTS
            ),
        )
        if sampling_covariance_by_trait is None:
            contrasts = contrast_result
            contrast_coefficients = None
        else:
            contrasts, contrast_coefficients, leaf_names = contrast_result
            leaf_names_by_trait[trait] = leaf_names
        for node in tree.traverse():
            if node.is_leaf:
                continue
            branch_id = branch_ids[node]
            branch_clade_id = clades.clade_id_for_node(node)
            reconciliation = (
                None
                if reconciliation_by_id is None
                else reconciliation_by_id[branch_clade_id]
            )
            if reconciliation is not None:
                if event_type != "all" and reconciliation["event_type"] != event_type:
                    continue
                if eligible_only and reconciliation["eligible"] != "yes":
                    continue
                if (
                    eligible_only
                    and speciation_coverage == "complete"
                    and reconciliation["coverage_status"] != "complete"
                ):
                    continue
            numerator, denominator = orientation_by_node[node]
            row = {
                "tree_id": tree_id,
                "branch_id": branch_id,
                "branch_clade_id": branch_clade_id,
                "node_class": get_node_class(node),
                "descendant_taxa": clades.csv_for_node(node),
                "num_taxa": clades.count_for_node(node),
                "trait": trait,
                "evolution_model": evolution_model,
                "evolution_parameter_name": evolution_spec.parameter_name or "",
                "evolution_parameter": (
                    evolution_parameter if evolution_parameter is not None else ""
                ),
                "branch_length_mode": (
                    branch_length
                    if evolution_spec.branch_lengths_used
                    else "not-applicable"
                ),
                "numerator_branch_id": branch_ids[numerator],
                "numerator_clade_id": clades.clade_id_for_node(numerator),
                "denominator_branch_id": branch_ids[denominator],
                "denominator_clade_id": clades.clade_id_for_node(denominator),
                "numerator_taxa": clades.csv_for_node(numerator),
                "denominator_taxa": clades.csv_for_node(denominator),
            }
            row.update(contrasts[node])
            if contrast_coefficients is not None:
                coefficient_by_trait_and_clade[(trait, branch_clade_id)] = (
                    contrast_coefficients[node]
                )
            if reconciliation is not None:
                for column in RECONCILIATION_CONTEXT_COLUMNS:
                    row[column] = reconciliation[column]
            rows.append(row)
    columns = list(BASE_CONTRAST_COLUMNS)
    if reconciliation_by_id is not None:
        columns = (
            BASE_CONTRAST_COLUMNS[:4]
            + RECONCILIATION_CONTEXT_COLUMNS
            + BASE_CONTRAST_COLUMNS[4:]
        )
    table = pd.DataFrame(rows, columns=columns)
    covariance_table = pd.DataFrame(columns=SAMPLING_COVARIANCE_COLUMNS)
    if sampling_covariance_by_trait is not None:
        if replicate_model_by_trait is None:
            raise ValueError(
                "Replicate model metadata is required with sampling covariance."
            )
        table, covariance_table = _sampling_covariance_table(
            table,
            coefficient_by_trait_and_clade,
            leaf_names_by_trait,
            sampling_covariance_by_trait,
            replicate_model_by_trait,
            tip_summary,
            compact_sampling_covariance,
        )
    if return_sampling_covariance:
        return table, covariance_table
    return table


def contrast_main(args):
    outputs = [
        ("--outfile", args.outfile),
        ("--sampling-covariance-out", getattr(args, "sampling_covariance_out", None)),
        ("--tip-summary-out", getattr(args, "tip_summary_out", None)),
    ]
    validate_distinct_output_paths(outputs)
    validate_outputs_do_not_replace_inputs(
        [
            ("--infile", args.infile),
            ("--trait", args.trait),
            ("--reconciliation", getattr(args, "reconciliation", None)),
        ],
        outputs,
        label="Contrast output",
    )
    tree = read_tree(args.infile, args.format, args.quoted_node_names)
    columns = _parse_columns(args.columns)
    _validate_contrast_tree(tree, args.branch_length)
    clades = CladeIndex(tree)
    reconciliation_by_id, reconciliation_tree_id = _read_reconciliation(
        args.reconciliation, clades
    )
    if reconciliation_by_id is None and args.event_type != "all":
        raise ValueError("'--event-type' requires '--reconciliation'.")
    requested_tree_id = getattr(args, "tree_id", "")
    if (
        requested_tree_id != ""
        and reconciliation_tree_id != ""
        and requested_tree_id != reconciliation_tree_id
    ):
        raise ValueError("'--tree-id' does not match --reconciliation tree_id.")
    tree_id = requested_tree_id or reconciliation_tree_id
    replicate_requested = _validate_replicate_options(args)
    replicate_estimates = None
    if replicate_requested:
        covariance_out = getattr(args, "sampling_covariance_out", None)
        if covariance_out in (None, ""):
            raise ValueError(
                "Replicate-aware contrasts require '--sampling-covariance-out'."
            )
        if covariance_out == "-":
            raise ValueError("'--sampling-covariance-out' cannot be STDOUT.")
        replicate_estimates = _read_replicate_traits(args, tree, columns, tree_id)
        values_by_trait = replicate_estimates.values_by_trait
    else:
        values_by_trait = _read_numeric_traits(args, tree, columns)
    eligible_only = getattr(args, "eligible_only", None)
    if eligible_only is None:
        eligible_only = args.event_type == "speciation"
    output = build_contrast_table(
        tree,
        values_by_trait,
        branch_length=args.branch_length,
        evolution_model=args.evolution_model,
        evolution_parameter=args.evolution_parameter,
        reconciliation_by_id=reconciliation_by_id,
        event_type=args.event_type,
        eligible_only=eligible_only,
        speciation_coverage=getattr(args, "speciation_coverage", "complete"),
        tree_id=tree_id,
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
        return_sampling_covariance=replicate_estimates is not None,
    )
    if replicate_estimates is None:
        table = output
        file_outputs = []
    else:
        table, covariance_table = output
        file_outputs = [(args.sampling_covariance_out, covariance_table)]
        tip_summary_out = getattr(args, "tip_summary_out", None)
        if tip_summary_out is not None:
            if tip_summary_out == "-":
                raise ValueError("'--tip-summary-out' cannot be STDOUT.")
            file_outputs.append((tip_summary_out, replicate_estimates.tip_summary))
    if args.outfile != "-":
        file_outputs.append((args.outfile, table))
    if file_outputs:
        from nwkit.pgls_pipeline import _write_dataframes_transactionally

        _write_dataframes_transactionally(file_outputs)
    if args.outfile == "-":
        print(table.to_csv(sep="\t", index=False), end="")
