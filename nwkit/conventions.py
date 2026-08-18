from typing import Any

DEFAULT_TABLE_MISSING_VALUES = ("", "NA", "NaN", "nan", "?", "missing", "unknown")
DEFAULT_TABLE_MISSING_VALUES_CSV = ",".join(DEFAULT_TABLE_MISSING_VALUES)

REGRESSION_BUNDLE_SUFFIXES = {
    "reconciliation_out": ".reconciliation.tsv",
    "gene_contrasts_out": ".gene-contrasts.tsv",
    "species_contrasts_out": ".species-contrasts.tsv",
    "response_sampling_covariance_out": ".response-sampling-covariance.tsv",
    "response_tip_summary_out": ".response-tip-summary.tsv",
    "predictor_sampling_covariance_out": ".predictor-sampling-covariance.tsv",
    "predictor_tip_summary_out": ".predictor-tip-summary.tsv",
    "random_effects_out": ".random-effects.tsv",
    "sensitivity_out": ".sensitivity.tsv",
    "trait_origins_out": ".trait-origins.tsv",
    "outfile": ".regression.tsv",
}
REGRESSION_BUNDLE_LOCK_SUFFIX = ".regression-bundle.lock"


def regression_bundle_paths(prefix: str) -> dict[str, str]:
    """Return deterministic output paths for an end-to-end regression bundle."""
    return {
        argument: "{}{}".format(prefix, suffix)
        for argument, suffix in REGRESSION_BUNDLE_SUFFIXES.items()
    }


def regression_bundle_lock_path(prefix: str) -> str:
    """Return the internal lock path protecting a regression bundle transaction."""
    return "{}{}".format(prefix, REGRESSION_BUNDLE_LOCK_SUFFIX)


STDIN_INPUT_DESTS = (
    "infile",
    "infile2",
    "data",
    "evolution_covariance",
    "tree",
    "seqin",
    "trait",
    "table",
    "taxid_tsv",
    "weight_tsv",
    "species_map_tsv",
    "species_name_tsv",
    "name_tsv",
    "species_list",
    "reference",
    "reconciliation",
    "predictor_contrasts",
    "expression",
    "gene_tree",
    "gene_tree_ensemble",
    "reconciliation_tree",
    "response_sampling_covariance",
    "predictor_sampling_covariance",
    "species_traits",
    "species_tree",
    "root_source",
    "name_source",
    "support_source",
    "length_source",
    "tip_image_manifest",
    "mcmctree_posterior",
    "densitree_trees",
)


def get_stdin_input_options(args: Any) -> list[tuple[str, str]]:
    options = [
        (dest, "--{}".format(dest.replace("_", "-")))
        for dest in STDIN_INPUT_DESTS
        if getattr(args, dest, None) == "-"
    ]
    options.extend(
        ("property_source", "--property-source")
        for value in (getattr(args, "property_source", None) or [])
        if str(value).rsplit("@", 1)[-1] == "-"
    )
    return options
