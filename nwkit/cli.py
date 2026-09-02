import argparse
import math
import sys

from nwkit import __version__
from nwkit.asr_models import model_names
from nwkit.conventions import DEFAULT_TABLE_MISSING_VALUES_CSV, get_stdin_input_options
from nwkit.model_specs import (
    CONTRAST_EVOLUTION_MODELS,
    EVOLUTION_MODELS,
    RESPONSE_FAMILIES,
)
from nwkit.species_parser import (
    DEFAULT_SPECIES_PARSER,
    DEFAULT_SPECIES_REGEX,
    SUPPORTED_SPECIES_PARSERS,
)


def _canonical_long_option(action):
    return next(
        (
            option
            for option in action.option_strings
            if option.startswith("--") and "_" not in option
        ),
        None,
    )


class CanonicalHelpFormatter(argparse.HelpFormatter):
    """Display canonical options while keeping compatibility aliases parseable."""

    def _format_action_invocation(self, action):
        if not action.option_strings:
            return super()._format_action_invocation(action)
        short_options = [
            option for option in action.option_strings if not option.startswith("--")
        ]
        canonical_long_option = _canonical_long_option(action)
        option_strings = short_options
        if canonical_long_option is not None:
            option_strings.append(canonical_long_option)
        if not option_strings:
            option_strings = list(action.option_strings)
        if action.nargs == 0:
            return ", ".join(option_strings)
        default = self._get_default_metavar_for_optional(action)
        args_string = self._format_args(action, default)
        return ", ".join(
            "{} {}".format(option, args_string) for option in option_strings
        )


class NwkitArgumentParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault("formatter_class", CanonicalHelpFormatter)
        super().__init__(*args, **kwargs)


def strtobool(val):
    val = val.lower()
    if val in ("y", "yes", "t", "true", "on", "1"):
        return True
    elif val in ("n", "no", "f", "false", "off", "0"):
        return False
    else:
        raise ValueError(f"Invalid truth value: {val}")


def strtoautobool(val):
    if val.lower() == "auto":
        return None
    return strtobool(val)


def finite_float(value):
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise argparse.ArgumentTypeError(
            "'{}' is not a valid floating-point value.".format(value)
        ) from exc
    if not math.isfinite(number):
        raise argparse.ArgumentTypeError("Floating-point values must be finite.")
    return number


def nonnegative_finite_float(value):
    number = finite_float(value)
    if number < 0:
        raise argparse.ArgumentTypeError("Floating-point values must be non-negative.")
    return number


def auto_or_finite_float(value):
    if isinstance(value, str) and value.lower() == "auto":
        return "auto"
    return finite_float(value)


def unit_interval_float(value):
    number = finite_float(value)
    if not 0.0 <= number <= 1.0:
        raise argparse.ArgumentTypeError(
            "Value must be between 0.0 and 1.0, inclusive."
        )
    return number


# Main parser
parser = NwkitArgumentParser(
    prog="nwkit",
    description=(
        "A toolkit for processing and visualizing phylogenetic trees and for "
        "phylogeny-aware comparative analysis. See `nwkit SUBCOMMAND -h` for "
        "usage (e.g., nwkit constrain -h)"
    ),
)
parser.add_argument(
    "--version", action="version", version="%(prog)s {}".format(__version__)
)
parser.add_argument(
    "--debug", action="store_true", help="Show a Python traceback when a command fails."
)
subparsers = parser.add_subparsers(dest="command")

# Parent parser for shared options
p_audit = NwkitArgumentParser(add_help=False)
p_audit.add_argument(
    "--audit",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="Append a JSONL provenance record containing arguments, hashes, warnings, and runtime metadata.",
)
p_audit.add_argument(
    "--debug",
    action="store_true",
    default=argparse.SUPPRESS,
    help="Show a Python traceback when a command fails.",
)

p_parent = NwkitArgumentParser(add_help=False, parents=[p_audit])
p_parent.add_argument(
    "-i",
    "--infile",
    metavar="PATH",
    default="-",
    type=str,
    required=False,
    action="store",
    help='default=%(default)s: Input Newick/NHX or supported PAML/MCMCtree tree container. Use "-" for STDIN.',
)
p_parent.add_argument(
    "-o",
    "--outfile",
    metavar="PATH",
    default="-",
    type=str,
    required=False,
    action="store",
    help='default=%(default)s: Output newick file. Use "-" for STDOUT.',
)
p_parent.add_argument(
    "-f",
    "--format",
    metavar="auto|auto-strict|INT",
    default="auto",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: ETE tree format. "
    '"auto-strict" fails on ambiguous unquoted numeric internal labels. '
    "See https://etetoolkit.github.io/ete/tutorial/tutorial_trees.html",
)
p_parent.add_argument(
    "-of",
    "--outformat",
    metavar="auto|INT",
    default="auto",
    type=str,
    required=False,
    action="store",
    help='ETE tree format for --outfile. "auto" indicates the same format as --format.',
)
p_parent.add_argument(
    "--quoted-node-names",
    "--quoted_node_names",
    dest="quoted_node_names",
    metavar="yes|no",
    default="yes",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Whether node names are quoted in the input file.",
)
p_preserve_properties = NwkitArgumentParser(add_help=False)
p_preserve_properties.add_argument(
    "--preserve-properties",
    "--preserve_properties",
    dest="preserve_properties",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Preserve custom node properties such as NHX D/H/S annotations in output.",
)

p_download = NwkitArgumentParser(add_help=False)
p_download.add_argument(
    "--download-dir",
    "--download_dir",
    dest="download_dir",
    metavar="PATH",
    default="auto",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Shared download/cache directory for external resources such as the ETE4 "
    'NCBI taxonomy database. "auto" uses the ETE4 default cache location. '
    '"inferred" stores downloads under <outfile_dir>/downloads, or ./downloads when writing to STDOUT.',
)
p_download.add_argument(
    "--taxonomy-cache-max-age-days",
    "--taxonomy_cache_max_age_days",
    dest="taxonomy_cache_max_age_days",
    metavar="FLOAT",
    default=30.0,
    type=finite_float,
    required=False,
    action="store",
    help="default=%(default)s: Recheck the NCBI taxonomy archive after this many days. Use 0 to check every run.",
)
p_download.add_argument(
    "--refresh-taxonomy-cache",
    "--refresh_taxonomy_cache",
    dest="refresh_taxonomy_cache",
    action="store_true",
    help="Force a checksum check and rebuild the ETE4 taxonomy cache when NCBI has a newer archive.",
)

p_tree_input = NwkitArgumentParser(add_help=False, parents=[p_audit])
p_tree_input.add_argument(
    "-i",
    "--infile",
    metavar="PATH",
    default="-",
    type=str,
    required=False,
    action="store",
    help='default=%(default)s: Input Newick/NHX or supported PAML/MCMCtree tree container. Use "-" for STDIN.',
)
p_tree_input.add_argument(
    "-f",
    "--format",
    metavar="auto|auto-strict|INT",
    default="auto",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: ETE tree format. "
    '"auto-strict" fails on ambiguous unquoted numeric internal labels. '
    "See https://etetoolkit.github.io/ete/tutorial/tutorial_trees.html",
)
p_tree_input.add_argument(
    "--quoted-node-names",
    "--quoted_node_names",
    dest="quoted_node_names",
    metavar="yes|no",
    default="yes",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Whether node names are quoted in the input file.",
)

p_table_output = NwkitArgumentParser(add_help=False)
p_table_output.add_argument(
    "-o",
    "--outfile",
    metavar="PATH",
    default="-",
    type=str,
    required=False,
    action="store",
    help='default=%(default)s: Output table file. Use "-" for STDOUT.',
)

p_text_output = NwkitArgumentParser(add_help=False)
p_text_output.add_argument(
    "-o",
    "--outfile",
    metavar="PATH",
    default="-",
    type=str,
    required=False,
    action="store",
    help='default=%(default)s: Output text file. Use "-" for STDOUT.',
)

p_table_input = NwkitArgumentParser(add_help=False, parents=[p_audit])
p_table_input.add_argument(
    "-i",
    "--infile",
    metavar="PATH",
    default="-",
    type=str,
    required=False,
    action="store",
    help='default=%(default)s: Input table file. Use "-" for STDIN.',
)
p_table_input.add_argument(
    "-o",
    "--outfile",
    metavar="PATH",
    default="-",
    type=str,
    required=False,
    action="store",
    help='default=%(default)s: Output newick file. Use "-" for STDOUT.',
)
p_table_input.add_argument(
    "-of",
    "--outformat",
    metavar="auto|INT",
    default="auto",
    type=str,
    required=False,
    action="store",
    help='ETE tree format for --outfile. "auto" infers a suitable format from the input table.',
)

p_species = NwkitArgumentParser(add_help=False)
p_species.add_argument(
    "--species-parser",
    "--species_parser",
    dest="species_parser",
    metavar="PRESET",
    default=DEFAULT_SPECIES_PARSER,
    required=False,
    type=str,
    choices=list(SUPPORTED_SPECIES_PARSERS),
    help="default=%(default)s: Species parser preset. "
    "legacy keeps the historical GENUS_SPECIES[_...] behavior. "
    "taxonomic keeps qualifier-aware species labels such as GENUS_SPECIES_cf or GENUS_sp_LABEL "
    "for internal identity while simplifying taxonomy queries.",
)
p_species.add_argument(
    "--species-regex",
    "--species_regex",
    dest="species_regex",
    metavar="REGEX",
    default=DEFAULT_SPECIES_REGEX,
    required=False,
    type=str,
    help="default=%(default)s: Regular expression for extracting species IDs when --species-parser legacy is used. "
    "If the regex contains capture groups, non-empty captured groups are joined with underscores; "
    "otherwise, the full regex match is used. "
    "For GENUS_SPECIES and GENUS_SPECIES_GENEID labels, the default works as-is.",
)
p_species.add_argument(
    "--species-map-tsv",
    "--species_map_tsv",
    dest="species_map_tsv",
    metavar="PATH",
    default=None,
    required=False,
    type=str,
    help="default=%(default)s: Optional TSV overriding parsed species labels and/or taxonomy queries. "
    'The file must contain a "leaf_name" column and at least one of "species_label" or "taxonomy_query". '
    "Mapped rows override the selected parser preset and regex.",
)


def _tip_table_policy_parent(*, unmatched_default="warn"):
    parent = NwkitArgumentParser(add_help=False)
    parent.add_argument(
        "--missing-values",
        "--missing_values",
        dest="missing_values",
        metavar="CSV",
        default=DEFAULT_TABLE_MISSING_VALUES_CSV,
        type=str,
        required=False,
        action="store",
        help="default=%(default)s: Comma-separated table values treated as missing.",
    )
    parent.add_argument(
        "--unmatched",
        metavar="warn|error|ignore",
        default=unmatched_default,
        type=str,
        required=False,
        action="store",
        choices=["warn", "error", "ignore"],
        help="default=%(default)s: Policy for table rows or tree tips without a counterpart.",
    )
    return parent


p_tip_table_policy = _tip_table_policy_parent()
p_contrast_tip_table_policy = _tip_table_policy_parent(unmatched_default="error")


def command_annotate(args):
    from nwkit.annotate import annotate_main

    annotate_main(args)


pannotate = subparsers.add_parser(
    "annotate",
    help="Attach tip-table values and aggregate them as Newick properties",
    parents=[p_parent, p_tip_table_policy],
)
pannotate.add_argument(
    "--table",
    metavar="PATH",
    default=None,
    type=str,
    required=True,
    action="store",
    help="TSV containing a leaf_name column and one or more annotation columns.",
)
pannotate.add_argument(
    "--columns",
    metavar="COL1,COL2,...",
    default=None,
    type=str,
    required=False,
    action="store",
    help="Columns attached to matching tips. By default all non-key columns are attached.",
)
pannotate.add_argument(
    "--property-map",
    "--property_map",
    dest="property_map",
    metavar="COLUMN=PROPERTY",
    default=[],
    type=str,
    required=False,
    action="append",
    help="Attach a table column under a different Newick property name. May be repeated.",
)
pannotate.add_argument(
    "--aggregate",
    metavar="COLUMN:METHOD[:PROPERTY]",
    default=[],
    type=str,
    required=False,
    action="append",
    help="Aggregate descendant-tip values onto internal nodes. Methods: unique, mode, count, mean, sum, min, max, list. May be repeated.",
)
pannotate.add_argument(
    "--report",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="Optional TSV audit of attached, aggregated, missing, and unmatched values.",
)
pannotate.set_defaults(handler=command_annotate)


def command_asr(args):
    from nwkit.asr import asr_main

    asr_main(args)


pasr = subparsers.add_parser(
    "asr",
    help="Infer ancestral traits and impute missing tips under discrete Mk or continuous Gaussian models",
    parents=[p_tree_input, p_table_output, p_tip_table_policy],
)
pasr.add_argument(
    "--trait",
    metavar="PATH",
    default=None,
    type=str,
    required=True,
    action="store",
    help='TSV containing "leaf_name" and trait columns.',
)
pasr.add_argument(
    "--state-column",
    "--state_column",
    dest="state_column",
    metavar="STR",
    default=None,
    type=str,
    required=True,
    action="store",
    help="Trait column in --trait; MV-BM/MV-OU require two or more comma-separated columns.",
)
pasr.add_argument(
    "--trait-type",
    "--trait_type",
    dest="trait_type",
    metavar="auto|discrete|continuous",
    default="auto",
    choices=["auto", "discrete", "continuous"],
    help="default=%(default)s: Infer the type from non-missing values on tree tips. "
    "Numeric values select continuous; other values select discrete. "
    "Use discrete explicitly for numeric category codes.",
)
pasr.add_argument(
    "--states",
    metavar="STATE1,STATE2,...",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Optional comma-separated state order. Unlisted observed states are rejected; "
    "required with --transition-graph ordered.",
)
pasr.add_argument(
    "--model",
    metavar="MODEL",
    default=None,
    type=str,
    required=False,
    action="store",
    choices=model_names(),
    help="default=ER for discrete, BM for continuous: Discrete ER/SYM/ARD/F81/GTR/"
    "MK-REGIME/HRM/COVARION/MK-MIXTURE/THRESHOLD/CUSTOM or continuous BM/BMS/LAMBDA/KAPPA/DELTA/EB/ACDC/"
    "BMS-DRIFT/BM-DRIFT/MV-BM/MV-OU/OU/OUM/OUMA/OUMV/OUMVA.",
)
pasr.add_argument(
    "--rate",
    metavar="FLOAT",
    default=None,
    type=finite_float,
    required=False,
    action="store",
    help="default=%(default)s: Fixed ER off-diagonal rate; optimizer starting rate for other fitted discrete models.",
)
pasr.add_argument(
    "--rate-bounds",
    "--rate_bounds",
    dest="rate_bounds",
    metavar="MIN,MAX",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=1e-9,1e3: Positive bounds for fitted Mk rates and GTR exchangeabilities.",
)
pasr.add_argument(
    "--transition-graph",
    "--transition_graph",
    dest="transition_graph",
    metavar="complete|ordered|PATH",
    default=None,
    type=str,
    help="ER/SYM/ARD, MK-REGIME, HRM, or COVARION: Allowed observed transitions. 'ordered' connects "
    "adjacent explicit --states bidirectionally; PATH is a TSV from_state/to_state edge list. "
    "F81/GTR require complete.",
)
pasr.add_argument(
    "--regime-map",
    "--regime_map",
    dest="regime_map",
    metavar="PATH",
    default=None,
    type=str,
    help="Regime models only: TSV assigning every branch_id, including root 0, to a regime.",
)
pasr.add_argument(
    "--regime-model",
    "--regime_model",
    dest="regime_model",
    metavar="ER|SYM|ARD|F81|GTR",
    default=None,
    choices=["ER", "SYM", "ARD", "F81", "GTR"],
    help="default=ER for MK-REGIME: Rate-matrix structure fitted independently in each regime.",
)
pasr.add_argument(
    "--regime-parameters",
    "--regime_parameters",
    dest="regime_parameters",
    metavar="PATH",
    default=None,
    type=str,
    help="Continuous regime models: Optional complete TSV of fixed regime parameters. "
    "BMS uses sigma2; BMS-DRIFT uses sigma2/drift; OUM uses theta; "
    "OUMA adds alpha; OUMV adds sigma2; OUMVA uses all three.",
)
pasr.add_argument(
    "--hidden-categories",
    "--hidden_categories",
    dest="hidden_categories",
    metavar="INT",
    default=None,
    type=int,
    help="default=2 for HRM/COVARION: Number of latent rate classes; must be at least two.",
)
pasr.add_argument(
    "--mixture-model",
    "--mixture_model",
    dest="mixture_model",
    metavar="ER|SYM|ARD|F81|GTR",
    default=None,
    choices=["ER", "SYM", "ARD", "F81", "GTR"],
    help="default=ER for MK-MIXTURE: Shared base rate-matrix parameterization.",
)
pasr.add_argument(
    "--rate-mixture",
    "--rate_mixture",
    dest="rate_mixture",
    metavar="gamma|free",
    default=None,
    choices=["gamma", "free"],
    help="default=gamma for MK-MIXTURE: Across-character rate distribution.",
)
pasr.add_argument(
    "--rate-categories",
    "--rate_categories",
    dest="rate_categories",
    metavar="INT",
    default=None,
    type=int,
    help="default=4 for MK-MIXTURE: Number of discrete-gamma or ordered free-rate categories.",
)
pasr.add_argument(
    "--rate-matrix",
    "--rate_matrix",
    dest="rate_matrix",
    metavar="PATH",
    default=None,
    type=str,
    help="Discrete CUSTOM only: Labelled TSV Q matrix. The first column is state; "
    "zero diagonals are replaced by negative off-diagonal row sums.",
)
pasr.add_argument(
    "--thresholds",
    metavar="T1,T2,...",
    default=None,
    type=str,
    help="THRESHOLD only: Fixed increasing liability thresholds. If omitted, binary uses zero; ordinal fixes the first at zero and samples the rest.",
)
pasr.add_argument(
    "--liability-samples",
    "--liability_samples",
    dest="liability_samples",
    metavar="INT",
    default=None,
    type=int,
    help="default=1000 for THRESHOLD: Retained liability draws per MCMC chain.",
)
pasr.add_argument(
    "--liability-burnin",
    "--liability_burnin",
    dest="liability_burnin",
    metavar="INT",
    default=None,
    type=int,
    help="default=500 for THRESHOLD: Burn-in sweeps per chain.",
)
pasr.add_argument(
    "--liability-thin",
    "--liability_thin",
    dest="liability_thin",
    metavar="INT",
    default=None,
    type=int,
    help="default=1 for THRESHOLD: Sweeps between retained draws.",
)
pasr.add_argument(
    "--liability-chains",
    "--liability_chains",
    dest="liability_chains",
    metavar="INT",
    default=None,
    type=int,
    help="default=4 for THRESHOLD: Independent MCMC chains used for R-hat and ESS diagnostics.",
)
pasr.add_argument(
    "--liability-out",
    "--liability_out",
    dest="liability_out",
    metavar="PATH",
    default=None,
    type=str,
    help="THRESHOLD only: Optional TSV of posterior latent-liability moments.",
)
pasr.add_argument(
    "--root-prior",
    "--root_prior",
    dest="root_prior",
    metavar="equal|empirical|stationary|flat|fixed|gaussian",
    default=None,
    type=str,
    required=False,
    action="store",
    choices=["equal", "empirical", "stationary", "flat", "fixed", "gaussian"],
    help="model-specific default: Discrete equal (F81/GTR stationary), flat for BM-family "
    "models, and stationary for OU/OUM. OU also supports fixed or Gaussian roots. "
    "Stationary uses the fitted/fixed process equilibrium; "
    "independent of --input-rooted.",
)
pasr.add_argument(
    "--root-mean",
    "--root_mean",
    dest="root_mean",
    metavar="FLOAT",
    default=None,
    type=finite_float,
    help="OU with --root-prior fixed/gaussian: Required fixed root-state mean.",
)
pasr.add_argument(
    "--root-variance",
    "--root_variance",
    dest="root_variance",
    metavar="FLOAT",
    default=None,
    type=nonnegative_finite_float,
    help="OU with --root-prior gaussian: Required positive root-state variance.",
)
pasr.add_argument(
    "--ambiguous-separator",
    "--ambiguous_separator",
    dest="ambiguous_separator",
    metavar="STR",
    default=None,
    type=str,
    required=False,
    action="store",
    help='default=| for discrete: Separator for ambiguous or polymorphic states such as "A|B".',
)
pasr.add_argument(
    "--target",
    metavar="all|intnode,missing-leaf|leaf",
    default="all",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Comma-separated node classes to report. "
    'Use "missing-leaf" to report only imputed leaves. '
    'Legacy "tip" and "missing_tip" values remain accepted.',
)
pasr.add_argument(
    "--output",
    metavar="probabilities|map|summary",
    default=None,
    type=str,
    required=False,
    action="store",
    choices=["probabilities", "map", "summary"],
    help="default=probabilities for discrete, summary for continuous: "
    "Report discrete probabilities/MAP states or continuous means, variances, and intervals.",
)
pasr.add_argument(
    "--sigma2",
    metavar="FLOAT",
    default=None,
    type=nonnegative_finite_float,
    help="Scalar continuous models: Fixed diffusion variance rate per branch-length unit. "
    "If omitted, estimate it; MV-BM estimates a covariance matrix instead.",
)
pasr.add_argument(
    "--evolution-parameter",
    "--evolution_parameter",
    dest="evolution_parameter",
    metavar="FLOAT",
    default=None,
    type=finite_float,
    help="LAMBDA/KAPPA/DELTA/EB/ACDC only: Fix the model's branch-transform parameter; otherwise estimate it.",
)
pasr.add_argument(
    "--evolution-parameter-bounds",
    "--evolution_parameter_bounds",
    dest="evolution_parameter_bounds",
    metavar="MIN,MAX",
    default=None,
    type=str,
    help="LAMBDA/KAPPA/DELTA/EB/ACDC only: Bounds in natural parameter units for estimation.",
)
pasr.add_argument(
    "--alpha",
    metavar="FLOAT",
    default=None,
    type=nonnegative_finite_float,
    help="OU regime models: Fixed attraction strength per branch-length unit. If omitted, estimate it.",
)
pasr.add_argument(
    "--alpha-bounds",
    "--alpha_bounds",
    dest="alpha_bounds",
    metavar="MIN,MAX",
    default=None,
    type=str,
    help="OU models: Positive bounds for estimating alpha. Default scales to tree depth.",
)
pasr.add_argument(
    "--theta",
    metavar="FLOAT",
    default=None,
    type=finite_float,
    help="OU models: Fixed shared process optimum. Regime models can instead read regime optima from a table.",
)
pasr.add_argument(
    "--eb-rate",
    "--eb_rate",
    dest="eb_rate",
    metavar="FLOAT",
    default=None,
    type=finite_float,
    help="EB/ACDC alias for --evolution-parameter: Fixed exponential change in diffusion rate per branch-length unit.",
)
pasr.add_argument(
    "--eb-rate-bounds",
    "--eb_rate_bounds",
    dest="eb_rate_bounds",
    metavar="MIN,MAX",
    default=None,
    type=str,
    help="EB/ACDC alias for --evolution-parameter-bounds. EB defaults to declining rates only; ACDC permits acceleration.",
)
pasr.add_argument(
    "--drift",
    metavar="FLOAT",
    default=None,
    type=finite_float,
    help="BM-DRIFT: Fixed directional trend. BMS-DRIFT: Fix one shared drift or "
    "use regime parameters; otherwise estimate each regime drift.",
)
pasr.add_argument(
    "--standard-error-column",
    "--standard_error_column",
    dest="standard_error_column",
    metavar="COLUMN",
    default=None,
    type=str,
    help="Known non-negative measurement SEs. Multivariate models require a comma-separated SE column per trait.",
)
pasr.add_argument(
    "--ci-level",
    "--ci_level",
    dest="ci_level",
    metavar="FLOAT",
    default=None,
    type=finite_float,
    help="default=0.95 for continuous: Interval coverage strictly between zero and one. "
    "Intervals condition on fitted/fixed model parameters and the input tree; "
    "parameter/tree uncertainty is excluded.",
)
pasr.add_argument(
    "--profile-ci-level",
    "--profile_ci_level",
    dest="profile_ci_level",
    metavar="FLOAT",
    default=None,
    type=finite_float,
    help="LAMBDA/KAPPA/DELTA/EB/ACDC: Likelihood-ratio profile interval level for the fitted shape parameter.",
)
pasr.add_argument(
    "--model-out",
    "--model_out",
    dest="model_out",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Optional TSV reporting the selected model, parameters, "
    "likelihood convention, optimizer diagnostics, and uncertainty contract.",
)
pasr.add_argument(
    "--compare-models",
    "--compare_models",
    dest="compare_models",
    metavar="MODEL1,MODEL2,...",
    default=None,
    type=str,
    help="Fit compatible ASR models and calculate AIC, AICc, BIC, and criterion weights.",
)
pasr.add_argument(
    "--model-comparison-out",
    "--model_comparison_out",
    dest="model_comparison_out",
    metavar="PATH",
    default=None,
    type=str,
    help="Model-comparison TSV; required with --compare-models.",
)
pasr.add_argument(
    "--covariance-out",
    "--covariance_out",
    dest="covariance_out",
    metavar="PATH",
    default=None,
    type=str,
    help="MV-BM/MV-OU only: Optional tidy TSV of conditional trait covariances for selected nodes.",
)
pasr.add_argument(
    "--posterior-samples-out",
    "--posterior_samples_out",
    dest="posterior_samples_out",
    metavar="PATH",
    default=None,
    type=str,
    help="Scalar continuous models: Optional long-form TSV of joint all-node posterior draws.",
)
pasr.add_argument(
    "--posterior-samples",
    "--posterior_samples",
    dest="posterior_samples",
    metavar="INT",
    default=None,
    type=int,
    help="default=1000 with --posterior-samples-out: Number of joint posterior draws.",
)
pasr.add_argument(
    "--posterior-predictive-out",
    "--posterior_predictive_out",
    dest="posterior_predictive_out",
    metavar="PATH",
    default=None,
    type=str,
    help="Scalar continuous models: Optional posterior-predictive check TSV.",
)
pasr.add_argument(
    "--posterior-predictive-simulations",
    "--posterior_predictive_simulations",
    dest="posterior_predictive_simulations",
    metavar="INT",
    default=None,
    type=int,
    help="default=1000 with --posterior-predictive-out: Replicated datasets.",
)
pasr.add_argument(
    "--bootstrap-out",
    "--bootstrap_out",
    dest="bootstrap_out",
    metavar="PATH",
    default=None,
    type=str,
    help="Scalar continuous models: Optional parametric-bootstrap refit TSV.",
)
pasr.add_argument(
    "--bootstrap-simulations",
    "--bootstrap_simulations",
    dest="bootstrap_simulations",
    metavar="INT",
    default=None,
    type=int,
    help="default=100 with --bootstrap-out: Number of simulated refits.",
)
pasr.add_argument(
    "--tree-out",
    "--tree_out",
    dest="tree_out",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Optional Newick/NHX tree annotated with discrete states/probabilities "
    "or continuous means/uncertainties.",
)
pasr.add_argument(
    "--tree-outformat",
    "--tree_outformat",
    dest="tree_outformat",
    metavar="auto|INT",
    default="auto",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: ETE tree format for --tree-out.",
)
pasr.add_argument(
    "--tree-annotation",
    "--tree_annotation",
    dest="tree_annotation",
    metavar="state|probability|map|mean|summary|all",
    default=None,
    type=str,
    required=False,
    action="store",
    choices=["state", "probability", "map", "mean", "summary", "all"],
    help="default=map for discrete, summary for continuous: NHX properties written to --tree-out. "
    "Discrete accepts state/probability/map/all; continuous accepts mean/summary/all. "
    '"all" additionally includes observed values and, for discrete, per-state probabilities.',
)
pasr.add_argument(
    "--stochastic-map-out",
    "--stochastic_map_out",
    dest="stochastic_map_out",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Optional TSV summarizing sampled stochastic-map transition counts per branch.",
)
pasr.add_argument(
    "--n-sim",
    "--n_sim",
    dest="n_sim",
    metavar="INT",
    default=None,
    type=int,
    required=False,
    action="store",
    help="default=100 for discrete: Number of stochastic maps sampled when --stochastic-map-out is specified.",
)
pasr.add_argument(
    "--threads",
    metavar="INT",
    default=None,
    type=int,
    required=False,
    action="store",
    help="default=1 for discrete: Number of parallel workers used for stochastic mapping simulations.",
)
pasr.add_argument(
    "--seed",
    metavar="INT",
    default=None,
    type=int,
    required=False,
    action="store",
    help="default=%(default)s: Random seed for stochastic mapping, liability MCMC, or continuous simulation diagnostics.",
)
pasr.set_defaults(handler=command_asr)


def command_constrain(args):
    from nwkit.constrain import constrain_main

    constrain_main(args)


pconstrain = subparsers.add_parser(
    "constrain",
    help="Generate a species-tree-like Newick file for topological constraint",
    parents=[p_parent, p_download, p_species],
)
pconstrain.add_argument(
    "--species-list",
    "--species_list",
    dest="species_list",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Text file containing species names, one per line. "
    'Expected formats are "GENUS SPECIES", "GENUS_SPECIES", or "GENUS_SPECIES_OTHERINFO". '
    "Species-aware parsing uses the genus and species fields only, ignoring any remaining suffix. "
    'e.g., "Arabidopsis thaliana" and "Arabidopsis_thaliana_TAIR10"',
)
pconstrain.add_argument(
    "--taxid-tsv",
    "--taxid_tsv",
    dest="taxid_tsv",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help='default=%(default)s: TSV file containing species names in the "leaf_name" column and their NCBI Taxonomy IDs in the "taxid" column. '
    "When specified, the provided NCBI Taxonomy IDs are used instead of inferring them from species names. "
    "Either --species-list or --taxid-tsv must be specified, but not both. "
    "This option is currently compatible only with --backbone ncbi.",
)
pconstrain.add_argument(
    "--backbone",
    metavar="ncbi|ncbi_apgiv|ncbi_user|user",
    default="ncbi",
    type=str,
    required=False,
    action="store",
    choices=["ncbi", "ncbi_apgiv", "ncbi_user", "user"],
    help="default=%(default)s: The backbone for tree constraint. "
    '--infile is not required except for "user". '
    "ncbi: Infer NCBI Taxonomy ID from species name, and generate a tree based on the ranks. "
    "ncbi_apgiv: Infer NCBI Taxonomy ID from species name, and match it with the order-level angiosperm phylogeny in APG IV (https://doi.org/10.1111/boj.12385). "
    "ncbi_user: Infer NCBI Taxonomy ID from species name, and match the ranks with the labels of the user-provided tree. "
    "user: User-provided tree in --infile.",
)
pconstrain.add_argument(
    "--rank",
    metavar="no|species|genus|family|order|...",
    default="no",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Constrain at a particular taxonomic rank and above. "
    'For example, if "family" is specified, "genus" and "species" are not considered. '
    "This option is currently compatible only with --backbone ncbi",
)
pconstrain.add_argument(
    "--collapse",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help='default=%(default)s: For tip names of "GENUS_SPECIES_OTHERINFO", '
    "drop OTHERINFO and collapse clades if GENUS_SPECIES is identical. "
    "The output file may be used as a species tree for phylogeny reconciliation. ",
)
pconstrain.set_defaults(handler=command_constrain)


def command_collapse(args):
    from nwkit.collapse import collapse_main

    collapse_main(args)


pcollapse = subparsers.add_parser(
    "collapse",
    help="Collapse internal branches by support and/or branch length",
    parents=[p_parent],
)
pcollapse.add_argument(
    "--min-support",
    "--min_support",
    dest="min_support",
    metavar="FLOAT",
    default=None,
    type=finite_float,
    required=False,
    action="store",
    help="default=%(default)s: Collapse internal branches whose support is smaller than this threshold.",
)
pcollapse.add_argument(
    "--max-dist",
    "--max_dist",
    dest="max_dist",
    metavar="FLOAT",
    default=None,
    type=finite_float,
    required=False,
    action="store",
    help="default=%(default)s: Collapse internal branches whose branch length is at most this threshold.",
)
pcollapse.add_argument(
    "--preserve-branch-length",
    "--preserve_branch_length",
    dest="preserve_branch_length",
    metavar="yes|no",
    default="yes",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Add the deleted branch length to descendant branches when collapsing.",
)
pcollapse.set_defaults(handler=command_collapse)


def command_compose(args):
    from nwkit.compose import compose_main

    compose_main(args)


pcompose = subparsers.add_parser(
    "compose",
    help="Assemble compatible roots, values, and annotations from multiple trees",
    parents=[p_parent],
)
pcompose.add_argument(
    "--manifest",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="Optional JSON manifest defining root, name, support, length, and property sources.",
)
pcompose.add_argument(
    "--root-source",
    "--root_source",
    dest="root_source",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="Tree providing the root split.",
)
pcompose.add_argument(
    "--name-source",
    "--name_source",
    dest="name_source",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="Tree providing node names.",
)
pcompose.add_argument(
    "--support-source",
    "--support_source",
    dest="support_source",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="Tree providing internal-node support values.",
)
pcompose.add_argument(
    "--length-source",
    "--length_source",
    dest="length_source",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="Tree providing branch lengths.",
)
pcompose.add_argument(
    "--property-source",
    "--property_source",
    dest="property_source",
    metavar="SOURCE=TARGET@PATH",
    default=[],
    type=str,
    required=False,
    action="append",
    help="Tree providing an arbitrary NHX property. SOURCE@PATH preserves its name. May be repeated.",
)
pcompose.add_argument(
    "--source-format",
    "--source_format",
    dest="source_format",
    metavar="auto|auto-strict|INT",
    default="auto",
    type=str,
    required=False,
    action="store",
    help="Default ETE parser format for all source trees.",
)
pcompose.add_argument(
    "--taxon-mode",
    "--taxon_mode",
    dest="taxon_mode",
    metavar="exact|intersection",
    default="exact",
    type=str,
    required=False,
    action="store",
    choices=["exact", "intersection"],
    help="Require identical tip sets or map unique clades projected onto shared tips.",
)
pcompose.add_argument(
    "--match-basis",
    "--match_basis",
    dest="match_basis",
    metavar="clade|split",
    default="clade",
    type=str,
    required=False,
    action="store",
    choices=["clade", "split"],
    help="Match rooted descendant clades or root-independent canonical edge splits.",
)
pcompose.add_argument(
    "--root-edge-policy",
    "--root_edge_policy",
    dest="root_edge_policy",
    metavar="TARGET_PROPERTY=POLICY",
    default=[],
    type=str,
    required=False,
    action="append",
    help="Resolve root-edge ambiguity per target property. Policies: auto, skip, equal-only, matching-side, mean, min, max, edge-total. May be repeated.",
)
pcompose.add_argument(
    "--allow-projected-values",
    "--allow_projected_values",
    dest="allow_projected_values",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Permit support and branch-length transfer through projected matches. Use only when shared-tip equivalence is sufficient.",
)
pcompose.add_argument(
    "--policy",
    metavar="compatible-only|strict",
    default="compatible-only",
    type=str,
    required=False,
    action="store",
    choices=["compatible-only", "strict"],
    help="Skip incompatible values, or require exact (not projected) matches for every requested value.",
)
pcompose.add_argument(
    "--report",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="Optional per-clade TSV recording sources, matches, transferred values, and conflicts.",
)
pcompose.set_defaults(handler=command_compose)


def command_cladefreq(args):
    from nwkit.cladefreq import cladefreq_main

    cladefreq_main(args)


pcladefreq = subparsers.add_parser(
    "cladefreq",
    help="Summarize clade frequencies across a tree collection",
    parents=[p_tree_input, p_table_output],
)
pcladefreq.add_argument(
    "-r",
    "--reference",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Optional reference tree used to flag clades present in a named topology.",
)
pcladefreq.add_argument(
    "-rf",
    "--reference-format",
    "--reference_format",
    dest="reference_format",
    metavar="auto|auto-strict|INT",
    default="auto",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: ETE tree format for --reference.",
)
pcladefreq.add_argument(
    "--weight-tsv",
    "--weight_tsv",
    dest="weight_tsv",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Optional TSV assigning a positive weight to each input tree. "
    'The file must contain a "weight" column and may contain a 1-based "tree_id" column.',
)
pcladefreq.add_argument(
    "--support-scale",
    "--support_scale",
    dest="support_scale",
    metavar="percent|proportion",
    default="percent",
    type=str,
    required=False,
    action="store",
    choices=["percent", "proportion"],
    help="default=%(default)s: Scale used for the reported clade frequencies.",
)
pcladefreq.add_argument(
    "--threads",
    metavar="INT",
    default=1,
    type=int,
    required=False,
    action="store",
    help="default=%(default)s: Number of parallel workers used to parse and summarize input trees.",
)
pcladefreq.set_defaults(handler=command_cladefreq)


def command_consensus(args):
    from nwkit.consensus import consensus_main

    consensus_main(args)


pconsensus = subparsers.add_parser(
    "consensus",
    help="Generate a consensus tree or transfer consensus support to a reference tree",
    parents=[p_parent],
)
pconsensus.add_argument(
    "-r",
    "--reference",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Optional reference tree that receives consensus support values instead of building a de novo consensus tree.",
)
pconsensus.add_argument(
    "-rf",
    "--reference-format",
    "--reference_format",
    dest="reference_format",
    metavar="auto|auto-strict|INT",
    default="auto",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: ETE tree format for --reference.",
)
pconsensus.add_argument(
    "--method",
    metavar="greedy|majority|strict",
    default="greedy",
    type=str,
    required=False,
    action="store",
    choices=["greedy", "majority", "strict"],
    help="default=%(default)s: Consensus-clade selection rule. "
    "greedy uses --min-freq, majority keeps clades with frequency > 0.5, and strict keeps only clades present in all trees.",
)
pconsensus.add_argument(
    "--comparison",
    metavar="rooted|unrooted",
    default="rooted",
    type=str,
    choices=["rooted", "unrooted"],
    help="default=%(default)s: Treat input branches as rooted clades or root-independent unrooted splits.",
)
pconsensus.add_argument(
    "--min-freq",
    "--min_freq",
    dest="min_freq",
    metavar="0.0<=FLOAT<=1.0",
    default=0.5,
    type=unit_interval_float,
    required=False,
    action="store",
    help="default=%(default)s: Minimum clade frequency retained when --method greedy is used.",
)
pconsensus.add_argument(
    "--branch-length",
    "--branch_length",
    dest="branch_length",
    metavar="none|mean|median",
    default="none",
    type=str,
    required=False,
    action="store",
    choices=["none", "mean", "median"],
    help="default=%(default)s: How branch lengths are assigned in the de novo consensus tree.",
)
pconsensus.add_argument(
    "--weight-tsv",
    "--weight_tsv",
    dest="weight_tsv",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Optional TSV assigning a positive weight to each input tree. "
    'The file must contain a "weight" column and may contain a 1-based "tree_id" column.',
)
pconsensus.add_argument(
    "--support-scale",
    "--support_scale",
    dest="support_scale",
    metavar="percent|proportion",
    default="percent",
    type=str,
    required=False,
    action="store",
    choices=["percent", "proportion"],
    help="default=%(default)s: Scale used for consensus support values.",
)
pconsensus.add_argument(
    "--threads",
    metavar="INT",
    default=1,
    type=int,
    required=False,
    action="store",
    help="default=%(default)s: Number of parallel workers used to parse and summarize input trees.",
)
pconsensus.set_defaults(handler=command_consensus)


def command_contrast(args):
    from nwkit.contrast import contrast_main

    contrast_main(args)


pcontrast = subparsers.add_parser(
    "contrast",
    help="Calculate continuous-trait phylogenetic independent contrasts",
    parents=[p_tree_input, p_table_output, p_contrast_tip_table_policy],
)
pcontrast.add_argument(
    "--trait",
    metavar="PATH",
    default=None,
    type=str,
    required=True,
    action="store",
    help='TSV file containing a "leaf_name" column and numeric trait columns.',
)
pcontrast.add_argument(
    "--columns",
    metavar="COLUMN1,COLUMN2,...",
    default=None,
    type=str,
    required=True,
    action="store",
    help="Comma-separated numeric columns in --trait for which contrasts are calculated.",
)
pcontrast.add_argument(
    "--reconciliation",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Optional TSV from 'nwkit reconcile'. It supplies event annotations, species-branch mappings, lineage IDs, and contrast orientation.",
)
pcontrast.add_argument(
    "--tree-id",
    "--tree_id",
    dest="tree_id",
    metavar="TEXT",
    default="",
    type=str,
    required=False,
    action="store",
    help="default=empty: Stable gene-family/tree identifier. Required for unambiguous multi-tree aggregation.",
)
pcontrast.add_argument(
    "--event-type",
    "--event_type",
    dest="event_type",
    metavar="all|speciation|duplication|transfer|unresolved",
    default="all",
    type=str,
    required=False,
    action="store",
    choices=["all", "speciation", "duplication", "transfer", "unresolved"],
    help="default=%(default)s: With --reconciliation, report only this event type.",
)
pcontrast.add_argument(
    "--eligible-only",
    "--eligible_only",
    dest="eligible_only",
    metavar="auto|yes|no",
    default=None,
    type=strtoautobool,
    required=False,
    action="store",
    help="default=auto: With --reconciliation, retain only eligible rows. Auto means yes for --event-type speciation and no for all other event selections.",
)
pcontrast.add_argument(
    "--speciation-coverage",
    "--speciation_coverage",
    dest="speciation_coverage",
    metavar="complete|any",
    default="complete",
    type=str,
    required=False,
    action="store",
    choices=["complete", "any"],
    help="default=%(default)s: When eligible rows are requested, require complete sampling of both species-tree daughter clades or allow explicitly reported partial coverage.",
)
pcontrast.add_argument(
    "--branch-length",
    "--branch_length",
    dest="branch_length",
    metavar="original|unit",
    default="original",
    type=str,
    required=False,
    action="store",
    choices=["original", "unit"],
    help="default=%(default)s: Use positive input branch lengths or replace every non-root branch length with one.",
)
pcontrast.add_argument(
    "--evolution-model",
    "--evolution_model",
    dest="evolution_model",
    metavar="MODEL",
    default="brownian",
    type=str,
    choices=list(CONTRAST_EVOLUTION_MODELS),
    help="default=%(default)s: Evolutionary model represented as Brownian motion on a transformed tree.",
)
pcontrast.add_argument(
    "--evolution-parameter",
    "--evolution_parameter",
    dest="evolution_parameter",
    metavar="FLOAT",
    default=None,
    type=finite_float,
    help="Fixed shape parameter required by parameterized contrast evolution models.",
)
pcontrast.add_argument(
    "--biological-id",
    "--biological_id",
    dest="biological_id",
    metavar="COLUMN",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Column identifying independent biological observations. Its presence enables replicate-aware contrasts.",
)
pcontrast.add_argument(
    "--technical-id",
    "--technical_id",
    dest="technical_id",
    metavar="COLUMN",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Optional technical-replicate identifier nested within leaf and biological observation.",
)
pcontrast.add_argument(
    "--technical-aggregation",
    "--technical_aggregation",
    dest="technical_aggregation",
    metavar="error|mean",
    default="error",
    type=str,
    required=False,
    action="store",
    choices=["error", "mean"],
    help="default=%(default)s: Reject technical replicates or explicitly average already transformed continuous values within each biological observation.",
)
pcontrast.add_argument(
    "--batch",
    metavar="COLUMN",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Optional categorical batch column fitted as a fixed observation-level effect; confounded designs are rejected.",
)
pcontrast.add_argument(
    "--within-variance",
    "--within_variance",
    dest="within_variance",
    metavar="pooled|leaf|known-se",
    default="pooled",
    type=str,
    required=False,
    action="store",
    choices=["pooled", "leaf", "known-se"],
    help="default=%(default)s: Estimate a pooled or leaf-specific biological variance, or read known standard errors.",
)
pcontrast.add_argument(
    "--standard-error-columns",
    "--standard_error_columns",
    dest="standard_error_columns",
    metavar="COLUMN1,COLUMN2,...",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=<TRAIT>_se: Comma-separated known-SE columns matching --columns when --within-variance known-se is selected.",
)
pcontrast.add_argument(
    "--sample-size-columns",
    "--sample_size_columns",
    dest="sample_size_columns",
    metavar="COLUMN1,COLUMN2,...",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Optional positive-integer sample-size columns matching --columns for known-SE input and audit output.",
)
pcontrast.add_argument(
    "--sampling-covariance-out",
    "--sampling_covariance_out",
    dest="sampling_covariance_out",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Required replicate-aware long-form TSV containing propagated contrast sampling covariance.",
)
pcontrast.add_argument(
    "--tip-summary-out",
    "--tip_summary_out",
    dest="tip_summary_out",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Optional audit TSV of leaf means, biological sample sizes, within-leaf SDs, and standard errors.",
)
pcontrast.set_defaults(handler=command_contrast)


def command_diff(args):
    from nwkit.diff import diff_main

    diff_main(args)


pdiff = subparsers.add_parser(
    "diff",
    help="Report interpretable clade, root, value, and annotation differences between trees",
    parents=[p_tree_input, p_table_output],
)
pdiff.add_argument(
    "-i2",
    "--infile2",
    metavar="PATH",
    default=None,
    type=str,
    required=True,
    action="store",
    help="Second Newick tree.",
)
pdiff.add_argument(
    "-f2",
    "--format2",
    metavar="auto|auto-strict|INT",
    default="auto",
    type=str,
    required=False,
    action="store",
    help="ETE parser format for --infile2.",
)
pdiff.add_argument(
    "--taxon-mode",
    "--taxon_mode",
    dest="taxon_mode",
    metavar="exact|intersection",
    default="exact",
    type=str,
    required=False,
    action="store",
    choices=["exact", "intersection"],
    help="Require identical tip sets or compare unique projections onto shared tips.",
)
pdiff.add_argument(
    "--comparison",
    metavar="rooted|unrooted",
    default="rooted",
    type=str,
    required=False,
    action="store",
    choices=["rooted", "unrooted"],
    help="Compare rooted descendant clades or root-independent edge splits.",
)
pdiff.add_argument(
    "--target",
    metavar="all|root|leaf|intnode",
    default="intnode",
    type=str,
    required=False,
    action="store",
    choices=["all", "root", "leaf", "intnode"],
    help="Node class included in detailed rows.",
)
pdiff.add_argument(
    "--property",
    metavar="KEY",
    default=[],
    type=str,
    required=False,
    action="append",
    help="Additional NHX property compared and serialized as JSON. May be repeated.",
)
pdiff.add_argument(
    "--fail-on-difference",
    "--fail_on_difference",
    dest="fail_on_difference",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="Exit nonzero after writing the report if any difference is detected.",
)
pdiff.set_defaults(handler=command_diff)


def command_dist(args):
    from nwkit.dist import dist_main

    dist_main(args)


pdist = subparsers.add_parser(
    "dist",
    help="Calculate topology and branch-length distances between two trees",
    parents=[p_tree_input, p_table_output],
)
pdist.add_argument(
    "-i2",
    "--infile2",
    metavar="PATH",
    default=None,
    type=str,
    required=True,
    action="store",
    help="default=%(default)s: Input newick file 2.",
)
pdist.add_argument(
    "-f2",
    "--format2",
    metavar="auto|auto-strict|INT",
    default="auto",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: ETE tree format for --infile2.",
)
pdist.add_argument(
    "--metric",
    metavar="METRIC[,METRIC...]",
    default=None,
    type=str,
    required=False,
    action="append",
    help="default=all: Distance metric. May be repeated or comma-separated. "
    "Choices: all, rf, normalized-rf, weighted-rf, branch-score, "
    "path-topological, path-length.",
)
pdist.add_argument(
    "--comparison",
    metavar="rooted|unrooted",
    default="rooted",
    type=str,
    required=False,
    action="store",
    choices=["rooted", "unrooted"],
    help="default=%(default)s: Compare rooted descendant clades or root-independent edge splits. "
    "Leaf-to-leaf path metrics are root-independent.",
)
pdist.add_argument(
    "-d",
    "--dist",
    dest="dist",
    metavar="STR",
    default=None,
    type=str,
    required=False,
    action="store",
    help=argparse.SUPPRESS,
)
pdist.set_defaults(handler=command_dist)


def command_draw(args):
    from nwkit.draw import draw_main

    draw_main(args)


pdraw = subparsers.add_parser(
    "draw",
    help="Draw a phylogenetic tree with tip images, support, categorical properties, and node-probability annotations",
    parents=[p_tree_input, p_species, p_tip_table_policy],
)
pdraw.add_argument(
    "-o",
    "--outfile",
    metavar="PATH",
    default="nwkit_draw.pdf",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Output image file.",
)
pdraw.add_argument(
    "--image-format",
    "--image_format",
    dest="image_format",
    metavar="auto|pdf|png|svg",
    default="auto",
    type=str,
    required=False,
    action="store",
    choices=["auto", "pdf", "png", "svg"],
    help='default=%(default)s: Output image format. "auto" infers from --outfile and otherwise falls back to pdf.',
)
pdraw.add_argument(
    "--layout",
    metavar="LAYOUT",
    default="rectangular",
    type=str,
    choices=[
        "rectangular",
        "slanted",
        "cladogram",
        "circular",
        "radial",
        "unrooted",
        "spiral",
        "fractal",
    ],
    help="default=%(default)s: Tree geometry; subtree placement and label-aware spacing are controlled independently.",
)
pdraw.add_argument(
    "--subtree-packing",
    "--subtree_packing",
    dest="subtree_packing",
    metavar="standard|tidy",
    default="standard",
    type=str,
    choices=["standard", "tidy"],
    help='default=%(default)s: Subtree placement strategy. "tidy" compacts rectangular, circular, and spiral layouts while preserving tip order.',
)
pdraw.add_argument(
    "--spiral-turns",
    "--spiral_turns",
    dest="spiral_turns",
    metavar="FLOAT",
    default=None,
    type=finite_float,
    help="default=auto: Number of turns used by --layout spiral.",
)
pdraw.add_argument(
    "--angular-span",
    "--angular_span",
    dest="angular_span",
    metavar="DEGREES",
    default=360.0,
    type=finite_float,
    help="default=%(default)s: Angular span occupied by circular or radial geometry; 180 uses a semicircle.",
)
pdraw.add_argument(
    "--angular-center",
    "--angular_center",
    dest="angular_center",
    metavar="DEGREES",
    default=90.0,
    type=finite_float,
    help="default=%(default)s: Direction at the center of a circular or radial sector; 90 centers a 180-degree sector on the upper half-plane.",
)
pdraw.add_argument(
    "--unrooted-method",
    "--unrooted_method",
    dest="unrooted_method",
    metavar="equal-angle|equal-daylight",
    default="equal-angle",
    type=str,
    choices=["equal-angle", "equal-daylight"],
    help="default=%(default)s: Angular optimization used by --layout unrooted.",
)
pdraw.add_argument(
    "--daylight-iterations",
    "--daylight_iterations",
    dest="daylight_iterations",
    metavar="INT",
    default=5,
    type=int,
    help="default=%(default)s: Maximum equal-daylight refinement passes.",
)
pdraw.add_argument(
    "--max-visible-tips",
    "--max_visible_tips",
    dest="max_visible_tips",
    metavar="INT",
    default=None,
    type=int,
    help="default=None: Automatically collapse clades in a drawing-only copy until no more than INT tips remain visible.",
)
pdraw.add_argument(
    "--collapse-label",
    "--collapse_label",
    dest="collapse_label",
    metavar="TEXT",
    default=None,
    type=str,
    help='default="{first}…{last} (n={tips})": Label template for automatically collapsed clades; fields are clade, first, last, and tips.',
)
pdraw.add_argument(
    "--collapse-property-aggregation",
    "--collapse_property_aggregation",
    dest="collapse_property_aggregation",
    metavar="none|mean",
    default="none",
    type=str,
    choices=["none", "mean"],
    help="default=%(default)s: Preserve only complete constant properties on collapsed clades; mean explicitly averages complete numeric properties.",
)
pdraw.add_argument(
    "--species-overlap-node-plot",
    "--species_overlap_node_plot",
    dest="species_overlap_node_plot",
    metavar="yes|no|auto",
    default="auto",
    required=False,
    type=str,
    choices=["yes", "no", "auto"],
    help="default=%(default)s: Show speciation/duplication node markers. "
    '"yes" always tries to plot, "no" disables them, and "auto" enables them only when '
    "all tip labels are parseable by the configured species parser.",
)
pdraw.add_argument(
    "--ladderize",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Ladderize the input tree before drawing.",
)
pdraw.add_argument(
    "--trait",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help='default=%(default)s: Optional TSV containing "leaf_name" and one or more grouping columns for tip-label coloring.',
)
pdraw.add_argument(
    "--group-by",
    "--group_by",
    dest="group_by",
    metavar="STR",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Column name in --trait used to color tip labels.",
)
pdraw.add_argument(
    "--trait-palette",
    "--trait_palette",
    dest="trait_palette",
    metavar="STR",
    default="tab10",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Matplotlib colormap name used for trait-group colors.",
)
pdraw.add_argument(
    "--support-labels",
    "--support_labels",
    dest="support_labels",
    metavar="yes|no",
    default="yes",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Whether to draw support labels on branches.",
)
pdraw.add_argument(
    "--support-min",
    "--support_min",
    dest="support_min",
    metavar="FLOAT",
    default=None,
    type=finite_float,
    required=False,
    action="store",
    help="default=%(default)s: Show only support labels greater than or equal to this value.",
)
pdraw.add_argument(
    "--figure-width",
    "--figure_width",
    dest="figure_width",
    metavar="FLOAT",
    default=3.6,
    type=finite_float,
    required=False,
    action="store",
    help="default=%(default)s: Figure width in inches.",
)
pdraw.add_argument(
    "--figure-height",
    "--figure_height",
    dest="figure_height",
    metavar="FLOAT",
    default=None,
    type=finite_float,
    required=False,
    action="store",
    help="default=auto: Optional fixed figure height in inches. By default height follows the number of tips.",
)
pdraw.add_argument(
    "--label-panel-width",
    "--label_panel_width",
    dest="label_panel_width",
    metavar="FLOAT",
    default=None,
    type=finite_float,
    required=False,
    action="store",
    help="default=auto: Width in inches reserved for tip labels and badges.",
)
pdraw.add_argument(
    "--tip-image-manifest",
    "--tip_image_manifest",
    dest="tip_image_manifest",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=None: TSV produced by nwkit image, containing leaf_name and local_path columns. "
    "When a leaf has multiple rows, the first candidate is used. "
    "Matching images are drawn in an aligned column beside the tree.",
)
pdraw.add_argument(
    "--tip-image-root",
    "--tip_image_root",
    dest="tip_image_root",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=manifest directory: Base directory for relative local_path values in --tip-image-manifest. "
    "Required when the manifest is read from STDIN.",
)
pdraw.add_argument(
    "--tip-image-size",
    "--tip_image_size",
    dest="tip_image_size",
    metavar="FLOAT",
    default=18.0,
    type=finite_float,
    required=False,
    action="store",
    help="default=%(default)s: Maximum width or height of each tip image in points.",
)
pdraw.add_argument(
    "--tip-image-gap",
    "--tip_image_gap",
    dest="tip_image_gap",
    metavar="FLOAT",
    default=4.0,
    type=finite_float,
    required=False,
    action="store",
    help="default=%(default)s: Horizontal padding on each side of the aligned tip-image column in points.",
)
pdraw.add_argument(
    "--font-size",
    "--font_size",
    dest="font_size",
    metavar="FLOAT",
    default=8.0,
    type=finite_float,
    required=False,
    action="store",
    help="default=%(default)s: Font size in points.",
)
pdraw.add_argument(
    "--font-family",
    "--font_family",
    dest="font_family",
    metavar="STR",
    default="Helvetica",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Matplotlib font family used for labels and legends.",
)
pdraw.add_argument(
    "--branch-color",
    "--branch_color",
    dest="branch_color",
    metavar="COLOR",
    default="#000000",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Matplotlib color for branches.",
)
pdraw.add_argument(
    "--branch-width",
    "--branch_width",
    dest="branch_width",
    metavar="FLOAT",
    default=0.8,
    type=finite_float,
    help="default=%(default)s: Base branch width in points.",
)
pdraw.add_argument(
    "--terminal-branch-color",
    "--terminal_branch_color",
    dest="terminal_branch_color",
    metavar="COLOR",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=None: Optional color applied to terminal branches.",
)
pdraw.add_argument(
    "--branch-color-property",
    "--branch_color_property",
    dest="branch_color_property",
    metavar="PROPERTY",
    default=None,
    type=str,
    help="default=None: Map a categorical Newick/NHX branch property to color.",
)
pdraw.add_argument(
    "--branch-width-property",
    "--branch_width_property",
    dest="branch_width_property",
    metavar="PROPERTY",
    default=None,
    type=str,
    help="default=None: Map a numeric Newick/NHX branch property to width.",
)
pdraw.add_argument(
    "--branch-width-range",
    "--branch_width_range",
    dest="branch_width_range",
    metavar="MIN,MAX",
    default="0.4,2.5",
    type=str,
    help="default=%(default)s: Output width range for --branch-width-property.",
)
pdraw.add_argument(
    "--scale-bar",
    "--scale_bar",
    dest="scale_bar",
    metavar="none|auto|FLOAT",
    default="none",
    type=str,
    help="default=%(default)s: Add an exact scale to a layout with directly measurable branch-length segments.",
)
pdraw.add_argument(
    "--depth-guide",
    "--depth_guide",
    dest="depth_guide",
    metavar="none|auto|FLOAT",
    default="none",
    type=str,
    help="default=%(default)s: Add root-to-node depth guides to slanted, radial, or spiral; FLOAT sets the tick interval.",
)
pdraw.add_argument(
    "--branch-length-unit",
    "--branch_length_unit",
    dest="branch_length_unit",
    metavar="TEXT",
    default="",
    type=str,
    help='default="": Unit appended to scale-bar and depth-guide labels, for example substitutions/site or Ma.',
)
pdraw.add_argument(
    "--time-constraints",
    "--time_constraints",
    dest="time_constraints",
    metavar="auto|yes|no",
    default="auto",
    type=str,
    choices=["auto", "yes", "no"],
    help="default=%(default)s: Draw MCMCtree point, bounded, lower, and upper calibration glyphs when present.",
)
pdraw.add_argument(
    "--time-credible-intervals",
    "--time_credible_intervals",
    dest="time_credible_intervals",
    metavar="auto|yes|no",
    default="auto",
    type=str,
    choices=["auto", "yes", "no"],
    help="default=%(default)s: Draw node-age credible-interval whiskers from MCMCtree/FigTree annotations.",
)
pdraw.add_argument(
    "--mcmctree-posterior",
    "--mcmctree_posterior",
    dest="mcmctree_posterior",
    metavar="PATH",
    default=None,
    type=str,
    help="default=None: Read MCMCtree mcmc.txt ages on the input topology for dated-tree and DensiTree rendering.",
)
pdraw.add_argument(
    "--densitree-trees",
    "--densitree_trees",
    dest="densitree_trees",
    metavar="PATH",
    default=None,
    type=str,
    help="default=None: Overlay a posterior or bootstrap multi-Newick tree collection, including topology variation, against the input reference tree.",
)
pdraw.add_argument(
    "--posterior-point",
    "--posterior_point",
    dest="posterior_point",
    metavar="mean|median",
    default="mean",
    type=str,
    choices=["mean", "median"],
    help="default=%(default)s: Posterior point age used by --mcmctree-posterior.",
)
pdraw.add_argument(
    "--posterior-ci",
    "--posterior_ci",
    dest="posterior_ci",
    metavar="hpd|equal-tail",
    default="hpd",
    type=str,
    choices=["hpd", "equal-tail"],
    help="default=%(default)s: Node-age interval used by --mcmctree-posterior.",
)
pdraw.add_argument(
    "--posterior-ci-level",
    "--posterior_ci_level",
    dest="posterior_ci_level",
    metavar="0<FLOAT<1",
    default=0.95,
    type=finite_float,
    help="default=%(default)s: Credible mass for dated-tree node intervals.",
)
pdraw.add_argument(
    "--posterior-burnin",
    "--posterior_burnin",
    dest="posterior_burnin",
    metavar="INT",
    default=0,
    type=int,
    help="default=%(default)s: Leading MCMC rows or sampled trees to discard.",
)
pdraw.add_argument(
    "--posterior-thin",
    "--posterior_thin",
    dest="posterior_thin",
    metavar="INT",
    default=1,
    type=int,
    help="default=%(default)s: Keep every INT-th MCMC row or sampled tree.",
)
pdraw.add_argument(
    "--densitree",
    metavar="none|all|ci|both",
    default="none",
    type=str,
    choices=["none", "all", "ci", "both"],
    help="default=%(default)s: Draw retained trees, branchwise empirical path envelopes, or both from an MCMCtree age table or a multi-Newick tree collection; topology-varying samples are stratified before envelopes are calculated.",
)
pdraw.add_argument(
    "--densitree-alpha",
    "--densitree_alpha",
    dest="densitree_alpha",
    metavar="0<FLOAT<=1",
    default=0.035,
    type=finite_float,
    help="default=%(default)s: Opacity of each posterior tree in all/both mode.",
)
pdraw.add_argument(
    "--densitree-color",
    "--densitree_color",
    dest="densitree_color",
    metavar="COLOR",
    default="#0072B2",
    type=str,
    help="default=%(default)s: Color of posterior trees in all/both mode.",
)
pdraw.add_argument(
    "--densitree-ci-level",
    "--densitree_ci_level",
    dest="densitree_ci_level",
    metavar="0<FLOAT<1",
    default=0.95,
    type=finite_float,
    help="default=%(default)s: Central fraction of whole sampled paths retained in each within-topology branch envelope.",
)
pdraw.add_argument(
    "--densitree-ci-alpha",
    "--densitree_ci_alpha",
    dest="densitree_ci_alpha",
    metavar="0<FLOAT<=1",
    default=0.18,
    type=finite_float,
    help="default=%(default)s: Maximum opacity of branch-envelope polygons; topology groups are scaled by relative sample frequency.",
)
pdraw.add_argument(
    "--densitree-ci-color",
    "--densitree_ci_color",
    dest="densitree_ci_color",
    metavar="COLOR",
    default="#56B4E9",
    type=str,
    help="default=%(default)s: Color of branch-envelope polygons.",
)
pdraw.add_argument(
    "--tip-labels",
    "--tip_labels",
    dest="tip_labels",
    metavar="yes|no",
    default="yes",
    type=strtobool,
    help="default=%(default)s: Whether to draw tip labels. Disabling them is useful for dense overview layouts.",
)
pdraw.add_argument(
    "--tip-label-position",
    "--tip_label_position",
    dest="tip_label_position",
    metavar="auto|aligned|branch-end",
    default="auto",
    type=str,
    required=False,
    action="store",
    choices=["auto", "aligned", "branch-end"],
    help="default=%(default)s: Align labels at the right edge for rectangular layout and place them beside branch endpoints otherwise.",
)
pdraw.add_argument(
    "--tip-label-wrap",
    "--tip_label_wrap",
    dest="tip_label_wrap",
    metavar="none|auto|taxonomy|INT",
    default="none",
    type=str,
    help="default=%(default)s: Display-only wrapping; taxonomy preserves an underscore-delimited genus_species binomial on one line.",
)
pdraw.add_argument(
    "--tip-spacing",
    "--tip_spacing",
    dest="tip_spacing",
    metavar="uniform|label-aware",
    default="uniform",
    type=str,
    choices=["uniform", "label-aware"],
    help="default=%(default)s: Allocate tip rows or angular sectors uniformly, or from measured label and annotation heights.",
)
pdraw.add_argument(
    "--tip-label-font-style",
    "--tip_label_font_style",
    dest="tip_label_font_style",
    metavar="plain|italic|taxonomy",
    default="plain",
    type=str,
    choices=["plain", "italic", "taxonomy"],
    help="default=%(default)s: Typography for tip labels; taxonomy italicizes exact genus_species binomials conservatively.",
)
pdraw.add_argument(
    "--tip-track",
    "--tip_track",
    dest="tip_track",
    metavar="PROPERTY",
    default=[],
    type=str,
    action="append",
    help="default=[]: Add a categorical or continuous tip annotation track. May be repeated.",
)
pdraw.add_argument(
    "--tip-track-type",
    "--tip_track_type",
    dest="tip_track_type",
    metavar="auto|categorical|continuous",
    default="auto",
    type=str,
    choices=["auto", "categorical", "continuous"],
    help="default=%(default)s: Shared interpretation of values in --tip-track properties.",
)
pdraw.add_argument(
    "--tip-track-size",
    "--tip_track_size",
    dest="tip_track_size",
    metavar="FLOAT",
    default=5.0,
    type=finite_float,
    help="default=%(default)s: Width and height of each tip-track tile in points.",
)
pdraw.add_argument(
    "--tip-track-palette",
    "--tip_track_palette",
    dest="tip_track_palette",
    metavar="STR",
    default="viridis",
    type=str,
    help="default=%(default)s: Matplotlib colormap for continuous tip tracks; categorical tracks use --trait-palette.",
)
pdraw.add_argument(
    "--root-marker",
    "--root_marker",
    dest="root_marker",
    metavar="none|circle|diamond",
    default="none",
    type=str,
    required=False,
    action="store",
    choices=["none", "circle", "diamond"],
    help="default=%(default)s: Optional marker drawn at the displayed root.",
)
pdraw.add_argument(
    "--root-marker-color",
    "--root_marker_color",
    dest="root_marker_color",
    metavar="COLOR",
    default="#0072B2",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Matplotlib color for --root-marker.",
)
pdraw.add_argument(
    "--root-marker-size",
    "--root_marker_size",
    dest="root_marker_size",
    metavar="FLOAT",
    default=None,
    type=finite_float,
    required=False,
    action="store",
    help="default=font size: Diameter of --root-marker in points.",
)
pdraw.add_argument(
    "--tip-badge-property",
    "--tip_badge_property",
    dest="tip_badge_property",
    metavar="PROPERTY",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=None: Display the value of this Newick/NHX property as a badge beside each matching tip.",
)
pdraw.add_argument(
    "--tip-badge-missing-label",
    "--tip_badge_missing_label",
    dest="tip_badge_missing_label",
    metavar="TEXT",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=None: Badge text used when --tip-badge-property is absent or empty; missing badges are otherwise omitted.",
)
pdraw.add_argument(
    "--node-pie-properties",
    "--node_pie_properties",
    dest="node_pie_properties",
    metavar="PROP1,PROP2,...",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=None: Draw node pie markers from comma-separated numeric Newick/NHX properties.",
)
pdraw.add_argument(
    "--node-pie-target",
    "--node_pie_target",
    dest="node_pie_target",
    metavar="root,intnode,leaf|all",
    default="root,intnode",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Comma-separated node classes receiving --node-pie-properties.",
)
pdraw.add_argument(
    "--node-pie-leaf-filter",
    "--node_pie_leaf_filter",
    dest="node_pie_leaf_filter",
    metavar="PROPERTY:OP:VALUE",
    default=[],
    type=str,
    required=False,
    action="append",
    help="default=[]: Restrict leaf pies with OP in ge, gt, le, lt, eq, or ne; root and internal pies remain unaffected. May be repeated.",
)
pdraw.add_argument(
    "--node-label-property",
    "--node_label_property",
    dest="node_label_property",
    metavar="PROPERTY",
    default=None,
    type=str,
    required=False,
    action="store",
    help='default=None: Display this Newick/NHX property beside matching nodes; use "name" for ordinary Newick node labels.',
)
pdraw.add_argument(
    "--node-label-target",
    "--node_label_target",
    dest="node_label_target",
    metavar="root,intnode,leaf|all",
    default="intnode",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Comma-separated node classes eligible for --node-label-property.",
)
pdraw.add_argument(
    "--node-label-filter",
    "--node_label_filter",
    dest="node_label_filter",
    metavar="PROPERTY:OP:VALUE",
    default=[],
    type=str,
    required=False,
    action="append",
    help="default=[]: Filter property labels with OP in ge, gt, le, lt, eq, or ne. May be repeated.",
)
pdraw.add_argument(
    "--node-label-decimals",
    "--node_label_decimals",
    dest="node_label_decimals",
    metavar="INT",
    default=2,
    type=int,
    required=False,
    action="store",
    help="default=%(default)s: Decimal places used for numeric --node-label-property values.",
)
pdraw.add_argument(
    "--node-label-prefix",
    "--node_label_prefix",
    dest="node_label_prefix",
    metavar="TEXT",
    default="",
    type=str,
    required=False,
    action="store",
    help='default="": Text prepended to each displayed --node-label-property value.',
)
pdraw.add_argument(
    "--property-color",
    "--property_color",
    dest="property_color",
    metavar="VALUE=COLOR",
    default=[],
    type=str,
    required=False,
    action="append",
    help="default=[]: Color assigned to a tip-badge value or node-pie property/state. May be repeated.",
)
pdraw.add_argument(
    "--legend",
    metavar="yes|no",
    default="yes",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Whether to draw legends for enabled categorical layers.",
)
pdraw.add_argument(
    "--legend-columns",
    "--legend_columns",
    dest="legend_columns",
    metavar="auto|INT",
    default="auto",
    type=str,
    help="default=%(default)s: Number of legend columns when the legend is above the tree.",
)
pdraw.add_argument(
    "--legend-position",
    "--legend_position",
    dest="legend_position",
    metavar="auto|top|right",
    default="auto",
    type=str,
    choices=["auto", "top", "right"],
    help="default=%(default)s: Place dense legends to the right automatically or force top/right placement.",
)
pdraw.add_argument(
    "--collision-policy",
    "--collision_policy",
    dest="collision_policy",
    metavar="resolve|warn|error|ignore",
    default="resolve",
    type=str,
    choices=["resolve", "warn", "error", "ignore"],
    help="default=%(default)s: Resolve movable annotation collisions, report them, reject them, or leave them unchanged.",
)
pdraw.add_argument(
    "--layout-report",
    "--layout_report",
    dest="layout_report",
    metavar="PATH",
    default=None,
    type=str,
    help="default=None: Write a reproducible JSON report of layout choices, annotation collisions, branch crossings, wrapping, and collapsing; use - for STDOUT.",
)
pdraw.add_argument(
    "--transparent",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Save the figure with a transparent background.",
)
pdraw.set_defaults(handler=command_draw)


def command_drop(args):
    from nwkit.drop import drop_main

    drop_main(args)


pdrop = subparsers.add_parser(
    "drop",
    help="Remove node and branch information",
    parents=[p_parent, p_preserve_properties],
)
pdrop.add_argument(
    "-t",
    "--target",
    metavar="all|root|leaf|intnode",
    default="all",
    type=str,
    required=False,
    action="store",
    choices=["all", "root", "leaf", "intnode"],
    help="default=%(default)s: Nodes to be edited.",
)
pdrop.add_argument(
    "--name",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Drop node names.",
)
pdrop.add_argument(
    "--support",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Drop support values.",
)
pdrop.add_argument(
    "--length",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Drop branch length.",
)
pdrop.add_argument(
    "--fill",
    metavar="STR/NUMERIC",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Fill values instead of simply dropping them.",
)
pdrop.set_defaults(handler=command_drop)


def command_info(args):
    from nwkit.info import info_main

    info_main(args)


pinfo = subparsers.add_parser(
    "info",
    help="Print tree information",
    parents=[p_tree_input, p_text_output, p_species],
)
pinfo.set_defaults(handler=command_info)


def command_image(args):
    from nwkit.image import image_main

    image_main(args)


pimage = subparsers.add_parser(
    "image",
    help="Retrieve representative species images with license-aware filtering",
    parents=[p_tree_input, p_download, p_species],
)
pimage.add_argument(
    "--out-dir",
    "--out_dir",
    dest="out_dir",
    metavar="PATH",
    default=None,
    type=str,
    required=True,
    action="store",
    help="Output directory for images and metadata.",
)
pimage.add_argument(
    "--style",
    metavar="auto|photo|silhouette",
    default="auto",
    type=str,
    required=False,
    action="store",
    choices=["auto", "photo", "silhouette"],
    help="default=%(default)s: Image preference mode.",
)
pimage.add_argument(
    "--source",
    metavar="phylopic,bioicons,inaturalist,wikimedia,gbif,eol,idigbio,openverse,ncbi",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=style-dependent: Comma-separated provider priority list. "
    "Supported sources are phylopic, bioicons, inaturalist, wikimedia, gbif, eol, idigbio, openverse, and ncbi. "
    "When ncbi appears after other providers, it is used as a lazy fallback because it may download the NCBI taxonomy image dump.",
)
pimage.add_argument(
    "--license-max",
    "--license_max",
    dest="license_max",
    metavar="public-domain|cc-by|cc-by-sa|cc-by-nc|cc-by-nc-sa|any",
    default="cc-by-nc-sa",
    type=str,
    required=False,
    action="store",
    choices=["public-domain", "cc-by", "cc-by-sa", "cc-by-nc", "cc-by-nc-sa", "any"],
    help="default=%(default)s: Strongest restriction level allowed for candidate images.",
)
pimage.add_argument(
    "--allow-nd",
    "--allow_nd",
    dest="allow_nd",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Whether ND licenses are acceptable.",
)
pimage.add_argument(
    "--fallback-rank",
    "--fallback_rank",
    dest="fallback_rank",
    metavar="none|genus|family",
    default="none",
    type=str,
    required=False,
    action="store",
    choices=["none", "genus", "family"],
    help="default=%(default)s: Allow fallback matching above the species rank.",
)
pimage.add_argument(
    "--max-per-species",
    "--max_per_species",
    dest="max_per_species",
    metavar="INT",
    default=1,
    type=int,
    required=False,
    action="store",
    help="default=%(default)s: Number of images to keep per species.",
)
pimage.add_argument(
    "--max-download-bytes",
    "--max_download_bytes",
    dest="max_download_bytes",
    metavar="INT",
    default=104857600,
    type=int,
    required=False,
    action="store",
    help="default=%(default)s: Reject any individual media download larger than this many bytes.",
)
pimage.add_argument(
    "--query-cache-max-age-hours",
    "--query_cache_max_age_hours",
    dest="query_cache_max_age_hours",
    metavar="FLOAT",
    default=168.0,
    type=finite_float,
    required=False,
    action="store",
    help="default=%(default)s: Refresh provider and Bioicons query caches after this age; 0 disables expiration.",
)
pimage.add_argument(
    "--refresh-cache",
    "--refresh_cache",
    dest="refresh_cache",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Ignore existing provider and Bioicons query caches for this run.",
)
pimage.add_argument(
    "--species-name-tsv",
    "--species_name_tsv",
    dest="species_name_tsv",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help='default=%(default)s: Optional legacy-style TSV with "leaf_name" and "species_name" columns. Prefer --species-map-tsv.',
)
pimage.add_argument(
    "--name-tsv",
    "--name_tsv",
    dest="species_name_tsv",
    metavar="PATH",
    type=str,
    help=argparse.SUPPRESS,
)
pimage.add_argument(
    "--manifest-out",
    "--manifest_out",
    dest="manifest_out",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Optional custom output path for manifest.tsv.",
)
pimage.add_argument(
    "--manifest", dest="manifest_out", metavar="PATH", type=str, help=argparse.SUPPRESS
)
pimage.add_argument(
    "--attribution-out",
    "--attribution_out",
    dest="attribution_out",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Optional custom output path for ATTRIBUTION.md.",
)
pimage.add_argument(
    "--attribution",
    dest="attribution_out",
    metavar="PATH",
    type=str,
    help=argparse.SUPPRESS,
)
pimage.add_argument(
    "--fail-on-missing",
    "--fail_on_missing",
    dest="fail_on_missing",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Whether unresolved species should trigger a non-zero exit code.",
)
pimage.add_argument(
    "--output-format",
    "--output_format",
    dest="output_format",
    metavar="original|png|jpg",
    default="original",
    type=str,
    required=False,
    action="store",
    choices=["original", "png", "jpg"],
    help="default=%(default)s: Optional normalized output format. "
    'SVG rasterization requires the optional "nwkit[image]" dependency set.',
)
pimage.add_argument(
    "--max-edge",
    "--max_edge",
    dest="max_edge",
    metavar="INT",
    default=None,
    type=int,
    required=False,
    action="store",
    help="default=%(default)s: Downscale rasterized output so the longest edge is at most this value.",
)
pimage.add_argument(
    "--canvas",
    metavar="none|square",
    default="none",
    type=str,
    required=False,
    action="store",
    choices=["none", "square"],
    help="default=%(default)s: Add padding to a square canvas after trimming and resizing.",
)
pimage.add_argument(
    "--background",
    metavar="white|transparent",
    default="white",
    type=str,
    required=False,
    action="store",
    choices=["white", "transparent"],
    help="default=%(default)s: Background color used when --canvas square is requested.",
)
pimage.add_argument(
    "--trim",
    metavar="off|white|transparent|semantic",
    default="off",
    type=str,
    required=False,
    action="store",
    choices=["off", "white", "transparent", "semantic"],
    help="default=%(default)s: Trim uniform white margins, transparent margins, or an estimated main foreground region before resizing.",
)
pimage.add_argument(
    "--trim-shape",
    "--trim_shape",
    dest="trim_shape",
    metavar="bbox|square",
    default="bbox",
    type=str,
    required=False,
    action="store",
    choices=["bbox", "square"],
    help='default=%(default)s: Shape of the trimmed result. "square" center-crops the trimmed content to a square and may clip it.',
)
pimage.set_defaults(handler=command_image)


def command_intersection(args):
    from nwkit.intersection import intersection_main

    intersection_main(args)


pintersection = subparsers.add_parser(
    "intersection",
    help="Drop non-overlapping leaves/sequences in 2 trees or tree+alignment",
    parents=[p_parent],
)
pintersection.add_argument(
    "-i2",
    "--infile2",
    metavar="PATH",
    default="",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Input newick file 2. The intersected version of this file "
    "will not be generated, so if necessary, replace --infile and --infile2 and run again.",
)
pintersection.add_argument(
    "-f2",
    "--format2",
    metavar="auto|auto-strict|INT",
    default="auto",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: ETE tree format for --infile2.",
)
pintersection.add_argument(
    "-si",
    "--seqin",
    metavar="PATH",
    default="",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Input sequence file.",
)
pintersection.add_argument(
    "-so",
    "--seqout",
    metavar="PATH",
    default="",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Output sequence file.",
)
pintersection.add_argument(
    "-sf",
    "--seqformat",
    metavar="fasta",
    default="fasta",
    type=str,
    required=False,
    action="store",
    choices=["fasta"],
    help="default=%(default)s: Sequence format for --seqin and --seqout.",
)
pintersection.add_argument(
    "--match",
    metavar="complete|prefix|backward",
    default="complete",
    type=str,
    required=False,
    action="store",
    choices=["complete", "prefix", "backward"],
    help="default=%(default)s: Method for ID matching.",
)
pintersection.set_defaults(handler=command_intersection)


def command_label(args):
    from nwkit.label import label_main

    label_main(args)


plabel = subparsers.add_parser(
    "label", help="Add unique node labels", parents=[p_parent]
)
plabel.add_argument(
    "-t",
    "--target",
    metavar="all|root|leaf|intnode",
    default="all",
    type=str,
    required=False,
    action="store",
    choices=["all", "root", "leaf", "intnode"],
    help="default=%(default)s: Nodes to be edited.",
)
plabel.add_argument(
    "--prefix",
    metavar="STR",
    default="n",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Prefix for node labels.",
)
plabel.add_argument(
    "--force",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Whether to overwrite existing node names.",
)
plabel.set_defaults(handler=command_label)


def command_rename(args):
    from nwkit.rename import rename_main

    rename_main(args)


prename = subparsers.add_parser(
    "rename",
    help="Rename nodes using a TSV mapping or regular expression",
    parents=[p_parent],
)
prename.add_argument(
    "--name-tsv",
    "--name_tsv",
    dest="name_tsv",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help='default=%(default)s: TSV containing "old_name" and "new_name" columns.',
)
prename.add_argument(
    "--pattern",
    metavar="REGEX",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Optional regular expression used to rename target nodes.",
)
prename.add_argument(
    "--replacement",
    metavar="STR",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Replacement string used with --pattern.",
)
prename.add_argument(
    "-t",
    "--target",
    metavar="all|root|leaf|intnode",
    default="leaf",
    type=str,
    required=False,
    action="store",
    choices=["all", "root", "leaf", "intnode"],
    help="default=%(default)s: Nodes to be renamed.",
)
prename.add_argument(
    "--require-all-old-names",
    "--require_all_old_names",
    dest="require_all_old_names",
    metavar="yes|no",
    default="yes",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Whether every old_name in --name-tsv must match a target node.",
)
prename.add_argument(
    "--require-match",
    "--require_match",
    dest="require_match",
    metavar="yes|no",
    default="yes",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Whether regex mode requires at least one target node name to match.",
)
prename.add_argument(
    "--check-leaf-uniqueness",
    "--check_leaf_uniqueness",
    dest="check_leaf_uniqueness",
    metavar="yes|no",
    default="yes",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Whether to fail if the renamed tree has duplicated or empty leaf labels.",
)
prename.set_defaults(handler=command_rename)


def command_reconcile(args):
    from nwkit.reconcile import reconcile_main

    reconcile_main(args)


preconcile = subparsers.add_parser(
    "reconcile",
    help="Map a rooted gene tree onto a rooted species tree and annotate events",
    parents=[p_tree_input, p_table_output, p_species],
)
preconcile.add_argument(
    "--species-tree",
    "--species_tree",
    dest="species_tree",
    metavar="PATH",
    default=None,
    type=str,
    required=True,
    action="store",
    help="Rooted, strictly bifurcating species tree whose tip labels match parsed gene-tip species labels.",
)
preconcile.add_argument(
    "--tree-id",
    "--tree_id",
    dest="tree_id",
    metavar="TEXT",
    default="",
    type=str,
    required=False,
    action="store",
    help="default=empty: Stable gene-family/tree identifier. Required for unambiguous multi-tree aggregation.",
)
preconcile.add_argument(
    "--species-tree-format",
    "--species_tree_format",
    dest="species_tree_format",
    metavar="auto|auto-strict|INT",
    default="auto",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: ETE tree format for --species-tree.",
)
preconcile.add_argument(
    "--event-source",
    "--event_source",
    dest="event_source",
    metavar="lca|nhx|species-overlap",
    default="lca",
    type=str,
    required=False,
    action="store",
    choices=["lca", "nhx", "species-overlap"],
    help="default=%(default)s: Infer events by LCA reconciliation, read GeneRax-style NHX S/D/H properties, or use the species-overlap heuristic. Valid S placements are authoritative; invalid NHX annotations are never silently replaced.",
)
preconcile.add_argument(
    "--unmatched",
    metavar="error|warn|ignore",
    default="error",
    type=str,
    required=False,
    action="store",
    choices=["error", "warn", "ignore"],
    help="default=%(default)s: Policy for gene tips whose species label cannot be parsed or is absent from --species-tree.",
)
preconcile.set_defaults(handler=command_reconcile)


def command_regress(args):
    from nwkit.regress import regress_main

    regress_main(args)


pregress = subparsers.add_parser(
    "regress",
    help="Fit phylogeny-aware regressions, including Gaussian PGLS and non-Gaussian PGLMMs",
    parents=[p_audit, p_table_output],
    allow_abbrev=False,
    usage=(
        "%(prog)s --tree TREE --data TABLE --responses TRAITS --predictors TRAITS [options]\n"
        "       %(prog)s (--gene-tree TREE | --gene-tree-ensemble TREES) "
        "--species-tree TREE --expression TABLE --species-traits TABLE --tree-id ID "
        "--responses TRAITS --predictors TRAITS [options]\n"
        "       %(prog)s --response-contrasts TABLE --predictor-contrasts TABLE "
        "--responses TRAITS --predictors TRAITS [options]"
    ),
    description=(
        "Fit conventional tip-level regressions, end-to-end reconciled gene-tree "
        "regressions, or models from precomputed reconciled contrasts. The primary "
        "input paths select exactly one workflow; modeling options never select a "
        "workflow implicitly."
    ),
    epilog=(
        "Input modes: conventional uses --tree and --data; end-to-end reconciled "
        "analysis uses a gene tree or ensemble with --species-tree, --expression, "
        "--species-traits, and --tree-id; precomputed analysis uses "
        "--response-contrasts and --predictor-contrasts."
    ),
)
pregress_common = pregress.add_argument_group("common regression inputs")
pregress_ordinary = pregress.add_argument_group(
    "conventional tip-level regression inputs and evolution"
)
pregress_raw = pregress.add_argument_group("end-to-end reconciled inputs")
pregress_precomputed = pregress.add_argument_group("precomputed contrast inputs")
pregress_tree_parsing = pregress.add_argument_group("shared tree parsing")
pregress_response = pregress.add_argument_group("response specification")
pregress_predictor = pregress.add_argument_group("predictor specification")
pregress_tip_matching = pregress.add_argument_group("tip-table matching")
pregress_response_replicates = pregress.add_argument_group(
    "response replicates and known standard errors"
)
pregress_predictor_replicates = pregress.add_argument_group(
    "predictor replicates and known standard errors"
)
pregress_reconciled = pregress.add_argument_group("reconciled regression model")
pregress_inference = pregress.add_argument_group("inference and computation")
pregress_origin = pregress.add_argument_group(
    "end-to-end categorical-origin diagnostics"
)
pregress_diagnostics = pregress.add_argument_group(
    "reconciled diagnostics and auxiliary outputs"
)

pregress_precomputed.add_argument(
    "--response-contrasts",
    dest="response_contrasts",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help='Precomputed mode: reconciled gene-contrast TSV from "nwkit contrast". Use "-" for STDIN.',
)
pregress_precomputed.add_argument(
    "--predictor-contrasts",
    dest="predictor_contrasts",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="Precomputed mode: species-tree contrast TSV from 'nwkit contrast'.",
)
pregress_common.add_argument(
    "--responses",
    metavar="TRAIT1,TRAIT2,...",
    default=None,
    type=str,
    required=True,
    action="store",
    help="Comma-separated response columns in --data/--expression or trait names in --response-contrasts. One model is fitted per response and applicable tree_id.",
)
pregress_common.add_argument(
    "--predictors",
    metavar="TRAIT1,TRAIT2,...",
    default=None,
    type=str,
    required=True,
    action="store",
    help="Comma-separated predictor columns in --data/--species-traits or trait names in --predictor-contrasts.",
)
pregress_response.add_argument(
    "--categorical-responses",
    dest="categorical_responses",
    metavar="TRAIT1,TRAIT2,...",
    default=None,
    type=str,
    help="Responses explicitly treated as unordered categories; non-numeric responses are also detected automatically. Two levels use a logit model and three or more use a multinomial-logit model.",
)
pregress_response.add_argument(
    "--ordered-responses",
    dest="ordered_responses",
    metavar="TRAIT=LOW|MIDDLE|HIGH,...",
    default=None,
    type=str,
    help="Ordered responses and their complete low-to-high level order; a cumulative-logit phylogenetic mixed model is used.",
)
pregress_response.add_argument(
    "--response-reference",
    dest="response_reference",
    metavar="TRAIT=LEVEL,...",
    default=None,
    type=str,
    help="Reference level for each unordered categorical response.",
)
pregress_response.add_argument(
    "--response-family",
    dest="response_family",
    metavar="TRAIT=FAMILY,...",
    default=None,
    type=str,
    help=(
        "Explicit response likelihood by trait. Tree-structured non-Gaussian "
        "GLMMs use sparse Laplace calculations; fits above 5,000 tips are "
        "attempted with a warning. Dense/"
        "custom covariance fallbacks above 500 tips require "
        "'--allow-large-dense yes'. Sparse multinomial fits above 20,000 "
        "tip-by-level linear predictors are likewise attempted with a warning. "
        "Gaussian fits with general dense covariance accept at most 2,000 "
        "observations; diagonal plus low-rank Gaussian fits remain scalable. "
        "Supported families: {}."
    ).format(", ".join(sorted(RESPONSE_FAMILIES))),
)
pregress_response.add_argument(
    "--response-offset",
    dest="response_offset",
    metavar="TRAIT=COLUMN,...",
    default=None,
    type=str,
    help="Log-offset column for count-response traits (for example log library size).",
)
pregress_response.add_argument(
    "--response-trials",
    dest="response_trials",
    metavar="TRAIT=COLUMN,...",
    default=None,
    type=str,
    help="Positive integer trial-count column required by beta-binomial responses.",
)
pregress_response.add_argument(
    "--response-censor-lower",
    dest="response_censor_lower",
    metavar="TRAIT=COLUMN,...",
    default=None,
    type=str,
    help="Lower censor-bound column for censored-Gaussian responses; missing means no lower bound.",
)
pregress_response.add_argument(
    "--response-censor-upper",
    dest="response_censor_upper",
    metavar="TRAIT=COLUMN,...",
    default=None,
    type=str,
    help="Upper censor-bound column for censored-Gaussian responses; missing means no upper bound.",
)
pregress_response.add_argument(
    "--response-dispersion",
    dest="response_dispersion",
    metavar="TRAIT=FLOAT,...",
    default=None,
    type=str,
    help="Fix positive dispersion/shape/precision/SD parameters; omitted values are estimated.",
)
pregress_response.add_argument(
    "--response-zero-probability",
    dest="response_zero_probability",
    metavar="TRAIT=FLOAT,...",
    default=None,
    type=str,
    help="Fix structural-zero probability in (0,1); omitted zero components are estimated.",
)
pregress_response.add_argument(
    "--coefficient-penalty",
    dest="coefficient_penalty",
    metavar="none|gaussian|student-t",
    default=None,
    choices=["none", "gaussian", "student-t"],
    help="default=student-t: Weak coefficient regularization for non-Gaussian models to stabilize sparse or separated data.",
)
pregress_response.add_argument(
    "--coefficient-prior-sd",
    dest="coefficient_prior_sd",
    metavar="FLOAT",
    default=None,
    type=finite_float,
    help="default=2.5: Positive scale for Gaussian or Student-t coefficient regularization.",
)
pregress_response.add_argument(
    "--multivariate-responses",
    dest="multivariate_responses",
    metavar="yes|no",
    default="no",
    type=strtobool,
    help="default=%(default)s: Jointly fit continuous Gaussian responses and estimate their evolutionary covariance (tree-structured fits use sparse calculations and warn above 5,000 tips or 20,000 tip-trait cells; dense fallback supports 2,000 observed cells).",
)
pregress_response.add_argument(
    "--allow-missing-responses",
    dest="allow_missing_responses",
    metavar="yes|no",
    default="no",
    type=strtobool,
    help="default=%(default)s: Retain partially observed tips in a multivariate Gaussian likelihood.",
)
pregress_predictor.add_argument(
    "--categorical-predictors",
    dest="categorical_predictors",
    metavar="TRAIT1,TRAIT2,...",
    default=None,
    type=str,
    help="Predictors explicitly treated as unordered factors; non-numeric predictors are also detected automatically.",
)
pregress_predictor.add_argument(
    "--ordered-predictors",
    dest="ordered_predictors",
    metavar="TRAIT=LOW|MIDDLE|HIGH,...",
    default=None,
    type=str,
    help="Ordered predictors and their complete low-to-high level order; polynomial factor contrasts are used.",
)
pregress_predictor.add_argument(
    "--predictor-reference",
    dest="predictor_reference",
    metavar="TRAIT=LEVEL,...",
    default=None,
    type=str,
    help="Reference level for each unordered categorical predictor.",
)
pregress_predictor.add_argument(
    "--predictor-factor-coding",
    dest="predictor_factor_coding",
    metavar="treatment|sum",
    default=None,
    choices=["treatment", "sum"],
    help="default=treatment: Coding used for unordered categorical predictors.",
)
pregress_predictor.add_argument(
    "--predictor-categorical-replicate-policy",
    dest="predictor_categorical_replicate_policy",
    metavar="error|latent",
    default=None,
    choices=["error", "latent"],
    help="default=error: Require one predictor state per tip or propagate the empirical category mean with sample-size-scaled moment uncertainty.",
)
pregress_ordinary.add_argument(
    "--tree",
    metavar="PATH",
    default=None,
    type=str,
    help="Rooted species tree whose branch lengths define the tip covariance. Its presence selects conventional tip-level regression mode.",
)
pregress_ordinary.add_argument(
    "--data",
    metavar="PATH",
    default=None,
    type=str,
    help='Tip-level TSV with "leaf_name", --responses, --predictors, and optional replicate columns.',
)
pregress_ordinary.add_argument(
    "--tree-format",
    dest="tree_format",
    metavar="auto|auto-strict|INT",
    default=None,
    type=str,
    help="default=auto: ETE tree format for --tree.",
)
pregress_ordinary.add_argument(
    "--branch-length",
    dest="branch_length",
    metavar="original|unit",
    default=None,
    type=str,
    choices=["original", "unit"],
    help="default=original: Use positive tree branch lengths or unit lengths for conventional regression.",
)
pregress_ordinary.add_argument(
    "--evolution-model",
    dest="evolution_model",
    metavar="MODEL",
    default=None,
    type=str,
    choices=list(EVOLUTION_MODELS),
    help="default=brownian: Evolutionary residual covariance model. Shape parameters are estimated unless fixed below.",
)
pregress_ordinary.add_argument(
    "--evolution-parameter",
    dest="evolution_parameter",
    metavar="FLOAT",
    default=None,
    type=finite_float,
    help="Fix the selected model's shape parameter; otherwise parameterized models estimate it.",
)
pregress_ordinary.add_argument(
    "--evolution-covariance",
    dest="evolution_covariance",
    metavar="PATH",
    default=None,
    type=str,
    help="Wide named covariance TSV required by --evolution-model custom and optionally used in model comparison.",
)
pregress_ordinary.add_argument(
    "--compare-evolution-models",
    dest="compare_evolution_models",
    metavar="MODEL1,MODEL2,...",
    default=None,
    type=str,
    help="Fit the listed evolutionary models by ML and calculate AIC, AICc, BIC, and both information-criterion weights.",
)
pregress_ordinary.add_argument(
    "--model-comparison-out",
    dest="model_comparison_out",
    metavar="PATH",
    default=None,
    type=str,
    help="Model-comparison TSV; required with --compare-evolution-models.",
)
pregress_ordinary.add_argument(
    "--intercept",
    metavar="yes|no",
    default=None,
    type=strtobool,
    help="default=yes: Include an intercept in conventional tip-level regression.",
)
pregress_ordinary.add_argument(
    "--response-sampling-covariance-out",
    dest="response_sampling_covariance_out",
    metavar="PATH",
    default=None,
    type=str,
    help="Optional long-form species-mean response sampling covariance for replicate-aware conventional regression.",
)
pregress_ordinary.add_argument(
    "--response-tip-summary-out",
    dest="response_tip_summary_out",
    metavar="PATH",
    default=None,
    type=str,
    help="Optional per-species response mean, sample size, and uncertainty audit for conventional regression.",
)
pregress_ordinary.add_argument(
    "--predictor-evolution-model",
    dest="predictor_evolution_model",
    metavar="MODEL",
    default=None,
    type=str,
    choices=list(EVOLUTION_MODELS),
    help="default=--evolution-model: Evolutionary covariance model for latent conventional-regression predictors.",
)
pregress_ordinary.add_argument(
    "--predictor-evolution-parameter",
    dest="predictor_evolution_parameter",
    metavar="FLOAT",
    default=None,
    type=finite_float,
    help="Fix the latent predictor model's shape parameter; otherwise parameterized models estimate it.",
)
pregress_ordinary.add_argument(
    "--predictor-branch-length",
    dest="predictor_branch_length",
    metavar="original|unit",
    default=None,
    type=str,
    choices=["original", "unit"],
    help="default=--branch-length: Branch-length mode for latent conventional-regression predictors.",
)
pregress_ordinary.add_argument(
    "--predictor-sampling-covariance-out",
    dest="predictor_sampling_covariance_out",
    metavar="PATH",
    default=None,
    type=str,
    help="Optional long-form species-mean predictor sampling covariance.",
)
pregress_ordinary.add_argument(
    "--predictor-tip-summary-out",
    dest="predictor_tip_summary_out",
    metavar="PATH",
    default=None,
    type=str,
    help="Optional per-species predictor replicate and uncertainty audit.",
)
pregress_raw.add_argument(
    "--gene-tree",
    dest="gene_tree",
    metavar="PATH",
    default=None,
    type=str,
    help="Dated rooted gene tree used for expression contrasts. Its presence selects the end-to-end reconciled workflow.",
)
pregress_raw.add_argument(
    "--gene-tree-ensemble",
    dest="gene_tree_ensemble",
    metavar="PATH",
    default=None,
    type=str,
    help="Multi-Newick posterior/bootstrap gene-tree sample. Fits every tree and combines coefficient uncertainty across trees and reconciliations.",
)
pregress_raw.add_argument(
    "--reconciliation-tree",
    dest="reconciliation_tree",
    metavar="PATH",
    default=None,
    type=str,
    help="Optional annotation-bearing gene tree used only for reconciliation; defaults to --gene-tree and must have the same rooted topology.",
)
pregress_raw.add_argument(
    "--species-tree",
    dest="species_tree",
    metavar="PATH",
    default=None,
    type=str,
    help="Rooted species tree used for reconciliation and predictor contrasts.",
)
pregress_raw.add_argument(
    "--expression",
    metavar="PATH",
    default=None,
    type=str,
    help='Expression TSV with "leaf_name" and the --responses columns.',
)
pregress_raw.add_argument(
    "--species-traits",
    dest="species_traits",
    metavar="PATH",
    default=None,
    type=str,
    help='Species-trait TSV with "leaf_name" and the --predictors columns.',
)
pregress_raw.add_argument(
    "--tree-id",
    dest="tree_id",
    metavar="TEXT",
    default=None,
    type=str,
    help="Required non-empty gene-family identifier recorded throughout the output bundle.",
)
pregress_raw.add_argument(
    "--out-prefix",
    dest="out_prefix",
    metavar="PREFIX",
    default=None,
    type=str,
    help="Write the final model and all inspectable intermediate tables under this prefix instead of --outfile.",
)
pregress_raw.add_argument(
    "--gene-tree-format",
    dest="gene_tree_format",
    metavar="auto|auto-strict|INT",
    default=None,
    type=str,
    help="default=auto: ETE tree format for --gene-tree.",
)
pregress_raw.add_argument(
    "--reconciliation-tree-format",
    dest="reconciliation_tree_format",
    metavar="auto|auto-strict|INT",
    default=None,
    type=str,
    help="default=--gene-tree-format: ETE tree format for --reconciliation-tree.",
)
pregress_raw.add_argument(
    "--species-tree-format",
    dest="species_tree_format",
    metavar="auto|auto-strict|INT",
    default=None,
    type=str,
    help="default=auto: ETE tree format for --species-tree.",
)
pregress_tree_parsing.add_argument(
    "--quoted-node-names",
    dest="quoted_node_names",
    metavar="yes|no",
    default=None,
    type=strtobool,
    help="default=yes: Whether node names are quoted in conventional and end-to-end input trees.",
)
pregress_raw.add_argument(
    "--event-source",
    dest="event_source",
    metavar="lca|nhx|species-overlap",
    default=None,
    type=str,
    choices=["lca", "nhx", "species-overlap"],
    help="default=lca: Reconciliation event source; NHX reads GeneRax S/D/H properties.",
)
pregress_raw.add_argument(
    "--species-parser",
    dest="species_parser",
    metavar="PRESET",
    default=None,
    type=str,
    choices=list(SUPPORTED_SPECIES_PARSERS),
    help="default=legacy: Species parser preset for mapping gene tips.",
)
pregress_raw.add_argument(
    "--species-regex",
    dest="species_regex",
    metavar="REGEX",
    default=None,
    type=str,
    help="default=the legacy species regex: Extraction regex for gene-tip species IDs.",
)
pregress_raw.add_argument(
    "--species-map-tsv",
    dest="species_map_tsv",
    metavar="PATH",
    default=None,
    type=str,
    help='Optional mapping TSV with "leaf_name" and "species_label" columns.',
)
pregress_tip_matching.add_argument(
    "--unmatched",
    metavar="error|warn|ignore",
    default=None,
    type=str,
    choices=["error", "warn", "ignore"],
    help="default=error: Policy for rows and tree tips that do not match in raw or conventional input.",
)
pregress_tip_matching.add_argument(
    "--missing-values",
    dest="missing_values",
    metavar="CSV",
    default=None,
    type=str,
    help="default=the NWKIT missing-value set: Values treated as missing in raw or conventional trait tables.",
)
pregress_raw.add_argument(
    "--gene-branch-length",
    dest="gene_branch_length",
    metavar="original|unit",
    default=None,
    type=str,
    choices=["original", "unit"],
    help="default=original: Branch lengths used for expression PICs.",
)
pregress_raw.add_argument(
    "--species-branch-length",
    dest="species_branch_length",
    metavar="original|unit",
    default=None,
    type=str,
    choices=["original", "unit"],
    help="default=original: Branch lengths used for species-trait PICs.",
)
pregress_raw.add_argument(
    "--gene-evolution-model",
    dest="gene_evolution_model",
    metavar="MODEL",
    default=None,
    type=str,
    choices=list(CONTRAST_EVOLUTION_MODELS),
    help="default=brownian: Evolutionary model used for gene-expression contrasts.",
)
pregress_raw.add_argument(
    "--gene-evolution-parameter",
    dest="gene_evolution_parameter",
    metavar="auto|FLOAT",
    default=None,
    type=auto_or_finite_float,
    help="default=auto: Estimate each response's shape parameter, or fix it to FLOAT.",
)
pregress_raw.add_argument(
    "--species-evolution-model",
    dest="species_evolution_model",
    metavar="MODEL",
    default=None,
    type=str,
    choices=list(CONTRAST_EVOLUTION_MODELS),
    help="default=brownian: Evolutionary model used for species-trait contrasts.",
)
pregress_raw.add_argument(
    "--species-evolution-parameter",
    dest="species_evolution_parameter",
    metavar="auto|FLOAT",
    default=None,
    type=auto_or_finite_float,
    help="default=auto: Estimate each predictor's shape parameter by marginal ML, or fix it to FLOAT.",
)
pregress_response_replicates.add_argument(
    "--response-biological-id",
    dest="response_biological_id",
    metavar="COLUMN",
    default=None,
    type=str,
    help="Response-table column identifying independent biological observations in --data or --expression.",
)
pregress_response_replicates.add_argument(
    "--response-technical-id",
    dest="response_technical_id",
    metavar="COLUMN",
    default=None,
    type=str,
    help="Optional technical-replicate identifier nested within a biological observation.",
)
pregress_response_replicates.add_argument(
    "--response-technical-aggregation",
    dest="response_technical_aggregation",
    metavar="error|mean",
    default=None,
    type=str,
    choices=["error", "mean"],
    help="default=error: Reject or explicitly average technical replicates.",
)
pregress_response_replicates.add_argument(
    "--response-batch",
    dest="response_batch",
    metavar="COLUMN",
    default=None,
    type=str,
    help="Optional categorical response-table batch column fitted as a fixed effect.",
)
pregress_response_replicates.add_argument(
    "--response-within-variance",
    dest="response_within_variance",
    metavar="pooled|leaf|known-se",
    default=None,
    type=str,
    choices=["pooled", "leaf", "known-se"],
    help="default=pooled: Biological replicate variance model or known-SE input.",
)
pregress_response_replicates.add_argument(
    "--response-standard-error-columns",
    dest="response_standard_error_columns",
    metavar="COLUMN1,COLUMN2,...",
    default=None,
    type=str,
    help="Known-SE columns corresponding to --responses.",
)
pregress_response_replicates.add_argument(
    "--response-sample-size-columns",
    dest="response_sample_size_columns",
    metavar="COLUMN1,COLUMN2,...",
    default=None,
    type=str,
    help="Optional known-SE sample-size columns corresponding to --responses.",
)
pregress_predictor_replicates.add_argument(
    "--predictor-biological-id",
    dest="predictor_biological_id",
    metavar="COLUMN",
    default=None,
    type=str,
    help="Predictor-table column identifying independent biological observations.",
)
pregress_predictor_replicates.add_argument(
    "--predictor-technical-id",
    dest="predictor_technical_id",
    metavar="COLUMN",
    default=None,
    type=str,
    help="Optional technical-replicate identifier for predictor observations.",
)
pregress_predictor_replicates.add_argument(
    "--predictor-technical-aggregation",
    dest="predictor_technical_aggregation",
    metavar="error|mean",
    default=None,
    type=str,
    choices=["error", "mean"],
    help="default=error: Reject or explicitly average predictor technical replicates.",
)
pregress_predictor_replicates.add_argument(
    "--predictor-batch",
    dest="predictor_batch",
    metavar="COLUMN",
    default=None,
    type=str,
    help="Optional categorical predictor-table batch column fitted as a fixed effect.",
)
pregress_predictor_replicates.add_argument(
    "--predictor-within-variance",
    dest="predictor_within_variance",
    metavar="pooled|leaf|known-se",
    default=None,
    type=str,
    choices=["pooled", "leaf", "known-se"],
    help="default=pooled: Predictor biological-replicate variance model or known-SE input.",
)
pregress_predictor_replicates.add_argument(
    "--predictor-standard-error-columns",
    dest="predictor_standard_error_columns",
    metavar="COLUMN1,COLUMN2,...",
    default=None,
    type=str,
    help="Known-SE columns corresponding to --predictors.",
)
pregress_predictor_replicates.add_argument(
    "--predictor-sample-size-columns",
    dest="predictor_sample_size_columns",
    metavar="COLUMN1,COLUMN2,...",
    default=None,
    type=str,
    help="Optional known-SE sample-size columns corresponding to --predictors.",
)
pregress_reconciled.add_argument(
    "--event-weighting",
    dest="event_weighting",
    metavar="event|contrast",
    default=None,
    type=str,
    required=False,
    action="store",
    choices=["event", "contrast"],
    help="default=event: Give each species event equal total weight, or give each gene contrast equal weight. Event weighting prevents copy-rich events from dominating.",
)
pregress_reconciled.add_argument(
    "--speciation-coverage",
    dest="speciation_coverage",
    metavar="complete|any",
    default=None,
    type=str,
    required=False,
    action="store",
    choices=["complete", "any"],
    help="default=complete: Require complete daughter-clade sampling or include explicitly reported partial coverage.",
)
pregress_inference.add_argument(
    "--confidence-level",
    dest="confidence_level",
    metavar="FLOAT",
    default=0.95,
    type=unit_interval_float,
    required=False,
    action="store",
    help="default=%(default)s: Two-sided confidence-interval level strictly between zero and one.",
)
pregress_reconciled.add_argument(
    "--reconciled-model",
    dest="reconciled_model",
    metavar="hierarchical|replicate-reml|cluster-hc1",
    default=None,
    type=str,
    required=False,
    action="store",
    choices=["hierarchical", "replicate-reml", "cluster-hc1"],
    help="default=hierarchical: Fit the replicate-aware hierarchical model, omit random effects, or run the earlier cluster-HC1 estimator for sensitivity analysis.",
)
pregress_precomputed.add_argument(
    "--response-sampling-covariance",
    dest="response_sampling_covariance",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="Long-form response-contrast sampling covariance from replicate-aware 'nwkit contrast'. A zero matrix is used when the response has no sampling-variance columns.",
)
pregress_precomputed.add_argument(
    "--predictor-sampling-covariance",
    dest="predictor_sampling_covariance",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="Precomputed mode: long-form predictor-contrast sampling covariance from replicate-aware 'nwkit contrast'.",
)
pregress_inference.add_argument(
    "--inference",
    metavar="wald|parametric-bootstrap|likelihood-ratio|profile-likelihood",
    default="wald",
    type=str,
    required=False,
    action="store",
    choices=[
        "wald",
        "parametric-bootstrap",
        "likelihood-ratio",
        "profile-likelihood",
    ],
    help="default=%(default)s: Wald, family-specific parametric bootstrap, likelihood-ratio, or profile-likelihood inference. Tree-structured bootstrap draws use the sparse backend at large tip counts.",
)
pregress_inference.add_argument(
    "--allow-large-dense",
    dest="allow_large_dense",
    metavar="yes|no",
    default=False,
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Explicitly permit large Gaussian or non-Gaussian fits that require a dense covariance representation; nwkit reports an estimated dense working-memory requirement before attempting them.",
)
pregress_inference.add_argument(
    "--bootstrap-replicates",
    dest="bootstrap_replicates",
    metavar="INT",
    default=1000,
    type=int,
    required=False,
    action="store",
    help="default=%(default)s: Number of simulations for parametric-bootstrap inference.",
)
pregress_inference.add_argument(
    "--seed",
    metavar="INT",
    default=1,
    type=int,
    required=False,
    action="store",
    help="default=%(default)s: Non-negative parametric-bootstrap random seed.",
)
pregress_inference.add_argument(
    "--reml",
    metavar="yes|no",
    default="yes",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Use restricted maximum likelihood for Gaussian variance components; predictor-dependent errors-in-variables covariance is always fitted by ML.",
)
pregress_reconciled.add_argument(
    "--event-random-effect",
    dest="event_random_effect",
    metavar="auto|yes|no",
    default=None,
    type=str,
    required=False,
    action="store",
    choices=["auto", "yes", "no"],
    help="default=auto: Include a shared species-event response when identifiable.",
)
pregress_reconciled.add_argument(
    "--lineage-random-slope",
    dest="lineage_random_slope",
    metavar="auto|yes|no",
    default=None,
    type=str,
    required=False,
    action="store",
    choices=["auto", "yes", "no"],
    help="default=auto: Include partially pooled lineage-specific trait slopes when identifiable.",
)
pregress_reconciled.add_argument(
    "--lineage-inference",
    dest="lineage_inference",
    metavar="none|likelihood-ratio|parametric-bootstrap",
    default=None,
    type=str,
    required=False,
    action="store",
    choices=["none", "likelihood-ratio", "parametric-bootstrap"],
    help="default=none: Test lineage-slope heterogeneity and average-plus-lineage joint effects; bootstrap gives calibrated joint-null P-values.",
)
pregress_reconciled.add_argument(
    "--lineage-leave-one-out",
    dest="lineage_leave_one_out",
    metavar="yes|no",
    default=None,
    type=strtobool,
    required=False,
    action="store",
    help="default=no: Refit after removing each reconciled gene lineage and report coefficient sensitivity.",
)
pregress_origin.add_argument(
    "--categorical-origin-diagnostics",
    dest="categorical_origin_diagnostics",
    metavar="none|stochastic-map",
    default=None,
    type=str,
    required=False,
    action="store",
    choices=["none", "stochastic-map"],
    help="default=none: Estimate categorical predictor gains/losses on species-tree branches with an ER Mk stochastic map.",
)
pregress_origin.add_argument(
    "--origin-map-replicates",
    dest="origin_map_replicates",
    metavar="INT",
    default=None,
    type=int,
    required=False,
    action="store",
    help="default=200: Number of stochastic maps for categorical trait-origin diagnostics.",
)
pregress_origin.add_argument(
    "--origin-map-threads",
    dest="origin_map_threads",
    metavar="INT",
    default=None,
    type=int,
    required=False,
    action="store",
    help="default=1: Parallel workers for categorical stochastic maps.",
)
pregress_origin.add_argument(
    "--origin-min-posterior",
    dest="origin_min_posterior",
    metavar="FLOAT",
    default=None,
    type=float,
    required=False,
    action="store",
    help="default=0.5: Minimum transition posterior frequency used for origin leave-one-out.",
)
pregress_origin.add_argument(
    "--origin-leave-one-out",
    dest="origin_leave_one_out",
    metavar="yes|no",
    default=None,
    type=strtobool,
    required=False,
    action="store",
    help="default=no: Refit after omitting events descended from each credible mapped trait origin.",
)
pregress_diagnostics.add_argument(
    "--random-effects-out",
    dest="random_effects_out",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Optional TSV of species-event modes plus lineage deviations, total slopes, intervals, and reliability.",
)
pregress_diagnostics.add_argument(
    "--sensitivity-out",
    dest="sensitivity_out",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Optional TSV of lineage/origin leave-one-out diagnostics.",
)
pregress_origin.add_argument(
    "--trait-origins-out",
    dest="trait_origins_out",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Optional TSV of categorical predictor transition-origin diagnostics.",
)
pregress.set_defaults(handler=command_regress)


def command_mark(args):
    from nwkit.mark import mark_main

    mark_main(args)


pmark = subparsers.add_parser(
    "mark", help="Add texts to node labels", parents=[p_parent]
)
pmark.add_argument(
    "-p",
    "--pattern",
    metavar="REGEX",
    default=".*",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Regular expression for label search.",
)
pmark.add_argument(
    "-t",
    "--target",
    metavar="mrca|clade|leaf",
    default="clade",
    type=str,
    required=False,
    action="store",
    choices=["mrca", "clade", "leaf"],
    help="default=%(default)s: Nodes to be marked.",
)
pmark.add_argument(
    "--target-only-clade",
    "--target_only_clade",
    dest="target_only_clade",
    metavar="yes|no",
    default="yes",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Mark the label of MRCA/clade whose clade contains only target leaves. "
    "Use with --target mrca/clade",
)
pmark.add_argument(
    "--insert-text",
    "--insert-txt",
    "--insert_txt",
    dest="insert_txt",
    metavar="STR",
    default=None,
    type=str,
    required=True,
    action="store",
    help="default=%(default)s: Label to insert to the target node labels.",
)
pmark.add_argument(
    "--insert-separator",
    "--insert-sep",
    "--insert_sep",
    dest="insert_sep",
    metavar="STR",
    default="",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Separator for --insert-text.",
)
pmark.add_argument(
    "--insert-position",
    "--insert-pos",
    "--insert_pos",
    dest="insert_pos",
    metavar="prefix|suffix",
    default="suffix",
    type=str,
    required=False,
    action="store",
    choices=["prefix", "suffix"],
    help="default=%(default)s: Place to insert --insert-text.",
)
pmark.set_defaults(handler=command_mark)


def command_mcmctree(args):
    from nwkit.mcmctree import mcmctree_main

    mcmctree_main(args)


pmcmctree = subparsers.add_parser(
    "mcmctree",
    help="Prepare MCMCtree calibrations or summarize its posterior node ages",
    parents=[p_tree_input, p_download, p_species],
)
pmcmctree.add_argument(
    "-o",
    "--outfile",
    metavar="PATH",
    default="-",
    type=str,
    required=False,
    action="store",
    help='default=%(default)s: Output MCMCtree calibration tree or posterior dated NHX tree. Use "-" for STDOUT.',
)
pmcmctree.add_argument(
    "--left-species",
    "--left_species",
    dest="left_species",
    metavar="STR",
    default=None,
    type=str,
    required=False,
    help="default=%(default)s: Any species in the left clade. "
    "If you want to set a bound on the node splitting Homo_sapiens and Mus_musculus, "
    "specify one of them (e.g., Homo_sapiens).",
)
pmcmctree.add_argument(
    "--right-species",
    "--right_species",
    dest="right_species",
    metavar="STR",
    default=None,
    type=str,
    required=False,
    help="default=%(default)s: Any species in the right clade deriving from the common ancestor. If "
    "you want to set a bound on the node splitting Homo_sapiens and Mus_musculus, specify the "
    "other one that is not used as the left species (e.g., Mus_musculus).",
)
pmcmctree.add_argument(
    "--lower-bound",
    "--lower_bound",
    dest="lower_bound",
    metavar="FLOAT",
    default=None,
    type=str,
    help="default=%(default)s: Lower bound of the calibration point.",
)
pmcmctree.add_argument(
    "--lower-offset",
    "--lower_offset",
    dest="lower_offset",
    metavar="FLOAT",
    default="0.1",
    type=str,
    help="default=%(default)s: ",
)
pmcmctree.add_argument(
    "--lower-scale",
    "--lower_scale",
    dest="lower_scale",
    metavar="FLOAT",
    default="1",
    type=str,
    help="default=%(default)s: ",
)
pmcmctree.add_argument(
    "--lower-tail-prob",
    "--lower_tail_prob",
    "--lower_tailProb",
    dest="lower_tail_prob",
    metavar="FLOAT",
    default="0.025",
    type=str,
    help="default=%(default)s: Lower tail probability. Use 1e-300 for hard bound. Default=0.025",
)
pmcmctree.add_argument(
    "--upper-bound",
    "--upper_bound",
    dest="upper_bound",
    metavar="FLOAT",
    default=None,
    type=str,
    help="default=%(default)s: Upper bound of the calibration point. A point estimate can be "
    "specified by setting the same age in both lower and upper bounds "
    "(e.g., --lower-bound 5.2 --upper-bound 5.2)",
)
pmcmctree.add_argument(
    "--upper-tail-prob",
    "--upper_tail_prob",
    "--upper_tailProb",
    dest="upper_tail_prob",
    metavar="FLOAT",
    default="0.025",
    type=str,
    help="default=%(default)s: Upper tail probability. Use 1e-300 for hard bound. Default=0.025",
)
pmcmctree.add_argument(
    "--add-header",
    "--add_header",
    dest="add_header",
    metavar="yes|no",
    nargs="?",
    const=True,
    default=False,
    type=strtobool,
    help="default=%(default)s: Add the header required for mcmctree.",
)
pmcmctree.add_argument(
    "--timetree",
    metavar="point|ci|no",
    default="no",
    type=str,
    required=False,
    action="store",
    choices=["point", "ci", "no"],
    help="default=%(default)s: Obtain the divergence time from timetree.org. "
    "For point/ci, tip labels are parsed with --species-parser/--species-regex/--species-map-tsv. "
    "Repeated species labels are allowed only when their tips are monophyletic. "
    "point: point estimate, ci: 95 percent confidence interval as upper and lower bounds. "
    "no: disable the function.",
)
pmcmctree.add_argument(
    "--min-clade-prop",
    "--min_clade_prop",
    dest="min_clade_prop",
    metavar="0.0<=FLOAT<=1.0",
    default=0,
    type=unit_interval_float,
    required=False,
    action="store",
    help="default=%(default)s: Minimum proportion of the clade size to the total number of species. "
    "If the clade proportion is smaller than this value, time constraints are removed.",
)
pmcmctree.add_argument(
    "--higher-rank-search",
    "--higher_rank_search",
    dest="higher_rank_search",
    metavar="yes|no",
    default="yes",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Attempt to obtain timetree data using the taxids for higher taxonomic ranks "
    "if the species-level search failed.",
)
pmcmctree.add_argument(
    "--threads",
    metavar="INT",
    default=1,
    type=int,
    required=False,
    action="store",
    help="default=%(default)s: Number of parallel workers used for TimeTree HTTP requests at the same taxonomic rank.",
)
pmcmctree.add_argument(
    "--posterior",
    metavar="PATH",
    default=None,
    type=str,
    help="default=None: Summarize an MCMCtree mcmc.txt file on the input topology and emit a standard NHX dated tree instead of adding calibrations.",
)
pmcmctree.add_argument(
    "--posterior-point",
    "--posterior_point",
    dest="posterior_point",
    metavar="mean|median",
    default="mean",
    type=str,
    choices=["mean", "median"],
    help="default=%(default)s: Point age used for dated-tree branch lengths with --posterior.",
)
pmcmctree.add_argument(
    "--posterior-ci",
    "--posterior_ci",
    dest="posterior_ci",
    metavar="hpd|equal-tail",
    default="hpd",
    type=str,
    choices=["hpd", "equal-tail"],
    help="default=%(default)s: Marginal credible interval stored on each internal node with --posterior.",
)
pmcmctree.add_argument(
    "--posterior-ci-level",
    "--posterior_ci_level",
    dest="posterior_ci_level",
    metavar="0<FLOAT<1",
    default=0.95,
    type=finite_float,
    help="default=%(default)s: Credible mass for --posterior intervals.",
)
pmcmctree.add_argument(
    "--posterior-burnin",
    "--posterior_burnin",
    dest="posterior_burnin",
    metavar="INT",
    default=0,
    type=int,
    help="default=%(default)s: Additional leading posterior rows to discard.",
)
pmcmctree.add_argument(
    "--posterior-thin",
    "--posterior_thin",
    dest="posterior_thin",
    metavar="INT",
    default=1,
    type=int,
    help="default=%(default)s: Keep every INT-th posterior row.",
)
pmcmctree.set_defaults(handler=command_mcmctree)


def command_monophyly(args):
    from nwkit.monophyly import monophyly_main

    monophyly_main(args)


pmonophyly = subparsers.add_parser(
    "monophyly",
    help="Assess whether species or trait-defined groups are monophyletic",
    parents=[p_tree_input, p_table_output, p_species, p_tip_table_policy],
)
pmonophyly.add_argument(
    "--trait",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help='default=%(default)s: Optional TSV containing a "leaf_name" column and one or more grouping columns.',
)
pmonophyly.add_argument(
    "--group-by",
    "--group_by",
    dest="group_by",
    metavar="STR",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Column name in --trait used to define groups. "
    "When omitted, groups are inferred from species labels parsed from the leaf names.",
)
pmonophyly.add_argument(
    "--unrooted",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Evaluate monophyly in unrooted mode.",
)
pmonophyly.add_argument(
    "--fail-on-non-monophyly",
    "--fail_on_non_monophyly",
    dest="fail_on_non_monophyly",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Exit with a non-zero status if any group is not monophyletic.",
)
pmonophyly.set_defaults(handler=command_monophyly)


def command_nhx2nwk(args):
    from nwkit.nhx2nwk import nhx2nwk_main

    nhx2nwk_main(args)


pnhx2nwk = subparsers.add_parser(
    "nhx2nwk", help="Generate Newick from NHX", parents=[p_parent]
)
pnhx2nwk.add_argument(
    "-p",
    "--node-label",
    "--node_label",
    dest="node_label",
    metavar="B|D|H|S|...",
    default="",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: NHX attribute to use as internal node labels.",
)
pnhx2nwk.set_defaults(handler=command_nhx2nwk)


def command_nwk2table(args):
    from nwkit.nwk2table import nwk2table_main

    nwk2table_main(args)


pnwk2table = subparsers.add_parser(
    "nwk2table",
    help="Convert a Newick tree into a parent-child table",
    parents=[p_tree_input, p_table_output],
)
pnwk2table.add_argument(
    "--age",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Add node ages for ultrametric trees.",
)
pnwk2table.add_argument(
    "--sister",
    metavar="yes|no",
    default="yes",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Add the branch_id of one sister node when available.",
)
pnwk2table.set_defaults(handler=command_nwk2table)


def command_printlabel(args):
    from nwkit.printlabel import printlabel_main

    printlabel_main(args)


pprintlabel = subparsers.add_parser(
    "printlabel",
    help="Search and print node labels",
    parents=[p_tree_input, p_text_output],
)
pprintlabel.add_argument(
    "-p",
    "--pattern",
    metavar="REGEX",
    default=".*",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Regular expression for label search.",
)
pprintlabel.add_argument(
    "-t",
    "--target",
    metavar="all|root|leaf|intnode",
    default="all",
    type=str,
    required=False,
    action="store",
    choices=["all", "root", "leaf", "intnode"],
    help="default=%(default)s: Nodes to be searched.",
)
pprintlabel.add_argument(
    "--sister",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Show labels of the sisters instead of targets.",
)
pprintlabel.set_defaults(handler=command_printlabel)


def command_prune(args):
    from nwkit.prune import prune_main

    prune_main(args)


pprune = subparsers.add_parser(
    "prune", help="Prune leaves", parents=[p_parent, p_preserve_properties]
)
pprune.add_argument(
    "-p",
    "--pattern",
    metavar="REGEX",
    default=".*",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Regular expression for label search.",
)
pprune.add_argument(
    "--invert-match",
    "--invert_match",
    dest="invert_match",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Prune unmatched leaves.",
)
pprune.set_defaults(handler=command_prune)


def command_rescale(args):
    from nwkit.rescale import rescale_main

    rescale_main(args)


prescale = subparsers.add_parser(
    "rescale", help="Rescale branch length with a given factor", parents=[p_parent]
)
prescale.add_argument(
    "-t",
    "--target",
    metavar="all|root|leaf|intnode",
    default="all",
    type=str,
    required=False,
    action="store",
    choices=["all", "root", "leaf", "intnode"],
    help="default=%(default)s: Nodes to be edited.",
)
prescale.add_argument(
    "--factor",
    metavar="FLOAT",
    default=None,
    type=finite_float,
    required=True,
    action="store",
    help="default=%(default)s: Rescaling factor of branch length.",
)
prescale.set_defaults(handler=command_rescale)


def command_root(args):
    from nwkit.root import root_main

    root_main(args)


proot = subparsers.add_parser(
    "root",
    help="Place or transfer the tree root",
    parents=[p_parent, p_download, p_species],
)
proot.add_argument(
    "-i2",
    "--infile2",
    metavar="PATH",
    default="",
    type=str,
    required=False,
    action="store",
    help='default=%(default)s: Input newick file 2. Used when --method "transfer". '
    "Leaf labels should be matched to those in --infile.",
)
proot.add_argument(
    "-f2",
    "--format2",
    metavar="auto|auto-strict|INT",
    default="auto",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: ETE tree format for --infile2.",
)
proot.add_argument(
    "--taxon-mode",
    "--taxon_mode",
    dest="taxon_mode",
    metavar="exact|intersection",
    default="exact",
    type=str,
    required=False,
    action="store",
    choices=["exact", "intersection"],
    help="For root transfer, require identical tips or use an unambiguous root split projected onto shared tips.",
)
proot.add_argument(
    "--method",
    metavar="STR",
    default="midpoint",
    type=str,
    required=False,
    action="store",
    choices=[
        "midpoint",
        "outgroup",
        "transfer",
        "mad",
        "mv",
        "reconciliation",
        "taxonomy",
    ],
    help="default=%(default)s: "
    "midpoint: Midpoint rooting. "
    "outgroup: Outgroup rooting with --outgroup. "
    "transfer: Transfer the root position from --infile2 to --infile. "
    "The two trees should have the same bipartitions at the root node. "
    "mad: Minimal Ancestor Deviation rooting (Tria et al. 2017). "
    "mv: Minimum Variance rooting (Mai et al. 2017). "
    "reconciliation: Root by minimizing the weighted LCA-reconciliation "
    "duplication/loss count over every gene-tree edge. "
    "taxonomy: Root by transferring a taxonomy-derived root split onto --infile.",
)
proot.add_argument(
    "--species-tree",
    "--species_tree",
    dest="species_tree",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="Rooted, strictly bifurcating species tree required by --method reconciliation. Its tip labels must match parsed gene-tip species labels.",
)
proot.add_argument(
    "--species-tree-format",
    "--species_tree_format",
    dest="species_tree_format",
    metavar="auto|auto-strict|INT",
    default="auto",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: ETE tree format for --species-tree.",
)
proot.add_argument(
    "--duplication-cost",
    "--duplication_cost",
    dest="duplication_cost",
    metavar="FLOAT",
    default=1.0,
    type=nonnegative_finite_float,
    required=False,
    action="store",
    help="default=%(default)s: Non-negative duplication weight for --method reconciliation.",
)
proot.add_argument(
    "--loss-cost",
    "--loss_cost",
    dest="loss_cost",
    metavar="FLOAT",
    default=1.0,
    type=nonnegative_finite_float,
    required=False,
    action="store",
    help="default=%(default)s: Non-negative implied-loss weight for --method reconciliation. Both reconciliation weights cannot be zero.",
)
proot.add_argument(
    "--outgroup",
    metavar="STR",
    default="",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: An outgroup label or a comma-separated list of outgroup labels. "
    "For the latter, the clade containing all specified labels are used as an outgroup, "
    "so all labels do not have to be specified for a large clade.",
)
proot.add_argument(
    "--taxonomy-source",
    "--taxonomy_source",
    dest="taxonomy_source",
    metavar="ncbi[,opentree,timetree,...]",
    default="ncbi,opentree,timetree",
    type=str,
    required=False,
    action="store",
    help='default=%(default)s: Comma-separated taxonomy/tree sources used when --method "taxonomy" is selected. '
    "Each source is tried in order until rooting succeeds. When taxids are inferred from leaf labels, "
    "taxonomy queries are controlled by --species-parser/--species-regex/--species-map-tsv. "
    "Repeated species labels are allowed for timetree/opentree only when their tips are monophyletic. "
    "Supported sources: ncbi, timetree, opentree.",
)
proot.add_argument(
    "--taxid-tsv",
    "--taxid_tsv",
    dest="taxid_tsv",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help='default=%(default)s: TSV file containing tree leaf labels in the "leaf_name" column and '
    'their taxonomy IDs in the "taxid" column. Used when --method "taxonomy". '
    "If omitted, taxonomy IDs are inferred from leaf labels. This option is used by the NCBI source.",
)
proot.add_argument(
    "--rank",
    metavar="no|species|genus|family|order|...",
    default="no",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Highest taxonomic rank retained when deriving the taxonomy tree for "
    '--method "taxonomy".',
)
proot.set_defaults(handler=command_root)


def command_rootcompare(args):
    from nwkit.root_compare import root_compare_main

    root_compare_main(args)


prootcompare = subparsers.add_parser(
    "rootcompare",
    help="Compare rooting methods in a TSV table and marked PDF tree",
    parents=[p_tree_input, p_table_output, p_download, p_species],
)
prootcompare.add_argument(
    "--methods",
    metavar="all|METHOD[,METHOD,...]",
    default="all",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Rooting methods to compare. 'all' always runs midpoint, MAD, and MV; it also runs outgroup, transfer, and reconciliation when their required inputs are supplied, and independently runs configured taxonomy sources when species names can be parsed (or NCBI taxids are supplied).",
)
prootcompare.add_argument(
    "--exclude-methods",
    "--exclude_methods",
    dest="exclude_methods",
    metavar="METHOD[,METHOD,...]",
    default=None,
    type=str,
    required=False,
    action="store",
    help="Exclude methods from the resolved method set. 'taxonomy' excludes every taxonomy source.",
)
prootcompare.add_argument(
    "--figure-out",
    "--figure_out",
    dest="figure_out",
    metavar="PATH.pdf",
    default=None,
    type=str,
    required=True,
    action="store",
    help="PDF showing every best rooting position, including ties, as a branch marker.",
)
prootcompare.add_argument(
    "-i2",
    "--infile2",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="Rooted reference tree for the transfer method. With --methods all, supplying this option enables transfer automatically.",
)
prootcompare.add_argument(
    "-f2",
    "--format2",
    metavar="auto|auto-strict|INT",
    default="auto",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: ETE tree format for --infile2.",
)
prootcompare.add_argument(
    "--taxon-mode",
    "--taxon_mode",
    dest="taxon_mode",
    metavar="exact|intersection",
    default="exact",
    type=str,
    required=False,
    action="store",
    choices=["exact", "intersection"],
    help="default=%(default)s: Require identical tips for transfer or project the reference root onto shared tips.",
)
prootcompare.add_argument(
    "--species-tree",
    "--species_tree",
    dest="species_tree",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="Rooted, strictly bifurcating species tree for reconciliation. With --methods all, supplying this option enables reconciliation automatically.",
)
prootcompare.add_argument(
    "--species-tree-format",
    "--species_tree_format",
    dest="species_tree_format",
    metavar="auto|auto-strict|INT",
    default="auto",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: ETE tree format for --species-tree.",
)
prootcompare.add_argument(
    "--duplication-cost",
    "--duplication_cost",
    dest="duplication_cost",
    metavar="FLOAT",
    default=1.0,
    type=nonnegative_finite_float,
    required=False,
    action="store",
    help="default=%(default)s: Non-negative reconciliation duplication weight.",
)
prootcompare.add_argument(
    "--loss-cost",
    "--loss_cost",
    dest="loss_cost",
    metavar="FLOAT",
    default=1.0,
    type=nonnegative_finite_float,
    required=False,
    action="store",
    help="default=%(default)s: Non-negative reconciliation implied-loss weight. Both reconciliation weights cannot be zero.",
)
prootcompare.add_argument(
    "--outgroup",
    metavar="STR",
    default=None,
    type=str,
    required=False,
    action="store",
    help="Outgroup label or comma-separated labels. With --methods all, supplying this option enables outgroup rooting automatically.",
)
prootcompare.add_argument(
    "--taxonomy-source",
    "--taxonomy_source",
    dest="taxonomy_source",
    metavar="ncbi[,opentree,timetree,...]",
    default="ncbi,opentree,timetree",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Taxonomy sources to evaluate independently. Under --methods all they are enabled when every tip yields a species name; NCBI is also enabled by --taxid-tsv.",
)
prootcompare.add_argument(
    "--taxid-tsv",
    "--taxid_tsv",
    dest="taxid_tsv",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help='TSV with "leaf_name" and "taxid" columns for the NCBI taxonomy method.',
)
prootcompare.add_argument(
    "--rank",
    metavar="no|species|genus|family|order|...",
    default="no",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Highest taxonomic rank retained for NCBI rooting.",
)
prootcompare.add_argument(
    "--figure-width",
    "--figure_width",
    dest="figure_width",
    metavar="INCHES",
    default=None,
    type=finite_float,
    required=False,
    action="store",
    help="default=auto: PDF figure width in inches. Auto expands unrooted trees according to tip-label density, from 7.2 up to 200 inches.",
)
prootcompare.add_argument(
    "--figure-height",
    "--figure_height",
    dest="figure_height",
    metavar="INCHES",
    default=None,
    type=finite_float,
    required=False,
    action="store",
    help="default=auto: PDF figure height in inches.",
)
prootcompare.add_argument(
    "--font-size",
    "--font_size",
    dest="font_size",
    metavar="POINTS",
    default=8.0,
    type=finite_float,
    required=False,
    action="store",
    help="default=%(default)s: Tip-label and legend font size in points.",
)
prootcompare.add_argument(
    "--tip-labels",
    "--tip_labels",
    dest="tip_labels",
    metavar="auto|yes|no",
    default="auto",
    type=str,
    required=False,
    choices=["auto", "yes", "no"],
    action="store",
    help="default=%(default)s: Show tip labels. Auto omits them above 200 tips to keep large comparison PDFs readable.",
)
prootcompare.add_argument(
    "--unrooted-method",
    "--unrooted_method",
    dest="unrooted_method",
    metavar="auto|equal-daylight|equal-angle",
    default="auto",
    type=str,
    required=False,
    choices=["auto", "equal-daylight", "equal-angle"],
    action="store",
    help="default=%(default)s: Unrooted layout algorithm. Auto uses equal-daylight through 2,000 displayed nodes and equal-angle for larger trees.",
)
prootcompare.add_argument(
    "--layout-report",
    "--layout_report",
    dest="layout_report",
    metavar="PATH.json",
    default=None,
    type=str,
    required=False,
    action="store",
    help="Optional JSON report containing resolved layout dimensions and collision counts.",
)
prootcompare.set_defaults(handler=command_rootcompare)


def command_sanitize(args):
    from nwkit.sanitize import sanitize_main

    sanitize_main(args)


psanitize = subparsers.add_parser(
    "sanitize",
    help="Eliminate non-standard Newick flavors",
    parents=[p_parent, p_preserve_properties],
)
psanitize.add_argument(
    "--remove-singleton",
    "--remove_singleton",
    dest="remove_singleton",
    metavar="yes|no",
    default="yes",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Remove singleton nodes represented as double brackets.",
)
psanitize.add_argument(
    "--resolve-polytomy",
    "--resolve_polytomy",
    dest="resolve_polytomy",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Resolve multifurcation (polytomy) nodes into dichotomies with zero-length branches.",
)
psanitize.add_argument(
    "--name-quote",
    "--name_quote",
    dest="name_quote",
    metavar="none|single|double",
    default="none",
    type=str,
    required=False,
    action="store",
    choices=["none", "single", "double"],
    help="default=%(default)s: Quotation of node and leaf names."
    "none = no quote, single = ', double = \" ",
)
psanitize.set_defaults(handler=command_sanitize)


def command_shuffle(args):
    from nwkit.shuffle import shuffle_main

    shuffle_main(args)


pshuffle = subparsers.add_parser(
    "shuffle", help="Shuffle branches and/or labels", parents=[p_parent]
)
pshuffle.add_argument(
    "--topology",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Randomize entire tree topology and branch length. "
    "Without --label yes, new topology preserve the leaf label orders and is not completely randomized.",
)
pshuffle.add_argument(
    "--branch-length",
    "--branch_length",
    dest="branch_length",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Shuffle branch length. Automatically activated when --topology yes.",
)
pshuffle.add_argument(
    "--label",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Shuffle leaf labels.",
)
pshuffle.add_argument(
    "--seed",
    metavar="INT",
    default=None,
    type=int,
    required=False,
    action="store",
    help="default=%(default)s: Random seed for reproducible topology, branch-length, and label shuffling.",
)
pshuffle.set_defaults(handler=command_shuffle)


def command_skim(args):
    from nwkit.skim import skim_main

    skim_main(args)


pskim = subparsers.add_parser(
    "skim",
    help="Sample leaves from clades with shared traits",
    parents=[p_parent, p_tip_table_policy],
)
pskim.add_argument(
    "--trait",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help='default=%(default)s: Path to a trait table (TSV) containing leaf names in the "leaf_name" column and one or more trait columns. '
    "If not specified, all leaves are treated as a single group and randomly sampled.",
)
pskim.add_argument(
    "--group-by",
    "--group_by",
    dest="group_by",
    metavar="STR",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Column name in the trait table used to group leaves before sampling. "
    "Leaves in a clade are treated as a group if all non-missing values in the specified column are identical. "
    "If not specified, all leaves are treated as a single group.",
)
pskim.add_argument(
    "--retain-per-clade",
    "--retain_per_clade",
    dest="retain_per_clade",
    metavar="INT",
    default=1,
    type=int,
    required=False,
    action="store",
    help="default=%(default)s: Number of leaves to retain per clade (group). "
    "If a clade (group) contains fewer leaves than this value, all of its leaves are retained.",
)
pskim.add_argument(
    "--prioritize-non-missing",
    "--prioritize_non_missing",
    dest="prioritize_non_missing",
    metavar="yes|no",
    default="yes",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Whether to prioritize leaves with non-missing trait values when sampling.",
)
pskim.add_argument(
    "--filter-by",
    "--filter_by",
    dest="filter_by",
    metavar="STR",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Column name in the trait table used to rank leaves within each group before sampling. "
    "If not specified, leaves are randomly sampled within each group.",
)
pskim.add_argument(
    "--filter-mode",
    "--filter_mode",
    dest="filter_mode",
    metavar="ascending|descending",
    default="ascending",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Sorting order for --filter-by values. "
    "If multiple leaves within a group have the same value, they are randomly sampled.",
)
pskim.add_argument(
    "--only-contrastive-clades",
    "--only_contrastive_clades",
    dest="only_contrastive_clades",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Whether to output a pruned tree retaining only minimal clades with multiple non-missing trait values.",
)
pskim.add_argument(
    "--group-table-prefix",
    "--group_table_prefix",
    dest="group_table_prefix",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Output prefix for group-assignment TSVs named <PATH>.all.tsv and <PATH>.sampled.tsv.",
)
pskim.add_argument(
    "--output-groupfile",
    "--output_groupfile",
    dest="output_groupfile",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help=argparse.SUPPRESS,
)
pskim.add_argument(
    "--seed",
    metavar="INT",
    default=None,
    type=int,
    required=False,
    action="store",
    help="default=%(default)s: Random seed for reproducible sampling and tie-breaking.",
)
pskim.set_defaults(handler=command_skim)


def command_sample(args):
    from nwkit.sample import sample_main

    sample_main(args)


psample = subparsers.add_parser(
    "sample",
    help="Select a representative leaf subset by filters, ranks, and sampling method",
    parents=[p_parent, p_tip_table_policy],
)
psample.add_argument(
    "-n",
    "--n",
    metavar="INT",
    default=1,
    type=int,
    required=False,
    action="store",
    help="default=%(default)s: Number of leaves to select.",
)
psample.add_argument(
    "--method",
    metavar="max-pd|ranked",
    default="max-pd",
    type=str,
    required=False,
    action="store",
    choices=["max-pd", "ranked"],
    help="default=%(default)s: Leaf selection method. "
    "max-pd greedily maximizes phylogenetic diversity; ranked selects the first N leaves after filters and ranking.",
)
psample.add_argument(
    "--trait",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help='default=%(default)s: Optional TSV containing leaf names in the "leaf_name" column and metadata columns used by --filter/--rank. '
    "If not specified, all leaves are sampled from the input tree.",
)
psample.add_argument(
    "--filter",
    metavar="COLUMN:OP:VALUE",
    default=[],
    type=str,
    required=False,
    action="append",
    help="default=%(default)s: Keep leaves whose metadata column satisfies a condition. "
    "Can be specified multiple times. Supported operators: ge, gt, le, lt, eq, ne. "
    "Example: --filter busco_complete_pct:ge:80 --filter num_seq:le:200000",
)
psample.add_argument(
    "--rank",
    metavar="COLUMN:asc|desc",
    default=[],
    type=str,
    required=False,
    action="append",
    help="default=%(default)s: Rank retained leaves before sampling. Can be specified multiple times. "
    "For max-pd, rank order is used as a deterministic tie-breaker. "
    "Example: --rank num_seq:asc --rank busco_complete_pct:desc",
)
psample.add_argument(
    "--allow-fewer",
    "--allow_fewer",
    dest="allow_fewer",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Allow fewer than --n leaves when filters leave too few candidates.",
)
psample.add_argument(
    "--report",
    "--output-table",
    "--output_table",
    dest="report",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Optional TSV report for selected leaves, including selection order and metadata columns.",
)
psample.set_defaults(handler=command_sample)


def command_subtree(args):
    from nwkit.subtree import subtree_main

    subtree_main(args)


psubtree = subparsers.add_parser(
    "subtree", help="Generate a subtree Newick file", parents=[p_parent, p_species]
)
psubtree.add_argument(
    "--left-leaf",
    "--left_leaf",
    dest="left_leaf",
    metavar="STR",
    default=None,
    type=str,
    required=False,
    help="default=%(default)s: Any leaf names in the left clade. For example, "
    "to extract the subtree with the root node splitting Homo_sapiens and Mus_musculus, "
    "specify one of them (e.g., Homo_sapiens).",
)
psubtree.add_argument(
    "--right-leaf",
    "--right_leaf",
    dest="right_leaf",
    metavar="STR",
    default=None,
    type=str,
    required=False,
    help="default=%(default)s: Any leaf names in the right clade. For example, "
    "to extract the subtree with the root node splitting Homo_sapiens and Mus_musculus, "
    "specify the other one that is not used as --left-leaf (e.g., Mus_musculus).",
)
psubtree.add_argument(
    "--leaves",
    metavar="leaf1,leaf2,leaf3,...",
    default=None,
    type=str,
    required=False,
    help="default=%(default)s: Comma-separated list of leaves. "
    "The output subtree has their most-recent common ancestor as the root. "
    "--left-leaf and --right-leaf are ignored if this option is specified. "
    "Single leaf name may be specified in combination with --orthogroup yes.",
)
psubtree.add_argument(
    "--orthogroup",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: The output subtree represents orthogroup(s) that contain all "
    "specified leaves. Species assignment is controlled by "
    "--species-parser/--species-regex/--species-map-tsv.",
)
psubtree.add_argument(
    "--dup-conf-score-threshold",
    "--dup_conf_score_threshold",
    dest="dup_conf_score_threshold",
    metavar="FLOAT",
    default=0,
    type=unit_interval_float,
    required=False,
    action="store",
    help="default=%(default)s: The threshold of duplication-confidence score for orthogroup delimitation. "
    "0 = most stringent, 1 = most relaxed. "
    "For the score, see https://www.ensembl.org/info/genome/compara/homology_types.html",
)
psubtree.set_defaults(handler=command_subtree)


def command_transfer(args):
    from nwkit.transfer import transfer_main

    transfer_main(args)


ptransfer = subparsers.add_parser(
    "transfer",
    help="Transfer values and arbitrary properties between compatible tree clades",
    parents=[p_parent],
)
ptransfer.add_argument(
    "-i2",
    "--infile2",
    metavar="PATH",
    default=None,
    type=str,
    required=True,
    action="store",
    help="Source Newick tree. Topologies may differ; taxon matching is controlled by --taxon-mode.",
)
ptransfer.add_argument(
    "-f2",
    "--format2",
    metavar="auto|auto-strict|INT",
    default="auto",
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: ETE tree format for --infile2.",
)
ptransfer.add_argument(
    "-t",
    "--target",
    metavar="all|root|leaf|intnode",
    default="all",
    type=str,
    required=False,
    action="store",
    choices=["all", "root", "leaf", "intnode"],
    help="default=%(default)s: Nodes to be edited.",
)
ptransfer.add_argument(
    "--name",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: transfer node names.",
)
ptransfer.add_argument(
    "--support",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: transfer support values.",
)
ptransfer.add_argument(
    "--length",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: transfer branch length.",
)
ptransfer.add_argument(
    "--property",
    metavar="KEY",
    default=[],
    type=str,
    required=False,
    action="append",
    help="Transfer an arbitrary NHX property without renaming it. May be repeated or comma-separated.",
)
ptransfer.add_argument(
    "--property-map",
    "--property_map",
    dest="property_map",
    metavar="SOURCE=TARGET",
    default=[],
    type=str,
    required=False,
    action="append",
    help="Transfer and rename an arbitrary NHX property. May be repeated.",
)
ptransfer.add_argument(
    "--fill",
    metavar="STR/NUMERIC",
    default=None,
    type=str,
    required=False,
    action="store",
    help="default=%(default)s: Fill values instead of leaving as is, if no corresponding node is found.",
)
ptransfer.add_argument(
    "--taxon-mode",
    "--taxon_mode",
    dest="taxon_mode",
    metavar="exact|intersection",
    default="exact",
    type=str,
    required=False,
    action="store",
    choices=["exact", "intersection"],
    help="Require identical tips or map only unique projections onto shared tips.",
)
ptransfer.add_argument(
    "--match-basis",
    "--match_basis",
    dest="match_basis",
    metavar="clade|split",
    default="clade",
    type=str,
    required=False,
    action="store",
    choices=["clade", "split"],
    help="Match rooted descendant clades or root-independent canonical edge splits.",
)
ptransfer.add_argument(
    "--root-edge-policy",
    "--root_edge_policy",
    dest="root_edge_policy",
    metavar="TARGET_PROPERTY=POLICY",
    default=[],
    type=str,
    required=False,
    action="append",
    help="Resolve root-edge ambiguity per target property. Policies: auto, skip, equal-only, matching-side, mean, min, max, edge-total. May be repeated.",
)
ptransfer.add_argument(
    "--allow-projected-values",
    "--allow_projected_values",
    dest="allow_projected_values",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Permit support and branch-length transfer through projected matches. Use only when shared-tip equivalence is sufficient.",
)
ptransfer.add_argument(
    "--policy",
    metavar="compatible-only|strict",
    default="compatible-only",
    type=str,
    required=False,
    action="store",
    choices=["compatible-only", "strict"],
    help="Skip incompatible requests, or require exact (not projected) matches for every requested value.",
)
ptransfer.add_argument(
    "--report",
    metavar="PATH",
    default=None,
    type=str,
    required=False,
    action="store",
    help="Optional per-node TSV with matches, values, skipped requests, and ambiguity reasons.",
)
ptransfer.set_defaults(handler=command_transfer)


def command_validate(args):
    from nwkit.validate import validate_main

    validate_main(args)


pvalidate = subparsers.add_parser(
    "validate",
    help="Validate one or more Newick trees and report structural issues",
    parents=[p_tree_input, p_table_output, p_species],
)
pvalidate.add_argument(
    "--check-species",
    "--check_species",
    dest="check_species",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Check whether leaf labels are parseable as species labels and species groups are monophyletic.",
)
pvalidate.add_argument(
    "--require-rooted",
    "--require_rooted",
    dest="require_rooted",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Mark unrooted or unknown-rooting trees as invalid. Declared rooted polytomies are accepted.",
)
pvalidate.add_argument(
    "--require-ultrametric",
    "--require_ultrametric",
    dest="require_ultrametric",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Mark non-ultrametric trees as invalid.",
)
pvalidate.add_argument(
    "--require-same-leaf-set",
    "--require_same_leaf_set",
    dest="require_same_leaf_set",
    metavar="yes|no",
    default="yes",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: For multiple input trees, mark trees whose leaf set differs from the first tree as invalid.",
)
pvalidate.add_argument(
    "--require-same-rooting",
    "--require_same_rooting",
    dest="require_same_rooting",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: For multiple input trees with the same leaf set, mark trees whose rooting status or partition into immediate root-child clades differs from the first parsed tree as invalid.",
)
pvalidate.add_argument(
    "--require-binary",
    "--require_binary",
    dest="require_binary",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Mark trees with singleton or multifurcating internal nodes as invalid.",
)
pvalidate.add_argument(
    "--require-all-support",
    "--require_all_support",
    dest="require_all_support",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Mark trees with missing internal support values as invalid.",
)
pvalidate.add_argument(
    "--require-unambiguous-format",
    "--require_unambiguous_format",
    dest="require_unambiguous_format",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Mark trees with auto-format ambiguity as invalid.",
)
pvalidate.add_argument(
    "--require-unquoted-names",
    "--require_unquoted_names",
    dest="require_unquoted_names",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Mark trees containing quoted node names as invalid.",
)
pvalidate.add_argument(
    "--fail-on-issue",
    "--fail_on_issue",
    dest="fail_on_issue",
    metavar="yes|no",
    default="no",
    type=strtobool,
    required=False,
    action="store",
    help="default=%(default)s: Exit with a non-zero status if any reported issue is found.",
)
pvalidate.set_defaults(handler=command_validate)


def command_table2nwk(args):
    from nwkit.table2nwk import table2nwk_main

    table2nwk_main(args)


ptable2nwk = subparsers.add_parser(
    "table2nwk",
    help="Convert a parent-child table into a Newick tree",
    parents=[p_table_input],
)
ptable2nwk.set_defaults(handler=command_table2nwk)


def command_help(args):
    print(parser.parse_args([args.command, "--help"]))


parser_help = subparsers.add_parser("help", help="Show help messages")
parser_help.add_argument("command", help="command name which help is shown")
parser_help.set_defaults(handler=command_help)


def _add_input_rooting_options():
    from nwkit.rooting_state import ROOTING_INPUT_OPTIONS

    for command, command_parser in subparsers.choices.items():
        if command in {"table2nwk", "help"}:
            continue
        input_names = {action.dest for action in command_parser._actions}
        if command == "regress":
            input_names.discard("infile")  # Precomputed contrasts are a table.
        options: dict[str, list[str]] = {}
        for role, option in ROOTING_INPUT_OPTIONS.items():
            if role in input_names:
                options.setdefault(option, []).append("--" + role.replace("_", "-"))
        for option, inputs in options.items():
            flags = ["--" + option]
            if command != "regress":
                flags.append("--" + option.replace("-", "_"))
            command_parser.add_argument(
                *flags,
                choices=("auto", "yes", "no"),
                default="auto",
                help="default=%(default)s: Rooting interpretation for {}. "
                "auto respects [&R]/[&U] and root NHX, otherwise assumes rooted only "
                "for roots with at most two children; yes/no override without rerooting or "
                "relaxing other validation.".format(", ".join(inputs)),
            )


_add_input_rooting_options()


def _validate_stdin_ownership(args):
    stdin_options = [option for _dest, option in get_stdin_input_options(args)]
    if len(stdin_options) > 1:
        raise ValueError(
            "STDIN can be assigned to only one input option at a time: {}".format(
                ", ".join(stdin_options)
            )
        )


LEGACY_OPTION_CANONICAL_OVERRIDES = {
    "dist": {
        "-d": "--metric",
        "--dist": "--metric",
    },
    "skim": {
        "--output-groupfile": "--group-table-prefix",
        "--output_groupfile": "--group-table-prefix",
    },
}


def _deprecated_option_aliases(command):
    """Return deprecated option aliases and their canonical replacements."""
    command_parser = subparsers.choices.get(command)
    if command_parser is None:
        return dict()
    canonical_by_dest = dict()
    for action in command_parser._actions:
        if action.help == argparse.SUPPRESS:
            continue
        canonical = _canonical_long_option(action)
        if canonical is not None:
            canonical_by_dest[action.dest] = canonical
    aliases = dict()
    for action in command_parser._actions:
        canonical = canonical_by_dest.get(action.dest)
        if canonical is None:
            continue
        for option in action.option_strings:
            if option.startswith("--") and option != canonical:
                aliases[option] = canonical
    aliases.update(LEGACY_OPTION_CANONICAL_OVERRIDES.get(command, {}))
    return aliases


def _used_deprecated_option_aliases(argv, command):
    aliases = _deprecated_option_aliases(command)
    used = list()
    seen = set()
    for token in argv:
        if token == "--":
            break
        option = str(token).split("=", 1)[0]
        if option in aliases and option not in seen:
            used.append((option, aliases[option]))
            seen.add(option)
    return used


def _warn_deprecated_option_aliases(argv, command):
    for alias, canonical in _used_deprecated_option_aliases(argv, command):
        sys.stderr.write(
            "Warning: option '{}' is deprecated; use '{}' instead.\n".format(
                alias,
                canonical,
            )
        )


def main(argv=None):
    is_console_invocation = argv is None
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    args = parser.parse_args(raw_argv)
    if hasattr(args, "handler"):
        try:
            if getattr(args, "audit", None) == "-":
                raise ValueError("'--audit' requires a file path, not '-'.")

            def invoke_handler(parsed_args):
                _warn_deprecated_option_aliases(raw_argv, parsed_args.command)
                _validate_stdin_ownership(parsed_args)
                return parsed_args.handler(parsed_args)

            if getattr(args, "audit", None):
                from nwkit.provenance import run_with_audit

                return run_with_audit(args=args, argv=raw_argv, handler=invoke_handler)
            return invoke_handler(args)
        except Exception as exc:
            if (not is_console_invocation) or getattr(args, "debug", False):
                raise
            sys.stderr.write("nwkit: error: {}\n".format(exc))
            return 2
    parser.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
