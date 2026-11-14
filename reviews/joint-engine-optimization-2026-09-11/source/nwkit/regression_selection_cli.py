"""Exploratory penalized phylogenetic regression CLI."""

SELECTION_SUFFIXES = (
    "coefficients.tsv",
    "path.tsv",
    "cv.tsv",
    "predictions.tsv",
    "stability.tsv",
    "metadata.json",
)


def selection_output_paths(prefix):
    return {suffix: f"{prefix}.{suffix}" for suffix in SELECTION_SUFFIXES}


def register_regression_selection(subparsers, audit):
    parser = subparsers.add_parser(
        "regress-select",
        parents=[audit],
        help="Select predictors with phylogenetic elastic net and nested group CV (no P-values).",
    )
    parser.add_argument("--tree", required=True, help="Branch-length Newick tree.")
    parser.add_argument("--data", required=True, help="Numeric TSV keyed by leaf_name.")
    parser.add_argument("--response", required=True)
    predictor_group = parser.add_mutually_exclusive_group(required=True)
    predictor_group.add_argument(
        "--predictors", help="Comma-separated numeric predictor columns."
    )
    predictor_group.add_argument(
        "--predictor-file",
        "--predictor_file",
        help="One predictor column name per line, for large matrices.",
    )
    parser.add_argument(
        "--unpenalized",
        default="",
        help="Predictor subset exempt from penalties; intercept is always exempt.",
    )
    parser.add_argument(
        "--family",
        choices=("gaussian", "binomial", "poisson", "negative-binomial"),
        default="gaussian",
    )
    parser.add_argument(
        "--evolution-model",
        "--evolution_model",
        choices=["brownian", "independent"],
        default="brownian",
    )
    parser.add_argument(
        "--folds",
        required=True,
        help="TSV with leaf_name and fold; at least three user-defined phylogenetic groups.",
    )
    parser.add_argument(
        "--strengths",
        default="1,0.1,0.01",
        help="Positive comma-separated penalty strengths.",
    )
    parser.add_argument(
        "--l1-ratios",
        "--l1_ratios",
        default="1,0.5",
        help="Values in (0,1]; 1=lasso, smaller=elastic net.",
    )
    parser.add_argument(
        "--prediction",
        choices=["conditional", "fixed"],
        default="conditional",
        help="Plug-in prediction from training random modes or fixed effects only; neither integrates latent uncertainty.",
    )
    parser.add_argument(
        "--out-prefix",
        "--out_prefix",
        required=True,
        help="Prefix for coefficient/path/CV/prediction/stability TSVs and metadata JSON.",
    )
    parser.set_defaults(handler=_command)


def _command(args):
    from nwkit.regression_selection import selection_main

    return selection_main(args)
