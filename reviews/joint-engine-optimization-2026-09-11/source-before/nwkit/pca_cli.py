"""CLI registration for phylogenetic principal components."""


def register_pca(subparsers, tree_input, table_output, table_policy):
    parser = subparsers.add_parser(
        "pca",
        parents=[tree_input, table_output, table_policy],
        help="Phylogenetic PCA with optional ancestral morphospace visualization.",
    )
    parser.add_argument("--trait", required=True, help="Tip-keyed numeric trait TSV.")
    parser.add_argument(
        "--columns",
        required=True,
        help="Comma-separated continuous trait columns (at least two).",
    )
    parser.add_argument("--model", choices=["BM", "LAMBDA"], default="BM")
    parser.add_argument(
        "--mode",
        choices=["cov", "corr"],
        default="cov",
        help="Evolutionary covariance or correlation PCA (default: cov).",
    )
    parser.add_argument(
        "--lambda-value",
        "--lambda_value",
        type=float,
        help="Fix lambda in [0,1]; requires --model LAMBDA. Otherwise estimate it jointly by ML.",
    )
    parser.add_argument(
        "--missing",
        choices=["error", "drop"],
        default="error",
        help="Error on missing traits (default), or drop incomplete tips jointly.",
    )
    for name, help_text in (
        (
            "loadings-out",
            "Optional long TSV of rotation coefficients, evolutionary loadings, centers and scales.",
        ),
        (
            "eigenvalues-out",
            "Optional component eigenvalues and explained-variance TSV.",
        ),
        (
            "model-out",
            "Optional JSON containing fitted PCA transformation and used/excluded tips.",
        ),
        (
            "ancestral-out",
            "Optional long TSV of conditional ancestral PC means and intervals.",
        ),
        (
            "figure-out",
            "Optional PDF, SVG or PNG phylomorphospace and loadings figure.",
        ),
    ):
        parser.add_argument("--" + name, "--" + name.replace("-", "_"), help=help_text)
    parser.add_argument(
        "--ci-level",
        "--ci_level",
        type=float,
        default=0.95,
        help="Conditional ancestor interval/ellipse level (default: 0.95).",
    )
    parser.add_argument(
        "--figure-components",
        "--figure_components",
        default="1,2",
        help="Two distinct component numbers for the figure (default: 1,2).",
    )
    parser.add_argument(
        "--figure-tip-labels",
        "--figure_tip_labels",
        choices=["yes", "no"],
        default="yes",
    )
    parser.set_defaults(handler=_command_pca)


def _command_pca(args):
    from nwkit.pca import pca_main

    return pca_main(args)
