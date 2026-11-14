"""CLI registration for disparity-through-time."""


def register_dtt(subparsers, tree_input, table_output, table_policy):
    parser = subparsers.add_parser(
        "dtt",
        parents=[tree_input, table_output, table_policy],
        help="Continuous-trait disparity through time with a fitted BM null and MDI.",
    )
    parser.add_argument("--trait", required=True, help="Tip-keyed numeric trait TSV.")
    parser.add_argument(
        "--columns",
        required=True,
        help="One or more comma-separated traits, analyzed jointly.",
    )
    parser.add_argument(
        "--scale",
        choices=["raw", "standardize"],
        default="raw",
        help="Raw Euclidean geometry (default) or observed-SD standardization per trait.",
    )
    parser.add_argument(
        "--missing",
        choices=["error", "drop"],
        default="error",
        help="Error on incomplete tips (default), or prune incomplete tips jointly.",
    )
    parser.add_argument(
        "--n-sim",
        "--n_sim",
        type=int,
        default=999,
        help="BM simulations: default 999; 0 for observed-only DTT, otherwise 2–10,000.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Worker processes (default 1; range 1–32).",
    )
    parser.add_argument(
        "--seed", type=int, default=1, help="Nonnegative random seed (default 1)."
    )
    parser.add_argument(
        "--ci-level",
        "--ci_level",
        type=float,
        default=0.95,
        help="Pointwise BM simulation envelope level (default .95), not a confidence interval.",
    )
    parser.add_argument(
        "--mdi-range",
        "--mdi_range",
        default="0,1",
        help="Relative-time interval for MDI integration (default 0,1).",
    )
    parser.add_argument(
        "--figure-layout",
        "--figure_layout",
        choices=["auto", "trees", "heatmap"],
        default="heatmap",
        help="Single tree plus heatmap (default; auto is an alias), or individual trees.",
    )
    parser.add_argument(
        "--figure-columns",
        "--figure_columns",
        help="Comma-separated subset of --columns to display, in this order; does not change DTT calculation.",
    )
    parser.add_argument(
        "--figure-scale",
        "--figure_scale",
        choices=["standardize", "raw"],
        default="standardize",
        help="Heatmap colors: per-trait observed z scores (default) or a shared raw scale; independent of --scale.",
    )
    for name, help_text in (
        ("summary-out", "Optional one-row TSV with MDI and BM-null area quantiles."),
        (
            "clades-out",
            "Optional node/clade relative-disparity TSV with original branch IDs.",
        ),
        ("simulations-out", "Optional long TSV of every simulated DTT curve."),
        (
            "model-out",
            "Optional JSON of BM rates, scaling, taxa, integration and uncertainty conventions.",
        ),
        (
            "figure-out",
            "Optional PNG, PDF or SVG of observed trait trees, DTT, BM envelope and MDI distribution.",
        ),
    ):
        parser.add_argument("--" + name, "--" + name.replace("-", "_"), help=help_text)
    parser.set_defaults(handler=_command_dtt)


def _command_dtt(args):
    from nwkit.dtt import dtt_main

    return dtt_main(args)
