"""CLI registration for native duplication/loss/genome-multiplication counts."""


def register_wgd_count(subparsers, tree_input, table_output, table_policy):
    parser = subparsers.add_parser(
        "wgd-count",
        parents=[tree_input, table_output, table_policy],
        help="Native gene-count DL/WGM candidate scan with search-wide null bootstrap.",
    )
    parser.add_argument(
        "--ascertainment",
        choices=["observed", "root-clades"],
        default="observed",
        help="Observation selection (default observed); root-clades requires copies in every root-child clade.",
    )
    parser.add_argument(
        "--counts",
        required=True,
        help="Wide TSV: family_id and exact species-tree tip columns; missing is not zero.",
    )
    for name, default, kind, help_text in (
        (
            "detection-tsv",
            None,
            str,
            "Optional leaf_name,detection_probability TSV for known observation probabilities.",
        ),
        (
            "rate-groups-tsv",
            None,
            str,
            "Optional branch_id,regime TSV covering all non-root input branches.",
        ),
        (
            "candidate-branches",
            None,
            str,
            "Comma-separated input branch IDs; default scans every positive-length non-root branch.",
        ),
        (
            "event-fractions",
            "0.25,0.5,0.75",
            str,
            "Finite candidate position grid measured from each branch's parent (default .25,.5,.75).",
        ),
        (
            "multiplicity",
            2,
            int,
            "Genome multiplication factor >=2 (default 2); factor is not estimated.",
        ),
        (
            "family-gamma-shape",
            "1",
            str,
            "Fixed family-rate gamma shape (default 1); 'none' for homogeneous families.",
        ),
        (
            "family-rate-categories",
            4,
            int,
            "Equal-weight gamma quantile categories (default 4).",
        ),
        (
            "bootstrap",
            0,
            int,
            "Null bootstrap datasets (default 0: exploratory scan with no p-values).",
        ),
        ("seed", 1, int, "Nonnegative seed for local simulation RNG (default 1)."),
        (
            "alpha",
            0.05,
            float,
            "Conditional bootstrap decision level (default .05); not posterior probability.",
        ),
        (
            "max-states",
            256,
            int,
            "Maximum count bound including the independent doubling check (default 256).",
        ),
        (
            "state-tolerance",
            1e-7,
            float,
            "Maximum per-family log-likelihood change on doubling the state bound.",
        ),
        (
            "max-iterations",
            200,
            int,
            "Iterations per optimization start (default 200).",
        ),
        (
            "model-out",
            None,
            str,
            "Optional JSON containing fit, conditioning, state checks and calibration metadata.",
        ),
    ):
        flags = ["--" + name]
        if "-" in name:
            flags.append("--" + name.replace("-", "_"))
        parser.add_argument(*flags, default=default, type=kind, help=help_text)
    parser.add_argument(
        "--rate-model",
        "--rate_model",
        choices=["homogeneous", "terminal-internal"],
        default="terminal-internal",
        help="Background rate groups (default terminal-internal); explicit rate-group TSV overrides this.",
    )
    parser.set_defaults(handler=_command_wgd_count)


def _command_wgd_count(args):
    from nwkit.wgd_count import wgd_count_main

    wgd_count_main(args)
