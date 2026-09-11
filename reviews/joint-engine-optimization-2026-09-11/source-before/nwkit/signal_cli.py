"""CLI registration for continuous-trait phylogenetic signal."""


def register_signal(subparsers, tree_input, table_output, table_policy):
    parser = subparsers.add_parser(
        "signal",
        parents=[tree_input, table_output, table_policy],
        help="Estimate and test Blomberg's K and Pagel's lambda for continuous traits.",
    )
    parser.add_argument("--trait", required=True, help="Tip-keyed trait TSV.")
    parser.add_argument(
        "--columns", required=True, help="Comma-separated continuous trait columns."
    )
    parser.add_argument(
        "--standard-error-column",
        "--standard_error_column",
        help="Comma-separated known SE columns, one per trait.",
    )
    parser.add_argument("--method", choices=["K", "lambda", "both"], default="both")
    parser.add_argument("--test", choices=["yes", "no"], default="yes")
    parser.add_argument(
        "--n-sim",
        "--n_sim",
        type=int,
        default=999,
        help="K permutations (default: 999).",
    )
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument(
        "--ci-level",
        "--ci_level",
        type=float,
        default=0.95,
        help="Lambda profile interval level (default: 0.95).",
    )
    parser.add_argument(
        "--p-adjust",
        "--p_adjust",
        choices=["bh", "none"],
        default="bh",
        help="Adjustment across traits separately for each method (default: bh).",
    )
    parser.set_defaults(handler=_command_signal)


def _command_signal(args):
    from nwkit.signal import signal_main

    return signal_main(args)
