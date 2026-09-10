"""CLI registration for calibrated and external IC OU shift inference."""


def register_shift(subparsers, tree_input, table_output):
    parser = subparsers.add_parser(
        "shift",
        parents=[tree_input, table_output],
        help="Search for continuous-trait OU shifts with bootstrap calibration.",
        description="Small-tree calibrated OU selection; legacy IC selection requires Rscript and kfl1ou >= 3.0.9.",
    )
    parser.add_argument(
        "--selection",
        choices=["calibrated", "ic"],
        default="calibrated",
        help="Selection method (default: calibrated; 4..16 tips, <=2 shifts). IC is the legacy external method.",
    )
    parser.add_argument(
        "--calibration-replicates",
        "--calibration_replicates",
        type=int,
        default=199,
        help="Complete-search parametric calibration replicates (default: 199).",
    )
    parser.add_argument(
        "--calibration-level",
        "--calibration_level",
        type=float,
        default=0.05,
        help="Sequential test level (default: 0.05; plug-in calibration, not an exact error guarantee).",
    )
    parser.add_argument(
        "--trait", required=True, help="TSV with leaf_name and a numeric trait column."
    )
    parser.add_argument(
        "--state-column",
        "--state_column",
        required=True,
        help="Continuous trait column (complete observations required).",
    )
    parser.add_argument(
        "--model-out",
        "--model_out",
        required=True,
        help="JSON model, branch mapping and search metadata (file path).",
    )
    parser.add_argument(
        "--fit-out",
        "--fit_out",
        help="Optional full kfl1ou fit as RDS (uses temporary tip tokens documented in model JSON).",
    )
    parser.add_argument(
        "--rscript", default="Rscript", help="Rscript executable name or path."
    )
    parser.add_argument(
        "--max-shifts",
        "--max_shifts",
        type=int,
        default=2,
        help="Maximum number of shifts (default: 2).",
    )
    parser.add_argument(
        "--criterion",
        choices=["pBIC", "pBICess", "mBIC", "BIC", "AICc"],
        default=None,
        help="IC-only score (default with --selection ic: pBIC).",
    )
    parser.add_argument(
        "--root-model",
        "--root_model",
        choices=["OUfixedRoot", "OUrandomRoot"],
        default="OUfixedRoot",
        help="kfl1ou root treatment; ASR OUM instead defaults to a stationary root.",
    )
    parser.add_argument(
        "--search-strategy",
        "--search_strategy",
        choices=["auto", "lasso", "ensemble", "exhaustive"],
        default="auto",
    )
    parser.add_argument(
        "--exhaustive-max-configurations",
        "--exhaustive_max_configurations",
        type=int,
        default=5000,
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1,
        help="Calibration or legacy ensemble seed (default: 1).",
    )
    parser.add_argument(
        "--bootstrap",
        type=int,
        default=0,
        help="Parametric bootstrap refits for selection support (default: 0, disabled).",
    )
    parser.add_argument(
        "--convergence",
        action="store_true",
        help="Include shared-effect candidates; calibrated search enumerates them jointly, IC uses backward convergence (AICc/BIC/pBIC only).",
    )
    parser.add_argument(
        "--bootstrap-seed",
        "--bootstrap_seed",
        type=int,
        default=1,
        help="Parametric bootstrap seed (default: 1).",
    )
    parser.add_argument(
        "--effects-out",
        "--effects_out",
        help="Optional TSV of shift mean and optimum effects.",
    )
    parser.add_argument(
        "--regime-parameters-out",
        "--regime_parameters_out",
        help="Optional TSV of regime optima under the root=baseline convention.",
    )
    parser.add_argument(
        "--tip-summary-out",
        "--tip_summary_out",
        help="Optional TSV of observations, predictions, residuals and optima.",
    )
    parser.add_argument(
        "--standard-error-column",
        "--standard_error_column",
        help="Known SE of each tip observation; finite nonnegative values in --trait. Variances are SE squared, not estimated.",
    )
    parser.set_defaults(handler=_command_shift)


def _command_shift(args):
    from nwkit.shift import shift_main

    return shift_main(args)
