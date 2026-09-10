"""CLI registration for calibrated and external IC OU shift inference."""


def register_shift(subparsers, tree_input, table_output):
    parser = subparsers.add_parser(
        "shift",
        parents=[tree_input, table_output],
        help="Search for continuous-trait OU shifts with bootstrap calibration.",
        description="Small-tree calibrated OU selection; legacy IC requires Rscript and kfl1ou >= 3.0.9. pBIC also requires the optimum-coordinate capability check.",
    )
    parser.add_argument(
        "--selection",
        choices=["calibrated", "ic", "native"],
        default="calibrated",
        help="Method (default: calibrated; 4..16 tips, <=2 shifts). Native fits multivariate fixed layouts or runs experimental calibrated search without R; IC is the legacy external method.",
    )
    parser.add_argument(
        "--global-null-gate",
        "--global_null_gate",
        action="store_true",
        help="Native AIC search only: full-search plug-in bootstrap gate for the no-shift null; does not control false branches when shifts exist.",
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
        help="Calibration test level (default: 0.05; plug-in calibration, not an exact error guarantee).",
    )
    parser.add_argument(
        "--trait", required=True, help="TSV with leaf_name and a numeric trait column."
    )
    parser.add_argument(
        "--state-column",
        "--state_column",
        required=True,
        help="Continuous trait column; native mode accepts comma-separated columns and missing coordinates.",
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
        choices=["pBIC", "pBICess", "mBIC", "BIC", "AIC", "AICc"],
        default=None,
        help="Information criterion (native: AIC/AICc/BIC/pBIC; native default: bootstrap; ic default: pBIC). Uncorrected pBIC backends are rejected before fitting.",
    )
    parser.add_argument(
        "--root-model",
        "--root_model",
        choices=["OUfixedRoot", "OUrandomRoot"],
        default="OUfixedRoot",
        help="IC root treatment; calibrated contrasts remove the common root component. ASR OUM defaults to a stationary root.",
    )
    parser.add_argument(
        "--search-strategy",
        "--search_strategy",
        choices=["auto", "lasso", "ensemble", "exhaustive", "native-path"],
        default="auto",
        help="Candidate search; native-path requires native AIC or AICc without convergence.",
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
        help="Include shared-effect candidates; calibrated search enumerates them jointly, IC uses backward convergence (AIC/AICc/BIC/pBIC only).",
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
    _register_native_options(parser)
    parser.set_defaults(handler=_command_shift)


def _register_native_options(parser):
    for name, default in [
        ("candidate-pool", 24),
        ("refit-budget", 48),
        ("screening-budget", 2000),
        ("beam-width", 2),
        ("lasso-iterations", 150),
        ("search-memory-mb", 512),
    ]:
        parser.add_argument(
            "--" + name,
            "--" + name.replace("-", "_"),
            type=int,
            default=default,
            help=f"Native search budget (default: {default}).",
        )
    parser.add_argument(
        "--resume-model",
        "--resume_model",
        help="Reuse a completed native JSON only when input, configuration and implementation fingerprints match.",
    )
    parser.add_argument(
        "--regime-map",
        "--regime_map",
        help="Native fixed layout: complete branch_id/regime TSV, including root 0.",
    )
    parser.add_argument(
        "--alpha",
        help="Native fixed alpha in original branch-time units; one value or comma-separated values per trait (0/inf limits supported).",
    )
    parser.add_argument(
        "--process-tip-variance",
        "--process_tip_variance",
        help="Native fixed process tip variance, in original trait units squared; one value or one per trait.",
    )
    parser.add_argument(
        "--measurement-variance",
        "--measurement_variance",
        help="Native fixed extra observation variance in addition to known SE squared; one value or one per trait.",
    )
    parser.add_argument(
        "--estimate-measurement-error",
        "--estimate_measurement_error",
        action="store_true",
        help="Native mode: estimate an additional observation variance per trait.",
    )
    parser.add_argument(
        "--optimizer-starts",
        "--optimizer_starts",
        type=int,
        default=3,
        help="Native deterministic variance-optimizer starts (default: 3).",
    )
    parser.add_argument(
        "--max-iterations",
        "--max_iterations",
        type=int,
        default=300,
        help="Native iterations per variance-optimizer start (default: 300).",
    )


def _command_shift(args):
    if args.global_null_gate and (
        args.selection != "native" or args.criterion != "AIC" or args.regime_map
    ):
        raise ValueError(
            "--global-null-gate requires native AIC search without --regime-map."
        )
    if args.selection == "native":
        from nwkit.shift_native_output import native_main

        return native_main(args)
    if (
        any(
            getattr(args, name) is not None
            for name in (
                "regime_map",
                "alpha",
                "process_tip_variance",
                "measurement_variance",
            )
        )
        or args.estimate_measurement_error
        or args.optimizer_starts != 3
        or args.max_iterations != 300
        or args.resume_model
        or any(
            getattr(args, name) != value
            for name, value in [
                ("candidate_pool", 24),
                ("refit_budget", 48),
                ("screening_budget", 2000),
                ("beam_width", 2),
                ("lasso_iterations", 150),
                ("search_memory_mb", 512),
            ]
        )
    ):
        raise ValueError("Native fitting options require --selection native.")
    from nwkit.shift import shift_main

    return shift_main(args)
