"""CLI registration for family-aware focal-lineage Ks correction."""


def register_ksrate(subparsers, tree_input, table_output):
    parser = subparsers.add_parser(
        "ksrate",
        parents=[tree_input, table_output],
        help="Focal-lineage Ks correction with multiple outgroups and family bootstrap.",
    )
    parser.add_argument(
        "--ks-tsv",
        "--ks_tsv",
        required=True,
        help="TSV: species_a,species_b,family_id,ks; one ortholog observation per pair/family.",
    )
    for name, default, kind, help_text in (
        ("focals", None, str, "Comma-separated focal species; default every tree tip."),
        (
            "bootstrap",
            199,
            int,
            "Shared-family bootstrap draws (default 199; 0 disables bootstrap diagnostics).",
        ),
        ("seed", 1, int, "Nonnegative random seed (default 1)."),
        (
            "ci-level",
            0.95,
            float,
            "Interval level (default .95), conditional on input tree/pairs.",
        ),
        ("trios-out", None, str, "Optional trio-level diagnostics TSV."),
        (
            "model-out",
            None,
            str,
            "Optional JSON of estimators, selection, seed and uncertainty meaning.",
        ),
    ):
        flags = ["--" + name]
        if "-" in name:
            flags.append("--" + name.replace("-", "_"))
        parser.add_argument(*flags, default=default, type=kind, help=help_text)
    parser.add_argument(
        "--ci-method",
        "--ci_method",
        choices=["family-bootstrap-percentile", "pair-median-bonferroni"],
        default="family-bootstrap-percentile",
        help="Primary interval: approximate shared-family percentile bootstrap (default), or conservative finite-sample simultaneous pair-median bounds.",
    )
    parser.add_argument(
        "--outgroup-policy",
        "--outgroup_policy",
        choices=["nearest", "all"],
        default="nearest",
        help="Outgroups from the immediate parent clade (default), or all outside the focal/sister clade.",
    )
    parser.set_defaults(handler=_command_ksrate)


def _command_ksrate(args):
    from nwkit.ksrate import ksrate_main

    ksrate_main(args)
