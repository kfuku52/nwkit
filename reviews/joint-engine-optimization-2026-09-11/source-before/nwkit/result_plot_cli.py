"""Lightweight result-report options shared by inference and saved-result drawing."""


def register_result_plot_options(draw, reconcile, radte, finite_float):
    mode = draw.add_mutually_exclusive_group()
    mode.add_argument(
        "--reconciliation",
        metavar="TSV",
        default=None,
        help="Render saved reconciliation events with --infile and --species-tree; no inference is run.",
    )
    mode.add_argument(
        "--radte-prefix",
        "--radte_prefix",
        default=None,
        metavar="PREFIX",
        help="Render a saved RADTE result bundle with --species-tree; omit --infile. No dating is run.",
    )
    draw.add_argument(
        "--species-tree",
        "--species_tree",
        default=None,
        metavar="PATH",
        help="Species tree for a saved reconciliation or RADTE report.",
    )
    draw.add_argument(
        "--species-tree-format",
        "--species_tree_format",
        default="auto",
        metavar="FORMAT",
        help="default=auto: Format of the species tree used in result reports.",
    )
    draw.set_defaults(figure_width=None)
    for action in draw._actions:
        if action.dest == "figure_width":
            action.help = "default=3.6 for a single tree; auto for reconciliation/RADTE reports. Width in inches."
    for parser in (reconcile, radte):
        parser.add_argument(
            "--figure-out",
            "--figure_out",
            default=None,
            metavar="PATH",
            help="Optional .pdf, .svg, or .png result report, published with the numerical outputs.",
        )
        parser.add_argument(
            "--figure-width",
            "--figure_width",
            default=None,
            type=finite_float,
            metavar="INCHES",
            help="default=auto: Result report width (at least 8 inches).",
        )
        parser.add_argument(
            "--figure-height",
            "--figure_height",
            default=None,
            type=finite_float,
            metavar="INCHES",
            help="default=auto: Result report height; must accommodate all tip labels and diagnostics.",
        )
        parser.add_argument(
            "--font-size",
            "--font_size",
            default=8.0,
            type=finite_float,
            metavar="POINTS",
            help="default=8: Result report label font size.",
        )
        parser.add_argument(
            "--branch-length-unit",
            "--branch_length_unit",
            default="",
            metavar="TEXT",
            help="Display unit for dated result figures, for example Ma. Default: input time units.",
        )


def result_figure_options(args):
    return dict(
        width=getattr(args, "figure_width", None),
        height=getattr(args, "figure_height", None),
        font_size=getattr(args, "font_size", 8.0),
        time_unit=getattr(args, "branch_length_unit", ""),
    )
