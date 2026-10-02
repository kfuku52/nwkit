"""CLI registration for exact MUL-tree reconciliation and hypothesis search."""


def register_mul_reconcile(subparsers, tree_input, table_output, species_parser):
    parser = subparsers.add_parser(
        "mul-reconcile",
        parents=[tree_input, table_output, species_parser],
        help="Exact MUL-tree duplication/loss parsimony and polyploid hypothesis search.",
    )
    parser.add_argument("--species-tree", "--species_tree", required=True)
    parser.add_argument(
        "--h1",
        default=None,
        help="Space-separated polyploid clades or postorder internal node numbers; default all.",
    )
    parser.add_argument(
        "--h2",
        default=None,
        help="Space-separated second-parent clades or node numbers; default all.",
    )
    parser.add_argument("--multree", choices=("yes", "no"), default="no")
    parser.add_argument(
        "--report",
        default=None,
        help="All optimal mappings for the first lowest-scoring hypothesis.",
    )
    parser.add_argument(
        "--check-out",
        "--check_out",
        default=None,
        help="Exact mapping/state counts; no candidate-dependent gene filtering.",
    )
    parser.add_argument(
        "--tree-out",
        "--tree_out",
        default=None,
        help="First minimum-score candidate; all tied candidates remain in the score table.",
    )
    parser.add_argument("--model-out", "--model_out", default=None)
    for name, default in (
        ("cpus", 1),
        ("max-candidates", 10000),
        ("max-state-pairs", 10000000),
        ("max-maps", 100000),
    ):
        options = (
            ("--" + name, "--" + name.replace("-", "_"))
            if "-" in name
            else ("--" + name,)
        )
        parser.add_argument(*options, type=int, default=default)
    parser.set_defaults(handler=_command)


def _command(args):
    from nwkit.mul_reconcile import mul_reconcile_main

    mul_reconcile_main(args)
