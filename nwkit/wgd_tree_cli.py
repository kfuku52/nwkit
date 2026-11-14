"""Lazy CLI registration for supplied-parameter gene topology inference."""


def register_wgd_tree(subparsers, tree_input, table_output, species_parser):
    parser = subparsers.add_parser(
        "wgd-tree",
        parents=[tree_input, table_output, species_parser],
        help="Fixed-gene-topology native DL/WGD likelihood and conditional node assignments.",
    )
    parser.add_argument("--species-tree", "--species_tree", required=True)
    parser.add_argument(
        "--count-model",
        "--count_model",
        required=True,
        help="Native wgd-count model JSON including full species-tree identity.",
    )
    parser.add_argument(
        "--event-ids",
        "--event_ids",
        default=None,
        help="Comma-separated species_event_id values; default all supplied count candidates.",
    )
    parser.add_argument("--tree-id", "--tree_id", default="")
    parser.add_argument(
        "--ascertainment", choices=("observed", "root-clades"), default="observed"
    )
    parser.add_argument("--max-tips", "--max_tips", type=int, default=128)
    parser.add_argument("--tolerance", type=float, default=1e-6)
    parser.add_argument(
        "--origin-tolerance", "--origin_tolerance", type=float, default=1e-5
    )
    parser.add_argument("--model-out", "--model_out", default=None)
    parser.set_defaults(handler=_command)


def _command(args):
    from nwkit.wgd_tree import wgd_tree_main

    wgd_tree_main(args)
