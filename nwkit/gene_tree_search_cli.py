"""Register targeted gene-tree search without eagerly importing the scorer."""


def register_gene_tree_search(subparsers, tree_input, table_output, species_parser):
    parser = subparsers.add_parser(
        "gene-tree-search",
        parents=[tree_input, table_output, species_parser],
        help="Detect coupled tip sets, propose complete-tip regrafts, and optionally fit GeneRax EVAL.",
    )

    def add(option, **kwargs):
        flags = [option]
        if "-" in option[2:]:
            flags.append("--" + option[2:].replace("-", "_"))
        return parser.add_argument(*flags, **kwargs)

    add("--species-tree", required=True)
    add(
        "--evaluation",
        choices=("none", "generax"),
        default="none",
        help="none proposes candidates only; generax compares independently refitted joint likelihoods.",
    )
    add(
        "--alignment",
        help="Same full-tip aligned FASTA (.gz supported) for every topology.",
    )
    add(
        "--subst-model",
        help="Explicit GeneRax substitution model; required for --evaluation generax.",
    )
    add("--rec-model", choices=("UndatedDL", "UndatedDTL"), default="UndatedDL")
    add(
        "--root-policy",
        choices=("keep", "optimize"),
        default="optimize",
        help="GeneRax EVAL: preserve each supplied candidate root, or optimize its root (default).",
    )
    add(
        "--generax-command",
        default="generax",
        help="Executable/launcher command, parsed without a shell, e.g. 'mpiexec -np 4 generax'.",
    )
    add(
        "--workdir",
        help="New directory for retained GeneRax inputs, outputs and log; required for evaluation.",
    )
    add(
        "--tree-out",
        help="Best evaluated topology; without evaluation retains the unchanged input.",
    )
    add(
        "--candidates-out",
        help="TSV with candidate_id and complete-tip Newick for every retained proposal.",
    )
    add(
        "--sets-out",
        help="TSV of automatically detected sets and their structural evidence.",
    )
    add(
        "--report-out",
        help="JSON with input hashes, model, coverage limits, and evaluation provenance.",
    )
    for name, default, help_text in (
        (
            "max-moved-tips",
            8,
            "Computational cap on set size; set size is otherwise inferred automatically.",
        ),
        (
            "max-proposals",
            64,
            "Maximum detected sets, prioritized by diagnostic D+L reduction per removed tip.",
        ),
        (
            "max-set-states",
            20000,
            "Overlap-cover enumeration budget; truncation is reported.",
        ),
        (
            "beam-width",
            16,
            "Partial joint-regraft states retained at each component insertion; no improvement gate.",
        ),
        (
            "max-candidates",
            128,
            "Maximum complete-tip candidate topologies including the unchanged baseline.",
        ),
        (
            "max-evaluations",
            64,
            "Maximum GeneRax evaluations including the independently refitted baseline.",
        ),
        ("seed", 12345, "GeneRax optimization seed."),
        (
            "eval-rounds",
            2,
            "Equal EVAL rounds per topology; retain and warm-start its best joint fit.",
        ),
        (
            "timeout",
            3600,
            "GeneRax subprocess time limit in seconds; failure is fatal and logs are retained.",
        ),
    ):
        add("--" + name, type=int, default=default, help=help_text)
    parser.set_defaults(handler=_main)


def _main(args):
    from nwkit.gene_tree_search import gene_tree_search_main

    return gene_tree_search_main(args)
