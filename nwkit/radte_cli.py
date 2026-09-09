"""Lightweight argument registration (no numerical imports during --help)."""


def register_radte(subparsers, audit_parent, species_parent, finite_float):
    parser = subparsers.add_parser(
        "radte",
        parents=[audit_parent, species_parent],
        help="Date reconciled gene trees with shared speciation ages and native relaxed clocks",
        description="Experimental gene-tree dating with one time parameter per species event. "
        "See RADTE.md for conditional inference, sequence models, and uncertainty limits.",
    )

    def add(name, **kwargs):
        flags = ["--" + name.replace("_", "-")]
        if "_" in name:
            flags.append("--" + name)
        parser.add_argument(*flags, **kwargs)

    add(
        "backend",
        choices=["native", "mcmctree"],
        default="native",
        help="Dating backend (default: native). MCMCTree is an optional external reference with soft priors.",
    )

    for name, help_text in [
        ("gene_tree", "Rooted gene tree with substitution branch lengths."),
        (
            "generax_nhx",
            "GeneRax NHX gene tree; requires S annotations and resolved D/L events.",
        ),
        ("notung_parsable", "Notung parsable reconciliation accompanying --gene-tree."),
        ("reconciliation", "Reusable TSV from nwkit reconcile or nwkit radte."),
        (
            "species_node_bounds_tsv",
            "Optional hard species age intervals: node, age_min, age_max.",
        ),
        (
            "species_node_intervals_tsv",
            "External species age intervals for display only: node or species_event_id, lower, upper, level, kind, source.",
        ),
        ("alignment", "Optional aligned FASTA with exactly the gene-tip names."),
        (
            "gene_tree_ensemble",
            "Optional Newick gene-tree samples for --uncertainty input-ensemble.",
        ),
        (
            "species_tree_ensemble",
            "Optional jointly dated species-tree samples; preserve the full chronogram per sample.",
        ),
        (
            "likelihood_summary",
            "Previously saved native quadratic sequence likelihood (.json).",
        ),
    ]:
        add(name, default=None, help=help_text)
    add(
        "species_tree",
        required=True,
        help="Rooted ultrametric species tree in time units.",
    )
    add(
        "out_prefix",
        required=True,
        help="Filesystem prefix for dated tree, tables and manifest.",
    )
    add(
        "gene_tree_format",
        default="auto",
        help="ETE format for the gene tree (default: auto).",
    )
    add(
        "species_tree_format",
        default="auto",
        help="ETE format for the species tree (default: auto).",
    )
    add(
        "reconcile",
        choices=["lca"],
        help="Calculate reconciliation internally; no intermediate file required.",
    )
    add(
        "max_age",
        type=finite_float,
        help="Required upper bound for duplications above the species root.",
    )
    add(
        "rate_correlation",
        type=finite_float,
        default=0.0,
        help="Log-rate AR(1) correlation per gene-tree edge; 0 is independent (default).",
    )
    add(
        "rate_sd",
        type=finite_float,
        default=None,
        help="Log-rate SD; default estimates it from the gene branches. Zero selects a strict sequence clock.",
    )
    add(
        "substitution_model",
        choices=["jc69", "hky", "gtr", "f81", "poisson", "lg", "lg-f"],
        default=None,
        help="Default detects DNA (GTR) or protein (LG); choose explicitly for ambiguous alphabets.",
    )
    add(
        "kappa",
        type=finite_float,
        default=None,
        help="Fix the HKY exchangeability ratio; default estimates it from the alignment.",
    )
    add(
        "gtr_exchangeabilities",
        default=None,
        help="Fix six GTR exchangeabilities: AC,AG,AT,CG,CT,GT. Default estimates five relative rates.",
    )
    add(
        "gamma_shape",
        type=finite_float,
        default=None,
        help="Fix the site-rate gamma shape; default estimates it when categories > 1.",
    )
    add(
        "gamma_categories",
        type=int,
        default=None,
        help="Gamma categories; 1 means homogeneous sites (default: 4).",
    )
    add(
        "inference",
        choices=["auto", "marginal", "joint-map"],
        default="auto",
        help="Sequence inference: auto marginalizes rates with a validated quadratic likelihood, "
        "otherwise uses exact conditional joint MAP; marginal forbids that fallback.",
    )
    add(
        "likelihood",
        choices=["auto", "exact", "quadratic"],
        default=None,
        help="Sequence likelihood; auto validates a quadratic approximation and uses exact likelihood if needed.",
    )
    add(
        "uncertainty",
        choices=["none", "laplace", "profile", "bootstrap", "input-ensemble"],
        default="none",
        help="Conditional interval method (default: none); diagnostics describe unavailable intervals.",
    )
    add(
        "ensemble_within_uncertainty",
        choices=["none", "profile", "bootstrap"],
        default="none",
        help="Evaluate conditional intervals separately within each input ensemble sample (default: none).",
    )
    add(
        "interval_level",
        type=finite_float,
        default=0.95,
        help="Interval probability level (default: 0.95).",
    )
    add(
        "bootstrap_replicates",
        type=int,
        default=100,
        help="Bootstrap refits when requested (default: 100).",
    )
    add(
        "starts",
        type=int,
        default=3,
        help="Deterministic-seed optimizer starts (default: 3).",
    )
    add(
        "maxiter",
        type=int,
        default=2000,
        help="Maximum optimizer iterations per start (default: 2000).",
    )
    add("seed", type=int, default=1, help="Seed for starts/bootstrap (default: 1).")
    add(
        "mcmctree_bin",
        default=None,
        help="MCMCTree executable (reference backend only; default: PATH).",
    )
    add(
        "mcmctree_likelihood",
        choices=["exact", "approximate"],
        default=None,
        help="PAML sequence likelihood (default: exact); approximate builds out.BV with BASEML once, then uses usedata=2.",
    )
    add(
        "mcmctree_clock",
        choices=[1, 2, 3],
        type=int,
        default=None,
        help="PAML clock: 1 strict, 2 independent, 3 correlated (default: 2).",
    )
    for name, default_value in [
        ("burnin", 2000),
        ("sampfreq", 10),
        ("samples", 20000),
        ("chains", 2),
    ]:
        add(
            "mcmctree_" + name,
            type=int,
            default=None,
            help=f"MCMCTree {name} (reference backend only; default: {default_value}).",
        )
    add(
        "mcmctree_timeout",
        type=finite_float,
        default=None,
        help="Timeout seconds per MCMCTree chain (default: 1800).",
    )
    add(
        "mcmctree_rate_prior",
        default=None,
        help="PAML rgene_gamma triple in normalized time units (default: '2 20 1').",
    )
    add(
        "mcmctree_variance_prior",
        default=None,
        help="PAML sigma2_gamma triple (default: '1 10 1').",
    )
    parser.set_defaults(handler=command_radte)
    return parser


def command_radte(args):
    from nwkit.radte import radte_main

    return radte_main(args)


def register_radte_compare(subparsers, audit_parent):
    parser = subparsers.add_parser(
        "radte-compare",
        parents=[audit_parent],
        help="Compare saved fixed, bounded and species-ensemble RADTE results",
    )
    for mode in ("fixed", "bounded", "ensemble"):
        parser.add_argument(
            "--" + mode + "-prefix",
            "--" + mode + "_prefix",
            required=True,
            help=f"Saved {mode} RADTE output prefix.",
        )
    parser.add_argument(
        "--species-tree",
        "--species_tree",
        required=True,
        help="Common reference species chronogram.",
    )
    parser.add_argument(
        "--out-prefix",
        "--out_prefix",
        required=True,
        help="Output prefix for PDF, comparison TSV, components TSV and manifest.",
    )
    parser.set_defaults(handler=command_radte_compare)


def command_radte_compare(args):
    from nwkit.radte_compare import compare_main

    return compare_main(args)
