"""CLI registration for MUL D+L, conditional MSC and experimental locus MC."""


def register_mul_reconcile(subparsers, tree_input, table_output, species_parser):
    parser = subparsers.add_parser(
        "mul-reconcile",
        parents=[tree_input, table_output, species_parser],
        help="Exact MUL-tree D+L, conditional MSC, or experimental locus DL+ILS comparison.",
    )
    parser.add_argument("--species-tree", "--species_tree", required=True)
    parser.add_argument(
        "--score-model",
        "--score_model",
        choices=("dl", "msc", "locus-mc"),
        default="dl",
        help="dl: exact D+L; msc: conditional MSC; locus-mc: experimental DL+ILS null/allopolyploid comparison.",
    )
    parser.add_argument(
        "--locus-model",
        "--locus_model",
        default=None,
        help="locus-mc only: explicit model, selection, simulation budget and finite parameter grid JSON.",
    )
    parser.add_argument(
        "--locus-bootstrap",
        "--locus_bootstrap",
        type=int,
        default=None,
        help="locus-mc only: full-search null Monte Carlo replicates per generating null point (default 0).",
    )
    parser.add_argument(
        "--locus-null-calibration",
        "--locus_null_calibration",
        choices=("plug-in", "grid-supremum"),
        default=None,
        help="locus-mc only: fitted-null plug-in (default) or maximum P-value across all supplied null grid points; replicates per point.",
    )
    parser.add_argument(
        "--locus-calibration-out",
        "--locus_calibration_out",
        default=None,
        help="locus-mc only: per-replicate null search/calibration TSV.",
    )
    parser.add_argument(
        "--species-time-unit",
        "--species_time_unit",
        choices=("generations", "coalescent"),
        default=None,
        help="MSC only: explicit unit of ultrametric species lengths and hybridization age.",
    )
    parser.add_argument(
        "--effective-population-size",
        "--effective_population_size",
        type=float,
        default=None,
        help="MSC generations only: fixed shared diploid/subgenome Ne; scale is 2*Ne generations.",
    )
    parser.add_argument(
        "--hybridization-age",
        "--hybridization_age",
        type=float,
        default=None,
        help="MSC only: direct second-parent attachment age, strictly inside both parental branches.",
    )
    parser.add_argument(
        "--max-coalescent-states",
        "--max_coalescent_states",
        type=int,
        default=None,
        help="MSC only: per-family/candidate work limit (default 100000); fail, never truncate.",
    )
    parser.add_argument(
        "--max-coalescent-assignments",
        "--max_coalescent_assignments",
        type=int,
        default=None,
        help="MSC only: assignments per family (default 10000); fail, never truncate.",
    )
    parser.add_argument(
        "--msc-fit",
        "--msc_fit",
        choices=("fixed", "age", "ne", "joint"),
        default=None,
        help="MSC: fixed parameters (default), or bounded age, shared Ne, or joint estimation.",
    )
    for name, help_text in (
        (
            "hybridization-age-bounds",
            "MSC age/joint fit: explicit lower/upper ages in input units.",
        ),
        (
            "population-size-bounds",
            "MSC ne/joint fit: explicit positive lower/upper diploid Ne; generations only.",
        ),
    ):
        parser.add_argument(
            "--" + name,
            "--" + name.replace("-", "_"),
            nargs=2,
            type=float,
            default=None,
            help=help_text,
        )
    for name, help_text in (
        ("msc-grid-points", "MSC fitting grid points per coordinate (default 5, >=3)."),
        (
            "msc-fit-starts",
            "MSC fitting: optimize from this many best grid starts (default 3).",
        ),
        ("msc-maxiter", "MSC fitting: maximum iterations per optimizer (default 200)."),
        (
            "msc-max-evaluations",
            "MSC fitting: unique likelihood evaluations per candidate (default 5000); fail, never truncate.",
        ),
    ):
        parser.add_argument(
            "--" + name,
            "--" + name.replace("-", "_"),
            type=int,
            default=None,
            help=help_text,
        )
    parser.add_argument(
        "--msc-profile-out",
        "--msc_profile_out",
        default=None,
        help="MSC fitting only: nuisance-refitted profile grid; not confidence intervals.",
    )
    parser.add_argument(
        "--h1",
        default=None,
        help="Polyploid clade/postorder selectors: DL defaults to all; MSC/locus-mc require one fixed non-root clade.",
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
        help="DL: optimal best-hypothesis maps; MSC: per-gene likelihoods; locus-mc: per-family MC probabilities.",
    )
    parser.add_argument(
        "--check-out",
        "--check_out",
        default=None,
        help="DL: mapping/state counts; MSC: assignments/work; locus-mc: MC probabilities. No family filtering.",
    )
    parser.add_argument(
        "--tree-out",
        "--tree_out",
        default=None,
        help="Best DL or identified MSC candidate; unsupported in locus-mc. Ties remain in the primary table.",
    )
    parser.add_argument("--model-out", "--model_out", default=None)
    parser.add_argument(
        "--node-out",
        "--node_out",
        default=None,
        help="DL only: complete node assignments for every tied best hypothesis; topology/clade IDs, not WGD/SSD calls.",
    )
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
