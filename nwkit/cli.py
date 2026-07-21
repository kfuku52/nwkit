import argparse
import sys

from nwkit import __version__
from nwkit.conventions import DEFAULT_TABLE_MISSING_VALUES_CSV, get_stdin_input_options
from nwkit.species_parser import (
    DEFAULT_SPECIES_PARSER,
    DEFAULT_SPECIES_REGEX,
    SUPPORTED_SPECIES_PARSERS,
)


def _canonical_long_option(action):
    return next(
        (
            option
            for option in action.option_strings
            if option.startswith('--') and '_' not in option
        ),
        None,
    )


class CanonicalHelpFormatter(argparse.HelpFormatter):
    """Display canonical options while keeping compatibility aliases parseable."""

    def _format_action_invocation(self, action):
        if not action.option_strings:
            return super()._format_action_invocation(action)
        short_options = [
            option for option in action.option_strings if not option.startswith('--')
        ]
        canonical_long_option = _canonical_long_option(action)
        option_strings = short_options
        if canonical_long_option is not None:
            option_strings.append(canonical_long_option)
        if not option_strings:
            option_strings = list(action.option_strings)
        if action.nargs == 0:
            return ', '.join(option_strings)
        default = self._get_default_metavar_for_optional(action)
        args_string = self._format_args(action, default)
        return ', '.join(
            '{} {}'.format(option, args_string)
            for option in option_strings
        )


class NwkitArgumentParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('formatter_class', CanonicalHelpFormatter)
        super().__init__(*args, **kwargs)


def strtobool(val):
    val = val.lower()
    if val in ('y', 'yes', 't', 'true', 'on', '1'):
        return True
    elif val in ('n', 'no', 'f', 'false', 'off', '0'):
        return False
    else:
        raise ValueError(f"Invalid truth value: {val}")

# Main parser
parser = NwkitArgumentParser(
    prog='nwkit',
    description='A toolkit for Newick trees. See `nwkit SUBCOMMAND -h` for usage (e.g., nwkit constrain -h)',
)
parser.add_argument('--version', action='version', version='%(prog)s {}'.format(__version__))
subparsers = parser.add_subparsers(dest='command')

# Parent parser for shared options
p_audit = argparse.ArgumentParser(add_help=False)
p_audit.add_argument('--audit', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='Append a JSONL provenance record containing arguments, hashes, warnings, and runtime metadata.')

p_parent = argparse.ArgumentParser(add_help=False, parents=[p_audit])
p_parent.add_argument('-i', '--infile', metavar='PATH', default='-', type=str, required=False, action='store',
                      help='default=%(default)s: Input newick file. Use "-" for STDIN.')
p_parent.add_argument('-o', '--outfile', metavar='PATH', default='-', type=str, required=False, action='store',
                      help='default=%(default)s: Output newick file. Use "-" for STDOUT.')
p_parent.add_argument('-f', '--format', metavar='auto|auto-strict|INT', default='auto', type=str, required=False, action='store',
                      help='default=%(default)s: ETE tree format. '
                           '"auto-strict" fails on ambiguous unquoted numeric internal labels. '
                           'See https://etetoolkit.github.io/ete/tutorial/tutorial_trees.html')
p_parent.add_argument('-of', '--outformat', metavar='auto|INT', default='auto', type=str, required=False, action='store',
                      help='ETE tree format for --outfile. "auto" indicates the same format as --format.')
p_parent.add_argument('--quoted-node-names', '--quoted_node_names', dest='quoted_node_names', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                        help='default=%(default)s: Whether node names are quoted in the input file.')

p_download = argparse.ArgumentParser(add_help=False)
p_download.add_argument('--download-dir', '--download_dir', dest='download_dir', metavar='PATH', default='auto', type=str, required=False, action='store',
                        help='default=%(default)s: Shared download/cache directory for external resources such as the ETE4 '
                             'NCBI taxonomy database. "auto" uses the ETE4 default cache location. '
                             '"inferred" stores downloads under <outfile_dir>/downloads, or ./downloads when writing to STDOUT.')

p_tree_input = argparse.ArgumentParser(add_help=False, parents=[p_audit])
p_tree_input.add_argument('-i', '--infile', metavar='PATH', default='-', type=str, required=False, action='store',
                          help='default=%(default)s: Input newick file. Use "-" for STDIN.')
p_tree_input.add_argument('-f', '--format', metavar='auto|auto-strict|INT', default='auto', type=str, required=False, action='store',
                          help='default=%(default)s: ETE tree format. '
                               '"auto-strict" fails on ambiguous unquoted numeric internal labels. '
                               'See https://etetoolkit.github.io/ete/tutorial/tutorial_trees.html')
p_tree_input.add_argument('--quoted-node-names', '--quoted_node_names', dest='quoted_node_names', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                          help='default=%(default)s: Whether node names are quoted in the input file.')

p_table_output = argparse.ArgumentParser(add_help=False)
p_table_output.add_argument('-o', '--outfile', metavar='PATH', default='-', type=str, required=False, action='store',
                            help='default=%(default)s: Output table file. Use "-" for STDOUT.')

p_text_output = argparse.ArgumentParser(add_help=False)
p_text_output.add_argument('-o', '--outfile', metavar='PATH', default='-', type=str, required=False, action='store',
                           help='default=%(default)s: Output text file. Use "-" for STDOUT.')

p_table_input = argparse.ArgumentParser(add_help=False, parents=[p_audit])
p_table_input.add_argument('-i', '--infile', metavar='PATH', default='-', type=str, required=False, action='store',
                           help='default=%(default)s: Input table file. Use "-" for STDIN.')
p_table_input.add_argument('-o', '--outfile', metavar='PATH', default='-', type=str, required=False, action='store',
                           help='default=%(default)s: Output newick file. Use "-" for STDOUT.')
p_table_input.add_argument('-of', '--outformat', metavar='auto|INT', default='auto', type=str, required=False, action='store',
                           help='ETE tree format for --outfile. "auto" infers a suitable format from the input table.')

p_species = argparse.ArgumentParser(add_help=False)
p_species.add_argument('--species-parser', '--species_parser', dest='species_parser', metavar='PRESET',
                       default=DEFAULT_SPECIES_PARSER, required=False, type=str,
                       choices=list(SUPPORTED_SPECIES_PARSERS),
                       help='default=%(default)s: Species parser preset. '
                            'legacy keeps the historical GENUS_SPECIES[_...] behavior. '
                            'taxonomic keeps qualifier-aware species labels such as GENUS_SPECIES_cf or GENUS_sp_LABEL '
                            'for internal identity while simplifying taxonomy queries.')
p_species.add_argument('--species-regex', '--species_regex', dest='species_regex', metavar='REGEX',
                       default=DEFAULT_SPECIES_REGEX, required=False, type=str,
                       help='default=%(default)s: Regular expression for extracting species IDs when --species-parser legacy is used. '
                            'If the regex contains capture groups, non-empty captured groups are joined with underscores; '
                            'otherwise, the full regex match is used. '
                            'For GENUS_SPECIES and GENUS_SPECIES_GENEID labels, the default works as-is.')
p_species.add_argument('--species-map-tsv', '--species_map_tsv', dest='species_map_tsv', metavar='PATH',
                       default=None, required=False, type=str,
                       help='default=%(default)s: Optional TSV overriding parsed species labels and/or taxonomy queries. '
                            'The file must contain a "leaf_name" column and at least one of "species_label" or "taxonomy_query". '
                            'Mapped rows override the selected parser preset and regex.')

p_tip_table_policy = argparse.ArgumentParser(add_help=False)
p_tip_table_policy.add_argument('--missing-values', '--missing_values', dest='missing_values', metavar='CSV',
                                default=DEFAULT_TABLE_MISSING_VALUES_CSV, type=str, required=False, action='store',
                                help='default=%(default)s: Comma-separated table values treated as missing.')
p_tip_table_policy.add_argument('--unmatched', metavar='warn|error|ignore', default='warn', type=str, required=False,
                                action='store', choices=['warn', 'error', 'ignore'],
                                help='default=%(default)s: Policy for table rows or tree tips without a counterpart.')

def command_annotate(args):
    from nwkit.annotate import annotate_main
    annotate_main(args)
pannotate = subparsers.add_parser('annotate', help='Attach tip-table values and aggregate them as Newick properties', parents=[p_parent, p_tip_table_policy])
pannotate.add_argument('--table', metavar='PATH', default=None, type=str, required=True, action='store',
                       help='TSV containing a leaf_name column and one or more annotation columns.')
pannotate.add_argument('--columns', metavar='COL1,COL2,...', default=None, type=str, required=False, action='store',
                       help='Columns attached to matching tips. By default all non-key columns are attached.')
pannotate.add_argument('--property-map', '--property_map', dest='property_map', metavar='COLUMN=PROPERTY', default=[], type=str, required=False, action='append',
                       help='Attach a table column under a different Newick property name. May be repeated.')
pannotate.add_argument('--aggregate', metavar='COLUMN:METHOD[:PROPERTY]', default=[], type=str, required=False, action='append',
                       help='Aggregate descendant-tip values onto internal nodes. Methods: unique, mode, count, mean, sum, min, max, list. May be repeated.')
pannotate.add_argument('--report', metavar='PATH', default=None, type=str, required=False, action='store',
                       help='Optional TSV audit of attached, aggregated, missing, and unmatched values.')
pannotate.set_defaults(handler=command_annotate)

def command_asr(args):
    from nwkit.asr import asr_main
    asr_main(args)
pasr = subparsers.add_parser('asr', help='Infer categorical ancestral states and impute missing tip states under an Mk model', parents=[p_tree_input, p_table_output, p_tip_table_policy])
pasr.add_argument('--trait', metavar='PATH', default=None, type=str, required=True, action='store',
                  help='TSV file containing a "leaf_name" column and a categorical state column.')
pasr.add_argument('--state-column', '--state_column', dest='state_column', metavar='STR', default=None, type=str, required=True, action='store',
                  help='Column name in --trait containing categorical states.')
pasr.add_argument('--states', metavar='STATE1,STATE2,...', default=None, type=str, required=False, action='store',
                  help='default=%(default)s: Optional comma-separated state order. Unlisted observed states are rejected.')
pasr.add_argument('--model', metavar='ER|SYM|ARD', default='ER', type=str, required=False, action='store', choices=['ER', 'SYM', 'ARD'],
                  help='default=%(default)s: Mk rate model. ER uses one shared off-diagonal rate, '
                       'SYM uses symmetric pairwise rates, and ARD uses all rates different.')
pasr.add_argument('--rate', metavar='FLOAT', default=None, type=float, required=False, action='store',
                  help='default=%(default)s: Fixed ER off-diagonal transition rate. For SYM/ARD, this value is used as the initial rate for ML optimization.')
pasr.add_argument('--rate-bounds', '--rate_bounds', dest='rate_bounds', metavar='MIN,MAX', default=None, type=str, required=False, action='store',
                  help='default=1e-9,1e3: Positive bounds used when estimating Mk rates.')
pasr.add_argument('--root-prior', '--root_prior', dest='root_prior', metavar='equal|empirical', default='equal', type=str, required=False, action='store',
                  choices=['equal', 'empirical'],
                  help='default=%(default)s: Prior state frequencies at the root.')
pasr.add_argument('--ambiguous-separator', '--ambiguous_separator', dest='ambiguous_separator', metavar='STR', default='|', type=str, required=False, action='store',
                  help='default=%(default)s: Separator for ambiguous or polymorphic states such as "A|B".')
pasr.add_argument('--target', metavar='all|intnode,missing-leaf|leaf', default='all', type=str, required=False, action='store',
                  help='default=%(default)s: Comma-separated node classes to report. '
                       'Use "missing-leaf" to report only imputed leaves. '
                       'Legacy "tip" and "missing_tip" values remain accepted.')
pasr.add_argument('--output', metavar='probabilities|map', default='probabilities', type=str, required=False, action='store',
                  choices=['probabilities', 'map'],
                  help='default=%(default)s: Report posterior probabilities or only the maximum a posteriori state.')
pasr.add_argument('--model-out', '--model_out', dest='model_out', metavar='PATH', default=None, type=str, required=False, action='store',
                  help='default=%(default)s: Optional TSV reporting fitted Mk model metadata, rates, and log-likelihood.')
pasr.add_argument('--tree-out', '--tree_out', dest='tree_out', metavar='PATH', default=None, type=str, required=False, action='store',
                  help='default=%(default)s: Optional Newick/NHX tree annotated with ASR state and probability properties.')
pasr.add_argument('--tree-outformat', '--tree_outformat', dest='tree_outformat', metavar='auto|INT', default='auto', type=str, required=False, action='store',
                  help='default=%(default)s: ETE tree format for --tree-out.')
pasr.add_argument('--tree-annotation', '--tree_annotation', dest='tree_annotation', metavar='state|probability|map|all', default='map', type=str, required=False, action='store',
                  choices=['state', 'probability', 'map', 'all'],
                  help='default=%(default)s: NHX properties written to --tree-out. '
                       'map writes ASR state and MAP probability; "all" also writes per-state posterior probabilities.')
pasr.add_argument('--stochastic-map-out', '--stochastic_map_out', dest='stochastic_map_out', metavar='PATH', default=None, type=str, required=False, action='store',
                  help='default=%(default)s: Optional TSV summarizing sampled stochastic-map transition counts per branch.')
pasr.add_argument('--n-sim', '--n_sim', dest='n_sim', metavar='INT', default=100, type=int, required=False, action='store',
                  help='default=%(default)s: Number of stochastic maps sampled when --stochastic-map-out is specified.')
pasr.add_argument('--threads', metavar='INT', default=1, type=int, required=False, action='store',
                  help='default=%(default)s: Number of parallel workers used for stochastic mapping simulations.')
pasr.add_argument('--seed', metavar='INT', default=None, type=int, required=False, action='store',
                  help='default=%(default)s: Random seed for stochastic mapping.')
pasr.set_defaults(handler=command_asr)

def command_constrain(args):
    from nwkit.constrain import constrain_main
    constrain_main(args)
pconstrain = subparsers.add_parser('constrain', help='Generate a species-tree-like Newick file for topological constraint', parents=[p_parent, p_download, p_species])
pconstrain.add_argument('--species-list', '--species_list', dest='species_list', metavar='PATH', default=None, type=str, required=False, action='store',
                        help='default=%(default)s: Text file containing species names, one per line. '
                             'Expected formats are "GENUS SPECIES", "GENUS_SPECIES", or "GENUS_SPECIES_OTHERINFO". '
                             'Species-aware parsing uses the genus and species fields only, ignoring any remaining suffix. '
                             'e.g., "Arabidopsis thaliana" and "Arabidopsis_thaliana_TAIR10"')
pconstrain.add_argument('--taxid-tsv', '--taxid_tsv', dest='taxid_tsv', metavar='PATH', default=None, type=str, required=False, action='store',
                        help='default=%(default)s: TSV file containing species names in the "leaf_name" column and their NCBI Taxonomy IDs in the "taxid" column. '
                             'When specified, the provided NCBI Taxonomy IDs are used instead of inferring them from species names. '
                             'Either --species-list or --taxid-tsv must be specified, but not both. '
                             'This option is currently compatible only with --backbone ncbi.')
pconstrain.add_argument('--backbone', metavar='ncbi|ncbi_apgiv|ncbi_user|user', default='ncbi', type=str, required=False,
                        action='store', choices=['ncbi', 'ncbi_apgiv', 'ncbi_user', 'user'],
                        help='default=%(default)s: The backbone for tree constraint. '
                             '--infile is not required except for "user". '
                             'ncbi: Infer NCBI Taxonomy ID from species name, and generate a tree based on the ranks. '
                             'ncbi_apgiv: Infer NCBI Taxonomy ID from species name, and match it with the order-level angiosperm phylogeny in APG IV (https://doi.org/10.1111/boj.12385). '
                             'ncbi_user: Infer NCBI Taxonomy ID from species name, and match the ranks with the labels of the user-provided tree. '
                             'user: User-provided tree in --infile.')
pconstrain.add_argument('--rank', metavar='no|species|genus|family|order|...', default='no', type=str, required=False,
                        action='store',
                        help='default=%(default)s: Constrain at a particular taxonomic rank and above. '
                             'For example, if "family" is specified, "genus" and "species" are not considered. '
                             'This option is currently compatible only with --backbone ncbi')
pconstrain.add_argument('--collapse', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                        help='default=%(default)s: For tip names of "GENUS_SPECIES_OTHERINFO", '
                             'drop OTHERINFO and collapse clades if GENUS_SPECIES is identical. '
                             'The output file may be used as a species tree for phylogeny reconciliation. ')
pconstrain.set_defaults(handler=command_constrain)

def command_collapse(args):
    from nwkit.collapse import collapse_main
    collapse_main(args)
pcollapse = subparsers.add_parser('collapse', help='Collapse internal branches by support and/or branch length', parents=[p_parent])
pcollapse.add_argument('--min-support', '--min_support', dest='min_support', metavar='FLOAT', default=None, type=float, required=False, action='store',
                       help='default=%(default)s: Collapse internal branches whose support is smaller than this threshold.')
pcollapse.add_argument('--max-dist', '--max_dist', dest='max_dist', metavar='FLOAT', default=None, type=float, required=False, action='store',
                       help='default=%(default)s: Collapse internal branches whose branch length is at most this threshold.')
pcollapse.add_argument('--preserve-branch-length', '--preserve_branch_length', dest='preserve_branch_length', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                       help='default=%(default)s: Add the deleted branch length to descendant branches when collapsing.')
pcollapse.set_defaults(handler=command_collapse)

def command_compose(args):
    from nwkit.compose import compose_main
    compose_main(args)
pcompose = subparsers.add_parser('compose', help='Assemble compatible roots, values, and annotations from multiple trees', parents=[p_parent])
pcompose.add_argument('--manifest', metavar='PATH', default=None, type=str, required=False, action='store',
                      help='Optional JSON manifest defining root, name, support, length, and property sources.')
pcompose.add_argument('--root-source', '--root_source', dest='root_source', metavar='PATH', default=None, type=str, required=False, action='store',
                      help='Tree providing the root split.')
pcompose.add_argument('--name-source', '--name_source', dest='name_source', metavar='PATH', default=None, type=str, required=False, action='store',
                      help='Tree providing node names.')
pcompose.add_argument('--support-source', '--support_source', dest='support_source', metavar='PATH', default=None, type=str, required=False, action='store',
                      help='Tree providing internal-node support values.')
pcompose.add_argument('--length-source', '--length_source', dest='length_source', metavar='PATH', default=None, type=str, required=False, action='store',
                      help='Tree providing branch lengths.')
pcompose.add_argument('--property-source', '--property_source', dest='property_source', metavar='SOURCE=TARGET@PATH', default=[], type=str, required=False, action='append',
                      help='Tree providing an arbitrary NHX property. SOURCE@PATH preserves its name. May be repeated.')
pcompose.add_argument('--source-format', '--source_format', dest='source_format', metavar='auto|auto-strict|INT', default='auto', type=str, required=False, action='store',
                      help='Default ETE parser format for all source trees.')
pcompose.add_argument('--taxon-mode', '--taxon_mode', dest='taxon_mode', metavar='exact|intersection', default='exact', type=str, required=False, action='store',
                      choices=['exact', 'intersection'],
                      help='Require identical tip sets or map unique clades projected onto shared tips.')
pcompose.add_argument('--match-basis', '--match_basis', dest='match_basis', metavar='clade|split', default='clade', type=str, required=False, action='store',
                      choices=['clade', 'split'],
                      help='Match rooted descendant clades or root-independent canonical edge splits.')
pcompose.add_argument('--root-edge-policy', '--root_edge_policy', dest='root_edge_policy', metavar='TARGET_PROPERTY=POLICY', default=[], type=str, required=False, action='append',
                      help='Resolve root-edge ambiguity per target property. Policies: auto, skip, equal-only, matching-side, mean, min, max, edge-total. May be repeated.')
pcompose.add_argument('--allow-projected-values', '--allow_projected_values', dest='allow_projected_values', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                      help='default=%(default)s: Permit support and branch-length transfer through projected matches. Use only when shared-tip equivalence is sufficient.')
pcompose.add_argument('--policy', metavar='compatible-only|strict', default='compatible-only', type=str, required=False, action='store',
                      choices=['compatible-only', 'strict'],
                      help='Skip incompatible values, or require exact (not projected) matches for every requested value.')
pcompose.add_argument('--report', metavar='PATH', default=None, type=str, required=False, action='store',
                      help='Optional per-clade TSV recording sources, matches, transferred values, and conflicts.')
pcompose.set_defaults(handler=command_compose)

def command_cladefreq(args):
    from nwkit.cladefreq import cladefreq_main
    cladefreq_main(args)
pcladefreq = subparsers.add_parser('cladefreq', help='Summarize clade frequencies across a tree collection', parents=[p_tree_input, p_table_output])
pcladefreq.add_argument('-r', '--reference', metavar='PATH', default=None, type=str, required=False, action='store',
                        help='default=%(default)s: Optional reference tree used to flag clades present in a named topology.')
pcladefreq.add_argument('-rf', '--reference-format', '--reference_format', dest='reference_format', metavar='auto|auto-strict|INT', default='auto', type=str, required=False, action='store',
                        help='default=%(default)s: ETE tree format for --reference.')
pcladefreq.add_argument('--weight-tsv', '--weight_tsv', dest='weight_tsv', metavar='PATH', default=None, type=str, required=False, action='store',
                        help='default=%(default)s: Optional TSV assigning a positive weight to each input tree. '
                             'The file must contain a "weight" column and may contain a 1-based "tree_id" column.')
pcladefreq.add_argument('--support-scale', '--support_scale', dest='support_scale', metavar='percent|proportion', default='percent', type=str, required=False, action='store',
                        choices=['percent', 'proportion'],
                        help='default=%(default)s: Scale used for the reported clade frequencies.')
pcladefreq.add_argument('--threads', metavar='INT', default=1, type=int, required=False, action='store',
                        help='default=%(default)s: Number of parallel workers used to parse and summarize input trees.')
pcladefreq.set_defaults(handler=command_cladefreq)

def command_consensus(args):
    from nwkit.consensus import consensus_main
    consensus_main(args)
pconsensus = subparsers.add_parser('consensus', help='Generate a consensus tree or transfer consensus support to a reference tree', parents=[p_parent])
pconsensus.add_argument('-r', '--reference', metavar='PATH', default=None, type=str, required=False, action='store',
                        help='default=%(default)s: Optional reference tree that receives consensus support values instead of building a de novo consensus tree.')
pconsensus.add_argument('-rf', '--reference-format', '--reference_format', dest='reference_format', metavar='auto|auto-strict|INT', default='auto', type=str, required=False, action='store',
                        help='default=%(default)s: ETE tree format for --reference.')
pconsensus.add_argument('--method', metavar='greedy|majority|strict', default='greedy', type=str, required=False, action='store',
                        choices=['greedy', 'majority', 'strict'],
                        help='default=%(default)s: Consensus-clade selection rule. '
                             'greedy uses --min-freq, majority keeps clades with frequency > 0.5, and strict keeps only clades present in all trees.')
pconsensus.add_argument('--min-freq', '--min_freq', dest='min_freq', metavar='0.0<=FLOAT<=1.0', default=0.5, type=float, required=False, action='store',
                        help='default=%(default)s: Minimum clade frequency retained when --method greedy is used.')
pconsensus.add_argument('--branch-length', '--branch_length', dest='branch_length', metavar='none|mean|median', default='none', type=str, required=False, action='store',
                        choices=['none', 'mean', 'median'],
                        help='default=%(default)s: How branch lengths are assigned in the de novo consensus tree.')
pconsensus.add_argument('--weight-tsv', '--weight_tsv', dest='weight_tsv', metavar='PATH', default=None, type=str, required=False, action='store',
                        help='default=%(default)s: Optional TSV assigning a positive weight to each input tree. '
                             'The file must contain a "weight" column and may contain a 1-based "tree_id" column.')
pconsensus.add_argument('--support-scale', '--support_scale', dest='support_scale', metavar='percent|proportion', default='percent', type=str, required=False, action='store',
                        choices=['percent', 'proportion'],
                        help='default=%(default)s: Scale used for consensus support values.')
pconsensus.add_argument('--threads', metavar='INT', default=1, type=int, required=False, action='store',
                        help='default=%(default)s: Number of parallel workers used to parse and summarize input trees.')
pconsensus.set_defaults(handler=command_consensus)

def command_diff(args):
    from nwkit.diff import diff_main
    diff_main(args)
pdiff = subparsers.add_parser('diff', help='Report interpretable clade, root, value, and annotation differences between trees', parents=[p_tree_input, p_table_output])
pdiff.add_argument('-i2', '--infile2', metavar='PATH', default=None, type=str, required=True, action='store',
                   help='Second Newick tree.')
pdiff.add_argument('-f2', '--format2', metavar='auto|auto-strict|INT', default='auto', type=str, required=False, action='store',
                   help='ETE parser format for --infile2.')
pdiff.add_argument('--taxon-mode', '--taxon_mode', dest='taxon_mode', metavar='exact|intersection', default='exact', type=str, required=False, action='store',
                   choices=['exact', 'intersection'],
                   help='Require identical tip sets or compare unique projections onto shared tips.')
pdiff.add_argument('--comparison', metavar='rooted|unrooted', default='rooted', type=str, required=False, action='store',
                   choices=['rooted', 'unrooted'],
                   help='Compare rooted descendant clades or root-independent edge splits.')
pdiff.add_argument('--target', metavar='all|root|leaf|intnode', default='intnode', type=str, required=False, action='store',
                   choices=['all', 'root', 'leaf', 'intnode'],
                   help='Node class included in detailed rows.')
pdiff.add_argument('--property', metavar='KEY', default=[], type=str, required=False, action='append',
                   help='Additional NHX property compared and serialized as JSON. May be repeated.')
pdiff.add_argument('--fail-on-difference', '--fail_on_difference', dest='fail_on_difference', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                   help='Exit nonzero after writing the report if any difference is detected.')
pdiff.set_defaults(handler=command_diff)

def command_dist(args):
    from nwkit.dist import dist_main
    dist_main(args)
pdist = subparsers.add_parser('dist', help='Calculate topological distance between two trees', parents=[p_tree_input, p_table_output])
pdist.add_argument('-i2', '--infile2', metavar='PATH', default=None, type=str, required=True, action='store',
                   help='default=%(default)s: Input newick file 2.')
pdist.add_argument('-f2', '--format2', metavar='auto|auto-strict|INT', default='auto', type=str, required=False, action='store',
                   help='default=%(default)s: ETE tree format for --infile2.')
pdist.add_argument('-d', '--dist', metavar='STR', default='RF', type=str, required=False, action='store',
                   help='default=%(default)s: Distance calculation method. RF=Robinson-Foulds')
pdist.set_defaults(handler=command_dist)

def command_draw(args):
    from nwkit.draw import draw_main
    draw_main(args)
pdraw = subparsers.add_parser('draw', help='Draw a phylogenetic tree with optional speciation/duplication node markers', parents=[p_tree_input, p_species, p_tip_table_policy])
pdraw.add_argument('-o', '--outfile', metavar='PATH', default='nwkit_draw.pdf', type=str, required=False, action='store',
                   help='default=%(default)s: Output image file.')
pdraw.add_argument('--image-format', '--image_format', dest='image_format', metavar='auto|pdf|png|svg', default='auto', type=str, required=False, action='store',
                   choices=['auto', 'pdf', 'png', 'svg'],
                   help='default=%(default)s: Output image format. "auto" infers from --outfile and otherwise falls back to pdf.')
pdraw.add_argument('--species-overlap-node-plot', '--species_overlap_node_plot', dest='species_overlap_node_plot', metavar='yes|no|auto', default='auto', required=False, type=str,
                   choices=['yes', 'no', 'auto'],
                   help='default=%(default)s: Show speciation/duplication node markers. '
                        '"yes" always tries to plot, "no" disables them, and "auto" enables them only when '
                        'all tip labels are parseable by the configured species parser.')
pdraw.add_argument('--ladderize', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                   help='default=%(default)s: Ladderize the input tree before drawing.')
pdraw.add_argument('--trait', metavar='PATH', default=None, type=str, required=False, action='store',
                   help='default=%(default)s: Optional TSV containing "leaf_name" and one or more grouping columns for tip-label coloring.')
pdraw.add_argument('--group-by', '--group_by', dest='group_by', metavar='STR', default=None, type=str, required=False, action='store',
                   help='default=%(default)s: Column name in --trait used to color tip labels.')
pdraw.add_argument('--trait-palette', '--trait_palette', dest='trait_palette', metavar='STR', default='tab10', type=str, required=False, action='store',
                   help='default=%(default)s: Matplotlib colormap name used for trait-group colors.')
pdraw.add_argument('--support-labels', '--support_labels', dest='support_labels', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                   help='default=%(default)s: Whether to draw support labels on branches.')
pdraw.add_argument('--support-min', '--support_min', dest='support_min', metavar='FLOAT', default=None, type=float, required=False, action='store',
                   help='default=%(default)s: Show only support labels greater than or equal to this value.')
pdraw.set_defaults(handler=command_draw)

def command_drop(args):
    from nwkit.drop import drop_main
    drop_main(args)
pdrop = subparsers.add_parser('drop', help='Remove node and branch information', parents=[p_parent])
pdrop.add_argument('-t', '--target', metavar='all|root|leaf|intnode', default='all', type=str, required=False,
                   action='store', choices=['all', 'root', 'leaf', 'intnode'],
                   help='default=%(default)s: Nodes to be edited.')
pdrop.add_argument('--name', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                   help='default=%(default)s: Drop node names.')
pdrop.add_argument('--support', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                   help='default=%(default)s: Drop support values.')
pdrop.add_argument('--length', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                   help='default=%(default)s: Drop branch length.')
pdrop.add_argument('--fill', metavar='STR/NUMERIC', default=None, type=str, required=False, action='store',
                   help='default=%(default)s: Fill values instead of simply dropping them.')
pdrop.set_defaults(handler=command_drop)

def command_info(args):
    from nwkit.info import info_main
    info_main(args)
pinfo = subparsers.add_parser('info', help='Print tree information', parents=[p_tree_input, p_text_output, p_species])
pinfo.set_defaults(handler=command_info)

def command_image(args):
    from nwkit.image import image_main
    image_main(args)
pimage = subparsers.add_parser('image', help='Retrieve representative species images with license-aware filtering', parents=[p_tree_input, p_download, p_species])
pimage.add_argument('--out-dir', '--out_dir', dest='out_dir', metavar='PATH', default=None, type=str, required=True, action='store',
                    help='Output directory for images and metadata.')
pimage.add_argument('--style', metavar='auto|photo|silhouette', default='auto', type=str, required=False, action='store',
                    choices=['auto', 'photo', 'silhouette'],
                    help='default=%(default)s: Image preference mode.')
pimage.add_argument('--source', metavar='phylopic,bioicons,inaturalist,wikimedia,gbif,eol,idigbio,openverse,ncbi', default=None, type=str, required=False, action='store',
                    help='default=style-dependent: Comma-separated provider priority list. '
                         'Supported sources are phylopic, bioicons, inaturalist, wikimedia, gbif, eol, idigbio, openverse, and ncbi. '
                         'When ncbi appears after other providers, it is used as a lazy fallback because it may download the NCBI taxonomy image dump.')
pimage.add_argument('--license-max', '--license_max', dest='license_max', metavar='public-domain|cc-by|cc-by-sa|cc-by-nc|cc-by-nc-sa|any', default='cc-by-nc-sa', type=str, required=False, action='store',
                    choices=['public-domain', 'cc-by', 'cc-by-sa', 'cc-by-nc', 'cc-by-nc-sa', 'any'],
                    help='default=%(default)s: Strongest restriction level allowed for candidate images.')
pimage.add_argument('--allow-nd', '--allow_nd', dest='allow_nd', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                    help='default=%(default)s: Whether ND licenses are acceptable.')
pimage.add_argument('--fallback-rank', '--fallback_rank', dest='fallback_rank', metavar='none|genus|family', default='none', type=str, required=False, action='store',
                    choices=['none', 'genus', 'family'],
                    help='default=%(default)s: Allow fallback matching above the species rank.')
pimage.add_argument('--max-per-species', '--max_per_species', dest='max_per_species', metavar='INT', default=1, type=int, required=False, action='store',
                    help='default=%(default)s: Number of images to keep per species.')
pimage.add_argument('--species-name-tsv', '--species_name_tsv', dest='species_name_tsv', metavar='PATH', default=None, type=str, required=False, action='store',
                    help='default=%(default)s: Optional legacy-style TSV with "leaf_name" and "species_name" columns. Prefer --species-map-tsv.')
pimage.add_argument('--name-tsv', '--name_tsv', dest='species_name_tsv', metavar='PATH', type=str, help=argparse.SUPPRESS)
pimage.add_argument('--manifest-out', '--manifest_out', dest='manifest_out', metavar='PATH', default=None, type=str, required=False, action='store',
                    help='default=%(default)s: Optional custom output path for manifest.tsv.')
pimage.add_argument('--manifest', dest='manifest_out', metavar='PATH', type=str, help=argparse.SUPPRESS)
pimage.add_argument('--attribution-out', '--attribution_out', dest='attribution_out', metavar='PATH', default=None, type=str, required=False, action='store',
                    help='default=%(default)s: Optional custom output path for ATTRIBUTION.md.')
pimage.add_argument('--attribution', dest='attribution_out', metavar='PATH', type=str, help=argparse.SUPPRESS)
pimage.add_argument('--fail-on-missing', '--fail_on_missing', dest='fail_on_missing', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                    help='default=%(default)s: Whether unresolved species should trigger a non-zero exit code.')
pimage.add_argument('--output-format', '--output_format', dest='output_format', metavar='original|png|jpg', default='original', type=str, required=False, action='store',
                    choices=['original', 'png', 'jpg'],
                    help='default=%(default)s: Optional normalized output format. '
                         'Image post-processing options require the optional "nwkit[image]" dependencies.')
pimage.add_argument('--max-edge', '--max_edge', dest='max_edge', metavar='INT', default=None, type=int, required=False, action='store',
                    help='default=%(default)s: Downscale rasterized output so the longest edge is at most this value.')
pimage.add_argument('--canvas', metavar='none|square', default='none', type=str, required=False, action='store',
                    choices=['none', 'square'],
                    help='default=%(default)s: Add padding to a square canvas after trimming and resizing.')
pimage.add_argument('--background', metavar='white|transparent', default='white', type=str, required=False, action='store',
                    choices=['white', 'transparent'],
                    help='default=%(default)s: Background color used when --canvas square is requested.')
pimage.add_argument('--trim', metavar='off|white|transparent|semantic', default='off', type=str, required=False, action='store',
                    choices=['off', 'white', 'transparent', 'semantic'],
                    help='default=%(default)s: Trim uniform white margins, transparent margins, or an estimated main foreground region before resizing.')
pimage.add_argument('--trim-shape', '--trim_shape', dest='trim_shape', metavar='bbox|square', default='bbox', type=str, required=False, action='store',
                    choices=['bbox', 'square'],
                    help='default=%(default)s: Shape of the trimmed result. "square" center-crops the trimmed content to a square and may clip it.')
pimage.set_defaults(handler=command_image)

def command_intersection(args):
    from nwkit.intersection import intersection_main
    intersection_main(args)
pintersection = subparsers.add_parser('intersection', help='Drop non-overlapping leaves/sequences in 2 trees or tree+alignment', parents=[p_parent])
pintersection.add_argument('-i2', '--infile2', metavar='PATH', default='', type=str, required=False, action='store',
                           help='default=%(default)s: Input newick file 2. The intersected version of this file '
                                'will not be generated, so if necessary, replace --infile and --infile2 and run again.')
pintersection.add_argument('-f2', '--format2', metavar='auto|auto-strict|INT', default='auto', type=str, required=False, action='store',
                           help='default=%(default)s: ETE tree format for --infile2.')
pintersection.add_argument('-si', '--seqin', metavar='PATH', default='', type=str, required=False, action='store',
                           help='default=%(default)s: Input sequence file.')
pintersection.add_argument('-so', '--seqout', metavar='PATH', default='', type=str, required=False, action='store',
                           help='default=%(default)s: Output sequence file.')
pintersection.add_argument('-sf', '--seqformat', metavar='fasta', default='fasta', type=str, required=False, action='store',
                           choices=['fasta'], help='default=%(default)s: Sequence format for --seqin and --seqout.')
pintersection.add_argument('--match', metavar='complete|prefix|backward', default='complete', type=str,
                           required=False, action='store', choices=['complete','prefix','backward'],
                           help='default=%(default)s: Method for ID matching.')
pintersection.set_defaults(handler=command_intersection)

def command_label(args):
    from nwkit.label import label_main
    label_main(args)
plabel = subparsers.add_parser('label', help='Add unique node labels', parents=[p_parent])
plabel.add_argument('-t', '--target', metavar='all|root|leaf|intnode', default='all', type=str, required=False,
                   action='store', choices=['all', 'root', 'leaf', 'intnode'],
                   help='default=%(default)s: Nodes to be edited.')
plabel.add_argument('--prefix', metavar='STR', default='n', type=str, required=False, action='store',
                   help='default=%(default)s: Prefix for node labels.')
plabel.add_argument('--force', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                   help='default=%(default)s: Whether to overwrite existing node names.')
plabel.set_defaults(handler=command_label)

def command_rename(args):
    from nwkit.rename import rename_main
    rename_main(args)
prename = subparsers.add_parser('rename', help='Rename nodes using a TSV mapping or regular expression', parents=[p_parent])
prename.add_argument('--name-tsv', '--name_tsv', dest='name_tsv', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: TSV containing "old_name" and "new_name" columns.')
prename.add_argument('--pattern', metavar='REGEX', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: Optional regular expression used to rename target nodes.')
prename.add_argument('--replacement', metavar='STR', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: Replacement string used with --pattern.')
prename.add_argument('-t', '--target', metavar='all|root|leaf|intnode', default='leaf', type=str, required=False,
                     action='store', choices=['all', 'root', 'leaf', 'intnode'],
                     help='default=%(default)s: Nodes to be renamed.')
prename.add_argument('--require-all-old-names', '--require_all_old_names', dest='require_all_old_names', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Whether every old_name in --name-tsv must match a target node.')
prename.add_argument('--require-match', '--require_match', dest='require_match', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Whether regex mode requires at least one target node name to match.')
prename.add_argument('--check-leaf-uniqueness', '--check_leaf_uniqueness', dest='check_leaf_uniqueness', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Whether to fail if the renamed tree has duplicated or empty leaf labels.')
prename.set_defaults(handler=command_rename)

def command_mark(args):
    from nwkit.mark import mark_main
    mark_main(args)
pmark = subparsers.add_parser('mark', help='Add texts to node labels', parents=[p_parent])
pmark.add_argument('-p', '--pattern', metavar='REGEX', default='.*', type=str, required=False, action='store',
                   help='default=%(default)s: Regular expression for label search.')
pmark.add_argument('-t', '--target', metavar='mrca|clade|leaf', default='clade', type=str, required=False,
                   action='store', choices=['mrca', 'clade', 'leaf'],
                   help='default=%(default)s: Nodes to be marked.')
pmark.add_argument('--target-only-clade', '--target_only_clade', dest='target_only_clade', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                   help='default=%(default)s: Mark the label of MRCA/clade whose clade contains only target leaves. '
                        'Use with --target mrca/clade')
pmark.add_argument('--insert-text', '--insert-txt', '--insert_txt', dest='insert_txt', metavar='STR', default=None, type=str, required=True, action='store',
                   help='default=%(default)s: Label to insert to the target node labels.')
pmark.add_argument('--insert-separator', '--insert-sep', '--insert_sep', dest='insert_sep', metavar='STR', default='', type=str, required=False, action='store',
                   help='default=%(default)s: Separator for --insert-text.')
pmark.add_argument('--insert-position', '--insert-pos', '--insert_pos', dest='insert_pos', metavar='prefix|suffix', default='suffix', type=str, required=False,
                   action='store', choices=['prefix', 'suffix'],
                   help='default=%(default)s: Place to insert --insert-text.')
pmark.set_defaults(handler=command_mark)

def command_mcmctree(args):
    from nwkit.mcmctree import mcmctree_main
    mcmctree_main(args)
pmcmctree = subparsers.add_parser('mcmctree', help='Introduce divergence time constraints for PAML\'s mcmctree', parents=[p_tree_input, p_download, p_species])
pmcmctree.add_argument('-o', '--outfile', metavar='PATH', default='-', type=str, required=False, action='store',
                       help='default=%(default)s: Output MCMCtree Newick file. Use "-" for STDOUT.')
pmcmctree.add_argument('--left-species', '--left_species', dest='left_species', metavar='STR', default=None, type=str, required=False,
                       help='default=%(default)s: Any species in the left clade. '
                            'If you want to set a bound on the node splitting Homo_sapiens and Mus_musculus, '
                            'specify one of them (e.g., Homo_sapiens).')
pmcmctree.add_argument('--right-species', '--right_species', dest='right_species', metavar='STR', default=None, type=str, required=False,
                       help='default=%(default)s: Any species in the right clade deriving from the common ancestor. If '
                            'you want to set a bound on the node splitting Homo_sapiens and Mus_musculus, specify the '
                            'other one that is not used as the left species (e.g., Mus_musculus).')
pmcmctree.add_argument('--lower-bound', '--lower_bound', dest='lower_bound', metavar='FLOAT', default=None, type=str,
                       help='default=%(default)s: Lower bound of the calibration point.')
pmcmctree.add_argument('--lower-offset', '--lower_offset', dest='lower_offset', metavar='FLOAT', default='0.1', type=str, help='default=%(default)s: ')
pmcmctree.add_argument('--lower-scale', '--lower_scale', dest='lower_scale', metavar='FLOAT', default='1', type=str, help='default=%(default)s: ')
pmcmctree.add_argument('--lower-tail-prob', '--lower_tail_prob', '--lower_tailProb', dest='lower_tail_prob', metavar='FLOAT', default='0.025', type=str,
                       help='default=%(default)s: Lower tail probability. Use 1e-300 for hard bound. Default=0.025')
pmcmctree.add_argument('--upper-bound', '--upper_bound', dest='upper_bound', metavar='FLOAT', default=None, type=str,
                       help='default=%(default)s: Upper bound of the calibration point. A point estimate can be '
                            'specified by setting the same age in both lower and upper bounds '
                            '(e.g., --lower-bound 5.2 --upper-bound 5.2)')
pmcmctree.add_argument('--upper-tail-prob', '--upper_tail_prob', '--upper_tailProb', dest='upper_tail_prob', metavar='FLOAT', default='0.025', type=str,
                       help='default=%(default)s: Upper tail probability. Use 1e-300 for hard bound. Default=0.025')
pmcmctree.add_argument('--add-header', '--add_header', dest='add_header', metavar='yes|no', nargs='?', const=True, default=False, type=strtobool,
                       help='default=%(default)s: Add the header required for mcmctree.')
pmcmctree.add_argument('--timetree', metavar='point|ci|no', default='no', type=str, required=False, action='store',
                       choices=['point', 'ci', 'no'],
                       help='default=%(default)s: Obtain the divergence time from timetree.org. '
                            'For point/ci, tip labels are parsed with --species-parser/--species-regex/--species-map-tsv. '
                            'Repeated species labels are allowed only when their tips are monophyletic. '
                            'point: point estimate, ci: 95 percent confidence interval as upper and lower bounds. '
                            'no: disable the function.')
pmcmctree.add_argument('--min-clade-prop', '--min_clade_prop', dest='min_clade_prop', metavar='0.0<=FLOAT<=1.0', default=0, type=float, required=False, action='store',
                       help='default=%(default)s: Minimum proportion of the clade size to the total number of species. '
                            'If the clade proportion is smaller than this value, time constraints are removed.')
pmcmctree.add_argument('--higher-rank-search', '--higher_rank_search', dest='higher_rank_search', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                       help='default=%(default)s: Attempt to obtain timetree data using the taxids for higher taxonomic ranks '                            
                            'if the species-level search failed.')
pmcmctree.add_argument('--threads', metavar='INT', default=1, type=int, required=False, action='store',
                       help='default=%(default)s: Number of parallel workers used for TimeTree HTTP requests at the same taxonomic rank.')
pmcmctree.set_defaults(handler=command_mcmctree)

def command_monophyly(args):
    from nwkit.monophyly import monophyly_main
    monophyly_main(args)
pmonophyly = subparsers.add_parser('monophyly', help='Assess whether species or trait-defined groups are monophyletic', parents=[p_tree_input, p_table_output, p_species, p_tip_table_policy])
pmonophyly.add_argument('--trait', metavar='PATH', default=None, type=str, required=False, action='store',
                        help='default=%(default)s: Optional TSV containing a "leaf_name" column and one or more grouping columns.')
pmonophyly.add_argument('--group-by', '--group_by', dest='group_by', metavar='STR', default=None, type=str, required=False, action='store',
                        help='default=%(default)s: Column name in --trait used to define groups. '
                             'When omitted, groups are inferred from species labels parsed from the leaf names.')
pmonophyly.add_argument('--unrooted', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                        help='default=%(default)s: Evaluate monophyly in unrooted mode.')
pmonophyly.add_argument('--fail-on-non-monophyly', '--fail_on_non_monophyly', dest='fail_on_non_monophyly', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                        help='default=%(default)s: Exit with a non-zero status if any group is not monophyletic.')
pmonophyly.set_defaults(handler=command_monophyly)

def command_nhx2nwk(args):
    from nwkit.nhx2nwk import nhx2nwk_main
    nhx2nwk_main(args)
pnhx2nwk = subparsers.add_parser('nhx2nwk', help='Generate Newick from NHX', parents=[p_parent])
pnhx2nwk.add_argument('-p', '--node-label', '--node_label', dest='node_label', metavar='B|D|H|S|...', default='', type=str, required=False, action='store',
                         help='default=%(default)s: NHX attribute to use as internal node labels.')
pnhx2nwk.set_defaults(handler=command_nhx2nwk)

def command_nwk2table(args):
    from nwkit.nwk2table import nwk2table_main
    nwk2table_main(args)
pnwk2table = subparsers.add_parser('nwk2table', help='Convert a Newick tree into a parent-child table', parents=[p_tree_input, p_table_output])
pnwk2table.add_argument('--age', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                        help='default=%(default)s: Add node ages for ultrametric trees.')
pnwk2table.add_argument('--sister', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                        help='default=%(default)s: Add the branch_id of one sister node when available.')
pnwk2table.set_defaults(handler=command_nwk2table)

def command_printlabel(args):
    from nwkit.printlabel import printlabel_main
    printlabel_main(args)
pprintlabel = subparsers.add_parser('printlabel', help='Search and print node labels', parents=[p_tree_input, p_text_output])
pprintlabel.add_argument('-p', '--pattern', metavar='REGEX', default='.*', type=str, required=False, action='store',
                         help='default=%(default)s: Regular expression for label search.')
pprintlabel.add_argument('-t', '--target', metavar='all|root|leaf|intnode', default='all', type=str, required=False,
                         action='store',
                         choices=['all', 'root', 'leaf', 'intnode'],
                         help='default=%(default)s: Nodes to be searched.')
pprintlabel.add_argument('--sister', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                         help='default=%(default)s: Show labels of the sisters instead of targets.')
pprintlabel.set_defaults(handler=command_printlabel)

def command_prune(args):
    from nwkit.prune import prune_main
    prune_main(args)
pprune = subparsers.add_parser('prune', help='Prune leaves', parents=[p_parent])
pprune.add_argument('-p', '--pattern', metavar='REGEX', default='.*', type=str, required=False, action='store',
                         help='default=%(default)s: Regular expression for label search.')
pprune.add_argument('--invert-match', '--invert_match', dest='invert_match', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                         help='default=%(default)s: Prune unmatched leaves.')
pprune.set_defaults(handler=command_prune)

def command_rescale(args):
    from nwkit.rescale import rescale_main
    rescale_main(args)
prescale = subparsers.add_parser('rescale', help='Rescale branch length with a given factor', parents=[p_parent])
prescale.add_argument('-t', '--target', metavar='all|root|leaf|intnode', default='all', type=str, required=False,
                   action='store',
                   choices=['all', 'root', 'leaf', 'intnode'],
                   help='default=%(default)s: Nodes to be edited.')
prescale.add_argument('--factor', metavar='FLOAT', default=None, type=float, required=True, action='store',
                   help='default=%(default)s: Rescaling factor of branch length.')
prescale.set_defaults(handler=command_rescale)

def command_root(args):
    from nwkit.root import root_main
    root_main(args)
proot = subparsers.add_parser('root', help='Place or transfer the tree root', parents=[p_parent, p_download, p_species])
proot.add_argument('-i2', '--infile2', metavar='PATH', default='', type=str, required=False, action='store',
                   help='default=%(default)s: Input newick file 2. Used when --method "transfer". '
                        'Leaf labels should be matched to those in --infile.')
proot.add_argument('-f2', '--format2', metavar='auto|auto-strict|INT', default='auto', type=str, required=False, action='store',
                   help='default=%(default)s: ETE tree format for --infile2.')
proot.add_argument('--taxon-mode', '--taxon_mode', dest='taxon_mode', metavar='exact|intersection', default='exact', type=str, required=False, action='store',
                   choices=['exact', 'intersection'],
                   help='For root transfer, require identical tips or use an unambiguous root split projected onto shared tips.')
proot.add_argument('--method', metavar='STR', default='midpoint', type=str, required=False, action='store',
                   choices=['midpoint','outgroup','transfer','mad','mv','taxonomy',],
                   help='default=%(default)s: '
                        'midpoint: Midpoint rooting. '
                        'outgroup: Outgroup rooting with --outgroup. '
                        'transfer: Transfer the root position from --infile2 to --infile. '
                        'The two trees should have the same bipartitions at the root node. '
                        'mad: Minimal Ancestor Deviation rooting (Tria et al. 2017). '
                        'mv: Minimum Variance rooting (Mai et al. 2017). '
                        'taxonomy: Root by transferring a taxonomy-derived root split onto --infile.')
proot.add_argument('--outgroup', metavar='STR', default='', type=str, required=False, action='store',
                   help='default=%(default)s: An outgroup label or a comma-separated list of outgroup labels. '
                        'For the latter, the clade containing all specified labels are used as an outgroup, '
                        'so all labels do not have to be specified for a large clade.')
proot.add_argument('--taxonomy-source', '--taxonomy_source', dest='taxonomy_source', metavar='ncbi[,opentree,timetree,...]', default='ncbi,opentree,timetree', type=str, required=False, action='store',
                   help='default=%(default)s: Comma-separated taxonomy/tree sources used when --method "taxonomy" is selected. '
                        'Each source is tried in order until rooting succeeds. When taxids are inferred from leaf labels, '
                        'taxonomy queries are controlled by --species-parser/--species-regex/--species-map-tsv. '
                        'Repeated species labels are allowed for timetree/opentree only when their tips are monophyletic. '
                        'Supported sources: ncbi, timetree, opentree.')
proot.add_argument('--taxid-tsv', '--taxid_tsv', dest='taxid_tsv', metavar='PATH', default=None, type=str, required=False, action='store',
                   help='default=%(default)s: TSV file containing tree leaf labels in the "leaf_name" column and '
                        'their taxonomy IDs in the "taxid" column. Used when --method "taxonomy". '
                        'If omitted, taxonomy IDs are inferred from leaf labels. This option is used by the NCBI source.')
proot.add_argument('--rank', metavar='no|species|genus|family|order|...', default='no', type=str, required=False,
                   action='store',
                   help='default=%(default)s: Highest taxonomic rank retained when deriving the taxonomy tree for '
                        '--method "taxonomy".')
proot.set_defaults(handler=command_root)

def command_sanitize(args):
    from nwkit.sanitize import sanitize_main
    sanitize_main(args)
psanitize = subparsers.add_parser('sanitize', help='Eliminate non-standard Newick flavors', parents=[p_parent])
psanitize.add_argument('--remove-singleton', '--remove_singleton', dest='remove_singleton', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                       help='default=%(default)s: Remove singleton nodes represented as double brackets.')
psanitize.add_argument('--resolve-polytomy', '--resolve_polytomy', dest='resolve_polytomy', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                       help='default=%(default)s: Resolve multifurcation (polytomy) nodes into dichotomies with zero-length branches.')
psanitize.add_argument('--name-quote', '--name_quote', dest='name_quote', metavar='none|single|double', default='none', type=str, required=False,
                       action='store', choices=['none','single','double'],
                       help='default=%(default)s: Quotation of node and leaf names.'
                            'none = no quote, single = \', double = \" ')
psanitize.set_defaults(handler=command_sanitize)

def command_shuffle(args):
    from nwkit.shuffle import shuffle_main
    shuffle_main(args)
pshuffle = subparsers.add_parser('shuffle', help='Shuffle branches and/or labels', parents=[p_parent])
pshuffle.add_argument('--topology', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                       help='default=%(default)s: Randomize entire tree topology and branch length. '
                            'Without --label yes, new topology preserve the leaf label orders and is not completely randomized.')
pshuffle.add_argument('--branch-length', '--branch_length', dest='branch_length', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                       help='default=%(default)s: Shuffle branch length. Automatically activated when --topology yes.')
pshuffle.add_argument('--label', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                       help='default=%(default)s: Shuffle leaf labels.')
pshuffle.add_argument('--seed', metavar='INT', default=None, type=int, required=False, action='store',
                      help='default=%(default)s: Random seed for reproducible topology, branch-length, and label shuffling.')
pshuffle.set_defaults(handler=command_shuffle)

def command_skim(args):
    from nwkit.skim import skim_main
    skim_main(args)
pskim = subparsers.add_parser('skim', help='Sample leaves from clades with shared traits', parents=[p_parent, p_tip_table_policy])
pskim.add_argument('--trait', metavar='PATH', default=None, type=str, required=False, action='store',
                       help='default=%(default)s: Path to a trait table (TSV) containing leaf names in the "leaf_name" column and one or more trait columns. '
                            'If not specified, all leaves are treated as a single group and randomly sampled.')
pskim.add_argument('--group-by', '--group_by', dest='group_by', metavar='STR', default=None, type=str, required=False, action='store',
                       help='default=%(default)s: Column name in the trait table used to group leaves before sampling. '
                            'Leaves in a clade are treated as a group if all non-missing values in the specified column are identical. '
                            'If not specified, all leaves are treated as a single group.')
pskim.add_argument('--retain-per-clade', '--retain_per_clade', dest='retain_per_clade', metavar='INT', default=1, type=int, required=False, action='store',
                       help='default=%(default)s: Number of leaves to retain per clade (group). '
                            'If a clade (group) contains fewer leaves than this value, all of its leaves are retained.')
pskim.add_argument('--prioritize-non-missing', '--prioritize_non_missing', dest='prioritize_non_missing', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                       help='default=%(default)s: Whether to prioritize leaves with non-missing trait values when sampling.')
pskim.add_argument('--filter-by', '--filter_by', dest='filter_by', metavar='STR', default=None, type=str, required=False, action='store',
                       help='default=%(default)s: Column name in the trait table used to rank leaves within each group before sampling. '
                            'If not specified, leaves are randomly sampled within each group.')
pskim.add_argument('--filter-mode', '--filter_mode', dest='filter_mode', metavar='ascending|descending', default='ascending', type=str, required=False, action='store',
                       help='default=%(default)s: Sorting order for --filter-by values. '
                            'If multiple leaves within a group have the same value, they are randomly sampled.')
pskim.add_argument('--only-contrastive-clades', '--only_contrastive_clades', dest='only_contrastive_clades', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                       help='default=%(default)s: Whether to output a pruned tree retaining only minimal clades with multiple non-missing trait values.')
pskim.add_argument('--group-table-prefix', '--group_table_prefix', dest='group_table_prefix', metavar='PATH', default=None, type=str, required=False, action='store',
                       help='default=%(default)s: Output prefix for group-assignment TSVs named <PATH>.all.tsv and <PATH>.sampled.tsv.')
pskim.add_argument('--output-groupfile', '--output_groupfile', dest='output_groupfile', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                       help=argparse.SUPPRESS)
pskim.add_argument('--seed', metavar='INT', default=None, type=int, required=False, action='store',
                   help='default=%(default)s: Random seed for reproducible sampling and tie-breaking.')
pskim.set_defaults(handler=command_skim)

def command_sample(args):
    from nwkit.sample import sample_main
    sample_main(args)
psample = subparsers.add_parser('sample', help='Select a representative leaf subset by filters, ranks, and sampling method', parents=[p_parent, p_tip_table_policy])
psample.add_argument('-n', '--n', metavar='INT', default=1, type=int, required=False, action='store',
                     help='default=%(default)s: Number of leaves to select.')
psample.add_argument('--method', metavar='max-pd|ranked', default='max-pd', type=str, required=False, action='store',
                     choices=['max-pd', 'ranked'],
                     help='default=%(default)s: Leaf selection method. '
                          'max-pd greedily maximizes phylogenetic diversity; ranked selects the first N leaves after filters and ranking.')
psample.add_argument('--trait', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: Optional TSV containing leaf names in the "leaf_name" column and metadata columns used by --filter/--rank. '
                          'If not specified, all leaves are sampled from the input tree.')
psample.add_argument('--filter', metavar='COLUMN:OP:VALUE', default=[], type=str, required=False, action='append',
                     help='default=%(default)s: Keep leaves whose metadata column satisfies a condition. '
                          'Can be specified multiple times. Supported operators: ge, gt, le, lt, eq, ne. '
                          'Example: --filter busco_complete_pct:ge:80 --filter num_seq:le:200000')
psample.add_argument('--rank', metavar='COLUMN:asc|desc', default=[], type=str, required=False, action='append',
                     help='default=%(default)s: Rank retained leaves before sampling. Can be specified multiple times. '
                          'For max-pd, rank order is used as a deterministic tie-breaker. '
                          'Example: --rank num_seq:asc --rank busco_complete_pct:desc')
psample.add_argument('--allow-fewer', '--allow_fewer', dest='allow_fewer', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                     help='default=%(default)s: Allow fewer than --n leaves when filters leave too few candidates.')
psample.add_argument('--report', '--output-table', '--output_table', dest='report', metavar='PATH', default=None, type=str, required=False, action='store',
                     help='default=%(default)s: Optional TSV report for selected leaves, including selection order and metadata columns.')
psample.set_defaults(handler=command_sample)

def command_subtree(args):
    from nwkit.subtree import subtree_main
    subtree_main(args)
psubtree = subparsers.add_parser('subtree', help='Generate a subtree Newick file', parents=[p_parent, p_species])
psubtree.add_argument('--left-leaf', '--left_leaf', dest='left_leaf', metavar='STR', default=None, type=str, required=False,
                      help='default=%(default)s: Any leaf names in the left clade. For example, '
                           'to extract the subtree with the root node splitting Homo_sapiens and Mus_musculus, '
                           'specify one of them (e.g., Homo_sapiens).')
psubtree.add_argument('--right-leaf', '--right_leaf', dest='right_leaf', metavar='STR', default=None, type=str, required=False,
                      help='default=%(default)s: Any leaf names in the right clade. For example, '
                           'to extract the subtree with the root node splitting Homo_sapiens and Mus_musculus, '
                           'specify the other one that is not used as --left-leaf (e.g., Mus_musculus).')
psubtree.add_argument('--leaves', metavar='leaf1,leaf2,leaf3,...', default=None, type=str, required=False,
                      help='default=%(default)s: Comma-separated list of leaves. '
                           'The output subtree has their most-recent common ancestor as the root. '
                           '--left-leaf and --right-leaf are ignored if this option is specified. '
                           'Single leaf name may be specified in combination with --orthogroup yes.')
psubtree.add_argument('--orthogroup', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                       help='default=%(default)s: The output subtree represents orthogroup(s) that contain all '
                            'specified leaves. Species assignment is controlled by '
                            '--species-parser/--species-regex/--species-map-tsv.')
psubtree.add_argument('--dup-conf-score-threshold', '--dup_conf_score_threshold', dest='dup_conf_score_threshold', metavar='FLOAT', default=0, type=float, required=False, action='store',
                       help='default=%(default)s: The threshold of duplication-confidence score for orthogroup delimitation. '
                            '0 = most stringent, 1 = most relaxed. '
                            'For the score, see https://www.ensembl.org/info/genome/compara/homology_types.html')
psubtree.set_defaults(handler=command_subtree)

def command_transfer(args):
    from nwkit.transfer import transfer_main
    transfer_main(args)
ptransfer = subparsers.add_parser('transfer', help='Transfer values and arbitrary properties between compatible tree clades', parents=[p_parent])
ptransfer.add_argument('-i2', '--infile2', metavar='PATH', default=None, type=str, required=True, action='store',
                       help='Source Newick tree. Topologies may differ; taxon matching is controlled by --taxon-mode.')
ptransfer.add_argument('-f2', '--format2', metavar='auto|auto-strict|INT', default='auto', type=str, required=False, action='store',
                       help='default=%(default)s: ETE tree format for --infile2.')
ptransfer.add_argument('-t', '--target', metavar='all|root|leaf|intnode', default='all', type=str, required=False,
                       action='store', choices=['all', 'root', 'leaf', 'intnode'],
                       help='default=%(default)s: Nodes to be edited.')
ptransfer.add_argument('--name', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                       help='default=%(default)s: transfer node names.')
ptransfer.add_argument('--support', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                       help='default=%(default)s: transfer support values.')
ptransfer.add_argument('--length', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                       help='default=%(default)s: transfer branch length.')
ptransfer.add_argument('--property', metavar='KEY', default=[], type=str, required=False, action='append',
                       help='Transfer an arbitrary NHX property without renaming it. May be repeated or comma-separated.')
ptransfer.add_argument('--property-map', '--property_map', dest='property_map', metavar='SOURCE=TARGET', default=[], type=str, required=False, action='append',
                       help='Transfer and rename an arbitrary NHX property. May be repeated.')
ptransfer.add_argument('--fill', metavar='STR/NUMERIC', default=None, type=str, required=False, action='store',
                       help='default=%(default)s: Fill values instead of leaving as is, if no corresponding node is found.')
ptransfer.add_argument('--taxon-mode', '--taxon_mode', dest='taxon_mode', metavar='exact|intersection', default='exact', type=str, required=False, action='store',
                       choices=['exact', 'intersection'],
                       help='Require identical tips or map only unique projections onto shared tips.')
ptransfer.add_argument('--match-basis', '--match_basis', dest='match_basis', metavar='clade|split', default='clade', type=str, required=False, action='store',
                       choices=['clade', 'split'],
                       help='Match rooted descendant clades or root-independent canonical edge splits.')
ptransfer.add_argument('--root-edge-policy', '--root_edge_policy', dest='root_edge_policy', metavar='TARGET_PROPERTY=POLICY', default=[], type=str, required=False, action='append',
                       help='Resolve root-edge ambiguity per target property. Policies: auto, skip, equal-only, matching-side, mean, min, max, edge-total. May be repeated.')
ptransfer.add_argument('--allow-projected-values', '--allow_projected_values', dest='allow_projected_values', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                       help='default=%(default)s: Permit support and branch-length transfer through projected matches. Use only when shared-tip equivalence is sufficient.')
ptransfer.add_argument('--policy', metavar='compatible-only|strict', default='compatible-only', type=str, required=False, action='store',
                       choices=['compatible-only', 'strict'],
                       help='Skip incompatible requests, or require exact (not projected) matches for every requested value.')
ptransfer.add_argument('--report', metavar='PATH', default=None, type=str, required=False, action='store',
                       help='Optional per-node TSV with matches, values, skipped requests, and ambiguity reasons.')
ptransfer.set_defaults(handler=command_transfer)

def command_validate(args):
    from nwkit.validate import validate_main
    validate_main(args)
pvalidate = subparsers.add_parser('validate', help='Validate one or more Newick trees and report structural issues', parents=[p_tree_input, p_table_output, p_species])
pvalidate.add_argument('--check-species', '--check_species', dest='check_species', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                       help='default=%(default)s: Check whether leaf labels are parseable as species labels and species groups are monophyletic.')
pvalidate.add_argument('--require-rooted', '--require_rooted', dest='require_rooted', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                       help='default=%(default)s: Mark unrooted trees as invalid.')
pvalidate.add_argument('--require-ultrametric', '--require_ultrametric', dest='require_ultrametric', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                       help='default=%(default)s: Mark non-ultrametric trees as invalid.')
pvalidate.add_argument('--require-same-leaf-set', '--require_same_leaf_set', dest='require_same_leaf_set', metavar='yes|no', default='yes', type=strtobool, required=False, action='store',
                       help='default=%(default)s: For multiple input trees, mark trees whose leaf set differs from the first tree as invalid.')
pvalidate.add_argument('--require-same-rooting', '--require_same_rooting', dest='require_same_rooting', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                       help='default=%(default)s: For multiple input trees with the same leaf set, mark trees whose rooted/unrooted status or root bipartition differs from the first parsed tree as invalid.')
pvalidate.add_argument('--require-binary', '--require_binary', dest='require_binary', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                       help='default=%(default)s: Mark trees with singleton or multifurcating internal nodes as invalid.')
pvalidate.add_argument('--require-all-support', '--require_all_support', dest='require_all_support', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                       help='default=%(default)s: Mark trees with missing internal support values as invalid.')
pvalidate.add_argument('--require-unambiguous-format', '--require_unambiguous_format', dest='require_unambiguous_format', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                       help='default=%(default)s: Mark trees with auto-format ambiguity as invalid.')
pvalidate.add_argument('--require-unquoted-names', '--require_unquoted_names', dest='require_unquoted_names', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                       help='default=%(default)s: Mark trees containing quoted node names as invalid.')
pvalidate.add_argument('--fail-on-issue', '--fail_on_issue', dest='fail_on_issue', metavar='yes|no', default='no', type=strtobool, required=False, action='store',
                       help='default=%(default)s: Exit with a non-zero status if any reported issue is found.')
pvalidate.set_defaults(handler=command_validate)

def command_table2nwk(args):
    from nwkit.table2nwk import table2nwk_main
    table2nwk_main(args)
ptable2nwk = subparsers.add_parser('table2nwk', help='Convert a parent-child table into a Newick tree', parents=[p_table_input])
ptable2nwk.set_defaults(handler=command_table2nwk)

def command_help(args):
    print(parser.parse_args([args.command, '--help']))

parser_help = subparsers.add_parser('help', help='Show help messages')
parser_help.add_argument('command', help='command name which help is shown')
parser_help.set_defaults(handler=command_help)


def _validate_stdin_ownership(args):
    stdin_options = [option for _dest, option in get_stdin_input_options(args)]
    if len(stdin_options) > 1:
        raise ValueError(
            'STDIN can be assigned to only one input option at a time: {}'.format(
                ', '.join(stdin_options)
            )
        )


LEGACY_OPTION_CANONICAL_OVERRIDES = {
    'skim': {
        '--output-groupfile': '--group-table-prefix',
        '--output_groupfile': '--group-table-prefix',
    },
}


def _deprecated_option_aliases(command):
    """Return deprecated long-option aliases and their canonical replacements."""
    command_parser = subparsers.choices.get(command)
    if command_parser is None:
        return dict()
    canonical_by_dest = dict()
    for action in command_parser._actions:
        if action.help == argparse.SUPPRESS:
            continue
        canonical = _canonical_long_option(action)
        if canonical is not None:
            canonical_by_dest[action.dest] = canonical
    aliases = dict()
    for action in command_parser._actions:
        canonical = canonical_by_dest.get(action.dest)
        if canonical is None:
            continue
        for option in action.option_strings:
            if option.startswith('--') and option != canonical:
                aliases[option] = canonical
    aliases.update(LEGACY_OPTION_CANONICAL_OVERRIDES.get(command, {}))
    return aliases


def _used_deprecated_option_aliases(argv, command):
    aliases = _deprecated_option_aliases(command)
    used = list()
    seen = set()
    for token in argv:
        if token == '--':
            break
        option = str(token).split('=', 1)[0]
        if option in aliases and option not in seen:
            used.append((option, aliases[option]))
            seen.add(option)
    return used


def _warn_deprecated_option_aliases(argv, command):
    for alias, canonical in _used_deprecated_option_aliases(argv, command):
        sys.stderr.write(
            "Warning: option '{}' is deprecated; use '{}' instead.\n".format(
                alias,
                canonical,
            )
        )


def main(argv=None):
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    args = parser.parse_args(raw_argv)
    if hasattr(args, 'handler'):
        if getattr(args, 'audit', None) == '-':
            raise ValueError("'--audit' requires a file path, not '-'.")
        def invoke_handler(parsed_args):
            _warn_deprecated_option_aliases(raw_argv, parsed_args.command)
            _validate_stdin_ownership(parsed_args)
            return parsed_args.handler(parsed_args)
        if getattr(args, 'audit', None):
            from nwkit.provenance import run_with_audit
            return run_with_audit(args=args, argv=raw_argv, handler=invoke_handler)
        return invoke_handler(args)
    parser.print_help()
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
