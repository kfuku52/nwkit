import argparse
import io
import json
import re
import sys
from argparse import Namespace

import pandas as pd
import pytest

from nwkit.annotate import ANNOTATION_REPORT_COLUMNS
from nwkit.asr import _parse_targets, asr_main
from nwkit.cladefreq import CLADEFREQ_COLUMNS, cladefreq_main
from nwkit.cli import _canonical_long_option, _deprecated_option_aliases, main, parser
from nwkit.consensus import _read_tree_weights
from nwkit.constrain import check_input_file, read_taxid_tsv
from nwkit.diff import DIFF_COLUMNS
from nwkit.image import extract_species_mapping, read_name_tsv as read_species_name_tsv
from nwkit.monophyly import MONOPHYLY_COLUMNS, monophyly_main
from nwkit.nwk2table import nwk2table_main
from nwkit.rename import read_name_tsv as read_rename_name_tsv
from nwkit.species_parser import read_species_map_tsv
from nwkit.table2nwk import table2nwk_main
from nwkit.transfer import REPORT_COLUMNS
from nwkit.util import assign_branch_ids, read_tip_table, read_tree
from tests.helpers import make_args


def _command_parsers():
    subparsers = next(
        action
        for action in parser._actions
        if isinstance(action, argparse._SubParsersAction)
    )
    return subparsers.choices


def test_every_visible_long_option_has_a_canonical_kebab_case_spelling():
    for command, command_parser in _command_parsers().items():
        for action in command_parser._actions:
            long_options = [
                option for option in action.option_strings if option.startswith('--')
            ]
            if not long_options or action.help == argparse.SUPPRESS:
                continue
            assert '_' not in long_options[0], (command, action.dest, long_options)
            expected = '--{}'.format(action.dest.replace('_', '-'))
            assert expected in long_options, (command, action.dest, long_options)
            if '_' in action.dest:
                assert '--{}'.format(action.dest) in long_options, (
                    command,
                    action.dest,
                    long_options,
                )


def test_help_hides_snake_case_compatibility_aliases(capsys):
    for command, command_parser in _command_parsers().items():
        assert re.search(r'--[A-Za-z0-9-]*_[A-Za-z0-9_-]*', command_parser.format_help()) is None, command

    with pytest.raises(SystemExit) as exc_info:
        main(['sample', '-h'])
    assert exc_info.value.code == 0
    help_text = capsys.readouterr().out
    assert '--allow-fewer' in help_text
    assert '--allow_fewer' not in help_text
    assert '--output-table' not in help_text
    assert '--output_table' not in help_text

    image_help = _command_parsers()['image'].format_help()
    assert '--species-name-tsv' in image_help
    assert '--name-tsv' not in image_help
    assert '--manifest-out' in image_help
    assert '--manifest PATH' not in image_help

    skim_help = _command_parsers()['skim'].format_help()
    assert '--group-table-prefix' in skim_help
    assert '--output-groupfile' not in skim_help

    dist_help = _command_parsers()['dist'].format_help()
    assert '--metric' in dist_help
    assert '--comparison' in dist_help
    assert '--dist' not in dist_help
    assert re.search(r'(^|[, ]\s)-d([, ]|$)', dist_help) is None


def test_every_noncanonical_long_option_is_registered_as_deprecated():
    for command, command_parser in _command_parsers().items():
        aliases = _deprecated_option_aliases(command)
        canonical_by_dest = {
            action.dest: _canonical_long_option(action)
            for action in command_parser._actions
            if action.help != argparse.SUPPRESS
            and _canonical_long_option(action) is not None
        }
        for action in command_parser._actions:
            canonical = canonical_by_dest.get(action.dest)
            for option in action.option_strings:
                if not option.startswith('--') or option == canonical:
                    continue
                if canonical is None:
                    assert option in aliases, (command, option)
                else:
                    assert aliases.get(option) == canonical, (
                        command,
                        option,
                        canonical,
                    )


def test_required_options_do_not_advertise_nonempty_defaults():
    for command, command_parser in _command_parsers().items():
        for action in command_parser._actions:
            if action.required:
                assert action.default in (None, argparse.SUPPRESS), (
                    command,
                    action.dest,
                    action.default,
                )


def test_non_newick_outputs_do_not_expose_outformat():
    commands = _command_parsers()
    for command in (
        'asr', 'cladefreq', 'diff', 'dist', 'draw', 'info', 'mcmctree',
        'monophyly', 'nwk2table', 'printlabel', 'validate',
    ):
        assert '--outformat' not in commands[command]._option_string_actions


def test_compatibility_aliases_resolve_to_canonical_destinations():
    image_args = parser.parse_args([
        'image', '--out_dir', 'images', '--name_tsv', 'names.tsv',
        '--manifest', 'manifest.tsv', '--attribution', 'ATTRIBUTION.md',
    ])
    assert image_args.out_dir == 'images'
    assert image_args.species_name_tsv == 'names.tsv'
    assert image_args.manifest_out == 'manifest.tsv'
    assert image_args.attribution_out == 'ATTRIBUTION.md'

    sample_args = parser.parse_args(['sample', '--output_table', 'sample.tsv'])
    assert sample_args.report == 'sample.tsv'

    mcmctree_args = parser.parse_args([
        'mcmctree', '--lower_tailProb', '0.1', '--upper_tail_prob', '0.2',
    ])
    assert mcmctree_args.lower_tail_prob == '0.1'
    assert mcmctree_args.upper_tail_prob == '0.2'
    assert parser.parse_args(['mcmctree', '--add-header']).add_header is True
    assert parser.parse_args(['mcmctree', '--add-header', 'yes']).add_header is True
    assert parser.parse_args(['mcmctree', '--add-header', 'no']).add_header is False

    assert _parse_targets('leaf,missing-leaf') == {'leaf', 'missing-leaf'}
    assert _parse_targets('tip,missing_tip') == {'leaf', 'missing-leaf'}

    dist_args = parser.parse_args([
        'dist', '--infile2', 'tree2.nwk', '--metric', 'rf,path-topological',
        '--metric', 'normalized-rf', '--comparison', 'unrooted',
    ])
    assert dist_args.metric == ['rf,path-topological', 'normalized-rf']
    assert dist_args.dist is None
    assert dist_args.comparison == 'unrooted'


@pytest.mark.parametrize('legacy_option', ['-d', '--dist'])
def test_dist_legacy_metric_option_warns_and_preserves_rf_output(
    tmp_path,
    capsys,
    legacy_option,
):
    tree1_path = tmp_path / 'tree1.nwk'
    tree2_path = tmp_path / 'tree2.nwk'
    tree1_path.write_text('((A:1,B:1):1,(C:1,D:1):1);')
    tree2_path.write_text('((A:1,C:1):1,(B:1,D:1):1);')

    main([
        'dist', '--infile', str(tree1_path), '--infile2', str(tree2_path),
        legacy_option, 'RF',
    ])

    captured = capsys.readouterr()
    assert captured.out == 'rf_dist\tmax_rf_dist\n4\t4\n'
    assert (
        "Warning: option '{}' is deprecated; use '--metric' instead.".format(
            legacy_option,
        )
    ) in captured.err


def test_deprecated_option_aliases_warn_with_canonical_replacements(tmp_path, capsys):
    image_aliases = _deprecated_option_aliases('image')
    assert image_aliases['--out_dir'] == '--out-dir'
    assert image_aliases['--name-tsv'] == '--species-name-tsv'
    assert image_aliases['--manifest'] == '--manifest-out'
    assert image_aliases['--attribution'] == '--attribution-out'
    assert _deprecated_option_aliases('sample')['--output-table'] == '--report'
    assert _deprecated_option_aliases('skim')['--output-groupfile'] == '--group-table-prefix'

    tree_path = tmp_path / 'tree.nwk'
    tree_path.write_text('(A:1,B:1);')
    main(['info', '-i', str(tree_path), '-f', '1', '--quoted_node_names=yes'])
    captured = capsys.readouterr()
    warning_lines = [line for line in captured.err.splitlines() if line.startswith('Warning:')]
    assert warning_lines == [
        "Warning: option '--quoted_node_names' is deprecated; "
        "use '--quoted-node-names' instead."
    ]

    main(['info', '-i', str(tree_path), '-f', '1', '--quoted-node-names', 'yes'])
    assert 'deprecated' not in capsys.readouterr().err


def test_deprecated_option_warning_is_recorded_by_audit(tmp_path, capsys):
    tree_path = tmp_path / 'tree.nwk'
    tree_path.write_text('(A:1,B:1);')
    audit_path = tmp_path / 'audit.jsonl'

    main([
        'info', '-i', str(tree_path), '-f', '1', '--quoted_node_names', 'yes',
        '--audit', str(audit_path),
    ])

    warning = next(
        line
        for line in capsys.readouterr().err.splitlines()
        if line.startswith('Warning:')
    )
    record = json.loads(audit_path.read_text())
    assert "option '--quoted_node_names' is deprecated" in warning
    assert warning in record['warnings']


def test_only_one_input_can_read_standard_input():
    with pytest.raises(ValueError, match='STDIN can be assigned to only one input'):
        main(['asr', '--trait', '-', '--state-column', 'state'])


def test_auxiliary_outputs_require_file_paths():
    with pytest.raises(ValueError, match="'--audit' requires a file path"):
        main(['info', '--audit', '-'])

    with pytest.raises(ValueError, match="'--report' must be a file path"):
        main(['sample', '--report', '-'])

    asr_args = make_args(
        trait='traits.tsv',
        state_column='state',
        model_out='-',
    )
    with pytest.raises(ValueError, match='Auxiliary outputs require file paths'):
        asr_main(asr_args)


def test_info_and_printlabel_honor_outfile(tmp_path, capsys):
    tree_path = tmp_path / 'tree.nwk'
    tree_path.write_text('((A:1,B:1)AB:1,C:2)ROOT;')
    info_path = tmp_path / 'info.txt'
    labels_path = tmp_path / 'labels.txt'

    main(['info', '-i', str(tree_path), '-f', '1', '-o', str(info_path)])
    main([
        'printlabel', '-i', str(tree_path), '-f', '1', '-o', str(labels_path),
        '--target', 'leaf',
    ])

    assert capsys.readouterr().out == ''
    assert 'Number of leaves: 3' in info_path.read_text()
    assert labels_path.read_text().splitlines() == ['A', 'B', 'C']


def test_sample_uses_shared_report_option(tmp_path):
    tree_path = tmp_path / 'tree.nwk'
    tree_path.write_text('((A:1,B:1):1,C:2);')
    outfile = tmp_path / 'sampled.nwk'
    report = tmp_path / 'sampled.tsv'

    main([
        'sample', '-i', str(tree_path), '-o', str(outfile), '-n', '1',
        '--report', str(report),
    ])

    assert outfile.exists()
    assert list(pd.read_csv(report, sep='\t').columns[:2]) == [
        'sample_order', 'leaf_name',
    ]


def test_skim_uses_explicit_group_table_prefix(tmp_path):
    tree_path = tmp_path / 'tree.nwk'
    tree_path.write_text('((A:1,B:1):1,(C:1,D:1):1);')
    trait_path = tmp_path / 'traits.tsv'
    trait_path.write_text('leaf_name\tgroup\nA\tx\nB\tx\nC\ty\nD\ty\n')
    outfile = tmp_path / 'skimmed.nwk'
    prefix = tmp_path / 'groups'

    main([
        'skim', '-i', str(tree_path), '-o', str(outfile),
        '--trait', str(trait_path), '--group-by', 'group',
        '--group-table-prefix', str(prefix),
    ])

    assert outfile.exists()
    assert (tmp_path / 'groups.all.tsv').exists()
    assert (tmp_path / 'groups.sampled.tsv').exists()


def test_tip_table_preserves_identifier_like_missing_markers(tmp_path):
    path = tmp_path / 'traits.tsv'
    path.write_text(
        'leaf_name\tgroup\n'
        '001\tNA\n'
        'NA\tx\n'
        'NaN\ty\n'
        'null\tnull\n'
    )

    table, table_only, tree_only = read_tip_table(
        str(path),
        tree_leaf_names=['001', 'NA', 'NaN', 'null'],
        unmatched='error',
    )

    assert table['leaf_name'].tolist() == ['001', 'NA', 'NaN', 'null']
    assert pd.isna(table.loc[0, 'group'])
    assert table.loc[3, 'group'] == 'null'
    assert table_only == []
    assert tree_only == []


def test_tip_table_unmatched_policy_is_shared(tmp_path, capsys):
    path = tmp_path / 'traits.tsv'
    path.write_text('leaf_name\tgroup\nA\tx\nZ\ty\n')

    _, table_only, tree_only = read_tip_table(
        str(path), tree_leaf_names=['A', 'B'], unmatched='warn'
    )
    assert table_only == ['Z']
    assert tree_only == ['B']
    warning = capsys.readouterr().err
    assert 'Rows in --trait not found in tree: Z' in warning
    assert 'Tree tips not found in --trait: B' in warning

    with pytest.raises(ValueError, match='--trait and tree tips differ'):
        read_tip_table(str(path), tree_leaf_names=['A', 'B'], unmatched='error')

    read_tip_table(str(path), tree_leaf_names=['A', 'B'], unmatched='ignore')
    assert capsys.readouterr().err == ''


def test_specialized_tsv_inputs_support_standard_input(monkeypatch, tmp_path):
    monkeypatch.setattr(sys, 'stdin', io.StringIO('leaf_name\tgroup\nA\tx\n'))
    tip_table, _, _ = read_tip_table('-', required_columns=('group',))
    assert tip_table.to_dict('records') == [{'leaf_name': 'A', 'group': 'x'}]

    monkeypatch.setattr(
        sys, 'stdin', io.StringIO('leaf_name\tspecies_name\nA\tApis mellifera\n')
    )
    assert read_species_name_tsv('-') == {'A': 'Apis mellifera'}

    monkeypatch.setattr(
        sys, 'stdin', io.StringIO('leaf_name\tspecies_label\nA\tApis_mellifera\n')
    )
    assert read_species_map_tsv('-')['A'].species_label == 'Apis_mellifera'

    monkeypatch.setattr(sys, 'stdin', io.StringIO('old_name\tnew_name\nA\tAlpha\n'))
    assert read_rename_name_tsv('-') == {'A': 'Alpha'}

    monkeypatch.setattr(sys, 'stdin', io.StringIO('tree_id\tweight\n1\t2\n2\t1\n'))
    assert _read_tree_weights('-', 2) == [2.0, 1.0]

    monkeypatch.setattr(sys, 'stdin', io.StringIO('leaf_name\ttaxid\nA\t7460\n'))
    taxids = read_taxid_tsv('-')
    assert taxids.to_dict('records') == [{'leaf_name': 'A', 'taxid': 7460}]

    monkeypatch.setattr(
        sys,
        'stdin',
        io.StringIO('branch_id\tparent\tname\n0\t-1\troot\n1\t0\tA\n2\t0\tB\n'),
    )
    outfile = tmp_path / 'stdin-table.nwk'
    table2nwk_main(make_args(infile='-', outfile=str(outfile), outformat='auto'))
    assert set(read_tree(str(outfile), '1', True, quiet=True).leaf_names()) == {'A', 'B'}


def test_constrain_validates_standard_input_only_once(monkeypatch):
    monkeypatch.setattr(sys, 'stdin', io.StringIO('A_a\nB_b\n'))
    labels, taxid_df = check_input_file(Namespace(
        species_list='-', taxid_tsv=None, backbone='user'
    ))
    assert labels == ['A_a', 'B_b']
    assert taxid_df is None

    monkeypatch.setattr(sys, 'stdin', io.StringIO('leaf_name\ttaxid\nA_a\t1\n'))
    labels, taxid_df = check_input_file(Namespace(
        species_list=None, taxid_tsv='-', backbone='ncbi'
    ))
    assert labels is None
    assert taxid_df.to_dict('records') == [{'leaf_name': 'A_a', 'taxid': 1}]


def test_image_species_mapping_honors_shared_species_map(tmp_path):
    tree_path = tmp_path / 'tree.nwk'
    tree_path.write_text("(sample_1:1,sample_2:1,' sample 3 ':1);")
    mapping_path = tmp_path / 'species.tsv'
    mapping_path.write_text(
        'leaf_name\tspecies_label\ttaxonomy_query\n'
        'sample_1\tApis_mellifera\tApis mellifera\n'
        'sample_2\tHomo_sapiens\tHomo sapiens\n'
        ' sample 3 \tMus_musculus\tMus musculus\n'
    )
    tree = read_tree(str(tree_path), '1', True, quiet=True)
    args = Namespace(
        species_parser='legacy',
        species_regex=r'^([^_]+_[^_]+)(?:_|$)',
        species_map_tsv=str(mapping_path),
    )

    leaf_to_species, unmatched = extract_species_mapping(tree, args=args)

    assert leaf_to_species == {
        'sample_1': 'Apis mellifera',
        'sample_2': 'Homo sapiens',
        ' sample 3 ': 'Mus musculus',
    }
    assert unmatched == []


def test_rename_mapping_preserves_exact_whitespace_in_node_names(tmp_path):
    tree_path = tmp_path / 'tree.nwk'
    tree_path.write_text("(' A ':1,B:1);")
    mapping_path = tmp_path / 'names.tsv'
    mapping_path.write_text('old_name\tnew_name\n A \t Alpha \n')
    outfile = tmp_path / 'renamed.nwk'

    main([
        'rename', '-i', str(tree_path), '-f', '1', '-o', str(outfile),
        '--name-tsv', str(mapping_path),
    ])

    tree = read_tree(str(outfile), '1', True, quiet=True)
    assert set(tree.leaf_names()) == {' Alpha ', 'B'}


def test_table2nwk_preserves_names_that_resemble_missing_values(tmp_path):
    table_path = tmp_path / 'tree.tsv'
    table_path.write_text(
        'branch_id\tparent\tname\tdist\tsupport\n'
        '0\t-1\tnull\t0\tNA\n'
        '1\t0\tNA\t1\tNaN\n'
        '2\t0\tNaN\t1\t\n'
    )
    outfile = tmp_path / 'tree.nwk'

    table2nwk_main(make_args(
        infile=str(table_path), outfile=str(outfile), outformat='auto'
    ))
    tree = read_tree(str(outfile), '1', True, quiet=True)

    assert tree.name == 'null'
    assert set(tree.leaf_names()) == {'NA', 'NaN'}


def test_empty_tabular_results_keep_stable_headers(tmp_path):
    trees_path = tmp_path / 'two_tips.nwk'
    trees_path.write_text('(A:1,B:1);\n')
    cladefreq_path = tmp_path / 'cladefreq.tsv'
    cladefreq_main(make_args(
        infile=str(trees_path),
        outfile=str(cladefreq_path),
        reference=None,
        reference_format='auto',
        weight_tsv=None,
        support_scale='percent',
        threads=1,
    ))
    assert tuple(pd.read_csv(cladefreq_path, sep='\t').columns) == CLADEFREQ_COLUMNS

    trait_path = tmp_path / 'empty_traits.tsv'
    trait_path.write_text('leaf_name\tgroup\n')
    monophyly_path = tmp_path / 'monophyly.tsv'
    monophyly_main(make_args(
        infile=str(trees_path),
        outfile=str(monophyly_path),
        trait=str(trait_path),
        group_by='group',
        unrooted=False,
        fail_on_non_monophyly=False,
        unmatched='warn',
    ))
    assert tuple(pd.read_csv(monophyly_path, sep='\t').columns) == MONOPHYLY_COLUMNS


def test_report_vocabulary_uses_branch_ids_node_classes_and_taxa():
    assert 'branch_id' in ANNOTATION_REPORT_COLUMNS
    assert 'node_id' not in ANNOTATION_REPORT_COLUMNS
    assert 'target_branch_id' in REPORT_COLUMNS
    assert 'target_node_id' not in REPORT_COLUMNS
    assert {'target_branch_id', 'source_branch_id'} <= set(DIFF_COLUMNS)
    assert {'target_node_id', 'source_node_id'}.isdisjoint(DIFF_COLUMNS)
    assert {
        'num_shared_taxa', 'num_target_only_taxa', 'num_source_only_taxa',
    } <= set(REPORT_COLUMNS)
    assert not any(column.endswith('_taxon_count') for column in REPORT_COLUMNS)


def test_assigning_report_ids_does_not_leak_newick_properties():
    tree = read_tree('((A:1,B:1):1,C:2);', '1', True, quiet=True)
    assert list(assign_branch_ids(tree).values()) == list(range(5))
    assert all('branch_id' not in node.props for node in tree.traverse())


def test_branch_ids_match_between_nwk2table_and_asr(tmp_path):
    tree_path = tmp_path / 'tree.nwk'
    tree_path.write_text('((A:1,B:1)AB:1,C:2)ROOT;')
    node_table_path = tmp_path / 'nodes.tsv'
    nwk2table_main(make_args(
        infile=str(tree_path),
        outfile=str(node_table_path),
        format='1',
        age=False,
        sister=False,
    ))
    trait_path = tmp_path / 'traits.tsv'
    trait_path.write_text('leaf_name\tstate\nA\tx\nB\tx\nC\ty\n')
    asr_path = tmp_path / 'asr.tsv'
    asr_main(make_args(
        infile=str(tree_path),
        outfile=str(asr_path),
        format='1',
        trait=str(trait_path),
        state_column='state',
        states='x,y',
        missing_values=None,
        model='ER',
        rate=0.1,
        rate_bounds=None,
        root_prior='equal',
        target='all',
        output='map',
    ))

    nodes = pd.read_csv(node_table_path, sep='\t').set_index('name')
    asr = pd.read_csv(asr_path, sep='\t').set_index('name')
    for name in ('ROOT', 'AB', 'A', 'B', 'C'):
        assert int(nodes.loc[name, 'branch_id']) == int(asr.loc[name, 'branch_id'])
    assert set(asr['node_class']) == {'root', 'intnode', 'leaf'}
