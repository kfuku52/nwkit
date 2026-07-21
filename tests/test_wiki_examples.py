import pandas as pd

from nwkit.asr import asr_main
from nwkit.cladefreq import cladefreq_main
from nwkit.collapse import collapse_main
from nwkit.consensus import consensus_main
from nwkit.monophyly import monophyly_main
from nwkit.nwk2table import nwk2table_main
from nwkit.rename import rename_main
from nwkit.table2nwk import table2nwk_main
from nwkit.validate import validate_main
from tests.helpers import make_args


def _write_tree_collection(tmp_path, trees, name='trees.nwk'):
    path = tmp_path / name
    path.write_text('\n'.join(trees) + '\n')
    return str(path)


class TestWikiExamples:
    def test_asr_example(self, tmp_path):
        infile = tmp_path / 'tree.nwk'
        infile.write_text('((A:1,B:1):1,(C:1,D:1):1);')
        trait = tmp_path / 'traits.tsv'
        pd.DataFrame(
            {
                'leaf_name': ['A', 'B', 'C', 'D'],
                'state': ['x', 'x', 'y', ''],
            }
        ).to_csv(trait, sep='\t', index=False)
        outfile = tmp_path / 'asr.tsv'
        args = make_args(
            infile=str(infile),
            outfile=str(outfile),
            trait=str(trait),
            state_column='state',
            states='x,y',
            missing_values=None,
            model='ER',
            rate=0.1,
            rate_bounds=None,
            root_prior='equal',
            target='intnode,missing_tip',
            output='map',
        )
        asr_main(args)
        table = pd.read_csv(outfile, sep='\t')
        assert {'branch_id', 'node_class', 'name', 'state', 'probability'}.issubset(table.columns)
        assert len(table.index) == 4
        assert table.loc[table['name'] == 'D', 'state'].iloc[0] == 'y'

    def test_collapse_example(self, tmp_path):
        infile = tmp_path / 'tree.nwk'
        infile.write_text('((A:1,B:1)40:1,(C:1,D:1)90:1);')
        outfile = tmp_path / 'collapsed.nwk'
        args = make_args(
            infile=str(infile),
            outfile=str(outfile),
            format='0',
            outformat='auto',
            min_support=70.0,
            max_dist=None,
            preserve_branch_length=True,
        )
        collapse_main(args)
        assert outfile.read_text().strip() == '((C:1,D:1)90:1,A:2,B:2);'

    def test_cladefreq_example(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            [
                '((A:1,B:1):1,(C:1,D:1):1);',
                '((A:1,B:1):1,(C:1,D:1):1);',
                '((A:1,C:1):1,(B:1,D:1):1);',
            ],
        )
        outfile = tmp_path / 'cladefreq.tsv'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            reference=None,
            reference_format='auto',
            weight_tsv=None,
            support_scale='percent',
        )
        cladefreq_main(args)
        table = pd.read_csv(outfile, sep='\t')
        table = table[['descendant_taxa', 'num_taxa', 'weight_sum', 'frequency']]
        expected = {
            'A,B': (2, 2.0, 200.0 / 3.0),
            'C,D': (2, 2.0, 200.0 / 3.0),
            'A,C': (2, 1.0, 100.0 / 3.0),
            'B,D': (2, 1.0, 100.0 / 3.0),
        }
        assert set(table['descendant_taxa']) == set(expected.keys())
        for _, row in table.iterrows():
            num_taxa, weight_sum, frequency = expected[row['descendant_taxa']]
            assert int(row['num_taxa']) == num_taxa
            assert abs(float(row['weight_sum']) - weight_sum) < 1e-9
            assert abs(float(row['frequency']) - frequency) < 1e-9

    def test_consensus_example(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            [
                '((A:1,B:1):1,(C:1,D:1):1);',
                '((A:1,B:1):1,(C:1,D:1):1);',
                '((A:1,C:1):1,(B:1,D:1):1);',
            ],
        )
        outfile = tmp_path / 'consensus.nwk'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            min_freq=0.5,
            reference=None,
            reference_format='auto',
            support_scale='percent',
            method='greedy',
            branch_length='none',
            weight_tsv=None,
        )
        consensus_main(args)
        assert outfile.read_text().strip() == '((A,B)66.6667,(C,D)66.6667);'

    def test_monophyly_example(self, tmp_path):
        infile = tmp_path / 'tree.nwk'
        infile.write_text('(((A:1,B:1):1,C:1):1,(D:1,E:1):1);')
        trait_path = tmp_path / 'traits.tsv'
        pd.DataFrame(
            {
                'leaf_name': ['A', 'B', 'C', 'D', 'E'],
                'group': ['x', 'x', 'y', 'z', 'y'],
            }
        ).to_csv(trait_path, sep='\t', index=False)
        outfile = tmp_path / 'monophyly.tsv'
        args = make_args(
            infile=str(infile),
            outfile=str(outfile),
            trait=str(trait_path),
            group_by='group',
            unrooted=False,
        )
        monophyly_main(args)
        table = pd.read_csv(outfile, sep='\t')
        by_group = table.set_index('group')
        assert by_group.loc['x', 'status'] == 'monophyletic'
        assert bool(by_group.loc['x', 'is_monophyletic']) is True
        assert by_group.loc['x', 'target_taxa'] == 'A,B'
        assert by_group.loc['y', 'status'] == 'polyphyletic'
        assert bool(by_group.loc['y', 'is_monophyletic']) is False
        assert by_group.loc['y', 'target_taxa'] == 'C,E'
        assert by_group.loc['y', 'intruder_taxa'] == 'A,B,D'
        assert by_group.loc['z', 'status'] == 'monophyletic'

    def test_rename_example(self, tmp_path):
        infile = tmp_path / 'tree.nwk'
        infile.write_text('((A:1,B:1):1,(C:1,D:1):1);')
        name_tsv = tmp_path / 'names.tsv'
        pd.DataFrame(
            {
                'old_name': ['A', 'C'],
                'new_name': ['Alpha', 'Charlie'],
            }
        ).to_csv(name_tsv, sep='\t', index=False)
        outfile = tmp_path / 'renamed.nwk'
        args = make_args(
            infile=str(infile),
            outfile=str(outfile),
            name_tsv=str(name_tsv),
            target='leaf',
            require_all_old_names=True,
            check_leaf_uniqueness=True,
        )
        rename_main(args)
        assert outfile.read_text().strip() == '((Alpha:1,B:1):1,(Charlie:1,D:1):1);'

    def test_validate_example(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            [
                '((A:-1,A:1):1,B:1);',
                '((A:1,B:1):1,(C:1,E:1):1);',
            ],
        )
        outfile = tmp_path / 'validate.tsv'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            require_rooted=False,
            require_ultrametric=False,
            require_same_leaf_set=True,
            require_same_rooting=False,
            require_binary=False,
            require_all_support=False,
            require_unambiguous_format=False,
            require_unquoted_names=False,
            check_species=False,
            fail_on_issue=False,
        )
        validate_main(args)
        table = pd.read_csv(outfile, sep='\t')
        first_row = table.iloc[0]
        second_row = table.iloc[1]
        assert first_row['status'] == 'invalid'
        assert int(first_row['num_duplicate_leaf_names']) == 1
        assert int(first_row['num_negative_branch_nodes']) == 1
        assert bool(first_row['leaf_set_matches_first']) is True
        assert first_row['issues'] == 'duplicate_leaf_names,negative_branch_length'
        assert second_row['status'] == 'invalid'
        assert int(second_row['num_duplicate_leaf_names']) == 0
        assert int(second_row['num_negative_branch_nodes']) == 0
        assert bool(second_row['leaf_set_matches_first']) is False
        assert second_row['issues'] == 'leaf_set_mismatch'

    def test_nwk2table_example(self, tmp_path):
        infile = tmp_path / 'tree.nwk'
        infile.write_text('((A:1,B:2)80:3,C:4);')
        outfile = tmp_path / 'tree.tsv'
        args = make_args(
            infile=str(infile),
            outfile=str(outfile),
            format='0',
            age=False,
            sister=False,
        )
        nwk2table_main(args)
        table = pd.read_csv(outfile, sep='\t')
        assert table['branch_id'].tolist() == [0, 1, 2, 3, 4]
        assert table['parent'].tolist() == [-1, 0, 0, 1, 1]
        assert pd.isna(table.loc[0, 'name'])
        assert pd.isna(table.loc[0, 'dist'])
        assert pd.isna(table.loc[0, 'support'])
        assert float(table.loc[1, 'dist']) == 3.0
        assert float(table.loc[1, 'support']) == 80.0
        assert table.loc[2, 'name'] == 'C'
        assert float(table.loc[2, 'dist']) == 4.0
        assert table.loc[3, 'name'] == 'A'
        assert float(table.loc[3, 'dist']) == 1.0
        assert table.loc[4, 'name'] == 'B'
        assert float(table.loc[4, 'dist']) == 2.0

    def test_table2nwk_example(self, tmp_path):
        infile = tmp_path / 'tree.tsv'
        pd.DataFrame(
            [
                {'branch_id': 0, 'parent': -1, 'name': '', 'dist': '', 'support': ''},
                {'branch_id': 1, 'parent': 0, 'name': '', 'dist': 3, 'support': 80},
                {'branch_id': 2, 'parent': 0, 'name': 'C', 'dist': 4, 'support': ''},
                {'branch_id': 3, 'parent': 1, 'name': 'A', 'dist': 1, 'support': ''},
                {'branch_id': 4, 'parent': 1, 'name': 'B', 'dist': 2, 'support': ''},
            ]
        ).to_csv(infile, sep='\t', index=False)
        outfile = tmp_path / 'tree.nwk'
        args = make_args(
            infile=str(infile),
            outfile=str(outfile),
            outformat='auto',
        )
        table2nwk_main(args)
        assert outfile.read_text().strip() == '((A:1,B:2)80:3,C:4);'
