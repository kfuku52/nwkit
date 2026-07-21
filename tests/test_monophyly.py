import pandas as pd
import pytest

from nwkit.monophyly import monophyly_main
from tests.helpers import make_args


class TestMonophylyMain:
    def test_trait_groups_report_monophyly(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('(((A:1,B:1):1,C:1):1,(D:1,E:1):1);', 'tree.nwk')
        trait_path = tmp_path / 'traits.tsv'
        pd.DataFrame(
            {
                'leaf_name': ['A', 'B', 'C', 'D', 'E'],
                'group': ['x', 'x', 'y', 'z', 'y'],
            }
        ).to_csv(trait_path, sep='\t', index=False)
        outfile = tmp_path / 'monophyly.tsv'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            trait=str(trait_path),
            group_by='group',
            unrooted=False,
            fail_on_non_monophyly=False,
        )
        monophyly_main(args)
        table = pd.read_csv(outfile, sep='\t')
        row_x = table.loc[table['group'] == 'x'].iloc[0]
        row_y = table.loc[table['group'] == 'y'].iloc[0]
        assert bool(row_x['is_monophyletic']) is True
        assert bool(row_y['is_monophyletic']) is False
        assert row_y['num_intruder_taxa'] > 0

    def test_species_mode_groups_duplicate_species(self, tmp_nwk, tmp_path):
        infile = tmp_nwk(
            '(((Homo_sapiens_gene1:1,Homo_sapiens_gene2:1):1,Pan_troglodytes_gene1:1):1,Mus_musculus_gene1:1);',
            'tree.nwk',
        )
        outfile = tmp_path / 'monophyly.tsv'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            trait=None,
            group_by=None,
            unrooted=False,
            fail_on_non_monophyly=False,
        )
        monophyly_main(args)
        table = pd.read_csv(outfile, sep='\t')
        homo_row = table.loc[table['group'] == 'Homo_sapiens'].iloc[0]
        assert bool(homo_row['is_monophyletic']) is True
        assert homo_row['num_target_taxa'] == 2

    def test_fail_on_non_monophyly_raises(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('(((A:1,B:1):1,C:1):1,(D:1,E:1):1);', 'tree.nwk')
        trait_path = tmp_path / 'traits.tsv'
        pd.DataFrame(
            {
                'leaf_name': ['A', 'B', 'C', 'D', 'E'],
                'group': ['x', 'x', 'y', 'z', 'y'],
            }
        ).to_csv(trait_path, sep='\t', index=False)
        args = make_args(
            infile=infile,
            outfile='-',
            trait=str(trait_path),
            group_by='group',
            unrooted=False,
            fail_on_non_monophyly=True,
        )
        with pytest.raises(ValueError, match='Non-monophyletic group'):
            monophyly_main(args)

    def test_trait_mode_rejects_missing_leaf_names(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,B:1):1,C:1);', 'tree.nwk')
        trait_path = tmp_path / 'traits.tsv'
        pd.DataFrame(
            {
                'leaf_name': ['A', 'Z'],
                'group': ['x', 'x'],
            }
        ).to_csv(trait_path, sep='\t', index=False)
        args = make_args(
            infile=infile,
            outfile='-',
            trait=str(trait_path),
            group_by='group',
            unrooted=False,
            fail_on_non_monophyly=False,
            unmatched='error',
        )
        with pytest.raises(ValueError, match='--trait and tree tips differ'):
            monophyly_main(args)

    def test_trait_mode_rejects_duplicate_tree_leaf_names(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,A:1):1,B:1);', 'tree.nwk')
        trait_path = tmp_path / 'traits.tsv'
        pd.DataFrame(
            {
                'leaf_name': ['A', 'B'],
                'group': ['x', 'y'],
            }
        ).to_csv(trait_path, sep='\t', index=False)
        args = make_args(
            infile=infile,
            outfile='-',
            trait=str(trait_path),
            group_by='group',
            unrooted=False,
            fail_on_non_monophyly=False,
        )
        with pytest.raises(ValueError, match='Duplicated leaf labels'):
            monophyly_main(args)

    def test_trait_mode_accepts_leading_zero_leaf_names(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((001:1,002:1):1,003:1);', 'tree.nwk')
        trait_path = tmp_path / 'traits.tsv'
        pd.DataFrame(
            {
                'leaf_name': ['001', '002', '003'],
                'group': ['x', 'x', 'y'],
            }
        ).to_csv(trait_path, sep='\t', index=False)
        outfile = tmp_path / 'monophyly.tsv'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            trait=str(trait_path),
            group_by='group',
            unrooted=False,
            fail_on_non_monophyly=False,
        )
        monophyly_main(args)
        table = pd.read_csv(outfile, sep='\t')
        row_x = table.loc[table['group'] == 'x'].iloc[0]
        assert bool(row_x['is_monophyletic']) is True
        assert row_x['target_taxa'] == '001,002'

    def test_trait_mode_accepts_na_literal_leaf_names(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((NA:1,B:1):1,C:1);', 'tree.nwk')
        trait_path = tmp_path / 'traits.tsv'
        pd.DataFrame(
            {
                'leaf_name': ['NA', 'B', 'C'],
                'group': ['x', 'x', 'y'],
            }
        ).to_csv(trait_path, sep='\t', index=False)
        outfile = tmp_path / 'monophyly.tsv'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            trait=str(trait_path),
            group_by='group',
            unrooted=False,
            fail_on_non_monophyly=False,
        )
        monophyly_main(args)
        table = pd.read_csv(outfile, sep='\t')
        row_x = table.loc[table['group'] == 'x'].iloc[0]
        assert bool(row_x['is_monophyletic']) is True
        assert row_x['target_taxa'] == 'B,NA'
