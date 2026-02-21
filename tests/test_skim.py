import numpy
import pandas
import pytest
from argparse import Namespace
from ete4 import Tree

from nwkit.skim import (
    read_trait,
    mark_traits_to_nodes,
    add_group_ids,
    add_contrastive_clade_ids,
    sample_from_groups,
    skim_main,
)
from tests.helpers import make_args


def make_skim_args(**kwargs):
    defaults = {
        'trait': None,
        'group_by': None,
        'retain_per_clade': 1,
        'prioritize_non_missing': True,
        'filter_by': None,
        'filter_mode': 'ascending',
        'only_contrastive_clades': False,
        'output_groupfile': False,
    }
    defaults.update(kwargs)
    return Namespace(**defaults)


class TestReadTrait:
    def test_read_trait_fills_missing_and_drops_unknown(self, tmp_path):
        tree = Tree('((A:1,B:1):1,C:1);', parser=1)
        trait_path = tmp_path / 'trait.tsv'
        pandas.DataFrame(
            {'leaf_name': ['A', 'B', 'D'], 'trait': ['x', 'x', 'z']}
        ).to_csv(trait_path, sep='\t', index=False)
        args = make_skim_args(trait=str(trait_path))
        out = read_trait(args, tree)
        assert set(out['leaf_name']) == {'A', 'B', 'C'}
        assert int((out['leaf_name'] == 'C').sum()) == 1

    def test_duplicate_leaf_names_raise(self, tmp_path):
        tree = Tree('((A:1,B:1):1,C:1);', parser=1)
        trait_path = tmp_path / 'trait.tsv'
        pandas.DataFrame(
            {'leaf_name': ['A', 'A', 'B'], 'trait': ['x', 'y', 'x']}
        ).to_csv(trait_path, sep='\t', index=False)
        args = make_skim_args(trait=str(trait_path))
        with pytest.raises(ValueError, match="Duplicated 'leaf_name'"):
            read_trait(args, tree)

    def test_duplicate_non_string_leaf_names_raise_clear_error(self, tmp_path):
        tree = Tree('(A:1,:1,:1);', parser=0)
        trait_path = tmp_path / 'trait.tsv'
        pandas.DataFrame(
            {'leaf_name': [float('nan'), float('nan'), 'A'], 'trait': ['x', 'y', 'z']}
        ).to_csv(trait_path, sep='\t', index=False)
        args = make_skim_args(trait=str(trait_path))
        with pytest.raises(ValueError, match="Duplicated 'leaf_name'"):
            read_trait(args, tree)

    def test_missing_leaf_name_column_raises(self, tmp_path):
        tree = Tree('((A:1,B:1):1,C:1);', parser=1)
        trait_path = tmp_path / 'trait.tsv'
        pandas.DataFrame(
            {'species': ['A', 'B'], 'trait': ['x', 'y']}
        ).to_csv(trait_path, sep='\t', index=False)
        args = make_skim_args(trait=str(trait_path))
        with pytest.raises(ValueError, match="Column 'leaf_name'"):
            read_trait(args, tree)


class TestGrouping:
    def test_add_group_ids_for_two_homogeneous_clades(self):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        args = make_skim_args(group_by='trait')
        trait_df = pandas.DataFrame(
            {'leaf_name': ['A', 'B', 'C', 'D'], 'trait': ['x', 'x', 'y', 'y']}
        )
        marked_tree = mark_traits_to_nodes(tree, trait_df, args)
        grouped = add_group_ids(trait_df.copy(), marked_tree)
        leaf_to_group = dict(zip(grouped['leaf_name'], grouped['group']))
        assert leaf_to_group['A'] == leaf_to_group['B']
        assert leaf_to_group['C'] == leaf_to_group['D']
        assert leaf_to_group['A'] != leaf_to_group['C']

    def test_add_contrastive_clade_ids_marks_minimal_mixed_clades(self):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        args = make_skim_args(group_by='trait')
        trait_df = pandas.DataFrame(
            {'leaf_name': ['A', 'B', 'C', 'D'], 'trait': ['x', 'y', 'z', 'z']}
        )
        marked_tree = mark_traits_to_nodes(tree, trait_df, args)
        contrastive = add_contrastive_clade_ids(trait_df.copy(), marked_tree)
        marked = contrastive[~contrastive['contrastive_clade'].isna()]
        assert set(marked['leaf_name']) == {'A', 'B'}
        assert len(set(marked['contrastive_clade'])) == 1

    def test_missing_group_by_column_raises(self):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        args = make_skim_args(group_by='trait')
        trait_df = pandas.DataFrame({'leaf_name': ['A', 'B', 'C', 'D']})
        with pytest.raises(ValueError, match='group_by'):
            mark_traits_to_nodes(tree, trait_df, args)

    def test_mark_traits_to_nodes_handles_unnamed_leaves_without_keyerror(self):
        tree = Tree('(A:1,:1,:1);', parser=0)
        args = make_skim_args(group_by=None)
        trait_df = pandas.DataFrame({'leaf_name': list(tree.leaf_names())})
        marked = mark_traits_to_nodes(tree, trait_df, args)
        assert set(marked.leaf_names()) == {'A', None}


class TestSampling:
    def test_sample_from_groups_prioritizes_non_missing(self):
        trait_df = pandas.DataFrame(
            {
                'leaf_name': ['A', 'B', 'C', 'D'],
                'group': [1, 1, 2, 2],
                'trait': ['x', numpy.nan, 'y', numpy.nan],
                'score': [10, 0, 20, 0],
            }
        )
        args = make_skim_args(
            group_by='trait',
            retain_per_clade=1,
            prioritize_non_missing=True,
            filter_by=None,
        )
        sampled = sample_from_groups(trait_df, args)
        assert set(sampled['leaf_name']) == {'A', 'C'}

    def test_missing_filter_by_column_raises(self):
        trait_df = pandas.DataFrame(
            {
                'leaf_name': ['A', 'B'],
                'group': [1, 1],
            }
        )
        args = make_skim_args(filter_by='score')
        with pytest.raises(ValueError, match='filter_by'):
            sample_from_groups(trait_df, args)


class TestSkimMainValidation:
    def test_only_contrastive_clades_with_no_contrastive_raises(self, tmp_path):
        nwk_path = tmp_path / 'tree.nwk'
        nwk_path.write_text('((A:1,B:1):1,(C:1,D:1):1);')
        trait_path = tmp_path / 'trait.tsv'
        pandas.DataFrame(
            {'leaf_name': ['A', 'B', 'C', 'D'], 'trait': ['x', 'x', 'x', 'x']}
        ).to_csv(trait_path, sep='\t', index=False)
        out_tree = tmp_path / 'out.nwk'
        args = make_args(
            infile=str(nwk_path),
            outfile=str(out_tree),
            format='auto',
            outformat='auto',
            quoted_node_names=True,
            trait=str(trait_path),
            group_by='trait',
            retain_per_clade=1,
            prioritize_non_missing=True,
            filter_by=None,
            filter_mode='ascending',
            only_contrastive_clades=True,
            output_groupfile=False,
        )
        with pytest.raises(ValueError, match='No leaves were selected'):
            skim_main(args)

    def test_retain_per_clade_must_be_positive(self, tmp_path):
        nwk_path = tmp_path / 'tree.nwk'
        nwk_path.write_text('((A:1,B:1):1,(C:1,D:1):1);')
        trait_path = tmp_path / 'trait.tsv'
        pandas.DataFrame(
            {'leaf_name': ['A', 'B', 'C', 'D'], 'trait': ['x', 'x', 'y', 'y']}
        ).to_csv(trait_path, sep='\t', index=False)
        out_tree = tmp_path / 'out.nwk'
        args = make_args(
            infile=str(nwk_path),
            outfile=str(out_tree),
            format='auto',
            outformat='auto',
            quoted_node_names=True,
            trait=str(trait_path),
            group_by='trait',
            retain_per_clade=0,
            prioritize_non_missing=True,
            filter_by=None,
            filter_mode='ascending',
            only_contrastive_clades=False,
            output_groupfile=False,
        )
        with pytest.raises(ValueError, match='retain_per_clade'):
            skim_main(args)

    def test_output_groupfile_requires_file_out(self, tmp_path):
        nwk_path = tmp_path / 'tree.nwk'
        nwk_path.write_text('((A:1,B:1):1,(C:1,D:1):1);')
        trait_path = tmp_path / 'trait.tsv'
        pandas.DataFrame(
            {'leaf_name': ['A', 'B', 'C', 'D'], 'trait': ['x', 'x', 'y', 'y']}
        ).to_csv(trait_path, sep='\t', index=False)
        args = make_args(
            infile=str(nwk_path),
            outfile='-',
            format='auto',
            outformat='auto',
            quoted_node_names=True,
            trait=str(trait_path),
            group_by='trait',
            retain_per_clade=1,
            prioritize_non_missing=True,
            filter_by=None,
            filter_mode='ascending',
            only_contrastive_clades=False,
            output_groupfile=True,
        )
        with pytest.raises(ValueError, match='output_groupfile'):
            skim_main(args)

    def test_empty_leaf_labels_raise_clear_error(self, tmp_path):
        nwk_path = tmp_path / 'tree.nwk'
        nwk_path.write_text('(A:1,:1,B:1);')
        out_tree = tmp_path / 'out.nwk'
        args = make_args(
            infile=str(nwk_path),
            outfile=str(out_tree),
            format='auto',
            outformat='auto',
            quoted_node_names=True,
            trait=None,
            group_by=None,
            retain_per_clade=1,
            prioritize_non_missing=True,
            filter_by=None,
            filter_mode='ascending',
            only_contrastive_clades=False,
            output_groupfile=False,
        )
        with pytest.raises(ValueError, match='Empty leaf labels'):
            skim_main(args)

    def test_duplicate_leaf_labels_raise_clear_error(self, tmp_path):
        nwk_path = tmp_path / 'tree.nwk'
        nwk_path.write_text('((A:1,A:1):1,B:1);')
        out_tree = tmp_path / 'out.nwk'
        args = make_args(
            infile=str(nwk_path),
            outfile=str(out_tree),
            format='auto',
            outformat='auto',
            quoted_node_names=True,
            trait=None,
            group_by=None,
            retain_per_clade=1,
            prioritize_non_missing=True,
            filter_by=None,
            filter_mode='ascending',
            only_contrastive_clades=False,
            output_groupfile=False,
        )
        with pytest.raises(ValueError, match='Duplicated leaf labels'):
            skim_main(args)
