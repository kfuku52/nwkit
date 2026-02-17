import numpy
import pandas
from argparse import Namespace
from ete4 import Tree

from nwkit.skim import (
    read_trait,
    mark_traits_to_nodes,
    add_group_ids,
    add_contrastive_clade_ids,
    sample_from_groups,
)


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
