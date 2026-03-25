import pandas as pd
import pytest

from nwkit.consensus import consensus_main
from nwkit.util import read_tree
from tests.helpers import make_args


def _write_tree_collection(tmp_path, trees, name='trees.nwk'):
    path = tmp_path / name
    path.write_text('\n'.join(trees) + '\n')
    return str(path)


class TestConsensusMain:
    def test_strict_consensus_removes_conflicting_clades(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            [
                '((A:1,B:1):1,(C:1,D:1):1);',
                '((A:1,C:1):1,(B:1,D:1):1);',
            ],
            name='strict.nwk',
        )
        outfile = str(tmp_path / 'strict_consensus.nwk')
        args = make_args(
            infile=infile,
            outfile=outfile,
            min_freq=0.5,
            reference=None,
            reference_format='auto',
            support_scale='percent',
            method='strict',
            branch_length='none',
            weight_tsv=None,
        )
        consensus_main(args)
        tree = read_tree(outfile, format='auto', quoted_node_names=True, quiet=True)
        assert {frozenset(child.leaf_names()) for child in tree.get_children()} == {
            frozenset({'A'}),
            frozenset({'B'}),
            frozenset({'C'}),
            frozenset({'D'}),
        }

    def test_builds_majority_consensus_tree(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            [
                '((A:1,B:1):1,(C:1,D:1):1);',
                '((A:1,B:1):1,(C:1,D:1):1);',
                '((A:1,C:1):1,(B:1,D:1):1);',
            ],
        )
        outfile = str(tmp_path / 'consensus.nwk')
        args = make_args(
            infile=infile,
            outfile=outfile,
            min_freq=0.5,
            reference=None,
            reference_format='auto',
            support_scale='percent',
            method='greedy',
            branch_length='none',
            weight_tsv=None,
        )
        consensus_main(args)
        tree = read_tree(outfile, format='auto', quoted_node_names=True, quiet=True)
        root_children = {frozenset(child.leaf_names()) for child in tree.get_children()}
        assert root_children == {frozenset({'A', 'B'}), frozenset({'C', 'D'})}
        assert abs(tree.common_ancestor(['A', 'B']).support - (200.0 / 3.0)) < 1e-4
        assert abs(tree.common_ancestor(['C', 'D']).support - (200.0 / 3.0)) < 1e-4

    def test_transfers_consensus_support_to_reference_tree(self, tmp_path, tmp_nwk):
        infile = _write_tree_collection(
            tmp_path,
            [
                '((A:1,B:1):1,(C:1,D:1):1);',
                '((A:1,B:1):1,(C:1,D:1):1);',
                '((A:1,C:1):1,(B:1,D:1):1);',
            ],
        )
        reference = tmp_nwk('((A:1,C:1):1,(B:1,D:1):1);', 'reference.nwk')
        outfile = str(tmp_path / 'reference_supported.nwk')
        args = make_args(
            infile=infile,
            outfile=outfile,
            min_freq=0.5,
            reference=reference,
            reference_format='auto',
            support_scale='percent',
            method='greedy',
            branch_length='none',
            weight_tsv=None,
        )
        consensus_main(args)
        tree = read_tree(outfile, format='auto', quoted_node_names=True, quiet=True)
        assert abs(tree.common_ancestor(['A', 'C']).support - (100.0 / 3.0)) < 1e-4
        assert abs(tree.common_ancestor(['B', 'D']).support - (100.0 / 3.0)) < 1e-4

    def test_rejects_mismatched_leaf_sets(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            [
                '((A:1,B:1):1,(C:1,D:1):1);',
                '((A:1,B:1):1,(C:1,E:1):1);',
            ],
        )
        args = make_args(
            infile=infile,
            outfile='-',
            min_freq=0.5,
            reference=None,
            reference_format='auto',
            support_scale='percent',
            method='greedy',
            branch_length='none',
            weight_tsv=None,
        )
        with pytest.raises(ValueError, match='Leaf labels must be identical'):
            consensus_main(args)

    def test_weighted_consensus_prefers_heavier_tree(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            [
                '((A:1,B:1):1,(C:1,D:1):1);',
                '((A:1,C:1):1,(B:1,D:1):1);',
            ],
            name='weighted.nwk',
        )
        weight_tsv = tmp_path / 'weights.tsv'
        pd.DataFrame({'weight': [3.0, 1.0]}).to_csv(weight_tsv, sep='\t', index=False)
        outfile = str(tmp_path / 'weighted_consensus.nwk')
        args = make_args(
            infile=infile,
            outfile=outfile,
            min_freq=0.5,
            reference=None,
            reference_format='auto',
            support_scale='percent',
            method='greedy',
            branch_length='none',
            weight_tsv=str(weight_tsv),
        )
        consensus_main(args)
        tree = read_tree(outfile, format='auto', quoted_node_names=True, quiet=True)
        assert {frozenset(child.leaf_names()) for child in tree.get_children()} == {
            frozenset({'A', 'B'}),
            frozenset({'C', 'D'}),
        }
        assert abs(tree.common_ancestor(['A', 'B']).support - 75.0) < 1e-6

    def test_consensus_can_average_branch_lengths(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            [
                '((A:1,B:2):1,(C:3,D:4):5);',
                '((A:5,B:6):3,(C:7,D:8):9);',
            ],
            name='branch_lengths.nwk',
        )
        outfile = str(tmp_path / 'branch_length_consensus.nwk')
        args = make_args(
            infile=infile,
            outfile=outfile,
            min_freq=0.5,
            reference=None,
            reference_format='auto',
            support_scale='percent',
            method='majority',
            branch_length='mean',
            weight_tsv=None,
        )
        consensus_main(args)
        tree = read_tree(outfile, format='auto', quoted_node_names=True, quiet=True)
        assert abs(tree.common_ancestor(['A', 'B']).dist - 2.0) < 1e-6
        assert abs(next(tree.search_nodes(name='A')).dist - 3.0) < 1e-6
