import os
import pytest
from ete3 import TreeNode

from nwkit.subtree import subtree_main
from nwkit.util import read_tree
from tests.helpers import make_args, DATA_DIR


class TestSubtreeMain:
    def test_extract_subtree_left_right(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('(((a:1,b:1):1,c:1):1,((d:1,e:1),f:1):1):0;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            left_leaf='a', right_leaf='c', leaves=None,
            orthogroup=False,
        )
        subtree_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        leaf_names = set(tree.get_leaf_names())
        assert leaf_names == {'a', 'b', 'c'}

    def test_extract_subtree_with_leaves(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('(((a:1,b:1):1,c:1):1,((d:1,e:1),f:1):1):0;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            left_leaf=None, right_leaf=None, leaves='a,b,c',
            orthogroup=False,
        )
        subtree_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.get_leaf_names()) == {'a', 'b', 'c'}

    def test_extract_subtree_two_leaves(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('(((a:1,b:1):1,c:1):1,((d:1,e:1),f:1):1):0;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            left_leaf=None, right_leaf=None, leaves='d,e',
            orthogroup=False,
        )
        subtree_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.get_leaf_names()) == {'d', 'e'}

    def test_single_leaf(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('(((a:1,b:1):1,c:1):1,((d:1,e:1),f:1):1):0;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            left_leaf=None, right_leaf=None, leaves='a',
            orthogroup=False,
        )
        subtree_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.get_leaf_names()) == {'a'}

    def test_leaf_not_found_raises(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((a:1,b:1):1,c:1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            left_leaf=None, right_leaf=None, leaves='z',
            orthogroup=False,
        )
        with pytest.raises(Exception, match='Specified leaf not found'):
            subtree_main(args)

    def test_with_data_file(self, tmp_outfile):
        infile = os.path.join(DATA_DIR, 'subtree2', 'input.nwk')
        if not os.path.exists(infile):
            pytest.skip('Test data not found')
        args = make_args(
            infile=infile, outfile=tmp_outfile,
            left_leaf='a', right_leaf='c', leaves=None,
            orthogroup=False,
        )
        subtree_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.get_leaf_names()) == {'a', 'b', 'c'}

    def test_orthogroup_mode(self, tmp_nwk, tmp_outfile):
        nwk = '((Homo_sapiens_G1:1,(Homo_sapiens_G2:1,Mus_musculus_G1:1):1):1,(Danio_rerio_G1:1,Xenopus_laevis_G1:1):1);'
        path = tmp_nwk(nwk)
        args = make_args(
            infile=path, outfile=tmp_outfile,
            left_leaf=None, right_leaf=None, leaves='Homo_sapiens_G1',
            orthogroup=True, dup_conf_score_threshold=0,
        )
        subtree_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert len(tree.get_leaf_names()) >= 1

    def test_wiki_exact_example(self, tmp_nwk, tmp_outfile):
        """Wiki example: nwkit subtree --left_leaf a --right_leaf c

        Input:  (((a:1,b:1):1,c:1):1,((d:1,e:1),f:1):1):0;
        Output: ((a:1,b:1):1,c:1):1;
        """
        path = tmp_nwk('(((a:1,b:1):1,c:1):1,((d:1,e:1),f:1):1):0;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            left_leaf='a', right_leaf='c', leaves=None,
            orthogroup=False,
        )
        subtree_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.get_leaf_names()) == {'a', 'b', 'c'}
        # Verify branch lengths are preserved
        a_leaf = [l for l in tree.iter_leaves() if l.name == 'a'][0]
        assert abs(a_leaf.dist - 1.0) < 1e-6
        c_leaf = [l for l in tree.iter_leaves() if l.name == 'c'][0]
        assert abs(c_leaf.dist - 1.0) < 1e-6

    def test_subtree_all_branch_lengths_preserved(self, tmp_nwk, tmp_outfile):
        """Subtree extraction must preserve all internal branch lengths."""
        path = tmp_nwk('(((a:2,b:3):4,c:5):6,((d:1,e:1):1,f:1):1):0;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            left_leaf='a', right_leaf='c', leaves=None,
            orthogroup=False,
        )
        subtree_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.get_leaf_names()) == {'a', 'b', 'c'}
        leaves = {l.name: l.dist for l in tree.iter_leaves()}
        assert abs(leaves['a'] - 2.0) < 1e-6
        assert abs(leaves['b'] - 3.0) < 1e-6
        assert abs(leaves['c'] - 5.0) < 1e-6
        # Internal branch (a,b) parent should have dist 4
        ab_parent = tree.get_common_ancestor(['a', 'b'])
        assert abs(ab_parent.dist - 4.0) < 1e-6

    def test_leaves_mode_exact_distances(self, tmp_nwk, tmp_outfile):
        """Wiki --leaves mode: verify pairwise distances in extracted subtree."""
        path = tmp_nwk('(((a:2,b:3):4,c:5):6,(d:1,e:1):1):0;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            left_leaf=None, right_leaf=None, leaves='a,b',
            orthogroup=False,
        )
        subtree_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.get_leaf_names()) == {'a', 'b'}
        # a-b pairwise distance should be 2+3 = 5
        ab_dist = tree.get_distance('a', 'b')
        assert abs(ab_dist - 5.0) < 1e-6
