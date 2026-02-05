import os
import pytest
from ete3 import TreeNode

from nwkit.prune import prune_main
from nwkit.util import read_tree
from tests.helpers import make_args, DATA_DIR


class TestPruneMain:
    def test_prune_by_exact_name(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('(((A1:2,(B1:1,B2:1):1):1,(A2:1,C1:1):2):1,C2:4):0.25;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            pattern='B1', invert_match=False,
        )
        prune_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        leaf_names = set(tree.get_leaf_names())
        assert 'B1' not in leaf_names
        assert 'A1' in leaf_names

    def test_prune_by_regex(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('(((A1:2,(B1:1,B2:1):1):1,(A2:1,C1:1):2):1,C2:4):0.25;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            pattern='B.*', invert_match=False,
        )
        prune_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        leaf_names = set(tree.get_leaf_names())
        assert 'B1' not in leaf_names
        assert 'B2' not in leaf_names
        assert 'A1' in leaf_names

    def test_prune_invert_match(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('(((A1:2,(B1:1,B2:1):1):1,(A2:1,C1:1):2):1,C2:4):0.25;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            pattern='A.*', invert_match=True,
        )
        prune_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        leaf_names = set(tree.get_leaf_names())
        # Invert match: keep only A-matching leaves, prune everything else
        assert all(name.startswith('A') for name in leaf_names)

    def test_prune_preserves_remaining(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            pattern='A', invert_match=False,
        )
        prune_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        leaf_names = set(tree.get_leaf_names())
        assert leaf_names == {'B', 'C', 'D'}

    def test_with_data_file(self, tmp_outfile):
        infile = os.path.join(DATA_DIR, 'prune2', 'tree.nwk')
        if not os.path.exists(infile):
            pytest.skip('Test data not found')
        args = make_args(
            infile=infile, outfile=tmp_outfile,
            pattern='B1', invert_match=False,
        )
        prune_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert 'B1' not in set(tree.get_leaf_names())

    def test_wiki_pipe_separated_pattern(self, tmp_nwk, tmp_outfile):
        """Wiki example: nwkit prune --pattern "B1|C2"

        Input:  (((A1:2.0,(B1:1.0,B2:1.0):1.0):1.0,(A2:1.0,C1:1.0):2.0):1.0,C2:4.0):0.25;
        Output: (((A1:2,B2:2)1:1,(A2:1,C1:1)1:2)1:1)1:0.25;
        """
        path = tmp_nwk('(((A1:2.0,(B1:1.0,B2:1.0):1.0):1.0,(A2:1.0,C1:1.0):2.0):1.0,C2:4.0):0.25;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            pattern='B1|C2', invert_match=False,
        )
        prune_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        leaf_names = set(tree.get_leaf_names())
        assert leaf_names == {'A1', 'B2', 'A2', 'C1'}
        assert 'B1' not in leaf_names
        assert 'C2' not in leaf_names

    def test_wiki_prune_branch_length_after_removal(self, tmp_nwk, tmp_outfile):
        """Wiki: when B1 is pruned, B2 absorbs the dissolved parent's branch length.

        Input:  (((A1:2.0,(B1:1.0,B2:1.0):1.0):1.0,(A2:1.0,C1:1.0):2.0):1.0,C2:4.0):0.25;
        After pruning B1|C2:
        Output: (((A1:2,B2:2):1,(A2:1,C1:1):2):1):0.25;

        B2 dist = 1.0 (original) + 1.0 (dissolved parent) = 2.0
        A1 dist = 2.0 (unchanged)
        """
        path = tmp_nwk('(((A1:2.0,(B1:1.0,B2:1.0):1.0):1.0,(A2:1.0,C1:1.0):2.0):1.0,C2:4.0):0.25;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            pattern='B1|C2', invert_match=False,
        )
        prune_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        leaves = {l.name: l.dist for l in tree.iter_leaves()}
        assert abs(leaves['B2'] - 2.0) < 1e-6
        assert abs(leaves['A1'] - 2.0) < 1e-6
        assert abs(leaves['A2'] - 1.0) < 1e-6
        assert abs(leaves['C1'] - 1.0) < 1e-6

    def test_wiki_invert_match_exact(self, tmp_nwk, tmp_outfile):
        """Wiki invert_match example: keep only A.* and B.* tips.

        Input:  (((A1:2.0,(B1:1.0,B2:1.0):1.0):1.0,(A2:1.0,C1:1.0):2.0):1.0,C2:4.0):0.25;
        Output: (((A1:2,(B1:1,B2:1):1):1,A2:3):1):0.25;

        A2 absorbs C1's dissolved parent: dist = 1.0 + 2.0 = 3.0
        """
        path = tmp_nwk('(((A1:2.0,(B1:1.0,B2:1.0):1.0):1.0,(A2:1.0,C1:1.0):2.0):1.0,C2:4.0):0.25;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            pattern='A.*|B.*', invert_match=True,
        )
        prune_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        leaf_names = set(tree.get_leaf_names())
        assert leaf_names == {'A1', 'A2', 'B1', 'B2'}
        leaves = {l.name: l.dist for l in tree.iter_leaves()}
        assert abs(leaves['A1'] - 2.0) < 1e-6
        assert abs(leaves['B1'] - 1.0) < 1e-6
        assert abs(leaves['B2'] - 1.0) < 1e-6
        assert abs(leaves['A2'] - 3.0) < 1e-6
