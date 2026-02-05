import pytest
from ete3 import TreeNode

from nwkit.rescale import rescale_main
from nwkit.util import read_tree
from tests.helpers import make_args


class TestRescaleMain:
    def test_rescale_all(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:2):3,(C:4,D:5):6);')
        args = make_args(infile=path, outfile=tmp_outfile, target='all', factor=2.0)
        rescale_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        leaves = {l.name: l.dist for l in tree.iter_leaves()}
        assert abs(leaves['A'] - 2.0) < 1e-6
        assert abs(leaves['B'] - 4.0) < 1e-6
        assert abs(leaves['C'] - 8.0) < 1e-6
        assert abs(leaves['D'] - 10.0) < 1e-6

    def test_rescale_leaf_only(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:2):3,(C:4,D:5):6);')
        args = make_args(infile=path, outfile=tmp_outfile, target='leaf', factor=10.0)
        rescale_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        leaves = {l.name: l.dist for l in tree.iter_leaves()}
        assert abs(leaves['A'] - 10.0) < 1e-6
        # Internal nodes should not be scaled
        for node in tree.traverse():
            if not node.is_leaf() and not node.is_root():
                assert abs(node.dist - 3.0) < 1e-6 or abs(node.dist - 6.0) < 1e-6

    def test_rescale_intnode_only(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:2):3,(C:4,D:5):6);')
        args = make_args(infile=path, outfile=tmp_outfile, target='intnode', factor=0.5)
        rescale_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        leaves = {l.name: l.dist for l in tree.iter_leaves()}
        # Leaves should remain unchanged
        assert abs(leaves['A'] - 1.0) < 1e-6
        assert abs(leaves['B'] - 2.0) < 1e-6

    def test_rescale_factor_zero(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:2):3,(C:4,D:5):6);')
        args = make_args(infile=path, outfile=tmp_outfile, target='all', factor=0.0)
        rescale_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        for node in tree.traverse():
            assert abs(node.dist) < 1e-6

    def test_rescale_preserves_topology(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:2):3,(C:4,D:5):6);')
        args = make_args(infile=path, outfile=tmp_outfile, target='all', factor=3.0)
        rescale_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.get_leaf_names()) == {'A', 'B', 'C', 'D'}
