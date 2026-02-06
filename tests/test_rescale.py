import pytest
from ete4 import Tree

from nwkit.rescale import rescale_main
from nwkit.util import read_tree
from tests.helpers import make_args


class TestRescaleMain:
    def test_rescale_all(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:2):3,(C:4,D:5):6);')
        args = make_args(infile=path, outfile=tmp_outfile, target='all', factor=2.0)
        rescale_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        leaves = {l.name: l.dist for l in tree.leaves()}
        assert abs(leaves['A'] - 2.0) < 1e-6
        assert abs(leaves['B'] - 4.0) < 1e-6
        assert abs(leaves['C'] - 8.0) < 1e-6
        assert abs(leaves['D'] - 10.0) < 1e-6

    def test_rescale_leaf_only(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:2):3,(C:4,D:5):6);')
        args = make_args(infile=path, outfile=tmp_outfile, target='leaf', factor=10.0)
        rescale_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        leaves = {l.name: l.dist for l in tree.leaves()}
        assert abs(leaves['A'] - 10.0) < 1e-6
        # Internal nodes should not be scaled
        for node in tree.traverse():
            if not node.is_leaf and not node.is_root:
                assert abs(node.dist - 3.0) < 1e-6 or abs(node.dist - 6.0) < 1e-6

    def test_rescale_intnode_only(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:2):3,(C:4,D:5):6);')
        args = make_args(infile=path, outfile=tmp_outfile, target='intnode', factor=0.5)
        rescale_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        leaves = {l.name: l.dist for l in tree.leaves()}
        # Leaves should remain unchanged
        assert abs(leaves['A'] - 1.0) < 1e-6
        assert abs(leaves['B'] - 2.0) < 1e-6

    def test_rescale_factor_zero(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:2):3,(C:4,D:5):6);')
        args = make_args(infile=path, outfile=tmp_outfile, target='all', factor=0.0)
        rescale_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        for node in tree.traverse():
            if node.dist is not None:
                assert abs(node.dist) < 1e-6

    def test_rescale_preserves_topology(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:2):3,(C:4,D:5):6);')
        args = make_args(infile=path, outfile=tmp_outfile, target='all', factor=3.0)
        rescale_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_wiki_all_factor_half(self, tmp_nwk, tmp_outfile):
        """Wiki example: rescale all branches by factor 0.5.

        Input:  ((A:2,B:4):6,(C:8,D:10):12);
        Output: ((A:1,B:2):3,(C:4,D:5):6);
        """
        path = tmp_nwk('((A:2,B:4):6,(C:8,D:10):12);')
        args = make_args(infile=path, outfile=tmp_outfile, target='all', factor=0.5)
        rescale_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        leaves = {l.name: l.dist for l in tree.leaves()}
        assert abs(leaves['A'] - 1.0) < 1e-6
        assert abs(leaves['B'] - 2.0) < 1e-6
        assert abs(leaves['C'] - 4.0) < 1e-6
        assert abs(leaves['D'] - 5.0) < 1e-6
        for node in tree.traverse():
            if not node.is_leaf and not node.is_root:
                assert abs(node.dist - 3.0) < 1e-6 or abs(node.dist - 6.0) < 1e-6

    def test_wiki_leaf_only_factor_tenth(self, tmp_nwk, tmp_outfile):
        """Wiki example: rescale only leaf branches by factor 0.1.

        Input:  ((A:2,B:4):6,(C:8,D:10):12);
        Output: ((A:0.2,B:0.4):6,(C:0.8,D:1):12);
        """
        path = tmp_nwk('((A:2,B:4):6,(C:8,D:10):12);')
        args = make_args(infile=path, outfile=tmp_outfile, target='leaf', factor=0.1)
        rescale_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        leaves = {l.name: l.dist for l in tree.leaves()}
        assert abs(leaves['A'] - 0.2) < 1e-6
        assert abs(leaves['B'] - 0.4) < 1e-6
        assert abs(leaves['C'] - 0.8) < 1e-6
        assert abs(leaves['D'] - 1.0) < 1e-6
        for node in tree.traverse():
            if not node.is_leaf and not node.is_root:
                assert abs(node.dist - 6.0) < 1e-6 or abs(node.dist - 12.0) < 1e-6

    def test_intnode_exact_values(self, tmp_nwk, tmp_outfile):
        """Rescale only internal branches: leaves unchanged, internals scaled."""
        path = tmp_nwk('((A:1,B:2):3,(C:4,D:5):6);')
        args = make_args(infile=path, outfile=tmp_outfile, target='intnode', factor=0.5)
        rescale_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        leaves = {l.name: l.dist for l in tree.leaves()}
        assert abs(leaves['A'] - 1.0) < 1e-6
        assert abs(leaves['B'] - 2.0) < 1e-6
        assert abs(leaves['C'] - 4.0) < 1e-6
        assert abs(leaves['D'] - 5.0) < 1e-6
        for node in tree.traverse():
            if not node.is_leaf and not node.is_root:
                assert abs(node.dist - 1.5) < 1e-6 or abs(node.dist - 3.0) < 1e-6
