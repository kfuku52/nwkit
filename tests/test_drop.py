import os
import pytest

from nwkit.drop import drop_main
from nwkit.util import read_tree
from tests.helpers import make_args


class TestDropMain:
    def test_drop_all_names(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:1)AB:1,(C:1,D:1)CD:1)root;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            target='all', name=True, support=False, length=False, fill=None,
        )
        drop_main(args)
        # When all names are dropped, the output won't contain original names
        with open(tmp_outfile) as f:
            content = f.read()
        assert 'AB' not in content
        assert 'CD' not in content
        assert 'root' not in content

    def test_drop_leaf_names(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:1)AB:1,(C:1,D:1)CD:1)root;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            target='leaf', name=True, support=False, length=False, fill=None,
        )
        drop_main(args)
        # Leaf names should be dropped, internal names preserved
        with open(tmp_outfile) as f:
            content = f.read()
        assert 'A' not in content.replace('AB', '')
        assert 'B' not in content.replace('AB', '')
        assert 'AB' in content

    def test_drop_intnode_names(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:1)AB:1,(C:1,D:1)CD:1)root;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            target='intnode', name=True, support=False, length=False, fill=None,
        )
        drop_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        # Leaf names should be preserved
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_drop_with_fill(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:1)AB:1,(C:1,D:1)CD:1)root;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            target='leaf', name=True, support=False, length=False, fill='UNKNOWN',
        )
        drop_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        for leaf in tree.leaves():
            assert leaf.name == 'UNKNOWN'

    def test_drop_support(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:1)90:1,(C:1,D:1)85:1);')
        args = make_args(
            infile=path, outfile=tmp_outfile, format='0',
            target='all', name=False, support=True, length=False, fill=None,
        )
        drop_main(args)
        # Should complete without error
        assert os.path.exists(tmp_outfile)

    def test_drop_length(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:2):3,(C:4,D:5):6);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            target='all', name=False, support=False, length=True, fill=None,
        )
        drop_main(args)
        assert os.path.exists(tmp_outfile)

    def test_issue10_drop_root_length_no_trailing_colon(self, tmp_nwk, tmp_outfile):
        """Regression test for GitHub issue #10.

        nwkit drop --target root --length yes used to produce output ending
        with ':;' instead of ';', which is malformed Newick.
        """
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1)root:0.5;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            target='root', name=False, support=False, length=True, fill=None,
        )
        drop_main(args)
        with open(tmp_outfile) as f:
            content = f.read().strip()
        # Output must not end with ':;'
        assert not content.endswith(':;'), f'Malformed Newick ending with ":;": {content}'
        assert content.endswith(';')

    def test_wiki_fill_unknown_exact(self, tmp_nwk, tmp_outfile):
        """Wiki example: drop leaf names with --fill unknown.

        Input:  ((A:1,B:1)AB:1,(C:1,D:1)CD:1)root;
        Output: ((unknown:1,unknown:1)AB:1,(unknown:1,unknown:1)CD:1)root:0;
        """
        path = tmp_nwk('((A:1,B:1)AB:1,(C:1,D:1)CD:1)root;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            target='leaf', name=True, support=False, length=False, fill='unknown',
        )
        drop_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        for leaf in tree.leaves():
            assert leaf.name == 'unknown'
        # Internal names preserved
        internal_names = [n.name for n in tree.traverse() if not n.is_leaf and n.name]
        assert 'AB' in internal_names
        assert 'CD' in internal_names
        assert 'root' in internal_names
        # Branch lengths preserved
        for leaf in tree.leaves():
            assert abs(leaf.dist - 1.0) < 1e-6

    def test_drop_length_with_fill_exact(self, tmp_nwk, tmp_outfile):
        """Drop leaf branch lengths with fill=0: internal lengths preserved."""
        path = tmp_nwk('((A:1,B:2):3,(C:4,D:5):6);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            target='leaf', name=False, support=False, length=True, fill='0',
        )
        drop_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D'}
        # Leaf branch lengths should be 0 (filled)
        for leaf in tree.leaves():
            assert abs(leaf.dist) < 1e-6
        # Internal branch lengths should be preserved
        for node in tree.traverse():
            if not node.is_leaf and not node.is_root:
                assert abs(node.dist - 3.0) < 1e-6 or abs(node.dist - 6.0) < 1e-6

    def test_drop_non_numeric_fill_for_support_raises(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            target='all', name=False, support=True, length=False, fill='NA',
        )
        with pytest.raises(ValueError, match='must be numeric'):
            drop_main(args)

    def test_drop_non_numeric_fill_for_length_raises(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            target='all', name=False, support=False, length=True, fill='NA',
        )
        with pytest.raises(ValueError, match='must be numeric'):
            drop_main(args)
