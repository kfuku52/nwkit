import pytest
from ete4 import Tree

from nwkit.label import label_main
from nwkit.util import read_tree
from tests.helpers import make_args


class TestLabelMain:
    def test_label_unnamed_nodes(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            target='all', prefix='n', force=False,
        )
        label_main(args)
        tree = read_tree(tmp_outfile, format='1', quoted_node_names=True, quiet=True)
        for node in tree.traverse():
            assert node.name != ''

    def test_label_with_custom_prefix(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            target='all', prefix='node_', force=False,
        )
        label_main(args)
        tree = read_tree(tmp_outfile, format='1', quoted_node_names=True, quiet=True)
        internal_names = [n.name for n in tree.traverse() if not n.is_leaf]
        assert all(name.startswith('node_') for name in internal_names)

    def test_label_leaves_only(self, tmp_nwk, tmp_outfile):
        # Use a tree where leaves have names but we force-overwrite them
        path = tmp_nwk('((X:1,Y:1):1,(Z:1,W:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            target='leaf', prefix='leaf', force=True,
        )
        label_main(args)
        tree = read_tree(tmp_outfile, format='1', quoted_node_names=True, quiet=True)
        for leaf in tree.leaves():
            assert leaf.name.startswith('leaf')

    def test_label_force_overwrite(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:1)AB:1,(C:1,D:1)CD:1)root;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            target='all', prefix='n', force=True,
        )
        label_main(args)
        tree = read_tree(tmp_outfile, format='1', quoted_node_names=True, quiet=True)
        names = [n.name for n in tree.traverse()]
        # With force, all names should be overwritten
        assert all(name.startswith('n') for name in names)

    def test_label_no_force_preserves_existing(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:1)AB:1,(C:1,D:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            target='all', prefix='n', force=False,
        )
        label_main(args)
        tree = read_tree(tmp_outfile, format='1', quoted_node_names=True, quiet=True)
        names = [n.name for n in tree.traverse()]
        # AB should be preserved, unnamed nodes should get new names
        assert 'AB' in names

    def test_label_intnode_only(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            target='intnode', prefix='int', force=False,
        )
        label_main(args)
        tree = read_tree(tmp_outfile, format='1', quoted_node_names=True, quiet=True)
        # Leaf names should be preserved
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D'}
        # Internal nodes should have new labels
        for node in tree.traverse():
            if not node.is_leaf:
                assert node.name.startswith('int')

    def test_wiki_intnode_exact_names(self, tmp_nwk, tmp_outfile):
        """Wiki example: label internal nodes with default prefix 'n'.

        Input:  ((A:1,B:1):1,(C:1,D:1):1);
        Output: ((A:1,B:1)n1:1,(C:1,D:1)n2:1)n0:0;

        All 3 internal nodes get sequential labels n0, n1, n2.
        """
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            target='intnode', prefix='n', force=False,
        )
        label_main(args)
        tree = read_tree(tmp_outfile, format='1', quoted_node_names=True, quiet=True)
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D'}
        internal_names = sorted([n.name for n in tree.traverse() if not n.is_leaf])
        assert set(internal_names) == {'n0', 'n1', 'n2'}
        # Branch lengths preserved
        for leaf in tree.leaves():
            assert abs(leaf.dist - 1.0) < 1e-6

    def test_wiki_preserve_existing_names(self, tmp_nwk, tmp_outfile):
        """Wiki example: preserve existing name AB, label unnamed nodes.

        Input:  ((A:1,B:1)AB:1,(C:1,D:1):1);
        Output: ((A:1,B:1)AB:1,(C:1,D:1)n1:1)n0:0;
        """
        path = tmp_nwk('((A:1,B:1)AB:1,(C:1,D:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            target='all', prefix='n', force=False,
        )
        label_main(args)
        tree = read_tree(tmp_outfile, format='1', quoted_node_names=True, quiet=True)
        all_names = [n.name for n in tree.traverse()]
        assert 'AB' in all_names
        new_names = [n for n in all_names if n.startswith('n')]
        assert len(new_names) == 2
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_wiki_force_overwrite_all(self, tmp_nwk, tmp_outfile):
        """Wiki example: force overwrite all names with custom prefix.

        Input:  ((A:1,B:1):1,(C:1,D:1):1);
        Output: ((node3:1,node4:1)node1:1,(node5:1,node6:1)node2:1)node0:0;

        All 7 nodes renamed with prefix 'node' followed by sequential number.
        """
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            target='all', prefix='node', force=True,
        )
        label_main(args)
        tree = read_tree(tmp_outfile, format='1', quoted_node_names=True, quiet=True)
        all_names = [n.name for n in tree.traverse()]
        assert len(all_names) == 7
        assert all(name.startswith('node') for name in all_names)
        assert len(set(all_names)) == 7
        nums = sorted([int(name.replace('node', '')) for name in all_names])
        assert nums == list(range(7))
