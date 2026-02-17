import os
import pytest
from ete4 import Tree

from nwkit.transfer import transfer_main
from nwkit.util import read_tree
from tests.helpers import make_args, DATA_DIR


class TestTransferMain:
    def test_transfer_names(self, tmp_nwk, tmp_outfile):
        # tree1 has leaf names but no internal names
        path1 = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);', 'tree1.nwk')
        # tree2 has internal names
        path2 = tmp_nwk('((A:1,B:1)AB:1,(C:1,D:1)CD:1)root;', 'tree2.nwk')
        args = make_args(
            infile=path1, infile2=path2, outfile=tmp_outfile,
            format2='auto', target='all',
            name=True, support=False, length=False, fill=None,
        )
        transfer_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_transfer_branch_lengths(self, tmp_nwk, tmp_outfile):
        path1 = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);', 'tree1.nwk')
        path2 = tmp_nwk('((A:10,B:20):30,(C:40,D:50):60);', 'tree2.nwk')
        args = make_args(
            infile=path1, infile2=path2, outfile=tmp_outfile,
            format2='auto', target='all',
            name=False, support=False, length=True, fill=None,
        )
        transfer_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        leaves = {l.name: l.dist for l in tree.leaves()}
        assert abs(leaves['A'] - 10) < 1e-6
        assert abs(leaves['B'] - 20) < 1e-6

    def test_transfer_leaf_only(self, tmp_nwk, tmp_outfile):
        path1 = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);', 'tree1.nwk')
        path2 = tmp_nwk('((A:10,B:20):30,(C:40,D:50):60);', 'tree2.nwk')
        args = make_args(
            infile=path1, infile2=path2, outfile=tmp_outfile,
            format2='auto', target='leaf',
            name=False, support=False, length=True, fill=None,
        )
        transfer_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        leaves = {l.name: l.dist for l in tree.leaves()}
        assert abs(leaves['A'] - 10) < 1e-6

    def test_mismatched_leaves_raises(self, tmp_nwk):
        path1 = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);', 'tree1.nwk')
        path2 = tmp_nwk('((A:1,B:1):1,(C:1,E:1):1);', 'tree2.nwk')
        args = make_args(
            infile=path1, infile2=path2, outfile='-',
            format2='auto', target='all',
            name=True, support=False, length=False, fill=None,
        )
        with pytest.raises(Exception, match='Leaf labels'):
            transfer_main(args)

    def test_transfer_with_fill(self, tmp_nwk, tmp_outfile):
        # Same topology but different internal structure (5 leaves, both rooted)
        # tree1 has ((A,B),(C,(D,E))), tree2 has ((A,B),(C,(D,E))) with extra node
        path1 = tmp_nwk('(((A:1,B:1):1,C:1):1,(D:1,E:1):1);', 'tree1.nwk')
        path2 = tmp_nwk('(((A:1,C:1):1,B:1):1,(D:1,E:1):1);', 'tree2.nwk')
        args = make_args(
            infile=path1, infile2=path2, outfile=tmp_outfile,
            format2='auto', target='intnode',
            name=True, support=False, length=False, fill='NA',
        )
        transfer_main(args)
        assert os.path.exists(tmp_outfile)

    def test_wiki_name_transfer_exact(self, tmp_nwk, tmp_outfile):
        """Wiki example: transfer internal node names.

        Input1: ((A:1,B:2):3,(C:4,D:5):6);
        Input2: ((A:0.1,B:0.2)clade_AB:0.3,(C:0.4,D:0.5)clade_CD:0.6)root;
        Output: ((A:1,B:2)clade_AB:3,(C:4,D:5)clade_CD:6)root:0;
        """
        path1 = tmp_nwk('((A:1,B:2):3,(C:4,D:5):6);', 'tree1.nwk')
        path2 = tmp_nwk('((A:0.1,B:0.2)clade_AB:0.3,(C:0.4,D:0.5)clade_CD:0.6)root;', 'tree2.nwk')
        args = make_args(
            infile=path1, infile2=path2, outfile=tmp_outfile,
            format2='auto', target='intnode',
            name=True, support=False, length=False, fill=None,
        )
        transfer_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        internal_names = [n.name for n in tree.traverse() if not n.is_leaf]
        assert 'clade_AB' in internal_names
        assert 'clade_CD' in internal_names
        assert 'root' in internal_names
        leaves = {l.name: l.dist for l in tree.leaves()}
        assert abs(leaves['A'] - 1.0) < 1e-6
        assert abs(leaves['B'] - 2.0) < 1e-6
        assert abs(leaves['C'] - 4.0) < 1e-6
        assert abs(leaves['D'] - 5.0) < 1e-6

    def test_wiki_length_transfer_exact(self, tmp_nwk, tmp_outfile):
        """Transfer branch lengths from tree2 to tree1: all 4 leaves exact."""
        path1 = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);', 'tree1.nwk')
        path2 = tmp_nwk('((A:10,B:20):30,(C:40,D:50):60);', 'tree2.nwk')
        args = make_args(
            infile=path1, infile2=path2, outfile=tmp_outfile,
            format2='auto', target='all',
            name=False, support=False, length=True, fill=None,
        )
        transfer_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        leaves = {l.name: l.dist for l in tree.leaves()}
        assert abs(leaves['A'] - 10.0) < 1e-6
        assert abs(leaves['B'] - 20.0) < 1e-6
        assert abs(leaves['C'] - 40.0) < 1e-6
        assert abs(leaves['D'] - 50.0) < 1e-6

    def test_support_transfer_root_target_handles_missing_root_support(self, tmp_nwk, tmp_outfile):
        path1 = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);', 'tree1.nwk')
        path2 = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1)root;', 'tree2.nwk')
        args = make_args(
            infile=path1, infile2=path2, outfile=tmp_outfile,
            format2='auto', target='root',
            name=False, support=True, length=False, fill=None,
        )
        transfer_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        root = [n for n in tree.traverse() if n.is_root][0]
        assert root.support is None

    def test_transfer_handles_numeric_root_support_during_reroot(self, tmp_nwk, tmp_outfile):
        path1 = tmp_nwk('((A:1,B:1)10:1,(C:1,D:1)20:1)30;', 'tree1.nwk')
        path2 = tmp_nwk('((A:10,B:20)40:1,(C:30,D:40)50:1)60;', 'tree2.nwk')
        args = make_args(
            infile=path1, infile2=path2, outfile=tmp_outfile,
            format2='auto', target='leaf',
            name=False, support=False, length=True, fill=None,
        )
        transfer_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        leaves = {l.name: l.dist for l in tree.leaves()}
        assert abs(leaves['A'] - 10.0) < 1e-6
        assert abs(leaves['B'] - 20.0) < 1e-6
        assert abs(leaves['C'] - 30.0) < 1e-6
        assert abs(leaves['D'] - 40.0) < 1e-6

    def test_support_transfer_root_target_with_numeric_root_support(self, tmp_nwk, tmp_outfile):
        path1 = tmp_nwk('((A:1,B:1)10:1,(C:1,D:1)20:1)30;', 'tree1.nwk')
        path2 = tmp_nwk('((A:1,B:1)40:1,(C:1,D:1)50:1)60;', 'tree2.nwk')
        args = make_args(
            infile=path1, infile2=path2, outfile=tmp_outfile,
            format2='auto', target='root',
            name=False, support=True, length=False, fill=None,
        )
        transfer_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        root = [n for n in tree.traverse() if n.is_root][0]
        assert abs(root.support - 60.0) < 1e-6
