import os
import pytest
from ete3 import TreeNode

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
        assert set(tree.get_leaf_names()) == {'A', 'B', 'C', 'D'}

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
        leaves = {l.name: l.dist for l in tree.iter_leaves()}
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
        leaves = {l.name: l.dist for l in tree.iter_leaves()}
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
