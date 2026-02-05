import os
import pytest
from ete3 import TreeNode

from nwkit.root import midpoint_rooting, outgroup_rooting, transfer_root, root_main
from nwkit.util import read_tree, is_rooted
from tests.helpers import make_args, DATA_DIR


class TestMidpointRooting:
    def test_basic(self):
        tree = TreeNode(newick='(A:1,B:3,(C:2,D:4):2);', format=1)
        rooted = midpoint_rooting(tree)
        assert is_rooted(rooted)

    def test_already_rooted(self):
        tree = TreeNode(newick='((A:1,B:1):1,(C:1,D:1):1);', format=1)
        rooted = midpoint_rooting(tree)
        assert is_rooted(rooted)
        assert set(rooted.get_leaf_names()) == {'A', 'B', 'C', 'D'}


class TestOutgroupRooting:
    def test_single_outgroup(self):
        tree = TreeNode(newick='(A:1,(B:1,(C:1,D:1):1):1);', format=1)
        rooted = outgroup_rooting(tree, 'A')
        assert is_rooted(rooted)
        # A should be one of the subroot children
        subroot_children = rooted.get_children()
        outgroup_leaves = set()
        for child in subroot_children:
            outgroup_leaves.update(child.get_leaf_names())
        assert 'A' in outgroup_leaves

    def test_multiple_outgroups(self):
        tree = TreeNode(newick='((A:1,B:1):1,(C:1,(D:1,E:1):1):1);', format=1)
        rooted = outgroup_rooting(tree, 'D,E')
        assert is_rooted(rooted)
        leaf_names = set(rooted.get_leaf_names())
        assert leaf_names == {'A', 'B', 'C', 'D', 'E'}

    def test_outgroup_not_found(self):
        tree = TreeNode(newick='((A:1,B:1):1,(C:1,D:1):1);', format=1)
        with pytest.raises(SystemExit):
            outgroup_rooting(tree, 'Z')


class TestTransferRoot:
    def test_transfer(self):
        # tree_from is rooted with (A,B) | (C,D)
        tree_from = TreeNode(newick='((A:1,B:1):1,(C:1,D:1):1);', format=1)
        # tree_to is unrooted
        tree_to = TreeNode(newick='(A:1,B:1,(C:1,D:1):1);', format=1)
        result = transfer_root(tree_to, tree_from)
        assert is_rooted(result)
        assert set(result.get_leaf_names()) == {'A', 'B', 'C', 'D'}


class TestRootMain:
    def test_midpoint(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('(A:1,B:5,(C:2,D:4):2);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            method='midpoint',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert is_rooted(tree)

    def test_outgroup(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('(A:1,(B:1,(C:1,D:1):1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            method='outgroup', outgroup='A',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert is_rooted(tree)
        assert set(tree.get_leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_transfer(self, tmp_nwk, tmp_outfile):
        path1 = tmp_nwk('(A:1,B:1,(C:1,D:1):1);', 'tree1.nwk')
        path2 = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);', 'tree2.nwk')
        args = make_args(
            infile=path1, outfile=tmp_outfile,
            method='transfer', infile2=path2, format2='auto',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert is_rooted(tree)

    def test_transfer_mismatched_raises(self, tmp_nwk):
        path1 = tmp_nwk('(A:1,B:1,(C:1,D:1):1);', 'tree1.nwk')
        path2 = tmp_nwk('((A:1,B:1):1,(C:1,E:1):1);', 'tree2.nwk')
        args = make_args(
            infile=path1, outfile='-',
            method='transfer', infile2=path2, format2='auto',
        )
        with pytest.raises(Exception, match='Leaf labels'):
            root_main(args)

    def test_wiki_outgroup_single(self, tmp_nwk, tmp_outfile):
        """Wiki example: nwkit root --method outgroup --outgroup a

        Input:  (((a:1,b:1):1,c:1):1,((d:1,e:1),f:1):1):0;
        Output: (a:0.5,(b:1,(c:1,((d:1,e:1):1,f:1):2):1):0.5):0;
        """
        path = tmp_nwk('(((a:1,b:1):1,c:1):1,((d:1,e:1),f:1):1):0;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            method='outgroup', outgroup='a',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert is_rooted(tree)
        assert set(tree.get_leaf_names()) == {'a', 'b', 'c', 'd', 'e', 'f'}
        # 'a' should be in one of the two subroot clades, by itself
        children = tree.get_children()
        child_leaf_sets = [set(c.get_leaf_names()) for c in children]
        assert {'a'} in child_leaf_sets or any('a' in s and len(s) == 1 for s in child_leaf_sets)

    def test_wiki_outgroup_multiple(self, tmp_nwk, tmp_outfile):
        """Wiki example: nwkit root --method outgroup --outgroup a,b

        Input:  (((a:1,b:1):1,c:1):1,((d:1,e:1),f:1):1):0;
        Output: ((a:1,b:1):0.5,(c:1,((d:1,e:1):1,f:1):2):0.5):0;
        """
        path = tmp_nwk('(((a:1,b:1):1,c:1):1,((d:1,e:1),f:1):1):0;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            method='outgroup', outgroup='a,b',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert is_rooted(tree)
        assert set(tree.get_leaf_names()) == {'a', 'b', 'c', 'd', 'e', 'f'}
        # {a,b} should be one of the two subroot clades
        children = tree.get_children()
        child_leaf_sets = [set(c.get_leaf_names()) for c in children]
        assert {'a', 'b'} in child_leaf_sets

    def test_wiki_root_transfer(self, tmp_nwk, tmp_outfile):
        """Wiki example: nwkit root --method transfer

        input1.nwk: (((a:1,b:1):1,c:1):1,((d:1,e:1),f:1):1):0;
        input2.nwk: ((((a:1,b:1):1,c:1),f:1):1,(d:3,e:3):1):0;
        Output: ((d:1,e:1):0.5,(f:1,(c:1,(b:1,a:1):1):2):0.5)Root:0;
        """
        path1 = tmp_nwk('(((a:1,b:1):1,c:1):1,((d:1,e:1),f:1):1):0;', 'tree1.nwk')
        path2 = tmp_nwk('((((a:1,b:1):1,c:1),f:1):1,(d:3,e:3):1):0;', 'tree2.nwk')
        args = make_args(
            infile=path1, outfile=tmp_outfile,
            method='transfer', infile2=path2, format2='auto',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert is_rooted(tree)
        assert set(tree.get_leaf_names()) == {'a', 'b', 'c', 'd', 'e', 'f'}
        # {d,e} should be in one subroot clade (like the source tree)
        children = tree.get_children()
        child_leaf_sets = [set(c.get_leaf_names()) for c in children]
        assert {'d', 'e'} in child_leaf_sets

    def test_wiki_midpoint_rooting(self, tmp_nwk, tmp_outfile):
        """Wiki example: nwkit root --method midpoint

        Input:  ((((a:5,b:1):1,c:3),f:1):1,(d:1,e:1):1):0;
        """
        path = tmp_nwk('((((a:5,b:1):1,c:3),f:1):1,(d:1,e:1):1):0;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            method='midpoint',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert is_rooted(tree)
        assert set(tree.get_leaf_names()) == {'a', 'b', 'c', 'd', 'e', 'f'}
        # The longest path goes through 'a', so 'a' should be near one subroot
        children = tree.get_children()
        # One side should contain 'a'
        sides_with_a = [c for c in children if 'a' in c.get_leaf_names()]
        assert len(sides_with_a) == 1
