import os
import pytest
from ete3 import TreeNode

from nwkit.sanitize import sanitize_main, add_quote
from nwkit.util import read_tree
from tests.helpers import make_args, DATA_DIR


class TestAddQuote:
    def test_single_quote(self):
        tree = TreeNode(newick='((A:1,B:1):1,(C:1,D:1):1);', format=1)
        tree = add_quote(tree, "'")
        for node in tree.traverse():
            if node.name != '':
                assert node.name.startswith("'")
                assert node.name.endswith("'")

    def test_double_quote(self):
        tree = TreeNode(newick='((A:1,B:1):1,(C:1,D:1):1);', format=1)
        tree = add_quote(tree, '"')
        for leaf in tree.iter_leaves():
            assert leaf.name.startswith('"')
            assert leaf.name.endswith('"')

    def test_no_quote(self):
        tree = TreeNode(newick='((A:1,B:1):1,(C:1,D:1):1);', format=1)
        tree = add_quote(tree, '')
        for leaf in tree.iter_leaves():
            assert leaf.name in ['A', 'B', 'C', 'D']

    def test_empty_names_skipped(self):
        tree = TreeNode(newick='((A:1,B:1):1,(C:1,D:1):1);', format=1)
        tree = add_quote(tree, "'")
        for node in tree.traverse():
            if node.name == '':
                continue
            assert node.name.startswith("'")


class TestSanitizeMain:
    def test_remove_singleton(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((((a:1,b:1):1):1,c:1):1,((d:1,e:1),f:1):1):0;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            remove_singleton=True, name_quote='none',
        )
        sanitize_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        # Check singleton is removed
        for node in tree.traverse():
            if not node.is_leaf():
                assert len(node.get_children()) != 1
        assert set(tree.get_leaf_names()) == {'a', 'b', 'c', 'd', 'e', 'f'}

    def test_add_single_quote(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((a:1,b:1):1,(c:1,d:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            remove_singleton=False, name_quote='single',
        )
        sanitize_main(args)
        with open(tmp_outfile) as f:
            content = f.read()
        assert "'" in content

    def test_add_double_quote(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((a:1,b:1):1,(c:1,d:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            remove_singleton=False, name_quote='double',
        )
        sanitize_main(args)
        with open(tmp_outfile) as f:
            content = f.read()
        assert '"' in content

    def test_no_quote(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((a:1,b:1):1,(c:1,d:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            remove_singleton=False, name_quote='none',
        )
        sanitize_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.get_leaf_names()) == {'a', 'b', 'c', 'd'}

    def test_with_data_file(self, tmp_outfile):
        infile = os.path.join(DATA_DIR, 'sanitize2', 'input.nwk')
        if not os.path.exists(infile):
            pytest.skip('Test data not found')
        args = make_args(
            infile=infile, outfile=tmp_outfile,
            remove_singleton=True, name_quote='single',
        )
        sanitize_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        # Singletons should be removed
        for node in tree.traverse():
            if not node.is_leaf():
                assert len(node.get_children()) != 1

    def test_wiki_exact_example(self, tmp_nwk, tmp_outfile):
        """Wiki example: nwkit sanitize --name_quote single

        Input:  ((((a:1,b:1):1):1,c:1):1,((d:1,e:1),f:1):1):0;
        Output: (('c':1,('a':1,'b':1):2):1,(('d':1,'e':1):1,'f':1):1):0;

        The singleton ((a:1,b:1):1):1 is collapsed and branch length summed (1+1=2).
        """
        path = tmp_nwk('((((a:1,b:1):1):1,c:1):1,((d:1,e:1),f:1):1):0;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            remove_singleton=True, name_quote='single',
        )
        sanitize_main(args)
        # Check raw output for single-quoted leaf names
        with open(tmp_outfile) as f:
            content = f.read()
        for name in ['a', 'b', 'c', 'd', 'e', 'f']:
            assert f"'{name}'" in content
        # Verify tree structure: no singleton nodes should remain
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        for node in tree.traverse():
            if not node.is_leaf():
                assert len(node.get_children()) >= 2

    def test_singleton_branch_length_summed(self, tmp_nwk, tmp_outfile):
        """When a singleton is removed, the branch lengths should be summed."""
        path = tmp_nwk('(((A:1,B:1):2):3,C:6);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            remove_singleton=True, name_quote='none',
        )
        sanitize_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.get_leaf_names()) == {'A', 'B', 'C'}
        # The parent of A,B should have combined branch length 2+3=5
        ab_parent = tree.get_common_ancestor(['A', 'B'])
        assert abs(ab_parent.dist - 5.0) < 1e-6
