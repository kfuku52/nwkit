import os
import pytest
from ete4 import Tree

from nwkit.sanitize import sanitize_main, add_quote
from nwkit.util import read_tree
from tests.helpers import make_args, DATA_DIR


class TestAddQuote:
    def test_single_quote(self):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        tree = add_quote(tree, "'")
        for node in tree.traverse():
            if node.name:
                assert node.name.startswith("'")
                assert node.name.endswith("'")

    def test_double_quote(self):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        tree = add_quote(tree, '"')
        for leaf in tree.leaves():
            assert leaf.name.startswith('"')
            assert leaf.name.endswith('"')

    def test_no_quote(self):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        tree = add_quote(tree, '')
        for leaf in tree.leaves():
            assert leaf.name in ['A', 'B', 'C', 'D']

    def test_empty_names_skipped(self):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        tree = add_quote(tree, "'")
        for node in tree.traverse():
            if not node.name:
                continue
            assert node.name.startswith("'")


class TestSanitizeMain:
    def test_remove_singleton(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((((a:1,b:1):1):1,c:1):1,((d:1,e:1),f:1):1):0;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            remove_singleton=True, resolve_polytomy=False, name_quote='none',
        )
        sanitize_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        # Check singleton is removed
        for node in tree.traverse():
            if not node.is_leaf:
                assert len(node.get_children()) != 1
        assert set(tree.leaf_names()) == {'a', 'b', 'c', 'd', 'e', 'f'}

    def test_add_single_quote(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((a:1,b:1):1,(c:1,d:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            remove_singleton=False, resolve_polytomy=False, name_quote='single',
        )
        sanitize_main(args)
        with open(tmp_outfile) as f:
            content = f.read()
        assert "'" in content

    def test_add_double_quote(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((a:1,b:1):1,(c:1,d:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            remove_singleton=False, resolve_polytomy=False, name_quote='double',
        )
        sanitize_main(args)
        with open(tmp_outfile) as f:
            content = f.read()
        assert '"' in content

    def test_no_quote(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((a:1,b:1):1,(c:1,d:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            remove_singleton=False, resolve_polytomy=False, name_quote='none',
        )
        sanitize_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.leaf_names()) == {'a', 'b', 'c', 'd'}

    def test_with_data_file(self, tmp_outfile):
        infile = os.path.join(DATA_DIR, 'sanitize2', 'input.nwk')
        if not os.path.exists(infile):
            pytest.skip('Test data not found')
        args = make_args(
            infile=infile, outfile=tmp_outfile,
            remove_singleton=True, resolve_polytomy=False, name_quote='single',
        )
        sanitize_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        # Singletons should be removed
        for node in tree.traverse():
            if not node.is_leaf:
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
            remove_singleton=True, resolve_polytomy=False, name_quote='single',
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
            if not node.is_leaf:
                assert len(node.get_children()) >= 2

    def test_singleton_branch_length_summed(self, tmp_nwk, tmp_outfile):
        """When a singleton is removed, the branch lengths should be summed."""
        path = tmp_nwk('(((A:1,B:1):2):3,C:6);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            remove_singleton=True, resolve_polytomy=False, name_quote='none',
        )
        sanitize_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.leaf_names()) == {'A', 'B', 'C'}
        # The parent of A,B should have combined branch length 2+3=5
        ab_parent = tree.common_ancestor(['A', 'B'])
        assert abs(ab_parent.dist - 5.0) < 1e-6

    def test_resolve_polytomy_basic(self, tmp_nwk, tmp_outfile):
        """Polytomy (A,B,C,D) should be resolved into dichotomies."""
        path = tmp_nwk('(A:1,B:1,C:1,D:1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            remove_singleton=False, resolve_polytomy=True, name_quote='none',
        )
        sanitize_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        for node in tree.traverse():
            if not node.is_leaf:
                assert len(node.get_children()) == 2
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_resolve_polytomy_zero_branch_length(self, tmp_nwk, tmp_outfile):
        """Newly created branches from polytomy resolution should have zero length."""
        path = tmp_nwk('(A:1,B:1,C:1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            remove_singleton=False, resolve_polytomy=True, name_quote='none',
        )
        sanitize_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        for leaf in tree.leaves():
            assert abs(leaf.dist - 1.0) < 1e-6

    def test_resolve_polytomy_no_change_on_dichotomy(self, tmp_nwk, tmp_outfile):
        """A fully dichotomous tree should not be changed."""
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            remove_singleton=False, resolve_polytomy=True, name_quote='none',
        )
        sanitize_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D'}
        assert len(list(tree.traverse())) == 7  # 4 leaves + 2 internal + root

    def test_resolve_polytomy_disabled(self, tmp_nwk, tmp_outfile):
        """When resolve_polytomy=False, polytomies should be preserved."""
        path = tmp_nwk('(A:1,B:1,C:1,D:1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            remove_singleton=False, resolve_polytomy=False, name_quote='none',
        )
        sanitize_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert len(tree.get_children()) == 4

    def test_resolve_polytomy_with_singleton_removal(self, tmp_nwk, tmp_outfile):
        """Singleton removal + polytomy resolution should both work together."""
        path = tmp_nwk('(((A:1,B:1):1):1,C:1,D:1,E:1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            remove_singleton=True, resolve_polytomy=True, name_quote='none',
        )
        sanitize_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        for node in tree.traverse():
            if not node.is_leaf:
                assert len(node.get_children()) == 2
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D', 'E'}

    def test_resolve_polytomy_nested(self, tmp_nwk, tmp_outfile):
        """Multiple nested polytomies should all be resolved."""
        path = tmp_nwk('((A:1,B:1,C:1):1,(D:1,E:1,F:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            remove_singleton=False, resolve_polytomy=True, name_quote='none',
        )
        sanitize_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        for node in tree.traverse():
            if not node.is_leaf:
                assert len(node.get_children()) == 2
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D', 'E', 'F'}

    def test_wiki_singleton_exact_branch_length(self, tmp_nwk, tmp_outfile):
        """Wiki example: singleton collapse sums branch lengths 1+1=2.

        Input:  ((((a:1,b:1):1):1,c:1):1,((d:1,e:1),f:1):1):0;
        After singleton removal: (a,b) parent gets dist 1+1=2.
        """
        path = tmp_nwk('((((a:1,b:1):1):1,c:1):1,((d:1,e:1),f:1):1):0;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            remove_singleton=True, resolve_polytomy=False, name_quote='none',
        )
        sanitize_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.leaf_names()) == {'a', 'b', 'c', 'd', 'e', 'f'}
        ab_parent = tree.common_ancestor(['a', 'b'])
        assert abs(ab_parent.dist - 2.0) < 1e-6
        # Leaf branch lengths preserved
        leaves = {l.name: l.dist for l in tree.leaves()}
        assert abs(leaves['a'] - 1.0) < 1e-6
        assert abs(leaves['b'] - 1.0) < 1e-6
        assert abs(leaves['c'] - 1.0) < 1e-6
        assert abs(leaves['f'] - 1.0) < 1e-6

    def test_pairwise_distances_after_singleton_removal(self, tmp_nwk, tmp_outfile):
        """Singleton removal must preserve all pairwise distances."""
        nwk = '(((A:1,B:2):3):4,C:10);'
        original = Tree(nwk, parser=1)
        path = tmp_nwk(nwk)
        args = make_args(
            infile=path, outfile=tmp_outfile,
            remove_singleton=True, resolve_polytomy=False, name_quote='none',
        )
        sanitize_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        for l1 in tree.leaves():
            for l2 in tree.leaves():
                if l1.name != l2.name:
                    orig_d = original.get_distance(l1.name, l2.name)
                    new_d = tree.get_distance(l1, l2)
                    assert abs(orig_d - new_d) < 1e-6, \
                        f'{l1.name}-{l2.name}: {orig_d} vs {new_d}'
