import os
import pytest
from ete4 import Tree

from nwkit.root import midpoint_rooting, outgroup_rooting, transfer_root, mad_rooting, mv_rooting, root_main
from nwkit.util import read_tree, is_rooted
from tests.helpers import make_args, DATA_DIR, safe_get_distance


class TestMidpointRooting:
    def test_basic(self):
        tree = Tree('(A:1,B:3,(C:2,D:4):2);', parser=1)
        rooted = midpoint_rooting(tree)
        assert is_rooted(rooted)

    def test_already_rooted(self):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        rooted = midpoint_rooting(tree)
        assert is_rooted(rooted)
        assert set(rooted.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_symmetric_clock_exact_distances(self):
        """On a symmetric clock tree, all root-to-tip distances are equal."""
        tree = Tree('((A:2,B:2):1,(C:2,D:2):1);', parser=1)
        rooted = midpoint_rooting(tree)
        dists = [safe_get_distance(rooted, rooted, l) for l in rooted.leaves()]
        assert all(abs(d - 3.0) < 1e-6 for d in dists)

    def test_pairwise_distance_preservation(self):
        """Midpoint rooting must preserve all pairwise distances."""
        nwk = '((A:3,B:1):2,(C:4,D:2):1);'
        original = Tree(nwk, parser=1)
        tree = Tree(nwk, parser=1)
        rooted = midpoint_rooting(tree)
        for l1 in rooted.leaves():
            for l2 in rooted.leaves():
                if l1.name != l2.name:
                    assert abs(original.get_distance(l1.name, l2.name) -
                               rooted.get_distance(l1, l2)) < 1e-6


class TestOutgroupRooting:
    def test_single_outgroup(self):
        tree = Tree('(A:1,(B:1,(C:1,D:1):1):1);', parser=1)
        rooted = outgroup_rooting(tree, 'A')
        assert is_rooted(rooted)
        # A should be one of the subroot children
        subroot_children = rooted.get_children()
        outgroup_leaves = set()
        for child in subroot_children:
            outgroup_leaves.update(child.leaf_names())
        assert 'A' in outgroup_leaves

    def test_multiple_outgroups(self):
        tree = Tree('((A:1,B:1):1,(C:1,(D:1,E:1):1):1);', parser=1)
        rooted = outgroup_rooting(tree, 'D,E')
        assert is_rooted(rooted)
        leaf_names = set(rooted.leaf_names())
        assert leaf_names == {'A', 'B', 'C', 'D', 'E'}

    def test_outgroup_not_found(self):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        with pytest.raises(SystemExit):
            outgroup_rooting(tree, 'Z')

    def test_pairwise_distance_preservation(self):
        """Outgroup rooting must preserve all pairwise distances."""
        nwk = '((A:2,B:3):1,(C:4,D:5):1);'
        original = Tree(nwk, parser=1)
        tree = Tree(nwk, parser=1)
        rooted = outgroup_rooting(tree, 'A')
        for l1 in rooted.leaves():
            for l2 in rooted.leaves():
                if l1.name != l2.name:
                    assert abs(original.get_distance(l1.name, l2.name) -
                               rooted.get_distance(l1, l2)) < 1e-6

    def test_multiple_outgroup_exact_bipartition(self):
        """Outgroup rooting with multiple tips creates correct bipartition."""
        tree = Tree('((A:1,B:1):1,(C:1,(D:1,E:1):1):1);', parser=1)
        rooted = outgroup_rooting(tree, 'D,E')
        children = rooted.get_children()
        child_leaf_sets = [set(c.leaf_names()) for c in children]
        assert {'D', 'E'} in child_leaf_sets
        assert {'A', 'B', 'C'} in child_leaf_sets


class TestMadRooting:
    def test_basic(self):
        tree = Tree('((A:1,B:2):1,(C:3,(D:1,E:2):1):1);', parser=1)
        rooted = mad_rooting(tree)
        assert is_rooted(rooted)

    def test_preserves_leaves(self):
        tree = Tree('((A:1,B:2):1,(C:3,(D:1,E:2):1):1);', parser=1)
        rooted = mad_rooting(tree)
        assert set(rooted.leaf_names()) == {'A', 'B', 'C', 'D', 'E'}

    def test_asymmetric_branch_lengths(self):
        tree = Tree('((A:10,B:1):1,(C:1,D:1):1);', parser=1)
        rooted = mad_rooting(tree)
        assert is_rooted(rooted)
        assert set(rooted.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_wiki_exact_root_split(self):
        """MAD rooting on wiki tree: verify root split and pairwise distances.

        Input: ((A:1,B:2):1,(C:3,(D:1,E:2):1):1);
        Wiki output: ((A:1,B:2):1.58461,(C:3,(D:1,E:2):1):0.415389);
        Total root branch ≈ 2.0
        """
        tree = Tree('((A:1,B:2):1,(C:3,(D:1,E:2):1):1);', parser=1)
        original = Tree('((A:1,B:2):1,(C:3,(D:1,E:2):1):1);', parser=1)
        rooted = mad_rooting(tree)
        # Total root branch should sum to ≈ 2.0
        children = rooted.get_children()
        total = sum(c.dist for c in children)
        assert abs(total - 2.0) < 0.01
        # Pairwise distances must be preserved
        for l1 in rooted.leaves():
            for l2 in rooted.leaves():
                if l1.name != l2.name:
                    orig_d = original.get_distance(l1.name, l2.name)
                    new_d = rooted.get_distance(l1, l2)
                    assert abs(orig_d - new_d) < 0.01


class TestMvRooting:
    def test_basic(self):
        tree = Tree('((A:1,B:2):1,(C:3,(D:1,E:2):1):1);', parser=1)
        rooted = mv_rooting(tree)
        assert is_rooted(rooted)

    def test_preserves_leaves(self):
        tree = Tree('((A:1,B:2):1,(C:3,(D:1,E:2):1):1);', parser=1)
        rooted = mv_rooting(tree)
        assert set(rooted.leaf_names()) == {'A', 'B', 'C', 'D', 'E'}

    def test_clock_like_tree(self):
        """On a clock-like tree, MV root should achieve near-zero variance."""
        import numpy as np
        tree = Tree('((A:2,B:2):1,(C:2,D:2):1);', parser=1)
        rooted = mv_rooting(tree)
        dists = [safe_get_distance(rooted, rooted, leaf) for leaf in rooted.leaves()]
        assert np.var(dists) < 1e-10

    def test_asymmetric_branch_lengths(self):
        tree = Tree('((A:10,B:1):1,(C:1,D:1):1);', parser=1)
        rooted = mv_rooting(tree)
        assert is_rooted(rooted)
        assert set(rooted.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_wiki_exact_root_position(self):
        """MV rooting on wiki tree: root at x=5/12 from (C,(D,E)) side.

        Input: ((A:1,B:2):1,(C:3,(D:1,E:2):1):1);
        Expected root-to-tip distances:
          C: 3 + 5/12 = 41/12, D: 1+1+5/12 = 29/12, E: 2+1+5/12 = 41/12,
          A: 1 + 19/12 = 31/12, B: 2 + 19/12 = 43/12
        """
        tree = Tree('((A:1,B:2):1,(C:3,(D:1,E:2):1):1);', parser=1)
        rooted = mv_rooting(tree)
        dists = {l.name: safe_get_distance(rooted, rooted, l) for l in rooted.leaves()}
        assert abs(dists['C'] - 41/12) < 1e-4
        assert abs(dists['D'] - 29/12) < 1e-4
        assert abs(dists['E'] - 41/12) < 1e-4
        assert abs(dists['A'] - 31/12) < 1e-4
        assert abs(dists['B'] - 43/12) < 1e-4
        # Root children: (C,(D,E)) side=5/12, (A,B) side=19/12
        children = rooted.get_children()
        for c in children:
            if 'C' in c.leaf_names():
                assert abs(c.dist - 5/12) < 1e-4
            else:
                assert abs(c.dist - 19/12) < 1e-4

    def test_pairwise_distance_preservation(self):
        """MV rooting must preserve all pairwise distances."""
        nwk = '((A:1,B:2):1,(C:3,(D:1,E:2):1):1);'
        original = Tree(nwk, parser=1)
        tree = Tree(nwk, parser=1)
        rooted = mv_rooting(tree)
        for l1 in rooted.leaves():
            for l2 in rooted.leaves():
                if l1.name != l2.name:
                    assert abs(original.get_distance(l1.name, l2.name) -
                               rooted.get_distance(l1, l2)) < 1e-4

    def test_three_taxa_exact(self):
        """MV rooting on 3-taxon tree: verify pairwise distances preserved."""
        tree = Tree('(A:1,B:3,C:5);', parser=1)
        rooted = mv_rooting(tree)
        for l1 in rooted.leaves():
            for l2 in rooted.leaves():
                if l1.name != l2.name:
                    d = rooted.get_distance(l1, l2)
                    if {l1.name, l2.name} == {'A', 'B'}:
                        assert abs(d - 4) < 1e-6
                    elif {l1.name, l2.name} == {'A', 'C'}:
                        assert abs(d - 6) < 1e-6
                    elif {l1.name, l2.name} == {'B', 'C'}:
                        assert abs(d - 8) < 1e-6


class TestTransferRoot:
    def test_transfer(self):
        # tree_from is rooted with (A,B) | (C,D)
        tree_from = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        # tree_to is unrooted
        tree_to = Tree('(A:1,B:1,(C:1,D:1):1);', parser=1)
        result = transfer_root(tree_to, tree_from)
        assert is_rooted(result)
        assert set(result.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_pairwise_distance_preservation(self):
        """Transfer root must preserve pairwise distances from tree_to."""
        tree_from = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        tree_to = Tree('(A:2,B:3,(C:4,D:5):1);', parser=1)
        original = Tree('(A:2,B:3,(C:4,D:5):1);', parser=1)
        result = transfer_root(tree_to, tree_from)
        for l1 in result.leaves():
            for l2 in result.leaves():
                if l1.name != l2.name:
                    orig_d = original.get_distance(l1.name, l2.name)
                    new_d = result.get_distance(l1, l2)
                    assert abs(orig_d - new_d) < 1e-6, \
                        f'{l1.name}-{l2.name}: {orig_d} vs {new_d}'


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
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D'}

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

    def test_mad(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:2):1,(C:3,(D:1,E:2):1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            method='mad',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert is_rooted(tree)
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D', 'E'}

    def test_mv(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:2):1,(C:3,(D:1,E:2):1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            method='mv',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert is_rooted(tree)
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D', 'E'}

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
        assert set(tree.leaf_names()) == {'a', 'b', 'c', 'd', 'e', 'f'}
        # 'a' should be in one of the two subroot clades, by itself
        children = tree.get_children()
        child_leaf_sets = [set(c.leaf_names()) for c in children]
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
        assert set(tree.leaf_names()) == {'a', 'b', 'c', 'd', 'e', 'f'}
        # {a,b} should be one of the two subroot clades
        children = tree.get_children()
        child_leaf_sets = [set(c.leaf_names()) for c in children]
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
        assert set(tree.leaf_names()) == {'a', 'b', 'c', 'd', 'e', 'f'}
        # {d,e} should be in one subroot clade (like the source tree)
        children = tree.get_children()
        child_leaf_sets = [set(c.leaf_names()) for c in children]
        assert {'d', 'e'} in child_leaf_sets

    def test_wiki_midpoint_rooting(self, tmp_nwk, tmp_outfile):
        """Wiki example: nwkit root --method midpoint

        Input:  ((((a:5,b:1):1,c:3):1,f:1):1,(d:1,e:1):1):0;
        Output: ((a:5,b:1):0.5,(c:3,(f:1,(d:1,e:1):2):1):0.5):0;
        """
        # Note: Added explicit :1 to internal node that was missing dist in original wiki example
        path = tmp_nwk('((((a:5,b:1):1,c:3):1,f:1):1,(d:1,e:1):1):0;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            method='midpoint',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert is_rooted(tree)
        assert set(tree.leaf_names()) == {'a', 'b', 'c', 'd', 'e', 'f'}
        # Verify exact root-to-tip distances from wiki output
        dists = {l.name: safe_get_distance(tree, tree, l) for l in tree.leaves()}
        assert abs(dists['a'] - 5.5) < 1e-6
        assert abs(dists['b'] - 1.5) < 1e-6
        assert abs(dists['c'] - 3.5) < 1e-6
        assert abs(dists['f'] - 2.5) < 1e-6
        assert abs(dists['d'] - 4.5) < 1e-6
        assert abs(dists['e'] - 4.5) < 1e-6
        # Root children should both have dist 0.5
        children = tree.get_children()
        for c in children:
            assert abs(c.dist - 0.5) < 1e-6

    def test_wiki_outgroup_single_exact_distances(self, tmp_nwk, tmp_outfile):
        """Wiki outgroup single: verify exact root-to-tip distances.

        Input:  (((a:1,b:1):1,c:1):1,((d:1,e:1):1,f:1):1):0;
        Output: (a:0.5,(b:1,(c:1,((d:1,e:1):1,f:1):2):1):0.5):0;
        """
        # Note: Added explicit :1 to internal node that was missing dist in original wiki example
        path = tmp_nwk('(((a:1,b:1):1,c:1):1,((d:1,e:1):1,f:1):1):0;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            method='outgroup', outgroup='a',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        dists = {l.name: safe_get_distance(tree, tree, l) for l in tree.leaves()}
        assert abs(dists['a'] - 0.5) < 1e-6
        assert abs(dists['b'] - 1.5) < 1e-6
        assert abs(dists['c'] - 2.5) < 1e-6
        assert abs(dists['f'] - 4.5) < 1e-6
        assert abs(dists['d'] - 5.5) < 1e-6
        assert abs(dists['e'] - 5.5) < 1e-6

    def test_wiki_outgroup_multiple_exact_distances(self, tmp_nwk, tmp_outfile):
        """Wiki outgroup multiple: verify exact root-to-tip distances.

        Input:  (((a:1,b:1):1,c:1):1,((d:1,e:1):1,f:1):1):0;
        Output: ((a:1,b:1):0.5,(c:1,((d:1,e:1):1,f:1):2):0.5):0;
        """
        # Note: Added explicit :1 to internal node that was missing dist in original wiki example
        path = tmp_nwk('(((a:1,b:1):1,c:1):1,((d:1,e:1):1,f:1):1):0;')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            method='outgroup', outgroup='a,b',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        dists = {l.name: safe_get_distance(tree, tree, l) for l in tree.leaves()}
        assert abs(dists['a'] - 1.5) < 1e-6
        assert abs(dists['b'] - 1.5) < 1e-6
        assert abs(dists['c'] - 1.5) < 1e-6
        assert abs(dists['f'] - 3.5) < 1e-6
        assert abs(dists['d'] - 4.5) < 1e-6
        assert abs(dists['e'] - 4.5) < 1e-6

    def test_wiki_transfer_exact_distances(self, tmp_nwk, tmp_outfile):
        """Wiki transfer: verify exact root-to-tip distances.

        input1: (((a:1,b:1):1,c:1):1,((d:1,e:1):1,f:1):1):0;
        input2: ((((a:1,b:1):1,c:1):1,f:1):1,(d:3,e:3):1):0;
        Output: ((d:1,e:1):0.5,(f:1,(c:1,(b:1,a:1):1):2):0.5)Root:0;
        """
        # Note: Added explicit :1 to internal nodes that were missing dist in original wiki example
        path1 = tmp_nwk('(((a:1,b:1):1,c:1):1,((d:1,e:1):1,f:1):1):0;', 'tree1.nwk')
        path2 = tmp_nwk('((((a:1,b:1):1,c:1):1,f:1):1,(d:3,e:3):1):0;', 'tree2.nwk')
        args = make_args(
            infile=path1, outfile=tmp_outfile,
            method='transfer', infile2=path2, format2='auto',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        dists = {l.name: safe_get_distance(tree, tree, l) for l in tree.leaves()}
        assert abs(dists['d'] - 1.5) < 1e-6
        assert abs(dists['e'] - 1.5) < 1e-6
        assert abs(dists['f'] - 1.5) < 1e-6
        assert abs(dists['c'] - 3.5) < 1e-6
        assert abs(dists['a'] - 4.5) < 1e-6
        assert abs(dists['b'] - 4.5) < 1e-6

    def test_wiki_mv_exact_distances(self, tmp_nwk, tmp_outfile):
        """Wiki MV rooting: verify exact root-to-tip distances.

        Input: ((A:1,B:2):1,(C:3,(D:1,E:2):1):1);
        Output: ((C:3,(D:1,E:2):1):0.416667,(A:1,B:2):1.58333);
        Root at x=5/12 from (C,(D,E)) side.
        """
        path = tmp_nwk('((A:1,B:2):1,(C:3,(D:1,E:2):1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            method='mv',
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        dists = {l.name: safe_get_distance(tree, tree, l) for l in tree.leaves()}
        assert abs(dists['C'] - 41/12) < 1e-4
        assert abs(dists['D'] - 29/12) < 1e-4
        assert abs(dists['E'] - 41/12) < 1e-4
        assert abs(dists['A'] - 31/12) < 1e-4
        assert abs(dists['B'] - 43/12) < 1e-4
