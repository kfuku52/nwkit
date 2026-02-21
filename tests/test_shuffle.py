import os
import random
import pytest
from ete4 import Tree

from nwkit.shuffle import get_shuffled_branch_lengths, print_rf_dist, shuffle_main
from nwkit.util import read_tree
from tests.helpers import make_args, DATA_DIR


class TestGetShuffledBranchLengths:
    def test_same_length(self):
        tree = Tree('((A:1,B:2):3,(C:4,D:5):6);', parser=1)
        nodes = list(tree.traverse())
        shuffled = get_shuffled_branch_lengths(nodes)
        assert len(shuffled) == len(nodes)

    def test_same_values(self):
        tree = Tree('((A:1,B:2):3,(C:4,D:5):6);', parser=1)
        nodes = list(tree.traverse())
        original = sorted([n.dist for n in nodes if n.dist is not None])
        shuffled = sorted([d for d in get_shuffled_branch_lengths(nodes) if d is not None])
        assert original == shuffled


class TestPrintRfDist:
    def test_same_tree(self, capsys):
        t1 = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        t2 = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        print_rf_dist(t1, t2)
        captured = capsys.readouterr()
        assert 'Robinson-Foulds distance' in captured.err
        assert '0' in captured.err

    def test_different_trees(self, capsys):
        t1 = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        t2 = Tree('((A:1,C:1):1,(B:1,D:1):1);', parser=1)
        print_rf_dist(t1, t2)
        captured = capsys.readouterr()
        assert 'Robinson-Foulds distance' in captured.err


class TestShuffleMain:
    def test_shuffle_branch_length(self, tmp_nwk, tmp_outfile):
        random.seed(42)
        path = tmp_nwk('((A:1,B:2):3,(C:4,D:5):6);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            topology=False, branch_length=True, label=False,
        )
        shuffle_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D'}
        # Branch lengths should be shuffled (same values, different assignment)
        branch_lengths = sorted([n.dist for n in tree.traverse() if not n.is_root])
        assert sorted(branch_lengths) == sorted([1, 2, 3, 4, 5, 6])

    def test_shuffle_labels(self, tmp_nwk, tmp_outfile):
        random.seed(42)
        path = tmp_nwk('((A:1,B:2):3,(C:4,D:5):6);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            topology=False, branch_length=False, label=True,
        )
        shuffle_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        # Same leaf names, but potentially different positions
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_shuffle_labels_reports_nonzero_rf_when_topology_changes(self, tmp_nwk, tmp_outfile, capsys):
        random.seed(1)
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            topology=False, branch_length=False, label=True,
        )
        shuffle_main(args)
        captured = capsys.readouterr()
        rf_line = [line for line in captured.err.splitlines() if 'Robinson-Foulds distance' in line][0]
        rf_value = int(rf_line.split('=')[1].split('(')[0].strip())
        assert rf_value > 0

    def test_shuffle_topology(self, tmp_nwk, tmp_outfile):
        random.seed(42)
        path = tmp_nwk('((A:1,B:2):3,(C:4,D:5):6);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            topology=True, branch_length=False, label=False,
        )
        shuffle_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_shuffle_preserves_leaf_count(self, tmp_nwk, tmp_outfile):
        random.seed(42)
        nwk = '(((A:1,B:1):1,C:1):1,((D:1,E:1):1,F:1):1);'
        path = tmp_nwk(nwk)
        args = make_args(
            infile=path, outfile=tmp_outfile,
            topology=True, branch_length=True, label=True,
        )
        shuffle_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert len(list(tree.leaf_names())) == 6

    def test_shuffle_unrooted_tree_does_not_crash(self, tmp_nwk, tmp_outfile):
        random.seed(7)
        path = tmp_nwk('(A:1,B:1,C:1,D:1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            topology=False, branch_length=False, label=True,
        )
        shuffle_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_shuffle_with_unnamed_leaf_skips_rf_but_writes_output(self, tmp_nwk, tmp_outfile, capsys):
        random.seed(9)
        path = tmp_nwk('((:1,B:1):1,(C:1,D:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            topology=False, branch_length=False, label=True,
        )
        shuffle_main(args)
        captured = capsys.readouterr()
        assert 'Skipping RF distance' in captured.err
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert len(list(tree.leaf_names())) == 4
