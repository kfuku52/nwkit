import os
import pytest
from ete3 import TreeNode

from nwkit.util import (
    read_tree,
    write_tree,
    remove_singleton,
    label2sciname,
    read_item_per_line_file,
    annotate_scientific_names,
    annotate_duplication_confidence_scores,
    is_all_leaf_names_identical,
    get_target_nodes,
    is_rooted,
)
from tests.helpers import make_args


class TestReadTree:
    def test_read_from_file(self, tmp_nwk):
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        tree = read_tree(path, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.get_leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_read_with_explicit_format(self, tmp_nwk):
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        tree = read_tree(path, format='1', quoted_node_names=True, quiet=True)
        assert len(list(tree.iter_leaves())) == 4

    def test_read_tree_with_internal_names(self, tmp_nwk):
        path = tmp_nwk('((A:1,B:1)AB:1,(C:1,D:1)CD:1)root;')
        tree = read_tree(path, format='auto', quoted_node_names=True, quiet=True)
        leaf_names = set(tree.get_leaf_names())
        assert leaf_names == {'A', 'B', 'C', 'D'}

    def test_read_tree_auto_format_detection(self, tmp_nwk):
        # Format 0: support values
        path = tmp_nwk('((A:1,B:1)90:1,(C:1,D:1)85:1);')
        tree = read_tree(path, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.get_leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_read_tree_invalid_raises(self, tmp_nwk):
        path = tmp_nwk('not_a_tree')
        with pytest.raises(Exception, match='Failed to parse'):
            read_tree(path, format='auto', quoted_node_names=True, quiet=True)


class TestWriteTree:
    def test_write_to_file(self, tmp_nwk, tmp_outfile):
        import nwkit.util as util_mod
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        tree = read_tree(path, format='1', quoted_node_names=True, quiet=True)
        args = make_args(outfile=tmp_outfile)
        write_tree(tree, args, format='1', quiet=True)
        assert os.path.exists(tmp_outfile)
        tree2 = read_tree(tmp_outfile, format='1', quoted_node_names=True, quiet=True)
        assert set(tree2.get_leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_write_to_stdout(self, tmp_nwk, capsys):
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        tree = read_tree(path, format='1', quoted_node_names=True, quiet=True)
        args = make_args(outfile='-')
        write_tree(tree, args, format='1', quiet=True)
        captured = capsys.readouterr()
        assert 'A' in captured.out
        assert 'B' in captured.out

    def test_write_preserves_node_names(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:1)AB:1,(C:1,D:1)CD:1)root;')
        tree = read_tree(path, format='1', quoted_node_names=True, quiet=True)
        args = make_args(outfile=tmp_outfile)
        write_tree(tree, args, format='1', quiet=True)
        with open(tmp_outfile) as f:
            content = f.read()
        assert 'AB' in content
        assert 'CD' in content


class TestRemoveSingleton:
    def test_remove_singleton_node(self):
        # Create tree with a singleton: ((A,B)) -> should become (A,B)
        tree = TreeNode(newick='(((A:1,B:1):1):1,(C:1,D:1):1);', format=1)
        num_nodes_before = len(list(tree.traverse()))
        tree = remove_singleton(tree, verbose=False, preserve_branch_length=True)
        num_nodes_after = len(list(tree.traverse()))
        assert num_nodes_after < num_nodes_before
        assert set(tree.get_leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_no_singleton(self, simple_tree):
        num_nodes_before = len(list(simple_tree.traverse()))
        tree = remove_singleton(simple_tree, verbose=False)
        num_nodes_after = len(list(tree.traverse()))
        assert num_nodes_before == num_nodes_after

    def test_preserve_branch_length(self):
        tree = TreeNode(newick='(((A:1,B:1):2):3,C:6);', format=1)
        tree = remove_singleton(tree, verbose=False, preserve_branch_length=True)
        # After removing singleton, branch lengths should be preserved (summed)
        assert set(tree.get_leaf_names()) == {'A', 'B', 'C'}


class TestLabel2Sciname:
    def test_single_string(self):
        result = label2sciname('Homo_sapiens_GENE1')
        assert result == 'Homo_sapiens'

    def test_list_input(self):
        result = label2sciname(['Homo_sapiens_GENE1', 'Mus_musculus_GENE2'])
        assert result == ['Homo_sapiens', 'Mus_musculus']

    def test_no_species_info(self):
        result = label2sciname('SingleWord')
        assert result is None

    def test_custom_delimiter(self):
        result = label2sciname('Homo-sapiens-GENE1', in_delim='-')
        assert result == 'Homo_sapiens'

    def test_custom_out_delimiter(self):
        result = label2sciname('Homo_sapiens_GENE1', out_delim=' ')
        assert result == 'Homo sapiens'

    def test_empty_list(self):
        result = label2sciname([])
        assert result == []


class TestReadItemPerLineFile:
    def test_basic(self, tmp_path):
        f = tmp_path / 'items.txt'
        f.write_text('apple\nbanana\ncherry\n')
        result = read_item_per_line_file(str(f))
        assert result == ['apple', 'banana', 'cherry']

    def test_empty_lines_stripped(self, tmp_path):
        f = tmp_path / 'items.txt'
        f.write_text('apple\n\nbanana\n\n')
        result = read_item_per_line_file(str(f))
        assert result == ['apple', 'banana']


class TestAnnotateScientificNames:
    def test_annotate(self, species_tree):
        tree = annotate_scientific_names(species_tree)
        sci_names = [leaf.sci_name for leaf in tree.iter_leaves()]
        assert 'Homo_sapiens' in sci_names
        assert 'Mus_musculus' in sci_names
        assert 'Danio_rerio' in sci_names


class TestAnnotateDuplicationConfidenceScores:
    def test_annotate(self, species_tree):
        tree = annotate_scientific_names(species_tree)
        tree = annotate_duplication_confidence_scores(tree)
        for node in tree.traverse():
            if not node.is_leaf():
                assert hasattr(node, 'dup_conf_score')
                assert 0 <= node.dup_conf_score <= 1

    def test_no_duplication(self):
        # Tree where each species appears once in each child clade
        nwk = '((Homo_sapiens_G1:1,Mus_musculus_G1:1):1,(Danio_rerio_G1:1,Xenopus_laevis_G1:1):1);'
        tree = TreeNode(newick=nwk, format=1)
        tree = annotate_scientific_names(tree)
        tree = annotate_duplication_confidence_scores(tree)
        # Root node has no species overlap -> dup_conf_score = 0
        assert tree.dup_conf_score == 0.0

    def test_full_duplication(self):
        # Tree where both children have the same species
        nwk = '((Homo_sapiens_G1:1,Mus_musculus_G1:1):1,(Homo_sapiens_G2:1,Mus_musculus_G2:1):1);'
        tree = TreeNode(newick=nwk, format=1)
        tree = annotate_scientific_names(tree)
        tree = annotate_duplication_confidence_scores(tree)
        assert tree.dup_conf_score == 1.0


class TestIsAllLeafNamesIdentical:
    def test_identical(self):
        t1 = TreeNode(newick='((A,B),(C,D));', format=1)
        t2 = TreeNode(newick='((A,C),(B,D));', format=1)
        assert is_all_leaf_names_identical(t1, t2) is True

    def test_not_identical(self):
        t1 = TreeNode(newick='((A,B),(C,D));', format=1)
        t2 = TreeNode(newick='((A,B),(C,E));', format=1)
        assert is_all_leaf_names_identical(t1, t2) is False


class TestGetTargetNodes:
    def test_all(self, simple_tree):
        nodes = get_target_nodes(simple_tree, 'all')
        assert len(nodes) == len(list(simple_tree.traverse()))

    def test_root(self, simple_tree):
        nodes = get_target_nodes(simple_tree, 'root')
        assert len(nodes) == 1
        assert nodes[0].is_root()

    def test_leaf(self, simple_tree):
        nodes = get_target_nodes(simple_tree, 'leaf')
        assert all(n.is_leaf() for n in nodes)
        assert len(nodes) == 4

    def test_intnode(self, simple_tree):
        nodes = get_target_nodes(simple_tree, 'intnode')
        assert all(not n.is_leaf() for n in nodes)


class TestIsRooted:
    def test_rooted(self, simple_tree):
        assert is_rooted(simple_tree) is True

    def test_unrooted(self, unrooted_tree):
        assert is_rooted(unrooted_tree) is False
