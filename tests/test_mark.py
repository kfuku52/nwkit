import os
import sys
import pytest
from argparse import Namespace
from ete3 import TreeNode

from nwkit.mark import annotate_tree_attr, get_insert_nodes, label_insert_nodes, mark_main
from nwkit.util import read_tree
from tests.helpers import make_args, DATA_DIR


def make_mark_args(**kwargs):
    defaults = {
        'infile': '-',
        'outfile': '-',
        'format': 'auto',
        'outformat': 'auto',
        'quoted_node_names': True,
        'pattern': '.*',
        'target': 'clade',
        'target_only_clade': True,
        'insert_txt': 'MARKED',
        'insert_sep': '',
        'insert_pos': 'suffix',
    }
    defaults.update(kwargs)
    return Namespace(**defaults)


class TestAnnotateTreeAttr:
    def test_annotate_all_leaves(self):
        tree = TreeNode(newick='(((A1:2,(B1:1,B2:1):1):1,(A2:1,C1:1):2):1,C2:4):0.25;', format=1)
        args = make_mark_args(pattern='.*')
        tree = annotate_tree_attr(tree, args)
        for leaf in tree.iter_leaves():
            assert leaf.is_target_leaf is True

    def test_annotate_pattern_match(self):
        tree = TreeNode(newick='(((A1:2,(B1:1,B2:1):1):1,(A2:1,C1:1):2):1,C2:4):0.25;', format=1)
        args = make_mark_args(pattern='A.*')
        tree = annotate_tree_attr(tree, args)
        for leaf in tree.iter_leaves():
            if leaf.name.startswith('A'):
                assert leaf.is_target_leaf is True
            else:
                assert leaf.is_target_leaf is False

    def test_mrca_annotation(self):
        tree = TreeNode(newick='(((A1:2,(B1:1,B2:1):1):1,(A2:1,C1:1):2):1,C2:4):0.25;', format=1)
        args = make_mark_args(pattern='B.*')
        tree = annotate_tree_attr(tree, args)
        # The MRCA of B1 and B2 should be marked
        mrca_nodes = [n for n in tree.traverse() if n.is_target_only_mrca]
        assert len(mrca_nodes) >= 1

    def test_no_match(self):
        tree = TreeNode(newick='((A:1,B:1):1,(C:1,D:1):1);', format=1)
        args = make_mark_args(pattern='Z.*')
        tree = annotate_tree_attr(tree, args)
        for leaf in tree.iter_leaves():
            assert leaf.is_target_leaf is False


class TestGetInsertNodes:
    def test_mrca_target(self):
        tree = TreeNode(newick='(((A1:2,(B1:1,B2:1):1):1,(A2:1,C1:1):2):1,C2:4):0.25;', format=1)
        args = make_mark_args(pattern='B.*', target='mrca', target_only_clade=True)
        tree = annotate_tree_attr(tree, args)
        nodes = get_insert_nodes(tree, args)
        assert len(nodes) >= 1

    def test_leaf_target(self):
        tree = TreeNode(newick='(((A1:2,(B1:1,B2:1):1):1,(A2:1,C1:1):2):1,C2:4):0.25;', format=1)
        args = make_mark_args(pattern='A.*', target='leaf')
        tree = annotate_tree_attr(tree, args)
        nodes = get_insert_nodes(tree, args)
        assert len(nodes) == 2  # A1, A2
        assert all(n.is_leaf() for n in nodes)

    def test_clade_target(self):
        tree = TreeNode(newick='(((A1:2,(B1:1,B2:1):1):1,(A2:1,C1:1):2):1,C2:4):0.25;', format=1)
        args = make_mark_args(pattern='B.*', target='clade', target_only_clade=True)
        tree = annotate_tree_attr(tree, args)
        nodes = get_insert_nodes(tree, args)
        assert len(nodes) >= 2  # B1, B2, and their MRCA


class TestLabelInsertNodes:
    def test_suffix(self):
        tree = TreeNode(newick='((A:1,B:1):1,(C:1,D:1):1);', format=1)
        args = make_mark_args(pattern='A', target='leaf', insert_txt='Fg', insert_pos='suffix')
        tree = annotate_tree_attr(tree, args)
        tree = label_insert_nodes(tree, args)
        a_leaf = [l for l in tree.iter_leaves() if 'A' in l.name][0]
        assert a_leaf.name.endswith('Fg')

    def test_prefix(self):
        tree = TreeNode(newick='((A:1,B:1):1,(C:1,D:1):1);', format=1)
        args = make_mark_args(pattern='A', target='leaf', insert_txt='Fg', insert_pos='prefix')
        tree = annotate_tree_attr(tree, args)
        tree = label_insert_nodes(tree, args)
        a_leaf = [l for l in tree.iter_leaves() if 'A' in l.name][0]
        assert a_leaf.name.startswith('Fg')

    def test_separator(self):
        tree = TreeNode(newick='((A:1,B:1):1,(C:1,D:1):1);', format=1)
        args = make_mark_args(pattern='A', target='leaf', insert_txt='Fg', insert_sep='_', insert_pos='suffix')
        tree = annotate_tree_attr(tree, args)
        tree = label_insert_nodes(tree, args)
        a_leaf = [l for l in tree.iter_leaves() if 'A' in l.name][0]
        assert a_leaf.name == 'A_Fg'


class TestMarkMain:
    def test_mark_main_suffix(self, tmp_nwk, tmp_outfile, monkeypatch):
        monkeypatch.setattr(sys, 'argv', ['nwkit', 'mark'])
        path = tmp_nwk('(((A1:2,(B1:1,B2:1):1):1,(A2:1,C1:1):2):1,C2:4):0.25;')
        args = make_mark_args(infile=path, outfile=tmp_outfile, pattern='A.*',
                              target='leaf', insert_txt='Foreground', insert_pos='suffix')
        mark_main(args)
        tree = read_tree(tmp_outfile, format='1', quoted_node_names=True, quiet=True)
        for leaf in tree.iter_leaves():
            if 'A1' in leaf.name or 'A2' in leaf.name:
                assert 'Foreground' in leaf.name

    def test_with_data_file(self, tmp_outfile, monkeypatch):
        monkeypatch.setattr(sys, 'argv', ['nwkit', 'mark'])
        infile = os.path.join(DATA_DIR, 'mark1', 'tree.nwk')
        if not os.path.exists(infile):
            pytest.skip('Test data not found')
        args = make_mark_args(
            infile=infile, outfile=tmp_outfile,
            pattern='A.*', target='leaf',
            insert_txt='Foreground', insert_pos='suffix',
        )
        mark_main(args)
        tree = read_tree(tmp_outfile, format='1', quoted_node_names=True, quiet=True)
        marked_leaves = [l for l in tree.iter_leaves() if 'Foreground' in l.name]
        assert len(marked_leaves) >= 1

    def test_wiki_codeml_clade_marking(self, tmp_nwk, tmp_outfile, monkeypatch):
        """Wiki example: mark clade with #1 for PAML codeml two-ratio mode.

        nwkit mark --infile input.nwk --pattern "A.*|B.*" --insert_txt "#1"
            --target clade --target_only_clade yes --outfile output.nwk

        Input:  (((A1:2.0,(B1:1.0,B2:1.0):1.0):1.0,(A2:1.0,C1:1.0):2.0):1.0,C2:4.0):0.25;
        Output: (((A1#1:2,(B1#1:1,B2#1:1)#1:1)#1:1,(A2#1:1,C1:1):2):1,C2:4):0.25;
        """
        monkeypatch.setattr(sys, 'argv', ['nwkit', 'mark'])
        path = tmp_nwk('(((A1:2.0,(B1:1.0,B2:1.0):1.0):1.0,(A2:1.0,C1:1.0):2.0):1.0,C2:4.0):0.25;')
        args = make_mark_args(
            infile=path, outfile=tmp_outfile,
            pattern='A.*|B.*', target='clade', target_only_clade=True,
            insert_txt='#1', insert_sep='', insert_pos='suffix',
        )
        mark_main(args)
        tree = read_tree(tmp_outfile, format='1', quoted_node_names=True, quiet=True)
        # A1, B1, B2 should be marked with #1 (they are in the target-only clade)
        for leaf in tree.iter_leaves():
            if leaf.name.startswith(('A1', 'B1', 'B2')):
                assert '#1' in leaf.name
        # C1, C2 should NOT be marked
        c_leaves = [l for l in tree.iter_leaves() if l.name.startswith('C')]
        for leaf in c_leaves:
            assert '#1' not in leaf.name

    def test_wiki_pipe_separated_regex(self, tmp_nwk, tmp_outfile, monkeypatch):
        """Test pipe-separated regex pattern matching multiple leaf groups."""
        monkeypatch.setattr(sys, 'argv', ['nwkit', 'mark'])
        path = tmp_nwk('(((A1:2,(B1:1,B2:1):1):1,(A2:1,C1:1):2):1,C2:4):0.25;')
        args = make_mark_args(
            infile=path, outfile=tmp_outfile,
            pattern='A.*|B.*', target='leaf',
            insert_txt='FG', insert_sep='', insert_pos='suffix',
        )
        mark_main(args)
        tree = read_tree(tmp_outfile, format='1', quoted_node_names=True, quiet=True)
        for leaf in tree.iter_leaves():
            if leaf.name.replace('FG', '').startswith(('A', 'B')):
                assert 'FG' in leaf.name
            else:
                assert 'FG' not in leaf.name

    def test_mark_all_mrca_clade(self, tmp_nwk, tmp_outfile, monkeypatch):
        """Test --target clade --target_only_clade no marks the entire clade
        containing all matched leaves (not just the sub-clade)."""
        monkeypatch.setattr(sys, 'argv', ['nwkit', 'mark'])
        path = tmp_nwk('(((A1:2,(B1:1,B2:1):1):1,(A2:1,C1:1):2):1,C2:4):0.25;')
        args = make_mark_args(
            infile=path, outfile=tmp_outfile,
            pattern='A.*', target='clade', target_only_clade=False,
            insert_txt='#1', insert_sep='', insert_pos='suffix',
        )
        mark_main(args)
        tree = read_tree(tmp_outfile, format='1', quoted_node_names=True, quiet=True)
        # A1 and A2 are not in a target-only clade, so the all_mrca_clade
        # should include all nodes between and up to their MRCA
        marked = [l for l in tree.iter_leaves() if '#1' in l.name]
        assert len(marked) >= 2

    def test_issue9_no_match_no_indexerror(self, tmp_nwk, tmp_outfile, monkeypatch):
        """Regression test for GitHub issue #9.

        nwkit mark used to crash with IndexError when the pattern matched
        no leaves in the tree. The fix guards target_leaves[0] access.
        """
        monkeypatch.setattr(sys, 'argv', ['nwkit', 'mark'])
        path = tmp_nwk('(((A1:2,(B1:1,B2:1):1):1,(A2:1,C1:1):2):1,C2:4):0.25;')
        args = make_mark_args(
            infile=path, outfile=tmp_outfile,
            pattern='NONEXISTENT_SPECIES', target='clade', target_only_clade=True,
            insert_txt='#1', insert_sep='', insert_pos='suffix',
        )
        # Should not raise IndexError
        mark_main(args)
        tree = read_tree(tmp_outfile, format='1', quoted_node_names=True, quiet=True)
        # No leaves should be marked
        for leaf in tree.iter_leaves():
            assert '#1' not in leaf.name

    def test_wiki_target_mrca_only(self, tmp_nwk, tmp_outfile, monkeypatch):
        """Wiki --target mrca: only the MRCA node of matched leaves is marked.

        Pattern B.*, target mrca, target_only_clade yes:
        Output: (((A1:2,(B1:1,B2:1)#1:1):1,(A2:1,C1:1):2):1,C2:4):0.25;

        Only the parent of B1/B2 gets #1; B1/B2 themselves are NOT marked.
        """
        monkeypatch.setattr(sys, 'argv', ['nwkit', 'mark'])
        path = tmp_nwk('(((A1:2,(B1:1,B2:1):1):1,(A2:1,C1:1):2):1,C2:4):0.25;')
        args = make_mark_args(
            infile=path, outfile=tmp_outfile,
            pattern='B.*', target='mrca', target_only_clade=True,
            insert_txt='#1', insert_sep='', insert_pos='suffix',
        )
        mark_main(args)
        tree = read_tree(tmp_outfile, format='1', quoted_node_names=True, quiet=True)
        # B1 and B2 leaves should NOT be marked
        for leaf in tree.iter_leaves():
            if leaf.name.startswith('B'):
                assert '#1' not in leaf.name
        # The MRCA internal node should be marked
        marked_internal = [n for n in tree.traverse()
                           if not n.is_leaf() and '#1' in n.name]
        assert len(marked_internal) >= 1

    def test_wiki_target_leaf_only(self, tmp_nwk, tmp_outfile, monkeypatch):
        """Wiki --target leaf: only matched leaves are marked, not internal nodes.

        Pattern B.*, target leaf:
        Output: (((A1:2,(B1#1:1,B2#1:1):1):1,(A2:1,C1:1):2):1,C2:4):0.25;
        """
        monkeypatch.setattr(sys, 'argv', ['nwkit', 'mark'])
        path = tmp_nwk('(((A1:2,(B1:1,B2:1):1):1,(A2:1,C1:1):2):1,C2:4):0.25;')
        args = make_mark_args(
            infile=path, outfile=tmp_outfile,
            pattern='B.*', target='leaf', target_only_clade=True,
            insert_txt='#1', insert_sep='', insert_pos='suffix',
        )
        mark_main(args)
        tree = read_tree(tmp_outfile, format='1', quoted_node_names=True, quiet=True)
        # B1, B2 should be marked
        b_leaves = [l for l in tree.iter_leaves() if l.name.startswith('B')]
        assert len(b_leaves) == 2
        for leaf in b_leaves:
            assert '#1' in leaf.name
        # Non-B leaves should NOT be marked
        other_leaves = [l for l in tree.iter_leaves() if not l.name.startswith('B')]
        for leaf in other_leaves:
            assert '#1' not in leaf.name
        # Internal nodes should NOT be marked
        for node in tree.traverse():
            if not node.is_leaf():
                assert '#1' not in node.name

    def test_wiki_clade_branch_lengths_preserved(self, tmp_nwk, tmp_outfile, monkeypatch):
        """Marking should not alter branch lengths."""
        monkeypatch.setattr(sys, 'argv', ['nwkit', 'mark'])
        path = tmp_nwk('(((A1:2,(B1:1,B2:1):1):1,(A2:1,C1:1):2):1,C2:4):0.25;')
        args = make_mark_args(
            infile=path, outfile=tmp_outfile,
            pattern='B.*', target='clade', target_only_clade=True,
            insert_txt='#1', insert_sep='', insert_pos='suffix',
        )
        mark_main(args)
        tree = read_tree(tmp_outfile, format='1', quoted_node_names=True, quiet=True)
        leaves = {}
        for l in tree.iter_leaves():
            clean_name = l.name.replace('#1', '')
            leaves[clean_name] = l.dist
        assert abs(leaves['A1'] - 2.0) < 1e-6
        assert abs(leaves['B1'] - 1.0) < 1e-6
        assert abs(leaves['B2'] - 1.0) < 1e-6
        assert abs(leaves['A2'] - 1.0) < 1e-6
        assert abs(leaves['C1'] - 1.0) < 1e-6
        assert abs(leaves['C2'] - 4.0) < 1e-6
