import os
import pytest
from ete3 import TreeNode

from nwkit.intersection import (
    get_leaf_names,
    get_seq_names,
    get_new_seqs,
    match_complete,
    match_prefix,
    match_backward,
    get_remove_names,
    intersection_main,
)
from nwkit.util import read_tree
from tests.helpers import make_args, DATA_DIR


class TestGetLeafNames:
    def test_unique_names(self):
        tree = TreeNode(newick='((A,B),(C,D));', format=1)
        names = get_leaf_names(tree)
        assert set(names) == {'A', 'B', 'C', 'D'}

    def test_duplicate_raises(self):
        tree = TreeNode(newick='((A,A),(C,D));', format=1)
        with pytest.raises(AssertionError, match='unique'):
            get_leaf_names(tree)


class TestMatchFunctions:
    def test_complete_match(self):
        assert match_complete('abc', 'abc') is True
        assert match_complete('abc', 'ab') is False

    def test_prefix_match(self):
        assert match_prefix('abcdef', 'abc') is True
        assert match_prefix('abc', 'abcdef') is True
        assert match_prefix('abc', 'xyz') is False

    def test_backward_match(self):
        assert match_backward('abcdef', 'def') is True
        assert match_backward('def', 'abcdef') is True
        assert match_backward('abc', 'xyz') is False


class TestGetRemoveNames:
    def test_complete(self):
        arr1 = ['A', 'B', 'C', 'D']
        arr2 = ['A', 'B', 'C']
        result = get_remove_names(arr1, arr2, 'complete')
        assert result == ['D']

    def test_all_matched(self):
        arr1 = ['A', 'B', 'C']
        arr2 = ['A', 'B', 'C']
        result = get_remove_names(arr1, arr2, 'complete')
        assert result == []

    def test_none_matched(self):
        arr1 = ['A', 'B', 'C']
        arr2 = ['D', 'E', 'F']
        result = get_remove_names(arr1, arr2, 'complete')
        assert result == ['A', 'B', 'C']

    def test_prefix_mode(self):
        arr1 = ['ABC_001', 'DEF_002', 'GHI_003']
        arr2 = ['ABC', 'DEF']
        result = get_remove_names(arr1, arr2, 'prefix')
        assert result == ['GHI_003']

    def test_backward_mode(self):
        arr1 = ['pre_ABC', 'pre_DEF', 'pre_GHI']
        arr2 = ['ABC', 'DEF']
        result = get_remove_names(arr1, arr2, 'backward')
        assert result == ['pre_GHI']


class TestIntersectionMain:
    def test_tree_tree_intersection(self, tmp_nwk, tmp_outfile):
        path1 = tmp_nwk('(((A:1,B:1):1,(C:1,D:1):1):1,(E:1,F:1):1):0;', 'tree1.nwk')
        path2 = tmp_nwk('(((A:1,B:1):1,C:1):1,G:1):0;', 'tree2.nwk')
        args = make_args(
            infile=path1, infile2=path2, outfile=tmp_outfile,
            format2='auto', seqin='', seqout='', seqformat='fasta',
            match='complete',
        )
        intersection_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        leaf_names = set(tree.get_leaf_names())
        assert leaf_names == {'A', 'B', 'C'}

    def test_no_overlap_raises(self, tmp_nwk, tmp_outfile):
        path1 = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);', 'tree1.nwk')
        path2 = tmp_nwk('((E:1,F:1):1,(G:1,H:1):1);', 'tree2.nwk')
        args = make_args(
            infile=path1, infile2=path2, outfile=tmp_outfile,
            format2='auto', seqin='', seqout='', seqformat='fasta',
            match='complete',
        )
        with pytest.raises(AssertionError, match='No overlap'):
            intersection_main(args)

    def test_tree_seq_intersection(self, tmp_path):
        nwk_path = tmp_path / 'tree.nwk'
        nwk_path.write_text('(((A:1,B:1):1,(C:1,D:1):1):1,(E:1,F:1):1):0;')
        seq_path = tmp_path / 'seq.fasta'
        seq_path.write_text('>A\nATG\n>C\nATG\n>D\nATG\n>F\nATG\n')
        out_tree = str(tmp_path / 'out.nwk')
        out_seq = str(tmp_path / 'out.fasta')
        args = make_args(
            infile=str(nwk_path), infile2='', outfile=out_tree,
            seqin=str(seq_path), seqout=out_seq, seqformat='fasta',
            format2='auto', match='complete',
        )
        intersection_main(args)
        tree = read_tree(out_tree, format='auto', quoted_node_names=True, quiet=True)
        leaf_names = set(tree.get_leaf_names())
        assert leaf_names == {'A', 'C', 'D', 'F'}
        assert os.path.exists(out_seq)

    def test_with_data_files(self, tmp_outfile):
        infile1 = os.path.join(DATA_DIR, 'intersection3', 'input1.nwk')
        infile2 = os.path.join(DATA_DIR, 'intersection3', 'input2.nwk')
        if not os.path.exists(infile1):
            pytest.skip('Test data not found')
        args = make_args(
            infile=infile1, infile2=infile2, outfile=tmp_outfile,
            format2='auto', seqin='', seqout='', seqformat='fasta',
            match='complete',
        )
        intersection_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        leaf_names = set(tree.get_leaf_names())
        assert leaf_names == {'A', 'B', 'C'}

    def test_wiki_tree_seq_intersection(self, tmp_path):
        """Wiki example: intersection between tree and FASTA alignment.

        input.nwk: (((A:1,B:1):1,(C:1,D:1):1):1,(E:1,F:1):1):0;
        input.fasta: >A ATGCAATAA >C ATGCATTAA >D ATGGATTAA >F ATGAGGTAA
        output.nwk: tree pruned to A, C, D, F
        output.fasta: same as input (all seqs are in tree)
        """
        nwk_path = tmp_path / 'input.nwk'
        nwk_path.write_text('(((A:1,B:1):1,(C:1,D:1):1):1,(E:1,F:1):1):0;')
        seq_path = tmp_path / 'input.fasta'
        seq_path.write_text('>A\nATGCAATAA\n>C\nATGCATTAA\n>D\nATGGATTAA\n>F\nATGAGGTAA\n')
        out_tree = str(tmp_path / 'output.nwk')
        out_seq = str(tmp_path / 'output.fasta')
        args = make_args(
            infile=str(nwk_path), infile2='', outfile=out_tree,
            seqin=str(seq_path), seqout=out_seq, seqformat='fasta',
            format2='auto', match='complete',
        )
        intersection_main(args)
        tree = read_tree(out_tree, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.get_leaf_names()) == {'A', 'C', 'D', 'F'}
        # Output FASTA should contain all 4 sequences (all are in the tree)
        import Bio.SeqIO
        out_records = list(Bio.SeqIO.parse(out_seq, 'fasta'))
        out_seq_names = {r.name for r in out_records}
        assert out_seq_names == {'A', 'C', 'D', 'F'}

    def test_wiki_data_intersection4(self, tmp_path):
        """Wiki example with data/intersection4 files."""
        infile = os.path.join(DATA_DIR, 'intersection4', 'input.nwk')
        seqin = os.path.join(DATA_DIR, 'intersection4', 'input.fasta')
        if not os.path.exists(infile):
            pytest.skip('Test data not found')
        out_tree = str(tmp_path / 'output.nwk')
        out_seq = str(tmp_path / 'output.fasta')
        args = make_args(
            infile=infile, infile2='', outfile=out_tree,
            seqin=seqin, seqout=out_seq, seqformat='fasta',
            format2='auto', match='complete',
        )
        intersection_main(args)
        tree = read_tree(out_tree, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.get_leaf_names()) == {'A', 'C', 'D', 'F'}
        import Bio.SeqIO
        out_records = list(Bio.SeqIO.parse(out_seq, 'fasta'))
        assert len(out_records) == 4
