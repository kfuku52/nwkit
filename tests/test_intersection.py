import os
import stat
import pytest
from ete4 import Tree

from nwkit.fasta import FastaRecord, parse_fasta
from nwkit import intersection as intersection_mod
from nwkit.intersection import (
    get_leaf_names,
    get_seq_names,
    match_complete,
    match_prefix,
    match_backward,
    get_remove_names,
    intersection_main,
)
from nwkit.util import read_tree
from tests.helpers import make_args


class TestGetLeafNames:
    def test_unique_names(self):
        tree = Tree('((A,B),(C,D));', parser=1)
        names = get_leaf_names(tree)
        assert set(names) == {'A', 'B', 'C', 'D'}

    def test_duplicate_raises(self):
        tree = Tree('((A,A),(C,D));', parser=1)
        with pytest.raises(ValueError, match='unique'):
            get_leaf_names(tree)

class TestGetSeqNames:
    def test_unique_names(self):
        seqs = [
            FastaRecord(name='A', raw='>A\nATG\n'),
            FastaRecord(name='B', raw='>B\nATG\n'),
        ]
        names = get_seq_names(seqs)
        assert names == ['A', 'B']

    def test_duplicate_raises(self):
        seqs = [
            FastaRecord(name='A', raw='>A\nATG\n'),
            FastaRecord(name='A', raw='>A\nATG\n'),
        ]
        with pytest.raises(ValueError, match='unique'):
            get_seq_names(seqs)


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

    def test_prefix_mode_with_empty_label_raises(self):
        arr1 = ['ABC_001', None, 'GHI_003']
        arr2 = ['ABC', 'GHI']
        with pytest.raises(ValueError, match='match prefix'):
            get_remove_names(arr1, arr2, 'prefix')

    def test_backward_mode_with_empty_label_raises(self):
        arr1 = ['pre_ABC', 'pre_DEF']
        arr2 = ['ABC', None]
        with pytest.raises(ValueError, match='match backward'):
            get_remove_names(arr1, arr2, 'backward')


class TestIntersectionMain:
    @pytest.mark.skipif(
        os.name == 'nt',
        reason='requires POSIX open-file replacement semantics',
    )
    def test_staging_writer_uses_open_descriptor_not_replaced_path(
        self,
        monkeypatch,
        tmp_path,
    ):
        victim = tmp_path / 'victim.txt'
        victim.write_text('unchanged\n')
        real_mkstemp = intersection_mod.tempfile.mkstemp

        def replace_created_path(*args, **kwargs):
            fd, path = real_mkstemp(*args, **kwargs)
            os.remove(path)
            os.symlink(victim, path)
            return fd, path

        monkeypatch.setattr(
            intersection_mod.tempfile,
            'mkstemp',
            replace_created_path,
        )

        with pytest.raises(RuntimeError, match='staging file was replaced'):
            intersection_mod._stage_file(
                '-',
                lambda handle: handle.write('staged output\n'),
            )

        assert victim.read_text() == 'unchanged\n'

    @pytest.mark.skipif(
        os.name == 'nt',
        reason='requires POSIX open-file replacement semantics',
    )
    def test_stdout_commit_reads_validated_descriptor_not_replaced_path(
        self,
        tmp_path,
        capsys,
    ):
        secret = tmp_path / 'secret.txt'
        secret.write_text('SECRET-CONTENT\n')
        staged = intersection_mod._stage_file(
            '-',
            lambda handle: handle.write('expected output\n'),
        )
        os.remove(staged['path'])
        os.symlink(secret, staged['path'])

        intersection_mod._commit_staged_outputs([('-', staged)])

        assert capsys.readouterr().out == 'expected output\n'
        assert secret.read_text() == 'SECRET-CONTENT\n'

    @pytest.mark.skipif(
        not hasattr(os, 'mkfifo'),
        reason='FIFOs are unavailable on this platform',
    )
    def test_existing_fifo_output_is_rejected_without_opening_it(self, tmp_path):
        nwk_path = tmp_path / 'tree.nwk'
        fifo_path = tmp_path / 'output.nwk'
        nwk_path.write_text('(A:1,B:1);')
        os.mkfifo(fifo_path)
        args = make_args(
            infile=str(nwk_path), infile2='(A:1,B:1);',
            outfile=str(fifo_path), format2='auto', seqin='',
            seqout='', seqformat='fasta', match='complete',
        )

        with pytest.raises(ValueError, match='regular file'):
            intersection_main(args)

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
        leaf_names = set(tree.leaf_names())
        assert leaf_names == {'A', 'B', 'C'}

    def test_requires_second_input_or_sequences(self, tmp_nwk, tmp_outfile):
        path1 = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);', 'tree1.nwk')
        args = make_args(
            infile=path1, infile2='', outfile=tmp_outfile,
            format2='auto', seqin='', seqout='', seqformat='fasta',
            match='complete',
        )
        with pytest.raises(ValueError, match='infile2'):
            intersection_main(args)

    def test_rejects_both_second_tree_and_sequences_before_writing(self, tmp_path):
        tree1 = tmp_path / 'tree1.nwk'
        tree2 = tmp_path / 'tree2.nwk'
        seqin = tmp_path / 'input.fasta'
        outfile = tmp_path / 'output.nwk'
        seqout = tmp_path / 'output.fasta'
        tree1.write_text('(A:1,B:1);')
        tree2.write_text('(A:1,C:1);')
        seqin.write_text('>A\nATG\n')
        outfile.write_text('original tree\n')
        seqout.write_text('original sequences\n')
        args = make_args(
            infile=str(tree1), infile2=str(tree2), outfile=str(outfile),
            format2='auto', seqin=str(seqin), seqout=str(seqout),
            seqformat='fasta', match='complete',
        )

        with pytest.raises(ValueError, match='exactly one'):
            intersection_main(args)

        assert outfile.read_text() == 'original tree\n'
        assert seqout.read_text() == 'original sequences\n'

    def test_no_overlap_raises(self, tmp_nwk, tmp_outfile):
        path1 = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);', 'tree1.nwk')
        path2 = tmp_nwk('((E:1,F:1):1,(G:1,H:1):1);', 'tree2.nwk')
        args = make_args(
            infile=path1, infile2=path2, outfile=tmp_outfile,
            format2='auto', seqin='', seqout='', seqformat='fasta',
            match='complete',
        )
        with pytest.raises(ValueError, match='No overlap'):
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
        leaf_names = set(tree.leaf_names())
        assert leaf_names == {'A', 'C', 'D', 'F'}
        assert os.path.exists(out_seq)

    def test_tree_seq_no_overlap_leaves_both_outputs_unchanged(self, tmp_path):
        nwk_path = tmp_path / 'tree.nwk'
        seq_path = tmp_path / 'seq.fasta'
        out_tree = tmp_path / 'out.nwk'
        out_seq = tmp_path / 'out.fasta'
        nwk_path.write_text('(A:1,B:1);')
        seq_path.write_text('>X\nATG\n>Y\nATG\n')
        out_tree.write_text('original tree\n')
        out_seq.write_text('original sequences\n')
        args = make_args(
            infile=str(nwk_path), infile2='', outfile=str(out_tree),
            seqin=str(seq_path), seqout=str(out_seq), seqformat='fasta',
            format2='auto', match='complete',
        )

        with pytest.raises(ValueError, match='No overlap'):
            intersection_main(args)

        assert out_tree.read_text() == 'original tree\n'
        assert out_seq.read_text() == 'original sequences\n'

    def test_tree_seq_commit_failure_restores_both_outputs(self, monkeypatch, tmp_path):
        nwk_path = tmp_path / 'tree.nwk'
        seq_path = tmp_path / 'seq.fasta'
        out_tree = tmp_path / 'out.nwk'
        out_seq = tmp_path / 'out.fasta'
        nwk_path.write_text('(A:1,B:1);')
        seq_path.write_text('>A\nATG\n>B\nATG\n')
        out_tree.write_text('original tree\n')
        out_seq.write_text('original sequences\n')
        args = make_args(
            infile=str(nwk_path), infile2='', outfile=str(out_tree),
            seqin=str(seq_path), seqout=str(out_seq), seqformat='fasta',
            format2='auto', match='complete',
        )
        real_replace = os.replace
        replace_calls = 0

        def fail_second_replace(source, target):
            nonlocal replace_calls
            replace_calls += 1
            if replace_calls == 2:
                raise OSError('simulated second-output commit failure')
            return real_replace(source, target)

        monkeypatch.setattr(intersection_mod.os, 'replace', fail_second_replace)
        with pytest.raises(OSError, match='second-output commit failure'):
            intersection_main(args)

        assert out_tree.read_text() == 'original tree\n'
        assert out_seq.read_text() == 'original sequences\n'
        leftovers = [
            path.name
            for path in tmp_path.iterdir()
            if path.name.startswith('.out.')
        ]
        assert leftovers == []

    @pytest.mark.skipif(
        os.name == 'nt',
        reason='POSIX file modes are unavailable on Windows',
    )
    def test_tree_seq_atomic_replace_preserves_existing_modes(self, tmp_path):
        nwk_path = tmp_path / 'tree.nwk'
        seq_path = tmp_path / 'seq.fasta'
        out_tree = tmp_path / 'out.nwk'
        out_seq = tmp_path / 'out.fasta'
        nwk_path.write_text('(A:1,B:1);')
        seq_path.write_text('>A\nATG\n>B\nATG\n')
        out_tree.write_text('old tree\n')
        out_seq.write_text('old sequences\n')
        out_tree.chmod(0o640)
        out_seq.chmod(0o604)
        args = make_args(
            infile=str(nwk_path), infile2='', outfile=str(out_tree),
            seqin=str(seq_path), seqout=str(out_seq), seqformat='fasta',
            format2='auto', match='complete',
        )

        intersection_main(args)

        assert stat.S_IMODE(out_tree.stat().st_mode) == 0o640
        assert stat.S_IMODE(out_seq.stat().st_mode) == 0o604

    @pytest.mark.skipif(
        os.name == 'nt',
        reason='POSIX file modes are unavailable on Windows',
    )
    def test_tree_seq_new_outputs_honor_process_umask(self, tmp_path):
        nwk_path = tmp_path / 'tree.nwk'
        seq_path = tmp_path / 'seq.fasta'
        out_tree = tmp_path / 'out.nwk'
        out_seq = tmp_path / 'out.fasta'
        nwk_path.write_text('(A:1,B:1);')
        seq_path.write_text('>A\nATG\n>B\nATG\n')
        args = make_args(
            infile=str(nwk_path), infile2='', outfile=str(out_tree),
            seqin=str(seq_path), seqout=str(out_seq), seqformat='fasta',
            format2='auto', match='complete',
        )
        previous_umask = os.umask(0o027)
        try:
            intersection_main(args)
        finally:
            os.umask(previous_umask)

        assert stat.S_IMODE(out_tree.stat().st_mode) == 0o640
        assert stat.S_IMODE(out_seq.stat().st_mode) == 0o640

    def test_tree_seq_outputs_follow_existing_symlinks(self, tmp_path):
        nwk_path = tmp_path / 'tree.nwk'
        seq_path = tmp_path / 'seq.fasta'
        tree_target = tmp_path / 'tree-target.nwk'
        seq_target = tmp_path / 'seq-target.fasta'
        out_tree = tmp_path / 'out.nwk'
        out_seq = tmp_path / 'out.fasta'
        nwk_path.write_text('(A:1,B:1);')
        seq_path.write_text('>A\nATG\n>B\nATG\n')
        tree_target.write_text('old tree\n')
        seq_target.write_text('old sequences\n')
        out_tree.symlink_to(tree_target.name)
        out_seq.symlink_to(seq_target.name)
        args = make_args(
            infile=str(nwk_path), infile2='', outfile=str(out_tree),
            seqin=str(seq_path), seqout=str(out_seq), seqformat='fasta',
            format2='auto', match='complete',
        )

        intersection_main(args)

        assert out_tree.is_symlink()
        assert out_seq.is_symlink()
        assert tree_target.read_text().endswith(';')
        assert seq_target.read_text() == '>A\nATG\n>B\nATG\n'

    def test_tree_seq_intersection_preserves_selected_fasta_records(self, tmp_path):
        nwk_path = tmp_path / 'tree.nwk'
        nwk_path.write_text('((A:1,B:1):1,(C:1,D:1):1);')
        seq_path = tmp_path / 'seq.fasta'
        seq_path.write_text(
            '>A retained description\nAC\nGT\n'
            '>X removed description\nNNNN\n'
            '>C another retained description\nTGCA\n'
        )
        out_tree = str(tmp_path / 'out.nwk')
        out_seq = str(tmp_path / 'out.fasta')
        args = make_args(
            infile=str(nwk_path), infile2='', outfile=out_tree,
            seqin=str(seq_path), seqout=out_seq, seqformat='fasta',
            format2='auto', match='complete',
        )
        intersection_main(args)
        assert (tmp_path / 'out.fasta').read_text() == (
            '>A retained description\nAC\nGT\n'
            '>C another retained description\nTGCA\n'
        )

    def test_tree_seq_intersection_with_empty_seqout_writes_to_stdout(self, tmp_path, capsys):
        nwk_path = tmp_path / 'tree.nwk'
        nwk_path.write_text('((A:1,B:1):1,(C:1,D:1):1);')
        seq_path = tmp_path / 'seq.fasta'
        seq_path.write_text('>A\nATG\n>B\nATG\n>C\nATG\n>D\nATG\n')
        out_tree = str(tmp_path / 'out.nwk')
        args = make_args(
            infile=str(nwk_path), infile2='', outfile=out_tree,
            seqin=str(seq_path), seqout='', seqformat='fasta',
            format2='auto', match='complete',
        )
        intersection_main(args)
        captured = capsys.readouterr()
        assert '>A' in captured.out
        assert os.path.exists(out_tree)

    def test_tree_and_sequence_outputs_cannot_both_be_stdout(self, tmp_path):
        nwk_path = tmp_path / 'tree.nwk'
        nwk_path.write_text('((A:1,B:1):1,(C:1,D:1):1);')
        seq_path = tmp_path / 'seq.fasta'
        seq_path.write_text('>A\nATG\n>B\nATG\n>C\nATG\n>D\nATG\n')
        args = make_args(
            infile=str(nwk_path), infile2='', outfile='-',
            seqin=str(seq_path), seqout='', seqformat='fasta',
            format2='auto', match='complete',
        )
        with pytest.raises(ValueError, match='both be written to stdout'):
            intersection_main(args)

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
        assert set(tree.leaf_names()) == {'A', 'C', 'D', 'F'}
        # Output FASTA should contain all 4 sequences (all are in the tree).
        with open(out_seq, newline='') as fh:
            out_seq_names = {record.name for record in parse_fasta(fh)}
        assert out_seq_names == {'A', 'C', 'D', 'F'}
