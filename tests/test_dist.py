import os
import pytest
from ete4 import Tree

from nwkit.dist import dist_main
from nwkit.util import read_tree
from tests.helpers import make_args, DATA_DIR


class TestDistMain:
    def test_identical_trees(self, tmp_nwk, tmp_outfile):
        nwk = '((A:1,B:1):1,(C:1,D:1):1);'
        path1 = tmp_nwk(nwk, 'tree1.nwk')
        path2 = tmp_nwk(nwk, 'tree2.nwk')
        args = make_args(
            infile=path1, infile2=path2, outfile=tmp_outfile,
            dist='RF', format2='auto',
        )
        dist_main(args)
        with open(tmp_outfile) as f:
            content = f.read()
        assert '0' in content

    def test_file_output_has_header_and_data_on_separate_lines(self, tmp_nwk, tmp_outfile):
        nwk = '((A:1,B:1):1,(C:1,D:1):1);'
        path1 = tmp_nwk(nwk, 'tree1.nwk')
        path2 = tmp_nwk(nwk, 'tree2.nwk')
        args = make_args(
            infile=path1, infile2=path2, outfile=tmp_outfile,
            dist='RF', format2='auto',
        )
        dist_main(args)
        with open(tmp_outfile) as f:
            lines = f.read().splitlines()
        assert lines[0] == 'rf_dist\tmax_rf_dist'
        assert lines[1] == '0\t4'

    def test_different_trees(self, tmp_nwk, capsys):
        path1 = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);', 'tree1.nwk')
        path2 = tmp_nwk('((A:1,C:1):1,(B:1,D:1):1);', 'tree2.nwk')
        args = make_args(
            infile=path1, infile2=path2, outfile='-',
            dist='RF', format2='auto',
        )
        dist_main(args)
        captured = capsys.readouterr()
        assert 'rf_dist' in captured.out
        # These trees have different topologies, RF > 0
        lines = captured.out.strip().split('\n')
        vals = lines[1].split('\t')
        rf_dist = int(vals[0])
        assert rf_dist > 0

    def test_mismatched_leaves_raises(self, tmp_nwk):
        path1 = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);', 'tree1.nwk')
        path2 = tmp_nwk('((A:1,B:1):1,(C:1,E:1):1);', 'tree2.nwk')
        args = make_args(
            infile=path1, infile2=path2, outfile='-',
            dist='RF', format2='auto',
        )
        with pytest.raises(ValueError, match='Leaf name'):
            dist_main(args)

    def test_unsupported_dist_raises(self, tmp_nwk):
        nwk = '((A:1,B:1):1,(C:1,D:1):1);'
        path1 = tmp_nwk(nwk, 'tree1.nwk')
        path2 = tmp_nwk(nwk, 'tree2.nwk')
        args = make_args(
            infile=path1, infile2=path2, outfile='-',
            dist='NJ', format2='auto',
        )
        with pytest.raises(ValueError):
            dist_main(args)

    def test_missing_infile2_raises(self, tmp_nwk):
        nwk = '((A:1,B:1):1,(C:1,D:1):1);'
        path1 = tmp_nwk(nwk, 'tree1.nwk')
        args = make_args(
            infile=path1, infile2='', outfile='-',
            dist='RF', format2='auto',
        )
        with pytest.raises(ValueError, match='infile2'):
            dist_main(args)

    def test_duplicate_leaf_names_raise_clear_error(self, tmp_nwk):
        path1 = tmp_nwk('((A:1,A:1):1,B:1);', 'tree1.nwk')
        path2 = tmp_nwk('((A:1,B:1):1,B:1);', 'tree2.nwk')
        args = make_args(
            infile=path1, infile2=path2, outfile='-',
            dist='RF', format2='auto',
        )
        with pytest.raises(ValueError, match='Duplicated leaf labels'):
            dist_main(args)

    def test_empty_leaf_labels_raise_clear_error(self, tmp_nwk):
        path1 = tmp_nwk('((A:1,:1):1,B:1);', 'tree1.nwk')
        path2 = tmp_nwk('((A:1,B:1):1,:1);', 'tree2.nwk')
        args = make_args(
            infile=path1, infile2=path2, outfile='-',
            dist='RF', format2='auto',
        )
        with pytest.raises(ValueError, match='Empty leaf labels'):
            dist_main(args)

    def test_with_data_files(self, capsys):
        infile1 = os.path.join(DATA_DIR, 'dist1', 'input1.nwk')
        infile2 = os.path.join(DATA_DIR, 'dist1', 'input2.nwk')
        if not os.path.exists(infile1):
            pytest.skip('Test data not found')
        args = make_args(
            infile=infile1, infile2=infile2, outfile='-',
            dist='RF', format2='auto',
        )
        dist_main(args)
        captured = capsys.readouterr()
        assert 'rf_dist' in captured.out

    def test_wiki_exact_rf_values(self, capsys):
        """Wiki example: RF distance between dist1 trees should be 4, max 54."""
        infile1 = os.path.join(DATA_DIR, 'dist1', 'input1.nwk')
        infile2 = os.path.join(DATA_DIR, 'dist1', 'input2.nwk')
        if not os.path.exists(infile1):
            pytest.skip('Test data not found')
        args = make_args(
            infile=infile1, infile2=infile2, outfile='-',
            dist='RF', format2='auto',
        )
        dist_main(args)
        captured = capsys.readouterr()
        lines = captured.out.strip().split('\n')
        assert lines[0] == 'rf_dist\tmax_rf_dist'
        vals = lines[1].split('\t')
        assert int(vals[0]) == 4
        assert int(vals[1]) == 54

    def test_rf_dist_zero_for_identical(self, tmp_nwk, capsys):
        """RF distance should be 0 for identical trees."""
        nwk = '(((A:1,B:1):1,C:1):1,(D:1,E:1):1);'
        path1 = tmp_nwk(nwk, 'tree1.nwk')
        path2 = tmp_nwk(nwk, 'tree2.nwk')
        args = make_args(
            infile=path1, infile2=path2, outfile='-',
            dist='RF', format2='auto',
        )
        dist_main(args)
        captured = capsys.readouterr()
        lines = captured.out.strip().split('\n')
        vals = lines[1].split('\t')
        assert int(vals[0]) == 0

    def test_unrooted_tree_raises_clear_error(self, tmp_nwk):
        path1 = tmp_nwk('(A:1,B:1,C:1);', 'tree1.nwk')
        path2 = tmp_nwk('(A:1,B:1,C:1);', 'tree2.nwk')
        args = make_args(
            infile=path1, infile2=path2, outfile='-',
            dist='RF', format2='auto',
        )
        with pytest.raises(ValueError, match='requires rooted trees'):
            dist_main(args)
