import csv
import math

import pytest

from nwkit.dist import dist_main
from tests.helpers import make_args


def _read_tsv(path):
    with open(path, newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


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


class TestDistanceMetrics:
    def test_all_is_default_and_uses_stable_long_form_rows(self, tmp_nwk, tmp_outfile):
        path1 = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);', 'tree1.nwk')
        path2 = tmp_nwk('((A:1,C:1):1,(B:1,D:1):1);', 'tree2.nwk')
        args = make_args(
            infile=path1,
            infile2=path2,
            outfile=tmp_outfile,
            format2='auto',
            metric=None,
            dist=None,
            comparison='rooted',
        )

        dist_main(args)

        rows = _read_tsv(tmp_outfile)
        assert list(rows[0]) == [
            'metric', 'comparison', 'num_taxa', 'distance', 'max_distance',
        ]
        assert [row['metric'] for row in rows] == [
            'rf',
            'normalized-rf',
            'weighted-rf',
            'branch-score',
            'path-topological',
            'path-length',
        ]
        assert [float(row['distance']) for row in rows] == pytest.approx([
            4.0, 1.0, 4.0, 2.0, 2.0, 4.0,
        ])
        assert [row['comparison'] for row in rows] == [
            'rooted',
            'rooted',
            'rooted',
            'rooted',
            'root-independent',
            'root-independent',
        ]
        assert {row['num_taxa'] for row in rows} == {'4'}
        assert [row['max_distance'] for row in rows] == [
            '4', '1.0', '', '', '', '',
        ]

    def test_branch_length_metrics_match_l1_l2_and_path_definitions(
        self,
        tmp_nwk,
        tmp_outfile,
    ):
        path1 = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);', 'tree1.nwk')
        path2 = tmp_nwk('((A:2,B:3):4,(C:5,D:6):7);', 'tree2.nwk')
        args = make_args(
            infile=path1,
            infile2=path2,
            outfile=tmp_outfile,
            format2='auto',
            metric=['all'],
            dist=None,
            comparison='rooted',
        )

        dist_main(args)

        rows = {row['metric']: row for row in _read_tsv(tmp_outfile)}
        assert float(rows['rf']['distance']) == 0.0
        assert float(rows['normalized-rf']['distance']) == 0.0
        assert float(rows['weighted-rf']['distance']) == pytest.approx(21.0)
        assert float(rows['branch-score']['distance']) == pytest.approx(math.sqrt(91))
        assert float(rows['path-topological']['distance']) == 0.0
        assert float(rows['path-length']['distance']) == pytest.approx(math.sqrt(992))

    def test_unrooted_comparison_suppresses_the_binary_root_edge(
        self,
        tmp_nwk,
        tmp_outfile,
    ):
        path1 = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);', 'tree1.nwk')
        path2 = tmp_nwk('((A:1,C:1):1,(B:1,D:1):1);', 'tree2.nwk')
        args = make_args(
            infile=path1,
            infile2=path2,
            outfile=tmp_outfile,
            format2='auto',
            metric=['rf,normalized-rf,weighted-rf,branch-score'],
            dist=None,
            comparison='unrooted',
        )

        dist_main(args)

        rows = {row['metric']: row for row in _read_tsv(tmp_outfile)}
        assert float(rows['rf']['distance']) == 2.0
        assert float(rows['rf']['max_distance']) == 2.0
        assert float(rows['normalized-rf']['distance']) == 1.0
        assert float(rows['weighted-rf']['distance']) == 4.0
        assert float(rows['branch-score']['distance']) == pytest.approx(math.sqrt(8))

    def test_unrooted_branch_metrics_are_invariant_to_root_position(
        self,
        tmp_nwk,
        tmp_outfile,
    ):
        path1 = tmp_nwk('((A:1,B:1):2,(C:1,D:1):3);', 'tree1.nwk')
        path2 = tmp_nwk('(A:0.4,(B:1,(C:1,D:1):5):0.6);', 'tree2.nwk')
        args = make_args(
            infile=path1,
            infile2=path2,
            outfile=tmp_outfile,
            format2='auto',
            metric=['all'],
            dist=None,
            comparison='unrooted',
        )

        dist_main(args)

        assert [float(row['distance']) for row in _read_tsv(tmp_outfile)] == pytest.approx([
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        ])

    def test_single_taxon_rf_maximum_is_zero(self, tmp_nwk, tmp_outfile):
        path1 = tmp_nwk('A;', 'tree1.nwk')
        path2 = tmp_nwk('A;', 'tree2.nwk')
        args = make_args(
            infile=path1,
            infile2=path2,
            outfile=tmp_outfile,
            format2='auto',
            metric=['rf,normalized-rf,path-topological'],
            dist=None,
            comparison='rooted',
        )

        dist_main(args)

        rows = {row['metric']: row for row in _read_tsv(tmp_outfile)}
        assert rows['rf']['distance'] == '0'
        assert rows['rf']['max_distance'] == '0'
        assert rows['normalized-rf']['distance'] == '0.0'

    def test_repeatable_metrics_are_deduplicated(self, tmp_nwk, tmp_outfile):
        nwk = '((A:1,B:1):1,(C:1,D:1):1);'
        path1 = tmp_nwk(nwk, 'tree1.nwk')
        path2 = tmp_nwk(nwk, 'tree2.nwk')
        args = make_args(
            infile=path1,
            infile2=path2,
            outfile=tmp_outfile,
            format2='auto',
            metric=['rf,path-topological', 'rf'],
            dist=None,
            comparison='rooted',
        )

        dist_main(args)

        assert [row['metric'] for row in _read_tsv(tmp_outfile)] == [
            'rf', 'path-topological',
        ]

    def test_path_topological_accepts_unrooted_trees_without_lengths(
        self,
        tmp_nwk,
        tmp_outfile,
    ):
        path1 = tmp_nwk('(A,B,C,D);', 'tree1.nwk')
        path2 = tmp_nwk('(A,B,C,D);', 'tree2.nwk')
        args = make_args(
            infile=path1,
            infile2=path2,
            outfile=tmp_outfile,
            format2='auto',
            metric=['path-topological'],
            dist=None,
            comparison='rooted',
        )

        dist_main(args)

        rows = _read_tsv(tmp_outfile)
        assert rows[0]['comparison'] == 'root-independent'
        assert float(rows[0]['distance']) == 0.0

    def test_default_all_rejects_missing_branch_lengths(self, tmp_nwk):
        path1 = tmp_nwk('((A,B),C);', 'tree1.nwk')
        path2 = tmp_nwk('((A,B),C);', 'tree2.nwk')
        args = make_args(
            infile=path1,
            infile2=path2,
            outfile='-',
            format2='auto',
            metric=None,
            dist=None,
            comparison='rooted',
        )

        with pytest.raises(ValueError, match='branch length on every non-root edge'):
            dist_main(args)

    @pytest.mark.parametrize('length', ['-1', 'nan', 'inf'])
    def test_branch_metrics_reject_invalid_lengths(self, tmp_nwk, length):
        path1 = tmp_nwk(
            '((A:{0},B:1):1,(C:1,D:1):1);'.format(length),
            'tree1.nwk',
        )
        path2 = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);', 'tree2.nwk')
        args = make_args(
            infile=path1,
            infile2=path2,
            outfile='-',
            format2='auto',
            metric=['path-length'],
            dist=None,
            comparison='rooted',
        )

        with pytest.raises(ValueError, match='finite, nonnegative branch lengths'):
            dist_main(args)

    def test_unsupported_metric_raises_clear_error(self, tmp_nwk):
        nwk = '((A:1,B:1):1,(C:1,D:1):1);'
        path1 = tmp_nwk(nwk, 'tree1.nwk')
        path2 = tmp_nwk(nwk, 'tree2.nwk')
        args = make_args(
            infile=path1,
            infile2=path2,
            outfile='-',
            format2='auto',
            metric=['quartet'],
            dist=None,
            comparison='rooted',
        )

        with pytest.raises(ValueError, match="Unsupported metric 'quartet'"):
            dist_main(args)

    def test_metric_and_deprecated_dist_cannot_be_combined(self):
        args = make_args(metric=['rf'], dist='RF')
        with pytest.raises(ValueError, match='either'):
            dist_main(args)
