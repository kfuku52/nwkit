import pandas as pd

from nwkit.diff import compare_trees, diff_main
from nwkit.util import read_tree
from tests.helpers import make_args


def test_rooted_diff_reports_common_and_tree_specific_clades():
    target = read_tree('(((A:1,B:1):1,C:1):1,D:1);', '1', True, quiet=True)
    source = read_tree('(((A:1,C:1):1,B:1):1,D:1);', '1', True, quiet=True)
    rows = compare_trees(target, source, comparison='rooted', target_class='intnode')
    statuses = {row['status'] for row in rows if row['record_type'] == 'node'}
    assert 'exact_match' in statuses
    assert 'unmatched' in statuses
    assert 'source_only' in statuses


def test_unrooted_diff_with_partial_taxon_overlap(tmp_nwk, tmp_path):
    target_path = tmp_nwk('((A:1,B:1):1,(C:1,(D:1,E:1):1):1);', 'target.nwk')
    source_path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);', 'source.nwk')
    outfile = tmp_path / 'diff.tsv'
    args = make_args(
        infile=target_path,
        infile2=source_path,
        outfile=str(outfile),
        format='1',
        format2='1',
        taxon_mode='intersection',
        comparison='unrooted',
        target='intnode',
        property=[],
        fail_on_difference=False,
    )
    diff_main(args)
    rows = pd.read_csv(outfile, sep='\t')
    leaf_summary = rows[rows['comparison'] == 'leaf_set'].iloc[0]
    assert leaf_summary['status'] == 'different'
    assert 'E' in leaf_summary['target_taxa']
    assert set(rows['comparison']) >= {'leaf_set', 'root_split', 'unrooted_split'}
