import pandas as pd
import pytest

from nwkit.annotate import annotate_main
from nwkit.util import read_tree
from tests.helpers import make_args


def test_annotate_tips_and_aggregate_internal_properties(tmp_nwk, tmp_path):
    tree_path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);', 'tree.nwk')
    table_path = tmp_path / 'traits.tsv'
    table_path.write_text(
        'leaf_name\tstate\tvalue\n'
        'A\tx\t1\n'
        'B\tx\t2\n'
        'C\ty\t3\n'
        'D\tz\t4\n'
    )
    outfile = tmp_path / 'annotated.nwk'
    report = tmp_path / 'annotation.tsv'
    args = make_args(
        infile=tree_path,
        outfile=str(outfile),
        format='1',
        outformat='1',
        table=str(table_path),
        columns='state,value',
        property_map=[],
        aggregate=['state:unique:state_consensus', 'value:mean:value_mean'],
        missing_values=',NA',
        unmatched='error',
        report=str(report),
    )
    annotate_main(args)
    output = read_tree(str(outfile), format='1', quoted_node_names=True, quiet=True)
    assert output['A'].props['state'] == 'x'
    assert float(output['A'].props['value']) == pytest.approx(1.0)
    ab = output.common_ancestor(['A', 'B'])
    cd = output.common_ancestor(['C', 'D'])
    assert ab.props['state_consensus'] == 'x'
    assert float(ab.props['value_mean']) == pytest.approx(1.5)
    assert 'state_consensus' not in cd.props
    rows = pd.read_csv(report, sep='\t')
    assert set(rows['status']) >= {'annotated', 'aggregated', 'not_aggregated'}


def test_annotate_unmatched_error(tmp_nwk, tmp_path):
    tree_path = tmp_nwk('(A:1,B:1);', 'tree.nwk')
    table_path = tmp_path / 'traits.tsv'
    table_path.write_text('leaf_name\tstate\nA\tx\nC\ty\n')
    args = make_args(
        infile=tree_path,
        outfile=str(tmp_path / 'unused.nwk'),
        format='1',
        outformat='1',
        table=str(table_path),
        columns='state',
        property_map=[],
        aggregate=[],
        missing_values=',NA',
        unmatched='error',
        report=None,
    )
    with pytest.raises(ValueError, match='tips differ'):
        annotate_main(args)


def test_annotate_unmatched_report_keeps_node_class_vocabulary(tmp_nwk, tmp_path):
    tree_path = tmp_nwk('(A:1,B:1);', 'tree.nwk')
    table_path = tmp_path / 'traits.tsv'
    table_path.write_text('leaf_name\tstate\nA\tx\nZ\ty\n')
    report_path = tmp_path / 'report.tsv'
    annotate_main(make_args(
        infile=tree_path,
        outfile=str(tmp_path / 'annotated.nwk'),
        format='1',
        outformat='1',
        table=str(table_path),
        columns='state',
        property_map=[],
        aggregate=[],
        missing_values=',NA',
        unmatched='ignore',
        report=str(report_path),
    ))

    report = pd.read_csv(report_path, sep='\t', keep_default_na=False)
    assert set(report['node_class']) <= {'', 'root', 'intnode', 'leaf'}
    table_only = report.loc[report['reason'] == 'table_row_absent_from_tree'].iloc[0]
    assert table_only['node_class'] == ''
