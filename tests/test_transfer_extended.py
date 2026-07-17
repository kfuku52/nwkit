import pandas as pd
import pytest

from nwkit.transfer import transfer_main
from nwkit.util import read_tree
from tests.helpers import make_args


def test_partial_taxon_property_transfer_and_report(tmp_nwk, tmp_path):
    target_path = tmp_nwk(
        '((A:1,B:1):1,(C:1,(D:1,E:1):1):1);',
        'target.nwk',
    )
    source_path = tmp_nwk(
        '((A:1,B:1)AB:1[&&NHX:color=red],(C:1,D:1)CD:1[&&NHX:color=blue])R;',
        'source.nwk',
    )
    outfile = tmp_path / 'output.nwk'
    report = tmp_path / 'transfer.tsv'
    args = make_args(
        infile=target_path,
        infile2=source_path,
        outfile=str(outfile),
        format='1',
        format2='1',
        outformat='1',
        target='intnode',
        name=False,
        support=False,
        length=False,
        property=['color'],
        property_map=[],
        fill=None,
        taxon_mode='intersection',
        policy='compatible-only',
        report=str(report),
    )
    transfer_main(args)
    output = read_tree(str(outfile), format='1', quoted_node_names=True, quiet=True)
    assert output.common_ancestor(['A', 'B']).props['color'] == 'red'
    assert output.common_ancestor(['C', 'D', 'E']).props['color'] == 'blue'
    rows = pd.read_csv(report, sep='\t')
    assert set(rows['status']) >= {'transferred', 'ambiguous'}
    assert rows['shared_taxon_count'].unique().tolist() == [4]


def test_property_map_renames_nhx_property(tmp_nwk, tmp_outfile):
    target_path = tmp_nwk('((A:1,B:1):1,C:1);', 'target.nwk')
    source_path = tmp_nwk('((A:1,B:1):1[&&NHX:score=0.8],C:1);', 'source.nwk')
    args = make_args(
        infile=target_path,
        infile2=source_path,
        outfile=tmp_outfile,
        format='1',
        format2='1',
        outformat='1',
        target='intnode',
        name=False,
        support=False,
        length=False,
        property=[],
        property_map=['score=confidence'],
        fill=None,
        taxon_mode='exact',
        policy='compatible-only',
        report=None,
    )
    transfer_main(args)
    output = read_tree(tmp_outfile, format='1', quoted_node_names=True, quiet=True)
    assert float(output.common_ancestor(['A', 'B']).props['confidence']) == pytest.approx(0.8)


def test_strict_transfer_fails_on_incompatible_clade(tmp_nwk, tmp_path):
    target_path = tmp_nwk('(((A:1,B:1):1,C:1):1,D:1);', 'target.nwk')
    source_path = tmp_nwk('(((A:1,C:1)AC:1,B:1):1,D:1);', 'source.nwk')
    report = tmp_path / 'strict.tsv'
    args = make_args(
        infile=target_path,
        infile2=source_path,
        outfile=str(tmp_path / 'unused.nwk'),
        format='1',
        format2='1',
        outformat='1',
        target='intnode',
        name=True,
        support=False,
        length=False,
        property=[],
        property_map=[],
        fill=None,
        taxon_mode='exact',
        policy='strict',
        report=str(report),
    )
    with pytest.raises(ValueError, match='Strict transfer failed'):
        transfer_main(args)
    assert report.exists()
