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
        match_basis='clade',
        allow_projected_values=False,
        report=str(report),
    )
    transfer_main(args)
    output = read_tree(str(outfile), format='1', quoted_node_names=True, quiet=True)
    assert output.common_ancestor(['A', 'B']).props['color'] == 'red'
    assert output.common_ancestor(['C', 'D', 'E']).props['color'] == 'blue'
    rows = pd.read_csv(report, sep='\t')
    assert set(rows['status']) >= {'transferred', 'ambiguous'}
    projected = rows[rows['status'] == 'transferred']
    assert 'projected_match' in set(projected['match_status'])
    assert projected['projection_only'].any()
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
        match_basis='clade',
        allow_projected_values=False,
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
        match_basis='clade',
        allow_projected_values=False,
        report=str(report),
    )
    with pytest.raises(ValueError, match='Strict transfer failed'):
        transfer_main(args)
    assert report.exists()


def test_projected_branch_length_requires_explicit_opt_in(tmp_nwk, tmp_path):
    target_path = tmp_nwk(
        '((A:1,(B:1,X:1):1):5,(C:1,D:1):1);',
        'target.nwk',
    )
    source_path = tmp_nwk(
        '((A:1,(B:1,Y:1):1):7,(C:1,D:1):1);',
        'source.nwk',
    )
    report = tmp_path / 'blocked.tsv'
    outfile = tmp_path / 'blocked.nwk'
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
        length=True,
        property=[],
        property_map=[],
        fill=None,
        taxon_mode='intersection',
        policy='compatible-only',
        match_basis='clade',
        allow_projected_values=False,
        report=str(report),
    )
    transfer_main(args)
    output = read_tree(str(outfile), format='1', quoted_node_names=True, quiet=True)
    assert output.common_ancestor(['A', 'B', 'X']).dist == pytest.approx(5.0)
    rows = pd.read_csv(report, sep='\t')
    rejected = rows[rows['target_taxa'] == 'A,B,X'].iloc[0]
    assert rejected['match_status'] == 'projected_match'
    assert rejected['status'] == 'projected_value_rejected'
    assert rejected['source_value'] == pytest.approx(7.0)

    args.outfile = str(tmp_path / 'allowed.nwk')
    args.report = str(tmp_path / 'allowed.tsv')
    args.allow_projected_values = True
    transfer_main(args)
    allowed = read_tree(args.outfile, format='1', quoted_node_names=True, quiet=True)
    assert allowed.common_ancestor(['A', 'B', 'X']).dist == pytest.approx(7.0)


def test_length_transfer_from_different_root_preserves_target_root_ratio(tmp_nwk, tmp_path):
    target_path = tmp_nwk(
        '((A:1,B:1):2,(C:1,(D:1,(E:1,F:1):1):1):8);',
        'target.nwk',
    )
    source_path = tmp_nwk(
        '((E:6,F:7):4,(D:5,(C:4,(A:2,B:3):6):3):4);',
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
        target='all',
        name=False,
        support=False,
        length=True,
        property=[],
        property_map=[],
        fill=None,
        taxon_mode='exact',
        policy='strict',
        match_basis='clade',
        allow_projected_values=False,
        report=str(report),
    )

    transfer_main(args)

    output = read_tree(str(outfile), format='1', quoted_node_names=True, quiet=True)
    root_length_by_side = {
        frozenset(child.leaf_names()): child.dist
        for child in output.get_children()
    }
    assert root_length_by_side[frozenset({'A', 'B'})] == pytest.approx(1.2)
    assert root_length_by_side[frozenset({'C', 'D', 'E', 'F'})] == pytest.approx(4.8)
    rows = pd.read_csv(report, sep='\t')
    root_edge_rows = rows[rows['reason'] == 'matching_aligned_root_edge_total']
    assert len(root_edge_rows) == 2


def test_partial_tip_identity_does_not_imply_exact_terminal_branch(tmp_nwk, tmp_path):
    target_path = tmp_nwk(
        '((A:5,X:1):1,(B:1,C:1):1);',
        'target.nwk',
    )
    source_path = tmp_nwk(
        '((A:7,Y:1):1,(B:1,C:1):1);',
        'source.nwk',
    )
    report = tmp_path / 'terminal-branch.tsv'
    outfile = tmp_path / 'terminal-branch.nwk'
    args = make_args(
        infile=target_path,
        infile2=source_path,
        outfile=str(outfile),
        format='1',
        format2='1',
        outformat='1',
        target='leaf',
        name=False,
        support=False,
        length=True,
        property=[],
        property_map=[],
        fill=None,
        taxon_mode='intersection',
        policy='compatible-only',
        match_basis='clade',
        allow_projected_values=False,
        report=str(report),
    )
    transfer_main(args)
    output = read_tree(str(outfile), format='1', quoted_node_names=True, quiet=True)
    assert output['A'].dist == pytest.approx(5.0)
    rows = pd.read_csv(report, sep='\t')
    a_row = rows[rows['target_taxa'] == 'A'].iloc[0]
    assert a_row['match_status'] == 'projected_match'
    assert a_row['status'] == 'projected_value_rejected'


def test_strict_policy_rejects_projected_annotation_match(tmp_nwk, tmp_path):
    target_path = tmp_nwk(
        '((A:1,(B:1,X:1):1):1,(C:1,D:1):1);',
        'target.nwk',
    )
    source_path = tmp_nwk(
        '((A:1,(B:1,Y:1):1):1[&&NHX:color=red],(C:1,D:1):1);',
        'source.nwk',
    )
    report = tmp_path / 'strict-projected.tsv'
    args = make_args(
        infile=target_path,
        infile2=source_path,
        outfile=str(tmp_path / 'unused.nwk'),
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
        policy='strict',
        match_basis='clade',
        allow_projected_values=False,
        report=str(report),
    )
    with pytest.raises(ValueError, match='Strict transfer failed'):
        transfer_main(args)
    rows = pd.read_csv(report, sep='\t')
    rejected = rows[rows['target_taxa'] == 'A,B,X'].iloc[0]
    assert rejected['match_status'] == 'projected_match'
    assert rejected['status'] == 'projected_match_rejected'
    assert rejected['source_value'] == 'red'


def test_split_basis_transfers_edge_name_across_rootings(tmp_nwk, tmp_path):
    target_path = tmp_nwk(
        '((((A:1,B:1):1,C:1):1,D:1):1,E:1);',
        'target.nwk',
    )
    source_path = tmp_nwk(
        '(A:1,(B:1,(C:1,(D:1,E:1):1)EdgeAB:1):1);',
        'source.nwk',
    )
    outfile = tmp_path / 'split.nwk'
    args = make_args(
        infile=target_path,
        infile2=source_path,
        outfile=str(outfile),
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
        policy='compatible-only',
        match_basis='split',
        allow_projected_values=False,
        report=None,
    )
    transfer_main(args)
    output = read_tree(str(outfile), format='1', quoted_node_names=True, quiet=True)
    assert output.common_ancestor(['A', 'B']).name == 'EdgeAB'
