import pandas as pd
import pytest

from nwkit.compose import compose_main
from nwkit.util import read_tree
from tests.helpers import make_args


def test_compose_combines_root_names_support_lengths_and_property(tmp_nwk, tmp_path):
    topology = tmp_nwk('(A:1,B:1,(C:1,D:1):1);', 'topology.nwk')
    root_source = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);', 'root.nwk')
    name_source = tmp_nwk('((A:1,B:1)AB:1,(C:1,D:1)CD:1)ROOT;', 'names.nwk')
    support_source = tmp_nwk('((A:1,B:1)95:1,(C:1,D:1)80:1);', 'support.nwk')
    length_source = tmp_nwk('((A:2,B:3):4,(C:5,D:6):7);', 'length.nwk')
    property_source = tmp_nwk(
        '((A:1,B:1):1[&&NHX:color=red],(C:1,D:1):1[&&NHX:color=blue]);',
        'properties.nwk',
    )
    outfile = tmp_path / 'composed.nwk'
    report = tmp_path / 'composition.tsv'
    args = make_args(
        infile=topology,
        outfile=str(outfile),
        format='1',
        outformat='auto',
        manifest=None,
        root_source=root_source,
        name_source=name_source,
        support_source=support_source,
        length_source=length_source,
        property_source=['color@{}'.format(property_source)],
        source_format='auto',
        taxon_mode='exact',
        policy='compatible-only',
        match_basis='clade',
        allow_projected_values=False,
        report=str(report),
    )
    compose_main(args)
    output = read_tree(str(outfile), format='1', quoted_node_names=True, quiet=True)
    ab = output.common_ancestor(['A', 'B'])
    cd = output.common_ancestor(['C', 'D'])
    assert output.name == 'ROOT'
    assert ab.name == 'AB'
    assert ab.support == pytest.approx(95.0)
    assert ab.dist == pytest.approx(4.0)
    assert ab.props['color'] == 'red'
    assert cd.props['color'] == 'blue'
    rows = pd.read_csv(report, sep='\t')
    assert set(rows['source_property']) >= {'root', 'name', 'support', 'length', 'color'}
    assert (rows['status'] == 'transferred').any()


def test_compose_json_manifest_resolves_relative_paths(tmp_path):
    topology = tmp_path / 'topology.nwk'
    names = tmp_path / 'names.nwk'
    manifest = tmp_path / 'compose.json'
    outfile = tmp_path / 'output.nwk'
    topology.write_text('((A:1,B:1):1,C:1);')
    names.write_text('((A:1,B:1)AB:1,C:1)ROOT;')
    manifest.write_text('{"name": "names.nwk"}')
    args = make_args(
        infile=str(topology),
        outfile=str(outfile),
        format='1',
        outformat='auto',
        manifest=str(manifest),
        root_source=None,
        name_source=None,
        support_source=None,
        length_source=None,
        property_source=[],
        source_format='auto',
        taxon_mode='exact',
        policy='compatible-only',
        match_basis='clade',
        allow_projected_values=False,
        report=None,
    )
    compose_main(args)
    output = read_tree(str(outfile), format='1', quoted_node_names=True, quiet=True)
    assert output.common_ancestor(['A', 'B']).name == 'AB'


def test_compose_blocks_projected_lengths_without_opt_in(tmp_nwk, tmp_path):
    topology = tmp_nwk(
        '((A:1,(B:1,X:1):1):5,(C:1,D:1):1);',
        'topology.nwk',
    )
    lengths = tmp_nwk(
        '((A:1,(B:1,Y:1):1):7,(C:1,D:1):1);',
        'lengths.nwk',
    )
    report = tmp_path / 'composition.tsv'
    args = make_args(
        infile=topology,
        outfile=str(tmp_path / 'blocked.nwk'),
        format='1',
        outformat='1',
        manifest=None,
        root_source=None,
        name_source=None,
        support_source=None,
        length_source=lengths,
        property_source=[],
        source_format='1',
        taxon_mode='intersection',
        policy='compatible-only',
        match_basis='clade',
        allow_projected_values=False,
        report=str(report),
    )
    compose_main(args)
    blocked = read_tree(args.outfile, format='1', quoted_node_names=True, quiet=True)
    assert blocked.common_ancestor(['A', 'B', 'X']).dist == pytest.approx(5.0)
    rows = pd.read_csv(report, sep='\t')
    rejected = rows[rows['target_taxa'] == 'A,B,X'].iloc[0]
    assert rejected['status'] == 'projected_value_rejected'
    assert rejected['match_status'] == 'projected_match'

    args.outfile = str(tmp_path / 'allowed.nwk')
    args.report = str(tmp_path / 'allowed.tsv')
    args.allow_projected_values = True
    compose_main(args)
    allowed = read_tree(args.outfile, format='1', quoted_node_names=True, quiet=True)
    assert allowed.common_ancestor(['A', 'B', 'X']).dist == pytest.approx(7.0)
