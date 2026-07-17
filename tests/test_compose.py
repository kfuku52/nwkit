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


def test_compose_preserves_support_and_lengths_across_different_rootings(
    tmp_nwk,
    tmp_path,
):
    topology = tmp_nwk(
        '((E:1,F:1):0.5,(D:1,(C:1,(A:1,B:1):1):1):0.5);',
        'topology.nwk',
    )
    root_source = tmp_nwk(
        '((A:1,B:1):0.2,(C:1,(D:1,(E:1,F:1):1):1):0.8);',
        'root.nwk',
    )
    support_source = tmp_nwk(
        '((E:1,F:1)40:0.5,'
        '(D:1,(C:1,(A:1,B:1)20:1)30:1)40:0.5);',
        'support.nwk',
    )
    length_source = tmp_nwk(
        '((E:6,F:7):4,'
        '(D:5,(C:4,(A:2,B:3):2):3):4);',
        'lengths.nwk',
    )
    report = tmp_path / 'composition.tsv'
    args = make_args(
        infile=topology,
        outfile=str(tmp_path / 'output.nwk'),
        format='1',
        outformat='0',
        manifest=None,
        root_source=root_source,
        name_source=None,
        support_source=support_source,
        length_source=length_source,
        property_source=[],
        source_format='auto',
        taxon_mode='exact',
        policy='strict',
        match_basis='clade',
        allow_projected_values=False,
        report=str(report),
    )

    compose_main(args)

    output = read_tree(args.outfile, format='0', quoted_node_names=True, quiet=True)
    assert output.common_ancestor(['A', 'B']).support == pytest.approx(20.0)
    assert output.common_ancestor(['D', 'E', 'F']).support == pytest.approx(30.0)
    assert output.common_ancestor(['E', 'F']).support == pytest.approx(40.0)
    assert output.common_ancestor(['D', 'E', 'F']).dist == pytest.approx(3.0)
    assert output.common_ancestor(['E', 'F']).dist == pytest.approx(8.0)
    assert output['A'].dist == pytest.approx(2.0)
    root_length_by_side = {
        frozenset(child.leaf_names()): child.dist
        for child in output.get_children()
    }
    assert root_length_by_side[frozenset({'A', 'B'})] == pytest.approx(0.4)
    assert root_length_by_side[frozenset({'C', 'D', 'E', 'F'})] == pytest.approx(1.6)
    complement = next(
        child for child in output.get_children()
        if set(child.leaf_names()) == {'C', 'D', 'E', 'F'}
    )
    assert complement.support == pytest.approx(20.0)
    rows = pd.read_csv(report, sep='\t')
    abc_split = rows[
        (rows['source_property'] == 'support')
        & (rows['target_taxa'] == 'D,E,F')
    ].iloc[0]
    assert abc_split['source_value'] == pytest.approx(30.0)
    assert abc_split['status'] == 'transferred'


def test_compose_preserves_names_and_nhx_properties_across_rootings(tmp_nwk, tmp_path):
    topology = tmp_nwk(
        '((E:1,F:1):0.5,(D:1,(C:1,(A:1,B:1):1):1):0.5);',
        'topology.nwk',
    )
    root_source = tmp_nwk(
        '((A:1,B:1):0.2,(C:1,(D:1,(E:1,F:1):1):1):0.8);',
        'root.nwk',
    )
    annotation_source = tmp_nwk(
        '((E:1,F:1)EDGE_EF:0.5[&&NHX:edge_id=EF],'
        '(D:1,(C:1,(A:1,B:1)EDGE_AB:1[&&NHX:edge_id=AB])'
        'EDGE_ABC:1[&&NHX:edge_id=ABC])EDGE_EF:0.5[&&NHX:edge_id=EF])ROOT;',
        'annotations.nwk',
    )
    args = make_args(
        infile=topology,
        outfile=str(tmp_path / 'output.nwk'),
        format='1',
        outformat='1',
        manifest=None,
        root_source=root_source,
        name_source=annotation_source,
        support_source=None,
        length_source=None,
        property_source=['edge_id@{}'.format(annotation_source)],
        source_format='1',
        taxon_mode='exact',
        policy='compatible-only',
        match_basis='clade',
        allow_projected_values=False,
        report=None,
    )

    compose_main(args)

    output = read_tree(args.outfile, format='1', quoted_node_names=True, quiet=True)
    abc_edge = output.common_ancestor(['D', 'E', 'F'])
    assert abc_edge.name == 'EDGE_ABC'
    assert abc_edge.props['edge_id'] == 'ABC'
    assert output.common_ancestor(['E', 'F']).name == 'EDGE_EF'
    assert output.common_ancestor(['A', 'B']).props['edge_id'] == 'AB'


def test_compose_all_component_sources_with_distinct_root_positions(tmp_nwk, tmp_path):
    topology = tmp_nwk(
        '((E:1,F:1):1,(D:1,(C:1,(A:1,B:1):1):1):1);',
        'topology.nwk',
    )
    root_source = tmp_nwk(
        '((A:1,B:1):0.5,(C:1,(D:1,(E:1,F:1):1):1):1.5);',
        'root.nwk',
    )
    name_source = tmp_nwk(
        '(C:1,((A:1,B:1)NAME_AB:1,'
        '(D:1,(E:1,F:1)NAME_EF:1)NAME_ABC:1):1)NAME_ROOT;',
        'names.nwk',
    )
    support_source = tmp_nwk(
        '(D:1,(((A:1,B:1)21:1,C:1)31:1,(E:1,F:1)41:1):1);',
        'supports.nwk',
    )
    length_source = tmp_nwk(
        '(F:3,(E:4,(D:5,(C:6,(A:7,B:8):10):20):30):7);',
        'lengths.nwk',
    )
    property_source = tmp_nwk(
        '(A:1,(B:1,(C:1,(D:1,'
        '(E:1,F:1):1[&&NHX:edge_tag=EF])'
        ':1[&&NHX:edge_tag=ABC])'
        ':1[&&NHX:edge_tag=AB]):1[&&NHX:edge_tag=PENDANT_A])'
        '[&&NHX:edge_tag=ROOT_PROP];',
        'properties.nwk',
    )
    args = make_args(
        infile=topology,
        outfile=str(tmp_path / 'output.nwk'),
        format='1',
        outformat='1',
        manifest=None,
        root_source=root_source,
        name_source=name_source,
        support_source=support_source,
        length_source=length_source,
        property_source=['edge_tag@{}'.format(property_source)],
        source_format='auto',
        taxon_mode='exact',
        policy='compatible-only',
        match_basis='clade',
        allow_projected_values=False,
        report=str(tmp_path / 'composition.tsv'),
    )

    compose_main(args)

    output = read_tree(args.outfile, format='1', quoted_node_names=True, quiet=True)
    all_taxa = frozenset(output.leaf_names())
    annotations_by_split = dict()
    for node in output.traverse():
        if node.is_root or node.is_leaf:
            continue
        side = frozenset(node.leaf_names())
        split = tuple(sorted(
            (side, all_taxa - side),
            key=lambda taxa: (len(taxa), tuple(sorted(taxa))),
        ))
        annotations_by_split.setdefault(split, set()).add(
            (node.name, node.support, node.props.get('edge_tag'))
        )
    expected = {
        tuple(sorted(
            (frozenset({'A', 'B'}), frozenset({'C', 'D', 'E', 'F'})),
            key=lambda taxa: (len(taxa), tuple(sorted(taxa))),
        )): ('NAME_AB', 21.0, 'AB'),
        tuple(sorted(
            (frozenset({'A', 'B', 'C'}), frozenset({'D', 'E', 'F'})),
            key=lambda taxa: (len(taxa), tuple(sorted(taxa))),
        )): ('NAME_ABC', 31.0, 'ABC'),
        tuple(sorted(
            (frozenset({'E', 'F'}), frozenset({'A', 'B', 'C', 'D'})),
            key=lambda taxa: (len(taxa), tuple(sorted(taxa))),
        )): ('NAME_EF', 41.0, 'EF'),
    }
    for split, annotation in expected.items():
        assert annotations_by_split[split] == {annotation}
    assert output.name == 'NAME_ROOT'
    root_length_by_side = {
        frozenset(child.leaf_names()): child.dist
        for child in output.get_children()
    }
    assert root_length_by_side[frozenset({'A', 'B'})] == pytest.approx(2.5)
    assert root_length_by_side[frozenset({'C', 'D', 'E', 'F'})] == pytest.approx(7.5)
    assert output.common_ancestor(['D', 'E', 'F']).dist == pytest.approx(20.0)
    assert output.common_ancestor(['E', 'F']).dist == pytest.approx(30.0)


def test_compose_without_root_source_does_not_adopt_annotation_source_root(
    tmp_nwk,
    tmp_path,
):
    topology = tmp_nwk('(A:1,B:1,(C:1,D:1):1);', 'topology.nwk')
    support_source = tmp_nwk(
        '(A:1,(B:1,(C:1,D:1)80:1)70:1);',
        'support.nwk',
    )
    args = make_args(
        infile=topology,
        outfile=str(tmp_path / 'output.nwk'),
        format='1',
        outformat='0',
        manifest=None,
        root_source=None,
        name_source=None,
        support_source=support_source,
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

    output = read_tree(args.outfile, format='0', quoted_node_names=True, quiet=True)
    assert len(output.get_children()) == 3
    assert output.common_ancestor(['C', 'D']).support == pytest.approx(80.0)


def test_compose_split_resolves_equal_support_on_root_edges(tmp_nwk, tmp_path):
    topology = tmp_nwk(
        '((E:1,F:1):0.5,(D:1,(C:1,(A:1,B:1):1):1):0.5);',
        'topology.nwk',
    )
    root_source = tmp_nwk(
        '((A:1,B:1):0.2,(C:1,(D:1,(E:1,F:1):1):1):0.8);',
        'root.nwk',
    )
    support_source = tmp_nwk(
        '((E:1,F:1)40:0.5,'
        '(D:1,(C:1,(A:1,B:1)20:1)30:1)40:0.5);',
        'support.nwk',
    )
    report = tmp_path / 'composition.tsv'
    args = make_args(
        infile=topology,
        outfile=str(tmp_path / 'output.nwk'),
        format='1',
        outformat='0',
        manifest=None,
        root_source=root_source,
        name_source=None,
        support_source=support_source,
        length_source=None,
        property_source=[],
        source_format='auto',
        taxon_mode='exact',
        policy='strict',
        match_basis='split',
        allow_projected_values=False,
        report=str(report),
    )

    compose_main(args)

    output = read_tree(args.outfile, format='0', quoted_node_names=True, quiet=True)
    assert output.common_ancestor(['A', 'B']).support == pytest.approx(20.0)
    assert output.common_ancestor(['D', 'E', 'F']).support == pytest.approx(30.0)
    assert output.common_ancestor(['E', 'F']).support == pytest.approx(40.0)
    rows = pd.read_csv(report, sep='\t')
    root_edge_rows = rows[
        (rows['source_property'] == 'support')
        & (rows['ambiguity_resolution'].isin(['unique_edge_value', 'equal_edge_values']))
    ]
    assert len(root_edge_rows) == 3
    assert set(root_edge_rows['status']) == {'transferred'}
    assert set(root_edge_rows['ambiguity_policy']) == {'auto'}


def test_compose_split_auto_uses_matching_side_for_conflicting_support(
    tmp_nwk,
    tmp_path,
):
    topology = tmp_nwk(
        '((A:1,B:1):1,(C:1,(D:1,(E:1,F:1):1):1):1);',
        'topology.nwk',
    )
    support_source = tmp_nwk(
        '((E:1,F:1)40:0.5,'
        '(D:1,(C:1,(A:1,B:1)20:1)30:1)50:0.5);',
        'support.nwk',
    )
    report = tmp_path / 'composition.tsv'
    args = make_args(
        infile=topology,
        outfile=str(tmp_path / 'output.nwk'),
        format='1',
        outformat='0',
        manifest=None,
        root_source=None,
        name_source=None,
        support_source=support_source,
        length_source=None,
        property_source=[],
        source_format='auto',
        taxon_mode='exact',
        policy='compatible-only',
        match_basis='split',
        allow_projected_values=False,
        report=str(report),
    )

    compose_main(args)

    output = read_tree(args.outfile, format='0', quoted_node_names=True, quiet=True)
    assert output.common_ancestor(['E', 'F']).support == pytest.approx(40.0)
    rows = pd.read_csv(report, sep='\t')
    ef_edge = rows[
        (rows['source_property'] == 'support')
        & (rows['target_taxa'] == 'E,F')
    ].iloc[0]
    assert ef_edge['status'] == 'transferred'
    assert ef_edge['reason'] == 'matching_root_edge_side'
    assert ef_edge['ambiguity_policy'] == 'auto'
    assert ef_edge['ambiguity_resolution'] == 'matching_side'
    assert '40.0' in ef_edge['source_candidate_values']
    assert '50.0' in ef_edge['source_candidate_values']


def test_compose_split_equal_only_rejects_conflicting_root_edge_support(
    tmp_nwk,
    tmp_path,
):
    topology = tmp_nwk(
        '((A:1,B:1):1,(C:1,(D:1,(E:1,F:1):1):1):1);',
        'topology.nwk',
    )
    support_source = tmp_nwk(
        '((E:1,F:1)40:0.5,'
        '(D:1,(C:1,(A:1,B:1)20:1)30:1)50:0.5);',
        'support.nwk',
    )
    report = tmp_path / 'composition.tsv'
    args = make_args(
        infile=topology,
        outfile=str(tmp_path / 'output.nwk'),
        format='1',
        outformat='0',
        manifest=None,
        root_source=None,
        name_source=None,
        support_source=support_source,
        length_source=None,
        property_source=[],
        source_format='auto',
        taxon_mode='exact',
        policy='compatible-only',
        match_basis='split',
        root_edge_policy=['support=equal-only'],
        allow_projected_values=False,
        report=str(report),
    )

    compose_main(args)

    rows = pd.read_csv(report, sep='\t')
    ef_edge = rows[
        (rows['source_property'] == 'support')
        & (rows['target_taxa'] == 'E,F')
    ].iloc[0]
    assert ef_edge['status'] == 'ambiguous'
    assert ef_edge['reason'] == 'conflicting_root_edge_values'
    assert ef_edge['ambiguity_policy'] == 'equal-only'
    assert ef_edge['ambiguity_resolution'] == 'equal_values_required'


def test_compose_split_mean_reduces_conflicting_root_edge_support(tmp_nwk, tmp_path):
    topology = tmp_nwk(
        '((A:1,B:1):1,(C:1,(D:1,(E:1,F:1):1):1):1);',
        'topology.nwk',
    )
    support_source = tmp_nwk(
        '((E:1,F:1)40:0.5,'
        '(D:1,(C:1,(A:1,B:1)20:1)30:1)50:0.5);',
        'support.nwk',
    )
    report = tmp_path / 'composition.tsv'
    args = make_args(
        infile=topology,
        outfile=str(tmp_path / 'output.nwk'),
        format='1',
        outformat='0',
        manifest=None,
        root_source=None,
        name_source=None,
        support_source=support_source,
        length_source=None,
        property_source=[],
        source_format='auto',
        taxon_mode='exact',
        policy='strict',
        match_basis='split',
        root_edge_policy=['support=mean'],
        allow_projected_values=False,
        report=str(report),
    )

    compose_main(args)

    output = read_tree(args.outfile, format='0', quoted_node_names=True, quiet=True)
    assert output.common_ancestor(['E', 'F']).support == pytest.approx(45.0)
    rows = pd.read_csv(report, sep='\t')
    ef_edge = rows[
        (rows['source_property'] == 'support')
        & (rows['target_taxa'] == 'E,F')
    ].iloc[0]
    assert ef_edge['status'] == 'transferred'
    assert ef_edge['ambiguity_policy'] == 'mean'
    assert ef_edge['ambiguity_resolution'] == 'mean'
    assert ef_edge['source_candidate_count'] == 2


def test_compose_split_edge_total_transfers_lengths_and_keeps_root_ratio(
    tmp_nwk,
    tmp_path,
):
    topology = tmp_nwk(
        '((E:1,F:1):0.5,(D:1,(C:1,(A:1,B:1):1):1):0.5);',
        'topology.nwk',
    )
    root_source = tmp_nwk(
        '((A:1,B:1):0.2,(C:1,(D:1,(E:1,F:1):1):1):0.8);',
        'root.nwk',
    )
    length_source = tmp_nwk(
        '((E:6,F:7):4,(D:5,(C:4,(A:2,B:3):2):3):4);',
        'lengths.nwk',
    )
    report = tmp_path / 'composition.tsv'
    args = make_args(
        infile=topology,
        outfile=str(tmp_path / 'output.nwk'),
        format='1',
        outformat='1',
        manifest=None,
        root_source=root_source,
        name_source=None,
        support_source=None,
        length_source=length_source,
        property_source=[],
        source_format='1',
        taxon_mode='exact',
        policy='strict',
        match_basis='split',
        allow_projected_values=False,
        report=str(report),
    )

    compose_main(args)

    output = read_tree(args.outfile, format='1', quoted_node_names=True, quiet=True)
    root_length_by_side = {
        frozenset(child.leaf_names()): child.dist
        for child in output.get_children()
    }
    assert root_length_by_side[frozenset({'A', 'B'})] == pytest.approx(0.4)
    assert root_length_by_side[frozenset({'C', 'D', 'E', 'F'})] == pytest.approx(1.6)
    assert output.common_ancestor(['E', 'F']).dist == pytest.approx(8.0)
    rows = pd.read_csv(report, sep='\t')
    total_rows = rows[rows['ambiguity_policy'] == 'edge-total']
    assert set(total_rows['ambiguity_resolution']) >= {
        'edge_total_target_ratio',
        'edge_total_single_branch',
    }


def test_compose_manifest_configures_root_edge_policy(tmp_path):
    topology = tmp_path / 'topology.nwk'
    supports = tmp_path / 'supports.nwk'
    manifest = tmp_path / 'compose.json'
    outfile = tmp_path / 'output.nwk'
    topology.write_text('((A:1,B:1):1,(C:1,(D:1,(E:1,F:1):1):1):1);')
    supports.write_text(
        '((E:1,F:1)40:0.5,'
        '(D:1,(C:1,(A:1,B:1)20:1)30:1)50:0.5);'
    )
    manifest.write_text(
        '{"support":"supports.nwk","root_edge_policies":{"support":"max"}}'
    )
    args = make_args(
        infile=str(topology),
        outfile=str(outfile),
        format='1',
        outformat='0',
        manifest=str(manifest),
        root_source=None,
        name_source=None,
        support_source=None,
        length_source=None,
        property_source=[],
        source_format='auto',
        taxon_mode='exact',
        policy='strict',
        match_basis='split',
        allow_projected_values=False,
        report=None,
    )

    compose_main(args)

    output = read_tree(str(outfile), format='0', quoted_node_names=True, quiet=True)
    assert output.common_ancestor(['E', 'F']).support == pytest.approx(50.0)


def test_compose_manifest_property_can_override_root_edge_policy(tmp_path):
    topology = tmp_path / 'topology.nwk'
    states = tmp_path / 'states.nwk'
    manifest = tmp_path / 'compose.json'
    outfile = tmp_path / 'output.nwk'
    topology.write_text(
        '((A:1,B:1):1,(C:1,(D:1,(E:1,F:1):1):1):1);'
    )
    states.write_text(
        '((E:1,F:1):0.5[&&NHX:state=left],'
        '(D:1,(C:1,(A:1,B:1):1):1):0.5[&&NHX:state=right]);'
    )
    manifest.write_text(
        '{'
        '"root_edge_policies":{"edge_state":"skip"},'
        '"properties":[{'
        '"path":"states.nwk",'
        '"source":"state",'
        '"target":"edge_state",'
        '"root_edge_policy":"matching-side"'
        '}]'
        '}'
    )
    args = make_args(
        infile=str(topology),
        outfile=str(outfile),
        format='1',
        outformat='1',
        manifest=str(manifest),
        root_source=None,
        name_source=None,
        support_source=None,
        length_source=None,
        property_source=[],
        source_format='1',
        taxon_mode='exact',
        policy='compatible-only',
        match_basis='split',
        allow_projected_values=False,
        report=None,
    )

    compose_main(args)

    output = read_tree(str(outfile), format='1', quoted_node_names=True, quiet=True)
    assert output.common_ancestor(['E', 'F']).props['edge_state'] == 'left'

    args.outfile = str(tmp_path / 'cli-output.nwk')
    args.report = str(tmp_path / 'cli-report.tsv')
    args.root_edge_policy = ['edge_state=skip']
    compose_main(args)

    cli_output = read_tree(
        args.outfile,
        format='1',
        quoted_node_names=True,
        quiet=True,
    )
    assert 'edge_state' not in cli_output.common_ancestor(['E', 'F']).props
    rows = pd.read_csv(args.report, sep='\t')
    ef_row = rows[
        (rows['target_property'] == 'edge_state')
        & (rows['target_taxa'] == 'E,F')
    ].iloc[0]
    assert ef_row['ambiguity_policy'] == 'skip'
    assert ef_row['ambiguity_resolution'] == 'skipped'


def test_compose_failed_root_transfer_leaves_target_root_unchanged(tmp_nwk, tmp_path):
    topology = tmp_nwk('(((A:1,B:1):1,C:1):1,D:1);', 'topology.nwk')
    incompatible_root = tmp_nwk(
        '((A:1,C:1):1,(B:1,D:1):1);',
        'incompatible-root.nwk',
    )
    report = tmp_path / 'composition.tsv'
    args = make_args(
        infile=topology,
        outfile=str(tmp_path / 'output.nwk'),
        format='1',
        outformat='1',
        manifest=None,
        root_source=incompatible_root,
        name_source=None,
        support_source=None,
        length_source=None,
        property_source=[],
        source_format='1',
        taxon_mode='exact',
        policy='compatible-only',
        match_basis='clade',
        allow_projected_values=False,
        report=str(report),
    )

    compose_main(args)

    output = read_tree(args.outfile, format='1', quoted_node_names=True, quiet=True)
    assert {
        frozenset(child.leaf_names())
        for child in output.get_children()
    } == {
        frozenset({'A', 'B', 'C'}),
        frozenset({'D'}),
    }
    rows = pd.read_csv(report, sep='\t')
    assert rows.loc[rows['source_property'] == 'root', 'status'].tolist() == ['unmatched']
