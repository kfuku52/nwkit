import pandas as pd
import pytest
from ete4 import Tree

from nwkit.clade_mapping import canonical_split
from nwkit.transfer import (
    parse_root_edge_policies,
    transfer_main,
    transfer_properties,
)
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
    assert rows['num_shared_taxa'].unique().tolist() == [4]


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


def test_root_edge_matching_side_transfers_string_property(tmp_nwk, tmp_path):
    target_path = tmp_nwk(
        '((A:1,B:1):1,(C:1,(D:1,(E:1,F:1):1):1):1);',
        'target.nwk',
    )
    source_path = tmp_nwk(
        '((E:1,F:1):0.5[&&NHX:state=left],'
        '(D:1,(C:1,(A:1,B:1):1):1):0.5[&&NHX:state=right]);',
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
        property=[],
        property_map=['state=edge_state'],
        fill=None,
        taxon_mode='exact',
        policy='compatible-only',
        match_basis='split',
        root_edge_policy=['edge_state=matching-side'],
        allow_projected_values=False,
        report=str(report),
    )

    transfer_main(args)

    output = read_tree(str(outfile), format='1', quoted_node_names=True, quiet=True)
    assert output.common_ancestor(['E', 'F']).props['edge_state'] == 'left'
    rows = pd.read_csv(report, sep='\t')
    ef_row = rows[
        (rows['target_property'] == 'edge_state')
        & (rows['target_taxa'] == 'E,F')
    ].iloc[0]
    assert ef_row['ambiguity_policy'] == 'matching-side'
    assert ef_row['ambiguity_resolution'] == 'matching_side'


def test_root_edge_matching_side_uses_projected_taxa(tmp_nwk, tmp_path):
    target_path = tmp_nwk(
        '((A:1,B:1,X:1):1,(C:1,(D:1,(E:1,F:1):1):1):1);',
        'target.nwk',
    )
    source_path = tmp_nwk(
        '((E:1,F:1,Y:1):0.5[&&NHX:state=left],'
        '(D:1,(C:1,(A:1,B:1):1):1):0.5[&&NHX:state=right]);',
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
        property=[],
        property_map=['state=edge_state'],
        fill=None,
        taxon_mode='intersection',
        policy='compatible-only',
        match_basis='split',
        root_edge_policy=['edge_state=matching-side'],
        allow_projected_values=False,
        report=str(report),
    )

    transfer_main(args)

    output = read_tree(str(outfile), format='1', quoted_node_names=True, quiet=True)
    assert output.common_ancestor(['E', 'F']).props['edge_state'] == 'left'
    rows = pd.read_csv(report, sep='\t')
    ef_row = rows[
        (rows['target_property'] == 'edge_state')
        & (rows['target_taxa'] == 'E,F')
    ].iloc[0]
    assert ef_row['status'] == 'transferred'
    assert ef_row['match_status'] == 'projected_match'
    assert ef_row['shared_descendant_taxa'] == 'E,F'
    assert ef_row['ambiguity_resolution'] == 'matching_side'
    assert '"taxa":["E","F","Y"]' in ef_row['source_candidate_values']


def test_root_edge_total_splits_evenly_when_target_root_total_is_zero(
    tmp_nwk,
    tmp_path,
):
    target_path = tmp_nwk(
        '((A:1,B:1):0,(C:1,(D:1,(E:1,F:1):1):1):0);',
        'target.nwk',
    )
    source_path = tmp_nwk(
        '((E:1,F:1):2,(D:1,(C:1,(A:1,B:1):6):1):4);',
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
        policy='compatible-only',
        match_basis='split',
        root_edge_policy=['length=edge-total'],
        allow_projected_values=False,
        report=str(report),
    )

    transfer_main(args)

    output = read_tree(str(outfile), format='1', quoted_node_names=True, quiet=True)
    assert [child.dist for child in output.get_children()] == pytest.approx([3.0, 3.0])
    rows = pd.read_csv(report, sep='\t')
    equal_split_rows = rows[
        rows['ambiguity_resolution'] == 'edge_total_equal_zero_ratio'
    ]
    assert len(equal_split_rows) == 2
    assert set(equal_split_rows['status']) == {'transferred'}


def test_root_edge_total_is_invariant_across_every_target_source_root_pair():
    newick = '(((A:1,B:1):1,C:1):1,(D:1,(E:1,F:1):1):1);'
    base = Tree(newick, parser=1)
    all_taxa = frozenset(base.leaf_names())

    def node_split(node):
        side = frozenset(node.leaf_names())
        return canonical_split(side, all_taxa - side)

    def split_key(split):
        return tuple(tuple(sorted(side)) for side in split)

    root_sides_by_split = dict()
    for node in base.traverse():
        if not node.is_root:
            root_sides_by_split.setdefault(node_split(node), frozenset(node.leaf_names()))
    assert len(root_sides_by_split) == 9

    def rooted_on(descendant_taxa):
        tree = Tree(newick, parser=1)
        candidates = [
            node
            for node in tree.traverse()
            if not node.is_root
            and frozenset(node.leaf_names()) == descendant_taxa
        ]
        if not candidates:
            candidates = [
                node
                for node in tree.traverse()
                if not node.is_root
                and frozenset(node.leaf_names()) == all_taxa - descendant_taxa
            ]
        candidate = candidates[0]
        if candidate not in tree.get_children():
            tree.dist = 0.0
            tree.set_outgroup(candidate)
        return tree

    def assign_edge_totals(tree, scale, root_ratio):
        nodes_by_split = dict()
        for node in tree.traverse():
            if not node.is_root:
                nodes_by_split.setdefault(node_split(node), []).append(node)
        for index, split in enumerate(sorted(nodes_by_split, key=split_key), start=1):
            nodes = nodes_by_split[split]
            total = float(index * scale)
            if len(nodes) == 2:
                nodes[0].dist = total * root_ratio
                nodes[1].dist = total * (1.0 - root_ratio)
            else:
                assert len(nodes) == 1
                nodes[0].dist = total
        tree.dist = 0.0

    def edge_totals(tree):
        totals = dict()
        for node in tree.traverse():
            if not node.is_root:
                split = node_split(node)
                totals[split] = totals.get(split, 0.0) + float(node.dist or 0.0)
        return totals

    rootings = [rooted_on(side) for side in root_sides_by_split.values()]
    checked_pairs = 0
    for target_rooting in rootings:
        for source_rooting in rootings:
            target = target_rooting.copy(method='deepcopy')
            source = source_rooting.copy(method='deepcopy')
            assign_edge_totals(target, scale=100, root_ratio=0.2)
            assign_edge_totals(source, scale=10, root_ratio=0.3)
            expected = edge_totals(source)

            result = transfer_properties(
                target=target,
                source=source,
                property_specs=[('length', 'length')],
                target_class='all',
                taxon_mode='exact',
                policy='compatible-only',
                align_roots=False,
                exclude_root=True,
                match_basis='split',
                root_edge_policies={'length': 'edge-total'},
            )

            assert edge_totals(result['tree']) == pytest.approx(expected)
            checked_pairs += 1
    assert checked_pairs == 81


def test_aligned_root_length_skip_policy_keeps_target_root_lengths(tmp_nwk, tmp_path):
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
        policy='compatible-only',
        match_basis='clade',
        root_edge_policy=['length=skip'],
        allow_projected_values=False,
        report=str(report),
    )

    transfer_main(args)

    output = read_tree(str(outfile), format='1', quoted_node_names=True, quiet=True)
    root_length_by_side = {
        frozenset(child.leaf_names()): child.dist
        for child in output.get_children()
    }
    assert root_length_by_side[frozenset({'A', 'B'})] == pytest.approx(2.0)
    assert root_length_by_side[frozenset({'C', 'D', 'E', 'F'})] == pytest.approx(8.0)
    rows = pd.read_csv(report, sep='\t')
    skipped = rows[rows['ambiguity_resolution'] == 'skipped_aligned_root_edge']
    assert len(skipped) == 2
    assert set(skipped['status']) == {'ambiguous'}


def test_parse_root_edge_policies_validates_property_policy_combinations():
    assert parse_root_edge_policies([
        'support=mean',
        'length=edge-total',
        'color=matching-side',
    ]) == {
        'support': 'mean',
        'length': 'edge-total',
        'color': 'matching-side',
    }
    with pytest.raises(ValueError, match='not valid for node names'):
        parse_root_edge_policies(['name=mean'])
    with pytest.raises(ValueError, match='not valid for branch lengths'):
        parse_root_edge_policies(['length=max'])
    with pytest.raises(ValueError, match='Unknown root-edge policy'):
        parse_root_edge_policies(['support=first'])


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
