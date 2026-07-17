import pytest
from ete4 import Tree

from nwkit.clade_mapping import build_clade_mapping, projected_root_split
from nwkit.root import transfer_root_with_taxon_mode


def _find_match(mapping, taxa):
    taxa = frozenset(taxa)
    return next(match for match in mapping.matches if match.target_taxa == taxa)


def test_exact_mapping_rejects_different_tip_sets():
    target = Tree('((A:1,B:1):1,C:1);', parser=1)
    source = Tree('((A:1,B:1):1,D:1);', parser=1)
    with pytest.raises(ValueError, match='target-only=C'):
        build_clade_mapping(target, source, taxon_mode='exact')


def test_intersection_maps_unique_projected_clade():
    target = Tree('((A:1,B:1):1,(C:1,(D:1,E:1):1):1);', parser=1)
    source = Tree('((A:1,B:1)AB:1,(C:1,D:1)CD:1)R;', parser=1)
    mapping = build_clade_mapping(target, source, taxon_mode='intersection')
    match = _find_match(mapping, {'C', 'D', 'E'})
    assert match.status == 'projected_match'
    assert match.projected_taxa == frozenset({'C', 'D'})
    assert match.source.name == 'CD'
    assert mapping.target_only_taxa == frozenset({'E'})


def test_unique_projection_is_not_reported_as_an_exact_branch_match():
    target = Tree('((A:1,(B:1,X:1):1):1,(C:1,D:1):1);', parser=1)
    source = Tree('((A:1,(B:1,Y:1):1):1,(C:1,D:1):1);', parser=1)
    mapping = build_clade_mapping(target, source, taxon_mode='intersection')
    match = _find_match(mapping, {'A', 'B', 'X'})
    assert match.status == 'projected_match'
    assert match.projected_taxa == frozenset({'A', 'B'})
    assert match.source_taxa == frozenset({'A', 'B', 'Y'})


def test_split_basis_matches_same_edge_across_different_rootings():
    target = Tree('((((A:1,B:1)AB:1,C:1)ABC:1,D:1)ABCD:1,E:1)R;', parser=1)
    source = Tree('(A:1,(B:1,(C:1,(D:1,E:1)DE:1)CDE:1)BCDE:1)R;', parser=1)
    clade_mapping = build_clade_mapping(target, source, match_basis='clade')
    split_mapping = build_clade_mapping(target, source, match_basis='split')
    assert _find_match(clade_mapping, {'A', 'B'}).status == 'unmatched'
    split_match = _find_match(split_mapping, {'A', 'B'})
    assert split_match.status == 'exact_match'
    assert split_match.source_taxa == frozenset({'C', 'D', 'E'})
    assert split_match.reason == 'matching_canonical_split'


def test_intersection_rejects_nonunique_projection():
    target = Tree('(((A:1,B:1):1,X:1):1,(C:1,D:1):1);', parser=1)
    source = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
    mapping = build_clade_mapping(target, source, taxon_mode='intersection')
    match = _find_match(mapping, {'A', 'B'})
    assert match.status == 'ambiguous'
    assert match.reason == 'target_projection_not_unique'


def test_intersection_rejects_internal_clade_defined_by_one_shared_tip():
    target = Tree('((A:1,B:1):1,(C:1,(D:1,X:1):1):1);', parser=1)
    source = Tree('((A:1,B:1):1,(C:1,(D:1,Y:1):1):1);', parser=1)
    mapping = build_clade_mapping(target, source, taxon_mode='intersection')
    match = _find_match(mapping, {'D', 'X'})
    assert match.status == 'ambiguous'
    assert match.reason == 'fewer_than_two_shared_descendant_taxa'


def test_partial_root_transfer_matches_projected_split_and_preserves_tips():
    target = Tree('(A:1,B:2,(C:3,D:4):1,X:2);', parser=1)
    source = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
    result = transfer_root_with_taxon_mode(
        target,
        source,
        taxon_mode='intersection',
    )
    shared = frozenset({'A', 'B', 'C', 'D'})
    assert projected_root_split(result, shared) == projected_root_split(source, shared)
    assert set(result.leaf_names()) == {'A', 'B', 'C', 'D', 'X'}
