from dataclasses import dataclass

from nwkit.util import get_node_class, validate_unique_named_leaves


SUPPORTED_TAXON_MODES = ('exact', 'intersection')
SUPPORTED_MATCH_BASES = ('clade', 'split')


@dataclass(frozen=True)
class CladeMatch:
    target: object
    source: object | None
    target_candidates: tuple
    source_candidates: tuple
    status: str
    reason: str
    match_basis: str
    projected_taxa: frozenset
    projected_split: tuple | None
    target_taxa: frozenset
    source_taxa: frozenset | None


@dataclass(frozen=True)
class CladeMapping:
    matches: tuple
    shared_taxa: frozenset
    target_only_taxa: frozenset
    source_only_taxa: frozenset
    unmatched_source_nodes: tuple
    match_basis: str


def _leaf_set(tree):
    return frozenset(str(name) for name in tree.leaf_names())


def _node_taxon_sets(tree):
    cached = tree.get_cached_content(prop='name')
    return {
        node: frozenset(str(name) for name in cached[node])
        for node in tree.traverse()
    }


def _format_taxa(taxa):
    return ','.join(sorted(str(taxon) for taxon in taxa))


def _projected_split_from_taxa(taxa, shared_taxa):
    side = frozenset(taxa & shared_taxa)
    complement = frozenset(shared_taxa - side)
    if (not side) or (not complement):
        return None
    return canonical_split(side, complement)


def _mapping_key(node_class, taxa, shared_taxa, match_basis):
    projected_taxa = frozenset(taxa & shared_taxa)
    if node_class in ('root', 'leaf') or match_basis == 'clade':
        return projected_taxa
    return _projected_split_from_taxa(taxa, shared_taxa)


def build_clade_mapping(target, source, taxon_mode='exact', match_basis='clade'):
    if taxon_mode not in SUPPORTED_TAXON_MODES:
        raise ValueError(
            "Unsupported taxon mode '{}'. Choose from: {}.".format(
                taxon_mode,
                ', '.join(SUPPORTED_TAXON_MODES),
            )
        )
    if match_basis not in SUPPORTED_MATCH_BASES:
        raise ValueError(
            "Unsupported match basis '{}'. Choose from: {}.".format(
                match_basis,
                ', '.join(SUPPORTED_MATCH_BASES),
            )
        )
    validate_unique_named_leaves(target, option_name='--infile', context=' for clade mapping')
    validate_unique_named_leaves(source, option_name='--infile2', context=' for clade mapping')
    target_leaf_set = _leaf_set(target)
    source_leaf_set = _leaf_set(source)
    if taxon_mode == 'exact' and target_leaf_set != source_leaf_set:
        target_only = target_leaf_set - source_leaf_set
        source_only = source_leaf_set - target_leaf_set
        details = list()
        if target_only:
            details.append('target-only={}'.format(_format_taxa(target_only)))
        if source_only:
            details.append('source-only={}'.format(_format_taxa(source_only)))
        raise ValueError(
            'Leaf labels must match exactly when --taxon-mode exact ({})'.format(
                '; '.join(details)
            )
        )
    shared_taxa = target_leaf_set & source_leaf_set
    if not shared_taxa:
        raise ValueError('The input trees do not share any uniquely named tips.')

    target_taxa_by_node = _node_taxon_sets(target)
    source_taxa_by_node = _node_taxon_sets(source)
    same_leaf_set = target_leaf_set == source_leaf_set
    target_groups = dict()
    source_groups = dict()
    for node, taxa in target_taxa_by_node.items():
        node_class = get_node_class(node)
        key = (node_class, _mapping_key(node_class, taxa, shared_taxa, match_basis))
        target_groups.setdefault(key, list()).append(node)
    for node, taxa in source_taxa_by_node.items():
        node_class = get_node_class(node)
        key = (node_class, _mapping_key(node_class, taxa, shared_taxa, match_basis))
        source_groups.setdefault(key, list()).append(node)

    matches = list()
    matched_source_ids = set()
    for target_node in target.traverse():
        node_class = get_node_class(target_node)
        target_taxa = target_taxa_by_node[target_node]
        projected_taxa = frozenset(target_taxa & shared_taxa)
        projected_split = _projected_split_from_taxa(target_taxa, shared_taxa)
        match_key = _mapping_key(node_class, target_taxa, shared_taxa, match_basis)
        key = (node_class, match_key)
        target_candidates = target_groups.get(key, [])
        source_candidates = source_groups.get(key, [])
        source_node = None
        if node_class == 'root':
            source_node = source
            status = 'exact_match' if same_leaf_set else 'projected_match'
            reason = 'root_to_root'
        elif not projected_taxa:
            status = 'unmatched'
            reason = 'no_shared_descendant_taxa'
        elif taxon_mode == 'intersection' and node_class == 'intnode' and len(projected_taxa) < 2:
            status = 'ambiguous'
            reason = 'fewer_than_two_shared_descendant_taxa'
        elif (
            node_class == 'intnode'
            and (
                match_key is None
                or (match_basis == 'clade' and projected_taxa == shared_taxa)
            )
        ):
            status = 'ambiguous'
            reason = (
                'projection_contains_all_shared_taxa'
                if match_basis == 'clade'
                else 'projection_does_not_define_usable_split'
            )
        elif len(target_candidates) > 1:
            status = 'ambiguous'
            reason = 'target_projection_not_unique'
        elif len(source_candidates) == 0:
            status = 'unmatched'
            reason = 'clade_absent_from_source'
        elif len(source_candidates) > 1:
            status = 'ambiguous'
            reason = 'source_projection_not_unique'
        else:
            source_node = source_candidates[0]
            status = (
                'exact_match'
                if same_leaf_set
                else 'projected_match'
            )
            reason = (
                'matching_descendant_taxa'
                if match_basis == 'clade'
                else 'matching_canonical_split'
            )
        source_taxa = None
        if source_node is not None:
            source_taxa = source_taxa_by_node[source_node]
            matched_source_ids.add(id(source_node))
        matches.append(
            CladeMatch(
                target=target_node,
                source=source_node,
                target_candidates=tuple(target_candidates),
                source_candidates=tuple(source_candidates),
                status=status,
                reason=reason,
                match_basis=match_basis,
                projected_taxa=projected_taxa,
                projected_split=projected_split,
                target_taxa=target_taxa,
                source_taxa=source_taxa,
            )
        )
    unmatched_source_nodes = tuple(
        node for node in source.traverse()
        if id(node) not in matched_source_ids
    )
    return CladeMapping(
        matches=tuple(matches),
        shared_taxa=shared_taxa,
        target_only_taxa=target_leaf_set - source_leaf_set,
        source_only_taxa=source_leaf_set - target_leaf_set,
        unmatched_source_nodes=unmatched_source_nodes,
        match_basis=match_basis,
    )


def projected_root_split(tree, shared_taxa):
    children = tree.get_children()
    if len(children) != 2:
        return None
    taxon_sets = _node_taxon_sets(tree)
    sides = [frozenset(taxon_sets[child] & shared_taxa) for child in children]
    if (not sides[0]) or (not sides[1]) or ((sides[0] | sides[1]) != shared_taxa):
        return None
    return canonical_split(sides[0], sides[1])


def canonical_split(side_a, side_b):
    key_a = tuple(sorted(side_a))
    key_b = tuple(sorted(side_b))
    if (len(side_a), key_a) <= (len(side_b), key_b):
        return frozenset(side_a), frozenset(side_b)
    return frozenset(side_b), frozenset(side_a)


def node_projected_split(node, shared_taxa, taxon_sets=None):
    if node.is_root:
        return None
    if taxon_sets is None:
        root = node
        while root.up is not None:
            root = root.up
        taxon_sets = _node_taxon_sets(root)
    side = frozenset(taxon_sets[node] & shared_taxa)
    complement = frozenset(shared_taxa - side)
    if (not side) or (not complement):
        return None
    return canonical_split(side, complement)


def find_root_split_candidates(target, source, taxon_mode='exact'):
    target_leaf_set = _leaf_set(target)
    source_leaf_set = _leaf_set(source)
    if taxon_mode == 'exact' and target_leaf_set != source_leaf_set:
        raise ValueError('Leaf labels must match exactly when --taxon-mode exact.')
    shared_taxa = target_leaf_set & source_leaf_set
    if len(shared_taxa) < 2:
        raise ValueError('Root transfer requires at least two shared tips.')
    source_split = projected_root_split(source, shared_taxa)
    if source_split is None:
        raise ValueError('The source root does not define a usable bipartition on shared tips.')
    target_taxa_by_node = _node_taxon_sets(target)
    candidates = list()
    for node in target.traverse():
        split = node_projected_split(node, shared_taxa, taxon_sets=target_taxa_by_node)
        if split == source_split:
            candidates.append(node)
    return shared_taxa, source_split, candidates
