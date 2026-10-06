"""Targeted complete-tip proposals; D+L is a search heuristic, not evidence of error."""

from dataclasses import dataclass, field
from itertools import combinations, product
from typing import Any

from nwkit.clade_index import CladeIndex, LcaIndex
from nwkit.reconcile import _validate_rooted_binary_tree
from nwkit.rooting_state import set_rooting_info
from nwkit.util import copy_tree_iteratively


def topology_key(tree, *, rooted=True):
    index = CladeIndex(tree)
    full = index.mask_by_node[tree]
    masks = set()
    for node, mask in index.mask_by_node.items():
        if node.is_root or node.is_leaf:
            continue
        if not rooted:
            mask = min(mask, full ^ mask)
            if mask.bit_count() <= 1 or (full ^ mask).bit_count() <= 1:
                continue
        masks.add(mask)
    return index.names, tuple(sorted(masks))


@dataclass(frozen=True)
class DLCost:
    duplications: int
    losses: int
    overlaps: int

    @property
    def total(self):
        return self.duplications + self.losses


class DLContext:
    def __init__(self, species_tree, mapping):
        _validate_rooted_binary_tree(species_tree, "--species-tree")
        self.tree = species_tree
        self.mapping = mapping
        self.leaves = {str(n.name): n for n in species_tree.leaves()}
        missing = sorted(set(mapping.values()) - self.leaves.keys())
        if missing:
            raise ValueError(f"Gene species absent from --species-tree: {missing}")
        self.lca = LcaIndex(species_tree)
        self.species_clades = CladeIndex(species_tree)

    def annotation(self, gene):
        mapped: dict[Any, Any] = {}
        taxa: dict[Any, frozenset[str]] = {}
        duplications = losses = overlaps = 0
        for node in gene.traverse("postorder"):
            if node.is_leaf:
                label = self.mapping[str(node.name)]
                mapped[node] = self.leaves[label]
                taxa[node] = frozenset((label,))
                continue
            left, right = node.children
            parent = self.lca.common_ancestor(mapped[left], mapped[right])
            mapped[node] = parent
            taxa[node] = taxa[left] | taxa[right]
            duplicate = parent in (mapped[left], mapped[right])
            duplications += int(duplicate)
            overlaps += bool(taxa[left] & taxa[right])
            depth = self.lca.depth[self.lca.index_by_node[parent]]
            distance = sum(
                self.lca.depth[self.lca.index_by_node[mapped[c]]] - depth
                for c in (left, right)
            )
            losses += distance if duplicate else distance - 2
        return DLCost(duplications, losses, overlaps), mapped, taxa

    def cost(self, gene):
        return self.annotation(gene)[0]


@dataclass
class Proposal:
    tips: tuple[str, ...]
    reasons: set[str] = field(default_factory=set)
    constraints: set[str] = field(default_factory=set)
    diagnostic_gain: int = 0


def _remove_and_suppress(tree, tips):
    """Detach maximal selected clades, retaining every selected tip for regrafting."""
    tree = copy_tree_iteratively(tree)
    index = CladeIndex(tree)
    selected = set(tips)
    components = []
    for node in list(tree.traverse("preorder")):
        names = set(index.names_for_node(node))
        if names <= selected and (
            node.up is None or not set(index.names_for_node(node.up)) <= selected
        ):
            components.append(node)
    for node in components:
        node.detach()
    for node in list(tree.traverse("postorder")):
        if len(node.children) == 1:
            child = node.children[0]
            child.dist = (child.dist or 0.0) + (node.dist or 0.0)
            child.detach()
            if node.up is None:
                tree = child
            else:
                parent = node.up
                node.detach()
                parent.add_child(child)
    set_rooting_info(tree, True)
    return tree, sorted(components, key=lambda n: tuple(sorted(n.leaf_names())))


def discover_proposals(
    gene, context, *, max_moved_tips=8, max_proposals=64, max_set_states=20000
):
    """Enumerate overlap side-cover sets and LCA cut sets, then coupled unions.

    Multiple tips in a cover are proposed together even if no singleton helps.
    Limits are explicit computational budgets, never a claim of exhaustiveness.
    """
    index = CladeIndex(gene)
    baseline, mapped, taxa = context.annotation(gene)
    proposals: dict[tuple[str, ...], Proposal] = {}
    coverage = {
        "set_states": 0,
        "set_states_truncated": False,
        "oversized_sets": 0,
        "insufficient_backbone_sets": 0,
        "proposal_sets_discovered": 0,
        "proposal_sets_truncated": 0,
    }

    def add(tips, reason, constraints):
        key = tuple(sorted(tips))
        if not key:
            return
        if len(key) > len(index.names) - 2:
            coverage["insufficient_backbone_sets"] += 1
            return
        if len(key) > max_moved_tips:
            coverage["oversized_sets"] += 1
            return
        proposal = proposals.setdefault(key, Proposal(key))
        proposal.reasons.add(reason)
        proposal.constraints.update(constraints)

    for node in sorted(
        gene.traverse("postorder"),
        key=lambda n: (index.count_for_node(n), index.names_for_node(n)),
    ):
        if node.is_leaf:
            continue
        children = sorted(node.children, key=index.names_for_node)
        left, right = children
        constraint = index.clade_id_for_node(node)
        overlap = sorted(taxa[left] & taxa[right])
        if overlap:
            choices = [
                [
                    tuple(n for n in index.names_for_node(c) if context.mapping[n] == s)
                    for c in children
                ]
                for s in overlap
            ]
            for sides in product(*choices):
                if coverage["set_states"] >= max_set_states:
                    coverage["set_states_truncated"] = True
                    break
                coverage["set_states"] += 1
                add(
                    {n for side in sides for n in side},
                    "species_overlap_cover",
                    (constraint,),
                )
        if mapped[node] in (mapped[left], mapped[right]):
            # Removing all tips on one species branch lowers the child's LCA.
            # This also detects non-overlap duplications (NAD), which SO misses.
            for child in children:
                if mapped[child] is not mapped[node]:
                    continue
                for branch in mapped[node].children:
                    branch_taxa = set(context.species_clades.names_for_node(branch))
                    cut = {
                        n
                        for n in index.names_for_node(child)
                        if context.mapping[n] not in branch_taxa
                    }
                    add(cut, "lca_cut", (constraint,))
                add(
                    index.names_for_node(child),
                    "duplication_child_clade",
                    (constraint,),
                )

    def rank(proposal):
        return (
            -proposal.diagnostic_gain / len(proposal.tips),
            len(proposal.tips),
            proposal.tips,
        )

    for proposal in proposals.values():
        backbone, _ = _remove_and_suppress(gene, proposal.tips)
        proposal.diagnostic_gain = baseline.total - context.cost(backbone).total
    # Couple related anomalies without requiring successful single-tip endpoints.
    singles = sorted(proposals.values(), key=rank)[:max_proposals]
    coverage["coupling_seed_sets_truncated"] = max(0, len(proposals) - len(singles))
    for first, second in combinations(singles, 2):
        related = bool(first.constraints & second.constraints) or bool(
            {context.mapping[n] for n in first.tips}
            & {context.mapping[n] for n in second.tips}
        )
        if (
            related
            and not set(first.tips) <= set(second.tips)
            and not set(second.tips) <= set(first.tips)
        ):
            add(
                set(first.tips) | set(second.tips),
                "coupled_constraints",
                first.constraints | second.constraints,
            )
    for proposal in proposals.values():
        if "coupled_constraints" in proposal.reasons:
            backbone, _ = _remove_and_suppress(gene, proposal.tips)
            proposal.diagnostic_gain = baseline.total - context.cost(backbone).total
    ordered = sorted(proposals.values(), key=rank)
    coverage["proposal_sets_discovered"] = len(ordered)
    coverage["proposal_sets_truncated"] = max(0, len(ordered) - max_proposals)
    return ordered[:max_proposals], coverage


def _graft(backbone, component, target_names):
    tree = copy_tree_iteratively(backbone)
    index = CladeIndex(tree)
    target = next(n for n in tree.traverse() if index.names_for_node(n) == target_names)
    moved = copy_tree_iteratively(component)
    junction = type(tree)()
    if target.up is None:
        junction.add_child(target)
        tree = junction
    else:
        parent = target.up
        length = target.dist or 0.0
        target.detach()
        junction.dist = length / 2
        target.dist = length / 2
        parent.add_child(junction)
        junction.add_child(target)
    junction.add_child(moved)
    set_rooting_info(tree, True)
    return tree


@dataclass
class Candidate:
    tree: Any
    cost: DLCost
    moved_tips: tuple[str, ...] = ()
    components: int = 0
    reasons: tuple[str, ...] = ()
    id: str = ""


def candidate_rank(candidate):
    return (
        candidate.cost.total,
        candidate.cost.duplications,
        len(candidate.moved_tips),
        topology_key(candidate.tree),
    )


def generate_candidates(
    gene, context, proposals, *, beam_width=16, max_candidates=128, rooted=True
):
    def key(tree):
        return topology_key(tree, rooted=rooted)

    baseline = Candidate(copy_tree_iteratively(gene), context.cost(gene), id="baseline")
    candidates = {key(gene): baseline}
    coverage: dict[str, Any] = {
        "regraft_states": 0,
        "beam_states_discarded": 0,
        "candidate_topologies_discovered": 1,
        "candidate_topologies_truncated": 0,
    }
    champions = []
    for proposal in proposals:
        backbone, components = _remove_and_suppress(gene, proposal.tips)
        beam = [backbone]
        for component in components:
            states: dict[Any, Any] = {}
            for partial in beam:
                index = CladeIndex(partial)
                for target in partial.traverse():
                    candidate = _graft(partial, component, index.names_for_node(target))
                    identity = key(candidate)
                    states.setdefault(identity, candidate)
                    coverage["regraft_states"] += 1
            ordered = sorted(
                states.values(),
                key=lambda t: (
                    context.cost(t).total,
                    context.cost(t).duplications,
                    topology_key(t),
                ),
            )
            coverage["beam_states_discarded"] += max(0, len(ordered) - beam_width)
            beam = ordered[:beam_width]
        endpoints = []
        for tree in beam:
            candidate = Candidate(
                tree,
                context.cost(tree),
                proposal.tips,
                len(components),
                tuple(sorted(proposal.reasons)),
            )
            identity = key(tree)
            if identity == key(gene):
                continue
            old = candidates.get(identity)
            if old is None or len(candidate.moved_tips) < len(old.moved_tips):
                candidates[identity] = candidate
            endpoints.append(identity)
        if endpoints:
            champions.append(
                min(endpoints, key=lambda k: candidate_rank(candidates[k]))
            )
    # Reserve a best endpoint for each detected set, then fill by D+L rank.
    prioritized = [
        *sorted((candidates[k] for k in champions), key=candidate_rank),
        *sorted(candidates.values(), key=candidate_rank),
    ]
    retained = {key(gene): baseline}
    for candidate in prioritized:
        if len(retained) >= max_candidates:
            break
        retained.setdefault(key(candidate.tree), candidate)
    ordered = [
        baseline,
        *sorted(
            (c for c in retained.values() if c is not baseline), key=candidate_rank
        ),
    ]
    for number, candidate in enumerate(ordered[1:], 1):
        candidate.id = f"candidate_{number:06d}"
    coverage["candidate_topologies_discovered"] = len(candidates)
    coverage["candidate_topologies_truncated"] = len(candidates) - len(ordered)
    coverage["topology_equivalence"] = "rooted" if rooted else "unrooted"
    coverage["component_internal_topology"] = "preserved"
    return ordered, coverage
