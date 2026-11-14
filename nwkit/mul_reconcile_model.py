"""Exact LCA duplication/loss reconciliation over ambiguous MUL-tree tips."""

import re
from dataclasses import dataclass
from functools import lru_cache
from itertools import product
from typing import Any

from ete4 import Tree

from nwkit.clade_index import LcaIndex
from nwkit.rooting_state import require_rooted
from nwkit.util import (
    _IMMUTABLE_TREE_PROP_TYPES,
    copy_tree_iteratively,
    get_subtree_leaf_index_extents,
    validate_unique_named_leaves,
)


def validate_binary(tree, label, *, unique=True):
    require_rooted(tree, f"{label} must be rooted.")
    if unique:
        validate_unique_named_leaves(tree, label)
    if len(list(tree.leaves())) < 2 or any(
        len(node.children) != 2 for node in tree.traverse() if not node.is_leaf
    ):
        raise ValueError(f"{label} must be strictly binary with at least two tips.")
    if any(not leaf.name for leaf in tree.leaves()):
        raise ValueError(f"{label} has an unnamed tip.")


def label_nodes(tree):
    """GRAMPA-compatible closing-parenthesis/postorder internal numbering."""
    count = 0
    for node in tree.traverse("postorder"):
        if not node.is_leaf:
            count += 1
            node.name = f"<{count}>"


def _clade_selector_lookup(tree, nodes):
    leaf_names = list(tree.leaf_names())
    tip_order = {name: index for index, name in enumerate(leaf_names)}
    if len(tip_order) != len(leaf_names):
        clade_sets = {frozenset(node.leaf_names()): node for node in nodes}
        return lambda members: clade_sets.get(frozenset(members))

    # With unique labels, a clade is an entire interval in tree leaf order.
    extents = get_subtree_leaf_index_extents(tree, tip_order)
    clade_ranges = {extents[node]: node for node in nodes}

    def lookup(members):
        try:
            positions = [tip_order[name] for name in members]
        except KeyError:
            return None
        return clade_ranges.get((min(positions), max(positions), len(positions)))

    return lookup


def select_nodes(tree, selectors):
    nodes = tuple(tree.traverse("postorder"))
    if selectors is None:
        return nodes
    lookup = {node.name: node for node in nodes}
    lookup_clade = None
    result = []
    seen = set()
    for selector in selectors.split():
        key = f"<{selector}>" if selector.isdigit() else selector
        if key in lookup:
            node = lookup[key]
        else:
            members = selector.split(",")
            if len(members) != len(set(members)):
                raise ValueError(
                    f"Duplicate species in hypothesis selector: {selector}"
                )
            if lookup_clade is None:
                lookup_clade = _clade_selector_lookup(tree, nodes)
            node = lookup_clade(members)
            if node is None:
                raise ValueError(f"Unknown or non-monophyletic hypothesis: {selector}")
        if node not in seen:
            seen.add(node)
            result.append(node)
    if not result:
        raise ValueError("Hypothesis selectors must not be empty.")
    return tuple(result)


@dataclass(frozen=True)
class Hypothesis:
    id: int
    tree: Any
    h1: str = "NA"
    h2: str = "NA"
    kind: str = "no-polyploidy"


def hypotheses(tree, h1=None, h2=None, *, multree=False, max_candidates=10000):
    validate_binary(tree, "Species tree", unique=not multree)
    # Keep ETE's copy protocol for subclasses, cached sizes and custom properties.
    plain_tree = all(
        type(node) is Tree
        and node.size == (0.0, 0.0)
        and type(node.props) is dict
        and all(
            type(key) is str and type(value) in _IMMUTABLE_TREE_PROP_TYPES
            for key, value in node.props.items()
        )
        for node in tree.traverse()
    )

    def clone(node):
        if plain_tree:
            return copy_tree_iteratively(node, inherit_rooting=False)
        return node.copy()

    original = clone(tree)
    for leaf in original.leaves():
        leaf.add_prop("mul_species", leaf.name)
    label_nodes(original)
    if multree:
        if h1 is not None or h2 is not None:
            raise ValueError("Supplied MUL-trees cannot also request H1/H2 search.")
        return (Hypothesis(0, original, kind="supplied-MUL-tree"),)
    if any(
        leaf.name.endswith(("+", "*"))
        or re.fullmatch(r"<\d+>", leaf.name)
        or leaf.name.isdigit()
        for leaf in original.leaves()
    ):
        raise ValueError(
            "Species labels cannot be numeric, internal <n> tags, or end in '+'/'*'."
        )
    selected_h1, selected_h2 = select_nodes(original, h1), select_nodes(original, h2)
    leaf_positions = {leaf.name: index for index, leaf in enumerate(original.leaves())}
    bounds = get_subtree_leaf_index_extents(original, leaf_positions)
    pairs: list[tuple[str, str]] = []
    for first in selected_h1:
        first_low, first_high, _ = bounds[first]
        for second in selected_h2:
            second_low, second_high, _ = bounds[second]
            # Clade spans nest exactly along an ancestor chain.
            if (
                first is not second
                and first_low <= second_low
                and second_high <= first_high
            ):
                continue
            if len(pairs) + 1 >= max_candidates:
                raise ValueError(
                    "MUL-tree search exceeds --max-candidates; restrict H1/H2."
                )
            pairs.append((first.name, second.name))
    if not pairs:
        raise ValueError("No admissible H1/H2 combinations: H2 cannot be below H1.")
    result = [Hypothesis(0, original)]
    for number, (first, second) in enumerate(pairs, 1):
        candidate = clone(original)
        lookup = {node.name: node for node in candidate.traverse()}
        copied = clone(lookup[first])
        copied.up = None
        for leaf in copied.leaves():
            leaf.name += "*"
        for leaf in lookup[first].leaves():
            leaf.name += "+"
        destination = lookup[second]
        parent = destination.up
        position = parent.children.index(destination) if parent is not None else 0
        if parent is not None:
            destination.detach()
        joint = Tree()
        joint.add_child(destination)
        joint.add_child(copied)
        if parent is None:
            candidate = joint
        else:
            parent.add_child(joint)
            parent.children.remove(joint)
            parent.children.insert(position, joint)
        label_nodes(candidate)
        result.append(
            Hypothesis(
                number,
                candidate,
                first,
                second,
                "autopolyploid" if first == second else "allopolyploid",
            )
        )
    return tuple(result)


@dataclass
class Cell:
    score: int
    count: int
    choices: list[tuple[int, int]] | None


_CONSUMED_STATE: dict[int, Cell] = {}
_CONSUMED_SCORE: dict[int, int] = {}
_LCA_CACHE_SIZE = 4096
_MAX_CACHED_SUBTREE_TIPS = 8


class _SpeciesIndex:
    def __init__(self, species):
        validate_binary(species, "MUL-tree", unique=False)
        self.lca = LcaIndex(species)
        self.common_ancestor_indices = lru_cache(maxsize=_LCA_CACHE_SIZE)(
            self.lca.common_ancestor_indices
        )
        self.labels = tuple(node.name for node in self.lca.nodes)
        occurrences: dict[str, list[int]] = {}
        for index, node in enumerate(self.lca.nodes):
            if node.is_leaf:
                name = node.props.get("mul_species", node.name)
                occurrences.setdefault(name, []).append(index)
        self.occurrences = occurrences


class _GeneIndex:
    def __init__(self, gene, parser):
        validate_binary(gene, "Gene tree")
        self.gene = gene
        self.nodes = tuple(gene.traverse("postorder"))
        self.indices = {node: index for index, node in enumerate(self.nodes)}
        self.taxa = {
            node: parser.parse(node.name).species_label
            for node in self.nodes
            if node.is_leaf
        }
        self.children = tuple(
            tuple(self.indices[child] for child in node.children) for node in self.nodes
        )
        self.score_keys = None


def prepare_score_cache(genes):
    """Intern frequent small species-colored subtrees once for the search."""
    ids: dict[Any, int] = {}
    sizes: list[int] = []
    counts: list[int] = []
    for gene in genes:
        keys: list[int] = []
        for node, children in zip(gene.nodes, gene.children, strict=True):
            if children:
                a, b = (keys[i] for i in children)
                if a < 0 or b < 0 or sizes[a] + sizes[b] > _MAX_CACHED_SUBTREE_TIPS:
                    keys.append(-1)
                    continue
                key, size = (min(a, b), max(a, b)), sizes[a] + sizes[b]
            else:
                key, size = gene.taxa[node], 1
            index = ids.setdefault(key, len(ids))
            if index == len(counts):
                sizes.append(size)
                counts.append(0)
            counts[index] += 1
            keys.append(index)
        gene.score_keys = keys
    for gene in genes:
        gene.score_keys = tuple(
            i if i >= 0 and counts[i] > 1 else -1 for i in gene.score_keys
        )


def reconcile_genes(genes, species, parser, max_state_pairs=10000000):
    """Reuse the fixed species index while evaluating gene trees in input order."""
    yield from _reconcile_indexed_genes(
        (_GeneIndex(gene, parser) for gene in genes), species, max_state_pairs
    )


def _reconcile_indexed_genes(genes, species, max_state_pairs, *, keep_choices=True):
    species_index = None
    cache = None if keep_choices else {}
    for gene in genes:
        if species_index is None:
            species_index = _SpeciesIndex(species)
        yield Reconciliation._from_indexes(
            gene, species_index, max_state_pairs, keep_choices=keep_choices, cache=cache
        )


def _score_indexed_genes(genes, species, max_state_pairs):
    """Compute exact minimum costs without counting or storing optimal maps."""
    species_index = None
    cached: dict[int, dict[int, int]] = {}
    for gene in genes:
        if species_index is None:
            species_index = _SpeciesIndex(species)
        depth = species_index.lca.depth
        lca = species_index.common_ancestor_indices
        states: list[dict[int, int]] = []
        state_pairs = 0
        for index, (node, children) in enumerate(
            zip(gene.nodes, gene.children, strict=True)
        ):
            key = -1 if gene.score_keys is None else gene.score_keys[index]
            if not children:
                combined = cached.get(key) if key >= 0 else None
                if combined is None:
                    options = species_index.occurrences.get(gene.taxa[node])
                    if options is None:
                        raise ValueError(
                            f"Gene tip does not match a species: {node.name}"
                        )
                    combined = dict.fromkeys(options, 0)
                    if key >= 0:
                        cached[key] = combined
                states.append(combined)
                continue
            left, right = states[children[0]], states[children[1]]
            state_pairs += len(left) * len(right)
            if state_pairs > max_state_pairs:
                raise ValueError("Reconciliation exceeds --max-state-pairs.")
            combined = cached.get(key) if key >= 0 else None
            if combined is None:
                combined = _combine_scores(left, right, depth, lca)
                if key >= 0:
                    cached[key] = combined
            states.append(combined)
            # Each gene subtree is consumed exactly once in postorder.
            states[children[0]] = states[children[1]] = _CONSUMED_SCORE
        yield min(score + depth[mapped] for mapped, score in states[-1].items())


def _combine_scores(left, right, depth, lca):
    combined: dict[int, int] = {}
    for first, score_left in left.items():
        partial = score_left + depth[first]
        for second, score_right in right.items():
            parent = lca(first, second)
            cost = partial + score_right + depth[second] - 2 * depth[parent]
            # A speciation adds -2 to the edge losses, a duplication +1.
            cost += 1 if parent in (first, second) else -2
            previous = combined.get(parent)
            if previous is None or cost < previous:
                combined[parent] = cost
    return combined


class Reconciliation:
    """State = gene node and its species-node map; costs compose through LCA."""

    def __init__(
        self, gene, species, parser, max_state_pairs=10000000, *, _species_index=None
    ):
        gene_index = _GeneIndex(gene, parser)
        if _species_index is None:
            _species_index = _SpeciesIndex(species)
        self._solve(gene_index, _species_index, max_state_pairs)

    @classmethod
    def _from_indexes(
        cls, gene, species, max_state_pairs, *, keep_choices=True, cache=None
    ):
        result = cls.__new__(cls)
        result._solve(
            gene, species, max_state_pairs, keep_choices=keep_choices, cache=cache
        )
        return result

    def _solve(self, gene, species, max_state_pairs, *, keep_choices=True, cache=None):
        self._keep_choices = keep_choices
        self._species_index = species
        self.gene = gene.gene
        self.nodes = gene.nodes
        self.gene_index = gene.indices
        self.lca = species.lca
        self.species = self.lca.nodes
        self.depth = self.lca.depth
        self.labels = species.labels
        self.states: list[dict[int, Cell]] = []
        self.combinations = 1
        self.ambiguous_tips = 0
        self.state_pairs = 0
        keys = gene.score_keys if cache is not None else None
        for index, (node, children) in enumerate(
            zip(self.nodes, gene.children, strict=True)
        ):
            if not children:
                options = species.occurrences.get(gene.taxa[node], ())
                if not options:
                    raise ValueError(f"Gene tip does not match a species: {node.name}")
                self.combinations *= len(options)
                self.ambiguous_tips += len(options) > 1
                self.states.append({index: Cell(0, 1, None) for index in options})
            else:
                left, right = self.states[children[0]], self.states[children[1]]
                self.state_pairs += len(left) * len(right)
                if self.state_pairs > max_state_pairs:
                    raise ValueError("Reconciliation exceeds --max-state-pairs.")
                key = -1 if keys is None else keys[index]
                combined = cache.get(key) if key >= 0 else None
                if combined is None:
                    combined = self._combine(left, right, keep_choices)
                    if key >= 0:
                        cache[key] = combined
                self.states.append(combined)
                if not keep_choices:
                    self.states[children[0]] = self.states[children[1]] = (
                        _CONSUMED_STATE
                    )
        # GRAMPA includes missing lineages above the mapped gene root.
        self.score = min(
            cell.score + self.depth[index] for index, cell in self.states[-1].items()
        )
        self.roots = tuple(
            index
            for index, cell in self.states[-1].items()
            if cell.score + self.depth[index] == self.score
        )
        self.num_maps = sum(self.states[-1][index].count for index in self.roots)

    def local_cost(self, parent, left, right):
        duplication = int(parent in (left, right))
        losses = self.depth[left] + self.depth[right] - 2 * self.depth[parent]
        losses += 2 * duplication - 2
        return duplication, losses

    def _combine(self, left, right, keep_choices):
        result: dict[int, Cell] = {}
        for first, second in product(left, right):
            parent = self._species_index.common_ancestor_indices(first, second)
            duplication, losses = self.local_cost(parent, first, second)
            score = left[first].score + right[second].score + duplication + losses
            count = left[first].count * right[second].count
            current = result.get(parent)
            if current is None or score < current.score:
                result[parent] = Cell(
                    score, count, [(first, second)] if keep_choices else None
                )
            elif score == current.score:
                current.count += count
                if current.choices is not None:
                    current.choices.append((first, second))
        return result

    def _trace(self, index, mapped, rank):
        # Unrank a backpointer path without recursion, including deep comb trees.
        mapping = {}
        pending = [(index, mapped, rank)]
        while pending:
            index, mapped, rank = pending.pop()
            mapping[index] = mapped
            node = self.nodes[index]
            if node.is_leaf:
                continue
            left, right = (self.gene_index[child] for child in node.children)
            choices = self.states[index][mapped].choices
            assert choices is not None
            for first, second in choices:
                right_count = self.states[right][second].count
                count = self.states[left][first].count * right_count
                if rank >= count:
                    rank -= count
                    continue
                left_rank, right_rank = divmod(rank, right_count)
                pending.append((right, second, right_rank))
                pending.append((left, first, left_rank))
                break
        return mapping

    def mappings(self, max_maps=100000):
        if not self._keep_choices:
            raise ValueError("Reconciliation did not retain optimal mappings.")
        if self.num_maps > max_maps:
            raise ValueError(
                "Optimal mapping output exceeds --max-maps; no maps were truncated."
            )
        for root in self.roots:
            for rank in range(self.states[-1][root].count):
                mapping = self._trace(len(self.nodes) - 1, root, rank)
                duplications, losses = 0, self.depth[root]
                node_rows = []
                for index, node in enumerate(self.nodes):
                    dup, loss = 0, 0
                    if not node.is_leaf:
                        first, second = (
                            mapping[self.gene_index[c]] for c in node.children
                        )
                        dup, loss = self.local_cost(mapping[index], first, second)
                        duplications += dup
                        losses += loss
                    node_rows.append(
                        {
                            "gene_node": index,
                            "gene_label": node.name,
                            "mul_node": mapping[index],
                            "mul_label": self.labels[mapping[index]],
                            "duplication": dup,
                            "child_edge_losses": loss,
                            "root_losses": self.depth[root] if node is self.gene else 0,
                        }
                    )
                yield duplications, losses, node_rows
