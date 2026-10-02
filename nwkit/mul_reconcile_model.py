"""Exact LCA duplication/loss reconciliation over ambiguous MUL-tree tips."""

import re
from dataclasses import dataclass
from itertools import product
from typing import Any

from nwkit.clade_index import LcaIndex
from nwkit.rooting_state import require_rooted
from nwkit.util import validate_unique_named_leaves


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


def select_nodes(tree, selectors):
    nodes = tuple(tree.traverse("postorder"))
    if selectors is None:
        return nodes
    lookup = {node.name: node for node in nodes}
    clades = {frozenset(node.leaf_names()): node for node in nodes}
    result = []
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
            node = clades.get(frozenset(members))
            if node is None:
                raise ValueError(f"Unknown or non-monophyletic hypothesis: {selector}")
        if node not in result:
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
    original = tree.copy()
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
    pairs: list[tuple[str, str]] = []
    for first in selected_h1:
        for second in selected_h2:
            if first is not second and first in second.ancestors():
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
        candidate = original.copy()
        lookup = {node.name: node for node in candidate.traverse()}
        copied = lookup[first].copy()
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
        from ete4 import Tree

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
    choices: list[tuple[int, int]]


class Reconciliation:
    """State = gene node and its species-node map; costs compose through LCA."""

    def __init__(self, gene, species, parser, max_state_pairs=10000000):
        validate_binary(gene, "Gene tree")
        validate_binary(species, "MUL-tree", unique=False)
        self.gene = gene
        self.nodes = tuple(gene.traverse("postorder"))
        self.gene_index = {node: index for index, node in enumerate(self.nodes)}
        self.lca = LcaIndex(species)
        self.species = self.lca.nodes
        self.depth = self.lca.depth
        self.labels = tuple(node.name for node in self.species)
        self.states: list[dict[int, Cell]] = []
        occurrences: dict[str, list[int]] = {}
        for index, node in enumerate(self.species):
            if node.is_leaf:
                name = node.props.get("mul_species", node.name)
                occurrences.setdefault(name, []).append(index)
        self.combinations = 1
        self.ambiguous_tips = 0
        self.state_pairs = 0
        for node in self.nodes:
            if node.is_leaf:
                name = parser.parse(node.name).species_label
                options = occurrences.get(name, ())
                if not options:
                    raise ValueError(f"Gene tip does not match a species: {node.name}")
                self.combinations *= len(options)
                self.ambiguous_tips += len(options) > 1
                self.states.append({index: Cell(0, 1, []) for index in options})
            else:
                left, right = (self.states[self.gene_index[c]] for c in node.children)
                self.state_pairs += len(left) * len(right)
                if self.state_pairs > max_state_pairs:
                    raise ValueError("Reconciliation exceeds --max-state-pairs.")
                self.states.append(self._combine(left, right))
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

    def _combine(self, left, right):
        result: dict[int, Cell] = {}
        for first, second in product(left, right):
            parent = self.lca.common_ancestor_indices(first, second)
            duplication, losses = self.local_cost(parent, first, second)
            score = left[first].score + right[second].score + duplication + losses
            count = left[first].count * right[second].count
            current = result.get(parent)
            if current is None or score < current.score:
                result[parent] = Cell(score, count, [(first, second)])
            elif score == current.score:
                current.count += count
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
            for first, second in self.states[index][mapped].choices:
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
