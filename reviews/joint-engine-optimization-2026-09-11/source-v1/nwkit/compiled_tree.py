"""Immutable array-oriented views of ETE trees.

The numerical engines in :mod:`nwkit` repeatedly need the same topology but
different model parameters.  ``CompiledTree`` records that topology once,
without copying or mutating the input tree.  Callers must discard the compiled
view after changing the tree.
"""

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True, slots=True)
class CompiledTree:
    """Preorder-indexed topology shared by Gaussian and CTMC calculations."""

    tree: Any
    nodes: tuple[Any, ...]
    index_by_node: dict[Any, int]
    parents: tuple[int, ...]
    children: tuple[tuple[int, ...], ...]
    postorder: tuple[int, ...]
    leaf_indices: tuple[int, ...]
    leaf_index_by_name: dict[str, int]

    @classmethod
    def from_tree(cls, tree) -> "CompiledTree":
        nodes = tuple(tree.traverse(strategy="preorder"))
        if not nodes or nodes[0] is not tree:
            raise ValueError("Failed to compile the input tree in preorder.")
        index_by_node = {node: index for index, node in enumerate(nodes)}
        if len(index_by_node) != len(nodes):
            raise ValueError("A compiled tree cannot contain repeated node objects.")

        parents = [-1] * len(nodes)
        child_lists: list[list[int]] = [[] for _ in nodes]
        leaf_indices = []
        leaf_index_by_name: dict[str, int] = {}
        for index, node in enumerate(nodes):
            if index:
                try:
                    parent = index_by_node[node.up]
                except KeyError as exc:
                    raise ValueError(
                        "A compiled tree contains a node whose parent is outside the tree."
                    ) from exc
                if parent >= index:
                    raise ValueError("Compiled-tree preorder must place parents first.")
                parents[index] = parent
                child_lists[parent].append(index)
            if node.is_leaf:
                name = str(node.name)
                if name in leaf_index_by_name:
                    raise ValueError(
                        "A compiled tree requires unique leaf names: {}.".format(name)
                    )
                leaf_indices.append(index)
                leaf_index_by_name[name] = index

        return cls(
            tree=tree,
            nodes=nodes,
            index_by_node=index_by_node,
            parents=tuple(parents),
            children=tuple(tuple(items) for items in child_lists),
            postorder=tuple(range(len(nodes) - 1, -1, -1)),
            leaf_indices=tuple(leaf_indices),
            leaf_index_by_name=leaf_index_by_name,
        )

    def require_tree(self, tree) -> None:
        if tree is not self.tree:
            raise ValueError("A compiled tree can only be reused with its source tree.")
