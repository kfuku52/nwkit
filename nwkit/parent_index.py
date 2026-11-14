"""Validated parent tables without renumbering caller-owned node identifiers."""

from collections.abc import Mapping

from ete4 import Tree

from nwkit.clade_index import LcaIndex


def traversal_order(children, root):
    """Preorder with explicit child order, rejecting repeated/cyclic nodes."""
    order, pending, seen = [], [root], set()
    while pending:
        node = pending.pop()
        if node in seen:
            raise ValueError(f"Tree contains a cycle or repeated child: {node}")
        seen.add(node)
        order.append(node)
        pending.extend(reversed(children.get(node, ())))
    return order


def validate_parent_links(parents, children, root):
    if root not in parents:
        raise ValueError("Tree root is absent from parent table")
    if parents[root] in parents:
        raise ValueError("Tree root must not have a parent in the table")
    for node, child_ids in children.items():
        if node not in parents:
            raise ValueError(f"Unknown parent: {node}")
        for child in child_ids:
            if child == root or child not in parents or parents[child] != node:
                raise ValueError(f"Child and parent links disagree: {node}, {child}")
    order = traversal_order(children, root)
    missing = set(parents) - set(order)
    if missing:
        raise ValueError(
            f"Tree contains branches disconnected from its root: {sorted(missing)}"
        )
    return order


class ParentIndex(Mapping):
    """Lazy ancestor chains and O(log n) LCA using the existing LcaIndex."""

    def __init__(self, parents, root=None, children=None):
        self.parents = dict(parents)
        if root is None:
            roots = [node for node, parent in parents.items() if parent not in parents]
            if len(roots) != 1:
                raise ValueError("Parent table must have exactly one root")
            root = roots[0]
        if children is None:
            children = {}
            for node, parent in parents.items():
                if node != root:
                    children.setdefault(parent, []).append(node)
        self.order = validate_parent_links(parents, children, root)
        self.nodes = {node: Tree() for node in self.order}
        for node in self.order:
            if node != root:
                self.nodes[parents[node]].add_child(self.nodes[node])
        self.ids = {node: key for key, node in self.nodes.items()}
        self.lca = LcaIndex(self.nodes[root])

    def __len__(self):
        return len(self.parents)

    def __iter__(self):
        return iter(self.parents)

    def __getitem__(self, node):
        if node not in self.parents:
            raise KeyError(node)
        chain = []
        while node in self.parents:
            chain.append(node)
            node = self.parents[node]
        return chain

    def common_ancestor(self, left, right):
        return self.ids[self.lca.common_ancestor(self.nodes[left], self.nodes[right])]
