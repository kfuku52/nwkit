from collections import defaultdict
from typing import Any

"""Shared topology features for native trees and external tree adapters."""


def tree_children(node):
    return node.clades if hasattr(node, "clades") else node.children


def _tree_order(tree):
    root = tree.root if hasattr(tree, "root") else tree
    order, pending = [], [root]
    while pending:
        node = pending.pop()
        order.append(node)
        pending.extend(reversed(tree_children(node)))
    return root, order


def descendant_sets(tree):
    root, order = _tree_order(tree)
    result: dict[Any, set] = {}
    for node in reversed(order):
        children = tree_children(node)
        result[node] = (
            set().union(*(result[c] for c in children)) if children else {node.name}
        )
    return root, order, result


def unrooted_splits(tree, taxa):
    taxa = set(taxa)
    _, _, descendants = descendant_sets(tree)
    result = set()
    for tips in descendants.values():
        left = frozenset(tips & taxa)
        right = frozenset(taxa - left)
        if min(len(left), len(right)) >= 2:
            result.add(min((left, right), key=lambda x: (len(x), sorted(x))))
    return result


def verify_constraint(tree, constraint, taxa, scope="complete"):
    if scope not in {"complete", "partial"}:
        raise ValueError("Constraint scope must be complete or partial")
    taxa = set(taxa)
    _, order = _tree_order(constraint)
    labels = _validate_tip_labels(order)
    _, result_order = _tree_order(tree)
    if _validate_tip_labels(result_order) != taxa:
        raise ValueError("Result tree and selected taxa differ")
    if labels - taxa or (scope == "complete" and labels != taxa):
        raise ValueError("Constraint and selected taxa differ")
    required = unrooted_splits(constraint, labels)
    if not required:
        raise ValueError(
            "Constraint must specify at least one nontrivial unrooted split"
        )
    missing = required - unrooted_splits(tree, labels)
    if missing:
        raise ValueError(
            f"Result tree violates {len(missing)} required unrooted constraint splits"
        )
    return required


def _validate_tip_labels(order):
    names = [node.name for node in order if not tree_children(node)]
    if any(not name for name in names) or len(names) != len(set(names)):
        raise ValueError("Tree tips must have unique nonempty names")
    return set(names)


def foreground_components(tree, foreground):
    foreground = set(foreground)
    root, order = _tree_order(tree)
    unknown = foreground - _validate_tip_labels(order)
    if unknown:
        raise ValueError(f"Foreground tips absent from tree: {sorted(unknown)}")
    all_foreground: dict[Any, bool] = {}
    for node in reversed(order):
        children = tree_children(node)
        all_foreground[node] = (
            all(all_foreground[child] for child in children)
            if children
            else node.name in foreground
        )
    components, pending = [], [root]
    while pending:
        node = pending.pop()
        if all_foreground[node]:
            names, clade = [], [node]
            while clade:
                child = clade.pop()
                children = tree_children(child)
                if children:
                    clade.extend(children)
                else:
                    names.append(child.name)
            components.append(tuple(sorted(names)))
        else:
            pending.extend(tree_children(node))
    return sorted(components)


def patristic_distances(tree):
    """All tip distances in O(N^2), without repeated whole-tree LCA scans."""
    graph = defaultdict(list)
    root = tree.root if hasattr(tree, "root") else tree
    order, pending = [], [root]
    while pending:
        parent = pending.pop()
        order.append(parent)
        pending.extend(reversed(tree_children(parent)))
        for child in tree_children(parent):
            length = (
                getattr(child, "branch_length", None)
                if hasattr(child, "clades")
                else child.dist
            ) or 0.0
            graph[parent].append((child, length))
            graph[child].append((parent, length))
    _validate_tip_labels(order)
    leaves = {
        tip: tip.name for tip in (node for node in order if not tree_children(node))
    }
    result = {}
    for tip, name in leaves.items():
        row, stack = {}, [(tip, None, 0.0)]
        while stack:
            node, previous, distance = stack.pop()
            if node in leaves:
                row[leaves[node]] = distance
            stack.extend(
                (child, node, distance + length)
                for child, length in graph[node]
                if child is not previous
            )
        result[name] = row
    return result


def projected_clade_candidates(wanted, records, shared_taxa, original_taxa):
    wanted, shared_taxa, original_taxa = (
        set(wanted),
        set(shared_taxa),
        set(original_taxa),
    )
    return [
        {**row, "added_descendants": sorted(set(row["descendants"]) - original_taxa)}
        for row in records
        if not row["is_root"] and set(row["descendants"]) & shared_taxa == wanted
    ]
