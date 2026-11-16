from collections import defaultdict
from itertools import chain
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


def _complement_intervals(anchor_path, tip_count):
    """Summarize the sides opposite the anchor, from the root to the tip."""
    lo, hi, size = tip_count, -1, 0
    count = 0
    intervals = set()
    for first, last, sibling_size, siblings in reversed(anchor_path):
        components = siblings + bool(size)
        lo = min(lo, first)
        hi = max(hi, last)
        size += sibling_size
        # A single component repeats the previous complement or a sibling.
        if components > 1 and size <= tip_count - 2:
            count += 1
            if hi - lo + 1 == size:
                intervals.add((lo, hi))
    return count, intervals


def _split_intervals(postorder, positions, rooted):
    """Count distinct splits, indexing those contiguous in the reference order.

    Every reference tip occurs once. Extra tips are ignored in the projection.
    """
    tip_count = len(positions)
    max_size = tip_count - (1 if rooted else 2)
    empty = (tip_count, -1, 0)
    extents: dict[Any, tuple[int, int, int] | None] = {}
    count = 0
    intervals = set()
    anchor_path = []
    for node, children in postorder:
        if not children:
            index = positions.get(node.name)
            if index is None:
                extents[node] = empty  # Tip outside the projected taxon set.
            elif index == 0 and not rooted:
                extents[node] = None  # Mark the path to the reference anchor.
            else:
                extents[node] = (index, index, 1)
            continue
        lo, hi, size = empty
        has_anchor = False
        components = len(children)
        for child in children:
            extent = extents.pop(child)
            if extent is None:
                has_anchor = True
                components -= 1
                continue
            first, last, child_size = extent
            if not child_size:
                components -= 1
                continue
            if first < lo:
                lo = first
            if last > hi:
                hi = last
            size += child_size
        if has_anchor:
            anchor_path.append((lo, hi, size, components))
            extents[node] = None
        else:
            extents[node] = (lo, hi, size)
            # A single projected child displays the same split one edge below.
            if components > 1 and size <= max_size:
                count += 1
                if hi - lo + 1 == size:
                    intervals.add((lo, hi))

    if anchor_path:
        extra_count, extra_intervals = _complement_intervals(anchor_path, tip_count)
        count += extra_count
        intervals.update(extra_intervals)
    return count, intervals


def _smaller_split_side(tip_names, lo, hi):
    inside = hi - lo + 1
    outside = len(tip_names) - inside
    if inside < outside:
        return frozenset(tip_names[lo : hi + 1])
    complement = frozenset(chain(tip_names[:lo], tip_names[hi + 1 :]))
    if outside < inside:
        return complement
    side = frozenset(tip_names[lo : hi + 1])
    return min(side, complement, key=sorted)


def verify_constraint(tree, constraint, taxa, scope="complete"):
    if scope not in {"complete", "partial"}:
        raise ValueError("Constraint scope must be complete or partial")
    taxa = set(taxa)
    _, order = _tree_order(constraint)
    tip_names = _validated_tip_names(order)
    labels = set(tip_names)
    _, result_order = _tree_order(tree)
    if set(_validated_tip_names(result_order)) != taxa:
        raise ValueError("Result tree and selected taxa differ")
    if labels - taxa or (scope == "complete" and labels != taxa):
        raise ValueError("Constraint and selected taxa differ")
    positions = {name: index for index, name in enumerate(tip_names)}
    _, required = _split_intervals(
        ((node, tree_children(node)) for node in reversed(order)), positions, False
    )
    if not required:
        raise ValueError(
            "Constraint must specify at least one nontrivial unrooted split"
        )
    _, displayed = _split_intervals(
        ((node, tree_children(node)) for node in reversed(result_order)),
        positions,
        False,
    )
    missing = required - displayed
    if missing:
        raise ValueError(
            f"Result tree violates {len(missing)} required unrooted constraint splits"
        )
    return {_smaller_split_side(tip_names, lo, hi) for lo, hi in required}


def _validated_tip_names(order):
    names = [node.name for node in order if not tree_children(node)]
    if any(not name for name in names) or len(names) != len(set(names)):
        raise ValueError("Tree tips must have unique nonempty names")
    return names


def foreground_components(tree, foreground):
    foreground = set(foreground)
    root, order = _tree_order(tree)
    unknown = foreground.difference(_validated_tip_names(order))
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
    _validated_tip_names(order)
    leaves = {
        tip: tip.name for tip in (node for node in order if not tree_children(node))
    }
    result = {}
    for tip, name in leaves.items():
        row, stack = {}, [(tip, None, 0.0)]
        while stack:
            node, previous, distance = stack.pop()
            leaf_name = leaves.get(node)
            if leaf_name is not None:
                row[leaf_name] = distance
                if previous is not None:
                    continue  # A reached tip has no onward edges.
            for child, length in graph[node]:
                if child is not previous:
                    stack.append((child, node, distance + length))
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
