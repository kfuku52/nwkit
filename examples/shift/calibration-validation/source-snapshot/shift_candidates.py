"""Enumerate identifiable small-tree shift locations and optimum equalities.

No likelihood or score formula is duplicated here. Branch IDs are NWKIT's
level-order IDs, including zero for the baseline optimum.
"""

import itertools
import math
from typing import Any

from nwkit.util import assign_branch_ids


def set_partitions(values):
    if not values:
        yield ()
        return
    first, *rest = values
    for groups in set_partitions(rest):
        yield ((first,), *groups)
        for index, group in enumerate(groups):
            yield (*groups[:index], (first, *group), *groups[index + 1 :])


def canonical_groups(groups):
    return tuple(sorted(tuple(sorted(group)) for group in groups))


def tip_groups(tree, selected, groups=None):
    ids = assign_branch_ids(tree)
    aliases = (
        {branch: branch for branch in (0, *selected)}
        if groups is None
        else {branch: min(group) for group in groups for branch in group}
    )
    labels: dict[Any, int] = {}
    tips: dict[int, list[str]] = {}
    for node in tree.traverse("preorder"):
        labels[node] = aliases[ids[node]] if ids[node] in aliases else labels[node.up]
        if node.is_leaf:
            tips.setdefault(labels[node], []).append(node.name)
    return tuple(sorted(tuple(sorted(names)) for names in tips.values()))


def has_redundant_change(tree, selected, groups):
    ids = assign_branch_ids(tree)
    aliases = {branch: min(group) for group in groups for branch in group}
    labels = {}
    for node in tree.traverse("preorder"):
        if node.is_root:
            labels[node] = aliases[0]
            continue
        labels[node] = aliases[ids[node]] if ids[node] in selected else labels[node.up]
        if ids[node] in selected and labels[node] == labels[node.up]:
            return True
    return False


def enumerate_candidates(tree, max_shifts=2, limit=5000):
    if max_shifts not in (0, 1, 2) or limit < 1:
        raise ValueError(
            "Joint reference supports 0..2 shifts and a positive candidate limit"
        )
    nodes = list(tree.traverse())
    if len(list(tree.leaves())) > 16:
        raise ValueError("Joint reference is limited to at most 16 tips")
    if any(len(node.children) != 2 for node in nodes if not node.is_leaf):
        raise ValueError("Joint reference requires a fully bifurcating tree")
    if any(
        not math.isfinite(node.dist) or node.dist <= 0
        for node in nodes
        if not node.is_root
    ):
        raise ValueError("Joint reference requires finite positive branches")
    heights = [tree.get_distance(tree, tip) for tip in tree.leaves()]
    if max(heights) - min(heights) > max(heights) * 1e-8:
        raise ValueError("Joint reference requires an ultrametric tree")
    ids = assign_branch_ids(tree)
    branches = sorted(branch for branch in ids.values() if branch)
    candidates: list[dict[str, Any]] = []
    counts = {
        "configurations_considered": 0,
        "unidentifiable_configurations": 0,
        "redundant_equalities": 0,
    }
    for size in range(max_shifts + 1):
        for selected in itertools.combinations(branches, size):
            counts["configurations_considered"] += 1
            # Empty ancestral regimes cause rank-deficient unconstrained designs.
            # With <=2 shifts this removes the pair of root-child edges.
            if len(tip_groups(tree, selected)) != size + 1:
                counts["unidentifiable_configurations"] += 1
                continue
            for groups in sorted(
                canonical_groups(p) for p in set_partitions((0, *selected))
            ):
                if has_redundant_change(tree, selected, groups):
                    counts["redundant_equalities"] += 1
                    continue
                candidates.append(
                    {
                        "candidate_id": len(candidates),
                        "shift_branch_ids": list(selected),
                        "groups": [list(group) for group in groups],
                    }
                )
                if len(candidates) > limit:
                    raise ValueError("Joint candidate limit exceeded before fitting")
    return candidates, counts
