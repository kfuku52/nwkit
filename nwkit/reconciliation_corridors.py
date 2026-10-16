"""Display geometry for species corridors and required terminal extinctions."""

from __future__ import annotations

import textwrap
from bisect import bisect_left, bisect_right
from collections import defaultdict


def corridor_lanes(recs, trees, spacing):
    """Share lanes across an ancestor and its descendant copy branches."""
    from nwkit.reconciliation_scene import _canonical_order

    ranks: dict[str, dict[str, int]] = {}
    ends: dict[str, dict[str, int]] = {}
    grouped: defaultdict[str, defaultdict[str, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )
    empty: defaultdict[str, defaultdict[str, list[dict]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for rec in recs:
        rid, tree = rec["id"], trees[rec["lower_tree"]]
        order = _canonical_order(tree)
        ranks[rid] = {nid: i for i, nid in enumerate(order)}
        ends[rid] = {}
        for nid in reversed(order):
            ends[rid][nid] = max(
                [ranks[rid][nid]] + [ends[rid][c] for c in tree.nodes[nid]["children"]]
            )
        for event in rec["events"]:
            grouped[event["upper_node"]][rid].add(event["lower_node"])
        for event, branch in marginalized_terminal_sides(rec, trees[rec["upper_tree"]]):
            grouped[branch][rid]  # Keep empty-side lanes separate from H ancestry.
            empty[branch][rid].append(event)
    lanes: dict[tuple[str, ...], float] = {}
    widths = {}
    lower = {r["id"]: trees[r["lower_tree"]] for r in recs}
    for branch, histories in grouped.items():
        frontiers = {}
        for rid, members in histories.items():
            ordered = sorted(members, key=ranks[rid].__getitem__)
            frontiers[rid] = [
                nid
                for i, nid in enumerate(ordered)
                if i == len(ordered) - 1 or ranks[rid][ordered[i + 1]] > ends[rid][nid]
            ]
        count = sum(
            len(nodes) + len(empty[branch][rid]) for rid, nodes in frontiers.items()
        )
        rows = max(
            max(1, len(textwrap.wrap(lower[rid].nodes[nid]["label"], 24)))
            for rid, members in histories.items()
            for nid in members | {e["lower_node"] for e in empty[branch][rid]}
        )
        widths[branch] = spacing * rows * count + 8
        base = 0
        for rid in sorted(histories):
            indices = [ranks[rid][nid] for nid in frontiers[rid]]
            for nid in histories[rid]:
                first = bisect_left(indices, ranks[rid][nid])
                last = bisect_right(indices, ends[rid][nid]) - 1
                center = base + (first + last) / 2 - (count - 1) / 2
                lanes[rid, branch, nid] = spacing * rows * center
            for index, event in enumerate(empty[branch][rid]):
                center = base + len(indices) + index - (count - 1) / 2
                lanes[rid, branch, event["id"], "empty"] = spacing * rows * center
            base += len(indices) + len(empty[branch][rid])
    return lanes, widths


def corridor_centers(tree, leaves, widths, spacing):
    """Pack subtree envelopes so neighboring corridors retain visible space."""
    rank = {nid: i for i, nid in enumerate(leaves)}
    ranges, heights, centers, offsets = {}, {}, {}, {}
    for nid in reversed(tree.order):
        children = tree.nodes[nid]["children"]
        if not children:
            ranges[nid] = (rank[nid], rank[nid], 1)
            lines = max(1, len(textwrap.wrap(tree.nodes[nid]["label"], 24)))
            heights[nid] = max(widths.get(nid, 12) + spacing, 2 * spacing * lines)
            centers[nid] = heights[nid] / 2
            continue
        children = sorted(children, key=lambda c: ranges[c][0])
        ranges[nid] = (
            min(ranges[c][0] for c in children),
            max(ranges[c][1] for c in children),
            sum(ranges[c][2] for c in children),
        )
        cursor = 0.0
        for child in children:
            offsets[child] = cursor
            cursor += heights[child] + spacing
        cursor -= spacing
        center = sum(offsets[c] + centers[c] for c in children) / len(children)
        half = widths.get(nid, 12) / 2 + spacing / 2
        low, high = min(0, center - half), max(cursor, center + half)
        for child in children:
            offsets[child] -= low
        centers[nid], heights[nid] = center - low, high - low
    if any(hi - lo + 1 != count for lo, hi, count in ranges.values()):
        # An explicit interleaved leaf order cannot be packed as planar subtrees.
        # Preserve that order with conservative clearance rather than reorder it.
        step = max(widths.values(), default=12) + 2 * spacing
        y = {nid: (i + 0.5) * step for i, nid in enumerate(leaves)}
        for nid in reversed(tree.order):
            children = tree.nodes[nid]["children"]
            if children:
                y[nid] = sum(y[c] for c in children) / len(children)
        return y, len(leaves) * step
    starts, y = {tree.root: 0.0}, {}
    for nid in tree.order:
        y[nid] = starts[nid] + centers[nid]
        for child in tree.nodes[nid]["children"]:
            starts[child] = starts[nid] + offsets[child]
    return y, heights[tree.root]


def required_terminal_losses(rec, upper):
    """Only interpret the documented aioi pure-loss terminal boundary case.

    General empty observations are not biological losses. Even here the stub
    represents required extinction; its midpoint is never an inferred time.
    """
    model = rec["model"]
    conditions = model.get("conditions", {})
    if (
        not isinstance(model.get("id"), str)
        or model.get("id") not in {"undated_DL", "undated_DTL", "undated_DTLX"}
        or rec["history_kind"] != "representative"
        or rec["selection"].get("measurement")
        != "partially_marginalized_ordered_derivation_v1"
        or not isinstance(conditions, dict)
        or any(
            isinstance(conditions.get(key), bool) or conditions.get(key) != 0
            for key in ("duplication", "transfer")
        )
    ):
        return []
    detection = conditions.get("detection", {})
    mapping = conditions.get("mapping", {})
    if not isinstance(detection, dict) or not isinstance(mapping, dict) or not mapping:
        return []
    result = []
    for event in rec["events"]:
        annotations = event.get("annotations", {})
        branch = annotations.get("marginalized_species")
        if (
            event["type"] != "speciation"
            or type(annotations.get("observed_arity")) is not int
            or annotations.get("observed_arity") != 1
        ):
            continue
        if not isinstance(branch, str) or branch not in upper.nodes:
            continue
        node = upper.nodes[branch]
        if (
            not node["children"]
            and upper.parents.get(branch) == event["upper_node"]
            and detection.get(node["label"]) == 1
            and not isinstance(detection.get(node["label"]), bool)
            and node["label"] not in mapping.values()
            and isinstance(conditions.get("loss"), (int, float))
            and not isinstance(conditions.get("loss"), bool)
            and conditions["loss"] > 0
        ):
            result.append((event, branch))
    return result


def marginalized_terminal_sides(rec, upper):
    if (
        rec["history_kind"] != "representative"
        or rec["selection"].get("measurement")
        != "partially_marginalized_ordered_derivation_v1"
    ):
        return []
    result = []
    for event in rec["events"]:
        annotations = event.get("annotations", {})
        branch = annotations.get("marginalized_species")
        if (
            event["type"] == "speciation"
            and type(annotations.get("observed_arity")) is int
            and annotations["observed_arity"] == 1
            and isinstance(branch, str)
            and branch in upper.nodes
            and not upper.nodes[branch]["children"]
            and upper.parents.get(branch) == event["upper_node"]
        ):
            result.append((event, branch))
    return result


def add_marginalized_side_stubs(scene, rec, upper, xy, positions, lineage_color, lanes):
    from nwkit.reconciliation_scene import EVENT_STYLE, Primitive, _route

    losses = required_terminal_losses(rec, upper)
    loss_ids = {event["id"] for event, _ in losses}
    for event, branch in marginalized_terminal_sides(rec, upper):
        start = positions[event["id"]]
        x, y = xy[branch]
        y += lanes[rec["id"], branch, event["id"], "empty"]
        parent_x = xy[upper.parents[branch]][0]
        end = ((parent_x + x) / 2, y)
        owner = (rec["id"], event["lower_node"])
        tree = rec["lower_tree"]
        is_loss = event["id"] in loss_ids
        prefix = "integrated-loss" if is_loss else "empty-side"
        name = f"{prefix}:{rec['id']}:{event['id']}:{branch}"
        note = (
            (
                "Display interpretation: required terminal extinction with recorded "
                "D=T=0 and detection=1. The nonobserved side is integrated in E. "
                "The midpoint is a layout convention, not an inferred loss time "
                "or an additional recorded event."
            )
            if is_loss
            else (
                "Display of the source-recorded marginalized terminal side E. "
                "Biological loss and incomplete detection histories remain integrated. "
                "The midpoint is a layout convention, not an inferred event time "
                "or an additional recorded event."
            )
        )
        scene.primitives.append(
            Primitive(
                "line",
                f"{name}:stub",
                _route(upper, xy, event["upper_node"], branch, start, end),
                color=lineage_color if is_loss else "#6b7280",
                owner=owner,
                tree=tree,
                annotation=note,
            )
        )
        kind = "loss" if is_loss else "empty_observation"
        color, marker = EVENT_STYLE[kind]
        scene.primitives.append(
            Primitive(
                "event",
                name,
                [end],
                text=kind,
                color=color,
                marker=marker,
                owner=owner,
                tree=tree,
                annotation=note,
            )
        )
    if losses:
        scene.notices.append(
            "Loss stubs show required terminal extinction only for recorded "
            "aioi D=T=0, detection=1 boundaries. Their midpoint is a display "
            "convention; the side history remains integrated in E and no loss "
            "time, additional event record or posterior support is inferred."
        )
