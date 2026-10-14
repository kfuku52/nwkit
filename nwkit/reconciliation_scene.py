"""Topological two-level geometry, independent of matplotlib and inference."""

from __future__ import annotations

import math
import textwrap
from collections import Counter, defaultdict
from dataclasses import dataclass, field

from nwkit.reconciliation_document import ReconciliationDocument
from nwkit.reconciliation_labels import model_label, source_label, summary_label


@dataclass
class Primitive:
    kind: str
    id: str
    points: list[tuple[float, float]]
    text: str = ""
    color: str = "#27364b"
    width: float = 1.2
    marker: str = "o"
    owner: tuple[str, str] | None = None
    panel: str | None = None
    tree: str | None = None


@dataclass
class Scene:
    width: float
    height: float
    primitives: list[Primitive] = field(default_factory=list)
    event_points: dict[tuple[str, ...], tuple[float, float]] = field(
        default_factory=dict
    )
    caption: str = "Topology layout; event positions do not represent times."
    notices: list[str] = field(default_factory=list)


EVENT_STYLE = {
    "leaf": ("#27364b", "o"),
    "speciation": ("#0072b2", "o"),
    "duplication": ("#d55e00", "s"),
    "loss": ("#b91c1c", "x"),
    "transfer": ("#9c4caa", ">"),
    "transfer_out": ("#9c4caa", ">"),
    "transfer_in": ("#9c4caa", "<"),
    "origin": ("#15803d", "*"),
    "unresolved": ("#6b7280", "X"),
    "empty_observation": ("#6b7280", "D"),
}


def _signatures(tree):
    result: dict[str, tuple[str, ...]] = {}
    for nid in reversed(tree.order):
        node = tree.nodes[nid]
        result[nid] = (
            tuple(sorted(label for c in node["children"] for label in result[c]))
            if node["children"]
            else (node["label"],)
        )
    return result


def _coordinates(tree, width, weights=None, spacing=12, leaf_order=None):
    signatures = _signatures(tree)
    depth, leaves, stack = {}, [], [(tree.root, 1)]
    while stack:
        nid, d = stack.pop()
        depth[nid] = d
        children = sorted(tree.nodes[nid]["children"], key=lambda c: (signatures[c], c))
        if not children:
            leaves.append(nid)
        stack.extend((c, d + 1) for c in reversed(children))
    if leaf_order is not None:
        leaves = list(leaf_order)
    y, current = {}, 0.0
    for nid in leaves:
        lines = max(1, len(textwrap.wrap(tree.nodes[nid]["label"], 24)))
        weight = max(lines, (weights or {}).get(nid, 1))
        y[nid] = current + spacing * weight
        current += 2 * spacing * weight
    for nid in reversed(tree.order):
        children = tree.nodes[nid]["children"]
        if children:
            y[nid] = sum(y[c] for c in children) / len(children)
    unit = (width - 160) / max(depth.values())
    xy = {nid: (30 + d * unit, y[nid]) for nid, d in depth.items()}
    return xy, unit, current


def _upper_paths(tree, xy, unit):
    result = {}
    for nid in tree.order:
        x, y = xy[nid]
        if nid in tree.parents:
            parent = xy[tree.parents[nid]]
            result[nid] = [parent, (parent[0], y), (x, y)]
        else:
            result[nid] = [(x - unit, y), (x, y)]
    return result


def _descendant_path(tree, start, end):
    path = [end]
    while path[-1] != start:
        if path[-1] not in tree.parents:
            raise ValueError(
                "Embedded placement changes upper branch without a recorded transfer."
            )
        path.append(tree.parents[path[-1]])
    return list(reversed(path))


def _route(upper, xy, start_branch, end_branch, start, end):
    if start_branch == end_branch:
        if end[0] < start[0] - 1e-9:
            raise ValueError(
                "Recorded event order runs backwards within an upper branch."
            )
        return [start, end]
    branches = _descendant_path(upper, start_branch, end_branch)
    start_offset = start[1] - xy[start_branch][1]
    end_offset = end[1] - xy[end_branch][1]
    points = [start, (xy[start_branch][0], start[1])]
    previous_offset = start_offset
    for index, branch in enumerate(branches[1:], 1):
        parent = xy[upper.parents[branch]]
        offset = start_offset + (end_offset - start_offset) * index / (
            len(branches) - 1
        )
        advance = xy[branch][0] - parent[0]
        shift = min(
            advance / 3,
            max(0, (previous_offset + max(abs(previous_offset), abs(offset))) / 2),
        )
        elbow = parent[0] + shift
        points.extend(
            [(elbow, parent[1] + previous_offset), (elbow, xy[branch][1] + offset)]
        )
        if branch != end_branch:
            points.append((xy[branch][0], xy[branch][1] + offset))
        previous_offset = offset
    points.append(end)
    return points


def _canonical_order(tree):
    signatures = _signatures(tree)
    result, stack = [], [tree.root]
    while stack:
        nid = stack.pop()
        result.append(nid)
        stack.extend(
            reversed(
                sorted(tree.nodes[nid]["children"], key=lambda c: (signatures[c], c))
            )
        )
    return result


def _event_geometry(rec, lower, upper, xy, unit, lanes):
    by_node, by_branch = defaultdict(list), defaultdict(list)
    rank = {nid: index for index, nid in enumerate(_canonical_order(lower))}
    event_rank = {e["id"]: index for index, e in enumerate(rec["events"])}
    for event in rec["events"]:
        if event["upper_node"] is None:
            raise ValueError(
                "Embedded view requires explicit placement; use the existing paired report for unmapped TSV events."
            )
        by_node[event["lower_node"]].append(event)
        by_branch[event["upper_node"]].append(event)
    if set(by_node) != set(lower.nodes):
        raise ValueError(
            "Embedded history requires at least one recorded event per lower node."
        )
    positions, branches = {}, {}
    for branch, events in by_branch.items():
        events.sort(key=lambda e: (rank[e["lower_node"]], event_rank[e["id"]]))
        for index, event in enumerate(events):
            x, y = xy[branch]
            fraction = (index + 1) / (len(events) + 1)
            if event["type"] == "leaf" or (
                event["type"] == "speciation"
                and event.get("annotations", {}).get("observed_arity", 2) == 2
            ):
                fraction = 1
            y += lanes[rec["id"], branch, event["lower_node"]]
            positions[event["id"]] = (x - unit + unit * fraction, y)
            branches[event["id"]] = branch
    return by_node, positions, branches


def _lower_primitives(scene, rec, lower, upper, xy, unit, offset, lanes):
    by_node, positions, branches = _event_geometry(rec, lower, upper, xy, unit, lanes)
    color = ("#287d8e", "#bb7432", "#5e5fa8", "#478443")[offset % 4]
    for nid in lower.order:
        events = by_node[nid]
        links = [(events[i - 1], events[i]) for i in range(1, len(events))]
        if nid in lower.parents:
            parent_events = by_node[lower.parents[nid]]
            links.insert(0, (parent_events[-1], events[0]))
        for before, after in links:
            start, end = positions[before["id"]], positions[after["id"]]
            transfer = before.get("recipient_lower_node") == nid and before["type"] in {
                "transfer",
                "transfer_out",
            }
            if transfer:
                if before["recipient_node"] != branches[after["id"]]:
                    raise ValueError(
                        "Transfer endpoint does not match its recorded recipient entry."
                    )
                scene.primitives.append(
                    Primitive(
                        "transfer",
                        f"transfer:{rec['id']}:{before['id']}",
                        [start, end],
                        color="#9c4caa",
                        owner=(rec["id"], nid),
                        tree=lower.id,
                    )
                )
            else:
                points = _route(
                    upper, xy, branches[before["id"]], branches[after["id"]], start, end
                )
                scene.primitives.append(
                    Primitive(
                        "line",
                        f"line:{rec['id']}:{before['id']}:{after['id']}",
                        points,
                        color=color,
                        owner=(rec["id"], nid),
                        tree=lower.id,
                    )
                )
        for event in events:
            point = positions[event["id"]]
            scene.event_points[rec["id"], event["id"]] = point
            event_color, marker = EVENT_STYLE[event["type"]]
            scene.primitives.append(
                Primitive(
                    "event",
                    f"event:{rec['id']}:{event['id']}",
                    [point],
                    text=event["type"],
                    color=event_color,
                    marker=marker,
                    owner=(rec["id"], nid),
                    tree=lower.id,
                )
            )
        if not lower.nodes[nid]["children"] and lower.nodes[nid]["label"]:
            point = positions[events[-1]["id"]]
            scene.primitives.append(
                Primitive(
                    "text",
                    f"lower-label:{rec['id']}:{nid}",
                    [(point[0] + 5, point[1])],
                    text=lower.nodes[nid]["label"],
                    color=color,
                    owner=(rec["id"], nid),
                    tree=lower.id,
                )
            )


def _standalone_scene(document, width, spacing, view=None):
    scene = Scene(width, 0)
    for tree in document.trees.values():
        xy, unit, height = _coordinates(
            tree,
            width,
            spacing=spacing,
            leaf_order=(view or {}).get("leaf_order", {}).get(tree.id),
        )
        shifted = {nid: (x, y + scene.height) for nid, (x, y) in xy.items()}
        for nid, points in _upper_paths(tree, shifted, unit).items():
            scene.primitives.append(
                Primitive(
                    "line",
                    f"tree:{tree.id}:{nid}",
                    points,
                    owner=(tree.id, nid),
                    tree=tree.id,
                )
            )
        for nid, node in tree.nodes.items():
            if node["label"]:
                x, y = shifted[nid]
                scene.primitives.append(
                    Primitive(
                        "text",
                        f"label:{tree.id}:{nid}",
                        [(x + 5, y)],
                        text=node["label"],
                        owner=(tree.id, nid),
                        tree=tree.id,
                    )
                )
        scene.height += height + 35
    return scene


def _display_selection(scene, document):
    selection = document.record["provenance"].get("display_selection")
    if not selection:
        return scene
    tree = document.trees[selection["tree"]]
    selected = {nid for nid in tree.nodes if tree.ancestral(selection["node"], nid)}
    rec_ids = {
        r["id"]
        for r in document.record["reconciliations"]
        if r["lower_tree"] == tree.id
    }
    owners = rec_ids or {tree.id}
    scene.primitives = [
        p
        for p in scene.primitives
        if p.owner is None or p.tree != tree.id or p.owner[1] in selected
    ]
    event_keys = {
        (r["id"], e["id"])
        for r in document.record["reconciliations"]
        if r["id"] in owners
        for e in r["events"]
        if e["lower_node"] in selected
    }
    scene.event_points = {
        key: xy
        for key, xy in scene.event_points.items()
        if key[0] not in rec_ids or key in event_keys
    }
    scene.notices.append(
        "Display-only clade selection; scientific trees/events are unchanged. Incoming transfer arrows retain their recorded external donor context."
    )
    return scene


def build_scene(
    document: ReconciliationDocument, *, width=720.0, font_size=8, view=None
) -> Scene:
    if not math.isfinite(width) or width < 360:
        raise ValueError("Document width must be at least 5 inches.")
    if not math.isfinite(font_size) or not 0 < font_size <= 72:
        raise ValueError("Font size must be positive and at most 72 points.")
    spacing = max(12, font_size * 1.5)
    recs = document.record["reconciliations"]
    upper_ids = {r["upper_tree"] for r in recs}
    lower_ids = {r["lower_tree"] for r in recs}
    if (view or {}).get("layout") == "standalone" or not recs:
        scene = _standalone_scene(document, width, spacing, view)
        scene.notices.append(source_label(document.record))
        return _display_selection(scene, document)
    if _needs_composition(document, upper_ids, lower_ids, view):
        from nwkit.reconciliation_composition import compose_document

        return _display_selection(
            compose_document(document, width=width, font_size=font_size, view=view),
            document,
        )
    if any(r["history_kind"] == "expectation" for r in recs):
        raise ValueError(
            "Expectation maps require a separate view; N2 renders histories, not expectation values as events."
        )
    upper = document.trees[next(iter(upper_ids))]
    occupancy: Counter = Counter()
    for rec in recs:
        lower = document.trees[rec["lower_tree"]]
        for event in rec["events"]:
            if event["type"] == "leaf":
                label = lower.nodes[event["lower_node"]]["label"]
                occupancy[event["upper_node"]] += max(1, len(textwrap.wrap(label, 24)))
    xy, unit, height = _coordinates(
        upper,
        width,
        {nid: max(1, count) for nid, count in occupancy.items()},
        spacing,
        leaf_order=(view or {}).get("leaf_order", {}).get(upper.id),
    )
    scene = Scene(width, height)
    by_branch = defaultdict(set)
    lower_by_rec = {r["id"]: document.trees[r["lower_tree"]] for r in recs}
    for rec in recs:
        for event in rec["events"]:
            by_branch[event["upper_node"]].add((rec["id"], event["lower_node"]))
    lanes, widths = {}, {}
    for branch, members in by_branch.items():
        rows = max(
            max(1, len(textwrap.wrap(lower_by_rec[rid].nodes[nid]["label"], 24)))
            for rid, nid in members
        )
        widths[branch] = spacing * rows * len(members) + 8
        for index, (rid, nid) in enumerate(sorted(members)):
            lanes[rid, branch, nid] = spacing * rows * (index - (len(members) - 1) / 2)
    for nid, points in _upper_paths(upper, xy, unit).items():
        scene.primitives.append(
            Primitive(
                "upper",
                f"upper:{upper.id}:{nid}",
                points,
                color="#e1e8ee",
                width=max(12, widths.get(nid, 12)),
                owner=(upper.id, nid),
                tree=upper.id,
            )
        )
        node = upper.nodes[nid]
        x, y = xy[nid]
        if node["label"]:
            scene.primitives.append(
                Primitive(
                    "text",
                    f"upper-label:{nid}",
                    [(x, y - max(9, widths.get(nid, 12) / 2 + 5))],
                    text=node["label"],
                    color="#596b80",
                    owner=(upper.id, nid),
                    tree=upper.id,
                )
            )
    for index, rec in enumerate(sorted(recs, key=lambda r: r["id"])):
        _lower_primitives(
            scene, rec, document.trees[rec["lower_tree"]], upper, xy, unit, index, lanes
        )
        scene.notices.append(model_label(rec))
        if (
            rec["selection"].get("measurement")
            == "partially_marginalized_ordered_derivation_v1"
        ):
            scene.notices.append(
                "Nonobserved side histories are integrated in E; this representative is not full R MAP or marginal H MAP, and its events are not posterior support."
            )
        if rec.get("summary"):
            scene.notices.append(summary_label(rec))
    ys = [y for p in scene.primitives for _, y in p.points]
    low = min(ys) if ys else 0
    for primitive in scene.primitives:
        primitive.points = [(x, y - low + 15) for x, y in primitive.points]
    scene.event_points = {
        key: (x, y - low + 15) for key, (x, y) in scene.event_points.items()
    }
    scene.height = (max(ys) - low + 35) if ys else height
    scene.notices.append(source_label(document.record))
    scene.notices.append(
        "Retained branch units (not used as layout time): "
        + ", ".join(
            f"{t.id}={t.record['branch_length_unit']}" for t in document.trees.values()
        )
    )
    return _display_selection(scene, document)


def _needs_composition(document, upper_ids, lower_ids, view):
    return bool(
        document.record["transformations"]
        or len(upper_ids) != 1
        or upper_ids & lower_ids
        or set(document.trees) != upper_ids | lower_ids
        or (view or {}).get("layout") == "paired"
    )
