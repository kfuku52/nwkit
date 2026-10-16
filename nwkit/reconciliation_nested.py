"""Three layers in recorded corridors; arc length is a display coordinate only."""

from __future__ import annotations

import copy
import math
from collections import defaultdict

from nwkit.reconciliation_composition import subset
from nwkit.reconciliation_corridors import (
    add_terminal_side_notice,
    corridor_lanes,
    terminal_side_annotations,
)
from nwkit.reconciliation_labels import model_label, source_label
from nwkit.reconciliation_scene import (
    EVENT_STYLE,
    Primitive,
    _canonical_order,
    _descendant_path,
    build_scene,
)


def nested_chain(document):
    recs = document.record["reconciliations"]
    if document.record["transformations"] or len(recs) < 2:
        return None
    lower = {r["lower_tree"] for r in recs}
    roots = [r for r in recs if r["upper_tree"] not in lower]
    if len(roots) != 1:
        return None
    middle = roots[0]
    bottom = [r for r in recs if r is not middle]
    if any(r["upper_tree"] != middle["lower_tree"] for r in bottom):
        return None
    if set(document.trees) != {middle["upper_tree"], middle["lower_tree"]} | {
        r["lower_tree"] for r in bottom
    }:
        return None
    return middle, bottom


def _curve(primitive):
    if primitive.kind != "transfer":
        return primitive.points
    (x, y), (xx, yy) = primitive.points[0], primitive.points[-1]
    cx, cy = (x + xx) / 2 - (yy - y) * 0.25, (y + yy) / 2 + (xx - x) * 0.25
    return [
        (
            (1 - t) ** 2 * x + 2 * (1 - t) * t * cx + t * t * xx,
            (1 - t) ** 2 * y + 2 * (1 - t) * t * cy + t * t * yy,
        )
        for t in (i / 32 for i in range(33))
    ]


def _point(path, fraction):
    lengths = [math.dist(a, b) for a, b in zip(path[:-1], path[1:], strict=True)]
    total = sum(lengths)
    target, traveled = fraction * total, 0.0
    for a, b, length in zip(path[:-1], path[1:], lengths, strict=True):
        if length and traveled + length >= target - 1e-9:
            t = min(1, max(0, (target - traveled) / length))
            return (a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t), (
                (b[1] - a[1]) / length,
                -(b[0] - a[0]) / length,
            )
        traveled += length
    return path[-1], (0.0, -1.0)


def _offset(path, amount):
    path = [p for i, p in enumerate(path) if i == 0 or p != path[i - 1]]
    if len(path) < 2:
        return [(path[0][0], path[0][1] - amount)] * 2
    normals = []
    for a, b in zip(path[:-1], path[1:], strict=True):
        length = math.dist(a, b)
        normals.append(((b[1] - a[1]) / length, -(b[0] - a[0]) / length))
    out: list[tuple[float, float]] = []
    for i, p in enumerate(path):
        left, right = normals[max(0, i - 1)], normals[min(i, len(normals) - 1)]
        denominator = 1 + left[0] * right[0] + left[1] * right[1]
        if denominator < 1e-8:
            out.extend(
                (p[0] + amount * n[0], p[1] + amount * n[1]) for n in (left, right)
            )
        else:
            out.append(
                (
                    p[0] + amount * (left[0] + right[0]) / denominator,
                    p[1] + amount * (left[1] + right[1]) / denominator,
                )
            )
    return out


def _segment(path, first, last):
    if last < first - 1e-9:
        raise ValueError("Recorded event order runs backwards in a middle branch.")
    a, _ = _point(path, first)
    b, _ = _point(path, last)
    lengths = [math.dist(x, y) for x, y in zip(path[:-1], path[1:], strict=True)]
    total = sum(lengths)
    distance = 0.0
    out = [a]
    for point, length in zip(path[1:], lengths, strict=True):
        distance += length
        if first * total < distance < last * total:
            out.append(point)
    out.append(b)
    return out


def _middle_paths(scene, rec, tree):
    paths: defaultdict[str, list[tuple[float, float]]] = defaultdict(list)
    for primitive in scene.primitives:
        if (
            primitive.owner
            and primitive.owner[0] == rec["id"]
            and primitive.kind in {"line", "transfer"}
            and not primitive.annotation
        ):
            points = _curve(primitive)
            if paths[primitive.owner[1]]:
                paths[primitive.owner[1]].extend(points[1:])
            else:
                paths[primitive.owner[1]].extend(points)
    by_node = defaultdict(list)
    for event in rec["events"]:
        by_node[event["lower_node"]].append(event)
    for node, events in by_node.items():
        if len(paths[node]) < 2:
            end = scene.event_points[rec["id"], events[-1]["id"]]
            paths[node] = [(end[0] - 35, end[1]), end]
    if set(paths) != set(tree.nodes):
        raise ValueError(
            "Middle tree lacks a recorded incoming branch for nested placement."
        )
    return paths


def _bottom_geometry(recs, document, paths, lanes, widths):
    by_node = {}
    positions = {}
    fractions = {}
    for rec in recs:
        grouped = defaultdict(list)
        tree = document.trees[rec["lower_tree"]]
        for event in rec["events"]:
            if event["upper_node"] not in paths:
                raise ValueError(
                    "Nested events require explicit middle-tree placement."
                )
            grouped[event["lower_node"]].append(event)
        if set(grouped) != set(tree.nodes):
            raise ValueError(
                "Nested history requires recorded events for every lower node."
            )
        by_node[rec["id"]] = grouped
    forks = _nested_forks(recs, paths, lanes, widths)
    for rec in recs:
        tree = document.trees[rec["lower_tree"]]
        rank = {nid: i for i, nid in enumerate(_canonical_order(tree))}
        by_branch = defaultdict(list)
        for event in rec["events"]:
            by_branch[event["upper_node"]].append(event)
        for branch, events in by_branch.items():
            events.sort(key=lambda e: (rank[e["lower_node"]], rec["events"].index(e)))
            interior_end = min(
                (
                    forks[rec["id"], e["id"]]
                    for e in events
                    if e["type"] == "speciation"
                ),
                default=1,
            )
            for i, event in enumerate(events):
                fraction = interior_end * (i + 1) / (len(events) + 1)
                if event["type"] == "leaf":
                    fraction = 1
                elif event["type"] == "speciation":
                    fraction = forks[rec["id"], event["id"]]
                point, normal = _point(paths[branch], fraction)
                amount = lanes[rec["id"], branch, event["lower_node"]]
                positions[rec["id"], event["id"]] = (
                    point[0] + amount * normal[0],
                    point[1] + amount * normal[1],
                )
                fractions[rec["id"], event["id"]] = fraction
    return by_node, positions, fractions


def _nested_forks(recs, paths, lanes, widths):
    groups = defaultdict(list)
    for rec in recs:
        for event in rec["events"]:
            if event["type"] == "speciation":
                branch = event["upper_node"]
                groups[branch].append(
                    (
                        lanes[rec["id"], branch, event["lower_node"]],
                        rec["id"],
                        event["id"],
                    )
                )
    result = {}
    for branch, events in groups.items():
        path = paths[branch]
        total = sum(math.dist(a, b) for a, b in zip(path[:-1], path[1:], strict=True))
        tail = next(
            (
                (a, b)
                for a, b in reversed(list(zip(path[:-1], path[1:], strict=True)))
                if a != b
            ),
            None,
        )
        span = min(
            12 * (len(events) - 1),
            widths[branch] / 2 - 3,
            math.dist(*tail) * 0.15 if tail and tail[0][1] == tail[1][1] else 0,
        )
        for i, (_, rid, eid) in enumerate(sorted(events)):
            result[rid, eid] = (
                1 - span * (1 - i / max(1, len(events) - 1)) / total if total else 1
            )
    return result


def _joined_parts(parts, start, end, speciation):
    if speciation and len(parts) > 1 and all(p[1] == parts[0][0][1] for p in parts[0]):
        parts = parts[1:]
        first = parts[0]
        if len(first) > 1 and first[0][0] == first[1][0]:
            first = list(first)
            first[0], first[1] = start, (start[0], first[1][1])
            parts[0] = first
    points = [p for part in parts for p in part]
    points[0], points[-1] = start, end
    # Mitered rectilinear hosts retain right angles at lane changes as well.
    if all(
        a[0] == b[0] or a[1] == b[1]
        for part in parts
        for a, b in zip(part[:-1], part[1:], strict=True)
    ):
        out = [points[0]]
        for point in points[1:]:
            if out[-1][0] != point[0] and out[-1][1] != point[1]:
                out.append((out[-1][0], point[1]))
            out.append(point)
        points = out
    return points


def _bottom_primitives(
    scene, rec, lower, middle, paths, grouped, positions, fractions, lane
):
    rid = rec["id"]
    for nid in lower.order:
        events = grouped[nid]
        links = list(zip(events[:-1], events[1:], strict=True))
        if nid in lower.parents:
            links.insert(0, (grouped[lower.parents[nid]][-1], events[0]))
        for before, after in links:
            start, end = positions[rid, before["id"]], positions[rid, after["id"]]
            transfer = (
                before["type"] in {"transfer", "transfer_out"}
                and before.get("recipient_lower_node") == nid
            )
            if transfer:
                if before["recipient_node"] != after["upper_node"]:
                    raise ValueError(
                        "Nested transfer entry does not match its recipient."
                    )
                scene.primitives.append(
                    Primitive(
                        "transfer",
                        f"nested-transfer:{rid}:{before['id']}",
                        [start, end],
                        color=EVENT_STYLE["transfer"][0],
                        owner=(rid, nid),
                        tree=lower.id,
                    )
                )
                continue
            a, b = before["upper_node"], after["upper_node"]
            branches = _descendant_path(middle, a, b)
            parts = []
            for j, branch in enumerate(branches):
                first = fractions[rid, before["id"]] if j == 0 else 0
                last = fractions[rid, after["id"]] if j == len(branches) - 1 else 1
                part = _segment(paths[branch], first, last)
                left = lane[rid, a, before["lower_node"]]
                right = lane[rid, b, nid]
                offset = left + (right - left) * j / max(1, len(branches) - 1)
                parts.append(_offset(part, offset))
            points = _joined_parts(parts, start, end, before["type"] == "speciation")
            scene.primitives.append(
                Primitive(
                    "line",
                    f"nested-line:{rid}:{before['id']}:{after['id']}",
                    points,
                    color="#176b73",
                    width=1.4,
                    owner=(rid, nid),
                    tree=lower.id,
                )
            )
        for event in events:
            point = positions[rid, event["id"]]
            color, marker = EVENT_STYLE[event["type"]]
            scene.primitives.append(
                Primitive(
                    "event",
                    f"event:{rid}:{event['id']}",
                    [point],
                    text=event["type"],
                    color=color,
                    marker=marker,
                    owner=(rid, nid),
                    tree=lower.id,
                )
            )
            scene.event_points[rid, event["id"]] = point
        if not lower.nodes[nid]["children"] and lower.nodes[nid]["label"]:
            point = positions[rid, events[-1]["id"]]
            scene.primitives.append(
                Primitive(
                    "text",
                    f"nested-label:{rid}:{nid}",
                    [(point[0] + 5, point[1])],
                    text=lower.nodes[nid]["label"],
                    owner=(rid, nid),
                    tree=lower.id,
                )
            )


def build_nested_scene(document, *, width, font_size, view=None):
    chain = nested_chain(document)
    if not chain:
        raise ValueError(
            "Nested view requires one recorded middle tree with adjacent lower histories."
        )
    rec, bottom = chain
    middle = document.trees[rec["lower_tree"]]
    part = subset(document, {rec["lower_tree"], rec["upper_tree"]}, [rec])
    scene = build_scene(part, width=width, font_size=font_size, view=view)
    spacing = max(6, font_size * 1.5)
    lanes, widths = corridor_lanes(bottom, document.trees, spacing)
    max_width = max(18, max(widths.values(), default=8) + 4)
    scale = max(1, max_width / max(12, font_size * 1.5))
    for primitive in scene.primitives:
        primitive.points = [(x, y * scale) for x, y in primitive.points]
    scene.event_points = {
        key: (x, y * scale) for key, (x, y) in scene.event_points.items()
    }
    scene.height *= scale
    paths = _middle_paths(scene, rec, middle)
    by_node, positions, fractions = _bottom_geometry(
        bottom, document, paths, lanes, widths
    )
    corridors = []
    for primitive in scene.primitives:
        if primitive.kind == "upper":
            primitive.width += max_width + 12
        if (
            primitive.owner
            and primitive.owner[0] == rec["id"]
            and primitive.kind in {"line", "transfer"}
            and not primitive.annotation
        ):
            corridor = copy.deepcopy(primitive)
            corridor.points = _curve(primitive)
            corridor.kind = "upper"
            corridor.color = "#efd9b1"
            corridor.width = max(18, widths.get(primitive.owner[1], 8) + 4)
            corridors.append(corridor)
    scene.primitives.extend(corridors)
    # Root incoming middle corridor can have no serialized edge primitive.
    if not any(
        p.owner == (rec["id"], middle.root) and p.kind == "upper" for p in corridors
    ):
        scene.primitives.append(
            Primitive(
                "upper",
                f"middle-stem:{middle.id}",
                paths[middle.root],
                color="#efd9b1",
                width=max_width,
                owner=(rec["id"], middle.root),
                tree=middle.id,
            )
        )
    _separate_layer_labels(scene, document, rec, middle, bottom, widths, font_size)
    for r in bottom:
        _bottom_primitives(
            scene,
            r,
            document.trees[r["lower_tree"]],
            middle,
            paths,
            by_node[r["id"]],
            positions,
            fractions,
            lanes,
        )
        _nested_side_stubs(scene, r, middle, paths, positions, lanes)
        scene.notices.append(model_label(r))
    low = min(y - p.width / 2 for p in scene.primitives for _, y in p.points)
    shift = max(0, 15 - low)
    for primitive in scene.primitives:
        primitive.points = [(x, y + shift) for x, y in primitive.points]
    scene.event_points = {
        key: (x, y + shift) for key, (x, y) in scene.event_points.items()
    }
    scene.height = (
        max(y + p.width / 2 for p in scene.primitives for _, y in p.points) + 20
    )
    scene.notices.append(
        "Three-layer nested topology: outer tree in gray, middle incoming branches in tan, lower histories in teal. Middle markers use short leaders at corridor edges; leaf labels use separate layer columns. Corridor distance and lane positions are display coordinates, with no physical time interpretation."
    )
    scene.notices.append(source_label(document.record))
    return scene


def _separate_layer_labels(scene, document, rec, middle, bottom, widths, font_size):
    lower_label_width = max(
        (
            min(24, len(node["label"])) * font_size + 4
            for r in bottom
            for node in document.trees[r["lower_tree"]].nodes.values()
            if not node["children"]
        ),
        default=40,
    )
    middle_label_width = max(
        (
            min(24, len(node["label"])) * font_size + 4
            for node in middle.nodes.values()
            if not node["children"]
        ),
        default=40,
    )
    middle_column = lower_label_width + 18
    outer_column = middle_column + middle_label_width + 18
    leaders = []
    for primitive in scene.primitives:
        if (
            primitive.owner
            and primitive.tree == middle.id
            and primitive.owner[0] == rec["id"]
            and not primitive.annotation
        ):
            if primitive.kind == "text":
                primitive.points = [(x + middle_column, y) for x, y in primitive.points]
            elif primitive.kind == "event":
                x, y = primitive.points[0]
                offset = max(18, widths.get(primitive.owner[1], 8) + 4) / 2 + 6
                primitive.points = [(x, y - offset)]
                leaders.append(
                    Primitive(
                        "line",
                        f"middle-leader:{primitive.id}",
                        [(x, y), (x, y - offset)],
                        color="#a16207",
                        width=0.6,
                        owner=primitive.owner,
                        tree=primitive.tree,
                    )
                )
                event_id = primitive.id.removeprefix("event:" + rec["id"] + ":")
                scene.event_points[rec["id"], event_id] = primitive.points[0]
        elif primitive.kind == "text" and primitive.id.startswith("upper-label:"):
            primitive.points = [(x + outer_column, y) for x, y in primitive.points]
    scene.primitives.extend(leaders)
    scene.width = max(
        scene.width, max(x for p in scene.primitives for x, _ in p.points) + 120
    )


def _nested_side_stubs(scene, rec, middle, paths, positions, lanes):
    for event, branch, kind, name, note in terminal_side_annotations(rec, middle):
        start = positions[rec["id"], event["id"]]
        part = _offset(
            _segment(paths[branch], 0, 0.5),
            lanes[rec["id"], branch, event["id"], "empty"],
        )
        points = _joined_parts([[start, start], part], start, part[-1], True)
        color, marker = EVENT_STYLE[kind]
        owner = (rec["id"], event["lower_node"])
        scene.primitives.extend(
            [
                Primitive(
                    "line",
                    f"{name}:stub",
                    points,
                    color="#176b73" if kind == "loss" else "#6b7280",
                    owner=owner,
                    tree=rec["lower_tree"],
                    annotation=note,
                ),
                Primitive(
                    "event",
                    name,
                    [points[-1]],
                    text=kind,
                    color=color,
                    marker=marker,
                    owner=owner,
                    tree=rec["lower_tree"],
                    annotation=note,
                ),
            ]
        )
    add_terminal_side_notice(scene, rec, middle)
