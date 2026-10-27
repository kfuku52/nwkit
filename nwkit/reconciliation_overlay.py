"""Compatible source event forks on one G; X-replaced facts stay in details."""

from __future__ import annotations

from dataclasses import dataclass

from nwkit.reconciliation_corridors import corridor_lanes
from nwkit.reconciliation_lineages import LineagePlan, lineage_plan
from nwkit.reconciliation_projection import ObservedContext
from nwkit.reconciliation_scene import EVENT_STYLE, Primitive, Scene


@dataclass
class HistoryLayout:
    context: ObservedContext
    lanes: dict
    widths: dict
    tips: dict
    plan: LineagePlan


def history_layout(context, spacing):
    rec, upper, history = context.rec, context.upper, context.history
    lanes, widths = corridor_lanes(
        [rec], {upper.id: upper, history.id: history}, spacing
    )
    return HistoryLayout(context, lanes, widths, {}, lineage_plan(context))


def _origin_point(event, layout, xy, unit, widths, paths):
    upper, plan = layout.context.upper, layout.plan
    branch = event["upper_node"]
    right, center = xy[branch]
    left = xy[upper.parents[branch]][0] if branch in upper.parents else right - unit
    candidates = []
    line = paths[plan.tree.root]
    for a, b in zip(line[:-1], line[1:], strict=True):
        lo, hi = max(left, min(a[0], b[0])), min(right, max(a[0], b[0]))
        if a[1] == b[1] and hi > lo and abs(a[1] - center) < widths[branch] / 2:
            candidates.append((hi - lo, ((lo + hi) / 2, a[1])))
    if not candidates:
        raise ValueError(
            "Recorded origin has no incoming G route in its species branch."
        )
    return max(candidates)[1]


def recorded_overlay(layout, xy, unit, widths, spacing, root_path, positions, paths):
    """Draw actual source forks compatible with G, without any H comparison."""
    context, plan = layout.context, layout.plan
    scene = Scene(0, 0)
    events = {e["id"]: e for e in context.rec["events"]}
    points = {eid: positions[node] for eid, node in plan.main.events.items()}
    points.update(
        (e["id"], _origin_point(e, layout, xy, unit, widths, paths))
        for e in plan.origins
    )
    for eid, point in points.items():
        event = events[eid]
        color, marker = EVENT_STYLE[event["type"]]
        scene.event_points[context.rec["id"], eid] = point
        scene.primitives.append(
            Primitive(
                "event",
                f"event:{context.rec['id']}:{eid}",
                [point],
                text=event["type"],
                color=color,
                marker=marker,
                owner=(context.rec["id"], event["lower_node"]),
                tree=context.history.id,
                annotation="Recorded source event compatible with G. X-replaced source events remain in the X details and scientific document.",
            )
        )
    for node, (event, kind, name, note) in plan.main.sides.items():
        color, marker = EVENT_STYLE[kind]
        scene.primitives.append(
            Primitive(
                "event",
                name,
                [positions[node]],
                text=kind,
                color=color,
                marker=marker,
                owner=(context.rec["id"], event["lower_node"]),
                tree=context.history.id,
                annotation=note,
            )
        )
    hidden = {eid for ids in plan.hidden.values() for eid in ids}
    if set(points) & hidden or set(points) | hidden != set(events):
        raise ValueError(
            "Single-G display did not account for every recorded source event."
        )
    return scene
