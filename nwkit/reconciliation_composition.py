"""Compose neutral adjacent-level and H/G views; no inference or time placement."""

from __future__ import annotations

import copy
from collections import Counter

from nwkit.reconciliation_document import validate_document
from nwkit.reconciliation_labels import runtime_label, source_label
from nwkit.reconciliation_scene import Primitive, Scene, _standalone_scene, build_scene
from nwkit.reconciliation_view import view_for_document


def append_scene(target, source, label):
    offset = target.height + 32
    target.primitives.append(
        Primitive(
            "text",
            f"panel:{len(target.primitives)}",
            [(20, target.height + 14)],
            text=label,
        )
    )
    for p in source.primitives:
        p = copy.deepcopy(p)
        p.points = [(x, y + offset) for x, y in p.points]
        target.primitives.append(p)
    target.event_points.update(
        {key: (x, y + offset) for key, (x, y) in source.event_points.items()}
    )
    target.height += source.height + 60
    target.notices.extend(source.notices)


def subset(document, trees, recs=()):
    record = copy.deepcopy(document.record)
    record["trees"] = [t for t in record["trees"] if t["id"] in trees]
    record["reconciliations"] = list(recs)
    record["transformations"] = []
    record["provenance"].pop("display_selection", None)
    return validate_document(record)


def compose_document(document, *, width, font_size, view=None):
    from nwkit.reconciliation_nested import build_nested_scene, nested_chain

    if (view or {}).get("layout") != "paired" and nested_chain(document):
        return build_nested_scene(document, width=width, font_size=font_size, view=view)
    scene = Scene(width, 0)
    displayed = set()
    recs = document.record["reconciliations"]
    if (view or {}).get("layout") == "paired":
        append_scene(
            scene,
            _standalone_scene(document, width, max(12, font_size * 1.5), view),
            "Recorded trees",
        )
        displayed.update(document.trees)
    else:
        groups: dict[str, list[dict]] = {}
        for rec in recs:
            groups.setdefault(rec["upper_tree"], []).append(rec)
        lower_ids = {r["lower_tree"] for r in recs}
        # Root level first, then stable ID order for adjacent levels.
        for upper, group in sorted(
            groups.items(), key=lambda x: (x[0] in lower_ids, x[0])
        ):
            tree_ids = {upper} | {r["lower_tree"] for r in group}
            part = subset(document, tree_ids, group)
            append_scene(
                scene,
                build_scene(part, width=width, font_size=font_size, view=view),
                f"Level: {', '.join(sorted(tree_ids - {upper}))} → {upper}",
            )
            displayed.update(tree_ids)
    for tid in sorted(set(document.trees) - displayed):
        part = subset(document, {tid})
        append_scene(
            scene,
            build_scene(part, width=width, font_size=font_size, view=view),
            f"Recorded tree: {tid}",
        )
    for transform in document.record["transformations"]:
        scene.notices.append(
            f"Transformation {transform['id']}: {transform['source_tree']} → {transform['target_tree']}; model={transform['model_id']}"
        )
        for block in transform["blocks"]:
            scene.notices.append(
                f"{block['id']}: two adjacent H nodes map to two adjacent {transform['target_tree']} nodes; resolution={block['resolution']}; three boundary clades retained in source JSON. Purple branches mark the consumed nodes. X1 has no physical event time or ILS interpretation."
            )
            for side, field in (
                ("source_tree", "source_nodes"),
                ("target_tree", "target_nodes"),
            ):
                tid = transform[side]
                owners = {tid} | {r["id"] for r in recs if r["lower_tree"] == tid}
                for p in scene.primitives:
                    if p.owner and p.owner[0] in owners and p.owner[1] in block[field]:
                        p.color = "#a03faa"
                        p.width = max(p.width, 2.4)
    if len({r["upper_tree"] for r in recs} & {r["lower_tree"] for r in recs}):
        scene.notices.append(
            "Multi-layer overview uses linked adjacent-level panels. Repeated middle-tree IDs identify the same scientific tree; coordinates are topological in each panel."
        )
    scene.notices.append(source_label(document.record))
    return scene


def apply_view(scene, document, view):
    if not view:
        return scene
    colors = view.get("event_colors", {})
    for p in scene.primitives:
        if p.kind == "event" and p.text in colors:
            p.color = colors[p.text]
    for item in view.get("collapsed_clades", []):
        tree = document.trees[item["tree"]]
        descendants = {
            nid
            for nid in tree.nodes
            if nid != item["node"] and tree.ancestral(item["node"], nid)
        }
        owners = {tree.id} | {
            r["id"]
            for r in document.record["reconciliations"]
            if r["lower_tree"] == tree.id
        }
        removed = [
            p
            for p in scene.primitives
            if p.owner and p.owner[0] in owners and p.owner[1] in descendants
        ]
        scene.primitives = [p for p in scene.primitives if p not in removed]
        anchors = [
            p.points[-1]
            for p in scene.primitives
            if p.owner and p.owner[0] in owners and p.owner[1] == item["node"]
        ]
        if anchors:
            scene.primitives.append(
                Primitive(
                    "text",
                    f"collapse:{tree.id}:{item['node']}",
                    [anchors[-1]],
                    text=f"collapsed ({len(descendants)} descendants)",
                    color="#6b7280",
                )
            )
        scene.notices.append(
            f"Display-only collapse {tree.id}/{item['node']}; original events and summaries retained."
        )
    return scene


def _event_key(document, rec, event):
    def leaves(tid, nid):
        tree = document.trees[tid]
        return tuple(
            sorted(
                n["label"]
                for key, n in tree.nodes.items()
                if not n["children"] and tree.ancestral(nid, key)
            )
        )

    return (
        event["type"],
        leaves(rec["lower_tree"], event["lower_node"]),
        leaves(rec["upper_tree"], event["upper_node"]) if event["upper_node"] else None,
        leaves(rec["upper_tree"], event["recipient_node"])
        if event.get("recipient_node")
        else None,
        leaves(rec["lower_tree"], event["recipient_lower_node"])
        if event.get("recipient_lower_node")
        else None,
    )


def comparison_scene(record, documents, *, width, font_size):
    scene = Scene(width, 0)
    reference = None
    for panel in record["panels"]:
        label = f"{panel['label']} — inference: {panel['inference_status']}; history: {panel['history_status']}"
        document = documents.get(panel["id"])
        if document and panel["history_status"] in {"available", "partial"}:
            part = build_scene(
                document, width=width, font_size=font_size, view=record["view"]
            )
            apply_view(part, document, view_for_document(record["view"], document))
            recs = document.record["reconciliations"]
            conditions = [(r["history_kind"], r["selection"], r["model"]) for r in recs]
            events = Counter(
                _event_key(document, r, e) for r in recs for e in r["events"]
            )
            if reference and reference[0] == conditions:
                scene.notices.append(
                    f"{panel['label']}: event multiset difference from first compatible panel: +{sum((events - reference[1]).values())}/-{sum((reference[1] - events).values())}. Compare recorded events only."
                )
            elif reference:
                scene.notices.append(
                    f"{panel['label']}: event difference unavailable because model conditions/history selection differ."
                )
            else:
                reference = conditions, events
        elif document and not document.record["reconciliations"]:
            part = build_scene(
                document,
                width=width,
                font_size=font_size,
                view={"layout": "standalone"},
            )
            part.notices.append(
                f"{panel['label']}: inferred tree retained; history {panel['history_status']}, so no event placement is displayed."
            )
        else:
            part = Scene(
                width,
                55,
                [
                    Primitive(
                        "text",
                        f"status:{panel['id']}",
                        [(30, 25)],
                        text=f"History {panel['history_status']}; no event geometry available.",
                    )
                ],
            )
        append_scene(scene, part, label)
        scene.notices.append(runtime_label(panel))
    scene.notices.append(
        f"Comparison mode={record['mode']}; {len(record['panels'])} panels, {len(documents)} source documents. Missing/failed/timeout panels are retained."
    )
    return scene
