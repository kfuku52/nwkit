"""Verified recorded-history context for schematic observed-G embedding."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from typing import Any

from nwkit.reconciliation_document import TreeRecord

LOSS_PROJECTION = "prune_recorded_losses_and_suppress_unary"


def _digest(label=None, children=()):
    value = ["leaf", label] if label is not None else ["fork", sorted(children)]
    return hashlib.sha256(json.dumps(value, ensure_ascii=True).encode()).hexdigest()


def projected_nodes(history, source, events, projection):
    """Prove the entire loss-pruned topology and derive, rather than trust, its map."""
    losses = {e["lower_node"] for e in events if e["type"] == "loss"}
    if any(
        history.nodes[n]["children"]
        or history.nodes[n].get("annotations", {}).get("observed", True)
        for n in losses
    ):
        raise ValueError(
            "Loss projection requires explicit unobserved terminal loss nodes."
        )
    shapes: dict[str, str | None] = {}
    retained: dict[str, str] = {}
    for nid in reversed(history.order):
        node = history.nodes[nid]
        children = [shapes[c] for c in node["children"] if shapes[c] is not None]
        if nid in losses:
            shapes[nid] = None
        elif not node["children"]:
            if not node.get("annotations", {}).get("observed", True):
                raise ValueError(
                    "Unobserved history tips require recorded loss events."
                )
            key = _digest(node["label"])
            shapes[nid] = key
            retained[key] = nid
        elif len(children) == 1:
            shapes[nid] = children[0]
        elif not children:
            shapes[nid] = None
        else:
            key = _digest(children=children)
            shapes[nid] = key
            retained[key] = nid
    target_shapes: dict[str, str] = {}
    mapping: dict[str, str] = {}
    for nid in reversed(source.order):
        node = source.nodes[nid]
        key = (
            _digest(children=[target_shapes[c] for c in node["children"]])
            if node["children"]
            else _digest(node["label"])
        )
        target_shapes[nid] = key
        if key not in retained:
            raise ValueError(
                "Loss projection does not reconstruct the whole X1 source tree."
            )
        mapping[nid] = retained[key]
    if target_shapes[source.root] != shapes[history.root]:
        raise ValueError(
            "Loss projection does not reconstruct the whole X1 source tree."
        )
    recorded = projection.get("annotations", {}).get("retained_nodes")
    if recorded is not None and recorded != mapping:
        raise ValueError(
            "Loss projection retained-node map differs from the verified topology."
        )
    return mapping


@dataclass
class ObservedContext:
    rec: dict[str, Any]
    transform: dict[str, Any]
    upper: TreeRecord
    history: TreeRecord
    source: TreeRecord
    target: TreeRecord
    mapping: dict[str, str]
    source_to_history: dict[str, str]


def observed_context(document):
    recs, transforms = (
        document.record["reconciliations"],
        document.record["transformations"],
    )
    x1 = [t for t in transforms if t["model_id"] == "X1"]
    if len(recs) != 1 or len(x1) != 1:
        raise ValueError(
            "Observed embedding requires one recorded history and one X1 transformation."
        )
    rec, transform = recs[0], x1[0]
    if rec["history_kind"] == "expectation":
        raise ValueError(
            "Observed embedding does not turn expectations into gene events."
        )
    if any(e["type"] in {"empty_observation", "unresolved"} for e in rec["events"]):
        raise ValueError(
            "Observed embedding requires resolved recorded histories; use --view embedded or paired."
        )
    upper, history, source, target = [
        document.trees[tid]
        for tid in (
            rec["upper_tree"],
            rec["lower_tree"],
            transform["source_tree"],
            transform["target_tree"],
        )
    ]
    if len({upper.id, source.id, target.id}) != 3:
        raise ValueError("Observed embedding requires distinct S, H and G trees.")
    if history.id == source.id:
        if len(transforms) != 1:
            raise ValueError(
                "Observed embedding supports one X1 transformation and an optional verified loss projection."
            )
        if any(e["type"] == "loss" for e in rec["events"]):
            raise ValueError(
                "Recorded losses require an explicit terminal-loss projection before X1."
            )
        mapping = {nid: nid for nid in source.nodes}
    else:
        projections = [
            t
            for t in transforms
            if t["source_tree"] == history.id
            and t["target_tree"] == source.id
            and t["model_id"] == LOSS_PROJECTION
            and not t["blocks"]
        ]
        if len(transforms) != 2 or len(projections) != 1:
            raise ValueError(
                "Observed embedding requires a verified loss projection connecting the recorded history to X1."
            )
        mapping = projected_nodes(history, source, rec["events"], projections[0])
    if set(document.trees) != {upper.id, history.id, source.id, target.id}:
        raise ValueError(
            "Observed embedding requires exactly S, the recorded history, X1 source H and target G."
        )
    placements = {}
    for event in rec["events"]:
        if event["type"] != "leaf":
            continue
        nid, branch = event["lower_node"], event["upper_node"]
        label = history.nodes[nid]["label"]
        if (
            history.nodes[nid]["children"]
            or label in placements
            or branch not in upper.nodes
            or upper.nodes[branch]["children"]
        ):
            raise ValueError(
                "Observed embedding requires one explicit species-tip placement per observed history leaf."
            )
        placements[label] = branch
    leaves = {n["label"] for n in target.nodes.values() if not n["children"]}
    if set(placements) != leaves:
        raise ValueError(
            "Observed embedding lacks a recorded species-tip placement for a G leaf."
        )
    return ObservedContext(
        rec, transform, upper, history, source, target, placements, mapping
    )


def history_x1_edges(context):
    """Trace the consumed H internal edges through suppressed history nodes."""
    result = set()
    for block in context.transform["blocks"]:
        parent, child = [context.source_to_history[n] for n in block["source_nodes"]]
        while child != parent:
            result.add(child)
            if child not in context.history.parents:
                raise ValueError("X1 history correspondence is not ancestral.")
            child = context.history.parents[child]
    return result
