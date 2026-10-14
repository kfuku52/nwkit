"""Versioned, file-only reconciliation documents; never performs inference."""

from __future__ import annotations

import copy
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

VERSION = "0.1.0"
EVENT_TYPES = {
    "leaf",
    "speciation",
    "duplication",
    "loss",
    "transfer",
    "transfer_out",
    "transfer_in",
    "origin",
    "unresolved",
    "empty_observation",
}
HISTORY_KINDS = {"representative", "sample", "expectation", "unspecified"}


@dataclass
class TreeRecord:
    record: dict[str, Any]
    nodes: dict[str, dict[str, Any]]
    parents: dict[str, str]
    order: list[str]

    @property
    def id(self) -> str:
        return self.record["id"]

    @property
    def root(self) -> str:
        return self.record["root"]

    def ancestral(self, ancestor: str, descendant: str) -> bool:
        while descendant != ancestor and descendant in self.parents:
            descendant = self.parents[descendant]
        return descendant == ancestor


@dataclass
class ReconciliationDocument:
    record: dict[str, Any]
    trees: dict[str, TreeRecord]


def _fields(value, required, optional, label):
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be an object.")
    missing, extra = (
        set(required) - value.keys(),
        value.keys() - set(required) - set(optional),
    )
    if missing or extra:
        raise ValueError(
            f"{label}: missing {sorted(missing)}; unknown {sorted(extra)}."
        )
    if "annotations" in value and not isinstance(value["annotations"], dict):
        raise ValueError(f"{label} annotations must be an object.")
    if "role" in value:
        _string(value["role"], f"{label} role")
    if "summary" in value and not isinstance(value["summary"], dict):
        raise ValueError(f"{label} summary must be an object.")


def _string(value, label):
    if not isinstance(value, str) or not value:
        raise ValueError(f"{label} must be a nonempty string.")


def _number(value, label, *, probability=False):
    if value is None:
        return
    if (
        isinstance(value, bool)
        or not isinstance(value, (float, int))
        or not math.isfinite(value)
    ):
        raise ValueError(f"{label} must be a finite number or null.")
    if value < 0 or (probability and value > 1):
        raise ValueError(f"{label} is out of range.")


def _unique(records, label):
    if not isinstance(records, list):
        raise ValueError(f"{label} must be an array.")
    result = {}
    for item in records:
        if not isinstance(item, dict):
            raise ValueError(f"{label} member must be an object.")
        identifier = item.get("id")
        _string(identifier, f"{label} ID")
        if identifier in result:
            raise ValueError(f"Duplicate {label} ID: {identifier}.")
        result[identifier] = item
    return result


def _tree(record):
    _fields(
        record,
        {"id", "root", "nodes", "branch_length_unit"},
        {"role", "annotations"},
        "tree",
    )
    _string(record["root"], "tree root")
    _string(record["branch_length_unit"], "branch_length_unit")
    nodes = _unique(record["nodes"], "node")
    if record["root"] not in nodes:
        raise ValueError("Tree root must reference a node.")
    parents = {}
    for nid, node in nodes.items():
        _fields(
            node, {"id", "label", "children", "branch_length"}, {"annotations"}, "node"
        )
        if not isinstance(node["label"], str) or not isinstance(node["children"], list):
            raise ValueError("Node label/children must be string/array.")
        _number(node["branch_length"], "branch_length")
        for child in node["children"]:
            if not isinstance(child, str) or child not in nodes or child in parents:
                raise ValueError("Unknown child or multiple parents in tree.")
            parents[child] = nid
    order, seen, pending = [], set(), [record["root"]]
    while pending:
        nid = pending.pop()
        if nid in seen:
            raise ValueError("Tree contains a cycle.")
        seen.add(nid)
        order.append(nid)
        pending.extend(reversed(nodes[nid]["children"]))
    if len(seen) != len(nodes) or record["root"] in parents:
        raise ValueError("Tree has a cycle or disconnected nodes.")
    labels = [
        n["label"]
        for n in nodes.values()
        if not n["children"] and n.get("annotations", {}).get("observed", True)
    ]
    if any(not label for label in labels) or len(labels) != len(set(labels)):
        raise ValueError(
            "Observed leaf labels must be nonempty and unique within each tree."
        )
    return TreeRecord(record, nodes, parents, order)


def _event(event, lower, upper):
    _fields(
        event,
        {"id", "type", "lower_node", "upper_node"},
        {
            "donor_node",
            "recipient_node",
            "recipient_lower_node",
            "time",
            "support",
            "annotations",
        },
        "event",
    )
    if event["type"] not in EVENT_TYPES or event["lower_node"] not in lower.nodes:
        raise ValueError("Unknown event type or lower node.")
    for field in ("upper_node", "donor_node", "recipient_node"):
        if event.get(field) is not None and event[field] not in upper.nodes:
            raise ValueError(f"Event {field} does not belong to the upper tree.")
    if (
        event.get("recipient_lower_node") is not None
        and event["recipient_lower_node"] not in lower.nodes
    ):
        raise ValueError("Unknown recipient_lower_node.")
    if event["type"] in {"transfer", "transfer_out"} and any(
        event.get(k) is None
        for k in ("donor_node", "recipient_node", "recipient_lower_node")
    ):
        raise ValueError(
            "Transfer requires explicit donor, recipient and recipient lower node."
        )
    if event["type"] in {"transfer", "transfer_out"} and (
        event["upper_node"] != event["donor_node"]
        or event["donor_node"] == event["recipient_node"]
    ):
        raise ValueError(
            "Transfer placement must be its donor and recipient must differ."
        )
    if event["type"] == "leaf" and (
        lower.nodes[event["lower_node"]]["children"]
        or (
            event["upper_node"] is not None
            and upper.nodes[event["upper_node"]]["children"]
        )
    ):
        raise ValueError("Leaf event requires lower and upper leaves.")
    if (
        event["type"] == "speciation"
        and event["upper_node"] is not None
        and not upper.nodes[event["upper_node"]]["children"]
    ):
        raise ValueError("Speciation cannot map to an upper leaf.")
    if event.get("time") is not None:
        _fields(event["time"], {"lower", "upper", "unit"}, set(), "event time")
        time = event["time"]
        _number(time["lower"], "time lower")
        _number(time["upper"], "time upper")
        _string(time["unit"], "time unit")
        if (
            time["lower"] is None
            or time["upper"] is None
            or time["lower"] > time["upper"]
        ):
            raise ValueError("Invalid event time interval.")
    if event.get("support") is not None:
        support = event["support"]
        _fields(support, {"kind", "value", "conditional_on"}, set(), "event support")
        _string(support["kind"], "support kind")
        _number(support["value"], "support value", probability=True)
        if support["value"] is not None:
            _string(support["conditional_on"], "support conditional_on")


def _reconciliation(rec, trees):
    _fields(
        rec,
        {
            "id",
            "lower_tree",
            "upper_tree",
            "history_kind",
            "events",
            "model",
            "selection",
        },
        {"summary", "annotations"},
        "reconciliation",
    )
    if rec["lower_tree"] == rec["upper_tree"] or any(
        rec[k] not in trees for k in ("lower_tree", "upper_tree")
    ):
        raise ValueError("Reconciliation must connect distinct known trees.")
    if (
        rec["history_kind"] not in HISTORY_KINDS
        or not isinstance(rec["model"], dict)
        or not isinstance(rec["selection"], dict)
    ):
        raise ValueError("Invalid history kind/model/selection.")
    if rec["history_kind"] == "representative":
        _string(rec["selection"].get("method"), "Representative selection method")
    lower, upper = trees[rec["lower_tree"]], trees[rec["upper_tree"]]
    events = _unique(rec["events"], "event")
    for event in events.values():
        _event(event, lower, upper)
    _transfer_entries(rec["events"], lower)
    # Hierarchy graph cycles are checked separately; event placement is not inferred.


def _transfer_entries(events, lower):
    by_node: dict[str, list[dict]] = {}
    for event in events:
        by_node.setdefault(event["lower_node"], []).append(event)
    paired = set()
    for event in events:
        if event["type"] not in {"transfer", "transfer_out"}:
            continue
        source, target = event["lower_node"], event["recipient_lower_node"]
        if target == source:
            sequence = by_node[source]
            index = sequence.index(event) + 1
            if index >= len(sequence):
                raise ValueError("Transfer lacks a recipient entry event.")
            entry = sequence[index]
        else:
            if lower.parents.get(target) != source or target not in by_node:
                raise ValueError(
                    "Transfer recipient entry must be the same lineage or an immediate child."
                )
            entry = by_node[target][0]
        if entry["upper_node"] != event["recipient_node"]:
            raise ValueError(
                "Transfer recipient entry placement differs from its recipient."
            )
        paired.add(entry["id"])
    if any(e["type"] == "transfer_in" and e["id"] not in paired for e in events):
        raise ValueError("Unpaired transfer_in event.")


def _x1_blocks(transform, trees):
    source, target = trees[transform["source_tree"]], trees[transform["target_tree"]]
    if any(
        len(n["children"]) not in {0, 2}
        or (not n["children"] and not n.get("annotations", {}).get("observed", True))
        for tree in (source, target)
        for n in tree.nodes.values()
    ):
        raise ValueError(
            "X1 requires binary trees of observed copies after unary suppression."
        )

    def leaves(tree, nid):
        stack, labels = [nid], []
        while stack:
            current = stack.pop()
            children = tree.nodes[current]["children"]
            if children:
                stack.extend(children)
            else:
                labels.append(tree.nodes[current]["label"])
        return sorted(labels)

    if leaves(source, source.root) != leaves(target, target.root):
        raise ValueError("X1 source/target leaf identities differ.")
    consumed: set[str] = set()
    for block in _unique(transform["blocks"], "X1 block").values():
        _fields(
            block,
            {"id", "source_nodes", "target_nodes", "boundaries", "resolution"},
            {"annotations"},
            "X1 block",
        )
        if block["resolution"] not in {1, 2} or isinstance(block["resolution"], bool):
            raise ValueError("X1 resolution must be 1 or 2.")
        for field, tree in (("source_nodes", source), ("target_nodes", target)):
            pair = block[field]
            if (
                not isinstance(pair, list)
                or len(pair) != 2
                or len(set(pair)) != 2
                or any(n not in tree.nodes for n in pair)
            ):
                raise ValueError(
                    "X1 consumed nodes must be a pair of known distinct nodes."
                )
            if tree.parents.get(pair[1]) != pair[0]:
                raise ValueError(
                    "X1 consumed nodes must be adjacent root/internal-child nodes."
                )
            if (
                len(tree.nodes[pair[0]]["children"]) != 2
                or len(tree.nodes[pair[1]]["children"]) != 2
            ):
                raise ValueError("X1 consumed nodes must be binary internal nodes.")
        if consumed.intersection(block["source_nodes"]):
            raise ValueError("X1 latent blocks overlap.")
        consumed.update(block["source_nodes"])
        if not isinstance(block["boundaries"], list) or len(block["boundaries"]) != 3:
            raise ValueError("X1 requires three boundary clades.")
        for boundary in block["boundaries"]:
            _fields(boundary, {"source_node", "target_node"}, set(), "X1 boundary")
            if (
                boundary["source_node"] not in source.nodes
                or boundary["target_node"] not in target.nodes
            ):
                raise ValueError("Unknown X1 boundary node.")
            if leaves(source, boundary["source_node"]) != leaves(
                target, boundary["target_node"]
            ):
                raise ValueError("X1 boundary leaf identities differ.")
        _validate_x1_resolution(block, source, target, leaves)


def _validate_x1_resolution(block, source, target, leaves):
    def actual(tree, pair):
        root, child = pair
        return set(tree.nodes[child]["children"]) | (
            set(tree.nodes[root]["children"]) - {child}
        )

    boundaries = block["boundaries"]
    for field, tree, pair in (
        ("source_node", source, block["source_nodes"]),
        ("target_node", target, block["target_nodes"]),
    ):
        if {b[field] for b in boundaries} != actual(tree, pair):
            raise ValueError("X1 boundaries must be the three actual boundary clades.")
    a, b, c = [boundary["source_node"] for boundary in boundaries]
    if {a, b} != set(source.nodes[block["source_nodes"][1]]["children"]):
        raise ValueError("X1 boundaries must be ordered A,B,C for ((A,B),C).")
    chosen = a if block["resolution"] == 1 else b
    expected = sorted(leaves(source, chosen) + leaves(source, c))
    if leaves(target, block["target_nodes"][1]) != expected:
        raise ValueError(
            "X1 target local split does not match the recorded resolution."
        )


def validate_document(record) -> ReconciliationDocument:
    """Validate structure and identity without changing model or missing values."""
    record = copy.deepcopy(record)
    _fields(
        record,
        {
            "schema_version",
            "kind",
            "trees",
            "reconciliations",
            "transformations",
            "provenance",
        },
        {"title", "annotations"},
        "document",
    )
    if (
        record["schema_version"] != VERSION
        or record["kind"] != "reconciliation_document"
    ):
        raise ValueError("Unsupported reconciliation document kind/schema_version.")
    # Also rejects nonfinite values inside arbitrary model/extension metadata.
    try:
        json.dumps(record, allow_nan=False)
    except (ValueError, TypeError) as exc:
        raise ValueError("Document must contain finite JSON values.") from exc
    if not isinstance(record["provenance"], dict):
        raise ValueError("Document provenance must be an object.")
    if "title" in record and not isinstance(record["title"], str):
        raise ValueError("Document title must be a string.")
    trees = {
        key: _tree(value) for key, value in _unique(record["trees"], "tree").items()
    }
    if not trees:
        raise ValueError("Document must contain at least one tree.")
    recs = _unique(record["reconciliations"], "reconciliation")
    for rec in recs.values():
        _reconciliation(rec, trees)
    _validate_display_selection(record, trees, recs)
    edges: dict[str, set[str]] = {}
    for rec in recs.values():
        edges.setdefault(rec["lower_tree"], set()).add(rec["upper_tree"])
    for start in edges:
        pending, seen = list(edges[start]), set()
        while pending:
            current = pending.pop()
            if current == start:
                raise ValueError("Reconciliation hierarchy contains a cycle.")
            if current not in seen:
                seen.add(current)
                pending.extend(edges.get(current, ()))
    for transform in _unique(record["transformations"], "transformation").values():
        _fields(
            transform,
            {"id", "source_tree", "target_tree", "model_id", "blocks"},
            {"annotations"},
            "transformation",
        )
        if transform["source_tree"] == transform["target_tree"] or any(
            transform[k] not in trees for k in ("source_tree", "target_tree")
        ):
            raise ValueError("Transformation requires distinct known trees.")
        _string(transform["model_id"], "transformation model_id")
        if not isinstance(transform["blocks"], list) or any(
            not isinstance(b, dict) for b in transform["blocks"]
        ):
            raise ValueError("Transformation blocks must be an array.")
        if transform["model_id"] == "X1":
            _x1_blocks(transform, trees)
    return ReconciliationDocument(record, trees)


def _validate_display_selection(record, trees, recs):
    selection = record["provenance"].get("display_selection")
    if selection is None:
        return
    _fields(
        selection,
        {"tree", "node", "inference_recomputed", "summaries"},
        set(),
        "display_selection",
    )
    _string(selection["tree"], "display_selection tree")
    _string(selection["node"], "display_selection node")
    _string(selection["summaries"], "display_selection summaries")
    if (
        selection["tree"] not in trees
        or selection["node"] not in trees[selection["tree"]].nodes
        or selection["inference_recomputed"] is not False
    ):
        raise ValueError("Invalid display_selection reference or recomputation flag.")
    if any(r["upper_tree"] == selection["tree"] for r in recs.values()):
        raise ValueError("display_selection must identify a standalone or lower clade.")


def document_json(document: ReconciliationDocument) -> str:
    return (
        json.dumps(
            document.record,
            ensure_ascii=False,
            allow_nan=False,
            sort_keys=True,
            indent=2,
        )
        + "\n"
    )


def select_document_clade(document, selection):
    try:
        tid, nid = selection.split("/", 1)
        document.trees[tid].nodes[nid]
    except (ValueError, KeyError) as exc:
        raise ValueError(
            "--document-clade must identify TREE/NODE in the document."
        ) from exc
    if any(r["upper_tree"] == tid for r in document.record["reconciliations"]):
        raise ValueError(
            "Select a lower clade; upper-tree pruning would alter the reconciliation context."
        )
    record = copy.deepcopy(document.record)
    record["provenance"]["display_selection"] = {
        "tree": tid,
        "node": nid,
        "inference_recomputed": False,
        "summaries": "original full-source summaries",
    }
    record["title"] = record.get("title", "") + " (selected clade)"
    return validate_document(record)


def read_document(path) -> ReconciliationDocument:
    from nwkit.reconciliation_import import import_tree_document, import_xml_document

    text = Path(path).read_text(encoding="utf-8-sig")
    if text.lstrip().startswith("{"):
        record = load_json(text)
    elif text.lstrip().startswith("<"):
        record = import_xml_document(text)
    else:
        record = import_tree_document(text)
    return validate_document(record)


def load_json(text):
    """Strict JSON transport shared by scientific documents and display state."""

    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError(f"Duplicate JSON key: {key}.")
            result[key] = value
        return result

    def nonfinite(value):
        raise ValueError(f"Nonfinite JSON value: {value}.")

    def finite_float(value):
        number = float(value)
        if not math.isfinite(number):
            raise ValueError("Nonfinite JSON number.")
        return number

    return json.loads(
        text,
        object_pairs_hook=pairs,
        parse_constant=nonfinite,
        parse_float=finite_float,
    )
