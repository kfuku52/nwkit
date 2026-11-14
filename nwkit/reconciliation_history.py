"""Preserve source event forks around an X1 edit without inventing G events."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field

from nwkit.reconciliation_corridors import terminal_side_annotations
from nwkit.reconciliation_document import TreeRecord
from nwkit.reconciliation_projection import _digest


@dataclass
class HistoryPart:
    tree: TreeRecord
    placed: dict[str, str]
    sources: dict[str, str]
    events: dict[str, str]
    jumps: dict = field(default_factory=dict)
    births: dict[str, str] = field(default_factory=dict)
    sides: dict = field(default_factory=dict)
    entries: dict[str, str] = field(default_factory=dict)


def tree_record(record, nodes, root):
    order, parents, pending = [], {}, [root]
    while pending:
        n = pending.pop()
        order.append(n)
        for child in nodes[n]["children"]:
            parents[child] = n
        pending.extend(reversed(nodes[n]["children"]))
    return TreeRecord(
        dict(record, root=root), {n: nodes[n] for n in order}, parents, order
    )


def _boundary_roots(context, block):
    """Keep chains outside the consumed internal edge on their boundary lineage."""
    source, history = context.source_to_history, context.history
    pair = {source[n] for n in block["source_nodes"]}
    inner = source[block["source_nodes"][1]]
    internal = set(pair)
    while inner != source[block["source_nodes"][0]]:
        internal.add(inner)
        inner = history.parents[inner]
    roots = {}
    for i, boundary in enumerate(block["boundaries"], 1):
        node = source[boundary["source_node"]]
        while history.parents[node] not in internal:
            node = history.parents[node]
        roots[node] = boundary["target_node"], f"{block['id']}.{i}"
    return internal, roots


def _source_event_ids(context, block, roots):
    """Keep X-replaced source facts in details without drawing another tree."""
    history = context.history
    root = context.source_to_history[block["source_nodes"][0]]
    included = set()
    pending = [root]
    while pending:
        n = pending.pop()
        if n in roots:
            continue
        included.add(n)
        pending.extend(history.nodes[n]["children"])
    rank = {n: i for i, n in enumerate(history.order)}
    events = [
        e
        for e in context.rec["events"]
        if e["lower_node"] in included and e["type"] != "origin"
    ]
    return tuple(e["id"] for e in sorted(events, key=lambda e: rank[e["lower_node"]]))


def expanded_history(context, retained):
    """Apply X only to topology; retain compatible full-history chains and losses."""
    history = context.history
    by_node = defaultdict(list)
    for event in context.rec["events"]:
        by_node[event["lower_node"]].append(event)
    identifiers = {n: retained.get(n, f"display-history:{n}") for n in history.nodes}
    if len(set(identifiers.values())) != len(identifiers):
        raise ValueError("History display-node IDs collide with recorded G nodes.")
    nodes = {
        identifiers[n]: dict(
            node, id=identifiers[n], children=[identifiers[c] for c in node["children"]]
        )
        for n, node in history.nodes.items()
    }
    hidden, consumed = {}, set()
    for block in context.transform["blocks"]:
        internal, roots = _boundary_roots(context, block)
        consumed.update(internal)
        hidden[block["id"]] = _source_event_ids(context, block, roots)
        boundaries = {g: identifiers[n] for n, (g, _) in roots.items()}
        pair = set(block["target_nodes"])
        for n in pair:
            nodes[n]["children"] = [
                c if c in pair else boundaries[c]
                for c in context.target.nodes[n]["children"]
            ]
    tree = tree_record(context.target.record, nodes, identifiers[history.root])
    sources = {key: n for n, key in identifiers.items() if key in tree.nodes}
    placed = {key: by_node[n][-1]["upper_node"] for key, n in sources.items()}
    events = {
        e["id"]: identifiers[e["lower_node"]]
        for e in context.rec["events"]
        if e["lower_node"] not in consumed
        and identifiers[e["lower_node"]] in tree.nodes
        and e["type"] != "origin"
    }
    main = HistoryPart(tree, placed, sources, events)
    return main, hidden


def add_side_ports(part, context):
    """A recorded marginalized side is a fork, without inventing a sampled loss."""
    nodes = {
        n: dict(v, children=list(v["children"])) for n, v in part.tree.nodes.items()
    }
    annotations = list(terminal_side_annotations(context.rec, context.upper))
    known = {e["id"] for e, *_ in annotations}
    for event in context.rec["events"]:
        recorded = event.get("annotations", {})
        branch = recorded.get("marginalized_species")
        if (
            event["id"] not in known
            and event["type"] in {"speciation", "transfer", "transfer_out"}
            and recorded.get("observed_arity") == 1
            and branch in context.upper.nodes
        ):
            name = f"empty-side:{context.rec['id']}:{event['id']}:{branch}"
            note = "Source-recorded unobserved side, integrated in E. The stub is display geometry, not an explicit loss, a loss count or an event time."
            annotations.append((event, branch, "empty_observation", name, note))
    for event, branch, kind, name, note in annotations:
        parent = part.events.get(event["id"])
        if parent is None:
            continue
        key = f"display-side:{name}"
        nodes[key] = dict(id=key, label="", children=[], branch_length=None)
        nodes[parent]["children"].append(key)
        part.placed[key] = branch
        part.sides[key] = event, kind, name, note
    part.tree = tree_record(part.tree.record, nodes, part.tree.root)


def add_transfers(part, context, fallback):
    """Each source transfer has its own fresh recipient entry, including dead sides."""
    nodes = {
        n: dict(v, children=list(v["children"])) for n, v in part.tree.nodes.items()
    }
    source_ids = part.entries
    for event in context.rec["events"]:
        if event["type"] not in {"transfer", "transfer_out"}:
            continue
        donor = part.events.get(event["id"])
        recipient = source_ids.get(event["recipient_lower_node"])
        if donor is not None and event["recipient_lower_node"] == event["lower_node"]:
            recipient = next(c for c in nodes[donor]["children"] if c not in part.sides)
        if donor is None:
            child = fallback.get(event["id"])
            if child is None:
                continue
            recipient = source_ids.get(event["recipient_lower_node"], child)
            donor = _transfer_donor(part, nodes, recipient, event)
        if donor is None or recipient is None:
            continue
        entry = f"display-transfer:{event['id']}:entry"
        if entry in nodes:
            raise ValueError("Transfer display-port ID collides with a recorded node.")
        nodes[entry] = dict(
            id=entry, label="", children=[recipient], branch_length=None
        )
        parent = next(
            n for n in nodes if n != entry and recipient in nodes[n]["children"]
        )
        nodes[parent]["children"] = [
            entry if c == recipient else c for c in nodes[parent]["children"]
        ]
        part.placed[entry] = event["recipient_node"]
        part.jumps[entry], part.births[event["id"]] = event, parent
        for e in context.rec["events"]:
            if (
                e["type"] == "transfer_in"
                and e["lower_node"] == event["recipient_lower_node"]
                and e["id"] in part.events
            ):
                part.events[e["id"]] = entry
    part.tree = tree_record(part.tree.record, nodes, part.tree.root)


def _transfer_donor(part, nodes, recipient, event):
    """Keep an incoming source lineage separate even beside a same-branch fork."""
    parent = next((n for n in nodes if recipient in nodes[n]["children"]), None)
    donor = f"display-transfer:{event['id']}:donor"
    if donor in nodes:
        raise ValueError("Transfer display-port ID collides with a recorded node.")
    nodes[donor] = dict(id=donor, label="", children=[recipient], branch_length=None)
    if parent is not None:
        nodes[parent]["children"] = [
            donor if c == recipient else c for c in nodes[parent]["children"]
        ]
    else:
        part.tree.record = dict(part.tree.record, root=donor)
    part.placed[donor] = event["donor_node"]
    return donor


def verify_observed_projection(part, target):
    """Prove that pruning unobserved sides and unary display ports recovers G."""
    shapes: dict[str, str | None] = {}
    for n in reversed(part.tree.order):
        children = [
            shapes[c] for c in part.tree.nodes[n]["children"] if shapes[c] is not None
        ]
        if len(children) == 1:
            shapes[n] = children[0]
        elif children:
            shapes[n] = _digest(children=children)
        elif n in target.nodes and not target.nodes[n]["children"]:
            shapes[n] = _digest(target.nodes[n]["label"])
        else:
            shapes[n] = None
    expected: dict[str, str] = {}
    for n in reversed(target.order):
        node = target.nodes[n]
        expected[n] = (
            _digest(children=[expected[c] for c in node["children"]])
            if node["children"]
            else _digest(node["label"])
        )
        if shapes.get(n) != expected[n]:
            raise ValueError(
                "Event-expanded display does not reconstruct exact G topology."
            )
    if shapes[part.tree.root] != expected[target.root]:
        raise ValueError("Event-expanded display has a different observed root.")
