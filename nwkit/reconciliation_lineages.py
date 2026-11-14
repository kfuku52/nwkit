"""Transfer-aware subdivisions of G, preserving its recorded branching topology."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

from nwkit.reconciliation_document import TreeRecord
from nwkit.reconciliation_event_chains import expand_event_chains
from nwkit.reconciliation_history import (
    HistoryPart,
    add_side_ports,
    add_transfers,
    expanded_history,
    verify_observed_projection,
)
from nwkit.reconciliation_scene import _signatures


@dataclass
class LineagePlan:
    tree: TreeRecord
    placed: dict[str, str]
    jumps: dict
    births: dict[str, str]
    owners: dict[str, str]
    retained: dict[str, str]
    descendants: dict[str, set[str]]
    main: HistoryPart
    hidden: dict[str, tuple[str, ...]]
    origins: list[dict]


def _correspondence(context):
    source, target = context.source, context.target
    signatures = _signatures(target)
    by_clade = {clade: n for n, clade in signatures.items()}
    mapping = {
        n: by_clade[clade]
        for n, clade in _signatures(source).items()
        if clade in by_clade
    }
    for block in context.transform["blocks"]:
        mapping.update(zip(block["source_nodes"], block["target_nodes"], strict=True))
    if set(mapping) != set(source.nodes):
        raise ValueError("X1 correspondence does not cover every source node.")
    retained = {context.source_to_history[n]: g for n, g in mapping.items()}
    descendants: dict[str, set[str]] = {}
    for n in reversed(context.history.order):
        descendants[n] = (
            {retained[n]}
            if n in retained
            else set().union(
                *(descendants[c] for c in context.history.nodes[n]["children"])
            )
        )
    return retained, descendants, signatures


def _transfers(context, descendants, signatures):
    transfers = defaultdict(list)
    history = context.history
    rank = {n: i for i, n in enumerate(history.order)}
    for event in context.rec["events"]:
        if event["type"] not in {"transfer", "transfer_out"}:
            continue
        if not event.get("recipient_lower_node") or not event.get("recipient_node"):
            raise ValueError(
                "Observed embedding requires an explicit transfer recipient lineage."
            )
        recipient = event["recipient_lower_node"]
        candidates = descendants[recipient]
        if not candidates:
            continue  # An entirely lost recipient is drawn as an explicit loss spur.
        if len(candidates) != 1:
            raise ValueError("Transfer recipient must identify one observed lineage.")
        child = next(iter(candidates))
        observed = []
        pending = [recipient]
        while pending:
            n = pending.pop()
            node = history.nodes[n]
            pending.extend(node["children"])
            if not node["children"] and node.get("annotations", {}).get(
                "observed", True
            ):
                observed.append(node["label"])
        if tuple(sorted(observed)) != signatures[child]:
            raise ValueError(
                "X1 crosses a transfer recipient boundary; a single-lineage G embedding is undefined. Use --view embedded or paired."
            )
        transfers[child].append(event)
    for events in transfers.values():
        events.sort(key=lambda e: rank[e["lower_node"]])
    return transfers


def _placements(context, retained, transfers):
    upper, target = context.upper, context.target
    depth = {upper.root: 0}
    for n in upper.order[1:]:
        depth[n] = depth[upper.parents[n]] + 1

    def lca(a, b):
        while a != b:
            if depth[a] >= depth[b]:
                a = upper.parents[a]
            else:
                b = upper.parents[b]
        return a

    consumed = {n for b in context.transform["blocks"] for n in b["target_nodes"]}
    fixed = {}
    for e in context.rec["events"]:
        if e["lower_node"] in retained and retained[e["lower_node"]] not in consumed:
            fixed[retained[e["lower_node"]]] = e["upper_node"]
    placed: dict[str, str] = {}
    for n in reversed(target.order):
        node = target.nodes[n]
        requirements = [
            transfers[c][0]["donor_node"] if c in transfers else placed[c]
            for c in node["children"]
        ]
        branch = context.mapping[node["label"]] if not requirements else requirements[0]
        for other in requirements[1:]:
            branch = lca(branch, other)
        placed[n] = fixed.get(n, branch)
        if not upper.ancestral(placed[n], branch):
            raise ValueError("Recorded placement conflicts with the X1 lineage routes.")
    return placed


def lineage_plan(context):
    """One event-expanded G; incompatible source forks stay in X details."""
    retained, descendants, signatures = _correspondence(context)
    transfers = _transfers(context, descendants, signatures)
    placed = _placements(context, retained, transfers)
    main, hidden = expanded_history(context, retained)
    main.placed.update(placed)
    fallback = {e["id"]: child for child, events in transfers.items() for e in events}
    expand_event_chains(main, context)
    add_side_ports(main, context)
    add_transfers(main, context, fallback)
    verify_observed_projection(main, context.target)
    owners = {n: n for n in context.target.nodes}
    for n in reversed(main.tree.order):
        if n in owners:
            continue
        children = main.tree.nodes[n]["children"]
        owners[n] = owners.get(children[0], n) if children else n
    return LineagePlan(
        main.tree,
        main.placed,
        main.jumps,
        main.births,
        owners,
        retained,
        descendants,
        main,
        hidden,
        [e for e in context.rec["events"] if e["type"] == "origin"],
    )
