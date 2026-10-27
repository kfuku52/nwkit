"""Expand ordered events sharing one recorded H node into real display forks."""

from collections import defaultdict

from nwkit.reconciliation_history import tree_record


def expand_event_chains(part, context):
    grouped = defaultdict(list)
    for event in context.rec["events"]:
        if event["id"] in part.events:
            grouped[event["lower_node"]].append(event)
    nodes = {
        n: dict(v, children=list(v["children"])) for n, v in part.tree.nodes.items()
    }
    part.entries = {source: n for n, source in part.sources.items()}
    root = part.tree.root
    for source, events in grouped.items():
        node = part.entries[source]
        origin = next((e for e in events if e["type"] == "origin"), None)
        biological = [e for e in events if e["type"] != "origin"]
        if not biological:
            continue
        # Same-branch origins annotate the incoming stem without an extra fork.
        chain = list(biological)
        if origin and origin["upper_node"] != biological[0]["upper_node"]:
            chain.insert(0, origin)
        prefixes = []
        for event in chain[:-1]:
            key = f"display-event:{part.tree.id}:{event['id']}:{node}"
            if key in nodes:
                raise ValueError("Event display-node ID collides with a recorded node.")
            nodes[key] = dict(id=key, label="", children=[], branch_length=None)
            part.placed[key] = event["upper_node"]
            part.sources[key] = source
            part.events[event["id"]] = key
            prefixes.append(key)
        part.events[chain[-1]["id"]] = node
        part.placed[node] = chain[-1]["upper_node"]
        if prefixes:
            first = prefixes[0]
            if node == root:
                root = first
            else:
                parent = part.tree.parents[node]
                nodes[parent]["children"] = [
                    first if c == node else c for c in nodes[parent]["children"]
                ]
            for before, after in zip(prefixes, [*prefixes[1:], node], strict=True):
                nodes[before]["children"] = [after]
            part.entries[source] = first
        if origin and origin not in chain:
            part.events[origin["id"]] = part.entries[source]
    part.tree = tree_record(part.tree.record, nodes, root)
