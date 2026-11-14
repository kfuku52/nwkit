"""Adapters for numbered phylogenetic edge tables (including ape)."""

from ete4 import Tree

from nwkit.clade_index import CladeIndex
from nwkit.util import validate_unique_named_leaves


def ape_tree(data):
    tips = data["tip.label"]
    nnode = data["Nnode"]
    if not isinstance(tips, list) or not tips or type(nnode) is not int or nnode < 1:
        raise ValueError("Expected ape tip labels and a positive internal node count")
    labels = data.get("node.label") or []
    if not isinstance(labels, list):
        raise ValueError("Expected an array of ape node labels")
    size = len(tips) + nnode
    edges = data["edge"]
    parents = {}
    children: dict[int, list[int]] = {}
    for edge in edges:
        if len(edge) != 2 or any(
            type(node) is not int or not 1 <= node <= size for node in edge
        ):
            raise ValueError("Invalid ape edge indices")
        parent, child = edge
        if parent <= len(tips) or child in parents:
            raise ValueError("Invalid ape parent or repeated child")
        parents[child] = parent
        children.setdefault(parent, []).append(child)
    roots = set(range(1, size + 1)) - parents.keys()
    if len(roots) != 1 or len(edges) != size - 1:
        raise ValueError("Expected one connected ape tree")
    root = roots.pop()
    order, pending, seen = [], [root], set()
    while pending:
        node = pending.pop()
        if node in seen:
            raise ValueError("Cycle in ape edges")
        seen.add(node)
        order.append(node)
        pending.extend(children.get(node, ()))
    if len(seen) != size or any(
        node not in children for node in range(len(tips) + 1, size + 1)
    ):
        raise ValueError("Disconnected or empty internal ape node")
    nodes = {node: Tree() for node in order}
    for node in order:
        value = (
            tips[node - 1]
            if node <= len(tips)
            else (
                labels[node - len(tips) - 1] if node - len(tips) <= len(labels) else ""
            )
        )
        nodes[node].name = "" if value is None else str(value)
        for child in children.get(node, ()):
            nodes[node].add_child(nodes[child])
    return nodes[root], [nodes[node] for node in range(len(tips) + 1, size + 1)]


def numbered_clade_keys(data, tip_set=None, *, compact=False):
    tree, internal = ape_tree(data)
    validate_unique_named_leaves(tree, "numbered tree")
    if compact:
        import hashlib
        import json

        catalog = sorted(set(data["tip.label"] if tip_set is None else tip_set))
        bit_by_name = {name: 1 << index for index, name in enumerate(catalog)}
        masks: dict[object, int] = {}
        for node in tree.traverse(strategy="postorder"):
            mask = bit_by_name.get(node.name, 0) if node.is_leaf else 0
            for child in node.children:
                mask |= masks[child]
            masks[node] = mask
        catalog_id = hashlib.sha256(
            json.dumps(catalog, ensure_ascii=False).encode()
        ).hexdigest()
        prefix = "clade-bits-v1:" + catalog_id + ":"
        return [None] * len(data["tip.label"]) + [
            prefix + format(masks[node], "x") if masks[node].bit_count() >= 2 else None
            for node in internal
        ]
    index = CladeIndex(tree)
    allowed = set(index.names if tip_set is None else tip_set)
    mask = sum(1 << i for i, name in enumerate(index.names) if name in allowed)
    result: list[str | None] = [None] * (len(data["tip.label"]) + data["Nnode"])
    for number, node in enumerate(internal, start=len(data["tip.label"])):
        projected = index.mask_by_node[node] & mask
        if projected.bit_count() >= 2:
            result[number] = "\t".join(index.names_for_mask(projected))
    return result
