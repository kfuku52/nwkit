"""Explicit adapters to the common document; preserves source event semantics."""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter
from dataclasses import dataclass, field
from typing import Any

from defusedxml.ElementTree import fromstring

from nwkit.reconciliation_xml import serialize_xml


@dataclass(eq=False)
class RawNode:
    label: str = ""
    length: float | None = None
    children: list[RawNode] = field(default_factory=list)
    events: list[dict] = field(default_factory=list)
    annotations: dict = field(default_factory=dict)


def empty_document(title=""):
    return {
        "schema_version": "0.1.0",
        "kind": "reconciliation_document",
        "title": title,
        "trees": [],
        "reconciliations": [],
        "transformations": [],
        "provenance": {},
    }


def flatten_tree(root, tree_id, role, unit="unspecified"):
    order, stack = [], [root]
    while stack:
        current = stack.pop()
        order.append(current)
        stack.extend(current.children)
    digests: dict[RawNode, str] = {}
    signatures: dict[str, str] = {}
    for node in reversed(order):
        signature = json.dumps(
            [
                node.label if not node.children else "",
                sorted(digests[c] for c in node.children),
            ],
            ensure_ascii=False,
        )
        digest = hashlib.sha256(signature.encode()).hexdigest()
        if digest in signatures and signatures[digest] != signature:
            raise ValueError("Structural node digest collision.")
        signatures[digest] = signature
        digests[node] = digest
    # Canonical order makes IDs independent of XML child order. Identical
    # unobserved subtrees have explicit occurrence IDs, never array indices.
    counts: Counter = Counter()
    ids, canonical, stack = {}, [], [root]
    while stack:
        node = stack.pop()
        digest = digests[node]
        counts[digest] += 1
        ids[node] = "n:" + digest + (f":{counts[digest]}" if counts[digest] > 1 else "")
        canonical.append(node)
        stack.extend(
            reversed(sorted(node.children, key=lambda c: (digests[c], c.label)))
        )
    nodes = [
        {
            "id": ids[n],
            "label": n.label,
            "children": [
                ids[c] for c in sorted(n.children, key=lambda c: (digests[c], c.label))
            ],
            "branch_length": n.length,
            "annotations": n.annotations,
        }
        for n in canonical
    ]
    return {
        "id": tree_id,
        "role": role,
        "root": ids[root],
        "branch_length_unit": unit,
        "nodes": nodes,
    }, ids


def _raw_ete(tree):
    originals = list(tree.traverse("postorder"))
    raw: dict[Any, RawNode] = {}
    for node in originals:
        raw[node] = RawNode(
            str(node.name or ""),
            node.dist,
            [raw[c] for c in node.children],
            annotations=dict(node.props),
        )
    return raw[tree], raw


def import_tree_document(text):
    from ete4 import Tree

    tree = Tree(text.strip(), parser=1)
    raw, _ = _raw_ete(tree)
    record, _ = flatten_tree(raw, "tree", "standalone")
    result = empty_document("Tree")
    result["trees"].append(record)
    result["provenance"] = {
        "format": "newick",
        "sha256": hashlib.sha256(text.encode()).hexdigest(),
    }
    return result


def _tag(element):
    return element.tag.rsplit("}", 1)[-1]


def _check_namespace(element):
    if element.tag.startswith("{"):
        namespace = element.tag[1:].split("}", 1)[0]
        if namespace not in {
            "http://www.phyloxml.org",
            "http://www.recgenetreexml.org",
            "http://www.recspeciestreexml.org",
            "http://www.recg.org",
        }:
            raise ValueError(f"Unsupported XML core namespace: {namespace}.")


def _children(element, name):
    children = [child for child in element if _tag(child) == name]
    for child in children:
        _check_namespace(child)
    return children


def _one(element, name, *, required=False):
    matches = _children(element, name)
    if len(matches) != 1:
        if matches or required:
            raise ValueError(f"Expected exactly one XML {name}.")
        return None
    return matches[0]


def _xml_tree(phylogeny):
    if phylogeny.attrib.get("rooted", "true") not in {"true", "1"}:
        raise ValueError("Reconciliation documents require a rooted XML tree.")
    root = _one(phylogeny, "clade", required=True)
    order, stack = [], [root]
    while stack:
        element = stack.pop()
        order.append(element)
        stack.extend(_children(element, "clade"))
    raw: dict[Any, RawNode] = {}
    allowed_events = {
        "leaf",
        "speciation",
        "duplication",
        "loss",
        "branchingOut",
        "transferBack",
    }
    for element in reversed(order):
        name = _one(element, "name")
        length = _one(element, "branch_length")
        value = element.attrib.get("branch_length")
        if length is not None:
            if value is not None:
                raise ValueError("Duplicate XML branch length representations.")
            value = length.text
        events = _one(element, "eventsRec")
        recorded = []
        if events is not None:
            for event in events:
                _check_namespace(event)
                tag = _tag(event)
                if tag not in allowed_events:
                    raise ValueError(f"Unsupported recPhyloXML event: {tag}.")
                recorded.append({"tag": tag, "attributes": dict(event.attrib)})
        label = (name.text or "") if name is not None else ""
        leaf_events = [e for e in recorded if e["tag"] == "leaf"]
        if leaf_events and leaf_events[-1]["attributes"].get("geneName"):
            gene_name = leaf_events[-1]["attributes"]["geneName"]
            if label and label != gene_name:
                raise ValueError("XML name and geneName disagree.")
            label = gene_name
        opaque = [
            serialize_xml(c)
            for c in element
            if _tag(c) not in {"name", "branch_length", "clade", "eventsRec"}
        ]
        annotations: dict[str, Any] = {"opaque_xml": opaque} if opaque else {}
        if recorded and recorded[-1]["tag"] == "loss":
            annotations["observed"] = False
        raw[element] = RawNode(
            label,
            None if value is None else float(value),
            [raw[c] for c in _children(element, "clade")],
            recorded,
            annotations,
        )
    return raw[root]


def _phylogenies(container):
    phylogenies = _children(container, "phylogeny")
    if not phylogenies:
        raise ValueError("XML container has no phylogeny.")
    return phylogenies


def _species_names(root, ids):
    names = {}
    for node, nid in ids.items():
        if node.label:
            if node.label in names:
                raise ValueError(
                    "recPhyloXML upper names must be unique for speciesLocation references."
                )
            names[node.label] = nid
    return names


def _xml_events(ids, upper_names):
    events: list[dict] = []
    by_node: dict[RawNode, list[dict]] = {}
    mapping = {"branchingOut": "transfer_out", "transferBack": "transfer_in"}
    for node, nid in ids.items():
        by_node[node] = []
        for index, raw in enumerate(node.events):
            attrs, tag = raw["attributes"], raw["tag"]
            location = (
                attrs.get("destinationSpecies")
                if tag == "transferBack"
                else attrs.get("speciesLocation")
            )
            if not location or location not in upper_names:
                raise ValueError(
                    f"Unknown or missing XML species reference: {location!r}."
                )
            event = {
                "id": f"{nid}:e{index}",
                "type": mapping.get(tag, tag),
                "lower_node": nid,
                "upper_node": upper_names[location],
                "annotations": {"recphyloxml": raw},
            }
            if "confidence" in attrs and not math.isfinite(float(attrs["confidence"])):
                raise ValueError("Nonfinite source XML confidence.")
            # Confidence has no declared conditioning or scale in this profile.
            # Preserve it verbatim as a source attribute, not as a probability.
            events.append(event)
            by_node[node].append(event)
    paired = set()
    for node in ids:
        outgoing = [e for e in by_node[node] if e["type"] == "transfer_out"]
        if outgoing:
            targets = [
                child
                for child in node.children
                if by_node[child] and by_node[child][0]["type"] == "transfer_in"
            ]
            if len(outgoing) != 1 or len(targets) != 1:
                raise ValueError(
                    "Unsupported unpaired/ambiguous XML transfer; no recipient is invented."
                )
            child = targets[0]
            outgoing[0].update(
                donor_node=outgoing[0]["upper_node"],
                recipient_node=by_node[child][0]["upper_node"],
                recipient_lower_node=ids[child],
            )
            paired.add(by_node[child][0]["id"])
    if any(e["type"] == "transfer_in" and e["id"] not in paired for e in events):
        raise ValueError("Unsupported unpaired XML transferBack.")
    return events


def import_xml_document(text):
    root = fromstring(text, forbid_dtd=True, forbid_entities=True, forbid_external=True)
    _check_namespace(root)
    result = empty_document("Reconciliation" if _tag(root) == "recPhylo" else "Tree")
    result["provenance"] = {
        "format": _tag(root),
        "sha256": hashlib.sha256(text.encode()).hexdigest(),
    }
    if _tag(root) == "phyloxml":
        for index, phy in enumerate(_phylogenies(root)):
            tree, _ = flatten_tree(
                _xml_tree(phy),
                f"tree:{index}",
                "standalone",
                phy.attrib.get("branch_length_unit", "unspecified"),
            )
            result["trees"].append(tree)
        return result
    if _tag(root) != "recPhylo":
        raise ValueError("Unsupported XML root; expected phyloxml or recPhylo.")
    if any(_tag(c) not in {"spTree", "recGeneTree"} for c in root):
        raise ValueError(
            "Unsupported recPhyloXML top-level element; refusing to omit it."
        )
    upper_container = _one(root, "spTree", required=True)
    upper_phy = _phylogenies(upper_container)
    if len(upper_phy) != 1:
        raise ValueError("recPhyloXML requires exactly one upper phylogeny.")
    upper_root = _xml_tree(upper_phy[0])
    upper, upper_ids = flatten_tree(
        upper_root,
        "upper",
        "upper",
        upper_phy[0].attrib.get("branch_length_unit", "unspecified"),
    )
    result["trees"].append(upper)
    names = _species_names(upper_root, upper_ids)
    for container in _children(root, "recGeneTree"):
        for phy in _phylogenies(container):
            index = len(result["reconciliations"])
            lower, ids = flatten_tree(
                _xml_tree(phy),
                f"lower:{index}",
                "lower",
                phy.attrib.get("branch_length_unit", "unspecified"),
            )
            result["trees"].append(lower)
            result["reconciliations"].append(
                {
                    "id": f"history:{index}",
                    "lower_tree": lower["id"],
                    "upper_tree": upper["id"],
                    "history_kind": "unspecified",
                    "selection": {},
                    "model": {"source_format": "recPhyloXML", "conditions": {}},
                    "events": _xml_events(ids, names),
                }
            )
    if not result["reconciliations"]:
        raise ValueError("recPhyloXML has no reconciled trees.")
    return result


def import_saved_tsv(args):
    from nwkit.clade_index import CladeIndex
    from nwkit.result_plot_data import read_result_table, reconciliation_plot_data
    from nwkit.util import read_tree

    gene = read_tree(
        args.infile, args.format, args.quoted_node_names, rooted=args.input_rooted
    )
    species = read_tree(
        args.species_tree,
        args.species_tree_format,
        True,
        rooted=args.species_tree_rooted,
    )
    table = read_result_table(args.reconciliation)
    data = reconciliation_plot_data(gene, species, table)
    result, indices = empty_document("Saved reconciliation"), {}
    for tid, tree in (("upper", species), ("lower", gene)):
        raw, nodes = _raw_ete(tree)
        record, ids = flatten_tree(raw, tid, tid)
        result["trees"].append(record)
        ci = CladeIndex(tree)
        indices[tid] = {ci.clade_id_for_node(n): ids[nodes[n]] for n in nodes}
    events, summaries = [], []
    for clade, row in data.events.items():
        event = {
            "id": clade,
            "type": row["event_type"],
            "lower_node": indices["lower"][clade],
            "upper_node": indices["upper"].get(row["species_event_id"]),
            "annotations": {"nwkit_tsv": row},
        }
        # Legacy tables report a mapped transfer node, not a donor/recipient pair.
        # This remains unresolved rather than becoming a fabricated full transfer.
        if event["type"] == "transfer":
            raise ValueError(
                "Legacy TSV transfer lacks a donor/recipient pair; use the existing paired report or an explicit document."
            )
        events.append(event)
        if row.get("implied_losses") not in (None, "", "NA"):
            summaries.append(
                {
                    "lower_node": indices["lower"][clade],
                    "implied_losses": int(row["implied_losses"]),
                }
            )
    sources = sorted(
        {str(r.get("event_source", "unreported")) for r in data.events.values()}
    )
    result["reconciliations"].append(
        {
            "id": "saved",
            "lower_tree": "lower",
            "upper_tree": "upper",
            "history_kind": "unspecified",
            "selection": {},
            "model": {"event_sources": sources, "conditions": {}},
            "summary": {"implied_losses": summaries},
            "events": events,
        }
    )
    return result
