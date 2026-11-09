"""Raw ordered recPhyloXML evidence, preserving unpaired source events.

Unlike the resolved reconciliation document, this profile does not invent or
require transfer pairing. Consumers may retain incomplete source evidence.
"""

from typing import Any

from defusedxml.ElementTree import fromstring

from nwkit.reconciliation_import import _one, _phylogenies, _tag, _xml_tree


def read_event_records(text, *, expected_clades=None):
    root = fromstring(text, forbid_dtd=True, forbid_entities=True, forbid_external=True)
    if _tag(root) != "recPhylo":
        raise ValueError("Expected recPhyloXML")
    upper_phy = _phylogenies(_one(root, "spTree", required=True))
    lower_phy = _phylogenies(_one(root, "recGeneTree", required=True))
    if len(upper_phy) != 1 or len(lower_phy) != 1:
        raise ValueError("Expected one species and one gene phylogeny")
    upper, lower = _xml_tree(upper_phy[0]), _xml_tree(lower_phy[0])
    species, order, stack = {}, [], [upper]
    content: dict[Any, list[str]] = {}
    while stack:
        node = stack.pop()
        if not node.label or node.label in species:
            raise ValueError("Duplicate XML species label: " + node.label)
        species[node.label] = node
        order.append(node)
        stack.extend(reversed(node.children))
    for node in reversed(order):
        tips = (
            [node.label]
            if not node.children
            else [tip for child in node.children for tip in content[child]]
        )
        if len(tips) != len(set(tips)):
            raise ValueError("Duplicate XML species tip")
        content[node] = tips
    species_taxa = {name: content[node] for name, node in species.items()}
    observed = [frozenset(tips) for tips in species_taxa.values()]
    if len(set(observed)) != len(observed):
        raise ValueError("Repeated XML species clade")
    if expected_clades is not None and set(observed) != {
        frozenset(tips) for tips in expected_clades
    }:
        raise ValueError("XML species topology/taxa differ from the reference tree")
    rows: list[dict[str, Any]] = []
    stack = [(lower, None)]
    seen_genes = set()
    while stack:
        node, parent = stack.pop()
        ordinal = len(rows)
        if not node.events:
            raise ValueError("Gene clade has no reconciliation event")
        events = []
        for event in node.events:
            kind, attrs = event["tag"], event["attributes"]
            field = (
                "destinationSpecies" if kind == "transferBack" else "speciesLocation"
            )
            location = attrs.get(field)
            if location not in species:
                raise ValueError(f"Unmapped {kind} species: {location}")
            events.append({"type": kind, "species": location, "attributes": attrs})
        leaves = [event for event in events if event["type"] == "leaf"]
        if leaves:
            if (
                node.children
                or len(leaves) != 1
                or node.label in seen_genes
                or len(species_taxa[leaves[0]["species"]]) != 1
            ):
                raise ValueError("Invalid or duplicate extant gene leaf: " + node.label)
            seen_genes.add(node.label)
        elif not node.children and not any(event["type"] == "loss" for event in events):
            raise ValueError("Terminal gene clade is neither a leaf nor a loss")
        rows.append(
            dict(ordinal=ordinal, parent=parent, label=node.label, events=events)
        )
        stack.extend((child, ordinal) for child in reversed(node.children))
    return species_taxa, rows
