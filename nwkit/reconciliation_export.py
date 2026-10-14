"""recPhyloXML round trip for the supported imported profile, without data loss."""

from xml.etree import ElementTree as ET

from defusedxml.ElementTree import fromstring

from nwkit.reconciliation_document import document_json


def _xml_clades(tree, events=None):
    elements = {}
    for nid in tree.order:
        node = tree.nodes[nid]
        element = ET.Element("clade")
        if node["label"]:
            ET.SubElement(element, "name").text = node["label"]
        if node["branch_length"] is not None:
            ET.SubElement(element, "branch_length").text = repr(node["branch_length"])
        extra = node.get("annotations", {})
        if set(extra) - {"observed", "opaque_xml"}:
            raise ValueError(
                "XML export cannot represent node extension annotations; retain canonical JSON."
            )
        for opaque in extra.get("opaque_xml", []):
            element.append(
                fromstring(
                    opaque, forbid_dtd=True, forbid_entities=True, forbid_external=True
                )
            )
        if events and nid in events:
            event_element = ET.SubElement(element, "eventsRec")
            for raw in events[nid]:
                ET.SubElement(event_element, raw["tag"], raw["attributes"])
        elements[nid] = element
    for nid in tree.order:
        for child in tree.nodes[nid]["children"]:
            elements[nid].append(elements[child])
    return elements[tree.root]


def recphyloxml_text(document):
    record = document.record
    recs = record["reconciliations"]
    if (
        not recs
        or record["transformations"]
        or record.get("annotations")
        or record["provenance"].get("display_selection")
    ):
        raise ValueError(
            "XML export requires complete two-level imported recPhyloXML histories without transformations/extensions."
        )
    upper_ids = {r["upper_tree"] for r in recs}
    if len(upper_ids) != 1 or set(document.trees) != upper_ids | {
        r["lower_tree"] for r in recs
    }:
        raise ValueError("XML export requires one upper tree and its lower histories.")
    upper = document.trees[next(iter(upper_ids))]
    root = ET.Element("recPhylo")
    species = ET.SubElement(ET.SubElement(root, "spTree"), "phylogeny", rooted="true")
    if upper.record["branch_length_unit"] != "unspecified":
        species.set("branch_length_unit", upper.record["branch_length_unit"])
    species.append(_xml_clades(upper))
    names = {nid: n["label"] for nid, n in upper.nodes.items()}
    for rec in recs:
        if (
            rec["history_kind"] != "unspecified"
            or rec["selection"]
            or rec.get("summary")
            or rec.get("annotations")
            or rec["model"] != {"source_format": "recPhyloXML", "conditions": {}}
        ):
            raise ValueError(
                "XML export would lose model/history/summary metadata; retain canonical JSON."
            )
        lower = document.trees[rec["lower_tree"]]
        events: dict[str, list[dict]] = {}
        for event in rec["events"]:
            annotations = event.get("annotations", {})
            if set(annotations) != {"recphyloxml"}:
                raise ValueError(
                    "XML export requires recorded source XML event semantics."
                )
            raw = annotations["recphyloxml"]
            tag = raw["tag"]
            expected = {
                "branchingOut": "transfer_out",
                "transferBack": "transfer_in",
            }.get(tag, tag)
            if event["type"] != expected or event.get("time") is not None:
                raise ValueError(
                    "XML event metadata differs from canonical event semantics."
                )
            attrs = dict(raw["attributes"])
            reference = (
                "destinationSpecies" if tag == "transferBack" else "speciesLocation"
            )
            if event["upper_node"] is None or not names[event["upper_node"]]:
                raise ValueError("XML species references require named upper nodes.")
            attrs[reference] = names[event["upper_node"]]
            if (
                "confidence" in attrs
                and event.get("support")
                and event["support"].get("value") != float(attrs["confidence"])
            ):
                raise ValueError("XML confidence and canonical support differ.")
            if event.get("support") and "confidence" not in attrs:
                raise ValueError("XML export would lose event support.")
            events.setdefault(event["lower_node"], []).append(
                {"tag": tag, "attributes": attrs}
            )
        phy = ET.SubElement(
            ET.SubElement(root, "recGeneTree"), "phylogeny", rooted="true"
        )
        if lower.record["branch_length_unit"] != "unspecified":
            phy.set("branch_length_unit", lower.record["branch_length_unit"])
        phy.append(_xml_clades(lower, events))
    return ET.tostring(root, encoding="unicode", xml_declaration=True) + "\n"


def document_text(document, path):
    if str(path).lower().endswith((".xml", ".recphyloxml")):
        return recphyloxml_text(document)
    return document_json(document)
