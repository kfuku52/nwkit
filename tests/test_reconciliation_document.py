"""File-contract semantics, compatibility and rejection without inference."""

import copy
import json
from pathlib import Path

import jsonschema
import pytest

from nwkit.reconciliation_document import (
    document_json,
    read_document,
    validate_document,
)
from nwkit.reconciliation_export import recphyloxml_text
from nwkit.reconciliation_import import import_xml_document

ROOT = Path(__file__).resolve().parents[1]
TOY = ROOT / "examples/reconciliation/toy.recphyloxml"


@pytest.fixture
def document():
    return read_document(TOY)


def test_schema_and_imported_event_contract(document):
    schema = json.loads(
        (ROOT / "nwkit/data_schema/reconciliation_document.schema.json").read_text()
    )
    jsonschema.Draft202012Validator.check_schema(schema)
    jsonschema.Draft202012Validator(schema).validate(document.record)
    events = document.record["reconciliations"][0]["events"]
    assert [e["type"] for e in events].count("loss") == 1
    assert [e["type"] for e in events].count("duplication") == 1
    donor = next(e for e in events if e["type"] == "transfer_out")
    recipient = next(e for e in events if e["type"] == "transfer_in")
    upper = document.trees["upper"]
    assert upper.nodes[donor["donor_node"]]["label"] == "C"
    assert upper.nodes[donor["recipient_node"]]["label"] == "A"
    assert donor["recipient_lower_node"] == recipient["lower_node"]
    assert all(e.get("support") is None for e in events)
    assert all(e.get("time") is None for e in events)
    assert document.record["reconciliations"][0]["history_kind"] == "unspecified"


@pytest.mark.parametrize("invalid", ["selection", "support", "display"])
def test_schema_rejects_missing_display_and_probability_meaning(document, invalid):
    record = copy.deepcopy(document.record)
    rec = record["reconciliations"][0]
    if invalid == "selection":
        rec.update(history_kind="representative", selection={"method": 1})
    elif invalid == "support":
        rec["events"][0]["support"] = {
            "kind": "posterior",
            "value": 0.9,
            "conditional_on": None,
        }
    else:
        record["provenance"]["display_selection"] = {"tree": "upper"}
    schema = json.loads(
        (ROOT / "nwkit/data_schema/reconciliation_document.schema.json").read_text()
    )
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.Draft202012Validator(schema).validate(record)
    with pytest.raises(ValueError):
        validate_document(record)


def test_loaded_display_selection_rejects_unknown_and_upper_references(document):
    record = copy.deepcopy(document.record)
    selection = {
        "tree": "upper",
        "node": document.trees["upper"].root,
        "inference_recomputed": False,
        "summaries": "full source",
    }
    record["provenance"]["display_selection"] = selection
    with pytest.raises(ValueError, match="standalone or lower"):
        validate_document(record)
    selection["node"] = "missing"
    with pytest.raises(ValueError, match="reference"):
        validate_document(record)


def test_schema_inventory_and_reserved_comparison_contract():
    import hashlib

    manifest = json.loads(
        (ROOT / "docs/guides/reconciliation_contract_manifest.json").read_text()
    )
    schemas = {}
    for filename, digest in manifest["schemas"].items():
        data = (ROOT / "nwkit/data_schema" / filename).read_bytes()
        assert hashlib.sha256(data).hexdigest() == digest
        schemas[filename] = json.loads(data)
        jsonschema.Draft202012Validator.check_schema(schemas[filename])
    runtime = {"seconds": None, "measured_at": None, "conditions": {}}
    comparison = {
        "schema_version": "0.1.0",
        "kind": "reconciliation_comparison",
        "mode": "general",
        "panels": [
            {
                "id": "any-method",
                "label": "Any method",
                "document": None,
                "document_sha256": None,
                "inference_status": "success",
                "history_status": "missing",
                "inference_runtime": runtime,
                "history_runtime": runtime,
            }
        ],
        "shared_context": {
            "upper_tree_sha256": None,
            "gene_tree_sha256": None,
            "mapping_sha256": None,
        },
        "view": {
            "schema_version": "0.1.0",
            "kind": "reconciliation_view",
            "layout": "topology_embedded",
        },
    }
    validator = jsonschema.Draft202012Validator(
        schemas["reconciliation_comparison.schema.json"]
    )
    validator.validate(comparison)
    comparison["panels"][0]["history_status"] = "available"
    with pytest.raises(jsonschema.ValidationError):
        validator.validate(comparison)


def test_json_and_xml_round_trips_preserve_histories(document, tmp_path):
    json_path, xml_path = tmp_path / "doc.json", tmp_path / "doc.xml"
    json_path.write_text(document_json(document))
    assert read_document(json_path).record == document.record
    xml_path.write_text(recphyloxml_text(document))
    reread = read_document(xml_path)
    assert reread.record["trees"] == document.record["trees"]
    assert reread.record["reconciliations"] == document.record["reconciliations"]


def test_xml_export_refuses_to_lose_model_or_history_meaning(document):
    record = copy.deepcopy(document.record)
    record["reconciliations"][0]["model"]["conditions"] = {"origin": "crown"}
    with pytest.raises(ValueError, match="lose model"):
        recphyloxml_text(validate_document(record))


def test_xml_confidence_scale_is_not_converted_to_probability():
    text = TOY.read_text().replace(
        '<duplication speciesLocation="AB"/>',
        '<duplication speciesLocation="AB" confidence="95"/>',
    )
    doc = validate_document(import_xml_document(text))
    event = next(
        e
        for e in doc.record["reconciliations"][0]["events"]
        if e["type"] == "duplication"
    )
    assert event.get("support") is None
    assert event["annotations"]["recphyloxml"]["attributes"]["confidence"] == "95"
    assert 'confidence="95"' in recphyloxml_text(doc)


def test_foreign_xml_event_namespace_is_not_reinterpreted():
    text = TOY.read_text().replace(
        '<duplication speciesLocation="AB"/>',
        '<foreign:duplication xmlns:foreign="urn:other-model" speciesLocation="AB"/>',
    )
    with pytest.raises(ValueError, match="namespace"):
        import_xml_document(text)


def test_child_order_and_internal_rename_do_not_change_structural_ids(document):
    from xml.etree import ElementTree as ET

    from defusedxml.ElementTree import fromstring

    root = fromstring(TOY.read_text())
    for element in root.iter("clade"):
        children = [c for c in element if c.tag == "clade"]
        for child in children:
            element.remove(child)
        element.extend(reversed(children))
    changed = validate_document(
        import_xml_document(ET.tostring(root, encoding="unicode"))
    )
    assert changed.record["trees"] == document.record["trees"]
    assert changed.record["reconciliations"] == document.record["reconciliations"]
    text = TOY.read_text().replace(
        "<name>s1</name>", "<name>different-internal-name</name>"
    )
    renamed = validate_document(import_xml_document(text))
    for tid in document.trees:
        assert set(renamed.trees[tid].nodes) == set(document.trees[tid].nodes)


@pytest.mark.parametrize(
    "mutation",
    [
        lambda r: r.update(schema_version="999"),
        lambda r: r["trees"][0]["nodes"].append(
            copy.deepcopy(r["trees"][0]["nodes"][0])
        ),
        lambda r: r["trees"][0]["nodes"][0].update(branch_length=float("nan")),
        lambda r: r["trees"][0]["nodes"][0]["children"].append(r["trees"][0]["root"]),
        lambda r: r["reconciliations"][0]["events"][0].update(upper_node="foreign"),
        lambda r: r["reconciliations"][0].update(history_kind="representative"),
        lambda r: r["reconciliations"][0]["events"][0].update(
            support={"kind": "posterior", "value": 0.5, "conditional_on": None}
        ),
        lambda r: r["reconciliations"][0]["events"][0].update(
            time={"lower": 2, "upper": 1, "unit": "Ma"}
        ),
    ],
)
def test_invalid_contracts_are_rejected(document, mutation):
    record = copy.deepcopy(document.record)
    mutation(record)
    with pytest.raises(ValueError):
        validate_document(record)


def test_hierarchy_cycles_rejected(document):
    record = copy.deepcopy(document.record)
    rec = copy.deepcopy(record["reconciliations"][0])
    rec.update(id="cycle", lower_tree="upper", upper_tree="lower:0", events=[])
    record["reconciliations"].append(rec)
    with pytest.raises(ValueError, match="hierarchy.*cycle"):
        validate_document(record)


@pytest.mark.parametrize(
    "text,match",
    [
        ('<!DOCTYPE x [<!ENTITY a "external">]><recPhylo/>', "DTDForbidden"),
        (
            TOY.read_text().replace(
                '<duplication speciesLocation="AB"/>', "<bifurcationOut/>"
            ),
            "Unsupported recPhyloXML event",
        ),
        (
            TOY.read_text().replace(
                'destinationSpecies="A"', 'destinationSpecies="missing"'
            ),
            "species reference",
        ),
        (
            TOY.read_text().replace(
                '<branchingOut speciesLocation="C"/>',
                '<duplication speciesLocation="C"/>',
            ),
            "unpaired XML transferBack",
        ),
    ],
)
def test_xml_limits_and_entity_rejection(text, match):
    with pytest.raises(Exception, match=match):
        import_xml_document(text)


def test_newick_phyloxml_labels_and_opaque_properties(tmp_path):
    path = tmp_path / "single.nwk"
    path.write_text("(NA:1,001:2)root;")
    tree = read_document(path).trees["tree"]
    assert {n["label"] for n in tree.nodes.values()} == {"NA", "001", "root"}
    path.write_text(
        '<phyloxml><phylogeny rooted="true"><clade><name>root</name><clade><name>NA</name><property ref="example">literal</property></clade></clade></phylogeny></phyloxml>'
    )
    doc = read_document(path)
    node = next(n for n in doc.trees["tree:0"].nodes.values() if n["label"] == "NA")
    assert "literal" in node["annotations"]["opaque_xml"][0]


def test_duplicate_json_keys_rejected(document, tmp_path):
    path = tmp_path / "bad.json"
    path.write_text(
        document_json(document).replace(
            '"schema_version": "0.1.0"',
            '"schema_version": "0.1.0", "schema_version": "0.1.0"',
        )
    )
    with pytest.raises(ValueError, match="Duplicate JSON key"):
        read_document(path)


def test_validation_detaches_caller_data(document):
    record = copy.deepcopy(document.record)
    checked = validate_document(record)
    record["trees"][0]["nodes"][0]["label"] = "changed"
    assert checked.trees["upper"].nodes[checked.trees["upper"].root]["label"] == "ABC"


def test_x1_contract_validates_boundaries_resolution_and_overlap():
    from nwkit.reconciliation_import import RawNode, empty_document, flatten_tree

    def triple(tid, paired):
        leaves = {label: RawNode(label) for label in "ABC"}
        inner = RawNode(children=[leaves[paired[0]], leaves[paired[1]]])
        outer = RawNode(
            children=[inner, leaves[next(l for l in "ABC" if l not in paired)]]
        )
        tree, ids = flatten_tree(outer, tid, "latent" if tid == "H" else "observed")
        return (
            tree,
            ids[outer],
            ids[inner],
            {label: ids[n] for label, n in leaves.items()},
        )

    h, hr, hi, hl = triple("H", "AB")
    g, gr, gi, gl = triple("G", "AC")
    record = empty_document("Synthetic X1 contract fixture")
    record["trees"] = [h, g]
    block = {
        "id": "x",
        "source_nodes": [hr, hi],
        "target_nodes": [gr, gi],
        "resolution": 1,
        "boundaries": [{"source_node": hl[l], "target_node": gl[l]} for l in "ABC"],
    }
    record["transformations"] = [
        {
            "id": "X",
            "source_tree": "H",
            "target_tree": "G",
            "model_id": "X1",
            "blocks": [block],
        }
    ]
    validate_document(record)
    bad = copy.deepcopy(record)
    bad["transformations"][0]["blocks"][0]["resolution"] = 2
    with pytest.raises(ValueError, match="target local split"):
        validate_document(bad)
    bad = copy.deepcopy(record)
    second = copy.deepcopy(block)
    second["id"] = "overlap"
    bad["transformations"][0]["blocks"].append(second)
    with pytest.raises(ValueError, match="overlap"):
        validate_document(bad)
