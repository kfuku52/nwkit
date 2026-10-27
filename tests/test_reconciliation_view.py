"""Generic comparison, display state, multi-layer layout and offline export."""

import base64
import copy
import hashlib
import io
import json
from pathlib import Path
from xml.etree import ElementTree as ET

import pytest
from pypdf import PdfReader

from nwkit.cli import main
from nwkit.reconciliation_composition import apply_view, comparison_scene
from nwkit.reconciliation_document import read_document, validate_document
from nwkit.reconciliation_html import _svg
from nwkit.reconciliation_scene import build_scene
from nwkit.reconciliation_view import read_comparison, validate_view

TOY = Path(__file__).parents[1] / "examples/reconciliation/toy.json"


def comparison(
    tmp_path, statuses=("available", "missing", "failed", "timeout", "unsupported")
):
    source = tmp_path / "source.json"
    source.write_bytes(TOY.read_bytes())
    runtime = {"seconds": None, "measured_at": None, "conditions": {}}
    panels = [
        {
            "id": str(i),
            "label": f"method {i}",
            "document": "source.json" if state == "available" else None,
            "document_sha256": hashlib.sha256(source.read_bytes()).hexdigest()
            if state == "available"
            else None,
            "inference_status": "success",
            "history_status": state,
            "inference_runtime": runtime,
            "history_runtime": runtime,
        }
        for i, state in enumerate(statuses)
    ]
    record = {
        "schema_version": "0.1.0",
        "kind": "reconciliation_comparison",
        "mode": "general",
        "panels": panels,
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
    path = tmp_path / "comparison.json"
    path.write_text(json.dumps(record))
    return path, record


def test_comparison_retains_missing_and_separate_timing(tmp_path):
    path, record = comparison(tmp_path)
    loaded, docs, inputs = read_comparison(path)
    assert len(loaded["panels"]) == 5 and len(inputs) == 1
    scene = comparison_scene(loaded, docs, width=720, font_size=8)
    assert sum(p.id.startswith("status:") for p in scene.primitives) == 4
    assert (
        sum(
            "inference runtime=" in n and "history runtime=" in n for n in scene.notices
        )
        == 5
    )
    main(
        [
            "draw",
            "--comparison",
            str(path),
            "--outfile",
            str(tmp_path / "comparison.pdf"),
        ]
    )
    assert PdfReader(tmp_path / "comparison.pdf").pages
    import zipfile

    outfile = tmp_path / "comparison.html"
    main(["draw", "--comparison", str(path), "--outfile", str(outfile)])
    encoded = (
        outfile.read_text().split('id="sourcebundle">', 1)[1].split("</script>", 1)[0]
    )
    with zipfile.ZipFile(io.BytesIO(base64.b64decode(encoded))) as archive:
        manifest = json.loads(archive.read("comparison.json"))
        for panel in manifest["panels"]:
            if panel["document"] is not None:
                assert (
                    hashlib.sha256(archive.read(panel["document"])).hexdigest()
                    == panel["document_sha256"]
                )

    source = tmp_path / "source.json"
    source.write_text(source.read_text() + " ")
    with pytest.raises(ValueError, match="hash mismatch"):
        read_comparison(path)


@pytest.mark.parametrize(
    "mutation",
    ["duplicate", "missing_hash", "fixed_context", "runtime", "unknown_view"],
)
def test_comparison_rejects_invalid_semantics(tmp_path, mutation):
    path, record = comparison(tmp_path)
    if mutation == "duplicate":
        record["panels"][1]["id"] = "0"
    elif mutation == "missing_hash":
        record["panels"][0]["document_sha256"] = None
    elif mutation == "fixed_context":
        record["mode"] = "fixed_gene_tree"
    elif mutation == "runtime":
        record["panels"][0]["inference_runtime"]["seconds"] = -1
    else:
        record["view"]["inference"] = True
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError):
        read_comparison(path)


def test_no_model_aware_diff_across_different_conditions(tmp_path):
    path, record = comparison(tmp_path, ("available", "available"))
    doc = read_document(TOY).record
    doc["reconciliations"][0]["model"]["conditions"]["eta"] = 0.3
    other = tmp_path / "other.json"
    other.write_text(json.dumps(doc))
    record["panels"][1].update(
        document="other.json",
        document_sha256=hashlib.sha256(other.read_bytes()).hexdigest(),
    )
    path.write_text(json.dumps(record))
    loaded, docs, _ = read_comparison(path)
    scene = comparison_scene(loaded, docs, width=720, font_size=8)
    assert any("conditions/history selection differ" in n for n in scene.notices)


def test_display_state_does_not_rewrite_science(tmp_path):
    doc = read_document(TOY)
    before = copy.deepcopy(doc.record)
    tree = doc.trees["lower:0"]
    leaves = [nid for nid, n in tree.nodes.items() if not n["children"]]
    view = {
        "schema_version": "0.1.0",
        "kind": "reconciliation_view",
        "layout": "paired",
        "leaf_order": {tree.id: list(reversed(leaves))},
        "collapsed_clades": [{"tree": tree.id, "node": tree.root}],
        "event_colors": {"duplication": "#123456"},
    }
    validate_view(view, doc)
    full = build_scene(doc, view=view)
    partial = apply_view(copy.deepcopy(full), doc, view)
    assert len(partial.primitives) < len(full.primitives)
    assert doc.record == before
    bad = copy.deepcopy(view)
    bad["leaf_order"][tree.id] = leaves[:-1]
    with pytest.raises(ValueError, match="every leaf"):
        validate_view(bad, doc)


def test_three_levels_preserve_middle_ids():
    doc = read_document(TOY).record
    middle = copy.deepcopy(doc["trees"][0])
    middle["id"] = "middle"
    doc["trees"].append(middle)
    doc["reconciliations"][0]["upper_tree"] = "middle"
    doc["reconciliations"].append(
        {
            "id": "middle_to_upper",
            "lower_tree": "middle",
            "upper_tree": "upper",
            "history_kind": "unspecified",
            "model": {},
            "selection": {},
            "events": [
                {
                    "id": n["id"],
                    "type": "speciation" if n["children"] else "leaf",
                    "lower_node": n["id"],
                    "upper_node": n["id"],
                }
                for n in middle["nodes"]
            ],
        }
    )
    validated = validate_document(doc)
    scene = build_scene(validated)
    assert any("Three-layer nested topology" in n for n in scene.notices)
    assert len(scene.event_points) == sum(
        len(r["events"]) for r in doc["reconciliations"]
    )
    assert validated.record == doc


def test_offline_html_and_vector_pdf_are_self_contained(tmp_path):
    doc = read_document(TOY).record
    doc["title"] = '</script><img src="https://invalid.test" onerror="alert(1)">'
    source = tmp_path / "source.json"
    source.write_text(json.dumps(doc))
    outfile, view = tmp_path / "offline.html", tmp_path / "view.json"
    main(
        [
            "draw",
            "--document",
            str(source),
            "--outfile",
            str(outfile),
            "--view-out",
            str(view),
        ]
    )
    text = outfile.read_text()
    assert "<script src=" not in text and "<img src=" not in text
    assert "fetch(" not in text and "XMLHttpRequest" not in text
    encoded = text.split('id="pdfbytes">', 1)[1].split("</script>", 1)[0]
    reader = PdfReader(io.BytesIO(base64.b64decode(encoded)))
    assert reader.pages and not reader.pages[0].images
    assert "Download vector PDF" in text and "Save view JSON" in text
    validate_view(json.loads(view.read_text()), read_document(source))
    ET.fromstring(_svg(build_scene(read_document(source)), 8))


def test_comparison_output_cannot_alias_referenced_document(tmp_path):
    path, _ = comparison(tmp_path)
    old = (tmp_path / "source.json").read_bytes()
    with pytest.raises(ValueError):
        main(
            [
                "draw",
                "--comparison",
                str(path),
                "--outfile",
                str(tmp_path / "source.json"),
            ]
        )
    assert (tmp_path / "source.json").read_bytes() == old


def test_three_layers_follow_recorded_middle_transfer_curves(tmp_path):
    import math

    doc = read_document(TOY).record
    middle = doc["trees"][1]
    bottom = copy.deepcopy(middle)
    bottom["id"] = "bottom"
    for n in bottom["nodes"]:
        if not n["children"]:
            n["label"] = "nested_" + n["label"]
    doc["trees"].append(bottom)
    doc["reconciliations"].append(
        {
            "id": "bottom_to_middle",
            "lower_tree": "bottom",
            "upper_tree": middle["id"],
            "history_kind": "unspecified",
            "selection": {},
            "model": {"source": "synthetic geometry only"},
            "events": [
                {
                    "id": n["id"],
                    "type": "speciation" if n["children"] else "leaf",
                    "lower_node": n["id"],
                    "upper_node": n["id"],
                }
                for n in bottom["nodes"]
            ],
        }
    )
    validated = validate_document(doc)
    scene = build_scene(validated)
    assert any(
        p.kind == "upper" and len(p.points) >= 33 and p.owner for p in scene.primitives
    )
    assert all(
        math.isfinite(v) for p in scene.primitives for xy in p.points for v in xy
    )
    assert len(scene.event_points) == sum(
        len(r["events"]) for r in doc["reconciliations"]
    )
    assert min(y for p in scene.primitives for _, y in p.points) >= 0
    source = tmp_path / "three.json"
    source.write_text(json.dumps(doc))
    main(["draw", "--document", str(source), "--outfile", str(tmp_path / "three.pdf")])
    assert PdfReader(tmp_path / "three.pdf").pages


def test_comparison_view_targets_only_existing_panel_clades(tmp_path):
    from nwkit.reconciliation_view import validate_comparison_view

    path, record = comparison(tmp_path)
    doc = read_document(TOY)
    rec = copy.deepcopy(doc.record)
    rec["trees"] = [t for t in rec["trees"] if t["id"] == "upper"]
    rec["reconciliations"] = []
    other = validate_document(rec)
    view = dict(
        record["view"],
        layout="paired",
        collapsed_clades=[{"tree": "lower:0", "node": doc.trees["lower:0"].root}],
    )
    validate_comparison_view(view, {"a": doc, "b": other})
    scene = comparison_scene(
        dict(record, view=view), {"0": doc}, width=720, font_size=8
    )
    assert any("Display-only collapse" in n for n in scene.notices)


@pytest.mark.parametrize(
    "field,value", [("leaf_order", []), ("event_colors", []), ("collapsed_clades", {})]
)
def test_display_state_rejects_bad_container_types(field, value):
    view = {
        "schema_version": "0.1.0",
        "kind": "reconciliation_view",
        "layout": "paired",
        field: value,
    }
    with pytest.raises(ValueError):
        validate_view(view)


def test_clade_selection_preserves_transformation_context():
    from nwkit.reconciliation_document import select_document_clade

    record = read_document(TOY).record
    counterpart = copy.deepcopy(record["trees"][1])
    counterpart["id"] = "G"
    record["trees"].append(counterpart)
    record["transformations"] = [
        {
            "id": "identity",
            "model_id": "identity",
            "source_tree": "lower:0",
            "target_tree": "G",
            "blocks": [],
        }
    ]
    document = validate_document(record)
    node = next(n["id"] for n in record["trees"][1]["nodes"] if n["label"] == "dup")
    selected = select_document_clade(document, f"lower:0/{node}")
    scene = build_scene(selected)
    assert selected.record["trees"] == document.record["trees"]
    assert selected.record["transformations"] == document.record["transformations"]
    assert any(p.owner and p.owner[0] == "G" for p in scene.primitives)
    assert all(
        not p.owner
        or p.owner[0] != "history:0"
        or selected.trees["lower:0"].ancestral(node, p.owner[1])
        for p in scene.primitives
    )


@pytest.mark.parametrize("replacement", [[], {}])
def test_view_scalar_containers_reject_cleanly(replacement):
    base = {
        "schema_version": "0.1.0",
        "kind": "reconciliation_view",
        "layout": "paired",
    }
    for extra in (
        {"layout": replacement},
        {"collapsed_clades": [{"tree": replacement, "node": "x"}]},
        {"output": {"format": replacement}},
    ):
        with pytest.raises(ValueError):
            validate_view(dict(base, **extra))


def test_shared_collapse_rejects_opaque_id_collision_between_different_clades(tmp_path):
    from nwkit.reconciliation_view import validate_comparison_view

    path, record = comparison(tmp_path)
    doc = read_document(TOY)
    other = copy.deepcopy(doc.record)
    for node in other["trees"][1]["nodes"]:
        if not node["children"]:
            node["label"] += "_other"
    view = dict(
        record["view"],
        layout="paired",
        collapsed_clades=[{"tree": "lower:0", "node": doc.trees["lower:0"].root}],
    )
    with pytest.raises(ValueError, match="different clades"):
        validate_comparison_view(view, {"a": doc, "b": validate_document(other)})


def test_html_initial_collapse_keeps_expandable_geometry_and_pdf_view(tmp_path):
    doc = read_document(TOY)
    tree = doc.trees["lower:0"]
    view = {
        "schema_version": "0.1.0",
        "kind": "reconciliation_view",
        "layout": "paired",
        "collapsed_clades": [{"tree": tree.id, "node": tree.root}],
        "event_colors": {"duplication": "#123456"},
    }
    config = tmp_path / "view.json"
    config.write_text(json.dumps(view))
    output = tmp_path / "collapsed.html"
    main(
        [
            "draw",
            "--document",
            str(TOY),
            "--view-config",
            str(config),
            "--outfile",
            str(output),
        ]
    )
    text = output.read_text()
    svg = ET.fromstring(text[text.index('<svg id="scene"') : text.index("</svg>") + 6])
    owners = [
        json.loads(e.attrib["data-detail"])["owner"]
        for e in svg.iter()
        if "data-detail" in e.attrib
    ]
    assert [
        tree.id,
        next(nid for nid, n in tree.nodes.items() if n["label"] == "gA1"),
    ] in owners
    encoded = text.split('id="pdfbytes">', 1)[1].split("</script>", 1)[0]
    pdf_text = "".join(
        page.extract_text()
        for page in PdfReader(io.BytesIO(base64.b64decode(encoded))).pages
    )
    assert "collapsed" in pdf_text and "gA1" not in pdf_text
    assert "#123456" in text
    # A layout override is revalidated before producing an invalid saved view.
    bad = dict(
        view, collapsed_clades=[{"tree": "upper", "node": doc.trees["upper"].root}]
    )
    config.write_text(json.dumps(bad))
    with pytest.raises(ValueError, match="Upper context"):
        main(
            [
                "draw",
                "--document",
                str(TOY),
                "--view-config",
                str(config),
                "--view",
                "embedded",
                "--outfile",
                str(tmp_path / "bad.pdf"),
            ]
        )


def test_comparison_missing_history_still_honors_tree_display_state(tmp_path):
    from nwkit.reconciliation_import import import_tree_document

    path, record = comparison(tmp_path, ("missing",))
    doc = validate_document(import_tree_document("(one,two);"))
    tree = next(iter(doc.trees.values()))
    record["view"] = dict(
        record["view"],
        layout="paired",
        collapsed_clades=[{"tree": tree.id, "node": tree.root}],
    )
    scene = comparison_scene(record, {"0": doc}, width=720, font_size=8)
    assert not any(p.text in {"one", "two"} for p in scene.primitives)
    assert any("collapsed" in p.text for p in scene.primitives)


def test_recorded_event_difference_preserves_hidden_side_and_arity(tmp_path):
    path, record = comparison(tmp_path, ("available", "available"))
    first = read_document(TOY).record
    rec = first["reconciliations"][0]
    upper = {n["label"]: n["id"] for n in first["trees"][0]["nodes"]}
    rec["events"][0]["annotations"].update(
        observed_arity=1, marginalized_species=upper["A"]
    )
    other = copy.deepcopy(first)
    other["reconciliations"][0]["events"][0]["annotations"]["marginalized_species"] = (
        upper["B"]
    )
    for i, doc in enumerate((first, other)):
        source = tmp_path / f"side-{i}.json"
        source.write_text(json.dumps(doc))
        record["panels"][i].update(
            document=source.name,
            document_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        )
    path.write_text(json.dumps(record))
    loaded, docs, _ = read_comparison(path)
    scene = comparison_scene(loaded, docs, width=720, font_size=8)
    assert any("difference" in text and "+1/-1" in text for text in scene.notices)


def test_comparison_hash_validates_exact_snapshot_used_for_render_and_zip(
    tmp_path, monkeypatch
):
    path, record = comparison(tmp_path)
    source = (tmp_path / "source.json").resolve()
    original = source.read_bytes()
    read = Path.read_bytes
    reads = []

    def capture_once(self):
        raw = read(self)
        if self == source:
            reads.append(self)
            # Simulate an external file edit after the verified snapshot is read.
            self.write_bytes(b"invalid external edit")
        return raw

    monkeypatch.setattr(Path, "read_bytes", capture_once)
    loaded, docs, _ = read_comparison(path)
    assert len(reads) == 1
    assert docs["0"].source_bytes == original
    assert docs["0"].record == read_document(TOY).record
    assert loaded == record


def test_comparison_event_positions_retain_each_panel_even_with_reused_ids(tmp_path):
    path, record = comparison(tmp_path, ("available", "available"))
    loaded, docs, _ = read_comparison(path)
    scene = comparison_scene(loaded, docs, width=720, font_size=8)
    expected = {
        (panel, rec["id"], e["id"])
        for panel, doc in docs.items()
        for rec in doc.record["reconciliations"]
        for e in rec["events"]
    }
    assert set(scene.event_points) == expected


def test_nested_comparison_retains_expanded_label_corridor_width(tmp_path):
    doc = read_document(TOY).record
    middle = doc["trees"][1]
    bottom = copy.deepcopy(middle)
    bottom["id"] = "bottom"
    for n in bottom["nodes"]:
        if not n["children"]:
            n["label"] = "nested_" + n["label"]
    doc["trees"].append(bottom)
    doc["reconciliations"].append(
        {
            "id": "bottom_to_middle",
            "lower_tree": "bottom",
            "upper_tree": middle["id"],
            "history_kind": "unspecified",
            "selection": {},
            "model": {},
            "events": [
                {
                    "id": n["id"],
                    "type": "speciation" if n["children"] else "leaf",
                    "lower_node": n["id"],
                    "upper_node": n["id"],
                }
                for n in bottom["nodes"]
            ],
        }
    )
    validated = validate_document(doc)
    nested = build_scene(validated, width=720, font_size=8)
    assert nested.width > 720
    path, record = comparison(tmp_path, ("available",))
    composed = comparison_scene(record, {"0": validated}, width=720, font_size=8)
    assert composed.width >= nested.width


def test_x1_highlight_preserves_biological_event_colors():
    from nwkit.reconciliation_import import import_tree_document
    from nwkit.reconciliation_scene import EVENT_STYLE

    source = import_tree_document("((A,B),C);")["trees"][0]
    target = import_tree_document("((A,C),B);")["trees"][0]
    upper = import_tree_document("species;")["trees"][0]
    source["id"], target["id"], upper["id"] = "H", "G", "S"
    a = {n["label"]: n["id"] for n in source["nodes"] if not n["children"]}
    b = {n["label"]: n["id"] for n in target["nodes"] if not n["children"]}
    hi = next(
        n["id"] for n in source["nodes"] if set(n["children"]) == {a["A"], a["B"]}
    )
    gi = next(
        n["id"] for n in target["nodes"] if set(n["children"]) == {b["A"], b["C"]}
    )
    record = {
        "schema_version": "0.1.0",
        "kind": "reconciliation_document",
        "provenance": {"description": "Synthetic color fixture; no inference"},
        "trees": [source, target, upper],
        "reconciliations": [
            {
                "id": "rec",
                "lower_tree": "H",
                "upper_tree": "S",
                "history_kind": "unspecified",
                "model": {},
                "selection": {},
                "events": [
                    {
                        "id": n["id"],
                        "type": "duplication" if n["children"] else "leaf",
                        "lower_node": n["id"],
                        "upper_node": upper["root"],
                    }
                    for n in source["nodes"]
                ],
            }
        ],
        "transformations": [
            {
                "id": "X",
                "source_tree": "H",
                "target_tree": "G",
                "model_id": "X1",
                "blocks": [
                    {
                        "id": "x",
                        "source_nodes": [source["root"], hi],
                        "target_nodes": [target["root"], gi],
                        "resolution": 1,
                        "boundaries": [
                            {"source_node": a[x], "target_node": b[x]} for x in "ABC"
                        ],
                    }
                ],
            }
        ],
    }
    scene = build_scene(validate_document(record))
    assert any(p.kind == "line" and p.color == "#a03faa" for p in scene.primitives)
    assert all(
        p.color == EVENT_STYLE[p.text][0]
        for p in scene.primitives
        if p.kind == "event" and p.id.startswith("event:")
    )


def test_collapse_event_positions_match_only_displayed_events():
    doc = read_document(TOY)
    tree = doc.trees["lower:0"]
    scene = apply_view(
        build_scene(doc),
        doc,
        {"collapsed_clades": [{"tree": tree.id, "node": tree.root}]},
    )
    rec = doc.record["reconciliations"][0]
    assert set(scene.event_points) == {
        (rec["id"], e["id"]) for e in rec["events"] if e["lower_node"] == tree.root
    }


def test_lower_clade_selection_keeps_middle_context_event_coordinates():
    from nwkit.reconciliation_document import select_document_clade

    record = read_document(TOY).record
    bottom = copy.deepcopy(record["trees"][1])
    bottom["id"] = "bottom"
    record["trees"].append(bottom)
    record["reconciliations"].append(
        {
            "id": "bottom_to_middle",
            "lower_tree": "bottom",
            "upper_tree": "lower:0",
            "history_kind": "unspecified",
            "model": {},
            "selection": {},
            "events": [
                {
                    "id": n["id"],
                    "type": "speciation" if n["children"] else "leaf",
                    "lower_node": n["id"],
                    "upper_node": n["id"],
                }
                for n in bottom["nodes"]
            ],
        }
    )
    doc = validate_document(record)
    leaf = next(n["id"] for n in bottom["nodes"] if n["label"] == "gA1")
    selected = select_document_clade(doc, f"bottom/{leaf}")
    scene = build_scene(selected)
    assert all(
        ("history:0", e["id"]) in scene.event_points
        for e in record["reconciliations"][0]["events"]
    )
    assert ("bottom_to_middle", leaf) in scene.event_points
