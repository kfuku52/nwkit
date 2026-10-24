"""Meaningful geometry and real vector-output/rollback checks for N2."""

import copy
from pathlib import Path

import pytest
from pypdf import PdfReader

from nwkit.cli import main
from nwkit.reconciliation_document import (
    read_document,
    select_document_clade,
    validate_document,
)
from nwkit.reconciliation_render import write_scene
from nwkit.reconciliation_scene import Primitive, Scene, build_scene

ROOT = Path(__file__).resolve().parents[1]
TOY = ROOT / "examples/reconciliation/toy.recphyloxml"


def test_every_recorded_event_and_transfer_endpoint_is_in_scene():
    doc = read_document(TOY)
    scene = build_scene(doc)
    rec = doc.record["reconciliations"][0]
    assert set(scene.event_points) == {(rec["id"], e["id"]) for e in rec["events"]}
    donor = next(e for e in rec["events"] if e["type"] == "transfer_out")
    receiver = next(e for e in rec["events"] if e["type"] == "transfer_in")
    arrow = next(p for p in scene.primitives if p.kind == "transfer")
    assert arrow.points == [
        scene.event_points[rec["id"], donor["id"]],
        scene.event_points[rec["id"], receiver["id"]],
    ]
    assert "do not represent times" in scene.caption


def test_multiple_lower_trees_have_separate_lanes():
    doc = read_document(TOY)
    record = copy.deepcopy(doc.record)
    tree = copy.deepcopy(record["trees"][1])
    tree["id"] = "lower:second"
    rec = copy.deepcopy(record["reconciliations"][0])
    rec.update(id="second", lower_tree=tree["id"])
    record["trees"].append(tree)
    record["reconciliations"].append(rec)
    scene = build_scene(validate_document(record))
    eid = rec["events"][0]["id"]
    assert (
        scene.event_points["history:0", eid][1] != scene.event_points["second", eid][1]
    )


def test_child_order_does_not_change_scene():
    doc = read_document(TOY)
    record = copy.deepcopy(doc.record)
    for tree in record["trees"]:
        tree["nodes"].reverse()
        for node in tree["nodes"]:
            node["children"].reverse()
    assert (
        build_scene(validate_document(record)).event_points
        == build_scene(doc).event_points
    )


@pytest.mark.parametrize("extension", ["pdf", "svg"])
def test_cli_real_vector_outputs_and_normalized_document(
    tmp_path, extension, monkeypatch
):
    import nwkit.reconcile

    def forbidden(*args, **kwargs):
        raise AssertionError("Drawing must not infer reconciliation")

    monkeypatch.setattr(nwkit.reconcile, "build_reconciliation_table", forbidden)
    output, normalized = tmp_path / f"figure.{extension}", tmp_path / "document.json"
    main(
        [
            "draw",
            "--document",
            str(TOY),
            "--outfile",
            str(output),
            "--document-out",
            str(normalized),
        ]
    )
    assert read_document(normalized).record == read_document(TOY).record
    if extension == "svg":
        text = output.read_text()
        assert (
            "event:history:0:" in text
            and "gA1" in text
            and "transfer:history:0:" in text
        )
        assert "<image" not in text
    else:
        pdf = PdfReader(output)
        assert "gA1" in "".join(page.extract_text() for page in pdf.pages)
        fonts = pdf.pages[0]["/Resources"]["/Font"].values()
        assert any(
            "/FontFile2"
            in f.get_object()["/DescendantFonts"][0].get_object()["/FontDescriptor"]
            for f in fonts
            if "/DescendantFonts" in f.get_object()
        )


def test_tiled_pdf_repeats_caption_and_marks_cross_page_transfer(tmp_path):
    scene = Scene(
        720,
        1000,
        [
            Primitive(
                "transfer", "transfer:test", [(100, 100), (400, 900)], color="#9c4caa"
            )
        ],
    )
    path = tmp_path / "tiled.pdf"
    write_scene(scene, path, image_format="pdf", pdf_pages="tiled")
    pages = PdfReader(path).pages
    assert len(pages) == 3
    assert all(
        "event positions do not represent times" in p.extract_text() for p in pages
    )
    assert any("Transfer continues: p" in p.extract_text() for p in pages)


@pytest.mark.parametrize(
    "points,crossings",
    [([(100, 100), (400, 900)], 1), ([(100, 400), (650, 400)], 2)],
)
def test_page_labels_follow_actual_rendered_curve_and_branch(points, crossings):
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.path import Path as ArtistPath

    from nwkit.reconciliation_render import _figure

    primitives = [
        Primitive("transfer", "curved", points),
        Primitive("line", "bent", [(60, 50), (60, 600), (450, 600)]),
    ]
    scene = Scene(720, 1000, primitives)
    edge = 455
    figure = _figure(scene, primitives, "QA", 8, 842, 595, 0, edge, 0, 3, edge, False)
    try:
        figure.canvas.draw()
        ax = figure.axes[0]
        curve = ax.patches[0].get_path()
        samples = np.concatenate(
            [
                segment(np.linspace(0, 1, 20001))
                for segment, code in curve.iter_bezier()
                if code == ArtistPath.CURVE3
            ]
        )
        boundary_x = samples[np.abs(samples[:, 1] - edge) < 0.1, 0]
        labels = [t for t in ax.texts if t.get_text() == "Transfer continues: p2"]
        assert len(labels) == crossings
        for label in labels:
            assert np.min(np.abs(boundary_x - label.get_position()[0])) < 2
        branch = next(t for t in ax.texts if t.get_text() == "Branch continues: p2")
        assert branch.get_position()[0] == 60
    finally:
        plt.close(figure)


def test_clade_selection_retains_upper_context_and_marks_display_only():
    doc = read_document(TOY)
    tree = doc.trees["lower:0"]
    nid = next(nid for nid, n in tree.nodes.items() if n["label"] == "dup")
    selected = select_document_clade(doc, f"lower:0/{nid}")
    assert selected.trees["upper"].record == doc.trees["upper"].record
    assert selected.trees["lower:0"].record == doc.trees["lower:0"].record
    assert (
        selected.record["provenance"]["display_selection"]["inference_recomputed"]
        is False
    )
    expected = [
        e
        for e in selected.record["reconciliations"][0]["events"]
        if tree.ancestral(nid, e["lower_node"])
    ]
    assert len(build_scene(selected).event_points) == len(expected)


def test_transferred_clade_selection_retains_incoming_context():
    doc = read_document(TOY)
    tree = doc.trees["lower:0"]
    nid = next(nid for nid, node in tree.nodes.items() if node["label"] == "gA3")
    selected = select_document_clade(doc, f"lower:0/{nid}")
    assert selected.record["trees"] == doc.record["trees"]
    scene = build_scene(selected)
    assert len(scene.event_points) == 2
    assert len([p for p in scene.primitives if p.kind == "transfer"]) == 1


def test_input_aliases_and_unsupported_options_leave_files_intact(tmp_path):
    source = tmp_path / "source.json"
    source.write_text((ROOT / "examples/reconciliation/toy.json").read_text())
    before = source.read_bytes()
    with pytest.raises(ValueError, match="overwrite input"):
        main(["draw", "--document", str(source), "--outfile", str(source)])
    assert source.read_bytes() == before
    with pytest.raises(ValueError, match="layout.*not supported"):
        main(
            [
                "draw",
                "--document",
                str(source),
                "--layout",
                "radial",
                "--outfile",
                str(tmp_path / "f.svg"),
            ]
        )
    with pytest.raises(ValueError, match="require --document"):
        main(["draw", "--pdf-pages", "tiled"])


def test_failed_related_output_preserves_previous_bundle(tmp_path, monkeypatch):
    import nwkit.output_transaction as transaction

    figure, document = tmp_path / "old.svg", tmp_path / "old.json"
    figure.write_text("previous figure")
    document.write_text("previous document")
    original = transaction._StagedOutputs.write_text

    def fail(self, target, writer):
        if target == str(document):
            raise OSError("injected document write failure")
        return original(self, target, writer)

    monkeypatch.setattr(transaction._StagedOutputs, "write_text", fail)
    with pytest.raises(OSError, match="injected"):
        main(
            [
                "draw",
                "--document",
                str(TOY),
                "--outfile",
                str(figure),
                "--document-out",
                str(document),
            ]
        )
    assert figure.read_text() == "previous figure"
    assert document.read_text() == "previous document"


def test_audit_lists_document_and_does_not_consume_stdin(tmp_path, monkeypatch):
    import sys

    class ForbiddenInput:
        def read(self, *args):
            raise AssertionError("Must not consume implicit STDIN")

    monkeypatch.setattr(sys, "stdin", ForbiddenInput())
    audit = tmp_path / "audit.json"
    main(
        [
            "draw",
            "--document",
            str(TOY),
            "--outfile",
            str(tmp_path / "audit.svg"),
            "--audit",
            str(audit),
        ]
    )
    assert str(TOY.resolve()) in audit.read_text()
    assert '"argument": "document"' in audit.read_text()


def test_expectations_and_missing_placements_are_explicit():
    record = copy.deepcopy(read_document(TOY).record)
    record["reconciliations"][0]["history_kind"] = "expectation"
    with pytest.raises(ValueError, match="Expectation"):
        build_scene(validate_document(record))
    record["reconciliations"][0]["history_kind"] = "unspecified"
    record["reconciliations"][0]["events"][0]["upper_node"] = None
    with pytest.raises(ValueError, match="explicit placement"):
        build_scene(validate_document(record))


def test_long_labels_have_distinct_copy_lanes():
    record = copy.deepcopy(read_document(TOY).record)
    lower = record["trees"][1]
    for n in lower["nodes"]:
        if not n["children"] and n["label"].startswith("gA"):
            n["label"] += "_" + "long_identifier_" * 5
    scene = build_scene(validate_document(record))
    labels = [
        p for p in scene.primitives if p.kind == "text" and p.text.startswith("gA")
    ]
    ys = sorted(p.points[0][1] for p in labels)
    assert min(b - a for a, b in zip(ys, ys[1:], strict=False)) >= 48


def test_thousand_species_and_genes_geometry_is_finite_and_complete():
    import math

    from nwkit.reconciliation_import import RawNode, empty_document, flatten_tree

    pairs = [(RawNode(f"S{i:04d}"), RawNode(f"g{i:04d}")) for i in range(1000)]
    correspondence = {lower: upper for upper, lower in pairs}
    while len(pairs) > 1:
        level = []
        for i in range(0, len(pairs), 2):
            if i + 1 == len(pairs):
                level.append(pairs[i])
            else:
                upper = RawNode(children=[pairs[i][0], pairs[i + 1][0]])
                lower = RawNode(children=[pairs[i][1], pairs[i + 1][1]])
                correspondence[lower] = upper
                level.append((upper, lower))
        pairs = level
    upper, uids = flatten_tree(pairs[0][0], "S", "upper")
    lower, gids = flatten_tree(pairs[0][1], "G", "lower")
    record = empty_document("Synthetic 1000-tip geometry fixture")
    record["trees"] = [upper, lower]
    record["reconciliations"] = [
        {
            "id": "history",
            "lower_tree": "G",
            "upper_tree": "S",
            "history_kind": "unspecified",
            "selection": {},
            "model": {"source": "synthetic"},
            "events": [
                {
                    "id": nid,
                    "type": "speciation" if node.children else "leaf",
                    "lower_node": nid,
                    "upper_node": uids[correspondence[node]],
                }
                for node, nid in gids.items()
            ],
        }
    ]
    scene = build_scene(validate_document(record))
    assert len(scene.event_points) == 1999
    assert all(
        math.isfinite(value)
        for p in scene.primitives
        for xy in p.points
        for value in xy
    )


def test_saved_tsv_keeps_loss_counts_as_summaries(tmp_path):
    gene, species, events = (
        tmp_path / "gene.nwk",
        tmp_path / "species.nwk",
        tmp_path / "events.tsv",
    )
    gene.write_text("((A_1:1,B_1:1):1,A_2:1);")
    species.write_text("(A:1,B:1);")
    mapping = tmp_path / "mapping.tsv"
    mapping.write_text("leaf_name\tspecies_label\nA_1\tA\nB_1\tB\nA_2\tA\n")
    main(
        [
            "reconcile",
            "--infile",
            str(gene),
            "--species-tree",
            str(species),
            "--species-map-tsv",
            str(mapping),
            "--outfile",
            str(events),
        ]
    )
    normalized = tmp_path / "saved.json"
    main(
        [
            "draw",
            "--infile",
            str(gene),
            "--species-tree",
            str(species),
            "--reconciliation",
            str(events),
            "--view",
            "embedded",
            "--outfile",
            str(tmp_path / "saved.svg"),
            "--document-out",
            str(normalized),
        ]
    )
    rec = read_document(normalized).record["reconciliations"][0]
    assert rec["summary"]["implied_losses"]
    assert all(e["type"] != "loss" for e in rec["events"])
