"""Mixed DTL/X histories, verified loss pruning, default views and export fidelity."""

import copy
import json

import pytest

from nwkit.cli import main
from nwkit.reconciliation_document import validate_document
from nwkit.reconciliation_import import RawNode, empty_document, flatten_tree
from nwkit.reconciliation_projection import LOSS_PROJECTION, observed_context
from nwkit.reconciliation_scene import EVENT_STYLE, build_scene
from nwkit.reconciliation_view import default_view


def mixed_document():
    a, b, c, foreign = [RawNode(label) for label in ("a", "b", "c", "b_from_C")]
    sa, sb, sc = [RawNode(label) for label in "ABC"]
    sab = RawNode(children=[sa, sb])
    sroot = RawNode(children=[sab, sc])
    species, sids = flatten_tree(sroot, "S", "species")
    ab = RawNode(children=[a, b])
    lost = RawNode(annotations={"observed": False})
    dup = RawNode(children=[ab, lost])
    transfer = RawNode(children=[c, foreign])
    full_root = RawNode(children=[dup, transfer])
    hroot = RawNode(children=[ab, transfer])
    inner = RawNode(children=[a, transfer])
    groot = RawNode(children=[inner, b])
    full, fids = flatten_tree(full_root, "F", "complete_history")
    source, hids = flatten_tree(hroot, "H", "latent")
    target, gids = flatten_tree(groot, "G", "observed")
    events = []
    for raw, kind, upper in (
        (full_root, "origin", sroot),
        (full_root, "speciation", sroot),
        (dup, "duplication", sab),
        (lost, "loss", sab),
        (ab, "speciation", sab),
        (transfer, "transfer", sc),
        (a, "leaf", sa),
        (b, "leaf", sb),
        (c, "leaf", sc),
        (foreign, "leaf", sb),
    ):
        event = {
            "id": f"e{len(events)}",
            "type": kind,
            "lower_node": fids[raw],
            "upper_node": sids[upper],
        }
        if kind == "transfer":
            event.update(
                donor_node=sids[sc],
                recipient_node=sids[sb],
                recipient_lower_node=fids[foreign],
            )
        events.append(event)
    record = empty_document("Synthetic mixed DTLX geometry fixture")
    record["trees"] = [species, full, source, target]
    record["reconciliations"] = [
        {
            "id": "R",
            "lower_tree": "F",
            "upper_tree": "S",
            "history_kind": "unspecified",
            "selection": {},
            "model": {},
            "events": events,
        }
    ]
    record["transformations"] = [
        {
            "id": "projection",
            "source_tree": "F",
            "target_tree": "H",
            "model_id": LOSS_PROJECTION,
            "blocks": [],
        },
        {
            "id": "X",
            "source_tree": "H",
            "target_tree": "G",
            "model_id": "X1",
            "blocks": [
                {
                    "id": "x",
                    "source_nodes": [hids[hroot], hids[ab]],
                    "target_nodes": [gids[groot], gids[inner]],
                    "resolution": 1,
                    "boundaries": [
                        {"source_node": hids[n], "target_node": gids[n]}
                        for n in (a, b, transfer)
                    ],
                }
            ],
        },
    ]
    return validate_document(record)


def accounted_events(scene):
    return set(scene.event_points) | {
        key for p in scene.primitives for key in p.source_event_ids
    }


def test_default_mixed_x1_view_keeps_exact_g_and_every_recorded_event():
    doc = mixed_document()
    original = copy.deepcopy(doc.record)
    view = default_view(doc)
    assert view["schema_version"] == "0.2.0" and view["layout"] == "observed_embedded"
    context = observed_context(doc)
    assert len(context.source_to_history) == len(doc.trees["H"].nodes)
    scene = build_scene(doc)
    assert {
        p.owner[1] for p in scene.primitives if p.id.startswith("observed-node:")
    } == set(doc.trees["G"].nodes)
    events = [p for p in scene.primitives if p.id.startswith("event:R:")]
    assert accounted_events(scene) == {
        ("R", e["id"]) for e in doc.record["reconciliations"][0]["events"]
    }
    assert all(p.color == EVENT_STYLE[p.text][0] for p in events)
    assert len(scene.event_points) == len(events)
    assert sum(p.kind == "upper" for p in scene.primitives) == len(doc.trees["S"].nodes)
    assert not any(p.id.startswith("panel:") for p in scene.primitives)
    assert any(p.kind == "transfer" for p in scene.primitives)
    assert not any(p.kind == "event" and p.text == "loss" for p in scene.primitives)
    assert not any(
        p.tree == "G" and p.text in {"duplication", "transfer", "loss", "origin"}
        for p in scene.primitives
    )
    assert doc.record == original


def test_unified_routes_share_observed_tips_and_retain_recorded_species_placements():
    doc = mixed_document()
    scene = build_scene(doc)
    history, species = doc.trees["F"], doc.trees["S"]
    tips = {
        doc.trees["G"].nodes[p.owner[1]]["label"]: p.points[0]
        for p in scene.primitives
        if p.id.startswith("observed-node:")
        and not doc.trees["G"].nodes[p.owner[1]]["children"]
    }
    pipes = {p.owner[1]: p for p in scene.primitives if p.kind == "upper"}
    for event in doc.record["reconciliations"][0]["events"]:
        if ("R", event["id"]) not in scene.event_points:
            continue
        point = scene.event_points["R", event["id"]]
        if event["type"] == "leaf":
            assert point == tips[history.nodes[event["lower_node"]]["label"]]
        if event["type"] in {"duplication", "loss", "transfer"}:
            branch = event["upper_node"]
            end = pipes[branch].points[-1]
            start_x = (
                pipes[species.parents[branch]].points[-1][0]
                if branch in species.parents
                else pipes[branch].points[0][0]
            )
            assert start_x <= point[0] <= end[0]
            assert abs(point[1] - end[1]) < pipes[branch].width / 2
    routes = [p for p in scene.primitives if p.kind == "line" and p.tree == "F"]
    assert not routes
    assert not any(p.id.startswith("context-") for p in scene.primitives)
    assert not any(p.dashed for p in scene.primitives)
    assert all(
        not p.dashed for p in scene.primitives if p.id.startswith("observed-edge:")
    )


def test_unified_svg_uses_square_species_junctions_and_only_solid_gene_routes():
    from xml.etree import ElementTree as ET

    from nwkit.reconciliation_html import _svg
    from tests.test_reconciliation_forks import compatible_duplication_document

    scene = build_scene(compatible_duplication_document())
    svg = ET.fromstring(_svg(scene, 8))
    species = []
    extinct = []
    for element in svg.iter():
        detail = (
            json.loads(element.attrib["data-detail"])
            if "data-detail" in element.attrib
            else {}
        )
        if detail.get("kind") == "upper":
            species.append(element)
            assert element.attrib["stroke-linecap"] == "butt"
            assert element.attrib["stroke-linejoin"] == "miter"
        if detail.get("kind") == "line" and detail.get("tree") == "F":
            extinct.append(element)
            assert "stroke-dasharray" not in element.attrib
            assert detail["id"].startswith("extinct:")
    assert species and extinct


def test_vector_renderer_uses_square_corridors_and_solid_extinct_sides():
    import matplotlib.pyplot as plt

    from nwkit.reconciliation_render import _primitive
    from tests.test_reconciliation_forks import compatible_duplication_document

    scene = build_scene(compatible_duplication_document())
    upper = next(p for p in scene.primitives if p.kind == "upper")
    extinct = next(p for p in scene.primitives if p.id.startswith("extinct:"))
    fig, ax = plt.subplots()
    try:
        _primitive(ax, upper, 8)
        assert len(ax.lines) == 2
        assert all(line.get_solid_capstyle() == "butt" for line in ax.lines)
        assert all(line.get_solid_joinstyle() == "miter" for line in ax.lines)
        _primitive(ax, extinct, 8)
        assert ax.lines[-1].get_linestyle() == "-"
    finally:
        plt.close(fig)


def test_species_tip_corridors_do_not_inherit_the_crowded_ancestral_width():
    doc = mixed_document()
    scene = build_scene(doc)
    species = doc.trees["S"]
    widths = {
        p.owner[1]: p.width
        for p in scene.primitives
        if p.id.startswith("observed-upper:")
    }
    c = next(n for n, node in species.nodes.items() if node["label"] == "C")
    sibling = next(n for n in species.nodes[species.root]["children"] if n != c)
    assert widths[c] == widths[sibling]
    a = next(n for n, node in species.nodes.items() if node["label"] == "A")
    b = next(n for n, node in species.nodes.items() if node["label"] == "B")
    # G has two lineages in AB and B; the discarded H comparison adds no lanes.
    assert widths[a] == widths[b] <= widths[c] < widths[species.root]


@pytest.mark.parametrize("font_size", [8, 16])
@pytest.mark.parametrize("reverse_species", [False, True])
def test_species_split_widths_match_without_squeezing_recorded_lanes(
    font_size, reverse_species
):
    from nwkit.reconciliation_overlay import history_layout

    doc = mixed_document()
    species = doc.trees["S"]
    tips = [n for n in species.order if not species.nodes[n]["children"]]
    view = default_view(doc)
    if reverse_species:
        view["leaf_order"] = {"S": list(reversed(tips))}
    scene = build_scene(doc, font_size=font_size, view=view)
    corridors = {p.owner[1]: p for p in scene.primitives if p.kind == "upper"}
    recorded = history_layout(observed_context(doc), max(12, font_size * 1.5))
    for branch, node in species.nodes.items():
        assert corridors[branch].width >= recorded.widths.get(branch, 0)
        children = node["children"]
        if children:
            assert len({corridors[child].width for child in children}) == 1
            # The two vertical segments meet at the same internal-node center.
            assert len({corridors[child].points[0] for child in children}) == 1
    leaf_events = [
        e for e in doc.record["reconciliations"][0]["events"] if e["type"] == "leaf"
    ]
    rec = doc.record["reconciliations"][0]
    for event in leaf_events:
        x, y = scene.event_points[rec["id"], event["id"]]
        corridor = corridors[event["upper_node"]]
        assert abs(y - corridor.points[-1][1]) <= corridor.width / 2 - 4
        assert corridor.points[0][0] < x <= corridor.points[-1][0]


@pytest.mark.parametrize(
    "damage", ["topology", "map", "unknown_projection", "missing_loss", "internal_loss"]
)
def test_invalid_projection_is_rejected_before_display(damage):
    record = copy.deepcopy(mixed_document().record)
    rec, projection = record["reconciliations"][0], record["transformations"][0]
    full = record["trees"][1]
    if damage == "topology":
        for node in full["nodes"]:
            if node["label"] in {"a", "c"}:
                node["label"] = {"a": "c", "c": "a"}[node["label"]]
    elif damage == "map":
        projection["annotations"] = {"retained_nodes": {"invented": "invented"}}
    elif damage == "unknown_projection":
        projection["model_id"] = "unknown"
    elif damage == "missing_loss":
        rec["events"] = [e for e in rec["events"] if e["type"] != "loss"]
    else:
        next(e for e in rec["events"] if e["type"] == "loss")["lower_node"] = full[
            "root"
        ]
    doc = validate_document(record)
    before = copy.deepcopy(doc.record)
    with pytest.raises(
        ValueError, match="projection|unobserved history|Unobserved history"
    ):
        build_scene(doc)
    assert doc.record == before


def test_cli_default_and_explicit_saved_legacy_view_both_round_trip(tmp_path):
    doc = mixed_document()
    path = tmp_path / "source.json"
    path.write_text(json.dumps(doc.record))
    saved, normalized = tmp_path / "saved.json", tmp_path / "normalized.json"
    for extension in ("svg", "pdf", "html"):
        main(
            [
                "draw",
                "--document",
                str(path),
                "--outfile",
                str(tmp_path / f"figure.{extension}"),
                "--view-out",
                str(saved),
                "--document-out",
                str(normalized),
            ]
        )
        assert json.loads(saved.read_text())["layout"] == "observed_embedded"
        assert json.loads(normalized.read_text()) == doc.record
    assert "observed-node:G:" in (tmp_path / "figure.svg").read_text()
    saved.write_text(
        json.dumps(
            {
                "schema_version": "0.1.0",
                "kind": "reconciliation_view",
                "layout": "topology_embedded",
            }
        )
    )
    main(
        [
            "draw",
            "--document",
            str(path),
            "--view-config",
            str(saved),
            "--outfile",
            str(tmp_path / "legacy.svg"),
        ]
    )
    assert "observed-node:G:" not in (tmp_path / "legacy.svg").read_text()
    main(
        [
            "draw",
            "--document",
            str(path),
            "--view",
            "paired",
            "--outfile",
            str(tmp_path / "paired.svg"),
        ]
    )
    assert "observed-node:G:" not in (tmp_path / "paired.svg").read_text()
