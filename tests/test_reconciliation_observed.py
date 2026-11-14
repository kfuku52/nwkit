"""Observed topology, schematic placement, source fidelity and export contracts."""

import copy
import json
from pathlib import Path

import jsonschema
import pytest
from pypdf import PdfReader

from nwkit.cli import main
from nwkit.reconciliation_document import validate_document
from nwkit.reconciliation_import import import_tree_document
from nwkit.reconciliation_observed import (
    _context,
    _crossing_gaps,
    _observed_geometry,
    _placements,
    observed_topology_scene,
)
from nwkit.reconciliation_ordering import MAX_LAYOUTS, choose_orientation, path_contacts
from nwkit.reconciliation_scene import build_scene
from nwkit.reconciliation_view import validate_view

VIEW = {
    "schema_version": "0.2.0",
    "kind": "reconciliation_view",
    "layout": "observed_embedded",
}


def fixture_document(resolution=1, multiple=False):
    h = "(((a,b),c),((d,e),f));" if multiple else "((a,b),c);"
    g = (
        "(((a,c),b),((d,f),e));"
        if multiple
        else ("((a,c),b);" if resolution == 1 else "((b,c),a);")
    )
    trees = [import_tree_document(n)["trees"][0] for n in (h.upper(), h, g)]
    for tree, tid in zip(trees, ("S", "H", "G"), strict=True):
        tree["id"] = tid

    def clades(tree):
        nodes = {n["id"]: n for n in tree["nodes"]}

        def labels(n):
            return (
                frozenset().union(*(labels(nodes[c]) for c in n["children"]))
                if n["children"]
                else frozenset([n["label"].lower()])
            )

        return {labels(n): n["id"] for n in nodes.values()}

    s, hmap, gmap = [clades(t) for t in trees]
    events = [
        {
            "id": "origin",
            "type": "origin",
            "lower_node": trees[1]["root"],
            "upper_node": trees[0]["root"],
        }
    ]
    for clade, nid in hmap.items():
        events.append(
            {
                "id": f"event:{nid}",
                "type": "leaf" if len(clade) == 1 else "speciation",
                "lower_node": nid,
                "upper_node": s[clade],
            }
        )
    blocks = []
    for labels in ["abc", "def"] if multiple else ["abc"]:
        a, b, c = [frozenset([x]) for x in labels]
        pair = a | c if resolution == 1 else b | c
        blocks.append(
            {
                "id": labels,
                "source_nodes": [hmap[a | b | c], hmap[a | b]],
                "target_nodes": [gmap[a | b | c], gmap[pair]],
                "resolution": resolution,
                "boundaries": [
                    {"source_node": hmap[x], "target_node": gmap[x]} for x in [a, b, c]
                ],
            }
        )
    return validate_document(
        {
            "schema_version": "0.1.0",
            "kind": "reconciliation_document",
            "provenance": {"description": "Synthetic geometry fixture; no inference"},
            "trees": trees,
            "reconciliations": [
                {
                    "id": "R",
                    "upper_tree": "S",
                    "lower_tree": "H",
                    "history_kind": "unspecified",
                    "selection": {},
                    "model": {},
                    "events": events,
                }
            ],
            "transformations": [
                {
                    "id": "X",
                    "source_tree": "H",
                    "target_tree": "G",
                    "model_id": "X1",
                    "blocks": blocks,
                }
            ],
        }
    )


@pytest.mark.parametrize("resolution,multiple", [(1, False), (2, False), (1, True)])
def test_observed_edges_retain_exact_g_topology_and_older_parents(resolution, multiple):
    doc = fixture_document(resolution, multiple)
    original = copy.deepcopy(doc.record)
    scene = observed_topology_scene(doc)
    target = doc.trees["G"]
    nodes = {
        p.owner[1]: p.points[0]
        for p in scene.primitives
        if p.id.startswith("observed-node:")
    }
    edges = [p for p in scene.primitives if p.id.startswith("observed-edge:")]
    assert {p.owner[1] for p in edges} == set(target.nodes)
    for child, parent in target.parents.items():
        fragments = [p for p in edges if p.owner[1] == child]
        assert fragments[0].points[0] == nodes[parent]
        assert fragments[-1].points[-1] == nodes[child]
        assert nodes[parent][0] < nodes[child][0]
    root = [p for p in edges if p.owner[1] == target.root]
    species_stem = next(
        p for p in scene.primitives if p.id == f"observed-upper:S:{doc.trees['S'].root}"
    )
    assert root[0].points[0][0] == species_stem.points[0][0]
    assert (
        abs(root[0].points[0][1] - species_stem.points[0][1]) < species_stem.width / 2
    )
    assert root[-1].points[-1] == nodes[target.root]
    assert all(p.color == "#287d8e" for p in root)
    from tests.test_reconciliation_projection import accounted_events

    assert accounted_events(scene) == {
        ("R", e["id"]) for e in original["reconciliations"][0]["events"]
    }
    assert all(
        p.text == "gene_node"
        for p in scene.primitives
        if p.kind == "event" and p.tree == "G"
    )
    assert doc.record == original
    full = build_scene(doc, view=VIEW)
    assert accounted_events(full) == {
        ("R", e["id"]) for e in original["reconciliations"][0]["events"]
    }


def test_discordant_join_is_in_ancestral_species_corridor():
    scene = observed_topology_scene(fixture_document())
    species_root = next(
        p
        for p in scene.primitives
        if p.id.startswith("observed-upper:") and len(p.points) == 2
    )
    joins = [
        p.points[0][0]
        for p in scene.primitives
        if p.kind == "event" and p.color == "#a03faa"
    ]
    assert len(joins) == 2 and max(joins) < species_root.points[-1][0]


def test_crossing_is_a_gap_instead_of_an_unmarked_join():
    paths = {"a": [(0, 10), (30, 10)], "b": [(15, 0), (15, 20)]}
    fragments = list(_crossing_gaps(paths))
    assert [p for n, p in fragments if n == "a"] == [
        [(0, 10), (12, 10)],
        [(18, 10), (30, 10)],
    ]
    assert [p for n, p in fragments if n == "b"] == [paths["b"]]


def test_stem_split_rotates_and_removes_two_unnecessary_crossings():
    doc = fixture_document()
    original = copy.deepcopy(doc.record)
    _, _, upper, target, mapping = _context(doc)
    placed = _placements(upper, target, mapping)
    xy, _, _, _, positions, paths, search = _observed_geometry(
        upper, target, placed, 12, 380, None, ()
    )
    assert search["previous_score"][:2] == (0, 3)
    assert search["score"][:3] == (0, 1, 0)
    assert search["evaluated"] == 4 and search["exhaustive"]
    assert search["flipped"] == (target.root,)
    b = next(nid for nid, node in target.nodes.items() if node["label"] == "b")
    other = next(nid for nid in target.nodes[target.root]["children"] if nid != b)
    # B now leaves above the AC branch in the stem; the one crossing is after S's root.
    assert paths[b][1][1] < positions[target.root][1] < positions[other][1]
    assert all(x > xy[upper.root][0] for _, _, x, _ in path_contacts(paths)[0])
    assert doc.record == original


@pytest.mark.parametrize("width,font_size", [(380, 8), (720, 8), (720, 16)])
@pytest.mark.parametrize("reverse_species", [False, True])
def test_join_pitch_and_centered_species_junction_bundles(
    width, font_size, reverse_species
):
    doc = fixture_document()
    _, _, upper, target, mapping = _context(doc)
    tips = sorted(
        (n for n, node in upper.nodes.items() if not node["children"]),
        key=lambda n: upper.nodes[n]["label"],
    )
    view = {"leaf_order": {upper.id: list(reversed(tips))}} if reverse_species else None
    xy, unit, _, _, positions, paths, search = _observed_geometry(
        upper,
        target,
        _placements(upper, target, mapping),
        max(12, font_size * 1.5),
        width,
        view,
        (),
    )
    pitch = search["horizontal_pitch"]
    joins = sorted(
        positions[n][0] for n, node in target.nodes.items() if node["children"]
    )
    assert joins[1] - joins[0] == pytest.approx(pitch)
    ab = next(n for n in upper.nodes if upper.nodes[n]["children"] and n != upper.root)
    mid_y = (xy[upper.root][1] + xy[ab][1]) / 2
    upward = sorted(
        {
            a[0]
            for points in paths.values()
            for a, b in zip(points[:-1], points[1:], strict=True)
            if a[0] == b[0]
            and abs(a[0] - xy[upper.root][0]) < unit / 2
            and min(a[1], b[1]) < mid_y < max(a[1], b[1])
        }
    )
    assert upward == pytest.approx(
        [xy[upper.root][0] - pitch / 2, xy[upper.root][0] + pitch / 2]
    )
    assert search["score"][:3] == (0, 1, 0)


def test_gene_stem_reaches_species_root_when_gene_lca_is_below_it():
    record = copy.deepcopy(fixture_document().record)
    species = import_tree_document("(((A,B),C),D);")["trees"][0]
    species["id"] = "S"
    record["trees"][0] = species
    document = validate_document(record)
    _, _, upper, target, mapping = _context(document)
    placed = _placements(upper, target, mapping)
    assert placed[target.root] != upper.root
    xy, unit, _, _, positions, paths, _ = _observed_geometry(
        upper, target, placed, 12, 720, None, ()
    )
    assert paths[target.root][0] == (xy[upper.root][0] - unit, xy[upper.root][1])
    assert paths[target.root][-1] == positions[target.root]
    assert len(paths[target.root]) > 2
    assert document.record == record


@pytest.mark.parametrize(
    "resolution,multiple,reverse_species",
    [(1, False, False), (2, False, False), (1, True, False), (1, False, True)],
)
def test_crossing_search_never_worsens_previous_geometry(
    resolution, multiple, reverse_species
):
    doc = fixture_document(resolution, multiple)
    _, _, upper, target, mapping = _context(doc)
    view = None
    if reverse_species:
        tips = [nid for nid, node in upper.nodes.items() if not node["children"]]
        view = {"leaf_order": {upper.id: list(reversed(tips))}}
    *_, search = _observed_geometry(
        upper, target, _placements(upper, target, mapping), 12, 380, view, ()
    )
    assert search["score"] <= search["previous_score"]


def test_contact_counts_keep_shared_nodes_but_detect_false_joins_and_overlaps():
    x = 1 / 3
    assert path_contacts({"a": [(x, 1), (x, 3)], "b": [(x, -1), (x, 1)]}) == (
        set(),
        set(),
    )
    contacts, overlaps = path_contacts(
        {"a": [(-10, 0), (10, 0)], "b": [(0, 0), (0, 10)]}
    )
    assert len(contacts) == 1 and not overlaps
    contacts, overlaps = path_contacts({"a": [(0, 0), (10, 0)], "b": [(5, 0), (15, 0)]})
    assert not contacts and len(overlaps) == 1


def test_large_orientation_search_has_a_fixed_candidate_limit():
    newick = "g000"
    for i in range(1, 100):
        newick = f"({newick},g{i:03d})"
    doc = validate_document(import_tree_document(newick + ";"))
    tree = next(iter(doc.trees.values()))
    original = copy.deepcopy(tree.record)
    seen = []

    def evaluate(order):
        seen.append(tuple(order))
        return (0, int(order[1] == seen[0][1]), 0, 0), tuple(order)

    chosen, search = choose_orientation(tree, evaluate)
    assert not search["exhaustive"] and len(seen) <= MAX_LAYOUTS
    assert search["score"][1] == 0 and chosen[1] != seen[0][1]
    assert tree.record == original


def test_branch_storage_order_does_not_change_crossing_reduction():
    doc = fixture_document(multiple=True)
    rotated = copy.deepcopy(doc.record)
    for tree in rotated["trees"]:
        for node in tree["nodes"]:
            node["children"].reverse()
    other = validate_document(rotated)

    def geometry(document):
        result = observed_topology_scene(document)
        return {p.id: p.points for p in result.primitives}

    assert geometry(doc) == geometry(other)


@pytest.mark.parametrize("kind", ["empty_observation", "unresolved"])
def test_unresolved_history_is_explicitly_unsupported_without_changing_source(kind):
    doc = fixture_document()
    doc.record["reconciliations"][0]["events"][1]["type"] = kind
    original = copy.deepcopy(doc.record)
    with pytest.raises(ValueError, match="resolved recorded histories"):
        build_scene(doc, view=VIEW)
    assert doc.record == original


def test_view_version_and_upper_collapse_contract():
    schema = json.loads(
        (
            Path(__file__).parents[1]
            / "nwkit/data_schema/reconciliation_view.schema.json"
        ).read_text()
    )
    validate_view(VIEW, fixture_document())
    jsonschema.validate(VIEW, schema)
    old = dict(VIEW, schema_version="0.1.0")
    with pytest.raises(ValueError, match="version/layout"):
        validate_view(old)
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(old, schema)
    with pytest.raises(ValueError, match="Upper context"):
        validate_view(
            dict(
                VIEW,
                collapsed_clades=[
                    {"tree": "S", "node": fixture_document().trees["S"].root}
                ],
            ),
            fixture_document(),
        )


def test_independent_g_tip_order_is_not_silently_ignored():
    doc = fixture_document()
    leaves = [nid for nid, node in doc.trees["G"].nodes.items() if not node["children"]]
    with pytest.raises(ValueError, match="set leaf_order on S"):
        build_scene(doc, view=dict(VIEW, leaf_order={"G": leaves}))


def test_unsupported_x1_cli_preserves_previous_output_bundle(tmp_path):
    doc = fixture_document()
    doc.record["reconciliations"][0]["events"][1]["type"] = "unresolved"
    source = tmp_path / "history.json"
    source.write_text(json.dumps(doc.record))
    paths = [
        tmp_path / name for name in ("old.pdf", "old.view.json", "old.document.json")
    ]
    for path in paths:
        path.write_bytes(b"previous output")
    with pytest.raises(ValueError, match="resolved recorded histories"):
        main(
            [
                "draw",
                "--document",
                str(source),
                "--view",
                "observed-embedded",
                "--outfile",
                str(paths[0]),
                "--view-out",
                str(paths[1]),
                "--document-out",
                str(paths[2]),
            ]
        )
    assert all(path.read_bytes() == b"previous output" for path in paths)


def test_cli_exports_and_saved_view_round_trip_without_scientific_changes(tmp_path):
    doc = fixture_document()
    source = tmp_path / "history.json"
    source.write_text(json.dumps(doc.record))
    view_path, normalized = tmp_path / "view.json", tmp_path / "normalized.json"
    for extension in ("svg", "pdf", "html"):
        main(
            [
                "draw",
                "--document",
                str(source),
                "--view",
                "observed-embedded",
                "--outfile",
                str(tmp_path / f"figure.{extension}"),
                "--view-out",
                str(view_path),
                "--document-out",
                str(normalized),
            ]
        )
        assert json.loads(normalized.read_text()) == doc.record
        assert json.loads(view_path.read_text())["layout"] == "observed_embedded"
    assert PdfReader(tmp_path / "figure.pdf").pages
    html = (tmp_path / "figure.html").read_text()
    assert "observed-ancestor-label:" in html and "observed-node:G:" in html
    main(
        [
            "draw",
            "--document",
            str(source),
            "--view-config",
            str(view_path),
            "--outfile",
            str(tmp_path / "again.svg"),
        ]
    )
    assert (tmp_path / "figure.svg").read_text().count("observed-node:G:") == (
        tmp_path / "again.svg"
    ).read_text().count("observed-node:G:")
