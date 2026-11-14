"""Readable species corridors and honest display of integrated terminal loss."""

import copy
import json
from itertools import combinations
from pathlib import Path
from xml.etree import ElementTree as ET

import pytest

from nwkit.cli import main
from nwkit.reconciliation_composition import apply_view
from nwkit.reconciliation_corridors import corridor_centers, corridor_lanes
from nwkit.reconciliation_document import read_document, validate_document
from nwkit.reconciliation_html import _svg
from nwkit.reconciliation_import import import_tree_document
from nwkit.reconciliation_scene import EVENT_STYLE, _route, build_scene


def loss_document():
    source = import_tree_document("((A,B),C);")
    upper = source["trees"][0]
    lower = import_tree_document("(a,c);")["trees"][0]
    upper["id"], lower["id"] = "S", "H"
    s = {n["label"]: n["id"] for n in upper["nodes"] if not n["children"]}
    h = {n["label"]: n["id"] for n in lower["nodes"] if not n["children"]}
    ab = next(n["id"] for n in upper["nodes"] if set(n["children"]) == {s["A"], s["B"]})
    events = [
        {
            "id": "root",
            "type": "speciation",
            "lower_node": lower["root"],
            "upper_node": upper["root"],
        },
        {
            "id": "unary",
            "type": "speciation",
            "lower_node": h["a"],
            "upper_node": ab,
            "annotations": {"observed_arity": 1, "marginalized_species": s["B"]},
        },
        {"id": "a", "type": "leaf", "lower_node": h["a"], "upper_node": s["A"]},
        {"id": "c", "type": "leaf", "lower_node": h["c"], "upper_node": s["C"]},
    ]
    source.update(
        trees=[upper, lower],
        transformations=[],
        reconciliations=[
            {
                "id": "rec",
                "lower_tree": "H",
                "upper_tree": "S",
                "history_kind": "representative",
                "selection": {
                    "method": "viterbi_derivation",
                    "measurement": "partially_marginalized_ordered_derivation_v1",
                },
                "model": {
                    "id": "undated_DL",
                    "conditions": {
                        "duplication": 0,
                        "transfer": 0,
                        "loss": 0.1,
                        "detection": {"A": 1, "B": 1, "C": 1},
                        "mapping": {"a": "A", "c": "C"},
                    },
                },
                "events": events,
            }
        ],
    )
    return validate_document(source)


def test_species_corridors_retain_space_and_aligned_tips():
    doc = loss_document()
    scene = build_scene(doc, font_size=12)
    upper = doc.trees["S"]
    paths = {p.owner[1]: p for p in scene.primitives if p.kind == "upper"}
    tips = [paths[nid] for nid, node in upper.nodes.items() if not node["children"]]
    assert len({p.points[-1][0] for p in tips}) == 1
    for a, b in zip(tips[:-1], tips[1:], strict=True):
        assert abs(a.points[-1][1] - b.points[-1][1]) > (a.width + b.width) / 2 + 10
    assert all(p.outline_color for p in tips)
    assert len({p.width for p in tips}) == 1  # Empty B remains a species corridor.


def test_loss_stub_ends_inside_B_and_preserves_every_source_record():
    doc = loss_document()
    before = copy.deepcopy(doc.record)
    scene = build_scene(doc)
    stub = next(
        p
        for p in scene.primitives
        if p.id.startswith("integrated-loss:") and p.kind == "line"
    )
    mark = next(
        p
        for p in scene.primitives
        if p.id.startswith("integrated-loss:") and p.kind == "event"
    )
    branch = next(
        p
        for p in scene.primitives
        if p.kind == "upper" and doc.trees["S"].nodes[p.owner[1]]["label"] == "B"
    )
    assert stub.points[-1] == mark.points[0]
    x, y = mark.points[0]
    assert branch.points[0][0] < x < branch.points[-1][0]
    assert y == branch.points[-1][1]
    assert stub.points[0] == scene.event_points["rec", "unary"]
    assert mark.marker == "x" and mark.text == "loss"
    assert "not an inferred loss time" in mark.annotation
    assert set(scene.event_points) == {
        ("rec", e["id"]) for e in before["reconciliations"][0]["events"]
    }
    assert doc.record == before
    assert not any(
        e["type"] == "loss" for e in doc.record["reconciliations"][0]["events"]
    )


@pytest.mark.parametrize(
    "change",
    [
        "sampling",
        "duplication",
        "transfer",
        "unknown_model",
        "unknown_selection",
        "missing_detection",
        "zero_loss",
        "observed_B",
        "recorded_leaf_in_B",
        "boolean_detection",
        "malformed_reference",
    ],
)
def test_empty_observation_is_not_automatically_labeled_loss(change):
    record = loss_document().record
    rec = record["reconciliations"][0]
    conditions = rec["model"]["conditions"]
    if change == "sampling":
        conditions["detection"]["B"] = 0.8
    elif change in {"duplication", "transfer"}:
        conditions[change] = 0.1
    elif change == "unknown_model":
        rec["model"]["id"] = ["unknown"]
    elif change == "unknown_selection":
        rec["selection"].pop("measurement")
    elif change == "missing_detection":
        conditions.pop("detection")
    elif change == "zero_loss":
        conditions["loss"] = 0
    elif change == "observed_B":
        conditions["mapping"]["a"] = "B"
    elif change == "recorded_leaf_in_B":
        rec["events"][3]["upper_node"] = rec["events"][1]["annotations"][
            "marginalized_species"
        ]
    elif change == "boolean_detection":
        conditions["detection"]["B"] = True
    else:
        rec["events"][1]["annotations"]["marginalized_species"] = []
    scene = build_scene(validate_document(record))
    assert not any(p.id.startswith("integrated-loss:") for p in scene.primitives)


def test_collapse_hides_the_derived_stub_with_its_source_lineage():
    doc = loss_document()
    scene = build_scene(doc)
    view = {"collapsed_clades": [{"tree": "H", "node": doc.trees["H"].root}]}
    apply_view(scene, doc, view)
    assert not any(p.id.startswith("integrated-loss:") for p in scene.primitives)


def test_two_histories_keep_separate_loss_endpoints():
    record = loss_document().record
    second = copy.deepcopy(record["reconciliations"][0])
    second["id"] = "second"
    record["reconciliations"].append(second)
    scene = build_scene(validate_document(record))
    marks = [p for p in scene.primitives if p.kind == "event" and p.annotation]
    assert len(marks) == 2
    assert marks[0].points[0][1] != marks[1].points[0][1]


def test_three_layers_do_not_nest_an_extinct_side_as_an_observed_branch():
    record = loss_document().record
    middle = record["trees"][1]
    bottom = copy.deepcopy(middle)
    bottom["id"] = "bottom"
    record["trees"].append(bottom)
    record["reconciliations"].append(
        {
            "id": "lower",
            "lower_tree": "bottom",
            "upper_tree": "H",
            "history_kind": "unspecified",
            "selection": {},
            "model": {},
            "events": [
                {
                    "id": n["id"],
                    "lower_node": n["id"],
                    "upper_node": n["id"],
                    "type": "speciation" if n["children"] else "leaf",
                }
                for n in bottom["nodes"]
            ],
        }
    )
    doc = validate_document(record)
    scene = build_scene(doc)
    assert len(scene.event_points) == sum(
        len(r["events"]) for r in record["reconciliations"]
    )
    assert not any(p.kind == "upper" and p.annotation for p in scene.primitives)
    mark = next(p for p in scene.primitives if p.kind == "event" and p.annotation)
    stub = next(p for p in scene.primitives if p.kind == "line" and p.annotation)
    assert stub.points[-1] == mark.points[0]


def test_canonical_duplicate_lanes_share_the_parent_lane():
    source = import_tree_document("((a1,a2),b);")
    source["trees"][0]["id"] = "H"
    doc = validate_document(source)
    tree = doc.trees["H"]
    rec = {
        "id": "rec",
        "lower_tree": "H",
        "upper_tree": "H",
        "model": {},
        "selection": {},
        "history_kind": "unspecified",
        "events": [
            {"upper_node": "A", "lower_node": nid}
            for nid in tree.nodes
            if tree.nodes[nid]["label"] != "b" and nid != tree.root
        ],
    }
    lanes, widths = corridor_lanes([rec], doc.trees, 12)
    ids = {n["label"]: nid for nid, n in tree.nodes.items()}
    parent = tree.parents[ids["a1"]]
    assert widths["A"] == 2 * 12 + 8  # Two extant lanes, not three serialized nodes.
    assert (
        lanes["rec", "A", ids["a1"]]
        < lanes["rec", "A", parent]
        < lanes["rec", "A", ids["a2"]]
    )


def test_interleaved_leaf_order_stays_explicit():
    doc = validate_document(import_tree_document("((A,B),C);"))
    tree = next(iter(doc.trees.values()))
    ids = {n["label"]: nid for nid, n in tree.nodes.items()}
    leaves = [ids[s] for s in ["B", "C", "A"]]
    y, _ = corridor_centers(tree, leaves, {nid: 30 for nid in tree.nodes}, 12)
    assert [nid for nid in sorted(leaves, key=y.__getitem__)] == leaves


def test_offset_gene_daughters_fork_at_the_same_parent_x():
    record = loss_document().record
    second = copy.deepcopy(record["reconciliations"][0])
    second["id"] = "second"
    record["reconciliations"].append(second)
    scene = build_scene(validate_document(record))
    assert (
        scene.event_points["rec", "root"][0] != scene.event_points["second", "root"][0]
    )
    for rid in ("rec", "second"):
        parent_x, _ = scene.event_points[rid, "root"]
        daughters = [
            p for p in scene.primitives if p.id.startswith(f"line:{rid}:root:")
        ]
        assert len(daughters) == 2
        for p in daughters:
            fork_x = next(
                a[0]
                for a, b in zip(p.points[:-1], p.points[1:], strict=True)
                if a[0] == b[0] and a[1] != b[1]
            )
            assert fork_x == parent_x


def duplicated_document(gene_tree):
    source = import_tree_document("(A,B);")
    upper = source["trees"][0]
    lower = import_tree_document(gene_tree)["trees"][0]
    upper["id"], lower["id"] = "S", "H"
    species = {n["label"]: n["id"] for n in upper["nodes"]}
    nodes = {n["id"]: n for n in lower["nodes"]}
    events = []
    for node in lower["nodes"]:
        children = node["children"]
        kind = (
            "speciation"
            if children and all(not nodes[c]["children"] for c in children)
            else "duplication"
        )
        events.append(
            {
                "id": node["id"],
                "lower_node": node["id"],
                "type": kind if children else "leaf",
                "upper_node": upper["root"]
                if children
                else species[node["label"][0].upper()],
            }
        )
    source.update(
        trees=[upper, lower],
        reconciliations=[
            {
                "id": "copies",
                "upper_tree": "S",
                "lower_tree": "H",
                "history_kind": "unspecified",
                "selection": {},
                "model": {},
                "events": events,
            }
        ],
    )
    return validate_document(source)


@pytest.mark.parametrize(
    "gene_tree,copies", [("((a1,b1),(a2,b2));", 2), ("(((a1,b1),(a2,b2)),(a3,b3));", 3)]
)
def test_duplicated_lineages_keep_separate_rectangular_speciation_forks(
    gene_tree, copies
):
    doc = duplicated_document(gene_tree)
    rec = doc.record["reconciliations"][0]
    events = rec["events"]
    upper = doc.trees["S"].record
    original = copy.deepcopy(doc.record)
    scene = build_scene(doc)
    speciation = [e for e in events if e["type"] == "speciation"]
    assert len(speciation) == copies
    fork_xs = {scene.event_points["copies", e["id"]][0] for e in speciation}
    assert len(fork_xs) == copies
    junction = next(
        p for p in scene.primitives if p.kind == "upper" and p.owner[1] == upper["root"]
    )
    assert all(abs(x - junction.points[-1][0]) < junction.width / 2 for x in fork_xs)
    for event in speciation:
        point = scene.event_points["copies", event["id"]]
        daughters = [
            p
            for p in scene.primitives
            if p.id.startswith(f"line:copies:{event['id']}:")
        ]
        assert len(daughters) == 2
        for daughter in daughters:
            assert daughter.points[0] == point
            turns = [
                (a, b)
                for a, b in zip(daughter.points[:-1], daughter.points[1:], strict=True)
                if a[1] != b[1]
            ]
            assert turns and all(a[0] == b[0] == point[0] for a, b in turns)
    assert doc.record == original


@pytest.mark.parametrize(
    "gene_tree", ["((a1,b1),(a2,b2));", "(((a1,b1),(a2,b2)),(a3,b3));"]
)
def test_recorded_duplications_and_species_forks_use_the_same_horizontal_pitch(
    gene_tree,
):
    doc = duplicated_document(gene_tree)
    rec = doc.record["reconciliations"][0]
    scene = build_scene(doc, width=380)
    xs = sorted(
        {
            scene.event_points[rec["id"], e["id"]][0]
            for e in rec["events"]
            if e["type"] in {"duplication", "speciation"}
        }
    )
    assert len(xs) >= 3
    assert [b - a for a, b in zip(xs[:-1], xs[1:], strict=True)] == pytest.approx(
        [6] * (len(xs) - 1)
    )


def test_equal_copy_counts_keep_the_same_x_pitch_in_unequal_species_branches():
    source = import_tree_document("((A,B),C);")
    upper = source["trees"][0]
    lower = import_tree_document("(((a1,a2),a3),((c1,c2),c3));")["trees"][0]
    upper["id"], lower["id"] = "S", "H"
    species = {n["label"]: n["id"] for n in upper["nodes"] if not n["children"]}
    ab = next(
        n["id"] for n in upper["nodes"] if n["children"] and n["id"] != upper["root"]
    )
    nodes = {n["id"]: n for n in lower["nodes"]}
    events = []
    for node in lower["nodes"]:
        leaf = node
        while leaf["children"]:
            leaf = nodes[leaf["children"][0]]
        branch = species[leaf["label"][0].upper()]
        if node["id"] == lower["root"]:
            kind, branch = "speciation", upper["root"]
        else:
            kind = "duplication" if node["children"] else "leaf"
            if (
                node["id"] in nodes[lower["root"]]["children"]
                and branch == species["A"]
            ):
                events.append(
                    {
                        "id": "pass-AB",
                        "type": "speciation",
                        "lower_node": node["id"],
                        "upper_node": ab,
                    }
                )
        events.append(
            {
                "id": node["id"],
                "type": kind,
                "lower_node": node["id"],
                "upper_node": branch,
            }
        )
    source.update(
        trees=[upper, lower],
        reconciliations=[
            {
                "id": "copies",
                "upper_tree": "S",
                "lower_tree": "H",
                "history_kind": "unspecified",
                "selection": {},
                "model": {},
                "events": events,
            }
        ],
    )
    doc = validate_document(source)
    scene = build_scene(doc, width=380)
    lengths = {}
    for label in ("A", "C"):
        corridor = next(
            p
            for p in scene.primitives
            if p.kind == "upper" and p.owner[1] == species[label]
        )
        lengths[label] = corridor.points[-1][0] - corridor.points[0][0]
        xs = sorted(
            scene.event_points["copies", e["id"]][0]
            for e in events
            if e["type"] == "duplication" and e["upper_node"] == species[label]
        )
        assert len(xs) == 2 and xs[1] - xs[0] == pytest.approx(6)
    assert lengths["C"] == pytest.approx(2 * lengths["A"])


def with_outer_layer(doc):
    record = copy.deepcopy(doc.record)
    middle = record["trees"][0]
    outer = copy.deepcopy(middle)
    outer["id"] = "outer"
    record["trees"].append(outer)
    record["reconciliations"].append(
        {
            "id": "outer",
            "lower_tree": middle["id"],
            "upper_tree": "outer",
            "history_kind": "unspecified",
            "selection": {},
            "model": {},
            "events": [
                {
                    "id": n["id"],
                    "lower_node": n["id"],
                    "upper_node": n["id"],
                    "type": "speciation" if n["children"] else "leaf",
                }
                for n in middle["nodes"]
            ],
        }
    )
    return validate_document(record)


@pytest.mark.parametrize("copies", [2, 3, 8])
def test_nested_copies_share_lanes_and_keep_distinct_rectangular_forks(copies):
    nodes = [f"(a{i},b{i})" for i in range(copies)]
    while len(nodes) > 1:
        nodes = [
            f"({nodes[i]},{nodes[i + 1]})" if i + 1 < len(nodes) else nodes[i]
            for i in range(0, len(nodes), 2)
        ]
    doc = with_outer_layer(duplicated_document(nodes[0] + ";"))
    original = copy.deepcopy(doc.record)
    scene = build_scene(doc)
    rec = doc.record["reconciliations"][0]
    fork_events = [e for e in rec["events"] if e["type"] == "speciation"]
    assert (
        len({scene.event_points["copies", e["id"]][0] for e in fork_events}) == copies
    )
    root = doc.trees["S"].root
    corridor = next(
        p
        for p in scene.primitives
        if p.kind == "upper" and p.tree == "S" and p.owner == ("outer", root)
    )
    assert corridor.width == copies * 12 + 12
    for event in fork_events:
        start = scene.event_points["copies", event["id"]]
        daughters = [
            p
            for p in scene.primitives
            if p.id.startswith(f"nested-line:copies:{event['id']}:")
        ]
        assert len(daughters) == 2
        for daughter in daughters:
            assert daughter.points[0] == start
            assert all(
                a[0] == b[0] or a[1] == b[1]
                for a, b in zip(daughter.points[:-1], daughter.points[1:], strict=True)
            )
            turns = [
                a[0]
                for a, b in zip(daughter.points[:-1], daughter.points[1:], strict=True)
                if a[1] != b[1]
            ]
            assert turns and turns[0] == start[0]
    assert doc.record == original


@pytest.mark.parametrize(
    "change,kind",
    [
        ("none", "loss"),
        ("sampling", "empty_observation"),
        ("duplication", "empty_observation"),
    ],
)
def test_nested_terminal_side_interpretation_matches_two_layers(change, kind):
    record = copy.deepcopy(loss_document().record)
    if change == "sampling":
        record["reconciliations"][0]["model"]["conditions"]["detection"]["B"] = 0.8
    elif change == "duplication":
        record["reconciliations"][0]["model"]["conditions"]["duplication"] = 0.1
    doc = with_outer_layer(validate_document(record))
    original = copy.deepcopy(doc.record)
    scene = build_scene(doc)
    marks = [p for p in scene.primitives if p.kind == "event" and p.annotation]
    assert len(marks) == 1 and marks[0].text == kind
    mark = marks[0]
    stub = next(p for p in scene.primitives if p.id == mark.id + ":stub")
    assert (
        stub.points[0] == scene.event_points["rec", "unary"]
        and stub.points[-1] == mark.points[0]
    )
    assert "not an inferred" in mark.annotation
    assert len(scene.event_points) == sum(
        len(r["events"]) for r in doc.record["reconciliations"]
    )
    assert doc.record == original


def test_nested_wide_labels_keep_layer_columns_separate(tmp_path, monkeypatch):
    import nwkit.reconciliation_render as render

    record = duplicated_document("((a1,b1),(a2,b2));").record
    for tree in record["trees"]:
        for node in tree["nodes"]:
            if not node["children"]:
                node["label"] += "W" * 70
    doc = with_outer_layer(validate_document(record))
    before = copy.deepcopy(doc.record)
    scene = build_scene(doc, font_size=8)
    original = render._figure
    checked = []

    def inspect(*args, **kwargs):
        fig = original(*args, **kwargs)
        fig.canvas.draw()
        ax = fig.axes[0]
        renderer = fig.canvas.get_renderer()
        labels = [t for t in ax.texts if "WW" in t.get_text()]
        assert len(labels) == 8
        boxes = [t.get_window_extent(renderer) for t in labels]
        assert all(not a.overlaps(b) for a, b in combinations(boxes, 2))
        assert all(b.x0 >= ax.bbox.x0 - 1 and b.x1 <= ax.bbox.x1 + 1 for b in boxes)
        checked.append(True)
        return fig

    monkeypatch.setattr(render, "_figure", inspect)
    render.write_scene(scene, tmp_path / "nested-wide.svg", font_size=8)
    assert checked and doc.record == before


def test_species_width_is_uniform_even_with_parallel_copies():
    doc = read_document(
        Path(__file__).parents[1] / "examples/reconciliation/toy.recphyloxml"
    )
    scene = build_scene(doc)
    assert len({p.width for p in scene.primitives if p.kind == "upper"}) == 1


def test_marginalized_side_is_separate_from_a_transferred_observed_copy():
    record = loss_document().record
    rec = record["reconciliations"][0]
    rec["model"]["conditions"]["transfer"] = 0.08
    branch = rec["events"][1]["annotations"]["marginalized_species"]
    rec["events"][3]["upper_node"] = branch
    scene = build_scene(validate_document(record))
    empty = next(
        p
        for p in scene.primitives
        if p.kind == "event" and p.id.startswith("empty-side:")
    )
    assert empty.marker == "D" and empty.text == "empty_observation"
    assert empty.points[0][1] != scene.event_points["rec", "c"][1]
    assert len(scene.event_points) == 4


def test_renderer_uses_requested_size_for_all_figure_text():
    import matplotlib.pyplot as plt

    from nwkit.reconciliation_render import _append_notes, _figure

    scene = build_scene(loss_document())
    primitives, height = _append_notes(scene, 8)
    fig = _figure(
        scene, primitives, "Title", 8, 720, height + 140, 0, height, 0, 1, height, True
    )
    try:
        assert fig.texts and all(t.get_fontsize() == 8 for t in fig.texts)
        assert all(t.get_fontsize() == 8 for t in fig.legends[0].get_texts())
    finally:
        plt.close(fig)


def test_duplication_daughters_are_rectangular_with_one_fork_x():
    start = (50, 30)
    for end in [(100, 20), (100, 40)]:
        points = _route(None, {}, "A", "A", start, end)
        assert points[0] == start and points[-1] == end
        assert points[1] == (start[0], end[1])
        assert all(
            a[0] == b[0] or a[1] == b[1]
            for a, b in zip(points[:-1], points[1:], strict=True)
        )


def test_transfer_has_a_distinct_hollow_source_and_semantic_palette():
    doc = read_document(
        Path(__file__).parents[1] / "examples/reconciliation/toy.recphyloxml"
    )
    scene = build_scene(doc)
    donor = next(
        p for p in scene.primitives if p.kind == "event" and p.text == "transfer_out"
    )
    arrow = next(p for p in scene.primitives if p.kind == "transfer")
    assert donor.marker == "O" and donor.points[0] == arrow.points[0]
    assert donor.color == arrow.color == "#e67e00"
    assert EVENT_STYLE["duplication"] == ("#d62728", "s")
    assert EVENT_STYLE["speciation"] == ("#0072b2", "o")
    assert EVENT_STYLE["loss"] == ("#d62728", "x")
    svg = ET.fromstring(_svg(scene, 8))
    circles = [
        e for e in svg.iter() if e.tag.endswith("circle") and e.get("fill") == "white"
    ]
    assert circles and any(e.get("stroke") == donor.color for e in circles)


def test_html_and_vector_exports_retain_display_interpretation(tmp_path):
    doc = loss_document()
    scene = build_scene(doc)
    svg = ET.fromstring(_svg(scene, 8))
    ns = {"s": "http://www.w3.org/2000/svg"}
    assert any(p.get("aria-hidden") == "true" for p in svg.findall("s:polyline", ns))
    details = [
        json.loads(e.attrib["data-detail"])
        for e in svg.iter()
        if "data-detail" in e.attrib
    ]
    assert any(
        d["text"] == "loss" and "Display interpretation" in d["annotation"]
        for d in details
    )
    source = tmp_path / "source.json"
    source.write_text(json.dumps(doc.record))
    normalized = tmp_path / "normalized.json"
    main(
        [
            "draw",
            "--document",
            str(source),
            "--outfile",
            str(tmp_path / "figure.svg"),
            "--document-out",
            str(normalized),
        ]
    )
    assert json.loads(normalized.read_text()) == doc.record
    assert "integrated-loss:" in (tmp_path / "figure.svg").read_text()
