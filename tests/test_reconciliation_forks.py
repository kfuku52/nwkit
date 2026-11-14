"""Source S/D are real forks, including chains suppressed by loss pruning at X."""

import copy

import pytest

from nwkit.reconciliation_document import validate_document
from nwkit.reconciliation_observed import _observed_geometry
from nwkit.reconciliation_overlay import history_layout
from nwkit.reconciliation_projection import observed_context
from nwkit.reconciliation_scene import build_scene
from tests.test_reconciliation_observed import fixture_document
from tests.test_reconciliation_projection import mixed_document


def _geometry(doc, width=720, font=8, reverse=False):
    context = observed_context(doc)
    spacing = max(12, font * 1.5)
    layout = history_layout(context, spacing)
    tips = [n for n in context.upper.order if not context.upper.nodes[n]["children"]]
    view = {"leaf_order": {context.upper.id: tips[::-1]}} if reverse else None
    geometry = _observed_geometry(
        context.upper, context.target, {}, spacing, width, view, (), layout
    )
    return context, layout, geometry


def _assert_fork(part, node, positions, paths):
    start = positions[node]
    directions = []
    for child in part.tree.nodes[node]["children"]:
        assert paths[child][0] == start
        directions.append(
            next((p[1] - start[1] for p in paths[child] if p[1] != start[1]), 0)
        )
    assert min(directions) < 0 < max(directions)


def assert_source_accounted(scene, doc):
    visible = set(scene.event_points)
    hidden = {key for p in scene.primitives for key in p.source_event_ids}
    expected = {
        (r["id"], e["id"]) for r in doc.record["reconciliations"] for e in r["events"]
    }
    assert not visible & hidden
    assert visible | hidden == expected


@pytest.mark.parametrize("width,font", [(380, 8), (720, 8), (720, 16)])
@pytest.mark.parametrize("reverse", [False, True])
def test_x_replaced_source_forks_stay_in_details_without_a_second_gene_tree(
    width, font, reverse
):
    doc = mixed_document()
    original = copy.deepcopy(doc.record)
    context, layout, geometry = _geometry(doc, width, font, reverse)
    hidden = {eid for ids in layout.plan.hidden.values() for eid in ids}
    for event in context.rec["events"]:
        if event["type"] in {"duplication", "loss", "speciation"}:
            assert event["id"] in hidden
    ordered = next(iter(layout.plan.hidden.values()))
    duplicate = next(e for e in context.rec["events"] if e["type"] == "duplication")
    child_species = next(
        e
        for e in context.rec["events"]
        if e["type"] == "speciation"
        and context.history.parents.get(e["lower_node"]) == duplicate["lower_node"]
    )
    assert ordered.index(duplicate["id"]) < ordered.index(child_species["id"])
    assert set(geometry[4]) == set(layout.plan.tree.nodes)
    scene = build_scene(doc, width=width, font_size=font)
    assert_source_accounted(scene, doc)
    assert not any(
        p.id.startswith(("history-context:", "history-transfer:", "context-"))
        for p in scene.primitives
    )
    assert not any(p.dashed for p in scene.primitives)
    assert doc.record == original


@pytest.mark.parametrize("resolution,multiple", [(1, False), (2, False), (1, True)])
def test_one_g_has_no_h_comparison_for_either_resolution_or_multiple_blocks(
    resolution, multiple
):
    doc = fixture_document(resolution, multiple)
    context, layout, geometry = _geometry(doc)
    scene = build_scene(doc)
    assert_source_accounted(scene, doc)
    assert len(layout.plan.hidden) == len(context.transform["blocks"])
    assert not any(
        p.id.startswith("context-") or p.color == "#4976a6" for p in scene.primitives
    )
    tips = [p for p in scene.primitives if p.id.startswith("observed-label:")]
    assert sorted(p.text for p in tips) == sorted(
        n["label"] for n in context.target.nodes.values() if not n["children"]
    )
    for block in context.transform["blocks"]:
        nodes = [
            p
            for p in scene.primitives
            if p.id.startswith("observed-node:") and p.owner[1] in block["target_nodes"]
        ]
        assert len(nodes) == 2 and all(p.source_event_ids for p in nodes)


def compatible_duplication_document():
    record = copy.deepcopy(mixed_document().record)
    rec, full = record["reconciliations"][0], record["trees"][1]
    transfer = next(e for e in rec["events"] if e["type"] == "transfer")
    for node in full["nodes"]:
        node["children"] = [
            "compatible-D" if n == transfer["lower_node"] else n
            for n in node["children"]
        ]
    full["nodes"].extend(
        [
            dict(
                id="compatible-D",
                label="",
                children=[transfer["lower_node"], "compatible-L"],
                branch_length=None,
            ),
            dict(
                id="compatible-L",
                label="",
                children=[],
                branch_length=None,
                annotations={"observed": False},
            ),
        ]
    )
    for eid, kind in [("compatible-D", "duplication"), ("compatible-L", "loss")]:
        rec["events"].append(
            dict(id=eid, type=kind, lower_node=eid, upper_node=transfer["upper_node"])
        )
    return validate_document(record)


def _projected_shape(part, target):
    result = {}
    for n in reversed(part.tree.order):
        children = [
            result[c] for c in part.tree.nodes[n]["children"] if result[c] is not None
        ]
        result[n] = (
            children[0]
            if len(children) == 1
            else tuple(sorted(children, key=repr))
            if children
            else target.nodes[n]["label"]
            if n in target.nodes and not target.nodes[n]["children"]
            else None
        )
    return result[part.tree.root]


def test_compatible_duplication_is_an_actual_fork_on_g_without_duplicating_the_history():
    doc = compatible_duplication_document()
    context, layout, geometry = _geometry(doc)
    part = layout.plan.main
    _assert_fork(part, part.events["compatible-D"], geometry[4], geometry[5])
    target = context.target
    from nwkit.reconciliation_history import HistoryPart

    expected = _projected_shape(HistoryPart(target, {}, {}, {}), target)
    assert _projected_shape(part, target) == expected
    assert "compatible-D" not in {
        eid for ids in layout.plan.hidden.values() for eid in ids
    }
    assert_source_accounted(build_scene(doc), doc)


def test_native_style_hidden_speciation_on_one_h_node_forks_in_its_recorded_branch():
    record = copy.deepcopy(fixture_document().record)
    species, history = record["trees"][:2]
    rec = record["reconciliations"][0]
    branch = next(n for n in species["nodes"] if n["label"] == "A")
    branch["children"] = ["A1", "A0"]
    species["nodes"].extend(
        dict(id=n, label=n, children=[], branch_length=None) for n in ("A1", "A0")
    )
    gene = next(n["id"] for n in history["nodes"] if n["label"] == "a")
    leaf = next(e for e in rec["events"] if e["lower_node"] == gene)
    leaf["upper_node"] = "A1"
    rec["events"].insert(
        rec["events"].index(leaf),
        dict(
            id="hidden-S",
            type="speciation",
            lower_node=gene,
            upper_node=branch["id"],
            annotations=dict(observed_arity=1, marginalized_species="A0"),
        ),
    )
    doc = validate_document(record)
    context, layout, geometry = _geometry(doc)
    part = layout.plan.main
    node = part.events["hidden-S"]
    assert node != part.events[leaf["id"]]
    assert part.placed[node] == branch["id"]
    _assert_fork(part, node, geometry[4], geometry[5])
    scene = build_scene(doc)
    assert any(p.text == "empty_observation" for p in scene.primitives)
    assert not any(p.text == "loss" for p in scene.primitives)
    assert_source_accounted(scene, doc)


def test_native_style_hidden_transfer_has_a_new_recipient_and_integrated_donor_side():
    record = copy.deepcopy(fixture_document().record)
    species, history = record["trees"][:2]
    rec = record["reconciliations"][0]
    locations = {n["label"]: n["id"] for n in species["nodes"] if n["label"]}
    gene = next(n["id"] for n in history["nodes"] if n["label"] == "c")
    leaf = next(e for e in rec["events"] if e["lower_node"] == gene)
    leaf["upper_node"] = locations["A"]
    rec["events"].insert(
        rec["events"].index(leaf),
        dict(
            id="hidden-T",
            type="transfer",
            lower_node=gene,
            upper_node=locations["C"],
            donor_node=locations["C"],
            recipient_node=locations["A"],
            recipient_lower_node=gene,
            annotations=dict(observed_arity=1, marginalized_species=locations["C"]),
        ),
    )
    doc = validate_document(record)
    context, layout, geometry = _geometry(doc)
    part = layout.plan.main
    entry = next(n for n, e in part.jumps.items() if e["id"] == "hidden-T")
    assert geometry[5][entry][0] == geometry[4][part.events["hidden-T"]]
    assert part.placed[entry] == locations["A"]
    assert part.placed[part.events["hidden-T"]] == locations["C"]
    assert any(
        event["id"] == "hidden-T" and kind == "empty_observation"
        for event, kind, *_ in part.sides.values()
    )
    scene = build_scene(doc)
    assert not any(p.text == "loss" for p in scene.primitives)
    assert_source_accounted(scene, doc)


def test_x_context_keeps_incoming_hidden_transfer_and_daughters_of_its_parent_duplication():
    record = copy.deepcopy(fixture_document().record)
    species, history, target = record["trees"]
    rec = record["reconciliations"][0]
    old_s, old_h, old_g = species["root"], history["root"], target["root"]
    for tree, root, leaf, label in [
        (species, "s-top", "s-D", "D"),
        (history, "h-top", "h-d", "d"),
        (target, "g-top", "g-d", "d"),
    ]:
        tree["nodes"].extend(
            [
                dict(
                    id=root, label="", children=[tree["root"], leaf], branch_length=None
                ),
                dict(id=leaf, label=label, children=[], branch_length=None),
            ]
        )
        tree["root"] = root
    origin = next(e for e in rec["events"] if e["type"] == "origin")
    origin.update(lower_node="h-top", upper_node="s-D")
    first = next(e for e in rec["events"] if e["lower_node"] == old_h)
    rec["events"].insert(
        rec["events"].index(first),
        dict(
            id="incoming-T",
            type="transfer",
            lower_node=old_h,
            upper_node="s-D",
            donor_node="s-D",
            recipient_node=old_s,
            recipient_lower_node=old_h,
            annotations=dict(observed_arity=1, marginalized_species="s-D"),
        ),
    )
    rec["events"].extend(
        [
            dict(
                id="parent-D", type="duplication", lower_node="h-top", upper_node="s-D"
            ),
            dict(id="d-leaf", type="leaf", lower_node="h-d", upper_node="s-D"),
        ]
    )
    doc = validate_document(record)
    context, layout, geometry = _geometry(doc)
    main = layout.plan.main
    _assert_fork(main, main.events["parent-D"], geometry[4], geometry[5])
    assert "incoming-T" in {eid for ids in layout.plan.hidden.values() for eid in ids}
    assert "incoming-T" not in main.events
    assert len(main.jumps) == 1
    assert layout.plan.placed[old_g] == old_s
    scene = build_scene(doc)
    pipe = next(
        p for p in scene.primitives if p.kind == "upper" and p.owner[1] == "s-D"
    )
    point = scene.event_points["R", "origin"]
    assert pipe.points[0][0] < point[0] < pipe.points[-1][0]
    assert abs(point[1] - pipe.points[-1][1]) < pipe.width / 2
    assert_source_accounted(scene, doc)
