"""A transfer subdivides one G edge and starts a fresh recipient lineage."""

import copy

import pytest

from nwkit.reconciliation_document import validate_document
from nwkit.reconciliation_observed import _observed_geometry
from nwkit.reconciliation_overlay import history_layout
from nwkit.reconciliation_projection import observed_context
from nwkit.reconciliation_scene import build_scene
from nwkit.reconciliation_view import default_view
from tests.test_reconciliation_projection import accounted_events, mixed_document


def _geometry(doc, width=720, font=8, reverse=False):
    context = observed_context(doc)
    spacing = max(12, font * 1.5)
    layout = history_layout(context, spacing)
    view = default_view(doc)
    if reverse:
        view["leaf_order"] = {
            "S": list(
                reversed(
                    [
                        n
                        for n in context.upper.order
                        if not context.upper.nodes[n]["children"]
                    ]
                )
            )
        }
    geometry = _observed_geometry(
        context.upper, context.target, {}, spacing, width, view, (), layout
    )
    return context, layout.plan, geometry


def _on_segment(point, a, b):
    x, y = point
    return (a[0] == b[0] == x and min(a[1], b[1]) <= y <= max(a[1], b[1])) or (
        a[1] == b[1] == y and min(a[0], b[0]) <= x <= max(a[0], b[0])
    )


@pytest.mark.parametrize("width,font", [(380, 8), (720, 8), (720, 16)])
@pytest.mark.parametrize("reverse", [False, True])
def test_transfer_arrives_on_a_new_branch_not_an_existing_lineage(width, font, reverse):
    doc = mixed_document()
    original = copy.deepcopy(doc.record)
    context, plan, geometry = _geometry(doc, width, font, reverse)
    xy, _, _, _, positions, paths, _ = geometry
    assert len(plan.jumps) == 1
    for entry, event in plan.jumps.items():
        arrow = paths[entry]
        assert arrow[0] == positions[plan.births[event["id"]]]
        assert arrow[-1] == positions[entry]
        branch = event["recipient_node"]
        left = xy[context.upper.parents[branch]][0]
        assert left < arrow[-1][0] < xy[branch][0]
        continuation = plan.tree.nodes[entry]["children"][0]
        assert paths[continuation][0] == arrow[-1]
        assert positions[continuation][0] > arrow[-1][0]
        assert plan.owners[entry] == continuation
        for n, points in paths.items():
            if (
                n not in plan.tree.nodes
                or n in plan.jumps
                or plan.owners[n] == continuation
            ):
                continue
            assert not any(
                _on_segment(arrow[-1], a, b)
                for a, b in zip(points[:-1], points[1:], strict=True)
            )
    # Suppressing the fresh-entry ports gives exactly G, including the X1 order.
    for child, expected_parent in context.target.parents.items():
        parent = plan.tree.parents[child]
        while parent not in context.target.nodes:
            assert len(plan.tree.nodes[parent]["children"]) == 1
            parent = plan.tree.parents[parent]
        assert parent == expected_parent
    scene = build_scene(doc, width=width, font_size=font)
    assert not any(p.dashed for p in scene.primitives)
    assert doc.record == original


def extinct_transfer_document():
    from tests.test_reconciliation_forks import compatible_duplication_document

    record = copy.deepcopy(compatible_duplication_document().record)
    rec = record["reconciliations"][0]
    duplicate = next(e for e in rec["events"] if e["id"] == "compatible-D")
    lost = next(e for e in rec["events"] if e["id"] == "compatible-L")
    species_c = next(n["id"] for n in record["trees"][0]["nodes"] if n["label"] == "A")
    duplicate.update(
        type="transfer",
        donor_node=duplicate["upper_node"],
        recipient_node=species_c,
        recipient_lower_node=lost["lower_node"],
    )
    lost["upper_node"] = species_c
    return validate_document(record)


def test_entirely_lost_transfer_still_has_a_new_entry_and_terminal_loss_branch():
    doc = extinct_transfer_document()
    scene = build_scene(doc)
    assert sum(p.kind == "transfer" for p in scene.primitives) == 2
    loss = next(
        e
        for e in doc.record["reconciliations"][0]["events"]
        if e["id"] == "compatible-L"
    )
    donor = next(
        e
        for e in doc.record["reconciliations"][0]["events"]
        if e.get("recipient_lower_node") == loss["lower_node"]
    )
    arrow = next(p for p in scene.primitives if p.id == f"transfer:R:{donor['id']}")
    spur = next(
        p
        for p in scene.primitives
        if p.kind == "line"
        and p.owner == ("R", loss["lower_node"])
        and p.points[-1] == scene.event_points["R", loss["id"]]
    )
    assert arrow.points[-1] == spur.points[0]
    assert spur.points[-1] == scene.event_points["R", loss["id"]]
    assert spur.points[0][0] < spur.points[-1][0]


def test_extinct_side_terminal_does_not_touch_an_existing_observed_lineage():
    from tests.test_reconciliation_forks import compatible_duplication_document

    scene = build_scene(compatible_duplication_document())
    loss = next(
        p.points[0] for p in scene.primitives if p.kind == "event" and p.text == "loss"
    )
    for p in scene.primitives:
        if p.id.startswith("observed-edge:"):
            assert not any(
                _on_segment(loss, a, b)
                for a, b in zip(p.points[:-1], p.points[1:], strict=True)
            )


def test_explicit_transfer_in_marker_is_at_the_fresh_entry_not_at_the_leaf():
    record = copy.deepcopy(mixed_document().record)
    rec = record["reconciliations"][0]
    transfer = next(e for e in rec["events"] if e["type"] == "transfer")
    transfer["type"] = "transfer_out"
    leaf = next(
        e for e in rec["events"] if e["lower_node"] == transfer["recipient_lower_node"]
    )
    rec["events"].insert(
        rec["events"].index(leaf),
        dict(
            id="incoming",
            type="transfer_in",
            lower_node=leaf["lower_node"],
            upper_node=leaf["upper_node"],
        ),
    )
    scene = build_scene(validate_document(record))
    arrow = next(p for p in scene.primitives if p.kind == "transfer")
    assert arrow.points[-1] == scene.event_points["R", "incoming"]
    assert arrow.points[-1] != scene.event_points["R", leaf["id"]]


def chained_transfer_document():
    record = copy.deepcopy(mixed_document().record)
    rec, full = record["reconciliations"][0], record["trees"][1]
    transfer = next(e for e in rec["events"] if e["type"] == "transfer")
    foreign = transfer["recipient_lower_node"]
    a = next(n["id"] for n in record["trees"][0]["nodes"] if n["label"] == "A")
    b = transfer["recipient_node"]
    parent = next(n for n in full["nodes"] if n["id"] == transfer["lower_node"])
    parent["children"] = ["chain-T" if n == foreign else n for n in parent["children"]]
    full["nodes"].extend(
        [
            dict(
                id="chain-T",
                label="",
                children=[foreign, "chain-loss"],
                branch_length=None,
            ),
            dict(
                id="chain-loss",
                label="",
                children=[],
                branch_length=None,
                annotations={"observed": False},
            ),
        ]
    )
    transfer.update(recipient_node=a, recipient_lower_node="chain-T")
    rec["events"].extend(
        [
            dict(
                id="chain-event",
                type="transfer",
                lower_node="chain-T",
                upper_node=a,
                donor_node=a,
                recipient_node=b,
                recipient_lower_node=foreign,
            ),
            dict(
                id="chain-loss-event",
                type="loss",
                lower_node="chain-loss",
                upper_node=a,
            ),
        ]
    )
    return validate_document(record)


def test_multiple_transfers_along_one_observed_edge_keep_their_order():
    doc = chained_transfer_document()
    _, plan, geometry = _geometry(doc)
    positions, paths = geometry[4:6]
    assert len(plan.jumps) == 2
    first, second = sorted(plan.jumps, key=lambda n: plan.tree.order.index(n))
    assert plan.jumps[second]["id"] == "chain-event"
    chain = []
    n = plan.tree.parents[second]
    while n != first:
        chain.append(n)
        n = plan.tree.parents[n]
    assert chain
    assert paths[second][0] == positions[plan.births["chain-event"]]
    assert len(accounted_events(build_scene(doc))) == 12


def test_internal_recipient_has_a_fresh_entry_before_its_species_descendants():
    record = copy.deepcopy(mixed_document().record)
    species, full = record["trees"][:2]
    rec = record["reconciliations"][0]
    branch = next(n for n in species["nodes"] if n["label"] == "B")
    branch["children"] = ["B1", "B2"]
    species["nodes"].extend(
        [dict(id=n, label=n, children=[], branch_length=None) for n in ("B1", "B2")]
    )
    transfer = next(e for e in rec["events"] if e["type"] == "transfer")
    foreign = transfer["recipient_lower_node"]
    donor = next(n for n in full["nodes"] if n["id"] == transfer["lower_node"])
    donor["children"] = [
        "entry-species-split" if n == foreign else n for n in donor["children"]
    ]
    full["nodes"].append(
        dict(id="entry-species-split", label="", children=[foreign], branch_length=None)
    )
    transfer["recipient_lower_node"] = "entry-species-split"
    for e in rec["events"]:
        if e["type"] == "leaf" and e["upper_node"] == branch["id"]:
            e["upper_node"] = "B1" if e["lower_node"] == foreign else "B2"
    rec["events"].append(
        dict(
            id="entry-speciation",
            type="speciation",
            lower_node="entry-species-split",
            upper_node=branch["id"],
        )
    )
    doc = validate_document(record)
    context, plan, geometry = _geometry(doc)
    entry = next(iter(plan.jumps))
    assert context.upper.nodes[plan.placed[entry]]["children"] == ["B1", "B2"]
    child = plan.tree.nodes[entry]["children"][0]
    # Keep the recorded unary species event on the recipient branch before B1.
    assert plan.placed[child] == branch["id"]
    descendant = plan.tree.nodes[child]["children"][0]
    assert plan.placed[descendant] == "B1"
    assert geometry[5][child][0] == geometry[4][entry]
    assert len(accounted_events(build_scene(doc))) == 11


def test_x1_crossing_a_transfer_recipient_boundary_fails_before_output_replacement(
    tmp_path,
):
    import json

    from nwkit.cli import main
    from tests.test_reconciliation_observed import fixture_document

    record = copy.deepcopy(fixture_document().record)
    species, history = record["trees"][:2]
    inner = next(
        n["id"]
        for n in history["nodes"]
        if len(n["children"]) == 2 and n["id"] != history["root"]
    )
    branch = next(
        e["upper_node"]
        for e in record["reconciliations"][0]["events"]
        if e["lower_node"] == inner
    )
    donor = next(n["id"] for n in species["nodes"] if n["label"] == "C")
    event = next(
        e
        for e in record["reconciliations"][0]["events"]
        if e["lower_node"] == history["root"] and e["type"] == "speciation"
    )
    event.update(
        type="transfer",
        upper_node=donor,
        donor_node=donor,
        recipient_node=branch,
        recipient_lower_node=inner,
    )
    doc = validate_document(record)
    with pytest.raises(ValueError, match="crosses a transfer recipient boundary"):
        build_scene(doc)
    source, output = tmp_path / "source.json", tmp_path / "figure.svg"
    source.write_text(json.dumps(record))
    output.write_text("preserve existing output")
    with pytest.raises(ValueError, match="crosses a transfer recipient boundary"):
        main(["draw", "--document", str(source), "--outfile", str(output)])
    assert output.read_text() == "preserve existing output"
