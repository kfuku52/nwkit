"""Schematic G topology inside S; placements are display data, never events."""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections import defaultdict

from nwkit.reconciliation_colors import GENE_COLORS, X1_COLOR, x1_internal_children
from nwkit.reconciliation_corridors import sibling_corridor_widths
from nwkit.reconciliation_ordering import EPSILON, choose_orientation, layout_score
from nwkit.reconciliation_overlay import history_layout, recorded_overlay
from nwkit.reconciliation_projection import observed_context
from nwkit.reconciliation_scene import (
    EVENT_STYLE,
    Primitive,
    Scene,
    _canonical_order,
    _coordinates,
    _descendant_path,
    _upper_paths,
)
from nwkit.reconciliation_spacing import horizontal_pitch, horizontal_tracks

NOTICE = (
    "Observed G topology in S: all joins and ancestral spacing are schematic. "
    "They are not inferred coalescences, event times, a stem origin, or evidence "
    "of ILS. Compatible source events and origin retain their recorded species "
    "placements. Source forks replaced by X remain in the X details and JSON; "
    "no H comparison is drawn. Transfer arrows create separate recipient lineages."
)


def _context(document):
    context = observed_context(document)
    return (
        context.rec,
        context.transform,
        context.upper,
        context.target,
        context.mapping,
    )


def _placements(upper, target, mapping):
    depth = {upper.root: 0}
    for nid in upper.order[1:]:
        depth[nid] = depth[upper.parents[nid]] + 1

    def lca(a, b):
        while a != b:
            if depth[a] >= depth[b]:
                a = upper.parents[a]
            else:
                b = upper.parents[b]
        return a

    placed = {}
    for nid in reversed(target.order):
        node = target.nodes[nid]
        if not node["children"]:
            placed[nid] = mapping[node["label"]]
        else:
            branch = placed[node["children"][0]]
            for child in node["children"][1:]:
                branch = lca(branch, placed[child])
            placed[nid] = branch
    return placed


def _lanes(
    upper,
    target,
    placed,
    spacing,
    order=None,
    jumps=(),
    *,
    members=None,
    root_branch=None,
):
    order = _canonical_order(target) if order is None else order
    rank = {nid: i for i, nid in enumerate(order)}
    ends: dict[str, int] = {}
    for nid in reversed(order):
        ends[nid] = max([rank[nid]] + [ends[c] for c in target.nodes[nid]["children"]])
    if members is None:
        members = defaultdict(set)
        for nid in order:
            start = (
                placed[nid]
                if nid in jumps
                else placed[target.parents[nid]]
                if nid in target.parents
                else root_branch or upper.root
            )
            for branch in _descendant_path(upper, start, placed[nid]):
                members[branch].add(nid)
    lanes, widths = {}, dict.fromkeys(upper.nodes, spacing + 8)
    for branch, nodes in members.items():
        ordered = sorted(nodes, key=rank.__getitem__)
        frontier = [
            n
            for i, n in enumerate(ordered)
            if i == len(ordered) - 1 or rank[ordered[i + 1]] > ends[n]
        ]
        indices = [rank[n] for n in frontier]
        widths[branch] = spacing * len(frontier) + 8
        for nid in nodes:
            first, last = (
                bisect_left(indices, rank[nid]),
                bisect_right(indices, ends[nid]) - 1,
            )
            lanes[branch, nid] = spacing * (
                (first + last) / 2 - (len(frontier) - 1) / 2
            )
    return lanes, members, widths


def _junction_routes(upper, target, placed, jumps=(), *, root_branch=None):
    routes, groups = {}, defaultdict(list)
    for nid in target.order:
        start = (
            placed[nid]
            if nid in jumps
            else placed[target.parents[nid]]
            if nid in target.parents
            else root_branch or upper.root
        )
        routes[nid] = _descendant_path(upper, start, placed[nid])
        for before, after in zip(routes[nid][:-1], routes[nid][1:], strict=True):
            groups[before, after].append(nid)
    return routes, groups


def _junction_forks(xy, lanes, groups, order, pitch, *, canonical=False):
    """Center each vertical bundle, keeping overlapping strokes on separate tracks."""
    rank = {n: i for i, n in enumerate(order)}
    spans, directions = defaultdict(list), {}
    for (branch, after), nodes in groups.items():
        for nid in nodes:
            incoming = xy[branch][1] + lanes[branch, nid]
            outgoing = xy[after][1] + lanes[after, nid]
            spans[branch].append(
                (min(incoming, outgoing), max(incoming, outgoing), nid)
            )
            down = outgoing > incoming
            directions[branch, nid] = down, -incoming if down else incoming, rank[nid]
    forks: dict[tuple[str, str], float] = {}
    for branch, strokes in spans.items():
        bundles: list[list[str]] = []
        high = float("-inf")
        for low, end, nid in sorted(strokes):
            if low > high + EPSILON:
                bundles.append([])
            bundles[-1].append(nid)
            high = max(high, end)
        for nodes in bundles:
            nodes.sort(
                key=rank.__getitem__
                if canonical
                else lambda n, branch=branch: directions[branch, n]
            )
            tracks = horizontal_tracks(
                xy[branch][0] + pitch * (len(nodes) - 1) / 2, len(nodes), pitch
            )
            if canonical:
                tracks.reverse()
            forks.update(
                ((branch, nid), x) for nid, x in zip(nodes, tracks, strict=True)
            )
    return forks


def _positions(upper, target, placed, xy, unit, lanes, order=None, *, pitch, forks):
    groups = defaultdict(list)
    for nid in _canonical_order(target) if order is None else order:
        groups[placed[nid]].append(nid)
    first_turn: dict[str, float] = {}
    for (branch, _), x in forks.items():
        first_turn[branch] = min(first_turn.get(branch, x), x)
    positions = {}
    for branch, nodes in groups.items():
        right, center = xy[branch]
        left = xy[upper.parents[branch]][0] if branch in upper.parents else right - unit
        joins = [n for n in nodes if target.nodes[n]["children"]]
        join_right = (
            first_turn[branch] - pitch
            if branch in first_turn
            else (left + right) / 2 + pitch * (len(joins) - 1) / 2
        )
        join_x = dict(
            zip(joins, horizontal_tracks(join_right, len(joins), pitch), strict=True)
        )
        for nid in nodes:
            x = join_x[nid] if target.nodes[nid]["children"] else right
            positions[nid] = x, center + lanes[branch, nid]
    return positions


def _edges(
    upper,
    target,
    xy,
    unit,
    positions,
    lanes,
    routes,
    forks,
    order=None,
    jumps=(),
    *,
    root_branch=None,
):
    paths = {}
    order = _canonical_order(target) if order is None else order
    for nid in order:
        if nid in jumps:
            paths[nid] = [positions[target.parents[nid]], positions[nid]]
            continue
        branches = routes[nid]
        if nid == target.root:
            branch = root_branch or upper.root
            left = (
                xy[upper.parents[branch]][0]
                if branch in upper.parents
                else xy[branch][0] - unit
            )
            start = (left, xy[branch][1] + lanes[branch, nid])
        else:
            start = positions[target.parents[nid]]
        points = [start, (start[0], xy[branches[0]][1] + lanes[branches[0], nid])]
        for before, after in zip(branches[:-1], branches[1:], strict=True):
            x = forks[before, nid]
            points.extend(
                [
                    (x, xy[before][1] + lanes[before, nid]),
                    (x, xy[after][1] + lanes[after, nid]),
                ]
            )
        end = positions[nid]
        points.extend([(end[0], xy[branches[-1]][1] + lanes[branches[-1], nid]), end])
        paths[nid] = [p for i, p in enumerate(points) if i == 0 or p != points[i - 1]]
    return paths


def _entry_lanes(target, placed, jumps, lanes):
    """Keep each recipient entry on its own continuation's reserved lane."""
    for n in reversed(target.order):
        if n in jumps:
            child = target.nodes[n]["children"][0]
            if (placed[n], child) in lanes:
                lanes[placed[n], n] = lanes[placed[n], child]


def _observed_geometry(
    upper, target, placed, spacing, width, view, priority, recorded=None
):
    jumps: dict = {}
    if recorded is not None:
        target, placed, jumps = (
            recorded.plan.tree,
            recorded.plan.placed,
            recorded.plan.jumps,
        )
    lanes, members, widths = _lanes(upper, target, placed, spacing, jumps=jumps)
    widths = sibling_corridor_widths(upper, widths)
    corridor_width = max(widths.values())
    xy, unit, height = _coordinates(
        upper,
        width,
        spacing=spacing,
        leaf_order=(view or {}).get("leaf_order", {}).get(upper.id),
        branch_widths=widths,
    )
    joins: defaultdict[str, int] = defaultdict(int)
    for nid in target.order:
        if target.nodes[nid]["children"]:
            joins[placed[nid]] += 1
    routes, groups = _junction_routes(upper, target, placed, jumps)
    turns: defaultdict[str, int] = defaultdict(int)
    for (branch, _), nodes in groups.items():
        turns[branch] += len(nodes)
    pitch = horizontal_pitch(
        unit,
        [joins[b] + turns[b] for b in members],
        [
            (count, widths[b], max(abs(lanes[b, n]) for n in members[b]))
            for b, count in turns.items()
        ],
        spacing,
    )

    def evaluate(order):
        lanes, _, _ = _lanes(
            upper, target, placed, spacing, order, jumps, members=members
        )
        _entry_lanes(target, placed, jumps, lanes)
        forks = _junction_forks(xy, lanes, groups, order, pitch)
        positions = _positions(
            upper, target, placed, xy, unit, lanes, order, pitch=pitch, forks=forks
        )
        paths = _edges(
            upper,
            target,
            xy,
            unit,
            positions,
            lanes,
            routes,
            forks,
            order,
            jumps,
        )
        return layout_score(
            {n: p for n, p in paths.items() if n not in jumps},
            xy[upper.root],
            widths[upper.root],
        ), (
            positions,
            paths,
        )

    geometry, search = choose_orientation(target, evaluate, priority)
    order = _canonical_order(target)
    lanes, _, _ = _lanes(upper, target, placed, spacing, order, jumps, members=members)
    _entry_lanes(target, placed, jumps, lanes)
    forks = _junction_forks(xy, lanes, groups, order, pitch, canonical=True)
    previous_positions = _positions(
        upper, target, placed, xy, unit, lanes, order, pitch=pitch, forks=forks
    )
    previous_paths = _edges(
        upper,
        target,
        xy,
        unit,
        previous_positions,
        lanes,
        routes,
        forks,
        order,
        jumps,
    )
    previous_score = layout_score(
        {n: p for n, p in previous_paths.items() if n not in jumps},
        xy[upper.root],
        widths[upper.root],
    )
    search["previous_score"] = previous_score
    search["horizontal_pitch"] = pitch
    search["corridor_widths"] = widths
    if previous_score < search["score"]:
        geometry = previous_positions, previous_paths
        search.update(score=previous_score, flipped=(), previous_retained=True)
    positions, paths = geometry
    return xy, unit, height, corridor_width, positions, paths, search


def _crossing_gaps(paths, gap=3):
    """Break horizontal strokes at interior crossings; no white background mask."""
    horizontal = sorted(
        (a[1], min(a[0], b[0]), max(a[0], b[0]), nid, i)
        for nid, points in paths.items()
        for i, (a, b) in enumerate(zip(points[:-1], points[1:], strict=True))
        if a[1] == b[1] and a[0] != b[0]
    )
    heights = [h[0] for h in horizontal]
    cuts = defaultdict(list)
    for nid, points in paths.items():
        for a, b in zip(points[:-1], points[1:], strict=True):
            if a[0] != b[0] or a[1] == b[1]:
                continue
            low, high = sorted((a[1], b[1]))
            for _, left, right, other, index in horizontal[
                bisect_right(heights, low) : bisect_left(heights, high)
            ]:
                if other != nid and left + gap < a[0] < right - gap:
                    cuts[other, index].append(a[0])
    for nid, points in paths.items():
        part = [points[0]]
        for i, (a, b) in enumerate(zip(points[:-1], points[1:], strict=True)):
            direction = 1 if b[0] >= a[0] else -1
            cursor = a[0]
            for x in sorted(cuts[nid, i], reverse=direction < 0):
                if direction * (x - gap * direction - cursor) > 0:
                    part.append((x - gap * direction, a[1]))
                    yield nid, part
                cursor = x + gap * direction
                part = [(cursor, a[1])]
            part.append(b)
        yield nid, part


def _source_ids(recorded):
    context = recorded.context
    return {
        node: tuple(
            (context.rec["id"], eid) for eid in recorded.plan.hidden[block["id"]]
        )
        for block in context.transform["blocks"]
        for node in block["target_nodes"]
    }


def _history_paths(scene, recorded, paths, internal):
    context, target, plan = recorded.context, recorded.context.target, recorded.plan
    source_ids = _source_ids(recorded)
    for index, (nid, points) in enumerate(
        _crossing_gaps({n: p for n, p in paths.items() if n not in plan.jumps})
    ):
        owner = plan.owners[nid]
        observed = owner in target.nodes
        prefix = "observed-edge" if observed else "extinct"
        tree = target.id if observed else context.history.id
        source = owner if observed else plan.main.sources.get(nid, nid)
        scene.primitives.append(
            Primitive(
                "line",
                f"{prefix}:{tree}:{source}:{index}",
                points,
                color=X1_COLOR if nid in internal else GENE_COLORS[0],
                width=2.4 if nid in internal else 1.5,
                owner=(target.id, owner) if observed else (context.rec["id"], source),
                tree=tree,
                annotation=NOTICE,
                source_event_ids=source_ids.get(nid, ()),
            )
        )
    for entry, event in plan.jumps.items():
        scene.primitives.append(
            Primitive(
                "transfer",
                f"transfer:{context.rec['id']}:{event['id']}",
                paths[entry],
                color=EVENT_STYLE["transfer"][0],
                owner=(context.rec["id"], event["recipient_lower_node"]),
                tree=context.history.id,
                annotation="This G edge enters a fresh recipient lineage and continues to its descendants; it does not merge into an existing gene branch.",
            )
        )


def observed_topology_scene(document, *, width=720, font_size=8, view=None):
    """One exact G in S, with compatible events and X-only source details."""
    context = observed_context(document)
    transform, upper, target, mapping = (
        context.transform,
        context.upper,
        context.target,
        context.mapping,
    )
    if target.id in (view or {}).get("leaf_order", {}):
        raise ValueError(
            "Observed G tips follow their species placements; set leaf_order on S instead of G."
        )
    placed = _placements(upper, target, mapping)
    spacing = max(12, font_size * 1.5)
    recorded = history_layout(context, spacing)
    consumed = {n for block in transform["blocks"] for n in block["target_nodes"]}
    xy, unit, height, corridor_width, positions, paths, search = _observed_geometry(
        upper, target, placed, spacing, width, view, consumed, recorded
    )
    widths = search["corridor_widths"]
    scene = Scene(
        width,
        height,
        caption="One G in S with compatible source events and X1; schematic geometry, no inferred times.",
    )
    for nid, points in _upper_paths(upper, xy, unit).items():
        scene.primitives.append(
            Primitive(
                "upper",
                f"observed-upper:{upper.id}:{nid}",
                points,
                color="#e1e8ee",
                width=widths[nid],
                owner=(upper.id, nid),
                tree=upper.id,
                outline_color="#a6b6c4",
            )
        )
        if upper.nodes[nid]["label"]:
            x, y = xy[nid]
            scene.primitives.append(
                Primitive(
                    "text",
                    f"observed-species:{nid}",
                    [(x, y - widths[nid] / 2 - 8)],
                    text=upper.nodes[nid]["label"],
                    color="#596b80",
                    owner=(upper.id, nid),
                    tree=upper.id,
                )
            )
    overlay = recorded_overlay(
        recorded, xy, unit, widths, spacing, paths[target.root], positions, paths
    )
    internal = x1_internal_children(target, transform["blocks"], "target_nodes")
    _history_paths(scene, recorded, paths, internal)
    source_ids = _source_ids(recorded)
    for nid, point in positions.items():
        if nid not in target.nodes:
            continue
        scene.primitives.append(
            Primitive(
                "event",
                f"observed-node:{target.id}:{nid}",
                [point],
                text="gene_node",
                color=X1_COLOR if nid in consumed else GENE_COLORS[0],
                marker="O",
                owner=(target.id, nid),
                tree=target.id,
                annotation=NOTICE,
                source_event_ids=source_ids.get(nid, ()),
            )
        )
        if not target.nodes[nid]["children"]:
            scene.primitives.append(
                Primitive(
                    "text",
                    f"observed-label:{nid}",
                    [(point[0] + 6, point[1])],
                    text=target.nodes[nid]["label"],
                    owner=(target.id, nid),
                    tree=target.id,
                )
            )
    scene.primitives.append(
        Primitive(
            "text",
            "observed-ancestor-label:0",
            [(10, xy[upper.root][1] - widths[upper.root] / 2 - 16)],
            text="Ancestor",
            color="#596b80",
            annotation=NOTICE,
        )
    )
    scene.primitives.extend(p for p in overlay.primitives if p.kind != "line")
    scene.event_points.update(overlay.event_points)
    low = (
        min(
            y - (p.width / 2 if p.kind == "upper" else 8)
            for p in scene.primitives
            for _, y in p.points
        )
        - 12
    )
    for primitive in scene.primitives:
        primitive.points = [(x, y - low) for x, y in primitive.points]
    scene.event_points = {
        key: (x, y - low) for key, (x, y) in scene.event_points.items()
    }
    scene.height = (
        max(
            y + (p.width / 2 if p.kind == "upper" else 8)
            for p in scene.primitives
            for _, y in p.points
        )
        + 15
    )
    scene.notices.append(NOTICE)
    scene.notices.append(
        "One set of teal gene routes preserves G after pruning recorded extinct sides and suppressing display subdivisions. X-replaced source events remain in X details and the scientific JSON. No H comparison or dashed gene routes are drawn."
    )
    scene.notices.append(
        "Gene branches use the usual teal. Purple marks only the two X1 G nodes and their connecting internal branch; boundary branches keep their usual color. Open circles are G nodes, not additional recorded events."
    )
    scene.notices.append(
        f"Display crossing reduction: {search['evaluated']} orientations evaluated; {search['score'][1]} intersections, {search['score'][2]} in the ancestral stem, {search['score'][0]} overlapping segments. Candidate search is bounded; this is not a global crossing-minimum certificate."
    )
    return scene
