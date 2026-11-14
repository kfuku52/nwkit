"""Bounded, deterministic branch rotations scored on actual orthogonal paths."""

from __future__ import annotations

import math
from bisect import bisect_left, bisect_right
from collections import defaultdict
from itertools import product

from nwkit.reconciliation_scene import _canonical_order

EPSILON = 1e-8
MAX_LAYOUTS = 128


def path_contacts(paths):
    """Count crossings/touches and overlaps, excluding genuine shared endpoints."""
    horizontal, vertical = [], []
    endpoints = {
        nid: {(round(x, 8), round(y, 8)) for x, y in (points[0], points[-1])}
        for nid, points in paths.items()
    }
    for nid, points in paths.items():
        for a, b in zip(points[:-1], points[1:], strict=True):
            if a == b:
                continue
            if a[1] == b[1]:
                horizontal.append((a[1], min(a[0], b[0]), max(a[0], b[0]), nid))
            elif a[0] == b[0]:
                vertical.append((a[0], min(a[1], b[1]), max(a[1], b[1]), nid))
            else:
                raise ValueError("Crossing reduction requires orthogonal G paths.")
    horizontal.sort()
    heights = [h[0] for h in horizontal]
    contacts = set()
    for x, low, high, nid in vertical:
        for y, left, right, other in horizontal[
            bisect_left(heights, low - EPSILON) : bisect_right(heights, high + EPSILON)
        ]:
            if nid == other or not left - EPSILON <= x <= right + EPSILON:
                continue
            if (round(x, 8), round(y, 8)) in endpoints[nid] & endpoints[other]:
                continue
            contacts.add((*sorted((nid, other)), round(x, 8), round(y, 8)))
    overlaps = set()
    for axis, segments in enumerate((horizontal, vertical)):
        groups = defaultdict(list)
        for position, low, high, nid in segments:
            groups[round(position, 8)].append((low, high, nid))
        for position, members in groups.items():
            if len(members) == 1:
                continue
            active: list[tuple[float, str]] = []
            for low, high, nid in sorted(members):
                active = [(end, owner) for end, owner in active if end >= low - EPSILON]
                for end, other in active:
                    if other == nid:
                        continue
                    if min(end, high) > low + EPSILON:
                        overlaps.add(
                            (
                                *sorted((nid, other)),
                                axis,
                                position,
                                round(low, 8),
                                round(min(end, high), 8),
                            )
                        )
                    else:
                        point = (low, position) if axis == 0 else (position, low)
                        if (round(point[0], 8), round(point[1], 8)) not in endpoints[
                            nid
                        ] & endpoints[other]:
                            contacts.add(
                                (
                                    *sorted((nid, other)),
                                    round(point[0], 8),
                                    round(point[1], 8),
                                )
                            )
                active.append((high, nid))
    return contacts, overlaps


def layout_score(paths, stem_end, corridor_width):
    contacts, overlaps = path_contacts(paths)
    stem = sum(
        x <= stem_end[0] + EPSILON
        and abs(y - stem_end[1]) <= corridor_width / 2 + EPSILON
        for _, _, x, y in contacts
    )
    length = sum(
        math.dist(a, b)
        for points in paths.values()
        for a, b in zip(points[:-1], points[1:], strict=True)
    )
    return len(overlaps), len(contacts), stem, round(length, 8)


def choose_orientation(tree, evaluate, priority=()):
    """Enumerate small trees; bounded local rotations for larger trees, retaining baseline."""
    canonical = _canonical_order(tree)
    rank = {nid: i for i, nid in enumerate(canonical)}
    children = {
        nid: sorted(tree.nodes[nid]["children"], key=rank.__getitem__)
        for nid in canonical
    }
    axes = [nid for nid in canonical if children[nid]]
    seen: set[frozenset[str]] = set()
    best = None
    best_flips: frozenset[str] = frozenset()
    best_geometry = None

    def consider(flips):
        nonlocal best, best_flips, best_geometry
        if flips in seen or len(seen) >= MAX_LAYOUTS:
            return
        seen.add(flips)
        order, pending = [], [tree.root]
        while pending:
            nid = pending.pop()
            order.append(nid)
            child_order = children[nid][::-1] if nid in flips else children[nid]
            pending.extend(reversed(child_order))
        score, geometry = evaluate(order)
        key = (*score, len(flips), tuple(rank[n] for n in canonical if n in flips))
        if best is None or key < best:
            best, best_flips, best_geometry = key, flips, geometry

    consider(frozenset())
    if len(axes) <= 7:
        for bits in product((False, True), repeat=len(axes)):
            consider(
                frozenset(n for n, flipped in zip(axes, bits, strict=True) if flipped)
            )
    else:
        preferred = set(priority)
        axes.sort(key=lambda n: (n not in preferred, rank[n]))
        for _ in range(2):
            previous = best_flips
            for nid in axes:
                consider(best_flips ^ {nid})
                if len(seen) >= MAX_LAYOUTS:
                    break
            if previous == best_flips or len(seen) >= MAX_LAYOUTS:
                break
    assert best is not None
    return best_geometry, {
        "evaluated": len(seen),
        "score": best[:4],
        "flipped": tuple(n for n in canonical if n in best_flips),
        "exhaustive": len(axes) <= 7,
    }
