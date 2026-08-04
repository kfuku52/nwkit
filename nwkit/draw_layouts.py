"""Coordinate generators for :mod:`nwkit.draw`.

The renderers in ``draw.py`` deliberately consume explicit edge paths rather
than inferring a drawing style from node coordinates.  This keeps annotations
attached to the same nodes while allowing orthogonal, warped, and recursive
layouts to coexist.
"""

from __future__ import annotations

import bisect
from dataclasses import dataclass, field
import math


@dataclass
class TreeDrawingLayout:
    """Geometry needed to render one tree layout."""

    name: str
    xcoord: dict
    ycoord: dict
    leaf_order: list
    edge_paths: dict = field(default_factory=dict)
    support_anchors: dict = field(default_factory=dict)
    support_angles: dict = field(default_factory=dict)
    label_angles: dict = field(default_factory=dict)
    root_path: list = field(default_factory=list)
    equal_aspect: bool = False
    spatial: bool = False
    metadata: dict = field(default_factory=dict)

    @property
    def bounds(self):
        points = [
            (float(x), float(y))
            for path in self.edge_paths.values()
            for x, y in path
        ]
        points.extend((float(x), float(y)) for x, y in self.root_path)
        points.extend(
            (float(self.xcoord[node]), float(self.ycoord[node]))
            for node in self.xcoord
        )
        if not points:
            return (0.0, 1.0, 0.0, 1.0)
        xs, ys = zip(*points)
        return (min(xs), max(xs), min(ys), max(ys))


def _normalized_leaf_weights(tree, leaf_weight_by_leaf=None):
    leaves = list(tree.leaves())
    leaf_weight_by_leaf = leaf_weight_by_leaf or {}
    raw = {
        leaf: max(float(leaf_weight_by_leaf.get(leaf, 1.0)), 1e-9)
        for leaf in leaves
    }
    mean = sum(raw.values()) / max(len(raw), 1)
    return leaves, {
        leaf: value / max(mean, 1e-12)
        for leaf, value in raw.items()
    }


def _weighted_leaf_positions(leaves, normalized_weight):
    if not leaves:
        return {}
    positions = {leaves[0]: 0.0}
    for previous, leaf in zip(leaves, leaves[1:]):
        positions[leaf] = positions[previous] + (
            (normalized_weight[previous] + normalized_weight[leaf]) / 2.0
        )
    return positions


def get_rectangular_coordinates(
    tree,
    use_topology_depth=False,
    leaf_weight_by_leaf=None,
):
    """Return conventional phylogram coordinates and the displayed tip order."""

    xcoord = {}
    ycoord = {}
    stack = [(tree, 0.0)]
    while stack:
        node, depth = stack.pop()
        xcoord[node] = float(depth)
        pending = []
        for child in node.get_children():
            child_dist = (
                1.0
                if use_topology_depth
                else (0.0 if child.dist is None else float(child.dist))
            )
            pending.append((child, depth + child_dist))
        stack.extend(reversed(pending))

    leaf_order, normalized_weight = _normalized_leaf_weights(
        tree,
        leaf_weight_by_leaf,
    )
    ycoord.update(_weighted_leaf_positions(leaf_order, normalized_weight))
    for node in tree.traverse(strategy='postorder'):
        if node in ycoord:
            continue
        child_ys = [ycoord[child] for child in node.get_children()]
        ycoord[node] = (
            float(sum(child_ys) / len(child_ys))
            if child_ys
            else float(len(ycoord))
        )
    return xcoord, ycoord, leaf_order


def _root_stub_from_largest_gap(tree, xcoord, ycoord, stub):
    """Point a root stub into the largest empty sector around the root."""

    child_angles = sorted(
        math.atan2(ycoord[child] - ycoord[tree], xcoord[child] - xcoord[tree])
        % (2.0 * math.pi)
        for child in tree.get_children()
        if math.hypot(
            xcoord[child] - xcoord[tree],
            ycoord[child] - ycoord[tree],
        ) > 1e-12
    )
    if child_angles:
        cyclic = child_angles + [child_angles[0] + (2.0 * math.pi)]
        gap_start, gap_end = max(
            zip(cyclic, cyclic[1:]),
            key=lambda pair: pair[1] - pair[0],
        )
        direction = (gap_start + gap_end) / 2.0
    else:
        direction = math.pi
    return [
        (
            xcoord[tree] + (float(stub) * math.cos(direction)),
            ycoord[tree] + (float(stub) * math.sin(direction)),
        ),
        (xcoord[tree], ycoord[tree]),
    ]


def _orthogonal_paths(tree, xcoord, ycoord):
    edge_paths = {}
    support_anchors = {}
    support_angles = {}
    for node in tree.traverse():
        if node.is_root:
            continue
        parent = node.up
        edge_paths[node] = [
            (xcoord[parent], ycoord[parent]),
            (xcoord[parent], ycoord[node]),
            (xcoord[node], ycoord[node]),
        ]
        support_anchors[node] = (
            xcoord[parent] + ((xcoord[node] - xcoord[parent]) * 0.5),
            ycoord[node],
        )
        support_angles[node] = 0.0
    x_values = list(xcoord.values())
    span = max(max(x_values, default=1.0) - min(x_values, default=0.0), 1.0)
    root_stub = max(span * 0.03, 0.03)
    root_path = [
        (xcoord[tree] - root_stub, ycoord[tree]),
        (xcoord[tree], ycoord[tree]),
    ]
    return edge_paths, support_anchors, support_angles, root_path


def rectangular_layout(
    tree,
    use_topology_depth=False,
    leaf_weight_by_leaf=None,
):
    xcoord, ycoord, leaf_order = get_rectangular_coordinates(
        tree,
        use_topology_depth=use_topology_depth,
        leaf_weight_by_leaf=leaf_weight_by_leaf,
    )
    edge_paths, support_anchors, support_angles, root_path = _orthogonal_paths(
        tree,
        xcoord,
        ycoord,
    )
    return TreeDrawingLayout(
        name='rectangular',
        xcoord=xcoord,
        ycoord=ycoord,
        leaf_order=leaf_order,
        edge_paths=edge_paths,
        support_anchors=support_anchors,
        support_angles=support_angles,
        label_angles={leaf: 0.0 for leaf in leaf_order},
        root_path=root_path,
    )


def _straight_paths(tree, xcoord, ycoord):
    edge_paths = {}
    support_anchors = {}
    support_angles = {}
    for node in tree.traverse():
        if node.is_root:
            continue
        parent = node.up
        start = (xcoord[parent], ycoord[parent])
        end = (xcoord[node], ycoord[node])
        edge_paths[node] = [start, end]
        support_anchors[node], support_angles[node] = _path_midpoint([start, end])
    x_values = list(xcoord.values())
    span = max(max(x_values, default=1.0) - min(x_values, default=0.0), 1.0)
    root_stub = max(span * 0.03, 0.03)
    root_path = [
        (xcoord[tree] - root_stub, ycoord[tree]),
        (xcoord[tree], ycoord[tree]),
    ]
    return edge_paths, support_anchors, support_angles, root_path


def slanted_layout(tree, use_topology_depth=False, leaf_weight_by_leaf=None):
    """Return a rooted phylogram with straight parent--child branches."""

    xcoord, ycoord, leaf_order = get_rectangular_coordinates(
        tree,
        use_topology_depth=use_topology_depth,
        leaf_weight_by_leaf=leaf_weight_by_leaf,
    )
    edge_paths, support_anchors, support_angles, root_path = _straight_paths(
        tree,
        xcoord,
        ycoord,
    )
    return TreeDrawingLayout(
        name='slanted',
        xcoord=xcoord,
        ycoord=ycoord,
        leaf_order=leaf_order,
        edge_paths=edge_paths,
        support_anchors=support_anchors,
        support_angles=support_angles,
        label_angles={leaf: 0.0 for leaf in leaf_order},
        root_path=root_path,
    )


def cladogram_layout(tree, leaf_weight_by_leaf=None):
    """Return an aligned-tip, topology-only slanted cladogram."""

    _, ycoord, leaf_order = get_rectangular_coordinates(
        tree,
        use_topology_depth=True,
        leaf_weight_by_leaf=leaf_weight_by_leaf,
    )
    height = {}
    for node in tree.traverse(strategy='postorder'):
        children = node.get_children()
        height[node] = 0.0 if not children else 1.0 + max(height[child] for child in children)
    total_height = height.get(tree, 0.0)
    xcoord = {
        node: total_height - height[node]
        for node in tree.traverse()
    }
    edge_paths, support_anchors, support_angles, root_path = _straight_paths(
        tree,
        xcoord,
        ycoord,
    )
    return TreeDrawingLayout(
        name='cladogram',
        xcoord=xcoord,
        ycoord=ycoord,
        leaf_order=leaf_order,
        edge_paths=edge_paths,
        support_anchors=support_anchors,
        support_angles=support_angles,
        label_angles={leaf: 0.0 for leaf in leaf_order},
        root_path=root_path,
    )


def _polar_coordinates(
    tree,
    use_topology_depth,
    angular_span_degrees=360.0,
    angular_center_degrees=90.0,
    leaf_weight_by_leaf=None,
    base_layout=None,
):
    if base_layout is None:
        base_layout = rectangular_layout(
            tree,
            use_topology_depth=use_topology_depth,
            leaf_weight_by_leaf=leaf_weight_by_leaf,
        )
    radius = base_layout.xcoord
    linear_y = base_layout.ycoord
    leaf_order = base_layout.leaf_order
    angular_span_degrees = float(angular_span_degrees)
    if (
        not math.isfinite(angular_span_degrees)
        or angular_span_degrees <= 0.0
        or angular_span_degrees > 360.0
    ):
        raise ValueError('--angular-span must be greater than zero and no greater than 360.')
    angular_center_degrees = float(angular_center_degrees)
    if not math.isfinite(angular_center_degrees):
        raise ValueError('--angular-center must be a finite number.')
    angular_center_degrees %= 360.0
    count = len(leaf_order)
    available = math.radians(angular_span_degrees)
    center = math.radians(angular_center_degrees)
    if count <= 1:
        start = center
        denominator = 1.0
        linear_origin = 0.0
    elif angular_span_degrees < 360.0:
        start = center - (available / 2.0)
        leaf_positions = [linear_y[leaf] for leaf in leaf_order]
        linear_origin = min(leaf_positions)
        denominator = max(max(leaf_positions) - linear_origin, 1e-12)
    else:
        start = center - math.pi
        _, normalized_weight = _normalized_leaf_weights(
            tree,
            leaf_weight_by_leaf,
        )
        leaf_positions = {leaf: linear_y[leaf] for leaf in leaf_order}
        first_at_seam = min(leaf_order, key=leaf_positions.get)
        last_at_seam = max(leaf_order, key=leaf_positions.get)
        linear_origin = leaf_positions[first_at_seam]
        linear_span = max(leaf_positions.values()) - linear_origin
        # A complete circle also needs clearance across its seam. This formula
        # reproduces the conventional weighted layout.  Tidy packing can move
        # traversal-adjacent leaves past one another in the perpendicular
        # coordinate, so the seam must be defined by the physical extrema,
        # not by the first and last leaves in traversal order.  Keeping the
        # occupied angular interval strictly below one turn makes the polar
        # transform injective away from the common root.
        seam_gap = (
            normalized_weight[first_at_seam]
            + normalized_weight[last_at_seam]
        ) / 2.0
        denominator = max(linear_span + seam_gap, 1e-12)
    angle = {
        node: start + (
            available * (linear_y[node] - linear_origin) / denominator
        )
        for node in tree.traverse()
    }
    xcoord = {
        node: radius[node] * math.cos(angle[node])
        for node in tree.traverse()
    }
    ycoord = {
        node: radius[node] * math.sin(angle[node])
        for node in tree.traverse()
    }
    return radius, angle, xcoord, ycoord, leaf_order


def _sample_arc(radius, start_angle, end_angle):
    if radius <= 1e-12 or abs(end_angle - start_angle) <= 1e-12:
        return [
            (radius * math.cos(start_angle), radius * math.sin(start_angle)),
        ]
    sample_count = max(
        2,
        int(math.ceil(abs(end_angle - start_angle) / math.radians(3.0))) + 1,
    )
    return [
        (
            radius * math.cos(
                start_angle + ((end_angle - start_angle) * index / (sample_count - 1))
            ),
            radius * math.sin(
                start_angle + ((end_angle - start_angle) * index / (sample_count - 1))
            ),
        )
        for index in range(sample_count)
    ]


def polar_layout(
    tree,
    mode,
    use_topology_depth=False,
    angular_span=360.0,
    angular_center=90.0,
    leaf_weight_by_leaf=None,
    base_layout=None,
    subtree_packing='standard',
):
    """Return rooted circular or straight radial geometry in a sector."""

    if mode not in {'circular', 'radial'}:
        raise ValueError('Unsupported polar layout: {}'.format(mode))
    radius, angle, xcoord, ycoord, leaf_order = _polar_coordinates(
        tree,
        use_topology_depth=use_topology_depth,
        angular_span_degrees=angular_span,
        angular_center_degrees=angular_center,
        leaf_weight_by_leaf=leaf_weight_by_leaf,
        base_layout=base_layout,
    )
    edge_paths = {}
    support_anchors = {}
    support_angles = {}
    for node in tree.traverse():
        if node.is_root:
            continue
        parent = node.up
        if mode == 'radial':
            path = [
                (xcoord[parent], ycoord[parent]),
                (xcoord[node], ycoord[node]),
            ]
            support_anchors[node], support_angles[node] = _path_midpoint(path)
        else:
            path = _sample_arc(radius[parent], angle[parent], angle[node])
            radial_end = (xcoord[node], ycoord[node])
            if path[-1] != radial_end:
                path.append(radial_end)
            support_radius = radius[parent] + ((radius[node] - radius[parent]) * 0.5)
            support_anchors[node] = (
                support_radius * math.cos(angle[node]),
                support_radius * math.sin(angle[node]),
            )
            support_angles[node] = math.degrees(angle[node])
        edge_paths[node] = path
    radius_values = list(radius.values())
    radius_span = max(max(radius_values, default=1.0), 1.0)
    root_path = _root_stub_from_largest_gap(
        tree,
        xcoord,
        ycoord,
        radius_span * 0.03,
    )
    return TreeDrawingLayout(
        name=mode,
        xcoord=xcoord,
        ycoord=ycoord,
        leaf_order=leaf_order,
        edge_paths=edge_paths,
        support_anchors=support_anchors,
        support_angles=support_angles,
        label_angles={
            leaf: math.degrees(angle[leaf])
            for leaf in leaf_order
        },
        root_path=root_path,
        equal_aspect=True,
        spatial=True,
        metadata={
            'subtree_packing': subtree_packing,
            'angular_span_degrees': float(angular_span),
            'angular_center_degrees': float(angular_center) % 360.0,
        },
    )


def _add_graph_edge(adjacency, first, second, length):
    adjacency[first].append((second, float(length)))
    adjacency[second].append((first, float(length)))


def _unrooted_graph(tree, use_topology_depth):
    nodes = list(tree.traverse())
    suppressed_root = tree if len(tree.get_children()) == 2 else None
    kept_nodes = [node for node in nodes if node is not suppressed_root]
    adjacency = {node: [] for node in kept_nodes}

    def edge_length(child):
        if use_topology_depth:
            return 1.0
        return 0.0 if child.dist is None else float(child.dist)

    for node in nodes:
        if node.is_root or node is suppressed_root:
            continue
        if suppressed_root is not None and node.up is suppressed_root:
            continue
        _add_graph_edge(adjacency, node.up, node, edge_length(node))
    if suppressed_root is not None:
        first, second = suppressed_root.get_children()
        _add_graph_edge(
            adjacency,
            first,
            second,
            edge_length(first) + edge_length(second),
        )
    return adjacency, suppressed_root


def _root_graph(adjacency, root):
    parent = {root: None}
    order = []
    stack = [root]
    while stack:
        node = stack.pop()
        order.append(node)
        for neighbor, _ in reversed(adjacency[node]):
            if neighbor is parent[node]:
                continue
            parent[neighbor] = node
            stack.append(neighbor)
    return parent, order


def _unrooted_component_weight(adjacency, leaf_weight_by_leaf=None):
    leaf_weight_by_leaf = leaf_weight_by_leaf or {}
    arbitrary_root = next(iter(adjacency))
    parent, order = _root_graph(adjacency, arbitrary_root)
    descendant_weight = {}
    for node in reversed(order):
        own_weight = (
            max(float(leaf_weight_by_leaf.get(node, 1.0)), 1e-9)
            if len(adjacency[node]) <= 1
            else 0.0
        )
        descendant_weight[node] = own_weight + sum(
            descendant_weight[neighbor]
            for neighbor, _ in adjacency[node]
            if parent.get(neighbor) is node
        )
    total_weight = max(descendant_weight[arbitrary_root], 1e-9)

    def component_weight(node, neighbor):
        if parent.get(neighbor) is node:
            return descendant_weight[neighbor]
        return total_weight - descendant_weight[node]

    return order, component_weight


def _unrooted_center(adjacency):
    order, component_weight = _unrooted_component_weight(adjacency)

    candidates = [node for node in order if len(adjacency[node]) > 1] or order
    order_index = {node: index for index, node in enumerate(order)}
    center = min(
        candidates,
        key=lambda node: (
            max(
                (
                    component_weight(node, neighbor)
                    for neighbor, _ in adjacency[node]
                ),
                default=0,
            ),
            order_index[node],
        ),
    )
    return center


def _component_nodes(adjacency, origin, neighbor):
    nodes = []
    stack = [(neighbor, origin)]
    while stack:
        node, previous = stack.pop()
        nodes.append(node)
        for next_node, _ in adjacency[node]:
            if next_node is previous:
                continue
            stack.append((next_node, node))
    return nodes


def _covering_arc(origin, nodes, xcoord, ycoord):
    angles = sorted(
        math.atan2(ycoord[node] - origin[1], xcoord[node] - origin[0])
        % (2.0 * math.pi)
        for node in nodes
        if math.hypot(
            xcoord[node] - origin[0],
            ycoord[node] - origin[1],
        ) > 1e-12
    )
    if not angles:
        return 0.0, 0.0
    if len(angles) == 1:
        return angles[0], 0.0
    cyclic = angles + [angles[0] + (2.0 * math.pi)]
    gap_index = max(
        range(len(angles)),
        key=lambda index: cyclic[index + 1] - cyclic[index],
    )
    start = cyclic[gap_index + 1] % (2.0 * math.pi)
    width = (cyclic[gap_index] + (2.0 * math.pi) - cyclic[gap_index + 1])
    return (start + (width / 2.0)) % (2.0 * math.pi), max(width, 0.0)


def _rotate_component(nodes, pivot, angle, xcoord, ycoord):
    cosine = math.cos(angle)
    sine = math.sin(angle)
    pivot_x = xcoord[pivot]
    pivot_y = ycoord[pivot]
    for node in nodes:
        delta_x = xcoord[node] - pivot_x
        delta_y = ycoord[node] - pivot_y
        xcoord[node] = pivot_x + (cosine * delta_x) - (sine * delta_y)
        ycoord[node] = pivot_y + (sine * delta_x) + (cosine * delta_y)


def _orientation(first, second, third):
    return (
        (second[0] - first[0]) * (third[1] - first[1])
        - (second[1] - first[1]) * (third[0] - first[0])
    )


def _proper_segment_intersection(first, second, third, fourth, tolerance=1e-12):
    if (
        max(first[0], second[0]) < min(third[0], fourth[0]) - tolerance
        or max(third[0], fourth[0]) < min(first[0], second[0]) - tolerance
        or max(first[1], second[1]) < min(third[1], fourth[1]) - tolerance
        or max(third[1], fourth[1]) < min(first[1], second[1]) - tolerance
    ):
        return False
    first_side = _orientation(first, second, third)
    second_side = _orientation(first, second, fourth)
    third_side = _orientation(third, fourth, first)
    fourth_side = _orientation(third, fourth, second)
    if (
        first_side * second_side < -(tolerance ** 2)
        and third_side * fourth_side < -(tolerance ** 2)
    ):
        return True

    def lies_on_segment(point, start, end):
        return (
            min(start[0], end[0]) - tolerance
            <= point[0]
            <= max(start[0], end[0]) + tolerance
            and min(start[1], end[1]) - tolerance
            <= point[1]
            <= max(start[1], end[1]) + tolerance
        )

    return (
        (abs(first_side) <= tolerance and lies_on_segment(third, first, second))
        or (abs(second_side) <= tolerance and lies_on_segment(fourth, first, second))
        or (abs(third_side) <= tolerance and lies_on_segment(first, third, fourth))
        or (abs(fourth_side) <= tolerance and lies_on_segment(second, third, fourth))
    )


def _undirected_graph_edges(adjacency):
    seen = set()
    edges = []
    for first, neighbors in adjacency.items():
        for second, _ in neighbors:
            key = frozenset((id(first), id(second)))
            if key in seen:
                continue
            seen.add(key)
            edges.append((first, second))
    return edges


def _graph_has_crossing(edges, xcoord, ycoord):
    records = []
    for first, second in edges:
        start = (xcoord[first], ycoord[first])
        end = (xcoord[second], ycoord[second])
        records.append({
            'nodes': (first, second),
            'segment': (start, end),
            'x0': min(start[0], end[0]),
            'x1': max(start[0], end[0]),
            'y0': min(start[1], end[1]),
            'y1': max(start[1], end[1]),
        })
    if not records:
        return False
    x_span = max(record['x1'] for record in records) - min(
        record['x0'] for record in records
    )
    y_span = max(record['y1'] for record in records) - min(
        record['y0'] for record in records
    )
    x_density = (
        math.inf
        if x_span <= 1e-12
        else sum(record['x1'] - record['x0'] for record in records) / x_span
    )
    y_density = (
        math.inf
        if y_span <= 1e-12
        else sum(record['y1'] - record['y0'] for record in records) / y_span
    )
    axis = 'x' if x_density <= y_density else 'y'
    ordered = sorted(records, key=lambda record: record[axis + '0'])
    active = []
    for record in ordered:
        active = [
            other
            for other in active
            if other[axis + '1'] >= record[axis + '0']
        ]
        first, second = record['nodes']
        for other in active:
            third, fourth = other['nodes']
            if first in (third, fourth) or second in (third, fourth):
                continue
            if _proper_segment_intersection(
                *record['segment'],
                *other['segment'],
            ):
                return True
        active.append(record)
    return False


def _equal_daylight_adjust(
    adjacency,
    xcoord,
    ycoord,
    max_iterations=5,
    tolerance=0.05,
    rotation_modifier=0.8,
):
    """Iteratively equalize visible angular gaps around internal nodes."""

    internal_nodes = [node for node in adjacency if len(adjacency[node]) > 2]
    if not internal_nodes:
        return {'iterations': 0, 'accepted_rotations': 0, 'rejected_rotations': 0}
    edges = _undirected_graph_edges(adjacency)
    accepted_rotations = 0
    rejected_rotations = 0

    def perform_pass(modifier):
        total_change = 0.0
        adjustments = 0
        for pivot in internal_nodes:
            components = [
                _component_nodes(adjacency, pivot, neighbor)
                for neighbor, _ in adjacency[pivot]
            ]
            origin = (xcoord[pivot], ycoord[pivot])
            arcs = []
            for component in components:
                midpoint, width = _covering_arc(origin, component, xcoord, ycoord)
                arcs.append([midpoint, width, component])
            arcs.sort(key=lambda item: item[0])
            total_daylight = (2.0 * math.pi) - sum(item[1] for item in arcs)
            if total_daylight <= 1e-9:
                continue
            daylight = total_daylight / len(arcs)
            desired_midpoint = arcs[0][0]
            previous_width = arcs[0][1]
            for index in range(1, len(arcs)):
                actual_midpoint, width, component = arcs[index]
                while actual_midpoint < desired_midpoint - math.pi:
                    actual_midpoint += 2.0 * math.pi
                desired_midpoint += (
                    (previous_width / 2.0) + daylight + (width / 2.0)
                )
                while actual_midpoint > desired_midpoint + math.pi:
                    actual_midpoint -= 2.0 * math.pi
                adjustment = (desired_midpoint - actual_midpoint) * modifier
                adjustment = max(
                    min(adjustment, math.pi / 3.0),
                    -math.pi / 3.0,
                )
                if abs(adjustment) > 1e-12:
                    _rotate_component(
                        component,
                        pivot,
                        adjustment,
                        xcoord,
                        ycoord,
                    )
                    total_change += abs(adjustment)
                    adjustments += 1
                previous_width = width
        return total_change, adjustments

    for iteration in range(1, max(int(max_iterations), 1) + 1):
        original = {
            node: (xcoord[node], ycoord[node])
            for node in adjacency
        }
        total_change = 0.0
        adjustments = 0
        accepted = False
        for attempt in range(10):
            for node, (x, y) in original.items():
                xcoord[node] = x
                ycoord[node] = y
            total_change, adjustments = perform_pass(
                float(rotation_modifier) * (0.5 ** attempt)
            )
            if not _graph_has_crossing(edges, xcoord, ycoord):
                accepted = True
                accepted_rotations += adjustments
                break
            rejected_rotations += adjustments
        if not accepted:
            for node, (x, y) in original.items():
                xcoord[node] = x
                ycoord[node] = y
            return {
                'iterations': iteration,
                'accepted_rotations': accepted_rotations,
                'rejected_rotations': rejected_rotations,
            }
        average_change = total_change / max(adjustments, 1)
        if average_change <= float(tolerance):
            return {
                'iterations': iteration,
                'accepted_rotations': accepted_rotations,
                'rejected_rotations': rejected_rotations,
            }
    return {
        'iterations': max(int(max_iterations), 1),
        'accepted_rotations': accepted_rotations,
        'rejected_rotations': rejected_rotations,
    }


def unrooted_layout(
    tree,
    use_topology_depth=False,
    method='equal-angle',
    daylight_iterations=5,
    leaf_weight_by_leaf=None,
):
    """Return an equal-angle or equal-daylight unrooted tree."""

    method = str(method).strip().lower()
    if method not in {'equal-angle', 'equal-daylight'}:
        raise ValueError("'--unrooted-method' must be equal-angle or equal-daylight.")

    adjacency, suppressed_root = _unrooted_graph(tree, use_topology_depth)
    if not adjacency:
        return TreeDrawingLayout(
            name='unrooted',
            xcoord={tree: 0.0},
            ycoord={tree: 0.0},
            leaf_order=list(tree.leaves()),
            equal_aspect=True,
            spatial=True,
        )
    center = _unrooted_center(adjacency)
    _, component_weight = _unrooted_component_weight(
        adjacency,
        leaf_weight_by_leaf=leaf_weight_by_leaf,
    )
    xcoord = {center: 0.0}
    ycoord = {center: 0.0}
    incoming_angle = {center: 0.0}
    stack = [(center, None, -math.pi, math.pi)]
    while stack:
        node, previous, sector_start, sector_end = stack.pop()
        neighbors = [
            (neighbor, length)
            for neighbor, length in adjacency[node]
            if neighbor is not previous
        ]
        total = float(sum(
            component_weight(node, neighbor)
            for neighbor, _ in neighbors
        ))
        if total <= 0.0:
            continue
        cursor = sector_start
        pending = []
        for neighbor, length in neighbors:
            fraction = component_weight(node, neighbor) / total
            child_start = cursor
            child_end = cursor + ((sector_end - sector_start) * fraction)
            direction = (child_start + child_end) / 2.0
            xcoord[neighbor] = xcoord[node] + (length * math.cos(direction))
            ycoord[neighbor] = ycoord[node] + (length * math.sin(direction))
            incoming_angle[neighbor] = direction
            pending.append((neighbor, node, child_start, child_end))
            cursor = child_end
        stack.extend(reversed(pending))

    daylight_metadata = {
        'daylight_iterations_actual': 0,
        'daylight_rotations_accepted': 0,
        'daylight_rotations_rejected': 0,
    }
    if method == 'equal-daylight':
        if len(adjacency) > 2000:
            raise ValueError(
                "Equal-daylight refinement is limited to 2,000 displayed nodes; "
                "use '--unrooted-method equal-angle' or collapse the drawing."
            )
        daylight_result = _equal_daylight_adjust(
            adjacency,
            xcoord,
            ycoord,
            max_iterations=daylight_iterations,
        )
        daylight_metadata = {
            'daylight_iterations_actual': daylight_result['iterations'],
            'daylight_rotations_accepted': daylight_result['accepted_rotations'],
            'daylight_rotations_rejected': daylight_result['rejected_rotations'],
        }

    if suppressed_root is not None:
        first, second = suppressed_root.get_children()
        first_length = 1.0 if use_topology_depth else float(first.dist or 0.0)
        second_length = 1.0 if use_topology_depth else float(second.dist or 0.0)
        total_length = first_length + second_length
        fraction = 0.5 if total_length <= 0.0 else first_length / total_length
        xcoord[suppressed_root] = xcoord[first] + ((xcoord[second] - xcoord[first]) * fraction)
        ycoord[suppressed_root] = ycoord[first] + ((ycoord[second] - ycoord[first]) * fraction)

    edge_paths = {}
    support_anchors = {}
    support_angles = {}
    for node in tree.traverse():
        if node.is_root:
            continue
        path = [
            (xcoord[node.up], ycoord[node.up]),
            (xcoord[node], ycoord[node]),
        ]
        edge_paths[node] = path
        support_anchors[node], support_angles[node] = _path_midpoint(path)
    leaves = list(tree.leaves())
    label_angles = {}
    for leaf in leaves:
        neighbors = adjacency.get(leaf, [])
        if neighbors:
            neighbor = neighbors[0][0]
            label_angles[leaf] = math.degrees(
                math.atan2(
                    ycoord[leaf] - ycoord[neighbor],
                    xcoord[leaf] - xcoord[neighbor],
                )
            )
        else:
            label_angles[leaf] = 0.0
    leaf_order = sorted(
        leaves,
        key=lambda leaf: math.atan2(ycoord[leaf], xcoord[leaf]),
    )
    return TreeDrawingLayout(
        name='unrooted' if method == 'equal-angle' else 'unrooted-daylight',
        xcoord=xcoord,
        ycoord=ycoord,
        leaf_order=leaf_order,
        edge_paths=edge_paths,
        support_anchors=support_anchors,
        support_angles=support_angles,
        label_angles=label_angles,
        equal_aspect=True,
        spatial=True,
        metadata=daylight_metadata,
    )


class _TidyBox:
    """Internal state for van der Ploeg's non-layered tidy-tree algorithm."""

    __slots__ = (
        'node', 'w', 'h', 'fixed', 'x', 'prelim', 'mod', 'shift', 'change',
        'tl', 'tr', 'el', 'er', 'msel', 'mser', 'children',
    )

    def __init__(self, node, width, height, fixed, children):
        self.node = node
        self.w = float(width)
        self.h = float(height)
        self.fixed = float(fixed)
        self.x = 0.0
        self.prelim = 0.0
        self.mod = 0.0
        self.shift = 0.0
        self.change = 0.0
        self.tl = None
        self.tr = None
        self.el = None
        self.er = None
        self.msel = 0.0
        self.mser = 0.0
        self.children = children


class _IndexYList:
    __slots__ = ('low', 'index', 'next')

    def __init__(self, low, index, next_item):
        self.low = float(low)
        self.index = int(index)
        self.next = next_item


def _tidy_bottom(box):
    return box.fixed + box.h


def _tidy_set_extremes(box):
    if not box.children:
        box.el = box
        box.er = box
        box.msel = 0.0
        box.mser = 0.0
        return
    box.el = box.children[0].el
    box.msel = box.children[0].msel
    box.er = box.children[-1].er
    box.mser = box.children[-1].mser


def _tidy_next_left_contour(box):
    return box.tl if not box.children else box.children[0]


def _tidy_next_right_contour(box):
    return box.tr if not box.children else box.children[-1]


def _tidy_distribute_extra(box, index, sibling_index, distance):
    if sibling_index == index - 1:
        return
    count = index - sibling_index
    box.children[sibling_index + 1].shift += distance / count
    box.children[index].shift -= distance / count
    box.children[index].change -= distance - (distance / count)


def _tidy_move_subtree(box, index, sibling_index, distance):
    child = box.children[index]
    child.mod += distance
    child.msel += distance
    child.mser += distance
    _tidy_distribute_extra(box, index, sibling_index, distance)


def _tidy_set_left_thread(box, index, contour, modifier_sum):
    left_extreme = box.children[0].el
    left_extreme.tl = contour
    difference = (modifier_sum - contour.mod) - box.children[0].msel
    left_extreme.mod += difference
    left_extreme.prelim -= difference
    box.children[0].el = box.children[index].el
    box.children[0].msel = box.children[index].msel


def _tidy_set_right_thread(box, index, contour, modifier_sum):
    right_extreme = box.children[index].er
    right_extreme.tr = contour
    difference = (modifier_sum - contour.mod) - box.children[index].mser
    right_extreme.mod += difference
    right_extreme.prelim -= difference
    box.children[index].er = box.children[index - 1].er
    box.children[index].mser = box.children[index - 1].mser


def _tidy_separate(box, index, index_y_list):
    right = box.children[index - 1]
    right_modifiers = right.mod
    left = box.children[index]
    left_modifiers = left.mod
    first_pair = True
    while right is not None and left is not None:
        while index_y_list is not None and _tidy_bottom(right) > index_y_list.low:
            index_y_list = index_y_list.next
        if index_y_list is None:
            break
        distance = (
            right_modifiers + right.prelim + right.w
        ) - (
            left_modifiers + left.prelim
        )
        # The first-pair condition incorporates the correction published with
        # the reference implementation after van der Ploeg (2014).
        if (first_pair and distance < 0.0) or distance > 0.0:
            left_modifiers += distance
            _tidy_move_subtree(box, index, index_y_list.index, distance)
            first_pair = False
        right_bottom = _tidy_bottom(right)
        left_bottom = _tidy_bottom(left)
        if right_bottom <= left_bottom:
            right = _tidy_next_right_contour(right)
            if right is not None:
                right_modifiers += right.mod
        if right_bottom >= left_bottom:
            left = _tidy_next_left_contour(left)
            if left is not None:
                left_modifiers += left.mod
    if right is None and left is not None:
        _tidy_set_left_thread(box, index, left, left_modifiers)
    elif right is not None and left is None:
        _tidy_set_right_thread(box, index, right, right_modifiers)


def _tidy_update_index_y(low, index, index_y_list):
    while index_y_list is not None and low >= index_y_list.low:
        index_y_list = index_y_list.next
    return _IndexYList(low, index, index_y_list)


def _tidy_position_root(box):
    first = box.children[0]
    last = box.children[-1]
    box.prelim = (
        first.prelim + first.mod + last.mod + last.prelim + last.w
    ) / 2.0 - (box.w / 2.0)


def _tidy_first_walk(root_box):
    boxes = []
    stack = [root_box]
    while stack:
        box = stack.pop()
        boxes.append(box)
        stack.extend(box.children)
    for box in reversed(boxes):
        if not box.children:
            _tidy_set_extremes(box)
            continue
        index_y_list = _tidy_update_index_y(
            _tidy_bottom(box.children[0].el),
            0,
            None,
        )
        for index in range(1, len(box.children)):
            child = box.children[index]
            lowest = _tidy_bottom(child.er)
            _tidy_separate(box, index, index_y_list)
            index_y_list = _tidy_update_index_y(lowest, index, index_y_list)
        _tidy_position_root(box)
        _tidy_set_extremes(box)


def _tidy_add_child_spacing(box):
    delta = 0.0
    modifier_delta = 0.0
    for child in box.children:
        delta += child.shift
        modifier_delta += delta + child.change
        child.mod += modifier_delta


def _tidy_second_walk(root_box, modifier_sum, coordinates):
    minimum = math.inf
    stack = [(root_box, float(modifier_sum))]
    while stack:
        box, inherited = stack.pop()
        current = inherited + box.mod
        box.x = box.prelim + current
        if box.node is not None:
            coordinates[box.node] = box.x + (box.w / 2.0)
        minimum = min(minimum, box.x)
        _tidy_add_child_spacing(box)
        stack.extend((child, current) for child in reversed(box.children))
    return minimum if math.isfinite(minimum) else 0.0


def _tidy_third_walk(root_box, shift, coordinates):
    stack = [root_box]
    while stack:
        box = stack.pop()
        box.x += shift
        if box.node is not None:
            coordinates[box.node] = box.x + (box.w / 2.0)
        stack.extend(reversed(box.children))


def tidy_layout(
    tree,
    use_topology_depth=False,
    terminal_extent_by_leaf=None,
    leaf_weight_by_leaf=None,
):
    """Draw a branch-length-aware non-layered tidy phylogram.

    This is an orientation-swapped implementation of van der Ploeg's
    public-domain linear-time algorithm.  A node box spans its incoming branch
    along the fixed axis; the algorithm calculates the compact perpendicular
    coordinate.  Optional terminal extents make branch-end labels part of the
    collision geometry.
    """

    terminal_extent_by_leaf = terminal_extent_by_leaf or {}
    xcoord, _, leaf_order = get_rectangular_coordinates(
        tree,
        use_topology_depth=use_topology_depth,
        leaf_weight_by_leaf=leaf_weight_by_leaf,
    )
    _, normalized_weight = _normalized_leaf_weights(
        tree,
        leaf_weight_by_leaf,
    )

    box_by_node = {}
    label_aware = leaf_weight_by_leaf is not None
    for node in tree.traverse(strategy='postorder'):
        children = [box_by_node[child] for child in node.get_children()]
        if node.is_root:
            fixed = xcoord[node]
            height = 0.0
        else:
            fixed = xcoord[node.up]
            height = max(xcoord[node] - fixed, 0.0)
        terminal_extent = (
            max(float(terminal_extent_by_leaf.get(node, 0.0)), 0.0)
            if node.is_leaf
            else 0.0
        )
        if node.is_leaf and label_aware and terminal_extent > 0.0:
            # A multiline label occupies only the terminal part of the fixed
            # axis, not the leaf's entire incoming branch. Represent it as a
            # one-child terminal box so the contour algorithm can let labels
            # and branches at disjoint horizontal positions pass vertically.
            children = [_TidyBox(
                node=None,
                width=normalized_weight[node],
                height=terminal_extent,
                fixed=xcoord[node],
                children=[],
            )]
            width = 1.0
        else:
            height += terminal_extent
            width = normalized_weight[node] if node.is_leaf else 1.0
        box_by_node[node] = _TidyBox(
            node=node,
            width=width,
            height=height,
            fixed=fixed,
            children=children,
        )
    root_box = box_by_node[tree]
    _tidy_first_walk(root_box)
    ycoord = {}
    minimum = _tidy_second_walk(root_box, 0.0, ycoord)
    if minimum != 0.0:
        _tidy_third_walk(root_box, -minimum, ycoord)
    if ycoord:
        origin = min(ycoord.values())
        ycoord = {node: value - origin for node, value in ycoord.items()}
    # The tidy algorithm allocates non-overlapping subtree boxes. Its box
    # midpoint is not necessarily the midpoint of the child nodes when leaf
    # boxes have different label-aware widths. Keep the compacted leaf
    # positions, then restore the phylogenetic drawing convention that every
    # parent joins its direct children at their arithmetic mean height.
    for node in tree.traverse(strategy='postorder'):
        children = node.get_children()
        if children:
            ycoord[node] = sum(ycoord[child] for child in children) / len(children)
    edge_paths, support_anchors, support_angles, root_path = _orthogonal_paths(
        tree,
        xcoord,
        ycoord,
    )
    return TreeDrawingLayout(
        name='rectangular',
        xcoord=xcoord,
        ycoord=ycoord,
        leaf_order=leaf_order,
        edge_paths=edge_paths,
        support_anchors=support_anchors,
        support_angles=support_angles,
        label_angles={leaf: 0.0 for leaf in leaf_order},
        root_path=root_path,
        metadata={'subtree_packing': 'tidy'},
    )


def _sample_segment(start, end, samples):
    samples = max(int(samples), 2)
    return [
        (
            start[0] + ((end[0] - start[0]) * index / (samples - 1)),
            start[1] + ((end[1] - start[1]) * index / (samples - 1)),
        )
        for index in range(samples)
    ]


def _warp_path(path, transform, vertical_sampler):
    warped = []
    for index, (start, end) in enumerate(zip(path, path[1:])):
        if abs(end[1] - start[1]) > 1e-12:
            source_segment = vertical_sampler(start, end)
        else:
            source_segment = _sample_segment(start, end, 8)
        segment = [transform(x, y) for x, y in source_segment]
        if index:
            segment = segment[1:]
        warped.extend(segment)
    return warped


def _path_midpoint(path):
    if not path:
        return (0.0, 0.0), 0.0
    if len(path) == 1:
        return path[0], 0.0
    lengths = []
    total = 0.0
    for start, end in zip(path, path[1:]):
        length = math.hypot(end[0] - start[0], end[1] - start[1])
        lengths.append(length)
        total += length
    target = total / 2.0
    travelled = 0.0
    for (start, end), length in zip(zip(path, path[1:]), lengths):
        if length <= 0.0:
            continue
        if travelled + length >= target:
            fraction = (target - travelled) / length
            point = (
                start[0] + ((end[0] - start[0]) * fraction),
                start[1] + ((end[1] - start[1]) * fraction),
            )
            angle = math.degrees(math.atan2(end[1] - start[1], end[0] - start[0]))
            return point, angle
        travelled += length
    start, end = path[-2], path[-1]
    angle = math.degrees(math.atan2(end[1] - start[1], end[0] - start[0]))
    return path[-1], angle


def spiral_layout(
    tree,
    use_topology_depth=False,
    turns=None,
    aspect_ratio=1.0,
    leaf_weight_by_leaf=None,
    terminal_extent_by_leaf=None,
    subtree_packing='standard',
):
    """Warp an orthogonal phylogram into an Archimedean spiral track."""

    if subtree_packing == 'tidy':
        base = tidy_layout(
            tree,
            use_topology_depth=use_topology_depth,
            terminal_extent_by_leaf=terminal_extent_by_leaf,
            leaf_weight_by_leaf=leaf_weight_by_leaf,
        )
    else:
        base = rectangular_layout(
            tree,
            use_topology_depth=use_topology_depth,
            leaf_weight_by_leaf=leaf_weight_by_leaf,
        )
    leaf_count = max(len(base.leaf_order), 1)
    if turns is None:
        turns = max(1.5, min(32.0, math.sqrt(leaf_count) / 2.0))
    turns = float(turns)
    if not math.isfinite(turns) or turns <= 0.0:
        raise ValueError('--spiral-turns must be a finite number greater than zero.')
    y_values = [base.ycoord[leaf] for leaf in base.leaf_order]
    y_min = min(y_values, default=0.0)
    y_span = max(max(y_values, default=1.0) - y_min, 1.0)
    x_values = list(base.xcoord.values())
    x_min = min(x_values, default=0.0)
    x_span = max(max(x_values, default=1.0) - x_min, 1.0)
    # Parameterize an injective spiral strip.  The centreline advances by one
    # radial pitch per turn, while tree depth occupies only 70% of that pitch.
    # Therefore two points whose angles differ by a complete turn cannot share
    # a radius.  The positive inner margin also prevents the negative-radius
    # fold that affected the former normal-offset construction.
    outer_radius = 1.0
    radial_pitch = outer_radius / (turns + 1.0)
    inner_radius = radial_pitch
    band_width = radial_pitch * 0.70
    horizontal_scale = max(float(aspect_ratio), 1e-6)

    def transform(x, y):
        fraction = (float(y) - y_min) / y_span
        theta = (-math.pi / 2.0) + (2.0 * math.pi * turns * fraction)
        centre_radius = inner_radius + (turns * radial_pitch * fraction)
        depth = (float(x) - x_min) / x_span
        radius = centre_radius + (band_width * (depth - 1.0))
        return (
            horizontal_scale * radius * math.cos(theta),
            radius * math.sin(theta),
        )

    angular_divisions = min(1024, max(48, int(math.ceil(turns * 32.0))))
    global_fractions = sorted({
        0.0,
        1.0,
        *(index / angular_divisions for index in range(angular_divisions + 1)),
        *(
            (float(base.ycoord[node]) - y_min) / y_span
            for node in tree.traverse()
        ),
    })

    def vertical_sampler(start, end):
        """Sample tracks at shared angles and every radial branch junction."""

        start_fraction = (float(start[1]) - y_min) / y_span
        end_fraction = (float(end[1]) - y_min) / y_span
        lower = min(start_fraction, end_fraction)
        upper = max(start_fraction, end_fraction)
        first_index = bisect.bisect_right(global_fractions, lower + 1e-12)
        last_index = bisect.bisect_left(global_fractions, upper - 1e-12)
        fractions = [
            start_fraction,
            *global_fractions[first_index:last_index],
            end_fraction,
        ]
        if end_fraction < start_fraction:
            fractions = list(reversed(sorted(fractions)))
        else:
            fractions = sorted(fractions)
        return [
            (float(start[0]), y_min + (fraction * y_span))
            for fraction in fractions
        ]

    edge_paths = {
        node: _warp_path(path, transform, vertical_sampler)
        for node, path in base.edge_paths.items()
    }
    xcoord = {}
    ycoord = {}
    for node in tree.traverse():
        xcoord[node], ycoord[node] = transform(base.xcoord[node], base.ycoord[node])
    support_anchors = {}
    support_angles = {}
    for node, path in edge_paths.items():
        # Support belongs on the terminal, branch-length-bearing portion of the
        # warped orthogonal edge rather than halfway around a connector arc.
        base_parent = node.up
        base_mid = (
            base.xcoord[base_parent] + ((base.xcoord[node] - base.xcoord[base_parent]) * 0.5),
            base.ycoord[node],
        )
        support_anchors[node] = transform(*base_mid)
        near_a = transform(base_mid[0] - (x_span * 1e-5), base_mid[1])
        near_b = transform(base_mid[0] + (x_span * 1e-5), base_mid[1])
        support_angles[node] = math.degrees(
            math.atan2(near_b[1] - near_a[1], near_b[0] - near_a[0])
        )
    label_angles = {}
    for leaf in base.leaf_order:
        x = base.xcoord[leaf]
        y = base.ycoord[leaf]
        at_leaf = transform(x, y)
        outside = transform(x + (x_span * 0.02), y)
        label_angles[leaf] = math.degrees(
            math.atan2(outside[1] - at_leaf[1], outside[0] - at_leaf[0])
        )
    root_path = _warp_path(base.root_path, transform, vertical_sampler)
    return TreeDrawingLayout(
        name='spiral',
        xcoord=xcoord,
        ycoord=ycoord,
        leaf_order=base.leaf_order,
        edge_paths=edge_paths,
        support_anchors=support_anchors,
        support_angles=support_angles,
        label_angles=label_angles,
        root_path=root_path,
        equal_aspect=True,
        spatial=True,
        metadata={'subtree_packing': subtree_packing},
    )


def fractal_layout(
    tree,
    aspect_ratio=1.0,
    leaf_weight_by_leaf=None,
    layout_name='fractal',
):
    """Return a rectangle-constrained, radial-fractal layout.

    Root clades divide a full turn, and every descendant clade recursively
    receives a nested angular sector in proportion to its tip count. Branch
    increments shrink with clade size, producing self-similar detail without
    branch crossings. The final affine fit uses the available rectangle.
    Geometry intentionally represents topology and balance, not branch length.
    """

    leaf_order = list(tree.leaves())
    leaf_weight_by_leaf = leaf_weight_by_leaf or {}
    leaf_counts = {}
    for node in tree.traverse(strategy='postorder'):
        if node.is_leaf:
            leaf_counts[node] = max(
                float(leaf_weight_by_leaf.get(node, 1.0)),
                1e-6,
            )
        else:
            leaf_counts[node] = max(
                1,
                sum(leaf_counts[child] for child in node.get_children()),
            )
    raw_x = {tree: 0.0}
    raw_y = {tree: 0.0}
    raw_radius = {tree: 0.0}
    incoming_angle = {tree: math.pi}

    stack = [(tree, 0.0, 2.0 * math.pi, 0.0, 1.0, True)]
    while stack:
        node, sector_start, sector_end, radius, scale, is_root = stack.pop()
        children = node.get_children()
        if not children:
            continue
        total = float(sum(leaf_counts[child] for child in children))
        if is_root:
            available_start = -math.pi / 2.0
            available_end = available_start + (2.0 * math.pi)
        else:
            sector_span = sector_end - sector_start
            inset = sector_span * 0.035 if len(children) > 1 else 0.0
            available_start = sector_start + inset
            available_end = sector_end - inset
        available_span = available_end - available_start
        cursor = available_start
        pending = []
        for child in children:
            weight = leaf_counts[child] / total
            child_sector_start = cursor
            child_sector_end = cursor + (available_span * weight)
            child_direction = (child_sector_start + child_sector_end) / 2.0
            if is_root:
                child_scale = scale * (0.45 + (0.10 * math.sqrt(weight)))
            else:
                child_scale = scale * max(0.80, 0.98 * math.sqrt(weight))
            child_radius = radius + child_scale
            raw_x[child] = child_radius * math.cos(child_direction)
            raw_y[child] = child_radius * math.sin(child_direction)
            raw_radius[child] = child_radius
            incoming_angle[child] = child_direction
            pending.append((
                child,
                child_sector_start,
                child_sector_end,
                child_radius,
                child_scale,
                False,
            ))
            cursor = child_sector_end
        stack.extend(reversed(pending))

    raw_x_values = list(raw_x.values())
    raw_y_values = list(raw_y.values())
    raw_x_min = min(raw_x_values, default=-1.0)
    raw_x_max = max(raw_x_values, default=1.0)
    raw_y_min = min(raw_y_values, default=-1.0)
    raw_y_max = max(raw_y_values, default=1.0)
    raw_x_span = max(raw_x_max - raw_x_min, 1e-9)
    raw_y_span = max(raw_y_max - raw_y_min, 1e-9)
    raw_center_x = (raw_x_min + raw_x_max) / 2.0
    raw_center_y = (raw_y_min + raw_y_max) / 2.0
    horizontal_scale = max(float(aspect_ratio), 1e-6)
    fit_x = (2.0 * horizontal_scale * 0.94) / raw_x_span
    fit_y = (2.0 * 0.94) / raw_y_span
    xcoord = {
        node: (raw_x[node] - raw_center_x) * fit_x
        for node in raw_x
    }
    ycoord = {
        node: (raw_y[node] - raw_center_y) * fit_y
        for node in raw_y
    }
    edge_paths = {}
    support_anchors = {}
    support_angles = {}
    label_angles = {}

    def fitted_polar(radius, angle):
        return (
            ((radius * math.cos(angle)) - raw_center_x) * fit_x,
            ((radius * math.sin(angle)) - raw_center_y) * fit_y,
        )

    for node in tree.traverse():
        if node.is_root:
            continue
        parent = node.up
        parent_radius = raw_radius[parent]
        child_radius = raw_radius[node]
        parent_angle = incoming_angle[parent]
        child_angle = incoming_angle[node]
        if parent.is_root:
            polar_path = [(0.0, child_angle)]
        else:
            arc_samples = max(
                2,
                int(math.ceil(abs(child_angle - parent_angle) / math.radians(4.0))) + 1,
            )
            polar_path = [
                (
                    parent_radius,
                    parent_angle + ((child_angle - parent_angle) * index / (arc_samples - 1)),
                )
                for index in range(arc_samples)
            ]
        polar_path.append((child_radius, child_angle))
        edge_paths[node] = [fitted_polar(radius, angle) for radius, angle in polar_path]
        support_radius = parent_radius + ((child_radius - parent_radius) / 2.0)
        support_anchors[node] = fitted_polar(support_radius, child_angle)
        support_start = fitted_polar(support_radius - 1e-5, child_angle)
        support_end = fitted_polar(support_radius + 1e-5, child_angle)
        angle = math.degrees(
            math.atan2(
                support_end[1] - support_start[1],
                support_end[0] - support_start[0],
            )
        )
        support_angles[node] = angle
        if node.is_leaf:
            label_angles[node] = angle

    root_children = tree.get_children()
    if root_children:
        child_angles = sorted(
            incoming_angle.get(child, 0.0) % (2.0 * math.pi)
            for child in root_children
        )
        cyclic = child_angles + [child_angles[0] + (2.0 * math.pi)]
        gap_start, gap_end = max(
            zip(cyclic, cyclic[1:]),
            key=lambda pair: pair[1] - pair[0],
        )
        root_direction = (gap_start + gap_end) / 2.0
    else:
        root_direction = math.pi
    stub = 0.035 * min(2.0 * horizontal_scale, 2.0)
    root_path = [
        (
            xcoord[tree] + (stub * math.cos(root_direction)),
            ycoord[tree] + (stub * math.sin(root_direction)),
        ),
        (xcoord[tree], ycoord[tree]),
    ]
    return TreeDrawingLayout(
        name=layout_name,
        xcoord=xcoord,
        ycoord=ycoord,
        leaf_order=leaf_order,
        edge_paths=edge_paths,
        support_anchors=support_anchors,
        support_angles=support_angles,
        label_angles=label_angles,
        root_path=root_path,
        equal_aspect=True,
        spatial=True,
    )


def make_tree_layout(
    tree,
    layout='rectangular',
    use_topology_depth=False,
    aspect_ratio=1.0,
    spiral_turns=None,
    angular_span=360.0,
    angular_center=90.0,
    terminal_extent_by_leaf=None,
    label_size_by_leaf=None,
    tip_spacing='uniform',
    subtree_packing='standard',
    unrooted_method='equal-angle',
    daylight_iterations=5,
):
    """Dispatch to a validated layout implementation."""

    layout = str(layout).strip().lower()
    tip_spacing = str(tip_spacing).strip().lower()
    subtree_packing = str(subtree_packing).strip().lower()
    supported_layouts = {
        'rectangular', 'slanted', 'cladogram', 'circular', 'radial',
        'unrooted', 'spiral', 'fractal',
    }
    if layout not in supported_layouts:
        raise ValueError("Unsupported '--layout': {}".format(layout))
    if tip_spacing not in {'uniform', 'label-aware'}:
        raise ValueError("'--tip-spacing' must be uniform or label-aware.")
    if subtree_packing not in {'standard', 'tidy'}:
        raise ValueError("'--subtree-packing' must be standard or tidy.")
    tidy_layouts = {'rectangular', 'circular', 'spiral'}
    if subtree_packing == 'tidy' and layout not in tidy_layouts:
        raise ValueError(
            "'--subtree-packing tidy' is supported only with rectangular, "
            'circular, and spiral layouts.'
        )
    angular_span = float(angular_span)
    angular_center = float(angular_center)
    if not math.isfinite(angular_span) or angular_span <= 0.0 or angular_span > 360.0:
        raise ValueError('--angular-span must be greater than zero and no greater than 360.')
    if not math.isfinite(angular_center):
        raise ValueError('--angular-center must be a finite number.')
    if layout not in {'circular', 'radial'}:
        if angular_span != 360.0:
            raise ValueError(
                "'--angular-span' is supported only with circular and radial layouts."
            )
        if angular_center % 360.0 != 90.0:
            raise ValueError(
                "'--angular-center' is supported only with circular and radial layouts."
            )
    label_size_by_leaf = label_size_by_leaf or {}
    leaf_weight_by_leaf = None
    if tip_spacing == 'label-aware':
        leaf_weight_by_leaf = {
            leaf: max(
                float(label_size_by_leaf.get(leaf, (0.0, 0.0))[1]),
                0.0,
            ) + (4.0 / 72.0)
            for leaf in tree.leaves()
        }
    if layout == 'rectangular':
        if subtree_packing == 'tidy':
            return tidy_layout(
                tree,
                use_topology_depth=use_topology_depth,
                terminal_extent_by_leaf=terminal_extent_by_leaf,
                leaf_weight_by_leaf=leaf_weight_by_leaf,
            )
        drawing = rectangular_layout(
            tree,
            use_topology_depth=use_topology_depth,
            leaf_weight_by_leaf=leaf_weight_by_leaf,
        )
        drawing.metadata['subtree_packing'] = subtree_packing
        return drawing
    if layout == 'slanted':
        return slanted_layout(
            tree,
            use_topology_depth=use_topology_depth,
            leaf_weight_by_leaf=leaf_weight_by_leaf,
        )
    if layout == 'cladogram':
        return cladogram_layout(
            tree,
            leaf_weight_by_leaf=leaf_weight_by_leaf,
        )
    if layout == 'spiral':
        return spiral_layout(
            tree,
            use_topology_depth=use_topology_depth,
            turns=spiral_turns,
            aspect_ratio=aspect_ratio,
            leaf_weight_by_leaf=leaf_weight_by_leaf,
            terminal_extent_by_leaf=terminal_extent_by_leaf,
            subtree_packing=subtree_packing,
        )
    if layout == 'fractal':
        return fractal_layout(
            tree,
            aspect_ratio=aspect_ratio,
            leaf_weight_by_leaf=leaf_weight_by_leaf,
        )
    if layout in {'circular', 'radial'}:
        base_layout = None
        if subtree_packing == 'tidy':
            base_layout = tidy_layout(
                tree,
                use_topology_depth=use_topology_depth,
                terminal_extent_by_leaf=terminal_extent_by_leaf,
                leaf_weight_by_leaf=leaf_weight_by_leaf,
            )
        return polar_layout(
            tree,
            mode=layout,
            use_topology_depth=use_topology_depth,
            angular_span=angular_span,
            angular_center=angular_center,
            leaf_weight_by_leaf=leaf_weight_by_leaf,
            base_layout=base_layout,
            subtree_packing=subtree_packing,
        )
    if layout == 'unrooted':
        return unrooted_layout(
            tree,
            use_topology_depth=use_topology_depth,
            method=unrooted_method,
            daylight_iterations=daylight_iterations,
            leaf_weight_by_leaf=leaf_weight_by_leaf,
        )
    raise ValueError("Unsupported '--layout': {}".format(layout))
