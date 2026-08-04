"""Post-render collision checks and deterministic layout reports."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import json
import math
import os
import sys


@dataclass
class DrawingArtist:
    artist: object
    kind: str
    priority: int = 50
    movable: bool = False
    shift_direction: tuple = (0.0, 1.0)
    owner: object = None
    total_shift: float = 0.0

    def shift(self, points):
        if not self.movable or not hasattr(self.artist, 'get_position'):
            return False
        dx, dy = self.shift_direction
        norm = math.hypot(dx, dy)
        if norm <= 1e-12:
            dx, dy, norm = 0.0, 1.0, 1.0
        x, y = self.artist.get_position()
        self.artist.set_position((
            float(x) + (float(points) * dx / norm),
            float(y) + (float(points) * dy / norm),
        ))
        self.total_shift += abs(float(points))
        return True


def _bbox_overlap(first, second, tolerance=0.75):
    overlap_x = min(first.x1, second.x1) - max(first.x0, second.x0)
    overlap_y = min(first.y1, second.y1) - max(first.y0, second.y0)
    return overlap_x > tolerance and overlap_y > tolerance


def _segment_intersects_bbox(start, end, bbox, padding=0.75):
    x0 = bbox.x0 - padding
    x1 = bbox.x1 + padding
    y0 = bbox.y0 - padding
    y1 = bbox.y1 + padding
    if max(start[0], end[0]) < x0 or min(start[0], end[0]) > x1:
        return False
    if max(start[1], end[1]) < y0 or min(start[1], end[1]) > y1:
        return False
    if x0 <= start[0] <= x1 and y0 <= start[1] <= y1:
        return True
    if x0 <= end[0] <= x1 and y0 <= end[1] <= y1:
        return True

    dx = end[0] - start[0]
    dy = end[1] - start[1]
    candidates = []
    if abs(dx) > 1e-12:
        candidates.extend((
            (x0 - start[0]) / dx,
            (x1 - start[0]) / dx,
        ))
    if abs(dy) > 1e-12:
        candidates.extend((
            (y0 - start[1]) / dy,
            (y1 - start[1]) / dy,
        ))
    for fraction in candidates:
        if 0.0 <= fraction <= 1.0:
            x = start[0] + (dx * fraction)
            y = start[1] + (dy * fraction)
            if x0 <= x <= x1 and y0 <= y <= y1:
                return True
    return False


def _line_segments(line):
    transformed = line.get_path().transformed(line.get_transform())
    vertices = transformed.vertices
    return list(zip(vertices, vertices[1:]))


def _proper_segment_intersection(first, second):
    """Return whether two finite line segments cross away from endpoints."""

    first_start, first_end = first
    second_start, second_end = second

    def orientation(start, end, point):
        return (
            ((end[0] - start[0]) * (point[1] - start[1]))
            - ((end[1] - start[1]) * (point[0] - start[0]))
        )

    values = (
        orientation(first_start, first_end, second_start),
        orientation(first_start, first_end, second_end),
        orientation(second_start, second_end, first_start),
        orientation(second_start, second_end, first_end),
    )
    coordinate_scale = max(
        1.0,
        abs(float(first_end[0]) - float(first_start[0])),
        abs(float(first_end[1]) - float(first_start[1])),
        abs(float(second_end[0]) - float(second_start[0])),
        abs(float(second_end[1]) - float(second_start[1])),
    )
    tolerance = 1e-10 * coordinate_scale * coordinate_scale
    return (
        ((values[0] > tolerance and values[1] < -tolerance)
         or (values[0] < -tolerance and values[1] > tolerance))
        and ((values[2] > tolerance and values[3] < -tolerance)
             or (values[2] < -tolerance and values[3] > tolerance))
    )


def _count_branch_crossings(branch_lines, segment_limit=8000):
    """Count crossing branch pairs with a bounded spatial-index audit."""

    branch_segments = [
        (owner, segment)
        for owner, line in branch_lines
        for segment in _line_segments(line)
    ]
    if len(branch_segments) > int(segment_limit):
        return 0, False
    cell_size = 36.0
    grid = {}
    tested = set()
    crossing_owners = set()
    for index, (owner, (start, end)) in enumerate(branch_segments):
        x0, x1 = sorted((float(start[0]), float(end[0])))
        y0, y1 = sorted((float(start[1]), float(end[1])))
        cells = [
            (x_cell, y_cell)
            for x_cell in range(
                math.floor(x0 / cell_size),
                math.floor(x1 / cell_size) + 1,
            )
            for y_cell in range(
                math.floor(y0 / cell_size),
                math.floor(y1 / cell_size) + 1,
            )
        ]
        candidates = {
            previous
            for cell in cells
            for previous in grid.get(cell, ())
        }
        for previous in candidates:
            pair = (previous, index)
            if pair in tested:
                continue
            tested.add(pair)
            other_owner, other_segment = branch_segments[previous]
            if owner is other_owner:
                continue
            if _proper_segment_intersection(other_segment, (start, end)):
                crossing_owners.add(tuple(sorted((id(owner), id(other_owner)))))
        for cell in cells:
            grid.setdefault(cell, []).append(index)
    return len(crossing_owners), True


def _artist_bounds(figure, artists):
    figure.canvas.draw()
    renderer = figure.canvas.get_renderer()
    bounds = {}
    for item in artists:
        if not item.artist.get_visible():
            continue
        if item.kind == 'tip_track' and hasattr(item.artist, 'offsetbox'):
            bounds[id(item)] = item.artist.offsetbox.get_window_extent(renderer)
        else:
            bounds[id(item)] = item.artist.get_window_extent(renderer=renderer)
    return bounds


def _find_collisions(figure, artists, branch_lines):
    visible = [item for item in artists if item.artist.get_visible()]
    bounds = _artist_bounds(figure, visible)
    collisions = []
    if bounds:
        overall_x = max(bound.x1 for bound in bounds.values()) - min(
            bound.x0 for bound in bounds.values()
        )
        overall_y = max(bound.y1 for bound in bounds.values()) - min(
            bound.y0 for bound in bounds.values()
        )
        x_density = sum(bound.width for bound in bounds.values()) / max(overall_x, 1.0)
        y_density = sum(bound.height for bound in bounds.values()) / max(overall_y, 1.0)
    else:
        x_density = y_density = 0.0
    sweep_x = x_density <= y_density
    lower = (lambda bound: bound.x0) if sweep_x else (lambda bound: bound.y0)
    upper = (lambda bound: bound.x1) if sweep_x else (lambda bound: bound.y1)
    ordered = sorted(visible, key=lambda item: lower(bounds[id(item)]))
    active = []
    for first in ordered:
        first_bounds = bounds[id(first)]
        active = [
            second
            for second in active
            if upper(bounds[id(second)]) >= lower(first_bounds)
        ]
        for second in active:
            if (
                first.owner is not None
                and first.owner is second.owner
                and {first.kind, second.kind}.issubset(
                    {'tip_label', 'tip_track', 'tip_badge'}
                )
            ):
                continue
            if _bbox_overlap(first_bounds, bounds[id(second)]):
                collisions.append(('artist', first, second))
        active.append(first)
    branch_segments = [
        (owner, segment)
        for owner, line in branch_lines
        for segment in _line_segments(line)
    ]
    cell_size = 36.0
    segment_grid = {}
    broad_segments = []
    for index, (_, (start, end)) in enumerate(branch_segments):
        x0, x1 = sorted((float(start[0]), float(end[0])))
        y0, y1 = sorted((float(start[1]), float(end[1])))
        x_cells = range(math.floor(x0 / cell_size), math.floor(x1 / cell_size) + 1)
        y_cells = range(math.floor(y0 / cell_size), math.floor(y1 / cell_size) + 1)
        cell_count = len(x_cells) * len(y_cells)
        if cell_count > 4096:
            broad_segments.append(index)
            continue
        for x_cell in x_cells:
            for y_cell in y_cells:
                segment_grid.setdefault((x_cell, y_cell), []).append(index)
    for item in visible:
        if item.kind == 'legend':
            continue
        item_bounds = bounds[id(item)]
        candidates = set(broad_segments)
        for x_cell in range(
            math.floor(item_bounds.x0 / cell_size),
            math.floor(item_bounds.x1 / cell_size) + 1,
        ):
            for y_cell in range(
                math.floor(item_bounds.y0 / cell_size),
                math.floor(item_bounds.y1 / cell_size) + 1,
            ):
                candidates.update(segment_grid.get((x_cell, y_cell), ()))
        for index in candidates:
            owner, (start, end) = branch_segments[index]
            if (
                owner is item.owner
                or (
                    item.kind == 'node_pie'
                    and item.owner is not None
                    and getattr(owner, 'up', None) is item.owner
                )
            ):
                continue
            if _segment_intersects_bbox(start, end, item_bounds):
                collisions.append(('branch', item, owner))
                break
    return collisions, bounds, True


def _collision_summary(collisions):
    counter = Counter()
    for collision_type, first, second in collisions:
        if collision_type == 'artist':
            pair = sorted((first.kind, second.kind))
            counter['{}:{}'.format(*pair)] += 1
        else:
            counter['{}:branch'.format(first.kind)] += 1
    return dict(sorted(counter.items()))


def _rectangle_union_area(rectangles):
    if not rectangles:
        return 0.0
    y_values = sorted({value for rect in rectangles for value in (rect.y0, rect.y1)})
    if len(y_values) < 2:
        return 0.0
    y_index = {value: index for index, value in enumerate(y_values)}
    segment_count = len(y_values) - 1
    cover = [0] * (segment_count * 4 + 8)
    length = [0.0] * (segment_count * 4 + 8)

    def update(index, left, right, query_left, query_right, delta):
        if query_right <= left or right <= query_left:
            return
        if query_left <= left and right <= query_right:
            cover[index] += delta
        else:
            middle = (left + right) // 2
            update(index * 2, left, middle, query_left, query_right, delta)
            update(index * 2 + 1, middle, right, query_left, query_right, delta)
        if cover[index] > 0:
            length[index] = y_values[right] - y_values[left]
        elif right - left == 1:
            length[index] = 0.0
        else:
            length[index] = length[index * 2] + length[index * 2 + 1]

    events = []
    for rect in rectangles:
        if rect.x1 <= rect.x0 or rect.y1 <= rect.y0:
            continue
        events.append((rect.x0, 1, y_index[rect.y0], y_index[rect.y1]))
        events.append((rect.x1, -1, y_index[rect.y0], y_index[rect.y1]))
    events.sort(key=lambda event: event[0])
    if not events:
        return 0.0
    area = 0.0
    previous_x = events[0][0]
    cursor = 0
    while cursor < len(events):
        x = events[cursor][0]
        area += (x - previous_x) * length[1]
        while cursor < len(events) and events[cursor][0] == x:
            _, delta, lower, upper = events[cursor]
            update(1, 0, segment_count, lower, upper, delta)
            cursor += 1
        previous_x = x
    return area


def evaluate_drawing(
    figure,
    artists,
    branch_lines,
    policy='resolve',
    max_iterations=24,
    emit_warning=True,
):
    """Resolve movable text collisions and return a serializable report."""

    policy = str(policy).strip().lower()
    if policy not in {'resolve', 'warn', 'error', 'ignore'}:
        raise ValueError("'--collision-policy' must be resolve, warn, error, or ignore.")
    initial, _, branch_check_complete = _find_collisions(
        figure, artists, branch_lines
    )
    iterations = 0
    if policy == 'resolve' and initial:
        movable = [
            item
            for item in artists
            if item.movable and hasattr(item.artist, 'get_position')
        ]

        def snapshot():
            return [
                (item, tuple(item.artist.get_position()), item.total_shift)
                for item in movable
            ]

        def restore(state):
            for item, position, total_shift in state:
                item.artist.set_position(position)
                item.total_shift = total_shift

        best_state = snapshot()
        best_score = (len(initial), 0.0)
        for iterations in range(1, int(max_iterations) + 1):
            collisions, bounds, _ = _find_collisions(
                figure, artists, branch_lines
            )
            score = (
                len(collisions),
                sum(item.total_shift for item in movable),
            )
            if score < best_score:
                best_state = snapshot()
                best_score = score
            moved = False
            for collision_type, first, second in collisions:
                if collision_type == 'artist':
                    candidates = sorted(
                        (first, second),
                        key=lambda item: (item.priority, item.total_shift),
                    )
                else:
                    candidates = [first]
                target = next(
                    (
                        item
                        for item in candidates
                        if item.movable and item.total_shift < 36.0
                    ),
                    None,
                )
                if target is not None:
                    direction = 1.0
                    if collision_type == 'artist':
                        other = second if target is first else first
                        target_bounds = bounds[id(target)]
                        other_bounds = bounds[id(other)]
                        delta = (
                            target_bounds.x0 + target_bounds.x1
                            - other_bounds.x0 - other_bounds.x1,
                            target_bounds.y0 + target_bounds.y1
                            - other_bounds.y0 - other_bounds.y1,
                        )
                        projection = (
                            delta[0] * target.shift_direction[0]
                            + delta[1] * target.shift_direction[1]
                        )
                        if projection < 0.0:
                            direction = -1.0
                    moved = target.shift(direction * 1.75) or moved
            if not collisions or not moved:
                break
        candidate, _, _ = _find_collisions(figure, artists, branch_lines)
        candidate_score = (
            len(candidate),
            sum(item.total_shift for item in movable),
        )
        if candidate_score > best_score:
            restore(best_state)
    final, bounds, branch_check_complete = _find_collisions(
        figure, artists, branch_lines
    )
    branch_crossing_count, branch_crossing_check_complete = (
        _count_branch_crossings(branch_lines)
    )
    if final and policy in {'resolve', 'warn'} and emit_warning:
        sys.stderr.write(
            'Drawing contains {} unresolved collision(s).\n'.format(len(final))
        )
    if final and policy == 'error':
        raise ValueError(
            'Drawing contains {} unresolved collision(s).'.format(len(final))
        )
    figure_area = max(float(figure.bbox.width * figure.bbox.height), 1.0)
    occupied_area = _rectangle_union_area(list(bounds.values()))
    return {
        'collision_policy': policy,
        'initial_collision_count': len(initial),
        'initial_collisions_by_kind': _collision_summary(initial),
        'final_collision_count': len(final),
        'final_collisions_by_kind': _collision_summary(final),
        'resolution_iterations': iterations,
        'moved_artist_count': sum(item.total_shift > 0.0 for item in artists),
        'maximum_artist_shift_points': max(
            (item.total_shift for item in artists),
            default=0.0,
        ),
        'annotation_occupied_fraction': min(occupied_area / figure_area, 1.0),
        'branch_collision_check_complete': branch_check_complete,
        'branch_crossing_count': branch_crossing_count,
        'branch_crossing_check_complete': branch_crossing_check_complete,
        'artist_counts': dict(sorted(Counter(item.kind for item in artists).items())),
    }


def write_layout_report(path, report):
    if path in (None, ''):
        return
    payload = json.dumps(report, sort_keys=True, indent=2) + '\n'
    if path == '-':
        sys.stdout.write(payload)
        return
    destination = os.path.realpath(path)
    os.makedirs(os.path.dirname(destination), exist_ok=True)
    with open(destination, 'w', encoding='utf-8') as handle:
        handle.write(payload)
