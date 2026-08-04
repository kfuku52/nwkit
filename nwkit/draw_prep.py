"""Tree preprocessing used only for rendering."""

from __future__ import annotations

import math

from ete4 import Tree


def _aggregate_leaf_properties(leaves, numeric_mode='none'):
    numeric_mode = str(numeric_mode).strip().lower()
    if numeric_mode not in {'none', 'mean'}:
        raise ValueError("'--collapse-property-aggregation' must be none or mean.")
    keys = sorted(
        {
            key
            for leaf in leaves
            for key in leaf.props
            if key not in {'name', 'dist', 'support'}
        }
    )
    aggregated = {}
    status = {}
    for key in keys:
        values = [leaf.props.get(key) for leaf in leaves]
        present = [value for value in values if value not in (None, '')]
        if not present:
            status[key] = 'missing'
            continue
        if len(present) != len(values):
            aggregated[key] = 'partial'
            status[key] = 'partial'
            continue
        if all(str(value) == str(present[0]) for value in values):
            aggregated[key] = present[0]
            status[key] = 'constant'
            continue
        try:
            numeric = [float(value) for value in values]
        except (TypeError, ValueError):
            aggregated[key] = 'mixed'
            status[key] = 'mixed'
        else:
            finite = [value for value in numeric if math.isfinite(value)]
            if len(finite) == len(values) and numeric_mode == 'mean':
                aggregated[key] = sum(finite) / len(finite)
                status[key] = 'mean'
            else:
                aggregated[key] = 'mixed'
                status[key] = 'mixed'
    return aggregated, status


def _copy_tree_iteratively(tree):
    """Copy an ETE tree without serializing or recursing through its shape."""

    drawing_tree = Tree()
    clone_by_node = {tree: drawing_tree}
    stack = [tree]
    while stack:
        node = stack.pop()
        clone = clone_by_node[node]
        for key, value in node.props.items():
            clone.props[key] = value
        children = node.get_children()
        for child in children:
            child_clone = clone.add_child()
            clone_by_node[child] = child_clone
        stack.extend(reversed(children))
    return drawing_tree


def collapse_tree_for_drawing(
    tree,
    max_visible_tips=None,
    label_template=None,
    property_aggregation='none',
):
    """Return a copied tree whose frontier contains at most ``max_visible_tips``."""

    if max_visible_tips in (None, ''):
        return tree, []
    maximum = int(max_visible_tips)
    if maximum < 2:
        raise ValueError("'--max-visible-tips' must be at least 2.")
    if len(list(tree.leaves())) <= maximum:
        return tree, []
    drawing_tree = _copy_tree_iteratively(tree)
    leaf_count = {}
    for node in drawing_tree.traverse(strategy='postorder'):
        leaf_count[node] = 1 if node.is_leaf else sum(
            leaf_count[child] for child in node.get_children()
        )
    frontier = [drawing_tree]
    while True:
        candidates = []
        for index, node in enumerate(frontier):
            children = node.get_children()
            increase = len(children) - 1
            if children and len(frontier) + increase <= maximum:
                candidates.append((leaf_count[node], -index, index, node))
        if not candidates:
            break
        _, _, index, node = max(candidates)
        frontier[index:index + 1] = node.get_children()
    collapsed = [node for node in frontier if not node.is_leaf]
    template = label_template or '{first}…{last} (n={tips})'
    report = []
    for node in collapsed:
        leaves = list(node.leaves())
        first = str(leaves[0].name or '?')
        last = str(leaves[-1].name or '?')
        clade = str(node.name or first)
        label = str(template).format(
            clade=clade,
            first=first,
            last=last,
            tips=len(leaves),
        )
        aggregated, property_status = _aggregate_leaf_properties(
            leaves,
            numeric_mode=property_aggregation,
        )
        for child in list(node.get_children()):
            child.detach()
        node.name = label
        for key, value in aggregated.items():
            node.props[key] = value
        node.props['nwkit_collapsed'] = 'true'
        node.props['nwkit_collapsed_tips'] = len(leaves)
        report.append({
            'label': label,
            'first_tip': first,
            'last_tip': last,
            'tip_count': len(leaves),
            'property_summary': {
                key: {
                    'status': property_status[key],
                    'value': aggregated.get(key),
                }
                for key in sorted(property_status)
            },
        })
    return drawing_tree, report
