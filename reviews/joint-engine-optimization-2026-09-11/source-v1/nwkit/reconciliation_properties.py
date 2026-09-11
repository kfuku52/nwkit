COLLAPSED_EVENT_BOUNDARY_PROP = "NWKIT_COLLAPSED_EVENT_BOUNDARY"


def property_is_true(value):
    return str(value).strip().upper() in {"Y", "YES", "TRUE", "1"}


def node_marks_non_speciation_boundary(node):
    if property_is_true(node.props.get(COLLAPSED_EVENT_BOUNDARY_PROP, "N")):
        return True
    transfer = str(node.props.get("H", "N")).strip().upper()
    if transfer == "Y" or transfer.startswith("Y@"):
        return True
    return property_is_true(node.props.get("D", "N"))


def preserve_collapsed_event_boundaries(tree):
    for node in tree.traverse(strategy="preorder"):
        if node.is_leaf or len(node.children) != 1:
            continue
        if node_marks_non_speciation_boundary(node):
            node.children[0].props[COLLAPSED_EVENT_BOUNDARY_PROP] = "Y"


def preserve_pruning_event_boundaries(tree, prune_flags):
    """Mark non-speciation events that pruning will collapse into one lineage."""
    for node in tree.traverse(strategy="preorder"):
        if node.is_leaf or prune_flags.get(node, False):
            continue
        retained_children = [
            child for child in node.children if not prune_flags.get(child, False)
        ]
        if len(retained_children) == 1 and node_marks_non_speciation_boundary(node):
            retained_children[0].props[COLLAPSED_EVENT_BOUNDARY_PROP] = "Y"
