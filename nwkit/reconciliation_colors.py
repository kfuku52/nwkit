"""Shared branch colors and the internal edges changed by an X1 block."""

GENE_COLORS = ("#287d8e", "#bb7432", "#5e5fa8", "#478443")
X1_COLOR = "#a03faa"


def x1_internal_children(tree, blocks, field):
    """Identify edges within each consumed pair, excluding incoming/boundary edges."""
    return {
        node
        for block in blocks
        for node in block[field]
        if tree.parents.get(node) in block[field]
    }
