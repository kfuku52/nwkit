import os
from argparse import Namespace

# Project root and data directory
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PROJECT_ROOT, 'data')


def safe_get_distance(tree, node1, node2):
    """Compute distance between nodes, treating None dist as 0.

    ete4's get_distance fails if any node along the path has None dist.
    This helper manually sums distances, treating None as 0.
    """
    # Find path from node1 to node2 via common ancestor
    ancestors1 = {node1}
    n = node1
    while n.up:
        n = n.up
        ancestors1.add(n)

    # Find common ancestor
    n = node2
    while n not in ancestors1:
        n = n.up
    common_anc = n

    # Sum distances from node1 to common ancestor
    dist = 0.0
    n = node1
    while n != common_anc:
        dist += n.dist if n.dist is not None else 0.0
        n = n.up

    # Sum distances from node2 to common ancestor
    n = node2
    while n != common_anc:
        dist += n.dist if n.dist is not None else 0.0
        n = n.up

    return dist


def make_args(**kwargs):
    """Create an argparse.Namespace with common defaults."""
    defaults = {
        'infile': '-',
        'outfile': '-',
        'format': 'auto',
        'outformat': 'auto',
        'quoted_node_names': True,
    }
    defaults.update(kwargs)
    return Namespace(**defaults)
