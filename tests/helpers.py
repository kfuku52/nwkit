from argparse import Namespace
from ete4 import Tree


def make_deep_ladder_tree(tip_count):
    root = Tree()
    current = root
    for index in range(tip_count - 2):
        current.add_child(name='T{}'.format(index), dist=1)
        current = current.add_child(dist=1)
    current.add_child(name='T{}'.format(tip_count - 2), dist=1)
    current.add_child(name='T{}'.format(tip_count - 1), dist=1)
    return root

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
        'download_dir': 'auto',
        'species_parser': 'legacy',
        'species_regex': r'^([^_]+_[^_]+)(?:_|$)',
        'species_map_tsv': None,
    }
    defaults.update(kwargs)
    return Namespace(**defaults)


def write_tree_collection(tmp_path, trees, name='trees.nwk'):
    """Write one Newick tree per line and return the resulting path."""
    path = tmp_path / name
    path.write_text('\n'.join(trees) + '\n')
    return str(path)
