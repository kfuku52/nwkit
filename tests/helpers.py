import os
from argparse import Namespace

# Project root and data directory
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PROJECT_ROOT, 'data')


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
