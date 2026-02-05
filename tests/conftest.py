import pytest
from ete3 import TreeNode
from tests.helpers import DATA_DIR


@pytest.fixture
def simple_tree():
    """A simple rooted tree: ((A:1,B:1):1,(C:1,D:1):1);"""
    return TreeNode(newick='((A:1,B:1):1,(C:1,D:1):1);', format=1)


@pytest.fixture
def six_leaf_tree():
    """Tree with 6 leaves: (((A:1,B:1):1,(C:1,D:1):1):1,(E:1,F:1):1):0;"""
    return TreeNode(newick='(((A:1,B:1):1,(C:1,D:1):1):1,(E:1,F:1):1):0;', format=1)


@pytest.fixture
def named_internal_tree():
    """Tree with named internal nodes."""
    return TreeNode(newick='(((A1:2,(B1:1,B2:1):1):1,(A2:1,C1:1):2):1,C2:4):0.25;', format=1)


@pytest.fixture
def unrooted_tree():
    """An unrooted tree (3 children at root): (A:1,B:1,C:1);"""
    return TreeNode(newick='(A:1,B:1,C:1);', format=1)


@pytest.fixture
def species_tree():
    """Tree with GENUS_SPECIES_GENEID leaf names."""
    nwk = '((Homo_sapiens_GENE1:1,Homo_sapiens_GENE2:1):1,(Mus_musculus_GENE1:1,Danio_rerio_GENE1:1):1);'
    return TreeNode(newick=nwk, format=1)


@pytest.fixture
def data_dir():
    """Path to the test data directory."""
    return DATA_DIR


@pytest.fixture
def tmp_nwk(tmp_path):
    """Helper to create a temporary Newick file. Returns a factory function."""
    def _make(nwk_str, name='tree.nwk'):
        p = tmp_path / name
        p.write_text(nwk_str)
        return str(p)
    return _make


@pytest.fixture
def tmp_outfile(tmp_path):
    """Return a path for a temporary output file."""
    return str(tmp_path / 'output.nwk')
