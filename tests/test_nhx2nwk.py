import os
import pytest
from ete3 import TreeNode

from nwkit.nhx2nwk import nhx2nwk_main
from nwkit.util import read_tree
from tests.helpers import make_args, DATA_DIR


class TestNhx2nwkMain:
    def test_no_node_label(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            node_label='',
        )
        nhx2nwk_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.get_leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_with_node_label_name_attr(self, tmp_nwk, tmp_outfile):
        # Test with the 'name' attribute (string type, avoids float issue with 'support')
        path = tmp_nwk('((A:1,B:1)AB:1,(C:1,D:1)CD:1)root;')
        args = make_args(
            infile=path, outfile=tmp_outfile, format='1',
            node_label='name',
        )
        nhx2nwk_main(args)
        assert os.path.exists(tmp_outfile)

    def test_with_data_file(self, tmp_outfile):
        infile = os.path.join(DATA_DIR, 'nhx2nwk_nodelabel', 'generax.nhx')
        if not os.path.exists(infile):
            pytest.skip('Test data not found')
        args = make_args(
            infile=infile, outfile=tmp_outfile,
            node_label='',
        )
        nhx2nwk_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert len(tree.get_leaf_names()) > 0

    def test_preserves_leaf_names(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((species_A:1,species_B:1):1,(species_C:1,species_D:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            node_label='',
        )
        nhx2nwk_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.get_leaf_names()) == {'species_A', 'species_B', 'species_C', 'species_D'}
