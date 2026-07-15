import os

from nwkit.nhx2nwk import nhx2nwk_main
from nwkit.util import read_tree
from tests.helpers import make_args


class TestNhx2nwkMain:
    def test_no_node_label(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            node_label='',
        )
        nhx2nwk_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_with_node_label_name_attr(self, tmp_nwk, tmp_outfile):
        # Test with the 'name' attribute (string type, avoids float issue with 'support')
        path = tmp_nwk('((A:1,B:1)AB:1,(C:1,D:1)CD:1)root;')
        args = make_args(
            infile=path, outfile=tmp_outfile, format='1',
            node_label='name',
        )
        nhx2nwk_main(args)
        assert os.path.exists(tmp_outfile)

    def test_preserves_leaf_names(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((species_A:1,species_B:1):1,(species_C:1,species_D:1):1);')
        args = make_args(
            infile=path, outfile=tmp_outfile,
            node_label='',
        )
        nhx2nwk_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.leaf_names()) == {'species_A', 'species_B', 'species_C', 'species_D'}

    def test_node_label_with_auto_outformat_preserves_internal_names(self, tmp_nwk, tmp_outfile):
        nhx = '((A:1[&&NHX:S=spA],B:1[&&NHX:S=spB]):1[&&NHX:S=AB],(C:1[&&NHX:S=spC],D:1[&&NHX:S=spD]):1[&&NHX:S=CD])[&&NHX:S=root];'
        path = tmp_nwk(nhx)
        args = make_args(
            infile=path, outfile=tmp_outfile,
            node_label='S',
            format='auto', outformat='auto',
        )
        nhx2nwk_main(args)
        tree = read_tree(tmp_outfile, format='1', quoted_node_names=True, quiet=True)
        internal_names = {n.name for n in tree.traverse() if not n.is_leaf}
        assert {'AB', 'CD', 'root'}.issubset(internal_names)
