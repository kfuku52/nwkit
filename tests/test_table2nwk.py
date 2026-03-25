import pandas as pd
import pytest

from nwkit.nwk2table import nwk2table_main
from nwkit.table2nwk import table2nwk_main
from nwkit.util import read_tree
from tests.helpers import make_args


class TestTable2NwkMain:
    def test_roundtrip_support_tree(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,B:2)80:3,C:4);', 'tree.nwk')
        table_path = tmp_path / 'tree.tsv'
        outfile = tmp_path / 'roundtrip.nwk'
        nwk2table_main(make_args(infile=infile, outfile=str(table_path), age=False, sister=True))
        table2nwk_main(make_args(infile=str(table_path), outfile=str(outfile), outformat='auto'))
        tree = read_tree(str(outfile), format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.leaf_names()) == {'A', 'B', 'C'}
        assert abs(next(tree.search_nodes(name='A')).dist - 1.0) < 1e-6
        assert abs(next(tree.search_nodes(name='B')).dist - 2.0) < 1e-6
        assert abs(tree.common_ancestor(['A', 'B']).support - 80.0) < 1e-6

    def test_roundtrip_named_internal_nodes(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,B:2)AB:3,C:4)root;', 'tree.nwk')
        table_path = tmp_path / 'tree.tsv'
        outfile = tmp_path / 'roundtrip_named.nwk'
        nwk2table_main(make_args(infile=infile, outfile=str(table_path), age=False, sister=True))
        table2nwk_main(make_args(infile=str(table_path), outfile=str(outfile), outformat='auto'))
        tree = read_tree(str(outfile), format='auto', quoted_node_names=True, quiet=True)
        assert tree.name == 'root'
        assert tree.common_ancestor(['A', 'B']).name == 'AB'
        assert abs(next(tree.search_nodes(name='A')).dist - 1.0) < 1e-6

    def test_rejects_multiple_roots(self, tmp_path):
        table_path = tmp_path / 'invalid.tsv'
        pd.DataFrame(
            [
                {'branch_id': 0, 'parent': -1, 'name': 'root1'},
                {'branch_id': 1, 'parent': -1, 'name': 'root2'},
            ]
        ).to_csv(table_path, sep='\t', index=False)
        args = make_args(infile=str(table_path), outfile='-', outformat='auto')
        with pytest.raises(ValueError, match='exactly one root'):
            table2nwk_main(args)
