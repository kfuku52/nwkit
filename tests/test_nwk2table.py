import pandas as pd
import pytest

from nwkit.nwk2table import nwk2table_main
from tests.helpers import make_args


class TestNwk2TableMain:
    def test_writes_parent_child_table(self, tmp_nwk, tmp_path):
        path = tmp_nwk('((A:1,B:2)80:3,C:4);', 'tree.nwk')
        outfile = tmp_path / 'tree.tsv'
        args = make_args(
            infile=path,
            outfile=str(outfile),
            age=False,
            sister=True,
        )
        nwk2table_main(args)
        table = pd.read_csv(outfile, sep='\t')
        assert {'branch_id', 'parent', 'name', 'dist', 'support', 'sister'}.issubset(table.columns)
        assert len(table.index) == 5
        assert len(table.loc[table['parent'] == -1].index) == 1
        a_row = table.loc[table['name'] == 'A'].iloc[0]
        assert abs(a_row['dist'] - 1.0) < 1e-6
        support_values = table['support'].dropna().tolist()
        assert any(abs(value - 80.0) < 1e-6 for value in support_values)

    def test_adds_ages_for_ultrametric_tree(self, tmp_nwk, tmp_path):
        path = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);', 'tree.nwk')
        outfile = tmp_path / 'tree.tsv'
        args = make_args(
            infile=path,
            outfile=str(outfile),
            age=True,
            sister=False,
        )
        nwk2table_main(args)
        table = pd.read_csv(outfile, sep='\t')
        assert 'age' in table.columns
        root_age = table.loc[table['parent'] == -1, 'age'].iloc[0]
        assert abs(root_age - 2.0) < 1e-6

    def test_age_requires_ultrametric_tree(self, tmp_nwk, tmp_path):
        path = tmp_nwk('((A:1,B:1):1,C:1);', 'tree.nwk')
        outfile = tmp_path / 'tree.tsv'
        args = make_args(
            infile=path,
            outfile=str(outfile),
            age=True,
            sister=False,
        )
        with pytest.raises(ValueError, match='ultrametric'):
            nwk2table_main(args)
