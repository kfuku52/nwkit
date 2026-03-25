import pytest

from nwkit.collapse import collapse_main
from nwkit.util import read_tree
from tests.helpers import make_args


class TestCollapseMain:
    def test_collapses_by_support(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:1)40:1,(C:1,D:1)90:1);', 'tree.nwk')
        args = make_args(
            infile=path,
            outfile=tmp_outfile,
            min_support=50.0,
            max_dist=None,
            preserve_branch_length=True,
        )
        collapse_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        child_leaf_sets = {frozenset(child.leaf_names()) for child in tree.get_children()}
        assert child_leaf_sets == {frozenset({'A'}), frozenset({'B'}), frozenset({'C', 'D'})}
        assert abs(next(tree.search_nodes(name='A')).dist - 2.0) < 1e-6
        assert abs(next(tree.search_nodes(name='B')).dist - 2.0) < 1e-6

    def test_collapses_by_branch_length(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((A:1,B:1):0.1,C:1);', 'tree.nwk')
        args = make_args(
            infile=path,
            outfile=tmp_outfile,
            min_support=None,
            max_dist=0.1,
            preserve_branch_length=True,
        )
        collapse_main(args)
        tree = read_tree(tmp_outfile, format='auto', quoted_node_names=True, quiet=True)
        child_leaf_sets = {frozenset(child.leaf_names()) for child in tree.get_children()}
        assert child_leaf_sets == {frozenset({'A'}), frozenset({'B'}), frozenset({'C'})}
        assert abs(next(tree.search_nodes(name='A')).dist - 1.1) < 1e-6
        assert abs(next(tree.search_nodes(name='B')).dist - 1.1) < 1e-6

    def test_requires_at_least_one_threshold(self, tmp_nwk):
        path = tmp_nwk('((A:1,B:1):1,C:1);', 'tree.nwk')
        args = make_args(
            infile=path,
            outfile='-',
            min_support=None,
            max_dist=None,
            preserve_branch_length=True,
        )
        with pytest.raises(ValueError, match='Specify at least one'):
            collapse_main(args)

    def test_support_threshold_requires_meaningful_support(self, tmp_nwk):
        path = tmp_nwk('((A:1,B:1)AB:1,C:1)root;', 'tree.nwk')
        args = make_args(
            infile=path,
            outfile='-',
            min_support=50.0,
            max_dist=None,
            preserve_branch_length=True,
        )
        with pytest.raises(ValueError, match='No meaningful support values'):
            collapse_main(args)
