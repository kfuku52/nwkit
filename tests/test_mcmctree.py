import os
import pytest
from argparse import Namespace
from ete3 import TreeNode

from nwkit.mcmctree import (
    mcmctree_main,
    add_common_anc_constraint,
    is_mrca_clade_root,
    remove_constraint_equal_upper,
    apply_min_clade_prop,
)
from nwkit.util import read_tree
from tests.helpers import make_args, DATA_DIR


def make_mcmctree_args(**kwargs):
    defaults = {
        'infile': '-',
        'outfile': '-',
        'format': 'auto',
        'outformat': 'auto',
        'quoted_node_names': True,
        'timetree': 'no',
        'left_species': None,
        'right_species': None,
        'lower_bound': None,
        'upper_bound': None,
        'lower_tailProb': '0.025',
        'upper_tailProb': '0.025',
        'lower_offset': '0.1',
        'lower_scale': '1',
        'add_header': False,
        'min_clade_prop': 0,
        'higher_rank_search': True,
    }
    defaults.update(kwargs)
    return Namespace(**defaults)


class TestAddCommonAncConstraint:
    def test_bound_constraint(self):
        tree = TreeNode(newick='((a:1,b:1):1,(c:1,d:1):1);', format=1)
        args = make_mcmctree_args(
            left_species='a', right_species='b',
            lower_bound='10.0', upper_bound='20.0',
        )
        tree = add_common_anc_constraint(tree, args)
        common_anc = tree.get_common_ancestor('a', 'b')
        assert 'B(10.0, 20.0' in common_anc.name

    def test_point_constraint(self):
        tree = TreeNode(newick='((a:1,b:1):1,(c:1,d:1):1);', format=1)
        args = make_mcmctree_args(
            left_species='a', right_species='b',
            lower_bound='15.0', upper_bound='15.0',
        )
        tree = add_common_anc_constraint(tree, args)
        common_anc = tree.get_common_ancestor('a', 'b')
        assert '@15.0' in common_anc.name


class TestIsMrcaCladeRoot:
    def test_no_missing_ids_in_result(self):
        """If 'missing_ids' not in result, should return True."""
        tree = TreeNode(newick='((a:1,b:1):1,(c:1,d:1):1);', format=1)
        node = tree.get_common_ancestor('a', 'b')
        result = 'some_api_response_data'
        assert is_mrca_clade_root(node, result, ncbi=None) is True

    def test_empty_missing_ids(self):
        """If missing_ids list is empty, should return True."""
        tree = TreeNode(newick='((a:1,b:1):1,(c:1,d:1):1);', format=1)
        node = tree.get_common_ancestor('a', 'b')
        result = 'missing_ids:[]'
        assert is_mrca_clade_root(node, result, ncbi=None) is True

    def test_issue11_null_missing_ids_raises(self):
        """Regression test for GitHub issue #11.

        When the TimeTree API returns 'null' for missing_ids, the function
        crashes with ValueError: invalid literal for int() with base 10: 'null'.
        The actual fix is in add_timetree_constraint() which catches
        'No TimeTree study info available for this MRCA' before calling
        is_mrca_clade_root. This test documents the known behavior.
        """
        tree = TreeNode(newick='((a:1,b:1):1,(c:1,d:1):1);', format=1)
        node = tree.get_common_ancestor('a', 'b')
        result = 'missing_ids:[null]'
        # The function itself still raises ValueError on 'null' input;
        # the fix is that callers filter out such responses before reaching here
        with pytest.raises(ValueError, match="invalid literal for int"):
            is_mrca_clade_root(node, result, ncbi=None)


class TestRemoveConstraintEqualUpper:
    def test_removes_duplicate_constraints(self):
        tree = TreeNode(newick='(((a:1,b:1):1,c:1):1,(d:1,e:1):1);', format=1)
        # Set same constraint on parent and child
        ab_node = tree.get_common_ancestor('a', 'b')
        abc_node = tree.get_common_ancestor('a', 'b', 'c')
        ab_node.name = "'B(10.0, 20.0, 0.025, 0.025)'"
        abc_node.name = "'B(10.0, 20.0, 0.025, 0.025)'"
        tree = remove_constraint_equal_upper(tree)
        # Child constraint should be removed (it matches parent)
        assert ab_node.name == ''
        # Parent constraint should remain
        assert abc_node.name == "'B(10.0, 20.0, 0.025, 0.025)'"


class TestApplyMinCladeProp:
    def test_removes_small_clades(self):
        tree = TreeNode(newick='(((a:1,b:1):1,c:1):1,(d:1,e:1):1);', format=1)
        ab_node = tree.get_common_ancestor('a', 'b')
        ab_node.name = "'B(10.0, 20.0)'"
        # min_clade_prop=0.5 means clade must have >= 2.5 leaves (50% of 5)
        tree = apply_min_clade_prop(tree, min_clade_prop=0.5)
        # ab_node has only 2 leaves (< 2.5), so its constraint should be removed
        assert ab_node.name == ''

    def test_keeps_large_clades(self):
        tree = TreeNode(newick='(((a:1,b:1):1,c:1):1,(d:1,e:1):1);', format=1)
        abc_node = tree.get_common_ancestor('a', 'b', 'c')
        abc_node.name = "'B(10.0, 20.0)'"
        # min_clade_prop=0.5 means clade must have >= 2.5 leaves
        tree = apply_min_clade_prop(tree, min_clade_prop=0.5)
        # abc_node has 3 leaves (>= 2.5), so its constraint should remain
        assert abc_node.name == "'B(10.0, 20.0)'"


class TestMcmctreeMain:
    def test_basic_constraint(self, tmp_nwk, tmp_outfile):
        """Test mcmctree_main with --timetree no and manual bounds."""
        path = tmp_nwk('((a:1,b:1):1,(c:1,d:1):1);')
        args = make_mcmctree_args(
            infile=path, outfile=tmp_outfile,
            left_species='a', right_species='b',
            lower_bound='10.0', upper_bound='20.0',
        )
        mcmctree_main(args)
        with open(tmp_outfile) as f:
            content = f.read()
        assert 'B(10.0, 20.0' in content
        assert content.strip().endswith(';')

    def test_point_estimate(self, tmp_nwk, tmp_outfile):
        """Test mcmctree_main with point estimate (lower == upper)."""
        path = tmp_nwk('((a:1,b:1):1,(c:1,d:1):1);')
        args = make_mcmctree_args(
            infile=path, outfile=tmp_outfile,
            left_species='a', right_species='b',
            lower_bound='15.0', upper_bound='15.0',
        )
        mcmctree_main(args)
        with open(tmp_outfile) as f:
            content = f.read()
        assert '@15.0' in content

    def test_add_header(self, tmp_nwk, tmp_outfile):
        """Test mcmctree_main with --add_header."""
        path = tmp_nwk('((a:1,b:1):1,(c:1,d:1):1);')
        args = make_mcmctree_args(
            infile=path, outfile=tmp_outfile,
            left_species='a', right_species='b',
            lower_bound='10.0', upper_bound='20.0',
            add_header=True,
        )
        mcmctree_main(args)
        with open(tmp_outfile) as f:
            content = f.read()
        lines = content.strip().split('\n')
        # First line should be "4 1" (4 leaves, 1 tree)
        assert lines[0].strip() == '4 1'

    def test_issue7_quoted_node_names_in_output(self, tmp_nwk, tmp_outfile):
        """Regression test for GitHub issue #7.

        When mcmctree adds constraint annotations like 'B(152.3, 236.2)',
        the output must be valid and parseable. The quoted node names must
        be properly handled.
        """
        path = tmp_nwk('((((agl:1,(cma:1,dvi:1):1):1,sor:1):1,((tmo:1,zmo:1):1,tca:1):1):1,ota:1);')
        args = make_mcmctree_args(
            infile=path, outfile=tmp_outfile,
            left_species='tca', right_species='sor',
            lower_bound='152.3', upper_bound='236.2',
            lower_tailProb='1e-300', upper_tailProb='1e-300',
        )
        mcmctree_main(args)
        with open(tmp_outfile) as f:
            content = f.read().strip()
        # Output should contain the constraint
        assert 'B(152.3, 236.2, 1e-300, 1e-300)' in content
        # Output should be valid (end with ;)
        assert content.endswith(';')
        # Output should not contain 'NoName' (these should be cleaned up)
        assert 'NoName' not in content

    def test_issue7_chained_constraints(self, tmp_nwk, tmp_outfile):
        """Regression test for GitHub issue #7 (chained mcmctree commands).

        Simulate the piping scenario: first add one constraint, then read
        the output and add another constraint. Both should be present.
        """
        import tempfile
        path = tmp_nwk('((((agl:1,(cma:1,dvi:1):1):1,sor:1):1,((tmo:1,zmo:1):1,tca:1):1):1,ota:1);')

        # First constraint
        with tempfile.NamedTemporaryFile(mode='w', suffix='.nwk', delete=False) as f:
            intermediate = f.name
        try:
            args1 = make_mcmctree_args(
                infile=path, outfile=intermediate,
                left_species='tca', right_species='sor',
                lower_bound='152.3', upper_bound='236.2',
                lower_tailProb='1e-300', upper_tailProb='1e-300',
            )
            mcmctree_main(args1)

            # Second constraint on the intermediate output
            args2 = make_mcmctree_args(
                infile=intermediate, outfile=tmp_outfile,
                left_species='sor', right_species='agl',
                lower_bound='154.0', upper_bound='245.8',
            )
            mcmctree_main(args2)
        finally:
            os.unlink(intermediate)

        with open(tmp_outfile) as f:
            content = f.read().strip()
        assert content.endswith(';')
        # Both constraints should be present
        assert 'B(152.3, 236.2, 1e-300, 1e-300)' in content
        assert 'B(154.0, 245.8, 0.025, 0.025)' in content

    def test_unrooted_tree_raises(self, tmp_nwk):
        """mcmctree requires a rooted tree (2 children at root)."""
        path = tmp_nwk('(a:1,b:1,c:1);')
        args = make_mcmctree_args(
            infile=path, outfile='-',
            left_species='a', right_species='b',
            lower_bound='10.0', upper_bound='20.0',
        )
        with pytest.raises(AssertionError, match='rooted'):
            mcmctree_main(args)

    def test_with_data_files(self, tmp_outfile):
        """Test with data files if available."""
        infile = os.path.join(DATA_DIR, 'mcmctree1', 'input.nwk')
        if not os.path.exists(infile):
            pytest.skip('Test data not found')
        args = make_mcmctree_args(
            infile=infile, outfile=tmp_outfile,
            left_species=None, right_species=None,
            timetree='no',
        )
        # Without left/right species and timetree='no', this should fail
        # because add_common_anc_constraint needs left/right species
        # Just verify the tree is read correctly
        tree = read_tree(infile, 'auto', True, quiet=True)
        assert len(tree.get_leaf_names()) > 0
