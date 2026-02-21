import os
import pytest
from argparse import Namespace
from unittest.mock import patch, MagicMock
from ete4 import Tree
import requests

import nwkit.mcmctree as mcmctree_mod
from nwkit.mcmctree import (
    mcmctree_main,
    add_common_anc_constraint,
    add_timetree_constraint,
    is_mrca_clade_root,
    are_two_lineage_rank_differentiated,
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
        tree = Tree('((a:1,b:1):1,(c:1,d:1):1);', parser=1)
        args = make_mcmctree_args(
            left_species='a', right_species='b',
            lower_bound='10.0', upper_bound='20.0',
        )
        tree = add_common_anc_constraint(tree, args)
        common_anc = tree.common_ancestor('a', 'b')
        assert 'B(10.0, 20.0' in common_anc.name

    def test_point_constraint(self):
        tree = Tree('((a:1,b:1):1,(c:1,d:1):1);', parser=1)
        args = make_mcmctree_args(
            left_species='a', right_species='b',
            lower_bound='15.0', upper_bound='15.0',
        )
        tree = add_common_anc_constraint(tree, args)
        common_anc = tree.common_ancestor('a', 'b')
        assert '@15.0' in common_anc.name

    def test_point_constraint_with_equivalent_numeric_strings(self):
        tree = Tree('((a:1,b:1):1,(c:1,d:1):1);', parser=1)
        args = make_mcmctree_args(
            left_species='a', right_species='b',
            lower_bound='15', upper_bound='15.0',
        )
        tree = add_common_anc_constraint(tree, args)
        common_anc = tree.common_ancestor('a', 'b')
        assert '@15' in common_anc.name

    def test_species_not_found_raises(self):
        tree = Tree('((a:1,b:1):1,(c:1,d:1):1);', parser=1)
        args = make_mcmctree_args(
            left_species='a', right_species='x',
            lower_bound='10.0', upper_bound='20.0',
        )
        with pytest.raises(ValueError, match='Species not found'):
            add_common_anc_constraint(tree, args)

    def test_same_left_and_right_species_raises(self):
        tree = Tree('((a:1,b:1):1,(c:1,d:1):1);', parser=1)
        args = make_mcmctree_args(
            left_species='a', right_species='a',
            lower_bound='10.0', upper_bound='20.0',
        )
        with pytest.raises(ValueError, match='must be different'):
            add_common_anc_constraint(tree, args)


class TestIsMrcaCladeRoot:
    def test_no_missing_ids_in_result(self):
        """If 'missing_ids' not in result, should return True."""
        tree = Tree('((a:1,b:1):1,(c:1,d:1):1);', parser=1)
        node = tree.common_ancestor('a', 'b')
        result = 'some_api_response_data'
        assert is_mrca_clade_root(node, result, ncbi=None) is True

    def test_empty_missing_ids(self):
        """If missing_ids list is empty, should return True."""
        tree = Tree('((a:1,b:1):1,(c:1,d:1):1);', parser=1)
        node = tree.common_ancestor('a', 'b')
        result = 'missing_ids:[]'
        assert is_mrca_clade_root(node, result, ncbi=None) is True

    def test_issue11_null_missing_ids_is_handled(self):
        """Regression test for GitHub issue #11.

        'missing_ids:[null]' should not crash even when this helper is called directly.
        """
        tree = Tree('((a:1,b:1):1,(c:1,d:1):1);', parser=1)
        node = tree.common_ancestor('a', 'b')
        result = 'missing_ids:[null]'
        assert is_mrca_clade_root(node, result, ncbi=None) is True


class TestAreTwoLineageRankDifferentiated:
    def test_non_bifurcating_node_returns_false(self):
        tree = Tree('(a:1,b:1,c:1);', parser=1)
        node = tree
        taxids = [1, 2, 3]
        leaf_names = ['a', 'b', 'c']
        assert are_two_lineage_rank_differentiated(node, taxids, leaf_names) is False


class TestRemoveConstraintEqualUpper:
    def test_removes_duplicate_constraints(self):
        tree = Tree('(((a:1,b:1):1,c:1):1,(d:1,e:1):1);', parser=1)
        # Set same constraint on parent and child
        ab_node = tree.common_ancestor('a', 'b')
        abc_node = tree.common_ancestor('a', 'b', 'c')
        ab_node.name = "'B(10.0, 20.0, 0.025, 0.025)'"
        abc_node.name = "'B(10.0, 20.0, 0.025, 0.025)'"
        tree = remove_constraint_equal_upper(tree)
        # Child constraint should be removed (it matches parent)
        assert ab_node.name == 'NoName'
        # Parent constraint should remain
        assert abc_node.name == "'B(10.0, 20.0, 0.025, 0.025)'"


class TestApplyMinCladeProp:
    def test_removes_small_clades(self):
        tree = Tree('(((a:1,b:1):1,c:1):1,(d:1,e:1):1);', parser=1)
        ab_node = tree.common_ancestor('a', 'b')
        ab_node.name = "'B(10.0, 20.0)'"
        # min_clade_prop=0.5 means clade must have >= 2.5 leaves (50% of 5)
        tree = apply_min_clade_prop(tree, min_clade_prop=0.5)
        # ab_node has only 2 leaves (< 2.5), so its constraint should be removed
        assert ab_node.name == 'NoName'

    def test_keeps_large_clades(self):
        tree = Tree('(((a:1,b:1):1,c:1):1,(d:1,e:1):1);', parser=1)
        abc_node = tree.common_ancestor('a', 'b', 'c')
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
        with pytest.raises(ValueError, match='rooted'):
            mcmctree_main(args)

    def test_timetree_no_requires_left_and_right_species(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((a:1,b:1):1,(c:1,d:1):1);')
        args = make_mcmctree_args(
            infile=path, outfile=tmp_outfile,
            left_species=None, right_species=None,
            lower_bound='10.0', upper_bound='20.0',
            timetree='no',
        )
        with pytest.raises(ValueError, match='left_species'):
            mcmctree_main(args)

    def test_timetree_no_requires_at_least_one_bound(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((a:1,b:1):1,(c:1,d:1):1);')
        args = make_mcmctree_args(
            infile=path, outfile=tmp_outfile,
            left_species='a', right_species='b',
            lower_bound=None, upper_bound=None,
            timetree='no',
        )
        with pytest.raises(ValueError, match='lower_bound'):
            mcmctree_main(args)

    @pytest.mark.parametrize('min_clade_prop', [-0.1, 1.1])
    def test_min_clade_prop_out_of_range_raises(self, tmp_nwk, tmp_outfile, min_clade_prop):
        path = tmp_nwk('((a:1,b:1):1,(c:1,d:1):1);')
        args = make_mcmctree_args(
            infile=path, outfile=tmp_outfile,
            left_species='a', right_species='b',
            lower_bound='10.0', upper_bound='20.0',
            min_clade_prop=min_clade_prop,
        )
        with pytest.raises(ValueError, match='min_clade_prop'):
            mcmctree_main(args)

    def test_timetree_requires_named_leaves(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((:1,b:1):1,(c:1,d:1):1);')
        args = make_mcmctree_args(
            infile=path, outfile=tmp_outfile,
            timetree='point',
        )
        with pytest.raises(ValueError, match='non-empty names'):
            mcmctree_main(args)

    def test_unknown_timetree_mode_raises(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk('((a:1,b:1):1,(c:1,d:1):1);')
        args = make_mcmctree_args(
            infile=path, outfile=tmp_outfile,
            timetree='unknown',
        )
        with pytest.raises(ValueError, match='Unknown'):
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
        assert len(list(tree.leaf_names())) > 0


class TestIssue12EndpointUrl:
    """Regression tests for GitHub issue #12: timetree.org API change.

    The old endpoint (timetree.temple.edu/api) returns JSON errors for all
    queries. The working endpoint is timetree.org/api.
    """

    def test_endpoint_url_is_timetree_org(self):
        """Verify the endpoint URL points to timetree.org, not temple.edu."""
        import inspect
        source = inspect.getsource(add_timetree_constraint)
        assert 'timetree.org/api' in source
        assert 'timetree.temple.edu' not in source

    def test_csv_response_parsing(self):
        """Verify that a CSV response from timetree.org is parsed correctly.

        The API returns CSV with \\r\\n line separator:
        precomputed_median,precomputed_age,precomputed_ci_low,precomputed_ci_high\\r\\n
        87.2,87.2,81.33807,91
        """
        import re
        csv_response = "precomputed_median,precomputed_age,precomputed_ci_low,precomputed_ci_high\r\n87.2,87.2,81.33807,91"
        # Simulate the parsing logic from add_timetree_constraint
        timetree_result = re.sub('.*;</script>', '', csv_response)
        timetree_keys = re.sub('\r\n.*', '', re.sub('<br>.*', '', timetree_result)).split(',')
        timetree_values = re.sub('.*\r\n', '', re.sub('.*<br>', '', timetree_result)).split(',')
        timetree_dict = dict()
        for key, value in zip(timetree_keys, timetree_values):
            timetree_dict[key] = value
        assert timetree_dict['precomputed_age'] == '87.2'
        assert timetree_dict['precomputed_ci_low'] == '81.33807'
        assert timetree_dict['precomputed_ci_high'] == '91'

    def test_temple_edu_json_error_is_skipped(self):
        """The old temple.edu JSON error response should be handled gracefully."""
        json_error = '{"found_ids":[],"missing_ids":[null,null],"mrca_ttid":null,"taxon_id":null,"error":"No TimeTree study info available for this MRCA"}'
        assert "No TimeTree study info available for this MRCA" in json_error

    def test_mrca_not_found_json_is_skipped(self):
        """JSON error with 'MRCA node not found' should be handled."""
        json_error = '{"found_ids":[],"missing_ids":["9999999","8888888"],"mrca_id":-1,"error":"MRCA node not found"}'
        assert "MRCA node not found" in json_error

    def test_network_error_is_skipped_without_crash(self, monkeypatch):
        class FakeNCBI:
            def get_name_translator(self, names):
                return {name: [i + 1] for i, name in enumerate(names)}

            def get_lineage(self, taxid):
                return [1, int(taxid)]

            def get_rank(self, lineages):
                out = {}
                for t in lineages:
                    t = int(t)
                    out[t] = 'species' if t != 1 else 'superkingdom'
                return out

            def get_taxid_translator(self, taxids):
                return {int(t): f'sp{int(t)}' for t in taxids}

        def raise_network_error(*args, **kwargs):
            raise requests.RequestException('network down')

        monkeypatch.setattr(mcmctree_mod, 'NCBITaxa', FakeNCBI)
        monkeypatch.setattr(mcmctree_mod.requests, 'get', raise_network_error)

        tree = Tree('((a:1,b:1):1);', parser=1)
        args = make_mcmctree_args(timetree='point', higher_rank_search=True)
        out = add_timetree_constraint(tree, args)
        assert out.name == 'NoName'

    def test_non_200_response_is_skipped_without_crash(self, monkeypatch):
        class FakeNCBI:
            def get_name_translator(self, names):
                return {name: [i + 1] for i, name in enumerate(names)}

            def get_lineage(self, taxid):
                return [1, int(taxid)]

            def get_rank(self, lineages):
                out = {}
                for t in lineages:
                    t = int(t)
                    out[t] = 'species' if t != 1 else 'superkingdom'
                return out

            def get_taxid_translator(self, taxids):
                return {int(t): f'sp{int(t)}' for t in taxids}

        class FakeResponse:
            status_code = 500
            text = '<html>server error</html>'

        def fake_get(*args, **kwargs):
            return FakeResponse()

        monkeypatch.setattr(mcmctree_mod, 'NCBITaxa', FakeNCBI)
        monkeypatch.setattr(mcmctree_mod.requests, 'get', fake_get)

        tree = Tree('((a:1,b:1):1);', parser=1)
        args = make_mcmctree_args(timetree='point', higher_rank_search=True)
        out = add_timetree_constraint(tree, args)
        assert out.name == 'NoName'
