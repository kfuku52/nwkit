import pandas as pd

from nwkit.validate import validate_main
from tests.helpers import make_args


def _write_tree_collection(tmp_path, trees, name='trees.nwk'):
    path = tmp_path / name
    path.write_text('\n'.join(trees) + '\n')
    return str(path)


class TestValidateMain:
    def test_reports_duplicate_leaf_and_negative_branch_length(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:-1,A:1):1,B:1);', 'tree.nwk')
        outfile = tmp_path / 'validate.tsv'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            check_species=False,
            require_rooted=False,
            require_ultrametric=False,
            require_same_leaf_set=True,
            fail_on_issue=False,
        )
        validate_main(args)
        table = pd.read_csv(outfile, sep='\t')
        assert table.loc[0, 'status'] == 'invalid'
        assert 'duplicate_leaf_names' in table.loc[0, 'issues']
        assert 'negative_branch_length' in table.loc[0, 'issues']

    def test_reports_leaf_set_mismatch_across_tree_collection(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            [
                '((A:1,B:1):1,(C:1,D:1):1);',
                '((A:1,B:1):1,(C:1,E:1):1);',
            ],
        )
        outfile = tmp_path / 'validate.tsv'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            check_species=False,
            require_rooted=False,
            require_ultrametric=False,
            require_same_leaf_set=True,
            fail_on_issue=False,
        )
        validate_main(args)
        table = pd.read_csv(outfile, sep='\t')
        assert table.loc[0, 'status'] == 'ok'
        assert table.loc[1, 'status'] == 'invalid'
        assert 'leaf_set_mismatch' in table.loc[1, 'issues']

    def test_require_same_rooting_compares_root_bipartition(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            [
                '((A:1,B:1):1,(C:1,(D:1,E:1):1):1);',
                '(A:0.5,(B:1,(C:1,(D:1,E:1):1):2):0.5);',
                '(((E:1,D:1):1,C:1):1,(B:1,A:1):1);',
            ],
        )
        outfile = tmp_path / 'validate.tsv'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            check_species=False,
            require_rooted=False,
            require_ultrametric=False,
            require_same_leaf_set=True,
            require_same_rooting=True,
            fail_on_issue=False,
        )

        validate_main(args)

        table = pd.read_csv(outfile, sep='\t')
        assert bool(table.loc[0, 'rooting_matches_first']) is True
        assert bool(table.loc[1, 'rooting_matches_first']) is False
        assert table.loc[1, 'status'] == 'invalid'
        assert 'rooting_mismatch' in table.loc[1, 'issues']
        assert bool(table.loc[2, 'rooting_matches_first']) is True
        assert table.loc[2, 'status'] == 'ok'

    def test_species_check_reports_non_monophyletic_species(self, tmp_nwk, tmp_path):
        infile = tmp_nwk(
            '((Homo_sapiens_gene1:1,Pan_troglodytes_gene1:1):1,(Homo_sapiens_gene2:1,Mus_musculus_gene1:1):1);',
            'tree.nwk',
        )
        outfile = tmp_path / 'validate.tsv'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            check_species=True,
            require_rooted=False,
            require_ultrametric=False,
            require_same_leaf_set=True,
            fail_on_issue=False,
        )
        validate_main(args)
        table = pd.read_csv(outfile, sep='\t')
        assert bool(table.loc[0, 'species_parseable']) is True
        assert bool(table.loc[0, 'species_groups_monophyletic']) is False
        assert 'species_not_monophyletic' in table.loc[0, 'issues']

    def test_can_report_format_ambiguity_and_support_requirements(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,B:1)40:1,(C:1,D:1):1);', 'tree.nwk')
        outfile = tmp_path / 'validate.tsv'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            check_species=False,
            require_rooted=False,
            require_ultrametric=False,
            require_same_leaf_set=True,
            require_same_rooting=False,
            require_binary=False,
            require_all_support=True,
            require_unambiguous_format=True,
            require_unquoted_names=False,
            fail_on_issue=False,
        )
        validate_main(args)
        table = pd.read_csv(outfile, sep='\t')
        assert bool(table.loc[0, 'format_ambiguous']) is True
        assert 'format_ambiguous' in table.loc[0, 'issues']
        assert 'missing_support' in table.loc[0, 'issues']

    def test_parse_errors_are_reported_instead_of_aborting(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            [
                '((A:1,B:1):1,(C:1,D:1):1);',
                '((A:1,B:1):1,(C:1,D:1):1));',
            ],
            name='invalid_collection.nwk',
        )
        outfile = tmp_path / 'validate.tsv'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            check_species=False,
            require_rooted=False,
            require_ultrametric=False,
            require_same_leaf_set=True,
            require_same_rooting=False,
            require_binary=False,
            require_all_support=False,
            require_unambiguous_format=False,
            require_unquoted_names=False,
            fail_on_issue=False,
        )
        validate_main(args)
        table = pd.read_csv(outfile, sep='\t')
        assert table.loc[0, 'status'] == 'ok'
        assert table.loc[1, 'status'] == 'invalid'
        assert bool(table.loc[1, 'parse_ok']) is False
        assert 'parse_error' in table.loc[1, 'issues']

    def test_malformed_first_tree_leaves_first_tree_comparisons_blank(self, tmp_path):
        infile = _write_tree_collection(
            tmp_path,
            [
                '((A:1,B:1):1,(C:1,D:1):1));',
                '((A:1,B:1):1,(C:1,D:1):1);',
            ],
            name='first_invalid_collection.nwk',
        )
        outfile = tmp_path / 'validate.tsv'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            check_species=False,
            require_rooted=False,
            require_ultrametric=False,
            require_same_leaf_set=True,
            require_same_rooting=True,
            require_binary=False,
            require_all_support=False,
            require_unambiguous_format=False,
            require_unquoted_names=False,
            fail_on_issue=False,
        )
        validate_main(args)
        table = pd.read_csv(outfile, sep='\t')
        assert table.loc[0, 'status'] == 'invalid'
        assert table.loc[1, 'status'] == 'ok'
        assert pd.isna(table.loc[1, 'leaf_set_matches_first'])
        assert pd.isna(table.loc[1, 'rooting_matches_first'])
        issues = '' if pd.isna(table.loc[1, 'issues']) else str(table.loc[1, 'issues'])
        assert 'leaf_set_mismatch' not in issues
        assert 'rooting_mismatch' not in issues
