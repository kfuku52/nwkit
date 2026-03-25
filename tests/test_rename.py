import pandas as pd
import pytest

from nwkit.rename import rename_main
from nwkit.util import read_tree
from tests.helpers import make_args


class TestRenameMain:
    def test_renames_leaf_nodes_from_tsv(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,B:1):1,(C:1,D:1):1);', 'tree.nwk')
        mapping_tsv = tmp_path / 'names.tsv'
        mapping_tsv.write_text('old_name\tnew_name\nA\tAlpha\nC\tCharlie\n')
        outfile = tmp_path / 'renamed.nwk'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            name_tsv=str(mapping_tsv),
            target='leaf',
            require_all_old_names=True,
            check_leaf_uniqueness=True,
        )
        rename_main(args)
        tree = read_tree(str(outfile), format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.leaf_names()) == {'Alpha', 'B', 'Charlie', 'D'}

    def test_can_rename_internal_nodes(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,B:1)AB:1,(C:1,D:1)CD:1)root;', 'tree.nwk')
        mapping_tsv = tmp_path / 'names.tsv'
        mapping_tsv.write_text('old_name\tnew_name\nAB\tclade_ab\nroot\tnew_root\n')
        outfile = tmp_path / 'renamed.nwk'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            name_tsv=str(mapping_tsv),
            target='intnode',
            require_all_old_names=True,
            check_leaf_uniqueness=True,
        )
        rename_main(args)
        tree = read_tree(str(outfile), format='auto', quoted_node_names=True, quiet=True)
        assert tree.name == 'new_root'
        assert tree.common_ancestor(['A', 'B']).name == 'clade_ab'

    def test_unused_old_name_raises_by_default(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,B:1):1,C:1);', 'tree.nwk')
        mapping_tsv = tmp_path / 'names.tsv'
        mapping_tsv.write_text('old_name\tnew_name\nZ\tZeta\n')
        args = make_args(
            infile=infile,
            outfile='-',
            name_tsv=str(mapping_tsv),
            target='leaf',
            require_all_old_names=True,
            check_leaf_uniqueness=True,
        )
        with pytest.raises(ValueError, match='were not found'):
            rename_main(args)

    def test_duplicate_leaf_names_after_rename_raise(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,B:1):1,C:1);', 'tree.nwk')
        mapping_tsv = tmp_path / 'names.tsv'
        mapping_tsv.write_text('old_name\tnew_name\nA\tX\nB\tX\n')
        args = make_args(
            infile=infile,
            outfile='-',
            name_tsv=str(mapping_tsv),
            target='leaf',
            require_all_old_names=True,
            check_leaf_uniqueness=True,
        )
        with pytest.raises(ValueError, match='Duplicated leaf labels'):
            rename_main(args)

    def test_duplicate_old_name_in_target_is_ambiguous(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,A:2):1,B:1);', 'tree.nwk')
        mapping_tsv = tmp_path / 'names.tsv'
        mapping_tsv.write_text('old_name\tnew_name\nA\tX\n')
        args = make_args(
            infile=infile,
            outfile='-',
            name_tsv=str(mapping_tsv),
            target='leaf',
            require_all_old_names=True,
            check_leaf_uniqueness=False,
        )
        with pytest.raises(ValueError, match='Matched multiple target nodes'):
            rename_main(args)

    def test_regex_rename_updates_matching_nodes(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((sample_A:1,sample_B:1):1,C:1);', 'tree.nwk')
        outfile = tmp_path / 'renamed_regex.nwk'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            name_tsv=None,
            pattern=r'^sample_',
            replacement='',
            target='leaf',
            require_match=True,
            check_leaf_uniqueness=True,
        )
        rename_main(args)
        tree = read_tree(str(outfile), format='auto', quoted_node_names=True, quiet=True)
        assert set(tree.leaf_names()) == {'A', 'B', 'C'}

    def test_regex_rename_requires_a_match_by_default(self, tmp_nwk):
        infile = tmp_nwk('((A:1,B:1):1,C:1);', 'tree.nwk')
        args = make_args(
            infile=infile,
            outfile='-',
            name_tsv=None,
            pattern=r'^sample_',
            replacement='',
            target='leaf',
            require_match=True,
            check_leaf_uniqueness=True,
        )
        with pytest.raises(ValueError, match='No target node names matched'):
            rename_main(args)

    def test_leaf_rename_preserves_support_values_in_auto_mode(self, tmp_nwk, tmp_path):
        infile = tmp_nwk('((A:1,B:1)0.95:1,(C:1,D:1)80:1);', 'tree.nwk')
        outfile = tmp_path / 'renamed_with_support.nwk'
        args = make_args(
            infile=infile,
            outfile=str(outfile),
            name_tsv=None,
            pattern=r'^A$',
            replacement='A1',
            target='leaf',
            require_match=True,
            check_leaf_uniqueness=True,
            outformat='auto',
            format='0',
        )
        rename_main(args)
        tree = read_tree(str(outfile), format='0', quoted_node_names=True, quiet=True)
        assert abs(tree.common_ancestor(['A1', 'B']).support - 0.95) < 1e-9
        assert abs(tree.common_ancestor(['C', 'D']).support - 80.0) < 1e-9
