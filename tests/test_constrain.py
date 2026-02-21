import pandas
import pytest
from argparse import Namespace
from ete4 import Tree

from nwkit.constrain import (
    get_max_ancestor_overlap_node,
    check_input_file,
    constrain_main,
    collapse_genes,
    get_taxid_counts,
    initialize_tree,
    match_taxa,
    taxid2tree,
)
from tests.helpers import make_args


class TestGetMaxAncestorOverlapNode:
    def test_returns_none_when_no_overlap(self):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        for node in tree.traverse():
            node.add_props(ancestors=[1, 2])
        assert get_max_ancestor_overlap_node(tree, [9, 10]) is None

    def test_uses_node_level_ancestors_for_max_overlap(self):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        for node in tree.traverse():
            node.add_props(ancestors=[])
        left = tree.common_ancestor(['A', 'B'])
        right = tree.common_ancestor(['C', 'D'])
        left.add_props(ancestors=[10, 11])
        right.add_props(ancestors=[10, 11, 12])
        result = get_max_ancestor_overlap_node(tree, [10, 12])
        assert result is right


class TestTaxid2Tree:
    def test_single_species_returns_single_leaf_tree(self):
        lineages = {'Homo_sapiens_gene1': [1, 2759, 7711]}
        taxid_counts = get_taxid_counts(lineages)
        tree = taxid2tree(lineages, taxid_counts)
        assert set(tree.leaf_names()) == {'Homo_sapiens_gene1'}

    def test_empty_lineages_raise(self):
        with pytest.raises(ValueError, match='No valid taxa'):
            taxid2tree({}, get_taxid_counts({}))


class TestCheckInputFile:
    def test_empty_species_list_raises(self, tmp_path):
        species_path = tmp_path / 'species.txt'
        species_path.write_text('\n')
        args = Namespace(species_list=str(species_path), taxid_tsv=None, backbone='ncbi')
        with pytest.raises(ValueError, match='species_list is empty'):
            check_input_file(args)

    def test_duplicate_species_list_entries_raise(self, tmp_path):
        species_path = tmp_path / 'species.txt'
        species_path.write_text('A\nA\n')
        args = Namespace(species_list=str(species_path), taxid_tsv=None, backbone='ncbi')
        with pytest.raises(ValueError, match='Duplicate entries'):
            check_input_file(args)

    def test_taxid_tsv_duplicate_leaf_name_raises(self, tmp_path):
        tsv_path = tmp_path / 'taxid.tsv'
        pandas.DataFrame(
            {'leaf_name': ['A', 'A'], 'taxid': [9606, 9606]}
        ).to_csv(tsv_path, sep='\t', index=False)
        args = Namespace(species_list=None, taxid_tsv=str(tsv_path), backbone='ncbi')
        with pytest.raises(ValueError, match='Duplicate values'):
            check_input_file(args)

    def test_taxid_tsv_missing_taxid_raises(self, tmp_path):
        tsv_path = tmp_path / 'taxid.tsv'
        pandas.DataFrame(
            {'leaf_name': ['A', 'B'], 'taxid': [9606, None]}
        ).to_csv(tsv_path, sep='\t', index=False)
        args = Namespace(species_list=None, taxid_tsv=str(tsv_path), backbone='ncbi')
        with pytest.raises(ValueError, match='missing values'):
            check_input_file(args)

    def test_taxid_tsv_non_numeric_taxid_raises(self, tmp_path):
        tsv_path = tmp_path / 'taxid.tsv'
        pandas.DataFrame(
            {'leaf_name': ['A', 'B'], 'taxid': ['9606', 'abc']}
        ).to_csv(tsv_path, sep='\t', index=False)
        args = Namespace(species_list=None, taxid_tsv=str(tsv_path), backbone='ncbi')
        with pytest.raises(ValueError, match='non-numeric'):
            check_input_file(args)

    def test_taxid_tsv_non_integer_taxid_raises(self, tmp_path):
        tsv_path = tmp_path / 'taxid.tsv'
        pandas.DataFrame(
            {'leaf_name': ['A', 'B'], 'taxid': [9606, 123.5]}
        ).to_csv(tsv_path, sep='\t', index=False)
        args = Namespace(species_list=None, taxid_tsv=str(tsv_path), backbone='ncbi')
        with pytest.raises(ValueError, match='non-integer'):
            check_input_file(args)


class TestCollapseGenes:
    def test_invalid_leaf_name_format_raises(self):
        tree = Tree('((A:1,Homo_sapiens_gene1:1):1,B:1);', parser=1)
        with pytest.raises(ValueError, match='GENUS_SPECIES'):
            collapse_genes(tree)

    def test_unnamed_leaf_raises_value_error(self):
        tree = Tree('((:1,Homo_sapiens_gene1:1):1,Mus_musculus_gene1:1);', parser=1)
        with pytest.raises(ValueError, match='GENUS_SPECIES'):
            collapse_genes(tree)


class TestMatchTaxa:
    def test_mixed_label_formats_do_not_depend_on_first_label(self):
        tree = Tree('(Homo sapiens:1,Mus musculus:1);', parser=1)
        tree = initialize_tree(tree)
        labels = ['UnknownLabel', 'Homo_sapiens_gene1']
        out = match_taxa(tree=tree, labels=labels, backbone_method='user')
        homo = [leaf for leaf in out.leaves() if leaf.name == 'Homo sapiens'][0]
        mus = [leaf for leaf in out.leaves() if leaf.name == 'Mus musculus'][0]
        assert homo.props.get('has_taxon') is True
        assert homo.props.get('taxon_names') == ['Homo_sapiens_gene1']
        assert mus.props.get('has_taxon') is False


class TestConstrainMain:
    def test_user_backbone_no_match_raises_clear_error(self, tmp_path):
        species_path = tmp_path / 'species.txt'
        species_path.write_text('Pan_troglodytes_gene1\n')
        tree_path = tmp_path / 'tree.nwk'
        tree_path.write_text('(Homo_sapiens:1,Mus_musculus:1);\n')
        args = make_args(
            infile=str(tree_path),
            outfile='-',
            species_list=str(species_path),
            taxid_tsv=None,
            backbone='user',
            rank='no',
            collapse=False,
        )
        with pytest.raises(ValueError, match='No taxa from --species_list matched'):
            constrain_main(args)

    def test_user_backbone_duplicate_leaf_labels_raise_clear_error(self, tmp_path):
        species_path = tmp_path / 'species.txt'
        species_path.write_text('A\n')
        tree_path = tmp_path / 'tree.nwk'
        tree_path.write_text('((A:1,A:1):1,B:1);\n')
        args = make_args(
            infile=str(tree_path),
            outfile='-',
            species_list=str(species_path),
            taxid_tsv=None,
            backbone='user',
            rank='no',
            collapse=False,
        )
        with pytest.raises(ValueError, match='Duplicated leaf labels'):
            constrain_main(args)

    def test_user_backbone_empty_leaf_labels_raise_clear_error(self, tmp_path):
        species_path = tmp_path / 'species.txt'
        species_path.write_text('A\n')
        tree_path = tmp_path / 'tree.nwk'
        tree_path.write_text('(A:1,:1,B:1);\n')
        args = make_args(
            infile=str(tree_path),
            outfile='-',
            species_list=str(species_path),
            taxid_tsv=None,
            backbone='user',
            rank='no',
            collapse=False,
        )
        with pytest.raises(ValueError, match='Empty leaf labels'):
            constrain_main(args)
