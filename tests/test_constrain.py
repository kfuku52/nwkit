import pandas as pd
import pytest
import numpy as np
from argparse import Namespace
from ete4 import Tree

from nwkit.constrain import (
    get_max_ancestor_overlap_node,
    check_input_file,
    constrain_main,
    collapse_genes,
    get_lineages,
    get_lineages_from_taxid,
    get_mrca_taxid,
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
        pd.DataFrame(
            {'leaf_name': ['A', 'A'], 'taxid': [9606, 9606]}
        ).to_csv(tsv_path, sep='\t', index=False)
        args = Namespace(species_list=None, taxid_tsv=str(tsv_path), backbone='ncbi')
        with pytest.raises(ValueError, match='Duplicate values'):
            check_input_file(args)

    def test_taxid_tsv_missing_taxid_raises(self, tmp_path):
        tsv_path = tmp_path / 'taxid.tsv'
        pd.DataFrame(
            {'leaf_name': ['A', 'B'], 'taxid': [9606, None]}
        ).to_csv(tsv_path, sep='\t', index=False)
        args = Namespace(species_list=None, taxid_tsv=str(tsv_path), backbone='ncbi')
        with pytest.raises(ValueError, match='missing values'):
            check_input_file(args)

    def test_taxid_tsv_non_numeric_taxid_raises(self, tmp_path):
        tsv_path = tmp_path / 'taxid.tsv'
        pd.DataFrame(
            {'leaf_name': ['A', 'B'], 'taxid': ['9606', 'abc']}
        ).to_csv(tsv_path, sep='\t', index=False)
        args = Namespace(species_list=None, taxid_tsv=str(tsv_path), backbone='ncbi')
        with pytest.raises(ValueError, match='non-numeric'):
            check_input_file(args)

    def test_taxid_tsv_non_integer_taxid_raises(self, tmp_path):
        tsv_path = tmp_path / 'taxid.tsv'
        pd.DataFrame(
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


class TestNcbiDownloadDirRouting:
    def test_get_lineages_uses_args_for_ncbi_db(self, monkeypatch, tmp_path):
        observed = dict()

        class FakeNCBI:
            def __init__(self):
                self.db = None

            def get_name_translator(self, names):
                return {'Homo sapiens': [9606]}

            def get_lineage(self, taxid):
                return [1, 9606]

        def fake_get_ete_ncbitaxa(args=None):
            observed['download_dir'] = getattr(args, 'download_dir', None)
            return FakeNCBI()

        monkeypatch.setattr('nwkit.constrain.get_ete_ncbitaxa', fake_get_ete_ncbitaxa)
        args = make_args(download_dir=str(tmp_path / 'cache'))
        lineages = get_lineages(['Homo_sapiens_gene1'], rank='no', args=args)
        assert observed['download_dir'] == str(tmp_path / 'cache')
        assert lineages == {'Homo_sapiens_gene1': [1, 9606]}

    def test_get_lineages_from_taxid_uses_args_for_ncbi_db(self, monkeypatch, tmp_path):
        observed = dict()

        class FakeNCBI:
            def __init__(self):
                self.db = None

            def get_lineage(self, taxid):
                return [1, int(taxid)]

            def get_rank(self, lineage):
                return {taxid: 'no_rank' for taxid in lineage}

        def fake_get_ete_ncbitaxa(args=None):
            observed['download_dir'] = getattr(args, 'download_dir', None)
            return FakeNCBI()

        monkeypatch.setattr('nwkit.constrain.get_ete_ncbitaxa', fake_get_ete_ncbitaxa)
        taxid_df = pd.DataFrame({'leaf_name': ['A'], 'taxid': [9606]})
        args = make_args(download_dir=str(tmp_path / 'cache'))
        lineages = get_lineages_from_taxid(taxid_df, rank='no', args=args)
        assert observed['download_dir'] == str(tmp_path / 'cache')
        assert lineages == {'A': [1, 9606]}

    def test_match_taxa_ncbi_uses_args_for_ncbi_db(self, monkeypatch, tmp_path):
        observed = dict()

        class FakeNCBI:
            def __init__(self):
                self.db = None

            def get_name_translator(self, names):
                return {'Homo sapiens': [9606]}

            def get_lineage(self, taxid):
                return [1, 10, 9606]

            def get_taxid_translator(self, taxids):
                return {1: 'root', 10: 'Primates', 9606: 'Homo sapiens'}

        def fake_get_ete_ncbitaxa(args=None):
            observed['download_dir'] = getattr(args, 'download_dir', None)
            return FakeNCBI()

        monkeypatch.setattr('nwkit.constrain.get_ete_ncbitaxa', fake_get_ete_ncbitaxa)
        tree = initialize_tree(Tree('(Primates:1,Plants:1);', parser=1))
        args = make_args(download_dir=str(tmp_path / 'cache'))
        out = match_taxa(tree=tree, labels=['Homo_sapiens'], backbone_method='ncbi_user', args=args)
        assert observed['download_dir'] == str(tmp_path / 'cache')
        primates = [leaf for leaf in out.leaves() if leaf.name == 'Primates'][0]
        assert primates.props.get('has_taxon') is True
        assert primates.props.get('taxon_names') == ['Homo_sapiens']

    def test_taxid2tree_uses_args_for_ncbi_db(self, monkeypatch, tmp_path):
        observed = dict()

        class FakeNCBI:
            def __init__(self):
                self.db = None

            def get_lineage(self, taxid):
                return [1, int(taxid)]

        def fake_get_ete_ncbitaxa(args=None):
            observed['download_dir'] = getattr(args, 'download_dir', None)
            return FakeNCBI()

        monkeypatch.setattr('nwkit.constrain.get_ete_ncbitaxa', fake_get_ete_ncbitaxa)
        args = make_args(download_dir=str(tmp_path / 'cache'))
        lineages = {'A': [1, 10], 'B': [1, 10], 'C': [1, 20], 'D': [1, 20]}
        tree = taxid2tree(lineages, get_taxid_counts(lineages), args=args)
        assert observed['download_dir'] == str(tmp_path / 'cache')
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_get_mrca_taxid_uses_args_for_ncbi_db(self, monkeypatch, tmp_path):
        observed = dict()

        class FakeNCBI:
            def __init__(self):
                self.db = None

            def get_lineage(self, taxid):
                taxid = int(taxid)
                if taxid == 10:
                    return [1, 10]
                return [1, 10, taxid]

        def fake_get_ete_ncbitaxa(args=None):
            observed['download_dir'] = getattr(args, 'download_dir', None)
            return FakeNCBI()

        monkeypatch.setattr('nwkit.constrain.get_ete_ncbitaxa', fake_get_ete_ncbitaxa)
        args = make_args(download_dir=str(tmp_path / 'cache'))
        multi_counts = np.array([[10, 2], [9606, 2]])
        assert get_mrca_taxid(multi_counts, args=args) == 9606
        assert observed['download_dir'] == str(tmp_path / 'cache')


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
