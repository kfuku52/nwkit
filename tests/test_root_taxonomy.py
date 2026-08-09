import pandas as pd
import pytest
from ete4 import Tree

import nwkit.root as root_mod
from nwkit.root import DEFAULT_TAXONOMY_SOURCE_CHAIN, taxonomy_rooting
from tests.helpers import make_args
from tests.root_test_support import (
    install_fake_ncbi,
    install_fake_opentree,
    install_fake_timetree,
)


class TestTaxonomyHttpBounds:
    def test_streaming_reader_stops_at_decoded_body_limit(self, monkeypatch):
        monkeypatch.setattr(root_mod, 'TAXONOMY_HTTP_MAX_BYTES', 8)
        monkeypatch.setattr(root_mod, 'TAXONOMY_HTTP_CHUNK_BYTES', 4)

        class StreamingResponse:
            headers = {}
            encoding = 'utf-8'

            def __init__(self):
                self.num_chunks_yielded = 0
                self.closed = False

            @property
            def text(self):
                raise AssertionError('Streaming responses must not access .text')

            def iter_content(self, chunk_size):
                assert chunk_size == 4
                for chunk in (b'1234', b'5678', b'9abc', b'defg'):
                    self.num_chunks_yielded += 1
                    yield chunk

            def close(self):
                self.closed = True

        response = StreamingResponse()
        with pytest.raises(ValueError, match='exceeds the size limit'):
            root_mod._bounded_response_text(response, 'test source')

        assert response.num_chunks_yielded == 3
        assert response.closed is True


class TestTaxonomyRooting:
    def test_malformed_opentree_json_falls_back_to_next_source(self, monkeypatch):
        install_fake_opentree(
            monkeypatch,
            tnrs_json=[],
            induced_subtree_json={'newick': ''},
        )
        reference = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        monkeypatch.setattr(
            root_mod,
            '_build_ncbi_reference_tree',
            lambda **kwargs: (reference, {'A', 'B', 'C', 'D'}, set()),
        )
        tree = Tree('((A:1,B:1):1,C:1,D:1);', parser=1)
        rooted = taxonomy_rooting(tree, taxonomy_source='opentree,ncbi', rank='no')
        assert {'A', 'B'} in [set(child.leaf_names()) for child in rooted.get_children()]

    def test_default_source_chain_constant(self):
        assert DEFAULT_TAXONOMY_SOURCE_CHAIN == 'ncbi,opentree,timetree'

    def test_ncbi_source_passes_args_to_ncbi_helpers(self, monkeypatch, tmp_path):
        observed = dict()

        def fake_resolve_ncbi_lineages(tree, taxid_tsv=None, rank='no', args=None, verbose=False):
            observed['lineages_download_dir'] = getattr(args, 'download_dir', None)
            return ({
                'A': [1, 10],
                'B': [1, 10],
                'C': [1, 20],
                'D': [1, 20],
            }, {})

        def fake_taxid2tree(lineages, taxid_counts, args=None):
            observed['tree_download_dir'] = getattr(args, 'download_dir', None)
            return Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)

        monkeypatch.setattr(root_mod, '_resolve_ncbi_lineages', fake_resolve_ncbi_lineages)
        monkeypatch.setattr(root_mod, 'taxid2tree', fake_taxid2tree)

        tree = Tree('(A:1,B:1,(C:1,D:1):1);', parser=1)
        args = make_args(download_dir=str(tmp_path / 'cache'))
        rooted = taxonomy_rooting(tree, taxonomy_source='ncbi', rank='no', args=args)
        assert observed['lineages_download_dir'] == str(tmp_path / 'cache')
        assert observed['tree_download_dir'] == str(tmp_path / 'cache')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'A', 'B'} in child_leaf_sets
        assert {'C', 'D'} in child_leaf_sets

    def test_infers_ncbi_lineages_from_leaf_labels(self, monkeypatch):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={
                'Homo sapiens': 9606,
                'Pan troglodytes': 9598,
                'Arabidopsis thaliana': 3702,
                'Oryza sativa': 4530,
            },
            lineage_by_taxid={
                1: [1],
                10: [1, 10],
                20: [1, 20],
                9606: [1, 10, 9606],
                9598: [1, 10, 9598],
                3702: [1, 20, 3702],
                4530: [1, 20, 4530],
            },
        )
        original = Tree('(Homo_sapiens:1,Pan_troglodytes:1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);', parser=1)
        tree = Tree('(Homo_sapiens:1,Pan_troglodytes:1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);', parser=1)
        rooted = taxonomy_rooting(tree, taxonomy_source='ncbi', taxid_tsv=None, rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'Homo_sapiens', 'Pan_troglodytes'} in child_leaf_sets
        assert {'Arabidopsis_thaliana', 'Oryza_sativa'} in child_leaf_sets
        for l1 in rooted.leaves():
            for l2 in rooted.leaves():
                if l1.name != l2.name:
                    assert abs(original.get_distance(l1.name, l2.name) - rooted.get_distance(l1, l2)) < 1e-6

    def test_ncbi_taxonomic_uses_taxonomy_query_fallbacks(self, monkeypatch):
        observed = dict(queries=[])

        class FakeNCBI:
            def __init__(self, *args, **kwargs):
                self.db = None

            def get_name_translator(self, names):
                observed['queries'].extend(list(names))
                mapping = dict()
                for name in names:
                    if name == 'Dictyostelium discoideum':
                        mapping[name] = [101]
                    elif name == 'Amoeba':
                        mapping[name] = [102]
                    elif name == 'Homo sapiens':
                        mapping[name] = [201]
                    elif name == 'Pan troglodytes':
                        mapping[name] = [202]
                return mapping

            def get_lineage(self, taxid):
                lineage_by_taxid = {
                    1: [1],
                    10: [1, 10],
                    20: [1, 20],
                    101: [1, 10, 101],
                    102: [1, 10, 102],
                    201: [1, 20, 201],
                    202: [1, 20, 202],
                }
                return lineage_by_taxid[int(taxid)]

            def get_taxid_translator(self, taxids):
                names = {
                    1: 'root',
                    10: 'Amoebozoa',
                    20: 'Primates',
                    101: 'Dictyostelium discoideum',
                    102: 'Amoeba',
                    201: 'Homo sapiens',
                    202: 'Pan troglodytes',
                }
                return {int(taxid): names[int(taxid)] for taxid in taxids if int(taxid) in names}

            def get_rank(self, taxids):
                ranks = {
                    1: 'no rank',
                    10: 'phylum',
                    20: 'order',
                    101: 'species',
                    102: 'genus',
                    201: 'species',
                    202: 'species',
                }
                return {int(taxid): ranks.get(int(taxid), 'no rank') for taxid in taxids}

        monkeypatch.setattr(root_mod, 'get_ete_ncbitaxa', lambda args=None: FakeNCBI())

        tree = Tree(
            '(Dictyostelium_cf_discoideum:1,Amoeba_sp_JDSRuffled:1,(Homo_sapiens:1,Pan_troglodytes:1):1);',
            parser=1,
        )
        args = make_args(species_parser='taxonomic')
        rooted = taxonomy_rooting(tree, taxonomy_source='ncbi', taxid_tsv=None, rank='no', args=args)
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'Dictyostelium_cf_discoideum', 'Amoeba_sp_JDSRuffled'} in child_leaf_sets
        assert {'Homo_sapiens', 'Pan_troglodytes'} in child_leaf_sets
        assert 'Dictyostelium discoideum' in observed['queries']
        assert 'Amoeba' in observed['queries']
        assert 'Dictyostelium cf discoideum' not in observed['queries']
        assert 'Amoeba sp JDSRuffled' not in observed['queries']

    def test_taxid_tsv_requires_exact_leaf_label_match(self, monkeypatch, tmp_path):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={},
            lineage_by_taxid={
                1: [1],
                10: [1, 10],
                20: [1, 20],
                100: [1, 10, 100],
                101: [1, 10, 101],
                200: [1, 20, 200],
                201: [1, 20, 201],
            },
        )
        tsv_path = tmp_path / 'taxid.tsv'
        pd.DataFrame(
            {'leaf_name': ['A', 'B', 'C', 'X'], 'taxid': [100, 101, 200, 201]}
        ).to_csv(tsv_path, sep='\t', index=False)
        tree = Tree('(A:1,B:1,(C:1,D:1):1);', parser=1)
        with pytest.raises(ValueError, match='match the leaf labels'):
            taxonomy_rooting(tree, taxonomy_source='ncbi', taxid_tsv=str(tsv_path), rank='no')

    def test_taxid_tsv_numeric_leaf_label_mismatch_raises_valueerror(self, monkeypatch, tmp_path):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={},
            lineage_by_taxid={
                1: [1],
                10: [1, 10],
                20: [1, 20],
                100: [1, 10, 100],
                200: [1, 20, 200],
            },
        )
        tsv_path = tmp_path / 'taxid.tsv'
        pd.DataFrame(
            {'leaf_name': [1, 3], 'taxid': [100, 200]}
        ).to_csv(tsv_path, sep='\t', index=False)
        tree = Tree('(1:1,2:1);', parser=1)
        with pytest.raises(ValueError, match='match the leaf labels'):
            taxonomy_rooting(tree, taxonomy_source='ncbi', taxid_tsv=str(tsv_path), rank='no')

    def test_taxid_tsv_numeric_leaf_labels_are_accepted(self, monkeypatch, tmp_path):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={},
            lineage_by_taxid={
                1: [1],
                10: [1, 10],
                20: [1, 20],
                100: [1, 10, 100],
                101: [1, 10, 101],
                200: [1, 20, 200],
            },
        )
        tsv_path = tmp_path / 'taxid.tsv'
        pd.DataFrame(
            {'leaf_name': [1, 2, 3], 'taxid': [100, 101, 200]}
        ).to_csv(tsv_path, sep='\t', index=False)
        tree = Tree('(1:1,2:1,3:1);', parser=1)
        rooted = taxonomy_rooting(tree, taxonomy_source='ncbi', taxid_tsv=str(tsv_path), rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'1', '2'} in child_leaf_sets
        assert {'3'} in child_leaf_sets

    def test_taxid_tsv_leading_zero_leaf_labels_are_accepted(self, monkeypatch, tmp_path):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={},
            lineage_by_taxid={
                1: [1],
                10: [1, 10],
                20: [1, 20],
                100: [1, 10, 100],
                101: [1, 10, 101],
                200: [1, 20, 200],
            },
        )
        tsv_path = tmp_path / 'taxid.tsv'
        pd.DataFrame(
            {'leaf_name': ['001', '002', '003'], 'taxid': [100, 101, 200]}
        ).to_csv(tsv_path, sep='\t', index=False)
        tree = Tree('(001:1,002:1,003:1);', parser=1)
        rooted = taxonomy_rooting(tree, taxonomy_source='ncbi', taxid_tsv=str(tsv_path), rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'001', '002'} in child_leaf_sets
        assert {'003'} in child_leaf_sets

    def test_taxid_tsv_na_literal_leaf_labels_are_accepted(self, monkeypatch, tmp_path):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={},
            lineage_by_taxid={
                1: [1],
                10: [1, 10],
                20: [1, 20],
                100: [1, 10, 100],
                101: [1, 10, 101],
                200: [1, 20, 200],
            },
        )
        tsv_path = tmp_path / 'taxid.tsv'
        pd.DataFrame(
            {'leaf_name': ['NA', 'B', 'C'], 'taxid': [100, 101, 200]}
        ).to_csv(tsv_path, sep='\t', index=False)
        tree = Tree('(NA:1,B:1,C:1);', parser=1)
        rooted = taxonomy_rooting(tree, taxonomy_source='ncbi', taxid_tsv=str(tsv_path), rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'NA', 'B'} in child_leaf_sets
        assert {'C'} in child_leaf_sets

    def test_ambiguous_taxonomy_root_raises(self, monkeypatch, tmp_path):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={},
            lineage_by_taxid={
                1: [1],
                100: [1, 100],
                200: [1, 200],
                300: [1, 300],
            },
        )
        tsv_path = tmp_path / 'taxid.tsv'
        pd.DataFrame(
            {'leaf_name': ['A', 'B', 'C'], 'taxid': [100, 200, 300]}
        ).to_csv(tsv_path, sep='\t', index=False)
        tree = Tree('(A:1,B:1,C:1);', parser=1)
        with pytest.raises(ValueError, match='ambiguous'):
            taxonomy_rooting(tree, taxonomy_source='ncbi', taxid_tsv=str(tsv_path), rank='no')

    def test_incompatible_taxonomy_bipartition_raises(self, monkeypatch):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={
                'Homo sapiens': 9606,
                'Pan troglodytes': 9598,
                'Arabidopsis thaliana': 3702,
                'Oryza sativa': 4530,
            },
            lineage_by_taxid={
                1: [1],
                10: [1, 10],
                20: [1, 20],
                9606: [1, 10, 9606],
                9598: [1, 10, 9598],
                3702: [1, 20, 3702],
                4530: [1, 20, 4530],
            },
        )
        tree = Tree('((Homo_sapiens:1,Arabidopsis_thaliana:1):1,(Pan_troglodytes:1,Oryza_sativa:1):1);', parser=1)
        with pytest.raises(ValueError, match='ncbi-derived root bipartition'):
            taxonomy_rooting(tree, taxonomy_source='ncbi', taxid_tsv=None, rank='no')

    def test_ncbi_placeholder_taxa_are_excluded_when_not_interfering(self, monkeypatch):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={
                'Homo sapiens': 9606,
                'Pan troglodytes': 9598,
                'Arabidopsis thaliana': 3702,
                'Oryza sativa': 4530,
                'Unknown': 32644,
            },
            lineage_by_taxid={
                1: [1],
                10: [1, 10],
                20: [1, 20],
                2787823: [1, 2787823],
                12908: [1, 2787823, 12908],
                9606: [1, 10, 9606],
                9598: [1, 10, 9598],
                3702: [1, 20, 3702],
                4530: [1, 20, 4530],
                32644: [1, 2787823, 12908, 32644],
            },
            taxid_to_name={
                1: 'root',
                10: 'Primates',
                20: 'Plants',
                2787823: 'unclassified entries',
                12908: 'unclassified sequences',
                9606: 'Homo sapiens',
                9598: 'Pan troglodytes',
                3702: 'Arabidopsis thaliana',
                4530: 'Oryza sativa',
                32644: 'unidentified',
            },
            rank_by_taxid={32644: 'species'},
        )
        tree = Tree(
            '((Homo_sapiens:1,(Pan_troglodytes:1,(Unknown_amoeba:1,Unknown_amoeboid:1):1):1):1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);',
            parser=1,
        )
        rooted = taxonomy_rooting(tree, taxonomy_source='ncbi', taxid_tsv=None, rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'Arabidopsis_thaliana', 'Oryza_sativa'} in child_leaf_sets
        assert {'Homo_sapiens', 'Pan_troglodytes', 'Unknown_amoeba', 'Unknown_amoeboid'} in child_leaf_sets

    def test_ncbi_placeholder_taxa_on_root_path_raise(self, monkeypatch):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={
                'Homo sapiens': 9606,
                'Pan troglodytes': 9598,
                'Arabidopsis thaliana': 3702,
                'Oryza sativa': 4530,
                'Unknown': 32644,
            },
            lineage_by_taxid={
                1: [1],
                10: [1, 10],
                20: [1, 20],
                2787823: [1, 2787823],
                12908: [1, 2787823, 12908],
                9606: [1, 10, 9606],
                9598: [1, 10, 9598],
                3702: [1, 20, 3702],
                4530: [1, 20, 4530],
                32644: [1, 2787823, 12908, 32644],
            },
            taxid_to_name={
                1: 'root',
                10: 'Primates',
                20: 'Plants',
                2787823: 'unclassified entries',
                12908: 'unclassified sequences',
                9606: 'Homo sapiens',
                9598: 'Pan troglodytes',
                3702: 'Arabidopsis thaliana',
                4530: 'Oryza sativa',
                32644: 'unidentified',
            },
            rank_by_taxid={32644: 'species'},
        )
        tree = Tree(
            '(((Homo_sapiens:1,Pan_troglodytes:1):1,(Unknown_amoeba:1,Unknown_amoeboid:1):1):1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);',
            parser=1,
        )
        with pytest.raises(ValueError, match='interfere'):
            taxonomy_rooting(tree, taxonomy_source='ncbi', taxid_tsv=None, rank='no')

    def test_timetree_maps_species_names_back_to_original_labels(self, monkeypatch):
        calls = install_fake_timetree(
            monkeypatch,
            upload_html='<div id="prunetree-msg-box"></div>',
            newick_text='((Homo_sapiens:1,Pan_troglodytes:1):10,(Arabidopsis_thaliana:1,Oryza_sativa:1):20);',
        )
        tree = Tree(
            '(Homo_sapiens_gene1:1,Pan_troglodytes_gene1:1,(Arabidopsis_thaliana_gene1:1,Oryza_sativa_gene1:1):1);',
            parser=1,
        )
        rooted = taxonomy_rooting(tree, taxonomy_source='timetree', rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'Homo_sapiens_gene1', 'Pan_troglodytes_gene1'} in child_leaf_sets
        assert {'Arabidopsis_thaliana_gene1', 'Oryza_sativa_gene1'} in child_leaf_sets
        assert len(calls) == 2
        assert all(stream is True for _, stream in calls)

    def test_timetree_allows_monophyletic_duplicate_species_labels(self, monkeypatch):
        install_fake_timetree(
            monkeypatch,
            upload_html='<div id="prunetree-msg-box"></div>',
            newick_text='((Homo_sapiens:1,Pan_troglodytes:1):10,((Arabidopsis_thaliana:1,Oryza_sativa:1):20,Saccharomyces_cerevisiae:1):30);',
        )
        tree = Tree(
            '((Arabidopsis_thaliana_gene1:1,Oryza_sativa_gene1:1):1,(Saccharomyces_cerevisiae_gene1:1,((Homo_sapiens_gene1:1,Homo_sapiens_gene2:1):1,Pan_troglodytes_gene1:1):1):1);',
            parser=1,
        )
        rooted = taxonomy_rooting(tree, taxonomy_source='timetree', rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'Homo_sapiens_gene1', 'Homo_sapiens_gene2', 'Pan_troglodytes_gene1'} in child_leaf_sets
        assert {'Arabidopsis_thaliana_gene1', 'Oryza_sativa_gene1', 'Saccharomyces_cerevisiae_gene1'} in child_leaf_sets

    def test_timetree_unresolved_names_are_excluded_when_not_interfering(self, monkeypatch):
        install_fake_timetree(
            monkeypatch,
            upload_html=(
                '<button id="prunetree-msg-btn" class="error-btn-enabled">Unresolved Names (1)</button>'
                '<div id="unresolved-names">Unknown amoeba (not found in NCBI taxonomy)<br/></div>'
            ),
            newick_text='((Homo_sapiens:1,Pan_troglodytes:1):10,(Arabidopsis_thaliana:1,Oryza_sativa:1):20);',
        )
        tree = Tree(
            '((Homo_sapiens:1,(Pan_troglodytes:1,Unknown_amoeba:1):1):1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);',
            parser=1,
        )
        rooted = taxonomy_rooting(tree, taxonomy_source='timetree', rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'Homo_sapiens', 'Pan_troglodytes', 'Unknown_amoeba'} in child_leaf_sets
        assert {'Arabidopsis_thaliana', 'Oryza_sativa'} in child_leaf_sets

    def test_source_chain_falls_back_to_timetree(self, monkeypatch):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={},
            lineage_by_taxid={},
        )
        install_fake_timetree(
            monkeypatch,
            upload_html='<div id="prunetree-msg-box"></div>',
            newick_text='((Homo_sapiens:1,Pan_troglodytes:1):10,(Arabidopsis_thaliana:1,Oryza_sativa:1):20);',
        )
        tree = Tree('(Homo_sapiens:1,Pan_troglodytes:1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);', parser=1)
        rooted = taxonomy_rooting(tree, taxonomy_source='ncbi,timetree', taxid_tsv=None, rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'Homo_sapiens', 'Pan_troglodytes'} in child_leaf_sets
        assert {'Arabidopsis_thaliana', 'Oryza_sativa'} in child_leaf_sets

    def test_opentree_maps_species_names_back_to_original_labels(self, monkeypatch):
        calls = install_fake_opentree(
            monkeypatch,
            tnrs_json={
                'results': [
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 1, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 2, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 3, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 4, 'is_suppressed_from_synth': False}}]},
                ],
            },
            induced_subtree_json={
                'broken': {},
                'newick': '((((Homo_sapiens,Pan_troglodytes)Homininae)Primates,(((Arabidopsis_thaliana,Oryza_sativa)Mesangiospermae)Embryophyta)Viridiplantae)Eukaryota);',
            },
        )
        tree = Tree(
            '(Homo_sapiens_gene1:1,Pan_troglodytes_gene1:1,(Arabidopsis_thaliana_gene1:1,Oryza_sativa_gene1:1):1);',
            parser=1,
        )
        rooted = taxonomy_rooting(tree, taxonomy_source='opentree', rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'Homo_sapiens_gene1', 'Pan_troglodytes_gene1'} in child_leaf_sets
        assert {'Arabidopsis_thaliana_gene1', 'Oryza_sativa_gene1'} in child_leaf_sets
        assert len(calls) == 2
        assert all(stream is True for _, stream in calls)

    def test_opentree_allows_monophyletic_duplicate_species_labels(self, monkeypatch):
        install_fake_opentree(
            monkeypatch,
            tnrs_json={
                'results': [
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 1, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 2, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 3, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 4, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 5, 'is_suppressed_from_synth': False}}]},
                ],
            },
            induced_subtree_json={
                'broken': {},
                'newick': '((Homo_sapiens,Pan_troglodytes)Primates,((Arabidopsis_thaliana,Oryza_sativa)Mesangiospermae,Saccharomyces_cerevisiae)Opisthokonta)Eukaryota;',
            },
        )
        tree = Tree(
            '((Arabidopsis_thaliana_gene1:1,Oryza_sativa_gene1:1):1,(Saccharomyces_cerevisiae_gene1:1,((Homo_sapiens_gene1:1,Homo_sapiens_gene2:1):1,Pan_troglodytes_gene1:1):1):1);',
            parser=1,
        )
        rooted = taxonomy_rooting(tree, taxonomy_source='opentree', rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'Homo_sapiens_gene1', 'Homo_sapiens_gene2', 'Pan_troglodytes_gene1'} in child_leaf_sets
        assert {'Arabidopsis_thaliana_gene1', 'Oryza_sativa_gene1', 'Saccharomyces_cerevisiae_gene1'} in child_leaf_sets

    def test_opentree_unresolved_name_is_excluded_when_not_interfering(self, monkeypatch):
        install_fake_opentree(
            monkeypatch,
            tnrs_json={
                'results': [
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 1, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 2, 'is_suppressed_from_synth': False}}]},
                    {'matches': []},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 3, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 4, 'is_suppressed_from_synth': False}}]},
                ],
            },
            induced_subtree_json={
                'broken': {},
                'newick': '((Homo_sapiens,Pan_troglodytes)Primates,(Arabidopsis_thaliana,Oryza_sativa)Mesangiospermae)Eukaryota;',
            },
        )
        tree = Tree(
            '((Homo_sapiens:1,(Pan_troglodytes:1,Batrachochytrium_dendrobatidis:1):1):1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);',
            parser=1,
        )
        rooted = taxonomy_rooting(tree, taxonomy_source='opentree', rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'Homo_sapiens', 'Pan_troglodytes', 'Batrachochytrium_dendrobatidis'} in child_leaf_sets
        assert {'Arabidopsis_thaliana', 'Oryza_sativa'} in child_leaf_sets

    def test_opentree_unresolved_name_on_root_path_raises(self, monkeypatch):
        install_fake_opentree(
            monkeypatch,
            tnrs_json={
                'results': [
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 1, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 2, 'is_suppressed_from_synth': False}}]},
                    {'matches': []},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 3, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 4, 'is_suppressed_from_synth': False}}]},
                ],
            },
            induced_subtree_json={
                'broken': {},
                'newick': '((Homo_sapiens,Pan_troglodytes)Primates,(Arabidopsis_thaliana,Oryza_sativa)Mesangiospermae)Eukaryota;',
            },
        )
        tree = Tree(
            '(((Homo_sapiens:1,Pan_troglodytes:1):1,Batrachochytrium_dendrobatidis:1):1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);',
            parser=1,
        )
        with pytest.raises(ValueError, match='interfere'):
            taxonomy_rooting(tree, taxonomy_source='opentree', rank='no')

    def test_source_chain_falls_back_to_opentree(self, monkeypatch):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={},
            lineage_by_taxid={},
        )
        install_fake_opentree(
            monkeypatch,
            tnrs_json={
                'results': [
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 1, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 2, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 3, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 4, 'is_suppressed_from_synth': False}}]},
                ],
            },
            induced_subtree_json={
                'broken': {},
                'newick': '((Homo_sapiens,Pan_troglodytes)Primates,(Arabidopsis_thaliana,Oryza_sativa)Mesangiospermae)Eukaryota;',
            },
        )
        tree = Tree('(Homo_sapiens:1,Pan_troglodytes:1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);', parser=1)
        rooted = taxonomy_rooting(tree, taxonomy_source='ncbi,opentree', taxid_tsv=None, rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'Homo_sapiens', 'Pan_troglodytes'} in child_leaf_sets
        assert {'Arabidopsis_thaliana', 'Oryza_sativa'} in child_leaf_sets

    def test_default_source_chain_falls_back_to_opentree(self, monkeypatch):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={},
            lineage_by_taxid={},
        )
        install_fake_opentree(
            monkeypatch,
            tnrs_json={
                'results': [
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 1, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 2, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 3, 'is_suppressed_from_synth': False}}]},
                    {'matches': [{'is_approximate_match': False, 'taxon': {'ott_id': 4, 'is_suppressed_from_synth': False}}]},
                ],
            },
            induced_subtree_json={
                'broken': {},
                'newick': '((Homo_sapiens,Pan_troglodytes)Primates,(Arabidopsis_thaliana,Oryza_sativa)Mesangiospermae)Eukaryota;',
            },
        )
        tree = Tree('(Homo_sapiens:1,Pan_troglodytes:1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);', parser=1)
        rooted = taxonomy_rooting(tree, rank='no')
        child_leaf_sets = [set(child.leaf_names()) for child in rooted.get_children()]
        assert {'Homo_sapiens', 'Pan_troglodytes'} in child_leaf_sets
        assert {'Arabidopsis_thaliana', 'Oryza_sativa'} in child_leaf_sets

    def test_unknown_taxonomy_source_raises(self):
        tree = Tree('(A:1,B:1,(C:1,D:1):1);', parser=1)
        with pytest.raises(ValueError, match='Unknown taxonomy source'):
            taxonomy_rooting(tree, taxonomy_source='ncbi,unknown', rank='no')
