import pytest
from ete4 import Tree

from nwkit.util import (
    annotate_duplication_confidence_scores,
    annotate_scientific_names,
    extract_species_label,
    extract_taxonomy_query,
    get_monophyletic_species_groups,
    get_subtree_leaf_bitmasks,
    get_subtree_leaf_name_sets,
    get_target_nodes,
    is_all_leaf_names_identical,
    is_rooted,
    label2sciname,
    read_item_per_line_file,
    remove_singleton,
    validate_unique_named_leaves,
)
from tests.helpers import make_args


class TestSpeciesGrouping:
    def test_monophyletic_duplicate_species_are_grouped(self):
        tree = Tree('(((Homo_sapiens_gene1:1,Homo_sapiens_gene2:1):1,Pan_troglodytes_gene1:1):1,Mus_musculus_gene1:1);', parser=1)
        leaf_name_to_sci_name, species_to_leaf_names = get_monophyletic_species_groups(tree)
        assert leaf_name_to_sci_name['Homo_sapiens_gene1'] == 'Homo_sapiens'
        assert leaf_name_to_sci_name['Homo_sapiens_gene2'] == 'Homo_sapiens'
        assert species_to_leaf_names['Homo_sapiens'] == ['Homo_sapiens_gene1', 'Homo_sapiens_gene2']

    def test_split_duplicate_species_raise(self):
        tree = Tree('((Homo_sapiens_gene1:1,Pan_troglodytes_gene1:1):1,(Homo_sapiens_gene2:1,Mus_musculus_gene1:1):1);', parser=1)
        with pytest.raises(ValueError, match='not monophyletic'):
            get_monophyletic_species_groups(tree)


class TestRemoveSingleton:
    def test_remove_singleton_node(self):
        # Create tree with a singleton: ((A,B)) -> should become (A,B)
        tree = Tree('(((A:1,B:1):1):1,(C:1,D:1):1);', parser=1)
        num_nodes_before = len(list(tree.traverse()))
        tree = remove_singleton(tree, verbose=False, preserve_branch_length=True)
        num_nodes_after = len(list(tree.traverse()))
        assert num_nodes_after < num_nodes_before
        assert set(tree.leaf_names()) == {'A', 'B', 'C', 'D'}

    def test_no_singleton(self, simple_tree):
        num_nodes_before = len(list(simple_tree.traverse()))
        tree = remove_singleton(simple_tree, verbose=False)
        num_nodes_after = len(list(tree.traverse()))
        assert num_nodes_before == num_nodes_after

    def test_preserve_branch_length(self):
        tree = Tree('(((A:1,B:1):2):3,C:6);', parser=1)
        tree = remove_singleton(tree, verbose=False, preserve_branch_length=True)
        # After removing singleton, branch lengths should be preserved (summed)
        assert set(tree.leaf_names()) == {'A', 'B', 'C'}

    def test_remove_singleton_root_wrapper(self):
        tree = Tree('((A:1,B:1):1);', parser=1)
        tree = remove_singleton(tree, verbose=False, preserve_branch_length=True)
        assert set(tree.leaf_names()) == {'A', 'B'}
        assert len(tree.get_children()) == 2
        assert abs(tree.get_distance('A', 'B') - 2.0) < 1e-6
        for node in tree.traverse():
            if not node.is_leaf:
                assert len(node.get_children()) != 1

    def test_remove_singleton_root_combines_stem_lengths(self):
        tree = Tree('((A:1,B:1)INNER:3)ROOT:4;', parser=1)
        tree = remove_singleton(tree, verbose=False, preserve_branch_length=True)
        assert tree.name == 'INNER'
        assert tree.dist == pytest.approx(7.0)

    @pytest.mark.parametrize(
        'newick',
        [
            '(((A:1,B:1):1e308):1e308,C:1);',
            '((A:1,B:1):1e308):1e308;',
        ],
    )
    def test_remove_singleton_rejects_overflow_before_mutating_tree(self, newick):
        tree = Tree(newick, parser=1)
        original_nodes = list(tree.traverse())
        original_state = [
            (
                node,
                node.up,
                tuple(node.get_children()),
                node.name,
                node.dist,
                node.support,
            )
            for node in original_nodes
        ]

        with pytest.raises(ValueError, match='finite branch length'):
            remove_singleton(
                tree,
                verbose=False,
                preserve_branch_length=True,
            )

        assert list(tree.traverse()) == original_nodes
        assert [
            (
                node,
                node.up,
                tuple(node.get_children()),
                node.name,
                node.dist,
                node.support,
            )
            for node in original_nodes
        ] == original_state


class TestLabel2Sciname:
    def test_single_string(self):
        result = label2sciname('Homo_sapiens_GENE1')
        assert result == 'Homo_sapiens'

    def test_list_input(self):
        result = label2sciname(['Homo_sapiens_GENE1', 'Mus_musculus_GENE2'])
        assert result == ['Homo_sapiens', 'Mus_musculus']

    def test_no_species_info(self):
        result = label2sciname('SingleWord')
        assert result is None

    def test_custom_delimiter(self):
        result = label2sciname('Homo-sapiens-GENE1', in_delim='-')
        assert result == 'Homo_sapiens'

    def test_custom_out_delimiter(self):
        result = label2sciname('Homo_sapiens_GENE1', out_delim=' ')
        assert result == 'Homo sapiens'

    def test_empty_list(self):
        result = label2sciname([])
        assert result == []

    def test_none_scalar(self):
        result = label2sciname(None)
        assert result is None

    def test_list_with_none(self):
        result = label2sciname(['Homo_sapiens_GENE1', None, 'Mus_musculus_GENE2'])
        assert result == ['Homo_sapiens', None, 'Mus_musculus']


class TestExtractSpeciesLabel:
    def test_default_regex_parses_with_suffix(self):
        assert extract_species_label('Homo_sapiens_GENE1') == 'Homo_sapiens'

    def test_default_regex_parses_exact_binomial(self):
        assert extract_species_label('Homo_sapiens') == 'Homo_sapiens'

    def test_custom_regex_uses_capture_groups(self):
        assert extract_species_label(
            'Homo.sapiens|GENE1',
            species_regex=r'^([A-Za-z]+)\.([A-Za-z]+)\|',
        ) == 'Homo_sapiens'

    def test_custom_regex_allows_space_output(self):
        assert extract_species_label(
            'Homo.sapiens|GENE1',
            species_regex=r'^([A-Za-z]+)\.([A-Za-z]+)\|',
            out_delim=' ',
        ) == 'Homo sapiens'

    def test_taxonomic_keeps_species_label_qualifier(self):
        args = make_args(species_parser='taxonomic')
        assert extract_species_label('Dictyostelium_cf_discoideum', args=args) == 'Dictyostelium_cf_discoideum'
        assert extract_species_label('Amoeba_sp_JDSRuffled', args=args) == 'Amoeba_sp_JDSRuffled'
        assert extract_species_label('Solanum_lycopersicum_cultivar_Heinz1706_gene1', args=args) == 'Solanum_lycopersicum_cultivar_Heinz1706'
        assert extract_species_label('Escherichia_coli_serovar_O157_gene1', args=args) == 'Escherichia_coli_serovar_O157'

    def test_taxonomic_taxonomy_query_falls_back(self):
        args = make_args(species_parser='taxonomic')
        assert extract_taxonomy_query('Dictyostelium_cf_discoideum', args=args) == 'Dictyostelium discoideum'
        assert extract_taxonomy_query('Amoeba_sp_JDSRuffled', args=args) == 'Amoeba'
        assert extract_taxonomy_query('Solanum_lycopersicum_cultivar_Heinz1706_gene1', args=args) == 'Solanum lycopersicum'
        assert extract_taxonomy_query('Escherichia_coli_serovar_O157_gene1', args=args) == 'Escherichia coli'

    def test_species_map_tsv_overrides_species_label_and_taxonomy_query(self, tmp_path):
        map_path = tmp_path / 'species_map.tsv'
        map_path.write_text(
            'leaf_name\tspecies_label\ttaxonomy_query\n'
            'Sample42\tMapped_species\tMapped species\n'
        )
        args = make_args(species_map_tsv=str(map_path))
        assert extract_species_label('Sample42', args=args) == 'Mapped_species'
        assert extract_taxonomy_query('Sample42', args=args) == 'Mapped species'


class TestReadItemPerLineFile:
    def test_basic(self, tmp_path):
        f = tmp_path / 'items.txt'
        f.write_text('apple\nbanana\ncherry\n')
        result = read_item_per_line_file(str(f))
        assert result == ['apple', 'banana', 'cherry']

    def test_empty_lines_stripped(self, tmp_path):
        f = tmp_path / 'items.txt'
        f.write_text('apple\n\nbanana\n\n')
        result = read_item_per_line_file(str(f))
        assert result == ['apple', 'banana']

    def test_crlf_and_whitespace_are_normalized(self, tmp_path):
        f = tmp_path / 'items.txt'
        f.write_bytes(b' apple \r\nbanana\r\n\r\n  cherry  \r\n')
        result = read_item_per_line_file(str(f))
        assert result == ['apple', 'banana', 'cherry']


class TestAnnotateScientificNames:
    def test_annotate(self, species_tree):
        tree = annotate_scientific_names(species_tree)
        sci_names = [leaf.props.get('sci_name') for leaf in tree.leaves()]
        assert 'Homo_sapiens' in sci_names
        assert 'Mus_musculus' in sci_names
        assert 'Danio_rerio' in sci_names

    def test_annotate_with_unnamed_leaf(self):
        tree = Tree('((:1,Homo_sapiens_gene1:1):1,Mus_musculus_gene1:1);', parser=1)
        tree = annotate_scientific_names(tree)
        unnamed = [leaf for leaf in tree.leaves() if not leaf.name][0]
        assert unnamed.props.get('sci_name') is None


class TestAnnotateDuplicationConfidenceScores:
    def test_annotate(self, species_tree):
        tree = annotate_scientific_names(species_tree)
        tree = annotate_duplication_confidence_scores(tree)
        for node in tree.traverse():
            if not node.is_leaf:
                assert 'dup_conf_score' in node.props
                assert 0 <= node.props.get('dup_conf_score') <= 1

    def test_no_duplication(self):
        # Tree where each species appears once in each child clade
        nwk = '((Homo_sapiens_G1:1,Mus_musculus_G1:1):1,(Danio_rerio_G1:1,Xenopus_laevis_G1:1):1);'
        tree = Tree(nwk, parser=1)
        tree = annotate_scientific_names(tree)
        tree = annotate_duplication_confidence_scores(tree)
        # Root node has no species overlap -> dup_conf_score = 0
        assert tree.props.get('dup_conf_score') == 0.0

    def test_full_duplication(self):
        # Tree where both children have the same species
        nwk = '((Homo_sapiens_G1:1,Mus_musculus_G1:1):1,(Homo_sapiens_G2:1,Mus_musculus_G2:1):1);'
        tree = Tree(nwk, parser=1)
        tree = annotate_scientific_names(tree)
        tree = annotate_duplication_confidence_scores(tree)
        assert tree.props.get('dup_conf_score') == 1.0

    def test_polytomy_internal_node_gets_zero_dup_conf_score(self):
        nwk = '(Homo_sapiens_G1:1,Homo_sapiens_G2:1,Mus_musculus_G1:1);'
        tree = Tree(nwk, parser=1)
        tree = annotate_scientific_names(tree)
        tree = annotate_duplication_confidence_scores(tree)
        assert tree.props.get('dup_conf_score') == 0.0


class TestGetSubtreeLeafNameSets:
    def test_sets_are_correct(self):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        sets = get_subtree_leaf_name_sets(tree)
        assert sets[tree] == {'A', 'B', 'C', 'D'}
        assert sets[tree.common_ancestor(['A', 'B'])] == {'A', 'B'}
        assert sets[tree.common_ancestor(['C', 'D'])] == {'C', 'D'}
        for leaf in tree.leaves():
            assert sets[leaf] == {leaf.name}


class TestGetSubtreeLeafBitmasks:
    def test_bitmasks_are_correct(self):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        leaf_name_to_bit = {'A': 0, 'B': 1, 'C': 2, 'D': 3}
        masks = get_subtree_leaf_bitmasks(tree, leaf_name_to_bit)
        assert masks[tree] == 0b1111
        assert masks[tree.common_ancestor(['A', 'B'])] == 0b0011
        assert masks[tree.common_ancestor(['C', 'D'])] == 0b1100
        for leaf in tree.leaves():
            assert masks[leaf] == (1 << leaf_name_to_bit[leaf.name])

    def test_missing_leaf_in_mapping_raises_clear_error(self):
        tree = Tree('((A:1,B:1):1,(C:1,D:1):1);', parser=1)
        leaf_name_to_bit = {'A': 0, 'B': 1, 'C': 2}
        with pytest.raises(ValueError, match='not found in reference mapping'):
            get_subtree_leaf_bitmasks(tree, leaf_name_to_bit)


class TestValidateUniqueNamedLeaves:
    def test_duplicate_leaf_names_raise(self):
        tree = Tree('((A:1,A:1):1,B:1);', parser=1)
        with pytest.raises(ValueError, match='Duplicated leaf labels'):
            validate_unique_named_leaves(tree, '--infile', " for 'transfer'")

    def test_empty_leaf_names_raise(self):
        tree = Tree('(A:1,:1,B:1);', parser=1)
        with pytest.raises(ValueError, match='Empty leaf labels'):
            validate_unique_named_leaves(tree, '--infile', " for 'transfer'")


class TestIsAllLeafNamesIdentical:
    def test_identical(self):
        t1 = Tree('((A,B),(C,D));', parser=1)
        t2 = Tree('((A,C),(B,D));', parser=1)
        assert is_all_leaf_names_identical(t1, t2) is True

    def test_not_identical(self):
        t1 = Tree('((A,B),(C,D));', parser=1)
        t2 = Tree('((A,B),(C,E));', parser=1)
        assert is_all_leaf_names_identical(t1, t2) is False

    def test_verbose_mode_handles_none_leaf_names(self):
        t1 = Tree('((A,B),(C,D));', parser=1)
        t2 = Tree('((A,B),(C,D));', parser=1)
        first_leaf = next(iter(t2.leaves()))
        first_leaf.name = None
        assert is_all_leaf_names_identical(t1, t2, verbose=True) is False


class TestGetTargetNodes:
    def test_all(self, simple_tree):
        nodes = get_target_nodes(simple_tree, 'all')
        assert len(nodes) == len(list(simple_tree.traverse()))

    def test_root(self, simple_tree):
        nodes = get_target_nodes(simple_tree, 'root')
        assert len(nodes) == 1
        assert nodes[0].is_root

    def test_leaf(self, simple_tree):
        nodes = get_target_nodes(simple_tree, 'leaf')
        assert all(n.is_leaf for n in nodes)
        assert len(nodes) == 4

    def test_intnode(self, simple_tree):
        nodes = get_target_nodes(simple_tree, 'intnode')
        assert all(not n.is_leaf for n in nodes)

    def test_invalid_target_raises(self, simple_tree):
        with pytest.raises(ValueError, match='Unknown target'):
            get_target_nodes(simple_tree, 'unknown_target')


class TestIsRooted:
    def test_rooted(self, simple_tree):
        assert is_rooted(simple_tree) is True

    def test_unrooted(self, unrooted_tree):
        assert is_rooted(unrooted_tree) is False

    def test_single_leaf_tree_is_rooted(self):
        tree = Tree('A;', parser=1)
        assert is_rooted(tree) is True

    def test_polytomy_root_tree_is_unrooted(self):
        tree = Tree('(A:1,B:1,C:1,D:1);', parser=1)
        assert is_rooted(tree) is False
