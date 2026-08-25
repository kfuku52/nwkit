import pytest

from nwkit.clade_mapping import canonical_split
from nwkit.root import root_main
from nwkit.util import is_rooted, read_tree
from tests.helpers import make_args, safe_get_distance
from tests.root_test_support import (
    install_fake_ncbi,
    install_fake_opentree,
    install_fake_timetree,
)


class TestRootMain:
    def test_writes_preserved_nhx_properties_after_rerooting(
        self, tmp_nwk, tmp_outfile
    ):
        path = tmp_nwk(
            "(((A:2[&&NHX:tip_tag=tip_A],B:3)AB:5[&&NHX:edge_tag=ab],C:7)"
            "ABC:11[&&NHX:root_edge_tag=same],"
            "(D:13,(E:17,F:19)EF:23[&&NHX:edge_tag=ef])"
            "DEF:29[&&NHX:root_edge_tag=same])ROOT:0[&&NHX:root_tag=root_value];"
        )
        args = make_args(
            infile=path,
            outfile=tmp_outfile,
            format="1",
            outformat="1",
            method="outgroup",
            outgroup="A",
        )

        root_main(args)

        tree = read_tree(tmp_outfile, format="1", quoted_node_names=True, quiet=True)
        assert tree.props["root_tag"] == "root_value"
        leaf_a = next(leaf for leaf in tree.leaves() if leaf.name == "A")
        assert leaf_a.props["tip_tag"] == "tip_A"
        all_taxa = frozenset(tree.leaf_names())
        internal_properties = {
            canonical_split(
                frozenset(node.leaf_names()),
                all_taxa - frozenset(node.leaf_names()),
            ): node.props.get("edge_tag")
            for node in tree.traverse()
            if not node.is_leaf and not node.is_root
        }
        assert (
            internal_properties[
                canonical_split(frozenset({"A", "B"}), frozenset({"C", "D", "E", "F"}))
            ]
            == "ab"
        )
        assert (
            internal_properties[
                canonical_split(frozenset({"E", "F"}), frozenset({"A", "B", "C", "D"}))
            ]
            == "ef"
        )

    def test_midpoint(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk("(A:1,B:5,(C:2,D:4):2);")
        args = make_args(
            infile=path,
            outfile=tmp_outfile,
            method="midpoint",
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format="auto", quoted_node_names=True, quiet=True)
        assert is_rooted(tree)

    def test_outgroup(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk("(A:1,(B:1,(C:1,D:1):1):1);")
        args = make_args(
            infile=path,
            outfile=tmp_outfile,
            method="outgroup",
            outgroup="A",
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format="auto", quoted_node_names=True, quiet=True)
        assert is_rooted(tree)
        assert set(tree.leaf_names()) == {"A", "B", "C", "D"}

    def test_outgroup_requires_label(self, tmp_nwk):
        path = tmp_nwk("(A:1,(B:1,(C:1,D:1):1):1);")
        args = make_args(
            infile=path,
            outfile="-",
            method="outgroup",
            outgroup="",
        )
        with pytest.raises(ValueError, match="outgroup"):
            root_main(args)

    def test_outgroup_requires_label_none(self, tmp_nwk):
        path = tmp_nwk("(A:1,(B:1,(C:1,D:1):1):1);")
        args = make_args(
            infile=path,
            outfile="-",
            method="outgroup",
            outgroup=None,
        )
        with pytest.raises(ValueError, match="outgroup"):
            root_main(args)

    def test_unknown_method_raises(self, tmp_nwk):
        path = tmp_nwk("(A:1,(B:1,(C:1,D:1):1):1);")
        args = make_args(
            infile=path,
            outfile="-",
            method="unknown",
        )
        with pytest.raises(ValueError, match="Unknown rooting method"):
            root_main(args)

    def test_transfer(self, tmp_nwk, tmp_outfile):
        path1 = tmp_nwk("(A:1,B:1,(C:1,D:1):1);", "tree1.nwk")
        path2 = tmp_nwk("((A:1,B:1):1,(C:1,D:1):1);", "tree2.nwk")
        args = make_args(
            infile=path1,
            outfile=tmp_outfile,
            method="transfer",
            infile2=path2,
            format2="auto",
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format="auto", quoted_node_names=True, quiet=True)
        assert is_rooted(tree)

    def test_transfer_single_leaf_tree(self, tmp_nwk, tmp_outfile):
        path1 = tmp_nwk("A;", "tree1.nwk")
        path2 = tmp_nwk("A;", "tree2.nwk")
        args = make_args(
            infile=path1,
            outfile=tmp_outfile,
            method="transfer",
            infile2=path2,
            format="1",
            format2="1",
            outformat="1",
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format="1", quoted_node_names=True, quiet=True)
        assert set(tree.leaf_names()) == {"A"}

    def test_transfer_mismatched_raises(self, tmp_nwk):
        path1 = tmp_nwk("(A:1,B:1,(C:1,D:1):1);", "tree1.nwk")
        path2 = tmp_nwk("((A:1,B:1):1,(C:1,E:1):1);", "tree2.nwk")
        args = make_args(
            infile=path1,
            outfile="-",
            method="transfer",
            infile2=path2,
            format2="auto",
        )
        with pytest.raises(Exception, match="Leaf labels"):
            root_main(args)

    def test_transfer_requires_rooted_infile2(self, tmp_nwk):
        path1 = tmp_nwk("(A:1,B:1,(C:1,D:1):1);", "tree1.nwk")
        path2 = tmp_nwk("(A:1,B:1,C:1,D:1);", "tree2.nwk")
        args = make_args(
            infile=path1,
            outfile="-",
            method="transfer",
            infile2=path2,
            format2="auto",
        )
        with pytest.raises(ValueError, match="must be rooted"):
            root_main(args)

    def test_transfer_requires_bifurcating_root_infile2(self, tmp_nwk):
        path1 = tmp_nwk("(A:1,B:1,(C:1,D:1):1);", "tree1.nwk")
        path2 = tmp_nwk("(((A:1,B:1):1,C:1):1);", "tree2.nwk")
        args = make_args(
            infile=path1,
            outfile="-",
            method="transfer",
            infile2=path2,
            format2="auto",
        )
        with pytest.raises(ValueError, match="exactly two children"):
            root_main(args)

    def test_transfer_incompatible_root_bipartition_raises(self, tmp_nwk):
        path1 = tmp_nwk("((A:1,B:1):1,(C:1,D:1):1);", "tree1.nwk")
        path2 = tmp_nwk("((A:1,C:1):1,(B:1,D:1):1);", "tree2.nwk")
        args = make_args(
            infile=path1,
            outfile="-",
            method="transfer",
            infile2=path2,
            format2="auto",
        )
        with pytest.raises(ValueError, match="No root bipartition"):
            root_main(args)

    def test_transfer_requires_infile2(self, tmp_nwk):
        path1 = tmp_nwk("(A:1,B:1,(C:1,D:1):1);", "tree1.nwk")
        args = make_args(
            infile=path1,
            outfile="-",
            method="transfer",
            infile2="",
            format2="auto",
        )
        with pytest.raises(ValueError, match="infile2"):
            root_main(args)

    def test_transfer_duplicate_leaf_names_raise(self, tmp_nwk):
        path1 = tmp_nwk("((A:1,A:2):1,B:1);", "tree1.nwk")
        path2 = tmp_nwk("((A:1,A:2):1,B:1);", "tree2.nwk")
        args = make_args(
            infile=path1,
            outfile="-",
            method="transfer",
            infile2=path2,
            format2="auto",
        )
        with pytest.raises(ValueError, match="Duplicated leaf labels"):
            root_main(args)

    def test_transfer_empty_leaf_labels_raise(self, tmp_nwk):
        path1 = tmp_nwk("(A:1,(:1,B:1):1);", "tree1.nwk")
        path2 = tmp_nwk("((A:1,B:1):1,:1);", "tree2.nwk")
        args = make_args(
            infile=path1,
            outfile="-",
            method="transfer",
            infile2=path2,
            format2="auto",
        )
        with pytest.raises(ValueError, match="Empty leaf labels"):
            root_main(args)

    def test_mad(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk("((A:1,B:2):1,(C:3,(D:1,E:2):1):1);")
        args = make_args(
            infile=path,
            outfile=tmp_outfile,
            method="mad",
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format="auto", quoted_node_names=True, quiet=True)
        assert is_rooted(tree)
        assert set(tree.leaf_names()) == {"A", "B", "C", "D", "E"}

    def test_mad_with_too_few_leaves_raises(self, tmp_nwk):
        path = tmp_nwk("(A:1,B:1);")
        args = make_args(
            infile=path,
            outfile="-",
            method="mad",
        )
        with pytest.raises(ValueError, match="at least 3 leaves"):
            root_main(args)

    def test_mv(self, tmp_nwk, tmp_outfile):
        path = tmp_nwk("((A:1,B:2):1,(C:3,(D:1,E:2):1):1);")
        args = make_args(
            infile=path,
            outfile=tmp_outfile,
            method="mv",
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format="auto", quoted_node_names=True, quiet=True)
        assert is_rooted(tree)
        assert set(tree.leaf_names()) == {"A", "B", "C", "D", "E"}

    def test_reconciliation(self, tmp_nwk, tmp_outfile):
        gene_path = tmp_nwk(
            "(A_a_g1:2,B_b_g1:3,(C_c_g1:4,D_d_g1:5):6);",
            "gene.nwk",
        )
        species_path = tmp_nwk(
            "((A_a:1,B_b:1):1,(C_c:1,D_d:1):1);",
            "species.nwk",
        )
        args = make_args(
            infile=gene_path,
            outfile=tmp_outfile,
            method="reconciliation",
            species_tree=species_path,
            species_tree_format="auto",
            duplication_cost=1.0,
            loss_cost=1.0,
        )

        root_main(args)

        tree = read_tree(tmp_outfile, format="auto", quoted_node_names=True, quiet=True)
        assert {frozenset(child.leaf_names()) for child in tree.get_children()} == {
            frozenset({"A_a_g1", "B_b_g1"}),
            frozenset({"C_c_g1", "D_d_g1"}),
        }

    def test_reconciliation_requires_species_tree(self, tmp_nwk):
        gene_path = tmp_nwk("(A_a_g1:1,B_b_g1:1);", "gene.nwk")
        args = make_args(
            infile=gene_path,
            outfile="-",
            method="reconciliation",
            species_tree=None,
        )

        with pytest.raises(ValueError, match="species-tree.*required"):
            root_main(args)

    def test_reconciliation_large_finite_score_still_writes_tree(
        self,
        tmp_nwk,
        tmp_outfile,
        capsys,
    ):
        gene_path = tmp_nwk(
            "(A_a_g1:1,A_a_g2:1,(A_a_g3:1,A_a_g4:1):1);",
            "gene.nwk",
        )
        species_path = tmp_nwk("(A_a:1,B_b:1);", "species.nwk")
        args = make_args(
            infile=gene_path,
            outfile=tmp_outfile,
            method="reconciliation",
            species_tree=species_path,
            species_tree_format="auto",
            duplication_cost=1e308,
            loss_cost=0,
        )

        root_main(args)

        rooted = read_tree(
            tmp_outfile,
            format="auto",
            quoted_node_names=True,
            quiet=True,
        )
        assert len(list(rooted.leaves())) == 4
        assert "score: 3e+308" in capsys.readouterr().err

    def test_taxonomy(self, monkeypatch, tmp_nwk, tmp_outfile):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={
                "Homo sapiens": 9606,
                "Pan troglodytes": 9598,
                "Arabidopsis thaliana": 3702,
                "Oryza sativa": 4530,
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
        path = tmp_nwk(
            "(Homo_sapiens:1,Pan_troglodytes:1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);"
        )
        args = make_args(
            infile=path,
            outfile=tmp_outfile,
            method="taxonomy",
            taxonomy_source="ncbi",
            taxid_tsv=None,
            rank="no",
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format="auto", quoted_node_names=True, quiet=True)
        child_leaf_sets = [set(child.leaf_names()) for child in tree.get_children()]
        assert {"Homo_sapiens", "Pan_troglodytes"} in child_leaf_sets
        assert {"Arabidopsis_thaliana", "Oryza_sativa"} in child_leaf_sets

    def test_taxonomy_source_chain(self, monkeypatch, tmp_nwk, tmp_outfile):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={},
            lineage_by_taxid={},
        )
        install_fake_timetree(
            monkeypatch,
            upload_html='<div id="prunetree-msg-box"></div>',
            newick_text="((Homo_sapiens:1,Pan_troglodytes:1):10,(Arabidopsis_thaliana:1,Oryza_sativa:1):20);",
        )
        path = tmp_nwk(
            "(Homo_sapiens:1,Pan_troglodytes:1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);"
        )
        args = make_args(
            infile=path,
            outfile=tmp_outfile,
            method="taxonomy",
            taxonomy_source="ncbi,timetree",
            taxid_tsv=None,
            rank="no",
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format="auto", quoted_node_names=True, quiet=True)
        child_leaf_sets = [set(child.leaf_names()) for child in tree.get_children()]
        assert {"Homo_sapiens", "Pan_troglodytes"} in child_leaf_sets
        assert {"Arabidopsis_thaliana", "Oryza_sativa"} in child_leaf_sets

    def test_taxonomy_opentree(self, monkeypatch, tmp_nwk, tmp_outfile):
        install_fake_opentree(
            monkeypatch,
            tnrs_json={
                "results": [
                    {
                        "matches": [
                            {
                                "is_approximate_match": False,
                                "taxon": {
                                    "ott_id": 1,
                                    "is_suppressed_from_synth": False,
                                },
                            }
                        ]
                    },
                    {
                        "matches": [
                            {
                                "is_approximate_match": False,
                                "taxon": {
                                    "ott_id": 2,
                                    "is_suppressed_from_synth": False,
                                },
                            }
                        ]
                    },
                    {
                        "matches": [
                            {
                                "is_approximate_match": False,
                                "taxon": {
                                    "ott_id": 3,
                                    "is_suppressed_from_synth": False,
                                },
                            }
                        ]
                    },
                    {
                        "matches": [
                            {
                                "is_approximate_match": False,
                                "taxon": {
                                    "ott_id": 4,
                                    "is_suppressed_from_synth": False,
                                },
                            }
                        ]
                    },
                ],
            },
            induced_subtree_json={
                "broken": {},
                "newick": "((Homo_sapiens,Pan_troglodytes)Primates,(Arabidopsis_thaliana,Oryza_sativa)Mesangiospermae)Eukaryota;",
            },
        )
        path = tmp_nwk(
            "(Homo_sapiens:1,Pan_troglodytes:1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);"
        )
        args = make_args(
            infile=path,
            outfile=tmp_outfile,
            method="taxonomy",
            taxonomy_source="opentree",
            taxid_tsv=None,
            rank="no",
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format="auto", quoted_node_names=True, quiet=True)
        child_leaf_sets = [set(child.leaf_names()) for child in tree.get_children()]
        assert {"Homo_sapiens", "Pan_troglodytes"} in child_leaf_sets
        assert {"Arabidopsis_thaliana", "Oryza_sativa"} in child_leaf_sets

    def test_taxonomy_defaults(self, monkeypatch, tmp_nwk, tmp_outfile):
        install_fake_ncbi(
            monkeypatch,
            name_to_taxid={},
            lineage_by_taxid={},
        )
        install_fake_opentree(
            monkeypatch,
            tnrs_json={
                "results": [
                    {
                        "matches": [
                            {
                                "is_approximate_match": False,
                                "taxon": {
                                    "ott_id": 1,
                                    "is_suppressed_from_synth": False,
                                },
                            }
                        ]
                    },
                    {
                        "matches": [
                            {
                                "is_approximate_match": False,
                                "taxon": {
                                    "ott_id": 2,
                                    "is_suppressed_from_synth": False,
                                },
                            }
                        ]
                    },
                    {
                        "matches": [
                            {
                                "is_approximate_match": False,
                                "taxon": {
                                    "ott_id": 3,
                                    "is_suppressed_from_synth": False,
                                },
                            }
                        ]
                    },
                    {
                        "matches": [
                            {
                                "is_approximate_match": False,
                                "taxon": {
                                    "ott_id": 4,
                                    "is_suppressed_from_synth": False,
                                },
                            }
                        ]
                    },
                ],
            },
            induced_subtree_json={
                "broken": {},
                "newick": "((Homo_sapiens,Pan_troglodytes)Primates,(Arabidopsis_thaliana,Oryza_sativa)Mesangiospermae)Eukaryota;",
            },
        )
        path = tmp_nwk(
            "(Homo_sapiens:1,Pan_troglodytes:1,(Arabidopsis_thaliana:1,Oryza_sativa:1):1);"
        )
        args = make_args(
            infile=path,
            outfile=tmp_outfile,
            method="taxonomy",
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format="auto", quoted_node_names=True, quiet=True)
        child_leaf_sets = [set(child.leaf_names()) for child in tree.get_children()]
        assert {"Homo_sapiens", "Pan_troglodytes"} in child_leaf_sets
        assert {"Arabidopsis_thaliana", "Oryza_sativa"} in child_leaf_sets

    def test_wiki_outgroup_single(self, tmp_nwk, tmp_outfile):
        """Wiki example: nwkit root --method outgroup --outgroup a

        Input:  (((a:1,b:1):1,c:1):1,((d:1,e:1),f:1):1):0;
        Output: (a:0.5,(b:1,(c:1,((d:1,e:1):1,f:1):2):1):0.5):0;
        """
        path = tmp_nwk("(((a:1,b:1):1,c:1):1,((d:1,e:1),f:1):1):0;")
        args = make_args(
            infile=path,
            outfile=tmp_outfile,
            method="outgroup",
            outgroup="a",
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format="auto", quoted_node_names=True, quiet=True)
        assert is_rooted(tree)
        assert set(tree.leaf_names()) == {"a", "b", "c", "d", "e", "f"}
        # 'a' should be in one of the two subroot clades, by itself
        children = tree.get_children()
        child_leaf_sets = [set(c.leaf_names()) for c in children]
        assert {"a"} in child_leaf_sets or any(
            "a" in s and len(s) == 1 for s in child_leaf_sets
        )

    def test_wiki_outgroup_multiple(self, tmp_nwk, tmp_outfile):
        """Wiki example: nwkit root --method outgroup --outgroup a,b

        Input:  (((a:1,b:1):1,c:1):1,((d:1,e:1),f:1):1):0;
        Output: ((a:1,b:1):0.5,(c:1,((d:1,e:1):1,f:1):2):0.5):0;
        """
        path = tmp_nwk("(((a:1,b:1):1,c:1):1,((d:1,e:1),f:1):1):0;")
        args = make_args(
            infile=path,
            outfile=tmp_outfile,
            method="outgroup",
            outgroup="a,b",
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format="auto", quoted_node_names=True, quiet=True)
        assert is_rooted(tree)
        assert set(tree.leaf_names()) == {"a", "b", "c", "d", "e", "f"}
        # {a,b} should be one of the two subroot clades
        children = tree.get_children()
        child_leaf_sets = [set(c.leaf_names()) for c in children]
        assert {"a", "b"} in child_leaf_sets

    def test_wiki_root_transfer(self, tmp_nwk, tmp_outfile):
        """Wiki example: nwkit root --method transfer

        input1.nwk: (((a:1,b:1):1,c:1):1,((d:1,e:1),f:1):1):0;
        input2.nwk: ((((a:1,b:1):1,c:1),f:1):1,(d:3,e:3):1):0;
        Output: ((d:1,e:1):0.5,(f:1,(c:1,(b:1,a:1):1):2):0.5)Root:0;
        """
        path1 = tmp_nwk("(((a:1,b:1):1,c:1):1,((d:1,e:1),f:1):1):0;", "tree1.nwk")
        path2 = tmp_nwk("((((a:1,b:1):1,c:1),f:1):1,(d:3,e:3):1):0;", "tree2.nwk")
        args = make_args(
            infile=path1,
            outfile=tmp_outfile,
            method="transfer",
            infile2=path2,
            format2="auto",
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format="auto", quoted_node_names=True, quiet=True)
        assert is_rooted(tree)
        assert set(tree.leaf_names()) == {"a", "b", "c", "d", "e", "f"}
        # {d,e} should be in one subroot clade (like the source tree)
        children = tree.get_children()
        child_leaf_sets = [set(c.leaf_names()) for c in children]
        assert {"d", "e"} in child_leaf_sets

    def test_wiki_midpoint_rooting(self, tmp_nwk, tmp_outfile):
        """Wiki example: nwkit root --method midpoint

        Input:  ((((a:5,b:1):1,c:3):1,f:1):1,(d:1,e:1):1):0;
        The exact diameter midpoint is the endpoint of the length-5 branch to a.
        """
        # Note: Added explicit :1 to internal node that was missing dist in original wiki example
        path = tmp_nwk("((((a:5,b:1):1,c:3):1,f:1):1,(d:1,e:1):1):0;")
        args = make_args(
            infile=path,
            outfile=tmp_outfile,
            method="midpoint",
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format="auto", quoted_node_names=True, quiet=True)
        assert is_rooted(tree)
        assert set(tree.leaf_names()) == {"a", "b", "c", "d", "e", "f"}
        # Verify exact root-to-tip distances from the diameter midpoint.
        dists = {l.name: safe_get_distance(tree, tree, l) for l in tree.leaves()}
        assert abs(dists["a"] - 5.0) < 1e-6
        assert abs(dists["b"] - 1.0) < 1e-6
        assert abs(dists["c"] - 4.0) < 1e-6
        assert abs(dists["f"] - 3.0) < 1e-6
        assert abs(dists["d"] - 5.0) < 1e-6
        assert abs(dists["e"] - 5.0) < 1e-6
        assert sorted(child.dist for child in tree.get_children()) == [0.0, 5.0]

    def test_wiki_outgroup_single_exact_distances(self, tmp_nwk, tmp_outfile):
        """Wiki outgroup single: verify exact root-to-tip distances.

        Input:  (((a:1,b:1):1,c:1):1,((d:1,e:1):1,f:1):1):0;
        Output: (a:0.5,(b:1,(c:1,((d:1,e:1):1,f:1):2):1):0.5):0;
        """
        # Note: Added explicit :1 to internal node that was missing dist in original wiki example
        path = tmp_nwk("(((a:1,b:1):1,c:1):1,((d:1,e:1):1,f:1):1):0;")
        args = make_args(
            infile=path,
            outfile=tmp_outfile,
            method="outgroup",
            outgroup="a",
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format="auto", quoted_node_names=True, quiet=True)
        dists = {l.name: safe_get_distance(tree, tree, l) for l in tree.leaves()}
        assert abs(dists["a"] - 0.5) < 1e-6
        assert abs(dists["b"] - 1.5) < 1e-6
        assert abs(dists["c"] - 2.5) < 1e-6
        assert abs(dists["f"] - 4.5) < 1e-6
        assert abs(dists["d"] - 5.5) < 1e-6
        assert abs(dists["e"] - 5.5) < 1e-6

    def test_wiki_outgroup_multiple_exact_distances(self, tmp_nwk, tmp_outfile):
        """Wiki outgroup multiple: verify exact root-to-tip distances.

        Input:  (((a:1,b:1):1,c:1):1,((d:1,e:1):1,f:1):1):0;
        Output: ((a:1,b:1):0.5,(c:1,((d:1,e:1):1,f:1):2):0.5):0;
        """
        # Note: Added explicit :1 to internal node that was missing dist in original wiki example
        path = tmp_nwk("(((a:1,b:1):1,c:1):1,((d:1,e:1):1,f:1):1):0;")
        args = make_args(
            infile=path,
            outfile=tmp_outfile,
            method="outgroup",
            outgroup="a,b",
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format="auto", quoted_node_names=True, quiet=True)
        dists = {l.name: safe_get_distance(tree, tree, l) for l in tree.leaves()}
        assert abs(dists["a"] - 1.5) < 1e-6
        assert abs(dists["b"] - 1.5) < 1e-6
        assert abs(dists["c"] - 1.5) < 1e-6
        assert abs(dists["f"] - 3.5) < 1e-6
        assert abs(dists["d"] - 4.5) < 1e-6
        assert abs(dists["e"] - 4.5) < 1e-6

    def test_wiki_transfer_exact_distances(self, tmp_nwk, tmp_outfile):
        """Wiki transfer: verify exact root-to-tip distances.

        input1: (((a:1,b:1):1,c:1):1,((d:1,e:1):1,f:1):1):0;
        input2: ((((a:1,b:1):1,c:1):1,f:1):1,(d:3,e:3):1):0;
        Output: ((d:1,e:1):0.5,(f:1,(c:1,(b:1,a:1):1):2):0.5)Root:0;
        """
        # Note: Added explicit :1 to internal nodes that were missing dist in original wiki example
        path1 = tmp_nwk("(((a:1,b:1):1,c:1):1,((d:1,e:1):1,f:1):1):0;", "tree1.nwk")
        path2 = tmp_nwk("((((a:1,b:1):1,c:1):1,f:1):1,(d:3,e:3):1):0;", "tree2.nwk")
        args = make_args(
            infile=path1,
            outfile=tmp_outfile,
            method="transfer",
            infile2=path2,
            format2="auto",
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format="auto", quoted_node_names=True, quiet=True)
        dists = {l.name: safe_get_distance(tree, tree, l) for l in tree.leaves()}
        assert abs(dists["d"] - 1.5) < 1e-6
        assert abs(dists["e"] - 1.5) < 1e-6
        assert abs(dists["f"] - 1.5) < 1e-6
        assert abs(dists["c"] - 3.5) < 1e-6
        assert abs(dists["a"] - 4.5) < 1e-6
        assert abs(dists["b"] - 4.5) < 1e-6

    def test_wiki_mv_exact_distances(self, tmp_nwk, tmp_outfile):
        """Wiki MV rooting: verify exact root-to-tip distances.

        Input: ((A:1,B:2):1,(C:3,(D:1,E:2):1):1);
        Output: ((C:3,(D:1,E:2):1):0.416667,(A:1,B:2):1.58333);
        Root at x=5/12 from (C,(D,E)) side.
        """
        path = tmp_nwk("((A:1,B:2):1,(C:3,(D:1,E:2):1):1);")
        args = make_args(
            infile=path,
            outfile=tmp_outfile,
            method="mv",
        )
        root_main(args)
        tree = read_tree(tmp_outfile, format="auto", quoted_node_names=True, quiet=True)
        dists = {l.name: safe_get_distance(tree, tree, l) for l in tree.leaves()}
        assert abs(dists["C"] - 41 / 12) < 1e-4
        assert abs(dists["D"] - 29 / 12) < 1e-4
        assert abs(dists["E"] - 41 / 12) < 1e-4
        assert abs(dists["A"] - 31 / 12) < 1e-4
        assert abs(dists["B"] - 43 / 12) < 1e-4
