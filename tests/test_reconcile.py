import pytest
from ete4 import Tree

from nwkit.reconcile import (
    _report_unmatched_species,
    build_reconciliation_table,
    lca_duplication_loss_contribution,
)

SPECIES_TREE = "((A_a:1,B_b:1):1,C_c:2);"
GENE_TREE_WITH_TWO_COPIES = (
    "(((A_a_g1:1,B_b_g1:1):1,C_c_g1:2):1,((A_a_g2:1,B_b_g2:1):1,C_c_g2:2):1);"
)


def _species_mapping(tree):
    return {
        str(leaf.name): "_".join(str(leaf.name).split("_")[:2])
        for leaf in tree.leaves()
    }


def test_lca_duplication_loss_contribution_counts_species_edges():
    species_tree = Tree("((A:1,B:1)AB:1,C:1)ABC;", parser=1)
    nodes = {str(node.name): node for node in species_tree.traverse()}
    depths = {
        species_tree: 0,
        nodes["AB"]: 1,
        nodes["A"]: 2,
        nodes["B"]: 2,
        nodes["C"]: 1,
    }

    assert lca_duplication_loss_contribution(
        species_tree,
        [nodes["A"], nodes["C"]],
        depths,
    ) == (0, 1)
    assert lca_duplication_loss_contribution(
        species_tree,
        [species_tree, nodes["A"]],
        depths,
    ) == (1, 2)
    assert lca_duplication_loss_contribution(
        species_tree,
        (node for node in [species_tree, nodes["A"]]),
        depths,
    ) == (1, 2)


def test_lca_reconciliation_keeps_parallel_speciation_events_as_separate_rows():
    gene_tree = Tree(GENE_TREE_WITH_TWO_COPIES, parser=1)
    species_tree = Tree(SPECIES_TREE, parser=1)
    table = build_reconciliation_table(
        gene_tree,
        species_tree,
        _species_mapping(gene_tree),
        event_source="lca",
    )

    internal = table[table["node_class"].isin(["root", "intnode"])]
    assert internal["event_type"].tolist().count("duplication") == 1
    speciation = internal[internal["event_type"] == "speciation"]
    assert len(speciation) == 4
    assert set(speciation["eligible"]) == {"yes"}
    assert speciation["species_branch_id"].value_counts().to_dict() == {0: 2, 1: 2}
    assert speciation["lineage_id"].nunique() == 2
    assert all(speciation.groupby("lineage_id").size() == 2)


def test_lca_reconciliation_allows_gene_loss_without_fabricating_an_event():
    gene_tree = Tree(
        "(((A_a_g1:1,B_b_g1:1):1,C_c_g1:2):1,(A_a_g2:1,C_c_g2:2):1);",
        parser=1,
    )
    species_tree = Tree(SPECIES_TREE, parser=1)
    table = build_reconciliation_table(
        gene_tree,
        species_tree,
        _species_mapping(gene_tree),
        event_source="lca",
    )

    speciation = table[
        (table["event_type"] == "speciation") & (table["eligible"] == "yes")
    ]
    assert speciation["species_branch_id"].value_counts().to_dict() == {0: 2, 1: 1}
    partial = speciation[speciation["coverage_status"] == "partial"]
    assert len(partial) == 1
    assert partial.iloc[0]["numerator_species_coverage"] == pytest.approx(0.5)


def test_nhx_events_are_used_without_silent_lca_replacement():
    gene_tree = Tree(
        "(((A_a_g1:1,B_b_g1:1):1[&&NHX:D=Y],C_c_g1:2):1[&&NHX:D=N],"
        "((A_a_g2:1,B_b_g2:1):1[&&NHX:D=N],C_c_g2:2):1[&&NHX:D=N])"
        "[&&NHX:D=N];",
        parser=1,
    )
    species_tree = Tree(SPECIES_TREE, parser=1)
    table = build_reconciliation_table(
        gene_tree,
        species_tree,
        _species_mapping(gene_tree),
        event_source="nhx",
    )

    root = table[table["gene_branch_id"] == 0].iloc[0]
    assert root["event_type"] == "speciation"
    assert root["eligible"] == "no"
    assert root["reason"] == "children_do_not_match_distinct_species_clades"
    annotated_duplication = table[table["descendant_taxa"] == "A_a_g1,B_b_g1"].iloc[0]
    assert annotated_duplication["event_type"] == "duplication"
    assert annotated_duplication["eligible"] == "no"


def test_missing_nhx_event_is_reported_as_unresolved():
    gene_tree = Tree("(A_a_g1:1,B_b_g1:1);", parser=1)
    species_tree = Tree("(A_a:1,B_b:1);", parser=1)
    table = build_reconciliation_table(
        gene_tree,
        species_tree,
        _species_mapping(gene_tree),
        event_source="nhx",
    )

    root = table.iloc[0]
    assert root["event_type"] == "unresolved"
    assert root["event_status"] == "unresolved"
    assert root["reason"] == "missing_nhx_D"


def test_generax_transfer_annotation_retains_source_and_destination():
    gene_tree = Tree("(A_a_g1:1,B_b_g1:1):1[&&NHX:D=N:H=Y@A_a@B_b];", parser=1)
    species_tree = Tree("(A_a:1,B_b:1);", parser=1)
    table = build_reconciliation_table(
        gene_tree,
        species_tree,
        _species_mapping(gene_tree),
        event_source="nhx",
    )

    root = table.iloc[0]
    assert root["event_type"] == "transfer"
    assert root["event_status"] == "resolved"
    assert root["transfer_source_species"] == "A_a"
    assert root["transfer_destination_species"] == "B_b"


def test_generax_S_annotation_controls_species_mapping_for_transfer():
    gene_tree = Tree("(A_g1:1,C_g1:1)[&&NHX:S=A:D=N:H=Y@A@C];", parser=1)
    species_tree = Tree("((A:1,B:1)AB:1,C:2)ABC;", parser=1)
    table = build_reconciliation_table(
        gene_tree,
        species_tree,
        {"A_g1": "A", "C_g1": "C"},
        event_source="nhx",
    )

    root = table.iloc[0]
    assert root["event_type"] == "transfer"
    assert root["species_name"] == "A"
    assert root["species_event_taxa"] == "A"


def test_invalid_generax_S_annotation_is_not_replaced_by_lca_mapping():
    gene_tree = Tree("(A_g1:1,B_g1:1)[&&NHX:S=missing:D=N:H=N];", parser=1)
    species_tree = Tree("(A:1,B:1)AB;", parser=1)
    table = build_reconciliation_table(
        gene_tree,
        species_tree,
        {"A_g1": "A", "B_g1": "B"},
        event_source="nhx",
    )

    root = table.iloc[0]
    assert root["mapping_status"] == "unmapped"
    assert root["species_event_id"] == ""
    assert root["reason"] == "nhx_S_not_in_species_tree"


def test_missing_generax_S_annotation_explicitly_falls_back_to_lca_mapping():
    gene_tree = Tree("(A_g1:1,B_g1:1)[&&NHX:D=N:H=N];", parser=1)
    species_tree = Tree("(A:1,B:1)AB;", parser=1)
    table = build_reconciliation_table(
        gene_tree,
        species_tree,
        {"A_g1": "A", "B_g1": "B"},
        event_source="nhx",
    )

    root = table.iloc[0]
    assert root["mapping_status"] == "mapped"
    assert root["species_name"] == "AB"
    assert root["eligible"] == "yes"


def test_malformed_generax_transfer_annotation_is_unresolved():
    gene_tree = Tree("(A_a_g1:1,B_b_g1:1)[&&NHX:D=N:H=Y@A_a];", parser=1)
    species_tree = Tree("(A_a:1,B_b:1);", parser=1)
    table = build_reconciliation_table(
        gene_tree,
        species_tree,
        _species_mapping(gene_tree),
        event_source="nhx",
    )
    assert table.iloc[0]["reason"] == "invalid_nhx_H"


def test_generax_transfer_endpoints_are_validated_against_species_tree():
    gene_tree = Tree("(A_g1:1,B_g1:1)[&&NHX:S=A:D=N:H=Y@A@missing];", parser=1)
    species_tree = Tree("(A:1,B:1)AB;", parser=1)
    table = build_reconciliation_table(
        gene_tree,
        species_tree,
        {"A_g1": "A", "B_g1": "B"},
        event_source="nhx",
    )

    root = table.iloc[0]
    assert root["event_type"] == "unresolved"
    assert root["reason"] == "nhx_H_endpoint_not_in_species_tree"
    assert root["transfer_source_species"] == ""
    assert root["transfer_destination_species"] == ""


def test_generax_transfer_source_must_match_S_annotation():
    gene_tree = Tree("(A_g1:1,B_g1:1)[&&NHX:S=A:D=N:H=Y@B@A];", parser=1)
    species_tree = Tree("(A:1,B:1)AB;", parser=1)
    table = build_reconciliation_table(
        gene_tree,
        species_tree,
        {"A_g1": "A", "B_g1": "B"},
        event_source="nhx",
    )

    root = table.iloc[0]
    assert root["event_type"] == "unresolved"
    assert root["reason"] == "nhx_H_source_conflicts_with_S"


def test_clade_and_species_event_ids_are_invariant_to_child_order():
    gene1 = Tree("((A_a_g1:1,B_b_g1:1):1,C_c_g1:2);", parser=1)
    gene2 = Tree("(C_c_g1:2,(B_b_g1:1,A_a_g1:1):1);", parser=1)
    species1 = Tree("((A_a:1,B_b:1):1,C_c:2);", parser=1)
    species2 = Tree("(C_c:2,(B_b:1,A_a:1):1);", parser=1)
    table1 = build_reconciliation_table(
        gene1, species1, _species_mapping(gene1), tree_id="OG1"
    ).set_index("descendant_taxa")
    table2 = build_reconciliation_table(
        gene2, species2, _species_mapping(gene2), tree_id="OG1"
    ).set_index("descendant_taxa")

    assert table1["gene_clade_id"].to_dict() == table2["gene_clade_id"].to_dict()
    assert table1["species_event_id"].to_dict() == table2["species_event_id"].to_dict()
    assert table1["tree_id"].unique().tolist() == ["OG1"]


def test_collapsed_non_speciation_boundary_starts_a_new_lineage():
    gene_tree = Tree("((A_a_g1:1,B_b_g1:1):1,C_c_g1:2);", parser=1)
    species_tree = Tree("((A_a:1,B_b:1):1,C_c:2);", parser=1)
    ab = gene_tree.common_ancestor(["A_a_g1", "B_b_g1"])
    ab.props["NWKIT_COLLAPSED_EVENT_BOUNDARY"] = "Y"
    table = build_reconciliation_table(
        gene_tree, species_tree, _species_mapping(gene_tree)
    )
    ab_row = table[table["descendant_taxa"] == "A_a_g1,B_b_g1"].iloc[0]
    root_row = table[table["node_class"] == "root"].iloc[0]
    assert ab_row["collapsed_event_boundary"] == "yes"
    assert ab_row["lineage_clade_id"] == ab_row["gene_clade_id"]
    assert ab_row["lineage_clade_id"] != root_row["lineage_clade_id"]


def test_unmatched_species_are_errors_by_default():
    species_tree = Tree("(A_a:1,B_b:1);", parser=1)
    with pytest.raises(ValueError, match="could not be mapped"):
        _report_unmatched_species(
            {"A_a_g1": "A_a", "X_x_g1": "X_x"},
            species_tree,
            policy="error",
        )


def test_reconciliation_species_mapping_must_exactly_cover_gene_tips():
    gene_tree = Tree("(A_g1:1,B_g1:1);", parser=1)
    species_tree = Tree("(A:1,B:1);", parser=1)
    with pytest.raises(ValueError, match="exactly cover"):
        build_reconciliation_table(
            gene_tree,
            species_tree,
            {"A_g1": "A", "extra": "B"},
        )


def test_reconciliation_rejects_mapping_keys_that_collide_as_strings():
    gene_tree = Tree("(1:1,2:1);", parser=1)
    species_tree = Tree("(A:1,B:1);", parser=1)
    with pytest.raises(ValueError, match="remain unique"):
        build_reconciliation_table(
            gene_tree,
            species_tree,
            {1: "A", "1": "A", 2: "B"},
        )
