import pandas as pd
import pytest
from ete4 import Tree

from nwkit.clade_index import CladeIndex
from nwkit.cli import main
from nwkit.contrast import (
    _read_reconciliation,
    build_contrast_table,
    calculate_contrasts,
)
from nwkit.reconcile import build_reconciliation_table


def test_calculate_contrasts_uses_standard_pic_branch_length_adjustment():
    tree = Tree("((A:1,B:1):1,C:2);", parser=1)
    contrasts = calculate_contrasts(tree, {"A": 1.0, "B": 3.0, "C": 4.0})

    ab = tree.common_ancestor(["A", "B"])
    assert contrasts[ab]["raw_contrast"] == pytest.approx(-2.0)
    assert contrasts[ab]["contrast_variance"] == pytest.approx(2.0)
    assert contrasts[ab]["standardized_contrast"] == pytest.approx(-(2.0**0.5))
    assert contrasts[tree]["raw_contrast"] == pytest.approx(-2.0)
    assert contrasts[tree]["contrast_variance"] == pytest.approx(3.5)
    assert contrasts[tree]["ancestral_estimate"] == pytest.approx(20.0 / 7.0)


def test_contrasts_support_transformed_tree_evolution_models():
    tree = Tree("((A:1,B:1):1,C:2);", parser=1)
    values = {"A": 1.0, "B": 3.0, "C": 4.0}
    contrasts = calculate_contrasts(
        tree,
        values,
        evolution_model="delta",
        evolution_parameter=2.0,
    )
    ab = tree.common_ancestor(["A", "B"])
    assert contrasts[ab]["contrast_variance"] == pytest.approx(3.0)
    table = build_contrast_table(
        tree,
        {"value": values},
        evolution_model="delta",
        evolution_parameter=2.0,
    )
    assert set(table["evolution_model"]) == {"delta"}
    assert set(table["evolution_parameter_name"]) == {"delta"}
    assert set(table["evolution_parameter"]) == {2.0}


def test_independent_contrasts_ignore_branch_lengths_and_report_not_applicable():
    tree = Tree("((A:0,B:1):1,C:2);", parser=1)
    table = build_contrast_table(
        tree,
        {"value": {"A": 1.0, "B": 3.0, "C": 4.0}},
        evolution_model="independent",
    )
    assert set(table["branch_length_mode"]) == {"not-applicable"}


def test_contrast_sign_is_invariant_to_newick_child_order():
    values = {"A": 1.0, "B": 3.0, "C": 4.0}
    tree1 = Tree("((A:1,B:1):1,C:2);", parser=1)
    tree2 = Tree("(C:2,(B:1,A:1):1);", parser=1)
    table1 = build_contrast_table(tree1, {"value": values})
    table2 = build_contrast_table(tree2, {"value": values})

    observed1 = table1.set_index("descendant_taxa")["standardized_contrast"]
    observed2 = table2.set_index("descendant_taxa")["standardized_contrast"]
    pd.testing.assert_series_equal(observed1.sort_index(), observed2.sort_index())


def test_reconciled_contrasts_keep_species_branch_and_lineage_clusters():
    species_tree = Tree("((A_a:1,B_b:1):1,C_c:2);", parser=1)
    gene_tree = Tree(
        "(((A_a_g1:1,B_b_g1:1):1,C_c_g1:2):1,((A_a_g2:1,B_b_g2:1):1,C_c_g2:2):1);",
        parser=1,
    )
    species_by_gene = {
        str(leaf.name): "_".join(str(leaf.name).split("_")[:2])
        for leaf in gene_tree.leaves()
    }
    reconciliation = build_reconciliation_table(
        gene_tree, species_tree, species_by_gene
    )
    reconciliation_by_id = {
        record["gene_clade_id"]: record for record in reconciliation.to_dict("records")
    }
    values = {
        "A_a_g1": 1.0,
        "B_b_g1": 2.0,
        "C_c_g1": 4.0,
        "A_a_g2": 2.0,
        "B_b_g2": 4.0,
        "C_c_g2": 8.0,
    }

    table = build_contrast_table(
        gene_tree,
        {"expression": values},
        reconciliation_by_id=reconciliation_by_id,
    )

    assert len(table) == 4
    assert set(table["event_type"]) == {"speciation"}
    assert set(table["eligible"]) == {"yes"}
    assert table["species_branch_id"].value_counts().to_dict() == {0: 2, 1: 2}
    assert table["lineage_id"].nunique() == 2


def test_reconcile_and_contrast_cli_end_to_end(tmp_path):
    species_tree = tmp_path / "species.nwk"
    gene_tree = tmp_path / "gene.nwk"
    trait = tmp_path / "expression.tsv"
    reconciliation = tmp_path / "reconciliation.tsv"
    contrasts = tmp_path / "contrasts.tsv"
    species_tree.write_text("((A_a:1,B_b:1):1,C_c:2);")
    gene_tree.write_text(
        "(((A_a_g1:1,B_b_g1:1):1,C_c_g1:2):1,((A_a_g2:1,B_b_g2:1):1,C_c_g2:2):1);"
    )
    trait.write_text(
        "leaf_name\texpression\n"
        "A_a_g1\t1\nB_b_g1\t2\nC_c_g1\t4\n"
        "A_a_g2\t2\nB_b_g2\t4\nC_c_g2\t8\n"
    )

    assert (
        main(
            [
                "reconcile",
                "--infile",
                str(gene_tree),
                "--species-tree",
                str(species_tree),
                "--outfile",
                str(reconciliation),
            ]
        )
        is None
    )
    assert (
        main(
            [
                "contrast",
                "--infile",
                str(gene_tree),
                "--trait",
                str(trait),
                "--columns",
                "expression",
                "--reconciliation",
                str(reconciliation),
                "--event-type",
                "speciation",
                "--outfile",
                str(contrasts),
            ]
        )
        is None
    )

    output = pd.read_csv(contrasts, sep="\t")
    assert len(output) == 4
    assert output["species_branch_id"].value_counts().to_dict() == {0: 2, 1: 2}

    duplication_contrasts = tmp_path / "duplication-contrasts.tsv"
    assert (
        main(
            [
                "contrast",
                "--infile",
                str(gene_tree),
                "--trait",
                str(trait),
                "--columns",
                "expression",
                "--reconciliation",
                str(reconciliation),
                "--event-type",
                "duplication",
                "--outfile",
                str(duplication_contrasts),
            ]
        )
        is None
    )
    duplication_output = pd.read_csv(duplication_contrasts, sep="\t")
    assert len(duplication_output) == 1
    assert duplication_output.iloc[0]["event_type"] == "duplication"


def test_reconcile_and_contrast_outputs_cannot_replace_inputs(tmp_path):
    species_tree = tmp_path / "species.nwk"
    gene_tree = tmp_path / "gene.nwk"
    trait = tmp_path / "trait.tsv"
    reconciliation = tmp_path / "reconciliation.tsv"
    species_tree.write_text("((A_a:1,B_b:1):1,C_c:2);")
    original_gene_tree = "((A_a_g1:1,B_b_g1:1):1,C_c_g1:2);"
    gene_tree.write_text(original_gene_tree)
    trait.write_text("leaf_name\texpression\nA_a_g1\t1\nB_b_g1\t2\nC_c_g1\t4\n")

    with pytest.raises(ValueError, match="must not overwrite input"):
        main(
            [
                "reconcile",
                "--infile",
                str(gene_tree),
                "--species-tree",
                str(species_tree),
                "--outfile",
                str(gene_tree),
            ]
        )
    assert gene_tree.read_text() == original_gene_tree

    main(
        [
            "reconcile",
            "--infile",
            str(gene_tree),
            "--species-tree",
            str(species_tree),
            "--outfile",
            str(reconciliation),
        ]
    )
    with pytest.raises(ValueError, match="must not overwrite input"):
        main(
            [
                "contrast",
                "--infile",
                str(gene_tree),
                "--trait",
                str(trait),
                "--columns",
                "expression",
                "--reconciliation",
                str(reconciliation),
                "--outfile",
                str(trait),
            ]
        )


def test_reconciliation_taxon_lists_round_trip_tip_names_with_commas(tmp_path):
    gene_tree = Tree("((A:1,B:1):1,C:2);", parser=1)
    species_tree = Tree("((X:1,Y:1):1,Z:2);", parser=1)
    gene_tree["A"].name = "gene,A"
    species_tree["X"].name = "species,X"
    species_by_gene = {"gene,A": "species,X", "B": "Y", "C": "Z"}
    reconciliation = build_reconciliation_table(
        gene_tree, species_tree, species_by_gene, event_source="lca"
    )
    path = tmp_path / "reconciliation.tsv"
    reconciliation.to_csv(path, sep="\t", index=False)

    parsed, _ = _read_reconciliation(path, CladeIndex(gene_tree))

    assert parsed is not None
    assert len(parsed) == 5
    root = max(parsed.values(), key=lambda record: int(record["num_taxa"]))
    assert '"gene,A"' in root["descendant_taxa"]


def test_contrast_rejects_nonpositive_branch_lengths():
    tree = Tree("((A:0,B:1):1,C:1);", parser=1)
    with pytest.raises(ValueError, match="positive finite branch lengths"):
        calculate_contrasts(tree, {"A": 1.0, "B": 2.0, "C": 3.0})


def test_unit_branch_lengths_are_an_explicit_fallback():
    tree = Tree("((A,B),C);", parser=1)
    contrasts = calculate_contrasts(
        tree,
        {"A": 1.0, "B": 2.0, "C": 3.0},
        branch_length="unit",
    )
    assert len(contrasts) == 2


def test_reconciliation_matches_reordered_gene_tree_by_stable_clade_id(tmp_path):
    species_tree = Tree("((A_a:1,B_b:1):1,C_c:2);", parser=1)
    original = Tree("((A_a_g1:1,B_b_g1:1):1,C_c_g1:2);", parser=1)
    reordered = Tree("(C_c_g1:2,(B_b_g1:1,A_a_g1:1):1);", parser=1)
    reconciliation = build_reconciliation_table(
        original,
        species_tree,
        {
            "A_a_g1": "A_a",
            "B_b_g1": "B_b",
            "C_c_g1": "C_c",
        },
        tree_id="OG1",
    )
    path = tmp_path / "reconciliation.tsv"
    reconciliation.to_csv(path, sep="\t", index=False)
    reconciliation_by_id, tree_id = _read_reconciliation(
        str(path), CladeIndex(reordered)
    )
    table = build_contrast_table(
        reordered,
        {"expression": {"A_a_g1": 1, "B_b_g1": 2, "C_c_g1": 4}},
        reconciliation_by_id=reconciliation_by_id,
        event_type="speciation",
        tree_id=tree_id,
    )
    assert len(table) == 2
    assert set(table["tree_id"]) == {"OG1"}
    assert set(table["gene_clade_id"]) == set(table["branch_clade_id"])
    species_contrasts = build_contrast_table(
        species_tree,
        {"trait": {"A_a": 1, "B_b": 2, "C_c": 4}},
    ).set_index("branch_clade_id")
    for _, row in table.iterrows():
        species_row = species_contrasts.loc[row["species_event_id"]]
        assert row["species_numerator_event_id"] == species_row["numerator_clade_id"]
        assert (
            row["species_denominator_event_id"] == species_row["denominator_clade_id"]
        )


def test_species_event_join_is_invariant_to_species_tree_child_order():
    species1 = Tree("((A:1,B:1):1,(C:1,D:1):1);", parser=1)
    species2 = Tree("((D:1,C:1):1,(B:1,A:1):1);", parser=1)
    values = {"A": 10, "B": 20, "C": 300, "D": 400}
    table1 = build_contrast_table(species1, {"trait": values})
    table2 = build_contrast_table(species2, {"trait": values})
    observed1 = table1.set_index("branch_clade_id")["raw_contrast"].sort_index()
    observed2 = table2.set_index("branch_clade_id")["raw_contrast"].sort_index()
    pd.testing.assert_series_equal(observed1, observed2)


def test_complete_coverage_is_safe_default_and_partial_is_explicit():
    species_tree = Tree("((A_a:1,B_b:1):1,C_c:2);", parser=1)
    gene_tree = Tree("(A_a_g1:1,C_c_g1:2);", parser=1)
    reconciliation = build_reconciliation_table(
        gene_tree,
        species_tree,
        {"A_a_g1": "A_a", "C_c_g1": "C_c"},
    )
    by_id = {
        record["gene_clade_id"]: record for record in reconciliation.to_dict("records")
    }
    values = {"A_a_g1": 1, "C_c_g1": 4}
    complete = build_contrast_table(
        gene_tree,
        {"expression": values},
        reconciliation_by_id=by_id,
        event_type="speciation",
    )
    partial = build_contrast_table(
        gene_tree,
        {"expression": values},
        reconciliation_by_id=by_id,
        event_type="speciation",
        speciation_coverage="any",
    )
    assert complete.empty
    assert len(partial) == 1
    assert partial.iloc[0]["coverage_status"] == "partial"


def test_extreme_finite_inputs_never_emit_nonfinite_outputs():
    tree = Tree("(A:1,B:1);", parser=1)
    with pytest.raises(ValueError, match="non-finite raw contrast"):
        calculate_contrasts(tree, {"A": 1e308, "B": -1e308})
    with pytest.raises(ValueError, match="non-finite standardized contrast"):
        calculate_contrasts(
            Tree("(A:1e-300,B:1e-300);", parser=1),
            {"A": 1e308, "B": 0},
        )


def test_stable_variance_formula_avoids_overflow_and_underflow():
    large = Tree("((A:1e200,B:1e200):1e200,C:1e200);", parser=1)
    assert len(calculate_contrasts(large, {"A": 1, "B": 2, "C": 3})) == 2
    tiny = Tree("((A:1e-300,B:1e-300):1e-300,C:1e-300);", parser=1)
    result = calculate_contrasts(tiny, {"A": 1, "B": 2, "C": 3})
    assert result[tiny]["contrast_variance"] == pytest.approx(
        2.5e-300, rel=1e-12, abs=0
    )


def test_reconciliation_enum_values_are_validated(tmp_path):
    tree = Tree("(A:1,B:1);", parser=1)
    species = Tree("(A:1,B:1);", parser=1)
    reconciliation = build_reconciliation_table(tree, species, {"A": "A", "B": "B"})
    reconciliation.loc[0, "eligible"] = "YES"
    path = tmp_path / "bad.tsv"
    reconciliation.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="invalid eligible"):
        _read_reconciliation(str(path), CladeIndex(tree))


def test_reconciliation_branch_and_lineage_identifiers_are_validated(tmp_path):
    tree = Tree("((A:1,B:1):1,C:2);", parser=1)
    reconciliation = build_reconciliation_table(
        tree, Tree("((A:1,B:1):1,C:2);", parser=1), {"A": "A", "B": "B", "C": "C"}
    )
    reconciliation.loc[1, "gene_branch_id"] = reconciliation.loc[0, "gene_branch_id"]
    path = tmp_path / "duplicate-branch.tsv"
    reconciliation.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="duplicated gene_branch_id"):
        _read_reconciliation(str(path), CladeIndex(tree))

    reconciliation = build_reconciliation_table(
        tree, Tree("((A:1,B:1):1,C:2);", parser=1), {"A": "A", "B": "B", "C": "C"}
    )
    root_index = reconciliation.index[reconciliation["node_class"] == "root"][0]
    leaf = reconciliation[reconciliation["node_class"] == "leaf"].iloc[0]
    reconciliation.loc[root_index, "lineage_clade_id"] = leaf["gene_clade_id"]
    reconciliation.loc[root_index, "lineage_id"] = leaf["gene_branch_id"]
    path = tmp_path / "bad-lineage.tsv"
    reconciliation.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="ancestor clade"):
        _read_reconciliation(str(path), CladeIndex(tree))


def test_reconciliation_orientation_and_coverage_are_validated(tmp_path):
    tree = Tree("((A:1,B:1):1,C:2);", parser=1)
    species = Tree("((A:1,B:1):1,C:2);", parser=1)
    reconciliation = build_reconciliation_table(
        tree, species, {"A": "A", "B": "B", "C": "C"}
    )
    eligible_index = reconciliation.index[reconciliation["eligible"] == "yes"][0]
    reconciliation.loc[eligible_index, "contrast_denominator_gene_clade_id"] = (
        reconciliation.loc[eligible_index, "contrast_numerator_gene_clade_id"]
    )
    path = tmp_path / "bad-orientation.tsv"
    reconciliation.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="immediate children"):
        _read_reconciliation(str(path), CladeIndex(tree))

    reconciliation = build_reconciliation_table(
        tree, species, {"A": "A", "B": "B", "C": "C"}
    )
    eligible_index = reconciliation.index[reconciliation["eligible"] == "yes"][0]
    reconciliation.loc[eligible_index, "numerator_species_coverage"] = 0.5
    path = tmp_path / "bad-coverage.tsv"
    reconciliation.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError, match="must match observed/full counts"):
        _read_reconciliation(str(path), CladeIndex(tree))
