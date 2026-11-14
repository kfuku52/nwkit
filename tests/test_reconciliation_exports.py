"""Reusable optimal roots and per-event losses for workflow consumers."""

import subprocess
import sys

import pandas as pd
import pytest
from ete4 import Tree

from nwkit.clade_mapping import projected_root_split
from nwkit.reconcile import build_reconciliation_table
from nwkit.root import reconciliation_rooting
from nwkit.util import read_tree, read_tree_strings


def cli(*args):
    return subprocess.run(
        [sys.executable, "-m", "nwkit", *map(str, args)], text=True, capture_output=True
    )


def test_losses_include_missing_lineages_but_not_above_gene_root():
    species = Tree("(A_a:1,(B_b:1,C_c:1):1);", parser=1)
    gene = Tree("((A_a_1:1,B_b_1:1):1,A_a_2:2);", parser=1)
    mapping = {name: "_".join(name.split("_")[:2]) for name in gene.leaf_names()}
    table = build_reconciliation_table(gene, species, mapping)
    assert table.implied_losses.sum() == 2
    assert table.loc[table.event_type.eq("duplication"), "implied_losses"].tolist() == [
        1
    ]
    assert table.loc[table.event_type.eq("speciation"), "implied_losses"].tolist() == [
        1
    ]
    assert table.loc[table.event_type.eq("leaf"), "implied_losses"].eq(0).all()
    # A family confined to B has no inferred ancestral losses outside its root.
    restricted = Tree("(B_b_1:1,B_b_2:1);", parser=1)
    table = build_reconciliation_table(
        restricted, species, {name: "B_b" for name in restricted.leaf_names()}
    )
    assert table.implied_losses.sum() == 0


def test_non_lca_and_unmapped_losses_are_undefined():
    species = Tree("(A_a:1,B_b:1);", parser=1)
    gene = Tree("(A_a_1:1,B_b_1:1);", parser=1)
    table = build_reconciliation_table(
        gene, species, {"A_a_1": "A_a", "B_b_1": "B_b"}, event_source="species-overlap"
    )
    assert table.implied_losses.isna().all()
    table = build_reconciliation_table(
        gene, species, {"A_a_1": "A_a", "B_b_1": "absent"}
    )
    assert pd.isna(table.iloc[0].implied_losses)


def test_all_optimal_roots_are_exported_without_distance_changes(tmp_path):
    source = tmp_path / "gene.nwk"
    source.write_text("((A_a_1:1,B_b_1:2):3,(A_a_2:4,B_b_2:5):6);")
    species_path = tmp_path / "species.nwk"
    species_path.write_text("(A_a:1,B_b:1);")
    selected, collection = tmp_path / "selected.nwk", tmp_path / "roots.nwk"
    result = cli(
        "root",
        "--method",
        "reconciliation",
        "--infile",
        source,
        "--species-tree",
        species_path,
        "--outfile",
        selected,
        "--candidates-out",
        collection,
        "--duplication-cost",
        "1.5",
    )
    assert result.returncode == 0, result.stderr
    original = read_tree(str(source), "auto", True)
    species = read_tree(str(species_path), "auto", True)
    names = frozenset(original.leaf_names())
    _, evaluation = reconciliation_rooting(
        original.copy(),
        species,
        {name: "_".join(name.split("_")[:2]) for name in names},
        duplication_cost=1.5,
        _return_evaluation=True,
    )
    trees = []
    for index, text in enumerate(read_tree_strings(str(collection))):
        path = tmp_path / f"candidate{index}.nwk"
        path.write_text(text)
        trees.append(read_tree(str(path), "auto", True))
    assert len(trees) == len(evaluation.candidates)
    assert {projected_root_split(tree, names) for tree in trees} == {
        c.split for c in evaluation.candidates
    }
    for tree in trees:
        for left in names:
            for right in names:
                assert tree.get_distance(left, right) == pytest.approx(
                    original.get_distance(left, right)
                )
    assert projected_root_split(
        read_tree(str(selected), "auto", True), names
    ) == projected_root_split(trees[0], names)


@pytest.mark.parametrize(
    "case", ["same_output", "input_overlap", "wrong_method", "bad_species"]
)
def test_candidate_export_rejects_invalid_request_before_replacing_files(
    tmp_path, case
):
    source, species = tmp_path / "gene.nwk", tmp_path / "species.nwk"
    selected, collection = tmp_path / "selected.nwk", tmp_path / "roots.nwk"
    source.write_text("(A_a_1:1,B_b_1:1);")
    species.write_text("(A_a:1,B_b:1);" if case != "bad_species" else "(A_a:1,C_c:1);")
    selected.write_text("old selected")
    collection.write_text("old collection")
    before = {
        path: path.read_bytes() for path in (source, species, selected, collection)
    }
    target = (
        selected
        if case == "same_output"
        else source
        if case == "input_overlap"
        else collection
    )
    result = cli(
        "root",
        "--method",
        "midpoint" if case == "wrong_method" else "reconciliation",
        "--infile",
        source,
        "--species-tree",
        species,
        "--outfile",
        selected,
        "--candidates-out",
        target,
    )
    assert result.returncode != 0
    assert {path: path.read_bytes() for path in before} == before


@pytest.mark.parametrize("target", ["", "-"])
def test_candidate_export_rejects_non_file_target(tmp_path, target):
    source, species, selected = (
        tmp_path / name for name in ("gene", "species", "selected")
    )
    source.write_text("(A_a_1:1,B_b_1:1);")
    species.write_text("(A_a:1,B_b:1);")
    selected.write_text("previous")
    result = cli(
        "root",
        "--method",
        "reconciliation",
        "--infile",
        source,
        "--species-tree",
        species,
        "--outfile",
        selected,
        "--candidates-out",
        target,
    )
    assert result.returncode != 0
    assert selected.read_text() == "previous"


@pytest.mark.parametrize("failure_call", [2, 6])
def test_candidate_write_failure_preserves_both_outputs(
    tmp_path, monkeypatch, failure_call
):
    from types import SimpleNamespace

    import nwkit.root as root_module

    gene = Tree("((A_a_1:1,A_a_2:2):3,(A_a_3:4,A_a_4:5):6);", parser=1)
    species = Tree("(A_a:1,B_b:1);", parser=1)
    rooted, evaluation = reconciliation_rooting(
        gene,
        species,
        {name: "A_a" for name in gene.leaf_names()},
        _return_evaluation=True,
    )
    selected, collection = tmp_path / "selected", tmp_path / "collection"
    selected.write_text("old selected")
    collection.write_text("old collection")
    original_write = root_module.write_tree
    calls = 0

    def fail_second_write(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == failure_call:
            raise OSError("injected disk failure")
        return original_write(*args, **kwargs)

    monkeypatch.setattr(root_module, "write_tree", fail_second_write)
    args = SimpleNamespace(
        outfile=str(selected), candidates_out=str(collection), outformat=1
    )
    with pytest.raises(OSError, match="injected disk failure"):
        root_module._write_reconciliation_candidates(rooted, evaluation, args, set())
    assert selected.read_text() == "old selected"
    assert collection.read_text() == "old collection"


@pytest.mark.parametrize("weights", [(1.5, 1), (0, 1), (1, 0), (0.1, 0.3)])
def test_optimal_root_scores_match_independent_lca_counts(weights):
    from nwkit.root import _root_by_outgroup_set

    species = Tree("((A_a:1,B_b:1):1,(C_c:1,D_d:1):1);", parser=1)
    gene = Tree("(((A_a_1:1,C_c_1:2):1,B_b_1:3):1,(A_a_2:4,D_d_1:5):1);", parser=1)
    mapping = {name: "_".join(name.split("_")[:2]) for name in gene.leaf_names()}
    rooted, evaluation = reconciliation_rooting(
        gene,
        species,
        mapping,
        duplication_cost=weights[0],
        loss_cost=weights[1],
        _return_evaluation=True,
    )
    # Enumerate every physical edge independently of the rooting optimizer.
    names = frozenset(gene.leaf_names())
    scored = {}
    for node in gene.traverse():
        if node.is_root:
            continue
        candidate = _root_by_outgroup_set(gene, set(node.leaf_names()))
        table = build_reconciliation_table(candidate, species, mapping)
        duplications = int(table.event_type.eq("duplication").sum())
        losses = int(table.implied_losses.sum())
        scored[projected_root_split(candidate, names)] = (duplications, losses)
    scores = {
        split: dup * weights[0] + loss * weights[1]
        for split, (dup, loss) in scored.items()
    }
    best = min(scores.values())
    assert {candidate.split for candidate in evaluation.candidates} == {
        split for split, score in scores.items() if score == pytest.approx(best)
    }
    for candidate in evaluation.candidates:
        assert (
            candidate.metrics["duplications"],
            candidate.metrics["losses"],
        ) == scored[candidate.split]
