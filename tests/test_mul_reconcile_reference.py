"""Optional executed comparisons against the original author's program/data."""

import itertools
import os
import random
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from nwkit.mul_reconcile import run_search, topology_text
from nwkit.mul_reconcile_model import Reconciliation, hypotheses
from tests.test_mul_reconcile import exhaustive, parser, tree

REFERENCE = os.environ.get("NWKIT_GRAMPA_REFERENCE")
pytestmark = pytest.mark.skipif(
    not REFERENCE,
    reason="Set NWKIT_GRAMPA_REFERENCE to the author's grampa.py checkout.",
)


def compare(tmp_path, species, gene_text, h1=None, expected_differences=None):
    genes = [tree(line) for line in gene_text.splitlines() if line.strip()]
    candidates = hypotheses(tree(species), h1)
    native = {candidate.id: candidate for candidate in candidates}
    native_scores = {
        (native[number].h1, native[number].h2): score
        for number, score, _ in run_search(candidates, genes, parser())
    }
    input_file = tmp_path / "genes.nwk"
    input_file.write_text(gene_text)
    args = [
        sys.executable,
        str(REFERENCE),
        "-s",
        species,
        "-g",
        str(input_file),
        "-o",
        str(tmp_path / "reference"),
        "-p",
        "1",
        "-v",
        "0",
        "--maps",
    ]
    if h1:
        args += ["-h1", h1]
    result = subprocess.run(args, capture_output=True, text=True, timeout=180)
    assert result.returncode == 0, result.stdout + result.stderr
    scores = pd.read_csv(
        tmp_path / "reference/grampa-scores.txt",
        sep="\t",
        comment="#",
        keep_default_na=False,
    )
    reference_scores = {
        (str(row["h1.node"]), str(row["h2.node"])): int(row["score"])
        for _, row in scores.iterrows()
    }
    differences = {
        pair: (reference_scores[pair], native_scores[pair])
        for pair in reference_scores
        if reference_scores[pair] != native_scores[pair]
    }
    assert reference_scores.keys() == native_scores.keys()
    assert differences == (expected_differences or {})
    checks = pd.read_csv(tmp_path / "reference/grampa-checknums.txt", sep="\t")
    assert not checks["over.cap.filtered"].astype(str).str.contains("Y").any()
    return len(genes), len(native_scores)


def test_original_program_small_null_auto_allo_missing_and_duplicate_cases(tmp_path):
    compare(
        tmp_path,
        "((A,X),B);",
        "((a_A,x1_X),(b_B,x2_X));\n((a_A,x_X),b_B);\n"
        "((x1_X,x2_X),(a_A,b_B));\n(a_A,b_B);\n(x1_X,x2_X);\n",
        expected_differences={("<2>", "<2>"): (15, 13)},
    )


def test_grouped_root_auto_difference_against_original_lca_and_full_enumeration(
    monkeypatch,
):
    """The reference grouping can exclude the optimal split of distinct species."""
    monkeypatch.syspath_prepend(str(Path(str(REFERENCE)).parent))
    from grampa_lib.mul_recon import reconLCA
    from grampa_lib.recontree import treeParse

    candidate = hypotheses(tree("((A,X),B);"), "2", "2")[1]
    species_text = topology_text(candidate.tree, internal_names=False).replace("+", "")
    species_info, _ = treeParse(species_text)
    observed = []
    gene_texts = [
        "((a_A,x1_X),(b_B,x2_X));",
        "((a_A,x_X),b_B);",
        "((x1_X,x2_X),(a_A,b_B));",
        "(a_A,b_B);",
        "(x1_X,x2_X);",
    ]
    for text in gene_texts:
        gene = tree(text)
        gene_info, _ = treeParse(text)
        tips = [label for label, info in gene_info.items() if info[2] == "tip"]
        scores = []
        for copies in itertools.product(("", "*"), repeat=len(tips)):
            assignment = {
                label: label.rsplit("_", 1)[1] + copy
                for label, copy in zip(tips, copies, strict=True)
            }
            maps = {
                label: [assignment[label]] if label in assignment else []
                for label in gene_info
            }
            scores.append(reconLCA(gene_info, species_info, maps))
        native = Reconciliation(gene, candidate.tree, parser())
        assert native.score == min(scores) == sum(exhaustive(gene, candidate.tree)[0])
        observed.append(native.score)
    assert sum(observed) == 13
    print(
        f"Whole-root auto: exact original-LCA enumeration={observed}, sum=13; grouped reference=15."
    )


def test_author_manual_dataset_all_candidate_scores(tmp_path):
    root = Path(str(REFERENCE)).parent
    species = (root / "data/manual_species_tree.tre").read_text().strip().rstrip(
        ";"
    ) + ";"
    # The author supplies one tree per line without terminal semicolons.
    genes = (
        "\n".join(
            line.strip().rstrip(";") + ";"
            for line in (root / "data/manual_gene_trees.txt").read_text().splitlines()
            if line.strip()
        )
        + "\n"
    )
    counts = compare(tmp_path, species, genes, "x,y,z")
    assert counts[0] == 25
    print(
        f"Author manual dataset: {counts[0]} trees, {counts[1]} hypotheses; every score matched."
    )


def test_seeded_every_candidate_against_original_lca_assignment_oracle(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(str(REFERENCE)).parent))
    from grampa_lib.mul_recon import reconLCA
    from grampa_lib.recontree import treeParse

    rng = random.Random(20261004)
    candidates = hypotheses(tree("((A,X),B);"))
    checked = 0
    for _ in range(12):
        subtrees = [
            f"g{i}_{rng.choice(('A', 'X', 'B'))}" for i in range(rng.randint(2, 6))
        ]
        while len(subtrees) > 1:
            first = subtrees.pop(rng.randrange(len(subtrees)))
            second = subtrees.pop(rng.randrange(len(subtrees)))
            subtrees.append(f"({first},{second})")
        text = subtrees[0] + ";"
        gene = tree(text)
        gene_info, _ = treeParse(text)
        tips = [label for label, info in gene_info.items() if info[2] == "tip"]
        for candidate in candidates:
            species_info, _ = treeParse(
                topology_text(candidate.tree, internal_names=False).replace("+", "")
            )
            species_tips = [
                label for label, info in species_info.items() if info[2] == "tip"
            ]
            choices = [
                [
                    label
                    for label in species_tips
                    if label.rstrip("*") == tip.rsplit("_", 1)[1]
                ]
                for tip in tips
            ]
            scores = []
            for assignment in itertools.product(*choices):
                selected = dict(zip(tips, assignment, strict=True))
                maps = {
                    label: [selected[label]] if label in selected else []
                    for label in gene_info
                }
                scores.append(reconLCA(gene_info, species_info, maps))
            native = Reconciliation(gene, candidate.tree, parser())
            assert native.score == min(scores)
            assert native.num_maps == scores.count(min(scores))
            checked += 1
    assert checked == 240
    print(f"Original-LCA exhaustive oracle: {checked} seeded gene/hypothesis pairs.")
