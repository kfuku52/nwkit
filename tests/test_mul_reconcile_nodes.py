import csv
import json
import subprocess
import sys

import pytest

from nwkit.clade_index import CladeIndex
from nwkit.mul_reconcile_nodes import rooted_topology_id
from nwkit.util import read_tree


def tree(text):
    return read_tree(text, "auto", True, quiet=True)


def run(tmp_path, *extra, gene="(x1_X,x2_X);", species="((A,X),B);"):
    (tmp_path / "gene.nwk").write_text(gene)
    (tmp_path / "species.nwk").write_text(species)
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "nwkit",
            "mul-reconcile",
            "--infile",
            str(tmp_path / "gene.nwk"),
            "--species-tree",
            str(tmp_path / "species.nwk"),
            "--species-regex",
            r".*_([^_]+)$",
            "--h1",
            "X",
            "--h2",
            "A B X",
            "--outfile",
            str(tmp_path / "scores.tsv"),
            "--node-out",
            str(tmp_path / "nodes.tsv"),
            "--model-out",
            str(tmp_path / "model.json"),
            *extra,
        ],
        text=True,
        capture_output=True,
    )


def test_rooted_identity_is_order_and_length_independent_but_topology_sensitive():
    first = tree("((a:1,b:2):3,c:4);")
    assert rooted_topology_id(first) == rooted_topology_id(tree("(c:99,(b:7,a:8):1);"))
    assert rooted_topology_id(first) != rooted_topology_id(tree("((a,c),b);"))
    with pytest.raises(ValueError, match="rooted"):
        rooted_topology_id(tree("[&U]((a,b),c);"))


def test_all_tied_best_candidates_and_complete_assignments_are_serial_parallel_identical(
    tmp_path,
):
    result = run(tmp_path, "--report", str(tmp_path / "legacy.tsv"))
    assert result.returncode == 0, result.stderr
    model = json.loads((tmp_path / "model.json").read_text())
    assert model["best_hypotheses"] == [1, 2, 3]
    with (tmp_path / "nodes.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {int(r["mul.tree"]) for r in rows} == set(model["best_hypotheses"])
    index = CladeIndex(tree("(x1_X,x2_X);"))
    ids = {index.clade_id_for_node(n) for n in index.mask_by_node}
    for candidate in model["best_hypotheses"]:
        subset = [r for r in rows if int(r["mul.tree"]) == candidate]
        count = int(subset[0]["optimal.mappings"])
        assert count >= 2
        assert len(subset) == count * 3
        for mapping in range(1, count + 1):
            group = [r for r in subset if int(r["mapping.id"]) == mapping]
            assert {r["gene_clade_id"] for r in group} == ids
            assert (
                sum(
                    int(r[k])
                    for r in group
                    for k in ("duplication", "child_edge_losses", "root_losses")
                )
                == 2
            )
    with (tmp_path / "legacy.tsv").open(newline="") as handle:
        legacy = list(csv.DictReader(handle, delimiter="\t"))
    assert {int(r["mul.tree"]) for r in legacy} == {model["reported_hypothesis"]}
    before = (tmp_path / "nodes.tsv").read_bytes()
    result = run(tmp_path, "--cpus", "2")
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "nodes.tsv").read_bytes() == before


def test_mapping_limit_preserves_the_whole_previous_bundle(tmp_path):
    targets = [tmp_path / name for name in ("scores.tsv", "nodes.tsv", "model.json")]
    for target in targets:
        target.write_text("previous " + target.name)
    before = {p: p.read_bytes() for p in targets}
    result = run(tmp_path, "--max-maps", "1")
    assert result.returncode != 0 and "no maps were truncated" in result.stderr
    assert {p: p.read_bytes() for p in targets} == before


def test_multiple_genes_use_global_best_candidates_and_per_gene_scores(tmp_path):
    result = run(tmp_path, gene="(x1_X,x2_X);\n(a_A,x_X);")
    assert result.returncode == 0, result.stderr
    model = json.loads((tmp_path / "model.json").read_text())
    assert model["num_gene_trees"] == 2 and model["best_hypotheses"] == [2]
    with (tmp_path / "nodes.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {r["mul.tree"] for r in rows} == {"2"}
    assert {r["gene.tree"] for r in rows} == {"1", "2"}
    for number, gene, score in (("1", "(x1_X,x2_X);", 2), ("2", "(a_A,x_X);", 1)):
        subset = [r for r in rows if r["gene.tree"] == number]
        assert {r["gene_topology_id"] for r in subset} == {
            rooted_topology_id(tree(gene))
        }
        assert {int(r["total.score"]) for r in subset} == {score}
        count = int(subset[0]["optimal.mappings"])
        assert len(subset) == 3 * count
        for mapping in range(1, count + 1):
            group = [r for r in subset if int(r["mapping.id"]) == mapping]
            assert (
                sum(
                    int(r[k])
                    for r in group
                    for k in ("duplication", "child_edge_losses", "root_losses")
                )
                == score
            )


@pytest.mark.parametrize("score_model", ["msc", "locus-mc"])
def test_node_output_is_not_silently_ignored_for_other_models(tmp_path, score_model):
    result = run(tmp_path, "--score-model", score_model)
    assert (
        result.returncode != 0
        and "--node-out requires --score-model dl" in result.stderr
    )
    assert not (tmp_path / "nodes.tsv").exists()


@pytest.mark.parametrize("target", ["gene.nwk", "scores.tsv", "-"])
def test_node_target_protects_inputs_and_other_outputs(tmp_path, target):
    result = run(
        tmp_path, "--node-out", str(tmp_path / target) if target != "-" else "-"
    )
    assert result.returncode != 0
    assert (tmp_path / "gene.nwk").read_text() == "(x1_X,x2_X);"


def test_duplicate_supplied_mul_occurrences_are_not_given_spurious_clade_identity(
    tmp_path,
):
    result = run(
        tmp_path, "--multree", "yes", "--h1", "", "--h2", "", species="((A,X),(B,X));"
    )
    assert result.returncode != 0
    # A supplied MUL input without selectors reaches the occurrence guard.
    command = [
        sys.executable,
        "-m",
        "nwkit",
        "mul-reconcile",
        "--infile",
        str(tmp_path / "gene.nwk"),
        "--species-tree",
        str(tmp_path / "species.nwk"),
        "--multree",
        "yes",
        "--species-regex",
        r".*_([^_]+)$",
        "--node-out",
        str(tmp_path / "nodes.tsv"),
    ]
    result = subprocess.run(command, text=True, capture_output=True)
    assert result.returncode != 0 and "Duplicated leaf labels" in result.stderr
