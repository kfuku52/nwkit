import itertools
import json
import random
import subprocess
import sys

import pandas as pd
import pytest

from nwkit.mul_reconcile import annotated_mapping_text, run_search
from nwkit.mul_reconcile_model import Reconciliation, hypotheses
from nwkit.species_parser import get_species_parser
from nwkit.util import read_tree


def tree(text):
    return read_tree(text, "auto", True, quiet=True)


def parser():
    return get_species_parser(species_parser="legacy", species_regex=r".*_([^_]+)$")


def exhaustive(gene, species, *, with_mappings=False):
    """Independent enumeration and ancestor-path scoring, with no DP/LCA index."""
    leaves = list(gene.leaves())
    choices = [
        [
            node
            for node in species.leaves()
            if node.props.get("mul_species", node.name) == leaf.name.rsplit("_", 1)[1]
        ]
        for leaf in leaves
    ]
    scored = []
    for assignment in itertools.product(*choices):
        mapped = dict(zip(leaves, assignment, strict=True))
        duplication = loss = 0
        for node in gene.traverse("postorder"):
            if node.is_leaf:
                continue
            first, second = (mapped[child] for child in node.children)
            second_path = [second, *second.ancestors()]
            parent = next(
                ancestor
                for ancestor in [first, *first.ancestors()]
                if ancestor in second_path
            )
            mapped[node] = parent
            dup = int(parent is first or parent is second)
            duplication += dup
            for child in (first, second):
                steps = 0
                while child is not parent:
                    child = child.up
                    steps += 1
                loss += steps - 1 + dup
        loss += len(list(mapped[gene].ancestors()))
        species_index = {
            node: index for index, node in enumerate(species.traverse("preorder"))
        }
        signature = tuple(
            species_index[mapped[node]] for node in gene.traverse("postorder")
        )
        scored.append((duplication + loss, duplication, loss, signature))
    best = min(row[0] for row in scored)
    return sorted(
        (dup, loss, signature) if with_mappings else (dup, loss)
        for score, dup, loss, signature in scored
        if score == best
    )


@pytest.mark.parametrize(
    "gene",
    [
        "(a_A,b_B);",
        "((a_A,x1_X),(b_B,x2_X));",
        "(((x1_X,x2_X),a_A),(x3_X,b_B));",
        "(x1_X,x2_X);",
        "(((a_A,b_B),x1_X),x2_X);",
    ],
)
def test_dynamic_program_matches_all_assignments_for_every_hypothesis(gene):
    genes = tree(gene)
    for candidate in hypotheses(tree("((A,X),B);")):
        result = Reconciliation(genes, candidate.tree, parser())
        expected = exhaustive(genes, candidate.tree)
        observed = sorted((dup, loss) for dup, loss, _ in result.mappings())
        assert observed == expected
        assert result.num_maps == len(expected)
        assert result.score == sum(expected[0])


def test_seeded_random_topologies_and_copy_counts_against_full_enumeration():
    rng = random.Random(20261002)
    candidates = hypotheses(tree("((A,X),(B,Y));"))
    checked = 0
    for _ in range(24):
        subtrees = [
            f"g{number}_{rng.choice(('A', 'X', 'B', 'Y'))}"
            for number in range(rng.randint(2, 7))
        ]
        while len(subtrees) > 1:
            first = subtrees.pop(rng.randrange(len(subtrees)))
            second = subtrees.pop(rng.randrange(len(subtrees)))
            subtrees.append(f"({first},{second})")
        gene = tree(subtrees[0] + ";")
        for candidate in candidates:
            result = Reconciliation(gene, candidate.tree, parser())
            expected = exhaustive(gene, candidate.tree)
            assert sorted((dup, loss) for dup, loss, _ in result.mappings()) == expected
            assert result.score == sum(expected[0])
            assert result.num_maps == len(expected)
            checked += 1
    assert checked == 24 * len(candidates)
    print(f"Seeded independent enumeration: {checked} gene/hypothesis combinations.")


@pytest.mark.parametrize(
    "species,gene,h1,h2",
    [
        ("((A,X),B);", "((a_A,x1_X),(b_B,x2_X));", "X", "B"),
        ("(X,(A,B));", "((x1_X,x2_X),(a_A,b_B));", "X", "X"),
    ],
)
def test_exact_allo_and_auto_positive_controls(species, gene, h1, h2):
    candidates = hypotheses(tree(species), h1, h2)
    scores = run_search(candidates, [tree(gene)], parser())
    assert scores[0][1] > 0
    assert scores[1][1] == 0
    assert candidates[1].kind == ("autopolyploid" if h1 == h2 else "allopolyploid")


def test_no_polyploid_null_and_root_loss_contract():
    species = tree("((A,X),B);")
    candidates = hypotheses(species)
    scores = run_search(candidates, [tree("((a_A,x_X),b_B);")], parser())
    assert scores[0][1] == 0
    assert all(score > 0 for _, score, _ in scores[1:])
    result = Reconciliation(tree("(x1_X,x2_X);"), candidates[0].tree, parser())
    assert result.score == 3  # One duplication, two unobserved ancestral lineages.
    for dup, loss, rows in result.mappings():
        assert dup == sum(row["duplication"] for row in rows)
        assert loss == sum(
            row["child_edge_losses"] + row["root_losses"] for row in rows
        )


def test_selectors_numbering_and_candidate_placement():
    candidates = hypotheses(tree("((A,X),B);"), "1", "B")
    assert candidates[1].h1 == "<1>"
    assert sorted(candidates[1].tree.leaf_names()) == ["A*", "A+", "B", "X*", "X+"]
    assert len(hypotheses(tree("((A,X),B);"), "A,X", "B")) == 2
    with pytest.raises(ValueError, match="admissible"):
        hypotheses(tree("((A,X),B);"), "A,X", "X")
    with pytest.raises(ValueError, match="non-monophyletic"):
        hypotheses(tree("((A,X),B);"), "A,B")
    with pytest.raises(ValueError, match="max-candidates"):
        hypotheses(tree("((A,X),B);"), max_candidates=1)


@pytest.mark.parametrize("label", ["1", "<1>", "A+", "X*"])
def test_reserved_species_labels_are_rejected_before_node_lookup(label):
    with pytest.raises(ValueError, match="Species labels"):
        hypotheses(tree(f"(({label},X),B);"))


def test_supplied_multree_and_tie_mapping_counts():
    candidate = hypotheses(tree("((A,X),(B,X));"), multree=True)[0]
    result = Reconciliation(tree("(x1_X,x2_X);"), candidate.tree, parser())
    assert sorted((dup, loss) for dup, loss, _ in result.mappings()) == exhaustive(
        result.gene, candidate.tree
    )
    assert result.num_maps >= 2
    with pytest.raises(ValueError, match="max-maps"):
        list(result.mappings(max_maps=1))


def test_supplied_multree_all_optimal_node_assignments_against_enumeration():
    rng = random.Random(20261003)
    candidate = hypotheses(tree("(((A,X),(B,X)),(A,X));"), multree=True)[0]
    for _ in range(30):
        subtrees = [
            f"g{i}_{rng.choice(('A', 'X', 'B'))}" for i in range(rng.randint(2, 7))
        ]
        while len(subtrees) > 1:
            first = subtrees.pop(rng.randrange(len(subtrees)))
            second = subtrees.pop(rng.randrange(len(subtrees)))
            subtrees.append(f"({first},{second})")
        gene = tree(subtrees[0] + ";")
        result = Reconciliation(gene, candidate.tree, parser())
        observed = sorted(
            (dup, loss, tuple(row["mul_node"] for row in rows))
            for dup, loss, rows in result.mappings()
        )
        assert observed == exhaustive(gene, candidate.tree, with_mappings=True)
        assert len({signature for _, _, signature in observed}) == result.num_maps


def test_deep_gene_tree_mapping_and_annotation_without_recursion():
    from ete4 import Tree

    gene = Tree()
    cursor = gene
    for index in range(1200):
        cursor.add_child(name=f"g{index}_A")
        cursor = cursor.add_child()
    cursor.name = "g1200_A"
    result = Reconciliation(gene, hypotheses(tree("(A,B);"))[0].tree, parser())
    dup, loss, rows = next(result.mappings())
    assert (dup, loss, result.score, result.num_maps) == (1200, 1, 1201, 1)
    assert len(rows) == 2401
    assert sum(row["child_edge_losses"] + row["root_losses"] for row in rows) == loss
    assert annotated_mapping_text(gene, rows).count("[A+-0]") == 1201


def test_child_order_and_branch_lengths_do_not_change_parsimony():
    gene = tree("((a_A,x1_X),(b_B,x2_X));")
    species = tree("((A,X),B);")
    before = [row[1] for row in run_search(hypotheses(species), [gene], parser())]
    for node in gene.traverse():
        node.dist = 987.65
        node.children.reverse()
    after = [row[1] for row in run_search(hypotheses(species), [gene], parser())]
    assert before == after


def test_unused_check_rows_are_not_collected():
    rows = run_search(
        hypotheses(tree("(A,B);")), [tree("(a_A,b_B);")], parser(), collect_checks=False
    )
    assert all(checks == [] for _, _, checks in rows)


def test_legacy_annotation_escapes_comment_metacharacters_but_json_retains_names():
    gene = tree("('a_A]x','b_B[x');")
    candidate = hypotheses(tree("('A]x','B[x');"))[0]
    result = Reconciliation(gene, candidate.tree, parser())
    _, _, mapping = next(result.mappings())
    annotated = annotated_mapping_text(gene, mapping)
    assert "[A%5Dx+-0]" in annotated
    assert "[B%5Bx+-0]" in annotated
    assert {row["mul_label"] for row in mapping} == {"A]x", "B[x", "<1>"}


@pytest.mark.parametrize(
    "species,gene,message",
    [
        ("[&U]((A,X),B);", "(a_A,b_B);", "rooted"),
        ("((A,X),B);", "[&U](a_A,b_B);", "rooted"),
        ("((A,X),B);", "(a_A,b_B,x_X);", "rooted|binary"),
        ("((A,X),B);", "(a_A,b_C);", "match a species"),
        ("((A,X),B);", "(a_A,a_A);", "unique|[Dd]uplicat"),
    ],
)
def test_invalid_inputs_fail(species, gene, message):
    with pytest.raises(ValueError, match=message):
        candidates = hypotheses(tree(species))
        run_search(candidates, [tree(gene)], parser())


def test_state_limit_is_failure_not_candidate_dependent_filter():
    candidate = hypotheses(tree("((A,X),B);"), "X", "B")[1]
    with pytest.raises(ValueError, match="max-state-pairs"):
        Reconciliation(tree("((a_A,x1_X),(b_B,x2_X));"), candidate.tree, parser(), 1)


def invoke(tmp_path, extra=(), text=None):
    species = tmp_path / "species.nwk"
    genes = tmp_path / "genes.nwk"
    species.write_text("((A,X),B);")
    genes.write_text(text or "((a_A,x1_X),(b_B,x2_X));\n(a_A,b_B);\n")
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "nwkit",
            "mul-reconcile",
            "-i",
            str(genes),
            "--species-tree",
            str(species),
            "--species-regex",
            r".*_([^_]+)$",
            "--h1",
            "X",
            "--h2",
            "B",
            "-o",
            str(tmp_path / "scores.tsv"),
            *extra,
        ],
        capture_output=True,
        text=True,
    )


def test_cli_staged_outputs_and_parallel_consistency(tmp_path):
    args = [
        "--report",
        str(tmp_path / "detail.tsv"),
        "--model-out",
        str(tmp_path / "model.json"),
        "--check-out",
        str(tmp_path / "check.tsv"),
        "--tree-out",
        str(tmp_path / "best.nwk"),
    ]
    result = invoke(tmp_path, args)
    assert result.returncode == 0, result.stderr
    initial = (tmp_path / "scores.tsv").read_bytes()
    metadata = json.loads((tmp_path / "model.json").read_text())
    details = pd.read_csv(tmp_path / "detail.tsv", sep="\t")
    scores = pd.read_csv(tmp_path / "scores.tsv", sep="\t")
    assert details["total.score"].sum() == scores.iloc[0]["score"]
    assert metadata["gene_filtering"].startswith("none")
    assert json.loads(details.iloc[0]["node.maps"])
    assert details.iloc[0]["maps.format"] == "grampa-annotated-newick"
    assert "x1_X[X+-0]" in details.iloc[0]["maps"]
    parallel = invoke(tmp_path, [*args, "--cpus", "2"])
    assert parallel.returncode == 0, parallel.stderr
    assert (tmp_path / "scores.tsv").read_bytes() == initial


def test_failure_preserves_all_existing_outputs(tmp_path):
    score = tmp_path / "scores.tsv"
    detail = tmp_path / "detail.tsv"
    score.write_text("old score")
    detail.write_text("old detail")
    result = invoke(tmp_path, ["--report", str(detail), "--max-state-pairs", "1"])
    assert result.returncode != 0
    assert score.read_text() == "old score"
    assert detail.read_text() == "old detail"
    result = invoke(tmp_path, ["--report", str(tmp_path / "genes.nwk")])
    assert result.returncode != 0
    assert "overwrite" in result.stderr.lower()


def test_mapping_cap_failure_preserves_outputs_without_truncation(tmp_path):
    score = tmp_path / "scores.tsv"
    detail = tmp_path / "detail.tsv"
    score.write_text("old score")
    detail.write_text("old detail")
    result = invoke(
        tmp_path, ["--report", str(detail), "--max-maps", "1"], text="(x1_X,x2_X);\n"
    )
    assert result.returncode != 0
    assert "no maps were truncated" in result.stderr
    assert score.read_text() == "old score"
    assert detail.read_text() == "old detail"


def test_primary_stdout_and_stdin_ownership(tmp_path):
    result = invoke(tmp_path, ["--outfile", "-"])
    assert result.returncode == 0, result.stderr
    assert result.stdout.startswith("mul.tree\th1.node\th2.node\tscore\t")
    assert "MUL reconciliation" not in result.stdout
    result = invoke(tmp_path, ["--report", "-"])
    assert result.returncode != 0
    result = invoke(tmp_path, ["--infile", "-", "--species-tree", "-"])
    assert result.returncode != 0
    assert "STDIN" in result.stderr


def test_empty_gene_collection_is_not_a_zero_score_analysis(tmp_path):
    result = invoke(tmp_path, text="\n")
    assert result.returncode != 0
