#!/usr/bin/env python3
"""Run deterministic output checks and the predefined input-case corpus.

The script deliberately uses different computational routes for the checks:
NWKIT's public CLI is compared with DendroPy for rooted RF distances and with
direct clade counting for tree-collection summaries. Mk marginal
probabilities are compared with exhaustive enumeration of all internal-state
assignments on four-tip trees.
"""

from __future__ import annotations

import argparse
import csv
import io
import itertools
import math
from pathlib import Path
import random
import subprocess
import tempfile

import dendropy
from dendropy.calculate import treecompare
import numpy as np
import pandas as pd
from ete4 import Tree

from nwkit.asr import compute_mk_marginals


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS_DIR = PROJECT_ROOT / "paper" / "results"
NWKIT = Path(__file__).resolve().parents[2] / ".not-used"


def run_cli(arguments: list[str], stdin: str | None = None) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["nwkit", *arguments],
        input=stdin,
        text=True,
        capture_output=True,
        check=True,
    )


def random_rooted_tree(labels: list[str], rng: random.Random) -> str:
    clusters = list(labels)
    while len(clusters) > 2:
        first_index, second_index = sorted(rng.sample(range(len(clusters)), 2), reverse=True)
        first = clusters.pop(first_index)
        second = clusters.pop(second_index)
        first_length = rng.uniform(0.05, 2.0)
        second_length = rng.uniform(0.05, 2.0)
        clusters.append(f"({first}:{first_length:.8f},{second}:{second_length:.8f})")
    left_length = rng.uniform(0.05, 2.0)
    right_length = rng.uniform(0.05, 2.0)
    return f"({clusters[0]}:{left_length:.8f},{clusters[1]}:{right_length:.8f});"


def dendropy_tree(newick: str, namespace: dendropy.TaxonNamespace) -> dendropy.Tree:
    return dendropy.Tree.get(
        data=newick,
        schema="newick",
        rooting="force-rooted",
        taxon_namespace=namespace,
        preserve_underscores=True,
    )


def dendropy_clades(tree: dendropy.Tree) -> set[frozenset[str]]:
    all_labels = {taxon.label for taxon in tree.taxon_namespace}
    clades: set[frozenset[str]] = set()
    for node in tree.postorder_node_iter():
        if node.is_leaf() or node.parent_node is None:
            continue
        labels = frozenset(leaf.taxon.label for leaf in node.leaf_iter())
        if 1 < len(labels) < len(all_labels):
            clades.add(labels)
    return clades


def ete_clades(tree: Tree) -> set[frozenset[str]]:
    all_labels = set(tree.leaf_names())
    clades: set[frozenset[str]] = set()
    for node in tree.traverse():
        if node.is_root or node.is_leaf:
            continue
        labels = frozenset(node.leaf_names())
        if 1 < len(labels) < len(all_labels):
            clades.add(labels)
    return clades


def parse_rf_table(text: str) -> int:
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    return int(rows[0]["rf_dist"])


def check_rooted_rf(rng: random.Random, repetitions: int = 50) -> dict[str, object]:
    labels = [f"T{index:02d}" for index in range(16)]
    namespace = dendropy.TaxonNamespace(labels)
    passed = 0
    differences: list[int] = []
    with tempfile.TemporaryDirectory(prefix="nwkit-rf-") as temporary_directory:
        directory = Path(temporary_directory)
        for index in range(repetitions):
            first = random_rooted_tree(labels, rng)
            second = random_rooted_tree(labels, rng)
            first_path = directory / f"first-{index}.nwk"
            second_path = directory / f"second-{index}.nwk"
            first_path.write_text(first)
            second_path.write_text(second)
            observed = parse_rf_table(
                run_cli(["dist", "-i", str(first_path), "-i2", str(second_path)]).stdout
            )
            first_tree = dendropy_tree(first, namespace)
            second_tree = dendropy_tree(second, namespace)
            first_tree.encode_bipartitions()
            second_tree.encode_bipartitions()
            expected = int(treecompare.symmetric_difference(first_tree, second_tree))
            difference = observed - expected
            differences.append(difference)
            passed += int(difference == 0)
    return {
        "check": "rooted_RF_vs_DendroPy",
        "n": repetitions,
        "passed": passed,
        "agreement": passed / repetitions,
        "max_abs_error": max(abs(value) for value in differences),
        "notes": "Rooted symmetric difference on random 16-tip binary trees",
    }


def check_clade_frequencies(rng: random.Random, repetitions: int = 20) -> dict[str, object]:
    labels = [f"T{index:02d}" for index in range(12)]
    namespace = dendropy.TaxonNamespace(labels)
    passed = 0
    maximum_error = 0.0
    with tempfile.TemporaryDirectory(prefix="nwkit-cladefreq-") as temporary_directory:
        directory = Path(temporary_directory)
        for replicate in range(repetitions):
            tree_strings = [random_rooted_tree(labels, rng) for _ in range(25)]
            tree_path = directory / f"trees-{replicate}.nwk"
            tree_path.write_text("\n".join(tree_strings) + "\n")
            observed_table = pd.read_csv(
                io.StringIO(
                    run_cli([
                        "cladefreq",
                        "-i",
                        str(tree_path),
                        "--support-scale",
                        "proportion",
                    ]).stdout
                ),
                sep="\t",
            )
            observed = {
                frozenset(str(row.leaf_set).split(",")): float(row.frequency)
                for row in observed_table.itertuples()
            }
            counts: dict[frozenset[str], int] = {}
            for tree_string in tree_strings:
                for clade in dendropy_clades(dendropy_tree(tree_string, namespace)):
                    counts[clade] = counts.get(clade, 0) + 1
            expected = {clade: count / len(tree_strings) for clade, count in counts.items()}
            keys = set(observed) | set(expected)
            error = max(abs(observed.get(key, 0.0) - expected.get(key, 0.0)) for key in keys)
            maximum_error = max(maximum_error, error)
            passed += int(error < 1e-12)
    return {
        "check": "clade_frequency_vs_direct_count",
        "n": repetitions,
        "passed": passed,
        "agreement": passed / repetitions,
        "max_abs_error": maximum_error,
        "notes": "Exact descendant-tip-set counts in collections of 25 random 12-tip trees",
    }


def check_majority_consensus(rng: random.Random, repetitions: int = 20) -> dict[str, object]:
    labels = [f"T{index:02d}" for index in range(12)]
    passed = 0
    with tempfile.TemporaryDirectory(prefix="nwkit-consensus-") as temporary_directory:
        directory = Path(temporary_directory)
        for replicate in range(repetitions):
            # Resampling a smaller topology pool makes majority clades likely.
            pool = [random_rooted_tree(labels, rng) for _ in range(4)]
            tree_strings = [pool[rng.randrange(len(pool))] for _ in range(31)]
            tree_path = directory / f"trees-{replicate}.nwk"
            output_path = directory / f"consensus-{replicate}.nwk"
            tree_path.write_text("\n".join(tree_strings) + "\n")
            run_cli([
                "consensus",
                "-i",
                str(tree_path),
                "-o",
                str(output_path),
                "--method",
                "majority",
                "--branch-length",
                "none",
            ])
            observed = ete_clades(Tree(output_path.read_text(), parser=0))
            counts: dict[frozenset[str], int] = {}
            namespace = dendropy.TaxonNamespace(labels)
            for tree_string in tree_strings:
                for clade in dendropy_clades(dendropy_tree(tree_string, namespace)):
                    counts[clade] = counts.get(clade, 0) + 1
            expected = {
                clade for clade, count in counts.items()
                if count / len(tree_strings) > 0.5
            }
            passed += int(observed == expected)
    return {
        "check": "majority_consensus_vs_direct_count",
        "n": repetitions,
        "passed": passed,
        "agreement": passed / repetitions,
        "max_abs_error": 0 if passed == repetitions else 1,
        "notes": "Clades occurring in more than half of 31 rooted trees",
    }


def transition_probability(parent_state: int, child_state: int, branch_length: float, rate: float) -> float:
    decay = math.exp(-2.0 * rate * branch_length)
    if parent_state == child_state:
        return 0.5 + 0.5 * decay
    return 0.5 - 0.5 * decay


def brute_force_marginals(tree: Tree, tip_states: dict[str, int], rate: float) -> dict[Tree, np.ndarray]:
    internal_nodes = [node for node in tree.traverse() if not node.is_leaf]
    weights: list[tuple[tuple[int, ...], float]] = []
    for assignment in itertools.product((0, 1), repeat=len(internal_nodes)):
        state_by_node = dict(zip(internal_nodes, assignment))
        probability = 0.5
        for node in tree.traverse():
            if node.is_root:
                continue
            parent_state = state_by_node[node.up]
            child_state = tip_states[node.name] if node.is_leaf else state_by_node[node]
            probability *= transition_probability(parent_state, child_state, float(node.dist), rate)
        weights.append((assignment, probability))
    total = sum(weight for _, weight in weights)
    marginals = {node: np.zeros(2, dtype=float) for node in internal_nodes}
    for assignment, weight in weights:
        for node, state in zip(internal_nodes, assignment):
            marginals[node][state] += weight / total
    return marginals


def check_mk_marginals(rng: random.Random, repetitions: int = 100) -> dict[str, object]:
    maximum_error = 0.0
    passed = 0
    states = ["0", "1"]
    for _ in range(repetitions):
        lengths = [rng.uniform(0.05, 2.0) for _ in range(6)]
        tree = Tree(
            "((A:{0},B:{1}):{2},(C:{3},D:{4}):{5});".format(*[f"{value:.10f}" for value in lengths]),
            parser=1,
        )
        tip_state_indices = {name: rng.randrange(2) for name in ("A", "B", "C", "D")}
        observed = {name: states[index] for name, index in tip_state_indices.items()}
        likelihood = {
            name: np.asarray([1.0, 0.0]) if index == 0 else np.asarray([0.0, 1.0])
            for name, index in tip_state_indices.items()
        }
        rate = rng.uniform(0.02, 1.0)
        nwkit_marginals, _ = compute_mk_marginals(
            tree=tree,
            states=states,
            observed_state_by_leaf=observed,
            likelihood_by_leaf=likelihood,
            model="ER",
            rate=rate,
            root_prior_mode="equal",
        )
        expected = brute_force_marginals(tree, tip_state_indices, rate)
        error = max(
            float(np.max(np.abs(nwkit_marginals[node] - expected[node])))
            for node in expected
        )
        maximum_error = max(maximum_error, error)
        passed += int(error < 1e-10)
    return {
        "check": "Mk_marginals_vs_exhaustive_enumeration",
        "n": repetitions,
        "passed": passed,
        "agreement": passed / repetitions,
        "max_abs_error": maximum_error,
        "notes": "Two-state ER model on random four-tip trees; all internal-state assignments enumerated",
    }


def check_seed_reproducibility() -> dict[str, object]:
    tree = "(((A:1,B:2):3,(C:4,D:5):6):7,((E:8,F:9):10,(G:11,H:12):13):14);"
    first = run_cli(["shuffle", "--label", "yes", "--seed", "20260716"], stdin=tree).stdout
    second = run_cli(["shuffle", "--label", "yes", "--seed", "20260716"], stdin=tree).stdout
    third = run_cli(["shuffle", "--label", "yes", "--seed", "20260717"], stdin=tree).stdout
    passed = int((first == second) and (first != third))
    return {
        "check": "shuffle_seed_reproducibility",
        "n": 1,
        "passed": passed,
        "agreement": float(passed),
        "max_abs_error": 0 if passed else 1,
        "notes": "Exact output equality for repeated seed and inequality for a second seed",
    }


VALIDATION_CASES = [
    {
        "case": "valid_default",
        "newick": "((A:1,B:1)X:1,(C:1,D:1)Y:1);",
        "options": [],
        "expected": "",
    },
    {
        "case": "duplicate_tip_label",
        "newick": "((A:1,A:1)X:1,(C:1,D:1)Y:1);",
        "options": [],
        "expected": "duplicate_leaf_names",
    },
    {
        "case": "empty_tip_label",
        "newick": "((A:1,:1)X:1,(C:1,D:1)Y:1);",
        "options": [],
        "expected": "empty_leaf_names",
    },
    {
        "case": "negative_branch_length",
        "newick": "((A:-1,B:1)X:1,(C:1,D:1)Y:1);",
        "options": [],
        "expected": "negative_branch_length",
    },
    {
        "case": "singleton_node",
        "newick": "(((A:1)Z:0,B:1)X:1,(C:1,D:1)Y:1);",
        "options": ["--require-binary", "yes"],
        "expected": "not_binary",
    },
    {
        "case": "multifurcation",
        "newick": "((A:1,B:1,C:1)X:1,D:2);",
        "options": ["--require-binary", "yes"],
        "expected": "not_binary",
    },
    {
        "case": "missing_support",
        "newick": "((A:1,B:1):1,(C:1,D:1)80:1);",
        "options": ["--require-all-support", "yes"],
        "expected": "missing_support",
    },
    {
        "case": "non_ultrametric",
        "newick": "((A:1,B:2)X:1,(C:1,D:1)Y:1);",
        "options": ["--require-ultrametric", "yes"],
        "expected": "not_ultrametric",
    },
    {
        "case": "unrooted",
        "newick": "(A:1,B:1,C:1,D:1);",
        "options": ["--require-rooted", "yes"],
        "expected": "not_rooted",
    },
    {
        "case": "ambiguous_numeric_internal_labels",
        "newick": "((A:1,B:1)90:1,(C:1,D:1)80:1);",
        "options": ["--require-unambiguous-format", "yes"],
        "expected": "format_ambiguous",
    },
    {
        "case": "quoted_node_name",
        "newick": "(('A a':1,B:1)X:1,(C:1,D:1)Y:1);",
        "options": ["--require-unquoted-names", "yes"],
        "expected": "quoted_node_names",
    },
    {
        "case": "leaf_set_mismatch",
        "newick": "((A:1,B:1)X:1,(C:1,D:1)Y:1);\n((A:1,B:1)X:1,(C:1,E:1)Y:1);",
        "options": ["--require-same-leaf-set", "yes"],
        "expected": "leaf_set_mismatch",
    },
]


def run_validation_corpus() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    with tempfile.TemporaryDirectory(prefix="nwkit-validation-") as temporary_directory:
        directory = Path(temporary_directory)
        for index, case in enumerate(VALIDATION_CASES):
            input_path = directory / f"case-{index}.nwk"
            input_path.write_text(str(case["newick"]) + "\n")
            result = run_cli(["validate", "-i", str(input_path), *case["options"]])
            table = pd.read_csv(io.StringIO(result.stdout), sep="\t", keep_default_na=False)
            issues = {
                issue
                for value in table["issues"].tolist()
                for issue in str(value).split(",")
                if issue
            }
            expected = str(case["expected"])
            outcome_matched = bool(
                (not expected and set(table["status"]) == {"ok"}) or expected in issues
            )
            rows.append({
                "case": case["case"],
                "expected_issue": expected or "none",
                "reported_issues": ",".join(sorted(issues)) or "none",
                "outcome_matched": outcome_matched,
                "num_trees": len(table.index),
            })
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--seed", type=int, default=20260716)
    args = parser.parse_args()
    args.results_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    checks = [
        check_rooted_rf(rng),
        check_clade_frequencies(rng),
        check_majority_consensus(rng),
        check_mk_marginals(rng),
        check_seed_reproducibility(),
    ]
    pd.DataFrame(checks).to_csv(args.results_dir / "correctness_checks.tsv", sep="\t", index=False)
    input_cases = run_validation_corpus()
    input_cases.to_csv(args.results_dir / "input_case_checks.tsv", sep="\t", index=False)

    if any(check["passed"] != check["n"] for check in checks):
        raise SystemExit("At least one correctness check failed")
    if not input_cases["outcome_matched"].all():
        raise SystemExit("At least one predefined input case did not match its expected outcome")


if __name__ == "__main__":
    main()
